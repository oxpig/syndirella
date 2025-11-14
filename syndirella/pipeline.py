#!/usr/bin/env python3
"""
syndirella.run_pipeline.py

This script contains the main pipeline for syndirella.
"""

import datetime
import logging
import time
import traceback
from typing import List, Tuple, Dict, Optional
from dataclasses import dataclass

import pandas as pd
from rdkit.Chem import inchi
from rdkit.Chem.EnumerateStereoisomers import EnumerateStereoisomers, StereoEnumerationOptions

import syndirella.utils.check_inputs as check_inputs
import syndirella.utils.fairy as fairy
import syndirella.utils.structure_outputs as structure_outputs
from syndirella.route.Cobbler import Cobbler
from syndirella.utils.error import *
from syndirella.route.CobblersWorkshop import CobblersWorkshop
from syndirella.slipper.Slipper import Slipper
from syndirella.slipper.SlipperFitter import SlipperFitter
from syndirella.constants import DatabaseSearchTool, RetrosynthesisTool, DEFAULT_DATABASE_SEARCH_TOOL, DEFAULT_RETROSYNTHESIS_TOOL

logger = logging.getLogger(__name__)


@dataclass
class PipelineConfig:
    """Configuration class for pipeline settings to reduce variable definitions."""
    # Required settings
    csv_path: str
    output_dir: str
    template_dir: str
    hits_path: str
    batch_num: int
    atom_diff_min: int
    atom_diff_max: int
    scaffold_place_num: int
    retro_tool: RetrosynthesisTool
    db_search_tool: DatabaseSearchTool
    long_code_column: str
    
    # Optional settings with defaults
    manual_routes: bool = False
    only_scaffold_place: bool = False
    scaffold_place: bool = True
    elab_single_reactant: bool = False
    additional_columns: List[str] = None
    reference_db: str = None
    metadata_path: str = None
    assert_scaffold_intra_geom_flatness: bool = True
    
    def __post_init__(self):
        if self.additional_columns is None:
            self.additional_columns = ['compound_set']
    
    @classmethod
    def from_settings(cls, settings: Dict) -> 'PipelineConfig':
        """Create PipelineConfig from settings dictionary."""
        try:
            return cls(
                csv_path=settings['input'],
                output_dir=settings['output'],
                template_dir=settings['templates'],
                hits_path=settings['hits_path'],
                metadata_path=settings.get('metadata', None),
                batch_num=settings['batch_num'],
                atom_diff_min=settings['atom_diff_min'],
                atom_diff_max=settings['atom_diff_max'],
                scaffold_place_num=settings['scaffold_place_num'],
                retro_tool=RetrosynthesisTool.from_string(settings.get('retro_tool', DEFAULT_RETROSYNTHESIS_TOOL.value)),
                db_search_tool=DatabaseSearchTool.from_string(settings.get('db_search_tool', DEFAULT_DATABASE_SEARCH_TOOL.value)),
                long_code_column=settings['long_code_column'],
                manual_routes=settings.get('manual', False),
                only_scaffold_place=settings.get('only_scaffold_place', False),
                scaffold_place=not settings.get('no_scaffold_place', False),
                elab_single_reactant=settings.get('elab_single_reactant', False),
                reference_db=settings.get('reference_db', None),
                assert_scaffold_intra_geom_flatness=not settings.get('no_assert_scaffold_intra_geom_flatness', False)
            )
        except KeyError as e:
            logger.critical(f"Missing critical argument to run pipeline: {e}")
            raise


def assert_scaffold_placement(scaffold: str,
                              template_path: str,
                              hits_path: str,
                              hits_names: List[str],
                              output_dir: str,
                              scaffold_place_num: int,
                              assert_intra_geom_flatness: bool = True,
                              ) -> Dict[Chem.Mol, str]:
    """
    Assert that the scaffold can be placed for any stereoisomers. If not, raise an error.
    """
    scaffold_mol = Chem.MolFromSmiles(scaffold)
    if not scaffold_mol:
        logger.critical(f"Could not create a molecule from the smiles {scaffold}.")
        raise MolError(smiles=scaffold)
    # enumerate stereoisomers
    opts = StereoEnumerationOptions(unique=True)
    isomers = list(EnumerateStereoisomers(scaffold_mol, options=opts))
    slipper_fitter = SlipperFitter(template_path, hits_path, hits_names, output_dir)
    placements: Dict[Chem.Mol, str | None] = {}
    for i, isomer in enumerate(isomers):
        scaffold_name: str = f'scaffold-{chr(65 + i)}'
        can_be_placed: str | None = slipper_fitter.check_scaffold(scaffold=isomer,
                                                                  scaffold_name=scaffold_name,
                                                                  scaffold_place_num=scaffold_place_num,
                                                                  assert_intra_geom_flatness=assert_intra_geom_flatness)  # path to scaffold if successful
        placements[isomer] = can_be_placed  # absolute path to minimised.mol scaffold, checked to exist
    if not any(placements.values()):
        logger.critical(f"Scaffold {scaffold} could not be placed successfully.")
        raise ScaffoldPlacementError(smiles=scaffold)
    return placements


def elaborate_from_cobbler_workshops(cobbler_workshops: List[CobblersWorkshop],
                                     template_path: str,
                                     hits_path: str,
                                     hits: List[str],
                                     batch_num: int,
                                     csv_path: str,
                                     output_dir: str,
                                     scaffold_placements: Dict[Chem.Mol, str],
                                     additional_info=None):
    """
    Does elaboration once the cobbler workshops are created.
    """
    if additional_info is None:
        additional_info = []
    for workshop in cobbler_workshops:
        try:
            slipper = None
            final_library = workshop.get_final_library()
            if final_library is None:
                logger.warning(f"Could not get the final library for compound {workshop.scaffold}. Skipping...")
                continue
            slipper = Slipper(library=final_library, template=template_path, hits_path=hits_path, hits_names=hits,
                              batch_num=batch_num, atoms_ids_expansion=None, additional_info=additional_info,
                              scaffold_placements=scaffold_placements)
            slipper.get_products()
            slipper.place_products()
            try:
                slipper.write_products_to_hippo()  # only write at the end after placement, to get correct route_uuid
            finally:
                slipper.clean_up_placements()
        except Exception as e:
            tb = traceback.format_exc()
            logger.critical(f"Error elaborating compound {workshop.scaffold}. {tb}")
            structure_outputs.structure_pipeline_outputs(error=e,
                                                         csv_path=csv_path,
                                                         output_dir=output_dir,
                                                         workshop=workshop if workshop is not None else None,
                                                         slipper=slipper if slipper is not None else None)
            continue
        structure_outputs.structure_pipeline_outputs(csv_path=csv_path,
                                                     output_dir=output_dir,
                                                     workshop=workshop,
                                                     slipper=slipper,
                                                     error=None)


def start_elaboration(product: str,
                      template_path: str,
                      hits_path: str,
                      hits: List[str],
                      output_dir: str,
                      scaffold_place_num: int,
                      scaffold_place: bool,
                      assert_scaffold_intra_geom_flatness: bool = True,) -> Tuple[float, Dict[Chem.Mol, str | None] | Dict]:
    """
    Start the elaboration process for a compound.
    """
    start_time = time.time()
    scaffold_placements = {}
    
    if scaffold_place:
        scaffold_placements = assert_scaffold_placement(
            scaffold=product,
            template_path=template_path,
            hits_path=hits_path,
            hits_names=hits,
            output_dir=output_dir,
            scaffold_place_num=scaffold_place_num,
            assert_intra_geom_flatness=assert_scaffold_intra_geom_flatness,
        )
    
    return start_time, scaffold_placements


def elaborate_compound_with_manual_routes(product: str,
                                          reactants: List[Tuple[str, str]],
                                          reaction_names: List[str],
                                          num_steps: int,
                                          hits: List[str],
                                          template_path: str,
                                          hits_path: str,
                                          batch_num: int,
                                          output_dir: str,
                                          csv_path: str,
                                          atom_diff_min: int,
                                          atom_diff_max: int,
                                          scaffold_place_num: int,
                                          only_scaffold_place: bool,
                                          scaffold_place: bool,
                                          elab_single_reactant: bool,
                                          retro_tool: RetrosynthesisTool,
                                          db_search_tool: DatabaseSearchTool,
                                          reference_db: str,
                                          assert_scaffold_intra_geom_flatness: bool = True,
                                          additional_info=None):
    """
    Elaborate compound using manual routes.
    """
    start_time, scaffold_placements = start_elaboration(
        product=product, 
        template_path=template_path,
        hits_path=hits_path,
        hits=hits, 
        output_dir=output_dir,
        scaffold_place_num=scaffold_place_num,
        scaffold_place=scaffold_place,
        assert_scaffold_intra_geom_flatness=assert_scaffold_intra_geom_flatness,
    )
    
    if not only_scaffold_place:

        # Create manual route workshop
        workshop = CobblersWorkshop(
            scaffold=product,
            reactants=reactants,
            reaction_names=reaction_names,
            num_steps=num_steps,
            output_dir=output_dir,
            atom_diff_min=atom_diff_min,
            atom_diff_max=atom_diff_max,
            elab_single_reactant=elab_single_reactant,
            db_search_tool=db_search_tool,
            retro_tool=retro_tool,
            id=fairy.generate_inchi_ID(product, isomeric=False),
            filter=False,
            reference_db=reference_db
        )
        
        elaborate_from_cobbler_workshops(
            cobbler_workshops=[workshop],
            template_path=template_path,
            hits_path=hits_path,
            hits=hits,
            batch_num=batch_num,
            additional_info=additional_info,
            csv_path=csv_path,
            output_dir=output_dir,
            scaffold_placements=scaffold_placements
        )
        
        end_time = time.time()
        elapsed_time = end_time - start_time
        logger.info(
            f"Finished Syndirella 👑pipeline for compound {product} | {inchi.MolToInchiKey(Chem.MolFromSmiles(product))} "
            f"after {datetime.timedelta(seconds=elapsed_time)}")
        logger.info("")


def elaborate_compound_full_auto(product: str,
                                 hits: List[str],
                                 template_path: str,
                                 hits_path: str,
                                 batch_num: int,
                                 output_dir: str,
                                 csv_path: str,
                                 atom_diff_min: int,
                                 atom_diff_max: int,
                                 scaffold_place_num: int,
                                 only_scaffold_place: bool,
                                 scaffold_place: bool,
                                 elab_single_reactant: bool,
                                 retro_tool: RetrosynthesisTool,
                                 db_search_tool: DatabaseSearchTool,
                                 reference_db: str | None = None,
                                 assert_scaffold_intra_geom_flatness: bool = True,
                                 additional_info=None):
    """
    Elaborate compound using full automatic retrosynthesis.
    """
    start_time, scaffold_placements = start_elaboration(
        product=product, 
        template_path=template_path,
        hits_path=hits_path,
        hits=hits, 
        output_dir=output_dir,
        scaffold_place_num=scaffold_place_num,
        scaffold_place=scaffold_place,
        assert_scaffold_intra_geom_flatness=assert_scaffold_intra_geom_flatness,
    )
    
    if not only_scaffold_place:
        cobbler = Cobbler(
            scaffold_compound=product,
            output_dir=output_dir, 
            atom_diff_min=atom_diff_min,
            atom_diff_max=atom_diff_max,
            elab_single_reactant=elab_single_reactant,
            retro_tool=retro_tool,
            db_search_tool=db_search_tool,
            reference_db=reference_db
        )
        cobbler_workshops: List[CobblersWorkshop] = cobbler.get_routes()
        elaborate_from_cobbler_workshops(
            cobbler_workshops=cobbler_workshops, 
            template_path=template_path,
            hits_path=hits_path, 
            hits=hits, 
            batch_num=batch_num,
            additional_info=additional_info, 
            csv_path=csv_path, 
            output_dir=output_dir,
            scaffold_placements=scaffold_placements
        )
        end_time = time.time()
        elapsed_time = end_time - start_time
        logger.info(
            f"Finished Syndirella 👑pipeline for compound {product} | {inchi.MolToInchiKey(Chem.MolFromSmiles(product))} "
            f"after {datetime.timedelta(seconds=elapsed_time)}")
        logger.info("")


def process_row(row: pd.Series, config: PipelineConfig):
    """
    Process a single row from the input CSV.
    """
    additional_info: dict = check_inputs.format_additional_info(row, config.additional_columns)
    template_path: str = check_inputs.get_template_path(
        template_dir=config.template_dir,
        template=row['template'],
        metadata_path=config.metadata_path
    )

    hits: List[str] = check_inputs.get_exact_hit_names(
        row=row, 
        metadata_path=config.metadata_path,
        hits_path=config.hits_path,
    )

    try:
        if config.manual_routes:
            reactants, reaction_names, num_steps = check_inputs.format_manual_route(row)
            elaborate_compound_with_manual_routes(
                product=row['smiles'],
                reactants=reactants,
                reaction_names=reaction_names,
                num_steps=num_steps,
                hits=hits,
                template_path=template_path,
                hits_path=config.hits_path,
                batch_num=config.batch_num,
                output_dir=config.output_dir,
                additional_info=additional_info,
                csv_path=config.csv_path,
                atom_diff_min=config.atom_diff_min,
                atom_diff_max=config.atom_diff_max,
                scaffold_place_num=config.scaffold_place_num,
                only_scaffold_place=config.only_scaffold_place,
                scaffold_place=config.scaffold_place,
                elab_single_reactant=config.elab_single_reactant,
                retro_tool=config.retro_tool,
                db_search_tool=config.db_search_tool,
                reference_db=config.reference_db,
                assert_scaffold_intra_geom_flatness=config.assert_scaffold_intra_geom_flatness,
            )
        else:
            elaborate_compound_full_auto(
                product=row['smiles'],
                hits=hits,
                template_path=template_path,
                hits_path=config.hits_path,
                batch_num=config.batch_num,
                output_dir=config.output_dir,
                additional_info=additional_info,
                csv_path=config.csv_path,
                atom_diff_min=config.atom_diff_min,
                atom_diff_max=config.atom_diff_max,
                scaffold_place_num=config.scaffold_place_num,
                only_scaffold_place=config.only_scaffold_place,
                scaffold_place=config.scaffold_place,
                elab_single_reactant=config.elab_single_reactant,
                retro_tool=config.retro_tool,
                db_search_tool=config.db_search_tool,
                reference_db=config.reference_db,
                assert_scaffold_intra_geom_flatness=config.assert_scaffold_intra_geom_flatness,
            )
    except Exception as e:
        tb = traceback.format_exc()
        logger.critical(f"Error elaborating compound {row['smiles']}. {tb}")
        structure_outputs.structure_pipeline_outputs(
            error=e,
            csv_path=config.csv_path,
            output_dir=config.output_dir,
            smiles=row['smiles'],
            template_path=template_path,
            hits=hits,
            additional_info=additional_info
        )


def run_pipeline(settings: Dict):
    """
    Run the whole syndirella pipeline! 👑
    """
    # Create configuration object
    config = PipelineConfig.from_settings(settings)
    
    # Log pipeline configuration
    if not config.only_scaffold_place:
        logger.info(f"Running the pipeline with {'manual' if config.manual_routes else 'full auto'} routes.")
        logger.info(f"Database search tool: {config.db_search_tool}")
        logger.info(f"Retrosynthesis tool: {config.retro_tool}")
    
    if config.only_scaffold_place:
        logger.info("Only placing scaffolds!")
    
    if not config.scaffold_place:
        logger.warning("Skipping initial scaffold placement! Immediately starting elaboration process.")
    
    if config.elab_single_reactant:
        logger.info("'--elab_single_reactant' set. Only elaborating a single reactant per elaboration series.")

    # Validate inputs
    check_inputs.check_pipeline_inputs(
        csv_path=config.csv_path,
        template_dir=config.template_dir,
        hits_path=config.hits_path,
        metadata_path=config.metadata_path,
        additional_columns=config.additional_columns,
        manual_routes=config.manual_routes,
        long_code_column=config.long_code_column
    )

    # Load data and process each row
    df = pd.read_csv(config.csv_path)
    for index, row in df.iterrows():
        process_row(row=row, config=config)

    logger.info("Pipeline complete.")

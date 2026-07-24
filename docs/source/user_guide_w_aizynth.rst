==========
User Guide
==========

Installation and setup
======================

1. Install Syndirella:
----------------------

On Mac OS and Linux you can install from PyPI using Conda.

.. code-block:: bash

    conda create -n syndirella python=3.10
    conda activate syndirella
    pip install syndirella
    pip install aizynthfinder
    python -c 'import pyrosetta_installer; pyrosetta_installer.install_pyrosetta()'  # run once, for Fragmenstein

.. note::

    PyRosetta is available for academic and non-commercial use (see `PyRosetta License <https://github.com/RosettaCommons/rosetta/blob/main/LICENSE.PyRosetta.md>`_).

**Troubleshooting:**

**cgrtools fails to install:**

.. code-block:: bash

    pip install "cython<3.2"
    conda install -c conda-forge c-compiler cxx-compiler
    pip install --no-build-isolation cgrtools

.. attention::

    Installation and usage have not been tested on Windows OS.

2. Setup AiZynthFinder:
-------------------------------------

Syndirella can use AiZynthFinder for retrosynthesis functionality (see Retrosynthesis for more information).

**Automatic Setup (Recommended):**
Syndirella will automatically handle AiZynthFinder setup when you first use it. Simply run:

.. code-block:: bash

    syndirella setup-aizynth

.. attention::

    ⚠️ This will download large model files (~750MB total) if run the first time (will take ~5 min). These are required to run AiZynthFinder. This will automatically download the required data to `[syndirella_package_path]/aizynth` directory and create the configuration file automatically. 

**Manual Setup (Alternative):**
If you prefer manual setup:

.. code-block:: bash

    cd [syndirella_package_path]/aizynth
    download_public_data .
    # Update config.yml (if you prefer)
    export AIZYNTH_CONFIG_FILE="path/to/syndirella/aizynth/config.yml"

Command Line Interface
======================

Syndirella provides a command-line interface with multiple subcommands. Get help with `-h` or `--help`.

**Available Commands:**

- **setup-aizynth**: Setup AiZynthFinder data and configuration
- **run**: Run the main Syndirella pipeline  
- **add-reaction**: Add a new reaction to the library

**Main Help Output:**

.. code-block:: bash

    usage: syndirella [-h] {setup-aizynth,run,add-reaction} ...

    Run the Syndirella pipeline with specified configurations.

    Available commands:
      setup-aizynth    Setup AiZynthFinder data and configuration
      run              Run the main Syndirella pipeline  
      add-reaction     Add a new reaction to the library

    Syndirella is installed at [path_to_installation]

**Run Command Help:**

.. code-block:: bash

    usage: syndirella run [-h] -i INPUT -o OUTPUT [-t TEMPLATES] [--hits_path HITS_PATH]
                      [--products PRODUCTS] [--batch_num BATCH_NUM] [--manual]
                      [--only_scaffold_place] [--scaffold_place_num SCAFFOLD_PLACE_NUM]
                      [--retro_tool {manifold,aizynthfinder}]
                      [--db_search_tool {manifold,arthor,hippo}] [--profile]
                      [--atom_diff_min ATOM_DIFF_MIN] [--atom_diff_max ATOM_DIFF_MAX]
                      [--just_retro] [--no_scaffold_place] [--elab_single_reactant]
                      [--reference_db REFERENCE_DB] [--no_assert_scaffold_intra_geom_flatness]
                      [--substructure_check_products]

    Run the full Syndirella pipeline with specified input files and parameters.

    options:
    -h, --help            show this help message and exit
    -i INPUT, --input INPUT
                            Input .csv file path for the pipeline.
    -o OUTPUT, --output OUTPUT
                            Output directory for the pipeline results.
    -t TEMPLATES, --templates TEMPLATES
                            Absolute path to a directory containing the template(s).
    --hits_path HITS_PATH
                            Absolute path to hits_path for placements (.sdf or .mol).
    --products PRODUCTS   Absolute path to products for placements.
    --batch_num BATCH_NUM
                            Batch number for processing.
    --manual              Use manual routes for processing.
    --only_scaffold_place
                            Only place scaffolds. Do not continue to elaborate.
    --scaffold_place_num SCAFFOLD_PLACE_NUM
                            Number of times to attempt scaffold placement.
    --retro_tool {manifold,aizynthfinder}
                            Retrosynthesis tool to use.
    --db_search_tool {manifold,arthor,hippo}
                            Database search tool to use.
    --profile             Run the pipeline with profiling.
    --atom_diff_min ATOM_DIFF_MIN
                            Minimum atom difference between elaborations and scaffold to keep.
    --atom_diff_max ATOM_DIFF_MAX
                            Maximum atom difference between elaborations and scaffold to keep.
    --just_retro          Only run retrosynthesis querying of scaffolds.
    --no_scaffold_place   Do not place scaffolds initially before elaborating.
    --elab_single_reactant
                            Only elaborate one reactant per elaboration series.
    --reference_db REFERENCE_DB
                            Path to reference HIPPO database file for superstructure search, must set
                            --db_search_tool to 'hippo'.
    --no_assert_scaffold_intra_geom_flatness
                            Don't check scaffold for intra geometry or flatness.
    --substructure_check_products
                            Only keep products that are superstructures of their scaffolds.

**Add Reaction Command Help:**

.. code-block:: bash

    usage: syndirella add-reaction [-h] --name NAME --smirks SMIRKS [--find_parent] [--fp_type {maccs_rxn_fp,morgan_rxn_fp}] [--threshold THRESHOLD] [--similarity_metric {tanimoto,dice,cosine}]

    Add a new reaction SMIRKS to the reaction library with optional parent finding.

    options:
      -h, --help            show this help message and exit
      --name NAME           Name of the new reaction.
      --smirks SMIRKS       SMIRKS string for the reaction.
      --find_parent         If True, treat as a child reaction and find parent based on similarity. (default: False)
      --fp_type {maccs_rxn_fp,morgan_rxn_fp}
                            Fingerprint type for similarity calculation. (default: maccs_rxn_fp)
      --threshold THRESHOLD
                            Similarity threshold for finding parent reaction. (default: 0.2)
      --similarity_metric {tanimoto,dice,cosine}
                            Similarity metric for finding parent reaction. (default: tanimoto)

**Setup AiZynthFinder Command Help:**

.. code-block:: bash

    usage: syndirella setup-aizynth [-h]

    Automatically download AiZynthFinder data and create configuration file.

    options:
      -h, --help  show this help message and exit

Default Tools
=============

Syndirella uses the following default tools:

**Default Retrosynthesis Tool**: ``aizynthfinder``
    - Alternative: ``manifold``
    - Set with: ``--retro_tool {aizynthfinder,manifold}``

**Default Database Search Tool**: ``arthor``  
    - Alternative: ``manifold``
    - Set with: ``--db_search_tool {arthor,manifold, or hippo if installed}``

.. note::

    If using Manifold for retrosynthesis or database search, you must set up your Manifold API credentials.
    See the :doc:`Retrosynthesis <retrosynthesis>` section for detailed setup instructions.

Basic Usage
===========

Elaborate a set of scaffolds using these steps:

1. Setup fragments and protein templates
----------------------------------------

Download the fragment hits from Fragalysis. In the download folder the important files are:

::

    target_name_combined.sdf # fragment poses with long code names
    /aligned_files/fragment_name/fragment_name_apo-desolv.pdb # apo pdb used for placement

.. attention::

    **IMPORTANT**: The template string in your CSV must **exactly match** the PDB filename (without extension). The hit names in your CSV must **exactly match** the molecule names in the SDF file.

2. Create input csv
-------------------

**Critical Requirements for Exact Matching:**

- **Template names**: Must exactly match the PDB filename (without .pdb extension)
- **Hit names**: Must exactly match the molecule names in the SDF file
- **No metadata file needed**: Direct matching eliminates the need for metadata.csv

Syndirella can be run either in *automatic* or *manual* mode.

**Automatic**:
    Scaffolds can be elaborated by routes automatically proposed by Manifold.
    An example template is at ``examples/run_syndirella_example/syndirella_input_example_automatic.csv``.

Required headers:

    ``smiles``:
        smiles string of :term:`scaffold`.
    ``hit1``:
        string that **exactly matches** the molecule name in the SDF file for 1 fragment inspiring hit.
    ``template``:
        string that **exactly matches** the PDB filename (without extension) to use for :term:`placement`.
    ``compound_set``:
        string or int identifier.

Not required headers:

    ``hitX``:
        string of short code of additional fragment inspiring hit.

.. note::

    Any number of fragment inspirations can be used. You just need to specify in a seperate header. Ex.
    ``hit1, hit2, hit3, hit4, hit5``.


**Manual**:
    You can set the exact route to elaborate the scaffold with the reaction names, exact reactants, and number of steps in the route.
    An example template is at ``examples/run_syndirella_example/syndirella_input_example_manual.csv``.

Required headers:

    ``smiles``:
        smiles string of scaffold.
    ``hit1``:
        string that **exactly matches** the molecule name in the SDF file for 1 fragment inspiring hit.
    ``template``:
        string that **exactly matches** the PDB filename (without extension) to use for :term:`placement`.
    ``compound_set``:
        string or int identifier.
    ``reaction_name_step1``:
        string of reaction name.
    ``reactant_step1``:
        smiles string of reactant.

Not required headers:

    ``reactant2_step1``:
        smiles string of second reactant in reaction step 1.
    ``product_stepX``:
        smiles string of product of step X. Only required for internal or first step to specify reactant for next step. Not required
        if step is final step of route (as the scaffold is the final product).
    ``reaction_name_stepX``:
        string of reaction name of step X.
    ``reactant_stepX``:
        smiles string of reactant that is *not* a product of previous step.
    ``hitX``:
        string of short code of additional fragment inspiring hit. Any number of hits can be used.

3. Run!
-------

**Important: Path Requirements**

All file paths must be **absolute paths** (not relative paths). This includes:
- Input CSV file path
- Output directory path  
- Template directory path
- Hits path (SDF/MOL file) - optional
- Metadata CSV file path - optional
- Products path - optional

Run pipeline in *automatic* mode:

.. code-block:: bash

    syndirella run --input [path_to_automatic.csv] --output [path_to_output_dir] --templates [path_to_templates_dir]
    --hits_path [path_to_fragments.sdf]


Run pipeline in *manual* mode:
    Add ``--manual`` flag.

4. Outputs
----------

**Output directory structure:**

🔑🔑🔑: Inchi key of flat scaffold (removed stereochemistry). Example: ``ZJENMQHSGLZNHL-UHFFFAOYSA-N``

.. code-block::

    output_dir
    ├── 🔑🔑🔑-scaffold-check # scaffold check directory per scaffold
    │   └── scaffold-check
    │       ├── scaffold-check.holo_minimised.pdb
    │       ├── scaffold-check.minimised.json
    │       └── scaffold-check.minimised.mol
    ├── 🔑🔑🔑 # directory per scaffold
    │   ├── extra
    │   │   ├── 🔑🔑🔑_[route_uuid]_[rxn_name]_r[reactant_num]_[step_num]of[total_steps].pkl.gz # reactants for step
    │   │   └── continued for all steps...
    │   ├── output
    │   │   ├── fstein_input.pkl.gz   # Fragmenstein placement input 
    │   │   ├── fstein_output.pkl.gz  # Fragmenstein placement output 
    │   │   ├── 🔑🔑🔑_[route_uuid]_[num]-[stereoisomer]
    │   │   │   ├── 🔑🔑🔑_[route_uuid]_[num]-[stereoisomer].mol
    │   │   │   ├── 🔑🔑🔑_[route_uuid]_[num]-[stereoisomer].json # energy values
    │   │   └── continued for all products...
    │   ├── 🔑🔑🔑_[route_uuid]_structured_output.pkl.gz # KEY OUTPUT FILE - full routes and placements
    │   ├── 🔑🔑🔑_[route_uuid]_[rxn_name]_products_[last_step]of[total_steps].pkl.gz & .csv # final products
    │   ├── 🔑🔑🔑_[route_uuid]_[rxn_name]_products_[last_step]of[total_steps]_placements.pkl.gz & .csv # merged placements with products info
    │   └── 🔑🔑🔑_[route_uuid]_fragmenstein_placements.pkl.gz & .csv # fragmenstein output
    ├── continued for all scaffolds...
    └── [input_csv]_output_YYYYMMDD_HHMM.csv # summary stats of all scaffolds

**Important output files:**

**[input_csv]_output_YYYYMMDD_HHMM.csv:**
    Summary stats of all scaffolds. Most columns are self-explanatory. The following columns might need clarification:

    ``total_num_products_enumstereo``:
        Total number of products enumerated with stereochemistry in the final step. This is counting the number of unique
        products with stereochemistry, so if a product with same stereochemistry is generated multiple times via different routes
        it will only be counted once.

    ``total_num_unique_products``:
        Total number of unique products without stereochemistry in the final step. If a product is generated multiple times
        by different routes it will only be counted once.

**🔑🔑🔑_[route_uuid]_[rxn_name]_products_[last_step]of[total_steps]_placements.pkl.gz & .csv:**
    Merged placements with products info.

**🔑🔑🔑_[route_uuid]_structured_output.pkl.gz:**
    **⭐ KEY OUTPUT FILE ⭐** - Contains complete synthesis routes and placement information. This is the primary file to read for detailed results including:
    
    - Full synthesis routes with reaction names, reactants, and products for each step
    - Placement information with energy values (ΔΔG, ΔG_bound, ΔG_unbound)
    - Structural quality metrics (comRMSD, intra-geometry checks)
    - Product stereochemistry and atom differences
    - Success flags and error information
    - Paths to molecular structure files
    
    This file contains all the information needed to reproduce and analyze the elaborations.

.. note::

    Placements of products are labeled succesful if:
        1. ΔΔG < 0.
        2. comRMSD < 2.0 Å.
        3. Pose of product passes `PoseBusters <https://github.com/maabuu/posebusters>`_ intrageometry checks:
            - Bond lengths: The bond lengths in the input molecule are within 0.75 of the lower and 1.25 of the upper bounds determined by distance geometry.
            - Bond angles: The angles in the input molecule are within 0.75 of the lower and 1.25 of the upper bounds determined by distance geometry.
            - Planar aromatic rings: All atoms in aromatic rings with 5 or 6 members are within 0.25 Å of the closest shared plane.
            - Planar double bonds: The two carbons of aliphatic carbon–carbon double bonds and their four neighbours are within 0.25 Å of the closest shared plane.
            - Internal steric clash: The interatomic distance between pairs of non-covalently bound atoms is above 60% of the lower bound distance apart determined by distance geometry.


Usage Option: Only Place Scaffolds (or Specifically Don't Place)
===============================================================

You can run Syndirella to only place scaffolds. It will not perform the full elaboration procedure.
A Fragmenstein placements CSV (``{inchi}-scaffold-check_fragmenstein_placements.csv``) is written in each scaffold-check directory.

.. code-block:: bash

    syndirella run --input [path_to_automatic.csv] --output [path_to_output_dir] --templates [path_to_templates_dir]
    --hits_path [path_to_fragments.sdf] --only_scaffold_place

You can also specify to not place the scaffold (most likely you confirmed placement using another method).

.. code-block:: bash

    syndirella run --input [path_to_automatic.csv] --output [path_to_output_dir] --templates [path_to_templates_dir]
    --hits_path [path_to_fragments.sdf] --no_scaffold_place


Usage Option: Only Get Retrosynthesis Routes of Scaffolds
========================================================

You can run Syndirella to find the Top 5 retrosynthesis routes of the scaffolds. It will identify the routes that contains
all reactions you have encoded in the RXN_SMIRKS_CONSTANTS.json file (a CAR route) and routes that don't contain those
reactions (non-CAR route).

.. code-block:: bash

    syndirella run --input [path_to_automatic.csv] --output [path_to_output_dir] --just_retro

**Output file:**
    - ``justretroquery_[retro_tool]_[input_csv_name].csv``: CSV file with all route information (e.g., ``justretroquery_aizynthfinder_input.csv`` or ``justretroquery_manifold_input.csv``)

.. note::

    The CSV file can be opened directly in Excel or any spreadsheet application, or read using pandas: ``pd.read_csv('justretroquery_[retro_tool]_[input_csv_name].csv')``

Structure of the important columns (where X is 0-4 for the top 5 routes):

    ``routeX``:
        List of dictionaries for each step in the route. Each dictionary contains:
        - ``name``: Reaction name
        - ``reactantSmiles``: Tuple of reactant SMILES strings
        - ``productSmiles``: Expected product SMILES from retrosynthesis
        - ``smirks_validated``: Boolean indicating if applying Syndirella's SMIRKS to the reactants produces the expected product (by InChI-key comparison)
        - ``actual_product_smiles``: SMILES of the product actually produced by the SMIRKS (if validation was performed)
    
    ``routeX_names``:
        List of reaction names in the route.
    
    ``routeX_CAR``:
        Boolean indicating if all reactions in the route are in RXN_SMIRKS_CONSTANTS.json (Chemically Accessible Reactions).
    
    ``routeX_non_CAR``:
        List of reaction names that are not in RXN_SMIRKS_CONSTANTS.json. Or None if all reactions are in RXN_SMIRKS_CONSTANTS.json.
    
    ``routeX_smirks_validated``:
        Boolean indicating if all reactions in the route passed SMIRKS validation (True if all validated, False if any failed, None if validation couldn't be performed).
    
    ``routeX_num_validated``:
        Number of reactions in the route that passed SMIRKS validation.
    
    ``routeX_num_failed_validation``:
        Number of reactions in the route that failed SMIRKS validation.
    
    ``routeX_CAR_and_validated``:
        **This is the most important column to check!** Boolean that is True only if:
        - All reactions are in Syndirella's reaction library (CAR = True), AND
        - All reactions passed SMIRKS validation (smirks_validated = True)
        
        Routes with ``routeX_CAR_and_validated = True`` are fully compatible with Syndirella and have been validated to work correctly.

.. attention::

    **Look for routes where ``route0_CAR_and_validated = True`` (or route1, route2, etc.).** These are the routes that:
    
    1. Use only reactions from Syndirella's reaction library (CAR routes)
    2. Have been validated to produce the expected products when applying Syndirella's SMIRKS patterns
    
    These routes are the most reliable for use in the full Syndirella pipeline.

If there are `NaN` values for all route columns, it means that there are no routes found for the scaffold.

Usage Option: Only Elaborate One Reactant per Series
========================================================

.. attention::

    This functionality is only provided for single step reactions.

You can have Syndirella output elaboration series for one reactant at a time. For example, if the route is a single step
amidation, there will be two elaboration series output: (1) only elaborating reactant 1 and (2) only elaborating reactant 2.

.. note::

    Each series per reactant will be handled as seperate, so they will have their own unique route uuids. If an
    alternative route is found for the original route, the alternative route will produce two seperate series as well
    for each reactant elaboration.

.. code-block:: bash

    syndirella run --input [path_to_input.csv] --output [path_to_output_dir] --templates [path_to_templates_dir]
    --hits_path [path_to_fragments.sdf] --elab_single_reactant




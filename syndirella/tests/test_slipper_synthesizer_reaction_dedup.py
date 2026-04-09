#!/usr/bin/env python3
"""Targeted tests for reaction output deduplication in SlipperSynthesizer."""

from types import SimpleNamespace

import pandas as pd
from rdkit import Chem

from syndirella.slipper.slipper_synthesizer.SlipperSynthesizer import SlipperSynthesizer


class DummyReactionPattern:
    def __init__(self, products):
        self._products = products

    def RunReactants(self, reactants):
        return self._products


def _make_synthesizer(products) -> SlipperSynthesizer:
    reaction = SimpleNamespace(
        reaction_pattern=DummyReactionPattern(products),
        scaffold=Chem.MolFromSmiles("CC"),
        reaction_name="dummy_reaction",
    )
    library = SimpleNamespace(
        route_uuid="route-1",
        current_step=1,
        num_steps=2,
        atom_diff_min=-999,
        atom_diff_max=999,
        reaction=reaction,
        id="lib-1",
        elab_single_reactant_int=None,
    )
    return SlipperSynthesizer(library=library, output_dir="/tmp")


def _mock_rdkit_products():
    return [
        (Chem.MolFromSmiles("CCO"), Chem.MolFromSmiles("NCC")),
        (Chem.MolFromSmiles("CCO"),),
        (None,),
    ]


def test_apply_reaction_deduplicates_and_sets_multi_product_flag():
    synth = _make_synthesizer(_mock_rdkit_products())
    row = pd.Series(
        {
            "r1_mol": Chem.MolFromSmiles("C"),
            "r2_mol": Chem.MolFromSmiles("N"),
            "flag": ("existing_flag",),
        }
    )

    out = synth.apply_reaction(row)

    expected = [
        (Chem.MolToSmiles(Chem.MolFromSmiles("CCO"), isomericSmiles=False), None),
        (Chem.MolToSmiles(Chem.MolFromSmiles("NCC"), isomericSmiles=False), None),
    ]
    assert out["combined"] == expected
    assert out["flag"] == ["existing_flag", "one_of_multiple_products"]


def test_apply_reaction_single_deduplicates_and_preserves_existing_string_flag():
    synth = _make_synthesizer(_mock_rdkit_products())
    row = pd.Series(
        {
            "r1_mol": Chem.MolFromSmiles("C"),
            "flag": "selectivity_issue",
        }
    )

    out = synth.apply_reaction_single(row)

    expected = [
        (Chem.MolToSmiles(Chem.MolFromSmiles("CCO"), isomericSmiles=False), None),
        (Chem.MolToSmiles(Chem.MolFromSmiles("NCC"), isomericSmiles=False), None),
    ]
    assert out["combined"] == expected
    assert out["flag"] == ["selectivity_issue", "one_of_multiple_products"]

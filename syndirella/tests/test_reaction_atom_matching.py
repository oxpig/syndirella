"""
test_reaction_atom_matching.py

Unit tests for Reaction atom/attachment matching. These run offline (no
retrosynthesis or data downloads) by constructing a Reaction directly from a
scaffold and its reactants.
"""
import logging
import unittest

from rdkit import Chem

from syndirella.route.Reaction import Reaction
from syndirella.route.SMARTSHandler import SMARTSHandler


class TestReactionAtomMatching(unittest.TestCase):
    def setUp(self):
        logging.disable(logging.CRITICAL)
        self.smarts_handler = SMARTSHandler()

    def tearDown(self):
        logging.disable(logging.NOTSET)

    def _build_reaction(self, scaffold_smiles, reactant_smiles, reaction_name):
        scaffold = Chem.MolFromSmiles(scaffold_smiles)
        reactants = [Chem.MolFromSmiles(s) for s in reactant_smiles]
        return Reaction(scaffold=scaffold,
                        reactants=reactants,
                        reaction_name=reaction_name,
                        smarts_handler=self.smarts_handler,
                        route_uuid='test')

    def test_issue_107_amidation_piperidine_and_pyridine(self):
        """
        Regression test for https://github.com/oxpig/syndirella/issues/107.

        The scaffold contains both a pyridine ring and a piperidine ring. The MCS
        of the piperidine amine reactant also matches the aromatic pyridine ring,
        so blindly taking the first substructure match mapped the amine onto the
        wrong ring and pushed the attachment index onto a carbon instead of the
        ring nitrogen, causing a spurious SMARTSError.
        """
        rxn = self._build_reaction(
            scaffold_smiles='O=C(Cc1ccccn1)N1CCCCC1',
            reactant_smiles=['O=C(O)Cc1ccccn1', 'C1CCNCC1'],
            reaction_name='Amidation')

        matched = {Chem.MolToSmiles(mol): set(ids)
                   for _, (mol, ids, _) in rxn.matched_smarts_index_to_reactant.items()}

        # The piperidine's matched atom must be its ring nitrogen (index 3).
        self.assertIn('C1CCNCC1', matched)
        self.assertEqual(matched['C1CCNCC1'], {3})

    def test_amidation_simple_case_still_matches(self):
        """A plain amidation with unambiguous rings should be unaffected."""
        rxn = self._build_reaction(
            scaffold_smiles='CC(=O)NCc1ccccc1',
            reactant_smiles=['CC(=O)O', 'NCc1ccccc1'],
            reaction_name='Amidation')
        self.assertEqual(len(rxn.matched_smarts_index_to_reactant), 2)


if __name__ == '__main__':
    unittest.main()

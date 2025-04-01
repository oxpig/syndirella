# python
import logging
import unittest

from syndirella.fairy import generate_inchi_ID
from syndirella.route.CobblersWorkshop import CobblersWorkshop


class TestAdditionalRouteGeneration(unittest.TestCase):
    # TODO: Need to check when getting AiZynthFinder functionality
    def setUp(self):
        self.product = 'CC(=O)Nc1cncc(CC(=O)NC(C)=O)c1'
        self.reactants = [('CCOC(=O)Cc1cncc(N)c1', 'CC(=O)Cl')]
        self.reaction_names = ['Amide_Schotten-Baumann_with_amine']
        self.num_steps = 1
        self.output_dir = '/Users/kate_fieseler/PycharmProjects/syndirella/syndirella/tests/additional_route'
        self.filter = False
        self.id = generate_inchi_ID(self.product)

    def test_get_additional_routes(self):
        logging.basicConfig(level=logging.INFO)
        workshop = CobblersWorkshop(self.product, self.reactants, self.reaction_names,
                                    self.num_steps, self.output_dir, self.filter, id=self.id)
        routes = workshop.get_additional_routes(edit_route=True)
        self.assertIsInstance(routes, list)


if __name__ == '__main__':
    unittest.main()

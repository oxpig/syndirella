from syndirella.utils.input import process_routes
import pandas as pd

df = pd.read_csv('/Users/kate_fieseler/PycharmProjects/retrievesynthesizable/NCS1/d2r/smiles_routes_1_step.csv')
process_routes(df, smiles_col_idx=0, reaction_name_col_idx=2, reactants_col_idx=3, results_dir="NCS1/d2r", output_filename="1_step_to_elaborate.csv")

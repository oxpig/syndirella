import os

cli_default_settings = dict(
    ranking=os.environ.get('SYNDIRELLA_RANKING', 'num_atom_diff'),
    cutoff=float(os.environ.get('SYNDIRELLA_CUTOFF', 5.0)),
    rxn_smarts_path=os.environ.get('SYNDIRELLA_RXN_SMARTS_PATH', "/Users/kate_fieseler/PycharmProjects/retrievesynthesizable/RXN_SMARTS_CONSTANTS.json"),
    supressed_exceptions=(Exception,),
)

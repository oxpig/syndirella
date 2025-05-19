# AiZynthFinder

AiZynthFinder is a tool for outputting retrosynthesis routes with purchasable reactants for a given molecule and I'm
using it to find routes limited to
the SMIRKS library I already have.

## Method:

Given a molecule it is broken down into precursors that can be purchased or cannot be broken down further. Uses a Monte
Carlo tree search, where each node in the tree is a set of molecules. At each iteration a leaf node is selected and
deemed
to be most promising. A neural network policy then selects a reaction template to apply and create new precursours.
Procedure is repeated until a purchasable precursor has been found or the tree has reached a maximum depth.

Understanding route and reactions:

```python
finder = AiZynthFinder()  # AiZynthFinder object
routes = finder.routes  # routes found after performing search

for i, route in enumerate(routes):
    rxn_tree = route['reaction_tree']
    for step, reaction in enumerate(rxn_tree.reactions()):
        reaction_data = {
            'smiles': finder.target_smiles,  # queried smiles
            'route_id': i,
            'route_score': route['score']['state score'],
            # 0.95 * N solved precursors / 0.05 * N reactions (high is good)
            'reaction_metadata': reaction.metadata,
            'step': step,
            'n_steps': len(list(rxn_tree.reactions())),
            'product': reaction.mol.smiles,  # product of step
            'reactants': tuple([umol.smiles for umol_tuple in reaction.reactants for umol in umol_tuple]),
            'reaction': reaction.smiles,  # SMIRKS of reaction
            'label': None
        }
```

### Configuration:

You need to set the stock, expansion policy, and filter policy.

Stock: Collection of compounds that serve as stop-conditions for the tree search. Default is `zinc_stock.hdf5` compounds
from ZINC (molecular weight up to 250 D and log P up to 3.5) that had reactivity labeled as “standard” or “reactive”,
resulting in 17,422,831 compounds.

Expansion policy: Neural network trained to score reaction templates given a molecule at a node. Default is
`uspto_model.onnx` which according to [this paper](https://pubs.rsc.org/en/content/articlelanding/2020/sc/c9sc04944d)
contains 302,282 templates.

Filter policy: Model to filter out predicted reactions that are unlikely to be chemically viable. Most likely "trained
on positive reactions and a virtually enumerated set of negative reactions.", as
described [in this paper](https://pubs.rsc.org/en/content/articlelanding/2020/sc/c9sc04944d).

### How is performance affected b
# Syndirella 👑

Syndirella (Synthetically Directed Elaborations) is a tool for generating and scoring synthetically accessible 
elaborations of molecules designed from fragment screens. 

## Table of Contents

- [Installation](#installation)
- [Usage](#usage)
- [Input](#input)
- [Features](#features)
- [Important Notes](#important-notes)
- [Documentation](#documentation)
- [Built With](#built-with)
- [Authors](#authors)
- [License](#license)
- [Acknowledgments](#acknowledgments)

## Installation
```bash
git clone https://github.com/kate-fie/syndirella
cd syndirella
make
```

## Usage
```bash
syndirella --help
```

## Input
- Compounds to elaborate (base compounds) should recapitulate at least **2** fragment hits. 
- Please use the [syndirella_input_template.csv](syndirella/syndirella_input_template.csv) as a guide for input file formatting.
  - The format is as follows:
    - 'smiles' - SMILES string of the base compound
    - 'reactants' - a list of tuples of SMILES strings of reactants in the order of first to last step. Please use quotes around each SMILES string.
    - 'reaction_names' - a list of reaction names in the order of first to last step. Please use quotes around each reaction name.
    - 'num_steps' - the number of steps in route
    - 'hits' - a space seperated list of fragment hits
    - 'compound_set' - a name for the compound set or any additional identifier you would like to use
## Features

## Important Notes
- If you change the `RXN_SMARTS_CONSTANTS.json` file you will need to recompile all outputs for your library and products <br />
to accurately reflect the changes.







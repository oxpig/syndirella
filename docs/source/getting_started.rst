Getting Started
===============

This guide will help you get started with Syndirella, covering the installation and basic usage.

.. contents::
   :local:
   :depth: 2


Installation
------------

**Using the Makefile**

   .. code-block:: bash

       git clone https://github.com/kate-fie/syndirella.git
       cd syndirella
       make

**Using pip** NOT YET IMPLEMENTED

   .. code-block:: bash

       pip install syndirella

Setting Up
----------

Once you have installed Syndirella, follow these steps to set it up:

1. **Setting environment variables**:

   .. code-block:: bash

       export SYNDIRELLA_BASE_PATH=/path/to/syndirella
       export MANIFOLD_API_KEY=[API_KEY]

You need to have a Manifold account and an API key to use Syndirella. You can find more information 'here <https://api.postera.ai/api/v1/docs/>'_.

2. **Configuration of reaction constants** (if applicable):
    1. Reactions
        Syndirella has 20 encoded reactions which you can find in 'syndirella/constants/RXN_SMARTS_CONSTANTS.json'. You can edit this file to add more reactions to your liking.
    2. Alternative Reactions
        You can provide alternative reactions to your preference if they are easier or more successful. Add them in 'syndirella/constants/ADDITIONAL_RXN_OPTIONS.json'.
    3. Alternative Reactants
        If there are reactants for certain reactions that you prefer to use, you can add them in 'syndirella/constants/RXN_SMARTS_CONSTANTS.json'.

Basic Usage
-----------

Elaborate a set of base compounds using these steps:

1. **Create CSV of Base Compounds**:
    The input CSV must contain base_compound SMILES and fragment hits names exactly how they appear in the SDF. You can either choose to run the pipeline with manual or automatic route finding.

    **Manual Routes**:
    You can select the exact route to elaborate the base_compound with the reaction names, exact reactants, and number of steps in the route. An example template is at 'syndirella/syndirella_input_template_manual.csv'.

    **Automatic Routes**:
    Base compounds can be elaborated by routes automatically proposed by Manifold. All you need in the CSV is the base_compound SMILES and hits names. An example template is at 'syndirella_input_template.csv'.

2. **Download protein structures and fragment hits from Fragalysis**:
    Go to 'Fragalysis <https://fragalysis.xchem.diamond.ac.uk/viewer/react/landing>'_ and download a dataset.

3. **Optional: Relax protein template**:
    Choose the holo protein template you would like to run placements with and relax with PyRosetta. Use '

    Manually check fragment hits are where you are expecting them!

3. **Run pipeline**:






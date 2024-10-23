
==========
User Guide
==========

Installation
============

On Mac OS and Linux you can install from PyPI using Conda.

.. code-block:: bash

    pip install --upgrade syndirella

.. attention::

    Installation and usage have not been tested on Windows OS.

Setting Up
==========

Once you have installed Syndirella, follow these steps to set it up:

1. Manifold API Key:
-------------------------------------

Syndirella uses `Postera's Manifold <https://app.postera.ai/>`_ API to perform retrosynthesis searches of input scaffolds
and superstructure searches of reactants. You can find more information on Manifold API access
`here <https://api.postera.ai/api/v1/docs/>`_.

.. attention::

    A Manifold Pro account is required for API access.

Set the Manifold API key as an environment variable.

.. code-block:: bash

   export MANIFOLD_API_KEY=[API_KEY]

2. Configuration of reaction constants (optional):
-------------------------------------------------------

All reactions and alternative reaction/reactant specifications are formatted with
`SMIRKS <https://www.daylight.com/dayhtml_tutorials/languages/smirks/index.html>`_ or
`SMARTS <https://www.daylight.com/dayhtml_tutorials/languages/smarts/index.html>`_. If you'd like to
edit or add to the libraries, please feel free to add directly to the following files or open an issue and I can help in
formatting SMARTS/SMIRKS strings to your needs as they can be quite finicky.

**Reactions**:
Syndirella has >20 encoded reactions which you can find in
`RXN_SMARTS_CONSTANTS.json <https://github.com/kate-fie/syndirella/blob/13f73d8beda750c023739729d2681f5939d29e29/syndirella/constants/RXN_SMARTS_CONSTANTS.json>`_.

**Alternative Reactions**:
You can provide alternatives for reactions based on preference. Add them in
`ADDITIONAL_RXN_OPTIONS.json <https://github.com/kate-fie/syndirella/blob/13f73d8beda750c023739729d2681f5939d29e29/syndirella/constants/ADDITIONAL_RXN_OPTIONS.json>`_.

**Alternative Reactants**:
If there are reactants for certain reactions that you prefer to use, you can add them in
`REACTANT_FILTER_CONSTANTS.json <https://github.com/kate-fie/syndirella/blob/13f73d8beda750c023739729d2681f5939d29e29/syndirella/constants/REACTANT_FILTER_CONSTANTS.json>`_.

Please find more details and specification on editing these files in Reaction Constants.

Basic Usage
===========

Elaborate a set of scaffolds using these steps:

1. Setup fragments and protein templates
----------------------------------------

Download the fragment hits from Fragalysis. In the download folder the important files are:

::

    metadata.csv # contains short and long codes
    target_name_combined.sdf # fragment poses with long code names
    /aligned_files/fragment_name/fragment_name_apo-desolv.pdb # apo pdb used for placement

2. Create input csv
-------------------

Syndirella can be run either in *automatic* or *manual* mode.

**Automatic**:
    Scaffolds can be elaborated by routes automatically proposed by Manifold.
    An example template is at ``syndirella/syndirella_input_template.csv``.

Required headers:

    ``smiles``:
        smiles string of :term:`scaffold`.
    ``hit1``:
        string of short code (or any substring of short code found in the metadata.csv) of 1 fragment inspiring hit.
    ``template``:
        path to apo protein template to use for :term:`placement`.
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
    An example template is at ``syndirella/syndirella_input_template_manual.csv``.

Required headers:

    ``smiles``:
        smiles string of scaffold.
    ``hit1``:
        string of short code (or any substring of short code found in the metadata.csv) of 1 fragment inspiring hit.
    ``template``:
        path to apo protein template to use for :term:`placement`.
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

Run pipeline in *automatic* mode:

.. code-block:: bash

    syndirella --input [path_to_automatic.csv] --output [path_to_output_dir] --templates [path_to_templates_dir]
    --hits_path [path_to_fragments.sdf] --metadata [path_to_metadata.csv]


Run pipeline in *manual* mode:
    Add ``--manual`` flag.

4. Outputs
----------

**Output directory structure:**

🔑🔑🔑: Inchi key of scaffold. Example: ``ZJENMQHSGLZNHL-UHFFFAOYSA-N``

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
    │   │   ├── 🔑🔑🔑_[route_uuid]_[num]-[stereoisomer]
    │   │   │   ├── 🔑🔑🔑_[route_uuid]_[num]-[stereoisomer].mol
    │   │   │   ├── 🔑🔑🔑_[route_uuid]_[num]-[stereoisomer].json # energy values
    │   │   └── continued for all products...
    │   ├── 🔑🔑🔑_[route_uuid]_[rxn_name]_products_[last_step]of[total_steps].pkl.gz & .csv # final products
    │   ├── 🔑🔑🔑_[route_uuid]_[rxn_name]_products_[last_step]of[total_steps]_placements.pkl.gz & .csv # merged placements with products info
    │   ├── 🔑🔑🔑_[route_uuid]_fragmenstein_placements.pkl.gz & .csv # fragmenstein output
    │   └── 🔑🔑🔑_[route_uuid]_to_hippo.pkl.gz # full routes and placements
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

**🔑🔑🔑_[route_uuid]_to_hippo.pkl.gz:**
    Full routes and placements.

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


Usage Option: Only Place Scaffolds
==================================

You can run Syndirella to only place scaffolds. It will not perform the full elaboration procedure.

.. code-block:: bash

    syndirella --input [path_to_automatic.csv] --output [path_to_output_dir] --templates [path_to_templates_dir]
    --hits_path [path_to_fragments.sdf] --metadata [path_to_metadata.csv] --scaffold_place


Command Line Interface
======================

.. code-block:: bash

    usage: syndirella [-h] -i INPUT -o OUTPUT -t TEMPLATES --hits_path HITS_PATH --metadata METADATA [--products PRODUCTS] [--batch_num BATCH_NUM] [--manual] [--scaffold_place]
                  [--scaffold_place_num SCAFFOLD_PLACE_NUM] [--profile] [--atom_diff_min ATOM_DIFF_MIN] [--atom_diff_max ATOM_DIFF_MAX] [--long_code_column LONG_CODE_COLUMN]

    Run the Syndirella pipeline with specified configurations.

    options:
      -h, --help            show this help message and exit
      -i INPUT, --input INPUT
                            Input .csv file path for the pipeline. (default: None)
      -o OUTPUT, --output OUTPUT
                            Output directory for the pipeline results. (default: None)
      -t TEMPLATES, --templates TEMPLATES
                            Absolute path to a directory containing the template(s). (default: None)
      --hits_path HITS_PATH
                            Absolute path to hits_path for placements (.sdf or .mol). (default: None)
      --metadata METADATA   Absolute path to metadata for placements. (default: None)
      --products PRODUCTS   Absolute path to products for placements. (default: None)
      --batch_num BATCH_NUM
                            Batch number for processing. (default: 10000)
      --manual              Use manual routes for processing. (default: False)
      --scaffold_place      Only place scaffolds. Do not continue to elaborate. (default: False)
      --scaffold_place_num SCAFFOLD_PLACE_NUM
                            Number of times to attempt scaffold placement. (default: 5)
      --profile             Run the pipeline with profiling. (default: False)
      --atom_diff_min ATOM_DIFF_MIN
                            Minimum atom difference between elaborations and scaffold to keep. (default: 0)
      --atom_diff_max ATOM_DIFF_MAX
                            Maximum atom difference between elaborations and scaffold to keep. (default: 10)
      --long_code_column LONG_CODE_COLUMN
                            Column name for long code in metadata csv to match to SDF name. The column can contain a substring for the
                            SDF name. (default: Long code)


    Syndirella is installed at [path_to_installation]




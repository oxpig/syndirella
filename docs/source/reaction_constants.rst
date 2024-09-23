
==================
Reaction Constants
==================

All reactions and alternative reaction/reactant specifications are formatted with
`SMIRKS <https://www.daylight.com/dayhtml_tutorials/languages/smirks/index.html>`_ or
`SMARTS <https://www.daylight.com/dayhtml_tutorials/languages/smarts/index.html>`_. If you'd like to
edit or add to the libraries, please feel free to add directly to the following files in your forked version or open an issue and I can help in
formatting SMARTS/SMIRKS strings to your needs as they can be quite finicky.

These are found at `syndirella/syndirella/constants <https://github.com/kate-fie/syndirella/tree/e563796e62c604d08aa9ee16beed26a9eee694c0/syndirella/constants>`_.

`RXN_SMIRKS_CONSTANTS.json <https://github.com/kate-fie/syndirella/blob/e563796e62c604d08aa9ee16beed26a9eee694c0/syndirella/constants/RXN_SMIRKS_CONSTANTS.json>`_
-------------------------

This file defines the SMIRKS for reactions.

The ``.json`` structure follows this format:

- **Key (e.g., "Amidation")**: The name of the reaction (with ``_`` as spaces)

- **Value (e.g., "[#6:1](=[#8:2])-[#8;H1].[$([N+0&H1,N+0&H2]);!$(NC=*);!$(NS);!$(N=*);!$(N-O);!$(N-o):3]>>[#6:1](=[#8:2])-[#7X3:3]")**: The reaction SMIRKS.


`REACTANT_FILTER_CONSTANTS.json <https://github.com/kate-fie/syndirella/blob/e563796e62c604d08aa9ee16beed26a9eee694c0/syndirella/constants/REACTANT_FILTER_CONSTANTS.json>`_
----------------------

This file defines options to transform reactants in specific reactants. Such as transforming a boronic acid to a boronate ester. The transformed
reactant is used as an additional reactant in the library of reactants. It used to structure inputs for ``rdkit.Chem.ReplaceSubstructs()``.

The ``.json`` structure follows this format:

- **Key (e.g., "Sp2-sp2_Suzuki_coupling")**: This indicates the specific reaction type being targeted.

- **SMARTS Pattern (e.g., "[#5](-[#8])(-[#8])")**: The SMARTS pattern to be modified for a reactant. Each reaction type can have multiple SMARTS patterns that describe specific substructures to look for in the reactants. These patterns are used to identify the appropriate atoms for modification.

- **Values**:
    - ``name``: Describes the action to take (e.g., `"add_boronate_ester"`). This is the modification to be applied when the SMARTS pattern is found.
    - ``to_add``: Specifies the SMARTS string of the atoms to be added to the structure. For instance, a boronate ester group ("[#6]-[#6]1(-[#6])-[#8]-[#5]-[#8]-[#6]-1(-[#6])-[#6]") might be added.
    - ``connecting_atom_id``: Specifies the index (0-based) of the atom in ``to_add`` that will connect to the existing structure. This number corresponds to the atom's position in the SMARTS string.

Example entry:

.. code-block:: json

    {
      "Sp2-sp2_Suzuki_coupling": {
        "[#5](-[#8])(-[#8])": {
          "name": "add_boronate_ester",
          "to_add": "[#6]-[#6]1(-[#6])-[#8]-[#5]-[#8]-[#6]-1(-[#6])-[#6]",
          "connecting_atom_id": 4
        },
        "[#9,#17,#53]": {
          "name": "add_bromine",
          "to_add": "[#35]",
          "connecting_atom_id": 0
        }
      }
    }


`ADDITIONAL_RXN_OPTIONS.json <https://github.com/kate-fie/syndirella/blob/e563796e62c604d08aa9ee16beed26a9eee694c0/syndirella/constants/ADDITIONAL_RXN_OPTIONS.json>`_
--------------------------------

This file specifies reactions that should be replaced by another reaction to make a new route for the pipeline. Similar
to ``REACTANT_FILTER_CONSTANTS.json`` as it used to structure inputs for ``rdkit.Chem.ReplaceSubstructs()``. But instead
of adding the edited reactant to the library of reactants of one reaction, it is creating a completely new route.

The ``.json`` structure follows this format:

- ``name``: The reaction name that triggers for an additional route to be created with this reaction replaced  (e.g., "Amide_schotten-baumann"). Must match the reaction name in the ``RXN_SMIRKS_CONSTANTS.json``.

- ``replace_with``: The new reaction to replace the original reaction in a new route. This is the reaction you're swapping in (e.g., "Amidation").

- ``reactant_smarts_to_replace``: A SMARTS string of the atoms in the reactant you are directly editing for use in the new reaction.

- ``reactant_id_to_replace``: Specifies which reactant (by index, starting from 1) is being replaced in the original reaction SMIRKS.

- ``reactant_smarts_to_replace_with``: A SMARTS string specifying the reactant pattern that will replace the atoms found in ``reactant_smarts_to_replace``.

- ``replacement_connecting_atom_id``: The index of the atom in the ``reactant_smarts_to_replace_with`` that connects to the rest of the reactant. The index is 0-based, meaning it refers to the atom index in SMARTS minus 1.

Example entry:

.. code-block:: json

    {
        "name": "Amide_schotten-baumann",
        "replace_with": "Amidation",
        "reactant_smarts_to_replace": "[#6X3;!$(C-N):1](=[OX1:2])-[#17,#9,#35:3]",
        "reactant_id_to_replace": 2,
        "reactant_smarts_to_replace_with": "[#6X3;!$(C-N):1](=[OX1:2])-[#8;H1:3]",
        "replacement_connecting_atom_id": 0
    }

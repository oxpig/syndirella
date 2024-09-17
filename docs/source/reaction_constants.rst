
==================
Reaction Constants
==================

Syndirella uses SMARTS to define reactions and alternative reaction options. You can edit these constants to suit your needs.

.. contents::
   :local:
   :depth: 2

`ADDITIONAL_RXN_OPTIONS_V2.json`
--------------------------------
name: Amide_schotten-baumann,
"replace_with": "Amidation",
"reactant_smarts_to_replace": "[#6X3;!$(C-N):1](=[OX1:2])-[#17,#9,#35:3]",
"reactant_id_to_replace": 2,
"reactant_smarts_to_replace_with": "[#6X3;!$(C-N):1](=[OX1:2])-[#8;H1:3]",
"replacement_connecting_atom_id": 0 # needs to be atom index in SMARTS minus 1. In this case, 0 corresponds to the carbon with index 1.

`
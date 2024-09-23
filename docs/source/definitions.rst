
===========
Definitions
===========

Classes
-------

.. glossary::

    Cobbler
        Represents a full elaboration for one scaffold. So it can contain multiple routes (multiple CobblersWorkshops).

    CobblersWorkshop
        Represents a single route for a scaffold. It contains a list of reactants and reaction names that describe
        each step in the route. Each step is represented by a CobblersBench object. If a step has a replacement for it,
        a new CobblersWorkshop is created with the replacement.

    CobblerBench
        Represents a single step in a route. Performs the whole process of finding analogues of reactants. Since the
        elaborated products are 'slippers' in this analogy, the CobblerBench is where the material to make the slippers
        (superstructure of reactants) are found.

    Reaction
        Represents a single reaction within the CobblerBench. Contains the reactants, products, and the name of the reaction.
        Assigns the reactants to the SMARTS pattern of the specific reactant in the full reaction SMIRKS.

    Library
        Contains the superstructures of reactant scaffolds found. Performs filtering for reactants before finding products.

    Slipper
        Handles finding the products of a Library and placing them. Is instantiated for each step of a route.

    SlipperSynthesizer
        Finds the products of a Library object with a cartesian multiplication of the superstructures of reactants after
        filtering.

    SlipperFitter
        'Fits' the final products into the protein template by performing constrained energy minimisation to the experimental
        fragment structures using the matched atoms. This placement procedure is carried out by
        `Fragmenstein <https://github.com/matteoferla/Fragmenstein>`_.

Modules
-------

.. glossary::

    Fairy
        This module provides functions to find similar cheaper reactants, filter out reactants based on simple filters,
        fingerprint generation, etc. Provides the 'fairy dust' to improve syntheses structured by ML tools.


General Terminology
-------------------

.. glossary::

    Scaffold
        Compound to elaborate from.

    Placement
        Constrained docking to the fragment inspirations.
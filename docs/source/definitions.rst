
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

    Slipper
        Represents a single product of a reaction. Contains the SMILES of the product and the reaction name that was
        used to generate it.





General Terminology
-------------------

.. glossary::

    Scaffold
        Compound to elaborate from.

    Placement
        Constrained docking to the fragment inspirations.
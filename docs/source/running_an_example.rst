Running an example
===============

This guide will walk you through an example to run Syndirella's pipeline.

.. contents::
   :local:
   :depth: 2

1. Make Input CSV

First we need to prepare our inputs into a .csv file. I usually do this using a combination of Postera's Manifold
(to draw structures in, label designs with their fragment inspiration, and download CSV) and excel to further edit.

Let's start with three designs that I've manually designed:

[2D pics]

I already have a route in mind for the first one but for the next two I'll have Manifold give me one. So that means we
need to structure two input csvs, a manual format for the first compound and an automatic format for the last two.

Manual route:

[CSV]

The CSV must contain the columns:

:code:`smiles`:

:code:`reactants`:

:code:`reaction_names`:

:code:`num_steps`:

:code:`hits`:

Optional columns:

:code:`compound_set`:

2. Prepare templates

3. Structure directories

4. Prepare job script

5. Run!




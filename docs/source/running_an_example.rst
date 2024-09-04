
==================
Running an example
==================

This guide will walk you through an example to run Syndirella's pipeline. All resources to run this example are found here.

Preparing the inputs
====================

I carry out these same steps in the example jupyter notebook.

Fragments and templates

Follow-up designs

Make input csv

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

Structure directories

Run
Using job script
Locally

Analysing the outputs
=====================

Structure of _output_YYYYMMDD_HHMM.csv

total_num_products_enumstereo:
    Total number of products enumerated with stereochemistry in the final step. This is counting the number of unique
    products with stereochemistry, so if a product with stereochemistry is generated multiple times via different routes
    it will only be counted once.

total_num_unique_products:
    Total number of unique products without stereochemistry in the final step. If a product is generated multiple times
    by different routes it will only be counted once.

Understanding structure of output directory


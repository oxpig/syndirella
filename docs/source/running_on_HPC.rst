
==============
Running on HPC
==============

Syndirella can be run easily via CLI on a HPC. While Syndirella itself
does not have the capability to run in parallel on multiple threads, I've
found running simultaneous jobs for each scaffold with low memory (3GB) to
significantly increase speed.

Preparing the input files to run on HPC with example scripts for SLURM is
found at `syndirella/hpc <https://github.com/kate-fie/syndirella/tree/e6f744415704e63f13608f7f80bcc7ed962410ed/hpc>`_.

Follow the example jupyter notebook to prepare everything, `prepare_csvs_jobs.ipynb <https://github.com/kate-fie/syndirella/blob/fc9c087c6ee276b404a7a226d1a076ed12e3e6a0/hpc/prep_csvs_jobs.ipynb>`_.



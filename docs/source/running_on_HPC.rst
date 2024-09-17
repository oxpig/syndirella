
==============
Running on HPC
==============

Syndirella can be run easily via CLI on a HPC. While Syndirella itself
does not have the capability to run in parallel on multiple threads, I've
found running simultaneous jobs for each scaffold with low memory (3GB) to
significantly increase speed.

Preparing the input files to run on HPC with example scripts for SLURM is
found at `syndirella/hpc <https://github.com/kate-fie/syndirella/tree/3aeabb1462f5f02cfbda41405e59e2e1efd6d7d1/syndirella/hpc>`_.

Follow the example jupyter notebook to prepare everything, `prepare_csvs_jobs.ipynb <https://github.com/kate-fie/syndirella/blob/3aeabb1462f5f02cfbda41405e59e2e1efd6d7d1/syndirella/hpc/prep_csvs_jobs.ipynb>`_.



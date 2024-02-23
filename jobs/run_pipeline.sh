#!/bin/bash
# For scripts required see: https://gist.github.com/matteoferla/e0496d5766c12a0ae1738b943b41a536
# A few things don't work in CentOS 7 due to GNU lib C (glibc) 2.17, so it has to be run in a Singularity container

#: << USAGE #########
# source run_job.env
# condor_submit /data/xchem-fragalysis/shared/target_script.condor
#USAGE #########

#: << INSTALLATION #############
#export PIP_NO_USER=1; # gets ignored.
#export PYTHONUSERBASE=$CONDA_PREFIX; # double tap
#conda install -y -n base -c conda-forge openbabel plip;
#pip install -q fragmenstein
#pip install -q pyrosetta-help
#PYROSETTA_USERNAME=levinthal PYROSETTA_PASSWORD=paradox install_pyrosetta
#INSTALLATION ############

export HOST=${HOST:-$(hostname)}
export USER=${USER:-$(users)}
export HOME=${HOME:-$_CONDOR_SCRATCH_DIR}
source /etc/os-release;
echo "Running script ${0} as $USER in $HOST which runs $PRETTY_NAME"
# ---------------------------------------------------------------

source /data/xchem-fragalysis/kfieseler/.bashrc
conda activate fragmenstein

cd /data/xchem-fragalysis/kfieseler/syndirella

pwd;
export INPUT="/data/xchem-fragalysis/kfieseler/syndirella/jobs/syndirella_input.csv"
export OUTPUT_DIR="/data/xchem-fragalysis/kfieseler/A71EV2A_run2/test";
export TEMPLATE="/data/xchem-fragalysis/kfieseler/syndirella/exps/fragments/x0310_template.pdb";
export HITS="/data/xchem-fragalysis/kfieseler/syndirella/exps/fragments/clean_hits.sdf";

echo "Running syndirella pipeline"
echo "INPUT: $INPUT"
echo "OUTPUT_DIR: $OUTPUT_DIR"
echo "TEMPLATE: $TEMPLATE"
echo "HITS: $HITS"

nice -19 python jobs/run_job.py \
--input $INPUT \
--output $OUTPUT_DIR \
--template $TEMPLATE \
--hits $HITS;

echo 'COMPLETE'
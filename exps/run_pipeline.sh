#! bin/bash

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

cd /data/xchem-fragalysis/kfieseler/repos/elaborate

pwd;
export HOME_DIR="/data/xchem-fragalysis/kfieseler/syndirella/exps/original";

echo "Running syndirella pipeline"
echo "HOME_DIR: $HOME_DIR"

nice -19 python exps/run_job.py \
--output $HOME_DIR;

echo 'COMPLETE'
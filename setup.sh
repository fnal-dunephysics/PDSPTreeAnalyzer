export PDSPAna_WD=`pwd`
export PDSPAna_LIB_PATH=$PDSPAna_WD/lib/
mkdir -p $PDSPAna_LIB_PATH
mkdir -p $PDSPAna_WD/tar

export PDSPAnaV="v1"
mkdir -p $PDSPAna_WD/data/$PDSPAnaV
export DATA_DIR=$PDSPAna_WD/data/$PDSPAnaV

export PDSPAnaGridOutDir="/pnfs/dune/scratch/users/$USER/PDSP_out"
mkdir -p $PDSPAnaGridOutDir

#### USER INFO ####
export PDSPAnaLogEmail='sungbino@fnal.gov'
export PDSPAnaLogWeb=''
export PDSPAnaLogWebDir=''

#### use cvmfs for root ####
#source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh
source /cvmfs/fermilab.opensciencegrid.org/packages/common/setup-env.sh ## -- For Alma 9
source /cvmfs/larsoft.opensciencegrid.org/spack-packages/setup-env.sh

#setup root v6_22_08d -q e20:p392:prof
spack load root@6.28.12

if [[ $HOSTNAME == *"dunegpvm"*"fnal.gov" ]]; then

  echo "@@@@ Working on dunegpvm"
  export PDSPAnaRunlogDir="/dune/app/users/$USER/Grid/PDSPAnaRunlog/"
  export PDSPAnaOutputDir="/dune/app/users/$USER/Grid/PDSPAnaOutput/"

fi

alias pdout="cd $PDSPAnaOutputDir/$PDSPAnaV/"

export MYBIN=$PDSPAna_WD/bin/
export PYTHONDIR=$PDSPAna_WD/python/
export PATH=${MYBIN}:${PYTHONDIR}:${PATH}

export ROOT_INCLUDE_PATH=/cvmfs/larsoft.opensciencegrid.org/spack-packages/opt/spack/linux-almalinux9-x86_64_v2/gcc-12.2.0/root-6.28.12-sfwfmqorvxttrxgfrfhoq5kplou2pddd/include/:$ROOT_INCLUDE_PATH
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PDSPAna_LIB_PATH:/cvmfs/larsoft.opensciencegrid.org/spack-packages/opt/spack/linux-almalinux9-x86_64_v2/gcc-12.2.0/root-6.28.12-sfwfmqorvxttrxgfrfhoq5kplou2pddd/lib/

source $PDSPAna_WD/bin/BashColorSets.sh

## submodules ##
#source bin/CheckSubmodules.sh

## Todo list ##
python python/PrintToDoLists.py
source $PDSPAna_WD/tmp/ToDoLists.sh
rm $PDSPAna_WD/tmp/ToDoLists.sh

## Log Dir ##
echo "* Your Log Directory Usage"
#du -sh $PDSPAnaRunlogDir
echo "-----------------------------------------------------------------"
CurrentGitBranch=`git branch | grep \* | cut -d ' ' -f2`
printf "> Current PDSPAnaAnalyzer branch : "${BRed}$CurrentGitBranch${Color_Off}"\n"

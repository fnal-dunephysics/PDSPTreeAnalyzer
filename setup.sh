export PDSPAna_WD=`pwd`
export PDSPAna_LIB_PATH=$PDSPAna_WD/lib/
mkdir -p $PDSPAna_LIB_PATH
mkdir -p $PDSPAna_WD/tar

export PDSPAnaV="v1"
mkdir -p $PDSPAna_WD/data/$PDSPAnaV
export DATA_DIR=$PDSPAna_WD/data/$PDSPAnaV

#### USER INFO ####
export PDSPAnaLogEmail='sungbino@fnal.gov'
export PDSPAnaLogWeb=''
export PDSPAnaLogWebDir=''

#### use cvmfs for root ####
export CMS_PATH=/cvmfs/cms.cern.ch
source $CMS_PATH/cmsset_default.sh
export SCRAM_ARCH=slc7_amd64_gcc10
export cmsswrel='cmssw-patch/CMSSW_12_6_0_patch1'
cd /cvmfs/cms.cern.ch/$SCRAM_ARCH/cms/$cmsswrel/src
echo "@@@@ SCRAM_ARCH = "$SCRAM_ARCH
echo "@@@@ cmsswrel = "$cmsswrel
echo "@@@@ scram..."
eval `scramv1 runtime -sh`
cd -
source /cvmfs/cms.cern.ch/$SCRAM_ARCH/cms/$cmsswrel/external/$SCRAM_ARCH/bin/thisroot.sh

if [[ $HOSTNAME == *"dunegpvm"*"fnal.gov" ]]; then

  echo "@@@@ Working on dunegpvm"
  export PDSPAnaRunlogDir="/dune/app/users/$USER/Grid/PDSPAnaRunlog/"
  export PDSPAnaOutputDir="/dune/app/users/$USER/Grid/PDSPAnaOutput/"

fi

alias pdout="cd $PDSPAnaOutputDir/$PDSPAnaV/"

export MYBIN=$PDSPAna_WD/bin/
export PYTHONDIR=$PDSPAna_WD/python/
export PATH=${MYBIN}:${PYTHONDIR}:${PATH}

export ROOT_INCLUDE_PATH=$ROOT_INCLUDE_PATH:$PDSPAna_WD/DataFormats/include/:$PDSPAna_WD/AnalyzerTools/include/:$PDSPAna_WD/Analyzers/include/
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$PDSPAna_LIB_PATH

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

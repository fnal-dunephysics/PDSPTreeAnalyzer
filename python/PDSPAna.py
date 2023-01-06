#!/usr/bin/env python

import os,sys,time
import argparse
import datetime
from CheckJobStatus import *
from TimeTools import *
import random
import subprocess

## Arguments

parser = argparse.ArgumentParser(description='PDSPAna Command')
parser.add_argument('-a', dest='Analyzer', default="")
parser.add_argument('-i', dest='InputSample', default="")
parser.add_argument('-l', dest='InputSampleList', default="")
parser.add_argument('-n', dest='NJobs', default=1, type=int)
parser.add_argument('-o', dest='Outputdir', default="")
parser.add_argument('-q', dest='Queue', default="fastq")
parser.add_argument('-P', dest='Momentum', default="1.0")
parser.add_argument('--no_exec', action='store_true')
parser.add_argument('--userflags', dest='Userflags', default="")
parser.add_argument('--nmax', dest='NMax', default=0, type=int)
parser.add_argument('--reduction', dest='Reduction', default=1, type=float)
parser.add_argument('--memory', dest='Memory', default=0, type=float)
parser.add_argument('--batchname',dest='BatchName', default="")
args = parser.parse_args()

## make userflags as a list
Userflags = []
if args.Userflags != "":
  Userflags = (args.Userflags).split(',')

## Add Abosolute path for outputdir
if args.Outputdir!='':
  if args.Outputdir[0]!='/':
    args.Outputdir = os.getcwd()+'/'+args.Outputdir

## TimeStamp

# 1) dir/file name style
JobStartTime = datetime.datetime.now()
timestamp =  JobStartTime.strftime('%Y_%m_%d_%H%M%S')
# 2) log style
JobStartTime = datetime.datetime.now()
string_JobStartTime =  JobStartTime.strftime('%Y-%m-%d %H:%M:%S')
string_ThisTime = ""

## Environment Variables

USER = os.environ['USER']
SCRAM_ARCH = os.environ['SCRAM_ARCH']
cmsswrel = os.environ['cmsswrel']
PDSPAna_WD = os.environ['PDSPAna_WD']
PDSPAnaV = os.environ['PDSPAnaV']
SAMPLE_DATA_DIR = PDSPAna_WD+'/data/'+PDSPAnaV+'/sample/'
PDSPAnaRunlogDir = os.environ['PDSPAnaRunlogDir']
PDSPAnaOutputDir = os.environ['PDSPAnaOutputDir']
PDSPAna_LIB_PATH = os.environ['PDSPAna_LIB_PATH']
UID = str(os.getuid())
HOSTNAME = os.environ['HOSTNAME']
SampleHOSTNAME = HOSTNAME

## Check hostname
IsDUNEgpvm = ("dunegpvm" in HOSTNAME)
if IsDUNEgpvm:
  HOSTNAME = "DUNEgpvm"
  SampleHOSTNAME = "FNAL"

## Make Sample List

InputSamples = []
StringForHash = ""

## When using txt file for input (i.e., -l option)

if args.InputSampleList is not "":
  lines = open(args.InputSampleList)
  for line in lines:
    if "#" in line:
      continue
    line = line.strip('\n')
    InputSamples.append(line)
    StringForHash += line
else:
  InputSamples.append(args.InputSample)
  StringForHash += args.InputSample
FileRangesForEachSample = []

## add flags to hash
for flag in Userflags:
  StringForHash += flag

## Define MasterJobDir
MasterJobDir = PDSPAnaRunlogDir+'/'+timestamp+'__'+args.Analyzer+'__'+'Momentum'+args.Momentum
for flag in Userflags:
  MasterJobDir += '__'+flag
MasterJobDir += '__'+HOSTNAME+'/'

## Copy libray
os.system('mkdir -p '+MasterJobDir+'/tar/lib/')
os.system('cp '+PDSPAna_LIB_PATH+'/* '+MasterJobDir+'/tar/lib')

## Loop over samples

# true or false for each sample
SampleFinishedForEachSample = []
PostJobFinishedForEachSample = []
BaseDirForEachSample = []
XsecForEachSample = []
for InputSample in InputSamples:

  NJobs = args.NJobs

  SampleFinishedForEachSample.append(False)
  PostJobFinishedForEachSample.append(False)

  ## Global Varialbes

  IsDATA = False
  DataPeriod = ""
  if "Data" in InputSample:
    IsDATA = True

  ## Prepare output

  base_rundir = MasterJobDir
  base_rundir = base_rundir+"/"

  os.system('mkdir -p '+base_rundir)

  ## Get Sample Path

  lines_files = []

  tmpfilepath = SAMPLE_DATA_DIR+'/For'+SampleHOSTNAME+'/'+InputSample+'.txt'
  lines_files = open(tmpfilepath).readlines()
  os.system('cp '+tmpfilepath+' '+base_rundir+'/input_filelist.txt')

  NTotalFiles = len(lines_files)
  ## Write run script

  if IsDUNEgpvm:
    commandsfilename = args.Analyzer+'_'+args.Momentum+'_'+InputSample
    run_commands = open(base_rundir+'/'+commandsfilename+'.sh','w')
    print>>run_commands,'''#!/bin/bash
# Setup grid submission

outDir=$1
echo "@@ outDir : ${{outDir}}"

outDir=/pnfs/dune/persistent/users/sungbino/PDSP_out

echo "@@ pwd"
pwd
gridBaseDir='pwd'

echo "@@ mkdir output"
mkdir output
thisOutputCreationDir='pwd'/output

echo "@@ ls"
ls

nProcess=$PROCESS
echo "@@ nProcess : "${{nProcess}}

#### use cvmfs for root ####
export CMS_PATH=/cvmfs/cms.cern.ch
source $CMS_PATH/cmsset_default.sh
export SCRAM_ARCH=slc7_amd64_gcc700
export cmsswrel='cmssw-patch/CMSSW_10_4_0_patch1'
cd /cvmfs/cms.cern.ch/$SCRAM_ARCH/cms/$cmsswrel/src
echo "@@@@ SCRAM_ARCH = "$SCRAM_ARCH
echo "@@@@ cmsswrel = "$cmsswrel
echo "@@@@ scram..."
eval `scramv1 runtime -sh`
cd -
source /cvmfs/cms.cern.ch/$SCRAM_ARCH/cms/$cmsswrel/external/$SCRAM_ARCH/bin/thisroot.sh

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${{gridBaseDir}}/lib/
echo ${{LD_LIBRARY_PATH}}

echo "@@ Env setup finished"

#root -l -b -q run.C

echo "@@ Finished!!"
ls -alg ${{thisOutputCreationDir}}/

echo "ifdh cp "${{thisOutputCreationDir}}/${{outName}} ${{outDir}}/${{outName}}
ifdh cp ${{thisOutputCreationDir}}/${{outName}} ${{outDir}}/${{outName}}
echo "@@ Done!"

'''.format(MasterJobDir, base_rundir, SCRAM_ARCH, cmsswrel)
    run_commands.close()

  thisjob_dir = base_rundir

  runCfileFullPath = ""
  runfunctionname = "run"
  runCfileFullPath = base_rundir+'/tar/run.C'
    
  IncludeLine = ''

  out = open(runCfileFullPath, 'w')
  print>>out,'''{2}

void {1}(){{

  {0} m;


'''.format(args.Analyzer, runfunctionname, IncludeLine)

  out.write('  m.LogEvery = 100\n')
  out.write('  m.MCSample = "'+InputSample+'";\n')
  out.write('  m.Beam_Momentum = "'+str(args.Momentum)+'";\n')
  out.write('  m.SetTreeName();\n')
  for it_file in lines_files:
    thisfilename = it_file.strip('\n')
    out.write('  m.AddFile("'+thisfilename+'");\n')

  out.write('  m.SetOutfilePath("hists.root");\n')
  if args.Reduction>1:
    out.write('  m.MaxEvent=m.fChain->GetEntries()/'+str(args.Reduction)+';\n')
  if args.NMax > 0:
    out.write('  m.MaxEvent='+str(args.NMax)+';\n')
  print>>out,'''  m.Init();
  m.initializeAnalyzer();
  m.initializeAnalyzerTools();
  m.SwitchToTempDir();
  m.Loop();

  m.WriteHist();

}'''


  if IsDUNEgpvm:

    cwd = os.getcwd()
    os.chdir(base_rundir)
    os.system('tar -cvf PDSPAna_build.tar tar')
    if not args.no_exec:
      job_submit = '''jobsub_submit -G dune --role=Analysis --resource-provides=\"usage_model=DEDICATED,OPPORTUNISTIC\" -l \'+SingularityImage=\\\"/cvmfs/singularity.opensciencegrid.org/fermilab/fnal-wn-sl7:latest\\\"\' --append_condor_requirements=\'(TARGET.HAS_SINGULARITY=?=true)\' --tar_file_name \"dropbox:///{0}/PDSPAna_build.tar\" --email-to sungbin.oh555@gmail.com -N 1 --expected-lifetime 2h --memory 30GB "file://$(pwd)/{1}.sh\"'''.format(base_rundir, commandsfilename)
      os.system(job_submit)
      print(job_submit)
    os.chdir(cwd)

if args.no_exec:
  exit()




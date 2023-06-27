#!/usr/bin/env python

import os,sys,time
import argparse
import datetime
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
parser.add_argument('-P', dest='Momentum', default="1.0")
parser.add_argument('--suffix', dest='Suffix', default="")
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
PDSPAna_WD = os.environ['PDSPAna_WD']
PDSPAnaV = os.environ['PDSPAnaV']
SAMPLE_DATA_DIR = PDSPAna_WD+'/data/'+PDSPAnaV+'/sample/'
PDSPAnaRunlogDir = os.environ['PDSPAnaRunlogDir']
PDSPAnaOutputDir = os.environ['PDSPAnaOutputDir']
PDSPAna_LIB_PATH = os.environ['PDSPAna_LIB_PATH']
PDSPAnaGridOutDir = os.environ['PDSPAnaGridOutDir']
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

if args.InputSampleList != "":
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

## Copy headers
os.system('mkdir -p '+MasterJobDir+'/tar/include/')
os.system('cp '+PDSPAna_WD+'/DataFormats/include/*.h '+MasterJobDir+'/tar/include')
os.system('cp '+PDSPAna_WD+'/AnalyzerTools/include/*.h '+MasterJobDir+'/tar/include')
os.system('cp '+PDSPAna_WD+'/Analyzers/include/*.h '+MasterJobDir+'/tar/include')
os.system('cp -r '+PDSPAna_WD+'/data '+MasterJobDir+'/tar/')

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

  ## Get number of total entries

  sample_info_path = SAMPLE_DATA_DIR+'/CommonSampleInfo/'+InputSample+'.txt'
  lines_sample_info = open(sample_info_path).readlines()
  Total_entries = 0
  for line_sample_info in lines_sample_info:
    if '#' in line_sample_info:
      continue
    this_split = line_sample_info.split('\t')
    Total_entries = int(this_split[1][0:-1])
    print(Total_entries)
    


  ## Loop for number of jobs
  N_jobs = args.NJobs
  ## -- Loop for different set of modified model parameter values
  alpha_start = 0.875
  alpha_end = 0.975
  alpha_step = 0.005
  N_alpha_steps = int((alpha_end - alpha_start + 0.001) / alpha_step)

  beta_start = 0.185
  beta_end = 0.265
  beta_step = 0.002
  N_beta_steps = int((beta_end - beta_start + 0.001) / beta_step)

  for i_a in range(0, N_alpha_steps):
    for i_b in range(0, N_beta_steps):
      this_alpha = "{:.3f}".format(alpha_start + alpha_step * float(i_a))
      this_beta = "{:.3f}".format(beta_start + beta_step * float(i_b))
      print('Run for alpha :' + this_alpha + ', beta : ' + this_beta)
      
      ## Write run script

      if IsDUNEgpvm:
        commandsfilename = args.Analyzer+'_'+args.Momentum+'_'+InputSample
        outname = commandsfilename + args.Suffix + '_' + str(i_a) + '_' + str(i_b) + '.root'
        #outname = commandsfilename + '_test.root'
        run_commands = open(base_rundir+'/'+commandsfilename+'_'+str(i_a)+'_'+ str(i_b)+'.sh','w')
        run_commands.write('''#!/bin/bash
# Setup grid submission

outDir=$1
echo "@@ outDir : ${{outDir}}"

outDir={3}

echo "@@ pwd"
pwd
gridBaseDir=`pwd`

echo "@@ mkdir output"
mkdir output
thisOutputCreationDir=`pwd`/output

echo "@@ ls"
ls

echo "@@ ls tar"
ls -lhs ${{INPUT_TAR_DIR_LOCAL}}/tar
ls -lhs ${{INPUT_TAR_DIR_LOCAL}}/tar/lib
ls -lhs ${{INPUT_TAR_DIR_LOCAL}}/tar/include/
ls -lhs ${{INPUT_TAR_DIR_LOCAL}}/tar/data/v1/Momentum_reweight


nProcess=$PROCESS
echo "@@ nProcess : "${{nProcess}}

#### use cvmfs for root ####
source /cvmfs/dune.opensciencegrid.org/products/dune/setup_dune.sh
setup root v6_22_08d -q e20:p392:prof

#### setup directories
export ROOT_INCLUDE_PATH=${{INPUT_TAR_DIR_LOCAL}}/tar:$ROOT_INCLUDE_PATH
echo ${{ROOT_INCLUDE_PATH}}

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${{INPUT_TAR_DIR_LOCAL}}/tar/lib/
echo ${{LD_LIBRARY_PATH}}

export PDSPAnaV="v1"
export DATA_DIR=${{INPUT_TAR_DIR_LOCAL}}/tar/data/${{PDSPAnaV}}

setup ifdhc
export IFDH_CP_MAXRETRIES=2

echo "@@ Env setup finished"

root -l -b -q ${{INPUT_TAR_DIR_LOCAL}}/tar/run_{1}_{2}.C

echo "@@ Finished!!"
ls -alg ${{thisOutputCreationDir}}/

echo "ifdh cp "${{thisOutputCreationDir}}/{0} ${{outDir}}/{0}
ifdh rm ${{outDir}}/{0}
ifdh cp ${{thisOutputCreationDir}}/{0} ${{outDir}}/{0}
echo "@@ Done!"

        '''.format(outname, str(i_a), str(i_b), PDSPAnaGridOutDir))
      run_commands.close()

      thisjob_dir = base_rundir

      runCfileFullPath = ""
      runfunctionname = 'run_'+str(i_a)+'_'+str(i_b)
      runCfileFullPath = base_rundir+'/tar/run_'+str(i_a)+'_'+str(i_b)+'.C'
    
      IncludeLine = ''

      out = open(runCfileFullPath, 'w')
      out.write('''
R__LOAD_LIBRARY(libPhysics.so)
R__LOAD_LIBRARY(libMathMore.so)
R__LOAD_LIBRARY(libTree.so)
R__LOAD_LIBRARY(libHist.so)
R__LOAD_LIBRARY(libGpad.so)
R__LOAD_LIBRARY(libGraf3d.so)
R__LOAD_LIBRARY(libDataFormats.so)
R__LOAD_LIBRARY(libAnalyzerTools.so)
R__LOAD_LIBRARY(libAnalyzers.so)
void {1}(){{

  {0} m;


'''.format(args.Analyzer, runfunctionname))

      out.write('  m.LogEvery = 100;\n')
      out.write('  m.MCSample = "'+InputSample+'";\n')
      out.write('  m.Beam_Momentum = '+str(args.Momentum)+';\n')
      out.write('  m.alpha_new = '+this_alpha+';\n')
      out.write('  m.beta_new = '+this_beta+';\n')
      out.write('  m.SetTreeName();\n')
      for it_file in lines_files:
        thisfilename = it_file.strip('\n')
      out.write('  m.AddFile("'+thisfilename+'");\n')

      out.write('  m.SetOutfilePath("./output/' + outname + '");\n')
      out.write('  m.MaxEvent=-1;\n')
      out.write('''  m.Init();
  m.initializeAnalyzer();
  m.initializeAnalyzerTools();
  m.SwitchToTempDir();
  m.Loop();
  m.Extract_MPVs();

}''')
      out.close()

  if IsDUNEgpvm:

    cwd = os.getcwd()
    os.chdir(base_rundir)
    os.system('tar -cvf PDSPAna_build.tar tar')


    if not args.no_exec:
      for i_a in range(0, N_alpha_steps):
        for i_b in range(0, N_beta_steps):
          commandsfilename = args.Analyzer+'_'+args.Momentum+'_'+InputSample+'_'+str(i_a)+'_'+str(i_b)
          job_submit = '''jobsub_submit -G dune --role=Analysis --resource-provides=\"usage_model=DEDICATED,OPPORTUNISTIC\" -l \'+SingularityImage=\\\"/cvmfs/singularity.opensciencegrid.org/fermilab/fnal-wn-sl7:latest\\\"\' --append_condor_requirements=\'(TARGET.HAS_SINGULARITY=?=true)\' --tar_file_name "dropbox:///{0}/PDSPAna_build.tar\" --email-to sungbin.oh555@gmail.com -N 1 --expected-lifetime 8h --memory 6GB "file://{0}/{1}.sh\"'''.format(base_rundir, commandsfilename)
          os.system(job_submit)
          print(job_submit)
          os.chdir(cwd)

if args.no_exec:
  exit()




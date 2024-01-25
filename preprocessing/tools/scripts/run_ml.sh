#!/bin/bash
#echo "Starting job on " `date` #Date/time of start of job
#echo "Running on: `uname -a`" #Condor job is running on this node
#echo "System software: `cat /etc/redhat-release`" #Operating System on that node

source /cvmfs/cms.cern.ch/cmsset_default.sh 
export SCRAM_ARCH=slc7_amd64_gcc10

mldir=$1
middir=$2
infiles=$3
outfile=$4
maxEvents=$5
isMC=$6

filename=$(basename "$outpath")

# Root can't write directly to Hadoop so we use an intermediate directory
midpath=$middir/$filename

echo "Running NanoAOD production for file $filename"
#echo "Inpath: $inpath"
#echo "midpath: $midpath"
#echo "Outpath: $outpath"


cd $mldir
eval `scramv1 runtime -sh` # cmsenv is an alias not on the workers

cd ml_photons_cmssw/ml_photons/python/ 
cmsRun Prod_FlatAOD.py inputFiles=$infiles outputFile=$midpath maxEvents=$maxEvents isMC=$isMC

echo Finished!
echo Moving Output file: $outpath
cp $midpath $outpath
rm $midpath

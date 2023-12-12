#!/bin/bash
#echo "Starting job on " `date` #Date/time of start of job
#echo "Running on: `uname -a`" #Condor job is running on this node
#echo "System software: `cat /etc/redhat-release`" #Operating System on that node

source /cvmfs/cms.cern.ch/cmsset_default.sh 
export SCRAM_ARCH=slc7_amd64_gcc10

inpath=$1
outpath=$2
maxEvents=$3

filename=$(basename "$outpath")
midpath=/scratch365/atownse2/tmp/${filename}

echo "Running NanoAOD production for file $filename"
#echo "Inpath: $inpath"
#echo "midpath: $midpath"
#echo "Outpath: $outpath"


cd /afs/crc.nd.edu/user/a/atownse2/Public/RSTriPhoton/preprocessing/processMiniAOD/CMSSW_10_6_19_patch2/src/
eval `scramv1 runtime -sh` # cmsenv is an alias not on the workers

cd ml_photons_cmssw/ml_photons/python/ 
cmsRun Prod_FlatAOD.py $inpath $midpath #$maxEvents

echo Finished!
echo Moving Output file: $outpath
cp $midpath $outpath
rm $midpath

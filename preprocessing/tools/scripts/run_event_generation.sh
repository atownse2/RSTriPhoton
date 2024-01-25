#Generates events for a given gridpack
source /cvmfs/cms.cern.ch/cmsset_default.sh 
export SCRAM_ARCH=slc7_amd64_gcc10

#### Arguments

releasedir=$1
tmpdir=$2
gridpath=$3
outpath=$4
year=$5
nevents=$6
saveAOD=$7

outputdir=$(dirname "$outpath")
if [ ! -d "$outputdir" ]
then
  echo "Output directory does not exist, exiting."
  exit 1
fi


outfilename=$(basename "$outpath")
jobname=${outfilename%.*}

echo "Running event generation for $jobname"
echo "Output file name: $outfilename"


wkdir=$tmpdir/$jobname
if [ -d $wkdir ]
then
  echo "Removing existing temporary directory."
  rm -rf $wkdir
  mkdir $wkdir
else
  echo "Creating temporary directory"
  mkdir $wkdir
fi

#### File I/O
wmLHEGENfilepath=${wkdir}/wmLHEGEN.root
SIMfilepath=${wkdir}/SIM.root
DIGIfilepath=${wkdir}/DIGIPremix.root
HLTfilepath=${wkdir}/HLT.root
RECOfilepath=${wkdir}/RECO.root
MINIfilepath=${wkdir}/MINI.root

GEN=$releasedir/wmLHEGEN__CMSSW_10_6_22/src
SIM=$releasedir/SIM__CMSSW_10_6_17_patch1/src
DIGI=$releasedir/DIGIPremix__CMSSW_10_6_17_patch1/src
HLT=$releasedir/HLT__CMSSW_10_2_16_UL/src
RECO=$releasedir/RECO__CMSSW_10_6_17_patch1/src
MINI=$releasedir/MiniAOD__CMSSW_10_6_17_patch1/src

echo "Starting GEN to MiniAOD steps"

echo "GEN step"
cd $GEN
eval `scramv1 runtime -sh`
cd $wkdir

cmsRun $GEN/wmLHEGEN_step.py ${gridpath} ${nevents} ${wmLHEGENfilepath}

# Check exit status
if [ $? -ne 0 ]
then
  echo "GEN step failed with exit code $?, exiting."
  exit 1
fi


echo "SIM step"
cd $SIM
eval `scramv1 runtime -sh`
cmsRun SIM_step.py $wmLHEGENfilepath $SIMfilepath
rm $wmLHEGENfilepath

echo "DIGI step"
cd $DIGI
eval `scramv1 runtime -sh`
cmsRun DIGIPremix_step.py $SIMfilepath $DIGIfilepath
rm $SIMfilepath

echo "HLT step"
cd $HLT
eval `scramv1 runtime -sh`
cmsRun HLT_step.py $DIGIfilepath $HLTfilepath
rm $DIGIfilepath

echo "RECO step"
cd $RECO
eval `scramv1 runtime -sh`
cmsRun RECO_step.py $HLTfilepath $RECOfilepath
rm $HLTfilepath

if [ $saveAOD == "True" ]
then
  # Replace MiniAOD in outpath with AOD
  outMinifilepath=${outpath//MiniAOD/AOD}
  cp $RECOfilepath $outMinifilepath
fi

echo "MiniAOD step"
cd $MINI
eval `scramv1 runtime -sh` 
cmsRun MiniAOD_step.py $RECOfilepath $MINIfilepath
rm $RECOfilepath


echo "Done generating events for $jobname"
echo "Moving output file to hadoop"
mv $MINIfilepath $outpath

echo "Cleaning up, removing temporary directory"
rm -rf $wkdir

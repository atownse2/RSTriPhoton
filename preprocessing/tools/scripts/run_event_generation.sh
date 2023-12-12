#Generates events for a given gridpack
source /cvmfs/cms.cern.ch/cmsset_default.sh 
export SCRAM_ARCH=slc7_amd64_gcc10

#### Arguments

gridpath=$1
outpath=$2
year=$3
nevents=$4
saveAOD=$5

outputdir=$(dirname "$outpath")
if [ ! -d "$outputdir" ]
then
  echo "Output directory does not exist, exiting."
  exit 1
fi


outfilename=$(basename "$outpath")
fragment=${outfilename%.*}

echo "Running event generation for fragment $fragment"
echo "Output file name: $outfilename"


tmpdir=/scratch365/atownse2/tmp/${fragment}

if [ -d $tmpdir ]
then
  echo "Removing existing temporary directory."
  rm -rf $tmpdir
  mkdir $tmpdir
else
  echo "Creating temporary directory"
  mkdir $tmpdir
fi

#### File I/O
wmLHEGENfilepath=${tmpdir}/wmLHEGEN.root
SIMfilepath=${tmpdir}/SIM.root
DIGIfilepath=${tmpdir}/DIGIPremix.root
HLTfilepath=${tmpdir}/HLT.root
RECOfilepath=${tmpdir}/RECO.root
MINIfilepath=${tmpdir}/MINI.root

#### Map era = 2018 to era_tag = RunIISummer20UL18
declare -A era_tags=(
  ["2018"]="RunIISummer20UL18"
  ["2017"]="RunIISummer20UL17"
  ["2016"]="RunIISummer20UL16"
  ["2016APV"]="RunIISummer20UL16APV"
)
era_tag=${era_tags[$year]}


#### Configuration
releasedir=/afs/crc.nd.edu/user/a/atownse2/Public/RSTriPhoton/preprocessing/tools/EXO-MCsampleProductions/FullSimulation/${era_tag} 

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
cd $tmpdir

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


echo "Done generating events for fragment $fragment"
echo "Moving output file to hadoop"
mv $MINIfilepath $outpath

echo "Cleaning up, removing temporary directory"
rm -rf $tmpdir

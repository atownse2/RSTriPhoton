#Script to set up gridpack generation scripts
M_BKK=$1
M_R=$2
MGDIR=$3
outpath=$4
condor=$5

cd $MGDIR

# Deactivate any existing conda environment
if [ "$condor" = "False" ]
then
  source deactivate
fi

# If source deactivate doesn't work, try:
# conda init
# exec bash
# conda deactivate

fragment="BkkToGRadionToGGG_M1-${M_BKK}_R0-${M_R}"

if [ -d $fragment ]
then
  echo "Removing existing fragment directories."
  rm -rf ${fragment}*
fi

#Create cards from templates
echo Creating cards
cardsdir=cards/production/2017/13TeV/BkkToGRadionToGGG
cd $cardsdir
python BkkToGRadionToGGG_M1_R0_gen_card.py $M_BKK $M_R

cd $MGDIR

#Make gridpacks
echo Making gridpacks
./gridpack_generation.sh ${fragment} ${cardsdir}/${fragment}
cp ${fragment}_slc7_amd64_gcc10_CMSSW_12_4_8_tarball.tar.xz $outpath

echo Cleaning up...
rm -rf ${fragment}*
rm -rf $cardsdir/${fragment}*

echo Done
import os
import math
import numpy as np
import argparse

import htcondor

import sys
sys.path.append('/afs/crc.nd.edu/user/a/atownse2/Public/RSTriPhoton/')
from analysis.utils.sample_info import mass_grid, get_signal_filetag

nEventsPerBatch = 1000

gridpackdir = lambda fragment: f'/hadoop/store/user/atownse2/RSTriPhoton/gridpacks/{fragment}_slc7_amd64_gcc10_CMSSW_12_4_8_tarball.tar.xz'
MiniAODdir = lambda era: f'/hadoop/store/user/atownse2/RSTriPhoton/signal/{era}/MiniAOD'
mgdir = '/afs/crc.nd.edu/user/a/atownse2/Public/RSTriPhoton/preprocessing/generateSignal/genproductions/bin/MadGraph5_aMCatNLO/'
condordir = '/scratch365/atownse2/condor'
testdir = '/afs/crc.nd.edu/user/a/atownse2/Public/RSTriPhoton/test/'

parser = argparse.ArgumentParser(description='Generate signal MiniAOD.')
parser.add_argument('--test', '-t', action='store_true', help='Run in test mode generating 10 events for a single point in the mass grid.')
parser.add_argument('--batch', '-b', action='store_true', help='Submit jobs to condor.')
parser.add_argument('--m_moe', nargs=2, type=float, metavar=('M_BKK', 'MOE'), help='Specify a single point in the mass grid to generate events for.')
parser.add_argument('--m_m' , nargs=2, type=float, metavar=('M_BKK', 'M_R'), help='Specify a single point in the mass grid to generate events for.')
parser.add_argument('--nEvents', '-n', type=int, default=10, help='Number of events to generate per mass point.')
parser.add_argument('--era', type=str, default='2018', help='Specify the era to use for the fragment. Default is 2018.')
parser.add_argument('--gridpack_only', '-g', action='store_true', help='Only generate gridpacks. Do not generate events.')
parser.add_argument('--saveAOD', action='store_true', help='Save AOD as well as MiniAOD')

args = parser.parse_args()

nEventsTotal = args.nEvents
era = args.era
saveAOD = 'True' if args.saveAOD else 'False'

if args.m_moe:
  M_BKK = args.m_moe[0]
  MOE = args.m_moe[1]
  Mass_BKK_R = [(M_BKK, np.round(MOE*(M_BKK/2)))]
elif args.m_m:
  Mass_BKK_R = [(args.m_m[0], args.m_m[1])]
else:
  Mass_BKK_R = mass_grid

outputdir = MiniAODdir(era)
test_tags=''
if args.test:
  Mass_BKK_R = Mass_BKK_R[:1]
  nEventsTotal = 10
  outputdir = testdir
  test_tags = '-t'

if args.batch:
  test_tags += '-b'


def submit_batch_events(gridpath, outpath, nevents):
  frag_tag = outpath.split('/')[-1].replace('.root','')
  make_events = htcondor.Submit({
      "executable": "run_event_generation.sh",
      "arguments": f"{gridpath} {outpath} {era} {nevents} {saveAOD}",
      "output": f"{condordir}/out/{frag_tag}.out",
      "error" : f"{condordir}/err/{frag_tag}.err",
      "log"   : f"{condordir}/log/{frag_tag}.log",              
      "request_cpus": "1",
      "request_memory": "128MB",
      "request_disk": "128MB",
  })
  print(f"Submitting job for fragment: {frag_tag} with {nevents} events")
  schedd = htcondor.Schedd()
  submit_result = schedd.submit(make_events)

def submit_local_events(gridpath, outpath, nevents):
  command = f'./run_event_generation.sh {gridpath} {outpath} {era} {nevents} {saveAOD}'
  os.system(command)


for M_BKK,M_R in Mass_BKK_R:
  fragment = get_signal_filetag((M_BKK, M_R))

  gridpath = gridpackdir(fragment)

  # If gridpach doesn't exist, make it (locally for now)
  if not os.path.exists(gridpath):
      print(f'Gridpack {gridpath} does not exist. Making it now.')
      if os.path.exists(f'{mgdir}/{fragment}*'):
        os.system(f'rm -rf {mgdir}/{fragment}*')
      run = f'./run_gridpack.sh {M_BKK} {M_R} {mgdir} {gridpath} > log/{fragment}_gridpack.log'
      os.system(run)
  
  if args.gridpack_only:
    continue

  #Find files in directory and append instead of overwriting
  existing_batches = [ int(fname.split('_')[-1].split('.')[0]) for fname in os.listdir(outputdir) if fragment in fname ]

  # Get lowest batch number that doesn't already exist
  nBatches = int(math.ceil(nEventsTotal/nEventsPerBatch))
  batch_nums = []
  i = 0
  while len(batch_nums) < nBatches:
    if i not in existing_batches:
      batch_nums.append(i)
    i += 1

  nLastJob = nEventsTotal - (nBatches-1)*nEventsPerBatch
  for batch_num in batch_nums:
    nevents = nEventsPerBatch if batch_num != batch_nums[-1] else nLastJob
    if args.test:
      outpath = f'{outputdir}/{fragment}_{era}_MiniAOD_{batch_num}_{test_tags}.root'
    else:
      outpath = f'{outputdir}/{fragment}_{era}_MiniAOD_{batch_num}.root'

    if args.batch:
      submit_batch_events(gridpath, outpath, nevents)
    else:
      submit_local_events(gridpath, outpath, nevents)

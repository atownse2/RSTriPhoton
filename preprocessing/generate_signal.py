import sys
import os

top_dir = '/afs/crc.nd.edu/user/a/atownse2/Public/RSTriPhoton'
sys.path.append(top_dir)

from analysis.utils import sample_info as si

from utils import submit_jobs as sj
from utils import clean_up_files as clean

"""
Use cases:
- Generate up to 10000 events for every point in the mass grid
    python get_n_events_grid.py -n 10000
- Generate up to 1000 events for a single mass point
    python get_n_events_grid.py --m_m 500 50 -n 1000
    python get_n_events_grid.py --m_moe 500 0.1 -n 1000
"""

# Define paths
test_dir = f'{top_dir}/test'
gridpackdir = '/hadoop/store/user/atownse2/RSTriPhoton/gridpacks/'
mgdir = f'{top_dir}/preprocessing/tools/genproductions/bin/MadGraph5_aMCatNLO/'
run_gridpack = f"{top_dir}/preprocessing/tools/scripts/run_gridpack.sh"
run_event_generation = f"{top_dir}/preprocessing/tools/scripts/run_event_generation.sh"

# Event generation
n_events_per_file = 1000

# Gridpacks
def get_gridpack(M_BKK, M_R, remake=False):
    gridpath = gridpackdir + si.get_signal_filetag((M_BKK, M_R)) + '_slc7_amd64_gcc10_CMSSW_12_4_8_tarball.tar.xz'
    fragment = si.get_signal_filetag((M_BKK, M_R))
    if not os.path.exists(gridpath) or remake:
        if os.path.exists(f'{mgdir}/{fragment}*'):
            os.system(f'rm -rf {mgdir}/{fragment}*')
        os.system(f"{run_gridpack} {M_BKK} {M_R} {mgdir} {gridpath}")
    return gridpath

def generate_signal_point(signal_point,
                          n_events_total,
                          year="2018",
                          gridpack_only=False,
                          remake_gridpacks=False,
                          saveAOD="False",
                          batch=False,
                          test=False):

    outdir = si.get_tuple_dir('signal', 'MiniAOD', year=year)
    if test:
        outdir = test_dir
    
    fragment = si.get_signal_filetag(signal_point)
    gridpack_path = get_gridpack(*signal_point, remake=remake_gridpacks)

    if gridpack_only:
        return

    existing_files = {fname: info for fname, info in clean.get_good_files('MiniAOD', 'signal').items() if fragment in fname}
    existing_batches = [fname.split('_')[-1].split('.')[0] for fname in existing_files]
    n_existing_events = sum([info['n_events'] for info in existing_files.values()])
    
    ibatch = 0
    n_to_generate = n_events_total - n_existing_events
    print(f"Generating {n_to_generate} events for {fragment}")
    while n_to_generate > 0:
        if ibatch in existing_batches:
            continue

        n = n_events_per_file
        if n_to_generate < n_events_per_file:
            n = n_to_generate

        outpath = f"{outdir}/{fragment}_{ibatch}.root"
        arguments = f"{gridpack_path} {outpath} {year} {n} {saveAOD}"

        print(arguments)

        if batch:
            result = sj.submit_batch(run_event_generation,
                                    arguments,
                                    f"{fragment}_{ibatch}")
        else:
            result = sj.submit_interactive(run_event_generation,
                                            arguments)

        n_to_generate -= n
        ibatch += 1

        if test:
            break

    return

if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser(description='Generate MiniAOD signal events')
    parser.add_argument('--n_total_events', '-n', type=int, default=10, help='Number of events to generate per mass point.')
    parser.add_argument('--m_m' , nargs=2, type=float, metavar=('M_BKK', 'M_R'), help='Specify a single point in the mass grid to generate events for.')
    parser.add_argument('--m_moe', nargs=2, type=float, metavar=('M_BKK', 'MOE'), help='Specify a single point in the mass grid to generate events for.')
    parser.add_argument('--year', type=str, default='2018', help='Specify the era to use for the fragment. Default is 2018.')
    parser.add_argument('--batch', '-b', action='store_true', help='Submit jobs to condor.')
    parser.add_argument('--gridpack_only', '-g', action='store_true', help='Only generate gridpacks. Do not generate events.')
    parser.add_argument('--remake_gridpacks', '-r', action='store_true', help='Remake gridpacks even if they already exist.')
    parser.add_argument('--saveAOD', action='store_true', help='Save AOD as well as MiniAOD')
    parser.add_argument('--test', '-t', action='store_true', help='Run in test mode generating 10 events for a single point in the mass grid.')

    args = parser.parse_args()

    mass_grid = si.get_mass_grid('new')
    if args.m_moe:
        m_m = si.m_moe_to_m_m(args.m_moe)
        mass_grid = [m_m]
    elif args.m_m:
        mass_grid = [args.m_m]

    # Clean up files
    clean.clean_up_files('signal', args.year, 'MiniAOD')

    for signal_point in mass_grid:
        generate_signal_point(signal_point,
                              args.n_total_events,
                              year=args.year,
                              gridpack_only=args.gridpack_only,
                              remake_gridpacks=args.remake_gridpacks,
                              saveAOD=args.saveAOD,
                              batch=args.batch,
                              test=args.test)
import sys
import os

top_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.append(top_dir)

from preprocessing.utils import signal_info as si
from preprocessing.utils import tuple_info as ti
from preprocessing.utils import submit_jobs as sj
from preprocessing.utils import clean_up_files as clean 

# Define paths
test_dir = f'{top_dir}/test/MiniAOD/'

# gridpackdir = '/hadoop/store/user/atownse2/RSTriPhoton/gridpacks/'
gridpackdir = '/hadoop/store/user/atownse2/RSTriPhoton/testgridpacks/'
mgdir = f'{top_dir}/preprocessing/tools/genproductions/bin/MadGraph5_aMCatNLO/'
run_gridpack = f"{top_dir}/preprocessing/tools/scripts/run_gridpack.sh"

# Temporary directory for input/output files (maybe not necessary...)
tmp_dir = f'/scratch365/{os.environ["USER"]}/tmp'
if not os.path.isdir(tmp_dir):
    os.makedirs(tmp_dir)

# For event simulation and reconstruction
era_tags = { "2018": "RunIISummer20UL18", "2017": "RunIISummer20UL17", "2016": "RunIISummer20UL16", "2016APV": "RunIISummer20UL16APV" }
simrecodir = lambda year : f'{top_dir}/preprocessing/tools/EXO-MCsampleProductions/FullSimulation/{era_tags[year]}'
run_event_generation = f"{top_dir}/preprocessing/tools/scripts/run_event_generation.sh"

# Gridpacks
def get_gridpack(signal_point, remake=False, batch=False):
    fragment = si.signal_point_tag(signal_point)
    gridpath = gridpackdir + fragment + '_slc7_amd64_gcc10_CMSSW_12_4_8_tarball.tar.xz'
    M_BKK = signal_point['M_BKK']
    M_R = signal_point['M_R']
    arguments = f"{M_BKK} {M_R} {mgdir} {gridpath} {batch}"
    if not os.path.exists(gridpath) or remake:
        if os.path.exists(f'{mgdir}/{fragment}*'):
            os.system(f'rm -rf {mgdir}/{fragment}*')
        if batch:
            result = sj.submit_condor(run_gridpack,
                                    arguments,
                                    f'{fragment}_gridpack')
            return None
        else:
            result = sj.submit_interactive(run_gridpack,
                                            arguments)
    return gridpath


def generate_signal_point(
        signal_point,
        year,
        n_events_total,
        n_events_per_file=1000,
        gridpack_only=False,
        remake_gridpacks=False,
        saveAOD="False",
        batch=False,
        test=False):

    outdir = ti.get_tuple_dir('signal', 'MiniAOD')
    if test:
        outdir = test_dir
        n_events_total = 10
    
    dataset = si.signal_dataset(signal_point, year)
    gridpack_path = get_gridpack(signal_point, remake=remake_gridpacks, batch=batch)

    if gridpack_only or gridpack_path is None:
        return

    if not test:
        cleaner = clean.FileCleaner()
        cleaner.clean_up_files('signal', 'MiniAOD', datasets=[dataset])

        MiniAOD_info = cleaner.good_files['signal'][dataset]['MiniAOD']

        existing_files = MiniAOD_info['files'].keys()
        existing_events = MiniAOD_info['n_events']
    else:
        existing_files = []
        existing_events = 0

    n_to_generate = n_events_total - existing_events
    ibatch = -1
    while n_to_generate > 0:
        ibatch += 1

        filename = f"{dataset}_MiniAOD_{ibatch}.root"
        if filename in existing_files:
            continue

        if n_to_generate < n_events_per_file:
            n = n_to_generate
        else:
            n = n_events_per_file
        
        outpath = f"{outdir}/{filename}"
        arguments = f"{simrecodir(year)} {tmp_dir} {gridpack_path} {outpath} {year} {n} {saveAOD}"

        if batch:
            result = sj.submit_condor(run_event_generation,
                                    arguments,
                                    f"{dataset}_{ibatch}")
        else:
            result = sj.submit_interactive(run_event_generation,
                                            arguments)

        n_to_generate -= n

        if test:
            break

    return

if __name__ == "__main__":

    import argparse

    parser = argparse.ArgumentParser(description='Generate MiniAOD signal events')
    parser.add_argument('--n_total_events', '-n', type=int, default=10, help='Number of events to generate per mass point.')
    parser.add_argument('--n_events_per_file', '-nf', type=int, default=1000, help='Number of events per file.')
    parser.add_argument('--m_m' , nargs=2, type=float, metavar=('M_BKK', 'M_R'), help='Specify a single point in the mass grid to generate events for.')
    parser.add_argument('--m_moe', nargs=2, type=float, metavar=('M_BKK', 'MOE'), help='Specify a single point in the mass grid to generate events for.')
    parser.add_argument('--mass_grid_version', '-mgv', type=str, default='new', help='Specify the mass grid version to use. Default is current.')
    parser.add_argument('--year', type=str, default='2018', help='Specify the era to use for the fragment. Default is 2018.')
    parser.add_argument('--batch', '-b', action='store_true', help='Submit jobs to condor.')
    parser.add_argument('--gridpack_only', '-g', action='store_true', help='Only generate gridpacks. Do not generate events.')
    parser.add_argument('--remake_gridpacks', '-r', action='store_true', help='Remake gridpacks even if they already exist.')
    parser.add_argument('--saveAOD', action='store_true', help='Save AOD as well as MiniAOD')
    parser.add_argument('--test', '-t', action='store_true', help='Run in test mode generating 10 events for a single point in the mass grid.')

    args = parser.parse_args()

    signal_point = None
    if args.m_moe:
        signal_point = si.m_moe_to_m_m({'M_BKK': args.m_moe[0], 'MOE': args.m_moe[1]})
    elif args.m_m:
        signal_point = {'M_BKK': args.m_m[0], 'M_R': args.m_m[1]}

    if signal_point is not None:
        mass_grid = [signal_point]
    else:
        mass_grid = si.get_mass_grid(args.mass_grid_version)

    for signal_point in mass_grid:
        generate_signal_point(signal_point,
                              args.year,
                              args.n_total_events,
                              n_events_per_file=args.n_events_per_file,
                              gridpack_only=args.gridpack_only,
                              remake_gridpacks=args.remake_gridpacks,
                              saveAOD=args.saveAOD,
                              batch=args.batch,
                              test=args.test)
        if args.test:
            break
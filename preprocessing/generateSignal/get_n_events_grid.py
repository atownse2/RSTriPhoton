import sys
import os

sys.path.append('/afs/crc.nd.edu/user/a/atownse2/Public/RSTriPhoton/')
from analysis.utils import sample_info as si

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('n_total', type=int, help='total number of events')

parser.add_argument('--year', '-y', type=str, default='2018', help='year')
parser.add_argument('--test', '-t', help='Just print commands', action='store_true', default=False)

args = parser.parse_args()
n_total = args.n_total
n_files = n_total/1000

# Count number of miniAOD files for each mass point in grid
miniAOD_dir = "/hadoop/store/user/atownse2/RSTriPhoton/signal/2018/MiniAOD/"

signal_tags = si.get_all_signal_filetags()

signal_dict = {}
for signal_tag in signal_tags:
    signal_dict[signal_tag] = [f for f in os.listdir(miniAOD_dir) if signal_tag in f]

# for signal_tag in signal_tags:
#     print(f"{signal_tag}: {len(signal_dict[signal_tag])} files")


for signal_tag, files in signal_dict.items():
    n_batches = n_files - len(files)
    signal_point = si.get_mass_point(signal_tag)
    if n_batches > 0:
        n_events = 1000*n_batches
        print(f"{signal_point}: {n_batches} batches missing, generating now...")
        command = f"python3 generateSignal.py --m_m {signal_point[0]} {signal_point[1]} -n {n_events} -b"
        if not args.test:
            os.system(command)

        
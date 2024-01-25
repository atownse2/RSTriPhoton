# Unpack gridpack in temporary directoy and read cross section from .log file
# import os
# import sys

# import json

# top_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
# sys.path.append(top_dir)

# import argparse

# parser = argparse.ArgumentParser()
# parser.add_argument('--test', '-t', action='store_true', help='Run on one gridpack')

# args = parser.parse_args()

# gridpack_dir = '/hadoop/store/user/atownse2/RSTriPhoton/gridpacks/'
# tmpdir = '/tmp/atownse2'

# # Clean up tmpdir
# os.system(f'rm -rf {tmpdir}/*')

# log_file = 'gridpack_generation.log'
# xs_file = f'{top}/analysis/metadata/json/signal_xs.json'

# # Load existing cross sections
# if os.path.exists(xs_file):
#     with open(xs_file, 'r') as f:
#         xs_dict = json.load(f)
# else:
#     xs_dict = {}

# for point in sample_info.mass_grid:

#     if point in xs_dict:
#         continue

#     filetag = sample_info.get_signal_filetag(point)

#     # Get gridpack name
#     gridpack = [f for f in os.listdir(gridpack_dir) if filetag in f][0]

#     # Make workspace
#     os.chdir(tmpdir)
#     os.system(f'mkdir {filetag}')
#     os.chdir(f'{tmpdir}/{filetag}')

#     # Unpack gridpack
#     os.system(f'tar -xf {gridpack_dir+gridpack}')

#     # Read cross section from .log file
#     with open(log_file, 'r') as f:
#         for line in f:
#             if 'Cross-section :' in line:
#                 split_line = line.split()
#                 xs = float(split_line[2])
#                 error = float(split_line[4])
#                 if xs > 0:
#                     print(line)
#                     print(f'{filetag}: {xs} +/- {error}')
#                     xs_dict[filetag] = {'xs': xs, 'error': error}
    
#     # Clean up
#     os.system(f'rm -rf {tmpdir}/{filetag}')

#     if args.test:
#         break


# # Save cross sections to yaml file
# with open(xs_file, 'w') as f:
#     json.dump(xs_dict, f, indent=4)


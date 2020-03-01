#!/usr/bin/env python3




import argparse
import h5py




parser = argparse.ArgumentParser()
parser.add_argument('filename',                              help="HDF5 file containing a run_config group")
parser.add_argument('--newlines', '-n', action='store_true', help="Print newlines between config items")
args = parser.parse_args()




decode = lambda v: v.decode('UTF-8') if isinstance(v, bytes) else v
h5f = h5py.File(args.filename, 'r')
cfg = h5f['run_config']
sep = '\n' if args.newlines else ' '



print(sep.join(['{}={}'.format(item, decode(cfg[item][()])) for item in cfg if cfg[item][()]]))

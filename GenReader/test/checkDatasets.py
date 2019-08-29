import sys
import os
import subprocess
import argparse

from EFTGenReader.GenReader.DatasetHelper import DatasetHelper

arg_parser = argparse.ArgumentParser(prog='checkDatasets.py',description='Utility script for displaying various info about known datasets')
arg_parser.add_argument('--list',action='store_true',help='Print out all known sample names from the JSON database and exit')
arg_parser.add_argument('--nfiles',metavar='N',default=0,type=int,help='print out the first %(metavar)s files for each dataset. -1 prints all files')
arg_parser.add_argument('datasets',metavar='NAME',nargs='*',help='specify multiple datasets to run over, one after the other')
args = arg_parser.parse_args()

GIT_REPO_DIR = subprocess.check_output(['git','rev-parse','--show-toplevel']).strip()

def main():
    ds_helper = DatasetHelper()
    json_fpath = os.path.join(GIT_REPO_DIR,'GenReader/data/JSON/datasets.json')
    ds_helper.load(json_fpath)
    if args.list:
        print "---Available Samples---"
        for name in sorted(ds_helper.list()):
            print "{name}".format(name=name)
        return

    if len(args.datasets) == 0:
        print "No datasets specified, exiting..."
        return

    for name in args.datasets:
        if not ds_helper.exists(name):
            print "Unknown dataset: {name}".format(name)
            continue
        ds_helper.dump(name)
        if args.nfiles:
            lst = ds_helper.getFiles(name)
            max_idx = len(lst) if args.nfiles < 0 else min(args.nfiles,len(lst))
            for idx,fn in enumerate(lst):
                if idx >= args.nfiles and args.nfiles > 0:
                    break
                print "[{0}/{1}] {fpath}".format(idx+1,max_idx,fpath=fn)

if __name__ == "__main__":
    main()
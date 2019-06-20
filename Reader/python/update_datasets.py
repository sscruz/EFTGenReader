import os
import subprocess
import argparse
from DatasetHelper import DatasetHelper

arg_parser = argparse.ArgumentParser(prog='update_datasets.py')
arg_parser.add_argument('-n','--dry-run',action='store_true',help='perform a trial run with no changes made')
arg_parser.add_argument('-l','--load',metavar='FILE',dest='inf',default='datasets.json',help='load datasets from %(metavar)s')
arg_parser.add_argument('-s','--save',metavar='FILE',dest='outf',default='test.json',help='save datasets to %(metavar)s')
arg_parser.add_argument('-d','--dir',metavar='DIR',default='Reader/data/JSON',help='load and save files from %(metavar)s; %(metavar)s should be relative to the top directory of the git repo')

args = arg_parser.parse_args()

GIT_REPO_DIR = subprocess.check_output(['git','rev-parse','--show-toplevel']).strip()
DATASET_FILE = os.path.join(GIT_REPO_DIR,args.dir,args.inf)
SAVE_FILE    = os.path.join(GIT_REPO_DIR,args.dir,args.outf)

# Various examples of how to add/update/remove samples
def main():
    ds_helper = DatasetHelper()
    ds_helper.load(DATASET_FILE)

    # Add a sample by specifying all the relevant info directly
    ds_helper.updateDataset('new_dataset',
        dataset='/my_dataset/somwhere/on_DAS',
        loc=None,
        on_das=True,
        is_eft=True,        # Note: Currently, none of our private EFT samples are on DAS, so this should really be false
        central_xsec=42.0
    )

    new_ds = {
        'dataset': '/some_other_dataset/located/on_DAS',
        'loc': None,
        'on_das': True,
        'is_eft': True,
        'central_xsec': 1337.0
    }

    # Add a new sample by passing a dictionary with all the relevant info
    ds_helper.updateDataset('another_dataset',**new_ds)

    # Modify a setting in an already existing sample
    ds_helper.updateDataset('central_tZq',is_eft=True)
    ds_helper.updateDataset('central_tZq',is_eft=False)

    # Remove a sample 
    ds.removeDataset('new_dataset')

    # Save the dataset to json file
    if not args.dry_run:
        ds_helper.save(SAVE_FILE)

if __name__ == "__main__":
    #main()
    pass
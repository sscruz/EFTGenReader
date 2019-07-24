import os
import subprocess
import argparse
from DatasetHelper import DatasetHelper

arg_parser = argparse.ArgumentParser(prog='update_datasets.py')
arg_parser.add_argument('-n','--dry-run',action='store_true',help='perform a trial run with no changes made')
arg_parser.add_argument('-l','--load',metavar='FILE',dest='inf',default='datasets.json',help='load datasets from %(metavar)s')
arg_parser.add_argument('-s','--save',metavar='FILE',dest='outf',default='test.json',help='save datasets to %(metavar)s')
arg_parser.add_argument('-d','--dir',metavar='DIR',default='GenReader/data/JSON',help='load/save files from/to %(metavar)s; %(metavar)s should be relative to the top directory of the git repo')
arg_parser.add_argument('--list',action='store_true',help='List all samples from the loaded JSON file')

args = arg_parser.parse_args()

GIT_REPO_DIR = subprocess.check_output(['git','rev-parse','--show-toplevel']).strip()
DATASET_FILE = os.path.join(GIT_REPO_DIR,args.dir,args.inf)
SAVE_FILE    = os.path.join(GIT_REPO_DIR,args.dir,args.outf)

ttZ_xsec = 0.244109
ttW_xsec = 0.341762
tZq_xsec = 0.07358
ttH_xsec = 0.5638

# Basic example of adding/modifying/removing datasets
def example1(ds_helper):
    # Add a sample by specifying all the relevant info directly
    ds_helper.newDataset('new_dataset',
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
    ds_helper.newDataset('another_dataset',**new_ds)

    # Modify a setting in an already existing sample
    ds_helper.updateDataset('central_tZq',is_eft=True)
    ds_helper.updateDataset('central_tZq',is_eft=False)

    # Remove a sample 
    ds.removeDataset('new_dataset')

# Example of adding multiple new datasets via loops
def example2(ds_helper):
    proc_lst = [
        #('tHq',0.0),
        #('tllq4f',tZq_xsec),
        #('tllq',tZq_xsec),
        #('ttH',ttH_xsec),
        #('ttllNoHiggs',ttZ_xsec),
        #('ttllNuNuNoHiggs',ttZ_xsec),
        #('ttlnu',ttW_xsec),
        ('ttHJet',ttH_xsec)
    ]
    tag_lst = [
        ('HanModelxqcut25',True),
        ('HanModelxqcut35',True),
        ('HanModelxqcut45',True)
    ]
    run_lst = ['run0']
    ver_lst = [
        ('v1','qCut19'),
        ('v2','qCut30'),
        ('v3','qCut45')
    ]

    for p,central_xsec in proc_lst:
        for ver,vtag in ver_lst:
            for gtag,is_eft in tag_lst:
                for run in run_lst:
                    hadoop_path = '/hadoop/store/user/awightma/postLHE_step/2019_04_19/ttHJet-xqcutStudies/{ver}/mAOD_step_{process}_{tag}_{run}'.format(ver=ver,process=p,tag=gtag,run=run)
                    new_name = '{process}_{grp_tag}-{ver_tag}'.format(process=p,grp_tag=gtag,ver_tag=vtag)
                    ds_helper.newDataset(new_name,force=False,
                        dataset='',
                        loc=hadoop_path,
                        is_eft=is_eft,
                        on_das=False,
                        central_xsec=central_xsec
                    )
                    print new_name
                    ds_helper.dump(new_name)

# Example of modifying content of existing datasets
def example3(ds_helper):
    for name in ds_helper.list():
        ds_helper.updateDataset(name,datatier='MINIAODSIM')
        if ds_helper.getData(name,'on_das'):
            das_loc = ds_helper.getData(name,'dataset')
            ds_helper.updateDataset(name,loc=das_loc)

# Explicitly add a single new dataset 'by hand'
def example4(ds_helper):
    ds_helper.newDataset('tllq4fNoSchanW_NoMatching-SM',
        dataset='',
        loc='/hadoop/store/user/awightma/postLHE_step/2019_04_19/tllq4f-NoDim6Diagrams-NoMatching/v1/mAOD_step_tllq4f_NoDim6NoSchanW_run0',
        datatier='MINIAODSIM',
        on_das=False,
        is_eft=False,
        central_xsec=tZq_xsec
    )


def main():
    ds_helper = DatasetHelper()
    ds_helper.load(DATASET_FILE)

    if args.list:
        print "---Available Samples---"
        for sample_name in sorted(ds_helper.list()):
            print "{name}".format(name=sample_name)
        return

    #example1(ds_helper)
    #example2(ds_helper)
    #example3(ds_helper)
    example4(ds_helper)

    # Save the dataset to json file
    if not args.dry_run:
        ds_helper.save(SAVE_FILE)

if __name__ == "__main__":
    main()
import datetime
import os
import subprocess

from EFTGenReader.GenReader.DatasetHelper import DatasetHelper

from lobster import cmssw
from lobster.core import AdvancedOptions, Category, Config, Dataset, StorageConfiguration, Workflow

# NOTE1: The cmsRun config likely doesn't correctly handle input samples which do not have any EFT reweighting
# NOTE2: The input samples must be in MAOD format

GIT_REPO_DIR = subprocess.check_output(['git','rev-parse','--show-toplevel']).strip()
timestamp_tag = datetime.datetime.now().strftime('%Y%m%d_%H%M')

# NOTE: All samples must be located somewhere under the directory specified by input_path
input_path = "/store/user/awightma/"

out_ver = "v1"
tag = 'NoTopLeptons-NoCuts'

master_label = 'EFT_LHE_{tstamp}'.format(tstamp=timestamp_tag)

RUN_MODE = 'testing'

#output_path = "/store/user/$USER/KinematicGenHists/{tag}/{ver}".format(tag=tag,ver=out_ver)
if RUN_MODE == 'testing':
    output_path  = "/store/user/$USER/tests/lobster_{tstamp}".format(tstamp=timestamp_tag)
    workdir_path = "/tmpscratch/users/$USER/tests/lobster_{tstamp}".format(tstamp=timestamp_tag)
    plotdir_path = "~/www/lobster/tests/lobster_{tstamp}".format(tstamp=timestamp_tag)
elif RUN_MODE == 'mg_studies':
    output_path  = "/store/user/$USER/KinematicGenHists/{tag}/{ver}".format(tag=tag,ver=out_ver)
    workdir_path = "/tmpscratch/users/$USER/KinematicGenHists/{tag}/{ver}".format(tag=tag,ver=out_ver)
    plotdir_path = "~/www/lobster/KinematicGenHists/{tag}/{ver}".format(tag=tag,ver=out_ver)

storage = StorageConfiguration(
    input=[
        "hdfs://eddie.crc.nd.edu:19000"  + input_path,
        "root://deepthought.crc.nd.edu/" + input_path,  # Note the extra slash after the hostname!
        "gsiftp://T3_US_NotreDame"       + input_path,
        "srm://T3_US_NotreDame"          + input_path,
    ],
    output=[
        "hdfs://eddie.crc.nd.edu:19000"  + output_path,
        "file:///hadoop"                 + output_path,
        # ND is not in the XrootD redirector, thus hardcode server.
        "root://deepthought.crc.nd.edu/" + output_path, # Note the extra slash after the hostname!
        "gsiftp://T3_US_NotreDame"       + output_path,
        "srm://T3_US_NotreDame"          + output_path,
    ],
)

processing = Category(
    name='processing',
    cores=1,
    memory=1200,
    disk=1000
    #mode='fixed'
)

wf = []

ds_helper = DatasetHelper()
ds_helper.load(os.path.join(GIT_REPO_DIR,"GenReader/data/JSON/datasets.json"))

width = 1
samples = [
    'ttHJet_HanModelxqcut25-qCut19',
    'ttHJet_HanModelxqcut25-qCut30',
    'ttHJet_HanModelxqcut25-qCut45',
]

for idx,sample_name in enumerate(samples):
    sample_loc = ds_helper.getData(sample_name,'loc')
    full_path = sample_loc.split("/hadoop")[1]
    rel_path = os.path.relpath(full_path,input_path)

    print "[{0:0>{w}}/{1:0>{w}}] Sample: {sample}".format(idx+1,len(samples),sample=sample_name,w=width)
    print "\tFullPath:  {path}".format(path=full_path)
    print "\tInputPath: {path}".format(path=input_path)
    print "\tRelPath:   {path}".format(path=rel_path)

    # The workflow label can't have any dashes (-) in it, so remove them
    safe_label_name = sample_name.replace('-','')
    output = Workflow(
        label='output_{label}'.format(label=safe_label_name),
        command='cmsRun lobsterized_EFTGenReader_cfg.py',
        cleanup_input=False,
        outputs=['output_tree.root'],
        merge_size='1.0G',  # This is set to a large value to make sure the final output is merged into a single file
        dataset=Dataset(
            files=rel_path,
            files_per_task=1,
            patterns=["*.root"]
        ),
        merge_command='hadd @outputfiles @inputfiles',
        category=processing
    )

    wf.extend([output])

config = Config(
    label=master_label,
    workdir=workdir_path,
    plotdir=plotdir_path,
    storage=storage,
    workflows=wf,
    advanced=AdvancedOptions(
        bad_exit_codes=[127, 160],
        log_level=1,
    )
)
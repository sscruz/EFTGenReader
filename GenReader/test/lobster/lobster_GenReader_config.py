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
input_path = "/store/user/"

out_ver = "v1"
tag = 'NoTopLeptons-NoCuts'

#master_label = 'EFT_LHE_{tstamp}'.format(tstamp=timestamp_tag)
master_label = 'EFT_T3_{tstamp}'.format(tstamp=timestamp_tag)

RUN_MODE = 'testing'
# RUN_MODE = 'mg_studies'

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
    mode='fixed',
    cores=1,
    memory=1200,
    disk=1000
)

wf = []

ds_helper = DatasetHelper()
ds_helper.load(os.path.join(GIT_REPO_DIR,"GenReader/data/JSON/datasets.json"))

width = 1
samples = [
    'tllq4f_SMNoSchanW'
]

das_mode = False
hadoop_mode = False

for idx,sample_name in enumerate(samples):
    if not ds_helper.exists(sample_name):
        continue
    if ds_helper.getData(sample_name,'on_das'):
        das_mode = True
    else:
        hadoop_mode = True

if das_mode and hadoop_mode:
    print "[ERROR] The list of samples contains some that are on /hadoop and others that are on DAS!"
    raise RuntimeError
elif das_mode:
    print "Mode: DAS"
elif hadoop_mode:
    print "Mode: Hadoop"
    storage.input = [
        "hdfs://eddie.crc.nd.edu:19000"  + input_path,
        "root://deepthought.crc.nd.edu/" + input_path,  # Note the extra slash after the hostname!
        "gsiftp://T3_US_NotreDame"       + input_path,
        "srm://T3_US_NotreDame"          + input_path,
    ]
else:
    print "[ERROR] Unknown mode. Did you include at least one valid sample?"
    raise RuntimeError

for idx,sample_name in enumerate(samples):
    if not ds_helper.exists(sample_name):
        print "[{0:0>{w}}/{1:0>{w}}] Skipping unknown sample: {sample}".format(idx+1,len(samples),sample=sample_name,w=width)
        continue
    print "[{0:0>{w}}/{1:0>{w}}] Sample: {sample}".format(idx+1,len(samples),sample=sample_name,w=width)

    sample_loc = ds_helper.getData(sample_name,'loc')
    if hadoop_mode:
        full_path = sample_loc.split("/hadoop")[1]
        rel_path = os.path.relpath(full_path,input_path)
        ds = Dataset(
            files=rel_path,
            files_per_task=5,
            patterns=["*.root"]
        )
        merge_size = '256M'    # EFT samples with many reweight points are O(25M)
        print "\tFullPath:  {path}".format(path=full_path)
        print "\tInputPath: {path}".format(path=input_path)
        print "\tRelPath:   {path}".format(path=rel_path)
    elif das_mode:
        ds = cmssw.Dataset(
            dataset=sample_loc,
            events_per_task=100000
        )
        merge_size = '512K'     # non-EFT sample sizes are O(40K)

    cms_cmd = ['cmsRun','lobsterized_EFTGenReader_cfg.py']
    cms_cmd.extend([
        'datatier={tier}'.format(tier=ds_helper.getData(sample_name,'datatier')),
        'minPtJet=30.0',
        'maxEtaJet=2.5'
    ])
    if not ds_helper.getData(sample_name,'is_eft'):
        cms_cmd.extend(['iseft=False'])

    print "\tCommand:   {cmd}".format(cmd=' '.join(cms_cmd))

    # The workflow label can't have any dashes (-) in it, so remove them
    safe_label_name = sample_name.replace('-','')
    output = Workflow(
        label='output_{label}'.format(label=safe_label_name),
        command=' '.join(cms_cmd),
        cleanup_input=False,
        outputs=['output_tree.root'],
        merge_size=merge_size,  # Note: Lobster takes a very long time trying to merge large numbers of small files for some reason
        dataset=ds,
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
        xrootd_servers=['ndcms.crc.nd.edu',
                       'cmsxrootd.fnal.gov',
                       'deepthought.crc.nd.edu']
    )
)

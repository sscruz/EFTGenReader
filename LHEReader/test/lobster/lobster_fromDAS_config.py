import datetime
import os
import subprocess

from EFTGenReader.GenReader.DatasetHelper import DatasetHelper

from lobster import cmssw
from lobster.core import AdvancedOptions, Category, Config, Dataset, StorageConfiguration, Workflow

GIT_REPO_DIR = subprocess.check_output(['git','rev-parse','--show-toplevel']).strip()
tstamp_tag = datetime.datetime.now().strftime('%Y%m%d_%H%M')

RUN_MODE = 'testing'
#RUN_MODE = 'mg_studies'

ver = "v1"
era = "fromJSON"#"2019_04_19"
tag = "ExampleTag"

master_label = 'EFT_T3_{tstamp}'.format(tstamp=tstamp_tag)
#master_label = 'EFT_LHE_{tstamp}'.format(tstamp=tstamp_tag)

if RUN_MODE == 'testing':
    output_path  = "/store/user/$USER/tests/lobster_{tstamp}".format(tstamp=tstamp_tag)
    workdir_path = "/tmpscratch/users/$USER/tests/lobster_{tstamp}".format(tstamp=tstamp_tag)
    plotdir_path = "~/www/lobster/tests/lobster_{tstamp}".format(tstamp=tstamp_tag)
elif RUN_MODE == 'mg_studies':
    output_path  = "/store/user/$USER/summaryTree_LHE/{era}/{tag}/{ver}".format(era=era,tag=tag,ver=ver)
    workdir_path = "/tmpscratch/users/$USER/summaryTree_LHE/{era}/{tag}/{ver}".format(era=era,tag=tag,ver=ver)
    plotdir_path = "~/www/lobster/summaryTree_LHE/{era}/{tag}/{ver}".format(era=era,tag=tag,ver=ver)
else:
    print "[ERROR] Unknown run setup, {setup}".format(setup=RUN_MODE)
    raise ValueError

storage = StorageConfiguration(
    output=[
        "hdfs://eddie.crc.nd.edu:19000"  + output_path,
        "file:///hadoop"                 + output_path,
        # ND is not in the XrootD redirector, thus hardcode server.
        "root://deepthought.crc.nd.edu/" + output_path, # Note the extra slash after the hostname!
        "gsiftp://T3_US_NotreDame"       + output_path,
        "srm://T3_US_NotreDame"          + output_path,
        "file:///hadoop"                 + output_path,
    ]
)

processing = Category(
    name='processing',
    cores=1,
    memory=1200,
    disk=1000
    mode='fixed'
)

wf = []

ds_helper = DatasetHelper()
ds_helper.load(os.path.join(GIT_REPO_DIR,"GenReader/data/JSON/datasets.json"))

width = 1
samples = [
    'central_tZq'
]

print "Generating workflows:"
for idx,sample_name in enumerate(samples):
    if not ds_helper.exists(sample_name):
        print "[{0:0>{w}}/{1:0>{w}}] Skipping unknown sample: {sample}".format(idx+1,len(samples),sample=sample_name,w=width)
        continue
    if not ds_helper.getData('on_das'):
        print "[{0:0>{w}}/{1:0>{w}}] Skipping non-DAS sample: {sample}".format(idx+1,len(samples),sample=sample_name,w=width)
        continue
    sample_loc = ds_helper.getData(sample_name,'loc')
    ds = cmssw.Dataset(
        dataset=sample_loc,
        events_per_task=30000
    )

    cms_cmd = ['cmsRun','EFTLHEReader_cfg.py']
    cms_cmd.extend(['datatier=MINIAODSIM'])

    print "\t[{n}/{tot}] mAOD Input: {dir}".format(n=idx+1,tot=len(samples),dir=sample_name)
    print "\tCommand: {cmd}".format(cmd=' '.join(cms_cmd))

    safe_label_name = sample_name.replace('-','')
    output = Workflow(
        label='output_{label}'.format(label=safe_label_name),
        command=' '.join(cms_cmd),
        cleanup_input=False,
        merge_size=-1,
        dataset=ds,
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
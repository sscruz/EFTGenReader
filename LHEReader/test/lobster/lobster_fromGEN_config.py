import datetime
import os

from lobster import cmssw
from lobster.core import AdvancedOptions, Category, Config, Dataset, StorageConfiguration, Workflow

# This lobster config runs over LHE tier EDM root files, produced using custom MadGraph gridpacks
#   It produces output trees similar to OSTwoLepAna.cc, but can run over LHE level EDM files. It stores
#   only the eftwgts and the original xsec wgt as well as basic run,LS,etc. info

timestamp_tag = datetime.datetime.now().strftime('%Y%m%d_%H%M')

username = "awightma"

RUN_SETUP = 'local'
RUN_SETUP = 'full_production'
RUN_SETUP = 'mg_studies'

input_version  = "v1"
output_version = "v1"
grp_tag  = "2019_04_19/ValidationHanModelPlusJet"
test_tag = "lobster_20180505_1440"      # If the input LHE files were also produced in 'local' running
prod_tag = "Round1/Batch1"

# Only run over gridpacks from specific processes/coeffs/runs (i.e. MG starting points)
process_whitelist = []
coeff_whitelist   = []
runs_whitelist    = []

#master_label = 'EFT_T3_%s' % (timestamp_tag)
master_label = 'EFT_LHE_%s' % (timestamp_tag)

if RUN_SETUP == 'local':
    # For quick generic lobster workflow testing
    input_path   = "/store/user/%s/tests/%s" % (username,test_tag)
    output_path  = "/store/user/$USER/tests/lobster_%s" % (timestamp_tag)
    workdir_path = "/tmpscratch/users/$USER/tests/lobster_%s" % (timestamp_tag)
    plotdir_path = "~/www/lobster/tests/lobster_%s" % (timestamp_tag)
elif RUN_SETUP == 'mg_studies':
    # For MadGraph test studies
    input_path   = "/store/user/%s/genOnly_step/%s/%s/" % (username,grp_tag,input_version)
    output_path  = "/store/user/$USER/summaryTree_LHE/%s-GEN/%s" % (grp_tag,output_version)
    workdir_path = "/tmpscratch/users/$USER/summaryTree_LHE/%s-GEN/%s" % (grp_tag,output_version)
    plotdir_path = "~/www/lobster/summaryTree_LHE/%s-GEN/%s" % (grp_tag,output_version)
elif RUN_SETUP == 'full_production':
    input_path   = "/store/user/%s/FullProduction/%s/genOnly_step/%s/" % (username,prod_tag,input_version)
    output_path  = "/store/user/$USER/summaryTree_LHE/FP/%s-GEN/%s" % (prod_tag,output_version)
    workdir_path = "/tmpscratch/users/$USER/summaryTree_LHE/FP/%s-GEN/%s" % (prod_tag,output_version)
    plotdir_path = "~/www/lobster/summaryTree_LHE/FP/%s-GEN/%s" % (prod_tag,output_version)
else:
    print "Unknown run setup, %s" % (RUN_SETUP)
    raise ValueError

input_path_full = "/hadoop" + input_path

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
    disk=1000,
    mode='fixed'
)

gen_dirs = []
for f in os.listdir(input_path_full):
    dir_path = os.path.join(input_path_full,f)
    if not os.path.isdir(dir_path):
        continue
    arr = f.split('_')
    if arr[0] != 'gen':
        continue
    elif len(os.listdir(dir_path)) == 0:
        print "[WARNING] Skipping empty directory, %s" % (f)
        continue
    p,c,r = arr[2],arr[3],arr[4]
    if len(process_whitelist) > 0 and not p in process_whitelist:
        continue
    elif len(coeff_whitelist) > 0 and not c in coeff_whitelist:
        continue
    elif len(runs_whitelist) > 0 and not r in runs_whitelist:
        continue
    gen_dirs.append(f)

wf = []

print "Generating workflows:"
for idx,gen_dir in enumerate(gen_dirs):
    arr = gen_dir.split('_')
    p,c,r = arr[2],arr[3],arr[4]
    print "\t[{n}/{tot}] GEN Input: {dir}".format(n=idx+1,tot=len(gen_dirs),dir=gen_dir)
    output = Workflow(
        label='output_{p}_{c}_{r}'.format(p=p,c=c,r=r),
        command='cmsRun EFTLHEReader_cfg.py',
        merge_size='1.0G',
        cleanup_input=False,
        dataset=Dataset(
            files=gen_dir,
            files_per_task=5,   # Remember that the GEN step already does 5-10 files per task
            patterns=["*.root"]
        ),
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
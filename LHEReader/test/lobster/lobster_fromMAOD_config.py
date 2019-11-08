import datetime
import os

from lobster import cmssw
from lobster.core import AdvancedOptions, Category, Config, Dataset, StorageConfiguration, Workflow

# This lobster config runs over LHE tier EDM root files, produced using custom MadGraph gridpacks
#   It produces output trees similar to OSTwoLepAna.cc, but can run over LHE level EDM files. It stores
#   only the eftwgts and the original xsec wgt as well as basic run,LS,etc. info

timestamp_tag = datetime.datetime.now().strftime('%Y%m%d_%H%M')

#username = "awightma"
username = "kmohrman"

#RUN_SETUP = 'local'
RUN_SETUP = 'full_production'
#RUN_SETUP = 'mg_studies'

input_version  = ""
output_version = "v1"

#grp_tag  = "2019_04_19/tllq4f_t-channelMatched_pythia-JetMax1_b2"
grp_tag  = ""
#out_tag  = "tllq4f_t-channelMatched_pythia-JetMax1_b2"
#out_tag  = "2019_08_14_addPtBranches/ttXJet_R5B1-HanV2ModelNOttggh-HanV4Model-Comp_analysisEtaCut"
#out_tag  = "2019_08_14_addPtBranches/ttX-ttXJet_R5B1-HanModel-0Jetvs1JetTests_analysisEtaCut"
out_tag  = "ttXjet"
test_tag = "lobster_20180505_1440"      # If the input LHE files were also produced in 'local' running
prod_tag = "Round6/Batch2"

# Only run over gridpacks from specific processes/coeffs/runs (i.e. MG starting points)
#process_whitelist = ["ttlnuJet","ttllNuNuJetNoHiggs","ttHJet"]
#process_whitelist = ["ttlnuJet","ttHJet"]
#process_whitelist = ["ttH","ttllNuNuNoHiggs","ttlnu","ttHJet","ttllNuNuJetNoHiggs","ttlnuJet"]
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
    input_path   = "/store/user/%s/postLHE_step/%s/%s/" % (username,grp_tag,input_version)
    output_path  = "/store/user/$USER/summaryTree_LHE/%s-mAOD/%s" % (out_tag,output_version)
    workdir_path = "/tmpscratch/users/$USER/summaryTree_LHE/%s-mAOD/%s" % (out_tag,output_version)
    plotdir_path = "~/www/lobster/summaryTree_LHE/%s-mAOD/%s" % (out_tag,output_version)
elif RUN_SETUP == 'full_production':
    input_path   = "/store/user/{user}/FullProduction/{prod}/postLHE_step/{ver}/".format(user=username,prod=prod_tag,ver=input_version)
    output_path  = "/store/user/$USER/summaryTree_LHE/FP/{prod}/{tag}-mAOD/{ver}".format(tag=out_tag,prod=prod_tag,ver=output_version)
    workdir_path = "/tmpscratch/users/$USER/summaryTree_LHE/FP/{prod}/{tag}-mAOD/{ver}".format(tag=out_tag,prod=prod_tag,ver=output_version)
    plotdir_path = "~/www/lobster/summaryTree_LHE/FP/{prod}/{tag}-mAOD/{ver}".format(tag=out_tag,prod=prod_tag,ver=output_version)
else:
    print "Unknown run setup, %s" % (RUN_SETUP)
    raise ValueError

input_path = "/store/user/"
input_path_full = "/hadoop" + input_path

dir_list = [
            #os.path.join(input_path_full,"awightma/FullProduction/Round5/Batch1/postLHE_step/v1"),
            #os.path.join(input_path_full,"kmohrman/postLHE_step/2019_04_19/ttXJetTests-HanV4Model-xqcut10qCut19/v3"),
            #os.path.join(input_path_full,"kmohrman/postLHE_step/2019_04_19/ttHJet-ttWJet-HanV2ModelNOttgghCheck-xqcut10qCut19/v1"),
            #os.path.join(input_path_full,"kmohrman/postLHE_step/2019_04_19/ttZJet-HanV2ModelNOttgghCheck-xqcut10qCut19/v1"),
            #os.path.join(input_path_full,"/hadoop/store/user/kmohrman/postLHE_step/2019_04_19/ttX-HanModel/v1"),
            #os.path.join(input_path_full,"kmohrman/FullProduction/Round6/Batch1/postLHE_step/v1"), # FP: R6B1 ttXjet
            os.path.join(input_path_full,"kmohrman/FullProduction/Round6/Batch2/postLHE_step/v1"), # FP: R6B2 ttXjet
        ]

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

maod_dirs = []
for path in dir_list:
    for f in os.listdir(path):
        if not os.path.isdir(path):
            continue
        arr = f.split('_')
        if arr[0] != 'mAOD':
            continue
        elif len(os.listdir(path)) == 0:
            print "[WARNING] Skipping empty directory, %s" % (f)
            continue
        p,c,r = arr[2],arr[3],arr[4]
        if len(process_whitelist) > 0 and not p in process_whitelist:
            continue
        elif len(coeff_whitelist) > 0 and not c in coeff_whitelist:
            continue
        elif len(runs_whitelist) > 0 and not r in runs_whitelist:
            continue
        relpath = os.path.relpath(path,input_path_full)
        maod_dirs.append(os.path.join(relpath,f))

wf = []

print "Generating workflows:"
for idx,maod_dir in enumerate(maod_dirs):
    #arr = maod_dir.split('_')
    head,tail = os.path.split(maod_dir)
    arr = tail.split('_')
    p,c,r = arr[2],arr[3],arr[4]

    cms_cmd = ['cmsRun','EFTLHEReader_cfg.py']
    cms_cmd.extend(['datatier=MINIAODSIM'])
    
    print "\t[{n}/{tot}] mAOD Input: {dir}".format(n=idx+1,tot=len(maod_dirs),dir=maod_dir)
    print "\tCommand: {cmd}".format(cmd=' '.join(cms_cmd))

    output = Workflow(
        label='output_{p}_{c}_{r}'.format(p=p,c=c,r=r),
        command=' '.join(cms_cmd),
        merge_size='1.0G',
        cleanup_input=False,
        dataset=Dataset(
            files=maod_dir,
            files_per_task=5,
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

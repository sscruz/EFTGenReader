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
#tag = 'ttX-ttXJet_HanV4_nJets3orMoreLep_lepPt10-jetPt30_lepEta-jetEta_genOnly'
tag = 'testing_nJets_plots_July27_2020/ttX_central_cutsComp-nJets_test_lepPt10-jetPt30_lepEta-jetEta'
#tag = 'tllq4fNoSchanWNoHiggs0p_7p5mil_checks_ptEtaCuts'

#master_label = 'EFT_LHE_{tstamp}'.format(tstamp=timestamp_tag)
master_label = 'EFT_T3_{tstamp}'.format(tstamp=timestamp_tag)
#master_label = 'EFT_other_T3_{tstamp}'.format(tstamp=timestamp_tag)

RUN_MODE = 'testing'
#RUN_MODE = 'mg_studies'

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
    #disk=1000
    #disk=2900
    disk=4000
)

wf = []

ds_helper = DatasetHelper()
ds_helper.load(os.path.join(GIT_REPO_DIR,"GenReader/data/JSON/datasets.json"))

width = 1
samples = [
    #'central_tZq',
    #'central_tZq_new_pmx_v2',
    #'central_ttH',
    #'central_ttW',
    #'central_ttZ',
    #'tllq4fMatchedNoSchanW_fromMAOD_nJetMax1',
    #'tllq4fMatchedNoSchanW_fromGEN_nJetMax1',
    #'tllq4fMatchedNoSchanW_fromGEN_nJetMax2',
    #'tllq4fMatchedNoSchanW_fromGEN_nJetMax2doublecheck',
    #'tllq4fMatchedNoSchanW_fromMAOD_nJetMax2',
    #'tllq4fMatchedNoHiggs_fromMAOD',
    #'ttHJet_HanV2Model-xqcut10qcut15',
    #'ttHJet_HanV2Model-xqcut10qcut19',
    #'ttHJet_HanV2Model-xqcut10qcut25',
    #'ttHJet_HanV4Model-xqcut10qcut15',
    'ttHJet_HanV4Model-xqcut10qcut19',
    #'ttHJet_HanV4Model-xqcut10qcut25',
    #'ttllNuNuJetNoHiggs_HanV4Model-xqcut10qcut15',
    'ttllNuNuJetNoHiggs_HanV4Model-xqcut10qcut19',
    #'ttllNuNuJetNoHiggs_HanV4Model-xqcut10qcut25',
    #'ttlnuJet_HanV4Model-xqcut10qcut15',
    'ttlnuJet_HanV4Model-xqcut10qcut19',
    #'ttlnuJet_HanV4Model-xqcut10qcut25',
    #'ttHJet_HanModel16DttllScanpointsxqcut10-qCut15',
    #'ttHJet_HanModel16DttllScanpointsxqcut10-qCut19',
    #'ttHJet_HanModel16DttllScanpointsxqcut10-qCut25',
    #'ttHJet_HanV4StartPtRun2-xqcut10qcut15',
    #'ttHJet_HanV4StartPtRun2-xqcut10qcut19',
    #'ttHJet_HanV4StartPtRun2-xqcut10qcut25',
    #'ttllNuNuJetNoHiggs_HanV4StartPtRun2-xqcut10qcut15',
    #'ttllNuNuJetNoHiggs_HanV4StartPtRun2-xqcut10qcut19',
    #'ttllNuNuJetNoHiggs_HanV4StartPtRun2-xqcut10qcut25',
    #'ttlnuJet_HanV4StartPtRun1-xqcut10qcut15',
    #'ttlnuJet_HanV4StartPtRun1-xqcut10qcut19',
    #'ttlnuJet_HanV4StartPtRun1-xqcut10qcut25',
    #'ttHJet_FP_R6B1',
    #'ttHJet_FP_R6B2',
    #'ttHJet_FP_R6B5',
    #'ttllNuNuJetNoHiggs_FP_R6B1',
    #'ttllNuNuJetNoHiggs_FP_R6B2',
    #'ttllNuNuJetNoHiggs_FP_R6B5',
    #'ttlnuJet_FP_R6B1',
    #'ttlnuJet_FP_R6B2',
    #'ttlnuJet_FP_R6B5',
    #'ttHJet_FP_R5B1',
    #'ttlnuJet_FP_R5B1',
    #'ttllNuNuJetNoHiggs_FP_R5B1',
    #'tllq4fNoSchanWNoHiggs0p_B1_fromMAOD', # Andrew's 500 k
    #'tllq4fNoSchanWNoHiggs0p_B2_fromMAOD', # Andrew's 2 mil
    #'tllq4fNoSchanWNoHiggs0p_B3_fromMAOD', # Andrew's 5 mil
    #'tllq4fNoSchanWNoHiggs0p_FP_R6B4',     # My 5 mil extra
    #'tHq4f_HanV4SMCheck-fromMAOD',                      # SM checks for tHq
    #'tHq4f_HanV4SMCheck-b2-fromMAOD',
    #'tllq4fNoSchanWNoHiggs0p_HanV4SMCheck-fromMAOD',    # SM checks for tllq
    #'tllq4fNoSchanWNoHiggs0p_HanV4SMCheck-b2-fromMAOD',
    #'ttHJet_HanV4SMCheck-fromMAOD',                     # dim6=0 check for ttH
    #'ttHJet_HanV4SMCheck-b2-fromMAOD',
    #'ttHJet_HanModelOrigSMCheck-fromMAOD',
    #'ttlnuJet_HanV4SMCheck-fromMAOD',                   # dim6=0 check for ttW
    #'ttlnuJet_HanV4SMCheck-b2-fromMAOD',
    #'ttlnuJet_HanModelOrigSMCheck-fromMAOD',
    #'ttllNuNuJetNoHiggs_HanV4SMCheck-fromMAOD',         # dim6=0 check for ttZ
    #'ttllNuNuJetNoHiggs_HanV4SMCheck-b2-fromMAOD',
    #'ttllNuNuJetNoHiggs_HanModelOrigSMCheck-fromMAOD',
    #'ttHJet_HanV4_startPtChecksRun0',                   # Start point checks
    #'ttHJet_HanV4_startPtChecksRun1',
    #'ttHJet_HanV4_startPtChecksRun2',
    #'ttHJet_HanV4_startPtChecksRun3',
    #'ttllNuNuJetNoHiggs_HanV4_startPtChecksRun2',
    #'ttllNuNuJetNoHiggs_HanV4_startPtChecksRun3',
    #'ttlnuJet_HanV4_startPtChecksRun0',
    #'ttlnuJet_HanV4_startPtChecksRun1',
    #'ttlnuJet_HanV4_startPtChecksRun2',
    #'ttlnuJet_HanV4_startPtChecksRun3',
    #'tllq4fNoSchanWNoHiggs0p_HanV4_startPtChecksRun0',
    #'tllq4fNoSchanWNoHiggs0p_HanV4_startPtChecksRun1',
    #'tllq4fNoSchanWNoHiggs0p_HanV4_startPtChecksRun2',
    #'tllq4fNoSchanWNoHiggs0p_HanV4_startPtChecksRun3',
    #'ttH_HanV4ModelNoJets',             # HanV4 no jet samples
    #'ttlnu_HanV4ModelNoJets',           # HanV4 no jet samples
    #'ttllNuNuNoHiggs_HanV4ModelNoJets', # HanV4 no jet samples
    #'ttH_SM-SMEFTNLO-fromReza', # Reza's SM NLO sampels
    #'ttW_SM-SMEFTNLO-fromReza', # Reza's SM NLO sampels
    #'ttZ_SM-SMEFTNLO-fromReza', # Reza's SM NLO sampels
    #'ttWJet_HanV4', # (Private ttWJet NOT ttlnu!)
    #'ttZJet_HanV4', # (Private ttZJet NOT ttll!)
    #'ttW_HanV4',    # (Private ttW NOT ttlnu!)
    #'ttZ_HanV4',    # (Private ttZ NOT ttll!)
    'tHq4f_HanV4_GEN',
    'tllq4fNoSchanWNoHiggs0p_HanV4_GEN',
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
    is_eft = ds_helper.getData(sample_name,'is_eft')
    if hadoop_mode:
        full_path = sample_loc.split("/hadoop")[1]
        rel_path = os.path.relpath(full_path,input_path)
        ds = Dataset(
            files=rel_path,
            #files_per_task=5,
            files_per_task=ds_helper.getData(sample_name,'files_per_task'),
            patterns=["*.root"]
        )
        if is_eft:
            merge_size = '256M'     # EFT samples with many reweight points are O(25M)
        else:
            merge_size = '512K'     # non-EFT samples are O(50-100k)
        print "\tFullPath:  {path}".format(path=full_path)
        print "\tInputPath: {path}".format(path=input_path)
        print "\tRelPath:   {path}".format(path=rel_path)
    elif das_mode:
        ds = cmssw.Dataset(
            dataset=sample_loc,
            #events_per_task=100000
            events_per_task=300000
        )
        merge_size = '512K'     # non-EFT sample sizes are O(40K)
        #merge_size = -1

    #cms_cmd = ['cmsRun','lobsterized_EFTGenReader_cfg.py'] # Uncomment this if want to run EFTGenReader not EFTGenHistsWithCuts
    cms_cmd = ['cmsRun','lobsterized_EFTGenHistsWithCuts_cfg.py'] # Uncomment this if want to run EFTGenHistsWithCuts not EFTGenReader
    cms_cmd.extend([
        'datatier={tier}'.format(tier=ds_helper.getData(sample_name,'datatier')),
        'minPtLep=10',
        'minPtJet=30.0',
        'maxEtaJet=2.5',
        'maxEtaLep=2.5'
    ])
    if not is_eft:
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
        dashboard=False,
        bad_exit_codes=[127, 160],
        log_level=1,
        #xrootd_servers=['ndcms.crc.nd.edu',
        #               'cmsxrootd.fnal.gov',
        #               'deepthought.crc.nd.edu']
    )
)

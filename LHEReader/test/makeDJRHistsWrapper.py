import datetime
import os
import subprocess
from EFTGenReader.GenReader.make_html import make_html
from EFTGenReader.GenReader.utils import clean_dir,get_files,move_files

USE_SMPT = True; USE_REFPT = False  # SM base pt
#USE_SMPT = False; USE_REFPT = True # Ref base pt

def make_plots(sample_loc,output_fn):
    if not os.path.exists(sample_loc):
        return

    files = []

    print 'Loc: {loc}'.format(loc=sample_loc)
    print 'Files:'
    for fn in os.listdir(sample_loc):
        fpath = os.path.join(sample_loc,fn)
        if not os.path.isfile(fpath):
            continue
        if fn.split('.')[-1] != 'root':
            continue
        files.append(fpath)
        #print "\t{fn}".format(fn=fn)

    inputs_fn = 'infpaths.txt'
    with open(inputs_fn,'w') as f:
        for fpath in files:
            l = "{fpath}\n".format(fpath=fpath)
            f.write(l)

    macro_str = "makeDJRHists.C(\"{infs}\",\"{outf}\",{smpt},{refpt})".format(infs=inputs_fn,outf=output_fn,smpt=int(USE_SMPT),refpt=int(USE_REFPT))
    cmd = ['root','-b','-l','-q',macro_str]

    print 'Root Command: {0}'.format(' '.join(cmd))
    subprocess.check_call(['root','-b','-l','-q',macro_str])
    print ""


def main():

    ### Directories for pheno paper ###

    # Nominal ttH, ttW, ttZ
    dir_loc_procs = '/hadoop/store/user/kmohrman/summaryTree_LHE/2020_03_03_addPSweights/ttX-ttXJet_HanV4_QED1-and-noConstraints_moreStats-goodStartPt-GEN/v1'
    dir_list_procs = [
        ("output_ttHJet_HanV4lModel16DttllScanpoints_run1" , "ttHJet_HanV4"),
        ("output_ttWJet_HanV4goodStartPt_run1" , "ttWJet_HanV4"),
        ("output_ttZJet_HanV4goodStartPt_run2" , "ttZJet_HanV4"),
    ]
    # ttH q cut variations
    dir_loc_ttHqCutScan = '/hadoop/store/user/kmohrman/summaryTree_LHE/2019_08_14_addPtBranches/ttXJet-HanV4Model-xqcut10qCutTests_analysisEtaCut-GEN/v1'
    dir_list_ttHqCutScan = [
        ("output_ttHJet_HanV4lModel16DttllScanpointsqCut15_run1", "ttHJet_xqcut10qCut15"),
        ("output_ttHJet_HanV4lModel16DttllScanpointsqCut19_run1", "ttHJet_xqcut10qCut19"),
        ("output_ttHJet_HanV4lModel16DttllScanpointsqCut25_run1", "ttHJet_xqcut10qCut25"),
    ]

    # ttHJet: xqCut=10 qCut Scan samples
    dir_loc_ttHxqcutScan = '/hadoop/store/user/kmohrman/summaryTree_LHE/2019_08_14_addPtBranches/ttHJet_HanV4xqcutTests_analysisEtaCut-GEN/v1'
    dir_list_ttHxqcutScan = [
        ("output_ttHJet_HanV4ttXjetxqcut5qCut19_run1" , "ttHJet_HanV4_xqcut05qCut19"),
        ("output_ttHJet_HanV4ttXjetxqcut10qCut19_run1", "ttHJet_HanV4_xqcut10qCut19"),
        ("output_ttHJet_HanV4ttXjetxqcut15qCut19_run1", "ttHJet_HanV4_xqcut15qCut19"),
    ]

    # ttH missing 5 particle vertex
    dir_loc_ttHmissing5PV = '/hadoop/store/user/kmohrman/summaryTree_LHE/2019_08_14_addPtBranches/ttHJet-xqcutStudies-xqcut10qCutTests_analysisEtaCut-GEN/v3/'
    dir_list_ttHmissing5PV = [
        ("output_ttHJet_HanModel16DttllScanpointsqCut19_run1" , "ttH_xqcut10qCut19_missing5PV"),
    ]

    #dir_loc , dir_list = dir_loc_procs, dir_list_procs
    #dir_loc , dir_list = dir_loc_ttHqCutScan , dir_list_ttHqCutScan
    #dir_loc , dir_list = dir_loc_ttHxqcutScan , dir_list_ttHxqcutScan
    dir_loc , dir_list = dir_loc_ttHmissing5PV , dir_list_ttHmissing5PV

    ##################################

    ## Han model orig, Han model v2 w/o 4 particel verts, Han model v4 comp
    #dir_loc = '/hadoop/store/user/kmohrman/summaryTree_LHE/2019_08_14_addPtBranches/ttHJet-ttWJet_R5B1-HanV4Model-Comp_analysisEtaCut-mAOD/v1'
    #dir_list = [
    #    ("output_ttHJet_HanModel16DttllScanpoints_run1",               "ttHJet_Original_Model"),
    #    ("output_ttHJet_HanV2ModelNOttggh16DttllScanpointsqCut19_run1","ttHJet_HanV2_NOttggh"),
    #    ("output_ttHJet_HanV4lModel16DttllScanpointsqCut19_run1",      "ttHJet_HanV4"),
    #]

    # Han model v4
    #dir_loc = '/hadoop/store/user/kmohrman/summaryTree_LHE/2019_08_14_addPtBranches/ttXJet-HanV4Model-xqcut10qCutTests_analysisEtaCut-GEN/v1'
    #dir_list = [
    #    ("output_ttHJet_HanV4lModel16DttllScanpointsqCut15_run1", "ttHJet_xqcut10qCut15"),
    #    ("output_ttHJet_HanV4lModel16DttllScanpointsqCut19_run1", "ttHJet_xqcut10qCut19"),
    #    ("output_ttHJet_HanV4lModel16DttllScanpointsqCut25_run1", "ttHJet_xqcut10qCut25"),
    #    ("output_ttllNuNuJetNoHiggs_HanV4lModel16DttllScanpointsqCut15_run1", "ttZJet_xqcut10qCut15"),
    #    ("output_ttllNuNuJetNoHiggs_HanV4lModel16DttllScanpointsqCut19_run1", "ttZJet_xqcut10qCut19"),
    #    ("output_ttllNuNuJetNoHiggs_HanV4lModel16DttllScanpointsqCut25_run1", "ttZJet_xqcut10qCut25"),
    #    ("output_ttlnuJet_HanV4lModel16DttllScanpointsqCut15_run1", "ttWJet_xqcut10qCut15"),
    #    ("output_ttlnuJet_HanV4lModel16DttllScanpointsqCut19_run1", "ttWJet_xqcut10qCut19"),
    #    ("output_ttlnuJet_HanV4lModel16DttllScanpointsqCut25_run1", "ttWJet_xqcut10qCut25"),
    #]

    # Han model v2
    #dir_loc = '/hadoop/store/user/kmohrman/summaryTree_LHE/2019_08_14_addPtBranches/ttHJetTests-HanV2Model-xqcut10qCut19_analysisEtaCut-GEN/v1'
    #dir_list = [
        #("output_ttHJet_HanV2Model16DttllScanpointsqCut19_run1"   , "ttHJet"),
        #("output_ttHJetgg_HanV2Model16DttllScanpointsqCut19_run1" , "ttHJetgg"),
        #("output_ttHJetgq_HanV2Model16DttllScanpointsqCut19_run1" , "ttHJetgq"),
        #("output_ttHJetqq_HanV2Model16DttllScanpointsqCut19_run1" , "ttHJetqq"),
    #]

    # New updated dim 6 UFO model (without the -1 QED couplings)
    #dir_loc = '/hadoop/store/user/kmohrman/summaryTree_LHE/2019_08_14_addPtBranches/ttHJet-NewHanModelxqcut10qcut19_analysisEtaCut-GEN/v1'
    #dir_list = [("output_ttHJet_NewHanModel16DttllScanpointsqCut19_run1" , "ttHJet_NewModel")]

    # ttbarJet sub process tests:
    #dir_loc ='/hadoop/store/user/kmohrman/summaryTree_LHE/2019_08_14_addPtBranches/ttHJetSubProcTestctGEffectsStudies-xqcut10_analysisEtaCut-GEN/v1'
    #dir_list = [
    #    ("output_ttH1Jetgq_HanModel16DttllScanpointsAdamsReqxqcut10qCut19_run1" , "ttHJetgq"),
    #    ("output_ttHJetgg_HanModel16DttllScanpointsAdamsReqxqcut10qCut19_run1" , "ttHJetgg"),
    #    #("output_ttHJetqq_HanModel16DttllScanpointsAdamsReqxqcut10qCut19_run1" , "ttHJetqq"),
    #]

    #ttbarJetgg and ttHJetqq xqcut=10-qCut=19 (Adam's request)
    #dir_loc = '/hadoop/store/user/kmohrman/summaryTree_LHE/2019_08_14_addPtBranches/HanModel16DttllScanpointsAdamsReqxqcut10_analysisEtaCut-GEN/v3'
    #dir_list = [
        #("output_ttbarJetgg_HanModel16DttllScanpointsAdamsReqxqcut10qCut19_run1" , "ttbarJetgg"),
        #("output_ttHJetqq_HanModel16DttllScanpointsAdamsReqxqcut10qCut19_run1"   , "ttHJetqq")
    #]

    #dir_loc = '/hadoop/store/user/kmohrman/summaryTree_LHE/tllq4f_t-channelMatched_pythia-JetMax1_b2-mAOD/v1'
    #dir_list = [("output_tllq4fMatchedNoSchanW_HanModel16DttllScanpointsXQCUT10_run1","tllq4fMatchedNoSchanW")]


    # tllq4f_t-channelMatched_pythia-JetMax1_b2-mAOD: nJetMax Tests 
    ##dir_loc = '/hadoop/store/user/kmohrman/summaryTree_LHE/tllq4f_t-channelMatched_pythia-JetMax1_b2-mAOD/v1' # For max n jet 1 (From Andrew's MAOD file)
    #dir_loc = '/hadoop/store/user/kmohrman/summaryTree_LHE/2019_08_14_addPtBranches/tllq4fMatchedNoSchanW_xqcut10qCut19-nJetMax2_analysisEtaCut-GEN/v1' # From my gen only file
    #dir_loc = '/hadoop/store/user/kmohrman/summaryTree_LHE/2019_08_14_addPtBranches/tllq4fMatchedNoSchanW_xqcut10qCut19-nJetMaxTests_analysisEtaCut-GEN/v1'
    #dir_list = [
        ##("output_tllq4fMatchedNoSchanW_HanModel16DttllScanpointsXQCUT10_run1","tllq4fMatchedNoSchanW") # Goes with first dir
        #("output_tllq4fMatchedNoSchanW_HanModel16DttllScanpointsXQCUT10qCut19nJetMax2_run1","tllq4fMatchedNoSchanW_nJetMax2") # Goes with second dir
        #("output_tllq4fMatchedNoSchanW_HanModel16DttllScanpointsXQCUT10qCut19nJetMax1_run1"       , "nJetMax1"),
        #("output_tllq4fMatchedNoSchanW_HanModel16DttllScanpointsXQCUT10qCut19nJetMax2_run1"       , "nJetMax2"),
        #("output_tllq4fMatchedNoSchanW_HanModel16DttllScanpointsXQCUT10qCut19nJetMaxDefault_run1" , "nJetMaxDefault"),
    #]

    # Central tZq (does not have DJR values)
    #dir_loc = '/hadoop/store/user/kmohrman/summaryTree_LHE/fromJSON/central_tZq/v1'
    #dir_list = [('output_central_tZq' , 'central_tZq')]


    # ttHJet: xqCut Scan qCut=25 samples (Reza's request)
    #dir_loc = '/hadoop/store/user/kmohrman/summaryTree_LHE/2019_08_14_addPtBranches/ttHJet-xqcutStudies-xqcutScan_analysisEtaCut-GEN/v1'
    #dir_list = [
    #    ("output_ttHJet_HanModel16DttllScanpointsqCut19_run1"        , "xqcut10qCut19"), # Include 10-19 for comparison
    #    ("output_ttHJet_HanModel16DttllScanpointsxqcut05qCut10_run1" , "xqcut05qCut10"), # Include 05-10 since the only other 05-10 we have has starting point of 4
    #    ("output_ttHJet_HanModel16DttllScanpointsxqcut05qCut25_run1" , "xqcut05qCut25"), # xqcutScan with qCut=25 
    #    ("output_ttHJet_HanModel16DttllScanpointsqCut25_run1"        , "xqcut10qCut25"), # xqcutScan with qCut=25
    #    ("output_ttHJet_HanModel16DttllScanpointsxqcut15qCut25_run1" , "xqcut15qCut25"), # xqcutScan with qCut=25
    #    ("output_ttHJet_HanModel16DttllScanpointsxqcut20qCut25_run1" , "xqcut20qCut25"), # xqcutScan with qCut=25
    #]

    # ttHJet: xqCut=10 qCut Scan samples
    #dir_loc = '/hadoop/store/user/kmohrman/summaryTree_LHE/2019_08_14_addPtBranches/ttHJet-xqcutStudies-xqcut10qCutTests_analysisEtaCut-GEN/v3/'
    #dir_list = [
    #    ("output_ttHJet_HanModel16DttllScanpointsqCut15_run1" , "xqcut10qCut15"),
    #    ("output_ttHJet_HanModel16DttllScanpointsqCut19_run1" , "xqcut10qCut19"),
    #    ("output_ttHJet_HanModel16DttllScanpointsqCut25_run1" , "xqcut10qCut25"),
    #]

    # Full production samples (all processes)
    #dir_loc = '/hadoop/store/user/kmohrman/summaryTree_LHE/FP/Round5/Batch1/allProcesses-mAOD/v1/'
    #dir_list = [
        #("output_tHq4fMatched_HanModel16DttllScanpoints_run1"         , "tHq4fMatched"),
        #("output_tllq4fMatchedNoHiggs_HanModel16DttllScanpoints_run1" , "tllq4fMatchedNoHiggs"),
        #("output_ttHJet_HanModel16DttllScanpoints_run1"               , "ttHJet"),
        #("output_ttllNuNuJetNoHiggs_HanModel16DttllScanpoints_run1"   , "ttllNuNuJetNoHiggs"),
        #("output_ttlnuJet_HanModel16DttllScanpoints_run1"             , "ttlnuJet"),
    #]

    for dir_name,info in dir_list:
        print dir_name, info
        full_path = os.path.join(dir_loc,dir_name)
        print "Full path: " , full_path
        #outf = 'FP_{info}_'.format(info=info) # Base pt and WC value are appended in .C file
        outf = '{info}_'.format(info=info) # Base pt and WC value are appended in .C file
        print "Name out output file:" , outf
        make_plots(full_path,outf)

        #return # Return here if just finding xsec
        #print "TEST, this should not print !!!"

        output_dir = ''
        #output_dir = '/afs/crc.nd.edu/user/k/kmohrman/www/EFT/DJRplots/WCtests/ctGscan/xqcut10qCut19'
        #output_dir = '/afs/crc.nd.edu/user/k/kmohrman/www/EFT/DJRplots/testing/FP_dir_test/{outf}WCvalsAN_noctGscan_djr'.format(outf=outf)
        #output_dir = '/afs/crc.nd.edu/user/k/kmohrman/www/EFT/DJRplots/ModelFileTests/dim6top_LO_UFO_HanV4_2/noCuts/{outf}ctGScan250'.format(outf=outf)
        #output_dir = '/afs/crc.nd.edu/user/k/kmohrman/www/EFT/DJRplots/forPhenoFinalDraft/ttH_xqcut_scan/{outf}'.format(outf=outf)

        if output_dir == '':
            timestamp_tag = datetime.datetime.now().strftime('%Y%m%d_%H%M')
            output_dir = '/afs/crc.nd.edu/user/k/kmohrman/www/EFT/DJRplots/testing/{time}_djr'.format(time=timestamp_tag)

        if USE_SMPT:
            #output_dir = '/afs/crc.nd.edu/user/k/kmohrman/www/EFT/DJRplots/testing/test_WCscan_SMbase'
            #output_dir = '/afs/crc.nd.edu/user/k/kmohrman/www/EFT/DJRplots/RwgtTests/basePtSM_xqcut10qCut19'
            #output_dir = '/afs/crc.nd.edu/user/k/kmohrman/www/EFT/DJRplots/RwgtTests/PtEtaCutsTests/PtCuts_EtaCut2p5/3lep/basePtSM_xqcut10qCut19_PtEtaCuts2lep'
            #output_dir = '/afs/crc.nd.edu/user/k/kmohrman/www/EFT/DJRplots/FPSamples/Round5Batch1/2lep_PtEtaCuts/basePtSM/{outf}WCvalsAN_noctGscan_djr'.format(outf=outf)
            print "Output dir:" , output_dir
        else:
            #output_dir = '/afs/crc.nd.edu/user/k/kmohrman/www/EFT/DJRplots/RwgtTests/basePtRef_xqcut10qCut19'
            #output_dir = '/afs/crc.nd.edu/user/k/kmohrman/www/EFT/DJRplots/RwgtTests/PtEtaCutsTests/PtCuts_EtaCut2p5/3lep/basePtRef_xqcut10qCut19_PtEtaCuts2lep'
            #output_dir = '/afs/crc.nd.edu/user/k/kmohrman/www/EFT/DJRplots/FPSamples/Round5Batch1/2lep_PtEtaCuts/basePtRef/{outf}WCvalsAN_noctGscan_djr'.format(outf=outf)
            print "Output dir:" , output_dir
            

        if not os.path.exists(output_dir):
            os.mkdir(output_dir);

        #imgs = get_files('.',targets=["^ttHJet_.*\.png$"])
        imgs = get_files('.',targets=["^.*\.png$"])
        imgs = imgs + get_files('.',targets=["^.*\.pdf$"])
        move_files(files=imgs,target=output_dir)
        make_html(output_dir)

    #imgs = get_files('.',targets=["^ttHJet_.*\.png$"])
    #move_files(files=imgs,target=output_dir)
    #make_html(output_dir)

if __name__ == "__main__":
    main()


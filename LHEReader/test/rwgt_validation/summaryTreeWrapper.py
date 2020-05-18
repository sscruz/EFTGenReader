# Wrapper to 'readOutputTree.C', which specifies the root files to run over
# This code reads output tress produced by EFTLHEReader.cc and groups them by process+coeff, then
#   processes them in 'readOutputTree.C'

import subprocess
import os
from utils import regex_match, getInfoByTag, getDirectories, groupByProcess, groupByCoefficient

#['ttH']
#['ctWctZctpcpQMAxisScan','cpQ3cptcptbcbWctGAxisScan','cptbcbWctGAxisScan']

ALL_INFO = [
    {
        'tag': '2018_05_06/AllProcessesAllCoeffsAxisScans',
        'grp_name': '1DScan',#'1DScan',
        'version': 'v1',
        'include': False,
        'p_wl': ['tllq'],#['ttH','ttlnu','tllq'],
        'c_wl': ["cbW"],#["ctW","ctp","cpQM","ctZ","ctG","cbW","cpQ3","cptb","cpt","ctl1","cQe1"],
        'r_wl': [],
    },
    {
        'tag': '2018_05_06/ttbar1DAxisScans',
        'grp_name': '',
        'version': 'v1',
        'include': False,
        'p_wl': [],
        'c_wl': ['ctW'],#['ctp','cbW','ctZ','ctW','ctG','cQe1'],
        'r_wl': [],
    },
    {
        'tag': '2018_05_06/4DAxisScans',
        'grp_name': '4DScan',#'4DAxisScans',
        'version': 'v1',
        'include': False,
        'p_wl': ['ttH'],
        'c_wl': ['ctWctZctpcpQMAxisScan'],
        'r_wl': [],
    },
    {
        'tag': '2018_05_06/9DSetStartAxisScans',
        'grp_name': '9DSetStartAxisScans',
        'version': 'v1',
        'include': False,
        'p_wl': [],
        'c_wl': [],
        'r_wl': [],
    },
    {
        'tag': '2018_05_06/9DSetStartFullScans',
        'grp_name': '9DSetStartFullScans',
        'version': 'v1',
        'include': False,
        'p_wl': [],
        'c_wl': [],
        'r_wl': [],
    },
    {
        'tag': '2018_05_06/16DRobertSetupFullScan',
        'grp_name': '16DFullScan',
        'version': 'v1',
        'include': False,
        'p_wl': ['ttH'],#['ttH','ttll','ttlnu'],
        'c_wl': [],
        'r_wl': [],
    },
    {
        'tag': '2018_05_06/9DSetStartLXPLUSFullScan',
        'grp_name': '9DSetStartLXPLUSFullScan',
        'version': 'v1',
        'include': False,
        'p_wl': [],
        'c_wl': [],
        'r_wl': [],
    },
    {
        'tag': '2018_05_06/16DBadStart',
        'grp_name': '16DBadStart',
        'version': 'v1',
        'include': False,
        'p_wl': [],
        'c_wl': [],
        'r_wl': [],
    },
    {
        'tag': '2018_05_06/16DBadStart500k',
        'grp_name': '16DBadStart',
        'version': 'v1',
        'include': False,
        'p_wl': [],
        'c_wl': [],
        'r_wl': [],
    },
    {
        'tag': '2018_05_06/1DSMStart500k',
        'grp_name': '',
        'version': 'v1',
        'include': False,
        'p_wl': ['ttH'],#['ttH','ttll','ttlnu'],
        'c_wl': ['ctW'],
        'r_wl': [],
    },
    {
        'tag': '2018_08_24/ttHPDFCompare200k',
        'grp_name': '',
        'version': 'v1',
        'include': False,
        'p_wl': [],
        'c_wl': [],
        'r_wl': [],
    },
    {# New Reference runs (missing ttllDecay)
        'tag': '2018_08_24/TopDecayDefaultPDFsReference',
        'grp_name': 'NewRef',
        'version': 'v1',
        'include': False,
        'p_wl': ['tllqDecay'],
        'c_wl': ['cbWTopDecaysDefaultPDFs'],
        'r_wl': [],
    },
    {# Robert starting values using new PDF and Top decays (with dim6)
        'tag': '2018_08_24/16DRobertStartTopDecayDefaultPDF',
        'grp_name': '',
        'version': 'v1',
        'include': False,
        'p_wl': ['ttHDecay'],
        'c_wl': [],
        'r_wl': [],
    },
    {# Remake of the TopDecayDefaultPDFsReference gridpacks to test stability
        'tag': '2018_08_24/tllq1DRedux',
        'grp_name': '1DRedo',
        'version': 'v1',
        'include': False,
        'p_wl': ['tllqDecay'],
        'c_wl': ['cbWTopDecaysDefaultPDFs'],
        'r_wl': [],#['run1','run2','run3'],
    },
    {# Identical gridpack settings to compare production on lxplus vs CMS connect
        'tag': '2018_08_24/tllqCrossCheck',
        'grp_name': '',
        'version': 'v1',
        'include': False,
        'p_wl': [],
        'c_wl': [],
        'r_wl': [],
    },
    {# Decays top quarks w/ dim6 with 15 starting points in range [-25,25]
        'tag': '2018_08_24/HighResWithDim6Decays',
        'grp_name': 'WithDim6Decays',
        'version': 'v2',
        'include': False,
        'p_wl': ['ttHDecay'],
        'c_wl': [],
        'r_wl': [],
    },
    {# Decays top quarks w/o dim6 with 15 starting points in range [-25,25]
        'tag': '2018_08_24/HighResNoDim6Decays',
        'grp_name': 'NoDim6Decays',
        'version': 'v3',
        'include': False,
        'p_wl': [],
        'c_wl': [],
        'r_wl': [],
    },
    {# Decays top quarks w/o dim6 and top width is not auto calculated
        'tag': '2018_08_24/NoDim6DecaysNoAutoDecays',
        'grp_name': 'NoAuto',
        'version': 'v1',
        'include': False,
        'p_wl': [],
        'c_wl': [],
        'r_wl': ['run0','run3','run7','run10','run14'],
    },
    {# Has cut_decays set to false (Note: These gridpacks also had no auto top width calculation)
        'tag': '2018_08_24/CutDecaysFalseNoAuto',
        'grp_name': 'CutFalseNoAuto',
        'version': 'v1',
        'include': False,
        'p_wl': [],
        'c_wl': [],
        'r_wl': [],
    },
    {# Has cut_decays set to true (Note: These gridpacks also had no auto top width calculation)
        'tag': '2018_08_24/CutDecaysTrueNoAuto',
        'grp_name': 'CutTrueNoAuto',
        'version': 'v1',
        'include': False,
        'p_wl': [],
        'c_wl': [],
        'r_wl': [],
    },
    {# Has cut_decays set to false (Note: These gridpacks include the auto top width calculation)
        'tag': '2018_08_24/CutDecaysFalseAutoTop',
        'grp_name': 'CutFalse',
        'version': 'v1',
        'include': False,
        'p_wl': [],
        'c_wl': [],
        'r_wl': [],
    },
    {# Has cut_decays set to true (Note: These gridpacks include the auto top width calculation)
        'tag': '2018_08_24/CutDecaysTrueAutoTop',
        'grp_name': 'CutTrue',
        'version': 'v1',
        'include': False,
        'p_wl': [],
        'c_wl': [],
        'r_wl': [],
    },
    {# A check to make sure the gridpacks produce sensible results (can switch between 'ttH','ttll','ttlnu', and 'tllq')
        'tag': '2018_08_24/ttHRefCheck',
        'grp_name': 'RefCheckV3',
        'version': 'v3',
        'include': False,
        'p_wl': [],
        'c_wl': [],
        'r_wl': [],
    },
    {# A check to make sure the gridpacks produce sensible results (can switch between 'ttH','ttll','ttlnu', and 'tllq')
        'tag': '2018_08_24/ttllRefCheck',
        'grp_name': 'RefCheckV2',
        'version': 'v2',
        'include': False,
        'p_wl': [],
        'c_wl': [],
        'r_wl': [],
    },
    {# A check to make sure the gridpacks produce sensible results (can switch between 'ttH','ttll','ttlnu', and 'tllq')
        'tag': '2018_08_24/ttlnuRefCheck',
        'grp_name': 'RefCheckV3',
        'version': 'v3',
        'include': False,
        'p_wl': [],
        'c_wl': [],
        'r_wl': [],
    },
    {# A check to make sure the gridpacks produce sensible results (can switch between 'ttH','ttll','ttlnu', and 'tllq')
        'tag': '2018_08_24/tllqRefCheck',
        'grp_name': 'RefCheckV3',
        'version': 'v3',
        'include': False,
        'p_wl': [],
        'c_wl': [],
        'r_wl': [],
    },
    {# A check to make sure the gridpacks produce sensible results (can switch between 'ttH','ttll','ttlnu', and 'tllq')
        'tag': '2018_08_24/ttbarRefCheck',
        'grp_name': 'RefCheck',
        'version': 'v1',
        'include': False,
        'p_wl': [],
        'c_wl': [],
        'r_wl': [],
    },
    {# Looks specifically at tllq cQlMi to try and determine validity of the MadGraph computation
        'tag': '2018_08_24/tllqCloserLook',
        'grp_name': '',
        'version': 'v1',
        'include': False,
        'p_wl': [],
        'c_wl': [],
        'r_wl': [],
    },
    {# 1D scans for all 16 WC in all signal processes + ttbar    ###### TEST #######
        'tag': '2018_08_24/FullReference',
        'grp_name': '1D-Ref',
        'version': 'v1',
        'include': False,
        #'include': True,
        'p_wl': [],
        'c_wl': [],
        'r_wl': [],
    },
    {# 200k events using the "16DOldLimitsAxisScan" gridpacks
        'tag': '2018_08_24/16DUniformStartCheck',
        'grp_name': 'RefValidation',
        'version': 'v3',
        'include': False,
        'p_wl': ['ttll'],
        'c_wl': [],
        'r_wl': ['run1'],
    },
    {# ttbar reference scan attempting to reweight using Alexander's setup
        'tag': '2018_08_24/AlexanderDecayChainCheck',
        'grp_name': 'DecayChainCheck',
        'version': 'v2',
        'include': False,
        'p_wl': [],
        'c_wl': [],
        'r_wl': [],
    },
    {# 1D scans including top decays with no Dim6 operators
        'tag': '2018_08_24/HighResNoDim6Decays',
        'grp_name': 'TopDecays',
        'version': 'v3',
        'include': False,
        'p_wl': [],
        'c_wl': [],
        'r_wl': [],
    },
    {# Using each_coupling_order
        'tag': '2018_08_24/KelciFirstLook',
        'grp_name': '',
        'version': 'v1',
        'include': False,
        'p_wl': [],
        'c_wl': [],
        'r_wl': [],
    },
    {# tHq 1-D reference scans
        'tag': '2018_08_24/tHqRefScans',
        'grp_name': '1D-Ref',
        'version': 'v3',
        'include': False,
        'p_wl': [],
        'c_wl': [],
        'r_wl': [],
    },
    {# tHq 16-D scans using the ttll 16D scanpoints from round4
        'tag': '2018_08_24/tHq16DttllScanpoints',
        'grp_name': '',
        'version': 'v1',
        'include': False,
        'p_wl': [],
        'c_wl': [],
        'r_wl': ['run1'],
    },
    {# tHq/tHlnu 1-D reference scans for all 4Quark operators
        'tag': '2018_08_24/tHqAndtHlnuRefAll4Quark',
        'grp_name': '1D-Ref',
        'version': 'v1',
        'include': False,
        'p_wl': [],
        'c_wl': [],
        'r_wl': [],
    },
    {# tHq FullProduction sample (5M events!), Note: Code won't be able to find the scanpoints file
        'tag': 'FP/Round4/Batch6',
        'grp_name': '',
        'version': 'v1',
        'include': False,
        'p_wl': ['tHq'],
        'c_wl': [],
        'r_wl': ['run1'],
    },
    {# tHlnu 1-D reference scans
        'tag': '2018_08_24/tHlnuRefScans',
        'grp_name': '1D-Ref',
        'version': 'v1',
        'include': False,
        'p_wl': [],
        'c_wl': [],
        'r_wl': [],
    },
    {# ttbar 1-D reference scans for the 4Hvy quark operators
        'tag': '2018_08_24/ttbarRef4HvyQuarkScans',
        'grp_name': '',
        'version': 'v2',
        'include': False,
        'p_wl': [],
        'c_wl': [],
        'r_wl': [],
    },
    {# ttbar 1-D reference scans for all 4Quark operators
        'tag': '2018_08_24/ttbarRefAll4Quark',
        'grp_name': '1D-Ref',
        'version': 'v2',
        'include': False,
        'p_wl': [],
        'c_wl': [],
        'r_wl': [],
    },
    {# ttH 1-D reference scans for all 4Quark operators
        'tag': '2018_08_24/ttHRefAll4Quark',
        'grp_name': '1D-Ref',
        'version': 'v1',
        'include': False,
        'p_wl': [],
        'c_wl': [],
        'r_wl': [],
    },
    {# ttlnu 1-D reference scans for all 4Quark operators
        'tag': '2018_08_24/ttlnuRefAll4Quark',
        'grp_name': '1D-Ref',
        'version': 'v1',
        'include': False,
        'p_wl': [],
        'c_wl': [],
        'r_wl': [],
    },
    {# ttll 1-D reference scans for all 4Quark operators
        'tag': '2018_08_24/ttllRefAll4Quark',
        'grp_name': '1D-Ref',
        'version': 'v1',
        'include': False,
        'p_wl': [],
        'c_wl': [],
        'r_wl': [],
    },
    {# tllq 1-D reference scans for all 4Quark operators
        'tag': '2018_08_24/tllqRefAll4Quark',
        'grp_name': '1D-Ref',
        'version': 'v1',
        'include': False,
        'p_wl': [],
        'c_wl': [],
        'r_wl': [],
    },
    {# tHq/tHlnu 1-D reference scans for all 4Quark operators
        'tag': '2018_08_24/tHqAndtHlnuRefAll4Quark',
        'grp_name': '1D-Ref',
        'version': 'v1',
        'include': False,
        'p_wl': [],
        'c_wl': [],
        'r_wl': [],
    },
    {# ttll 16D samples using central cuts in the MG card, uses exact same scanpoints as in ttll_16DOldLimits*_run1
        'tag': '2019_04_19/UsingCentralCuts',
        'grp_name': 'CentralCuts',
        'version': 'v2',    #v1 = ttll, v2 = tllq
        'include': False,
        'p_wl': [],
        'c_wl': [],
        'r_wl': [],
    },
    {# Validation set of events for samples with +1 parton using Han's model
        'tag': '2019_04_19/ValidationHanModelPlusJet',
        'grp_name': 'plus1Jet',
        'version': 'v1',
        'include': False,
        'p_wl': [],
        'c_wl': [],
        'r_wl': [],
    },
    {# Validation set of events for samples with +1 parton using Han's model
        'tag': '2019_04_19/ValidationHanModelPlusJet-mAOD',
        'grp_name': 'matched-1Jet',
        'version': 'v1',
        'include': False,
        'p_wl': ['ttHJet'],
        'c_wl': ['HanModel16DttllScanpoints'],
        'r_wl': [],
    },
    {# Validation set of events for samples with +1 parton using Han's model (xqcut/qCut tests)
        'tag': '2019_04_19/ttHJet-xqcutStudies-mAOD',
        'grp_name': 'v2-qcut30',
        'version': 'v2',
        'include': False,
        'p_wl': [],
        'c_wl': [],
        'r_wl': [],
    },
    {# Validation set of events for samples with +1 parton using Han's model (xqcut/qCut tests)
        'tag': '2019_04_19/ttHJet-xqcutStudies-mAOD',
        'grp_name': 'v3-qcut45',
        'version': 'v3',
        'include': False,
        'p_wl': [],
        'c_wl': [],
        'r_wl': [],
    },
    {# Validation set of events for samples with +1 parton using Han's model (xqcut/qCut tests)
        'tag': '2019_04_19/ttHJet-xqcutStudies-mAOD',
        'grp_name': 'v4-qcut19',
        'version': 'v4',
        'include': False,
        'p_wl': [],
        'c_wl': [],
        'r_wl': [],
    },
    {# Validation set of events for samples with +1 parton using Han's model (xqcut/qCut tests)
        'tag': '2019_04_19/ttHJet-xqcutStudies-mAOD',
        'grp_name': 'v5-qcut30',
        'version': 'v5',
        'include': False,
        'p_wl': [],
        'c_wl': [],
        'r_wl': [],
    },
    {# Validation set of events for samples with +1 parton using Han's model (xqcut/qCut tests)
        'tag': '2019_04_19/ttHJet-xqcutStudies-mAOD',
        'grp_name': 'v6-qcut45',
        'version': 'v6',
        'include': False,
        'p_wl': [],
        'c_wl': [],
        'r_wl': [],
    },
    {# Validation set of events for samples with +1 parton using Han's model (xqcut/qCut tests)
        'tag': 'ttHJet-xqcutStudies-mAOD',
        'grp_name': '',
        'version': 'vTotal',
        #'include': True,
        'include': False,
        'p_wl': [],
        #'c_wl': [],
        'c_wl': ["HanModelxqcut10qCut30","HanModelxqcut10qCut19"],
        'r_wl': [],
    },
    {# Validation set of events for samples with +1 parton using Han's model (xqcut/qCut tests)
        'tag': 'ttHJet-xqcutStudies-qCutScan-GEN',
        'grp_name': '',
        'version': 'v1',
        #'include': True,
        'include': False,
        'p_wl': [],
        #'c_wl': [
        #    "HanModelxqcut05qCut10",
        #    "HanModelxqcut05qCut19",
        #    "HanModelxqcut05qCut25",
        #    "HanModelxqcut05qCut40",
        #    "HanModelxqcut05qCut60",
        #    "HanModelxqcut10qCut19",
        # ],
        'c_wl': [],
        'r_wl': [],
    },
    {# Validation set of events for samples with +1 parton using Han's model (xqcut/qCut tests): Scanpoints file
        'tag': 'ttHJet-xqcutStudies-xqcut10qCutTests-GEN',
        'grp_name': '',
        'version': 'v2',
        #'include': True,
        'include': False,
        'p_wl': [],
        'c_wl': [
            "HanModel16DttllScanpointsqCut15",
            "HanModel16DttllScanpointsqCut19",
            "HanModel16DttllScanpointsqCut25",
         ],
        #'c_wl': [],
        'r_wl': [],
    },
    {# Validation set of events for samples with +1 parton using Han's model (xqcut/qCut tests): Scanpoints file (no ctG)
        'tag': 'ttHJet-HanModelNoctG16DttllScanpoints-GEN',
        'grp_name': '',
        'version': 'v1',
        #'include': True,
        'include': False,
        'p_wl': [],
        'c_wl': [],
        'r_wl': [],
    },
    {# TEST with different base paths
        'basepath': '/hadoop/store/user/kmohrman/summaryTree_LHE/2019_04_19/',
        'tag': 'ttHJet-xqcutStudies-qCutScan-GEN',
        'grp_name': '',
        'version': 'v1',
        #'include': True,
        'include': False,
        'p_wl': [],
        'c_wl': ["HanModelxqcut10qCut19"],
        'r_wl': [],
    },
    {# TEST with different base paths
        'basepath': '/afs/crc.nd.edu/user/k/kmohrman/CMSSW_Releases/CMSSW_9_4_6/src/EFTGenReader/LHEReader/test/rwgt_validation/root_files_from_hadoop/',
        'tag': 'ttHJet-xqcutStudies-mAOD',
        'grp_name': '',
        'version': 'vTotal',
        #'include': True,
        'include': False,
        'p_wl': [],
        'c_wl': ["HanModelxqcut10qCut19"],
        'r_wl': [],
    },
    {# Validation set of events for samples with +1 parton using Han's model
        'tag': 'reference_scans/HanModel_1jet/ttH-tllq4f-tHq4f-ttlnu-GEN',
        'grp_name': 'Jet1DRef',
        'version': 'v1',
        'include': False,
        'p_wl': [],
        'c_wl': [],
        'r_wl': [],
    },
    {# Validation set of events for samples with +1 parton using Han's model (xqcut/qCut tests): Scanpoints files xqcut tests
        'tag': 'ttHJet-xqcutStudies-xqcutScan_analysisEtaCut-GEN',
        'grp_name': '',
        'version': 'v1',
        #'include': True,
        'include': False,
        'p_wl': [],
        'c_wl': ["HanModel16DttllScanpointsxqcut05qCut10",
            "HanModel16DttllScanpointsxqcut05qCut25",
            "HanModel16DttllScanpointsqCut25",
            "HanModel16DttllScanpointsxqcut15qCut25",
            "HanModel16DttllScanpointsxqcut20qCut25",
        ],
        'r_wl': [],
    },
    {# FP samples: 
        'tag': 'allProcesses-mAOD',
        'grp_name': '',
        'version': 'v1',
        #'include': True,
        'include': False,
        'p_wl': [],
        'c_wl': [],
        'r_wl': [],
    },
    {# Han model v2 xqcut 10 qcut 15,19,25
        'tag': 'ttHJet-HanV2Model-xqcut10qCutTests_analysisEtaCut-GEN',
        'grp_name': '',
        'version': 'v1',
        #'include': True,
        'include': False,
        'p_wl': [],
        'c_wl': [],
        'r_wl': [],
    },
    {## Han model v4 xqcut 10 qcut 15,19,25
        'tag': 'ttXJet-HanV4Model-xqcut10qCutTests_analysisEtaCut-GEN',
        'grp_name': '',
        'version': 'v1',
        #'include': True,
        'include': False,
        'p_wl': ["ttHJet"],
        'c_wl': [],
        'r_wl': [],
    },
    {# Han model v4 0jet vs 1jet comp ####### (Use this for 0j vs 1j tests for all but ttZ)! #######
        'tag': 'ttX-ttXJet-HanV4Model-0Jetvs1JetTests_analysisEtaCut-GEN',
        'grp_name': '',
        'version': 'v1',
        #'include': True,
        'include': False,
        #'p_wl': ["ttH","ttHJet","ttlnuJet","ttlnu","ttllNuNuNoHiggs"],
        'p_wl': ["ttllNuNuNoHiggs"],
        'c_wl': [],
        'r_wl': [],
    },
    {# Han model v4 0jet vs 1jet comp ####### Has ttZJet runs 1,2,3  #######
        'tag': 'ttX-ttXJet_HanV4_0Jetvs1JetTests_with-ttZjetRun2Run3_analysisEtaCut-GEN',
        'grp_name': '',
        'version': 'v1',
        'include': True,
        #'include': False,
        #'p_wl': ["ttH","ttHJet","ttlnuJet","ttlnu","ttllNuNuJetNoHiggs","ttllNuNuNoHiggs"],
        'p_wl': ["ttH","ttHJet"],
        #'p_wl': ["ttH","ttHJet","ttllNuNuJetNoHiggs","ttllNuNuNoHiggs"],
        #'p_wl': ['ttllNuNuJetNoHiggs','ttllNuNuNoHiggs'],
        'c_wl': [],
        'r_wl': [],
    },
    {# Han model v4 vs han origional model comp
        'tag': 'ttXJet_R5B1-HanV4Model-Comp_analysisEtaCut-mAOD',
        'grp_name': '',
        'version': 'v1',
        #'include': True,
        'include': False,
        'p_wl': ["ttlnuJet","ttHJet","ttllNuNuJetNoHiggs"],
        'c_wl': ["HanModel16DttllScanpoints"],
        'r_wl': [],
    },
    {# ttH and ttW: Han model v4 vs Han model v2 (w/o missing 5 verts added in) vs han origional model comp
        'tag': 'ttHJet-ttWJet_R5B1-HanV4Model-Comp_analysisEtaCut-mAOD',
        'grp_name': '',
        'version': 'v1',
        #'include': True,
        'include': False,
        'p_wl': [],
        'c_wl': [],
        'r_wl': [],
    },
    {# ttZ (and also ttH, ttW): Han model v4 vs Han model v2 (w/o missing 5 verts added in) vs han origional model comp
        'tag': 'ttXJet_R5B1-HanV2ModelNOttggh-HanV4Model-Comp_analysisEtaCut-mAOD',
        'grp_name': '',
        'version': 'v1',
        #'include': True,
        'include': False,
        'p_wl': ["ttllNuNuJetNoHiggs"],
        'c_wl': [],
        'r_wl': [],
    },
    {# HanModel (original) 0jet vs 1jet (FP) comp
        'tag': 'ttX-ttXJet_R5B1-HanModel-0Jetvs1JetTests_analysisEtaCut-mAOD',
        'grp_name': '',
        'version': 'v1',
        #'include': True,
        'include': False,
        'p_wl': [],
        'c_wl': [],
        'r_wl': [],
    },
    {# HanModelV4 starting points comp for ttHJet (runs 0,2) and ttWJet (all runs)
        'tag': 'ttHJet-ttWJet_HanV4ttXJetStartPtChecks-xqcut10qCut19_analysisEtaCut-GEN',
        'grp_name': '',
        'version': 'v1',
        #'include': True,
        'include': False,
        'p_wl': [],
        'c_wl': [],
        'r_wl': [],
    },
    {# HanModelV4 starting points comp for ttHJet (all runs) and ttZJet (runs 2 and 3)
        'tag': 'ttHJet-ttZJet_HanV4ttXJetStartPtChecks-xqcut10qCut19_analysisEtaCut-GEN',
        'grp_name': '',
        'version': 'v1',
        #'include': True,
        'include': False,
        'p_wl': ["ttllNuNuJetNoHiggs"],
        'c_wl': [],
        'r_wl': ["run2"],
    },
    {# HanModelV4 R6B1 starting point ttXJet qCut tests
        'tag': 'ttXJet_HanV4ttXJetStartPtChecks-xqcut10qCutTests_analysisEtaCut-GEN',
        'grp_name': '',
        'version': 'v2',
        #'include': True,
        'include': False,
        'p_wl': ["ttHJet"], ##### Use this version of tthJet for the comparison to dedicated gridpacks #####
        'c_wl': ["HanV4ttXJetStartPtChecksqCut19"],
        'r_wl': [],
    },
    {# HanModelV4 starting point chekcs for tHq 100k events samples
        'tag': 'tHq4f_HanV4tHqStartPtChecks-allRunsMatchOff_analysisEtaCut-GEN',
        'grp_name': '',
        'version': 'v1',
        #'include': True,
        'include': False,
        'p_wl': [],
        'c_wl': [],
        'r_wl': [],
    },
    {# HanModelV4 starting point chekcs for tZq 100k events samples
        'tag': 'tllq4fNoSchanWNoHiggs0p_HanV4tZqStartPtChecks-allRunsMatchOff_analysisEtaCut-GEN',
        'grp_name': '',
        'version': 'v1',
        #'include': True,
        'include': False,
        'p_wl': [],
        'c_wl': [],
        'r_wl': ["run0"],
    },
    {# HanModelV4 Round 6, (batch depends on HADOOP_BASE_PATH) FP samples:
        #'tag': 'ttXjet-mAOD', # Batch1 (ttZ, ttW)
        #'tag': 'tHq4f-mAOD', # Batch3 (tHq)
        'tag': 'tZq4f-mAOD', # Batch4 (tZq)
        #'tag': 'ttHjet-mAOD', # Batch7 (ttH good start pt)
        'grp_name': '',
        'version': 'v1',
        #'include': True,
        'include': False,
        #'p_wl': ["ttllNuNuJetNoHiggs","ttlnuJet"], # ttH from batch 1 had the bad start pt
        'p_wl': [],
        'c_wl': [],
        'r_wl': [],
    },
    {# HanModelV4 ttHJet dedicated ctG=[-3,3] axis scan ####### AXIS SCAN #######
        'tag': 'ttHJet_HanV4ctGAxisScan_analysisEtaCut-GEN',
        'grp_name': '',
        'version': 'v1',
        #'include': True,
        'include': False,
        'p_wl': [],
        'c_wl': [],
        'r_wl': [],
    },
    {# HanModelV4 ttHJet dedicated ctW=[-4,4] axis scan (for smeft comp, qed1, qcd2, dim6=2) ####### AXIS SCAN #######
        'tag': 'ttHJet_HanV4_cbW-AxisScan-withRwgt_smeftComp_QED1_QCD2_DIM62-GEN',
        'grp_name': '',
        'version': 'v1',
        #'include': True,
        'include': False,
        'p_wl': [],
        'c_wl': [],
        'r_wl': [],
        'basepath' : "/hadoop/store/user/kmohrman/summaryTree_LHE/2020_03_03_addPSweights/",
    },
    {# HanModelV4 ttHJet xqcut and qcut scan (12 samples in all)
        'tag': 'ttHJet_HanV4xqcutTests_analysisEtaCut-GEN',
        'grp_name': '',
        'version': 'v1',
        #'include': True,
        'include': False,
        'p_wl': ["ttHJet"],
        'c_wl': ["HanV4ttXjetxqcut10qCut19"],
        'r_wl': [],
    },
    {# HanModelV4 ttXJet comp with SMEFT (in proc card: QED=1, QCD=2, DIM6=2, and ttZ not ttll, ttW not ttlnu)
        # NOTE: This is a bad comparison! DIM6 should NOT be 2 for this comparison (NP=2 in smeft NLO is more or less DIM6=1, not DIM6=2 in dim6Top)
        'tag': 'ttXJet_HanV4_semftComp_QED1_QCD2_DIM62-GEN',
        'grp_name': '',
        'version': 'v1',
        #'include': True,
        'include': False,
        #'p_wl': ["ttHJetSMEFTcomp","ttWJetSMEFTcomp","ttZJetSMEFTcomp"],
        'p_wl': ["ttHJetSMEFTcomp"],
        'c_wl': [],
        'r_wl': [],
        'basepath' : "/hadoop/store/user/kmohrman/summaryTree_LHE/2020_03_03_addPSweights/",
    },
    {# HanModelV4 ttHJet, comp with NLO, so QED=1, QCD=2 (and DIM6=1)
        'tag': 'ttHJet_HanV4_withRwgt_smeftComp_QED1_QCD2_DIM61-GEN',
        'grp_name': '',
        'version': 'v1',
        'include': True,
        #'include': False,
        'p_wl': [],
        'c_wl': [],
        'r_wl': [],
        'basepath' : "/hadoop/store/user/kmohrman/summaryTree_LHE/2020_03_03_addPSweights",
    },
]

# These are the tags to get the runs which should be used as reference points in the plotting
REF_TAGS = [
    #'ttbar1DAxisScans',
    #'AllProcessAllCoeffsAxisScans',
    #'HighResNoDim6Decays',
    #'ttHRefCheck',
    #'ttllRefCheck',
    #'ttlnuRefCheck',
    #'tllqRefCheck',
    #'ttbarRefCheck',
    #'tllqCloserLook',
    #'FullReference',
    #'tHqRefScans',
    #'tHlnuRefScans',
    #'ttbarRefAll4Quark',
    #'ttHRefAll4Quark',
    #'ttlnuRefAll4Quark',
    #'ttllRefAll4Quark',
    #'tllqRefAll4Quark',
    #'tHqAndtHlnuRefAll4Quark',
    #'AlexanderDecayChainCheck'
    #'HighResNoDim6Decays'
    #'ttH-tllq4f-tHq4f-ttlnu-GEN'

    # Axis scans (don't need to include as "True" in ALL_INFO)
    #'ttHJet_HanV4ctGAxisScan_analysisEtaCut-GEN'
    #'ttHJet_HanV4_cbW-AxisScan-withRwgt_smeftComp_QED1_QCD2_DIM62-GEN'
]

# Dictionary to map certain MG processes to another for use in finding reference samples
REF_PROCESS_MAP = {
    #'tllq4fMatchedNoHiggs': 'tllq',
    #'ttHJet': 'ttH',
    #'ttH': 'ttHJet',
    'ttH': 'ttHJetSMEFTcomp'
    #'ttlnuJet': 'ttlnu',
    #'tHq4fMatched': 'tHq',
    #'ttllNuNuJetNoHiggs': 'ttll'
}

#HADOOP_PATH = "/hadoop/store/user/awightma/summaryTree_LHE/%s/%s/" % (grp_tag,version)
#HADOOP_BASE_PATH = "/hadoop/store/user/awightma/summaryTree_LHE/"
#HADOOP_BASE_PATH = "/hadoop/store/user/kmohrman/summaryTree_LHE/"

#HADOOP_BASE_PATH = "/afs/crc.nd.edu/user/k/kmohrman/CMSSW_Releases/CMSSW_9_4_6/src/EFTGenReader/LHEReader/test/rwgt_validation/root_files_from_hadoop/"
#HADOOP_BASE_PATH = "/hadoop/store/user/kmohrman/summaryTree_LHE/2019_04_19/"
HADOOP_BASE_PATH = "/hadoop/store/user/kmohrman/summaryTree_LHE/2019_08_14_addPtBranches/" ###
#HADOOP_BASE_PATH = "/hadoop/store/user/kmohrman/summaryTree_LHE/2020_03_03_addPSweights/"

#HADOOP_BASE_PATH = "/hadoop/store/user/kmohrman/summaryTree_LHE/FP/Round5/Batch1/"
#HADOOP_BASE_PATH = "/hadoop/store/user/kmohrman/summaryTree_LHE/FP/Round6/Batch1/"
#HADOOP_BASE_PATH = "/hadoop/store/user/kmohrman/summaryTree_LHE/FP/Round6/Batch3/"
#HADOOP_BASE_PATH = "/hadoop/store/user/kmohrman/summaryTree_LHE/FP/Round6/Batch4/"
#HADOOP_BASE_PATH = "/hadoop/store/user/kmohrman/summaryTree_LHE/FP/Round6/Batch7/"

def runByProcess():
    #file_dirs = groupByProcess(HADOOP_PATH,grp_tag,process_whitelist,coeff_whitelist,run_whitelist)
    #file_dirs = groupByCoefficient(HADOOP_PATH,grp_tag,process_whitelist,coeff_whitelist,run_whitelist)

    #NOTE: The output name could be duplicated and overwrite a previous run
    file_dirs = {}
    all_grouped_file_dirs_dict = {}
    spacing = 0
    for idx,info in enumerate(ALL_INFO):
        if not info['include']:
            continue
        # TEST for multiple base paths!
        if "basepath" in info.keys():
            path = os.path.join(info['basepath'],info['tag'],info['version'])
        else:
            path = os.path.join(HADOOP_BASE_PATH,info['tag'],info['version'])
        #path = os.path.join(HADOOP_BASE_PATH,info['tag'],info['version'])
        #grouped_dirs = groupByProcess(path,info['tag'],info['grp_name'],info['p_wl'],info['c_wl'],info['r_wl'])
        grouped_dirs = groupByProcess(path,info['tag'],'',info['p_wl'],info['c_wl'],info['r_wl'])
        print "\nThe current grouped dirs are:" , grouped_dirs , "\n"
        for tup,dirs in grouped_dirs.iteritems():
            if "ttH" in tup[1]:
                proc = "ttH"
            elif "ttlnu" in tup[1] or "ttW" in tup[1]:
                proc = "ttW"
            elif "ttll" in tup[1] or "ttZ" in tup[1]:
                proc = "ttZ"
            else:
                print "\nError: Unknown process" , tup[1] , "exiting"
                raise BaseException
            if not all_grouped_file_dirs_dict.has_key(proc):
                all_grouped_file_dirs_dict[proc] = []
            all_grouped_file_dirs_dict[proc].extend(dirs)
            #if not file_dirs.has_key(tup):
            #    file_dirs[tup] = []
            #file_dirs[tup].extend(dirs)
        #file_dirs.update(grouped_dirs)
        #for k in file_dirs.keys():
        #    if len(k[2]) > 0:
        #        spacing = max([spacing,len("%s_%s" % (k[1],k[2]))])
        #    else:
        #        spacing = max([spacing,len("%s" % (k[1]))])

    # Run root macro once per process
    count = 0
    #for tup,fdirs in file_dirs.iteritems():
    #    if len(tup) != 3:
    #        print "[WARNING] Unknown file tuple,",tup
    #        continue
    #    process = tup[1]
    #    grp_name = tup[2]
    #    if len(grp_name) > 0:
    #        output_name = process + "_" + grp_name
    #    else:
    #        output_name = process
    print "ALL GRP DIRS !!!!" , all_grouped_file_dirs_dict
    for proc,fdirs in all_grouped_file_dirs_dict.iteritems():
        output_name = proc

        ref_dirs = []
        for rtag in REF_TAGS:
            info_list = getInfoByTag(ALL_INFO,rtag)
            for info in info_list:
                if "basepath" in info.keys():
                    ref_path = os.path.join(info['basepath'],info['tag'],info['version'])
                else:
                    ref_path = os.path.join(HADOOP_BASE_PATH,info['tag'],info['version'])
                #ref_path = os.path.join(HADOOP_BASE_PATH,info['tag'],info['version'])
                #search_proc = process
                search_proc = proc
                #if REF_PROCESS_MAP.has_key(process):
                print "proc: " , proc , REF_PROCESS_MAP.keys()
                if REF_PROCESS_MAP.has_key(proc):
                    #search_proc = REF_PROCESS_MAP[process]
                    search_proc = REF_PROCESS_MAP[proc]
                ref_dirs += getDirectories(ref_path,
                    p_wl=[search_proc],
                    c_wl=[],#'cQeiRefV1AxisScan'],
                    r_wl=[]
                )
        #ref_dirs = []
        print "REF DIRS !!! " , ref_dirs

        print tup
        for x in fdirs:
            print "\t{}".format(x)

        dir_inputs = 'dirpaths.txt'
        ref_inputs = 'refpaths.txt'
        with open(dir_inputs,'w') as f:
            # Write the path directories to all relevant coeffs/runs
            for fd in sorted(fdirs):
                l = "%s\n" % (fd)
                f.write(l)
        with open(ref_inputs,'w') as f:
            for fd in sorted(ref_dirs):
                l = "%s\n" % (fd)
                f.write(l)

        print "info given to rungridpackvalidation: " , "out name: " , output_name , "dir_inputs:", dir_inputs , "ref_inputs:" , ref_inputs
        #print "[%d/%d] %s (dirs %d, ref %d):" % (count+1,len(file_dirs.keys()),output_name.ljust(spacing),len(fdirs),len(ref_dirs))
        print "[%d/%d] %s (dirs %d, ref %d):" % (count+1,len(all_grouped_file_dirs_dict.keys()),output_name.ljust(spacing),len(fdirs),len(ref_dirs))
        #subprocess.check_call(['root','-b','-l','-q','readLHEOutputTree.C+(\"%s\",\"%s\",\"%s\")' % (output_name,dir_inputs,grp_name)])
        #subprocess.check_call(['root','-b','-l','-q','rwgtDistributions.C+(\"%s_limit_fit\",\"%s\",\"%s\")' % (process,dir_inputs,"")])
        #subprocess.check_call(['root','-b','-l','-q','runXsec.C+(\"%s\",\"%s\",\"%s\")' % (process,dir_inputs,grp_name)])
        #subprocess.check_call(['root','-b','-l','-q','runGridpackValidation.C+(\"%s\",\"%s\",\"%s\",\"%s\")' % (output_name,dir_inputs,ref_inputs,grp_name)])
        subprocess.check_call(['root','-b','-l','-q','runGridpackValidation.C+(\"%s\",\"%s\",\"%s\")' % (output_name,dir_inputs,ref_inputs)])
        #subprocess.check_call(['root','-b','-l','-q','run2DComparison.C+(\"%s\",\"%s\",\"%s\")' % (output_name,dir_inputs,grp_name)])
        count += 1

def runByCoeff(tags,runs):
    procs = ["ttH","ttll","ttlnu","tllq"]
    for idx,tag in enumerate(tags):
        print "[%d/%d] %s:" % (idx+1,len(tags),tag)
        fdirs = []  # fdirs can have runs from different processes/grp_tags/runs, so make sure we can handle that
        for idx,info in enumerate(ALL_INFO):
            if not info['include']:
                continue
            path = os.path.join(HADOOP_BASE_PATH,info['tag'],info['version'])
            fdirs += getDirectories(path,
                p_wl=procs,
                c_wl=tag,
                r_wl=runs
            )
        dir_inputs = 'dirpaths.txt'
        with open(dir_inputs,'w') as f:
            for fd in sorted(fdirs):
                l = "%s\n" % (fd)
                f.write(l)
        subprocess.check_call(['root','-b','-l','-q','runLayeredPlots.C+(\"%s\")' % (dir_inputs)])
    return



runByProcess()

#tags = ["ctGRefV1AxisScan","cpQMRefV1AxisScan"]; runs = ["run5"]
#runByCoeff(tags,runs)

#tags = ["ctpRefV1AxisScan"]; runs = ["run2"]
#runByCoeff(tags,runs)

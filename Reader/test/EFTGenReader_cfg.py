import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing
import sys
import os

from EFTGenReader.Reader.DatasetHelper import DatasetHelper

options = VarParsing.VarParsing('analysis')
#options.maxEvents = 100000#-1
options.maxEvents = 10

nd_redirect = "root://ndcms.crc.nd.edu/"
fnal_redirect = "root://cmsxrootd.fnal.gov/"
global_redirect = "root://cms-xrd-global.cern.ch/"

ds_helper = DatasetHelper()
ds_helper.load("../data/JSON/datasets.json")
ds_helper.root_redirect = nd_redirect

root_redirect = nd_redirect

central_ttZ_xsec = 0.244109
central_ttZ = [
    "/store/mc/RunIIFall17MiniAODv2/TTZToLLNuNu_M-10_TuneCP5_PSweights_13TeV-amcatnlo-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/20000/024736F4-0150-E811-BB33-3417EBE34CAB.root",
    "/store/mc/RunIIFall17MiniAODv2/TTZToLLNuNu_M-10_TuneCP5_PSweights_13TeV-amcatnlo-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/30000/FEA77880-4176-E811-9434-0CC47AFB8104.root",
    "/store/mc/RunIIFall17MiniAODv2/TTZToLLNuNu_M-10_TuneCP5_PSweights_13TeV-amcatnlo-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/30000/FE92B0F0-E683-E811-AF02-1866DAEB296C.root",
]
for i,s in enumerate(central_ttZ): central_ttZ[i] = root_redirect + central_ttZ[i]

central_ttW_xsec = 0.341762
central_ttW = [
    "/store/mc/RunIIFall17MiniAODv2/TTWJetsToLNu_TuneCP5_PSweights_13TeV-amcatnloFXFX-madspin-pythia8/MINIAODSIM/PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14-v1/70000/FE95B21E-35AD-E811-8F38-5065F37D9082.root",
    "/store/mc/RunIIFall17MiniAODv2/TTWJetsToLNu_TuneCP5_PSweights_13TeV-amcatnloFXFX-madspin-pythia8/MINIAODSIM/PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14-v1/70000/F88C7524-99B2-E811-83B4-6CC2173C4580.root",
    "/store/mc/RunIIFall17MiniAODv2/TTWJetsToLNu_TuneCP5_PSweights_13TeV-amcatnloFXFX-madspin-pythia8/MINIAODSIM/PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14-v1/70000/F6BA9F5D-3FAF-E811-B5CF-FA163EC5703C.root",
    "/store/mc/RunIIFall17MiniAODv2/TTWJetsToLNu_TuneCP5_PSweights_13TeV-amcatnloFXFX-madspin-pythia8/MINIAODSIM/PU2017_12Apr2018_new_pmx_94X_mc2017_realistic_v14-v1/70000/F60FE38F-52AE-E811-9B84-44A842CFD626.root",
]
for i,s in enumerate(central_ttW): central_ttW[i] = root_redirect + central_ttW[i]

central_tZq_xsec = 0.07473
central_tZq = [
    "/store/mc/RunIIFall17MiniAODv2/tZq_ll_4f_ckm_NLO_TuneCP5_PSweights_13TeV-amcatnlo-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/90000/24D78155-C843-E811-B2A8-7CD30AC030A0.root",
    "/store/mc/RunIIFall17MiniAODv2/tZq_ll_4f_ckm_NLO_TuneCP5_PSweights_13TeV-amcatnlo-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/50000/FEE79580-7742-E811-AA75-90B11C2CC96F.root",
    "/store/mc/RunIIFall17MiniAODv2/tZq_ll_4f_ckm_NLO_TuneCP5_PSweights_13TeV-amcatnlo-pythia8/MINIAODSIM/PU2017_12Apr2018_94X_mc2017_realistic_v14-v1/50000/E2273B78-4941-E811-9B8B-0025905C3DD8.root",
]
for i,s in enumerate(central_tZq): central_tZq[i] = root_redirect + central_tZq[i]

from Configuration.StandardSequences.Eras import eras
process = cms.Process("Demo", eras.Run2_2017)

process.load('FWCore.MessageService.MessageLogger_cfi')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents) # number of events
)

process.MessageLogger.cerr.FwkReport.reportEvery = 500

hadoop_protocol = "file://"
user_path = "/hadoop/store/user/awightma/"

loc_ttll_FP_R4B9       = "FullProduction/Round4/Batch9/postLHE_step/v1/mAOD_step_ttll_16DcentralCuts_run0/"
loc_ttllNoHiggs_SM     = "postLHE_step/2019_04_19/ttll-tllq-ProcessCardStudies/v2/mAOD_step_ttllNoHiggs_NoDim6_run0/"
loc_ttllNoHiggs_EFT    = "postLHE_step/2019_04_19/ttll-tllq-ProcessCardStudies/v2/mAOD_step_ttllNoHiggs_16DttllScanpoints_run1/"

loc_tllq_FR_R4B9       = "FullProduction/Round4/Batch9/postLHE_step/v1/mAOD_step_tllq_16DcentralCuts_run0/"
loc_tllq4f_SMNoSchanW  = "postLHE_step/2019_04_19/tllq4f-NoDim6Diagrams/v2/mAOD_step_tllq4f_NoDim6NoSchanW_run0/"
loc_tllq4fMatched_SM   = "postLHE_step/2019_04_19/ttll-tllq-ProcessCardStudies/v1/mAOD_step_tllq4fMatched_NoDim6_run0/"
loc_tllq4fMatched_EFT  = "postLHE_step/2019_04_19/ttll-tllq-ProcessCardStudies/v2/mAOD_step_tllq4fMatched_16DttllScanpoints_run1/"
loc_tllq4fNoHiggs_SM   = "postLHE_step/2019_04_19/ttll-tllq-ProcessCardStudies/v1/mAOD_step_tllq4fNoHiggs_NoDim6_run0/"
loc_tllq4fNoHiggs_EFT  = "postLHE_step/2019_04_19/ttll-tllq-ProcessCardStudies/v2/mAOD_step_tllq4fNoHiggs_16DttllScanpoints_run1/"

loc_ttlnu_FP_R4B9      = "FullProduction/Round4/Batch9/postLHE_step/v1/mAOD_step_ttlnu_16DcentralCuts_run0/"
loc_ttlnu_EFT          = "postLHE_step/2019_04_19/ttlnu-NoDim6Diagrams/v2/mAOD_step_ttlnu_NoDim6_run0/"
loc_ttlnuJet_EFT       = "postLHE_step/2019_04_19/ttlnu-NoDim6Diagrams/v2/mAOD_step_ttlnuJet_1JetNoDim6_run0/"
loc_ttlnu_NoPDFWeights = "postLHE_step/2019_04_19/ttlnuJet-NoPDFWeights/v1/mAOD_step_ttlnuJet_16DttlnuScanpointsAutoJetCuts_run1/"

priv_ttll_FP_R4B9    = [loc_ttll_FP_R4B9   ,"privateTTZOldCuts_NoTopLeptons_output_tree.root"       ,True ,central_ttZ_xsec]
priv_ttllNoHiggs_SM  = [loc_ttllNoHiggs_SM ,"privateTTZ_NoTopLeptons_output_tree.root"              ,False,central_ttZ_xsec]
priv_ttllNoHiggs_EFT = [loc_ttllNoHiggs_EFT,"privateTTZNoHiggs-NoDim6_NoTopLeptons_output_tree.root",True ,central_ttZ_xsec]

priv_tllq_FR_R4B9      = [loc_tllq_FR_R4B9     ,"privateTZQ_NoTopLeptons_output_tree.root"                  ,True ,central_tZq_xsec]
priv_tllq4f_SMNoSchanW = [loc_tllq4f_SMNoSchanW,"privateTZQ4f-NoDim6-NoSchanW_NoTopLeptons_output_tree.root",False,central_tZq_xsec]
priv_tllq4fMatched_SM  = [loc_tllq4fMatched_SM ,"privateTZQ4fMatched-NoDim6_NoTopLeptons_output_tree.root"  ,False,central_tZq_xsec]
priv_tllq4fMatched_EFT = [loc_tllq4fMatched_EFT,"privateTZQ4fMatched-EFT_NoTopLeptons_output_tree.root"     ,True ,central_tZq_xsec]
priv_tllq4fNoHiggs_SM  = [loc_tllq4fNoHiggs_SM ,"privateTZQ4fNoHiggs-NoDim6_NoTopLeptons_output_tree.root"  ,False,central_tZq_xsec]
priv_tllq4fNoHiggs_EFT = [loc_tllq4fNoHiggs_EFT,"privateTZQ4fNoHiggs-EFT_NoTopLeptons_output_tree.root"     ,True ,central_tZq_xsec]

priv_ttlnu_FP_R4B9      = [loc_ttlnu_FP_R4B9     ,"privateTTZNoHiggs-EFT_NoTopLeptons_output_tree.root" ,True,central_ttW_xsec]
priv_ttlnu_EFT          = [loc_ttlnu_EFT         ,"privateTTW_NoTopLeptons_output_tree.root"            ,True,central_ttW_xsec]
priv_ttlnuJet_EFT       = [loc_ttlnuJet_EFT      ,"privateTTW-NoDim6-1Jet_NoTopLeptons_output_tree.root",True,central_ttW_xsec]
priv_ttlnu_NoPDFWeights = [loc_ttlnu_NoPDFWeights,"privateTTW-NoDim6_NoTopLeptons_output_tree.root"     ,True,central_ttW_xsec]

priv_loc,out_fname,iseft,NLO_xsec_norm = priv_tllq4fMatched_EFT
#iseft = False
#out_fname = "TEST_output_tree.root"
#out_fname = "centralTTZ_NoTopLeptons_output_tree.root"
#out_fname = "centralTTW_NoTopLeptons_output_tree.root"
#out_fname = "centralTZQ_NoTopLeptons_output_tree.root"

private_files = []
private_path = os.path.join(user_path,priv_loc)
for fname in os.listdir(private_path):
    fpath = os.path.join(private_path,fname)
    if not ".root" in fname: continue
    private_files.append(hadoop_protocol + fpath)

ds_name = 'tllq4fMatched_SM'

files     = ds_helper.getFiles(ds_name)
is_eft    = ds_helper.getData('is_eft')
xsec_norm = ds_helper.getData('central_xsec')

out_fname = "%s_NoTopLeptons_output_tree.root" % (ds_name)
out_fname = "TEST_output_tree.root"

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        *files
        #*private_files
        #*central_ttZ
        #*central_ttW
        #*central_tZq
    )
)

#process.load("ttH-13TeVMultiLeptons.TemplateMakers.EFTGenReader_cfi")
process.load("EFTGenReader.Reader.EFTGenReader_cfi")

process.EFTGenReader.debug = False
process.EFTGenReader.iseft = is_eft
process.EFTGenReader.gp_events = 500
process.EFTGenReader.norm_type = 0      # 0 - No norm, 1 - unit norm, 2 - xsec norm
process.EFTGenReader.xsec_norm = xsec_norm
process.EFTGenReader.intg_lumi = 1.0

process.TFileService = cms.Service("TFileService",
    fileName = cms.string(os.path.join("output",out_fname))
)

process.p = cms.Path(process.EFTGenReader)

# summary
process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(False),
    SkipEvent = cms.untracked.vstring('ProductNotFound')
)
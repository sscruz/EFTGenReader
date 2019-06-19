import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing
import sys
import os

from EFTGenReader.Reader.DatasetHelper import DatasetHelper

options = VarParsing.VarParsing('analysis')
options.maxEvents = 100000#-1
#options.maxEvents = 10

nd_redirect = "root://ndcms.crc.nd.edu/"
fnal_redirect = "root://cmsxrootd.fnal.gov/"
global_redirect = "root://cms-xrd-global.cern.ch/"

cmssw_base_dir = os.environ['CMSSW_BASE']
dataset_fpath = os.path.join(cmssw_base_dir,"src/EFTGenReader/Reader/data/JSON/datasets.json")

ds_helper = DatasetHelper()
ds_helper.load(dataset_fpath)
ds_helper.root_redirect = nd_redirect

from Configuration.StandardSequences.Eras import eras
process = cms.Process("Demo", eras.Run2_2017)

process.load('FWCore.MessageService.MessageLogger_cfi')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents) # number of events
)

process.MessageLogger.cerr.FwkReport.reportEvery = 500

# central datasets
#   'central_ttZ'
#   'central_ttW'
#   'central_tZq'

# private ttll datasets
#   'ttll_FP_R4B9'
#   'ttllNoHiggs_SM'
#   'ttllNoHiggs_EFT'
#   'ttllnunuNoHiggs_SM'

# private tllq datasets
#   'tllq_FR_R4B9'
#   'tllq4f_SMNoSchanW'
#   'tllq4fMatched_SM'
#   'tllq4fMatched_EFT'
#   'tllq4fNoHiggs_SM'
#   'tllq4fNoHiggs_EFT'

# private ttlnu datasets
#   'ttlnu_FP_R4B9'
#   'ttlnu_EFT'
#   'ttlnuJet_EFT'
#   'ttlnu_NoPDFWeights'

ds_name = 'ttllnunuNoHiggs_SM'

files     = ds_helper.getFiles(ds_name)
is_eft    = ds_helper.getData(ds_name,'is_eft')
xsec_norm = ds_helper.getData(ds_name,'central_xsec')

out_fname = "%s_NoTopLeptons_output_tree.root" % (ds_name)
#out_fname = "TEST_output_tree.root"

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        *files
    )
)

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
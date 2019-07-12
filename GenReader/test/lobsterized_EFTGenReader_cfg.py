import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing
import sys
import os

options = VarParsing.VarParsing('analysis')

# Setup and register default options
options.maxEvents = 100
options.register("dataset","central_ttH",
    VarParsing.VarParsing.multiplicity.singleton,
    VarParsing.VarParsing.varType.string,"name of the dataset as it appears in the JSON file")
options.register("test",False,
    VarParsing.VarParsing.multiplicity.singleton,
    VarParsing.VarParsing.varType.bool, "changes the output name to a dummy value")
options.register("debug",False,
    VarParsing.VarParsing.multiplicity.singleton,
    VarParsing.VarParsing.varType.bool, "run in debug mode")
options.register("normType",1,
    VarParsing.VarParsing.multiplicity.singleton,
    VarParsing.VarParsing.varType.int,"how to normalize the histograms; 0 - no norm, 1 - unit norm (default), 2 - xsec norm")
options.register("intgLumi",1.0,
    VarParsing.VarParsing.multiplicity.singleton,
    VarParsing.VarParsing.varType.float,"intg. lumi to scale the histograms to (no effect for unit norm mode)")
options.register("fnSuffix","_NoTopLeptons_output_tree",
    VarParsing.VarParsing.multiplicity.singleton,
    VarParsing.VarParsing.varType.string,"string to append to the end of the output root file")
options.register("minPtJet",-1.0,
    VarParsing.VarParsing.multiplicity.singleton,
    VarParsing.VarParsing.varType.float,"max pt cut for genjets")
options.register("maxEtaJet",-1.0,
    VarParsing.VarParsing.multiplicity.singleton,
    VarParsing.VarParsing.varType.float,"max eta cut for genjets (-1 means no cut)")
options.register("minPtLep",-1.0,
    VarParsing.VarParsing.multiplicity.singleton,
    VarParsing.VarParsing.varType.float,"max pt cut for genleptons")
options.register("maxEtaLep",-1.0,
    VarParsing.VarParsing.multiplicity.singleton,
    VarParsing.VarParsing.varType.float,"max eta cut for genleptons (-1 means no cut)")

# Get and parse the command line arguments
options.parseArguments()

nd_redirect = "root://ndcms.crc.nd.edu/"
fnal_redirect = "root://cmsxrootd.fnal.gov/"
global_redirect = "root://cms-xrd-global.cern.ch/"

cmssw_base_dir = os.environ['CMSSW_BASE']

# TODO: Add check for if the ds sample exists and error if not

from Configuration.StandardSequences.Eras import eras
process = cms.Process("Demo", eras.Run2_2017)

process.load('FWCore.MessageService.MessageLogger_cfi')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents) # number of events
)

process.MessageLogger.cerr.FwkReport.reportEvery = 500

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        'file:///hadoop/store/user/awightma/postLHE_step/2019_04_19/ValidationHanModelPlusJet/v1/mAOD_step_ttHJet_HanModelFixedStartingPt_run0/HIG-RunIIFall17MiniAOD-00821ND_10257.root'
    )
)

process.load("EFTGenReader.GenReader.EFTGenReader_cfi")

process.EFTGenReader.debug     = options.debug
process.EFTGenReader.norm_type = options.normType      # 0 - No norm, 1 - unit norm, 2 - xsec norm
process.EFTGenReader.intg_lumi = options.intgLumi
process.EFTGenReader.gp_events = 500

# Cut settings
process.EFTGenReader.min_pt_jet = options.minPtJet
process.EFTGenReader.min_pt_lep = options.minPtLep
process.EFTGenReader.max_eta_jet = options.maxEtaJet
process.EFTGenReader.max_eta_lep = options.maxEtaLep

process.EFTGenReader.iseft = True

process.TFileService = cms.Service("TFileService",
    fileName = cms.string("output_tree.root")
)

process.p = cms.Path(process.EFTGenReader)

# summary
process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(False),
    SkipEvent = cms.untracked.vstring('ProductNotFound')
)
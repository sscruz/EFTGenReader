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
options.register("datatier","GEN",
    VarParsing.VarParsing.multiplicity.singleton,
    VarParsing.VarParsing.varType.string,"EDM datatier format of the input root files")
options.register("test",False,
    VarParsing.VarParsing.multiplicity.singleton,
    VarParsing.VarParsing.varType.bool, "changes the output name to a dummy value")
options.register("debug",False,
    VarParsing.VarParsing.multiplicity.singleton,
    VarParsing.VarParsing.varType.bool, "run in debug mode")
options.register("iseft",True,
    VarParsing.VarParsing.multiplicity.singleton,
    VarParsing.VarParsing.varType.bool, "the sample has EFT reweighting")
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

process.load("EFTGenReader.GenReader.EFTGenHistsWithCuts_cfi")

process.EFTGenHistsWithCuts.debug     = options.debug
process.EFTGenHistsWithCuts.norm_type = options.normType      # 0 - No norm, 1 - unit norm, 2 - xsec norm
process.EFTGenHistsWithCuts.intg_lumi = options.intgLumi
process.EFTGenHistsWithCuts.gp_events = 500

# Cut settings
process.EFTGenHistsWithCuts.min_pt_jet = options.minPtJet
process.EFTGenHistsWithCuts.min_pt_lep = options.minPtLep
process.EFTGenHistsWithCuts.max_eta_jet = options.maxEtaJet
process.EFTGenHistsWithCuts.max_eta_lep = options.maxEtaLep

process.EFTGenHistsWithCuts.iseft = options.iseft

if options.datatier == "GEN":
    process.EFTGenHistsWithCuts.GenParticles = cms.InputTag("genParticles")
    process.EFTGenHistsWithCuts.GenJets      = cms.InputTag("ak4GenJets")
elif options.datatier == "MINIAODSIM": # Note: The particel level stuff can ONLY be run over GEN samples, so this mini AOD option is currently broken
    process.EFTGenHistsWithCuts.GenParticles = cms.InputTag("prunedGenParticles")
    process.EFTGenHistsWithCuts.GenJets      = cms.InputTag("slimmedGenJets")
else:
    print "[ERROR] Unknown datatier: {}".format(options.datatier)
    raise RuntimeError

process.TFileService = cms.Service("TFileService",
    fileName = cms.string("output_tree.root")
)

####### Particle level stuff from Reza ######
process.load("GeneratorInterface.RivetInterface.mergedGenParticles_cfi")
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load("GeneratorInterface.RivetInterface.genParticles2HepMC_cfi")
process.load("GeneratorInterface.RivetInterface.particleLevel_cfi")
from PhysicsTools.PatAlgos.slimming.prunedGenParticles_cfi import *
from PhysicsTools.PatAlgos.slimming.packedGenParticles_cfi import *
process.prunedGenParticlesOne = prunedGenParticles.clone()
process.prunedGenParticlesOne.select = cms.vstring( "keep    *")
process.prunedGenParticlesTwo = prunedGenParticles.clone()
process.prunedGenParticlesTwo.select = cms.vstring( "keep    *")
process.prunedGenParticlesTwo.src =  cms.InputTag("prunedGenParticlesOne")
process.mergedGenParticles.inputPruned = cms.InputTag("prunedGenParticlesTwo")
process.mergedGenParticles.inputPacked = cms.InputTag("mypackedGenParticles")

# Particle level jets and leptons
process.EFTGenHistsWithCuts.ParticleLevelJets = cms.InputTag("particleLevel","jets")
process.EFTGenHistsWithCuts.ParticleLevelLeptons = cms.InputTag("particleLevel","leptons")

process.mypackedGenParticles = packedGenParticles.clone()
process.mypackedGenParticles.inputCollection = cms.InputTag("prunedGenParticlesOne")
process.mypackedGenParticles.map = cms.InputTag("prunedGenParticlesTwo")
process.mypackedGenParticles.inputOriginal = cms.InputTag("genParticles")

process.p1 = cms.Path(
    process.prunedGenParticlesOne *
    process.prunedGenParticlesTwo *
    process.mypackedGenParticles *
    process.mergedGenParticles *
    process.genParticles2HepMC *
    process.particleLevel *
    process.EFTGenHistsWithCuts
    )
####### End particle level stuff from Reza ######
#process.p = cms.Path(process.EFTGenHistsWithCuts) # This line is here if we are not including the particle level stuff

# summary
process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(False),
    SkipEvent = cms.untracked.vstring('ProductNotFound')
)

####### Particle level stuff from Reza ######
process.out = cms.OutputModule(
    "PoolOutputModule",
    fileName = cms.untracked.string("EDMiii.root")
    )

process.outpath = cms.EndPath(process.out)
###### End particle level stuff from Reza ######

import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing
import sys
import os

from EFTGenReader.GenReader.DatasetHelper import DatasetHelper

options = VarParsing.VarParsing('analysis')

# Setup and register default options
options.maxEvents = -1
options.register("fnSuffix","_output_tree",
    VarParsing.VarParsing.multiplicity.singleton,
    VarParsing.VarParsing.varType.string,"string to append to the end of the output root file")
options.register("minPtJet",30,
    VarParsing.VarParsing.multiplicity.singleton,
    VarParsing.VarParsing.varType.float,"max pt cut for genjets")
options.register("maxEtaJet",2.4,
    VarParsing.VarParsing.multiplicity.singleton,
    VarParsing.VarParsing.varType.float,"max eta cut for genjets (-1 means no cut)")
options.register("minPtLep",20,
    VarParsing.VarParsing.multiplicity.singleton,
    VarParsing.VarParsing.varType.float,"max pt cut for genleptons")
options.register("maxEtaLep",2.4,
    VarParsing.VarParsing.multiplicity.singleton,
    VarParsing.VarParsing.varType.float,"max eta cut for genleptons (-1 means no cut)")

# Get and parse the command line arguments
options.parseArguments()


from Configuration.StandardSequences.Eras import eras
process = cms.Process("Demo", eras.Run2_2017)

process.load('FWCore.MessageService.MessageLogger_cfi')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents) # number of events
)

process.MessageLogger.cerr.FwkReport.reportEvery = 1

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        options.inputFiles
    )
)

process.load("EFTGenReader.GenReader.EFTGenHistsWithCuts_cfi")


# Cut settings
process.EFTGenHistsWithCuts.min_pt_jet = options.minPtJet
process.EFTGenHistsWithCuts.min_pt_lep = options.minPtLep
process.EFTGenHistsWithCuts.max_eta_jet = options.maxEtaJet
process.EFTGenHistsWithCuts.max_eta_lep = options.maxEtaLep

process.EFTGenHistsWithCuts.iseft = True
process.EFTGenHistsWithCuts.xsec_norm = True

process.EFTGenHistsWithCuts.GenParticles = cms.InputTag("genParticles")
process.EFTGenHistsWithCuts.GenJets      = cms.InputTag("ak4GenJets")

process.TFileService = cms.Service("TFileService",
                                   fileName = cms.string('output.root')
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

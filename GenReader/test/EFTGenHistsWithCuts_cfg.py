import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing
import sys
import os

from EFTGenReader.GenReader.DatasetHelper import DatasetHelper

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
options.register("fnSuffix","_output_tree",
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
dataset_fpath = os.path.join(cmssw_base_dir,"src/EFTGenReader/GenReader/data/JSON/datasets.json")

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
#   'central_ttH'

# private ttll datasets
#   'ttll_FP_R4B9'
#   'ttll_SM'
#   'ttllNoHiggs_SM'
#   'ttllNoHiggs_EFT'
#   'ttllnunuNoHiggs_SM'

# private tllq datasets
#   'tllq_FR_R4B9'
#   'tllq_SM'
#   'tllq4f_SMNoSchanW'
#   'tllq4fMatched_SM'
#   'tllq4fMatched_EFT'
#   'tllq4fNoHiggs_SM'
#   'tllq4fNoHiggs_EFT'

# private ttlnu datasets
#   'ttlnu_FP_R4B9'
#   'ttlnu_SM'
#   'ttlnu_EFT'     --> Not really EFT and is probably just an equiv re-running of 'ttlnu_SM'
#   'ttlnuJet_EFT'  --> Also probably doesn't have any EFT diagrams
#   'ttlnu_NoPDFWeights'

# private ttH datasets
#   'ttH_SM'

ds_name   = options.dataset
files     = ds_helper.getFiles(ds_name)
is_eft    = ds_helper.getData(ds_name,'is_eft')
xsec_norm = ds_helper.getData(ds_name,'central_xsec')
datatier  = ds_helper.getData(ds_name,'datatier')

#out_fname = "%s_NoTopLeptons_output_tree.root" % (ds_name)
out_fname = "%s%s.root" % (ds_name,options.fnSuffix)
if options.test:
    out_fname = "TEST_output_tree.root"
out_path = os.path.join("output",out_fname)

print "Using Sample: %s" % (ds_name)
print "Save output to: %s" % (out_path)

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        *files
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

process.EFTGenHistsWithCuts.iseft = is_eft
process.EFTGenHistsWithCuts.xsec_norm = xsec_norm

if datatier == "GEN":
    process.EFTGenHistsWithCuts.GenParticles = cms.InputTag("genParticles")
    process.EFTGenHistsWithCuts.GenJets      = cms.InputTag("ak4GenJets")
elif datatier == "MINIAODSIM":
    process.EFTGenHistsWithCuts.GenParticles = cms.InputTag("prunedGenParticles")
    process.EFTGenHistsWithCuts.GenJets      = cms.InputTag("slimmedGenJets")
else:
    print "[ERROR] Unknown datatier: {}".format(datatier)
    raise RuntimeError

process.TFileService = cms.Service("TFileService",
    fileName = cms.string(out_path)
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

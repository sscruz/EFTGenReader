import FWCore.ParameterSet.Config as cms
import FWCore.ParameterSet.VarParsing as VarParsing
import sys
import os

options = VarParsing.VarParsing('analysis')
#options.maxEvents = -1
options.maxEvents = 500

from Configuration.StandardSequences.Eras import eras
process = cms.Process("Demo", eras.Run2_2017)

process.load('FWCore.MessageService.MessageLogger_cfi')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(options.maxEvents) # number of events
)

process.MessageLogger.cerr.FwkReport.reportEvery = 100

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        #"file:///hadoop/store/user/awightma/LHE_step/2018_04_17/v1/lhe_step_ttH_cblS1_run0/HIG-RunIIFall17wmLHE-00000ND_20123.root",
        #"file:///hadoop/store/user/awightma/LHE_step/2018_04_17/v1/lhe_step_ttH_cblS1_run0/HIG-RunIIFall17wmLHE-00000ND_20443.root"
        #"file:///hadoop/store/user/awightma/LHE_step/2019_04_19/ValidationHanModelPlusJet/v1/lhe_step_ttHJet_HanModel16DttllScanpoints_run1/HIG-RunIIFall17wmLHE-00000ND_7021.root",
        #"file:///hadoop/store/user/awightma/LHE_step/2019_04_19/ValidationHanModelPlusJet/v1/lhe_step_ttHJet_HanModel16DttllScanpoints_run1/HIG-RunIIFall17wmLHE-00000ND_7022.root",
        #"file:///hadoop/store/user/awightma/LHE_step/2019_04_19/ValidationHanModelPlusJet/v1/lhe_step_ttHJet_HanModel16DttllScanpoints_run1/HIG-RunIIFall17wmLHE-00000ND_7023.root"
        #"file:///hadoop/store/user/awightma/postLHE_step/2019_04_19/ValidationHanModelPlusJet/v1/mAOD_step_ttHJet_HanModel16DttllScanpoints_run1/HIG-RunIIFall17MiniAOD-00821ND_74026.root"
        "file:///hadoop/store/user/kmohrman/genOnly_step/2019_04_19/ttHJet-xqcutStudies-xqcut10qCutTests/v1/gen_step_ttHJet_HanModel16DttllScanpointsqCut19_run1/GEN-00000_967.root"
    )
)

process.load("EFTGenReader.LHEReader.EFTLHEReader_cfi")
process.EFTLHEReader.min_pt_jet  = cms.double(-1)
process.EFTLHEReader.min_pt_lep  = cms.double(-1)
process.EFTLHEReader.max_eta_jet = cms.double(2.5)
process.EFTLHEReader.max_eta_lep = cms.double(2.5)

process.TFileService = cms.Service("TFileService",
    fileName = cms.string("output_tree.root")
)

process.p = cms.Path(process.EFTLHEReader)

# summary
process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(False),
    SkipEvent = cms.untracked.vstring('ProductNotFound')
)

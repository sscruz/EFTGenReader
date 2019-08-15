import FWCore.ParameterSet.Config as cms

EFTLHEReader = cms.EDAnalyzer("EFTLHEReader",
    min_pt_jet   = cms.double(-1.0),
    min_pt_lep   = cms.double(-1.0),
    max_eta_jet  = cms.double(-1.0),
    max_eta_lep  = cms.double(-1.0),

    LHEInfo      = cms.InputTag("externalLHEProducer"),
    GENInfo      = cms.InputTag("generator"),
    GenParticles = cms.InputTag(""),
    GenJets      = cms.InputTag("")
)

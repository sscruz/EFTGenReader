import FWCore.ParameterSet.Config as cms

EFTGenHistsWithCuts = cms.EDAnalyzer("EFTGenHistsWithCuts",
    debug = cms.bool(False),
    iseft = cms.bool(False),
    gp_events = cms.int32(500),
    norm_type = cms.int32(0), # 0 - No norm, 1 - unit norm, 2 - xsec norm
    xsec_norm = cms.double(1.0),
    intg_lumi = cms.double(1.0),

    # Designed to be as close as possible to TOP-19-001
    min_pt_jet = cms.double(30),
    min_pt_lep = cms.double(10),
    max_eta_jet = cms.double(2.4),
    max_eta_lep = cms.double(2.4),
    staggered_pt_cuts_lep = cms.vdouble(25.0,15.0,10.0,10.0),
    min_njets_2lss = cms.int32(4),
    min_njets_3l = cms.int32(2),
    min_njets_4l = cms.int32(2),
    max_njet_bins_2lss = cms.int32(7),
    max_njet_bins_3l = cms.int32(5),
    max_njet_bins_4l = cms.int32(4),

    LHEInfo      = cms.InputTag("externalLHEProducer"),
    GENInfo      = cms.InputTag("generator"),
    GenParticles = cms.InputTag(""),
    GenJets      = cms.InputTag("")
)

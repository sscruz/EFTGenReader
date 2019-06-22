import FWCore.ParameterSet.Config as cms

EFTGenReader = cms.EDAnalyzer("EFTGenReader",
    debug = cms.bool(False),
    iseft = cms.bool(False),
    gp_events = cms.int32(500),
    norm_type = cms.int32(0), # 0 - No norm, 1 - unit norm, 2 - xsec norm
    xsec_norm = cms.double(1.0),
    intg_lumi = cms.double(1.0),
    min_pt_jet = cms.double(-1.0),
    min_pt_lep = cms.double(-1.0),
    max_eta_jet = cms.double(-1.0),
    max_eta_lep = cms.double(-1.0)
)
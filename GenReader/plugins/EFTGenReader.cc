// Created by Andrew Wightman

//#include "EFTGenReader/GenReader/interface/EFTGenReader.h"
#include "EFTGenReader.h"

EFTGenReader::EFTGenReader(const edm::ParameterSet& iConfig)
{
    debug = iConfig.getParameter<bool> ("debug");
    iseft = iConfig.getParameter<bool> ("iseft");

    min_pt_jet  = iConfig.getParameter<double> ("min_pt_jet");
    min_pt_lep  = iConfig.getParameter<double> ("min_pt_lep");
    max_eta_jet = iConfig.getParameter<double> ("max_eta_jet");
    max_eta_lep = iConfig.getParameter<double> ("max_eta_lep");
    
    parse_params(); // Currently doesn't do anything
    lheInfo_token_      = consumes<LHEEventProduct>(iConfig.getParameter<edm::InputTag>("LHEInfo"));
    genInfo_token_      = consumes<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("GENInfo"));
    genParticles_token_ = consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("GenParticles"));
    genJets_token_      = consumes<std::vector<reco::GenJet> >(iConfig.getParameter<edm::InputTag>("GenJets"));
}

EFTGenReader::~EFTGenReader(){}

void EFTGenReader::beginJob()
{
    total_ls = 0;
    total_events = 0;
    total_orig_xsec = 0.;
    total_sm_xsec = 0.;

    edm::Service<TFileService> newfs;

    int pdg_bins = 100;
    int njet_bins = 16;
    int pt_bins = 5;
    int eta_bins = 10;
    int invmass_bins = 30;
    int deltaR_bins = 10;

    // Book the histograms that we will fill in the event loop

    //h_lep_ptEFT = newfs->make<TH1EFT>("h_lep_ptEFT","h_lep_ptEFT",pt_bins,0,300);

    h_eventsumEFT = newfs->make<TH1EFT>("h_eventsumEFT","h_eventsumEFT",1,0,1);

    // lep pt
    h_lep_ptEFT = newfs->make<TH1EFT>("h_lep_ptEFT","h_lep_ptEFT",pt_bins,0,300);
    h_lep_ptSM = newfs->make<TH1D>("h_lep_ptSM","h_lep_ptSM",pt_bins,0,300);
    h_e_ptEFT = newfs->make<TH1EFT>("h_e_ptEFT","h_e_ptEFT",pt_bins,0,300);
    h_e_ptSM = newfs->make<TH1D>("h_e_ptSM","h_e_ptSM",pt_bins,0,300);
    h_mu_ptEFT = newfs->make<TH1EFT>("h_mu_ptEFT","h_mu_ptEFT",pt_bins,0,300);
    h_mu_ptSM = newfs->make<TH1D>("h_mu_ptSM","h_mu_ptSM",pt_bins,0,300);
    h_tau_ptEFT = newfs->make<TH1EFT>("h_tau_ptEFT","h_tau_ptEFT",pt_bins,0,300);
    h_tau_ptSM = newfs->make<TH1D>("h_tau_ptSM","h_tau_ptSM",pt_bins,0,300);
    
    h_lep1_ptEFT = newfs->make<TH1EFT>("h_lep1_ptEFT","h_lep1_ptEFT",pt_bins,0,300);
    h_lep1_ptSM = newfs->make<TH1D>("h_lep1_ptSM","h_lep1_ptSM",pt_bins,0,300);
    h_e1_ptEFT = newfs->make<TH1EFT>("h_e1_ptEFT","h_e1_ptEFT",pt_bins,0,300);
    h_e1_ptSM = newfs->make<TH1D>("h_e1_ptSM","h_e1_ptSM",pt_bins,0,300);
    h_mu1_ptEFT = newfs->make<TH1EFT>("h_mu1_ptEFT","h_mu1_ptEFT",pt_bins,0,300);
    h_mu1_ptSM = newfs->make<TH1D>("h_mu1_ptSM","h_mu1_ptSM",pt_bins,0,300);
    h_tau1_ptEFT = newfs->make<TH1EFT>("h_tau1_ptEFT","h_tau1_ptEFT",pt_bins,0,300);
    h_tau1_ptSM = newfs->make<TH1D>("h_tau1_ptSM","h_tau1_ptSM",pt_bins,0,300);
    
    h_lep2_ptEFT = newfs->make<TH1EFT>("h_lep2_ptEFT","h_lep2_ptEFT",pt_bins,0,300);
    h_lep2_ptSM = newfs->make<TH1D>("h_lep2_ptSM","h_lep2_ptSM",pt_bins,0,300);
    h_e2_ptEFT = newfs->make<TH1EFT>("h_e2_ptEFT","h_e2_ptEFT",pt_bins,0,300);
    h_e2_ptSM = newfs->make<TH1D>("h_e2_ptSM","h_e2_ptSM",pt_bins,0,300);
    h_mu2_ptEFT = newfs->make<TH1EFT>("h_mu2_ptEFT","h_mu2_ptEFT",pt_bins,0,300);
    h_mu2_ptSM = newfs->make<TH1D>("h_mu2_ptSM","h_mu2_ptSM",pt_bins,0,300);
    h_tau2_ptEFT = newfs->make<TH1EFT>("h_tau2_ptEFT","h_tau2_ptEFT",pt_bins,0,300);
    h_tau2_ptSM = newfs->make<TH1D>("h_tau2_ptSM","h_tau2_ptSM",pt_bins,0,300);
    
    h_lepSum_ptEFT = newfs->make<TH1EFT>("h_lepSum_ptEFT","h_lepSum_ptEFT",pt_bins,0,300);
    h_lepSum_ptSM = newfs->make<TH1D>("h_lepSum_ptSM","h_lepSum_ptSM",pt_bins,0,300);
    h_eSum_ptEFT = newfs->make<TH1EFT>("h_eSum_ptEFT","h_eSum_ptEFT",pt_bins,0,300);
    h_eSum_ptSM = newfs->make<TH1D>("h_eSum_ptSM","h_eSum_ptSM",pt_bins,0,300);
    h_muSum_ptEFT = newfs->make<TH1EFT>("h_muSum_ptEFT","h_muSum_ptEFT",pt_bins,0,300);
    h_muSum_ptSM = newfs->make<TH1D>("h_muSum_ptSM","h_muSum_ptSM",pt_bins,0,300);
    h_tauSum_ptEFT = newfs->make<TH1EFT>("h_tauSum_ptEFT","h_tauSum_ptEFT",pt_bins,0,300);
    h_tauSum_ptSM = newfs->make<TH1D>("h_tauSum_ptSM","h_tauSum_ptSM",pt_bins,0,300);

    // lep eta
    h_lep_etaEFT = newfs->make<TH1EFT>("h_lep_etaEFT","h_lep_etaEFT",eta_bins,-5.0,5.0);
    h_lep_etaSM = newfs->make<TH1D>("h_lep_etaSM","h_lep_etaSM",eta_bins,-5.0,5.0);
    h_e_etaEFT = newfs->make<TH1EFT>("h_e_etaEFT","h_e_etaEFT",eta_bins,-5.0,5.0);
    h_e_etaSM = newfs->make<TH1D>("h_e_etaSM","h_e_etaSM",eta_bins,-5.0,5.0);
    h_mu_etaEFT = newfs->make<TH1EFT>("h_mu_etaEFT","h_mu_etaEFT",eta_bins,-5.0,5.0);
    h_mu_etaSM = newfs->make<TH1D>("h_mu_etaSM","h_mu_etaSM",eta_bins,-5.0,5.0);
    h_tau_etaEFT = newfs->make<TH1EFT>("h_tau_etaEFT","h_tau_etaEFT",eta_bins,-5.0,5.0);
    h_tau_etaSM = newfs->make<TH1D>("h_tau_etaSM","h_tau_etaSM",eta_bins,-5.0,5.0);

    // dilep inv mass
    //h_mllEFT = newfs->make<TH1EFT>("h_mllEFT","h_mllEFT",pt_bins,0,300);
    //h_mllSM = newfs->make<TH1D>("h_mllSM","h_mllSM",pt_bins,0,300);
    h_mllEFT = newfs->make<TH1EFT>("h_mllEFT","h_mllEFT",invmass_bins,0,300);
    h_mllSM = newfs->make<TH1D>("h_mllSM","h_mllSM",invmass_bins,0,300);
    h_meeEFT = newfs->make<TH1EFT>("h_meeEFT","h_meeEFT",invmass_bins,0,300);
    h_meeSM = newfs->make<TH1D>("h_meeSM","h_meeSM",invmass_bins,0,300);
    h_mmumuEFT = newfs->make<TH1EFT>("h_mmumuEFT","h_mmumuEFT",invmass_bins,0,300);
    h_mmumuSM = newfs->make<TH1D>("h_mmumuSM","h_mmumuSM",invmass_bins,0,300);
    h_mtautauEFT = newfs->make<TH1EFT>("h_mtautauEFT","h_mtautauEFT",invmass_bins,0,300);
    h_mtautauSM = newfs->make<TH1D>("h_mtautauSM","h_mtautauSM",invmass_bins,0,300);

    // Jet histograms
    //h_nJetsEFT = newfs->make<TH1EFT>("h_njetsEFT","h_njetsEFT",16,0,15);
    //h_nJetsSM = newfs->make<TH1D>("h_njetsSM","h_njetsSM",16,0,15);
    h_nJetsEFT = newfs->make<TH1EFT>("h_njetsEFT","h_njetsEFT",njet_bins,0,njet_bins);
    h_nJetsSM = newfs->make<TH1D>("h_njetsSM","h_njetsSM",njet_bins,0,njet_bins);

    h_jet_ptEFT = newfs->make<TH1EFT>("h_jet_ptEFT","h_jet_ptEFT",pt_bins,0,300);
    h_jet_ptSM = newfs->make<TH1D>("h_jet_ptSM","h_jet_ptSM",pt_bins,0,300);

    h_jet1_ptEFT = newfs->make<TH1EFT>("h_jet1_ptEFT","h_jet1_ptEFT",pt_bins,0,300);
    h_jet1_ptSM = newfs->make<TH1D>("h_jet1_ptSM","h_jet1_ptSM",pt_bins,0,300);

    h_jet2_ptEFT = newfs->make<TH1EFT>("h_jet2_ptEFT","h_jet2_ptEFT",pt_bins,0,300);
    h_jet2_ptSM = newfs->make<TH1D>("h_jet2_ptSM","h_jet2_ptSM",pt_bins,0,300);

    h_jet3_ptEFT = newfs->make<TH1EFT>("h_jet3_ptEFT","h_jet2_ptEFT",pt_bins,0,300);
    h_jet3_ptSM = newfs->make<TH1D>("h_jet3_ptSM","h_jet2_ptSM",pt_bins,0,300);

    h_jet4_ptEFT = newfs->make<TH1EFT>("h_jet4_ptEFT","h_jet2_ptEFT",pt_bins,0,300);
    h_jet4_ptSM = newfs->make<TH1D>("h_jet4_ptSM","h_jet2_ptSM",pt_bins,0,300);


    // Jet Eta histograms
    h_jet_etaEFT = newfs->make<TH1EFT>("h_jet_etaEFT","h_jet_etaEFT",eta_bins,-5.0,5.0);
    h_jet_etaSM = newfs->make<TH1D>("h_jet_etaSM","h_jet_etaSM",eta_bins,-5.0,5.0);

    h_jet1_etaEFT = newfs->make<TH1EFT>("h_jet1_etaEFT","h_jet1_etaEFT",eta_bins,-5.0,5.0);
    h_jet1_etaSM = newfs->make<TH1D>("h_jet1_etaSM","h_jet1_etaSM",eta_bins,-5.0,5.0);

    h_jet2_etaEFT = newfs->make<TH1EFT>("h_jet2_etaEFT","h_jet2_etaEFT",eta_bins,-5.0,5.0);
    h_jet2_etaSM = newfs->make<TH1D>("h_jet2_etaSM","h_jet2_etaSM",eta_bins,-5.0,5.0);

    h_jet3_etaEFT = newfs->make<TH1EFT>("h_jet3_etaEFT","h_jet3_etaEFT",eta_bins,-5.0,5.0);
    h_jet3_etaSM = newfs->make<TH1D>("h_jet3_etaSM","h_jet3_etaSM",eta_bins,-5.0,5.0);

    h_jet4_etaEFT = newfs->make<TH1EFT>("h_jet4_etaEFT","h_jet4_etaEFT",eta_bins,-5.0,5.0);
    h_jet4_etaSM = newfs->make<TH1D>("h_jet4_etaSM","h_jet4_etaSM",eta_bins,-5.0,5.0);

    // pdgId histograms
    h_pdgIdEFT = newfs->make<TH1EFT>("h_pdgIdEFT","h_pdgIdEFT",2*pdg_bins,-pdg_bins+1,pdg_bins);
    h_pdgIdSM = newfs->make<TH1D>("h_pdgIdSM","h_pdgIdSM",2*pdg_bins,-pdg_bins+1,pdg_bins);

    h_pdgIdFromZEFT = newfs->make<TH1EFT>("h_pdgIdFromZEFT","h_pdgIdFromZEFT",2*pdg_bins,-pdg_bins+1,pdg_bins);
    h_pdgIdFromZSM = newfs->make<TH1D>("h_pdgIdFromZSM","h_pdgIdFromZSM",2*pdg_bins,-pdg_bins+1,pdg_bins);

    h_pdgIdFromHEFT = newfs->make<TH1EFT>("h_pdgIdFromHEFT","h_pdgIdFromHEFT",2*pdg_bins,-pdg_bins+1,pdg_bins);
    h_pdgIdFromHSM = newfs->make<TH1D>("h_pdgIdFromHSM","h_pdgIdFromHSM",2*pdg_bins,-pdg_bins+1,pdg_bins);

    h_pdgIdFromWEFT = newfs->make<TH1EFT>("h_pdgIdFromWEFT","h_pdgIdFromWEFT",2*pdg_bins,-pdg_bins+1,pdg_bins);
    h_pdgIdFromWSM = newfs->make<TH1D>("h_pdgIdFromWSM","h_pdgIdFromWSM",2*pdg_bins,-pdg_bins+1,pdg_bins);

    h_pdgIdFromPhotonEFT = newfs->make<TH1EFT>("h_pdgIdFromPhotonEFT","h_pdgIdFromPhotonEFT",2*pdg_bins,-pdg_bins+1,pdg_bins);
    h_pdgIdFromPhotonSM = newfs->make<TH1D>("h_pdgIdFromPhotonSM","h_pdgIdFromPhotonSM",2*pdg_bins,-pdg_bins+1,pdg_bins);

    h_pdgIdFromGluonEFT = newfs->make<TH1EFT>("h_pdgIdFromGluonEFT","h_pdgIdFromGluonEFT",2*pdg_bins,-pdg_bins+1,pdg_bins);
    h_pdgIdFromGluonSM = newfs->make<TH1D>("h_pdgIdFromGluonSM","h_pdgIdFromGluonSM",2*pdg_bins,-pdg_bins+1,pdg_bins);

    h_pdgIdFromQCDEFT = newfs->make<TH1EFT>("h_pdgIdFromQCDEFT","h_pdgIdFromQCDEFT",2*pdg_bins,-pdg_bins+1,pdg_bins);
    h_pdgIdFromQCDSM = newfs->make<TH1D>("h_pdgIdFromQCDSM","h_pdgIdFromQCDSM",2*pdg_bins,-pdg_bins+1,pdg_bins);

    h_pdgIdFromOtherEFT = newfs->make<TH1EFT>("h_pdgIdFromOtherEFT","h_pdgIdFromOtherEFT",2*pdg_bins,-pdg_bins+1,pdg_bins);
    h_pdgIdFromOtherSM = newfs->make<TH1D>("h_pdgIdFromOtherSM","h_pdgIdFromOtherSM",2*pdg_bins,-pdg_bins+1,pdg_bins);

    h_pdgIdLepMotherEFT = newfs->make<TH1EFT>("h_pdgIdLepMotherEFT","h_pdgIdLepMotherEFT",2*pdg_bins,-pdg_bins+1,pdg_bins);
    h_pdgIdLepMotherSM = newfs->make<TH1D>("h_pdgIdLepMotherSM","h_pdgIdLepMotherSM",2*pdg_bins,-pdg_bins+1,pdg_bins);
    h_pdgIdElectronMotherEFT = newfs->make<TH1EFT>("h_pdgIdElectronMotherEFT","h_pdgIdElectronMotherEFT",2*pdg_bins,-pdg_bins+1,pdg_bins);
    h_pdgIdElectronMotherSM = newfs->make<TH1D>("h_pdgIdElectronMotherSM","h_pdgIdElectronMotherSM",2*pdg_bins,-pdg_bins+1,pdg_bins);
    h_pdgIdMuonMotherEFT = newfs->make<TH1EFT>("h_pdgIdMuonMotherEFT","h_pdgIdMuonMotherEFT",2*pdg_bins,-pdg_bins+1,pdg_bins);
    h_pdgIdMuonMotherSM = newfs->make<TH1D>("h_pdgIdMuonMotherSM","h_pdgIdMuonMotherSM",2*pdg_bins,-pdg_bins+1,pdg_bins);
    h_pdgIdTauMotherEFT = newfs->make<TH1EFT>("h_pdgIdTauMotherEFT","h_pdgIdTauMotherEFT",2*pdg_bins,-pdg_bins+1,pdg_bins);
    h_pdgIdTauMotherSM = newfs->make<TH1D>("h_pdgIdTauMotherSM","h_pdgIdTauMotherSM",2*pdg_bins,-pdg_bins+1,pdg_bins);

    h_pdgIdLepGrMotherEFT = newfs->make<TH1EFT>("h_pdgIdLepGrMotherEFT","h_pdgIdLepGrMotherEFT",2*pdg_bins,-pdg_bins+1,pdg_bins);
    h_pdgIdLepGrMotherSM = newfs->make<TH1D>("h_pdgIdLepGrMotherSM","h_pdgIdLepGrMotherSM",2*pdg_bins,-pdg_bins+1,pdg_bins);
    h_pdgIdElectronGrMotherEFT = newfs->make<TH1EFT>("h_pdgIdElectronGrMotherEFT","h_pdgIdElectronGrMotherEFT",2*pdg_bins,-pdg_bins+1,pdg_bins);
    h_pdgIdElectronGrMotherSM = newfs->make<TH1D>("h_pdgIdElectronGrMotherSM","h_pdgIdElectronGrMotherSM",2*pdg_bins,-pdg_bins+1,pdg_bins);
    h_pdgIdMuonGrMotherEFT = newfs->make<TH1EFT>("h_pdgIdMuonGrMotherEFT","h_pdgIdMuonGrMotherEFT",2*pdg_bins,-pdg_bins+1,pdg_bins);
    h_pdgIdMuonGrMotherSM = newfs->make<TH1D>("h_pdgIdMuonGrMotherSM","h_pdgIdMuonGrMotherSM",2*pdg_bins,-pdg_bins+1,pdg_bins);
    h_pdgIdTauGrMotherEFT = newfs->make<TH1EFT>("h_pdgIdTauGrMotherEFT","h_pdgIdTauGrMotherEFT",2*pdg_bins,-pdg_bins+1,pdg_bins);
    h_pdgIdTauGrMotherSM = newfs->make<TH1D>("h_pdgIdTauGrMotherSM","h_pdgIdTauGrMotherSM",2*pdg_bins,-pdg_bins+1,pdg_bins);

    // misc. observables
    h_deltaREFT    = newfs->make<TH1EFT>("h_deltaREFT","h_deltaREFT",deltaR_bins,0,5);
    h_deltaRSM     = newfs->make<TH1D>("h_deltaRSM","h_deltaRSM",deltaR_bins,0,5);
    // NOTE: Should re write how these dR hists are filled to only fill with min dR; till then, leave commented
    //h_e_deltaREFT  = newfs->make<TH1EFT>("h_e_deltaREFT","h_e_deltaREFT",deltaR_bins,0,5);
    //h_e_deltaRSM   = newfs->make<TH1D>("h_e_deltaRSM","h_e_deltaRSM",deltaR_bins,0,5);
    //h_mu_deltaREFT = newfs->make<TH1EFT>("h_mu_deltaREFT","h_mu_deltaREFT",deltaR_bins,0,5);
    //h_mu_deltaRSM  = newfs->make<TH1D>("h_mu_deltaRSM","h_mu_deltaRSM",deltaR_bins,0,5);
    //h_tau_deltaREFT = newfs->make<TH1EFT>("h_tau_deltaREFT","h_tau_deltaREFT",deltaR_bins,0,5);
    //h_tau_deltaRSM  = newfs->make<TH1D>("h_tau_deltaRSM","h_tau_deltaRSM",deltaR_bins,0,5);
    
    h_prompt_leptonsEFT   = newfs->make<TH1EFT>("h_prompt_leptonsEFT","h_prompt_leptonsEFT",30,0,11);
    h_prompt_leptonsSM    = newfs->make<TH1D>("h_prompt_leptonsSM","h_prompt_leptonsSM",30,0,11);
    h_prompt_electronsEFT = newfs->make<TH1EFT>("h_prompt_electronsEFT","h_prompt_electronsEFT",30,0,11);
    h_prompt_electronsSM  = newfs->make<TH1D>("h_prompt_electronsSM","h_prompt_electronsSM",30,0,11);
    h_prompt_muonsEFT     = newfs->make<TH1EFT>("h_prompt_muonsEFT","h_prompt_muonsEFT",30,0,11);
    h_prompt_muonsSM      = newfs->make<TH1D>("h_prompt_muonsSM","h_prompt_muonsSM",30,0,11);
    h_prompt_tausEFT      = newfs->make<TH1EFT>("h_prompt_tausEFT","h_prompt_tausEFT",30,0,11);
    h_prompt_tausSM       = newfs->make<TH1D>("h_prompt_tausSM","h_prompt_tausSM",30,0,11);

    // Don't normalize these plots
    //h_SMwgt_norm = newfs->make<TH1D>("h_SMwgt_norm","h_SMwgt_norm",350,-0.1,2.0);
    h_SMwgt_norm = newfs->make<TH1D>("h_SMwgt_norm","h_SMwgt_norm",350,-8,1); // Log x scale
    binLogX(h_SMwgt_norm);

    summaryTree = newfs->make<TTree>("summaryTree","Summary Event Values");
    tree_add_branches();
}

void EFTGenReader::endJob()
{
    std::string delim = " | ";

    std::cout << "Total events processed: " << total_events << std::endl;
    std::cout << "Total LS: " << total_ls << std::endl;
    std::cout << "                " 
        << std::setw(13) << "xsec / LS" << delim
        << std::setw(13) << "xsec / Events" << delim
        << std::setw(13) << "xsec / 500"
        << std::endl;
    std::cout << "SM Xsec (pb):   "
        << std::setw(13) << total_sm_xsec / double(total_ls) << delim
        << std::setw(13) << total_sm_xsec / double(total_events) << delim
        << std::setw(13) << total_sm_xsec / 500.
        << std::endl;
    std::cout << "Orig Xsec (pb): "
        << std::setw(11) << total_orig_xsec / double(total_ls) << delim
        << std::setw(11) << total_orig_xsec / double(total_events) << delim
        << std::setw(11) << total_orig_xsec / 500.
        << std::endl;
}

void EFTGenReader::analyze(const edm::Event& event, const edm::EventSetup& evsetup)
{
    eventcount++;

    // set tree vars to default values
    initialize_variables();

    edm::Handle<LHEEventProduct> LHEInfo;
    edm::Handle<reco::GenParticleCollection> prunedParticles;
    edm::Handle<std::vector<reco::GenJet> > genJets;

    event.getByToken(lheInfo_token_,LHEInfo);
    event.getByToken(genParticles_token_,prunedParticles);
    event.getByToken(genJets_token_,genJets);

    if (debug){    
        if (eventcount == 3) dumpJets(*genJets);
    }

    reco::GenParticleCollection gen_leptons = GetGenLeptons(*prunedParticles);
    std::vector<reco::GenJet> gen_jets = GetGenJets(*genJets);

    originalXWGTUP_intree = LHEInfo->originalXWGTUP();  // original cross-section
    double sm_wgt = 0.;
    std::vector<WCPoint> wc_pts;
    if (iseft) {// Add EFT weights 
        for (auto wgt_info: LHEInfo->weights()) {
            auto LHEwgtstr = std::string(wgt_info.id);
            std::size_t foundstr = LHEwgtstr.find("EFTrwgt"); // only save our EFT weights
            if (foundstr!=std::string::npos) {
                WCPoint wc_pt(wgt_info.id,wgt_info.wgt);
                wc_pts.push_back(wc_pt);
                if (wc_pt.isSMPoint()) {
                    sm_wgt = wgt_info.wgt;
                }
            }
        }
    } else {
        sm_wgt = originalXWGTUP_intree;
        WCPoint wc_pt("smpt",sm_wgt);
        wc_pts.push_back(wc_pt);
    }

    WCFit eft_fit(wc_pts,"");

    if (debug) {
        for (uint i=0; i < wc_pts.size(); i++){
            WCPoint wc_pt = wc_pts.at(i);
            double pt_wgt = wc_pt.wgt;
            double fit_val = eft_fit.evalPoint(&wc_pt);
            std::cout << std::setw(3) << i << ": " << std::setw(12) << pt_wgt << " | " << std::setw(12) << fit_val << " | " << std::setw(12) << (pt_wgt-fit_val) << std::endl;
        }
    }

    h_eventsumEFT->Fill(0.5,1,eft_fit);

    total_sm_xsec += sm_wgt;
    total_orig_xsec += originalXWGTUP_intree;

    h_SMwgt_norm->Fill(sm_wgt);

    h_prompt_leptonsEFT->Fill(gen_leptons.size(),1.0,eft_fit);   h_prompt_leptonsSM->Fill(gen_leptons.size(),sm_wgt);


    double min_dR = -1;

    int n_prompt_electrons = 0;
    int n_prompt_muons = 0;
    int n_prompt_taus = 0;
    for (size_t i = 0; i < gen_leptons.size(); i++) {
        const reco::GenParticle& p1 = gen_leptons.at(i);
        int id = p1.pdgId();
        int st = p1.status();
        double pt = p1.p4().Pt();
        double eta = p1.p4().Eta();

        int direct_id = p1.mother()->pdgId();
        const reco::Candidate* p_mom = GetGenMotherNoFsr(&p1);
        const reco::Candidate* p_gmom = GetGenMotherNoFsr(p_mom);
        int mom_id = p_mom->pdgId();
        int gmom_id = p_gmom->pdgId();

        h_pdgIdLepGrMotherEFT->Fill(gmom_id,1.0,eft_fit); h_pdgIdLepGrMotherSM->Fill(gmom_id,sm_wgt);
        h_pdgIdLepMotherEFT->Fill(mom_id,1.0,eft_fit);    h_pdgIdLepMotherSM->Fill(mom_id,sm_wgt);
        h_lep_ptEFT->Fill(pt,1.0,eft_fit);                h_lep_ptSM->Fill(pt,sm_wgt);
        h_lep_etaEFT->Fill(eta,1.0,eft_fit);              h_lep_etaSM->Fill(eta,sm_wgt);
        if (abs(id) == 11) { // electron
            n_prompt_electrons += 1;
            h_pdgIdElectronGrMotherEFT->Fill(gmom_id,1.0,eft_fit); h_pdgIdElectronGrMotherSM->Fill(gmom_id,sm_wgt);
            h_pdgIdElectronMotherEFT->Fill(mom_id,1.0,eft_fit);    h_pdgIdElectronMotherSM->Fill(mom_id,sm_wgt);
            h_e_ptEFT->Fill(pt,1.0,eft_fit);                       h_e_ptSM->Fill(pt,sm_wgt);
            h_e_etaEFT->Fill(eta,1.0,eft_fit);                     h_e_etaSM->Fill(eta,sm_wgt);
        } else if (abs(id) == 13) { // muon
            n_prompt_muons += 1;
            h_pdgIdMuonGrMotherEFT->Fill(gmom_id,1.0,eft_fit);     h_pdgIdMuonGrMotherSM->Fill(gmom_id,sm_wgt);
            h_pdgIdMuonMotherEFT->Fill(mom_id,1.0,eft_fit);        h_pdgIdMuonMotherSM->Fill(mom_id,sm_wgt);
            h_mu_ptEFT->Fill(pt,1.0,eft_fit);                      h_mu_ptSM->Fill(pt,sm_wgt);
            h_mu_etaEFT->Fill(eta,1.0,eft_fit);                    h_mu_etaSM->Fill(eta,sm_wgt);
        //} else if (abs(id) == 13) { // tau
        } else if (abs(id) == 15) { // tau
            n_prompt_taus += 1;
            h_pdgIdTauGrMotherEFT->Fill(gmom_id,1.0,eft_fit);      h_pdgIdTauGrMotherSM->Fill(gmom_id,sm_wgt);
            h_pdgIdTauMotherEFT->Fill(mom_id,1.0,eft_fit);         h_pdgIdTauMotherSM->Fill(mom_id,sm_wgt);
            h_tau_ptEFT->Fill(pt,1.0,eft_fit);                     h_tau_ptSM->Fill(pt,sm_wgt);
            h_tau_etaEFT->Fill(eta,1.0,eft_fit);                   h_tau_etaSM->Fill(eta,sm_wgt);
        }

        // mom_id and direct_id should be the same thing for the gen leptons selected by GetGenLeptons()
        if (mom_id == 21) {
            h_pdgIdFromGluonEFT->Fill(id,1.0,eft_fit);
            h_pdgIdFromGluonSM->Fill(id,sm_wgt);
        } else if (mom_id >= 1 && mom_id <= 5) {
            h_pdgIdFromQCDEFT->Fill(id,1.0,eft_fit);
            h_pdgIdFromQCDSM->Fill(id,sm_wgt);  
        } else if (mom_id == 22) {
            h_pdgIdFromPhotonEFT->Fill(id,1.0,eft_fit);
            h_pdgIdFromPhotonSM->Fill(id,sm_wgt);
        } else if (mom_id == 23) {
            h_pdgIdFromZEFT->Fill(id,1.0,eft_fit);
            h_pdgIdFromZSM->Fill(id,sm_wgt);
        } else if (mom_id == -24 || mom_id == 24) {
            h_pdgIdFromWEFT->Fill(id,1.0,eft_fit);
            h_pdgIdFromWSM->Fill(id,sm_wgt);
        } else if (mom_id == 25) {
            h_pdgIdFromHEFT->Fill(id,1.0,eft_fit);
            h_pdgIdFromHSM->Fill(id,sm_wgt);
        } else {
            h_pdgIdFromOtherEFT->Fill(id,1.0,eft_fit);
            h_pdgIdFromOtherSM->Fill(id,sm_wgt);
        }

        if (i == 0) {
            h_lep1_ptEFT->Fill(pt,1.0,eft_fit);
            h_lep1_ptSM->Fill(pt,sm_wgt);
            if (abs(id) == 11) { // electron
                h_e1_ptEFT->Fill(pt,1.0,eft_fit);
                h_e1_ptSM->Fill(pt,sm_wgt);
            } else if (abs(id) == 13) {// muon
                h_mu1_ptEFT->Fill(pt,1.0,eft_fit);
                h_mu1_ptSM->Fill(pt,sm_wgt);
            } else if (abs(id) == 15) {// tau
                h_tau1_ptEFT->Fill(pt,1.0,eft_fit);
                h_tau1_ptSM->Fill(pt,sm_wgt);
            }
        } else if (i == 1) {
            h_lep2_ptEFT->Fill(pt,1.0,eft_fit);
            h_lep2_ptSM->Fill(pt,sm_wgt);
            if (abs(id) == 11) { // electron
                h_e2_ptEFT->Fill(pt,1.0,eft_fit);
                h_e2_ptSM->Fill(pt,sm_wgt);
            } else if (abs(id) == 13) { // muon
                h_mu2_ptEFT->Fill(pt,1.0,eft_fit);
                h_mu2_ptSM->Fill(pt,sm_wgt);
            } else if (abs(id) == 15) { // tau
                h_tau2_ptEFT->Fill(pt,1.0,eft_fit);
                h_tau2_ptSM->Fill(pt,sm_wgt);
            }
        }
        for (size_t j = 0; j < gen_leptons.size(); j++) {
            if (i >= j) continue;

            const reco::GenParticle& p2 = gen_leptons.at(j);
            double dR = getdR(p1,p2);
            double mll = getInvMass(p1,p2);
            auto p4vec = p1.p4() + p2.p4();
            double sum_pt = p4vec.Pt();

            int id1 = p1.pdgId();
            int id2 = p2.pdgId();

            if (min_dR < 0 or std::min(dR,min_dR) == dR){ // Add in this to fix the dR plots!!!
                min_dR = dR;
            }

            h_mllEFT->Fill(mll,1.0,eft_fit);              h_mllSM->Fill(mll,sm_wgt);
            //h_deltaREFT->Fill(dR,1.0,eft_fit);            h_deltaRSM->Fill(dR,sm_wgt);
            h_lepSum_ptEFT->Fill(sum_pt,1.0,eft_fit);     h_lepSum_ptSM->Fill(sum_pt,sm_wgt);
            if (abs(id1) == 11 and abs(id2) == 11){ // both are electrons
                h_meeEFT->Fill(mll,1.0,eft_fit);          h_meeSM->Fill(mll,sm_wgt);
                //h_e_deltaREFT->Fill(dR,1.0,eft_fit);      h_e_deltaRSM->Fill(dR,sm_wgt);
                h_eSum_ptEFT->Fill(sum_pt,1.0,eft_fit);   h_eSum_ptSM->Fill(sum_pt,sm_wgt);
            } else if (abs(id1) == 13 and abs(id2) == 13){ // both are muons
                h_mmumuEFT->Fill(mll,1.0,eft_fit);        h_mmumuSM->Fill(mll,sm_wgt);
                //h_mu_deltaREFT->Fill(dR,1.0,eft_fit);     h_mu_deltaRSM->Fill(dR,sm_wgt);
                h_muSum_ptEFT->Fill(sum_pt,1.0,eft_fit);  h_muSum_ptSM->Fill(sum_pt,sm_wgt);
            } else if (abs(id1) == 15 and abs(id2) == 15){ // both are taus
                h_mtautauEFT->Fill(mll,1.0,eft_fit);      h_mtautauSM->Fill(mll,sm_wgt);
                //h_tau_deltaREFT->Fill(dR,1.0,eft_fit);    h_tau_deltaRSM->Fill(dR,sm_wgt);
                h_tauSum_ptEFT->Fill(sum_pt,1.0,eft_fit); h_tauSum_ptSM->Fill(sum_pt,sm_wgt);
            }
        }
    }
    h_deltaREFT->Fill(min_dR,1.0,eft_fit); h_deltaRSM->Fill(min_dR,sm_wgt);

    h_prompt_electronsEFT->Fill(n_prompt_electrons,1.0,eft_fit); h_prompt_electronsSM->Fill(n_prompt_electrons,sm_wgt);
    h_prompt_muonsEFT->Fill(n_prompt_muons,1.0,eft_fit);         h_prompt_muonsSM->Fill(n_prompt_muons,sm_wgt);
    h_prompt_tausEFT->Fill(n_prompt_taus,1.0,eft_fit);           h_prompt_tausSM->Fill(n_prompt_taus,sm_wgt);

    h_nJetsEFT->Fill(gen_jets.size(),1.0,eft_fit);
    h_nJetsSM->Fill(gen_jets.size(),sm_wgt);

    for (size_t i = 0; i < gen_jets.size(); i++) {
        const reco::GenJet& p1 = gen_jets.at(i);
        double pt = p1.p4().Pt();
        double eta = p1.p4().Eta();

        h_jet_ptEFT->Fill(pt,1.0,eft_fit);   h_jet_ptSM->Fill(pt,sm_wgt);
        h_jet_etaEFT->Fill(eta,1.0,eft_fit); h_jet_etaSM->Fill(eta,sm_wgt);
        if (i == 0) {
            h_jet1_ptEFT->Fill(pt,1.0,eft_fit);
            h_jet1_ptSM->Fill(pt,sm_wgt);
            h_jet1_etaEFT->Fill(eta,1.0,eft_fit);
            h_jet1_etaSM->Fill(eta,sm_wgt);
        } else if (i == 1) {
            h_jet2_ptEFT->Fill(pt,1.0,eft_fit);
            h_jet2_ptSM->Fill(pt,sm_wgt);
            h_jet2_etaEFT->Fill(eta,1.0,eft_fit);
            h_jet2_etaSM->Fill(eta,sm_wgt);
        } else if (i == 2) {
            h_jet3_ptEFT->Fill(pt,1.0,eft_fit);
            h_jet3_ptSM->Fill(pt,sm_wgt);
            h_jet3_etaEFT->Fill(eta,1.0,eft_fit);
            h_jet3_etaSM->Fill(eta,sm_wgt);
        } else if (i == 3) {
            h_jet4_ptEFT->Fill(pt,1.0,eft_fit);
            h_jet4_ptSM->Fill(pt,sm_wgt);
            h_jet4_etaEFT->Fill(eta,1.0,eft_fit);
            h_jet4_etaSM->Fill(eta,sm_wgt);
        }
    }

    eventnum_intree = event.id().event();
    lumiBlock_intree = event.id().luminosityBlock();
    runNumber_intree = event.id().run();

}

void EFTGenReader::beginRun(edm::Run const& run, edm::EventSetup const& evsetup)
{
    eventcount = 0;
}

void EFTGenReader::endRun(edm::Run const& run, edm::EventSetup const& evsetup)
{
    total_events += eventcount;
}

void EFTGenReader::beginLuminosityBlock(edm::LuminosityBlock const& lumi, edm::EventSetup const& evsetup)
{
    total_ls += 1;
}

void EFTGenReader::endLuminosityBlock(edm::LuminosityBlock const& lumi, edm::EventSetup const& evsetup){}

DEFINE_FWK_MODULE(EFTGenReader);


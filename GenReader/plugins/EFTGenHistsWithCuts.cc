
#include "EFTGenHistsWithCuts.h"
using namespace reco;

EFTGenHistsWithCuts::EFTGenHistsWithCuts(const edm::ParameterSet& iConfig)
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

    // Particle level stuff
    particleLevelJetsToken_ = consumes<std::vector<reco::GenJet>>(iConfig.getParameter<edm::InputTag>("ParticleLevelJets"));
    particleLevelLeptonsToken_ = consumes<std::vector<reco::GenJet>>(iConfig.getParameter<edm::InputTag>("ParticleLevelLeptons"));
    
}

EFTGenHistsWithCuts::~EFTGenHistsWithCuts(){}


void EFTGenHistsWithCuts::beginJob()
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
    h_eventsumEFT = newfs->make<TH1EFT>("h_eventsumEFT","h_eventsumEFT",1,0,1);

    // Jet histograms
    h_nJetsEFT= newfs->make<TH1EFT>("h_njetsEFT","h_njetsEFT",njet_bins,0,njet_bins);
    h_nJetsSM = newfs->make<TH1D>("h_njetsSM","h_njetsSM",njet_bins,0,njet_bins);

    // These are njets histograms with some basic selections
    h_nJets_3orMoreLep_EFT = newfs->make<TH1EFT>("h_nJets_3orMoreLep_EFT","h_nJets_3orMoreLep_EFT",njet_bins,0,njet_bins);
    h_nJets_3orMoreLep_SM  = newfs->make<TH1D>("h_nJets_3orMoreLep_SM","h_nJets_3orMoreLep_SM",njet_bins,0,njet_bins);
    h_nJets_3Lep_EFT       = newfs->make<TH1EFT>("h_nJets_3Lep_EFT","h_nJets_3Lep_EFT",njet_bins,0,njet_bins);
    h_nJets_3Lep_SM        = newfs->make<TH1D>("h_nJets_3Lep_SM","h_nJets_3Lep_SM",njet_bins,0,njet_bins);
    h_nJets_cleanJets_3orMoreLep_EFT = newfs->make<TH1EFT>("h_nJets_cleanJets_3orMoreLep_EFT","h_nJets_cleanJets_3orMoreLep_EFT",njet_bins,0,njet_bins);
    h_nJets_cleanJets_3orMoreLep_SM  = newfs->make<TH1D>("h_nJets_cleanJets_3orMoreLep_SM","h_nJets_cleanJets_3orMoreLep_SM",njet_bins,0,njet_bins);
    h_nJets_cleanJets_3Lep_EFT       = newfs->make<TH1EFT>("h_nJets_cleanJets_3Lep_EFT","h_nJets_cleanJets_3Lep_EFT",njet_bins,0,njet_bins);
    h_nJets_cleanJets_3Lep_SM        = newfs->make<TH1D>("h_nJets_cleanJets_3Lep_SM","h_nJets_cleanJets_3Lep_SM",njet_bins,0,njet_bins);

    // Particle level hists
    h_pl_nJets_EFT = newfs->make<TH1EFT>("h_pl_njetsEFT","h_pl_njetsEFT",njet_bins,0,njet_bins);
    h_pl_nJets_SM  = newfs->make<TH1D>("h_pl_njetsSM","h_pl_njetsSM",njet_bins,0,njet_bins);
    h_pl_nJets_3Lep_EFT = newfs->make<TH1EFT>("h_pl_nJets_3Lep_EFT","h_pl_nJets_3Lep_EFT",njet_bins,0,njet_bins);
    h_pl_nJets_3Lep_SM  = newfs->make<TH1D>("h_pl_nJets_3Lep_SM","h_pl_nJets_3Lep_SM",njet_bins,0,njet_bins);

    // Don't normalize these plots
    h_SMwgt_norm = newfs->make<TH1D>("h_SMwgt_norm","h_SMwgt_norm",350,-0.1,2.0);
    h_SMwgt_norm = newfs->make<TH1D>("h_SMwgt_norm","h_SMwgt_norm",350,-8,1); // Log x scale
    binLogX(h_SMwgt_norm);

    summaryTree = newfs->make<TTree>("summaryTree","Summary Event Values");
    tree_add_branches();
}

void EFTGenHistsWithCuts::endJob()
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

void EFTGenHistsWithCuts::analyze(const edm::Event& event, const edm::EventSetup& evsetup)
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

    reco::GenParticleCollection gen_leptons = GetGenLeptons(*prunedParticles);
    reco::GenParticleCollection gen_b = GetGenParticlesSubset(*prunedParticles, 5);
    std::vector<reco::GenJet> gen_jets = GetGenJets(*genJets);

    // Particle level stuff:
    edm::Handle<std::vector<reco::GenJet>> particleLevelJetsHandle_;
    edm::Handle<std::vector<reco::GenJet>> particleLevelLeptonsHandle_;
    event.getByToken(particleLevelJetsToken_,particleLevelJetsHandle_);
    event.getByToken(particleLevelLeptonsToken_,particleLevelLeptonsHandle_);
    std::vector<reco::GenJet> pl_jets    = MakePtEtaCuts(*particleLevelJetsHandle_,"jet");
    std::vector<reco::GenJet> pl_leptons = MakePtEtaCuts(*particleLevelLeptonsHandle_,"lep");

    // Clean jets
    std::vector<reco::GenJet> gen_jets_clean = CleanGenJets(gen_jets,gen_leptons);

    // Make pt, eta cuts on jets (after doing jet cleaning)
    gen_leptons = MakePtEtaCuts_forGenParticlesCollection(gen_leptons,"lep");

    // Get just charged leptons (recall std::vector<reco::GenParticle>> is an alias for std::vector<reco::GenParticle>>)
    reco::GenParticleCollection gen_leptons_charged = GetChargedGenParticle(gen_leptons);

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

    total_sm_xsec += sm_wgt;
    total_orig_xsec += originalXWGTUP_intree;

    h_eventsumEFT->Fill(0.5,1,eft_fit);
    h_SMwgt_norm->Fill(sm_wgt);

    // njets histograms
    h_nJetsEFT->Fill(gen_jets.size(),1.0,eft_fit);
    h_nJetsSM->Fill(gen_jets.size(),sm_wgt);

    // Particle level njets hists
    h_pl_nJets_EFT->Fill(pl_jets.size(),1.0,eft_fit);
    h_pl_nJets_SM->Fill(pl_jets.size(),sm_wgt);

    // njets histograms with some basic slection criteria
    if (gen_leptons_charged.size() >= 3) {
        h_nJets_3orMoreLep_EFT->Fill(gen_jets.size(),1.0,eft_fit);
        h_nJets_3orMoreLep_SM->Fill(gen_jets.size(),sm_wgt);
        h_nJets_cleanJets_3orMoreLep_EFT->Fill(gen_jets_clean.size(),1.0,eft_fit);
        h_nJets_cleanJets_3orMoreLep_SM->Fill(gen_jets_clean.size(),sm_wgt);
    }
    if (gen_leptons_charged.size() == 3) {
        h_nJets_3Lep_EFT->Fill(gen_jets.size(),1.0,eft_fit);
        h_nJets_3Lep_SM->Fill(gen_jets.size(),sm_wgt);
        h_nJets_cleanJets_3Lep_EFT->Fill(gen_jets_clean.size(),1.0,eft_fit);
        h_nJets_cleanJets_3Lep_SM->Fill(gen_jets_clean.size(),sm_wgt);
        // Particle level
        h_pl_nJets_3Lep_EFT->Fill(pl_jets.size(),1.0,eft_fit);
        h_pl_nJets_3Lep_SM->Fill(pl_jets.size(),sm_wgt);
    }


    eventnum_intree = event.id().event();
    lumiBlock_intree = event.id().luminosityBlock();
    runNumber_intree = event.id().run();

    summaryTree->Fill();

}

void EFTGenHistsWithCuts::beginRun(edm::Run const& run, edm::EventSetup const& evsetup)
{
    eventcount = 0;
}

void EFTGenHistsWithCuts::endRun(edm::Run const& run, edm::EventSetup const& evsetup)
{
    total_events += eventcount;
}

void EFTGenHistsWithCuts::beginLuminosityBlock(edm::LuminosityBlock const& lumi, edm::EventSetup const& evsetup)
{
    total_ls += 1;
}

void EFTGenHistsWithCuts::endLuminosityBlock(edm::LuminosityBlock const& lumi, edm::EventSetup const& evsetup){}


DEFINE_FWK_MODULE(EFTGenHistsWithCuts);


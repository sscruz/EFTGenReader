
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

    // Automatically declare histograms, store in hist_dict
    for (size_t i=0; i<lep_cats_vect.size(); i++){
        TString lep_cat = lep_cats_vect.at(i);
        for (size_t j=0; j<hist_info_vec.size(); j++){
            TString multiplicity_type = hist_info_vec.at(j).h_multiplicity; // How many histogram per event, e.g. one hist for all event (like nJets), or a hist for pairs of evens (like dR histograms)
            TString h_type = hist_info_vec.at(j).h_type;
            int     h_bins = hist_info_vec.at(j).h_bins;
            int     h_min  = hist_info_vec.at(j).h_min;
            int     h_max  = hist_info_vec.at(j).h_max;
            size_t  h_no   = abs(hist_info_vec.at(j).h_no);
           if (multiplicity_type == "all"){
                //TString hist_name = constructHistName(lep_cat,h_type,-1);
                TString hist_name = constructHistName(lep_cat,h_type,{});
                hist_dict[hist_name] = newfs->make<TH1EFT>(hist_name,hist_name,h_bins,h_min,h_max);
            }
            else{
                for (size_t k=0; k<h_no; k++){
                    if (multiplicity_type == "single"){
                        //TString hist_name = constructHistName(lep_cat,h_type,k+1);
                        TString hist_name = constructHistName(lep_cat,h_type,{k+1});
                        //std::cout << hist_name << std::endl;
                        hist_dict[hist_name] = newfs->make<TH1EFT>(hist_name,hist_name,h_bins,h_min,h_max);
                    }
                    if (multiplicity_type == "pair"){
                        for (size_t l=0; l<h_no; l++){
                            if (k<l){
                                TString hist_name = constructHistName(lep_cat,h_type,{k+1,l+1});
                                hist_dict[hist_name] = newfs->make<TH1EFT>(hist_name,hist_name,h_bins,h_min,h_max);
                                //std::cout << hist_name << std::endl;
                            }
                        }
                    }
                }
            }
        }
    }
    // Testing analysis catetory yield hists (for PL vs RECO)
    for (size_t i=0; i<ana_cats_vct.size(); i++){
        TString ana_cat = ana_cats_vct.at(i);
        TString hist_name = constructHistName(ana_cat,"yield",{});
        hist_dict[hist_name] = newfs->make<TH1EFT>(hist_name,hist_name,1,0,1);
        TString h_pl_njet_name = constructHistName(ana_cat,"yield-njets",{});
        hist_dict[h_pl_njet_name] = newfs->make<TH1EFT>(h_pl_njet_name,h_pl_njet_name,12,0,12);
        // TH1Ds for event sums
        TString h_ana_cat_pass_name = constructHistName(ana_cat,"n-events-pass",{});
        hist_TH1D_dict[h_ana_cat_pass_name] = newfs->make<TH1D>(h_ana_cat_pass_name,h_ana_cat_pass_name,1,0,1);
    }

    //// Declare histograms by hand ////

    int pdg_bins = 100;
    int njet_bins = 16;
    int pt_bins = 5;
    int eta_bins = 10;
    int invmass_bins = 30;
    int deltaR_bins = 10;
 
    //important variables for hists
    double lep_pt_min = 0;
    double lep_pt_max = 300;
    double jet_pt_min = 30;
    double jet_pt_max = 300;
    double jetpt_bin_width = (jet_pt_max-jet_pt_min)/pt_bins;


    //bin sizes for 2D jet vs. bjet hist
    int njet_bins_jetbjet = 8; //number of jet bins 
    int nbjet_bins_jetbjet = 5; //number of bjet bins

    // TEST PLvsGEN //
    h_PLvsGEN_gen_njetsclean    = newfs->make<TH1EFT>("h_PLvsGEN_gen_njetsclean","h_PLvsGEN_gen_njetsclean",12,0,12);
    h_PLvsGEN_gen_nbjets0p1        = newfs->make<TH1EFT>("h_PLvsGEN_gen_nbjets0p1","h_PLvsGEN_gen_nbjets0p1",10,0,10);
    h_PLvsGEN_gen_nbjets0p01        = newfs->make<TH1EFT>("h_PLvsGEN_gen_nbjets0p01","h_PLvsGEN_gen_nbjets0p01",10,0,10);
    h_PLvsGEN_gen_nchleps       = newfs->make<TH1EFT>("h_PLvsGEN_gen_nchleps","h_PLvsGEN_gen_nchleps",8,0,8);
    h_PLvsGEN_gen_nchlepsnotau  = newfs->make<TH1EFT>("h_PLvsGEN_gen_nchlepsnotau","h_PLvsGEN_gen_nchlepsnotau",8,0,8);
    h_PLvsGEN_gen_lepIDs        = newfs->make<TH1EFT>("h_PLvsGEN_gen_lepIDs","h_PLvsGEN_gen_lepIDs",8,10,18);
    h_PLvsGEN_gen_lepIDsnotau   = newfs->make<TH1EFT>("h_PLvsGEN_gen_lepIDsnotau","h_PLvsGEN_gen_lepIDsnotau",8,10,18);
    h_PLvsGEN_gen_sumallpt      = newfs->make<TH1EFT>("h_PLvsGEN_gen_sumallpt","h_PLvsGEN_gen_sumallpt",20,0,1500); 
    h_PLvsGEN_gen_sumallptnotau = newfs->make<TH1EFT>("h_PLvsGEN_gen_sumallptnotau","h_PLvsGEN_gen_sumallptnotau",20,0,1500); 
    h_PLvsGEN_gen_sumallleppt      = newfs->make<TH1EFT>("h_PLvsGEN_gen_sumallleppt","h_PLvsGEN_gen_sumallleppt",10,0,600);
    h_PLvsGEN_gen_sumalllepptnotau = newfs->make<TH1EFT>("h_PLvsGEN_gen_sumalllepptnotau","h_PLvsGEN_gen_sumalllepptnotau",10,0,600);
    h_PLvsGEN_gen_ht               = newfs->make<TH1EFT>("h_PLvsGEN_gen_ht","h_PLvsGEN_gen_ht",20,0,1500);

    h_PLvsGEN_pl_njets      = newfs->make<TH1EFT>("h_PLvsGEN_pl_njets","h_PLvsGEN_pl_njets",12,0,12);
    h_PLvsGEN_pl_njetsclean = newfs->make<TH1EFT>("h_PLvsGEN_pl_njetsclean","h_PLvsGEN_pl_njetsclean",12,0,12);
    h_PLvsGEN_pl_nbjets     = newfs->make<TH1EFT>("h_PLvsGEN_pl_nbjets","h_PLvsGEN_pl_nbjets",10,0,10);
    h_PLvsGEN_pl_nbjetsclean = newfs->make<TH1EFT>("h_PLvsGEN_pl_nbjetsclean","h_PLvsGEN_pl_nbjetsclean",10,0,10);
    h_PLvsGEN_pl_nleps      = newfs->make<TH1EFT>("h_PLvsGEN_pl_nleps","h_PLvsGEN_pl_nleps",8,0,8);
    h_PLvsGEN_pl_lepIDs     = newfs->make<TH1EFT>("h_PLvsGEN_pl_lepIDs","h_PLvsGEN_pl_lepIDs",8,10,18);
    h_PLvsGEN_pl_sumallpt   = newfs->make<TH1EFT>("h_PLvsGEN_pl_sumallpt","h_PLvsGEN_pl_sumallpt",20,0,1500); 
    h_PLvsGEN_pl_sumallleppt= newfs->make<TH1EFT>("h_PLvsGEN_pl_sumallleppt","h_PLvsGEN_pl_sumallleppt",10,0,600);
    h_PLvsGEN_pl_ht         = newfs->make<TH1EFT>("h_PLvsGEN_pl_ht","h_PLvsGEN_pl_ht",20,0,1500);
    ////

    //2D jets vs. bjets hists for various leptons categories
    h_2lss_jetbjetEFT = newfs->make<TH2EFT>("h_2lss_jetbjetEFT","h_2lss_jetbjetEFT;N_{jets};N_{bjets}",njet_bins_jetbjet,0,njet_bins_jetbjet,nbjet_bins_jetbjet,0,nbjet_bins_jetbjet);
    h_2lss_jetbjetSM =	newfs->make<TH2D>("h_2lss_jetbjetSM","h_2lss_jetbjetSM;N_{jets};N_{bjets}",njet_bins_jetbjet,0,njet_bins_jetbjet,nbjet_bins_jetbjet,0,nbjet_bins_jetbjet);
    h_3l_jetbjetEFT =	newfs->make<TH2EFT>("h_3l_jetbjetEFT","h_3l_jetbjetEFT;N_{jets};N_{bjets}",njet_bins_jetbjet,0,njet_bins_jetbjet,nbjet_bins_jetbjet,0,nbjet_bins_jetbjet);
    h_3l_jetbjetSM =	newfs->make<TH2D>("h_3l_jetbjetSM","h_3l_jetbjetSM;N_{jets};N_{bjets}",njet_bins_jetbjet,0,njet_bins_jetbjet,nbjet_bins_jetbjet,0,nbjet_bins_jetbjet);
    h_3l_sfz_jetbjetEFT =	newfs->make<TH2EFT>("h_3l_sfz_jetbjetEFT","h_3l_sfz_jetbjetEFT;N_{jets};N_{bjets}",njet_bins_jetbjet,0,njet_bins_jetbjet,nbjet_bins_jetbjet,0,nbjet_bins_jetbjet);
    h_3l_sfz_jetbjetSM =	newfs->make<TH2D>("h_3l_sfz_jetbjetSM","h_3l_sfz_jetbjetSM;N_{jets};N_{bjets}",njet_bins_jetbjet,0,njet_bins_jetbjet,nbjet_bins_jetbjet,0,nbjet_bins_jetbjet);
    h_4l_jetbjetEFT =	newfs->make<TH2EFT>("h_4l_jetbjetEFT","h_4l_jetbjetEFT;N_{jets};N_{bjets}",njet_bins_jetbjet,0,njet_bins_jetbjet,nbjet_bins_jetbjet,0,nbjet_bins_jetbjet);
    h_4l_jetbjetSM =	newfs->make<TH2D>("h_4l_jetbjetSM","h_4l_jetbjetSM;N_{jets};N_{bjets}",njet_bins_jetbjet,0,njet_bins_jetbjet,nbjet_bins_jetbjet,0,nbjet_bins_jetbjet);

    // Book the histograms that we will fill in the event loop
    h_eventsumEFT = newfs->make<TH1EFT>("h_eventsumEFT","h_eventsumEFT",1,0,1);

    // Jet histograms
    //h_nJetsEFT= newfs->make<TH1EFT>("h_njetsEFT","h_njetsEFT",njet_bins,0,njet_bins);
    //h_nJetsSM = newfs->make<TH1D>("h_njetsSM","h_njetsSM",njet_bins,0,njet_bins);

    // These are njets histograms with some basic selections
    //h_nJets_3orMoreLep_EFT = newfs->make<TH1EFT>("h_nJets_3orMoreLep_EFT","h_nJets_3orMoreLep_EFT",njet_bins,0,njet_bins);
    //h_nJets_3orMoreLep_SM  = newfs->make<TH1D>("h_nJets_3orMoreLep_SM","h_nJets_3orMoreLep_SM",njet_bins,0,njet_bins);
    //h_nJets_3Lep_EFT       = newfs->make<TH1EFT>("h_nJets_3Lep_EFT","h_nJets_3Lep_EFT",njet_bins,0,njet_bins);
    //h_nJets_3Lep_SM        = newfs->make<TH1D>("h_nJets_3Lep_SM","h_nJets_3Lep_SM",njet_bins,0,njet_bins);
    //h_nJets_cleanJets_3orMoreLep_EFT = newfs->make<TH1EFT>("h_nJets_cleanJets_3orMoreLep_EFT","h_nJets_cleanJets_3orMoreLep_EFT",njet_bins,0,njet_bins);
    //h_nJets_cleanJets_3orMoreLep_SM  = newfs->make<TH1D>("h_nJets_cleanJets_3orMoreLep_SM","h_nJets_cleanJets_3orMoreLep_SM",njet_bins,0,njet_bins);
    //h_nJets_cleanJets_3Lep_EFT       = newfs->make<TH1EFT>("h_nJets_cleanJets_3Lep_EFT","h_nJets_cleanJets_3Lep_EFT",njet_bins,0,njet_bins);
    //h_nJets_cleanJets_3Lep_SM        = newfs->make<TH1D>("h_nJets_cleanJets_3Lep_SM","h_nJets_cleanJets_3Lep_SM",njet_bins,0,njet_bins);

    // Particle level hists
    //h_pl_nJets_EFT = newfs->make<TH1EFT>("h_pl_njetsEFT","h_pl_njetsEFT",njet_bins,0,njet_bins);
    //h_pl_nJets_SM  = newfs->make<TH1D>("h_pl_njetsSM","h_pl_njetsSM",njet_bins,0,njet_bins);
    h_pl_nJets_3Lep_EFT = newfs->make<TH1EFT>("h_pl_nJets_3Lep_EFT","h_pl_nJets_3Lep_EFT",njet_bins,0,njet_bins);
    h_pl_nJets_3Lep_SM  = newfs->make<TH1D>("h_pl_nJets_3Lep_SM","h_pl_nJets_3Lep_SM",njet_bins,0,njet_bins);
    h_pl_clean_nJets_3Lep_EFT = newfs->make<TH1EFT>("h_pl_clean_nJets_3Lep_EFT","h_pl_clean_nJets_3Lep_EFT",njet_bins,0,njet_bins);

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
    //std::vector<reco::GenJet> gen_bjets = GetGenBJets(*genJets); // Does not work

    // Particle level stuff:
    edm::Handle<std::vector<reco::GenJet>> particleLevelJetsHandle_;
    edm::Handle<std::vector<reco::GenJet>> particleLevelLeptonsHandle_;
    event.getByToken(particleLevelJetsToken_,particleLevelJetsHandle_);
    event.getByToken(particleLevelLeptonsToken_,particleLevelLeptonsHandle_);
    std::vector<reco::GenJet> pl_jets    = MakePtEtaCuts(*particleLevelJetsHandle_,min_pt_jet,max_eta_jet);
    std::vector<reco::GenJet> pl_leptons = MakePtEtaCuts(*particleLevelLeptonsHandle_,min_pt_lep,max_eta_lep);
    std::vector<reco::GenJet> pl_bjets = GetGenBJets(pl_jets);

    // Do not use!!! Should clean the PL leptons, not the PL jets
    std::vector<reco::GenJet> pl_jets_clean = CleanCollection(pl_jets,*particleLevelLeptonsHandle_,0.4); // Clean gen jets
    std::vector<reco::GenJet> pl_bjets_clean = GetGenBJets(pl_jets_clean);

    // Clean PL leptons
    pl_leptons = CleanCollection(pl_leptons,pl_jets,0.4);

    // Clean jets
    std::vector<reco::GenJet> gen_jets_clean = CleanCollection(gen_jets,gen_leptons,0.4);
    // Get gen b jets from cleaned gen jets
    //std::vector<reco::GenJet> gen_bjets = GetGenBJets(gen_jets_clean); // Does not work for gen, only pl
    // Find the b particles
    gen_b = MakePtEtaCuts(gen_b,min_pt_jet,max_eta_jet);
    std::vector<reco::GenJet> gen_bjets_fromDRtest0p1 = GetGenJetsFromDR(gen_jets_clean,gen_b,0.1);
    std::vector<reco::GenJet> gen_bjets_fromDRtest0p01 = GetGenJetsFromDR(gen_jets_clean,gen_b,0.01);

    // Make pt, eta cuts on leptons (after doing jet cleaning)
    gen_leptons = MakePtEtaCuts(gen_leptons,min_pt_lep,max_eta_lep);

    // Get just charged leptons (recall std::vector<reco::GenParticle>> is an alias for std::vector<reco::GenParticle>>)
    reco::GenParticleCollection gen_leptons_charged = getChargedParticles(gen_leptons);


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

    // Find what lepton category (if any) this even falls into
    TString lep_cat_name = getLepCat(gen_leptons);
    std::vector<TString> cats_vect;
    if (lep_cat_name != "none"){
        cats_vect = {lep_cat_name,"anyLepCat"};
    } else {
        cats_vect = {lep_cat_name};
    }

    for (auto lep_cat: cats_vect){
        // Loop over jets and fill jet hists automatically
        double ht=0;
        for (size_t i = 0; i < gen_jets_clean.size(); i++) {
            const reco::GenJet& p1= gen_jets.at(i);
            double pt = p1.p4().Pt();
            double eta = p1.p4().Eta();
            //std::cout << pt << std::endl;
            TString h_pt_name = constructHistName(lep_cat,"jet_pt",{i+1});
            TString h_eta_name = constructHistName(lep_cat,"jet_eta",{i+1});
            fillHistIfExists(h_pt_name,pt,eft_fit);
            fillHistIfExists(h_eta_name,eta,eft_fit);
            ht = ht + pt;
            for (size_t j = 0; j < gen_jets_clean.size(); j++) {
                const reco::GenJet& p2 = gen_jets.at(j);
                double dR = getdR(p1,p2);
                int hist_number = 10*(i+1)+(j+1);
                TString h_dR_name = constructHistName(lep_cat,"jet_dR",{i+1,j+1});
                fillHistIfExists(h_dR_name,dR,eft_fit);
            }
        }

        // Fill jet hists that include info for all jets in event
        TString h_ht_name = constructHistName(lep_cat,"ht",{});
        fillHistIfExists(h_ht_name,ht,eft_fit);
        TString h_njet_name = constructHistName(lep_cat,"njets",{});
        fillHistIfExists(h_njet_name,gen_jets_clean.size(),eft_fit);

        // Loop over leptonss and fill hists automatically
        for (size_t i = 0; i < gen_leptons_charged.size(); i++) {
            const reco::GenParticle& p1 = gen_leptons_charged.at(i);
            double pt = p1.p4().Pt();
            double eta = p1.p4().Eta();
            TString h_pt_name = constructHistName(lep_cat,"lep_pt",{i+1});
            TString h_eta_name = constructHistName(lep_cat,"lep_eta",{i+1});
            fillHistIfExists(h_pt_name,pt,eft_fit);
            fillHistIfExists(h_eta_name,eta,eft_fit);
            for (size_t j = 0; j < gen_leptons_charged.size(); j++) {
                const reco::GenParticle& p2 = gen_leptons_charged.at(j);
                double dR = getdR(p1,p2);
                double mll = getInvMass(p1,p2);
                int hist_number = 10*(i+1)+(j+1);
                TString h_dR_name = constructHistName(lep_cat,"lep_dR",{i+1,j+1});
                TString h_mll_name = constructHistName(lep_cat,"lep_mll",{i+1,j+1});
                fillHistIfExists(h_dR_name,dR,eft_fit);
                fillHistIfExists(h_mll_name,mll,eft_fit);
            }
        }
    }
    

    //// Filling histograms by hand ////

    /* // Probably we can delete these now
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
    */
    if (gen_leptons_charged.size() == 3) {
        // We not automatically fill the equivalent of thse, probably can delete
        //h_nJets_3Lep_EFT->Fill(gen_jets.size(),1.0,eft_fit);
        //h_nJets_3Lep_SM->Fill(gen_jets.size(),sm_wgt);
        //h_nJets_cleanJets_3Lep_EFT->Fill(gen_jets_clean.size(),1.0,eft_fit);
        //h_nJets_cleanJets_3Lep_SM->Fill(gen_jets_clean.size(),sm_wgt);
        // Particle level
        h_pl_nJets_3Lep_EFT->Fill(pl_jets.size(),1.0,eft_fit);
        h_pl_nJets_3Lep_SM->Fill(pl_jets.size(),sm_wgt);
    }
    if (gen_leptons_charged.size() == 3) {
        h_pl_clean_nJets_3Lep_EFT->Fill(pl_jets_clean.size(),1.0,eft_fit);
    }

    //////////////////////////////////
    ////////// TEST PLvsGEN //////////
    //////////////////////////////////

    /*
    std::cout << "\tgen_leptons: " << gen_leptons.size() << std::endl;
    std::cout << "\tgen_leptons_charged: " << gen_leptons_charged.size() << std::endl;
    for (auto l: gen_leptons_charged){
        std::cout << "\tcharge, id : " << l.charge() << " " << l.pdgId() << std::endl;
    }
    std::cout << "\tgen_jets: " << gen_jets.size() << std::endl;
    std::cout << "\tgen_jets_clean: " << gen_jets_clean.size() << std::endl;
    std::cout << "\tgen_b: " << gen_b.size() << std::endl;
    //std::cout << "\tgen_bjets: " << gen_bjets.size() << std::endl;

    std::cout << " " << std::endl;
    std::cout << "\tpl_leptons: " << pl_leptons.size() << std::endl;
    for (auto l: pl_leptons){
        std::cout << "\tcharge, id : " << l.charge() << " " << l.pdgId() << std::endl;
    }
    std::cout << "\tpl_jets: " << pl_jets.size() << std::endl;
    std::cout << "\tpl_jets_clean: " << pl_jets_clean.size() << std::endl;
    std::cout << "\tpl_bjets: " << pl_bjets.size() << std::endl;
    */

    //std::cout << "\t" << .size() << std::endl;
    //Fill->(,1.0,eft_fit);

    double gen_leppt = 0;
    double gen_leppt_notau = 0;
    double gen_ht = 0;
    double pl_leppt = 0;
    double pl_ht = 0;
    reco::GenParticleCollection gen_leptons_charged_notau;

    // GEN
    h_PLvsGEN_gen_njetsclean->Fill(gen_jets_clean.size(),1.0,eft_fit);
    h_PLvsGEN_gen_nchleps->Fill(gen_leptons_charged.size(),1.0,eft_fit);
    h_PLvsGEN_gen_nbjets0p1->Fill(gen_bjets_fromDRtest0p1.size(),1.0,eft_fit);
    h_PLvsGEN_gen_nbjets0p01->Fill(gen_bjets_fromDRtest0p01.size(),1.0,eft_fit);
    for (auto l: gen_leptons_charged){
        h_PLvsGEN_gen_lepIDs->Fill(abs(l.pdgId()),1.0,eft_fit);
        gen_leppt = gen_leppt + l.p4().Pt();
        if (abs(l.pdgId()) != 15){
            gen_leptons_charged_notau.push_back(l);
            h_PLvsGEN_gen_lepIDsnotau->Fill(abs(l.pdgId()),1.0,eft_fit);
            gen_leppt_notau = gen_leppt_notau + l.p4().Pt();
        }
    }
    h_PLvsGEN_gen_nchlepsnotau->Fill(gen_leptons_charged_notau.size(),1.0,eft_fit);
    for (auto j: gen_jets_clean){
        gen_ht = gen_ht + j.p4().Pt();
    }
    h_PLvsGEN_gen_sumallleppt->Fill(gen_leppt,1.0,eft_fit);
    h_PLvsGEN_gen_sumalllepptnotau->Fill(gen_leppt_notau,1.0,eft_fit);
    h_PLvsGEN_gen_ht->Fill(gen_ht,1.0,eft_fit);
    h_PLvsGEN_gen_sumallpt->Fill(gen_ht+gen_leppt,1.0,eft_fit);
    h_PLvsGEN_gen_sumallptnotau->Fill(gen_ht+gen_leppt_notau,1.0,eft_fit);

    // PL
    h_PLvsGEN_pl_njets->Fill(pl_jets.size(),1.0,eft_fit);
    h_PLvsGEN_pl_njetsclean->Fill(pl_jets_clean.size(),1.0,eft_fit);
    h_PLvsGEN_pl_nbjets->Fill(pl_bjets.size(),1.0,eft_fit);
    h_PLvsGEN_pl_nbjetsclean->Fill(pl_bjets_clean.size(),1.0,eft_fit);
    h_PLvsGEN_pl_nleps->Fill(pl_leptons.size(),1.0,eft_fit);
    for (auto l: pl_leptons){
        h_PLvsGEN_pl_lepIDs->Fill(abs(l.pdgId()),1.0,eft_fit);
        pl_leppt = pl_leppt + l.p4().Pt();
    }
    for (auto j: pl_jets_clean){
        pl_ht = pl_ht + j.p4().Pt();
    }
    h_PLvsGEN_pl_sumallleppt->Fill(pl_leppt,1.0,eft_fit);
    h_PLvsGEN_pl_ht->Fill(pl_ht,1.0,eft_fit);
    h_PLvsGEN_pl_sumallpt->Fill(pl_leppt+pl_ht,1.0,eft_fit);

    //////////////////////////////////////////

    // Testing analysis catetory yield hists (for PL vs RECO)
    TString ana_cat = getAnaCat(pl_leptons,pl_jets,pl_bjets);
    // Yield hist
    TString ana_cat_hist_name = constructHistName(ana_cat,"yield",{});
    fillHistIfExists(ana_cat_hist_name,0.5,eft_fit);
    // Yield njets hists
    TString h_pl_njet_name = constructHistName(ana_cat,"yield-njets",{});
    fillHistIfExists(h_pl_njet_name,GetNJetsForLepCat(pl_leptons,pl_jets),eft_fit);
    // N events passing hist
    TString h_ana_cat_pass_name = constructHistName(ana_cat,"n-events-pass",{});
    fillTH1DHistIfExists(h_ana_cat_pass_name,0.5);

    //////////////////////////////////////////
    //////////////////////////////////////////


    //Filling the 2D hists for jets vs. bjets (different lepton categories). Overblow bins taken care of. 

    double njet_max = 8.0;
    double nbjet_max = 5.0;

    //lepton categories
    if(pl_leptons.size() == 2)
    {
	const reco::GenParticle& p1 = pl_leptons.at(0);
        const reco::GenParticle& p2 = pl_leptons.at(1);
	double pt1 = p1.p4().Pt();
        double pt2 = p2.p4().Pt();
        int id1 = p1.pdgId();
        int id2 = p2.pdgId();
        int id_product = id1*id2;
        if(id_product > 0 && pt1 > 25 && pt2 > 15)
        {
	    if(pl_jets.size()<njet_max && pl_bjets.size()<nbjet_max)
            {
                h_2lss_jetbjetEFT->Fill(pl_jets.size(),pl_bjets.size(),1.0,eft_fit);
                h_2lss_jetbjetSM->Fill(pl_jets.size(),pl_bjets.size(),sm_wgt);
            }

	    else if(pl_jets.size()<njet_max && pl_bjets.size()>nbjet_max)
            {
                h_2lss_jetbjetEFT->Fill(pl_jets.size(),nbjet_max-1,1.0,eft_fit);
                h_2lss_jetbjetSM->Fill(pl_jets.size(),nbjet_max-1,sm_wgt);
            }

	    else if(pl_jets.size()>njet_max && pl_bjets.size()>nbjet_max)
            {
                h_2lss_jetbjetEFT->Fill(njet_max-1,pl_bjets.size(),1.0,eft_fit);
                h_2lss_jetbjetSM->Fill(njet_max-1,pl_bjets.size(),sm_wgt);
            }

	    else
	    {
                h_2lss_jetbjetEFT->Fill(njet_max-1,nbjet_max-1,1.0,eft_fit);
                h_2lss_jetbjetSM->Fill(njet_max-1,nbjet_max-1,sm_wgt);
            }
 	  

        }
    }

    if(pl_leptons.size() == 3)
    {
	const reco::GenParticle& p1 = pl_leptons.at(0);
        const reco::GenParticle& p2 = pl_leptons.at(1);
        const reco::GenParticle& p3 = pl_leptons.at(2);
	double pt1 = p1.p4().Pt();
        double pt2 = p2.p4().Pt();
        double pt3 = p3.p4().Pt();
	if(pt1 > 25 && pt2 > 15 && pt3 > 10)
	{
            bool isSFOSZ = false;
            //look for OS
            if(abs(p1.charge() + p2.charge() + p3.charge()) == 1) {
                //loop over all leptons
                for(size_t i = 0; i < pl_leptons.size(); i++) {
                    //loop over leptons not seen yet (prevents (1,2) (2,1) duplicates)
                    for(size_t j = i+1; j < pl_leptons.size(); j++) {
                        //look for same flavor only (Z->ee or Z->mumu)
                        if(abs(pl_leptons.at(i).pdgId()) != abs(pl_leptons.at(j).pdgId())) continue;
                        //look for opposite sign only (e+e- or mu+mu-)
                        if(pl_leptons.at(i).charge() * pl_leptons.at(j).charge() > 0) continue;
                        //look for |M(ll) - M(Z)| < 10 GeV
                        if(fabs((pl_leptons.at(i).p4() + pl_leptons.at(j).p4()).M() - 91.2) < 10) {
                            isSFOSZ = true;
                            break; //done searching
                        }
                    }
                    if(isSFOSZ) break; //no need to search any further
                }
            }
            if(pl_jets.size()<njet_max && pl_bjets.size()<nbjet_max)
            {
                if(isSFOSZ) {
                h_3l_sfz_jetbjetEFT->Fill(pl_jets.size(),pl_bjets.size(),1.0,eft_fit);
                h_3l_sfz_jetbjetSM->Fill(pl_jets.size(),pl_bjets.size(),sm_wgt);
                }
                else {
                h_3l_jetbjetEFT->Fill(pl_jets.size(),pl_bjets.size(),1.0,eft_fit);
                h_3l_jetbjetSM->Fill(pl_jets.size(),pl_bjets.size(),sm_wgt);
                }
            }

	    else if(pl_jets.size()<njet_max && pl_bjets.size()>nbjet_max)
            {
                if(isSFOSZ) {
                h_3l_sfz_jetbjetEFT->Fill(pl_jets.size(),nbjet_max-1,1.0,eft_fit);
                h_3l_sfz_jetbjetSM->Fill(pl_jets.size(),nbjet_max-1,sm_wgt);
                }
                else {
                h_3l_jetbjetEFT->Fill(pl_jets.size(),nbjet_max-1,1.0,eft_fit);
                h_3l_jetbjetSM->Fill(pl_jets.size(),nbjet_max-1,sm_wgt);
                }
            }

	    else if(pl_jets.size()>njet_max && pl_bjets.size()>nbjet_max)
            {
                if(isSFOSZ) {
                h_3l_sfz_jetbjetEFT->Fill(njet_max-1,pl_bjets.size(),1.0,eft_fit);
                h_3l_sfz_jetbjetSM->Fill(njet_max-1,pl_bjets.size(),sm_wgt);
                }
                else {
                h_3l_jetbjetEFT->Fill(njet_max-1,pl_bjets.size(),1.0,eft_fit);
                h_3l_jetbjetSM->Fill(njet_max-1,pl_bjets.size(),sm_wgt);
                }
            }

	    else
	    {
                if(isSFOSZ) {
                h_3l_sfz_jetbjetEFT->Fill(njet_max-1,nbjet_max-1,1.0,eft_fit);
                h_3l_sfz_jetbjetSM->Fill(njet_max-1,nbjet_max-1,sm_wgt);
                }
                else {
                h_3l_jetbjetEFT->Fill(njet_max-1,nbjet_max-1,1.0,eft_fit);
                h_3l_jetbjetSM->Fill(njet_max-1,nbjet_max-1,sm_wgt);
                }
            }
	}
    }

    if(pl_leptons.size() >= 4)
    {
	const reco::GenParticle& p1 = pl_leptons.at(0);
        const reco::GenParticle& p2 = pl_leptons.at(1);
        const reco::GenParticle& p3 = pl_leptons.at(2);
        const reco::GenParticle& p4 = pl_leptons.at(3);
    	double pt1 = p1.p4().Pt();
    	double pt2 = p2.p4().Pt();
    	double pt3 = p3.p4().Pt();
    	double pt4 = p4.p4().Pt();

	if(pt1 > 25 && pt2 > 15 && pt3 > 10 && pt4 > 10)
	{
            if(pl_jets.size()<njet_max && pl_bjets.size()<nbjet_max)
            {
                h_4l_jetbjetEFT->Fill(pl_jets.size(),pl_bjets.size(),1.0,eft_fit);
                h_4l_jetbjetSM->Fill(pl_jets.size(),pl_bjets.size(),sm_wgt);
            }

	    else if(pl_jets.size()<njet_max && pl_bjets.size()>nbjet_max)
            {
                h_4l_jetbjetEFT->Fill(pl_jets.size(),nbjet_max-1,1.0,eft_fit);
                h_4l_jetbjetSM->Fill(pl_jets.size(),nbjet_max-1,sm_wgt);
            }

	    else if(pl_jets.size()>njet_max && pl_bjets.size()>nbjet_max)
            {
                h_4l_jetbjetEFT->Fill(njet_max-1,pl_bjets.size(),1.0,eft_fit);
                h_4l_jetbjetSM->Fill(njet_max-1,pl_bjets.size(),sm_wgt);
            }

	    else
	    {
                h_4l_jetbjetEFT->Fill(njet_max-1,nbjet_max-1,1.0,eft_fit);
                h_4l_jetbjetSM->Fill(njet_max-1,nbjet_max-1,sm_wgt);
            }
	}
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



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
    staggered_pt_cuts_lep = iConfig.getParameter<std::vector<double>> ("staggered_pt_cuts_lep");
    min_njets_2lss  = iConfig.getParameter<int> ("min_njets_2lss");
    min_njets_3l    = iConfig.getParameter<int> ("min_njets_3l");
    min_njets_4l    = iConfig.getParameter<int> ("min_njets_4l");
    max_njet_bins_2lss = iConfig.getParameter<int> ("max_njet_bins_2lss");
    max_njet_bins_3l   = iConfig.getParameter<int> ("max_njet_bins_3l");
    max_njet_bins_4l   = iConfig.getParameter<int> ("max_njet_bins_4l");

    
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
    for (size_t i=0; i<ana_cats_vct.size(); i++){
        TString ana_cat = ana_cats_vct.at(i);
        for (size_t j=0; j<hist_info_vec.size(); j++){
            TString multiplicity_type = hist_info_vec.at(j).h_multiplicity; // How many histogram per event, e.g. one hist for all event (like nJets), or a hist for pairs of evens (like dR histograms)
            TString h_type = hist_info_vec.at(j).h_type;
            int     h_bins = hist_info_vec.at(j).h_bins;
            int     h_min  = hist_info_vec.at(j).h_min;
            int     h_max  = hist_info_vec.at(j).h_max;
            size_t  h_no   = abs(hist_info_vec.at(j).h_no);
           if (multiplicity_type == "all"){
                TString hist_name = ConstructHistName(ana_cat,h_type,{});
                hist_dict[hist_name] = newfs->make<TH1EFT>(hist_name,hist_name,h_bins,h_min,h_max);
            }
            else{
                for (size_t k=0; k<h_no; k++){
                    if (multiplicity_type == "single"){
                        TString hist_name = ConstructHistName(ana_cat,h_type,{k+1});
                        //std::cout << hist_name << std::endl;
                        hist_dict[hist_name] = newfs->make<TH1EFT>(hist_name,hist_name,h_bins,h_min,h_max);
                    }
                    if (multiplicity_type == "pair"){
                        for (size_t l=0; l<h_no; l++){
                            if (k<l){
                                TString hist_name = ConstructHistName(ana_cat,h_type,{k+1,l+1});
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
        TString hist_name = ConstructHistName(ana_cat,"yield",{});
        hist_dict[hist_name] = newfs->make<TH1EFT>(hist_name,hist_name,1,0,1);
        TString h_pl_njet_name = ConstructHistName(ana_cat,"yield-njets",{});
        hist_dict[h_pl_njet_name] = newfs->make<TH1EFT>(h_pl_njet_name,h_pl_njet_name,12,0,12);
        // TH1Ds for event sums
        TString h_ana_cat_pass_name = ConstructHistName(ana_cat,"n-events-pass",{});
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

    //2D jets vs. bjets hists for various leptons categories
    h_2lss_jetbjetEFT = newfs->make<TH2EFT>("h_2lss_jetbjetEFT","h_2lss_jetbjetEFT;N_{jets};N_{bjets}",njet_bins_jetbjet,0,njet_bins_jetbjet,nbjet_bins_jetbjet,0,nbjet_bins_jetbjet);
    h_2lss_jetbjetSM =  newfs->make<TH2D>("h_2lss_jetbjetSM","h_2lss_jetbjetSM;N_{jets};N_{bjets}",njet_bins_jetbjet,0,njet_bins_jetbjet,nbjet_bins_jetbjet,0,nbjet_bins_jetbjet);
    h_3l_jetbjetEFT =   newfs->make<TH2EFT>("h_3l_jetbjetEFT","h_3l_jetbjetEFT;N_{jets};N_{bjets}",njet_bins_jetbjet,0,njet_bins_jetbjet,nbjet_bins_jetbjet,0,nbjet_bins_jetbjet);
    h_3l_jetbjetSM =    newfs->make<TH2D>("h_3l_jetbjetSM","h_3l_jetbjetSM;N_{jets};N_{bjets}",njet_bins_jetbjet,0,njet_bins_jetbjet,nbjet_bins_jetbjet,0,nbjet_bins_jetbjet);
    h_3l_sfz_jetbjetEFT =       newfs->make<TH2EFT>("h_3l_sfz_jetbjetEFT","h_3l_sfz_jetbjetEFT;N_{jets};N_{bjets}",njet_bins_jetbjet,0,njet_bins_jetbjet,nbjet_bins_jetbjet,0,nbjet_bins_jetbjet);
    h_3l_sfz_jetbjetSM =        newfs->make<TH2D>("h_3l_sfz_jetbjetSM","h_3l_sfz_jetbjetSM;N_{jets};N_{bjets}",njet_bins_jetbjet,0,njet_bins_jetbjet,nbjet_bins_jetbjet,0,nbjet_bins_jetbjet);
    h_4l_jetbjetEFT =   newfs->make<TH2EFT>("h_4l_jetbjetEFT","h_4l_jetbjetEFT;N_{jets};N_{bjets}",njet_bins_jetbjet,0,njet_bins_jetbjet,nbjet_bins_jetbjet,0,nbjet_bins_jetbjet);
    h_4l_jetbjetSM =    newfs->make<TH2D>("h_4l_jetbjetSM","h_4l_jetbjetSM;N_{jets};N_{bjets}",njet_bins_jetbjet,0,njet_bins_jetbjet,nbjet_bins_jetbjet,0,nbjet_bins_jetbjet);

    //1D Z boson hists                                                                                                                                                                                     
    //h_3l_sfz_Zpt = newfs->make<TH1EFT>("h_3l_sfz_Zpt", "h_3l_sfz_Zpt", 10, lep_pt_min, 400);
    //h_3l_sfz_cos = newfs->make<TH1EFT>("h_3l_sfz_cos", "h_3l_sfz_cos", 5, -1, 1);

    // Book the histograms that we will fill in the event loop
    h_eventsumEFT = newfs->make<TH1EFT>("h_eventsumEFT","h_eventsumEFT",1,0,1);

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
    event.getByToken(lheInfo_token_,LHEInfo);

    /*
    // Gen (no longer used)
    edm::Handle<reco::GenParticleCollection> prunedParticles;
    edm::Handle<std::vector<reco::GenJet> > genJets;
    event.getByToken(genParticles_token_,prunedParticles);
    event.getByToken(genJets_token_,genJets);
    reco::GenParticleCollection gen_leptons = GetGenLeptons(*prunedParticles);
    //reco::GenParticleCollection gen_b = GetGenParticlesSubset(*prunedParticles, 5);
    std::vector<reco::GenJet> gen_jets = GetGenJets(*genJets);
    // Clean jets
    std::vector<reco::GenJet> gen_jets_clean = CleanCollection(gen_jets,gen_leptons,0.4);
    // Get gen b jets (has not really been tested)
    //gen_b = MakeBaselinePtEtaCuts(gen_b,min_pt_jet,max_eta_jet);
    //std::vector<reco::GenJet> gen_bjets_fromDRtest0p1 = GetGenJetsFromDR(gen_jets_clean,gen_b,0.1);
    // Make pt, eta cuts on leptons (after doing jet cleaning)
    gen_leptons = MakeBaselinePtEtaCuts(gen_leptons,min_pt_lep,max_eta_lep);
    // Get just charged gen leptons (recall std::vector<reco::GenParticle>> is an alias for std::vector<reco::GenParticle>>)
    reco::GenParticleCollection gen_leptons_charged = GetChargedParticles(gen_leptons);
    */

    // Particle level //
    edm::Handle<std::vector<reco::GenJet>> particleLevelJetsHandle_;
    edm::Handle<std::vector<reco::GenJet>> particleLevelLeptonsHandle_;
    event.getByToken(particleLevelJetsToken_,particleLevelJetsHandle_);
    event.getByToken(particleLevelLeptonsToken_,particleLevelLeptonsHandle_);
    std::vector<reco::GenJet> pl_jets    = MakeBaselinePtEtaCuts(*particleLevelJetsHandle_,min_pt_jet,max_eta_jet);
    std::vector<reco::GenJet> pl_bjets   = GetGenBJets(pl_jets);
    std::vector<reco::GenJet> pl_leptons = MakeBaselinePtEtaCuts(*particleLevelLeptonsHandle_,min_pt_lep,max_eta_lep);
    pl_leptons = MakeStaggeredPtCuts(pl_leptons,staggered_pt_cuts_lep,min_pt_lep);
    pl_leptons = CleanCollection(pl_leptons,pl_jets,0.4); // Clean PL leptons (Should clean the PL leptons, not the PL jets: https://twiki.cern.ch/twiki/bin/view/LHCPhysics/ParticleLevelTopDefinitions)


    // Get eft_fit
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

    // Keep track of total xsec
    total_sm_xsec += sm_wgt;
    total_orig_xsec += originalXWGTUP_intree;

    // Fill h_eventsumEFT and h_SMwgt_norm hists
    h_eventsumEFT->Fill(0.5,1,eft_fit);
    h_SMwgt_norm->Fill(sm_wgt);

    // Find what analysis categories (if any) this even falls into
    TString ana_cat_name = GetAnaCat(pl_leptons,pl_jets,pl_bjets);
    //TString lep_cat_name = GetLepCat(pl_leptons);
    //std::cout << "Lep cat: " << lep_cat_name << ", Ana cat: " << ana_cat_name << std::endl;
    std::vector<TString> cats_vect;
    if (ana_cat_name != "none"){
        cats_vect = {ana_cat_name,"anyAnaCat"};
    } else {
        cats_vect = {ana_cat_name};
    }

    // Loop over leptop categories
    for (auto ana_cat: cats_vect){

        // Loop over jets and fill jet hists automatically
        double ht=0;
        for (size_t i = 0; i < pl_jets.size(); i++) {
            const reco::GenJet& p1 = pl_jets.at(i);
            double pt = p1.p4().Pt();
            double eta = p1.p4().Eta();
            //std::cout << pt << std::endl;
            TString h_pt_name = ConstructHistName(ana_cat,"jet_pt",{i+1});
            TString h_eta_name = ConstructHistName(ana_cat,"jet_eta",{i+1});
            FillHistIfExists(h_pt_name,pt,eft_fit);
            FillHistIfExists(h_eta_name,eta,eft_fit);
            ht = ht + pt;
            for (size_t j = 0; j < pl_jets.size(); j++) {
                const reco::GenJet& p2 = pl_jets.at(j);
                double dR = getdR(p1,p2);
                int hist_number = 10*(i+1)+(j+1);
                TString h_dR_name = ConstructHistName(ana_cat,"jet_dR",{i+1,j+1});
                FillHistIfExists(h_dR_name,dR,eft_fit);
            }
        }

        // Fill jet hists that include info for all jets in event
        TString h_ht_name = ConstructHistName(ana_cat,"ht",{});
        FillHistIfExists(h_ht_name,ht,eft_fit);
        TString h_njet_name = ConstructHistName(ana_cat,"njets",{});
        FillHistIfExists(h_njet_name,pl_jets.size(),eft_fit);

        // Loop over leptonss and fill hists automatically
        for (size_t i = 0; i < pl_leptons.size(); i++) {
            const reco::GenJet& p1 = pl_leptons.at(i);
            double pt = p1.p4().Pt();
            double eta = p1.p4().Eta();
            TString h_pt_name = ConstructHistName(ana_cat,"lep_pt",{i+1});
            TString h_eta_name = ConstructHistName(ana_cat,"lep_eta",{i+1});
            FillHistIfExists(h_pt_name,pt,eft_fit);
            FillHistIfExists(h_eta_name,eta,eft_fit);
            for (size_t j = 0; j < pl_leptons.size(); j++) {
                const reco::GenJet& p2 = pl_leptons.at(j);
                double dR = getdR(p1,p2);
                double mll = GetInvMass(p1,p2);
                int hist_number = 10*(i+1)+(j+1);
                TString h_dR_name = ConstructHistName(ana_cat,"lep_dR",{i+1,j+1});
                TString h_mll_name = ConstructHistName(ana_cat,"lep_mll",{i+1,j+1});
                FillHistIfExists(h_dR_name,dR,eft_fit);
                FillHistIfExists(h_mll_name,mll,eft_fit);
            }
        }
        if(ana_cat == "3l-sfz-1b" || ana_cat == "3l-sfz-2b") {
           int lep1 = -1;
           int lep2 = -1;
           for(size_t i = 0; i < pl_leptons.size(); i++) {
              for(size_t j = 0; j < pl_leptons.size(); j++) {
                 if (fabs((pl_leptons.at(i).p4() + pl_leptons.at(j).p4()).M() - 91.2) < 10) {
                    lep1 = i;
                    lep2 = j;
                    break; //done searching
                 }                                                                                                                                                             
              }
           }
           TString h_3l_sfz_Zpt = ConstructHistName(ana_cat, "Zpt", {});
           FillHistIfExists(h_3l_sfz_Zpt, pl_leptons.at(lep1).p4().Pt()+pl_leptons.at(lep2).p4().Pt(), eft_fit);
           double s = GetCosThetaStar(pl_leptons, lep1, lep2);
           TString h_3l_sfz_cos = ConstructHistName(ana_cat, "cos", {});
           FillHistIfExists(h_3l_sfz_cos, s, eft_fit);
        }
            
    }
    
    //////////////////////////////////////////

    // These were some hists made for the PL vs RECO checks. Are they still useful?

    // Testing analysis catetory yield hists (for PL vs RECO)
    // Yield hist
    TString ana_cat_hist_name = ConstructHistName(ana_cat_name,"yield",{});
    FillHistIfExists(ana_cat_hist_name,0.5,eft_fit);
    // Yield njets hists
    TString h_pl_njet_name = ConstructHistName(ana_cat_name,"yield-njets",{});
    FillHistIfExists(h_pl_njet_name,GetNJetsForLepCat(pl_leptons,pl_jets),eft_fit);
    // N events passing hist
    TString h_ana_cat_pass_name = ConstructHistName(ana_cat_name,"n-events-pass",{});
    FillTH1DHistIfExists(h_ana_cat_pass_name,0.5);

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


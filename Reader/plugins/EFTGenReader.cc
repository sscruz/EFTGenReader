// Created by Andrew Wightman

#include "EFTGenReader/Reader/interface/EFTGenReader.h"

EFTGenReader::EFTGenReader(const edm::ParameterSet& constructparams)
{
    debug = constructparams.getParameter<bool> ("debug");
    iseft = constructparams.getParameter<bool> ("iseft");
    gp_events = constructparams.getParameter<int> ("gp_events");
    norm_type = constructparams.getParameter<int> ("norm_type");
    xsec_norm = constructparams.getParameter<double> ("xsec_norm");
    intg_lumi = constructparams.getParameter<double> ("intg_lumi");

    entire_pset = constructparams;
    parse_params(); // Currently doesn't do anything
    lheInfo_token_ = consumes<LHEEventProduct>(edm::InputTag("externalLHEProducer"));
    genInfo_token_ = consumes<GenEventInfoProduct>(edm::InputTag("generator"));
    genParticles_token_ = consumes<reco::GenParticleCollection>(edm::InputTag("prunedGenParticles"));
    slimmed_genJets_token_ = consumes<std::vector<reco::GenJet> >(edm::InputTag("slimmedGenJets"));
    slimmed_genJetsAK8_token_ = consumes<std::vector<reco::GenJet> >(edm::InputTag("slimmedGenJetsAK8"));
    slimmed_genJetsAK8SoftDropSubJets_token_ = consumes<std::vector<reco::GenJet> >(edm::InputTag("slimmedGenJetsAK8SoftDropSubJets"));
    //genPackedParticles_token_ = consumes<std::vector<pat::PackedGenParticleCollection> >(edm::InputTag("packedGenParticles"));
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
    int pt_bins = 200;
    int eta_bins = 100;

    // Book the histograms that we will fill in the event loop

    // lep pt
    h_lep_pt = newfs->make<TH1D>("h_lep_pt","h_lep_pt",pt_bins,0,300);
    h_lep_ptSM = newfs->make<TH1D>("h_lep_ptSM","h_lep_ptSM",pt_bins,0,300);
    th1d_hists.push_back(h_lep_pt);
    th1d_hists.push_back(h_lep_ptSM);
    
    h_lep1_pt = newfs->make<TH1D>("h_lep1_pt","h_lep1_pt",pt_bins,0,300);
    h_lep1_ptSM = newfs->make<TH1D>("h_lep1_ptSM","h_lep1_ptSM",pt_bins,0,300);
    th1d_hists.push_back(h_lep1_pt);
    th1d_hists.push_back(h_lep1_ptSM);
    
    h_lep2_pt = newfs->make<TH1D>("h_lep2_pt","h_lep2_pt",pt_bins,0,300);
    h_lep2_ptSM = newfs->make<TH1D>("h_lep2_ptSM","h_lep2_ptSM",pt_bins,0,300);
    th1d_hists.push_back(h_lep2_pt);
    th1d_hists.push_back(h_lep2_ptSM);
    
    h_lepSum_pt = newfs->make<TH1D>("h_lepSum_pt","h_lepSum_pt",pt_bins,0,300);
    h_lepSum_ptSM = newfs->make<TH1D>("h_lepSum_ptSM","h_lepSum_ptSM",pt_bins,0,300);
    th1d_hists.push_back(h_lepSum_pt);
    th1d_hists.push_back(h_lepSum_ptSM);

    // lep eta
    h_lep_eta = newfs->make<TH1D>("h_lep_eta","h_lep_eta",eta_bins,-5.0,5.0);
    h_lep_etaSM = newfs->make<TH1D>("h_lep_etaSM","h_lep_etaSM",eta_bins,-5.0,5.0);
    th1d_hists.push_back(h_lep_eta);
    th1d_hists.push_back(h_lep_etaSM);

    // dilep inv mass
    h_mll = newfs->make<TH1D>("h_mll","h_mll",pt_bins,0,300);
    h_mllSM = newfs->make<TH1D>("h_mllSM","h_mllSM",pt_bins,0,300);
    th1d_hists.push_back(h_mll);
    th1d_hists.push_back(h_mllSM);

    // Jet histograms
    h_nJets = newfs->make<TH1D>("h_njets","h_njets",16,0,15);
    h_nJetsSM = newfs->make<TH1D>("h_njetsSM","h_njetsSM",16,0,15);
    th1d_hists.push_back(h_nJets);
    th1d_hists.push_back(h_nJetsSM);

    h_jet_pt = newfs->make<TH1D>("h_jet_pt","h_jet_pt",pt_bins,0,300);
    h_jet_ptSM = newfs->make<TH1D>("h_jet_ptSM","h_jet_ptSM",pt_bins,0,300);
    th1d_hists.push_back(h_jet_pt);
    th1d_hists.push_back(h_jet_ptSM);

    h_jet1_pt = newfs->make<TH1D>("h_jet1_pt","h_jet1_pt",pt_bins,0,300);
    h_jet1_ptSM = newfs->make<TH1D>("h_jet1_ptSM","h_jet1_ptSM",pt_bins,0,300);
    th1d_hists.push_back(h_jet1_pt);
    th1d_hists.push_back(h_jet1_ptSM);

    h_jet2_pt = newfs->make<TH1D>("h_jet2_pt","h_jet2_pt",pt_bins,0,300);
    h_jet2_ptSM = newfs->make<TH1D>("h_jet2_ptSM","h_jet2_ptSM",pt_bins,0,300);
    th1d_hists.push_back(h_jet2_pt);
    th1d_hists.push_back(h_jet2_ptSM);

    h_jet_eta = newfs->make<TH1D>("h_jet_eta","h_jet_eta",eta_bins,-5.0,5.0);
    h_jet_etaSM = newfs->make<TH1D>("h_jet_etaSM","h_jet_etaSM",eta_bins,-5.0,5.0);
    th1d_hists.push_back(h_jet_eta);
    th1d_hists.push_back(h_jet_etaSM);

    // pdgId histograms
    h_pdgId = newfs->make<TH1D>("h_pdgId","h_pdgId",2*pdg_bins,-pdg_bins+1,pdg_bins);
    h_pdgIdSM = newfs->make<TH1D>("h_pdgIdSM","h_pdgIdSM",2*pdg_bins,-pdg_bins+1,pdg_bins);
    th1d_hists.push_back(h_pdgId);
    th1d_hists.push_back(h_pdgIdSM);

    h_pdgIdFromZ = newfs->make<TH1D>("h_pdgIdFromZ","h_pdgIdFromZ",2*pdg_bins,-pdg_bins+1,pdg_bins);
    h_pdgIdFromZSM = newfs->make<TH1D>("h_pdgIdFromZSM","h_pdgIdFromZSM",2*pdg_bins,-pdg_bins+1,pdg_bins);
    th1d_hists.push_back(h_pdgIdFromZ);
    th1d_hists.push_back(h_pdgIdFromZSM);

    h_pdgIdFromH = newfs->make<TH1D>("h_pdgIdFromH","h_pdgIdFromH",2*pdg_bins,-pdg_bins+1,pdg_bins);
    h_pdgIdFromHSM = newfs->make<TH1D>("h_pdgIdFromHSM","h_pdgIdFromHSM",2*pdg_bins,-pdg_bins+1,pdg_bins);
    th1d_hists.push_back(h_pdgIdFromH);
    th1d_hists.push_back(h_pdgIdFromHSM);

    h_pdgIdFromW = newfs->make<TH1D>("h_pdgIdFromW","h_pdgIdFromW",2*pdg_bins,-pdg_bins+1,pdg_bins);
    h_pdgIdFromWSM = newfs->make<TH1D>("h_pdgIdFromWSM","h_pdgIdFromWSM",2*pdg_bins,-pdg_bins+1,pdg_bins);
    th1d_hists.push_back(h_pdgIdFromW);
    th1d_hists.push_back(h_pdgIdFromWSM);

    h_pdgIdFromPhoton = newfs->make<TH1D>("h_pdgIdFromPhoton","h_pdgIdFromPhoton",2*pdg_bins,-pdg_bins+1,pdg_bins);
    h_pdgIdFromPhotonSM = newfs->make<TH1D>("h_pdgIdFromPhotonSM","h_pdgIdFromPhotonSM",2*pdg_bins,-pdg_bins+1,pdg_bins);
    th1d_hists.push_back(h_pdgIdFromPhoton);
    th1d_hists.push_back(h_pdgIdFromPhotonSM);

    h_pdgIdFromGluon = newfs->make<TH1D>("h_pdgIdFromGluon","h_pdgIdFromGluon",2*pdg_bins,-pdg_bins+1,pdg_bins);
    h_pdgIdFromGluonSM = newfs->make<TH1D>("h_pdgIdFromGluonSM","h_pdgIdFromGluonSM",2*pdg_bins,-pdg_bins+1,pdg_bins);
    th1d_hists.push_back(h_pdgIdFromGluon);
    th1d_hists.push_back(h_pdgIdFromGluonSM);

    h_pdgIdFromQCD = newfs->make<TH1D>("h_pdgIdFromQCD","h_pdgIdFromQCD",2*pdg_bins,-pdg_bins+1,pdg_bins);
    h_pdgIdFromQCDSM = newfs->make<TH1D>("h_pdgIdFromQCDSM","h_pdgIdFromQCDSM",2*pdg_bins,-pdg_bins+1,pdg_bins);
    th1d_hists.push_back(h_pdgIdFromQCD);
    th1d_hists.push_back(h_pdgIdFromQCDSM);

    h_pdgIdFromOther = newfs->make<TH1D>("h_pdgIdFromOther","h_pdgIdFromOther",2*pdg_bins,-pdg_bins+1,pdg_bins);
    h_pdgIdFromOtherSM = newfs->make<TH1D>("h_pdgIdFromOtherSM","h_pdgIdFromOtherSM",2*pdg_bins,-pdg_bins+1,pdg_bins);
    th1d_hists.push_back(h_pdgIdFromOther);
    th1d_hists.push_back(h_pdgIdFromOtherSM);

    h_pdgIdLepMother = newfs->make<TH1D>("h_pdgIdLepMother","h_pdgIdLepMother",2*pdg_bins,-pdg_bins+1,pdg_bins);
    h_pdgIdLepMotherSM = newfs->make<TH1D>("h_pdgIdLepMotherSM","h_pdgIdLepMotherSM",2*pdg_bins,-pdg_bins+1,pdg_bins);
    th1d_hists.push_back(h_pdgIdLepMother);
    th1d_hists.push_back(h_pdgIdLepMotherSM);

    h_pdgIdLepGrMother = newfs->make<TH1D>("h_pdgIdLepGrMother","h_pdgIdLepGrMother",2*pdg_bins,-pdg_bins+1,pdg_bins);
    h_pdgIdLepGrMotherSM = newfs->make<TH1D>("h_pdgIdLepGrMotherSM","h_pdgIdLepGrMotherSM",2*pdg_bins,-pdg_bins+1,pdg_bins);
    th1d_hists.push_back(h_pdgIdLepGrMother);
    th1d_hists.push_back(h_pdgIdLepGrMotherSM);

    // misc. observables
    h_deltaR = newfs->make<TH1D>("h_deltaR","h_deltaR",100,0,5);
    h_deltaRSM = newfs->make<TH1D>("h_deltaRSM","h_deltaRSM",100,0,5);
    th1d_hists.push_back(h_deltaR);
    th1d_hists.push_back(h_deltaRSM);
    
    h_prompt_leptons = newfs->make<TH1D>("h_prompt_leptons","h_prompt_leptons",30,0,11);
    h_prompt_leptonsSM = newfs->make<TH1D>("h_prompt_leptonsSM","h_prompt_leptonsSM",30,0,11);
    th1d_hists.push_back(h_prompt_leptons);
    th1d_hists.push_back(h_prompt_leptonsSM);
    
    // Don't normalize these plots
    h_SMwgt_norm = newfs->make<TH1D>("h_SMwgt_norm","h_SMwgt_norm",350,-0.1,2.0);

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

    int nbins; double value; double norm;
    //for (uint i = 0; i < th1d_hists.size(); i++) {
    //    nbins = th1d_hists.at(i)->
    //}

    // Normalize all the plots to unit area
    for (TH1D* h: th1d_hists) {
        nbins = h->GetNbinsX();
        value = h->Integral(0,nbins+1);
        if (value == 0) value = 1.;
        if (norm_type == 1) {// Unit normalize        
            h->Scale(1./value);
        } else if (norm_type == 2) {
            // xsec_norm should be the NLO SM xsec for the sample
            if (iseft) {
                norm = intg_lumi*xsec_norm*total_sm_xsec / total_sm_xsec;
            } else {
                //norm = intg_lumi*xsec_norm / total_events;
                norm = intg_lumi*total_sm_xsec / total_events;
            }
            //norm = norm / value;
            h->Scale(norm);
        }
    }
}

void EFTGenReader::analyze(const edm::Event& event, const edm::EventSetup& evsetup)
{
    //std::string sm_wgt_str1 = "EFTrwgt183_ctW_0.0_ctp_0.0_cpQM_0.0_ctli_0.0_cQei_0.0_ctZ_0.0_cQlMi_0.0_cQl3i_0.0_ctG_0.0_ctlTi_0.0_cbW_0.0_cpQ3_0.0_ctei_0.0_cpt_0.0_ctlSi_0.0_cptb_0.0";
    //std::string sm_wgt_str2 = "EFTrwgt183_ctW_0.0_ctp_0.0_cpQM_0.0_cpt_0.0_cQei_0.0_ctZ_0.0_cQlMi_0.0_cQl3i_0.0_ctG_0.0_ctlTi_0.0_cbW_0.0_cpQ3_0.0_ctei_0.0_ctli_0.0_ctlSi_0.0_cptb_0.0";
    eventcount++;

    //std::cout << "EVENTNUMBER: " << eventcount << std::endl;

    // set tree vars to default values
    initialize_variables();

    edm::Handle<LHEEventProduct> LHEInfo;
    edm::Handle<reco::GenParticleCollection> prunedParticles;
    //edm::Handle<pat::PackedGenParticle> packedParticles;

    edm::Handle<std::vector<reco::GenJet> > slimGenJets;
    edm::Handle<std::vector<reco::GenJet> > slimGenJetsAK8;
    edm::Handle<std::vector<reco::GenJet> > slimGenJetsAK8SoftDropSubJets;

    event.getByToken(lheInfo_token_,LHEInfo);
    event.getByToken(genParticles_token_,prunedParticles);
    //event.getByToken(genPackedParticles_token_,packedParticles);

    event.getByToken(slimmed_genJets_token_,slimGenJets);
    event.getByToken(slimmed_genJetsAK8_token_,slimGenJetsAK8);
    event.getByToken(slimmed_genJetsAK8SoftDropSubJets_token_,slimGenJetsAK8SoftDropSubJets);


    if (debug){    
        //if (eventcount == 3) dumpParticles(*prunedParticles);
        if (eventcount == 3) dumpJets(*slimGenJets);
    }

    originalXWGTUP_intree = LHEInfo->originalXWGTUP();  // original cross-section
    double sm_wgt = 0.;
    double norm_sm_wgt = -1.;
    if (iseft) {// Add EFT weights 
        for (auto wgt_info: LHEInfo->weights()) {
            auto LHEwgtstr = std::string(wgt_info.id);
            std::size_t foundstr = LHEwgtstr.find("EFTrwgt"); // only save our EFT weights
            if (foundstr!=std::string::npos) {
                //eftwgts_intree[wgt_info.id] = wgt_info.wgt;
                //if (eventcount == 1) {
                //    std::cout << wgt_info.id << std::endl;
                //}
                //if (LHEwgtstr == sm_wgt_str1 || LHEwgtstr == sm_wgt_str2) {
                //    norm_sm_wgt = wgt_info.wgt / originalXWGTUP_intree;
                //    sm_wgt = wgt_info.wgt;
                //}
                WCPoint wc_pt(wgt_info.id);
                if (wc_pt.isSMPoint()) {
                    norm_sm_wgt = wgt_info.wgt / originalXWGTUP_intree;
                    sm_wgt = wgt_info.wgt;
                    break;
                }
            }
        }
    } else {
        norm_sm_wgt = 1.0;
        sm_wgt = originalXWGTUP_intree;
    }

    //if (eventcount % 25 == 1) {
    //    std::cout << "EVENTNUMBER: " << eventcount << std::endl;
    //    std::cout << "\tsm_wgt:   " << sm_wgt << " (" << norm_sm_wgt << ")" << std::endl;
    //    std::cout << "\torig_wgt: " << originalXWGTUP_intree << std::endl;
    //}

    total_sm_xsec += sm_wgt;
    total_orig_xsec += originalXWGTUP_intree;

    //std::cout << "SM Wgt: " << sm_wgt << std::endl;
    //std::cout << "Orig Wgt: " << originalXWGTUP_intree << std::endl;
    //std::cout << "Norm Wgt: " << norm_sm_wgt << std::endl;
    //std::cout << std::endl;

    h_SMwgt_norm->Fill(norm_sm_wgt);

    reco::GenParticleCollection gen_leptons = GetGenLeptons(*prunedParticles);

    //if (gen_leptons.size() != 2) {
    //    std::cout << "EVENTNUMBER: " << eventcount << std::endl;
    //    dumpParticles(*prunedParticles);
    //}
    //dumpParticles(gen_leptons);

    h_prompt_leptons->Fill(gen_leptons.size());
    h_prompt_leptonsSM->Fill(gen_leptons.size(),norm_sm_wgt);

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

        //std::cout << "pdgId: " << id
        //          << "\n\tPt: " << pt << "\tEta: " << eta
        //          << "\n\tStatus: " << st
        //          << "\n\tDirect pdgId:  " << direct_id
        //          << "\n\tMother pdgId:  " << mom_id
        //          << "\n\tGMother pdgId: " << gmom_id << std::endl;

        h_pdgId->Fill(id);                 h_pdgIdSM->Fill(id,norm_sm_wgt);
        h_lep_pt->Fill(pt);                h_lep_ptSM->Fill(pt,norm_sm_wgt);
        h_lep_eta->Fill(eta);              h_lep_etaSM->Fill(eta,norm_sm_wgt);
        h_pdgIdLepMother->Fill(mom_id);    h_pdgIdLepMotherSM->Fill(mom_id,norm_sm_wgt);
        h_pdgIdLepGrMother->Fill(gmom_id); h_pdgIdLepGrMotherSM->Fill(gmom_id,norm_sm_wgt);

        // mom_id and direct_id should be the same thing for the gen leptons selected by GetGenLeptons()
        if (mom_id == 21) {
            h_pdgIdFromGluon->Fill(id);
            h_pdgIdFromGluonSM->Fill(id,norm_sm_wgt);
        } else if (mom_id >= 1 && mom_id <= 5) {
            h_pdgIdFromQCD->Fill(id);
            h_pdgIdFromQCDSM->Fill(id,norm_sm_wgt);  
        } else if (mom_id == 22) {
            h_pdgIdFromPhoton->Fill(id);
            h_pdgIdFromPhotonSM->Fill(id,norm_sm_wgt);
        } else if (mom_id == 23) {
            h_pdgIdFromZ->Fill(id);
            h_pdgIdFromZSM->Fill(id,norm_sm_wgt);
        } else if (mom_id == -24 || mom_id == 24) {
            h_pdgIdFromW->Fill(id);
            h_pdgIdFromWSM->Fill(id,norm_sm_wgt);
        } else if (mom_id == 25) {
            h_pdgIdFromH->Fill(id);
            h_pdgIdFromHSM->Fill(id,norm_sm_wgt);
        } else {
            h_pdgIdFromOther->Fill(id);
            h_pdgIdFromOtherSM->Fill(id,norm_sm_wgt);
        }

        if (i == 0) {
            h_lep1_pt->Fill(pt);
            h_lep1_ptSM->Fill(pt,norm_sm_wgt);
        } else if (i == 1) {
            h_lep2_pt->Fill(pt);
            h_lep2_ptSM->Fill(pt,norm_sm_wgt);
        }
        for (size_t j = 0; j < gen_leptons.size(); j++) {
            if (i >= j) continue;
            //std::cout << "\t(i,j): " << "(" << i << "," << j << ")" << std::endl;

            const reco::GenParticle& p2 = gen_leptons.at(j);
            double dR = getdR(p1,p2);
            double mll = getInvMass(p1,p2);
            auto p4vec = p1.p4() + p2.p4();
            double sum_pt = p4vec.Pt();

            //std::cout << "\t\tdR: " << dR
            //          << "\n\t\tMll: " << mll << std::endl;

            h_mll->Fill(mll);          h_mllSM->Fill(mll,norm_sm_wgt);
            h_deltaR->Fill(dR);        h_deltaRSM->Fill(dR,norm_sm_wgt);
            h_lepSum_pt->Fill(sum_pt); h_lepSum_ptSM->Fill(sum_pt,norm_sm_wgt);
        }
    }

    std::vector<reco::GenJet> gen_jets = GetGenJets(*slimGenJets);

    h_nJets->Fill(gen_jets.size());
    h_nJetsSM->Fill(gen_jets.size(),norm_sm_wgt);

    for (size_t i = 0; i < gen_jets.size(); i++) {
        const reco::GenJet& p1 = gen_jets.at(i);
        double pt = p1.p4().Pt();
        double eta = p1.p4().Eta();

        h_jet_pt->Fill(pt);   h_jet_ptSM->Fill(pt,norm_sm_wgt);
        h_jet_eta->Fill(eta); h_jet_etaSM->Fill(eta,norm_sm_wgt);
        if (i == 0) {
            h_jet1_pt->Fill(pt);
            h_jet1_ptSM->Fill(pt,norm_sm_wgt);
        } else if (i == 1) {
            h_jet2_pt->Fill(pt);
            h_jet2_ptSM->Fill(pt,norm_sm_wgt);
        }
    }

    eventnum_intree = event.id().event();
    lumiBlock_intree = event.id().luminosityBlock();
    runNumber_intree = event.id().run();

    pruned_genParticles_intree = *prunedParticles;

    slimmed_genJets_intree = *slimGenJets;
    slimmed_genJetsAK8_intree = *slimGenJetsAK8;
    slimmed_genJetsAK8SoftDropSubJets_intree = *slimGenJetsAK8SoftDropSubJets;

    //if (eventcount < 10) {
    //    summaryTree->Fill();
    //}

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
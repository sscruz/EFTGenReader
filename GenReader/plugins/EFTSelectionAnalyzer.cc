// Created by Andrew Wightman

#include "EFTSelectionAnalyzer.h"

EFTSelectionAnalyzer::EFTSelectionAnalyzer(const edm::ParameterSet& iConfig)
{
    debug = iConfig.getParameter<bool> ("debug");
    iseft = iConfig.getParameter<bool> ("iseft");

    min_pt_jet  = iConfig.getParameter<double> ("min_pt_jet");
    min_pt_lep  = iConfig.getParameter<double> ("min_pt_lep");
    max_eta_jet = iConfig.getParameter<double> ("max_eta_jet");
    max_eta_lep = iConfig.getParameter<double> ("max_eta_lep");
    
    lheInfo_token_      = consumes<LHEEventProduct>(iConfig.getParameter<edm::InputTag>("LHEInfo"));
    genInfo_token_      = consumes<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("GENInfo"));
    genParticles_token_ = consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("GenParticles"));
    genJets_token_      = consumes<std::vector<reco::GenJet> >(iConfig.getParameter<edm::InputTag>("GenJets"));
}

EFTSelectionAnalyzer::~EFTSelectionAnalyzer(){}

void EFTSelectionAnalyzer::beginJob()
{
    total_ls = 0;
    total_events = 0;
    total_orig_xsec = 0.;
    total_sm_xsec = 0.;

    fit_all_incl  = new WCFit;
    fit_2lss_incl = new WCFit;
    fit_2lss_ee   = new WCFit;
    fit_2lss_emu  = new WCFit;
    fit_2lss_mumu = new WCFit;
}

void EFTSelectionAnalyzer::endJob()
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
        << std::setw(13) << total_orig_xsec / double(total_ls) << delim
        << std::setw(13) << total_orig_xsec / double(total_events) << delim
        << std::setw(13) << total_orig_xsec / 500.
        << std::endl;

    WCPoint* rwgt_pt = new WCPoint("rwgt",0.0);

    fit_all_incl->setTag("All Events"); 
    fit_2lss_incl->setTag("2lss Events"); 
    fit_2lss_ee->setTag("2lss ee Events"); 
    fit_2lss_emu->setTag("2lss emu Events"); 
    fit_2lss_mumu->setTag("2lss mumu Events"); 

    std::cout << fit_all_incl->getTag()  << " SM Eval: " << fit_all_incl->evalPoint(rwgt_pt)  << std::endl;
    std::cout << fit_2lss_incl->getTag() << " SM Eval: " << fit_2lss_incl->evalPoint(rwgt_pt) << std::endl;
    std::cout << fit_2lss_ee->getTag()   << " SM Eval: " << fit_2lss_ee->evalPoint(rwgt_pt)   << std::endl;
    std::cout << fit_2lss_emu->getTag()  << " SM Eval: " << fit_2lss_emu->evalPoint(rwgt_pt)  << std::endl;
    std::cout << fit_2lss_mumu->getTag() << " SM Eval: " << fit_2lss_mumu->evalPoint(rwgt_pt) << std::endl;

    delete fit_all_incl;
    delete fit_2lss_incl;
    delete fit_2lss_ee;
    delete fit_2lss_emu;
    delete fit_2lss_mumu;
    delete rwgt_pt;
}

void EFTSelectionAnalyzer::analyze(const edm::Event& event, const edm::EventSetup& evsetup)
{
    eventcount++;

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
    //std::vector<reco::GenJet> gen_jets = GetGenJets(*genJets);

    double sm_wgt = 0.;
    double orig_wgt = LHEInfo->originalXWGTUP();  // original cross-section
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
        sm_wgt = orig_wgt;
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

    total_sm_xsec += sm_wgt;
    total_orig_xsec += orig_wgt;


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

        if (abs(id) == 11) { // electron
            n_prompt_electrons += 1;
        } else if (abs(id) == 13) { // muon
            n_prompt_muons += 1;
        } else if (abs(id) == 13) { // tau
            n_prompt_taus += 1;
        }
    }
    fit_all_incl->addFit(eft_fit);
    if (gen_leptons.size() == 2){
        fit_2lss_incl->addFit(eft_fit);
        if (n_prompt_electrons == 2) {
            fit_2lss_ee->addFit(eft_fit);
        } else if (n_prompt_electrons == 1 and n_prompt_muons ==1){
            fit_2lss_emu->addFit(eft_fit);
        } else if (n_prompt_muons == 2){
            fit_2lss_mumu->addFit(eft_fit);
        }
    }

}

void EFTSelectionAnalyzer::beginRun(edm::Run const& run, edm::EventSetup const& evsetup)
{
    eventcount = 0;
}

void EFTSelectionAnalyzer::endRun(edm::Run const& run, edm::EventSetup const& evsetup)
{
    total_events += eventcount;
}

void EFTSelectionAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const& lumi, edm::EventSetup const& evsetup)
{
    total_ls += 1;
}

void EFTSelectionAnalyzer::endLuminosityBlock(edm::LuminosityBlock const& lumi, edm::EventSetup const& evsetup){}

DEFINE_FWK_MODULE(EFTSelectionAnalyzer);


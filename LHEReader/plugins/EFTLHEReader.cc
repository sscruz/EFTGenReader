// Created by Andrew Wightman

//#include "EFTGenReader/LHEReader/interface/EFTLHEReader.h"
#include "EFTLHEReader.h"

EFTLHEReader::EFTLHEReader(const edm::ParameterSet& iConfig)
{
    min_pt_jet = iConfig.getParameter<double> ("min_pt_jet");
    min_pt_lep = iConfig.getParameter<double> ("min_pt_lep");
    max_eta_jet = iConfig.getParameter<double> ("max_eta_jet");
    max_eta_lep = iConfig.getParameter<double> ("max_eta_lep");

    //entire_pset = iConfig;
    
    parse_params(); // Currently doesn't do anything
    //lheInfo_token_      = consumes<LHEEventProduct>(edm::InputTag("externalLHEProducer"));
    //genInfo_token_      = consumes<GenEventInfoProduct>(edm::InputTag("generator")); // Associating the token with a moduel form the edm file
    //genParticles_token_ = consumes<std::vector<reco::GenParticle> >(edm::InputTag("genParticles"));
    //genJets_token_      = consumes<std::vector<reco::GenJet> >(edm::InputTag("ak4GenJets")); // or ak8GenJets?
    lheInfo_token_      = consumes<LHEEventProduct>(iConfig.getParameter<edm::InputTag>("LHEInfo"));
    genInfo_token_      = consumes<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("GENInfo"));
    genParticles_token_ = consumes<reco::GenParticleCollection>(iConfig.getParameter<edm::InputTag>("GenParticles"));
    genJets_token_      = consumes<std::vector<reco::GenJet> >(iConfig.getParameter<edm::InputTag>("GenJets"));
}

EFTLHEReader::~EFTLHEReader(){}

void EFTLHEReader::beginJob()
{
    edm::Service<TFileService> newfs;
    summaryTree = newfs->make<TTree>("summaryTree","Summary Event Values");

    tree_add_branches();
}

void EFTLHEReader::endJob(){}

void EFTLHEReader::analyze(const edm::Event& event, const edm::EventSetup& evsetup)
{
    eventcount++;

    // set tree vars to default values
    initialize_variables();

    edm::Handle<LHEEventProduct> LHEInfo;
    event.getByToken(lheInfo_token_,LHEInfo);

    edm::Handle<GenEventInfoProduct> GENInfo;
    event.getByToken(genInfo_token_,GENInfo); // GENInfo will store the DJR values

    edm::Handle<std::vector<reco::GenParticle> > genParticles;
    event.getByToken(genParticles_token_, genParticles);

    edm::Handle<std::vector<reco::GenJet> > genJets;
    event.getByToken(genJets_token_, genJets);

    originalXWGTUP_intree = LHEInfo->originalXWGTUP();  // original cross-section

    // Add EFT weights
    std::vector<WCPoint> wc_pts;
    for (auto wgt_info: LHEInfo->weights())
    {
        auto LHEwgtstr = std::string(wgt_info.id);
        std::size_t foundstr = LHEwgtstr.find("EFTrwgt"); // only save our EFT weights
        if (foundstr!=std::string::npos) {
            eftwgts_intree[wgt_info.id] = wgt_info.wgt;
            WCPoint wc_pt(wgt_info.id,wgt_info.wgt);
            wc_pts.push_back(wc_pt);
        }
    }

    WCFit event_fit(wc_pts,"");

    std::vector<reco::GenParticle> gen_leptons = GetGenLeptons(*genParticles); // Note GetGenLeptons sorts by pt
    for (size_t i = 0; i < gen_leptons.size(); i++) {
        const reco::GenParticle& p = gen_leptons.at(i);
        double pt = p.p4().Pt();

        if (i==0) {
           genLep_pt1_intree = pt;
        } else if (i==1) {
            genLep_pt2_intree = pt;
        } else if (i==2) {
            genLep_pt3_intree = pt;
        }
    }

    std::vector<reco::GenJet> gen_jets = GetGenJets(*genJets); //
    for (size_t i = 0; i < gen_jets.size(); i++) {
        const reco::GenJet& p = gen_jets.at(i);
        double pt = p.p4().Pt();

        if (i==0) {
           genJet_pt1_intree = pt;
        } else if (i==1) {
            genJet_pt2_intree = pt;
        } else if (i==2) {
            genJet_pt3_intree = pt;
        } else if (i==3) {
            genJet_pt4_intree = pt;
        }
    }

    eventnum_intree = event.id().event();
    lumiBlock_intree = event.id().luminosityBlock();
    runNumber_intree = event.id().run();
    wcFit_intree = event_fit;

    nMEpartons_intree         = GENInfo->nMEPartons();
    nMEpartonsFiltered_intree = GENInfo->nMEPartonsFiltered();
    genWgt_intree             = GENInfo->weight();

    djrvalues_intree.clear();
    for (auto djr_val: GENInfo->DJRValues()) {
        djrvalues_intree.push_back(djr_val); // Change to a double
    }

    summaryTree->Fill();
}

void EFTLHEReader::beginRun(edm::Run const& run, edm::EventSetup const& evsetup)
{
    eventcount = 0;
}

void EFTLHEReader::endRun(edm::Run const& run, edm::EventSetup const& evsetup)
{
    std::cout << "total events processed: " << eventcount << std::endl;
}

void EFTLHEReader::beginLuminosityBlock(edm::LuminosityBlock const& lumi, edm::EventSetup const& evsetup){}
void EFTLHEReader::endLuminosityBlock(edm::LuminosityBlock const& lumi, edm::EventSetup const& evsetup){}

DEFINE_FWK_MODULE(EFTLHEReader);

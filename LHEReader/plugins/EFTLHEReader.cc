// Created by Andrew Wightman

#include "EFTGenReader/LHEReader/interface/EFTLHEReader.h"

EFTLHEReader::EFTLHEReader(const edm::ParameterSet& constructparams)
{
    entire_pset = constructparams;
    parse_params(); // Currently doesn't do anything
    lheInfo_token_ = consumes<LHEEventProduct>(edm::InputTag("externalLHEProducer"));
    genInfo_token_ = consumes<GenEventInfoProduct>(edm::InputTag("generator")); // Associating the token with a moduel form the edm file
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

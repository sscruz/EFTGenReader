// Created by Andrew Wightman

//#include "EFTGenReader/LHEReader/interface/EFTLHEReader.h"
#include "EFTLHEReader.h"

EFTLHEReader::EFTLHEReader(const edm::ParameterSet& iConfig)
{
    min_pt_jet = iConfig.getParameter<double> ("min_pt_jet");
    min_pt_lep = iConfig.getParameter<double> ("min_pt_lep");
    max_eta_jet = iConfig.getParameter<double> ("max_eta_jet");
    max_eta_lep = iConfig.getParameter<double> ("max_eta_lep");

    is_4f_scheme = iConfig.getParameter<bool>("is4fScheme");

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

    consumes<LHERunInfoProduct, edm::InRun>(edm::InputTag("externalLHEProducer")); // Not sure what this does???
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

    // https://twiki.cern.ch/twiki/bin/viewauth/CMS/TopModGen#Event_Generation
    if (GENInfo->weights().size()==14) // test if PS weights present
    {
        std::vector<double> genwgtinfo = GENInfo->weights();
        preshowerISRweightUp_intree = genwgtinfo[2];
        preshowerFSRweightUp_intree = genwgtinfo[3];
        preshowerISRweightDown_intree = genwgtinfo[4];
        preshowerFSRweightDown_intree = genwgtinfo[5];
    }


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
        if ( LHEwgtstr.find("1006")!=std::string::npos ) muRWeightUp_intree      = wgt_info.wgt / originalXWGTUP_intree;
        if ( LHEwgtstr.find("1011")!=std::string::npos ) muRWeightDown_intree    = wgt_info.wgt / originalXWGTUP_intree;
        if ( LHEwgtstr.find("1016")!=std::string::npos ) muFWeightUp_intree      = wgt_info.wgt / originalXWGTUP_intree;
        if ( LHEwgtstr.find("1031")!=std::string::npos ) muFWeightDown_intree    = wgt_info.wgt / originalXWGTUP_intree;
        if ( LHEwgtstr.find("1021")!=std::string::npos ) muRmuFWeightUp_intree   = wgt_info.wgt / originalXWGTUP_intree;
        if ( LHEwgtstr.find("1041")!=std::string::npos ) muRmuFWeightDown_intree = wgt_info.wgt / originalXWGTUP_intree;
    }

    WCFit event_fit(wc_pts,"");

    //-------------
    // pdf weights and Q^2 weights
    //-------------
    std::vector<double> nnpdfWeights;
    // create the PDF
    TString pdf_set_name = "NNPDF31_nnlo_hessian_pdfas";
    if (is_4f_scheme) {
        pdf_set_name = "NNPDF31_nnlo_as_0118_nf_4";
    }
    LHAPDF::PDFSet nnpdfSet(pdf_set_name.Data());

    int pdfID_start = 306001;
    int pdfID_end = 306102;
    if (is_4f_scheme) {
        pdfID_start = 320901;
        pdfID_end = 321000;
    }

    // obtain weights
    auto& mcWeights = LHEInfo->weights();
    for (size_t i = 0; i < mcWeights.size(); i++)
    {
        // use the mapping to identify the weight
        auto wgtstr = std::string(mcWeights[i].id);
        std::size_t foundstr = wgtstr.find("rwgt");
        if ( foundstr!=std::string::npos ) continue; // skip EFT stuff here
        if (wgtstr=="dummy_point") continue;
        int idInt = stoi(mcWeights[i].id);
        if (pdfIdMap_.find(idInt) != pdfIdMap_.end())
        {
            int setId = pdfIdMap_[idInt];
            if (setId >= pdfID_start && setId <= pdfID_end) // NNPDF30_nlo_as_0118
            {
                // divide by original weight to get scale-factor-like number
                nnpdfWeights.push_back(mcWeights[i].wgt / originalXWGTUP_intree);
            }
        }
    }

    // create the combined up/down variations
    double weightUp = 1.;
    double weightDown = 1.;
    if (nnpdfWeights.size() > 0)
    {
        // in rare cases it might happen that not all weights are present, so in order to
        // use LHAPDF's uncertainty function, we fill up the vector
        // this is expected not to have a big impact
        while ((int)nnpdfWeights.size() < (pdfID_end - pdfID_start + 1))
        {
            nnpdfWeights.push_back(1.);
        }
        // first value must be the nominal value, i.e. 1 since we normalize by original weight
        nnpdfWeights.insert(nnpdfWeights.begin(), 1.);                                                                      // <<----------
        // calculate combined weights
        // std::cout << "nnpdfWeights size is: " << (int)nnpdfWeights.size() << std::endl;
        const LHAPDF::PDFUncertainty pdfUnc = nnpdfSet.uncertainty(nnpdfWeights, 68.268949);    // 2nd arg is the CL (in %) to rescale the uncertainties to, 68.268949 is the default
        weightUp = pdfUnc.central + pdfUnc.errplus;
        weightDown = pdfUnc.central - pdfUnc.errminus;
    }

    nnpdfWeightUp_intree = weightUp;
    nnpdfWeightDown_intree = weightDown;



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
    bool debug = false;
    eventcount = 0;
    /////////////////
    //// pdf weights mapping
    /////////////////
    edm::Handle<LHERunInfoProduct> runInfo;
    run.getByLabel("externalLHEProducer", runInfo);

    std::string weightTag,startStr,setStr,endStr;

    weightTag = "initrwgt";
    // The private samples have a slightly different formatting of the header
    startStr  = "&lt;weight id=";
    setStr    = " MUR=\"1.0\" MUF=\"1.0\" PDF=\"";
    endStr    = " &gt; PDF=306000"; // For 5f PDFs
    if (is_4f_scheme) {
        endStr    = " &gt; PDF=320900"; // For 4f PDFs
    }

    if (debug) std::cout << "before loop over pdf stuff in beginRun" << std::endl;
    for (std::vector<LHERunInfoProduct::Header>::const_iterator it = runInfo->headers_begin(); it != runInfo->headers_end(); it++) {
        if (it->tag() != weightTag) {
            //std::cout << it->tag() << std::endl;
            continue;
        }

        std::vector<std::string> lines = it->lines();
        for (size_t i = 0; i < lines.size(); i++) {
            //std::cout << lines[i];
            size_t startPos = lines[i].find(startStr);
            size_t setPos = lines[i].find(setStr);
            size_t endPos = lines[i].find(endStr);
            //std::cout << (startPos == std::string::npos) << std::endl;
            //std::cout << (setPos == std::string::npos)   << std::endl;
            //std::cout << (endPos == std::string::npos)   << std::endl;
            if (startPos == std::string::npos || setPos == std::string::npos || endPos == std::string::npos) {
                continue;
            }
            std::string weightId = lines[i].substr(startPos + startStr.size() + 1, setPos - startPos - startStr.size() - 2); // this has changed, modify it???
            std::string setId = lines[i].substr(setPos + setStr.size(), endPos - setPos - setStr.size() - 1); // this has changed, modify it???
            //std::cout << weightId << " : " << setId << std::endl;
            try {
                pdfIdMap_[stoi(weightId)] = stoi(setId);
            }
            catch (...) {
                std::cerr << "error while parsing the lhe run xml header: ";
                std::cerr << "cannot interpret as ints:" << weightId << " -> " << setId << std::endl;
            }
        }
    }

}

void EFTLHEReader::endRun(edm::Run const& run, edm::EventSetup const& evsetup)
{
    std::cout << "total events processed: " << eventcount << std::endl;
}

void EFTLHEReader::beginLuminosityBlock(edm::LuminosityBlock const& lumi, edm::EventSetup const& evsetup){}
void EFTLHEReader::endLuminosityBlock(edm::LuminosityBlock const& lumi, edm::EventSetup const& evsetup){}

DEFINE_FWK_MODULE(EFTLHEReader);

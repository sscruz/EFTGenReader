// created by Andrew Wightman

#include <cstdlib>
#include <memory>
#include <vector>
#include <string>
#include <map>
#include <unordered_map>
#include <boost/any.hpp>

#include <iostream>
#include <algorithm>
#include <exception>
#include <cmath> 
#include <iomanip>
#include <fstream>
#include <sstream>

#include "TROOT.h"
#include "TSystem.h"
#include "TStyle.h"
#include <TRandom3.h>
#include "TTree.h"
#include "TH1.h"
#include "TH2.h"
#include "TH2D.h"
#include "TVector.h"
#include "TLorentzVector.h"

// Framework
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "FWCore/ParameterSet/interface/ProcessDesc.h"
#include "FWCore/ParameterSet/interface/FileInPath.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
//#include "FWCore/PythonParameterSet/interface/PythonProcessDesc.h"
#include "FWCore/PythonParameterSet/interface/PyBind11ProcessDesc.h"

#include "FWCore/Framework/interface/Run.h"

// Physics
#include "LHAPDF/LHAPDF.h"
#include "Math/LorentzVector.h"
#include "SimDataFormats/GeneratorProducts/interface/LHERunInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"

//#include "EFTGenReader/GenReader/interface/WCPoint.h"
//#include "EFTGenReader/GenReader/interface/WCFit.h"
#include "EFTGenReader/EFTHelperUtilities/interface/WCPoint.h"
#include "EFTGenReader/EFTHelperUtilities/interface/WCFit.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/PatCandidates/interface/Jet.h"

// end includes
// -----------------------------------------------

class EFTLHEReader: public edm::EDAnalyzer
{
    private:
        // EDAnalyzer-specific:
        virtual void beginJob();
        virtual void analyze(const edm::Event&, const edm::EventSetup&);
        virtual void endJob();
        virtual void beginRun(edm::Run const&, edm::EventSetup const&);
        virtual void endRun(edm::Run const&, edm::EventSetup const&);
        virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);
        virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&);

        void tree_add_branches();
        void initialize_variables();
        void parse_params();

    public:
        explicit EFTLHEReader(const edm::ParameterSet&);
        ~EFTLHEReader();

        template <typename T> edm::Handle<T> get_collection(const edm::Event& event, const edm::EDGetTokenT<T>& token);
        std::vector<reco::GenJet> GetGenJets(const std::vector<reco::GenJet>& inputs);
        const reco::Candidate* GetGenMotherNoFsr(const reco::Candidate* p);
        reco::GenParticleCollection GetGenLeptons(const reco::GenParticleCollection& gen_particles);

        std::ofstream fout;
        FILE * ffout;

        std::string sampleName;
        int sampleNumber;
        int eventcount;

        //edm::ParameterSet entire_pset;

        edm::EDGetTokenT<LHEEventProduct> lheInfo_token_;
        edm::EDGetTokenT<GenEventInfoProduct> genInfo_token_; //GenEventInfoProduct name of the class (the token is of this type) (define token)
        edm::EDGetTokenT<std::vector<reco::GenParticle> > genParticles_token_;
        edm::EDGetTokenT<std::vector<reco::GenJet> > genJets_token_;

        std::map<int,int> pdfIdMap_;

        // declare the tree
        TTree * summaryTree;

        // tree branches
        std::unordered_map<std::string,double> eftwgts_intree;
        double originalXWGTUP_intree;
        int eventnum_intree;
        int lumiBlock_intree;
        int runNumber_intree;

        int nMEpartons_intree;
        int nMEpartonsFiltered_intree;
        double genWgt_intree;

        double preshowerISRweightUp_intree;
        double preshowerFSRweightUp_intree;
        double preshowerISRweightDown_intree;
        double preshowerFSRweightDown_intree;

        double muRWeightUp_intree;
        double muRWeightDown_intree;
        double muFWeightUp_intree;
        double muFWeightDown_intree;
        double muRmuFWeightUp_intree;
        double muRmuFWeightDown_intree;

        double genLep_pt1_intree;
        double genLep_pt2_intree;
        double genLep_pt3_intree;

        double genJet_pt1_intree;
        double genJet_pt2_intree;
        double genJet_pt3_intree;
        double genJet_pt4_intree;

        double nnpdfWeightUp_intree;
        double nnpdfWeightDown_intree;

        std::vector<double> djrvalues_intree; // Variabel to store the info we are getting from the module 

        WCFit wcFit_intree;

        // Kinematic and eta cuts
        double min_pt_jet;
        double min_pt_lep;
        double max_eta_jet;
        double max_eta_lep;

        bool is_4f_scheme;

};

void EFTLHEReader::tree_add_branches()
{
    summaryTree->Branch("eftwgts",&eftwgts_intree);
    summaryTree->Branch("originalXWGTUP",&originalXWGTUP_intree);

    summaryTree->Branch("eventnum",&eventnum_intree);
    summaryTree->Branch("lumiBlock",&lumiBlock_intree);
    summaryTree->Branch("runNumber",&runNumber_intree);

    summaryTree->Branch("nMEpartons",&nMEpartons_intree);
    summaryTree->Branch("nMEpartonsFiltered",&nMEpartonsFiltered_intree);
    summaryTree->Branch("genWgt",&genWgt_intree);
    summaryTree->Branch("DJRValues",&djrvalues_intree); // Add branch 

    summaryTree->Branch("wcFit",&wcFit_intree);

    summaryTree->Branch("genLep_pt1",&genLep_pt1_intree);
    summaryTree->Branch("genLep_pt2",&genLep_pt2_intree);
    summaryTree->Branch("genLep_pt3",&genLep_pt3_intree);

    summaryTree->Branch("genJet_pt1",&genJet_pt1_intree);
    summaryTree->Branch("genJet_pt2",&genJet_pt2_intree);
    summaryTree->Branch("genJet_pt3",&genJet_pt3_intree);
    summaryTree->Branch("genJet_pt4",&genJet_pt4_intree);

    summaryTree->Branch("psISRweightUp",&preshowerISRweightUp_intree);
    summaryTree->Branch("psFSRweightUp",&preshowerFSRweightUp_intree);
    summaryTree->Branch("psISRweightDown",&preshowerISRweightDown_intree);
    summaryTree->Branch("psFSRweightDown",&preshowerFSRweightDown_intree);

    summaryTree->Branch("muRWeightUp",&muRWeightUp_intree);
    summaryTree->Branch("muRWeightDown",&muRWeightDown_intree);
    summaryTree->Branch("muFWeightUp",&muFWeightUp_intree);
    summaryTree->Branch("muFWeightDown",&muFWeightDown_intree);
    summaryTree->Branch("muRmuFWeightUp",&muRmuFWeightUp_intree);
    summaryTree->Branch("muRmuFWeightDown",&muRmuFWeightDown_intree);

    summaryTree->Branch("nnpdfWeightUp",&nnpdfWeightUp_intree);
    summaryTree->Branch("nnpdfWeightDown",&nnpdfWeightDown_intree);

}

void EFTLHEReader::initialize_variables()
{
    eftwgts_intree.clear();
    originalXWGTUP_intree = -99;

    wcFit_intree.clear();

    genLep_pt1_intree = -999;
    genLep_pt2_intree = -999;
    genLep_pt3_intree = -999;
    genJet_pt1_intree = -999;
    genJet_pt2_intree = -999;
    genJet_pt3_intree = -999;
    genJet_pt4_intree = -999;
}

// Currently not setting any special parameters in config files
void EFTLHEReader::parse_params() {}

template<typename T>
edm::Handle<T>
EFTLHEReader::get_collection(const edm::Event& event, const edm::EDGetTokenT<T>& token)
{
    edm::Handle<T> handle;
    event.getByToken(token, handle);
    if (!handle.isValid())
        throw edm::Exception(edm::errors::InvalidReference, "Can't find a collection.");
    return handle;
}

// From EFTGenReader.h
reco::GenParticleCollection EFTLHEReader::GetGenLeptons(const reco::GenParticleCollection& gen_particles) {
    reco::GenParticleCollection gen_leptons;
    bool is_lepton;
    bool is_neutrino;
    for (size_t i = 0; i < gen_particles.size(); i++) {
        const reco::GenParticle& p = gen_particles.at(i);
        const reco::Candidate* p_mom_noFSR = GetGenMotherNoFsr(&p); // This will walk up the chain and so doesn't get direct mothers
        const reco::Candidate* p_gmom_noFSR = GetGenMotherNoFsr(p_mom_noFSR);
        const reco::Candidate* p_mom = p.mother();  // The direct mother
        int id = p.pdgId();
        int gmom_noFSR_id = p_gmom_noFSR->pdgId();
        is_lepton = (abs(id) == 11 || abs(id) == 13 || abs(id) == 15);
        is_neutrino = (abs(id) == 12 || abs(id) == 14 || abs(id) == 16);

        if (!is_lepton && !is_neutrino) {
            continue;
        }

        // Note: The kinematic cuts could realistically go anywhere in this loop
        if (p.p4().Pt() < min_pt_lep) {
            continue;
        } else if (max_eta_lep > 0.0 && fabs(p.eta()) >= max_eta_lep) {
            continue;
        }

        int mom_id = id;    // If no mother, set to own id
        if (p_mom) mom_id = p_mom->pdgId();

        bool is_fromEWBoson = (abs(mom_id) >= 22 && abs(mom_id) <= 25);
        bool is_fromTopSystem = (abs(mom_id) == 24 && abs(gmom_noFSR_id) == 6);

        //if (is_fromTopSystem) {
        //    continue;
        //}

        bool is_hard_process = p.isHardProcess();
        if (is_hard_process) {
            gen_leptons.push_back(p);
        }
        //bool is_fromHardGluon = (abs(mom_id) == 21 && p_mom->status() == 21);   // status == 21 corresponds to incoming hard particle
        //bool is_fromHardQCD = (abs(mom_id) >= 1 && abs(mom_id) <= 5 && p_mom->status() == 21);

        //if (is_fromEWBoson) {
        //    gen_leptons.push_back(p);
        //} else if (is_hard_process && is_fromHardGluon) {
        //    gen_leptons.push_back(p);
        //} else if (is_hard_process && is_fromHardQCD) {
        //    gen_leptons.push_back(p);
        //}
    }
    std::sort(gen_leptons.begin(),gen_leptons.end(), [] (reco::GenParticle a, reco::GenParticle b) { return a.p4().Pt() > b.p4().Pt();});
    return gen_leptons;
}

// From EFTGenReader.h
// This is ill-defined for particles with multiple mothers
const reco::Candidate* EFTLHEReader::GetGenMotherNoFsr(const reco::Candidate* p) {
    if (p->numberOfMothers()) {
        const reco::Candidate* mom = p->mother(0);
        if (mom->pdgId() != p->pdgId()) {
            return mom;
        } else {
            return GetGenMotherNoFsr(mom);
        }
    } else {
        return p;
    }
}

std::vector<reco::GenJet> EFTLHEReader::GetGenJets(const std::vector<reco::GenJet>& inputs) {
    std::vector<reco::GenJet> ret;
    for (size_t i = 0; i < inputs.size(); i++) {
        const reco::GenJet& j = inputs.at(i);
        if (j.p4().Pt() < min_pt_jet) {
            continue;
        } else if (max_eta_jet > 0.0 && fabs(j.eta()) >= max_eta_jet) {
            continue;
        }
        ret.push_back(j);
    }
    std::sort(ret.begin(),ret.end(), [] (reco::GenJet a, reco::GenJet b) { return a.p4().Pt() > b.p4().Pt();});
    return ret;
}


#ifndef EFTGENREADER_GENREADER_EFTGENHISTSWITHCUTS_h
#define EFTGENREADER_GENREADER_EFTGENHISTSWITHCUTS_h

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
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/PythonParameterSet/interface/PythonProcessDesc.h"
#include "DataFormats/Common/interface/Handle.h"

// Physics
#include "Math/LorentzVector.h"
#include "DataFormats/Math/interface/deltaR.h"
#include "SimDataFormats/GeneratorProducts/interface/LHEEventProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

#include "EFTGenReader/EFTHelperUtilities/interface/WCPoint.h"
#include "EFTGenReader/EFTHelperUtilities/interface/WCFit.h"
#include "EFTGenReader/EFTHelperUtilities/interface/TH1EFT.h"
#include "EFTGenReader/EFTHelperUtilities/interface/utilities.h"

// end includes
// -----------------------------------------------

class EFTGenHistsWithCuts: public edm::EDAnalyzer
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
        void dumpParticles(const reco::GenParticleCollection& particles);
        void dumpJets(const std::vector<reco::GenJet>& jets);
        const reco::Candidate* GetGenMotherNoFsr(const reco::Candidate* p);
        std::pair<const reco::Candidate*, const reco::Candidate*> GetGenDaughterNoFsr(const reco::Candidate* p);
		
        reco::GenParticleCollection GetGenLeptons(const reco::GenParticleCollection& gen_particles);
        reco::GenParticleCollection GetGenParticlesSubset(const reco::GenParticleCollection& gen_particles, int pdg_id);
        std::vector<reco::GenJet> GetGenJets(const std::vector<reco::GenJet>& inputs);
        std::vector<reco::GenJet> CleanGenJets(const std::vector<reco::GenJet>& gen_jets, const reco::GenParticleCollection& gen_leptons);
        std::vector<reco::GenJet> MakePtEtaCuts(const std::vector<reco::GenJet>& inputs, std::string input_type);

        ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> getSumTLV(reco::GenParticleCollection col);
        double getdPhi(reco::GenParticle p1, reco::GenParticle p2);
        double getdR(reco::GenParticle p1, reco::GenParticle p2);
        double getInvMass(reco::GenParticle p1, reco::GenParticle p2);

    public:
        explicit EFTGenHistsWithCuts(const edm::ParameterSet&);
        ~EFTGenHistsWithCuts();

        template <typename T> edm::Handle<T> get_collection(const edm::Event& event, const edm::EDGetTokenT<T>& token);

        std::ofstream fout;
        FILE * ffout;

        std::string sampleName;
        int sampleNumber;
        int eventcount;

        // runtime configurable parameters
        bool iseft;
        bool debug;
        double min_pt_jet;
        double min_pt_lep;
        double max_eta_jet;
        double max_eta_lep;

        edm::EDGetTokenT<LHEEventProduct> lheInfo_token_;
        edm::EDGetTokenT<GenEventInfoProduct> genInfo_token_;
        edm::EDGetTokenT<reco::GenParticleCollection> genParticles_token_;  // reco::GenParticlesCollection is an alias for std::vector<reco::GenParticle>>
        edm::EDGetTokenT<std::vector<reco::GenJet> > genJets_token_;

        // Particel level stuff:
        edm::EDGetTokenT<std::vector<reco::GenJet>> particleLevelJetsToken_;
        edm::EDGetTokenT<std::vector<reco::GenJet>> particleLevelLeptonsToken_;

        // Histograms
        TH1D* h_SMwgt_norm;
        TH1EFT* h_eventsumEFT;

        // Jet Histograms
        TH1EFT* h_nJetsEFT; TH1D* h_nJetsSM;

        // njets histograms with some selecitons
        TH1EFT* h_nJets_3orMoreLep_EFT;           TH1D* h_nJets_3orMoreLep_SM;
        TH1EFT* h_nJets_3Lep_EFT;                 TH1D* h_nJets_3Lep_SM;
        TH1EFT* h_nJets_cleanJets_3orMoreLep_EFT; TH1D* h_nJets_cleanJets_3orMoreLep_SM;
        TH1EFT* h_nJets_cleanJets_3Lep_EFT;       TH1D* h_nJets_cleanJets_3Lep_SM;

        // Particle level jets hists
        TH1EFT* h_pl_nJets_EFT; TH1D* h_pl_nJets_SM;
        TH1EFT* h_pl_nJets_3Lep_EFT; TH1D* h_pl_nJets_3Lep_SM;

        // declare the tree
        TTree * summaryTree;

        // tree branches
        double originalXWGTUP_intree;
        int eventnum_intree;
        int lumiBlock_intree;
        int runNumber_intree;

        WCFit eft_wgt_intree;
        double sm_wgt_intree;
        int n_jet_intree;
        int n_bjet_intree;

        // Misc. counters
        int total_ls;
        int total_events;
        double total_orig_xsec;
        double total_sm_xsec;
};

void EFTGenHistsWithCuts::tree_add_branches()
{
    summaryTree->Branch("originalXWGTUP",&originalXWGTUP_intree);
    summaryTree->Branch("eventnum",&eventnum_intree);
    summaryTree->Branch("lumiBlock",&lumiBlock_intree);
    summaryTree->Branch("runNumber",&runNumber_intree);

    summaryTree->Branch("eft_wgt",&eft_wgt_intree);
    summaryTree->Branch("sm_wgt",&sm_wgt_intree);
    summaryTree->Branch("n_jets",&n_jet_intree);
    summaryTree->Branch("n_bjets",&n_bjet_intree);
}

void EFTGenHistsWithCuts::initialize_variables()
{
    originalXWGTUP_intree = -99;
    eventnum_intree = -1;
    lumiBlock_intree = -1;
    runNumber_intree = -1;
}

// Currently not setting any special parameters in config files
void EFTGenHistsWithCuts::parse_params() {}

void EFTGenHistsWithCuts::dumpParticles(const reco::GenParticleCollection& particles) {
    for (size_t i = 0; i < particles.size(); i++) {
        const reco::GenParticle* p = &(particles.at(i));
        int id = p->pdgId(); int st = p->status();
        double pt = p->p4().Pt(); double eta = p->p4().Eta();
        double mass = p->mass();

        uint n_dau = p->numberOfDaughters();
        uint n_mom = p->numberOfMothers();

        const reco::Candidate* p_mom_noFSR = GetGenMotherNoFsr(p);
        const reco::Candidate* p_gmom_noFSR = GetGenMotherNoFsr(p_mom_noFSR);
        int mom_id_noFSR = p_mom_noFSR->pdgId();
        int gmom_id_noFSR = p_gmom_noFSR->pdgId();

        bool is_hard_process = p->isHardProcess();

        std::cout << "pdgId: " << id << std::endl;
        std::cout << "\tpt: "    << pt
                  << "   eta: "  << eta
                  << "   mass: " << mass << std::endl;

        std::cout << "\tStatus: " << st << std::endl;
        std::cout << "\tisHard: " << is_hard_process << std::endl;
        std::cout << "\tRefSelf: " << p << " (pdgId: " << id << ")" << std::endl;
        std::cout << "\tRefMotherNoISR: " << p_mom_noFSR << " (pdgId: " << mom_id_noFSR << ")" << std::endl;
        std::cout << "\tRefGMotherNoISR: " << p_gmom_noFSR << " (pdgId: " << gmom_id_noFSR << ")" << std::endl;
        for (size_t j=0; j < n_mom; j++) {
            const reco::Candidate * mom = p->mother(j);
            int mom_id = mom->pdgId();
            std::cout << "\tRefMother:  " << mom << " (pdgId: " << mom_id << ")" << std::endl;
        }
        for (size_t j=0; j < n_dau; j++) {
            const reco::Candidate * d = p->daughter(j);
            int dau_id = d->pdgId();
            std::cout << "\tRefDaughter: " << d << " (pdgId: " << dau_id << ")" << std::endl;
        }
    }
}

void EFTGenHistsWithCuts::dumpJets(const std::vector<reco::GenJet>& jets) {
    for (size_t i = 0; i < jets.size(); i++) {
        const reco::GenJet* p = &(jets.at(i));
        double pt = p->p4().Pt(); double eta = p->p4().Eta();
        double mass = p->mass();

        std::cout << "Jet: " << i << std::endl;
        std::cout << "\tpt: "    << pt
                  << "   eta: "  << eta
                  << "   mass: " << mass << std::endl;
    }
}

// This is ill-defined for particles with multiple mothers
const reco::Candidate* EFTGenHistsWithCuts::GetGenMotherNoFsr(const reco::Candidate* p) {
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

std::pair<const reco::Candidate*, const reco::Candidate*> EFTGenHistsWithCuts::GetGenDaughterNoFsr(const reco::Candidate* p) {
    const reco::Candidate* d0 = p;
    const reco::Candidate* d1 = p;

    if (p->numberOfDaughters()) {
        const reco::Candidate* d = p->daughter(0);
        if (d->pdgId() != p->pdgId()) {
            d0 = d;
        } else {
            d0 = GetGenDaughterNoFsr(d).first;
        }
        if (p->numberOfDaughters() > 1) {
            d = p->daughter(1);
            if (d->pdgId() != p->pdgId()) {
                d1 = d;
            } else {
                d1 = GetGenDaughterNoFsr(d).second;
            }
        }
    }
    std::pair<const reco::Candidate*, const reco::Candidate*> ret_pair(d0,d1);
    return ret_pair;
}

reco::GenParticleCollection EFTGenHistsWithCuts::GetGenLeptons(const reco::GenParticleCollection& gen_particles) {
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
        bool is_fromHardGluon = (abs(mom_id) == 21 && p_mom->status() == 21);   // status == 21 corresponds to incoming hard particle
        bool is_fromHardQCD = (abs(mom_id) >= 1 && abs(mom_id) <= 5 && p_mom->status() == 21);

        if (is_fromEWBoson) {
            gen_leptons.push_back(p);
        } else if (is_hard_process && is_fromHardGluon) {
            gen_leptons.push_back(p);
        } else if (is_hard_process && is_fromHardQCD) {
            gen_leptons.push_back(p);
        }
    }
    std::sort(gen_leptons.begin(),gen_leptons.end(), [] (reco::GenParticle a, reco::GenParticle b) { return a.p4().Pt() > b.p4().Pt();});
    return gen_leptons;
}

reco::GenParticleCollection EFTGenHistsWithCuts::GetGenParticlesSubset(const reco::GenParticleCollection& gen_particles, int pdg_id) {
    reco::GenParticleCollection gen_subset;
    for (size_t i = 0; i < gen_particles.size(); i++) {
        const reco::GenParticle& p = gen_particles.at(i);
        int id = p.pdgId();
        if (abs(id) == pdg_id){
            gen_subset.push_back(p);
        }
    }
    return gen_subset;
}

// Clean jets based on how close they are to leptons
std::vector<reco::GenJet> EFTGenHistsWithCuts::CleanGenJets(const std::vector<reco::GenJet>& gen_jets, const reco::GenParticleCollection& gen_leptons) {
    std::vector<reco::GenJet> cleanedJets;
    bool isClean;
    double coneSize = 0.4;
    for (size_t i = 0; i < gen_jets.size(); i++){
        isClean = true;
        const reco::GenJet& jet = gen_jets.at(i);
        for (size_t j = 0; j < gen_leptons.size(); j++){
            const reco::GenParticle& lep = gen_leptons.at(j);
            //double dR = getdR(jet,lep);
            if (getdR(jet,lep) <= coneSize){
                isClean = false;
                //break;
            }
        }
        if (isClean) {
            cleanedJets.push_back(jet);
        }
    }
    return cleanedJets;
}

std::vector<reco::GenJet> EFTGenHistsWithCuts::GetGenJets(const std::vector<reco::GenJet>& inputs) {
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

// Pt and eta cuts particles of type std::vector<reco::GenJet>, note if the type is not specified as "jet" or "lep", generic cuts of 2.5 and 10 are used
std::vector<reco::GenJet> EFTGenHistsWithCuts::MakePtEtaCuts(const std::vector<reco::GenJet>& inputs, std::string input_type) {
    std::vector<reco::GenJet> ret;
    double max_eta = 2.5;
    double min_pt = 10;
    if (input_type=="jet"){
        max_eta = max_eta_jet;
        min_pt = min_pt_jet;
        std::cout << "\nUsing jet pt and eta cuts:" << min_pt << " " << max_eta << "\n" << std::endl;
    } else if (input_type=="lep"){
        max_eta = max_eta_lep;
        min_pt = min_pt_lep;
        std::cout << "\nUsing lep pt and eta cuts:" << min_pt << " " << max_eta << "\n" << std::endl;
    }
    for (size_t i = 0; i < inputs.size(); i++) {
        const reco::GenJet& p = inputs.at(i);
        if (p.p4().Pt() < min_pt) {
            continue;
        } else if (max_eta > 0.0 && fabs(p.eta()) >= max_eta) {
            continue;
        }
        ret.push_back(p);
    }
    std::sort(ret.begin(),ret.end(), [] (reco::GenJet a, reco::GenJet b) { return a.p4().Pt() > b.p4().Pt();});
    return ret;
}

ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> EFTGenHistsWithCuts::getSumTLV(reco::GenParticleCollection col) {
    ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> p4vec;
    for (const auto it: col) p4vec += it.p4();
    return p4vec;
}

double EFTGenHistsWithCuts::getdPhi(reco::GenParticle p1, reco::GenParticle p2) {
    double dPhi = p2.p4().Phi() - p1.p4().Phi();
    double pi = 2.0*asin(1.0);
    if (dPhi>=pi) dPhi = dPhi - 2.0*pi;     // from TVector2
    if (dPhi<(-pi)) dPhi = dPhi + 2.0*pi;   // from TVector2
    return dPhi;
}

double EFTGenHistsWithCuts::getdR(reco::GenParticle p1, reco::GenParticle p2) {
    double dR = (p1.p4().Eta() - p2.p4().Eta())*(p1.p4().Eta() - p2.p4().Eta());
    dR += (getdPhi(p1,p2)*getdPhi(p1,p2));
    dR = sqrt(dR);
    return dR;
}

double EFTGenHistsWithCuts::getInvMass(reco::GenParticle p1, reco::GenParticle p2) {
    auto p4vec = p1.p4() + p2.p4();
    return p4vec.M();
}

template<typename T>
edm::Handle<T>
EFTGenHistsWithCuts::get_collection(const edm::Event& event, const edm::EDGetTokenT<T>& token)
{
    edm::Handle<T> handle;
    event.getByToken(token, handle);
    if (!handle.isValid())
        throw edm::Exception(edm::errors::InvalidReference, "Can't find a collection.");
    return handle;
}

#endif

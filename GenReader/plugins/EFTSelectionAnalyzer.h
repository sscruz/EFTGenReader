// created by Andrew Wightman and Kelci Mohrman
#ifndef EFTGENREADER_GENREADER_SELECTIONANALYZER_h
#define EFTGENREADER_GENREADER_SELECTIONANALYZER_h

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
//#include "FWCore/PythonParameterSet/interface/PythonProcessDesc.h"
#include "FWCore/PythonParameterSet/interface/PyBind11ProcessDesc.h"
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
//#include "EFTGenReader/EFTHelperUtilities/interface/TH1EFT.h"

#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Jet.h"

// end includes
// -----------------------------------------------

class EFTSelectionAnalyzer: public edm::EDAnalyzer
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

        void addLepCategory(std::string category_name, std::string fit_tag);

        void dumpParticles(const reco::GenParticleCollection& particles);
        void dumpJets(const std::vector<reco::GenJet>& jets);
        const reco::Candidate* GetGenMotherNoFsr(const reco::Candidate* p);
        std::pair<const reco::Candidate*, const reco::Candidate*> GetGenDaughterNoFsr(const reco::Candidate* p);

        std::vector<reco::GenJet>   MakeJetPtEtaCuts(const std::vector<reco::GenJet>& gen_jets, double pt_cut, double eta_cut);
        reco::GenParticleCollection MakePtEtaCuts(const reco::GenParticleCollection& gen_particles, double pt_cut, double eta_cut);
        reco::GenParticleCollection GetGenParticlesSubset(const reco::GenParticleCollection& gen_particles, int pdg_id);
        reco::GenParticleCollection GetGenElectrons(const reco::GenParticleCollection& gen_particles);
        reco::GenParticleCollection GetGenMuons(const reco::GenParticleCollection& gen_particles);
        reco::GenParticleCollection GetGenLeptons(const reco::GenParticleCollection& gen_particles);
        std::vector<reco::GenJet> GetGenJets(const std::vector<reco::GenJet>& inputs);
        std::vector<pat::Electron> MakePtEtaCutsPatElectrons(const std::vector<pat::Electron>& pat_electrons, double pt_cut, double eta_cut);
        std::vector<pat::Muon> MakePtEtaCutsPatMuons(const std::vector<pat::Muon>& pat_muon, double pt_cut, double eta_cut);
        std::vector<pat::Jet> MakePatJetPtEtaCuts(const std::vector<pat::Jet>& pat_jet, double pt_cut, double eta_cut);
        ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> getSumTLV(reco::GenParticleCollection col);
        double getdPhi(reco::GenParticle p1, reco::GenParticle p2);
        double getdR(reco::GenParticle p1, reco::GenParticle p2);
        double getInvMass(reco::GenParticle p1, reco::GenParticle p2);

    public:
        explicit EFTSelectionAnalyzer(const edm::ParameterSet&);
        ~EFTSelectionAnalyzer();

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

        edm::EDGetTokenT<std::vector<pat::Electron> > patElectrons_token_;
        edm::EDGetTokenT<std::vector<pat::Muon> > patMuons_token_;
        edm::EDGetTokenT<std::vector<pat::Jet> > patJets_token_;

        std::unordered_map<std::string,WCFit*> lep_fits_dict;
        std::vector<std::string> lep_category_names_vect; // keys in the dict (but do not have to be names of the fit tags)

        // Misc. counters
        int total_ls;
        int total_events;
        double total_orig_xsec;
        double total_sm_xsec;
};


void EFTSelectionAnalyzer::addLepCategory(std::string category_name, std::string fit_tag) {
    if (lep_fits_dict.count(category_name) > 0 ) {
        std::cout << "Warning: The category \"" << category_name << "\" already exists in dictionary, not adding." << std::endl;
        return;
    }
    WCFit* new_category_fit = new WCFit();
    new_category_fit->setTag(fit_tag);
    this->lep_fits_dict[category_name] = new_category_fit;
    this->lep_category_names_vect.push_back(category_name);
    return;
}

void EFTSelectionAnalyzer::dumpParticles(const reco::GenParticleCollection& particles) {
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

void EFTSelectionAnalyzer::dumpJets(const std::vector<reco::GenJet>& jets) {
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
const reco::Candidate* EFTSelectionAnalyzer::GetGenMotherNoFsr(const reco::Candidate* p) {
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

std::pair<const reco::Candidate*, const reco::Candidate*> EFTSelectionAnalyzer::GetGenDaughterNoFsr(const reco::Candidate* p) {
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

reco::GenParticleCollection EFTSelectionAnalyzer::GetGenParticlesSubset(const reco::GenParticleCollection& gen_particles, int pdg_id) {
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

std::vector<reco::GenJet> EFTSelectionAnalyzer::MakeJetPtEtaCuts(const std::vector<reco::GenJet>& gen_jets, double pt_cut, double eta_cut) {
    std::vector<reco::GenJet> gen_jets_cut;
    for (size_t i = 0; i < gen_jets.size(); i++) {
        const reco::GenJet& j = gen_jets.at(i);
        if (j.p4().Pt() < pt_cut) { // Do not include particles whose pt is less than the pt cut
            continue;
        } else if (eta_cut > 0.0 && fabs(j.eta()) >= eta_cut) { // Do not include particles whose eta is greater than the eta cut
            continue;
        } else {
            gen_jets_cut.push_back(j);
        }
    }
    return gen_jets_cut;
}
reco::GenParticleCollection EFTSelectionAnalyzer::MakePtEtaCuts(const reco::GenParticleCollection& gen_particles, double pt_cut, double eta_cut) {
    reco::GenParticleCollection gen_particles_cut;
    for (size_t i = 0; i < gen_particles.size(); i++) {
        const reco::GenParticle& p = gen_particles.at(i);
        if (p.p4().Pt() < pt_cut) { // Do not include particles whose pt is less than the pt cut
            continue;
        } else if (eta_cut > 0.0 && fabs(p.eta()) >= eta_cut) { // Do not include particles whose eta is greater than the eta cut
            continue;
        } else {
            gen_particles_cut.push_back(p);
        }
    }
    return gen_particles_cut;
}

// Reco level pt and eta cuts: 
std::vector<pat::Electron> EFTSelectionAnalyzer::MakePtEtaCutsPatElectrons(const std::vector<pat::Electron>& pat_electrons, double pt_cut, double eta_cut) {
    std::vector<pat::Electron> pat_electrons_cut;
    for (size_t i = 0; i < pat_electrons.size(); i++) {
        const pat::Electron& p = pat_electrons.at(i);
        if (p.p4().Pt() < pt_cut) { // Do not include particles whose pt is less than the pt cut
            continue;
        } else if (eta_cut > 0.0 && fabs(p.eta()) >= eta_cut) { // Do not include particles whose eta is greater than the eta cut
            continue;
        } else {
            pat_electrons_cut.push_back(p);
        }
    }
    return pat_electrons_cut;
}
std::vector<pat::Muon> EFTSelectionAnalyzer::MakePtEtaCutsPatMuons(const std::vector<pat::Muon>& pat_muons, double pt_cut, double eta_cut) {
    std::vector<pat::Muon> pat_muons_cut;
    for (size_t i = 0; i < pat_muons.size(); i++) {
        const pat::Muon& p = pat_muons.at(i);
        if (p.p4().Pt() < pt_cut) { // Do not include particles whose pt is less than the pt cut
            continue;
        } else if (eta_cut > 0.0 && fabs(p.eta()) >= eta_cut) { // Do not include particles whose eta is greater than the eta cut
            continue;
        } else {
            pat_muons_cut.push_back(p);
        }
    }
    return pat_muons_cut;
}
std::vector<pat::Jet> EFTSelectionAnalyzer::MakePatJetPtEtaCuts(const std::vector<pat::Jet>& pat_jets, double pt_cut, double eta_cut) {
    std::vector<pat::Jet> pat_jets_cut;
    for (size_t i = 0; i < pat_jets.size(); i++) {
        const pat::Jet& j = pat_jets.at(i);
        if (j.p4().Pt() < pt_cut) { // Do not include particles whose pt is less than the pt cut
            continue;
        } else if (eta_cut > 0.0 && fabs(j.eta()) >= eta_cut) { // Do not include particles whose eta is greater than the eta cut
            continue;
        } else {
            pat_jets_cut.push_back(j);
        }
    }
    return pat_jets_cut;
}

reco::GenParticleCollection EFTSelectionAnalyzer::GetGenLeptons(const reco::GenParticleCollection& gen_particles) {
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
        //is_neutrino = (abs(id) == 12 || abs(id) == 14 || abs(id) == 16);

        //if (!is_lepton && !is_neutrino) {
            //continue;
        //}
        if (!is_lepton) {
            continue;
        }

        //// Note: The kinematic cuts could realistically go anywhere in this loop
        //if (p.p4().Pt() < min_pt_lep) {
        //    continue;
        //} else if (max_eta_lep > 0.0 && fabs(p.eta()) >= max_eta_lep) {
        //    continue;
        //}

        int mom_id = id;    // If no mother, set to own id
        if (p_mom) mom_id = p_mom->pdgId();

        bool is_fromEWBoson = (abs(mom_id) >= 22 && abs(mom_id) <= 25);
        bool is_fromTopSystem = (abs(mom_id) == 24 && abs(gmom_noFSR_id) == 6);

        //if (is_fromTopSystem) {
            //continue;
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

std::vector<reco::GenJet> EFTSelectionAnalyzer::GetGenJets(const std::vector<reco::GenJet>& inputs) {
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

ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> EFTSelectionAnalyzer::getSumTLV(reco::GenParticleCollection col) {
    ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> p4vec;
    for (const auto it: col) p4vec += it.p4();
    return p4vec;
}

double EFTSelectionAnalyzer::getdPhi(reco::GenParticle p1, reco::GenParticle p2) {
    double dPhi = p2.p4().Phi() - p1.p4().Phi();
    double pi = 2.0*asin(1.0);
    if (dPhi>=pi) dPhi = dPhi - 2.0*pi;     // from TVector2
    if (dPhi<(-pi)) dPhi = dPhi + 2.0*pi;   // from TVector2
    return dPhi;
}

double EFTSelectionAnalyzer::getdR(reco::GenParticle p1, reco::GenParticle p2) {
    double dR = (p1.p4().Eta() - p2.p4().Eta())*(p1.p4().Eta() - p2.p4().Eta());
    dR += (getdPhi(p1,p2)*getdPhi(p1,p2));
    dR = sqrt(dR);
    return dR;
}

double EFTSelectionAnalyzer::getInvMass(reco::GenParticle p1, reco::GenParticle p2) {
    auto p4vec = p1.p4() + p2.p4();
    return p4vec.M();
}

template<typename T>
edm::Handle<T>
EFTSelectionAnalyzer::get_collection(const edm::Event& event, const edm::EDGetTokenT<T>& token)
{
    edm::Handle<T> handle;
    event.getByToken(token, handle);
    if (!handle.isValid())
        throw edm::Exception(edm::errors::InvalidReference, "Can't find a collection.");
    return handle;
}

#endif

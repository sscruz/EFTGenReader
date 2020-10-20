// created by Andrew Wightman
#ifndef EFTGENREADER_GENREADER_EFTGENREADER_h
#define EFTGENREADER_GENREADER_EFTGENREADER_h

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
//#include "​DataFormats/​JetReco/​interface/​GenJet.h"
#include "DataFormats/PatCandidates/interface/Jet.h"
#include "DataFormats/PatCandidates/interface/PackedGenParticle.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"

// Multilepton
//#include "ttH-13TeVMultiLeptons/TemplateMakers/interface/JobParameters.h"
//#include "ttH-13TeVMultiLeptons/TemplateMakers/src/classes.h"
//#include "ttH-13TeVMultiLeptons/TemplateMakers/interface/objectClasses.h"

//#include "ttH-13TeVMultiLeptons/TemplateMakers/interface/WCPoint.h"
//#include "EFTGenReader/GenReader/interface/WCPoint.h"
//#include "EFTGenReader/GenReader/interface/WCFit.h"
//#include "EFTGenReader/GenReader/interface/TH1EFT.h"
#include "EFTGenReader/EFTHelperUtilities/interface/WCPoint.h"
#include "EFTGenReader/EFTHelperUtilities/interface/WCFit.h"
#include "EFTGenReader/EFTHelperUtilities/interface/TH1EFT.h"
#include "EFTGenReader/EFTHelperUtilities/interface/utilities.h"

// end includes
// -----------------------------------------------

class EFTGenReader: public edm::EDAnalyzer
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
        ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> getSumTLV(reco::GenParticleCollection col);
        double getdPhi(reco::GenParticle p1, reco::GenParticle p2);
        double getdR(reco::GenParticle p1, reco::GenParticle p2);
        double getInvMass(reco::GenParticle p1, reco::GenParticle p2);

    public:
        explicit EFTGenReader(const edm::ParameterSet&);
        ~EFTGenReader();

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

        // Histograms
        TH1D* h_SMwgt_norm;

        TH1EFT* h_eventsumEFT;

        // pdgId Histograms
        TH1EFT* h_pdgIdEFT;            TH1D* h_pdgIdSM;
        TH1EFT* h_pdgIdFromZEFT;       TH1D* h_pdgIdFromZSM;
        TH1EFT* h_pdgIdFromHEFT;       TH1D* h_pdgIdFromHSM;
        TH1EFT* h_pdgIdFromWEFT;       TH1D* h_pdgIdFromWSM;
        TH1EFT* h_pdgIdFromPhotonEFT;  TH1D* h_pdgIdFromPhotonSM;
        TH1EFT* h_pdgIdFromGluonEFT;   TH1D* h_pdgIdFromGluonSM;
        TH1EFT* h_pdgIdFromQCDEFT;     TH1D* h_pdgIdFromQCDSM;
        TH1EFT* h_pdgIdFromOtherEFT;   TH1D* h_pdgIdFromOtherSM;
        TH1EFT* h_pdgIdLepMotherEFT;   TH1D* h_pdgIdLepMotherSM;
        TH1EFT* h_pdgIdLepGrMotherEFT; TH1D* h_pdgIdLepGrMotherSM;

        TH1EFT* h_pdgIdElectronMotherEFT;   TH1D* h_pdgIdElectronMotherSM;
        TH1EFT* h_pdgIdElectronGrMotherEFT; TH1D* h_pdgIdElectronGrMotherSM;
        TH1EFT* h_pdgIdMuonMotherEFT;       TH1D* h_pdgIdMuonMotherSM;
        TH1EFT* h_pdgIdMuonGrMotherEFT;     TH1D* h_pdgIdMuonGrMotherSM;
        TH1EFT* h_pdgIdTauMotherEFT;        TH1D* h_pdgIdTauMotherSM;
        TH1EFT* h_pdgIdTauGrMotherEFT;      TH1D* h_pdgIdTauGrMotherSM;

        
        // Lepton Histograms
        TH1EFT* h_prompt_leptonsEFT;   TH1D* h_prompt_leptonsSM;
        TH1EFT* h_lep_ptEFT;           TH1D* h_lep_ptSM;
        TH1EFT* h_lep_etaEFT;          TH1D* h_lep_etaSM;
        TH1EFT* h_lep1_ptEFT;          TH1D* h_lep1_ptSM;
        TH1EFT* h_lep2_ptEFT;          TH1D* h_lep2_ptSM;
        TH1EFT* h_mllEFT;              TH1D* h_mllSM;
        TH1EFT* h_deltaREFT;           TH1D* h_deltaRSM;
        TH1EFT* h_lepSum_ptEFT;        TH1D* h_lepSum_ptSM;

        // Lepton Histograms: electrons
        TH1EFT* h_prompt_electronsEFT; TH1D* h_prompt_electronsSM;
        TH1EFT* h_e_ptEFT;             TH1D* h_e_ptSM;
        TH1EFT* h_e_etaEFT;            TH1D* h_e_etaSM;
        TH1EFT* h_e1_ptEFT;            TH1D* h_e1_ptSM;
        TH1EFT* h_e2_ptEFT;            TH1D* h_e2_ptSM;
        TH1EFT* h_meeEFT;              TH1D* h_meeSM;
        //TH1EFT* h_e_deltaREFT;         TH1D* h_e_deltaRSM;
        TH1EFT* h_eSum_ptEFT;          TH1D* h_eSum_ptSM;

        // Lepton Histograms: muons
        TH1EFT* h_prompt_muonsEFT;     TH1D* h_prompt_muonsSM;
        TH1EFT* h_mu_ptEFT;            TH1D* h_mu_ptSM;
        TH1EFT* h_mu_etaEFT;           TH1D* h_mu_etaSM;
        TH1EFT* h_mu1_ptEFT;           TH1D* h_mu1_ptSM;
        TH1EFT* h_mu2_ptEFT;           TH1D* h_mu2_ptSM;
        TH1EFT* h_mmumuEFT;            TH1D* h_mmumuSM;
        //TH1EFT* h_mu_deltaREFT;        TH1D* h_mu_deltaRSM;
        TH1EFT* h_muSum_ptEFT;         TH1D* h_muSum_ptSM;

        // Lepton Histograms: Taus
        TH1EFT* h_prompt_tausEFT;      TH1D* h_prompt_tausSM;
        TH1EFT* h_tau_ptEFT;           TH1D* h_tau_ptSM;
        TH1EFT* h_tau_etaEFT;          TH1D* h_tau_etaSM;
        TH1EFT* h_tau1_ptEFT;          TH1D* h_tau1_ptSM;
        TH1EFT* h_tau2_ptEFT;          TH1D* h_tau2_ptSM;
        TH1EFT* h_mtautauEFT;          TH1D* h_mtautauSM;
        //TH1EFT* h_tau_deltaREFT;       TH1D* h_tau_deltaRSM;
        TH1EFT* h_tauSum_ptEFT;        TH1D* h_tauSum_ptSM;

        // Jet Histograms
        TH1EFT* h_nJetsEFT;   TH1D* h_nJetsSM;
        TH1EFT* h_jet_ptEFT;  TH1D* h_jet_ptSM;
        TH1EFT* h_jet1_ptEFT; TH1D* h_jet1_ptSM;
        TH1EFT* h_jet2_ptEFT; TH1D* h_jet2_ptSM;
        TH1EFT* h_jet3_ptEFT; TH1D* h_jet3_ptSM;
        TH1EFT* h_jet4_ptEFT; TH1D* h_jet4_ptSM;

        TH1EFT* h_jet_etaEFT;  TH1D* h_jet_etaSM;
        TH1EFT* h_jet1_etaEFT; TH1D* h_jet1_etaSM;
        TH1EFT* h_jet2_etaEFT; TH1D* h_jet2_etaSM;
        TH1EFT* h_jet3_etaEFT; TH1D* h_jet3_etaSM;
        TH1EFT* h_jet4_etaEFT; TH1D* h_jet4_etaSM;

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

void EFTGenReader::tree_add_branches()
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

void EFTGenReader::initialize_variables()
{
    originalXWGTUP_intree = -99;
    eventnum_intree = -1;
    lumiBlock_intree = -1;
    runNumber_intree = -1;
}

// Currently not setting any special parameters in config files
void EFTGenReader::parse_params() {}

void EFTGenReader::dumpParticles(const reco::GenParticleCollection& particles) {
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

void EFTGenReader::dumpJets(const std::vector<reco::GenJet>& jets) {
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
const reco::Candidate* EFTGenReader::GetGenMotherNoFsr(const reco::Candidate* p) {
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

std::pair<const reco::Candidate*, const reco::Candidate*> EFTGenReader::GetGenDaughterNoFsr(const reco::Candidate* p) {
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

reco::GenParticleCollection EFTGenReader::GetGenLeptons(const reco::GenParticleCollection& gen_particles) {
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

        if (is_fromTopSystem) {
            continue;
        }

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

reco::GenParticleCollection EFTGenReader::GetGenParticlesSubset(const reco::GenParticleCollection& gen_particles, int pdg_id) {
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

std::vector<reco::GenJet> EFTGenReader::GetGenJets(const std::vector<reco::GenJet>& inputs) {
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

ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> EFTGenReader::getSumTLV(reco::GenParticleCollection col) {
    ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> p4vec;
    for (const auto it: col) p4vec += it.p4();
    return p4vec;
}

double EFTGenReader::getdPhi(reco::GenParticle p1, reco::GenParticle p2) {
    double dPhi = p2.p4().Phi() - p1.p4().Phi();
    double pi = 2.0*asin(1.0);
    if (dPhi>=pi) dPhi = dPhi - 2.0*pi;     // from TVector2
    if (dPhi<(-pi)) dPhi = dPhi + 2.0*pi;   // from TVector2
    return dPhi;
}

double EFTGenReader::getdR(reco::GenParticle p1, reco::GenParticle p2) {
    double dR = (p1.p4().Eta() - p2.p4().Eta())*(p1.p4().Eta() - p2.p4().Eta());
    dR += (getdPhi(p1,p2)*getdPhi(p1,p2));
    dR = sqrt(dR);
    return dR;
}

double EFTGenReader::getInvMass(reco::GenParticle p1, reco::GenParticle p2) {
    auto p4vec = p1.p4() + p2.p4();
    return p4vec.M();
}

template<typename T>
edm::Handle<T>
EFTGenReader::get_collection(const edm::Event& event, const edm::EDGetTokenT<T>& token)
{
    edm::Handle<T> handle;
    event.getByToken(token, handle);
    if (!handle.isValid())
        throw edm::Exception(edm::errors::InvalidReference, "Can't find a collection.");
    return handle;
}

#endif

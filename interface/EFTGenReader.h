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
#include "FWCore/PythonParameterSet/interface/PythonProcessDesc.h"

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

#include "ttH-13TeVMultiLeptons/TemplateMakers/interface/WCPoint.h"

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

        // cfi params
        bool iseft;
        bool debug;
        int gp_events;
        int norm_type;
        double xsec_norm;
        double intg_lumi;

        std::vector<TH1D*> th1d_hists;

        TH1D* h_SMwgt_norm;

        // pdgId Histograms
        TH1D* h_pdgId;            TH1D* h_pdgIdSM;
        TH1D* h_pdgIdFromZ;       TH1D* h_pdgIdFromZSM;
        TH1D* h_pdgIdFromH;       TH1D* h_pdgIdFromHSM;
        TH1D* h_pdgIdFromW;       TH1D* h_pdgIdFromWSM;
        TH1D* h_pdgIdFromPhoton;  TH1D* h_pdgIdFromPhotonSM;
        TH1D* h_pdgIdFromGluon;   TH1D* h_pdgIdFromGluonSM;
        TH1D* h_pdgIdFromQCD;     TH1D* h_pdgIdFromQCDSM;
        TH1D* h_pdgIdFromOther;   TH1D* h_pdgIdFromOtherSM;
        TH1D* h_pdgIdLepMother;   TH1D* h_pdgIdLepMotherSM;
        TH1D* h_pdgIdLepGrMother; TH1D* h_pdgIdLepGrMotherSM;
        
        // Lepton Histograms
        TH1D* h_prompt_leptons; TH1D* h_prompt_leptonsSM;
        TH1D* h_lep_pt;         TH1D* h_lep_ptSM;
        TH1D* h_lep_eta;        TH1D* h_lep_etaSM;
        TH1D* h_lep1_pt;        TH1D* h_lep1_ptSM;
        TH1D* h_lep2_pt;        TH1D* h_lep2_ptSM;
        TH1D* h_mll;            TH1D* h_mllSM;
        TH1D* h_deltaR;         TH1D* h_deltaRSM;
        TH1D* h_lepSum_pt;      TH1D* h_lepSum_ptSM;

        // Jet Histograms
        TH1D* h_nJets;   TH1D* h_nJetsSM;
        TH1D* h_jet_pt;  TH1D* h_jet_ptSM;
        TH1D* h_jet1_pt; TH1D* h_jet1_ptSM;
        TH1D* h_jet2_pt; TH1D* h_jet2_ptSM;
        TH1D* h_jet_eta; TH1D* h_jet_etaSM;

        edm::ParameterSet entire_pset;

        // Note: reco::GenParticlesCollection is an alias for std::vector<reco::GenParticle>>
        edm::EDGetTokenT<LHEEventProduct> lheInfo_token_;
        edm::EDGetTokenT<GenEventInfoProduct> genInfo_token_;
        edm::EDGetTokenT<reco::GenParticleCollection> genParticles_token_;
        //edm::EDGetTokenT<std::vector<pat::PackedGenParticleCollection> > genPackedParticles_token_;
        edm::EDGetTokenT<std::vector<reco::GenJet> > slimmed_genJets_token_;
        edm::EDGetTokenT<std::vector<reco::GenJet> > slimmed_genJetsAK8_token_;
        edm::EDGetTokenT<std::vector<reco::GenJet> > slimmed_genJetsAK8SoftDropSubJets_token_;

        // declare the tree
        TTree * summaryTree;

        // tree branches
        //std::unordered_map<std::string,double> eftwgts_intree;
        double originalXWGTUP_intree;
        int eventnum_intree;
        int lumiBlock_intree;
        int runNumber_intree;

        int total_ls;
        int total_events;
        double total_orig_xsec;
        double total_sm_xsec;

        //std::vector<ttH::GenParticle> pruned_genParticles_intree;
        reco::GenParticleCollection pruned_genParticles_intree;
        std::vector<reco::GenJet> slimmed_genJets_intree;
        std::vector<reco::GenJet> slimmed_genJetsAK8_intree;
        std::vector<reco::GenJet> slimmed_genJetsAK8SoftDropSubJets_intree;
};

void EFTGenReader::tree_add_branches()
{
    //summaryTree->Branch("eftwgts",&eftwgts_intree);
    summaryTree->Branch("originalXWGTUP",&originalXWGTUP_intree);

    summaryTree->Branch("eventnum",&eventnum_intree);
    summaryTree->Branch("lumiBlock",&lumiBlock_intree);
    summaryTree->Branch("runNumber",&runNumber_intree);
    summaryTree->Branch("pruned_genParticles",&pruned_genParticles_intree);
    summaryTree->Branch("genJets",&slimmed_genJets_intree);
    summaryTree->Branch("genJets_AK8",&slimmed_genJetsAK8_intree);
    summaryTree->Branch("genJets_AK8SoftDrop",&slimmed_genJetsAK8SoftDropSubJets_intree);
}

void EFTGenReader::initialize_variables()
{
    //eftwgts_intree.clear();
    originalXWGTUP_intree = -99;
    eventnum_intree = -1;
    lumiBlock_intree = -1;
    runNumber_intree = -1;

    pruned_genParticles_intree.clear();
    slimmed_genJets_intree.clear();
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

std::vector<reco::GenJet> EFTGenReader::GetGenJets(const std::vector<reco::GenJet>& inputs) {
    std::vector<reco::GenJet> ret;
    for (size_t i = 0; i < inputs.size(); i++) {
        const reco::GenJet& j = inputs.at(i);

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
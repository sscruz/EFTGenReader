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
#include "EFTGenReader/EFTHelperUtilities/interface/TH2EFT.h"
#include "EFTGenReader/EFTHelperUtilities/interface/utilities.h"

#define IS_BHADRON_PDGID(id) ( ((abs(id)/100)%10 == 5) || (abs(id) >= 5000 && abs(id) <= 5999) )

// end includes
// -----------------------------------------------

// Structure for storing info for all of the hists we want to plot
struct HistInfo {
    TString h_multiplicity;
    TString h_type;
    int h_bins;
    int h_min;
    int h_max;
    int h_no;
};

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
        std::vector<reco::GenJet> GetGenBJets(const std::vector<reco::GenJet>& inputs);


        ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> getSumTLV(reco::GenParticleCollection col);
        double getdPhi(reco::GenParticle p1, reco::GenParticle p2);
        //double getdR(reco::GenParticle p1, reco::GenParticle p2);
        double getInvMass(reco::GenParticle p1, reco::GenParticle p2);

        //TString constructHistName(TString lep_cat, TString hist_type, int hist_no);
        TString constructHistName(TString lep_cat, TString hist_type, std::vector<size_t> hist_idx_vect);
        void fillHistIfExists(TString h_name, double val, WCFit eft_fit);


    public:
        explicit EFTGenHistsWithCuts(const edm::ParameterSet&);
        ~EFTGenHistsWithCuts();

        template <typename T> edm::Handle<T> get_collection(const edm::Event& event, const edm::EDGetTokenT<T>& token);
        template <typename T> std::vector<T> MakePtEtaCuts(const std::vector<T>& input, double min_pt, double max_eta);
        template <typename T1, typename T2> std::vector<T1> CleanGenJets(const std::vector<T1>& obj1_vector, const std::vector<T2>& obj2_vector, const double coneSize);
        template <typename T1, typename T2> double getdR(T1 p1, T2 p2);

        template <typename T> TString getLepCat(const std::vector<T>& leptons);
        template <typename T> int getChargeSum(const std::vector<T>& particles_vect);
        template <typename T> std::vector<T> getChargedParticles(const std::vector<T>& particles_vect);

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
        //TH1EFT* h_nJetsEFT; TH1D* h_nJetsSM;

        // njets histograms with some selecitons
        //TH1EFT* h_nJets_3orMoreLep_EFT;           TH1D* h_nJets_3orMoreLep_SM;
        //TH1EFT* h_nJets_3Lep_EFT;                 TH1D* h_nJets_3Lep_SM;
        //TH1EFT* h_nJets_cleanJets_3orMoreLep_EFT; TH1D* h_nJets_cleanJets_3orMoreLep_SM;
        //TH1EFT* h_nJets_cleanJets_3Lep_EFT;       TH1D* h_nJets_cleanJets_3Lep_SM;

        // Particle level jets hists
        //TH1EFT* h_pl_nJets_EFT; TH1D* h_pl_nJets_SM;
        TH1EFT* h_pl_nJets_3Lep_EFT; TH1D* h_pl_nJets_3Lep_SM;
        TH1EFT* h_pl_clean_nJets_3Lep_EFT;
        TH2EFT* h_bjet_jetEFT; TH2D* h_bjet_jetSM;
        TH2EFT* h_bjet_jet2lEFT; TH2D* h_bjet_jet2lSM;
        TH2EFT* h_bjet_jet3lEFT; TH2D* h_bjet_jet3lSM;
        TH2EFT* h_bjet_jet4lEFT; TH2D* h_bjet_jet4lSM;

	TH2EFT* h_2lss_jetbjetEFT; TH2D* h_2lss_jetbjetSM;
	TH2EFT* h_3l_jetbjetEFT;   TH2D* h_3l_jetbjetSM;
	TH2EFT* h_4l_jetbjetEFT;   TH2D* h_4l_jetbjetSM;

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

        // Set up the categories we want to consider and histogram types we want to make

        std::map<TString,TH1EFT*> hist_dict;
        std::vector<std::string> lep_cats_vect {"2lss","3l","4l"};

        int n_eta_bins = 12;
        int eta_min = -3;
        int eta_max = 3;

        int n_pt_bins = 12;
        int pt_min = 0;
        int pt_max = 500;
        int ht_max = 1000;

        int n_njet_bins = 14;
        int njet_min = 0;
        int njet_max = n_njet_bins;

        int n_dr_bins = 10;

        std::vector<HistInfo> hist_info_vec {
            //{"multiplicity type" , "name of hist type" , number of bins, min, max, number of hists of this type}
            {"single" , "jet_pt"  , n_pt_bins   , pt_min   , pt_max   , 4},
            {"single" , "jet_eta" , n_eta_bins  , eta_min  , eta_max  , 4},
            {"single" , "lep_pt"  , n_pt_bins   , pt_min   , pt_max   , 4},
            {"single" , "lep_eta" , n_eta_bins  , eta_min  , eta_max  , 4},
            {"all"    , "njets"   , n_njet_bins , njet_min , njet_max , 1},
            {"all"    , "ht"      , n_pt_bins   , pt_min   , ht_max   , 1},
            {"pair"   , "jet_dR"  , n_dr_bins   , 0 , 5 , 4},
            {"pair"   , "lep_dR"  , n_dr_bins   , 0 , 5 , 4},
            {"pair"   , "lep_mll"  , 30 , 0 , 500 , 4}

        };

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

        /*
        // Note: The kinematic cuts could realistically go anywhere in this loop
        if (p.p4().Pt() < min_pt_lep) {
            continue;
        } else if (max_eta_lep > 0.0 && fabs(p.eta()) >= max_eta_lep) {
            continue;
        }
        */

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

std::vector<reco::GenJet> EFTGenHistsWithCuts::GetGenBJets(const std::vector<reco::GenJet>& inputs) {
    std::vector<reco::GenJet> gen_bjets; //stores b-tagged jets
    for(auto &genJet : inputs) {
        std::vector< const reco::Candidate * > jconst=genJet.getJetConstituentsQuick(); //Loop over all particles in the jet
        for(size_t ijc=0; ijc <jconst.size(); ijc++) {
            const reco::Candidate *par=jconst[ijc];
            if(par->status() != 2) continue; //Status 2 is decayed or fragmented entry (i.e. decayed particle or parton produced in shower.) https://twiki.cern.ch/twiki/bin/view/CMSPublic/WorkBookGenParticleCandidate
            int absid = abs(par->pdgId()); //Store |pdgId|
            if(IS_BHADRON_PDGID(absid)) {
                gen_bjets.push_back(genJet); //Push jet to the gen_bjets vector
                break; //We found a B meson, meaning this is a b-jet! No need to search further
            }
        }
    }
    return gen_bjets;
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

double EFTGenHistsWithCuts::getInvMass(reco::GenParticle p1, reco::GenParticle p2) {
    auto p4vec = p1.p4() + p2.p4();
    return p4vec.M();
}

// Make a standardized hist name out of category info and hist type and hist number
/*
TString EFTGenHistsWithCuts::constructHistName(TString lep_cat, TString hist_type, int hist_no){
    TString ret_str;
    TString hist_no_str;
    if (hist_no == -1){
        hist_no_str = "";
    } else {
        hist_no_str = std::to_string(hist_no);
    }
    ret_str = "h_"+lep_cat+"_"+hist_type+"_"+hist_no_str;
    return ret_str;
}
*/
TString EFTGenHistsWithCuts::constructHistName(TString lep_cat, TString hist_type, std::vector<size_t> hist_idx_vect){
    TString ret_str;
    TString hist_no_str = "";
    int n_indices = hist_idx_vect.size();
    for (int i=0; i<n_indices; i++){
        if (i==0){ // If first index, add the "_" before first number
            hist_no_str = hist_no_str + "_";
        }
        hist_no_str = hist_no_str + std::to_string(hist_idx_vect.at(i));
        if (i<n_indices-1){ // If not last index, ad "-" after number
            hist_no_str = hist_no_str + "-";
        }
    }
    ret_str = "h_"+lep_cat+"_"+hist_type+hist_no_str;
    return ret_str;
}

// Fill TH1EFTs if they exist in the hist dictionary
void EFTGenHistsWithCuts::fillHistIfExists(TString h_name, double val, WCFit eft_fit){
    if (hist_dict.find(h_name) != hist_dict.end()){
        hist_dict[h_name]->Fill(val,1.0,eft_fit);
        if (h_name == "h_3l_njets_1"){
        }
    }
}

// Template function for making pt and eta cuts
template <typename T>
std::vector<T> EFTGenHistsWithCuts::MakePtEtaCuts(const std::vector<T>& input, double min_pt, double max_eta) {
    std::vector<T> ret;
    for (size_t i = 0; i < input.size(); i++) {
        const T& p = input.at(i);
        if (p.p4().Pt() < min_pt) {
            continue;
        } else if (max_eta > 0.0 && fabs(p.eta()) >= max_eta) {
            continue;
        }
        ret.push_back(p);
    }
    std::sort(ret.begin(),ret.end(), [] (T a, T b) { return a.p4().Pt() > b.p4().Pt();});
    return ret;
}

// Clean objects based on how close they are to other objects
template <typename T1, typename T2> 
std::vector<T1> EFTGenHistsWithCuts::CleanGenJets(const std::vector<T1>& obj1_vector, const std::vector<T2>& obj2_vector, const double coneSize) {
    std::vector<T1> cleaned_obj1_vector;
    bool isClean;
    for (size_t i = 0; i < obj1_vector.size(); i++){
        isClean = true;
        const T1& obj1 = obj1_vector.at(i);
        for (size_t j = 0; j < obj2_vector.size(); j++){
            const T2& obj2 = obj2_vector.at(j);
            if (getdR(obj1,obj2) <= coneSize){
                isClean = false;
            }
        }
        if (isClean) {
            cleaned_obj1_vector.push_back(obj1);
        }
    }
    return cleaned_obj1_vector;
}

// Get dR between two objects
template <typename T1, typename T2>
double EFTGenHistsWithCuts::getdR(T1 p1, T2 p2) {
    double dR = (p1.p4().Eta() - p2.p4().Eta())*(p1.p4().Eta() - p2.p4().Eta());
    dR += (getdPhi(p1,p2)*getdPhi(p1,p2));
    dR = sqrt(dR);
    return dR;
}

// Get lepton category of event
// NOTE: Right now just returns either 2lss, 3l, 4l (and 4l is for 4 or more l), or none
template <typename T>
TString EFTGenHistsWithCuts::getLepCat(const std::vector<T>& leptons) {
    TString lep_cat_name;
    const std::vector<T>& charged_leptons = getChargedParticles(leptons);
    int ch_sum = getChargeSum(leptons);
    if (charged_leptons.size() == 2 and ch_sum != 0){
        lep_cat_name = "2lss";
    } else if (charged_leptons.size() == 3) {
        lep_cat_name = "3l";
    } else if (charged_leptons.size() > 3) {
        lep_cat_name = "4l";
    } else {
        lep_cat_name = "none";
    }
    return lep_cat_name;
}

// Get charge sum of particles
template <typename T>
int EFTGenHistsWithCuts::getChargeSum(const std::vector<T>& particles_vect) {
    int ch_sum = 0;
    for (size_t i = 0; i < particles_vect.size(); i++) {
        const T& p = particles_vect.at(i);
        ch_sum = ch_sum + p.charge();
    }
    return ch_sum;
}

// Returns the charged particles from particles of type reco::GenParticleCollection
template <typename T>
std::vector<T> EFTGenHistsWithCuts::getChargedParticles(const std::vector<T>& particles_vect) {
    std::vector<T> ret;
    for (size_t i = 0; i < particles_vect.size(); i++) {
        const T& p =  particles_vect.at(i);
        if (p.charge() != 0){
            ret.push_back(p);
        }
    }
    return ret;
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

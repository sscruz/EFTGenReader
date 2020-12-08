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
        //std::vector<reco::GenJet> GetGenJetsFromDR(const std::vector<reco::GenJet>& gen_jets, const reco::GenParticleCollection& b_quarks, double dR_threshold);

        ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> GetSumTLV(reco::GenParticleCollection col);
        double GetdPhi(reco::GenParticle p1, reco::GenParticle p2);
        double GetInvMass(reco::GenParticle p1, reco::GenParticle p2);

        TString ConstructHistName(TString lep_cat, TString hist_type, std::vector<size_t> hist_idx_vect);
        void FillHistIfExists(TString h_name, double val, WCFit eft_fit);
        void FillTH1DHistIfExists(TString h_name, double val);

        bool has_substr(TString s, TString substr);


    public:
        explicit EFTGenHistsWithCuts(const edm::ParameterSet&);
        ~EFTGenHistsWithCuts();

        template <typename T> edm::Handle<T> get_collection(const edm::Event& event, const edm::EDGetTokenT<T>& token);
        template <typename T> std::vector<T> MakeBaselinePtEtaCuts(const std::vector<T>& input, double min_pt, double max_eta);
        template <typename T> std::vector<T> MakeStaggeredPtCuts(const std::vector<T>& particles, std::vector<double> pt_vector, double pt_min);
        template <typename T1, typename T2> std::vector<T1> CleanCollection(const std::vector<T1>& obj1_vector, const std::vector<T2>& obj2_vector, const double coneSize);
        template <typename T1, typename T2> double getdR(T1 p1, T2 p2);

        template <typename T> TString GetLepCat(const std::vector<T>& leptons);
        template <typename T> int GetChargeSum(const std::vector<T>& particles_vect);
        template <typename T> std::vector<T> GetChargedParticles(const std::vector<T>& particles_vect);
        template <typename T1, typename T2, typename T3> TString GetAnaCat(const std::vector<T1>& leptons, const std::vector<T2>& jets, const std::vector<T3>& bjets);
        template <typename T> bool IsSFOSZ(const std::vector<T>& leptons);
        template <typename T1, typename T2> int GetNJetsForLepCat(const std::vector<T1>& leptons, const std::vector<T2>& jets);

        std::ofstream fout;
        FILE * ffout;

        std::string sampleName;
        int sampleNumber;
        int eventcount;

        // runtime configurable parameters
        bool iseft;
        bool debug;

        // Cuts and parameters
        double min_pt_jet;
        double min_pt_lep;
        double max_eta_jet;
        double max_eta_lep;
        std::vector<double> staggered_pt_cuts_lep;
        int min_njets_2lss;
        int min_njets_3l;
        int min_njets_4l;
        int max_njet_bins_2lss;
        int max_njet_bins_3l;
        int max_njet_bins_4l;

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

        // Particle level jets hists
        TH2EFT* h_bjet_jetEFT; TH2D* h_bjet_jetSM;
        TH2EFT* h_bjet_jet2lEFT; TH2D* h_bjet_jet2lSM;
        TH2EFT* h_bjet_jet3lEFT; TH2D* h_bjet_jet3lSM;
        TH2EFT* h_bjet_jet4lEFT; TH2D* h_bjet_jet4lSM;

        TH2EFT* h_2lss_jetbjetEFT; TH2D* h_2lss_jetbjetSM;
        TH2EFT* h_3l_jetbjetEFT;   TH2D* h_3l_jetbjetSM;
        TH2EFT* h_3l_sfz_jetbjetEFT;   TH2D* h_3l_sfz_jetbjetSM;
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
        std::map<TString,TH1D*> hist_TH1D_dict;
        std::vector<std::string> lep_cats_vect {"2lss","3l","4l","anyLepCat"};

        std::vector<std::string> ana_cats_vct {"2lss-2b","3l-sfz-1b","3l-sfz-2b","3l-1b","3l-2b","4l-2b"}; // Sort of hard coded for now...

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

        int mom_id = id;    // If no mother, set to own id
        if (p_mom) mom_id = p_mom->pdgId();

        bool is_fromEWBoson = (abs(mom_id) >= 22 && abs(mom_id) <= 25);
        bool is_fromTopSystem = (abs(mom_id) == 24 && abs(gmom_noFSR_id) == 6);

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

/*
// TEST function for identifying GEN b jets with DR matching
// Not currently used, and if we eventually do want to use it, should probably test it further
// One thing in particular, shoul add check to make sure there is no double counting
std::vector<reco::GenJet> EFTGenHistsWithCuts::GetGenJetsFromDR(const std::vector<reco::GenJet>& gen_jets, const reco::GenParticleCollection& b_quarks, double dR_threshold) {
    //std::cout << "Threshold: " << dR_threshold << std::endl;
    std::vector<reco::GenJet> b_gen_jets;
    double min_dR;
    for(auto &j: gen_jets){
        //std::cout << "\nstarting new jet..." << std::endl;
        min_dR = 999990;
        for(auto &b: b_quarks){
            //std::cout << "\t dr is: " << getdR(j,b) << std::endl;
            if (getdR(j,b) < min_dR){
                min_dR = getdR(j,b);
                //std::cout << "\t\t NEW MIN dr is: " << getdR(j,b) << std::endl;
            }
        }
        //std::cout << "\t Min dr was: " << min_dR << std::endl;
        if (min_dR < dR_threshold){
            b_gen_jets.push_back(j);
            //std::cout << "\tPushing back this jet, id as b jet" << std::endl;
        }
    }
    return b_gen_jets;
}
*/

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

// For PL jets, does not work for GEN jets
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


ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> EFTGenHistsWithCuts::GetSumTLV(reco::GenParticleCollection col) {
    ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<double>> p4vec;
    for (const auto it: col) p4vec += it.p4();
    return p4vec;
}

double EFTGenHistsWithCuts::GetdPhi(reco::GenParticle p1, reco::GenParticle p2) {
    double dPhi = p2.p4().Phi() - p1.p4().Phi();
    double pi = 2.0*asin(1.0);
    if (dPhi>=pi) dPhi = dPhi - 2.0*pi;     // from TVector2
    if (dPhi<(-pi)) dPhi = dPhi + 2.0*pi;   // from TVector2
    return dPhi;
}

double EFTGenHistsWithCuts::GetInvMass(reco::GenParticle p1, reco::GenParticle p2) {
    auto p4vec = p1.p4() + p2.p4();
    return p4vec.M();
}

// Make a standardized hist name out of category info and hist type and hist number
TString EFTGenHistsWithCuts::ConstructHistName(TString lep_cat, TString hist_type, std::vector<size_t> hist_idx_vect){
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
void EFTGenHistsWithCuts::FillHistIfExists(TString h_name, double val, WCFit eft_fit){
    if (hist_dict.find(h_name) != hist_dict.end()){
        hist_dict[h_name]->Fill(val,1.0,eft_fit);
    }
}

// Fill TH1Ds if they exist in the hist dictionary
void EFTGenHistsWithCuts::FillTH1DHistIfExists(TString h_name, double val){
    if (hist_TH1D_dict.find(h_name) != hist_TH1D_dict.end()){
        hist_TH1D_dict[h_name]->Fill(val);
    }
}

// Check if string has given substring
bool EFTGenHistsWithCuts::has_substr(TString s, TString substr){
    return (s.Index(substr) != -1);
}

// Template function for making pt and eta cuts
template <typename T>
std::vector<T> EFTGenHistsWithCuts::MakeBaselinePtEtaCuts(const std::vector<T>& input, double min_pt, double max_eta) {
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

// Template function for making staggered pt cuts
template <typename T>
std::vector<T> EFTGenHistsWithCuts::MakeStaggeredPtCuts(const std::vector<T>& particles, std::vector<double> pt_vector, double pt_min) {
    std::vector<T> sorted_particles = particles;
    std::sort(sorted_particles.begin(),sorted_particles.end(), [] (T a, T b) { return a.p4().Pt() > b.p4().Pt();});
    std::vector<T> ret;
    for (size_t i = 0; i < sorted_particles.size(); i++) {
        const T& p = sorted_particles.at(i);
        double pt_threshold;
        if (i < pt_vector.size()){
            pt_threshold = pt_vector.at(i);
        } else {
            pt_threshold = pt_min;
        }
        if (p.p4().Pt() > pt_threshold) {
            ret.push_back(p);
        }
    }
    return ret;
}

// Clean objects based on how close they are to other objects
template <typename T1, typename T2> 
std::vector<T1> EFTGenHistsWithCuts::CleanCollection(const std::vector<T1>& obj1_vector, const std::vector<T2>& obj2_vector, const double coneSize) {
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
    dR += (GetdPhi(p1,p2)*GetdPhi(p1,p2));
    dR = sqrt(dR);
    return dR;
}

// Get lepton category of event
// NOTE: Right now just returns either 2lss, 3l, 4l (and 4l is for 4 or more l), or none
template <typename T>
TString EFTGenHistsWithCuts::GetLepCat(const std::vector<T>& leptons) {

    std::vector<T> leptons_ch;
    leptons_ch = GetChargedParticles(leptons);                                     // Make sure only looking at charged particles
    leptons_ch = MakeBaselinePtEtaCuts(leptons_ch,min_pt_lep,max_eta_lep);         // Make sure eta cuts have been made
    leptons_ch = MakeStaggeredPtCuts(leptons_ch,staggered_pt_cuts_lep,min_pt_lep); // Make sure appropriate analysis pt cuts have been made

    // Check which lepton category, if any, the event satisfies
    TString lep_cat_name = "none";
    int ch_sum = GetChargeSum(leptons_ch);
    if (leptons_ch.size() == 2 and ch_sum != 0){
        lep_cat_name = "2lss";
    } else if (leptons_ch.size() == 3) {
        lep_cat_name = "3l";
    } else if (leptons_ch.size() > 3) {
        lep_cat_name = "4l";
    }
    return lep_cat_name;
}

// Returns true if same flavor opposite sign pair in 3l event that's within z peak
template <typename T>
bool EFTGenHistsWithCuts::IsSFOSZ(const std::vector<T>& leptons) {
    bool is_sfoz = false;
    if (GetLepCat(leptons) == "3l"){
        const T& p1 = leptons.at(0);
        const T& p2 = leptons.at(1);
        const T& p3 = leptons.at(2);
        if(abs(p1.charge() + p2.charge() + p3.charge()) == 1) {
            //loop over all leptons
            for(size_t i = 0; i < leptons.size(); i++) {
                //loop over leptons not seen yet (prevents (1,2) (2,1) duplicates)
                for(size_t j = i+1; j < leptons.size(); j++) {
                    if (abs(leptons.at(i).pdgId()) != abs(leptons.at(j).pdgId())) continue; //look for same flavor only (Z->ee or Z->mumu)
                    if (leptons.at(i).charge() * leptons.at(j).charge() > 0) continue; //look for opposite sign only (e+e- or mu+mu-)
                    if (fabs((leptons.at(i).p4() + leptons.at(j).p4()).M() - 91.2) < 10) {
                        is_sfoz = true;
                        break; //done searching
                    }
                }
                if (is_sfoz) break; //no need to search any further
            }
        }
    }
    return is_sfoz;
}

// Get the analysis cateogry the event falls into
// Not sure if this is the right idea or not for how to set up this function...
template <typename T1, typename T2, typename T3>
TString EFTGenHistsWithCuts::GetAnaCat(const std::vector<T1>& leptons, const std::vector<T2>& jets, const std::vector<T3>& bjets){

    int njets = jets.size();
    int nbjets = bjets.size();
    TString cat = "none";
    TString lepcat = GetLepCat(leptons);

    // Get nbjets catetory
    TString bcat = "";
    if (nbjets == 1){
        bcat = "1b";
    } else if (nbjets > 1){
        bcat = "2b";
    }

    // Get the analysis category
    if (lepcat == "2lss"){
        if (njets >= min_njets_2lss){
            if (nbjets >=2){
                cat = lepcat+"-"+bcat;;
            }
        }
    } else if (lepcat == "3l"){
        if (njets >= min_njets_3l){
            if (nbjets >= 1){
                if (IsSFOSZ(leptons)){
                    cat = lepcat+"-"+"sfz"+"-"+bcat;
                } else {
                    cat = lepcat+"-"+bcat;
                }
            }
        }
    } else if (lepcat == "4l"){
        if (njets >= min_njets_4l){
            if (nbjets >= 2){
                cat = lepcat+"-"+bcat;
            }
        }
    }
    return cat;
}

// Returns number of jets, accounting for TOP-10-001 binning
template <typename T1, typename T2>
int EFTGenHistsWithCuts::GetNJetsForLepCat(const std::vector<T1>& leptons, const std::vector<T2>& jets){
    // These things should probably go somewhere that's not here:
    int njets = jets.size();
    int njets_ret;
    TString lepcat = GetLepCat(leptons);
    if ( (lepcat == "2lss") and (njets > max_njet_bins_2lss) ){
        njets_ret = max_njet_bins_2lss;
    } else if ( (lepcat == "3l") and (njets > max_njet_bins_3l) ){
        njets_ret = max_njet_bins_3l;
    } else if ( (lepcat == "4l") and (njets > max_njet_bins_4l) ){
        njets_ret = max_njet_bins_4l;
    } else {
        njets_ret = njets;
    }
    return njets_ret;
}


// Get charge sum of particles
template <typename T>
int EFTGenHistsWithCuts::GetChargeSum(const std::vector<T>& particles_vect) {
    int ch_sum = 0;
    for (size_t i = 0; i < particles_vect.size(); i++) {
        const T& p = particles_vect.at(i);
        ch_sum = ch_sum + p.charge();
    }
    return ch_sum;
}

// Returns the charged particles from particles of type reco::GenParticleCollection
template <typename T>
std::vector<T> EFTGenHistsWithCuts::GetChargedParticles(const std::vector<T>& particles_vect) {
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

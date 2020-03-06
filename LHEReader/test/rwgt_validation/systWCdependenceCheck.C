#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>
#include <unordered_map>
#include <set>

#include <math.h>

#include "TString.h"
#include "TSystemDirectory.h"
#include "TSystemFile.h"
#include "TFile.h"
#include "TList.h"
#include "TChain.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TF1.h"
#include "TGraph.h"
#include "TLatex.h"

#include "EFTGenReader/EFTHelperUtilities/interface/WCPoint.h"
#include "EFTGenReader/EFTHelperUtilities/interface/WCFit.h"
#include "EFTGenReader/EFTHelperUtilities/interface/TH1EFT.h"
#include "EFTGenReader/EFTHelperUtilities/interface/Stopwatch.h"
#include "EFTGenReader/EFTHelperUtilities/interface/split_string.h"
#include "makeEFTPlots.h"


void makePlot(TString sys, vector<WCPoint> pts_vect, vector<WCPoint> pts_vect_WgtU, vector<WCPoint> pts_vect_WgtD){
    WCFit* wc_fit = new WCFit(pts_vect,"test");
    WCFit* wc_fit_WgtU = new WCFit(pts_vect_WgtU,"test");
    WCFit* wc_fit_WgtD = new WCFit(pts_vect_WgtD,"test");

    std::string wc_name = "cpt"; // Should pass this to makePlot?
    TString save_name = sys+".png";
    //TString save_name = "TEST.png";
    TString plot_name = "";
    TString x_axis_name = wc_name+" Strength";
    TString y_axis_name = "\\sigma_{NP}/\\sigma_{SM}";
    int nom_clr = kBlack;
    int u_clr = kGreen;
    int d_clr = kBlue;
    int x_min = -20;
    int x_max = 20;
    float y_min = .9;

    TCanvas *c1 = new TCanvas("c1","",1200,800);
    c1->cd();
    c1->SetGrid(1,1);

    float left   = 0.32;
    float right  = 0.68;
    float top    = 0.88;
    float bottom = 0.8;
    TLegend *legend;
    legend = new TLegend(left,top,right,bottom);
    legend->SetNColumns(3);
    legend->SetBorderSize(0);

    // NOMINAL //
    double s0 = wc_fit->getCoefficient(wc_fit->kSMstr,wc_fit->kSMstr);
    double s1 = wc_fit->getCoefficient(wc_fit->kSMstr,wc_name);
    double s2 = wc_fit->getCoefficient(wc_name,wc_name);
    TF1* fit = new TF1("fit","pol2",x_min,x_max);
    fit->SetParameter(0,s0);
    fit->SetParameter(1,s1);
    fit->SetParameter(2,s2);
    fit->SetMinimum(y_min);
    fit->SetLineColor(nom_clr);
    fit->GetXaxis()->SetTitle(x_axis_name);
    fit->GetYaxis()->SetTitle(y_axis_name);
    fit->SetTitle(plot_name);
    fit->Draw();
    legend->AddEntry(fit,"Nominal","l");
    for (int i = 0; i < pts_vect.size(); i++) {
        WCPoint ref_pt = pts_vect.at(i);
        TGraph* ref_pt_gr = new TGraph(1);
        ref_pt_gr->SetPoint(0,ref_pt.getStrength(wc_name),ref_pt.wgt);
        ref_pt_gr->SetMarkerStyle(4);
        ref_pt_gr->SetMarkerColor(1);
        ref_pt_gr->Draw("P");
    }
    
    // UP //
    double s0_u = wc_fit_WgtU->getCoefficient(wc_fit_WgtU->kSMstr,wc_fit_WgtU->kSMstr);
    double s1_u = wc_fit_WgtU->getCoefficient(wc_fit_WgtU->kSMstr,wc_name);
    double s2_u = wc_fit_WgtU->getCoefficient(wc_name,wc_name);
    TF1* fit_u = new TF1("fit_u","pol2",-20,20);
    fit_u->SetParameter(0,s0_u);
    fit_u->SetParameter(1,s1_u);
    fit_u->SetParameter(2,s2_u);
    fit_u->SetLineColor(u_clr);
    float max_y_u = fit_u->Eval(x_max);
    float min_y_u = fit_u->Eval(0);
    fit_u->Draw("same");
    legend->AddEntry(fit_u,sys+" UP","l");
    for (int i = 0; i < pts_vect_WgtU.size(); i++) {
        WCPoint ref_pt_u = pts_vect_WgtU.at(i);
        TGraph* ref_pt_gr_u = new TGraph(1);
        ref_pt_gr_u->SetPoint(0,ref_pt_u.getStrength(wc_name),ref_pt_u.wgt);
        ref_pt_gr_u->SetMarkerStyle(4);
        ref_pt_gr_u->SetMarkerColor(u_clr);
        ref_pt_gr_u->Draw("P");
    }
    // DOWN //
    double s0_d = wc_fit_WgtD->getCoefficient(wc_fit_WgtD->kSMstr,wc_fit_WgtD->kSMstr);
    double s1_d = wc_fit_WgtD->getCoefficient(wc_fit_WgtD->kSMstr,wc_name);
    double s2_d = wc_fit_WgtD->getCoefficient(wc_name,wc_name);
    TF1* fit_d = new TF1("fit_d","pol2",-20,20);
    fit_d->SetParameter(0,s0_d);
    fit_d->SetParameter(1,s1_d);
    fit_d->SetParameter(2,s2_d);
    fit_d->SetLineColor(d_clr);
    float max_y_d = fit_d->Eval(x_max);
    float min_y_d = fit_d->Eval(0);
    fit_d->Draw("same");
    legend->AddEntry(fit_d,sys+" DOWN","l");
    for (int i = 0; i < pts_vect_WgtD.size(); i++) {
        WCPoint ref_pt_d = pts_vect_WgtD.at(i);
        TGraph* ref_pt_gr_d = new TGraph(1);
        ref_pt_gr_d->SetPoint(0,ref_pt_d.getStrength(wc_name),ref_pt_d.wgt);
        ref_pt_gr_d->SetMarkerStyle(4);
        ref_pt_gr_d->SetMarkerColor(d_clr);
        ref_pt_gr_d->Draw("P");
    }
    // Draw, save, delete
    float max_val = std::max(max_y_u,max_y_d);
    float min_val = std::min(min_y_u,min_y_d);
    //fit->SetMaximum(max_val);
    //fit->SetMinimum(min_val);
    legend->Draw();
    c1->Print(save_name,"png");
    //c1->Print("TEST.png","png");
    delete legend;
    delete c1;
}

std::map<string,string> make_wcpt_run_map() {
    // Set up the wc point string (depends a lot on the naming scheme!):
    // Right now this is set up for the 5 cpt axis scan files!!!
    std::map<string,string> ref_pts_dict;
    std::string wcname = "cpt";
    float range_max = 15;
    float npts = 5;
    float step = (range_max*2)/(npts-1);
    //std::cout << "STEP " << step << std::endl;
    int run = 0;
    for (float wcval=-range_max; wcval<=range_max; wcval=wcval+step){
        ref_pts_dict["run"+std::to_string(run)] = "wcpt_"+wcname+"_"+std::to_string(wcval);
        run = run + 1;
        std::cout << "Run: " << run << " wc val: " << wcval << std::endl;
    }
    return ref_pts_dict;
}

// Returns a vector of 3 maps (for nominal, up, down), each map contains a vector of the 5 WCPoints for each systematic
std::vector<std::map<std::string,std::vector<WCPoint>>> get_WCpt_syst_maps(TString run_dirs_file){

    std::vector<std::string> sys_names {"psISR","psFSR","muR","muF","muRmuF","nnpdf"};

    std::map<std::string,std::vector<WCPoint>> selection_pts_map; // Should only have key "nominal", could be a vect but want same type as selection_pts_WgtU/D_map
    std::map<std::string,std::vector<WCPoint>> selection_pts_WgtU_map;
    std::map<std::string,std::vector<WCPoint>> selection_pts_WgtD_map;
    WCPoint sm_pt("smpt");
    std::vector<TString> all_dirs;
    TString fdir;
    std::ifstream input_filenames(run_dirs_file);
    std::cout << "Using dirs:" << std::endl;
    while (input_filenames >> fdir) {
        all_dirs.push_back(fdir);
        std::cout << "\t" << fdir << std::endl;
    }
    input_filenames.close();
    int n_pts = all_dirs.size();
    std::cout << "Number of points to plot (number of files provied): " << n_pts << std::endl;

    std::map<string,string> ref_pts_dict = make_wcpt_run_map(); // Very bad, depends very much on naming scheme!)
    /*
    // Set up the wc point string (depends a lot on the naming scheme!):
    // Right now this is set up for the 5 cpt axis scan files!!!
    std::map<string,string> ref_pts_dict;
    std::string wcname = "cpt";
    float range_max = 15;
    float npts = 5;
    float step = (range_max*2)/(npts-1);
    std::cout << "STEP" << step << std::endl;
    int run = 0;
    for (float wcval=-range_max; wcval<=range_max; wcval=wcval+step){
        ref_pts_dict["run"+std::to_string(run)] = "wcpt_"+wcname+"_"+std::to_string(wcval);
        run = run + 1;
        std::cout << "Run: " << run << " wc val: " << wcval << std::endl;
    }
    */

    // Loop over files
    for (int idx=0; idx<all_dirs.size(); idx++){
        fdir = all_dirs.at(idx);
        std::string run_dir = getRunDirectory(fdir.Data());
        std::cout << "[" << (idx+1) << "/" << all_dirs.size() << "] Full Path: " << fdir << std::endl;

        std::vector<std::string> words;
        split_string(run_dir,words,"_");
        if (words.size() != 4) {
            std::cout << "[WARNING] Skipping invalid run directory!" << std::endl;
            continue;
        }

        // Get the tags from the file names
        std::string process   = words.at(1);
        std::string grp_tag   = words.at(2);
        std::string run_label = words.at(3);
        std::cout << "Process: " << process << " Grp tag: " << grp_tag << " Run label: " << run_label << std::endl;

        // Chain together all root files in the run directory
        TChain chain("EFTLHEReader/summaryTree");
        auto dir_files = getFiles(fdir);
        for (auto fn: dir_files) {
            //std::cout << "\tFiles: " << fn << std::endl;
            TString fname = fdir + "/" + fn;
            chain.Add(fname);
        }

        // Set up number of events in loop
        int chain_entries = chain.GetEntries();
        int last_entry = chain_entries;
        if (chain_entries > 100000) { // 100k != 10min for some reason
            std::cout << "Chain_entries: " << chain_entries << std::endl;
            last_entry = 100000;
        }
        //last_entry = 100; // For testing
        std::cout << "Last_entry: " << last_entry << std::endl;
        int first_entry = 0;

        WCFit* wcFit_intree = 0;
        WCFit inclusive_fit;
        WCFit selection_fit;
        WCFit selection_fit_WgtU;
        WCFit selection_fit_WgtD;

        double originalXWGTUP_intree = -1.;
        double genLep_pt3_intree;
        double genJet_pt4_intree;

        double psISRweightUp_intree;
        double psISRweightDown_intree;
        double psFSRweightUp_intree;
        double psFSRweightDown_intree;

        double muRWeightUp_intree;
        double muRWeightDown_intree;
        double muFWeightUp_intree;
        double muFWeightDown_intree;
        double muRmuFWeightUp_intree;
        double muRmuFWeightDown_intree;
        double nnpdfWeightUp_intree;
        double nnpdfWeightDown_intree;

        chain.SetBranchAddress("originalXWGTUP",&originalXWGTUP_intree);
        chain.SetBranchAddress("wcFit",&wcFit_intree);
        chain.SetBranchAddress("genLep_pt3",&genLep_pt3_intree);
        chain.SetBranchAddress("genJet_pt4",&genJet_pt4_intree);

        chain.SetBranchAddress("psISRweightUp",&psISRweightUp_intree);
        chain.SetBranchAddress("psFSRweightUp",&psFSRweightUp_intree);
        chain.SetBranchAddress("psISRweightDown",&psISRweightDown_intree);
        chain.SetBranchAddress("psFSRweightDown",&psFSRweightDown_intree);
        chain.SetBranchAddress("muRWeightUp",&muRWeightUp_intree);
        chain.SetBranchAddress("muRWeightDown",&muRWeightDown_intree);
        chain.SetBranchAddress("muFWeightUp",&muFWeightUp_intree);
        chain.SetBranchAddress("muFWeightDown",&muFWeightDown_intree);
        chain.SetBranchAddress("muRmuFWeightUp",&muRmuFWeightUp_intree);
        chain.SetBranchAddress("muRmuFWeightDown",&muRmuFWeightDown_intree);
        chain.SetBranchAddress("nnpdfWeightUp",&nnpdfWeightUp_intree);
        chain.SetBranchAddress("nnpdfWeightDown",&nnpdfWeightDown_intree);

        // WC Point 
        std::cout << "Run label: " << run_label << " , Dictionary entry: " << ref_pts_dict[run_label] << std::endl;
        std::string pt_str = ref_pts_dict[run_label];
        WCPoint selection_pt = WCPoint(pt_str); // Remember to only add wgt for events that pass selection
        std::map<std::string,WCFit> wcfit_systs_U_map;
        std::map<std::string,WCFit> wcfit_systs_D_map;
        std::map<std::string,WCPoint> wcpt_systs_U_map;
        std::map<std::string,WCPoint> wcpt_systs_D_map;
        for (auto sys: sys_names) {
            wcpt_systs_U_map[sys] = WCPoint(pt_str);
            wcpt_systs_D_map[sys] = WCPoint(pt_str);
        }

        // Event loop:
        int selection_events = 0;
        for (int i = first_entry; i < last_entry; i++) {
            chain.GetEntry(i);
            inclusive_fit.addFit(*wcFit_intree);
            // Define syst dict:
            std::map<std::string,std::map<std::string,double>> sys_map {
                {"psISR"  , {{"u",psISRweightUp_intree}  , {"d",psISRweightDown_intree}}},
                {"psFSR"  , {{"u",psFSRweightUp_intree}  , {"d",psFSRweightDown_intree}}},
                {"muR"    , {{"u",muRWeightUp_intree}    , {"d",muRWeightDown_intree}}},
                {"muF"    , {{"u",muFWeightUp_intree}    , {"d",muFWeightDown_intree}}},
                {"muRmuF" , {{"u",muRmuFWeightUp_intree} , {"d",muRmuFWeightDown_intree}}},
                {"nnpdf"  , {{"u",nnpdfWeightUp_intree}  , {"d",nnpdfWeightDown_intree}}}
            };
            if (genLep_pt3_intree > 10 and genJet_pt4_intree > 30){
                selection_events = selection_events + 1;
                selection_fit.addFit(*wcFit_intree);
                selection_pt.wgt += originalXWGTUP_intree;

                for (auto it = sys_map.begin(); it != sys_map.end(); it++){
                    std::string sys = it->first;

                    // Scale the up fit
                    wcFit_intree->scale(sys_map[sys]["u"]);       // Scale wcFit by sys UP weight
                    wcfit_systs_U_map[sys].addFit(*wcFit_intree); // Add scaled fit
                    wcFit_intree->scale(1/sys_map[sys]["u"]);     // Put wcFit back to orig value

                    // Scale the down fit
                    wcFit_intree->scale(sys_map[sys]["d"]);       // Scale wcFit by sys DOWN weight
                    wcfit_systs_D_map[sys].addFit(*wcFit_intree); // Add scaled fit
                    wcFit_intree->scale(1/sys_map[sys]["d"]);     // Put wcFit back to orig value

                    // Update the WC pts
                    wcpt_systs_U_map[sys].wgt += originalXWGTUP_intree*sys_map[sys]["u"];
                    wcpt_systs_D_map[sys].wgt += originalXWGTUP_intree*sys_map[sys]["d"];

                }
            }
        }
        std::cout << "\n Selected events over total: " << selection_events << "/" << last_entry << "->" << (float)selection_events/last_entry << "\n" << std::endl;

        // Normalize to SM
        double SM_xsec_incl = inclusive_fit.evalPoint(&sm_pt);
        double SM_xsec_sel = selection_fit.evalPoint(&sm_pt);
        selection_pt.scale(1.0/SM_xsec_sel);
        selection_pts_map["nominal"].push_back(selection_pt);

        for (auto sys: sys_names){
            double SM_xsec_sel_WgtU = wcfit_systs_U_map[sys].evalPoint(&sm_pt);
            double SM_xsec_sel_WgtD = wcfit_systs_D_map[sys].evalPoint(&sm_pt);
            inclusive_fit.scale(1.0/SM_xsec_incl);
            selection_fit.scale(1.0/SM_xsec_sel);

            //wcpt_systs_U_map[sys].scale(1.0/SM_xsec_sel);
            //wcpt_systs_D_map[sys].scale(1.0/SM_xsec_sel);
            wcpt_systs_U_map[sys].scale(1.0/SM_xsec_sel_WgtU);
            wcpt_systs_D_map[sys].scale(1.0/SM_xsec_sel_WgtD);

            // Fill the vector of WC points
            selection_pts_WgtU_map[sys].push_back(wcpt_systs_U_map[sys]);
            selection_pts_WgtD_map[sys].push_back(wcpt_systs_D_map[sys]);
        }
    }
    std::vector<std::map<std::string,std::vector<WCPoint>>> return_vect = {selection_pts_map,selection_pts_WgtU_map,selection_pts_WgtD_map};
    return return_vect;
}

// Returns a vector of 3 maps (nominal, up, down), each map contains a vector of the 5 WCPoints for qCut sys
std::vector<std::map<std::string,std::vector<WCPoint>>> get_WCpt_qCut_maps(TString run_dirs_file){

    std::map<std::string,std::vector<WCPoint>> selection_pts_map;
    std::map<std::string,std::vector<WCPoint>> selection_pts_WgtU_map;
    std::map<std::string,std::vector<WCPoint>> selection_pts_WgtD_map;
    WCPoint sm_pt("smpt");

    // Set up the wc point string (depends a lot on the naming scheme!):
    std::map<string,string> ref_pts_dict = make_wcpt_run_map();

    // Print dirs
    std::vector<TString> all_dirs;
    TString fdir;
    std::ifstream input_filenames(run_dirs_file);
    std::cout << "Using dirs:" << std::endl;
    while (input_filenames >> fdir) {
        all_dirs.push_back(fdir);
        std::cout << "\t" << fdir << std::endl;
    }
    input_filenames.close();

    // Loop over files
    for (int idx=0; idx<all_dirs.size(); idx++){

        fdir = all_dirs.at(idx);
        std::string run_dir = getRunDirectory(fdir.Data());
        std::cout << "[" << (idx+1) << "/" << all_dirs.size() << "] Full Path: " << fdir << std::endl;
        std::vector<std::string> words;
        split_string(run_dir,words,"_");
        if (words.size() != 4) {
            std::cout << "[WARNING] Skipping invalid run directory!" << std::endl;
            continue;
        }

        // Get the tags from the file names
        std::string process   = words.at(1);
        std::string grp_tag   = words.at(2);
        std::string run_label = words.at(3);
        std::cout << "Process: " << process << " Grp tag: " << grp_tag << " Run label: " << run_label << std::endl;

        // Check qCut type
        std::string f_type = "";
        if (grp_tag.find("qCut19") != std::string::npos){
            f_type = "qCut_nom";
            std::cout << "NOMINAL !!! " << grp_tag << std::endl;
        } else if (grp_tag.find("qCut25") != std::string::npos){
            f_type = "qCut_up";
            std::cout << "up !!! " << grp_tag << std::endl;
        } else if (grp_tag.find("qCut15") != std::string::npos){
            f_type = "qCut_down";
            std::cout << "down !!! " << grp_tag << std::endl;
        }

        // Chain together all root files in the run directory
        TChain chain("EFTLHEReader/summaryTree");
        auto dir_files = getFiles(fdir);
        for (auto fn: dir_files) {
            //std::cout << "\tFiles: " << fn << std::endl;
            TString fname = fdir + "/" + fn;
            chain.Add(fname);
        }

        // Set number of events in loop
        int chain_entries = chain.GetEntries();
        int last_entry = chain_entries;
        if (chain_entries > 100000) { // 100k != 10min for some reason
            std::cout << "Chain_entries: " << chain_entries << std::endl;
            last_entry = 100000;
        }
        //last_entry = 10; // For testing
        std::cout << "Last_entry: " << last_entry << std::endl;
        int first_entry = 0;

        WCFit* wcFit_intree = 0;
        WCFit selection_fit;
        WCFit selection_fit_WgtU;
        WCFit selection_fit_WgtD;

        double originalXWGTUP_intree = -1.;
        double genLep_pt3_intree;
        double genJet_pt4_intree;
        chain.SetBranchAddress("originalXWGTUP",&originalXWGTUP_intree);
        chain.SetBranchAddress("wcFit",&wcFit_intree);
        chain.SetBranchAddress("genLep_pt3",&genLep_pt3_intree);
        chain.SetBranchAddress("genJet_pt4",&genJet_pt4_intree);

        // Set up the WC Points
        //std::cout << "Run label: " << run_label << " , Dictionary entry: " << ref_pts_dict[run_label] << std::endl;
        std::string pt_str = ref_pts_dict[run_label];
        WCPoint selection_pt = WCPoint(pt_str); // Remember to only add wgt for events that pass selection
        std::map<std::string,WCFit> wcfit_systs_U_map;
        std::map<std::string,WCFit> wcfit_systs_D_map;
        std::map<std::string,WCPoint> wcpt_systs_U_map;
        std::map<std::string,WCPoint> wcpt_systs_D_map;
        wcpt_systs_U_map["qCut"] = WCPoint(pt_str);
        wcpt_systs_D_map["qCut"] = WCPoint(pt_str);

        // Event loop:
        int selection_events = 0;
        for (int i = first_entry; i < last_entry; i++) {

            chain.GetEntry(i);
            if (genLep_pt3_intree > 10 and genJet_pt4_intree > 30){
                selection_events = selection_events + 1;

                // Nominal
                //if (grp_tag.find("qCut19") != std::string::npos){
                if (f_type == "qCut_nom"){
                    selection_fit.addFit(*wcFit_intree);
                    selection_pt.wgt += originalXWGTUP_intree;
                }

                // Scale the up fit and update WC pt
                //if (grp_tag.find("qCut25") != std::string::npos){
                if (f_type == "qCut_up"){
                    wcfit_systs_U_map["qCut"].addFit(*wcFit_intree);
                    wcpt_systs_U_map["qCut"].wgt += originalXWGTUP_intree;
                }

                // Scale the down fit and update WC pt
                //if (grp_tag.find("qCut15") != std::string::npos){
                if (f_type == "qCut_down"){
                    wcfit_systs_D_map["qCut"].addFit(*wcFit_intree);
                    wcpt_systs_D_map["qCut"].wgt += originalXWGTUP_intree;
                }
            }
        }
        std::cout << "\n Selected events over total: " << selection_events << "/" << last_entry << "->" << (float)selection_events/last_entry << "\n" << std::endl;
        // Normalize to SM

        if (f_type == "qCut_nom"){
            double SM_xsec_sel = selection_fit.evalPoint(&sm_pt);
            selection_pt.scale(1.0/SM_xsec_sel);
            selection_pts_map["nominal"].push_back(selection_pt);
        }

        //selection_fit.scale(1.0/SM_xsec_sel);

        if (f_type == "qCut_up"){
            double SM_xsec_sel_WgtU = wcfit_systs_U_map["qCut"].evalPoint(&sm_pt);
            wcpt_systs_U_map["qCut"].scale(1.0/SM_xsec_sel_WgtU);
            //wcpt_systs_U_map["qCut"].scale(1.0/SM_xsec_sel); // Don't know how to do this norm since nominal was a different MG run?
            selection_pts_WgtU_map["qCut"].push_back(wcpt_systs_U_map["qCut"]);
        }
        if (f_type == "qCut_down"){
            double SM_xsec_sel_WgtD = wcfit_systs_D_map["qCut"].evalPoint(&sm_pt);
            wcpt_systs_D_map["qCut"].scale(1.0/SM_xsec_sel_WgtD);
            //wcpt_systs_D_map["qCut"].scale(1.0/SM_xsec_sel); // Don't know how to do this norm since nominal was a different MG run?
            selection_pts_WgtD_map["qCut"].push_back(wcpt_systs_D_map["qCut"]);
        }

    }
    std::vector<std::map<std::string,std::vector<WCPoint>>> return_vect = {selection_pts_map,selection_pts_WgtU_map,selection_pts_WgtD_map};
    return return_vect;
}

void systWCdependenceCheck(TString run_dirs_file, TString run_dirs_qCuts_file) {
    gStyle->SetOptStat(0);
    std::vector<std::string> sys_names {"psISR","psFSR","muR","muF","muRmuF","nnpdf"};

    ///*    
    std::vector<std::map<std::string,std::vector<WCPoint>>> systs_pts_vect;
    systs_pts_vect = get_WCpt_syst_maps(run_dirs_file);
    std::map<std::string,std::vector<WCPoint>> selection_pts_map      = systs_pts_vect.at(0);
    std::map<std::string,std::vector<WCPoint>> selection_pts_WgtU_map = systs_pts_vect.at(1);
    std::map<std::string,std::vector<WCPoint>> selection_pts_WgtD_map = systs_pts_vect.at(2);

    for(auto sys : sys_names){
        makePlot(sys,selection_pts_map["nominal"],selection_pts_WgtU_map[sys],selection_pts_WgtD_map[sys]);
    }
    //*/
    
    ///*
    std::cout << "\n ---------- Starting on the q cut calculations ----------\n" << std::endl;
    std::vector<std::map<std::string,std::vector<WCPoint>>> qCut_pts_vect;
    qCut_pts_vect = get_WCpt_qCut_maps(run_dirs_qCuts_file);
    //std::map<std::string,std::vector<WCPoint>> selection_pts_map = qCut_pts_vect.at(0);
    std::map<std::string,std::vector<WCPoint>> selection_pts_qCutU_map = qCut_pts_vect.at(1);
    std::map<std::string,std::vector<WCPoint>> selection_pts_qCutD_map = qCut_pts_vect.at(2);
    makePlot("qCut",selection_pts_map["nominal"],selection_pts_qCutU_map["qCut"],selection_pts_qCutD_map["qCut"]);
    //*/

}

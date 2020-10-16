#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TCanvas.h"
#include "TCut.h"
#include "TROOT.h"
#include "TChain.h"
#include "TLegend.h"

#include <iostream>
#include <fstream>
#include <vector>

#include "EFTGenReader/EFTHelperUtilities/interface/Stopwatch.h"
#include "EFTGenReader/EFTHelperUtilities/interface/WCFit.h"
#include "EFTGenReader/EFTHelperUtilities/interface/WCPoint.h"
#include "EFTGenReader/EFTHelperUtilities/interface/TH1EFT.h"

void printProgress(int current_index, int total_entries, int interval=20) {
    if (current_index % max(int(total_entries*interval/100.),interval) == 0) {
        float fraction = 100.*current_index/total_entries;
        std::cout << int(fraction) << " % processed " << std::endl;
    }
}

void setcanvas(TCanvas *c1, TPad **pad) {
    c1->SetLeftMargin(0.0);
    c1->SetTopMargin(0.00);
    c1->SetRightMargin(0.00);
    c1->SetBottomMargin(0.0);

    pad[0] = new TPad("pad0","pad",0,0.5,0.5,1.0);
    pad[1] = new TPad("pad1","pad",0.5,0.5,1.0,1.0);
    pad[2] = new TPad("pad2","pad",0,0,0.5,0.5);
    pad[3] = new TPad("pad3","pad",0.5,0.0,1.0,0.5);
    for(int k=0; k<4; k++) {
        pad[k]->Draw();
    }
    return;
}

void setlegend(TLegend *legend, TH1D *hall, TH1D *hmult0, TH1D *hmult1, TH1D *hmult2, TH1D *hmult3) {
    //legend->SetTextSize(0.050);
    legend->SetTextSize(0.055);
    legend->SetBorderSize(0);
    //legend->SetTextFont(62);
    legend->SetLineColor(0);
    legend->SetLineStyle(1);
    legend->SetLineWidth(1);
    legend->SetFillColor(0);
    legend->SetFillStyle(1001);

    legend->AddEntry(hmult0,"0 partons");
    legend->AddEntry(hmult1,"1 parton");
    legend->AddEntry(hall,"Total");
    //legend->AddEntry(hmult2,"2 partons");
    //legend->AddEntry(hmult3,"3 partons");
    return;
}

std::vector<TH1EFT*> makeTH1EFTs(const char *name, TChain *tree, int djr_idx, string rwgt_str, int nbins, double xlow, double xhigh) {
////std::vector<TH1EFT> makeTH1EFTs(const char *name, TChain *tree, int djr_idx, string rwgt_str, int nbins, double xlow, double xhigh) {

    //std::string sample_type = "central_sample";
    std::string sample_type = "EFT_sample";

    std::string cat = "no_cuts";
    //std::string cat = "2lep";
    //std::string cat = "3lep";

    // Set max number of events for event loop
    //int max_events = -1; // Run over all events
    //int max_events = 30000; // Debug
    //int max_events = 100000; 
    int max_events = 1000; 
    //int max_events = 300000; // For FP files (for ttH, ttll, ttlnu)
    //int max_events = 900000; // For FP files (for tHq, tllq)

    Stopwatch sw;
    std::set<int> unique_runs;

    TH1EFT *hall   = new TH1EFT(TString::Format("hall_%s",name),"",nbins,xlow,xhigh);
    TH1EFT *hmult0 = new TH1EFT(TString::Format("hmult0_%s",name),"",nbins,xlow,xhigh);
    TH1EFT *hmult1 = new TH1EFT(TString::Format("hmult1_%s",name),"",nbins,xlow,xhigh);
    TH1EFT *hmult2 = new TH1EFT(TString::Format("hmult2_%s",name),"",nbins,xlow,xhigh);
    TH1EFT *hmult3 = new TH1EFT(TString::Format("hmult3_%s",name),"",nbins,xlow,xhigh);

    int nMEpartons_intree;
    int nMEpartonsFiltered_intree;
    double genWgt_intree;
    std::vector<double> *djrvalues_intree = 0;
    WCFit* wcFit_intree = 0;
    double origxsec_intree;
    int lumiBlock_intree;

    double genLep_pt1_intree;
    double genLep_pt2_intree;
    double genLep_pt3_intree;
    double genJet_pt1_intree;
    double genJet_pt2_intree;
    double genJet_pt3_intree;
    double genJet_pt4_intree;

    tree->SetBranchAddress("nMEpartons",&nMEpartons_intree);
    tree->SetBranchAddress("nMEpartonsFiltered",&nMEpartonsFiltered_intree);
    tree->SetBranchAddress("genWgt",&genWgt_intree);
    tree->SetBranchAddress("DJRValues",&djrvalues_intree);
    tree->SetBranchAddress("wcFit",&wcFit_intree);
    tree->SetBranchAddress("originalXWGTUP",&origxsec_intree);
    tree->SetBranchAddress("lumiBlock",&lumiBlock_intree);

    tree->SetBranchAddress("genLep_pt1",&genLep_pt1_intree);
    tree->SetBranchAddress("genLep_pt2",&genLep_pt2_intree);
    tree->SetBranchAddress("genLep_pt3",&genLep_pt3_intree);
    tree->SetBranchAddress("genJet_pt1",&genJet_pt1_intree);
    tree->SetBranchAddress("genJet_pt2",&genJet_pt2_intree);
    tree->SetBranchAddress("genJet_pt3",&genJet_pt3_intree);
    tree->SetBranchAddress("genJet_pt4",&genJet_pt4_intree);

    WCFit inclusive_fit; 
    WCPoint* wc_pt = new WCPoint(rwgt_str,0);

    // Event loop:
    int tot_events = tree->GetEntries();
    sw.start("Full loop");
    int n_tot = 0;
    int n_count = 0;
    for (int i = 0; i < tot_events; i++) {
        n_tot = i;
        if (i == max_events){
            break;
        }
        sw.start("Event loop");
        tree->GetEntry(i);
        unique_runs.insert(lumiBlock_intree);
        if (sample_type == "EFT_sample"){
            if (djr_idx >= djrvalues_intree->size()) {
                continue;
            }
        }

        if (cat != "no_cuts"){
            // Cut on pt (this applies to both 2lep and 3lep catagories)
            if (genLep_pt1_intree < 25.0) {
                continue; // All events need at least 1 lepton with pt greater than 25
            } else if (genLep_pt2_intree < 15.0) {
                continue; // All events need at least 2 leptons (wiht pt greater than 15)
            }
            if (cat == "2lep") {
                // 2 leptons + 4 jets
                if (genLep_pt3_intree > 10.0) {
                    continue; // We want exactly 2 leptons with pt greater than 10 (skip if 3 leptons in event)
                } else if (genJet_pt4_intree < 30.0) {
                    continue; // We want at least 4 jets with pt greater than 30
                }
            } else if (cat == "3lep") { 
                // We won't split up 3 and 4 lep catagories
                if (genLep_pt3_intree < 10) {
                    continue; // We want at least 3 leptons with pt greater than 10
                } else if (genJet_pt2_intree < 30.0 ) {
                    continue; // We want at least 2 jets with pt greater than 30
                }
            } 
        }

        if (sample_type == "central_sample") { // For any sample w/o EFT rwgting
            //std::cout << origxsec_intree << std::endl;
            std::vector<WCPoint> tmp_pts;
            WCPoint tmp_SM_pt("smpt",origxsec_intree);
            tmp_pts.push_back(tmp_SM_pt);
            WCFit tmp_fit(tmp_pts,"");
            inclusive_fit.addFit(tmp_fit);
        }

        if (sample_type == "EFT_sample"){
            double djr_val = log10(djrvalues_intree->at(djr_idx)); // Remember to un comment when not doing central sampel!!!
            inclusive_fit.addFit(*wcFit_intree);

            printProgress(i,tot_events,10);

            double eft_weight = wcFit_intree->evalPoint(wc_pt);

            bool mult0 = (nMEpartonsFiltered_intree == 0);
            bool mult1 = (nMEpartonsFiltered_intree == 1);
            bool mult2 = (nMEpartonsFiltered_intree == 2);
            bool mult3 = (nMEpartonsFiltered_intree == 3);

            sw.start("Hist fill");
            
            wcFit_intree->scale(genWgt_intree);
            if (genWgt_intree != 0) {
                hall->Fill(djr_val,genWgt_intree*eft_weight,*wcFit_intree);
                if (mult0) {
                    hmult0->Fill(djr_val,genWgt_intree*eft_weight,*wcFit_intree);
                } else if (mult1) {
                    hmult1->Fill(djr_val,genWgt_intree*eft_weight,*wcFit_intree);
                } else if (mult2) {
                    hmult2->Fill(djr_val,genWgt_intree*eft_weight,*wcFit_intree);
                } else if (mult3) {
                    hmult3->Fill(djr_val,genWgt_intree*eft_weight,*wcFit_intree);
                }
            }
            sw.lap("Hist fill");
            sw.lap("Event loop");
        }
        n_count = n_count + 1;

    }
    sw.stop("Full loop");
    sw.readAllTimers(true,"");
    sw.readAllTimers(false,"");

    double norm_factor;
    if (sample_type == "EFT_sample"){
        norm_factor = unique_runs.size(); // number of lumi blocks
    } else if (sample_type == "central_sample"){
        norm_factor = max_events; // this one is correct?
        //norm_factor = n_count;
    }
    //int nlumi_blocks = unique_runs.size();
    // Divide by number of groups of 500 events (e.g. 1000000 events / 500 = 200 blocks)
    /*
    double incl_xsec = inclusive_fit.evalPoint(wc_pt) / nlumi_blocks; 
    double incl_xsec_err = inclusive_fit.evalPointError(wc_pt) / nlumi_blocks;
    */
    double incl_xsec = inclusive_fit.evalPoint(wc_pt) / norm_factor; 
    double incl_xsec_err = inclusive_fit.evalPointError(wc_pt) / norm_factor;

    std::cout << "\nXsec: " << incl_xsec << " +- " << incl_xsec_err 
              //<< " pb (events: " << n_count << "/" << max_events << " -> " << (float)n_count/max_events 
              << " pb (events: " << n_count << "/" << n_tot << " -> " << (float)n_count/n_tot 
              << ", type: " << sample_type 
              << ", cat: " << cat 
              << ", rwgt_str: \"" << rwgt_str << "\")\n"
              << std::endl;

    //hall->Scale(1.0/SM_norm)
    //hmult->0Scale(1.0/SM_norm);
    //hmult->1Scale(1.0/SM_norm);
    //hmult->2Scale(1.0/SM_norm);
    //hmult->3Scale(1.0/SM_norm);

    vector<TH1EFT*> hist_list; 
    hist_list.push_back(hall);
    hist_list.push_back(hmult0);
    hist_list.push_back(hmult1);
    hist_list.push_back(hmult2);
    hist_list.push_back(hmult3);
    //vector<TH1EFT> hist_list; 
    //hist_list.push_back(*hall);
    //hist_list.push_back(*hmult0);
    //hist_list.push_back(*hmult1);
    //hist_list.push_back(*hmult2);
    //hist_list.push_back(*hmult3);

    delete wc_pt;
    return hist_list;

}

//void makeplot(string rwgt_str,  const char *xlabel, std::vector<TH1EFT*> hist_list) {
void makeplot(WCPoint *wc_pt,  const char *xlabel, std::vector<TH1EFT*> hist_list) {
//void makeplot_Pointers(WCPoint *wc_pt,  const char *xlabel, std::vector<TH1EFT*> hist_list) {

    TH1EFT* hall   = hist_list[0];
    TH1EFT* hmult0 = hist_list[1];
    TH1EFT* hmult1 = hist_list[2];
    TH1EFT* hmult2 = hist_list[3];
    TH1EFT* hmult3 = hist_list[4];

    //WCPoint* wc_pt = new WCPoint(rwgt_str,0);
    for (Int_t bin_idx = 0; bin_idx <= hall->GetNbinsX()+1; bin_idx++) {
        double wcfit_bin_val = hall->GetBinFit(bin_idx).evalPoint(wc_pt);
        double wcfit_bin_err = hall->GetBinFit(bin_idx).evalPointError(wc_pt);
        hall->SetBinContent(bin_idx,wcfit_bin_val);
        hall->SetBinError(bin_idx,wcfit_bin_err);
    }
    for (Int_t bin_idx = 0; bin_idx <= hmult0->GetNbinsX()+1; bin_idx++) {
        double wcfit_bin_val = hmult0->GetBinFit(bin_idx).evalPoint(wc_pt);
        double wcfit_bin_err = hmult0->GetBinFit(bin_idx).evalPointError(wc_pt);
        hmult0->SetBinContent(bin_idx,wcfit_bin_val);
        hmult0->SetBinError(bin_idx,wcfit_bin_err);
    }
    for (Int_t bin_idx = 0; bin_idx <= hmult1->GetNbinsX()+1; bin_idx++) {
        double wcfit_bin_val = hmult1->GetBinFit(bin_idx).evalPoint(wc_pt);
        double wcfit_bin_err = hmult1->GetBinFit(bin_idx).evalPointError(wc_pt);
        hmult1->SetBinContent(bin_idx,wcfit_bin_val);
        hmult1->SetBinError(bin_idx,wcfit_bin_err);
    }
    for (Int_t bin_idx = 0; bin_idx <= hmult2->GetNbinsX()+1; bin_idx++) {
        double wcfit_bin_val = hmult2->GetBinFit(bin_idx).evalPoint(wc_pt);
        double wcfit_bin_err = hmult2->GetBinFit(bin_idx).evalPointError(wc_pt);
        hmult2->SetBinContent(bin_idx,wcfit_bin_val);
        hmult2->SetBinError(bin_idx,wcfit_bin_err);
    }
    for (Int_t bin_idx = 0; bin_idx <= hmult3->GetNbinsX()+1; bin_idx++) {
        double wcfit_bin_val = hmult3->GetBinFit(bin_idx).evalPoint(wc_pt);
        double wcfit_bin_err = hmult3->GetBinFit(bin_idx).evalPointError(wc_pt);
        hmult3->SetBinContent(bin_idx,wcfit_bin_val);
        hmult3->SetBinError(bin_idx,wcfit_bin_err);
    }

    //std::cout << "VALUE AT MAX BIN hall: " << hall->GetMaximumBin() << std::endl;
    //std::cout << "VALUE AT MAX : " << hall->GetMaximum() << std::endl;
    //std::cout << "VALUE AT MAX hmult1: " << hmult1->GetMaximum() << std::endl;
    //std::cout << "VALUE AT MAX BIN for hmult1: " << hmult1->GetMaximumBin() << std::endl;
    double max_y = hall->GetBinContent(hall->GetMaximumBin());
    hall->GetYaxis()->SetRangeUser(0.01,50*(hall->GetBinContent(hall->GetMaximumBin()))); // Do not use GetMaximum() here since it will always be whatever we set it to first

    //hall->SetLineColor(921);
    //hmult0->SetLineColor(600);
    //hmult1->SetLineColor(629);
    hall->SetLineColor(14);
    hmult0->SetLineColor(4);
    hmult1->SetLineColor(2);
    hmult2->SetLineColor(419);
    hmult3->SetLineColor(810);

    hall->SetLineWidth(4);
    hmult0->SetLineWidth(2);
    hmult1->SetLineWidth(2);
    hmult0->SetLineStyle(2);
    hmult1->SetLineStyle(2);
    hmult2->SetLineStyle(2);
    hmult3->SetLineStyle(2);

    hall->GetXaxis()->SetTitle(xlabel);

    // Adjust sizes
    hall->GetXaxis()->SetTitleSize(0.08);
    hall->GetXaxis()->SetTitleOffset(0.8);
    hall->GetXaxis()->SetLabelSize(0.055);
    hall->GetYaxis()->SetLabelSize(0.055);
    hall->GetXaxis()->SetLabelOffset(0.009);
    hall->GetYaxis()->SetLabelOffset(0.009);

    //TLegend *legend=new TLegend(0.67,0.87-4*0.06,0.87,0.87);
    //TLegend *legend=new TLegend(0.14,0.87-3*0.06,0.32,0.87); // left,top,right,bottom
    TLegend *legend=new TLegend(0.14,0.89,0.32,0.9-0.07*3); // left,top,right,bottom
    setlegend(legend, hall, hmult0, hmult1, hmult2, hmult3);

    hall->Draw("EHIST");
    hmult0->Draw("EHISTSAME");
    hmult1->Draw("EHISTSAME");
    //hmult2->Draw("EHISTSAME");
    //hmult3->Draw("EHISTSAME");
    gStyle->SetOptStat(0);
    gPad->SetLogy(1);

    legend->Draw();
    return;
}

// Version of this function with object instead of pointer
void makeplot_obj(WCPoint *wc_pt,  const char *xlabel, std::vector<TH1EFT*> hist_list) {

    TH1EFT hall   = *hist_list[0];
    TH1EFT hmult0 = *hist_list[1];
    TH1EFT hmult1 = *hist_list[2];
    TH1EFT hmult2 = *hist_list[3];
    TH1EFT hmult3 = *hist_list[4];

    //WCPoint* wc_pt = new WCPoint(rwgt_str,0);
    for (Int_t bin_idx = 0; bin_idx <= hall.GetNbinsX()+1; bin_idx++) {
        double wcfit_bin_val = hall.GetBinFit(bin_idx).evalPoint(wc_pt);
        double wcfit_bin_err = hall.GetBinFit(bin_idx).evalPointError(wc_pt);
        //hall->SetBinContent(bin_idx,wcfit_bin_val);
        //hall->SetBinError(bin_idx,wcfit_bin_err);
        hall.SetBinContent(bin_idx,wcfit_bin_val);
        hall.SetBinError(bin_idx,wcfit_bin_err);
    }
    for (Int_t bin_idx = 0; bin_idx <= hmult0.GetNbinsX()+1; bin_idx++) {
        double wcfit_bin_val = hmult0.GetBinFit(bin_idx).evalPoint(wc_pt);
        double wcfit_bin_err = hmult0.GetBinFit(bin_idx).evalPointError(wc_pt);
        //hmult0->SetBinContent(bin_idx,wcfit_bin_val);
        //hmult0->SetBinError(bin_idx,wcfit_bin_err);
        hmult0.SetBinContent(bin_idx,wcfit_bin_val);
        hmult0.SetBinError(bin_idx,wcfit_bin_err);
    }
    for (Int_t bin_idx = 0; bin_idx <= hmult1.GetNbinsX()+1; bin_idx++) {
        double wcfit_bin_val = hmult1.GetBinFit(bin_idx).evalPoint(wc_pt);
        double wcfit_bin_err = hmult1.GetBinFit(bin_idx).evalPointError(wc_pt);
        //hmult1->SetBinContent(bin_idx,wcfit_bin_val);
        //hmult1->SetBinError(bin_idx,wcfit_bin_err);
        hmult1.SetBinContent(bin_idx,wcfit_bin_val);
        hmult1.SetBinError(bin_idx,wcfit_bin_err);
    }
    for (Int_t bin_idx = 0; bin_idx <= hmult2.GetNbinsX()+1; bin_idx++) {
        double wcfit_bin_val = hmult2.GetBinFit(bin_idx).evalPoint(wc_pt);
        double wcfit_bin_err = hmult2.GetBinFit(bin_idx).evalPointError(wc_pt);
        //hmult2->SetBinContent(bin_idx,wcfit_bin_val);
        //hmult2->SetBinError(bin_idx,wcfit_bin_err);
        hmult2.SetBinContent(bin_idx,wcfit_bin_val);
        hmult2.SetBinError(bin_idx,wcfit_bin_err);
    }
    for (Int_t bin_idx = 0; bin_idx <= hmult3.GetNbinsX()+1; bin_idx++) {
        double wcfit_bin_val = hmult3.GetBinFit(bin_idx).evalPoint(wc_pt);
        double wcfit_bin_err = hmult3.GetBinFit(bin_idx).evalPointError(wc_pt);
        //hmult3->SetBinContent(bin_idx,wcfit_bin_val);
        //hmult3->SetBinError(bin_idx,wcfit_bin_err);
        hmult3.SetBinContent(bin_idx,wcfit_bin_val);
        hmult3.SetBinError(bin_idx,wcfit_bin_err);
    }

    std::cout << "VALUE AT MAX BIN hall: " << hall.GetMaximumBin() << std::endl;
    std::cout << "VALUE AT MAX : " << hall.GetMaximum() << std::endl;
    std::cout << "VALUE AT MAX hmult1: " << hmult1.GetMaximum() << std::endl;
    std::cout << "VALUE AT MAX BIN for hmult1: " << hmult1.GetMaximumBin() << std::endl;
    hall.GetYaxis()->SetRangeUser(0.001,hall.GetMaximum());
    hall.Draw("EHIST");
    std::cout << hall.GetBinContent(23) << std::endl;

    /*
    hall.SetLineColor(921);
    hmult0.SetLineColor(600);
    hmult1.SetLineColor(629);
    hmult2.SetLineColor(419);
    hmult3.SetLineColor(810);

    hall.SetLineWidth(2);
    hmult0.SetLineStyle(2);
    hmult1.SetLineStyle(2);
    hmult2.SetLineStyle(2);
    hmult3.SetLineStyle(2);

    hall.GetXaxis()->SetTitle(xlabel);

    // Adjust sizes
    hall.GetXaxis()->SetTitleSize(0.09);
    hall.GetXaxis()->SetTitleOffset(0.8);
    hall.GetXaxis()->SetLabelSize(0.06);
    hall.GetYaxis()->SetLabelSize(0.06);

    TLegend *legend=new TLegend(0.67,0.87-4*0.06,0.87,0.87);
    setlegend(legend, &hall, &hmult0, &hmult1, &hmult2, &hmult3);

    hall.Draw("EHIST");
    hmult0.Draw("EHISTSAME");
    hmult1.Draw("EHISTSAME");
    //hmult2.Draw("EHISTSAME");
    //hmult3.Draw("EHISTSAME");
    gStyle->SetOptStat(0);
    gPad->SetLogy(1);

    legend->Draw();
    */
    return;
}

void makeDJRHists(const TString & infile_spec, const TString & outfile, bool basePtSM, bool basePtRefPt) {

    // Configure the scan
    bool ctGscan     = true;  // Fill the dictionary if we are doing a ctG scan 
    bool WCscan      = true;   // Fill the dictionary if we are testing all WC at their ref values

    if ( basePtSM == 1 and basePtRefPt == 1 ) {
        std::cout << "\nERROR: Cannot set both basePtSM and basePtRefPt to true - please choose one base point to vary the WC with respect to\n" << std::endl;
        return; 
    } else if ( basePtSM == 0 and basePtRefPt == 0 ){
        std::cout << "\nERROR: No base point chosen - please choose one base point to vary the WC with respect to\n" << std::endl;
        return; 
    } else if  ( basePtSM ){
        std::cout << "\nBase point selected: SM point\n" << std::endl;
    } else if ( basePtRefPt ) {
        std::cout << "\nBase point selected: Ref point\n" << std::endl;
    }

    // Get the tree
    TH1::SetDefaultSumw2();
    TChain *tree = new TChain("EFTLHEReader/summaryTree");
    std::ifstream infiles(infile_spec);
    TString fn;
    while (infiles >> fn) {
        tree->Add(fn);
    }

    int nbins = 50.;
    double djrmin = -0.5;
    double djrmax = 3.;


    // Make the TH1EFts:
    std::vector<TH1EFT*> hist_list0 = makeTH1EFTs("djr0",tree,0,"",nbins,djrmin,djrmax);
    ////std::vector<TH1EFT> hist_list0 =  makeTH1EFTs("djr0",tree,0,"",nbins,djrmin,djrmax);
    //std::vector<TH1EFT*> hist_list0 = makeTH1EFTs("djr0",tree,0,"rwgt_ctG_250.0",nbins,djrmin,djrmax);
    //std::vector<TH1EFT*> hist_list1 = makeTH1EFTs("djr1",tree,1,"",nbins,djrmin,djrmax);
    //std::vector<TH1EFT*> hist_list2 = makeTH1EFTs("djr2",tree,2,"",nbins,djrmin,djrmax);
    //std::vector<TH1EFT*> hist_list3 = makeTH1EFTs("djr3",tree,3,"",nbins,djrmin,djrmax);

    // Return here when just finding xsec
    //return;

    // Declare the dictionary (for reweighting) and always include SM and RefPt in the dictionary. Also include a point from the lims in the AN: 
    std::map<string,std::pair<string,double> > rwgt_dict; 
    rwgt_dict["SM"] = std::make_pair("",0);
    rwgt_dict["RefPt"] = std::make_pair("",0);
    rwgt_dict["ANPt"] = std::make_pair("",0);
    rwgt_dict["LitPt"] = std::make_pair("",0);

    if ( ctGscan ) {
        string wc_name = "ctG";
        //double min = -3.5;
        //double max = 1.5;
        //double delta = 0.1;
        //double min = -250.0;
        //double max = 250.0;
        //double delta = 50.0;
        double min = -5.0;
        double max = 5.0;
        double delta = 1.0;
        double nsteps = (max-min)/delta;
        double wc_val;
        for (int i = 0; i <= nsteps; i++) {
            wc_val = min + i*delta;
            //std::cout << i << " " << wc_val << std::endl;
            string dict_key = wc_name + "_" + std::to_string(wc_val);
            rwgt_dict[dict_key] = std::make_pair(wc_name,wc_val);
        }
    }

    if ( WCscan ) {
        std::vector<string> wc_names_vect{ "ctW", "ctp", "cpQM", "ctei", "ctli", "cQei", "ctZ", "cQlMi", "cQl3i", "ctG", "ctlTi", "cbW", "cpQ3", "cptb", "cpt", "ctlSi"};
        //std::vector<double> ref_pt_vals_vect{ -8.303849, 64.337172, 45.883907, 24.328689, 24.43011, 23.757944, -6.093077, 23.951426, 21.540499, -3.609446, 21.809598, 49.595354, -51.106621, 136.133729, -43.552406, -20.005026}; // These are the values of the WC at the ref point
        std::vector<double> ref_pt_vals_vect{-4.0 , 41.0 , 29.0 , 7.0 , 8.0 , 7.0 , -4.0 , 7.0 , -8.0 , -2.0 , -2.0 , -5.0 , -10.0 , -18.0 , -25.0 , -9.0}; // These numbers are from table 18 of the AN
        for (int i = 0; i < wc_names_vect.size(); i++) {
            //cout << wc_names_vect.at(i) << " " << ref_pt_vals_vect.at(i) << endl;
            if (basePtSM ) {
                string dict_key = wc_names_vect.at(i) + "_" + std::to_string(ref_pt_vals_vect.at(i));
                rwgt_dict[dict_key] = std::make_pair(wc_names_vect.at(i),ref_pt_vals_vect.at(i));
            } else if (basePtRefPt ) {
                string dict_key = wc_names_vect.at(i) + "_0";
                rwgt_dict[dict_key] = std::make_pair(wc_names_vect.at(i),0);
            }
        }
    }

    // Print the dictionary: 
    std::cout << "\nRwgt dictionary:" << std::endl;
    for (auto it =rwgt_dict.begin(); it !=rwgt_dict.end(); it++ ) {
        std::pair<string,double> p = it->second;
        std::cout << "    " << it->first << ": " << p.first << " " << p.second << std::endl ;
    }
    std::cout << "\n" << std::endl;

    // Set the WCPoint points: 
    string sm_pt_str = "";
    string ref_pt_str = "rwgt_ctW_-8.303849_ctp_64.337172_cpQM_45.883907_ctei_24.328689_ctli_24.43011_cQei_23.757944_ctZ_-6.093077_cQlMi_23.951426_cQl3i_21.540499_ctG_-3.609446_ctlTi_21.809598_cbW_49.595354_cpQ3_-51.106621_cptb_136.133729_cpt_-43.552406_ctlSi_-20.005026";
    string an_pt_str = "rwgt_ctW_-4.0_ctp_41.0_cpQM_29.0_ctei_7.0_ctli_8.0_cQei_7.0_ctZ_-4.0_cQlMi_7.0_cQl3i_-8.0_ctG_-2.0_ctlTi_-2.0_cbW_-5.0_cpQ3_-10.0_cptb_-18.0_cpt_-25.0_ctlSi_-9.0";
    string arxiv1901_pt_str = "rwgt_ctG_0.4_ctW_-1.8_cbW_3.1_ctZ_4.0_cptb_-27_cpQ3_5.8_cpQM_-3.5_cpt_18_ctp_-60";
    WCPoint* sm_pt = new WCPoint(sm_pt_str);
    WCPoint* ref_pt = new WCPoint(ref_pt_str); 
    WCPoint* an_pt = new WCPoint(an_pt_str);
    WCPoint* lit_pt = new WCPoint(arxiv1901_pt_str);
    WCPoint* tmp_pt;

    // Set the base point: This is the point our scan changes values **with respect to** (either SM point, or ref point)
    WCPoint* base_pt;
    if ( basePtSM ){
        base_pt = sm_pt;
    } else if ( basePtRefPt ){
        base_pt = ref_pt;
    } 

    // Loop over the dictionary and rwgt the hists: 
    for (auto i = rwgt_dict.begin(); i != rwgt_dict.end(); i++ ) {

        // Append rwgt info to png name:
        TString outfile_TString = (class TString)outfile;
        TString outfile_TString_pdf = (class TString)outfile;
        TString rwgt_string_key = i->first;
        outfile_TString.Append(rwgt_string_key);
        outfile_TString_pdf.Append(rwgt_string_key);
        if ( basePtSM ){
            outfile_TString.Append("_basePtSM");
            outfile_TString_pdf.Append("_basePtSM");
        } else if (basePtRefPt ){
            outfile_TString.Append("_basePtRefPt");
            outfile_TString_pdf.Append("_basePtRefPt");
        }
        outfile_TString.Append(".png");
        outfile_TString_pdf.Append(".pdf");

        tmp_pt = base_pt;
        //double wc_val = 0;
        double orig_val;

        std::pair<string,double> rwgt_pair = i->second;

        // Set the rwgt pt that gets passed to makeplot:
        if (rwgt_string_key == "SM") {
            orig_val = sm_pt->getStrength(rwgt_pair.first);
            tmp_pt = sm_pt;
            std::cout << "    SM point" << std::endl;
            //wc_val = orig_val;
        }
        else if (rwgt_string_key == "RefPt") {
            orig_val = ref_pt->getStrength(rwgt_pair.first);
            tmp_pt = ref_pt;
            std::cout << "    Ref point" << std::endl;
            //wc_val = orig_val;
        } else if (rwgt_string_key == "ANPt") {
            orig_val = ref_pt->getStrength(rwgt_pair.first);
            tmp_pt = an_pt;
            std::cout << "    AN point" << std::endl;
        } else if (rwgt_string_key == "LitPt") {
            orig_val = ref_pt->getStrength(rwgt_pair.first);
            tmp_pt = lit_pt;
            std::cout << "    Literature point" << std::endl;
        } else {
            orig_val = base_pt->getStrength(rwgt_pair.first); 
            tmp_pt = base_pt;
            tmp_pt->setStrength(rwgt_pair.first,rwgt_pair.second);       
            std::cout << "    " << rwgt_string_key << ", " "rwgt point: " << rwgt_pair.first << ", " << rwgt_pair.second << std::endl;
        }

        //tmp_pt->setStrength(rwgt_pair.first,wc_val); // Set the base object to a new value

        std::cout << "\tOrig: " << orig_val << std::endl;
        std::cout << "\tNew: " << tmp_pt->getStrength(rwgt_pair.first) << std::endl;

        //TCanvas *c1 = new TCanvas("c1", "c1", 800, 600);
        TCanvas *c1 = new TCanvas("c1", "c1", 800, 700);

        // For just the first plot (0->1): 
        c1->SetLeftMargin(0.0);
        c1->SetTopMargin(0.00);
        c1->SetRightMargin(0.00);
        c1->SetBottomMargin(0.00);
        TPad* pad[1];
        pad[0] = new TPad("pad0","pad",0,0,1.0,1.0);
        pad[0]->Draw();
        pad[0]->cd();
        pad[0]->SetBottomMargin(.17);
        //makeplot(tmp_pt,"DJR 0->1",hist_list0);
        makeplot(tmp_pt,"DJR 0#rightarrow1",hist_list0);

        // For all 4 plots: 
        //TPad *pad[4];
        //setcanvas(c1,pad);
        //pad[0]->cd();
        //makeplot(tmp_pt,"DJR 0->1",hist_list0);
        //pad[1]->cd();
        //makeplot(tmp_pt,"DJR 1->2",hist_list1);
        //pad[2]->cd();
        //makeplot(tmp_pt,"DJR 2->3",hist_list2);
        //pad[3]->cd();
        //makeplot(tmp_pt,"DJR 3->4",hist_list3);

        c1->Print(outfile_TString);
        c1->Print(outfile_TString_pdf);
        delete c1;

        // Reset the base object to the original state:
        tmp_pt->setStrength(rwgt_pair.first,orig_val);

        std::cout << "\tReset: " << "tmp: " << tmp_pt->getStrength(rwgt_pair.first) << ", " << "base: " << base_pt->getStrength(rwgt_pair.first) << std::endl;
    }

    return;

}



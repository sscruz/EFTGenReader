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
#include "TGraphAsymmErrors.h"
#include "TLatex.h"
#include "TH1D.h"

#include "TGaxis.h"

#include "EFTGenReader/EFTHelperUtilities/interface/WCPoint.h"
#include "EFTGenReader/EFTHelperUtilities/interface/WCFit.h"
#include "EFTGenReader/EFTHelperUtilities/interface/TH1EFT.h"
#include "EFTGenReader/EFTHelperUtilities/interface/Stopwatch.h"
#include "EFTGenReader/EFTHelperUtilities/interface/split_string.h"
#include "makeEFTPlots.h"

// // // //
// NOTE the following: 
//    - if we ever look at a different set of WC values, would need to modify how we associate runs with WC strengths!
// // // //

// Checks if two vectors of WCPoints are the same size
void checkVectSizes(std::vector<WCPoint> v1, std::vector<WCPoint> v2){
    if ( v1.size() != v2.size() ){
        std::cout << "\nError: Vectors are not the same length! First is " << v1.size() << " while second is " << v2.size() << ". Exiting...\n" << std::endl;
        throw std::exception();
    }
}

// Returns a map (for "inc" and "dec") of vectors of WCpoints that are the max and min muR, muF, muRmuF
std::map<std::string,std::vector<WCPoint>> muRmuFenvelope(std::string wc_name, std::vector<WCPoint> v_nom, std::map<std::string,std::vector<WCPoint>> m1, std::map<std::string,std::vector<WCPoint>> m2){

    std::map<std::string,std::vector<WCPoint>> return_map;
    std::vector<std::string> s_lst {"muR","muF","muRmuF"};
    double wgt1, wgt2, wgt_min, wgt_max;
    for (int i=0; i<v_nom.size(); i++){
        WCPoint wcpt_inc;
        WCPoint wcpt_dec;
        wcpt_inc.setStrength(wc_name,v_nom.at(i).getStrength(wc_name));
        wcpt_dec.setStrength(wc_name,v_nom.at(i).getStrength(wc_name));
        wgt_min = 1000000;
        wgt_max = 0;
        for (auto sys : s_lst){
            wgt1 = m1[sys].at(i).wgt;
            wgt2 = m2[sys].at(i).wgt;
            //std::cout << "w1 and w2: " << wgt1 << " " << wgt2 << std::endl;
            if (wgt1 < wgt_min) {
                wgt_min = wgt1;
            }
            if (wgt1 > wgt_max) {
                wgt_max = wgt1;
            }
            if (wgt2 < wgt_min) {
                wgt_min = wgt2;
            }
            if (wgt2 > wgt_max) {
                wgt_max = wgt2;
            }
        }
        wcpt_inc.wgt = wgt_max;
        wcpt_dec.wgt = wgt_min;
        //std::cout << "w min  and w max: " << wgt_min << " " << wgt_max << std::endl;
        return_map["inc"].push_back(wcpt_inc);
        return_map["dec"].push_back(wcpt_dec);
    }
    return return_map;
}


// Returns a map of two vectors of WCPoints, one with all of the systs that inc the xsec added in quad, one wiht all the systs that dec the xsec added in quad
std::map<std::string,std::vector<WCPoint>> getQuadSums(std::string wc_name, std::vector<std::string> s_names, std::map<std::string,std::vector<WCPoint>> nom_map, std::map<std::string,std::vector<WCPoint>> up_map, std::map<std::string,std::vector<WCPoint>> down_map){

    map<std::string,std::vector<WCPoint>> return_map;
    std::map<std::string,std::vector<WCPoint>> return_map_inc;
    std::map<std::string,std::vector<WCPoint>> return_map_dec;
    double u , d , n , delta_inc , delta_dec , quad_sum_inc , quad_sum_dec;

    for (int i=0; i<nom_map["nominal"].size(); i++){
        //std::cout << "\npoint: " << i << "\n" << std::endl;

        quad_sum_inc = 0;
        quad_sum_dec = 0;
        WCPoint wcpt_inc;
        WCPoint wcpt_dec;
        wcpt_inc.setStrength(wc_name,nom_map["nominal"].at(i).getStrength(wc_name));
        wcpt_dec.setStrength(wc_name,nom_map["nominal"].at(i).getStrength(wc_name));

        n = nom_map["nominal"].at(i).wgt;
        for(auto sys : s_names){
            //std::cout << "sys: " << sys << std::endl;
            //std::cout << "WC Strength: " << nom_map["nominal"].at(i).getStrength("cpt") << " " << up_map[sys].at(i).getStrength("cpt") << " " << down_map[sys].at(i).getStrength("cpt") << std::endl;

            u = up_map[sys].at(i).wgt;
            d = down_map[sys].at(i).wgt;
            //std::cout << "weights nom, up, down: " << n << " " << u << " " << d << std::endl;

            if ( (u-n>0 and d-n>0) or (u-n<0 and d-n<0) ){ // Both u and d are in same dir
                //std::cout << "wgts in same dir" << std::endl;
                delta_inc = max( abs(u-n),abs(d-n) );
                delta_dec = max( abs(u-n),abs(d-n) );
            } else { // u and d are in opposite directions
                //std::cout << "wgts NOT in same dir" << std::endl;
                if (u>n) {
                    delta_inc = abs(u-n);
                    delta_dec = abs(d-n);
                } else if (u<n) {
                    delta_inc = abs(d-n);
                    delta_dec = abs(u-n);
                } else if ( u==n and d==n ){
                    delta_inc = 0;
                    delta_dec = 0;
                } else {
                    if (u==n){
                        std::cout << "up equals nom: " << u << " " << n << std::endl;
                    }
                    if (d==n){
                        std::cout << "down equals nom: " << d << " " << n << std::endl;
                    }
                    std::cout << "This should not happen! Fix the if statements to account for this possiblity. Exiting..." << std::endl;
                    throw std::exception();
                }
            }
            //std::cout << "delta inc and dec: " << delta_inc << " " << delta_dec << "\n" << std::endl;
            quad_sum_inc = quad_sum_inc + delta_inc*delta_inc;
            quad_sum_dec = quad_sum_dec + delta_dec*delta_dec;
        }
        //std::cout << "quad sum inc and dec: " << quad_sum_inc << " " << quad_sum_dec << std::endl;
        wcpt_inc.wgt = n+sqrt(quad_sum_inc);
        wcpt_dec.wgt = n-sqrt(quad_sum_dec);
        return_map["inc"].push_back(wcpt_inc);
        return_map["dec"].push_back(wcpt_dec);
    }
    return return_map;
}

// Returns a TH1D whose bins are f1/f2 in the range xmin to xmax
TH1D* divideTF1s(TString name, TF1* f1, TF1* f2, float xmin, float xmax){
    int nsteps = 1000;
    float step = (xmax - xmin)/nsteps;
    float xcoord;
    float ratio;
    TH1D* ratio_hst = new TH1D(name,"",nsteps,xmin,xmax);
    for (int i=0; i<nsteps; i++){
        xcoord = xmin + (step*i);
        ratio = f1->Eval(xcoord)/f2->Eval(xcoord);
        ratio_hst->SetBinContent(i+1,ratio);
    }
    return ratio_hst;
}

// Retrurns (f1(wcpt) - f2(wcpt))/f2(wcpt)
double getFracDiff(WCFit *f1, WCFit *f2, WCPoint wcpt_eval){
    double frac_diff = (f1->evalPoint(&wcpt_eval) - f2->evalPoint(&wcpt_eval) ) / f2->evalPoint(&wcpt_eval);
    //std::cout << frac_diff << " : " << f1->evalPoint(&wcpt_eval) << " " << f2->evalPoint(&wcpt_eval) << std::endl;
    return frac_diff;
}

// Print uncertainty size at specific points
void printFracUncty(std::string wc_name, TString sys, vector<WCPoint> pts_vect_nom, vector<WCPoint> pts_vect_WgtU, vector<WCPoint> pts_vect_WgtD){
    std::map<string,std::pair<double,double>> smefit_lims_dict {
        {"ctG"  , {-0.4, 0.4} },
        {"ctW"  , {-1.8, 0.9} },
        {"cbW"  , {-2.6, 3.1} },
        {"ctZ"  , {-2.1, 4.0} },
        {"cptb" , {-27, 8.7}  },
        {"cpQ3" , {-5.5, 5.8} },
        {"cpQM" , {-3.5, 3}   },
        {"cpt"  , {-13, 18}   },
        {"ctp"  , {-60, 10}   }
    };
    double wc_m = smefit_lims_dict[wc_name].first;
    double wc_p = smefit_lims_dict[wc_name].second;
    //double wc_m = -15;
    //double wc_p = 15;
    WCPoint wcpt_m_lim; // Minum limit
    WCPoint wcpt_p_lim; // Plus limit
    WCFit* wc_fit_nom  = new WCFit(pts_vect_nom,"test");
    WCFit* wc_fit_WgtU = new WCFit(pts_vect_WgtU,"test");
    WCFit* wc_fit_WgtD = new WCFit(pts_vect_WgtD,"test");

    wcpt_p_lim.setStrength(wc_name,wc_p);
    wcpt_m_lim.setStrength(wc_name,wc_m);
    double diff_p_avg = ( abs(getFracDiff(wc_fit_WgtU,wc_fit_nom,wcpt_p_lim)) + abs(getFracDiff(wc_fit_WgtD,wc_fit_nom,wcpt_p_lim)) ) /2;
    double diff_m_avg = ( abs(getFracDiff(wc_fit_WgtU,wc_fit_nom,wcpt_m_lim)) + abs(getFracDiff(wc_fit_WgtD,wc_fit_nom,wcpt_m_lim)) ) /2;
    std::cout << sys << " " << wc_p << " " << diff_p_avg << " " << wc_m << " " << diff_m_avg << " " << std::endl;
    //std::cout << "The fit vals at this point, nom, u, d: " << wc_fit_nom->evalPoint(&wcpt_p_lim) << " " << wc_fit_WgtU->evalPoint(&wcpt_p_lim) << " " << wc_fit_WgtD->evalPoint(&wcpt_p_lim) << std::endl;
    //std::cout << "The fit vals at this n point, nom, u, d: " << wc_fit_nom->evalPoint(&wcpt_m_lim) << " " << wc_fit_WgtU->evalPoint(&wcpt_m_lim) << " " << wc_fit_WgtD->evalPoint(&wcpt_m_lim) << std::endl;
}

// Get an error band between two TF1s (Not tested or used!!!)
TGraphAsymmErrors getErrorBandFromTF1s(TF1 *fit, TF1 *fit_u, TF1 *fit_d, double x_min, double x_max){
    const Int_t n = 100;
    double xcoord = 0;
    double step_size = (x_max-x_min)/float(n-1);
    Double_t x[n], y[n], exl[n], eyl[n], exh[n], eyh[n];
    for (int i=0; i<n; i++){
        xcoord = x_min + i*step_size;
        x[i] = xcoord;
        y[i] = fit->Eval(xcoord);
        // The errror band only makes sense in the same line is either always above or always below nominal
        // If this is not the case, these if statements do not cover that situation, we'll just get some negative error bars
        // Probably we should not be trying to draw error bands between u and d fits if they are crossing nominal anyway
        if (fit_u->Eval(xcoord) > fit->Eval(xcoord)){
            eyh[i] = fit_u->Eval(xcoord) - fit->Eval(xcoord);
            eyl[i] = fit->Eval(xcoord) - fit_d->Eval(xcoord);
        } else {
            eyh[i] = fit_d->Eval(xcoord) - fit->Eval(xcoord);
            eyl[i] = fit->Eval(xcoord) - fit_u->Eval(xcoord);
        }
        exl[i] = 0;
        exh[i] = 0;
    }
    auto error_band = new TGraphAsymmErrors(n,x,y,exl,exh,eyl,eyh);
    return *error_band;
}

// Get an error band between two TH1s
TGraphAsymmErrors* getErrorBandFromTH1s(TH1 *h, TH1 *h_u, TH1 *h_d){
    float xcoord;
    //const Int_t n = h->size();
    int bins = h->GetXaxis()->GetNbins();
    const Int_t n = bins;
    Double_t low = h->GetXaxis()->GetXmin();
    Double_t high = h->GetXaxis()->GetXmax();

    Double_t x[n], y[n], exl[n], eyl[n], exh[n], eyh[n];
    //for (int i = 1; i <= bins; i++) {
    for (int i = 0; i <= bins; i++) {
        xcoord = h->GetBinCenter(i);
        x[i] = xcoord;
        y[i] = h->GetBinContent(i);
        if (h_u->GetBinContent(i) > h->GetBinContent(i)){
            eyh[i] = h_u->GetBinContent(i) - h->GetBinContent(i);
            eyl[i] = h->GetBinContent(i) - h_d->GetBinContent(i);
        } else {
            eyh[i] = h_d->GetBinContent(i) - h->GetBinContent(i);
            eyl[i] = h->GetBinContent(i) - h_u->GetBinContent(i);
        }
        exl[i] = 0;
        exh[i] = 0;
        //std::cout << "i,x,y " << i << " " << x[i] << " " << y[i] << std::endl;
    }
    auto error_band  = new TGraphAsymmErrors(bins,x,y,exl,exh,eyl,eyh);
    return error_band;
}


// Make plots of the up down and nominal
//void makePlot(std::string wc_name, TString sys, vector<WCPoint> pts_vect, vector<WCPoint> pts_vect_WgtU, vector<WCPoint> pts_vect_WgtD){
void makePlot(std::string wc_name, TString sys, vector<vector<WCPoint>> pts_vect_vect, vector<vector<WCPoint>> pts_vect_WgtU_vect, vector<vector<WCPoint>> pts_vect_WgtD_vect, vector<TString> tag_vect={"",""}){
    std::cout << "\t" << sys << ":" << std::endl;

    bool include_ratio = true;
    bool draw_errorband = true;
    bool draw_updown = false;
    bool draw_points = false;

    TString save_name = sys;
    //TString save_name = sys+".png";
    //TString save_name_pdf = sys+".pdf";
    TString plot_name = "";
    TString x_axis_name = wc_name+" Strength";
    TString y_axis_name = "\\sigma_{NP}/\\sigma_{SM}";
    int nom_clr = kBlack;
    int u_clr = kGreen;
    int d_clr = kBlue;
    float x_min, x_max, y_min, y_max;
    if (wc_name=="cpt"){
        //x_min = -20;
        //x_max = 20;
        //y_min = .9;
        x_min = -13;
        x_max = 18;
        y_min = .95;
        y_max = 1.35;
    } else if (wc_name=="ctG"){
        //x_min = -2.5;
        //x_max = 2.5;
        //y_min = 0.0;
        //y_max = 15;
        x_min = -.4;
        x_max = 0.4;
        y_min = 0.6;
        y_max = 2.20;
    } else if (wc_name=="cpQ3"){
        x_min = -5.5;
        x_max = 5.8;
        //y_min = .95;
        //y_max = 1.18;
        y_min = 0.95;
        y_max = 1.25;
    } else {
        std::cout << "\nError: Unknown WC: " << wc_name << " , exiting.." << std::endl;
        throw std::exception();
    }

    int n_samples_to_plot = pts_vect_vect.size();
    if (n_samples_to_plot>1){
        include_ratio = false;       
        std::cout << "Not including ratio plot, too many samples." << std::endl;
    }

    //TCanvas *c1 = new TCanvas("c1","",1200,800);
    //TCanvas *c1 = new TCanvas("c1","",1000,900);
    TCanvas *c1 = new TCanvas("c1","",1000,1000);
    //TPad *pad1 = new TPad("pad1", "pad1", 0, 0.25, 1, 1.0); // xlow, ylow, xup, yup
    //TPad *pad2 = new TPad("pad2", "pad2", 0, 0.05, 1, 0.25); // xlow, ylow, xup, yup
    //TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0); // xlow, ylow, xup, yup
    //TPad *pad2 = new TPad("pad2", "pad2", 0, 0.00, 1, 0.3); // xlow, ylow, xup, yup
    double r = 0.40;
    double epsilon = 0.03;
    double pad_left_edge;
    double pad_right_edge;
    TPad *pad1 = new TPad("pad1", "pad1", 0, r-epsilon, 1, 1);
    TPad *pad2 = new TPad("pad2", "pad2", 0, 0, 1, r*(1-epsilon));
    if (include_ratio){
        pad1->SetLeftMargin(.17);  // Room for the legend
        //pad1->SetBottomMargin(.17); // Upper and lower plot are not joined
        pad1->SetBottomMargin(0.0);
        pad1->SetGrid(1,1);        // Grid
        pad1->Draw();              // Draw the upper pad: pad1
        c1->cd();                  // Go back to the main canvas before defining pad2
        pad2->SetLeftMargin(.17);  // Room for the legend
        //pad2->SetTopMargin(0.1);
        pad2->SetTopMargin(0.0);
        //pad2->SetBottomMargin(0.25);
        pad2->SetBottomMargin(0.45);
        pad2->SetGrid(1,1);
        pad2->Draw();
        pad1->cd();
        pad_left_edge = pad1->GetLeftMargin();
        pad_right_edge = 1-pad1->GetRightMargin();
    } else {
        c1->SetBottomMargin(.20);
        c1->SetLeftMargin(.20);
        c1->SetGrid(1,1);
        pad_left_edge = c1->GetLeftMargin();
        pad_right_edge = 1-c1->GetRightMargin();
    }

    // Set up the legend
    double ledgend_width = 0.5;
    double top, bottom, left, right;
    if (include_ratio){
        top    = 0.78;
        bottom = 0.62;
    } else {
        top = 0.85;
        bottom = 0.75;
    }
    left  = (pad_right_edge-pad_left_edge)/2 + pad_left_edge - ledgend_width/2; // Center the legend on the pad, not canvas
    right = (pad_right_edge-pad_left_edge)/2 + pad_left_edge + ledgend_width/2; // Center the legend on the pad, not canvas
    TLegend *legend;
    legend = new TLegend(left,top,right,bottom);
    if (draw_updown){
        legend->SetNColumns(3);
    } else {
        legend->SetNColumns(n_samples_to_plot);
    }
    legend->SetBorderSize(0);

    // Loop over all samples provided (up till now is just 1, or maybe 2 (if doing 0p and 0+1p)
    if (n_samples_to_plot != pts_vect_WgtU_vect.size() or n_samples_to_plot != pts_vect_WgtD_vect.size()){
        std::cout << "\nError: Mismatch is vector sizes passed to plotting function!" << std::endl;
        throw std::exception();
    }
    for (int sample_idx=0; sample_idx<n_samples_to_plot; sample_idx++){
        TString tag = tag_vect.at(sample_idx);
        std::cout << "Tag: " << tag << std::endl;
        vector<WCPoint> pts_vect = pts_vect_vect.at(sample_idx);            // Set of nominal points (for a given sample)
        vector<WCPoint> pts_vect_WgtU = pts_vect_WgtU_vect.at(sample_idx);  // Set of up variations (for a given sample)
        vector<WCPoint> pts_vect_WgtD = pts_vect_WgtD_vect.at(sample_idx);  // Set of down variations (for a given sample)
    
        WCFit* wc_fit = new WCFit(pts_vect,"test");
        WCFit* wc_fit_WgtU = new WCFit(pts_vect_WgtU,"test");
        WCFit* wc_fit_WgtD = new WCFit(pts_vect_WgtD,"test");

        // NOMINAL //
        double s0 = wc_fit->getCoefficient(wc_fit->kSMstr,wc_fit->kSMstr);
        double s1 = wc_fit->getCoefficient(wc_fit->kSMstr,wc_name);
        double s2 = wc_fit->getCoefficient(wc_name,wc_name);
        TF1* fit = new TF1("fit","pol2",x_min,x_max);
        fit->SetParameter(0,s0);
        fit->SetParameter(1,s1);
        fit->SetParameter(2,s2);
        fit->SetMinimum(y_min);
        fit->SetMaximum(y_max);
        fit->SetLineColor(nom_clr);
        fit->GetXaxis()->SetTitle(x_axis_name);
        fit->GetYaxis()->SetTitle(y_axis_name);
        if (wc_name == "ctG"){
            fit->GetXaxis()->SetNdivisions(6,2,0,true);
        }
        // Size and offset
        fit->GetXaxis()->SetTitleSize(0.060);
        fit->GetYaxis()->SetTitleSize(0.12);
        fit->GetYaxis()->SetTitleOffset(0.70);
        fit->GetXaxis()->SetTitleOffset(0.9);
        fit->GetYaxis()->SetLabelSize(0.065);
        fit->GetXaxis()->SetLabelSize(0.065);
        fit->GetYaxis()->SetLabelOffset(0.015);
        fit->GetXaxis()->SetLabelOffset(0.015);
        if (not include_ratio){
            fit->GetXaxis()->SetLabelSize(0.050);
            fit->GetYaxis()->SetLabelSize(0.050);
            fit->GetXaxis()->SetTitleOffset(1.1);
            fit->GetYaxis()->SetTitleOffset(1.0);
            fit->GetYaxis()->SetTitleSize(0.10);
        }
        fit->SetTitle(plot_name);
        if (sample_idx==0){
            fit->Draw();
        } else {
            fit->Draw("same");
        }
        if (tag != ""){
            legend->AddEntry(fit,tag,"l");
            if (tag.Index("0+1p") != -1){
                std::cout << "Custom 0+1p legend: " << tag << std::endl;
                fit->SetLineColor(kBlue);
                fit->SetLineWidth(3);
                fit->SetLineStyle(7); // Plot the 0+1p lines as dashed
            } else if (tag.Index("0p") != -1){
                std::cout << "Custom 0p legend: " << tag << std::endl;
                fit->SetLineColor(kBlue);
                fit->SetLineWidth(2);
            }
        } else {
            std::cout << "default legend: " << std::endl;
            legend->AddEntry(fit,"Nominal","l");
        }
        std::cout << "pts_vect.size() size: " << pts_vect.size() << std::endl;
        for (int i = 0; i < pts_vect.size(); i++) {
            WCPoint ref_pt = pts_vect.at(i);
            TGraph* ref_pt_gr = new TGraph(1);
            ref_pt_gr->SetPoint(0,ref_pt.getStrength(wc_name),ref_pt.wgt);
            std::cout << "NOMINAL VAL x:" << ref_pt.getStrength(wc_name) << " y: " << ref_pt.wgt << std::endl;
            ref_pt_gr->SetMarkerStyle(4);
            ref_pt_gr->SetMarkerColor(1);
            if (draw_points){
                ref_pt_gr->Draw("P");
            }
        }
        
        // UP //
        double s0_u = wc_fit_WgtU->getCoefficient(wc_fit_WgtU->kSMstr,wc_fit_WgtU->kSMstr);
        double s1_u = wc_fit_WgtU->getCoefficient(wc_fit_WgtU->kSMstr,wc_name);
        double s2_u = wc_fit_WgtU->getCoefficient(wc_name,wc_name);
        TF1* fit_u = new TF1("fit_u","pol2",x_min,x_max);
        fit_u->SetParameter(0,s0_u);
        fit_u->SetParameter(1,s1_u);
        fit_u->SetParameter(2,s2_u);
        fit_u->SetLineColor(u_clr);
        float max_y_u = fit_u->Eval(x_max);
        float min_y_u = fit_u->Eval(0);
        if (draw_updown){
            fit_u->Draw("same");
            if ( sys=="muRmuF_evnelope" or sys=="quad_sum"){
                legend->AddEntry(fit_u,sys+" HIGH","l");
            } else {
                legend->AddEntry(fit_u,sys+" UP","l");
            }
        }
        for (int i = 0; i < pts_vect_WgtU.size(); i++) {
            WCPoint ref_pt_u = pts_vect_WgtU.at(i);
            TGraph* ref_pt_gr_u = new TGraph(1);
            std::cout << "UP VAL x:" << ref_pt_u.getStrength(wc_name) << " y: " << ref_pt_u.wgt << std::endl;
            ref_pt_gr_u->SetPoint(0,ref_pt_u.getStrength(wc_name),ref_pt_u.wgt);
            ref_pt_gr_u->SetMarkerStyle(4);
            ref_pt_gr_u->SetMarkerColor(u_clr);
            if (draw_updown and draw_points) {
                ref_pt_gr_u->Draw("P");
            }
        }
        // DOWN //
        double s0_d = wc_fit_WgtD->getCoefficient(wc_fit_WgtD->kSMstr,wc_fit_WgtD->kSMstr);
        double s1_d = wc_fit_WgtD->getCoefficient(wc_fit_WgtD->kSMstr,wc_name);
        double s2_d = wc_fit_WgtD->getCoefficient(wc_name,wc_name);
        TF1* fit_d = new TF1("fit_d","pol2",x_min,x_max);
        fit_d->SetParameter(0,s0_d);
        fit_d->SetParameter(1,s1_d);
        fit_d->SetParameter(2,s2_d);
        fit_d->SetLineColor(d_clr);
        float max_y_d = fit_d->Eval(x_max);
        float min_y_d = fit_d->Eval(0);
        if (draw_updown){
            fit_d->Draw("same");
            if ( sys=="muRmuF_evnelope" or sys=="quad_sum"){
                legend->AddEntry(fit_d,sys+" LOW","l");
            } else {
                legend->AddEntry(fit_d,sys+" DOWN","l");
            }
        }
        for (int i = 0; i < pts_vect_WgtD.size(); i++) {
            WCPoint ref_pt_d = pts_vect_WgtD.at(i);
            TGraph* ref_pt_gr_d = new TGraph(1);
            std::cout << "DOWN VAL x:" << ref_pt_d.getStrength(wc_name) << " y: " << ref_pt_d.wgt << std::endl;
            ref_pt_gr_d->SetPoint(0,ref_pt_d.getStrength(wc_name),ref_pt_d.wgt);
            ref_pt_gr_d->SetMarkerStyle(4);
            ref_pt_gr_d->SetMarkerColor(d_clr);
            if (draw_updown and draw_points) {
                ref_pt_gr_d->Draw("P");
            }
        }

        if (include_ratio){
            pad2->cd();               // pad2 becomes the current pad
            TH1D* nom_ratio_hist;
            nom_ratio_hist = divideTF1s("n",fit,fit,x_min,x_max);
            nom_ratio_hist->SetLineColor(nom_clr);
            nom_ratio_hist->SetLineWidth(2); 
            if (wc_name == "ctG"){
                nom_ratio_hist->GetXaxis()->SetNdivisions(6,2,0,true);
            }
            nom_ratio_hist->GetYaxis()->SetNdivisions(5,2,1,true);
            nom_ratio_hist->GetYaxis()->SetTitleOffset(0.50);
            nom_ratio_hist->GetXaxis()->SetTitleSize(0.16);
            nom_ratio_hist->GetYaxis()->SetTitleSize(0.18);
            nom_ratio_hist->GetXaxis()->SetLabelSize(0.15);
            nom_ratio_hist->GetYaxis()->SetLabelSize(0.10);
            nom_ratio_hist->GetXaxis()->SetLabelOffset(0.015);
            nom_ratio_hist->GetYaxis()->SetLabelOffset(0.015);
            nom_ratio_hist->GetXaxis()->SetTickLength(0.09);
            nom_ratio_hist->GetYaxis()->SetTitle("Ratio");
            nom_ratio_hist->GetXaxis()->SetTitle(x_axis_name);
            if (sample_idx==0){
                nom_ratio_hist->Draw();
            } else {
                nom_ratio_hist->Draw("same");
            }
            TH1D* u_ratio_hist;
            u_ratio_hist = divideTF1s("u",fit_u,fit,x_min,x_max);
            u_ratio_hist->SetLineColor(u_clr);
            u_ratio_hist->SetLineWidth(2);
            TH1D* d_ratio_hist;
            d_ratio_hist = divideTF1s("d",fit_d,fit,x_min,x_max);
            d_ratio_hist->SetLineColor(d_clr);
            d_ratio_hist->SetLineWidth(2);
            if (draw_updown){
                u_ratio_hist->Draw("same");
                d_ratio_hist->Draw("same");
            }

            float ratio_y_lim = max( abs(u_ratio_hist->GetBinContent(1)-1), abs(d_ratio_hist->GetBinContent(1)-1) );
            //u_ratio_hist->SetMaximum(1.05*(1+ratio_y_lim));
            //u_ratio_hist->SetMinimum(0.95*(1-ratio_y_lim));
            nom_ratio_hist->SetMaximum(1.025); // Hard code (good for cpt)
            nom_ratio_hist->SetMinimum(0.975); // Hard code (good for cpt)
            if (draw_errorband){
                auto ratio_band = getErrorBandFromTH1s(nom_ratio_hist,u_ratio_hist,d_ratio_hist);
                ratio_band->SetFillStyle(3002);
                ratio_band->SetFillColor(fit->GetLineColor());
                ratio_band->Draw("3");
            }
            pad1->cd();
        }

        if (draw_errorband){
            const Int_t n = 100;
            double xcoord = 0;
            double step_size = (x_max-x_min)/float(n-1);
            Double_t x[n], y[n], exl[n], eyl[n], exh[n], eyh[n];
            for (int i=0; i<n; i++){
                xcoord = x_min + i*step_size;
                x[i] = xcoord;
                y[i] = fit->Eval(xcoord);
                // The errror band only makes sense in the same line is either always above or always below nominal
                // If this is not the case, these if statements do not cover that situation, we'll just get some negative error bars
                // Probably we should not be trying to draw error bands between u and d fits if they are crossing nominal anyway
                if (fit_u->Eval(xcoord) > fit->Eval(xcoord)){
                    eyh[i] = fit_u->Eval(xcoord) - fit->Eval(xcoord);
                    eyl[i] = fit->Eval(xcoord) - fit_d->Eval(xcoord);
                } else {
                    eyh[i] = fit_d->Eval(xcoord) - fit->Eval(xcoord);
                    eyl[i] = fit->Eval(xcoord) - fit_u->Eval(xcoord);
                }
                exl[i] = 0;
                exh[i] = 0;
            }
            auto error_band = new TGraphAsymmErrors(n,x,y,exl,exh,eyl,eyh);
            error_band->SetFillStyle(3002);
            error_band->SetFillColor(fit->GetLineColor());
            error_band->Draw("3");
        }

        float max_val = std::max(max_y_u,max_y_d);
        float min_val = std::min(min_y_u,min_y_d);
        //fit->SetMaximum(max_val);
        //fit->SetMinimum(min_val*.95); // Hard code (good for cpt)

    } // This closes the loop over samples

    // Draw, save, delete
    legend->Draw();
    c1->Print(save_name+".png","png");
    c1->Print(save_name+".pdf","pdf");
    delete legend;
    delete c1;

}

// Very dependent on naming schemes, retruns a dictionary that maps runs to wc values
// Right now this is set up to only handel the 5 cpt axis scan files or the 7 ctG axis scan files!!!
std::map<string,string> make_wcpt_run_map(TString wc_name) {
    std::map<string,string> ref_pts_dict;
    TString wcname = wc_name;
    float range_max; // Largest WC val in the scan (smallest should be negative of this) 
    float npts;      // Number of files in the scan
    int   run;      // Run number of first file in the scan
    std::cout << wc_name << std::endl;
    if (wc_name=="cpt"){
        range_max = 15;
        npts = 5;
        run = 0;
    } else if (wc_name=="ctG"){
        range_max = 2;
        npts = 5;
        run = 1;
    } else if (wc_name=="cpQ3"){
        range_max = 4;
        npts = 5;
        run = 0;
    } else {
        std::cout << "\nError: Do not know which runs correspond to which wc values for the given wc (" << wc_name << "), exiting.." << std::endl;
        throw std::exception();
    }
    float step = (range_max*2)/(npts-1);
    //std::cout << "STEP " << step << std::endl;
    for (float wcval=-range_max; wcval<=range_max; wcval=wcval+step){
        ref_pts_dict["run"+std::to_string(run)] = "wcpt_"+wcname+"_"+std::to_string(wcval);
        std::cout << "Run: " << run << ", wc val: " << wcval << ", Dict entry: " << "run"+std::to_string(run) << " : " << ref_pts_dict["run"+std::to_string(run)] << std::endl;
        run = run + 1;
    }
    std::cout << "ENTRY: " << ref_pts_dict["run1"] << std::endl;
    return ref_pts_dict;
}

// Returns a vector of 3 maps (for nominal, up, down), each map contains a vector of the 5 WCPoints for each systematic
std::vector<std::map<std::string,std::vector<WCPoint>>> get_WCpt_syst_maps(TString wc_name, TString run_dirs_file){

    std::vector<std::string> sys_names {"psISR","psFSR","muR","muF","muRmuF","nnpdf"};

    std::map<std::string,std::vector<WCPoint>> selection_pts_map; // Should only have key "nominal", could be a vect but want same type as selection_pts_WgtU/D_map
    std::map<std::string,std::vector<WCPoint>> selection_pts_WgtU_map;
    std::map<std::string,std::vector<WCPoint>> selection_pts_WgtD_map;

    TString sm_run;
    if (wc_name=="cpt"){
        sm_run = "run2";
    } else if (wc_name=="ctG"){
        sm_run = "run3";
    } else if (wc_name=="cpQ3"){
        sm_run = "run2";
    } else {
        std::cout << "\nError: Do not know which runs correspond to which wc values for the given wc (" << wc_name << "), exiting.." << std::endl;
        throw std::exception();
    }
    float sm_norm_NOM = 0;
    std::map<std::string,float> sm_norm_UP;
    std::map<std::string,float> sm_norm_DOWN;
    for (auto sys: sys_names){
        sm_norm_UP[sys] = 0;
        sm_norm_DOWN[sys] = 0;
    }

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

    std::map<string,string> ref_pts_dict = make_wcpt_run_map(wc_name); // Very bad, depends very much on naming scheme!)

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
            std::cout << "qCut variation: nominal. Group tag: " << grp_tag << std::endl;
        } else if (grp_tag.find("qCut25") != std::string::npos){
            f_type = "qCut_up";
            std::cout << "qCut variation: up. Group tag: " << grp_tag << std::endl;
        } else if (grp_tag.find("qCut15") != std::string::npos){
            f_type = "qCut_down";
            std::cout << "qCut variation: down. Group tag: " << grp_tag << std::endl;
        } else {
            f_type = "qCut_nom";
            std::cout << "qCut variation: NO qCut in file name, so assuming nominal. Group tag: " << grp_tag << std::endl;
        }

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
            //last_entry = 100000;
        }
        //last_entry = 100; // For testing
        std::cout << "Last_entry: " << last_entry << std::endl;
        int first_entry = 0;

        WCFit* wcFit_intree = 0;
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
        wcpt_systs_U_map["qCut"] = WCPoint(pt_str);
        wcpt_systs_D_map["qCut"] = WCPoint(pt_str);

        // Event loop:
        int selection_events = 0;
        for (int i = first_entry; i < last_entry; i++) {
            chain.GetEntry(i);

            // Define syst dict:
            std::map<std::string,std::map<std::string,double>> sys_map {
                {"psISR"  , {{"u",psISRweightUp_intree}  , {"d",psISRweightDown_intree}}},
                {"psFSR"  , {{"u",psFSRweightUp_intree}  , {"d",psFSRweightDown_intree}}},
                {"muR"    , {{"u",muRWeightUp_intree}    , {"d",muRWeightDown_intree}}},
                {"muF"    , {{"u",muFWeightUp_intree}    , {"d",muFWeightDown_intree}}},
                {"muRmuF" , {{"u",muRmuFWeightUp_intree} , {"d",muRmuFWeightDown_intree}}},
                {"nnpdf"  , {{"u",nnpdfWeightUp_intree}  , {"d",nnpdfWeightDown_intree}}}
            };

            //if (genLep_pt3_intree > 10 and genJet_pt4_intree > 30){
            if (true){ // Comment out selection cuts

                selection_events = selection_events + 1;
                if (f_type == "qCut_nom"){
                    selection_pt.wgt += originalXWGTUP_intree;
                    for (auto it = sys_map.begin(); it != sys_map.end(); it++){
                        std::string sys = it->first;
                        wcpt_systs_U_map[sys].wgt += originalXWGTUP_intree*sys_map[sys]["u"];
                        wcpt_systs_D_map[sys].wgt += originalXWGTUP_intree*sys_map[sys]["d"];
                    }
                } else if (f_type == "qCut_up"){
                    wcpt_systs_U_map["qCut"].wgt += originalXWGTUP_intree;
                } else if (f_type == "qCut_down"){
                    wcpt_systs_D_map["qCut"].wgt += originalXWGTUP_intree;
                }
            }
        }
        std::cout << "\n Selected events over total: " << selection_events << "/" << last_entry << "->" << (float)selection_events/last_entry << "\n" << std::endl;

        // Normalize and add to map
        if (f_type == "qCut_nom"){
            // Find the nominal SM value to normalize to
            if (run_label == sm_run){ 
                sm_norm_NOM = selection_pt.wgt;
                std::cout << "norm for nominal: " << sm_norm_NOM << " run : " << run_label << "grp tag: " << grp_tag << "\n" << std::endl;
            } else if (sm_norm_NOM == 0){
                std::cout << "\nError: No SM value to normalize to! Check order of files in input file (SM run needs to be first).\n" << std::endl;
                throw std::exception();
            }

            std::cout << "\n Weight of selection_pts_map[nominal] before norm!: " << selection_pt.wgt << "\n" << std::endl;

            // Normalize nominal to SM
            selection_pt.scale(1/sm_norm_NOM);
            selection_pts_map["nominal"].push_back(selection_pt);

            std::cout << "\n Weight of selection_pts_map[nominal] after norm!: " << selection_pt.wgt << "\n" << std::endl;

            // Loop over systematics, normalize and add to map
            for (auto sys: sys_names){

                // Find the up and down SM value to normalize to
                if (run_label == sm_run){ 
                    sm_norm_UP[sys] = wcpt_systs_U_map[sys].wgt;
                    sm_norm_DOWN[sys] = wcpt_systs_D_map[sys].wgt;
                    //std::cout << "\nnorm for up and down!!!, run : " << run_label << "\n" << std::endl;
                } else if ( sm_norm_UP[sys]==0 or sm_norm_DOWN[sys]==0){
                    std::cout << "\nError: No SM value to normalize to! Check order of files in input file (SM run needs to be first).\n" << std::endl;
                    throw std::exception();
                }

                // Normalize
                wcpt_systs_U_map[sys].scale(1/sm_norm_UP[sys]);   // Norm to up to up at 0
                wcpt_systs_D_map[sys].scale(1/sm_norm_DOWN[sys]); // Norm down to down at 0
                //wcpt_systs_U_map[sys].scale(1/sm_norm_NOM); // Norm up to nominal at 0
                //wcpt_systs_D_map[sys].scale(1/sm_norm_NOM); // Norm down to nominal at 0

                // Fill the vector of WC points
                selection_pts_WgtU_map[sys].push_back(wcpt_systs_U_map[sys]);
                selection_pts_WgtD_map[sys].push_back(wcpt_systs_D_map[sys]);
            }
        } 

        // Normalize and add qCut to map
        else if (f_type == "qCut_up"){
            if (run_label == sm_run){
                sm_norm_UP["qCut"] = wcpt_systs_U_map["qCut"].wgt;
            } else if (sm_norm_UP["qCut"]==0){
                std::cout << "\nError: No SM value to normalize to for qCut up! Check order of files in input file (SM run needs to be first).\n" << std::endl;
                throw std::exception();
            }
            //std::cout << " sm_norm_UP[qCut] " << sm_norm_UP["qCut"] << std::endl;
            wcpt_systs_U_map["qCut"].scale(1/sm_norm_UP["qCut"]); // Norm up to up at 0
            //wcpt_systs_U_map["qCut"].scale(1/sm_norm_NOM); // Norm up to nominal at 0
            selection_pts_WgtU_map["qCut"].push_back(wcpt_systs_U_map["qCut"]);
        } 
        else if (f_type == "qCut_down"){
            if (run_label == sm_run){
                sm_norm_DOWN["qCut"] = wcpt_systs_D_map["qCut"].wgt;
            } else if (sm_norm_DOWN["qCut"]==0){
                std::cout << "\nError: No SM value to normalize to for qCut down! Check order of files in input file (SM run needs to be first).\n" << std::endl;
                throw std::exception();
            }
            //std::cout << " sm_norm_DOWN[qCut] " << sm_norm_DOWN["qCut"] << std::endl;
            wcpt_systs_D_map["qCut"].scale(1/sm_norm_DOWN["qCut"]); // Norm down to down at 0
            //wcpt_systs_D_map["qCut"].scale(1/sm_norm_NOM); // Norm down to nominal at 0
            selection_pts_WgtD_map["qCut"].push_back(wcpt_systs_D_map["qCut"]);
        }

    }
    std::vector<std::map<std::string,std::vector<WCPoint>>> return_vect = {selection_pts_map,selection_pts_WgtU_map,selection_pts_WgtD_map};
    return return_vect;
}

void systWCdependenceCheck(TString proc_name, TString wc_name, TString run_dirs_file, TString run_dirs_qCuts_file, TString run_dirs_0p_file) {
    gStyle->SetOptStat(0);
    std::cout << "\nThe WC we scan over in these files: " << wc_name << "\n" << std::endl;

    // Fill the dictionaries of WC points for each syst, Note for the systs_pts_vect: 0 is nominal, 1 is up, 2 is down
    std::vector<std::string> sys_names {"psISR","psFSR","muR","muF","muRmuF","nnpdf","qCut"};
    std::vector<std::map<std::string,std::vector<WCPoint>>> systs_pts_vect;
    std::vector<std::map<std::string,std::vector<WCPoint>>> systs_pts_vect_p1p;
    std::vector<std::map<std::string,std::vector<WCPoint>>> systs_pts_vect_p0p;
    //systs_pts_vect = get_WCpt_syst_maps(wc_name,run_dirs_file);         // plus 1p (No qCut variations)
    systs_pts_vect_p1p = get_WCpt_syst_maps(wc_name,run_dirs_qCuts_file); // plus 1p (With qCut variations)
    systs_pts_vect_p0p = get_WCpt_syst_maps(wc_name,run_dirs_0p_file);    // plus 0p file

    // Make plots for all of the systematics individually
    //sys_names.push_back("qCut"); // If including qCut
    for(auto sys : sys_names){
        makePlot(string(wc_name),sys,{systs_pts_vect_p1p.at(0)["nominal"]},{systs_pts_vect_p1p.at(1)[sys]},{systs_pts_vect_p1p.at(2)[sys]});
        if (sys != "qCut"){
            std::cout << "Not quct: " << sys << std::endl;
            makePlot(string(wc_name),sys+"_0p",{systs_pts_vect_p0p.at(0)["nominal"]},{systs_pts_vect_p0p.at(1)[sys]},{systs_pts_vect_p0p.at(2)[sys]});
        }
    }

    // Get the envelope of muR and muF
    std::map<std::string,std::vector<WCPoint>> muRmuFenvelope_map_p1p;
    std::map<std::string,std::vector<WCPoint>> muRmuFenvelope_map_p0p;
    muRmuFenvelope_map_p1p = muRmuFenvelope(string(wc_name),systs_pts_vect_p1p.at(0)["nominal"],systs_pts_vect_p1p.at(1),systs_pts_vect_p1p.at(2));
    muRmuFenvelope_map_p0p = muRmuFenvelope(string(wc_name),systs_pts_vect_p0p.at(0)["nominal"],systs_pts_vect_p0p.at(1),systs_pts_vect_p0p.at(2));

    makePlot(string(wc_name),"muRmuF_evnelope",{systs_pts_vect_p1p.at(0)["nominal"]},{muRmuFenvelope_map_p1p["inc"]},{muRmuFenvelope_map_p1p["dec"]});
    systs_pts_vect_p1p.at(1)["muRmuFenv"] = muRmuFenvelope_map_p1p["inc"]; // Add to map, though up and down do not hold meaning here
    systs_pts_vect_p1p.at(2)["muRmuFenv"] = muRmuFenvelope_map_p1p["dec"]; // Add to map, though up and down do not hold meaning here
    makePlot(string(wc_name),"muRmuF_evnelope_0p",{systs_pts_vect_p0p.at(0)["nominal"]},{muRmuFenvelope_map_p0p["inc"]},{muRmuFenvelope_map_p0p["dec"]});
    systs_pts_vect_p0p.at(1)["muRmuFenv"] = muRmuFenvelope_map_p0p["inc"]; // Add to map, though up and down do not hold meaning here
    systs_pts_vect_p0p.at(2)["muRmuFenv"] = muRmuFenvelope_map_p0p["dec"]; // Add to map, though up and down do not hold meaning here
    
    // Add the systematics in quadrature
    std::vector<std::string> s_to_add_in_quad {"psISR","psFSR","nnpdf","qCut","muRmuFenv"}; // with qCut
    std::vector<std::string> s_to_add_in_quad_noqCut {"psISR","psFSR","nnpdf","muRmuFenv"}; // No qCut
    std::map<std::string,std::vector<WCPoint>> quad_sum_map_p1p;
    std::map<std::string,std::vector<WCPoint>> quad_sum_map_p0p;
    quad_sum_map_p1p = getQuadSums(string(wc_name),s_to_add_in_quad,systs_pts_vect_p1p.at(0),systs_pts_vect_p1p.at(1),systs_pts_vect_p1p.at(2));
    quad_sum_map_p0p = getQuadSums(string(wc_name),s_to_add_in_quad_noqCut,systs_pts_vect_p0p.at(0),systs_pts_vect_p0p.at(1),systs_pts_vect_p0p.at(2));
    makePlot(string(wc_name),"quad_sum",{systs_pts_vect_p1p.at(0)["nominal"]},{quad_sum_map_p1p["inc"]},{quad_sum_map_p1p["dec"]},{proc_name+" 0+1p "});
    makePlot(string(wc_name),"quad_sum_0p",{systs_pts_vect_p0p.at(0)["nominal"]},{quad_sum_map_p0p["inc"]},{quad_sum_map_p0p["dec"]},{proc_name+" 0+1p "});

    makePlot(
        string(wc_name),
        "quad_sum_0p_and_0p1p",
        {systs_pts_vect_p1p.at(0)["nominal"],systs_pts_vect_p0p.at(0)["nominal"]},
        {quad_sum_map_p1p["inc"],quad_sum_map_p0p["inc"]},
        {quad_sum_map_p1p["dec"],quad_sum_map_p0p["dec"]},
        {proc_name+" 0+1p ",proc_name+" 0p "} // Optional argument, each tag should correspond to the samples passed
    );

    // TEST print frac uncty
    for(auto sys : s_to_add_in_quad){
        printFracUncty(string(wc_name) , sys , systs_pts_vect_p1p.at(0)["nominal"] , systs_pts_vect_p1p.at(1)[sys], systs_pts_vect_p1p.at(2)[sys]);
    }
    printFracUncty(string(wc_name) , "quad sum" , systs_pts_vect_p1p.at(0)["nominal"] , quad_sum_map_p1p["inc"] , quad_sum_map_p1p["dec"]);


}

//{
//    //TString sample = "privateTTW";
//    //TString sample = "privateTTW-NoDim6";
//    TString sample = "privateTTW-NoDim6-1Jet";
//    //TString sample = "privateTZQ4f-NoDim6-NoSchanW";
//    //TString sample = "centralTTZ";
//    TFile* f = TFile::Open(sample+"_NoTopLeptons_output_tree.root");
//    TDirectory* td = f->GetDirectory("EFTGenReader");
//    TIter next(td->GetListOfKeys());
//    TKey* key;
//    while ((key = (TKey*)next())) {
//        TString cmp = "summaryTree";
//        TString s = key->GetName();
//        if (s == cmp) continue;
//        s += "_" + sample + ".png";
//        std::cout << s << std::endl;
//        TH1D* h = (TH1D*)td->Get(key->GetName());
//
//        TCanvas c1("c1","",1280,720);
//        c1.cd();
//        h->Draw();
//        c1.Print(s,"png");
//    }
//    f->Close();
//}

#include "EFTGenReader/EFTHelperUtilities/interface/TH1EFT.h"
#include "EFTGenReader/EFTHelperUtilities/interface/WCPoint.h"

TCanvas* findCanvas(TString search,std::vector<TCanvas*> canvs) {
    //TCanvas* canv = nullptr;
    for (auto c: canvs) {
        if (c->GetTitle() == search) {
            return c;
            //canv = c;
            //break;
        }
    }
    return nullptr;
}

int findCanvasIndex(TString search,std::vector<TCanvas*> canvs) {
    int idx = 0;
    for (auto c: canvs) {
        if (c->GetTitle() == search) {
            return idx;
        }
        idx++;
    }
    return -1;
}

// Returns a histogram that is the ratio of h1/h2
template<typename T>
TH1D* makeRatioHistogram(TString name,T* h1,T* h2) {
    if (h1->GetXaxis()->GetNbins() != h2->GetXaxis()->GetNbins()) {
        std::cout << "[Error] makeRatioHistogram() - bin mismatch between ratio of hists" << std::endl;
        throw;
    }
    TAxis* xaxis = h1->GetXaxis();
    int bins = xaxis->GetNbins();
    Double_t low = xaxis->GetXmin();
    Double_t high = xaxis->GetXmax();

    TAxis* yaxis = h1->GetYaxis();
    Double_t yaxis_sz = yaxis->GetLabelSize();

    // std::cout << "Ratio Name: " << name << std::endl;
    TH1D* h_ratio = new TH1D(name,"",bins,low,high);    // Make sure the title is empty!
    h_ratio->GetYaxis()->SetLabelSize(yaxis_sz*1.5);
    h_ratio->GetXaxis()->SetLabelSize(yaxis_sz*2.5);    // Note: If this is done after alphanumeric, it does nothing
    for (int i = 1; i <= bins; i++) {
        if (xaxis->IsAlphanumeric()) {
            TString bin_label = xaxis->GetBinLabel(i);
            h_ratio->GetXaxis()->SetBinLabel(i,bin_label);
        }
        double val1 = h1->GetBinContent(i);
        double err1 = h1->GetBinError(i);

        double val2 = h2->GetBinContent(i);
        double err2 = h2->GetBinError(i);

        double ratio;
        double ratio_err;
        if (val2 != 0) {
            ratio = abs(val1 / val2);
        } else if (val1 == 0 && val2 == 0) {
            ratio = 1.0;
        } else {
            ratio = -1;
        }
        ratio_err = 0.0;    // Ignore the error for now
        h_ratio->SetBinContent(i,ratio);
        h_ratio->SetBinError(i,ratio_err);
    }
    //h_ratio->GetYaxis()->SetRangeUser(0.0,2.3);
    return h_ratio;
}


//void makeEFTGenPlots(TString output_fname,std::vector<TString> input_fnames) {
void makeEFTGenPlots(std::vector<TString> input_fnames, TString wc_string) {

    // Set up the type of plots we want to make
    bool only_njets = false;
    bool only_jetPt_lepPt = false;
    bool only_SM = false;
    bool include_ratio = false;
    bool draw_sm = true; // Only really makes sense if we are only drawing one file currently
    std::string norm_type = "unit_norm"; //"SM_rel_norm";
    std::string plot_type = "0p_vs_1p_comp";

    std::vector<TFile*> files;
    TH1::SetDefaultSumw2();

    std::vector<int> clrs {kBlue,kRed,kGreen};
    // TLegend parameters
    double left,right,top,bottom,scale_factor,minimum;
    if (include_ratio) {
        //left         = 0.75-0.13; // SetRightMargin is 0.13
        //right        = 0.95-0.13; // SetRightMargin is 0.13
        left         = 0.74;
        right        = 0.94;
        scale_factor = 0.08;
        top          = 0.89;
    } else {
        left         = 0.81;
        right        = 0.98;
        scale_factor = 0.05;
        top          = 0.9;
    }
    minimum      = 0.1;

    //gStyle->SetPalette(kCLRS,palette);
    //if (input_fnames.size() <= 3) {
    //    gStyle->SetPalette(kCLRS2,pal2);
    //}
    gStyle->SetOptStat(0);
    gStyle->SetPadRightMargin(0.2);

    //std::cout << "Output File: " << output_fname << std::endl;
    for (auto s: input_fnames) {
        std::cout << "  Input File: " << s << std::endl;
        TFile* f = TFile::Open(s);
        if (f) {
            //f->Print();
            files.push_back(f);
        }
    }

    if (!files.size()) {
        return;
    }

    std::cout << "" << std::endl;

    std::vector<TCanvas*> canvs;
    std::vector<TLegend*> legs;
    std::map<std::string,TH1EFT*> empty_hists_dict; // emptytest
    std::map<std::string,double> max_yvals_dict; // emptytest
    //WCPoint* smpt = new WCPoint();
    WCPoint* wc_pt = new WCPoint(wc_string.Data());
    WCPoint* sm_pt = new WCPoint("smpt");
    int filecounter = 0;


    // Loop over files and histograms to  find largest y values:
    std::vector<TFile*> files1;
    for (auto s: input_fnames) {
        TFile* f1 = TFile::Open(s);
        if (f1) {
            files1.push_back(f1);
        }
    }
    /*
    std::map<string,float> maxYvals_dict;
    std::map<std::string,TH1D*> histSM_dict;
    for (auto f1: files1) {
        TDirectory* td1;
        td1 = f1->GetDirectory("EFTGenReader");
        if (td1 == 0){
            td1 = f1->GetDirectory("EFTGenHistsWithCuts"); // In case we are using files produced by the  EFTGenHistsWithCuts analyzer
        }
        TKey* key1;
        TIter next(td1->GetListOfKeys());
        while ((key1 = (TKey*)next())) {
            TString s1 = key1->GetName();
            TH1EFT* h_eventsumEFT1 = (TH1EFT*)td1->Get("h_eventsumEFT");
            double eventsum_SM1 = h_eventsumEFT1->GetBinFit(1).evalPoint(sm_pt);
            TH1EFT* h1 = (TH1EFT*)td1->Get(key1->GetName());
            //if (s1.Index("EFT") != -1) {
            if (s1.Index("SM") == -1 and s1.Index("summaryTree") == -1) { // Not all EFT hists have EFT in their name, but we assume all SM hists have SM in the name
                for (Int_t bin_idx = 0; bin_idx <= h1->GetNbinsX()+1; bin_idx++) {
                    double wcfit_bin_val = h1->GetBinFit(bin_idx).evalPoint(wc_pt);
                    h1->SetBinContent(bin_idx,wcfit_bin_val);
                }
            }
            // SM rel
            if (TString(norm_type).Contains("SM_rel")) {
                if (s1 != "h_SMwgt_norm" && s1.Contains("EFT")) {
                    auto tmps = s1;
                    std::cout << tmps << " " << tmps.ReplaceAll("EFT","SM") << std::endl;
                    auto hsm = (TH1D*)f1->Get("EFTGenReader/"+tmps);
                    std::cout << "EFTGenReader/"+tmps << std::endl;
                    if(hsm != nullptr) {
                    hsm->SetDirectory(0);
                    histSM_dict[tmps.Data()]=hsm;
                    }
                }
            }
            // SM norm
            if (norm_type == "SM_norm"){
                if (s1 != "h_SMwgt_norm") {
                    h1->Scale(1.0/eventsum_SM1);
                }
            }
            // Unit norm
            if (norm_type == "unit_norm") {
                if (s1 != "h_SMwgt_norm") {
                    Int_t nbins = h1->GetNbinsX();
                    Double_t intg = h1->Integral(0,nbins+1);
                    if (intg > 1.0) {
                        h1->Scale(1./intg);
                        h1->Scale(1./h1->GetBinWidth(1)); // Scale by bin width
                    }
                }
            }

            // Fill the dictionary:
            if (maxYvals_dict.count(string(s1)) == 0) {
                maxYvals_dict[string(s1)] = h1->GetMaximum();
            } else if (maxYvals_dict[string(s1)] < h1->GetMaximum()){
                //std::cout << "New largest value for " << s1 << ": " << h1->GetMaximum() << " (Old largest val: " << maxYvals_dict[string(s1)] << ")" << std::endl;
                maxYvals_dict[string(s1)] = h1->GetMaximum();
              }
        }
    }
    for (auto f1: files1) {
        f1->Close();
    }
    // End of loop over files to find largest y value
    */


    std::map<std::string,std::vector<TH1D*>> hist_dict;
    int clr_idx = 0;

    for (auto f: files) {
        f->GetList()->Clear(); // ???
        f->Print();

        TString fname = f->GetName();
        //cout << "name of the file !!! " << fname << std::endl;
        TString sub_str;
        Ssiz_t idx = fname.First('/');
        //cout << "idx !!! " << idx << std::endl;
        int idx_begin = 0;
        int idx_end = 0;
        if (idx != TString::kNPOS) {
            sub_str = fname(idx+1,fname.Length());
        }

        /*
        //TString marker = "output_";
        TString marker = "HanV4Model";
        //TString marker = "0p_";
        idx_begin = sub_str.Index(marker)+marker.Length();
        //idx_end = sub_str.Index("fromMAOD.root");
        idx_end = sub_str.Index(".root");
        sub_str = sub_str(idx_begin,idx_end-idx_begin);
        //cout << " this is the sub str now !!!!!!!!!! " << sub_str << std::endl;

        // Hard code for  making no jets vs jets plots!
        if (plot_type == "0p_vs_1p_comp") {
            if (sub_str.Index("NoJets") != -1){
                std::cout << "the no jets sub str! " << sub_str << std::endl;
                sub_str = "ttH 0 partons";
            } else {
                std::cout << "the plus jets sub str! " << sub_str << std::endl;
                sub_str = "ttH 0+1 partons";
            }
        }
        */

        TString marker = "output_";
        //TString marker = "HanV4Model";
        //TString marker = "0p_";
        //idx_begin = sub_str.Index(marker)+marker.Length();
        //idx_end = sub_str.Index("fromMAOD.root");
        //idx_end = sub_str.Index(".root");
        //sub_str = sub_str(idx_begin,idx_end-idx_begin);
        //cout << "\ntThis is the sub str now !!!!!!!!!!\n" << sub_str << std::endl;

        // This legend string stuff is very specific to the names of the files and is just leftover 
        // from the last thing I used this code for, i.e. this is very ad hoc and will probably not 
        // work great for situations besides the specific one I was using it for (sorry)
        //sub_str = "ttH";
        std::string leg_str;
        if (sub_str.Index("ttH") != -1){
            leg_str = "ttH ";
        } else if (sub_str.Index("ttW") != -1){
            leg_str = "ttW ";
        } else if (sub_str.Index("ttWZ") != -1){
            leg_str = "ttZ ";
        } else if (sub_str.Index("ttlnu") != -1){
            leg_str = "ttlnu ";
        } else if (sub_str.Index("ttll") != -1){
            leg_str = "ttll ";
        }
        if (plot_type == "0p_vs_1p_comp") {
            if (sub_str.Index("NoJets") != -1){
                //std::cout << "the no jets sub str! " << sub_str << std::endl;
                sub_str = leg_str+"0 partons";
            } else if (sub_str.Index("NLO") != -1){
                //std::cout << "the NLO sub str! " << sub_str << std::endl;
                sub_str = leg_str+"NLO";
            } else if (sub_str.Index("central") != -1){
                //std::cout << "The NLO sub str! " << sub_str << std::endl;
                sub_str = leg_str+"central";
            } else if (sub_str.Index("Jet") != -1) {
                //std::cout << "the plus jets sub str! " << sub_str << std::endl;
                sub_str = leg_str+"0+1 partons";
            } else {
                //std::cout << "the 0 jets sub str also! " << sub_str << std::endl;
                sub_str = leg_str+"0 partons";
            }
        }

        TDirectory* td;
        td = f->GetDirectory("EFTGenReader");
        if (td == 0){
            td = f->GetDirectory("EFTGenHistsWithCuts"); // In case we are using files produced by the  EFTGenHistsWithCuts analyzer
        }

        TH1EFT* h_eventsumEFT = (TH1EFT*)td->Get("h_eventsumEFT");
        double eventsum_SM = h_eventsumEFT->GetBinFit(1).evalPoint(sm_pt);

        TIter next(td->GetListOfKeys());
        TKey* key; TCanvas* c; TLegend* leg;
        while ((key = (TKey*)next())) {
            TString cmp = "summaryTree";
            TString s = key->GetName();
            if (s == cmp) continue;

            if(only_SM) { // For not plotting EFT plots
                if (s.Index("EFT") != -1) {
                    continue;
                }
            }
            if(only_njets){ // For only plotting njets plots
                if (s.Index("njets") == -1 and s.Index("nJets") == -1) { // Apparently I am inconsistent in njets naming convention
                    continue;
                }
            }
            if(only_jetPt_lepPt){
                if (s.Index("lep1_pt")==-1 and s.Index("lep2_pt")==-1 and s.Index("jet1_pt")==-1 and s.Index("jet2_pt")==-1){
                    continue;
                }
            }

            //TH1D* h = (TH1D*)td->Get(key->GetName());
            TH1EFT* h = (TH1EFT*)td->Get(key->GetName());
            h->SetMarkerStyle(kFullCircle);
            h->SetMarkerSize(0.25);
            h->SetOption("E");
            h->SetMarkerColor(clrs.at(clr_idx));
            if (only_njets){
                h->SetTitle("");
            }

            TH1EFT* h_sm_rwgt = (TH1EFT*)h->Clone(s+"_sm_rwgt");
            h_sm_rwgt->SetMarkerStyle(kFullCircle);
            h_sm_rwgt->SetMarkerSize(0.25);
            h_sm_rwgt->SetOption("E");
            h_sm_rwgt->SetMarkerColor(kBlack);
            //h_sm_rwgt->SetTitle("");

            //std::cout << "h bins: " << h->GetNbinsX() << std::endl;
            //std::cout << "h fit size: " << h->hist_fits.size() << std::endl;
            //std::cout << "h sm bins: " << h_sm_rwgt->GetNbinsX() << std::endl;
            //std::cout << "h sm fit size: " << h_sm_rwgt->hist_fits.size() << std::endl;
            //std::cout << "" << std::endl;

            //Int_t nbins = h->GetNbinsX();
            //Double_t intg = h->Integral(0,nbins+1);
            //if (intg > 1.0) {
            //    h->Scale(1./intg);
            //}

            //// Rebin the pt distrubutions
            //if (s.Index("_pt") != -1) {
            //    //cout << "this is a pt plot!" << s << "\n";
            //	h->Rebin(10);
            //}

            c = findCanvas(key->GetName(),canvs);
            int c_idx = findCanvasIndex(key->GetName(),canvs);

            TString canv_str = (TString)f->GetName() + "-" + (TString)key->GetName();
            //std::cout << "\nthe canvas info!!! " << canv_str << "\n" << std::endl;

            bool is_new = !c;
            if (is_new) {
                //std::cout << "is new !!!!" << std::endl;
                //c = new TCanvas(canv_str,key->GetName(),1280,720);
                c = new TCanvas(canv_str,key->GetName(),1100,800);

                if (include_ratio) {
                    Float_t small = .04;
                    const float padding = 1e-5;
                    const float ydiv = 0.3;
                    c->Divide(1,2,small,small);
                    c->GetPad(1)->SetPad(padding,ydiv+padding,1-padding,1-padding);
                    c->GetPad(1)->SetLeftMargin(.05);
                    c->GetPad(1)->SetRightMargin(.05);
                    c->GetPad(1)->SetBottomMargin(.3);
                    c->GetPad(1)->Modified();

                    c->GetPad(2)->SetLeftMargin(.05);
                    c->GetPad(2)->SetRightMargin(.05);
                    c->GetPad(2)->SetBottomMargin(.3);
                    c->GetPad(2)->SetPad(padding,padding,1-padding,ydiv-padding);
                    c->GetPad(2)->Modified();

                    c->cd(1);
                    gPad->SetBottomMargin(small);
                    gPad->Modified();

                    c->cd(2);
                    gPad->SetTopMargin(small);
                    gPad->SetTickx();
                    gPad->Modified();
                }

                canvs.push_back(c);
                bottom = std::max(top - scale_factor*(files.size()+1),minimum);
                leg = new TLegend(left,top,right,bottom);
                legs.push_back(leg);

                // emptytest This hist will be empty, we'll draw it first and use it for y axis range
                TH1EFT* h_empty = (TH1EFT*)h->Clone(key->GetName());
                for (Int_t bin_idx = 0; bin_idx <= h_empty->GetNbinsX()+1; bin_idx++) {
                    h_empty->SetBinContent(bin_idx,0);
                }
                empty_hists_dict[key->GetName()] = h_empty;
                h_empty->Draw();
                max_yvals_dict[key->GetName()] = 0;

            } else {
                leg = legs.at(c_idx);
            }

            // Log log scale for wgt hist
            TPad* canv_pad;
            if (include_ratio) {
                canv_pad = (TPad*)c->cd(1);
            } else {
                canv_pad = (TPad*)c->cd();
            }
            if (s == "h_SMwgt_norm") {
                canv_pad->SetLogy(1);
                canv_pad->SetLogx(1);
            }

            //if (s.Index("EFT") != -1) {
            if (s.Index("SM") == -1 and s.Index("summaryTree") == -1) { // Not all EFT hists have EFT in their name, but we assume all SM hists have SM in the name
                //std::cout << h->GetName() << std::endl;
                for (Int_t bin_idx = 0; bin_idx <= h->GetNbinsX()+1; bin_idx++) {
                    double wcfit_bin_val = h->GetBinFit(bin_idx).evalPoint(wc_pt);
                    double wcfit_bin_err = h->GetBinFit(bin_idx).evalPointError(wc_pt);
                    //std::cout << "this is the wcfit_bin_val: " << wcfit_bin_val << std::endl;
                    h->SetBinContent(bin_idx,wcfit_bin_val);
                    h->SetBinError(bin_idx,wcfit_bin_err);
                }
                for (Int_t bin_idx = 0; bin_idx <= h_sm_rwgt->GetNbinsX()+1; bin_idx++) {
                    double wcfit_bin_val_sm = h_sm_rwgt->GetBinFit(bin_idx).evalPoint(sm_pt);
                    double wcfit_bin_err_sm = h_sm_rwgt->GetBinFit(bin_idx).evalPointError(sm_pt);
                    //std::cout << "this is the wcfit_bin_val_sm: " << wcfit_bin_val_sm<< std::endl;
                    h_sm_rwgt->SetBinContent(bin_idx,wcfit_bin_val_sm);
                    h_sm_rwgt->SetBinError(bin_idx,wcfit_bin_err_sm);

                }
            }
            /*
            // SM rel
            if (TString(norm_type).Contains("SM_rel")){
                if (s.Index("EFT") != -1) {
                    auto tmps = s;
                    tmps.ReplaceAll("EFT","SM");
                    TH1D* h_sm = (TH1D*)histSM_dict[tmps.Data()];
                    double normH = 1.;
                    double normHsm = 1.;
                    if(TString(norm_type).Contains("norm")) {
                        double normHsm = 1./h_sm->Integral();
                        h->Scale(normH);
                        h_sm->Scale(normHsm);
                    }
                    for (Int_t bin_idx = 1; bin_idx <= h->GetNbinsX(); bin_idx++) {
                        double sm_bin_val = h_sm->GetBinContent(bin_idx);
                        double sm_bin_err = h_sm->GetBinContent(bin_idx);
                        double wcfit_bin_val = h->GetBinFit(bin_idx).evalPoint(wc_pt);
                        double wcfit_bin_err = h->GetBinFit(bin_idx).evalPointError(wc_pt);
                        if(bin_idx == h->GetNbinsX()) {
                            sm_bin_val += h_sm->GetBinContent(bin_idx+1);
                            sm_bin_err = sqrt(pow(sm_bin_err,2) + pow(h_sm->GetBinError(bin_idx+1),2));
                            wcfit_bin_val += h->GetBinFit(bin_idx+1).evalPoint(wc_pt);
                            wcfit_bin_err = sqrt(pow(wcfit_bin_err,2) + pow(h->GetBinFit(bin_idx+1).evalPointError(wc_pt),2));
                        }
                        if(TString(norm_type).Contains("norm")) {
                            wcfit_bin_val *= 1 / h->Integral();
                            wcfit_bin_err *= 1 / h->Integral();
                        }
                        wcfit_bin_err = wcfit_bin_val/sm_bin_val * sqrt(wcfit_bin_err/wcfit_bin_val + sm_bin_err/sm_bin_val);
                        wcfit_bin_val *= 1/h_sm->GetBinContent(bin_idx);
                        h->SetBinContent(bin_idx,wcfit_bin_val);
                        h->SetBinError(bin_idx,wcfit_bin_err);
                    }
                    if(TString(norm_type).Contains("norm")) tmps.ReplaceAll("SM","SM_norm");
                    tmps.ReplaceAll("SM","EFToSM");
                    h->SetName(tmps);
                    h->SetTitle(tmps);
                    h->Fit("pol1", "FSMEQW");
                    float chi = 0.;
                    for (Int_t bin_idx = 1; bin_idx <= h->GetNbinsX()+1; bin_idx++) {
                        double bin = h->GetBinCenter(bin_idx);
                        double bin_val = h->GetBinContent(bin_idx);
                        chi += pow(h->GetFunction("pol1")->Eval(bin) - bin_val, 2) / bin_val;
                    }
                    //std::cout << "Fit chi^2 for " << tmps << " = " << chi << std::endl;
                }
            }
            */
            // Unit norm
            if (norm_type == "unit_norm"){
                if (s != "h_SMwgt_norm") {
                    Int_t nbins = h->GetNbinsX();
                    Double_t intg = h->Integral(0,nbins+1);
                    if (intg > 1.0) {
                        h->Scale(1./intg);
                        h->Scale(1./h->GetBinWidth(1)); // Scale by bin width
                        //h->GetYaxis()->SetRangeUser(0.0,1.2*maxYvals_dict[string(s)]);
                        // SM histogram
                        Double_t intg_sm_rwgt = h_sm_rwgt->Integral(0,nbins+1); // Should this have its own "if intg_sm_rwgt>0" statement?
                        h_sm_rwgt->Scale(1./intg_sm_rwgt);
                        h_sm_rwgt->Scale(1./h_sm_rwgt->GetBinWidth(1)); // Scale by bin width
                    }
                } else {
                    std::cout << "\nNot normalizing this histogram: " << s << "\n" << std::endl;
                    //h->GetYaxis()->SetRangeUser(0.1,2*maxYvals_dict[string(s)]);
                }
            }

            // emptytest
            //std::cout << "s: " << s << " max_yvals_dict[string(s)]: " << max_yvals_dict[string(s)] << std::endl;
            if (h->GetMaximum() > max_yvals_dict[string(s)]){
                max_yvals_dict[string(s)] = h->GetMaximum();
                empty_hists_dict[string(s)]->GetYaxis()->SetRangeUser(0.0,1.2*max_yvals_dict[string(s)]);
                //std::cout << "\th->GetMaximum() is new max: " << h->GetMaximum() << std::endl;
            }
            /*
            if (is_new) {
                //h->Draw("E PLC PMC");
                h->Draw("E");
            } else {
                //h->Draw("E SAME PLC PMC");
                h->Draw("E SAME");
            }
            */
            h->Draw("E SAME"); // emptytest 
            if (draw_sm) {
                h_sm_rwgt->SetLineColor(kBlack);
                h_sm_rwgt->Draw("E SAME");
                leg->AddEntry(h_sm_rwgt,sub_str+"_SM","l");
                // Check if max val is larger than previous max y value, and set range of empty_hists_dict accordingly
                if (h_sm_rwgt->GetMaximum() > max_yvals_dict[string(s)]){
                    max_yvals_dict[string(s)] = h_sm_rwgt->GetMaximum();
                    empty_hists_dict[string(s)]->GetYaxis()->SetRangeUser(0.0,1.2*max_yvals_dict[string(s)]);
                }
            }
            //std::cout << "s: " << s << " max val: " << max_yvals_dict[string(s)] << std::endl;;

            //std::cout << "clr_idx: " << clr_idx << " clrs.at(clr_idx) " << clrs.at(clr_idx) << std::endl;
            h->SetLineColor(clrs.at(clr_idx));
            leg->AddEntry(h,sub_str,"l");
            if (include_ratio) {
                leg->SetBorderSize(0);
            }
            c->Update();
            hist_dict[string(s)].push_back(h);
        }
        clr_idx = clr_idx + 1;
        filecounter++;
    }

    // Get ratio hists and put them in a dictionary (format of hist_dict: {"hist name": [h1,h2]})
    TH1D* ratio_hist;
    std::map<std::string,std::vector<TH1D*>> ratio_hist_dict;
    //std::cout << "\n the length !!! " << hist_dict.size() << "\n" << std::endl;
    for (auto it = hist_dict.begin(); it != hist_dict.end(); it++ ){
        //std::cout << " it->first " << it->first << std:: endl;
        //std::cout << "hist_dict[it->first].size()" << hist_dict[it->first].size() << std::endl;
        for(int i=0; i<hist_dict[it->first].size(); i++){
            // Note: I added "i" to the name we pass to makeRatioHistogram to get rid of "Potential memory leak" issue, next 
            // time we want to actually include a ratio hist, should check if this causes any problems
            ratio_hist = makeRatioHistogram(it->first+to_string(i),hist_dict[it->first].at(i), hist_dict[it->first].at(0)); 
            ratio_hist_dict[string(it->first)].push_back(ratio_hist);
            // Print the hists
            for(int k=0; k<=hist_dict[it->first].at(i)->GetNbinsX(); k++){
                //std::cout << "hist " << i << " bin " << k << ": " << hist_dict[it->first].at(i)->GetBinContent(k) << std::endl;
            }
            for(int k=0; k<=ratio_hist->GetNbinsX(); k++){
                //std::cout << "h" << i << "/h0 bin " << k << ": " << ratio_hist->GetBinContent(k) << std::endl;
            }
        }
    }


    for (auto c: canvs) {

        c->cd();
        int idx = findCanvasIndex(c->GetTitle(),canvs);
        c->cd(1);
        legs.at(idx)->Draw();

        // Draw the ratio hists
        if (include_ratio) {
            c->cd(2);
            // Check the max and min ratio values to scale the y axis
            double max_ratio = 0;
            double min_ratio = 1;
            double ratio_high = 0;
            double ratio_low = 0;
            for (int i=0; i<ratio_hist_dict[c->GetTitle()].size(); i++){
                ratio_high = ratio_hist_dict[c->GetTitle()].at(i)->GetBinContent(ratio_hist_dict[c->GetTitle()].at(i)->GetMaximumBin());
                ratio_low  = ratio_hist_dict[c->GetTitle()].at(i)->GetBinContent(ratio_hist_dict[c->GetTitle()].at(i)->GetMinimumBin());
                if (ratio_high > max_ratio){
                    max_ratio = ratio_high;
                } 
                if (ratio_low < min_ratio){
                    min_ratio = ratio_low;
                }
            }
            //std::cout << "max for " << c->GetTitle() << " " << max_ratio << std::endl;
            //std::cout << "min for " << c->GetTitle() << " " << min_ratio << std::endl;

            // Plot the hists
            clr_idx = 0;
            for (int i=0; i<ratio_hist_dict[c->GetTitle()].size(); i++){
                TH1D* r_hist = ratio_hist_dict[c->GetTitle()].at(i);
                r_hist->GetYaxis()->SetNdivisions(305,true);
                if (i==0){
                    r_hist->Draw("e2");
                    //r_hist->SetMaximum(1.2*max_ratio);
                    //r_hist->SetMinimum(0.4*min_ratio);
                    // Hard code max,min: 0p vs 1p: 4,0; qCut comp ttH_ctG2: 1.75,.25
                    r_hist->SetMaximum(1.75); // Hard code max
                    r_hist->SetMinimum(.25);  // Hard code min
                    if (only_njets) {
                        r_hist->GetXaxis()->SetTitle("N jets");
                    }
                    r_hist->GetXaxis()->CenterTitle();
                    r_hist->GetXaxis()->SetTitleOffset(1.0);
                    r_hist->GetXaxis()->SetTitleSize(0.14);
                    r_hist->GetXaxis()->SetTitleFont(12);
                } else {
                    r_hist->Draw("e2 same");
                }
                //std::cout << "color idx " << clr_idx << std::endl;
                r_hist->SetLineColor(clrs.at(clr_idx));
                clr_idx = clr_idx + 1;
            }
            c->GetPad(1)->RedrawAxis();
            c->GetPad(2)->RedrawAxis();
            c->Update();
        }


        //TString s = (TString)c->GetTitle() + "_" + output_fname + ".png";
        // No need to uniquely name the images, since they will be placed in a dedicated directory
        TString s = (TString)c->GetTitle() + ".png";
        if(norm_type == "SM_rel") s.ReplaceAll("EFT","EFToSM");
        if(norm_type == "SM_rel_norm") s.ReplaceAll("EFT","EFToSMnorm");
        c->Print(s,"png");
        c->Print(s.ReplaceAll("png","pdf"),"pdf");
    }

    //// Print maxYvals dictionary:
    //for (auto it = maxYvals_dict.begin(); it != maxYvals_dict.end(); it++ ) {
    //    std::cout << it->first  // string (key)
    //              << ": "
    //              << it->second   // string's value
    //              << std::endl ;
    //}


    delete sm_pt;
    delete wc_pt;
    for (auto f: files) {
        f->Close();
    }
}


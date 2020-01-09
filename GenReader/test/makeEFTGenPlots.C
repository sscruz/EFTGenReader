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

//void makeEFTGenPlots(TString output_fname,std::vector<TString> input_fnames) {
void makeEFTGenPlots(std::vector<TString> input_fnames, TString wc_string) {
    std::vector<TFile*> files;

    TH1::SetDefaultSumw2();

    //const Int_t kCLRS = 9;
    const Int_t kCLRS = 5;
    Int_t palette[kCLRS];
    //palette[0] = kBlack;
    //palette[1] = kBlue;
    //palette[2] = kRed;
    //palette[3] = kGreen;
    //palette[4] = kMagenta;
    //palette[5] = kCyan;
    //palette[6] = kCyan+3;
    //palette[7] = kGreen+3;
    //palette[8] = kRed+3;
    //palette[9] = kMagenta+3;

    string colorScheme = "standard";
    if (colorScheme == "standard") {
        palette[0] = kBlack;
        palette[1] = kCyan;
        palette[2] = kBlue;
        palette[3] = kGreen;
        palette[4] = kRed;
    }

    // TLegend parameters
    double left,right,top,bottom,scale_factor,minimum;
    left         = 0.81;
    right        = 0.98;
    top          = 0.9;
    scale_factor = 0.05;
    minimum      = 0.1;

    gStyle->SetPalette(kCLRS,palette);
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
    std::map<string,float> maxYvals_dict;
    for (auto f1: files1) {
        TDirectory* td1 = f1->GetDirectory("EFTGenReader");
        TKey* key1;
        TIter next(td1->GetListOfKeys());
        while ((key1 = (TKey*)next())) {
            TString s1 = key1->GetName();
            TH1EFT* h_eventsumEFT1 = (TH1EFT*)td1->Get("h_eventsumEFT");
            double eventsum_SM1 = h_eventsumEFT1->GetBinFit(1).evalPoint(sm_pt);
            TH1EFT* h1 = (TH1EFT*)td1->Get(key1->GetName());
            if (s1.Index("EFT") != -1) {
                for (Int_t bin_idx = 0; bin_idx <= h1->GetNbinsX()+1; bin_idx++) {
                    double wcfit_bin_val = h1->GetBinFit(bin_idx).evalPoint(wc_pt);
                    h1->SetBinContent(bin_idx,wcfit_bin_val);
                }
            }
            //h1->Scale(1.0/eventsum_SM1);
            if (s1 != "h_SMwgt_norm") {
                h1->Scale(1.0/eventsum_SM1);
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

    for (auto f: files) {
        f->Print();

        TString fname = f->GetName();
        cout << "name of the file !!! " << fname << std::endl;
        TString sub_str;
        Ssiz_t idx = fname.First('/');
        //cout << "idx !!! " << idx << std::endl;
        int idx_begin = 0;
        int idx_end = 0;
        if (idx != TString::kNPOS) {
            sub_str = fname(idx+1,fname.Length());
        }
        TString marker = "output_";
        //TString marker = "0p_";
        idx_begin = sub_str.Index(marker)+marker.Length();
        //idx_end = sub_str.Index("fromMAOD.root");
        idx_end = sub_str.Index(".root");
        sub_str = sub_str(idx_begin,idx_end-idx_begin);
        //cout << " this is the sub str now !!!!!!!!!! " << sub_str << std::endl;

        TDirectory* td = f->GetDirectory("EFTGenReader");

        TH1EFT* h_eventsumEFT = (TH1EFT*)td->Get("h_eventsumEFT");
        double eventsum_SM = h_eventsumEFT->GetBinFit(1).evalPoint(sm_pt);

        TIter next(td->GetListOfKeys());
        TKey* key; TCanvas* c; TLegend* leg;
        while ((key = (TKey*)next())) {
            TString cmp = "summaryTree";
            TString s = key->GetName();
            if (s == cmp) continue;

            /*
            // FOR NOT PLOTTING EFT PLOTS
            if (s.Index("EFT") != -1) {
                continue;
            }
            */

            //TH1D* h = (TH1D*)td->Get(key->GetName());
            TH1EFT* h = (TH1EFT*)td->Get(key->GetName());
            h->SetMarkerStyle(kFullCircle);
            h->SetMarkerSize(0.25);
            h->SetOption("E");

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

            bool is_new = !c;
            if (is_new) {
                c = new TCanvas(canv_str,key->GetName(),1280,720);

                canvs.push_back(c);
                
                bottom = std::max(top - scale_factor*(files.size()+1),minimum);
                leg = new TLegend(left,top,right,bottom);
                legs.push_back(leg);
            } else {
                leg = legs.at(c_idx);
            }

            //c->cd();
            // Log log scale for wgt hist
            TPad* canv_pad = (TPad*)c->cd();
            if (s == "h_SMwgt_norm") {
                canv_pad->SetLogy(1);
                canv_pad->SetLogx(1);
            }

            if (s.Index("EFT") != -1) {
                //std::cout << h->GetName() << std::endl;
                for (Int_t bin_idx = 0; bin_idx <= h->GetNbinsX()+1; bin_idx++) {
                    double wcfit_bin_val = h->GetBinFit(bin_idx).evalPoint(wc_pt);
                    double wcfit_bin_err = h->GetBinFit(bin_idx).evalPointError(wc_pt);
                    h->SetBinContent(bin_idx,wcfit_bin_val);
                    h->SetBinError(bin_idx,wcfit_bin_err);
                }
            }
            // SM norm
            //h->Scale(1.0/eventsum_SM);
            //h->GetYaxis()->SetRangeUser(0.0,1.2*maxYvals_dict[string(s)]);
            if (s != "h_SMwgt_norm") {
                h->Scale(1.0/eventsum_SM);
                h->GetYaxis()->SetRangeUser(0.0,1.2*maxYvals_dict[string(s)]);
            } else {
                std::cout << "\nNot normalizing this histogram: " << s << "\n" << std::endl;
                h->GetYaxis()->SetRangeUser(0.1,2*maxYvals_dict[string(s)]);
            }

            /*
            // Unit norm
            Int_t nbins = h->GetNbinsX();
            Double_t intg = h->Integral(0,nbins+1);
            if (intg > 1.0) {
               h->Scale(1./intg);
            }
            */

            if (is_new) {
                h->Draw("E PLC PMC");
            } else {
                h->Draw("E SAME PLC PMC");
            }
            leg->AddEntry(h,sub_str,"l");
            c->Update();

        }

        filecounter++;
    }

    for (auto c: canvs) {

        c->cd();
        int idx = findCanvasIndex(c->GetTitle(),canvs);

        legs.at(idx)->Draw();

        //TString s = (TString)c->GetTitle() + "_" + output_fname + ".png";
        // No need to uniquely name the images, since they will be placed in a dedicated directory
        TString s = (TString)c->GetTitle() + ".png";
        c->Print(s,"png");
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

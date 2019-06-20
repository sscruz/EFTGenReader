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
void makeEFTGenPlots(std::vector<TString> input_fnames) {
    std::vector<TFile*> files;

    const Int_t kCLRS = 6;
    Int_t palette[kCLRS];
    palette[0] = kBlack;
    palette[1] = kBlue;
    palette[2] = kRed;
    palette[3] = kGreen;
    palette[4] = kMagenta;
    palette[5] = kCyan;

    //const Int_t kCLRS = 5;
    //Int_t palette[kCLRS];
    //palette[kCLRS-1] = kBlack;
    //palette[kCLRS-2] = kBlue;
    //palette[kCLRS-3] = kRed;
    //palette[kCLRS-4] = kGreen+2;
    //palette[kCLRS-5] = kMagenta+2;

    // TLegend parameters
    double left,right,top,bottom,scale_factor,minimum;
    left         = 0.81;
    right        = 0.98;
    top          = 0.9;
    scale_factor = 0.05;
    minimum      = 0.1;

    gStyle->SetPalette(kCLRS,palette);
    gStyle->SetOptStat(0);
    gStyle->SetPadRightMargin(0.2);

    //std::cout << "Output File: " << output_fname << std::endl;
    for (auto s: input_fnames) {
        std::cout << "  Input File: " << s << std::endl;
        TFile* f = TFile::Open(s);
        if (f) {
            f->Print();
            files.push_back(f);
        }
    }

    if (!files.size()) {
        return;
    }

    std::cout << "" << std::endl;

    std::vector<TCanvas*> canvs;
    std::vector<TLegend*> legs;
    int filecounter = 0;
    for (auto f: files) {
        f->Print();

        TString fname = f->GetName();

        Ssiz_t idx = fname.First('_');
        TString sub_str = fname(0,idx);

        TDirectory* td = f->GetDirectory("EFTGenReader");
        TIter next(td->GetListOfKeys());
        TKey* key; TCanvas* c; TLegend* leg;
        while ((key = (TKey*)next())) {
            TString cmp = "summaryTree";
            TString s = key->GetName();
            if (s == cmp) continue;

            TH1D* h = (TH1D*)td->Get(key->GetName());
            h->SetMarkerStyle(kFullCircle);

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

            c->cd();

            if (is_new) {
                h->Draw("PLC PMC");
            } else {
                h->Draw("SAME PLC PMC");
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

    for (auto f: files) {
        f->Close();
    }
}
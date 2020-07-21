#ifndef TH2EFT_H_
#define TH2EFT_H_

#include "TH2D.h"
#include <vector>
#include "TClass.h"
#include "WCFit.h"
#include "WCPoint.h"

class TH2EFT : public TH2D
{
    public:
    
        // ROOT needs these:
        TH2EFT();
        ~TH2EFT();
        
        // usual constructor:
        TH2EFT(const char *name, const char *title, Int_t nbinsx, Double_t xlow, Double_t xup, Int_t nbinsy, Double_t ylow, Double_t yup);
        
        std::vector<WCFit> hist_fits;
        std::vector<WCFit> hist_fitsx;
        std::vector<WCFit> hist_fitsy;
        //TODO(maybe?): Add over/underflow bin fit functions and update Fill to use them accordingly
        WCFit overflow_fit;
        WCFit underflow_fit;

        using TH2D::Fill;           // Bring the TH2D Fill fcts into scope
        using TH2D::GetBinContent;  // Bring the TH2D GetBinContent fcts into scope
        using TH2D::Scale;          // Bring the TH2D Scale fcts into scope (likely not needed)

        Int_t Fill(Double_t x, Double_t y, Double_t w, WCFit fit);
        WCFit GetBinFit(Int_t bin);
        WCFit GetSumFit(Int_t axis);
        Double_t GetBinContent(Int_t binx, Int_t biny, WCPoint wc_pt);
        //TH2EFT* Scale(WCPoint wc_pt);
        void Scale(WCPoint wc_pt);
        void ScaleFits(double amt, Int_t axis);
        void DumpFits();
        
        void SetBins(Int_t nx, Double_t xmin, Double_t xmax, Int_t ny, Double_t ymin, Double_t ymax);  // overriding virtual function from TH2
        Bool_t Add(const TH2 *h1, Double_t c1=1); // overriding virtual function from TH2
        Long64_t Merge(TCollection* list);

        ClassDef(TH2EFT,1); // ROOT needs this here
        //TODO(maybe?): Add member function to return specifically fit coeffs (rather then entire WCFit object)
};

// ROOT needs this here:
ClassImp(TH2EFT);
TH2EFT::TH2EFT() {}
TH2EFT::~TH2EFT() {}


TH2EFT::TH2EFT(const char *name, const char *title, Int_t nbinsx, Double_t xlow, Double_t xup, Int_t nbinsy, Double_t ylow, Double_t yup)
 : TH2D (name, title, nbinsx, xlow, xup, nbinsy, ylow, yup) 
{
    // Create/Initialize a fit function for each bin in the histogram
    WCFit new_fit;
    for (Int_t i = 0; i < this->FindFixBin(xup,yup); i++) {
        this->hist_fits.push_back(new_fit);
    }
    for (Int_t i = 0; i < nbinsx; i++) {
        this->hist_fitsx.push_back(new_fit);
    }
    for (Int_t i = 0; i < nbinsy; i++) {
        this->hist_fitsy.push_back(new_fit);
    }
}
void TH2EFT::SetBins(Int_t nx, Double_t xmin, Double_t xmax, Int_t ny, Double_t ymin, Double_t ymax)
{
    // Use this function with care! Non-over/underflow bins are simply
    // erased and hist_fits re-sized with empty fits.
    
    hist_fits.clear();
    hist_fitsx.clear();
    hist_fitsy.clear();
    WCFit new_fit;
    for (Int_t i = 0; i < nx; i++) {
        this->hist_fitsx.push_back(new_fit);
    }
    for (Int_t i = 0; i < ny; i++) {
        this->hist_fitsy.push_back(new_fit);
    }
    for (Int_t i = 0; i < nx+ny; i++) {
        this->hist_fits.push_back(new_fit);
    }
    
    TH2::SetBins(nx, xmin, xmax, ny, ymin, ymax);
}
Bool_t TH2EFT::Add(const TH2 *h1, Double_t c1)
{
    // check whether the object pointed to inherits from (or is a) TH2EFT:
    if (h1->IsA()->InheritsFrom(TH2EFT::Class())) {
        if ((this->hist_fitsx.size() == ((TH2EFT*)h1)->hist_fitsx.size()) && (this->hist_fitsy.size() == ((TH2EFT*)h1)->hist_fitsy.size())) {
            for (unsigned int i = 0; i < this->hist_fitsx.size(); i++) {
                // assumes this hist and the one whose fits we're adding have the same bins!
                this->hist_fitsx[i].addFit( ((TH2EFT*)h1)->hist_fitsx[i] );
            }
            for (unsigned int i = 0; i < this->hist_fitsy.size(); i++) {
                // assumes this hist and the one whose fits we're adding have the same bins!
                this->hist_fitsy[i].addFit( ((TH2EFT*)h1)->hist_fitsy[i] );
            }
        } else { 
            std::cout << "Attempt to add 2 TH2EFTs with different # of fits!" << std::endl;
            std::cout << this->hist_fitsx.size() << ", " << ((TH2EFT*)h1)->hist_fitsx.size() << std::endl;
            std::cout << this->hist_fitsy.size() << ", " << ((TH2EFT*)h1)->hist_fitsy.size() << std::endl;
        }
        this->overflow_fit.addFit( ((TH2EFT*)h1)->overflow_fit );
        this->underflow_fit.addFit( ((TH2EFT*)h1)->underflow_fit );
    }
    
    return TH2::Add(h1,c1); // I think this should work
}

// Custom merge function for using hadd
Long64_t TH2EFT::Merge(TCollection* list)
{
    TIter nexthist(list);
    TH2EFT *hist;
    while ((hist = (TH2EFT*)nexthist.Next())) {
        if (this->hist_fitsx.size() != hist->hist_fitsx.size()) {
            std::cout << "[WARNING] Skipping histogram with different # of fits" << std::endl;
            continue;
        }
        if (this->hist_fitsy.size() != hist->hist_fitsy.size()) {
            std::cout << "[WARNING] Skipping histogram with different # of fits" << std::endl;
            continue;
        }
        for (unsigned int i = 0; i < this->hist_fitsx.size(); i++) {
            this->hist_fitsx.at(i).addFit(hist->hist_fitsx.at(i));
        }
        for (unsigned int i = 0; i < this->hist_fitsy.size(); i++) {
            this->hist_fitsy.at(i).addFit(hist->hist_fitsy.at(i));
        }
        this->overflow_fit.addFit(hist->overflow_fit);
        this->underflow_fit.addFit(hist->underflow_fit);
    }

    return TH2::Merge(list);
}

Int_t TH2EFT::Fill(Double_t x, Double_t y, Double_t w, WCFit fit)
{
    Int_t bin_id = this->FindFixBin(x,y) - 1;
    Int_t nhists  = this->hist_fits.size();
    if (bin_id >= nhists) {
        // For now ignore events which enter overflow bin
        this->overflow_fit.addFit(fit);
        return Fill(x,y,w);
    } else if (bin_id < 0) {
        // For now ignore events which enter underflow bin
        this->underflow_fit.addFit(fit);
        return Fill(x,y,w);
    }
    this->hist_fits.at(bin_id).addFit(fit);
    return Fill(x,y,w); // the original TH2D member function
}

// Returns a fit function for a particular bin (no checks are made if the bin is an over/underflow bin)
WCFit TH2EFT::GetBinFit(Int_t bin)
{
    Int_t nhists = this->hist_fits.size();
    if (bin <= 0) {
        return this->underflow_fit;
    } else if (bin > nhists) {
        return this->overflow_fit;
    }
    else return this->hist_fits.at(bin - 1);
}

// Returns a WCFit whose structure constants are determined by summing structure constants from all bins
WCFit TH2EFT::GetSumFit(Int_t axis=0)
{
    WCFit summed_fit;
    if(axis == 1) {
        for (unsigned int i = 0; i < this->hist_fitsy.size(); i++) {
            summed_fit.addFit(this->hist_fitsy.at(i));
        }
        return summed_fit;
    }
    else {
        for (unsigned int i = 0; i < this->hist_fitsx.size(); i++) {
            summed_fit.addFit(this->hist_fitsx.at(i));
        }
        return summed_fit;
    }
}

// Returns a bin scaled by the the corresponding fit evaluated at a particular WC point
Double_t TH2EFT::GetBinContent(Int_t binx, Int_t biny, WCPoint wc_pt)
{
    Int_t bin = this->FindFixBin(binx-1,biny-1);
    if (this->GetBinFit(bin).getDim() <= 0) {
        // We don't have a fit for this bin, return regular bin contents
        return GetBinContent(binx, biny);
    }

    double scale_value = this->GetBinFit(bin).evalPoint(&wc_pt);
    Double_t num_events = GetBinContent(binx,biny);
    if (num_events == 0) {
        return 0.0;
    }

    return scale_value;
}
void TH2EFT::Scale(WCPoint wc_pt)
{
    // Warning: calling GetEntries after a call to this function will return a 
    // non-zero value, even if the histogram was never filled.
    
    for (Int_t i = 1; i <= this->GetNbinsX(); i++) {
        for (Int_t j = 1; j <= this->GetNbinsY(); j++) {
            Int_t bin = this->FindFixBin(i-1,j-1);
            Double_t new_content = this->GetBinContent(i,j,wc_pt);
            Double_t new_error = (GetBinFit(bin)).evalPointError(&wc_pt);
            this->SetBinContent(bin,new_content);
            this->SetBinError(bin,new_error);
        }
    }
    
}
// evalPointError disabled:
/*
void TH2EFT::Scale(WCPoint wc_pt)
{
    // Warning: calling GetEntries after a call to this function will return a 
    // non-zero value, even if the histogram was never filled.
    
    for (Int_t i = 1; i <= this->GetNbinsX(); i++) {
        Double_t old_content = this->GetBinContent(i);
        Double_t new_content = this->GetBinContent(i,wc_pt);
        Double_t old_error = this->GetBinError(i);
        this->SetBinContent(i,new_content);
        this->SetBinError(i,old_error*new_content/old_content);
    }
    
}
*/
// Uniformly scale all fits by amt
void TH2EFT::ScaleFits(double amt, Int_t axis=0)
{
    if(axis == 0){
        for (uint i = 0; i < this->hist_fitsx.size(); i++) {
            this->hist_fitsx.at(i).scale(amt);
        }
    }
    if(axis == 1){
        for (uint i = 0; i < this->hist_fitsy.size(); i++) {
            this->hist_fitsy.at(i).scale(amt);
        }
    }
}

// Display the fit parameters for all bins
void TH2EFT::DumpFits()
{
    for (uint i = 0; i < this->hist_fitsx.size(); i++) {
        this->hist_fitsx.at(i).dump();
    }
    for (uint i = 0; i < this->hist_fitsy.size(); i++) {
        this->hist_fitsy.at(i).dump();
    }
}
#endif

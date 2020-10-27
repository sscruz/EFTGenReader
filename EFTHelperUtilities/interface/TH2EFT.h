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
        TH2EFT(const char *name, const char *title, Int_t nbinsx, Double_t *x, Int_t nbinsy, Double_t* y);
        
        std::vector<WCFit> hist_fits;
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
        Double_t GetBinContent(Int_t bin, WCPoint wc_pt);
        TH2EFT* Rebin(Int_t nbinsx, Double_t *x, Int_t nbinsy, Double_t* y);
        //TH2EFT* Scale(WCPoint wc_pt);
        void Scale(WCPoint wc_pt);
        void ScaleFits(double amt, Int_t axis);
        void AddBinFit(Int_t bin, WCFit &fit);
        void DumpFits();
        void SetBinFit(Int_t bin, WCFit &fits);
        void Copy(TObject &obj) const;
        
        void SetBins(Int_t nx, Double_t xmin, Double_t xmax, Int_t ny, Double_t ymin, Double_t ymax);  // overriding virtual function from TH2
        Bool_t Add(const TH2 *h1, Double_t c1=1); // overriding virtual function from TH2
        Bool_t NormalizeTo(const TH2D *h1, Double_t c1=0.5);
        Long64_t Merge(TCollection* list);

        ClassDef(TH2EFT,1); // ROOT needs this here
        //TODO(maybe?): Add member function to return specifically fit coeffs (rather then entire WCFit object)
};

// ROOT needs this here:
//ClassImp(TH2EFT);

#endif

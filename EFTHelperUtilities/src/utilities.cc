//#include "../interface/utilities.h"
#include "EFTGenReader/EFTHelperUtilities/interface/utilities.h"

void binLogX(TH1D* hist) {
    TAxis* axis = hist->GetXaxis();
    int bins = axis->GetNbins();

    Double_t fro = axis->GetXmin();
    Double_t to  = axis->GetXmax();
    Double_t width = (to - fro) / bins;
    Axis_t *new_bins = new Axis_t[bins+1];

    for (int i = 0; i <= bins; i++) {
        new_bins[i] = TMath::Power(10,fro + i*width);
    }
    axis->Set(bins,new_bins);
    delete[] new_bins;
}

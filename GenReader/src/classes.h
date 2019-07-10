#include "EFTGenReader/GenReader/interface/WCPoint.h"
#include "EFTGenReader/GenReader/interface/WCFit.h"
#include "EFTGenReader/GenReader/interface/TH1EFT.h"

#ifndef EFTGENREADER_DUMMY_OBJS
#define EFTGENREADER_DUMMY_OBJS

namespace {
    struct WC_Objects {
        WCPoint dummy_wcpoint;
        WCFit dummy_wcfit;
        TH1EFT dummy_th1eft;
    };
}

#endif
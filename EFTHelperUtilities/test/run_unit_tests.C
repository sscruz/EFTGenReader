#include "EFTGenReader/EFTHelperUtilities/interface/TH1EFT.h"
#include "EFTGenReader/EFTHelperUtilities/interface/WCPoint.h"
#include "EFTGenReader/EFTHelperUtilities/interface/WCFit.h"

// run with: root -l -b -q run_unit_tests.C

bool test_wcfit() {
    std::string chk_str;

    bool unit_chk;
    int all_chks,units;

    // The structure constants
    double s00 = 1.0;
    double s10 = 1.5;
    double s11 = 1.25;

    std::vector<WCPoint> pts;
    std::vector<double> vals {-1.0,1.25,0.5,2.5,4};
    for (auto x: vals) {
        double y = s00*1.0 + s10*x + s11*x*x;
        pts.push_back(WCPoint(TString::Format("wgt_ctG_%.2f",x).Data(),y));
    }

    double chk_x = 1.5;
    double chk_y = s00*1.0 + s10*chk_x + s11*chk_x*chk_x;
    WCPoint chk_pt(TString::Format("wgt_ctG_%.2f",chk_x).Data(),0.0);

    std::cout << "Running unit tests for WCFit class" << std::endl;
    all_chks = 0;
    units = 0;

    WCFit fit_base = WCFit(pts,"base");
    unit_chk = (fit_base.evalPoint(&chk_pt) == chk_y);
    all_chks += unit_chk;
    units += 1;

    chk_str = unit_chk ? "Passed" : "Failed";
    std::cout << "--- UNIT 1 ---" << std::endl;
    std::cout << "chk_x    : " << chk_x << std::endl;
    std::cout << "chk_y    : " << chk_y << std::endl;
    std::cout << "evalPoint: " << fit_base.evalPoint(&chk_pt) << std::endl;
    std::cout << "test: " << chk_str << std::endl;
    std::cout << "--------------\n" <<std::endl;

    WCFit fit_new = WCFit();
    fit_new.setTag("new");
    fit_new.addFit(fit_base);
    unit_chk = (fit_base.evalPoint(&chk_pt) == chk_y);
    all_chks += unit_chk;
    units += 1;

    chk_str = unit_chk ? "Passed" : "Failed";
    std::cout << "--- UNIT 2 ---" << std::endl;
    std::cout << "chk_x    : " << chk_x << std::endl;
    std::cout << "chk_y    : " << chk_y << std::endl;
    std::cout << "evalPoint: " << fit_new.evalPoint(&chk_pt) << std::endl;
    std::cout << "test: " << chk_str << std::endl;
    std::cout << "--------------\n" <<std::endl;

    fit_new.addFit(fit_base);
    unit_chk = (fit_new.evalPoint(&chk_pt) == 2*chk_y);
    all_chks += unit_chk;
    units += 1;

    chk_str = unit_chk ? "Passed" : "Failed";
    std::cout << "--- UNIT 3 ---" << std::endl;
    std::cout << "chk_x    : " << chk_x << std::endl;
    std::cout << "chk_y    : " << 2*chk_y << std::endl;
    std::cout << "evalPoint: " << fit_new.evalPoint(&chk_pt) << std::endl;
    std::cout << "test: " << chk_str << std::endl;
    std::cout << "--------------\n" <<std::endl;

    unit_chk = (fit_base.evalPoint(&chk_pt) == chk_y);
    all_chks += unit_chk;
    units += 1;

    chk_str = unit_chk ? "Passed" : "Failed";
    std::cout << "--- UNIT 4 ---" << std::endl;
    std::cout << "chk_x    : " << chk_x << std::endl;
    std::cout << "chk_y    : " << chk_y << std::endl;
    std::cout << "evalPoint: " << fit_base.evalPoint(&chk_pt) << std::endl;
    std::cout << "test: " << chk_str << std::endl;
    std::cout << "--------------\n" <<std::endl;

    std::cout << "Passed Checks: " << TString::Format("%d/%d",all_chks,units) << std::endl;
    return (all_chks == units);
}


bool test_th1eft() {
    std::string chk_str;
    bool unit_chk;
    int all_chks,units;
    double result,expected,diff,tolerance;

    // The structure constants
    double s00 = 1.0;
    double s10 = 1.5;
    double s11 = 1.25;

    // A dummy WC name to use
    std::string wc_name = "ctG";

    std::vector<WCPoint> pts;
    std::vector<double> vals {-1.0,1.25,0.5,2.5,4};
    for (auto x: vals) {
        double y = s00*1.0 + s10*x + s11*x*x;
        pts.push_back(WCPoint(TString::Format("wgt_%s_%.2f",wc_name.c_str(),x).Data(),y));
    }

    WCFit fit_1(pts,"f1");
    WCFit fit_2 = WCFit();
    fit_2.setTag("f2");

    fit_2.addFit(fit_1);
    fit_2.addFit(fit_1);

    double chk_x = 1.5;
    double chk_y = s00*1.0 + s10*chk_x + s11*chk_x*chk_x;
    WCPoint chk_pt(TString::Format("wgt_%s_%.2f",wc_name.c_str(),chk_x).Data(),0.0);

    std::cout << "Running unit tests for TH1EFT class" << std::endl;
    all_chks = 0;
    units = 0;

    TH1EFT* h_base = new TH1EFT("h_base","h_base",1,0,1);

    h_base->Fill(0.5,1.0,fit_1);

    expected = 1.0;
    result = h_base->GetBinContent(1);

    unit_chk = (result == expected);
    all_chks += unit_chk;
    units += 1;

    chk_str = unit_chk ? "Passed" : "Failed";
    std::cout << "--- UNIT 1 ---" << std::endl;
    std::cout << "expected     : " << expected << std::endl;
    std::cout << "GetBinContent: " << result << std::endl;
    std::cout << "test: " << chk_str << std::endl;
    std::cout << "--------------\n" <<std::endl;

    //////////////////////////////////////////////////////

    expected = chk_y;
    result = h_base->GetBinContent(1,chk_pt);

    unit_chk = (result == expected);
    all_chks += unit_chk;
    units += 1;

    chk_str = unit_chk ? "Passed" : "Failed";
    std::cout << "--- UNIT 2 ---" << std::endl;
    std::cout << "expected     : " << expected << std::endl;
    std::cout << "GetBinContent: " << result << std::endl;
    std::cout << "test: " << chk_str << std::endl;
    std::cout << "--------------\n" <<std::endl;

    //////////////////////////////////////////////////////

    h_base->Scale(chk_pt);

    expected = fit_1.evalPoint(&chk_pt);
    result = h_base->GetBinContent(1);

    unit_chk = (result == expected);
    all_chks += unit_chk;
    units += 1;

    chk_str = unit_chk ? "Passed" : "Failed";
    std::cout << "--- UNIT 3 ---" << std::endl;
    std::cout << "chk_x        : " << chk_pt.getStrength(wc_name) << std::endl;
    std::cout << "expected     : " << expected << std::endl;
    std::cout << "GetBinContent: " << result << std::endl;
    std::cout << "test: " << chk_str << std::endl;
    std::cout << "--------------\n" <<std::endl;

    //////////////////////////////////////////////////////

    chk_x = 0.75;
    chk_y = s00*1.0 + s10*chk_x + s11*chk_x*chk_x;
    chk_pt.setStrength(wc_name,chk_x);
    h_base->Scale(chk_pt);

    expected = fit_1.evalPoint(&chk_pt);
    result = h_base->GetBinContent(1);

    unit_chk = (result == expected);
    all_chks += unit_chk;
    units += 1;

    chk_str = unit_chk ? "Passed" : "Failed";
    std::cout << "--- UNIT 4 ---" << std::endl;
    std::cout << "chk_x        : " << chk_pt.getStrength(wc_name) << std::endl;
    std::cout << "expected     : " << expected << std::endl;
    std::cout << "GetBinContent: " << result << std::endl;
    std::cout << "test: " << chk_str << std::endl;
    std::cout << "--------------\n" <<std::endl;

    //////////////////////////////////////////////////////

    h_base->Fill(0.5,1.0,fit_2);
    h_base->Scale(chk_pt);

    // First make sure the original WCFits weren't messed with
    expected = chk_y + 2*chk_y;
    result = fit_1.evalPoint(&chk_pt) + fit_2.evalPoint(&chk_pt);
    diff = fabs(expected - result);
    tolerance = 1e-10;

    unit_chk = (diff < tolerance);
    all_chks += unit_chk;
    units += 1;

    chk_str = unit_chk ? "Passed" : "Failed";
    std::cout << "--- UNIT 5 ---" << std::endl;
    std::cout << "chk_x        : " << chk_pt.getStrength(wc_name) << std::endl;
    std::cout << "expected     : " << expected << std::endl;
    std::cout << "fit_1 + fit_2: " << result << std::endl;
    std::cout << "difference   : " << diff << std::endl;
    std::cout << "tolerance    : " << tolerance << std::endl;
    std::cout << "test: " << chk_str << std::endl;
    std::cout << "--------------\n" <<std::endl;

    // Now check that the TH1EFT actually worked
    expected = fit_1.evalPoint(&chk_pt) + fit_2.evalPoint(&chk_pt);
    result = h_base->GetBinContent(1);

    unit_chk = (result == expected);
    all_chks += unit_chk;
    units += 1;

    chk_str = unit_chk ? "Passed" : "Failed";
    std::cout << "--- UNIT 6 ---" << std::endl;
    std::cout << "chk_x        : " << chk_pt.getStrength(wc_name) << std::endl;
    std::cout << "expected     : " << expected << std::endl;
    std::cout << "GetBinContent: " << result << std::endl;
    std::cout << "test: " << chk_str << std::endl;
    std::cout << "--------------\n" <<std::endl;

    //////////////////////////////////////////////////////

    TH1EFT* h_new = (TH1EFT*)h_base->Clone("h_new");

    chk_x = 0.975;
    chk_y = s00*1.0 + s10*chk_x + s11*chk_x*chk_x;
    chk_pt.setStrength(wc_name,chk_x);
    
    h_new->Scale(chk_pt);

    // First check that h_new has the right value
    expected = fit_1.evalPoint(&chk_pt) + fit_2.evalPoint(&chk_pt);
    result = h_new->GetBinContent(1);

    unit_chk = (result == expected);
    all_chks += unit_chk;
    units += 1;

    chk_str = unit_chk ? "Passed" : "Failed";
    std::cout << "--- UNIT 7 ---" << std::endl;
    std::cout << "chk_x        : " << chk_pt.getStrength(wc_name) << std::endl;
    std::cout << "expected     : " << expected << std::endl;
    std::cout << "GetBinContent: " << result << std::endl;
    std::cout << "test: " << chk_str << std::endl;
    std::cout << "--------------\n" <<std::endl;
    
    chk_x = 0.75;    // Needs to be w/e chk_x was before UNIT 7
    chk_y = s00*1.0 + s10*chk_x + s11*chk_x*chk_x;
    chk_pt.setStrength(wc_name,chk_x);

    // Next check that the h_base was unaffected when we scaled h_new
    expected = fit_1.evalPoint(&chk_pt) + fit_2.evalPoint(&chk_pt);
    result = h_base->GetBinContent(1);

    unit_chk = (result == expected);
    all_chks += unit_chk;
    units += 1;

    chk_str = unit_chk ? "Passed" : "Failed";
    std::cout << "--- UNIT 8 ---" << std::endl;
    std::cout << "chk_x        : " << chk_pt.getStrength(wc_name) << std::endl;
    std::cout << "expected     : " << expected << std::endl;
    std::cout << "GetBinContent: " << result << std::endl;
    std::cout << "test: " << chk_str << std::endl;
    std::cout << "--------------\n" <<std::endl;

    //////////////////////////////////////////////////////

    std::cout << "Passed Checks: " << TString::Format("%d/%d",all_chks,units) << std::endl;
    return (all_chks == units);
}

void run_unit_tests() {
    bool all_chks = true;

    all_chks = test_wcfit() && all_chks;
    std::cout << std::endl;

    all_chks = test_th1eft() && all_chks;
    std::cout << std::endl;

    if (all_chks) {
        std::cout << "All unit tests completed successfully!" << std::endl;
    } else {
        std::cout << "Some unit tests failed!" << std::endl;
    }

    return;
}

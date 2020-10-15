#include "EFTGenReader/EFTHelperUtilities/interface/TH1EFT.h"
#include "EFTGenReader/EFTHelperUtilities/interface/WCPoint.h"
#include "EFTGenReader/EFTHelperUtilities/interface/WCFit.h"

// run with: root -l -b -q run_unit_tests.C

double fval(std::vector<double> xvals, std::vector<double> svals) {
    // Ordering convention for the structure constants:
    // Dim=0 (0,0)
    // Dim=1 (0,0) (1,0) (1,1)
    // Dim=2 (0,0) (1,0) (1,1) (2,0) (2,1) (2,2)
    double y = 0.0;
    int idx = 0;
    for (int i=0; i < xvals.size(); i++) {
        for (int j=0; j <= i ; j++) {
            double c1 = xvals.at(i);
            double c2 = xvals.at(j);
            double s  = svals.at(idx);
            y += s*c1*c2;
            //std::cout << TString::Format("(%d,%d) ",i,j);
            idx++;
        }
    }
    //std::cout << std::endl;
    return y;
}

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

bool test_stats() {
    std::string chk_str;
    bool unit_chk;
    int all_chks,units;
    double result,expected,diff,tolerance;

    // Basically the SM 'strength'
    double x0 = 1.0;

    // Dummy WC names to use (needs to match dimension of pt
    std::vector<std::string> wc_names {"sm","ctG","ctZ"};

    // The structure constants, need to match dimension of pt
    std::vector<double> svals {
        1.15, // (00)
        1.35,1.25, // (10) (11)
        0.25,0.75,1.00, // (20) (21) (22)
    };
    // Make sure there are enough pts to fully determine the fit!
    std::vector<std::vector<double> > pts {
        {x0,-1.00, 0.00},
        {x0,-0.50, 0.25},
        {x0, 0.00, 0.35},
        {x0, 0.25, 0.05},
        {x0, 0.50,-0.05},
        {x0, 0.75, 0.25},
        {x0, 1.00,-0.35},
    };

    std::vector<WCPoint> wc_pts;
    for (auto pt: pts) {
        double y = fval(pt,svals);
        TString str = "wgt";
        for (int i=1; i < pt.size(); i++) { // NOTE: pt better not be size 0!!
            std::string wc_str = wc_names.at(i);
            str += TString::Format("_%s_%.2f",wc_str.c_str(),pt.at(i));
        }
        wc_pts.push_back(WCPoint(str.Data(),y));
    }

    WCFit fit_1(wc_pts,"f1");
    WCFit fit_2 = WCFit();
    fit_2.setTag("f2");

    int nevents = 5000;
    for (int i=0; i < nevents; i++) {
        fit_2.addFit(fit_1);
    }

    //////////////////////////////////////////////////////

    std::cout << "Running unit tests for stats unc." << std::endl;
    all_chks = 0;
    units = 0;

    // Needs to be the same size as wc_names
    std::vector<double> chk_x {x0,1.2,0.4};
    double chk_y = 0.0;
    double chk_e = 0.0;
    for (int i=0; i < nevents; i++) {
        double v = fval(chk_x,svals);
        chk_y += v;
        chk_e += v*v;
    }   
    chk_e = sqrt(chk_e);

    TString chk_wcstr("wgt");
    int sidx = 0;
    for (int i=0; i < wc_names.size(); i++) {
        if (i) { // Need to skip first entry since that's the SM 'strength'
            chk_wcstr += TString::Format("_%s_%.2f",wc_names.at(i).c_str(),chk_x.at(i));
        }
        for (int j=0; j <= i; j++) {
            double v = svals.at(sidx);
            std::cout << TString::Format("s%d%d: %.2f",i,j,v) << std::endl;
            sidx++;
        }
    }
    std::cout << std::endl;
    WCPoint chk_pt(chk_wcstr.Data(),0.0);

    //////////////////////////////////////////////////////

    // Basic check for proper adding of quadratic structure constants
    // Note: We expect the diff to grow with increased number of events due to the numeric precison
    expected = chk_y;
    result = fit_2.evalPoint(&chk_pt);
    diff = abs(expected - result);
    tolerance = 1e-4;

    unit_chk = (diff < tolerance); all_chks += unit_chk; units += 1;
    chk_str = unit_chk ? "Passed" : "Failed";
    std::cout << "--- UNIT 1 ---" << std::endl;
    std::cout << "evts     : " << nevents << std::endl;
    std::cout << "chk_wcstr: " << chk_wcstr << std::endl;
    std::cout << "expected : " << expected << std::endl;
    std::cout << "result   : " << result << std::endl;
    std::cout << "diff     : " << diff << std::endl;
    std::cout << "tolerance: " << tolerance << std::endl;
    std::cout << "test: " << chk_str << std::endl;
    std::cout << "--------------\n" <<std::endl;


    // Check the error calculation
    // Note: We expect the diff to grow with increased number of events due to the numeric precison
    expected = chk_e;
    result = fit_2.evalPointError(&chk_pt);
    diff = abs(expected - result);
    tolerance = 1e-05*sqrt(10*nevents);

    unit_chk = (diff < tolerance); all_chks += unit_chk; units += 1;
    chk_str = unit_chk ? "Passed" : "Failed";
    std::cout << "--- UNIT 2 ---" << std::endl;
    std::cout << "evts     : " << nevents << std::endl;
    std::cout << "chk_wcstr: " << chk_wcstr << std::endl;
    std::cout << "expected : " << expected << std::endl;
    std::cout << "result   : " << result << std::endl;
    std::cout << "diff     : " << diff << std::endl;
    std::cout << "tolerance: " << tolerance << std::endl;
    std::cout << "test: " << chk_str << std::endl;
    std::cout << "--------------\n" <<std::endl;

    // Now do the percent error
    // Note: The diff here also appears to grow apparently due to numeric precison, but much more slowly (it is still kind of concerning)
    expected = chk_e / chk_y;
    result = fit_2.evalPointError(&chk_pt) / fit_2.evalPoint(&chk_pt);
    diff = abs(expected - result);
    tolerance = 1e-04;

    unit_chk = (diff < tolerance); all_chks += unit_chk; units += 1;
    chk_str = unit_chk ? "Passed" : "Failed";
    std::cout << "--- UNIT 3 ---" << std::endl;
    std::cout << "evts     : " << nevents << std::endl;
    std::cout << "chk_wcstr: " << chk_wcstr << std::endl;
    std::cout << "expected : " << expected << std::endl;
    std::cout << "result   : " << result << std::endl;
    std::cout << "diff     : " << diff << std::endl;
    std::cout << "tolerance: " << tolerance << std::endl;
    std::cout << "test: " << chk_str << std::endl;
    std::cout << "--------------\n" <<std::endl;

    //////////////////////////////////////////////////////

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

    all_chks = test_stats() && all_chks;
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

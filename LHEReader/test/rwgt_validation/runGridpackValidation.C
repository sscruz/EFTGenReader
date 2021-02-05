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

#include "TRandom3.h"
#include "EFTGenReader/EFTHelperUtilities/interface/WCPoint.h"
#include "EFTGenReader/EFTHelperUtilities/interface/WCFit.h"
#include "EFTGenReader/EFTHelperUtilities/interface/TH1EFT.h"
#include "EFTGenReader/EFTHelperUtilities/interface/Stopwatch.h"
#include "EFTGenReader/EFTHelperUtilities/interface/split_string.h"
//#include "makeEFTPlots.h"
#include "make1DXsecPlots.h"

/*
    This script does the following:
        - Splits the plots by WC
        - 1D Inclusive Xsec WC scaling plots
        - Overlays multiple runs on the same plot
        - Overlays dedicated reference points on the corresponding plots
*/

//Ex: /hadoop/store/user/awightma/gridpack_scans/2018_05_06/scanpoints/ttH_2HeavyScan10kPilot_run0_scanpoints.txt

const std::string kMGStart   = "MGStart";   // The tag we use to designate MadGraph starting point in the scanpoints file
const std::string kOrig      = "original";
const std::string kOutputDir = "read_lhe_outputs";

/*
struct customPointInfo {
    TString proc_name;
    TString type;
    TString wc_name;
    double x;
    double y;
    double err;
};
*/

/*
// Print info about customPointsInfo object
void print_customPointInfo_object(std::vector<customPointInfo> points_info){
    for (int i=0; i<points_info.size(); i++){
        customPointInfo s = points_info.at(i);
        std::cout << "\nPrinting info about customPointInfo struct:" << std::endl;
        std::cout << "\tProc name: " << s.proc_name << ", Type: " << s.type << ", WC: " << s.wc_name << ", x val: " << s.x << ", y val: " << s.y << ", error: " << s.err << std::endl;
    }
}
*/

// Normalize a point (in the "customPointInfo" struct format) to a given value, propagate errors
std::vector<customPointInfo> normCustomPoint(std::vector<customPointInfo> points_info, double norm_pt_val, double norm_pt_err){
    std::cout << "norm pt: " << norm_pt_val << " norm err: " << norm_pt_err << std::endl;
    std::vector<customPointInfo> return_vector;
    for (int i=0; i<points_info.size(); i++){
        customPointInfo point_info = points_info.at(i);
        TString wc_name = point_info.wc_name;
        double y = point_info.y;
        double e = point_info.err;
        double y_norm = y/norm_pt_val;
        double e_norm = sqrt((e/y)*(e/y) + (norm_pt_err/norm_pt_val)*(norm_pt_err/norm_pt_val)) * (y/norm_pt_val);
        //std::cout << "\nBefore norm: " << wc_name << " y: " << y << " e: " << e << std::endl;
        //std::cout << "After norm: " << wc_name << " y: " << y_norm << " e: " << e_norm << "\n" << std::endl;
        point_info.y = y_norm;
        point_info.err = e_norm;
        return_vector.push_back(point_info);
    }
    return return_vector;
}

// Takes a vector of (x,y) points for a WC and returns a WCFit (where x = WC value, y = cross section at that value)
WCFit make_WCFit(std::string wc_name, std::string tag, std::vector<std::pair<double,double>> pts_pairs_vect, double norm_factor=1) {
    std::vector<WCPoint> pts_to_fit;
    for (int i; i<pts_pairs_vect.size(); i++){
        WCPoint wcpt;
        wcpt.setStrength(wc_name,pts_pairs_vect.at(i).first);
        wcpt.wgt = pts_pairs_vect.at(i).second / norm_factor;
        pts_to_fit.push_back(wcpt);
    }
    WCFit* wcfit = new WCFit(pts_to_fit,tag);
    return *wcfit;
}

// TestPrint
void print_xsec(WCFit wcfit){

    std::cout << "\nTag: " << wcfit.getTag() << std::endl;

    // Best fit val from AN
    std::map<std::string,float> bestfit_dict {
        {"ctW"   , -0.58},
        {"ctZ"   , -0.63},
        {"ctp"   , 25.50},
        {"cpQM"  , -1.07},
        {"ctG"   , -0.85},
        {"cbW"   , 3.17 },
        {"cpQ3"  , -1.81},
        {"cptb"  , 0.13 },
        {"cpt"   , -3.25},
        {"cQl3i" , -4.20},
        {"cQlMi" , 0.74 },
        {"cQei"  , -0.27},
        {"ctli"  , 0.33 },
        {"ctei"  , 0.33 },
        {"ctlSi" , -0.07},
        {"ctlTi" , -0.01}
    };

    WCPoint bestfit_wcpt;
    WCPoint sm_pt("smpt");
    WCPoint ctG2_wcpt;
    ctG2_wcpt.setStrength("ctG",2);

    ////////////////////
    /*
    // For Comp with 1601
    WCPoint comp_1601_pt;
    double ctW_val = -1.5;
    double ctZ_val = ctW_val*0.8768010;
    std::cout << "The ctZ and ctW values: " << ctW_val << " , " << ctZ_val << std::endl;
    comp_1601_pt.setStrength("ctW",ctW_val);
    std::cout << "At ctW_val " << ctW_val << " : " << wcfit.evalPoint(&comp_1601_pt) << " +- " << wcfit.evalPointError(&comp_1601_pt) << std::endl;
    comp_1601_pt.setStrength("ctZ",ctZ_val);
    std::cout << "At the 1601 comp point: " << wcfit.evalPoint(&comp_1601_pt) << " +- " << wcfit.evalPointError(&comp_1601_pt) << std::endl;
    std:: cout << " " << std::endl;

    WCPoint comp_1601_pt2;
    double ctW_val2 = -3;
    double ctZ_val2 = ctW_val2*0.8768010;
    std::cout << "The ctZ and ctW values: " << ctW_val2 << " , " << ctZ_val2 << std::endl;
    comp_1601_pt2.setStrength("ctW",ctW_val2);
    std::cout << "At ctW_val2 " << ctW_val2 << " : " << wcfit.evalPoint(&comp_1601_pt2) << " +- " << wcfit.evalPointError(&comp_1601_pt2) << std::endl;
    comp_1601_pt2.setStrength("ctZ",ctZ_val2);
    std::cout << "At the 1601 comp point: " << wcfit.evalPoint(&comp_1601_pt2) << " +- " << wcfit.evalPointError(&comp_1601_pt2) << std::endl;
    std:: cout << " " << std::endl;

    WCPoint comp_1601_pt3;
    double ctW_val3 = 3;
    double ctZ_val3 = ctW_val3*0.8768010;
    std::cout << "The ctZ and ctW values: " << ctW_val3 << " , " << ctZ_val3 << std::endl;
    comp_1601_pt3.setStrength("ctW",ctW_val3);
    std::cout << "At ctW_val3 " << ctW_val3 << " : " << wcfit.evalPoint(&comp_1601_pt3) << " +- " << wcfit.evalPointError(&comp_1601_pt3) << std::endl;
    comp_1601_pt3.setStrength("ctZ",ctZ_val3);
    std::cout << "At the 1601 comp point: " << wcfit.evalPoint(&comp_1601_pt3) << " +- " << wcfit.evalPointError(&comp_1601_pt3) << std::endl;
    std:: cout << " " << std::endl;
    */
    ////////////////////

    /*
    std::string wc;
    float wc_bestfit_val;
    float bestfit_val;
    float bestfit_val_err;
    for (auto it = bestfit_dict.begin(); it != bestfit_dict.end(); it++){
        wc = it->first;
        wc_bestfit_val = it->second; 
        std::cout << wc << " " << wc_bestfit_val << std::endl;
        bestfit_wcpt.setStrength(wc,wc_bestfit_val);
    }
    bestfit_val = wcfit.evalPoint(&bestfit_wcpt);
    bestfit_val_err = wcfit.evalPointError(&bestfit_wcpt);
    std::cout << "Best fit val is: " << bestfit_val << " +- " << bestfit_val_err << std::endl;
    std::cout << "SM val: " << wcfit.evalPoint(&sm_pt)<< " +- " << wcfit.evalPointError(&sm_pt)<<std::endl;
    std::cout << "ctG=2: " << wcfit.evalPoint(&ctG2_wcpt) << " +- " << wcfit.evalPointError(&ctG2_wcpt)<<std::endl;
    */

}

void runit(TString output_name,TString input_rundirs_spec,TString ref_rundirs_spec,TString grp_name="") {
    int run_idx = 0;

    ////////////////// Setting up the infor for the points from Reza //////////////////
    std::vector<customPointInfo> ttH_ctG_RezaPointsSMEFT;
    std::vector<customPointInfo> ttH_ctp_RezaPointsSMEFT;
    std::vector<customPointInfo> ttW_ctG_RezaPointsSMEFT;
    std::vector<customPointInfo> ttZ_ctG_RezaPointsSMEFT;
    std::vector<customPointInfo> ttZ_cpQM_RezaPointsSMEFT;
    std::vector<customPointInfo> ttZ_cpt_RezaPointsSMEFT;
    double ttH_NLO_SM = 0.4398;
    double ttH_NLO_SM_err = 0.001484;
    double ttW_NLO_SM = 0.5613;
    double ttW_NLO_SM_err = 0.0027;
    double ttZ_NLO_SM = 0.7333;
    double ttZ_NLO_SM_err = 0.003;
    ttH_ctG_RezaPointsSMEFT = {
        {"tth","nlo","ctG", -.3,0.3586,0.0013},
        {"tth","nlo","ctG", 0,ttH_NLO_SM,ttH_NLO_SM_err},
        {"tth","nlo","ctG", .3,0.6215,0.0023},
    };
    ttH_ctp_RezaPointsSMEFT = {
        {"tth","nlo","ctp", -50,7.935,0.024},
        {"tth","nlo","ctp", -30,3.835,0.01},
        {"tth","nlo","ctp", 0,ttH_NLO_SM,ttH_NLO_SM_err},
        {"tth","nlo","ctp", 5,0.2009,0.00068},
    };
    ttW_ctG_RezaPointsSMEFT = {
        {"ttW","nlo","ctG", -.3,0.5166,0.0022},
        {"ttW","nlo","ctG", 0,ttW_NLO_SM,ttW_NLO_SM_err},
        {"ttW","nlo","ctG", .3,0.621,0.0036},
    };
    ttZ_ctG_RezaPointsSMEFT = {
        {"ttZ","nlo","ctG", -.3,0.669,0.0027},
        {"ttZ","nlo","ctG", 0,ttZ_NLO_SM,ttZ_NLO_SM_err},
        {"ttZ","nlo","ctG", .3,0.8229,0.0033},
    };
    ttZ_cpQM_RezaPointsSMEFT = {
        {"ttZ","nlo","cpQM", -3,0.9815,0.0034},
        {"ttZ","nlo","cpQM", 0,ttZ_NLO_SM,ttZ_NLO_SM_err},
        {"ttZ","nlo","cpQM", 2,0.591,0.003523},
    };
    ttZ_cpt_RezaPointsSMEFT = {
        {"ttZ","nlo","cpt", -10,0.5501,0.0015},
        {"ttZ","nlo","cpt", 0,ttZ_NLO_SM,ttZ_NLO_SM_err},
        {"ttZ","nlo","cpt", 10,1.497,0.0045},
    };

    // Normalize the points to the SM
    ttH_ctG_RezaPointsSMEFT = normCustomPoint(ttH_ctG_RezaPointsSMEFT,ttH_NLO_SM,ttH_NLO_SM_err);
    ttH_ctp_RezaPointsSMEFT = normCustomPoint(ttH_ctp_RezaPointsSMEFT,ttH_NLO_SM,ttH_NLO_SM_err);
    ttW_ctG_RezaPointsSMEFT = normCustomPoint(ttW_ctG_RezaPointsSMEFT,ttW_NLO_SM,ttW_NLO_SM_err);
    ttZ_ctG_RezaPointsSMEFT = normCustomPoint(ttZ_ctG_RezaPointsSMEFT,ttZ_NLO_SM,ttZ_NLO_SM_err);
    ttZ_cpQM_RezaPointsSMEFT = normCustomPoint(ttZ_cpQM_RezaPointsSMEFT,ttZ_NLO_SM,ttZ_NLO_SM_err);
    ttZ_cpt_RezaPointsSMEFT = normCustomPoint(ttZ_cpt_RezaPointsSMEFT,ttZ_NLO_SM,ttZ_NLO_SM_err);

    std::vector<std::vector<customPointInfo>> reza_pts_info_vect = {ttH_ctG_RezaPointsSMEFT, ttH_ctp_RezaPointsSMEFT, ttW_ctG_RezaPointsSMEFT, ttZ_ctG_RezaPointsSMEFT, ttZ_cpQM_RezaPointsSMEFT, ttZ_cpt_RezaPointsSMEFT};
    for (auto points: reza_pts_info_vect){
        print_customPointInfo_object(points);
    };
    ////////////////////////////////////

    //return ;

    Stopwatch sw;

    std::vector<WCFit> target_fits;
    std::vector<WCFit> ref_fits;
    std::vector<WCPoint> ref_pts;

    // This will get filled with the params from arxiv 1607.05330 and 1601.08193 and passed to the plotting script
    // to be plotted for comparison with our 0j and 0+1j fits. The form of the vector is more or less:
    // [ ('LO name',[constant,lin,quad]) , ('NLO name',[constant,lin,quad]) ]
    std::vector<std::pair<std::string,std::vector<double>>> arxiv_fit_comps_vect;

    // Make custom WCFits to pass to the plot maker
    std::map<std::string,std::map<std::string,std::vector<WCFit>>> custom_fits_map; // e.g. {"ttH":{"ctG":[wcfit_lo,wcfit_nlo]}}
    std::map<std::string,std::map<std::string,std::vector<std::pair<double,double>>>> ptsMap; // e.g. {"ttH":{"ctG_lo":[(0,1),(1,2)]}}
    std::vector<std::string> custom_fits_wc_lst {"ctG","ctp","cpQM","cpt","ctZ","ctW"};

    double sm_ttH_nlo, sm_ttW_nlo, sm_ttZ_nlo, sm_ttZ_lo;
    /* // For PDF 2.3 (NOT what was used for the LO gridpacks):
    // ttH
    //double sm_ttH_nlo = 0.4521; // PDF 2.3
    //ptsMap["ttH"]["ctG_nlo"] = { {-1,0.5731} , {0,sm_ttH_nlo} , {1,1.482} };  // PDF 2.3
    //ptsMap["ttH"]["ctp_nlo"] = { {-5,0.7925} , {0,sm_ttH_nlo} , {5,0.2051} }; // PDF 2.3
    // ttW
    //double sm_ttW_nlo = 0.5591; // PDF 2.3
    //ptsMap["ttW"]["ctG_nlo"] = { {-1,0.4192} , {0,sm_ttW_nlo} , {5,0.7272} }; // PDF 2.3
    // ttZ
    //double sm_ttZ_nlo = 0.7441; // PDF 2.3
    //ptsMap["ttZ"]["ctG_nlo"]   = { {-1,0.6725} , {0,sm_ttZ_nlo} , {1,1.206} };  // PDF 2.3
    //ptsMap["ttZ"]["cpQM_nlo"]  = { {-5,1.208}  , {0,sm_ttZ_nlo} , {5,0.4392} }; // PDF 2.3
    //ptsMap["ttZ"]["cpt_nlo"]   = { {-5,0.5757} , {0,sm_ttZ_nlo} , {5,1.068} };  // PDF 2.3
    //ptsMap["ttZ"]["ctZ_nlo"]   = { {-5,2.261}  , {0,sm_ttZ_nlo} , {5,2.238} };  // PDF 2.3

    //double sm_ttZ_lo = 0.5841; //PDF 2.3
    //ptsMap["ttZ"]["cpt_lo"]   = { {-5,0.4395}  , {0,sm_ttZ_lo} , {5,0.8523} }; // PDF 2.3
    //ptsMap["ttZ"]["ctp_lo"]   = { {-5,0.5855}  , {0,sm_ttZ_lo} , {5,0.586} };  // PDF 2.3
    //ptsMap["ttZ"]["ctZ_lo"]   = { {-5,1.903}  , {0,sm_ttZ_lo} , {5,1.896} };   // PDF 2.3
    //ptsMap["ttZ"]["ctW_lo"]   = { {-5,0.5842}  , {0,sm_ttZ_lo} , {5,0.5846} }; // PDF 2.3
    */
    // NNPDF31_nnlo_hessian_pdfas PDF which has LHAPDF ID 306000
    //sm_ttZ_nlo = .7365;
    //ptsMap["ttZ"]["ctZ_nlo"]   = { {-2,.9671}  , {0,sm_ttZ_nlo} , {3,1.244} };
    //sm_ttZ_lo = .5219205;
    //ptsMap["ttZ"]["ctZ_lo"]   = { {-2,.7085581}  , {0,sm_ttZ_lo} , {3,.9369637} };

    //////////////////////////
    // Reza NLO SMEFT QED2, All previous had been QED1 (Note: QED=2 does not really seem to actually do anything, still just get QED=1 contributions)
    sm_ttZ_nlo = .732;
    ptsMap["ttZ"]["ctW_nlo"]   = { {-1.5,.7314}  , {0,sm_ttZ_nlo} , {.5,.7343} };
    sm_ttZ_lo = 6.471;
    ptsMap["ttZ"]["ctW_lo"]   = { {-1.5,6.482}  , {0,sm_ttZ_lo} , {.5,6.453} };

    sm_ttH_nlo = .4529;
    ptsMap["ttH"]["cpt_nlo"]   = { {-5,.4564}  , {0,sm_ttH_nlo} , {5,.449} };
    ptsMap["ttH"]["ctZ_nlo"]   = { {-5,.4526}  , {0,sm_ttH_nlo} , {5,.4536} };
    //////////////////////////

    // Make the WCFits
    for (auto wc : custom_fits_wc_lst){
        if ( ptsMap["ttH"].find(wc+"_nlo") != ptsMap["ttH"].end() ){
            custom_fits_map["ttH"][wc].push_back(make_WCFit(wc,"Reza smeft NLO",ptsMap["ttH"][wc+"_nlo"],sm_ttH_nlo));
        }
        if ( ptsMap["ttW"].find(wc+"_nlo") != ptsMap["ttW"].end() ){
            custom_fits_map["ttW"][wc].push_back(make_WCFit(wc,"Reza smeft NLO",ptsMap["ttW"][wc+"_nlo"],sm_ttW_nlo));
        }
        if ( ptsMap["ttZ"].find(wc+"_nlo") != ptsMap["ttZ"].end() ){
            custom_fits_map["ttZ"][wc].push_back(make_WCFit(wc,"Reza smeft NLO",ptsMap["ttZ"][wc+"_nlo"],sm_ttZ_nlo));
        }
        if (ptsMap["ttZ"].find(wc+"_lo") != ptsMap["ttZ"].end() ){
            custom_fits_map["ttZ"][wc].push_back(make_WCFit(wc,"Reza smeft LO",ptsMap["ttZ"][wc+"_lo"],sm_ttZ_lo));
        }
    }
    //ptsMap["ttH"]["_lo"]  = { {-10,} , {-5,} , {0,sm_ttH_lo} , {5,}, {10,} };
    //ptsMap["ttH"]["_nlo"] = { {-5,} , {0,sm_ttH_nlo} , {5,} };


    WCPoint sm_pt("smpt");

    std::string curr_process;

    PlotOptions xsec_plt_ops_1d;
    xsec_plt_ops_1d.output_dir = kOutputDir;
    xsec_plt_ops_1d.setXLimits(0.0,0.0);
    xsec_plt_ops_1d.setYLimits(0.0,1.3);

    std::vector<TString> all_dirs,tar_dirs,ref_dirs;
    TString fdir;
    std::ifstream input_filenames(input_rundirs_spec);
    while (input_filenames >> fdir) {
        tar_dirs.push_back(fdir);
        all_dirs.push_back(fdir);
    }
    input_filenames.close();
    std::ifstream ref_filenames(ref_rundirs_spec);
    while (ref_filenames >> fdir) {
        ref_dirs.push_back(fdir);
        if (!hasElement(all_dirs,fdir)) all_dirs.push_back(fdir);
    }
    ref_filenames.close();
    for (uint line_idx = 0; line_idx < all_dirs.size(); line_idx++) {
        fdir = all_dirs.at(line_idx);  // fdir will be a path to a hadoop directory with root files for a particluar run
        std::string run_dir = getRunDirectory(fdir.Data());
        std::cout << "[" << (line_idx+1) << "/" << all_dirs.size() << "] Full Path: " << fdir << std::endl;
        //std::cout << "\tRun Dir: " << run_dir << std::endl;

        bool is_tar = hasElement(tar_dirs,fdir);
        bool is_ref = hasElement(ref_dirs,fdir);

        std::vector<std::string> words;
        split_string(run_dir,words,"_");
        if (words.size() != 4) {
            std::cout << "[WARNING] Skipping invalid run directory!" << std::endl;
            continue;
        }

        if (run_idx == 0) {
            curr_process = words.at(1);
        }

        // Chain together all root files in the run directory
        TChain chain("EFTLHEReader/summaryTree");
        auto dir_files = getFiles(fdir);
        for (auto fn: dir_files) {
            //std::cout << "\tFiles: " << fn << std::endl;
            TString fname = fdir + "/" + fn;
            chain.Add(fname);
        }

        std::string process   = words.at(1);
        std::string grp_tag   = words.at(2);
        std::string run_label = words.at(3);

        /* This is all related to reading the scanpoints files, which is highly dependent on the naming scheme used for the samples

        std::string scanpoints_dir = getScanPointsDirectory(fdir.Data());
        std::string scanpoints_fpath = scanpoints_dir + process + "_" + grp_tag + "_" + run_label + "_scanpoints.txt";
        std::vector<WCPoint> scan_pts = parseScanPointsFile(scanpoints_fpath);

        WCPoint start_pt;
        bool found_start = false;
        for (auto& wc_pt: scan_pts) {
            if (wc_pt.tag == kMGStart) {
                start_pt = wc_pt;
                start_pt.wgt = 0.0;
                found_start = true;
                break;
            }
        }

        if (!found_start) {
            // Likely b/c the gridpack was produced before the scanpoint file was implemented, or it got deleted/moved
            std::cout << "[ERROR] Unable to find starting point in scanpoints file!" << std::endl;
            continue;
        }
        */

        int chain_entries = chain.GetEntries();
        int last_entry = chain_entries;
        if (chain_entries > 100000) { // 100k takes about 10min
            std::cout << "Chain_entries: " << chain_entries << std::endl;
            last_entry = 100000;
        }
        //last_entry = 100000; // For testing
        last_entry = 10; // For testing
        std::cout << "Last_entry: " << last_entry << std::endl;

        // Test LS stuff
        //int num_LS = 200;
        int num_LS = 10; // Testing
        std::set<int> unique_runs;

        int first_entry = 0;

        WCFit* wcFit_intree = 0;
        double originalXWGTUP_intree = -1.;
        int lumiBlock_intree;
        unordered_map<std::string,double>* eftwgts_intree = 0; // ResidCheck, note the "= 0" so that we do not get a seg fault

        chain.SetBranchAddress("originalXWGTUP",&originalXWGTUP_intree);
        chain.SetBranchAddress("wcFit",&wcFit_intree);
        chain.SetBranchAddress("lumiBlock",&lumiBlock_intree);

        chain.SetBranchAddress("eftwgts",&eftwgts_intree); // ResidCheck

        bool skip_progress = false;

        chain.GetEntry(first_entry);
        if (is_ref) {
            skip_progress = true;
        }

        WCFit inclusive_fit;
        //WCFit selection_fit;

        //int selection_events = 0;
        double genLep_pt1_intree;
        double genLep_pt2_intree;
        double genLep_pt3_intree;
        double genJet_pt4_intree;
        chain.SetBranchAddress("genLep_pt1",&genLep_pt1_intree);
        chain.SetBranchAddress("genLep_pt2",&genLep_pt2_intree);
        chain.SetBranchAddress("genLep_pt3",&genLep_pt3_intree);
        chain.SetBranchAddress("genJet_pt4",&genJet_pt4_intree);

        // Variables for 
        int lumiBlock_first;

        // Set up the wc point string (depends a lot on the naming scheme)
        std::map<string,string> ref_pts_dict;
        //std::string wcname = "ctG";
        //int range_max = 3; // for ctG
        std::string wcname = "cbW";
        float range_max = 4; // for cbW
        float npts = 5;
        float step = (range_max*2)/(npts-1);
        int run = 0;
        for (int wcval=-range_max; wcval<=range_max; wcval=wcval+step){
            //std::cout << "WCVAL : " << wcval << " LABEL: " << "run"+std::to_string(run) << " dict entry: " << "wcpt_"+wcname+"_"+std::to_string(wcval) << std::endl;
            ref_pts_dict["run"+std::to_string(run)] = "wcpt_"+wcname+"_"+std::to_string(wcval);
            run = run + 1;
        }
        std::cout << "Run label: " << run_label << " , Dictionary entry: " << ref_pts_dict[run_label] << std::endl;
        std::string pt_str = ref_pts_dict[run_label];
        WCPoint ref_fit_pt = WCPoint(pt_str);

        ofstream resids_file; // ResidTest
        resids_file.open("fit_resids_test_dir/test_resids_"+process+".txt"); // ResidTest

        sw.start("Full Loop");
        for (int i = first_entry; i < last_entry; i++) {
            sw.start("Event Loop");
            if (is_tar && !skip_progress) {
                printProgress(i - first_entry,last_entry - first_entry);
            }
            sw.start("Get Entry");
            chain.GetEntry(i);
            sw.lap("Get Entry");

            /*
            // Check for LS info for just running over one block:
            ////std::cout << lumiBlock_intree << std::endl;
            ////if (i == first_entry){
                ////lumiBlock_first = lumiBlock_intree;
            ////}
            ////if (lumiBlock_intree != lumiBlock_first){
                ////std::cout << "Breanking the loop!" << std::endl;
                ////break;
            ////}
            unique_runs.insert(lumiBlock_intree);
            if (unique_runs.size()+1 > num_LS) {
                break;
            }
            */

            //// ResidCheck ////
            //int counter = 0;
            for (auto it = eftwgts_intree->begin(); it != eftwgts_intree->end(); it++ ){
                std::string rwgt_str = it->first;
                WCPoint rwgt_wcpt = WCPoint(rwgt_str);
                double rwgt_wgt = it->second;
                double eval_wgt = wcFit_intree->evalPoint(&rwgt_wcpt);
                double resid = rwgt_wgt - eval_wgt;
                //std::cout << rwgt_str << "\n" << rwgt_wgt << " , " << eval_wgt << " , " << rwgt_wgt-eval_wgt << "\n" << std::endl;
                //std::cout << counter << " " << resid << std::endl;
                resids_file << resid << " " ;
                //counter = counter+1;
            }

            sw.start("Add Fit");
            inclusive_fit.addFit(*wcFit_intree);
            sw.lap("Add Fit");
            sw.lap("Event Loop");

            ref_fit_pt.wgt += originalXWGTUP_intree;

            //if (genLep_pt3_intree > 10 and genJet_pt4_intree > 100){
            //if (genLep_pt3_intree > 10 and genJet_pt4_intree > 30){
                //selection_events = selection_events + 1;
                //selection_fit.addFit(*wcFit_intree);
            //} else{
                //std::cout << "Skip; pt: " << genLep_pt3_intree << " , " << genJet_pt4_intree << std::endl;
            //}


        }
        sw.stop("Full Loop");
        resids_file.close(); // ResidsTest

        //std::cout << "\n Selected events over total: " << selection_events << "/" << last_entry << "->" << (float)selection_events/last_entry << "\n" << std::endl;

        // Normalize to SM
        //std::cout << "\nBefore all norm!!! incl SM xsec: " << inclusive_fit.evalPoint(&sm_pt) << " selection SM xsec: " << selection_fit.evalPoint(&sm_pt) << std::endl;

        double SM_xsec_incl = inclusive_fit.evalPoint(&sm_pt);
        //double SM_xsec_sel = selection_fit.evalPoint(&sm_pt);
        inclusive_fit.scale(1.0/SM_xsec_incl);
        /*
        // TMP LS stuff
        //inclusive_fit.scale(1.0/SM_xsec_incl); // TMP!!! do not norm to SM, LS stuff
        inclusive_fit.scale(1.0/unique_runs.size()); // TMP!! LS stuff
        std::cout << "\n\nnormalizing to LS!!! " << unique_runs.size() << "\n\n" << std::endl;
        selection_fit.scale(1.0/SM_xsec_sel);
        */

        //std::cout << "\nAfter all norm!!! incl SM xsec: " << inclusive_fit.evalPoint(&sm_pt) << " selection SM xsec: " << selection_fit.evalPoint(&sm_pt) << std::endl;

        ///* // Dump the fit functions
        //std::vector<std::string> list_of_WC = {"ctG","ctW"};
        std::vector<std::string> list_of_WC = {"ctp","cpQM","ctW","ctZ","ctG","cbW","cpQ3","cptb","cpt","cQl3i","cQlMi","cQei","ctli","ctei","ctlSi","ctlTi"};
        std::cout << " " << std:: endl;
        for (std::string WC : list_of_WC){ 
            inclusive_fit.dump(false,153,WC);
            //selection_fit.dump(false,153,WC);
            std::cout << " " << std:: endl;
        }
        //*/
        

        // Normalize ref pt and add to list
        std::cout << "Is ref? " << is_ref << std::endl;
        if (is_ref) {
            //ref_fits.push_back(inclusive_fit); ???
            ref_fit_pt.scale(1.0/SM_xsec_incl); // CHECK SM_xsec_incl vs SM_xsec_sel if making ref pts!!!!!!!
            ref_pts.push_back(ref_fit_pt);
        }

        if (is_tar) {

            std::cout << "Group tag: " << grp_tag << std::endl;
            std::string leg_tag;
            TString process_TStr = process;
            TString tmp_tag = grp_tag;

            // This will set up the names in the legend
            //TString comp_type = "0p1pComp"; // ResultsSec
            TString comp_type = "tag";
            //TString comp_type = "qCutScan";
            //TString comp_type = "matchScaleScan";
            //TString comp_type = "startPtComp";
            //TString comp_type = "ttWttZchecks"; // NLOCompSec

            // Get cleaned up version of process name
            if (process_TStr.Index("ttH") != -1){
                //leg_tag = "tth";
                leg_tag = "t#bar{t}h";
            // Note, do we still want ttll and ttlnu to show up at ttW and ttZ? Should re evalueate next time we want to plot ttll an ttlnu
            } else if (process_TStr.Index("ttZ") != -1 or process_TStr.Index("ttll") != -1) { 
                //leg_tag = "ttZ";
                leg_tag = "t#bar{t}Z";
            } else if (process_TStr.Index("ttW") != -1 or process_TStr.Index("ttlnu") != -1) {
                //leg_tag = "ttW";
                leg_tag = "t#bar{t}W";
            } else {
                std::cout << "Note: process " << process << " is not ttH, ttll, or ttlnu. Not cleaning up process name." << std::endl;
            }
            if (comp_type == "0p1pComp"){
                if (tmp_tag.Index("NoJets") != -1 or process_TStr.Index("Jet") == -1) {
                    //leg_tag = leg_tag + " 0p ";
                    leg_tag = leg_tag + " LO ";
                } else {
                    //leg_tag = leg_tag + " 0+1p  ";
                    leg_tag = leg_tag + "+j LO ";
                }
            } else if (comp_type == "tag"){
                leg_tag = leg_tag + grp_tag;
            } else if (comp_type == "qCutScan") {
                leg_tag = leg_tag + " 0+1p: xqcut10, " + tmp_tag(tmp_tag.Index("qCut"), tmp_tag.Length());
            } else if (comp_type == "matchScaleScan") {
                leg_tag = leg_tag + " 0+1p: " + tmp_tag(tmp_tag.Index("xqcut"),100);
            } else if (comp_type == "startPtComp"){
                leg_tag = tmp_tag + " " + run_label;
            } else if (comp_type == "ttWttZchecks"){
                leg_tag = process;
                cout << "tmp_tag" << tmp_tag << std::endl;
                if (tmp_tag.Index("QED1QCD2") != -1){
                    leg_tag = leg_tag + " qed=1 qcd=2";
                }
                if (tmp_tag.Index("QED1QCD3") != -1){
                    leg_tag = leg_tag + " qed=1";
                }
                if (tmp_tag.Index("QED2") != -1){
                    leg_tag = leg_tag + " qed=2";
                }
                //leg_tag = leg_tag + " " + run_label;
            }

            /* Misc legend settings
            std::cout << "TMP TAG: " << tmp_tag << std::endl;
            //if (tmp_tag.Index("SMEFTcomp") != -1 or tmp_tag.Index("QED1QCD2") != -1){
            if (tmp_tag.Index("dim6TopvMay2020") != -1){
                leg_tag = leg_tag + "dim6TopvMay2020";
            }
            if (tmp_tag.Index("HanV4") != -1){
                leg_tag = leg_tag + "HanV4";
            }
            */

            inclusive_fit.setTag(leg_tag); // If passing inclusive fit!!!
            target_fits.push_back(inclusive_fit);
            //selection_fit.setTag(leg_tag); // If passing selection fit!!!
            //target_fits.push_back(selection_fit);

            /*
            // Make a custom WCFit that is made from points were cpQM is -cpQ3
            std::vector<WCPoint> fitPts_vect;
            for (int i=-5; i<=5; i=i+2){
                std::cout << "FIRST LOOP!!! " << i << std::endl;
                WCPoint evalPt;
                WCPoint fitPt;
                evalPt.setStrength("cpQ3",i);
                evalPt.setStrength("cpQM",-i);
                fitPt.setStrength("cpQ3",i);
                fitPt.wgt = inclusive_fit.evalPoint(&evalPt);
                fitPts_vect.push_back(fitPt);
            }
            for (int i=0; i<fitPts_vect.size(); i++){
                std::cout << "SECOND LOOP: " << i << " " << fitPts_vect.at(i).getDim() << " " << fitPts_vect.at(i).getStrength("cpQ3") << " " << fitPts_vect.at(i).getStrength("cpQM") << " " << fitPts_vect.at(i).wgt << std::endl;
            }
            WCFit* customFit = new WCFit(fitPts_vect,leg_tag+"( cpQM at -cpQ3)");
            target_fits.push_back(*customFit);
            //supplemental_fits.push_back(customFit);
            */
        }

        std::cout << std::endl;
        run_idx++;
    }

    /*
    // TestPrint
    print_xsec(target_fits.at(0));
    std::cout << "\nThe size of the target fits: " << target_fits.size() << "\n" << std::endl;
    for (size_t i=0; i<target_fits.size(); i++){
        print_xsec(target_fits.at(i));
    }
    */

    bool print_stopwatch = 0;
    if (print_stopwatch) {    
        bool use_avg; std::string norm_name;
        use_avg = 1; norm_name = "";
        sw.readAllTimers(use_avg,norm_name);
        use_avg = 0; norm_name = "";
        sw.readAllTimers(use_avg,norm_name);
    }

    // Dynamically figure out which WC are present across all WCFits
    std::vector<std::string> wc_names;
    if (wc_names.size() == 0) {
        std::set<std::string> wc_set;
        for (uint i = 0; i < target_fits.size(); i++) {
            WCFit fit = target_fits.at(i);
            std::vector<std::string> names = fit.getNames();
            for (uint j = 0; j < names.size(); j++) {
                std::string s = names.at(j);
                if (s != "sm" && wc_set.count(s) == 0) {
                    wc_set.insert(s);
                    wc_names.push_back(s);
                }
            }
        }
    }


    // Plot specific 1-D fits
    for (auto& wc_name: wc_names) {
        arxiv_fit_comps_vect = {};
        std::vector<WCFit> subset_fits; // These are the fits we are actually going to plot
        std::vector<customPointInfo> points_to_plot_with_errorbars; // This is what we'll pass to the plotting script
        for (uint i = 0; i < target_fits.size(); i++) {
            if (target_fits.at(i).hasCoefficient(wc_name)) {
                // For comparing 1D fits to each other
                subset_fits.push_back(target_fits.at(i));
                //std::cout << "SIZE !!! " << target_fits.at(i).getDim() << std::endl; 
            }
        }
        if (subset_fits.size() == 0) {
            // No fit has the specified WC (shouldn't really be possible)
            continue;
        }

        xsec_plt_ops_1d.tag = output_name.Data();   // This becomes the save name for the plot
        xsec_plt_ops_1d.tag += "_" + wc_name;

        // Make a string for the cleaned up process name (e.g. ttllNuNuJetNoHiggs -> ttW)
        TString curr_process_clean = curr_process;
        if (curr_process_clean.Index("ttH") != -1){
            curr_process_clean = "tth";
        } else if (curr_process_clean.Index("ttll") != -1) {
            curr_process_clean = "ttZ";
        } else if (curr_process_clean.Index("ttlnu") != -1) {
            curr_process_clean = "ttW";
        }

        // Strip the 'i' from certain WC names (e.g. ctei --> cte) in the plot titles
        if (wc_name.back() == 'i') {
            int len = wc_name.size();
            xsec_plt_ops_1d.title = curr_process + " " + wc_name.substr(0,len-1);
            //xsec_plt_ops_1d.title = curr_process_clean + " " + wc_name.substr(0,len-1); // TMP!!! (for e.g. ttllNuNuJetNoHiggs -> ttW)
        } else {
            xsec_plt_ops_1d.title = curr_process + " " + wc_name;
            //xsec_plt_ops_1d.title = curr_process_clean + " " + wc_name; // TMP!!!
        }
        
        bool save_fits = false; // Don't save the fits for now
        if (save_fits) {
            std::string fitparams_fpath = kOutputDir + "/" + "fitparams_" + curr_process + "_" + wc_name + ".txt";
            make_fitparams_file(fitparams_fpath,subset_fits);
        }

        // Fill the vector of specific quadratic fits to pass to the plot maker (right now just used for NLO fit comps)
        std::pair<std::string,std::vector<double>> LO_pair;
        std::pair<std::string,std::vector<double>> NLO_pair;
        if (curr_process.find("ttH") != std::string::npos){
            std::cout << "Current proc is ttH: " << curr_process << std::endl;
            LO_pair.first  = "1607.05330 LO";
            NLO_pair.first = "1607.05330 NLO";
            if (wc_name=="ctG"){
                LO_pair.second  = {1,1.014,1.390};
                NLO_pair.second = {1,0.991,1.328};
            } else if (wc_name=="ctp"){
                LO_pair.second  = {1,-0.119,0.0035};
                NLO_pair.second = {1,-0.123,0.0037};
            }
            for (auto wc : custom_fits_wc_lst){
                if (wc == wc_name and custom_fits_map["ttH"].find(wc_name) != custom_fits_map["ttH"].end() ){
                    for (auto fit : custom_fits_map["ttH"][wc_name]){
                        //subset_fits.push_back(fit);
                    }
                }
            }
            // Reza smeft
            for (auto points_vects: reza_pts_info_vect){
                for (auto point: points_vects){
                    if (point.proc_name == "tth" and point.wc_name == wc_name){
                        //std::cout << "point.wc: " << point.wc_name << std::endl;
                        points_to_plot_with_errorbars.push_back(point);
                    }
                }
            }

        } else if (curr_process.find("ttlnu") != std::string::npos or curr_process.find("ttW") != std::string::npos){
            std::cout << "Current proc is ttW:" << curr_process << std::endl;
            for (auto wc : custom_fits_wc_lst){
                //if (wc == wc_name and custom_fits_map["ttW"].find(wc_name) != custom_fits_map["ttW"].end() ){
                if (wc == wc_name ){
                    for (auto fit : custom_fits_map["ttW"][wc_name]){
                        //subset_fits.push_back(fit);
                    }
                }
            }
            // Reza smeft
            for (auto points_vects: reza_pts_info_vect){
                for (auto point: points_vects){
                    if (point.proc_name == "ttW" and point.wc_name == wc_name){
                        //std::cout << "point.wc: " << point.wc_name << std::endl;
                        points_to_plot_with_errorbars.push_back(point);
                    }
                }
            }
        } else if (curr_process.find("ttll") != std::string::npos or curr_process.find("ttZ") != std::string::npos){
            std::cout << "Current proc is ttZ: " << curr_process << std::endl;
            LO_pair.first  = "1601.08193 LO";
            NLO_pair.first = "1601.08193 NLO";
            if (wc_name=="ctG"){
                //LO_pair.second  = {1,0.376,0.339}; // ttZ (Tab 3)
                //NLO_pair.second = {1,0.353,0.278}; // ttZ (Tab 3)
                LO_pair.second  = {1,0.37634549750590707, 0.33932790758729325}; // ttZ more precision (Tab 3)
                NLO_pair.second = {1,0.3532423208191126, 0.2781569965870307}; // ttZ more precision (Tab 3)
            } else if (wc_name=="cpQ3"){
                //LO_pair.second  = {1,0.103,0.00368}; // ttZ (Tab 3)
                //NLO_pair.second = {1,0.103,0.00432}; // ttZ (Tab 3)
                LO_pair.second  = {1,0.10278288264636387, 0.0036755053819900237}; // ttZ more precision (Tab 3)
                NLO_pair.second = {1,0.1030716723549488, 0.004323094425483504}; // ttZ more precision (Tab 3)
            } else if (wc_name=="cpt"){
                //LO_pair.second  = {1,0.0677,0.00381}; // ttZ (Tab 3)
                //NLO_pair.second = {1,0.0654,0.00444}; // ttZ (Tab 3)
                LO_pair.second  = {1,0.06773431346810187, 0.00380677343134681}; // ttZ more precision (Tab 3)
                NLO_pair.second = {1,0.06541524459613197, 0.004436860068259386}; // ttZ more precision (Tab 3)
            } else if (wc_name=="ctW"){
                //LO_pair.second   = {1,-0.000263,0.0274}; // ttZ (Tab 3)
                //NLO_pair.second  = {1,-0.00193,0.0275};  // ttZ (Tab 3)
                LO_pair.second  = {1,-0.00026253609871357313, 0.02743502231556839}; // ttZ more precision (Tab 3)
                NLO_pair.second = {1,-0.0019340159271899885, 0.027531285551763367}; // ttZ more precision (Tab 3)
            }
            for (auto wc : custom_fits_wc_lst){
                if (wc == wc_name and custom_fits_map["ttZ"].find(wc_name) != custom_fits_map["ttZ"].end()){
                    for (auto fit : custom_fits_map["ttZ"][wc_name]){
                        //subset_fits.push_back(fit);
                    }
                }
            }
            // Reza smeft
            for (auto points_vects: reza_pts_info_vect){
                for (auto point: points_vects){
                    if (point.proc_name == "ttZ" and point.wc_name == wc_name){
                        //std::cout << "point.wc: " << point.wc_name << std::endl;
                        points_to_plot_with_errorbars.push_back(point);
                    }
                }
            }

        }
        //arxiv_fit_comps_vect.push_back(NLO_pair); // Comment out to not plot arxiv comps
        //arxiv_fit_comps_vect.push_back(LO_pair);  // Comment out to not plot arxiv comps

        ///*
        make_1d_xsec_plot(
            xsec_plt_ops_1d,
            wc_name,
            subset_fits,
            ref_pts,
            arxiv_fit_comps_vect,
            points_to_plot_with_errorbars,
            curr_process
        );
        //*/
    }
    std::cout << "Finished!" << std::endl;
}

void runGridpackValidation(TString output_name,TString input_rundirs_spec,TString ref_rundirs_spec,TString grp_name="") {
    //gStyle->SetPadRightMargin(0.2);
    gStyle->SetOptStat(0);

    runit(output_name,input_rundirs_spec,ref_rundirs_spec,grp_name);
}

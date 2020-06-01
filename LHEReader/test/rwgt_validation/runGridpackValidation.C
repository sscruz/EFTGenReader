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
#include "makeEFTPlots.h"

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

void print_xsec(WCFit wcfit){

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

}

void runit(TString output_name,TString input_rundirs_spec,TString ref_rundirs_spec,TString grp_name="") {
    int run_idx = 0;

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

    // ttH
    double sm_ttH_nlo = 0.4521;
    ptsMap["ttH"]["ctG_nlo"] = { {-1,0.5731} , {0,sm_ttH_nlo} , {1,1.482} };
    ptsMap["ttH"]["ctp_nlo"] = { {-5,0.7925} , {0,sm_ttH_nlo} , {5,0.2051} };
    // ttW
    double sm_ttW_nlo = 0.5591;
    ptsMap["ttW"]["ctG_nlo"] = { {-1,0.4192} , {0,sm_ttW_nlo} , {5,0.7272} };
    // ttZ
    double sm_ttZ_nlo = 0.7441;
    ptsMap["ttZ"]["ctG_nlo"]   = { {-1,0.6725} , {0,sm_ttZ_nlo} , {1,1.206} };
    ptsMap["ttZ"]["cpQM_nlo"]  = { {-5,1.208}  , {0,sm_ttZ_nlo} , {5,0.4392} };
    ptsMap["ttZ"]["cpt_nlo"]   = { {-5,0.5757} , {0,sm_ttZ_nlo} , {5,1.068} };
    ptsMap["ttZ"]["ctZ_nlo"]   = { {-5,2.261}  , {0,sm_ttZ_nlo} , {5,2.238} };

    double sm_ttZ_lo = 0.5841;
    ptsMap["ttZ"]["cpt_lo"]   = { {-5,0.4395}  , {0,sm_ttZ_lo} , {5,0.8523} };
    ptsMap["ttZ"]["ctp_lo"]   = { {-5,0.5855}  , {0,sm_ttZ_lo} , {5,0.586} };
    ptsMap["ttZ"]["ctZ_lo"]   = { {-5,1.903}  , {0,sm_ttZ_lo} , {5,1.896} };
    ptsMap["ttZ"]["ctW_lo"]   = { {-5,0.5842}  , {0,sm_ttZ_lo} , {5,0.5846} };

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
        last_entry = 100; // For testing
        std::cout << "Last_entry: " << last_entry << std::endl;

        int first_entry = 0;

        WCFit* wcFit_intree = 0;
        double originalXWGTUP_intree = -1.;

        chain.SetBranchAddress("originalXWGTUP",&originalXWGTUP_intree);
        chain.SetBranchAddress("wcFit",&wcFit_intree);

        bool skip_progress = false;

        chain.GetEntry(first_entry);
        if (is_ref) {
            skip_progress = true;
        }

        WCFit inclusive_fit;
        WCFit selection_fit;

        int selection_events = 0;
        double genLep_pt1_intree;
        double genLep_pt2_intree;
        double genLep_pt3_intree;
        double genJet_pt4_intree;
        chain.SetBranchAddress("genLep_pt1",&genLep_pt1_intree);
        chain.SetBranchAddress("genLep_pt2",&genLep_pt2_intree);
        chain.SetBranchAddress("genLep_pt3",&genLep_pt3_intree);
        chain.SetBranchAddress("genJet_pt4",&genJet_pt4_intree);

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

        sw.start("Full Loop");
        for (int i = first_entry; i < last_entry; i++) {
            sw.start("Event Loop");
            if (is_tar && !skip_progress) {
                printProgress(i - first_entry,last_entry - first_entry);
            }
            sw.start("Get Entry");
            chain.GetEntry(i);
            sw.lap("Get Entry");

            sw.start("Add Fit");
            inclusive_fit.addFit(*wcFit_intree);
            sw.lap("Add Fit");
            sw.lap("Event Loop");

            ref_fit_pt.wgt += originalXWGTUP_intree;

            //if (genLep_pt3_intree > 10 and genJet_pt4_intree > 100){
            if (genLep_pt3_intree > 10 and genJet_pt4_intree > 30){
                selection_events = selection_events + 1;
                selection_fit.addFit(*wcFit_intree);
            } else{
                //std::cout << "Skip; pt: " << genLep_pt3_intree << " , " << genJet_pt4_intree << std::endl;
            }

        }
        sw.stop("Full Loop");

        std::cout << "\n Selected events over total: " << selection_events << "/" << last_entry << "->" << (float)selection_events/last_entry << "\n" << std::endl;

        // Normalize to SM
        //std::cout << "\nBefore all norm!!! incl SM xsec: " << inclusive_fit.evalPoint(&sm_pt) << " selection SM xsec: " << selection_fit.evalPoint(&sm_pt) << std::endl;

        double SM_xsec_incl = inclusive_fit.evalPoint(&sm_pt);
        double SM_xsec_sel = selection_fit.evalPoint(&sm_pt);
        inclusive_fit.scale(1.0/SM_xsec_incl);
        selection_fit.scale(1.0/SM_xsec_sel);

        //std::cout << "\nAfter all norm!!! incl SM xsec: " << inclusive_fit.evalPoint(&sm_pt) << " selection SM xsec: " << selection_fit.evalPoint(&sm_pt) << std::endl;

        /* // Dump the fit functions
        //std::vector<std::string> list_of_WC = {"ctG","ctW"};
        std::vector<std::string> list_of_WC = {"ctp","cpQM","ctW","ctZ","ctG","cbW","cpQ3","cptb","cpt","cQl3i","cQlMi","cQei","ctli","ctei","ctlSi","ctlTi"};
        std::cout << " " << std:: endl;
        for (std::string WC : list_of_WC){ 
            inclusive_fit.dump(false,153,WC);
            //selection_fit.dump(false,153,WC);
            std::cout << " " << std:: endl;
        }
        */

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
            //TString comp_type = "0p1pComp";
            //TString comp_type = "qCutScan";
            //TString comp_type = "matchScaleScan";
            //TString comp_type = "startPtComp";
            TString comp_type = "ttWttZchecks";

            // Get cleaned up version of process name
            if (process_TStr.Index("ttH") != -1){
                leg_tag = "tth";
            } else if (process_TStr.Index("ttll") != -1) {
                leg_tag = "ttZ";
            } else if (process_TStr.Index("ttlnu") != -1) {
                leg_tag = "ttW";
            } else {
                std::cout << "Note: process " << process << " is not ttH, ttll, or ttlnu. Not cleaning up process name." << std::endl;
            }
            if (comp_type == "0p1pComp"){
                if (tmp_tag.Index("NoJets") != -1 or process_TStr.Index("Jet") == -1) {
                    leg_tag = leg_tag + " 0p ";
                } else {
                    leg_tag = leg_tag + " 0+1p  ";
                }
            } else if (comp_type == "qCutScan") {
                leg_tag = leg_tag + " 0+1p: xqcut10, " + tmp_tag(tmp_tag.Index("qCut"), tmp_tag.Length());
            } else if (comp_type == "matchScaleScan") {
                leg_tag = leg_tag + " 0+1p: " + tmp_tag(tmp_tag.Index("xqcut"),100);
            } else if (comp_type == "startPtComp"){
                leg_tag = tmp_tag + " " + run_label;
            } else if (comp_type == "ttWttZchecks"){
                leg_tag = process;
                if (tmp_tag.Index("QED1QCD2") != -1){
                    leg_tag = leg_tag + " qed=1 qcd=2";
                }
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

    //print_xsec(target_fits.at(0));

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
                        subset_fits.push_back(fit);
                    }
                }
            }
        } else if (curr_process.find("ttlnu") != std::string::npos or curr_process.find("ttW") != std::string::npos){
            std::cout << "Current proc is ttW:" << curr_process << std::endl;
            for (auto wc : custom_fits_wc_lst){
                //if (wc == wc_name and custom_fits_map["ttW"].find(wc_name) != custom_fits_map["ttW"].end() ){
                if (wc == wc_name ){
                    for (auto fit : custom_fits_map["ttW"][wc_name]){
                        subset_fits.push_back(fit);
                    }
                }
            }
        } else if (curr_process.find("ttll") != std::string::npos or curr_process.find("ttZ") != std::string::npos){
            std::cout << "Current proc is ttZ: " << curr_process << std::endl;
            LO_pair.first  = "1601.08193 LO";
            NLO_pair.first = "1601.08193 NLO";
            //LO_pair.first  = "1601.08193 LO (mumu)";
            //NLO_pair.first = "1601.08193 NLO (mumu)";
            if (wc_name=="ctG"){
                LO_pair.second  = {1,0.376,0.339}; // ttZ (Tab 3)
                NLO_pair.second = {1,0.353,0.278}; // ttZ (Tab 3)
                //LO_pair.second  = {1,0.356,0.283}; // ttmu+mu- (Tab 6)
                //NLO_pair.second = {1,0.335,0.226}; // ttmu+mu- (Tab 6)
            } else if (wc_name=="cpQ3"){
                LO_pair.second  = {1,0.103,0.00368}; // ttZ (Tab 3)
                NLO_pair.second = {1,0.103,0.00432}; // ttZ (Tab 3)
                //LO_pair.second  = {1,0.0816,0.00319}; // ttmu+mu- (Tab 6)
                //NLO_pair.second = {1,0.0793,0.00311}; // ttmu+mu- (Tab 6)
            } else if (wc_name=="cpt"){
                LO_pair.second  = {1,0.0677,0.00381}; // ttZ (Tab 3)
                NLO_pair.second = {1,0.0654,0.00444}; // ttZ (Tab 3)
                //LO_pair.second  = {1,0.0537,0.00315}; // ttmu+mu- (Tab 6)
                //NLO_pair.second = {1,0.0504,0.00299}; // ttmu+mu- (Tab 6)
            } else if (wc_name=="ctW"){
                LO_pair.second   = {1,-0.000263,0.0274}; // ttZ (Tab 3)
                NLO_pair.second  = {1,-0.00193,0.0275};  // ttZ (Tab 3)
                //LO_pair.second  = {1,0.000789,0.0235}; // ttmu+mu- (Tab 6)
                //NLO_pair.second = {1,-0.00112,0.0227}; // ttmu+mu- (Tab 6)
            }
            for (auto wc : custom_fits_wc_lst){
                if (wc == wc_name and custom_fits_map["ttZ"].find(wc_name) != custom_fits_map["ttZ"].end()){
                    for (auto fit : custom_fits_map["ttZ"][wc_name]){
                        subset_fits.push_back(fit);
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
            arxiv_fit_comps_vect
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

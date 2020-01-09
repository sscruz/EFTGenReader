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

//#include "EFTGenReader/GenReader/interface/WCPoint.h"
//#include "EFTGenReader/GenReader/interface/WCFit.h"
//#include "EFTGenReader/GenReader/interface/TH1EFT.h"
//#include "EFTGenReader/GenReader/interface/Stopwatch.h"
//#include "EFTGenReader/GenReader/interface/split_string.h"
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

void runit(TString output_name,TString input_rundirs_spec,TString ref_rundirs_spec,TString grp_name="") {
    int run_idx = 0;

    Stopwatch sw;

    std::vector<WCFit> target_fits;
    std::vector<WCFit> ref_fits;
    std::vector<WCPoint> ref_pts;

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

        // Set up the wc point string (depends a lot on the naming scheme)
        std::map<string,string> ref_pts_dict;
        std::string wcname = "ctG";
        int range_max = 3;
        int run = 0;
        for (int wcval=-range_max; wcval<=range_max; wcval++){
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
        }
        sw.stop("Full Loop");

        // Normalize to SM
        double SM_xsec = inclusive_fit.evalPoint(&sm_pt);
        inclusive_fit.scale(1.0/SM_xsec);

        // Normalize ref pt and add to list
        std::cout << "Is ref? " << is_ref << std::endl;
        if (is_ref) {
            //ref_fits.push_back(inclusive_fit); ???
            ref_fit_pt.scale(1.0/SM_xsec);
            ref_pts.push_back(ref_fit_pt);
        }

        if (is_tar) {

            //std::string fit_tag;
            //fit_tag = grp_tag;
            //fit_tag = process + " " + fit_tag; // TMP!!!

            std::cout << "group tag: " << grp_tag << std::endl;
            std::string leg_tag;
            TString process_TStr = process;
            TString tmp_tag = grp_tag;

            // Get cleaned up version of process name
            //TString comp_type = "0p1pComp";
            //TString comp_type = "qCutScan";
            TString comp_type = "matchScaleScan";
            //TString comp_type = "startPtComp";
            if (process_TStr.Index("ttH") != -1){
                leg_tag = "tth";
            } else if (process_TStr.Index("ttll") != -1) {
                leg_tag = "ttZ";
            } else if (process_TStr.Index("ttlnu") != -1) {
                leg_tag = "ttW";
            } else {
                std::cout << "Warning: process " << process << " is not ttH, ttll, or ttlnu" << std::endl;
            }
            if (comp_type == "0p1pComp"){
                if (tmp_tag.Index("NoJets") != -1) {
                    leg_tag = leg_tag + " 0p";
                } else {
                    leg_tag = leg_tag + " 0+1p  ";
                }
            } else if (comp_type == "qCutScan") {
                leg_tag = leg_tag + " 0+1p: xqcut10, " + tmp_tag(tmp_tag.Index("qCut"), tmp_tag.Length());
            } else if (comp_type == "matchScaleScan") {
                leg_tag = leg_tag + " 0+1p: " + tmp_tag(tmp_tag.Index("xqcut"),100);
            } else if (comp_type == "startPtComp"){
                leg_tag = tmp_tag + " " + run_label;
            }
            inclusive_fit.setTag(leg_tag);

            //fit_tag.erase(0,25);
            //if (fit_tag[0] != 'x') {
            //    fit_tag = "xqcut10" + fit_tag;
            //}
            //fit_tag = "ttHJet_" + fit_tag;
            //inclusive_fit.setTag(fit_tag);

            target_fits.push_back(inclusive_fit);
        }

        std::cout << std::endl;
        run_idx++;
    }


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
        std::vector<WCFit> subset_fits; // These are the fits we are actually going to plot
        for (uint i = 0; i < target_fits.size(); i++) {
            if (target_fits.at(i).hasCoefficient(wc_name)) {
                // For comparing 1D fits to each other
                subset_fits.push_back(target_fits.at(i));
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
        //xsec_plt_ops_1d.title = ""; // NO TITLE!!!
        
        bool save_fits = false; // Don't save the fits for now
        if (save_fits) {
            std::string fitparams_fpath = kOutputDir + "/" + "fitparams_" + curr_process + "_" + wc_name + ".txt";
            make_fitparams_file(fitparams_fpath,subset_fits);
        }

        make_1d_xsec_plot(
            xsec_plt_ops_1d,
            wc_name,
            subset_fits,
            ref_pts
        );
    }

    std::cout << "Finished!" << std::endl;
}

void runGridpackValidation(TString output_name,TString input_rundirs_spec,TString ref_rundirs_spec,TString grp_name="") {
    //gStyle->SetPadRightMargin(0.2);
    gStyle->SetOptStat(0);

    runit(output_name,input_rundirs_spec,ref_rundirs_spec,grp_name);
}

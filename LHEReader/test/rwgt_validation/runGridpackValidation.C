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

#include "EFTGenReader/GenReader/interface/WCPoint.h"
#include "EFTGenReader/GenReader/interface/WCFit.h"
#include "EFTGenReader/GenReader/interface/TH1EFT.h"
#include "EFTGenReader/GenReader/interface/split_string.h"
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

    std::vector<WCFit> target_fits;
    std::vector<WCPoint> ref_pts;

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

        std::string process   = words.at(1);
        std::string grp_tag   = words.at(2);
        std::string run_label = words.at(3);

        // Chain together all root files in the run directory
        TChain chain("EFTLHEReader/summaryTree");
        auto dir_files = getFiles(fdir);
        for (auto fn: dir_files) {
            //std::cout << "\tFiles: " << fn << std::endl;
            TString fname = fdir + "/" + fn;
            chain.Add(fname);
        }

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

        int chain_entries = chain.GetEntries();
        int last_entry  = chain_entries;
        int first_entry = 0;

        std::unordered_map<std::string,double> *eftwgts_intree = 0;
        double originalXWGTUP_intree = -1.;
        int lumiBlock_intree = -1;

        chain.SetBranchAddress("eftwgts",&eftwgts_intree);
        chain.SetBranchAddress("originalXWGTUP",&originalXWGTUP_intree);
        chain.SetBranchAddress("lumiBlock",&lumiBlock_intree);

        bool skip_progress = false;

        chain.GetEntry(first_entry);
        if (eftwgts_intree->size() < 20) {
            skip_progress = true;
        }

        std::set<int> unique_runs;
        std::unordered_map<std::string,WCPoint> inclusive_xsec_wgts;    // Keys are the rwgt string ids
        inclusive_xsec_wgts[kOrig] = start_pt;
        for (int i = first_entry; i < last_entry; i++) {
            if (is_tar && !skip_progress) {
                printProgress(i - first_entry,last_entry - first_entry);
            }
            chain.GetEntry(i);
            unique_runs.insert(lumiBlock_intree);
            inclusive_xsec_wgts[kOrig].wgt += originalXWGTUP_intree;
            // Note: Commenting this out can make generating plots with reference points really slow
            //       will need to be commented out if we want to normalize reference points to SM
            if (!is_tar) {
                // This gridpack run is not being validated (i.e. We don't care about its rwgt points)
                continue;
            }
            for (auto& kv: *eftwgts_intree) {
                if (inclusive_xsec_wgts.find(kv.first) == inclusive_xsec_wgts.end()) {
                    inclusive_xsec_wgts[kv.first] = WCPoint(kv.first,0.0);
                }
                inclusive_xsec_wgts[kv.first].wgt += kv.second;
            }
        }

        // For now skip normalizing to SM
        double SM_xsec = -1.0;
        //for (auto& kv: inclusive_xsec_wgts) {
        //    if (inclusive_xsec_wgts[kv.first].isSMPoint()) {
        //        // Need to normalize, otherwise we will
        //        SM_xsec = inclusive_xsec_wgts[kv.first].wgt;
        //        break;
        //    }
        //}

        double gridpack_scale = 1.0 / unique_runs.size();
        double xsec_norm = 1.0 / SM_xsec;
        std::vector<WCPoint> fit_pts;
        for (auto& kv: inclusive_xsec_wgts) {
            if (SM_xsec < 0.0) {
                inclusive_xsec_wgts[kv.first].scale(gridpack_scale);
            } else {
                inclusive_xsec_wgts[kv.first].scale(xsec_norm);
            }
            fit_pts.push_back(inclusive_xsec_wgts[kv.first]);
        }

        std::cout << "\tDist: " << start_pt.getEuclideanDistance() << std::endl;
        std::cout << "\tkOrig wgt: " << inclusive_xsec_wgts[kOrig].wgt << std::endl;
        if (SM_xsec < 0.0) {
            std::cout << "\tSM Xsec: " << SM_xsec << std::endl;
        } else {
            std::cout << "\tSM Xsec: " << SM_xsec / unique_runs.size() << std::endl;
        }
        //for (auto& kv: inclusive_xsec_wgts) {
        //    if (kv.first == kOrig) {
        //        continue;
        //    }
        //    std::cout << "\t" << kv.first << ": " << kv.second.wgt << std::endl;
        //}

        if (is_ref) {
            // Only include specific runs to use as reference
            ref_pts.push_back(inclusive_xsec_wgts[kOrig]);
        }

        if (is_tar) {
            std::string fit_tag;    // The plots will have the process name in the title
            if (grp_name.Length() > 0) {
                // Overwrite the group tag with a custom one
                fit_tag = grp_name.Data();
                fit_tag += "_" + run_label;
                //fit_tag = process + "_" + grp_name.Data() + "_" + run_label;
            } else {
                // Default behaviour
                fit_tag = grp_tag + "_" + run_label;
            }

            if (fit_pts.size() <= 1) {
                std::cout << "Error: Not enough fit points for " << fit_tag << std::endl;
                continue;
            }
            WCFit inclusive_fit(fit_pts,fit_tag);
            inclusive_fit.setStart(inclusive_xsec_wgts[kOrig]);
            target_fits.push_back(inclusive_fit);
        }

        std::cout << std::endl;
        run_idx++;
    }

    // Dynamically figure out which WC are present    
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
            std::vector<std::string> words;
            split_string(target_fits.at(i).getTag(),words,"_");
            //if (target_fits.at(i).hasCoefficient(wc_name) && target_fits.at(i).getDim() > 1) {
            //    // We only want to plot the multi-dim fits, 1-D fits are assumed to be for reference points
            //    subset_fits.push_back(target_fits.at(i));
            //}
            if (target_fits.at(i).hasCoefficient(wc_name)) {
                // For comparing 1D fits to each other
                subset_fits.push_back(target_fits.at(i));
            }
        }
        if (subset_fits.size() == 0) {
            continue;
        }

        //for (uint j = 0; j < subset_fits.size(); j++) {
        //    std::string save_name = kOutputDir + "/" + "fitparams_" + curr_process + "_" + wc_name + ".txt";
        //    subset_fits.at(j).save(save_name,j);
        //}

        xsec_plt_ops_1d.tag = output_name.Data();   // This becomes the save name for the plot
        xsec_plt_ops_1d.tag += "_" + wc_name;
        //xsec_plt_ops_1d.title = xsec_plt_ops_1d.tag;
        //xsec_plt_ops_1d.title = curr_process + " " + wc_name + " Scan";
        // Strip the 'i' from certain WC names (e.g. ctei --> cte)
        if (wc_name.back() == 'i') {
            int len = wc_name.size();
            xsec_plt_ops_1d.title = curr_process + " " + wc_name.substr(0,len-1);
        } else {
            xsec_plt_ops_1d.title = curr_process + " " + wc_name;
        }
        
        std::string fitparams_fpath = kOutputDir + "/" + "fitparams_" + curr_process + "_" + wc_name + ".txt";
        //make_fitparams_file(fitparams_fpath,subset_fits);

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
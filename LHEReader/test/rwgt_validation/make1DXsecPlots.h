#ifndef MAKEEFTPLOTS_H_
#define MAKEEFTPLOTS_H_

#include <string>
#include <vector>
#include <algorithm>

#include "EFTGenReader/EFTHelperUtilities/interface/WCPoint.h"
#include "EFTGenReader/EFTHelperUtilities/interface/WCFit.h"

#include "TString.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TH1D.h"
#include "TMath.h"
#include "TF1.h"
#include "TF2.h"
#include "TGraph.h"
#include "TLatex.h"

#include "TGraphErrors.h"

// For storing custom points and fits, why not just use WCPoint?
struct customPointInfo {
    TString proc_name;
    TString type;
    TString wc_name;
    double x;
    double y;
    double err;
};

// Print info about customPointsInfo object
void print_customPointInfo_object(std::vector<customPointInfo> points_info){
    for (int i=0; i<points_info.size(); i++){
        customPointInfo s = points_info.at(i);
        std::cout << "\nPrinting info about customPointInfo struct:" << std::endl;
        std::cout << "\tProc name: " << s.proc_name << ", Type: " << s.type << ", WC: " << s.wc_name << ", x val: " << s.x << ", y val: " << s.y << ", error: " << s.err << std::endl;
    }
}


class PlotOptions
{
private:
    const int kNColors = 8;
    //const int clr_map[8] = {12,46,9,30,41,4,6,8};
    //const int clr_map[8] = {46,12,9,30,41,4,6,8};
    //const int clr_map[8] = {46,4,13,3,46,13}; // Was a good scheme for the cpQM cpQ3 plot NLO comp plot
    const int clr_map[8] = {46,46,8,8,4,4}; // Want private 0p and 0+1p to have same color, smefit LO and NLO to have same color etc
    int clr_idx = 0;
public:
    PlotOptions(){};
    ~PlotOptions(){};

    std::string tag;        // Names the plot options and is used for save file name
    std::string output_dir;
    std::string x_name,y_name,z_name;
    std::string title;
    double x_min,x_max;
    double y_min,y_max;
    double z_min,z_max;
    int nbins;

    // Returns the next color in the color map based on internal clr_idx
    int nextColor() {
        int clr = this->clr_map[this->clr_idx % this->kNColors];
        clr_idx++;
        return clr;
    }

    // Returns the color corresponding to the current internal clr_idx
    int currColor() {
        return this->clr_map[clr_idx % this->kNColors];
    }

    // Returns the clr_idx corresponding to the supplied color
    int findColor(int clr) {
        for (int idx = 0; idx < this->kNColors; idx++) {
            if (this->clr_map[idx] == clr) {
                return idx;
            }
        }
        return -1;
    }

    // Returns the color at the specified color index
    int getColor(int _idx) {
        return this->clr_map[_idx % this->kNColors];
    }

    void updateXLimits(double _low,double _high) {
        this->x_min = std::min(this->x_min,_low);
        this->x_max = std::max(this->x_max,_high);
    }

    void updateYLimits(double _low,double _high) {
        this->y_min = std::min(this->y_min,_low);
        this->y_max = std::max(this->y_max,_high);
    }

    void setXLimits(double _low,double _high) {
        this->x_min = _low;
        this->x_max = _high;
    }

    void setYLimits(double _low,double _high) {
        this->y_min = _low;
        this->y_max = _high;
    }
};

void printProgress(int current_index, int total_entries, int interval=20) {
    if (current_index % max(int(total_entries*interval/100.),interval) == 0) {
        float fraction = 100.*current_index/total_entries;
        std::cout << int(fraction) << " % processed " << std::endl;
    }
}

// Returns a histogram that is the ratio of h1/h2
template<typename T>
TH1D* makeRatioHistogram(TString name,T* h1,T* h2) {
    if (h1->GetXaxis()->GetNbins() != h2->GetXaxis()->GetNbins()) {
        std::cout << "[Error] makeRatioHistogram() - bin mismatch between ratio of hists" << std::endl;
        throw;
    }
    TAxis* xaxis = h1->GetXaxis();
    int bins = xaxis->GetNbins();
    Double_t low = xaxis->GetXmin();
    Double_t high = xaxis->GetXmax();

    TAxis* yaxis = h1->GetYaxis();
    Double_t yaxis_sz = yaxis->GetLabelSize();

    // std::cout << "Ratio Name: " << name << std::endl;
    TH1D* h_ratio = new TH1D(name,"",bins,low,high);    // Make sure the title is empty!
    h_ratio->GetYaxis()->SetLabelSize(yaxis_sz*1.5);
    h_ratio->GetXaxis()->SetLabelSize(yaxis_sz*2.5);    // Note: If this is done after alphanumeric, it does nothing
    for (int i = 1; i <= bins; i++) {
        if (xaxis->IsAlphanumeric()) {
            TString bin_label = xaxis->GetBinLabel(i);
            h_ratio->GetXaxis()->SetBinLabel(i,bin_label);
        }
        double val1 = h1->GetBinContent(i);
        double err1 = h1->GetBinError(i);

        double val2 = h2->GetBinContent(i);
        double err2 = h2->GetBinError(i);

        double ratio;
        double ratio_err;
        if (val2 != 0) {
            ratio = abs(val1 / val2);
        } else if (val1 == 0 && val2 == 0) {
            ratio = 1.0;
        } else {
            ratio = -1;
        }
        ratio_err = 0.0;    // Ignore the error for now
        h_ratio->SetBinContent(i,ratio);
        h_ratio->SetBinError(i,ratio_err);
    }
    //h_ratio->GetYaxis()->SetRangeUser(0.0,2.3);
    return h_ratio;
}

// Returns a vector of file names inside of the specified directory
vector<TString> getFiles(const char *dirname=".", const char *ext=".root") {
    vector<TString> dir_files;
    TSystemDirectory dir(dirname,dirname);
    TList *files = dir.GetListOfFiles();
    if (files) {
        TSystemFile *file;
        TString fname;
        TIter next(files);
        while ((file=(TSystemFile*)next())) {
            fname = file->GetName();
            if (!file->IsDirectory() && fname.EndsWith(ext)) {
                dir_files.push_back(fname);
            }
        }
    }
    delete files;
    return dir_files;
}

// Returns the directory name for the coefficient scan
// e.g: /hadoop/store/user/awightma/summaryTree_LHE/v1/output_ttH_ctlSI1_run2 --> output_ttH_ctlSI1_run2
std::string getRunDirectory(std::string str) {
    std::vector<std::string> words;
    split_string(str,words,"/");
    if (words.size() > 0) {
        return words.back();
    } else {
        return "";
    }
}

// Attempts to construct path to the scanpoints directory for a particular gridpack run
std::string getScanPointsDirectory(std::string summary_tree_path) {
    //NOTE: Assumes that the gridpack sub-directory is encoded into the summary tree directory path
    //NOTE: This won't work if the directory structure underneath 'grp_tag' directory is changed
    // Example Input: /hadoop/store/user/awightma/summaryTree_LHE/2018_05_06/ctW1dim/v1/output_ttH_ctW1dim_run6
    // Example Input: /hadoop/store/user/awightma/summaryTree_LHE/reference_scans/HanModel_1jet/ttH-tllq4f-tHq4f-ttlnu-GEN/v1/output_tllq4fMatchedNoHiggs_cpQ3HanModel1DRef_run0
    std::string base_dir = "/hadoop/store/user/awightma/gridpack_scans/";

    std::vector<std::string> words;
    split_string(summary_tree_path,words,"/");
    if (words.size() < 6) {
        std::cout << "[ERROR] Unable to parse summary tree path: " << summary_tree_path << std::endl;
        return "";
    } else if (words.at(5) != "summaryTree_LHE") {
        std::cout << "[ERROR] Invalid summary tree path (bad name): " << summary_tree_path << std::endl;
        return "";
    }

    uint start_idx = 0;
    uint end_idx = words.size() - 4;
    for (uint i=0; i < words.size(); i++) {
        if (words.at(i) == "summaryTree_LHE") {
            start_idx = i+1;
            break;
        }
    }

    if (start_idx > end_idx) {
        std::cout << "[ERROR] Invalid summary tree path (bad idx): " << summary_tree_path << std::endl;
    }

    std::string sub_path = "";
    for (uint i=start_idx; i <= end_idx; i++) {
        sub_path += words.at(i);
        if (i != end_idx) {
            sub_path += "/";
        }
    }

    //std::string scanpoints_dir = base_dir + words.at(6) + "/scanpoints/";
    std::string scanpoints_dir = base_dir + sub_path + "/scanpoints/";
    return scanpoints_dir;
}

// Reads a scanpoints file and returns a vector of WC points
std::vector<WCPoint> parseScanPointsFile(std::string fpath) {
    std::string line,header;
    std::vector<std::string> coeffs,words;
    std::vector<WCPoint> wc_pts;
    
    ifstream inf(fpath);
    if (!inf) {
        std::cout << "[ERROR] Unable to open file: " << fpath << std::endl;
        return wc_pts;
    }

    std::getline(inf,header);
    split_string(header,words," ");
    for (auto s: words) {
        if (s.size() == 0) {
            continue;
        }
        coeffs.push_back(s);
    }

    while (!inf.eof()) {
        std::getline(inf,line);
        words.clear();
        split_string(line,words," ");
        std::string pt_name;
        std::vector<double> strengths;
        for (uint i = 0; i < words.size(); i++) {
            std::string s = words.at(i);
            if (s.size() == 0) {
                continue;
            } else if (i == 0) {
                pt_name = words.at(i);
                continue;
            }
            strengths.push_back(std::stod(s));
        }
        if (strengths.size() != coeffs.size()) {
            std::cout << "[WANRING] Failed to parse line in scanpoints file, " << pt_name << std::endl;
            continue;
        }
        WCPoint wc_pt;
        wc_pt.tag = pt_name;
        for (uint i = 0; i < coeffs.size(); i++) {
            wc_pt.setStrength(coeffs.at(i),strengths.at(i));
        }
        wc_pts.push_back(wc_pt);
    }
    inf.close();
    return wc_pts;
}

bool hasElement(std::vector<TString> v, TString chk) {
    bool exists = false;
    for (auto& e: v) {
        exists = (e == chk);
        if (exists) break;
    }
    return exists;
}

// Reformats a 1-D histogram to have equal bin width on log scale
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

void padR(std::string& str,std::size_t pwidth,const char pchar=' ') {
    if (str.size() >= pwidth) {
        int i1 = pwidth - 3;
        int len = str.size() - i1;
        str.replace(i1,len,"...");
        return;
    }
    str.append(pwidth - str.size(),pchar);
}

std::vector<WCPoint> sortByStrength(vector<WCPoint> pts,std::string wc_name) {
    vector<WCPoint> sorted_pts(pts.begin(),pts.end());
    std::sort(sorted_pts.begin(),sorted_pts.end(),[&wc_name] (WCPoint a, WCPoint b) {return a.getStrength(wc_name) < b.getStrength(wc_name);});
    return sorted_pts;
}

std::vector<WCPoint> sortByWeight(vector<WCPoint> pts) {
    vector<WCPoint> sorted_pts(pts.begin(),pts.end());
    std::sort(sorted_pts.begin(),sorted_pts.end(),[] (WCPoint a, WCPoint b) {return a.wgt < b.wgt;});
    return sorted_pts;
}

// Returns the minimum of the WCFit for a particular WC (setting all others to zero)
double get1DMinimum(WCFit fit,std::string wc_name) {
    double c1 = fit.getCoefficient("sm",wc_name);
    double c2 = fit.getCoefficient(wc_name,wc_name);
    if (c2 == 0.0) {
        // The fit doesn't have a minimum w.r.t this WC
        return 0.0;
    }
    return -c1/(2.0*c2);
}

// Plots the inclusive cross-section as function of a particular WC parameter strength
// Note: quad_fit_vect is a vector of maps of the form [('name',[s0,s1,s2])]
void make_1d_xsec_plot(
    PlotOptions plt_ops,
    std::string wc_name,
    std::vector<WCFit> wc_fits,
    std::vector<WCPoint> ref_pts = {},   // orig_pts
    std::vector<std::pair<std::string,std::vector<double>>> quad_fit_vect = {},
    std::vector<customPointInfo> points_to_plot_with_errorbars = {},
    std::string proc_name = ""
) {

    //// Setup the plots ////

    // Colors
    // TMP: for making the private/arxiv/smeft LO vs NLO plots
    // comment out the plot ops colors for this and also note the colors would be messed up if we tried to plot ratio plots
    // Also note the error bands colors would be messed up so don't plot those either
    //int arxiv_color = 46;          // Dark red
    int arxiv_color = 6;          // Pink
    int private_color = 4;         // Blue
    int smeft_color = 2;           // Red
    int QEDorder_comp_color = 3;   // Bright green
    int QEDorder_2_comp_color = 7;   // Cyan
    int dim6TopvMay2020_color = 2; // Red
    int ttWttZ_color = 4;          // Blue

    // Misc
    //double y_axis_min = .8;
    //double y_axis_min = 0;
    bool color_auto = false;

    // For ratio plots
    TString ratio_plot_title = "(0+1p)/0p";
    double ratio_max_y = 1.5;
    double ratio_min_y = 0.8;
    //TString ratio_plot_title = "ratio to qCut19";
    //double ratio_max_y = 1.03;
    //double ratio_min_y = .97;

    // Plot features
    bool include_legend = true;
    bool include_title = false;
    bool include_ratio = true;
    bool include_error_bands = false;
    bool legend_centered = false;
    bool plotting_NLO_comp_fits = true;

    // X axis lims
    bool use_smefit_lims = true;
    bool use_asimov_lims = false;
    // Y axis lims
    double universal_ymin = 0;
    bool use_hardcode_ylims = true;

    // Skip certain plots or fits
    bool skip_4f = true;
    bool skip_0p = false;

    // SPECIFIC TO PHENO RESULTS SECTION!!!
    if (points_to_plot_with_errorbars.size() == 0){
        include_ratio = true;
        skip_0p = false;
    } else {
        include_ratio = false;
        skip_0p = true;
    }

    //// End setup ////

    if (skip_4f and wc_name.find("i") != std::string::npos){
        return;
    }

    if (wc_fits.size() < 2 and include_ratio == true) {
        std::cout << "\nOnly one fit, not making ratio plot!\n" << std::endl;
        include_ratio = false;
    }

    //plt_ops.x_name = "NP Strength";
    // Name the x axis
    TString x_axis_title = wc_name;
    if (wc_name.back() == 'i') {
        int len = wc_name.size();
        x_axis_title = wc_name.substr(0,len-1);
    }
    plt_ops.x_name = x_axis_title + " Strength";
    //plt_ops.y_name = "\\sigma_{NP}/\\sigma_{SM}";
    plt_ops.y_name = "\\mu = \\sigma/\\sigma_{SM}";

    TString file_type1 = "pdf";
    TString file_type2 = "png";

    TString x_axis_name = plt_ops.x_name;
    TString y_axis_name = plt_ops.y_name;
    TString plot_name = plt_ops.title;

    int fitpts_marker_style  = 7;
    int origpts_marker_style = 29;
    double fitpts_marker_size  = 0.7;
    double origpts_marker_size = 2.0;//2.0;//3.0;
    bool include_fitpts = true;

    plt_ops.setXLimits(0.0,0.0);
    plt_ops.setYLimits(0.0,1.3);
    //plt_ops.setYLimits(-0.01,0.1);

    // Setup the low and high limits for the plot
    std::vector<WCPoint> sorted_pts;
    double x_low,x_high,y_val;
    int pts_sz;
    for (auto& wc_fit: wc_fits) {
        // Adjust plot axis for the rwgt points
        sorted_pts = wc_fit.getFitPoints();
        pts_sz = sorted_pts.size();

        if (wc_fit.getDim() <= 1) {
            // For multi-dim fits the fit points aren't normally along the axis being displayed
            sorted_pts = sortByWeight(sorted_pts);
            plt_ops.updateYLimits(sorted_pts.at(0).wgt,sorted_pts.at(pts_sz - 1).wgt);
        }

        // Adjust plot axis for the ref points
        //for (uint j = 0; j < ref_pts.size(); j++) {
        for (auto& ref_pt: ref_pts) {
            if (!ref_pt.isSMPoint()) {
                // Always include a SM point
                if (ref_pt.getDim() != 1) {
                    // Don't try to plot pts which are in n-Dim WC phase space
                    continue;
                } else if (ref_pt.getStrength(wc_name) == 0.0) {
                    // The point is 1-D, but not for the WC we are plotting
                    continue;
                }
            }
            //WCPoint orig_pt = ref_pts.at(j);
            y_val = wc_fit.evalPoint(&ref_pt);

            plt_ops.updateXLimits(1.2*ref_pt.getStrength(wc_name),1.2*ref_pt.getStrength(wc_name));
            plt_ops.updateYLimits(y_val,y_val);
            plt_ops.updateYLimits(ref_pt.wgt,ref_pt.wgt);
        }
    }

    // Fill dict with x axis ranges from AN
    std::map<string,double> xlims_dict;
    std::vector<string> wc_names_vect{ "ctW", "ctp", "cpQM", "ctei", "ctli", "cQei", "ctZ", "cQlMi", "cQl3i", "ctG", "ctlTi", "cbW", "cpQ3", "cptb", "cpt", "ctlSi"};
    std::vector<double> an_lims_vals_vect{-4.0 , 41.0 , 29.0 , 7.0 , 8.0 , 7.0 , -4.0 , 7.0 , -8.0 , -2.0 , -2.0 , -5.0 , -10.0 , -18.0 , -25.0 , -9.0};
    for (int i = 0; i < wc_names_vect.size(); i++) {
        xlims_dict[wc_names_vect.at(i)] = abs(an_lims_vals_vect.at(i));
    }

    // Hard code some y lims for nice loking plots for pheno paper
    std::map<string,std::pair<double,double>> hard_code_y_lims;
    if (proc_name.find("ttH") != std::string::npos){
        hard_code_y_lims = {
            {"ctG"  , {0.6, 2.20} },
            {"ctW"  , {0.9, 1.45} },
            {"ctZ"  , {0.8, 2.30} },
            {"cpt"  , {0.9, 1.45} },
            {"ctp"  , {0.0, 32.0} }
        };
    } else if (proc_name.find("ttW") != std::string::npos){
        hard_code_y_lims = {
            {"ctG"  , {0.8, 1.35} },
            {"cpQ3" , {0.9, 1.40} },
            {"cpt"  , {0.8, 1.80} },
            {"ctW"  , {0.9, 1.20} },
            {"ctp"  , {0.9, 1.42} },
        };
    } else if (proc_name.find("ttZ") != std::string::npos){
        hard_code_y_lims = {
            {"cpQM" , {0.6, 1.90} },
            {"cpt"  , {0.5, 4.50} },
            {"ctG"  , {0.82, 1.40} },
            {"ctZ"  , {0.7, 3.70} },
        };
    }

    std::map<string,std::pair<double,double>> smefit_lims_dict {
        {"ctG"  , {-0.4, 0.4} },
        {"ctW"  , {-1.8, 0.9} },
        {"cbW"  , {-2.6, 3.1} },
        {"ctZ"  , {-2.1, 4.0} },
        {"cptb" , {-27, 8.7}  },
        {"cpQ3" , {-5.5, 5.8} },
        {"cpQM" , {-3.5, 3}   },
        {"cpt"  , {-13, 18}   },
        {"ctp"  , {-60, 10}   }
    };

    std::map<string,std::pair<double,double>> asimov_lims_dict {
        {"ctW"   , {-1.75, 2.58}   },
        {"ctZ"   , {-3.31, 3.29}   },
        {"ctp"   , {-10.14, 39.02} },
        {"cpQM"  , {-8.04, 26.49}  },
        {"ctG"   , {-1.02, 0.92}   },
        {"cbW"   , {-4.06, 4.03}   },
        {"cpQ3"  , {-8.12, 3.60}   },
        {"cptb"  , {-14.74, 14.57} },
        {"cpt"   , {-22.34, 13.86} },
        {"cQl3i" , {-7.93, 7.28}   },
        {"cQlMi" , {-4.30, 5.39}   },
        {"cQei"  , {-4.76, 4.84}   },
        {"ctli"  , {-4.65, 5.29}   },
        {"ctei"  , {-4.69, 5.27}   },
        {"ctlSi" , {-6.84, 6.84}   },
        {"ctlTi" , {-0.90, 0.90}   }
    };

    // Re-try setting axis ranges w/o using reference pts
    if (plt_ops.x_min == plt_ops.x_max) {
        for (auto& wc_fit: wc_fits) {
            sorted_pts = wc_fit.getFitPoints();
            pts_sz = sorted_pts.size();

            sorted_pts = sortByStrength(sorted_pts,wc_name);
            x_low  = sorted_pts.at(0).getStrength(wc_name);
            x_high = sorted_pts.at(pts_sz - 1).getStrength(wc_name);

            x_low = std::max(x_low,-25.0);
            x_high = std::min(x_high,25.0);

            // Set x axis range to range AN range
            //x_low = -xlims_dict[wc_name];
            //x_high = xlims_dict[wc_name];
            if ( smefit_lims_dict.find(wc_name) != smefit_lims_dict.end() and use_smefit_lims){
                //std::cout << "in smefit_lims_dict ! " << wc_name << std::endl;
                //std::cout << smefit_lims_dict[wc_name].first << " , " << smefit_lims_dict[wc_name].first << std::endl;
                x_low  = smefit_lims_dict[wc_name].first;
                x_high = smefit_lims_dict[wc_name].second;
            } else if (asimov_lims_dict.find(wc_name) != asimov_lims_dict.end() and use_asimov_lims){
                //std::cout << "in asimov_lims_dict! " << wc_name << std::endl;
                //std::cout << asimov_lims_dict[wc_name].first << " , " << asimov_lims_dict[wc_name].first << std::endl;
                x_low  = asimov_lims_dict[wc_name].first;
                x_high = asimov_lims_dict[wc_name].second;
            }

            plt_ops.updateXLimits(x_low,x_high);

            int n_chk_pts = 10;
            double delta = (x_high - x_low) / n_chk_pts;
            for (int i=0; i <= n_chk_pts; i++) {
                y_val = wc_fit.evalPoint(wc_name,x_low+delta*i);
                plt_ops.updateYLimits(y_val,y_val);
            }
        }
    } else {
        std::cout << "Not setting AN lims!!!" << std::endl;
    }

    TCanvas *c1 = new TCanvas("c1","",1000,1000);

    /*
    double pad_left_edge;
    double pad_right_edge;
    //Float_t small = .04;
    //Float_t small = .08;
    //Float_t small = .3;
    //const float padding = 1e-5;
    if (include_ratio) {
        const float ydiv = 0.28;
        //c1->Divide(1,2,small,small);
        c1->Divide(1,2);

        c1->GetPad(1)->SetPad(0,ydiv,1,1); // xlow, ylow, xup, yup
        //c1->GetPad(1)->SetPad(padding,ydiv+padding,1-padding,1-padding);
        //c1->GetPad(1)->SetLeftMargin(.08);
        //c1->GetPad(1)->SetRightMargin(.08);
        c1->GetPad(1)->SetRightMargin(0.06);
        c1->GetPad(1)->SetLeftMargin(0.15);
        //c1->GetPad(1)->SetBottomMargin(.3);
        c1->GetPad(1)->SetBottomMargin(.14);
        c1->GetPad(1)->SetGrid(1,1);
        c1->GetPad(1)->Modified();
        
        //c1->GetPad(2)->SetLeftMargin(.08);
        //c1->GetPad(2)->SetRightMargin(.08);
        c1->GetPad(2)->SetRightMargin(0.06);
        c1->GetPad(2)->SetLeftMargin(0.15);
        //c1->GetPad(2)->SetBottomMargin(.3);
        c1->GetPad(2)->SetBottomMargin(.18);
        //c1->GetPad(2)->SetPad(padding,padding,1-padding,ydiv-padding);
        c1->GetPad(2)->SetPad(0,0,1,ydiv); // xlow, ylow, xup, yup
        c1->GetPad(2)->SetTopMargin(.1);
        c1->GetPad(2)->SetGrid(1,1);
        c1->GetPad(2)->Modified();

        c1->cd(2);
        //gPad->SetTopMargin(small);
        //gPad->SetTopMargin(.1);
        gPad->SetTickx();
        //gPad->Modified();
        c1->cd(1);
        //gPad->SetBottomMargin(small);
        //gPad->SetBottomMargin(.1);
        //gPad->Modified();
        

        pad_left_edge = c1->GetPad(1)->GetLeftMargin();
        pad_right_edge = 1-c1->GetPad(1)->GetRightMargin();
    } else {
        c1->ToggleEventStatus();
        c1->cd();
        c1->SetGrid(1,1);
        c1->GetPad(0)->SetRightMargin(0.06);
        c1->GetPad(0)->SetLeftMargin(0.15);
        c1->GetPad(0)->SetBottomMargin(0.14);
        pad_left_edge = c1->GetPad(0)->GetLeftMargin();
        pad_right_edge = 1-c1->GetPad(0)->GetRightMargin();
    }
    */

    double r = 0.40;
    //double epsilon = 0.03;
    //double epsilon = 0.05;
    double epsilon = 0.07;
    double pad_left_edge;
    double pad_right_edge;
    double pad_top_edge;
    TPad *pad1 = new TPad("pad1", "pad1", 0, r-epsilon, 1, 1);
    TPad *pad2 = new TPad("pad2", "pad2", 0, 0, 1, r*(1-epsilon));
    if (include_ratio){
        pad1->SetLeftMargin(.17);  // Room for the legend
        pad1->SetBottomMargin(0.0);
        pad1->SetGrid(1,1);        // Grid
        pad1->Draw();              // Draw the upper pad: pad1
        c1->cd();                  // Go back to the main canvas before defining pad2
        pad2->SetLeftMargin(.17);  // Room for the legend
        pad2->SetTopMargin(0.0);
        pad2->SetBottomMargin(0.45);
        pad2->SetGrid(1,1);
        pad2->Draw();
        pad1->cd();
        pad_left_edge = pad1->GetLeftMargin();
        pad_right_edge = 1-pad1->GetRightMargin();
        pad_top_edge = 1-pad1->GetTopMargin();
    } else {
        c1->SetBottomMargin(.20);
        c1->SetLeftMargin(.20);
        c1->SetGrid(1,1);
        pad_left_edge = c1->GetLeftMargin();
        pad_right_edge = 1-c1->GetRightMargin();
        pad_top_edge = 1-c1->GetTopMargin();
    }

    
    double left,right,top,bottom,scale_factor,minimum;
    // The legend
    double ledgend_width = 0.65;
    TLegend *legend;
    //if (include_ratio or legend_centered) {
    if (legend_centered) {
        left  = (pad_right_edge-pad_left_edge)/2 + pad_left_edge - ledgend_width/2; // Center the legend on the pad, not canvas
        right = (pad_right_edge-pad_left_edge)/2 + pad_left_edge + ledgend_width/2; // Center the legend on the pad, not canvas
        top    = 0.88;
        bottom = 0.7;
        legend = new TLegend(left,top,right,bottom);
        legend->SetNColumns(wc_fits.size());
        legend->SetBorderSize(0);
    } else {
        /*
        left         = 0.81;
        right        = 0.98;
        top          = 0.9;
        scale_factor = 0.05;
        minimum      = 0.1;
        bottom = std::max(top - scale_factor*(wc_fits.size()+1),minimum);
        legend = new TLegend(left,top,right,bottom);
        */
        // Top right
        left   = pad_right_edge-ledgend_width;
        right  = pad_right_edge-0.003;
        top    = pad_top_edge-0.003;
        bottom = 0.71;
        legend = new TLegend(left,top,right,bottom);
        legend->SetNColumns(wc_fits.size());
        legend->SetBorderSize(0);        
    }

    bool include_orig_pts = (wc_fits.size() == ref_pts.size());

    std::vector<TH1D*> hist_vect;
    std::vector<int> fit_color_vect;
    std::vector<int> fit_style_vect;
    for (uint i = 0; i < wc_fits.size(); i++) {

        WCFit wc_fit = wc_fits.at(i);
        double s0 = wc_fit.getCoefficient(wc_fit.kSMstr,wc_fit.kSMstr);
        double s1 = wc_fit.getCoefficient(wc_fit.kSMstr,wc_name);
        double s2 = wc_fit.getCoefficient(wc_name,wc_name);

        // TMP, usefule for making pheno paper plots, make 4.1 and 4.2 plots at same time...
        if (skip_0p and (string(wc_fit.getTag()).find("0p") != std::string::npos) ){
            if (i==0){
                std::cout << "\nCannot skip fit " << wc_fit.getTag() << " as it is the first and would mess up the drawing...\n" << std::endl;
            }
            std::cout << "\nSkipping fit " << wc_fit.getTag() << " ...\n" << std::endl;
            continue;
        }

        TF1* fit = new TF1("fit","pol2",plt_ops.x_min,plt_ops.x_max);
        fit->SetParameter(0,s0);
        fit->SetParameter(1,s1);
        fit->SetParameter(2,s2);
        // Set axes name
        fit->GetXaxis()->SetTitle(x_axis_name);
        fit->GetYaxis()->SetTitle(y_axis_name);
        /*
        // Set size of axis lables and offset from axis
        fit->GetXaxis()->SetTitleSize(0.060);
        fit->GetYaxis()->SetTitleSize(0.075);
        fit->GetYaxis()->SetTitleOffset(0.98);
        fit->GetXaxis()->SetTitleOffset(0.98);
        // Set the sizes of the numbers on the axes:
        fit->GetYaxis()->SetLabelSize(0.045);
        fit->GetXaxis()->SetLabelSize(0.045);
        fit->GetYaxis()->SetLabelOffset(0.015);
        fit->GetXaxis()->SetLabelOffset(0.015);
        */
        fit->GetXaxis()->SetTitleSize(0.060);
        fit->GetYaxis()->SetTitleSize(0.12);
        fit->GetYaxis()->SetTitleOffset(0.68);
        fit->GetXaxis()->SetTitleOffset(0.9);
        fit->GetYaxis()->SetLabelSize(0.065);
        fit->GetXaxis()->SetLabelSize(0.065);
        fit->GetYaxis()->SetLabelOffset(0.015);
        fit->GetXaxis()->SetLabelOffset(0.015);
        if (not include_ratio){
            /*
            fit->GetXaxis()->SetTitleOffset(1.0);
            fit->GetYaxis()->SetTitleOffset(0.8);
            fit->GetXaxis()->SetLabelSize(0.061);
            fit->GetYaxis()->SetLabelSize(0.061);
            fit->GetXaxis()->SetTitleSize(0.065);
            fit->GetYaxis()->SetTitleSize(0.092);
            */
            fit->GetXaxis()->SetLabelSize(0.050);
            fit->GetYaxis()->SetLabelSize(0.050);
            fit->GetXaxis()->SetTitleOffset(1.1);
            fit->GetYaxis()->SetTitleOffset(1.0);
            fit->GetYaxis()->SetTitleSize(0.10);
            fit->GetXaxis()->SetNdivisions(7,5,0,true);
        }
        if (include_title){
            fit->SetTitle(plot_name);
        } else {
            fit->SetTitle(""); // no title
        }


        // Error bands: 
        const int npts = 100;
        Float_t x_vals[npts];
        Float_t y_vals[npts];
        Float_t xerr_vals[npts];
        Float_t yerr_vals[npts];

        float x_min = plt_ops.x_min;
        float x_max = plt_ops.x_max;
        float x_range = x_max - x_min;
        float step_size = x_range/(npts-1);

        TH1D* hist = new TH1D("hist","",npts,x_min,x_max); // Hist for getting ratio plot

        WCPoint wc_pt;
        float x_coord;
        for (int idx = 0; idx < npts; idx++){
            //x_coord = x_min + idx*step_size;
            x_coord = hist->GetBinCenter(idx+1);
            wc_pt.setStrength(wc_name,x_coord);
            x_vals[idx] = x_coord;
            y_vals[idx] = wc_fit.evalPoint(&wc_pt);
            xerr_vals[idx] = 0;
            yerr_vals[idx] = wc_fit.evalPointError(&wc_pt);

            hist->SetBinContent(idx+1,y_vals[idx]);
        }

        hist_vect.push_back(hist);
        TGraphErrors *err_graph = new TGraphErrors(npts, x_vals, y_vals, xerr_vals, yerr_vals);

        //fit->SetLineStyle(9);
        fit->SetLineWidth(4);

        std::string tmp_str = wc_fit.getTag();
        TString leg_str = tmp_str;

        // This is where figure out what color to draw the fit as. It's a mess.
        //std::cout << "\nThe leg str is: " << leg_str << "\n" << std::endl;
        if (string(leg_str).find("Jet") != std::string::npos or string(leg_str).find("NLO") != std::string::npos or string(tmp_str).find("0+1p") != std::string::npos){
           fit->SetLineStyle(7); // Plot the 0+1p lines as dashed 
        }
        if (string(leg_str).find("0p") != std::string::npos or string(tmp_str).find("0+1p") != std::string::npos){ // Private samples usually indicated this way
            fit->SetLineColor(private_color);
        //} else if (string(leg_str).find("ttlnu") != std::string::npos or string(leg_str).find("ttll") != std::string::npos or string(leg_str).find("ttH") != std::string::npos){ // If it's ttll or ttlnu, probably private
        } else if (string(leg_str).find("ttlnu") != std::string::npos or string(leg_str).find("ttll") != std::string::npos ){ // If it's ttll or ttlnu, probably private
            fit->SetLineColor(private_color);
        } else if (string(leg_str).find("smeft") != std::string::npos){ // Reza's samples
            fit->SetLineColor(smeft_color);
        } else if (string(leg_str).find("ttW") != std::string::npos or string(leg_str).find("ttZ") != std::string::npos or string(leg_str).find("ttH") != std::string::npos){ // New ttW, ttZ private samples for comp with Reza's NLO
            fit->SetLineColor(ttWttZ_color);
        }
        if (string(leg_str).find("qed=1") != std::string::npos){ // If it's a QED order constrained sample, override any previous color and color specifically for this
            fit->SetLineColor(QEDorder_comp_color);
        }
        if (string(leg_str).find("qed=2") != std::string::npos){ // If it's a QED order constrained sample, override any previous color and color specifically for this
            fit->SetLineColor(QEDorder_2_comp_color);
        }
        if (string(leg_str).find("dim6TopvMay2020") != std::string::npos){ // The newly updated dim6TOp modle
            fit->SetLineColor(dim6TopvMay2020_color);
        }

        // Auto color the lines
        if (color_auto){
            fit->SetLineColor(i+1);
        }

        if (i == 0) {
            fit->Draw();
            //hist->Draw("SAME");
            legend->AddEntry(fit,leg_str,"l");
        } else {
            fit->Draw("LSAME");
            //hist->Draw("SAME");
            legend->AddEntry(fit,leg_str,"l");
        }

        // Set y limits for plots
        double tmp_ymax = max( fit->Eval(fit->GetXmin()) , fit->Eval(fit->GetXmax()) );
        if (use_hardcode_ylims and hard_code_y_lims.find(wc_name) !=hard_code_y_lims.end()){ // We have a hard coded limit
            fit->GetYaxis()->SetRangeUser(hard_code_y_lims[wc_name].first,hard_code_y_lims[wc_name].second);
            std::cout << "\nSetting hard coded y lims to: " << hard_code_y_lims[wc_name].first << " , " << hard_code_y_lims[wc_name].second << std::endl;
        } else {
            fit->SetMaximum(1.5*tmp_ymax);
            fit->SetMinimum(universal_ymin);
        }

        if (include_error_bands){
            //err_graph->SetFillColor(plt_ops.getColor(i));
            err_graph->SetFillColor(fit->GetLineColor());
            //err_graph->SetFillStyle(3003);
            err_graph->SetFillStyle(3002);
            if (string(leg_str).find("cpQM at -cpQ3") == std::string::npos and string(leg_str).find("smeft") == std::string::npos){ // Don't plot the error band for the weird cpQ3 - cpQM plot or the smeft plots from Reza
                err_graph->Draw("3");
                //err_graph->Draw("SAME");
            }
        }
        fit_color_vect.push_back(fit->GetLineColor()); // Keep track of this so can match ratio plots to correct fits
        fit_style_vect.push_back(fit->GetLineStyle()); // Keep track of this so can match ratio plots to correct fits



        // Account for points y values, if not already using hard coded y lims
        int npts_to_plot = points_to_plot_with_errorbars.size();
        double tmp_max_nlo_pt = 0;
        for (int npt=0; npt<npts_to_plot; npt++){
            customPointInfo points_struct = points_to_plot_with_errorbars[npt];
            if (points_struct.y > tmp_max_nlo_pt){
                tmp_max_nlo_pt = points_struct.y;
            }
        }
        if (not use_hardcode_ylims and tmp_max_nlo_pt*1.5 > fit->GetMaximum()){
            fit->SetMaximum(tmp_max_nlo_pt*1.5);
        }

    }


    if (plotting_NLO_comp_fits and quad_fit_vect.size()!=0){
        std::cout << "Number of extra fits to plot: " << quad_fit_vect.size() << std::endl;
        double quad_fit_s0, quad_fit_s1, quad_fit_s2;
        for (int i=0; i<quad_fit_vect.size(); i++){
            if (quad_fit_vect[i].second.size() != 0){
                //std::cout << quad_fit_vect[i].first << std::endl;
                //std::cout << quad_fit_vect[i].second[0] << " " << quad_fit_vect[i].second[1] << " " << quad_fit_vect[i].second[2] << std::endl;
                TF1* comp_fit = new TF1("fit","pol2",plt_ops.x_min,plt_ops.x_max);
                comp_fit->SetParameter(0,quad_fit_vect[i].second[0]);
                comp_fit->SetParameter(1,quad_fit_vect[i].second[1]);
                comp_fit->SetParameter(2,quad_fit_vect[i].second[2]);
                //comp_fit->SetLineColor(plt_ops.getColor(i+wc_fits.size()));
                //comp_fit->SetLineStyle(7);
                TString leg_str_for_comp_fit = quad_fit_vect[i].first;
                legend->AddEntry(comp_fit,leg_str_for_comp_fit,"l");
                comp_fit->SetLineColor(arxiv_color);
                if (string(leg_str_for_comp_fit).find("NLO") != std::string::npos){
                    comp_fit->SetLineStyle(7); // NLO as dashed
                }
                /*
                // Specific to plotting NLO and LO comps (depends on knowing if 0p or 0+1p is plotted first):
                if (quad_fit_vect[i].first.find("NLO") != std::string::npos){
                   comp_fit->SetLineColor(plt_ops.getColor(0));                    
                } else {
                    comp_fit->SetLineColor(plt_ops.getColor(1));
                }
                */
                comp_fit->Draw("same");
            }
        }
    }
    // For the custom points with errorbars (points, not fits) (Reza smeftnlo in this case)
    if (plotting_NLO_comp_fits and points_to_plot_with_errorbars.size() !=0){
        int npts_to_plot = points_to_plot_with_errorbars.size();
        std::vector<double> x_vec;
        std::vector<double> y_vec;
        std::vector<double> x_e_vec;
        std::vector<double> y_e_vec;
        Float_t x_point_vals[npts_to_plot];
        Float_t y_point_vals[npts_to_plot];
        Float_t x_err_vals[npts_to_plot];
        Float_t y_err_vals[npts_to_plot];
        TString points_with_ebars_leg_str;
        for (int npt=0; npt<npts_to_plot; npt++){
            customPointInfo points_struct = points_to_plot_with_errorbars[npt];
            print_customPointInfo_object({points_struct});
            x_point_vals[npt]=points_struct.x;
            y_point_vals[npt]=points_struct.y;
            x_err_vals[npt]=0;
            y_err_vals[npt]=points_struct.err;
            points_with_ebars_leg_str = points_struct.proc_name + " " + points_struct.type; // Should be the same for each point in stuct?
            //if (y_point_vals[npt]*1.5 > fit->GetMaximum()){
                //fit->SetMaximum(1.5*y_point_vals[npt]);
            //}
        }
        TGraphErrors *points_with_ebars = new TGraphErrors(points_to_plot_with_errorbars.size(),x_point_vals,y_point_vals,x_err_vals,y_err_vals);
        points_with_ebars->SetMarkerStyle(4); // Open circles
        points_with_ebars->SetMarkerSize(2);
        points_with_ebars->SetMarkerColor(smeft_color);
        points_with_ebars->SetLineColor(smeft_color);
        points_with_ebars->SetLineWidth(5);
        points_with_ebars->Draw("same P");
        legend->AddEntry(points_with_ebars,points_with_ebars_leg_str,"p"); // Legend for points
        //legend->AddEntry(points_with_ebars,points_with_ebars_leg_str,"l"); // Legend for points

    }


    // Calculate the ratio hists
    if (include_ratio){
        //c1->cd(2);
        pad2->cd(); // pad2 becomes the current pad
        TH1D* ratio_hist;
        //std::cout << hist_vect.size() << std::endl;
        for (int i=0; i<hist_vect.size(); i++){
            ratio_hist = makeRatioHistogram("rhist",hist_vect.at(i),hist_vect.at(1)); // Need to choose which hist to divide w.r.t.
            //ratio_hist = makeRatioHistogram("rhist",hist_vect.at(i),hist_vect.at(0)); // Need to choose which hist to divide w.r.t.
            //ratio_hist->GetYaxis()->SetNdivisions(010,true);
            //ratio_hist->GetYaxis()->SetNdivisions(4,1,0,true);
            ratio_hist->GetYaxis()->SetNdivisions(4,5,0,true);
            ratio_hist->GetYaxis()->SetTitle(ratio_plot_title);
            /*
            ratio_hist->GetYaxis()->SetTitleOffset(0.4);
            //ratio_hist->GetYaxis()->SetTitleSize(0.12);
            ratio_hist->GetYaxis()->SetTitleSize(0.19);
            ratio_hist->GetYaxis()->SetLabelSize(0.15);
            ratio_hist->GetXaxis()->SetLabelSize(0.15);
            ratio_hist->GetYaxis()->SetLabelOffset(0.015);
            ratio_hist->GetXaxis()->SetLabelOffset(0.015);
            */

            ratio_hist->GetYaxis()->SetTitleOffset(0.50);
            ratio_hist->GetXaxis()->SetTitleSize(0.16);
            //ratio_hist->GetYaxis()->SetTitleSize(0.18);
            ratio_hist->GetYaxis()->SetTitleSize(0.16);
            ratio_hist->GetXaxis()->SetLabelSize(0.15);
            ratio_hist->GetYaxis()->SetLabelSize(0.10);
            ratio_hist->GetXaxis()->SetLabelOffset(0.015);
            ratio_hist->GetYaxis()->SetLabelOffset(0.015);
            ratio_hist->GetXaxis()->SetTickLength(0.09);
            ratio_hist->GetXaxis()->SetTitle(x_axis_name);

            ratio_hist->SetLineWidth(3);
            ratio_hist->SetLineColor(fit_color_vect.at(i));
            ratio_hist->SetLineStyle(fit_style_vect.at(i));
            //ratio_hist->SetLineColor(plt_ops.getColor(i));
            ratio_hist->SetMaximum(ratio_max_y);
            ratio_hist->SetMinimum(ratio_min_y);
            ratio_hist->Draw("SAME C");
        }
    }
    

    // Draw MadGraph reference points (if they happen to land on this particular 1-D axis scan)
    for (uint i = 0; i < ref_pts.size(); i++) {
        WCPoint ref_pt = ref_pts.at(i);
        if (!ref_pt.hasWC(wc_name)) {
            // Only include points relevant to this 1D scan
            continue;
        }

        //std::cout << "iteration: " << i << " wc name: " << wc_name << " val: " << ref_pt.getStrength(wc_name) << std::endl;
        if (!ref_pt.isSMPoint()) {
            // Need to explicitly include the SM point
            if (ref_pt.getDim() != 1) {
                // The ref point is n-Dim WC phase space (we probably messed up in selecting our reference gridpacks)
                continue;
            } else if (ref_pt.getStrength(wc_name) == 0.0) {
                // The ref point has no dependence on this particular 1D plot
                continue;
            }
        }

        TGraph* ref_pt_gr = new TGraph(1);
        ref_pt_gr->SetPoint(0,ref_pt.getStrength(wc_name),ref_pt.wgt);
        ref_pt_gr->SetMarkerStyle(origpts_marker_style);
        ref_pt_gr->SetMarkerSize(origpts_marker_size*1.5);
        ref_pt_gr->SetMarkerColor(1);  // Paint it black!
        ref_pt_gr->Draw("P");
    }

    legend->SetFillColor(0);
    legend->SetNColumns(2); // For two columns

    if (include_legend) {
        if (include_ratio){
            //c1->cd(1);
            pad1->cd();
            legend->Draw();
        } else {
            legend->Draw();
        }
    }

    double txt_x = 0.9;
    double txt_y = 0.93;
    double txt_sz = 0.02;
    double txt_xoff = 0.0;
    double txt_yoff = -0.02;

    TLatex CMS_text = TLatex(txt_x,txt_y,"CMS Preliminary Simulation");
    TLatex Lumi_text = TLatex(txt_x + txt_xoff,txt_y + txt_yoff,"Luminosity = 41.29 fb^{-1}");
    
    CMS_text.SetNDC(1);
    CMS_text.SetTextSize(txt_sz);
    CMS_text.SetTextAlign(30);

    Lumi_text.SetNDC(1);
    Lumi_text.SetTextSize(txt_sz);
    Lumi_text.SetTextAlign(30);

    //CMS_text.Draw("same");
    //Lumi_text.Draw("same");

    c1->Update();


    TString save_name;
    
    save_name = plt_ops.output_dir + "/" + plt_ops.tag + "." + file_type1;
    c1->Print(save_name,file_type1);

    save_name = plt_ops.output_dir + "/" + plt_ops.tag + "." + file_type2;
    c1->Print(save_name,file_type2);

    delete legend;
    delete c1;
}


void make_binned_weight_plot(
    TH1D *hist,
    std::string output_dir,
    std::string fname
)
{
    TCanvas *c1 = new TCanvas("c1","",1280,720);
    c1->ToggleEventStatus();
    c1->cd();
    c1->SetLogy();
    c1->SetLogx();
    c1->SetGrid(1,1);

    hist->Draw();

    TString save_name = output_dir + "/" + fname + ".pdf";
    c1->Print(save_name,"pdf");

    delete c1;
}

// Saves the specified list of fits to a txt file
void make_fitparams_file(
    std::string fpath,
    std::vector<WCFit> wc_fits
)
{
    for (uint i = 0; i < wc_fits.size(); i++) {
        WCFit fit = wc_fits.at(i);
        if (i == 0) {
            //fit.dump(false);
            fit.save(fpath,false);
        } else {
            //fit.dump(true);
            fit.save(fpath,true);
        }
    }
}

/*
// For 1-D fits, generate a fit using only MadGraph starting points
void make_dedicated_fits(
    std::string process,
    std::vector<std::string> wc_names,
    std::vector<WCPoint> orig_pts
)
{
    std::string output_dir = "read_lhe_outputs";
    for (auto& wc_name: wc_names) {
        std::vector<WCPoint> fit_pts;
        for (uint i = 0; i < orig_pts.size(); i++) {
            WCPoint wc_pt = orig_pts.at(i);
            if (wc_pt.hasWC(wc_name) && wc_pt.isSMPoint()) {
                // Always add SM points to the list
                fit_pts.push_back(wc_pt);
            } else if (wc_pt.hasWC(wc_name) && wc_pt.getDim() == 1 && wc_pt.getStrength(wc_name) != 0.0) {
                // Add 1-D points which are non-zero for the WC of interest
                fit_pts.push_back(wc_pt);
            }
        }

        if (fit_pts.size() < 3) {
            std::cout << "[ERROR] Not enough fit points for " << wc_name << ", skipping..." << std::endl;
            continue;
        }

        std::string fit_tag = process + "_" + wc_name + "_" + "orig";
        WCFit fit(fit_pts,fit_tag);
        std::string save_path = output_dir + "/" + "fitparams_" + process + "_" + wc_name + ".txt";
        fit.save(save_path);
        for (uint i = 0; i < fit_pts.size(); i++) {
            WCPoint wc_pt = fit_pts.at(i);
            wc_pt.dump(wc_name,false);
        }
        std::cout << std::endl;
    }
}
*/



#endif
/* MAKEEFTPLOTS */

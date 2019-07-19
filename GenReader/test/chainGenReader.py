import sys
import os
import subprocess
import argparse

from EFTGenReader.GenReader.DatasetHelper import DatasetHelper

# Wrapper script to run over multiple samples back-to-back

arg_parser = argparse.ArgumentParser(prog='chainGenReader.py',description='Wrapper script to run multiple samples to EFTGenReader_cfg in series')
arg_parser.add_argument('-n','--nevents',metavar='N',default=10,type=int,help='number of events to run over')
arg_parser.add_argument('--test',action='store_true',help='will pass the test option on through to the underlying cmsRun call')
arg_parser.add_argument('--list',action='store_true',help='Print out all known sample names from the JSON database and exit')
arg_parser.add_argument('--out-suffix',metavar='NAME',default='_NoTopLeptons_output_tree',help='add a different suffix for the output root files')
arg_parser.add_argument('--minptj',metavar='VAL',type=float,help='min pt cut for genjets')
arg_parser.add_argument('--minptl',metavar='VAL',type=float,help='min pt cut for genleptons')
arg_parser.add_argument('--maxetaj',metavar='VAL',type=float,help='max eta cut for genjets')
arg_parser.add_argument('--maxetal',metavar='VAL',type=float,help='max eta cut for genleptons')
arg_parser.add_argument('datasets',metavar='NAME',nargs='*',help='specify multiple datasets to run over, one after the other')
args = arg_parser.parse_args()

GIT_REPO_DIR = subprocess.check_output(['git','rev-parse','--show-toplevel']).strip()

# All samples
# central_ttZ central_ttW central_tZq central_ttH ttll_FP_R4B9 ttll_SM ttllNoHiggs_SM ttllNoHiggs_EFT ttllnunuNoHiggs_SM tllq_FR_R4B9 tllq_SM tllq4f_SMNoSchanW tllq4fMatched_SM tllq4fMatched_EFT tllq4fNoHiggs_SM tllq4fNoHiggs_EFT ttlnu_FP_R4B9 ttlnu_SM ttlnu_EFT ttlnuJet_EFT ttlnu_NoPDFWeights ttH_SM

# Central samples
# central_ttZ central_ttW central_tZq central_ttH

# private ttll samples
# ttll_FP_R4B9 ttll_SM ttllNoHiggs_SM ttllNoHiggs_EFT ttllnunuNoHiggs_SM

# private tllq samples
# tllq_FR_R4B9 tllq_SM tllq4f_SMNoSchanW tllq4fMatched_SM tllq4fMatched_EFT tllq4fNoHiggs_SM tllq4fNoHiggs_EFT

# private ttlnu samples
# ttlnu_FP_R4B9 ttlnu_SM ttlnu_EFT ttlnuJet_EFT ttlnu_NoPDFWeights

# private ttH samples
# ttH_SM

def main():
    if args.list:
        ds_helper = DatasetHelper()
        json_fpath = os.path.join(GIT_REPO_DIR,'GenReader/data/JSON/datasets.json')
        ds_helper.load(json_fpath)
        print "---Available Samples---"
        for sample_name in sorted(ds_helper.list()):
            print "{name}".format(name=sample_name)
        return

    if len(arg.datasets) == 0:
        print "No samples specified, exiting..."
        return

    lst = args.datasets
    if args.test and len(lst) > 1:
        print "WARNING: Only running over the first sample when in 'test' mode"
        lst = [lst[0]]
    max_events = args.nevents
    norm_type = 1
    intg_lumi = 1.0
    base_cmd = ["cmsRun","EFTGenReader_cfg.py"]
    width = 1
    if len(lst) >= 10:
        width = 2
    for idx,ds_name in enumerate(lst):
        print "\nProcessing dataset {0}... [{1:0>{w}}/{2}]".format(ds_name,idx+1,len(lst),w=width)
        full_cmd = [x for x in base_cmd]
        full_cmd.append("maxEvents={}".format(max_events))
        full_cmd.append("normType={}".format(norm_type))
        full_cmd.append("intgLumi={}".format(intg_lumi))
        full_cmd.append("fnSuffix={}".format(args.out_suffix))
        if args.test:
            full_cmd.append("test=True")

        # Add kinematic cut options
        if args.minptj:
            full_cmd.append("minPtJet={}".format(args.minptj))
        if args.minptl:
            full_cmd.append("minPtLep={}".format(args.minptl))
        if args.maxetaj:
            full_cmd.append("maxEtaJet={}".format(args.maxetaj))
        if args.maxetal:
            full_cmd.append("maxEtaLep={}".format(args.maxetal))

        full_cmd.append("dataset={}".format(ds_name))
        print "Full Command: {}".format(" ".join(full_cmd))
        subprocess.check_call(full_cmd)

if __name__ == "__main__":
    main()

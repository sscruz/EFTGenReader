import sys
import os
import subprocess
import argparse
import datetime
from EFTGenReader.GenReader.utils import regex_match,run_process,get_files,move_files,clean_dir
from EFTGenReader.GenReader.make_html import make_html

USER_DIR = os.path.expanduser('~')

tstamp = datetime.datetime.now().strftime('%Y%m%d-%H%M')

arg_parser = argparse.ArgumentParser(prog='wrapperEFTGenPlots.py')
arg_parser.add_argument('--web-dir',metavar='DIR',default=os.path.join(USER_DIR,'www'),help='location to save all sub-directories to; %(metavar)s can either be an absolute or relative path')
arg_parser.add_argument('-d','--dir',metavar='DIR',default='testing01',help='name of sub-directory to save plots in')
arg_parser.add_argument('-f','--force',action='store_true',help='remove all .png and .html files from output directory before moving plots over')
arg_parser.add_argument('--wcStr',action='append',metavar='NAME:VAL',help='specify WC value that the TH1EFT should reweight to. Syntax should be %(metavar)s')
arg_parser.add_argument('--from-file',metavar='FILE',help='adds files or directories to be processed read from %(metavar)s, one entry per line')
arg_parser.add_argument('--no-timestamp',action='store_true',help='when merging, do not include timestamp in the name of the merged file')
arg_parser.add_argument('infiles',metavar='FILE',nargs='*',help='%(metavar)s should be a path to root file (or directory with files) produced by EFTGenReader_cfg')
args = arg_parser.parse_args()

GIT_REPO_DIR = subprocess.check_output(['git','rev-parse','--show-toplevel']).strip()
TEST_DIR = os.path.join(GIT_REPO_DIR,'GenReader/test')

# Returns the mode we are in based on the first input
def get_input_mode(inputs):
    first = inputs[0]
    mode = None
    if os.path.isdir(first):
        mode = 'dir_mode'
    elif os.path.isfile(first):
        mode = 'file_mode'
    else:
        print "ERROR: Unknown input {path}".format(path=first)
        raise RuntimeError
    return mode

# Check that the inputs are either all files or all directories
def check_inputs(inputs):
    mode = get_input_mode(inputs)
    for idx,fpath in enumerate(inputs):
        if not os.path.exists(fpath):
            print "ERROR: Unknown fpath {path}".format(path=fpath)
            raise RuntimeError
        elif mode == 'dir_mode' and os.path.isfile(fpath):
            print "ERROR: Mixed file and directory inputs!"
            raise RuntimeError
        elif mode == 'file_mode' and os.path.isdir(fpath):
            print "ERROR: Mixed file and directory inputs!"
            raise RuntimeError
    return mode

def main():
    if not os.path.exists(args.web_dir):
        print "Unknown directory: %s" % (args.web_dir)
        return

    cur_dir = os.getcwd()

    wc_string = "wcpoint"
    if args.wcStr:
        for s in args.wcStr:
            wc,v = s.split(':')
            wc_string += "_{wc}_{val}".format(wc=wc,val=v)
    wc_string = '"{str}"'.format(str=wc_string)
    print "Reweight String: {str}".format(str=wc_string)

    inputs = []
    if args.from_file:
        if not os.path.exists(args.from_file):
            print "ERROR: Unknown file: {fpath}".format(fpath=args.from_file)
            raise RuntimeError
        with open(args.from_file) as f:
            for l in f:
                l = l.strip().split('#')[0]
                l = l.split('#')[0]
                if l.strip():
                    inputs.append(l.strip())
    inputs.extend(args.infiles)

    if len(inputs) == 0:
        print "No inputs specified, exiting..."
        return

    mode = check_inputs(inputs)

    file_lst = []
    for idx,fpath in enumerate(inputs):
        if mode == 'file_mode':
            file_lst.append(fpath)
        elif mode == 'dir_mode':
            tmp_str = os.path.normpath(fpath)
            tmp_str = os.path.split(tmp_str)[1]
            merge_name = "merged_{tstamp}_{name}.root".format(tstamp=tstamp,name=tmp_str)
            if args.no_timestamp:
                merge_name = "merged_{name}.root".format(name=tmp_str)
            merged_output = os.path.join(TEST_DIR,'output',merge_name)
            merged_output = os.path.relpath(merged_output,cur_dir)
            if os.path.exists(merged_output):
                print "Skipping {name}, since it already exists in {dir}!".format(name=merge_name,dir='output')
                continue
            to_merge = []
            for fn in os.listdir(fpath):
                if not ".root" in fn:
                    continue
                to_merge.append(os.path.join(fpath,fn))
            if len(to_merge):
                hadd_cmd = ['hadd',merged_output] + to_merge
                print "Merge command: {s1}\n\t{s2}".format(s1=' '.join(hadd_cmd[:2]),s2='\n\t'.join(hadd_cmd[2:]))
                run_process(hadd_cmd)
                file_lst.append(merged_output)

    s = ",".join('"{fn}"'.format(fn=fn) for fn in file_lst)
    subprocess.check_call(["root", "-b", "-l", "-q","makeEFTGenPlots.C(%s,%s)" % ("{{{}}}".format(s),wc_string)])

    print "Output Dir: {}".format(args.dir)

    output_dir = os.path.join(args.web_dir,args.dir)
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    if args.force:
        clean_dir(output_dir,["^.*\.png$","^index.html$"])

    imgs = get_files('.',targets=["^h_.*\.png$"])
    move_files(files=imgs,target=output_dir)
    imgs = get_files('.',targets=["^h_.*\.pdf$"])
    move_files(files=imgs,target=output_dir)
    make_html(output_dir)

if __name__ == "__main__":
    main()

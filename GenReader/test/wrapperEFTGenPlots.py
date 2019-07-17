import sys
import os
import subprocess
import argparse
from EFTGenReader.GenReader.utils import regex_match
from EFTGenReader.GenReader.make_html import make_html

USER_DIR = os.path.expanduser('~')

arg_parser = argparse.ArgumentParser(prog='wrapperEFTGenPlots.py')
arg_parser.add_argument('--web-dir',metavar='WEBDIR',default=os.path.join(USER_DIR,'www'),help='location to save all sub-directories to; %(metavar)s can either be an absolute or relative path')
arg_parser.add_argument('-d','--dir',metavar='DIR',default='testing01',help='name of sub-directory to save plots in')
arg_parser.add_argument('-f','--force',action='store_true',help='Remove all .png and .html files from output directory before moving plots over')
arg_parser.add_argument('--wcStr',action='append',metavar='NAME:VAL',help='Specify WC value that the TH1EFT should reweight to. Syntax should be %(metavar)s')
arg_parser.add_argument('infiles',metavar='FILE',nargs='+',help='%(metavar)s should be a path to root file produced by EFTGenReader_cfg')
args = arg_parser.parse_args()

#WEB_AREA_DIR = "/afs/crc.nd.edu/user/a/awightma/www"
#GEN_PLOTS_DIR = os.path.join(WEB_AREA_DIR,"eft_stuff/misc/gen_plots")

def get_files(tdir):
    if not os.path.exists(tdir): return []
    return [x for x in os.listdir(tdir) if os.path.isfile(x)]

def move_files(files,target):
    width = len(max(files,key=len))
    for src in files:
        dst = os.path.join(target,src)
        #print "[{0:0>{w1}}] {1:<{w2}} --> {2}".format(len(src),src,dst,w1=2,w2=width)
        os.rename(src,dst)

# Removes files from tdir which match any of the regex in targets list
def clean_dir(tdir,targets):
    fnames = regex_match(get_files(tdir),targets)
    if len(fnames) == 0: return
    print "Removing files from: {}".format(tdir)
    print "\tTargets: {}".format(targets)
    for fn in fnames:
        fpath = os.path.join(tdir,fn)
        print "\tRemoving {}".format(fn)
        #os.remove(fpath)


def main():
    if not os.path.exists(args.web_dir):
        print "Unknown directory: %s" % (args.web_dir)
        return

    cur_dir = os.getcwd()

    wc_string = "wcpoint"
    for s in args.wcStr:
        wc,v = s.split(':')
        wc_string += "_{wc}_{val}".format(wc=wc,val=v)
    wc_string = '"{str}"'.format(str=wc_string)

    if len(args.wcStr):
        print "Reweight String: {str}".format(str=wc_string)

    s = ",".join('"{fn}"'.format(fn=fn) for fn in args.infiles)
    subprocess.check_call(["root", "-b", "-l", "-q","makeEFTGenPlots.C(%s,%s)" % ("{{{}}}".format(s),wc_string)])

    print "Output Dir: {}".format(args.dir)

    output_dir = os.path.join(args.web_dir,args.dir)
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    if args.force:
        clean_dir(output_dir,["^.*\.png$","^index.html$"])

    imgs = regex_match(get_files('.'),["^h_.*\.png$"])
    move_files(files=imgs,target=output_dir)
    make_html(output_dir)

if __name__ == "__main__":
    main()
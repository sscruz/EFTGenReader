import sys
import os
import subprocess
from EFTGenReader.Reader.utils import regex_match
from EFTGenReader.Reader.make_html import make_html

WEB_AREA_DIR = "/afs/crc.nd.edu/user/a/awightma/www"
GEN_PLOTS_DIR = os.path.join(WEB_AREA_DIR,"eft_stuff/misc/gen_plots")
ROOT_OUTPUT_DIR = "output"

def get_files(tdir):
    if not os.path.exists(tdir): return []
    return [x for x in os.listdir(tdir) if os.path.isfile(x)]

def move_files(files,target):
    width = len(max(files,key=len))
    for src in files:
        dst = os.path.join(target,src)
        #print "[{0:0>{w1}}] {1:<{w2}} --> {2}".format(len(src),src,dst,w1=2,w2=width)
        os.rename(src,dst)

def main():
    if len(sys.argv) < 3:
        print "Not enough arguments!",sys.argv
        return
    if not os.path.exists(GEN_PLOTS_DIR):
        print "Unknown directory: %s" % (GEN_PLOTS_DIR)
        return

    cur_dir = os.getcwd()

    #outf = sys.argv[1]
    out_dir = sys.argv[1]
    infiles = sys.argv[2:]

    s = ",".join('"{}"'.format(os.path.join(ROOT_OUTPUT_DIR,fn)) for fn in infiles)
    #print "String: %s" % (s)
    #subprocess.check_call(["root", "-b", "-l", "-q","makeEFTGenPlots.C(\"%s\", %s)" % (outf,"{{{}}}".format(s))])
    subprocess.check_call(["root", "-b", "-l", "-q","makeEFTGenPlots.C(%s)" % ("{{{}}}".format(s))])

    print "Output Dir: %s" % (out_dir)

    #output_dir = os.path.join(GEN_PLOTS_DIR,"overlayTTW_NoTopLeptons")
    output_dir = os.path.join(GEN_PLOTS_DIR,out_dir)
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    imgs = regex_match(get_files('.'),["^h_.*\.png$"])
    move_files(files=imgs,target=output_dir)
    make_html(output_dir)

if __name__ == "__main__":
    main()
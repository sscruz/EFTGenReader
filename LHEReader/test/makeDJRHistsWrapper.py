import os
import subprocess
from EFTGenReader.GenReader.make_html import make_html
from EFTGenReader.GenReader.utils import clean_dir,get_files,move_files

#USE_SMPT = True; USE_REFPT = False  # SM base pt
USE_SMPT = False; USE_REFPT = True # Ref base pt

def make_plots(sample_loc,output_fn):
    if not os.path.exists(sample_loc):
        return

    files = []

    print 'Loc: {loc}'.format(loc=sample_loc)
    print 'Files:'
    for fn in os.listdir(sample_loc):
        fpath = os.path.join(sample_loc,fn)
        if not os.path.isfile(fpath):
            continue
        if fn.split('.')[-1] != 'root':
            continue
        files.append(fpath)
        print "\t{fn}".format(fn=fn)

    inputs_fn = 'infpaths.txt'
    with open(inputs_fn,'w') as f:
        for fpath in files:
            l = "{fpath}\n".format(fpath=fpath)
            f.write(l)

    macro_str = "makeDJRHists.C(\"{infs}\",\"{outf}\",{smpt},{refpt})".format(infs=inputs_fn,outf=output_fn,smpt=int(USE_SMPT),refpt=int(USE_REFPT))
    cmd = ['root','-b','-l','-q',macro_str]

    print 'Root Command: {0}'.format(' '.join(cmd))
    subprocess.check_call(['root','-b','-l','-q',macro_str])
    print ""


def main():

    dir_loc = '/hadoop/store/user/kmohrman/summaryTree_LHE/2019_08_14_addPtBranches/ttHJet-xqcutStudies-xqcut10qCutTests_analysisEtaCut-GEN/v3/'
    dir_list = [
        #("output_ttHJet_HanModel16DttllScanpointsqCut15_run1" , "xqcut10qCut15"),
        ("output_ttHJet_HanModel16DttllScanpointsqCut19_run1" , "xqcut10qCut19"),
        #("output_ttHJet_HanModel16DttllScanpointsqCut25_run1" , "xqcut10qCut25"),
    ]

    for dir_name,Cuts in dir_list:
        print dir_name, Cuts
        full_path = os.path.join(dir_loc,dir_name)
        print full_path
        outf = 'ttHJet_{Cuts}_'.format(Cuts=Cuts)
        make_plots(full_path,outf)

        #output_dir = '/afs/crc.nd.edu/user/k/kmohrman/www/EFT/DJRplots/WCtests/ctGscan/xqcut10qCut19'

        if USE_SMPT:
            #output_dir = '/afs/crc.nd.edu/user/k/kmohrman/www/EFT/DJRplots/testing/test_WCscan_SMbase'
            #output_dir = '/afs/crc.nd.edu/user/k/kmohrman/www/EFT/DJRplots/RwgtTests/basePtSM_xqcut10qCut19'
            output_dir = '/afs/crc.nd.edu/user/k/kmohrman/www/EFT/DJRplots/RwgtTests/PtEtaCutsTests/PtCuts_EtaCut2p5/3lep/basePtSM_xqcut10qCut19_PtEtaCuts2lep'
        else:
            #output_dir = '/afs/crc.nd.edu/user/k/kmohrman/www/EFT/DJRplots/testing/test_WCscan_Refbase'
            #output_dir = '/afs/crc.nd.edu/user/k/kmohrman/www/EFT/DJRplots/RwgtTests/basePtRef_xqcut10qCut19'
            output_dir = '/afs/crc.nd.edu/user/k/kmohrman/www/EFT/DJRplots/RwgtTests/PtEtaCutsTests/PtCuts_EtaCut2p5/3lep/basePtRef_xqcut10qCut19_PtEtaCuts2lep'

        if not os.path.exists(output_dir):
            os.mkdir(output_dir);

    imgs = get_files('.',targets=["^ttHJet_.*\.png$"])
    move_files(files=imgs,target=output_dir)
    make_html(output_dir)

if __name__ == "__main__":
    main()


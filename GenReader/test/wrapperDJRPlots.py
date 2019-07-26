import os
import subprocess
from EFTGenReader.GenReader.make_html import make_html
from EFTGenReader.GenReader.utils import clean_dir,get_files,move_files

# Example dir: "/hadoop/store/user/awightma/postLHE_step/2019_04_19/ttHJet-xqcutStudies/v3/mAOD_step_ttHJet_HanModelxqcut25_run0"

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

    macro_str = "plotdjr.C(\"{infs}\",\"{outf}\")".format(infs=inputs_fn,outf=output_fn)
    cmd = ['root','-b','-l','-q',macro_str]

    print 'Root Command: {0}'.format(' '.join(cmd))
    subprocess.check_call(['root','-b','-l','-q',macro_str])
    print ""

def main():
    xqcut_map = {
        'xqcut10': 'HanModel16DttllScanpoints',
        'xqcut25': 'HanModelxqcut25',
        'xqcut35': 'HanModelxqcut35',
        'xqcut45': 'HanModelxqcut45'
    }

    qCut_map = {
        'qCut19': 'v1',
        'qCut30': 'v2',
        'qCut45': 'v3'
    }

    hpath1 = '/hadoop/store/user/awightma/postLHE_step/2019_04_19/ttHJet-xqcutStudies'
    run1 = 'run0'
    lst1 = [
        #('xqcut25','qCut19'),
        ('xqcut25','qCut30'),
        ('xqcut25','qCut45'),

        #('xqcut35','qCut19'),
        #('xqcut35','qCut30'),
        ('xqcut35','qCut45'),

        #('xqcut45','qCut19'),
        #('xqcut45','qCut30'),
        ('xqcut45','qCut45'),
    ]

    hpath2 = '/hadoop/store/user/kmohrman/postLHE_step/2019_04_19/ttHJet-xqcutStudies'
    run2 = 'run1'
    lst2 = [
        ('xqcut10','qCut19'),
        ('xqcut10','qCut30'),
        ('xqcut10','qCut45')
    ]

    hpath = hpath1
    run   = run1
    lst   = lst1
    for xqcut,qCut in lst:
        v = qCut_map[qCut]
        d = 'mAOD_step_ttHJet_{grp}_{run}'.format(grp=xqcut_map[xqcut],run=run)
        loc = os.path.join(hpath,v,d)
        outf = 'ttHJet_{xqcut}{qCut}_djrplot.png'.format(qCut=qCut,xqcut=xqcut)
        make_plots(loc,outf)

    hpath = hpath2
    run   = run2
    lst   = lst2
    for xqcut,qCut in lst:
        v = qCut_map[qCut]
        d = 'mAOD_step_ttHJet_{grp}_{run}'.format(grp=xqcut_map[xqcut],run=run)
        loc = os.path.join(hpath,v,d)
        outf = 'ttHJet_{xqcut}{qCut}_djrplot.png'.format(qCut=qCut,xqcut=xqcut)
        make_plots(loc,outf)

    output_dir = '/afs/crc.nd.edu/user/a/awightma/www/eft_stuff/misc/djr_plots/ttHJet_2019-07-26_1319'
    imgs = get_files('.',targets=["^ttHJet_.*\.png$"])
    move_files(files=imgs,target=output_dir)
    make_html(output_dir)


if __name__ == "__main__":
    main()
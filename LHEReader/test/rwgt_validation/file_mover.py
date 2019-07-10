import os

from EFTGenReader.GenReader.utils import regex_match
from EFTGenReader.GenReader.make_html import make_html

USER_DIR = os.path.expanduser('~')

def get_files(tdir,ext=''):
    ret = []
    if not os.path.exists(tdir):
        return ret
    for fn in os.listdir(tdir):
        fpath = os.path.join(tdir,fn)
        if not os.path.isfile(fpath):
            continue
        elif len(fn.split('.')) < 2:
            continue
        elif len(ext) and fn.split('.')[-1] != ext:
            continue
        ret.append(fn)
    return ret

def move_files(files,target,indent=''):
    width = len(max(files,key=len))
    for fpath in files:
        src_path,fn = os.path.split(fpath)
        src = os.path.join(src_path,fn)
        dst = os.path.join(target,fn)
        #print "{indent}[{0:0>{w1}}] {1:<{w2}} --> {2}".format(len(src),src,dst,w1=2,w2=width,indent=indent)
        print "{indent}Moving: {file}".format(indent=indent,file=fn)
        os.rename(src,dst)

# Removes files from tdir which match any of the regex in targets list
def clean_dir(tdir,targets,indent=''):
    fnames = regex_match(get_files(tdir),targets)
    if len(fnames) == 0: return
    print "{indent}Removing files from: {dir}".format(dir=tdir,indent=indent)
    print "{indent}\tTargets: {tar}".format(tar=targets,indent=indent)
    for fn in fnames:
        fpath = os.path.join(tdir,fn)
        print "{indent}\tRemoving: {file}".format(file=fn,indent=indent)
        #os.remove(fpath)

def group_by(fn_lst,typ):
    # Returns: {'p1': [], 'p2': [], ...}
    ret = {}
    for fn in fn_lst:
        p,c,r = fn.split('_')[:3]
        if typ == 'process':
            grp = p
        elif typ == 'tag':
            grp = c
        elif typ == 'run':
            grp = r
        else:
            raise RuntimeError("Invalid group type: {type}".format(type=typ))
        if not ret.has_key(grp):
            ret[grp] = []
        ret[grp].append(fn)
    return ret

def move_rwgt_plots():
    force = False
    indent= ' '*4
    src_dir = "read_lhe_outputs"
    dst_dir = os.path.join(USER_DIR,'www/eft_stuff/misc/rwgt_validation_plots/plus_1jet-GEN_withHanModel1DRef')

    if not os.path.exists(dst_dir):
        print "Unknown destination directory: {dir}".format(dir=dst_dir)
        return


    imgs = get_files(src_dir,ext='png')
    grps = group_by(imgs,'process')
    for g,lst in grps.iteritems():
        fpath_lst = [os.path.join(src_dir,x) for x in lst]
        output_dir = os.path.join(dst_dir,g)

        print "Group: {grp}".format(grp=g)
        print "{indent}Output directory: {dir}".format(dir=output_dir,indent=indent)
        if not os.path.exists(output_dir):
            print "{indent}Making empty directory: {dir}".format(dir=output_dir,indent=indent)
            os.mkdir(output_dir)
        if force:
            clean_dir(output_dir,["^.*\.png$","^index.html$"],indent=indent)
        move_files(fpath_lst,output_dir,indent=indent)
        make_html(output_dir)


def main():
    move_rwgt_plots()

if __name__ == "__main__":
    main()
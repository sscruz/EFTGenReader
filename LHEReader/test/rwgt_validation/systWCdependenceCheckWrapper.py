import subprocess

def main():

    # Right now systWCdependenceCheck.C assumes these dirs have the 5 cpt axis scan files in it!!!
    wc_name = "cpt"
    f_name = "dirpaths_systCheck.txt" 
    qCut_f_name = "dirpaths_systCheck_qCutDirs.txt"

    subprocess.check_call(['root','-b','-l','-q','systWCdependenceCheck.C+(\"%s\",\"%s",\"%s")' % (wc_name,f_name,qCut_f_name)])

main()

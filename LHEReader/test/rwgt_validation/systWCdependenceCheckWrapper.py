import subprocess

def main():

    # If specify cpt, systWCdependenceCheck.C assumes these dirs have the 5 cpt axis scan files in it
    wc_name = "cpt"
    f_name = "dirpaths_systCheck.txt" 
    qCut_f_name = "dirpaths_systCheck_qCutDirs.txt"

    # If specify cpt, systWCdependenceCheck.C assumes these dirs have the 7 cpt axis scan files in it
    #wc_name = "ctG"
    #qCut_f_name = "dirpaths_systCheck_qCutDirs_ctG.txt"
    #f_name = qCut_f_name

    subprocess.check_call(['root','-b','-l','-q','systWCdependenceCheck.C+(\"%s\",\"%s",\"%s")' % (wc_name,f_name,qCut_f_name)])

main()

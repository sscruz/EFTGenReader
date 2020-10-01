import subprocess

def main():

    # Specify ttW cpQ3
    #proc_name = "ttW"
    #wc_name = "cpQ3"
    #f_name = "dirpaths_systCheck_qCutDirs_ttW-cqQ3.txt"
    #f_name_qCut = "dirpaths_systCheck_qCutDirs_ttW-cqQ3.txt"
    #f_name_0p = "dirpaths_systCheck_0p_ttW-cpQ3.txt"

    # If specify cpt, systWCdependenceCheck.C assumes these dirs have the 5 cpt axis scan files in it
    #proc_name = "tth"
    #wc_name = "cpt"
    #f_name = "dirpaths_systCheck.txt" 
    #f_name_qCut = "dirpaths_systCheck_qCutDirs_cpt.txt" 
    #f_name_0p = "dirpaths_systCheck_qCutDirs_cpt_0p.txt"

    # If specify cpt, systWCdependenceCheck.C assumes these dirs have the 7 cpt axis scan files in it
    proc_name = "tth"
    wc_name = "ctG"
    f_name_qCut = "dirpaths_systCheck_qCutDirs_ctG.txt"
    f_name = f_name_qCut
    f_name_0p = ""

    #subprocess.check_call(['root','-b','-l','-q','systWCdependenceCheck.C+(\"%s\",\"%s",\"%s",\"%s")' % (wc_name,f_name,f_name_qCut,f_name_0p)])
    subprocess.check_call(['root','-b','-l','-q','systWCdependenceCheck.C+(\"%s\",\"%s\",\"%s",\"%s",\"%s")' % (proc_name,wc_name,f_name,f_name_qCut,f_name_0p)])

main()

import subprocess

def main():

    f_name = "dirpaths_systCheck.txt" # Right now systWCdependenceCheck.C assumes this dir has the 5 cpt axis scan files in it!!!
    subprocess.check_call(['root','-b','-l','-q','systWCdependenceCheck.C+(\"%s\")' % (f_name)])

main()

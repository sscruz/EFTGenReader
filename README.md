# Set up Repo
To install this package:

    cd $CMSSW_BASE/src/
    cmsenv
    git clone https://github.com/Andrew42/EFTGenReader.git EFTGenReader
    scram b

# Run the code
To run the code:

    cd $CMSSW_BASE/src/EFTGenReader/Reader/test
    cmsRun EFTGenReader_cfg.py
    python wrapperEFTGenPlots.py dst_dir_name f1.root f2.root f3.root
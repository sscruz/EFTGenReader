# Set up Repo
To install this package:

    cd $CMSSW_BASE/src/
    cmsenv
    git clone https://github.com/Andrew42/EFTGenReader.git EFTGenReader
    scram b
In the `wrapperEFTGenPlots.py` file, make sure the `GEN_PLOTS_DIR` variable points to the directory that you want the output directory to be created.

# Run the code
If reading files from the grid (e.g. for central samples) make sure to initialize your voms proxy
    
    voms-proxy-init -voms cms -valid 192:00

To run the code:

    cd $CMSSW_BASE/src/EFTGenReader/Reader/test
    cmsRun EFTGenReader_cfg.py
    python wrapperEFTGenPlots.py output_dir_name output/f1.root output/f2.root output/f3.root
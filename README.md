### Set up the Repo
To install this package:

    cd $CMSSW_BASE/src/
    cmsenv
    git clone https://github.com/Andrew42/EFTGenReader.git EFTGenReader
    scram b
**Note:** In the `wrapperEFTGenPlots.py` file, make sure the `GEN_PLOTS_DIR` variable points to the directory that you want the output directory to be created.

### Run the code
If reading files from the grid (e.g. for central samples) make sure to initialize your voms proxy first
    
    voms-proxy-init -voms cms -valid 192:00

To run the code:

    cd $CMSSW_BASE/src/EFTGenReader/GenReader/test
    cmsRun EFTGenReader_cfg.py json_sample_name
    python wrapperEFTGenPlots.py --web-dir=/absolute/or/relative/path/to/some/output/directory --dir=sub-directory_name fpath1 fpath2 etc...

Where `json_sample_name` is the name of an entry in the `datasets.json` file. The non-DAS samples are only reachable from earth. `fpath1`, `fpath2`, etc. are the (relative) file paths to the root files produced by running `EFTGenReader_cfg`. Each root file passed to `wrapperEFTGenPlots` will have the same named histograms overlayed with one another on a single plot.

### cmsRun Options
The `EFTGenReader_cfg` supports a number of command-line options for configuring the cmsRun job, which are implemented using the CMSSW `VarParsing` package:
```
Command-line Options:
    dataset=Name
        Name of the dataset as it appears in the JSON file
    test=True
        Indicates that the output root file should use the dummy TEST name
    debug=True
        Run in debug mode
    normType
        How to normalize the histograms; 0 - no norm, 1 - unit norm (default), 2 - xsec norm
    intgLumi
        Value to scale all histograms by (no effect for unit norm mode)
Examples:
    # basic running
    cmsRun EFTGenReader_cfg.py dataset=central_ttH
    # run in test mode
    cmsRun EFTGenReader_cfg.py test=True dataset=central_ttH
    # perform no histogram normalizations
    cmsRun EFTGenReader_cfg.py normType=0 dataset=central_ttH
```

### chainGenReader
The python script `chainGenReader.py` can be used to do multiple `cmsRun` commands in series. All of the config options available to `EFTGenReader_cfg` can be specified at the command-line for `chainGenReader`, which will then be passed on to the underlying `cmsRun` call. All samples in the chain will be ran with the same command-line options. Type `python chainGenReader.py --help` for a complete list of available options.

### Misc. Other
The information used to find the samples is stored in the `datasets.json` file located in `GenReader/data/JSON`. Currently, the code is only able to run on samples that are located locally on hadoop, or on DAS. Files are found by looking in the directory specified by the `loc` entry. DAS files are found by querying via the command line utility `dasgoclient`, which should be available from in any recent `CMSSW` release.

Additional datasets can be added/modified/removed from the `datasets.json` file either by hand (not recommended) or by making use of the `update_datasets.py` script. The script shows some examples of how to add/modify/remove entries from the `dataset.json` file via the `DatasetHelper` class. Currently, there is no way to edit the JSON file entirely from the command-line.

**Note:** The versions of the `TH1EFT` and `WCFit` classes in this repo have been modified relative to the ones in the `ttH-Multilepton` repository. Here we have re-enabled the storing of the fit errors as well as switched back to the alternative `TH1EFT::Scale` member function.
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
    cmsRun EFTGenReader_cfg.py
    python wrapperEFTGenPlots.py output_dir_name fpath1 fpath2 etc...

Where `fpath1`, `fpath2`, etc. are the (relative) file paths to the root files produced by running `EFTGenReader_cfg`. Each root file passed to `wrapperEFTGenPlots` will have the same named histograms overlayed with one another on a single plot.

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

### Misc. Other
The information used to find the samples is stored in the `datasets.json` file located in `Reader/data/JSON`. Currently, the code is only able to run on samples that are located locally on hadoop, or on DAS. Local files are found by looking in the directory specified by the `loc` entry. DAS files are found by querying via the command line utility `dasgoclient`, which should be available from in any recent `CMSSW` release. The DAS query is specified via the `dataset` entry in the JSON file.

**Note:** Currently, all DAS entries have a value of `null` for the `loc` entry and all local entries have an empty string value for the `dataset` entry. In the future, this may be changed to consolidate the two into a single entry in the JSON file.

Additional datasets can be added/modified/removed from the `datasets.json` file either by hand or by making use of the `update_datasets.py` script. The script shows some examples of how to add/modify/remove entries from the `dataset.json` file via the `DatasetHelper` class. Currently, there is no way to edit the JSON file entirely from the command-line.
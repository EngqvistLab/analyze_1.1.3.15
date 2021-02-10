# Experimental exploration of the enzyme class EC 1.1.3.15
Code used in the study outlined in our pre-print "Experimental investigation of enzyme functional annotations reveals extensive annotation error" (https://www.biorxiv.org/content/10.1101/2020.12.18.423474v1). The code is free to use and modify under GPLv3, but we do ask that you cite our article if you make use of it.

![Figure 2](/results/figures/figure2_ela.png.png)

## Requirements and usage
The scripts in this repository are in the form of Jupyter Notebooks and rely on Python (3.7). A full list of required packages, and their versions, are included in the `./code/environment.yml` file. This file can be used as a reference to manually install the libraries using pip, but it is by far easier to make use of Miniconda (https://docs.conda.io) to create an environment with all the libraries inside.

```bash
conda env create -f envrionment.yml
conda env activate annotation_error
```

The Notebooks are abundantly annotated and should hopefully be self-explanatory. All data files needed to produce the publication figures (`./code/2_anlyze_data/1_analyze_1.1.3.15.ipynb`) are included in the respository. This notebook is computationally light. To re-run the full bioinformatic analsis on BRENDA (`./code/1_prepare_data/`) is computationally expensive and takes about one week on a machine with 12 cores. Additionally, the full analysis requires the downloading of some data files from the Zenodo data repository (). The downloaded should be carried out automatically when running the notebooks. *A pre-requisite for this* is that your system has `wget` and `unzip` installed (typically available in Unix systems). An alternative is to manually download the files and extract them.


Finally, the analysis relies on accessory code from some of our other repositories (https://github.com/EngqvistLab/orgtools, https://github.com/EngqvistLab/brenparse, https://github.com/EngqvistLab/UniRep50, https://github.com/mengqvist/wsvg). These will be automatically installed if using Miniconda and the environment.yml file.


### Project Organization
    ├── LICENSE
    ├── README.md                          <- The top-level README for developers using this project.
    │
    ├── code
    │   ├── 1_prepare_data                                    <- Folder holding the main data processing scripts
    │   │   ├── server_scripts                                <- Folder holding server scripts which were used to paralellize the analysis
    │   │   ├── 1_parse_brenda_html_and_fasta.ipynb           <- Notebook for obtaining UniProt identifiers, cluster sequences, and domain of life information
    │   │   ├── 2_compute_entire_brenda_identities.ipynb      <- Notebook for computing sequence identity between cluster representatives and closest characterized
    │   │   ├── 3_get_brenda_domains.ipynb                    <- Notebook for obtaining Pfam domains for cluster representatives
    │   │   └── 4_get_domain_and_identity_info_for_seqs.ipynb <- Notebook for obtaining Pfam domains for all sequences in 1.1.3.15, from two database versions
    │   │
    │   ├── 2_analyze_data                                    <- Folder containing data analysis script
    │   │   └── 1_analyze_1.1.3.15.ipynb                      <- Notebook for analyzing the prepared data files and make publication figures
    │   │
    │   └── environment.yml                                   <- File specifying the Miniconda environment used
    │
    ├── data
    │   ├── final                           <- Folder holding output files of processed data
    │   │   ├── brenda_2017_1               <- Folder holding data relating to the 2017.1 version of the BRENDA database
    │   │   └── brenda_2019_2               <- Folder holding data relating to the 2019.2 version of the BRENDA database
    │   │
    │   ├── intermediate                    <- Folder holding intermediate data that has been transformed
    │   │   ├── brenda_2017_1               <- Folder holding data relating to the 2017.1 version of the BRENDA database
    │   │   ├── brenda_2019_2               <- Folder holding data relating to the 2019.2 version of the BRENDA database
    │   │   └── BRENDA_for_paper            <- Folder holding data relating to the bioinformatic similarity analysis
    │   │
    │   ├── raw_external                    <- Folder holding raw unmodified external data
    │   │   ├── brenda_2017_1               <- Folder holding data relating to the 2017.1 version of the BRENDA database
    │   │   └── brenda_2019_2               <- Folder holding data relating to the 2019.2 version of the BRENDA database
    │   │
    │   └── raw_internal                    <- Folder holding raw unmodified experimental data
    │       └── experiments                 <- Folder 
    │           ├── alternative_activities  <- Folder holding enzyme activity data for the enzymes with non-canonical domains
    │           ├── batch_1                 <- Folder holding enzyme activity data from the screen for canonical EC 1.1.3.15 activity
    │           ├── batch_2                 <- Folder holding enzyme activity data from the screen for canonical EC 1.1.3.15 activity
    │           └── detailed                <- Folder holding information about enzyme expression, solubility, and activity
    │
    ├── doc                                 <- Folder holding the pre-print article
    │
    └── reports            
        └── figures                         <- Folder holding graphics and figures to be used in reporting


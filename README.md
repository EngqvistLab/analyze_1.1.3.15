# Experimental exploration of the enzyme class EC 1.1.3.15
Code used in the study outlined in our pre-print "Experimental investigation of enzyme functional annotations reveals extensive annotation error" (https://www.biorxiv.org/content/10.1101/2020.12.18.423474v1). The code is free to use and modify under GPLv3, but we do ask that you cite our article if you make use of it.

## xyz



### Project Organization
    ├── LICENSE
    ├── README.md                    <- The top-level README for developers using this project.
    │
    ├── code
    │   ├── 1_prepare_data                                    <- Folder
    │   │   ├── server_scripts                                <- Folder
    │   │   ├── 1_parse_brenda_html_and_fasta.ipynb           <- Notebook
    │   │   ├── 2_compute_entire_brenda_identities.ipynb      <- Notebook
    │   │   ├── 3_get_brenda_domains.ipynb                    <- Notebook
    │   │   └── 4_get_domain_and_identity_info_for_seqs.ipynb <- Notebook
    │   │
    │   ├── 2_analyze_data                                    <- Folder
    │   │   ├── 1_analyze_1.1.3.15.ipynb                      <- Notebook
    │   │   └── 2_plot_molecules.ipynb                        <- Notebook 
    │   │
    │   └── environment.yml                                   <- File
    │
    ├── data
    │   ├── final                    <- Folder holding output files of processed data
    │   │   ├── brenda_2017_1        <- Folder
    │   │   ├── brenda_2019_2        <- Folder
    │   │   ├── experiments          <- Folder
    │   │   ├── first_selection      <- Folder
    │   │   └── second_selection     <- Folder
    │   │
    │   ├── intermediate             <- Folder holding intermediate data that has been transformed
    │   │   ├── brenda_2017_1        <- Folder
    │   │   ├── brenda_2019_2        <- Folder
    │   │   └── BRENDA_for_paper     <- Folder
    │   │
    │   ├── raw_external             <- Folder holding raw unmodified external data
    │   │   ├── brenda_2017_1        <- Folder
    │   │   ├── brenda_2019_2        <- Folder
    │   │   ├── missing_sequences    <- Folder
    │   │   └── pfam                 <- Folder
    │   │
    │   └── raw_internal             <- Folder holding raw unmodified experimental data
    │       ├── experiments          <- Folder
    │       │   ├── batch_1          <- Folder containing enzyme activity data
    │       │   ├── batch_2          <- Folder containing enzyme activity data
    │       │   └── detailed         <- Folder containing enzyme activity data
    │       │
    │       ├── first_selection      <- Folder containing the first selection of proteins
    │       └── second_selection     <- Folder containing the second selection of proteins
    │
    ├── doc                          <- Folder containing the pre-print article
    │
    └── reports            
        └── figures                  <- Folder holding graphics and figures to be used in reporting


### Microfluidic enzyme assays
The Fluidigm Biomark HD microfluidic qPCR system was adapted to carrying out enzyme kinetic studies using a amplex-red based assay. This repository contains data, code and figures used in the study (https://www.biorxiv.org/content/10.1101/2020.09.18.303248v1).

![Graphical abstract](/results/figures/TOC.png)


### Project Organization
    ├── LICENSE
    ├── README.md                <- The top-level README for developers using this project.
    ├── code
    │   ├── flex_six_first_test  <- Code to analyze the FlexSix data.
    │   ├── flex_six_km_test     <- Code to analyze the FlexSix data.
    │   ├── plate_km_test        <- Code to analyze the plate data.
    │   ├── pub_figs             <- Code to generate color grid in Figure 1.
    │   └── resorufin_std_curves <- Code to analyze standard curves in FlexSix chip and microplate.
    │
    ├── data
    │   ├── final                    <- Final output files of processed data.
    │   │   ├── flex_six_first_test  <- The initial test, exploring detergents and BSA concentrations.
    │   │   ├── flex_six_km_test     <- Michaelis-Menten kinetics carried out on the FlexSix chip.
    │   │   ├── plate_km_test        <- Michaelis-Menten kinetics carried out in microplate.
    │   │   └── resorufin_std_curves <- Standard curves obtained using resorufin.
    │   │
    │   ├── intermediate             <- Intermediate data that has been transformed.
    │   │   ├── flex_six_first_test  <- The initial test, exploring detergents and BSA concentrations.
    │   │   ├── flex_six_km_test     <- Michaelis-Menten kinetics carried out on the FlexSix chip.
    │   │   ├── plate_km_test        <- Michaelis-Menten kinetics carried out in microplate.
    │   │   └── resorufin_std_curves <- Standard curves obtained using resorufin.
    │   │
    │   └── raw_internal             <- The raw unmodified data.
    │       ├── flex_six_first_test  <- The initial test, exploring detergents and BSA concentrations.
    │       ├── flex_six_km_test     <- Michaelis-Menten kinetics carried out on the FlexSix chip.
    │       └── resorufin_std_curves <- Standard curves obtained using resorufin.
    │
    ├── doc                          <- Folder containing the pre-print article.
    │
    └── reports            
        └── figures                  <- Generated graphics and figures to be used in reporting

### Code
This folder contains jupyter-notebooks to run the bioinformatic analysis and generating figures for the article. The simplest way to ensure that all required libraries are installed is to make use of Miniconda (https://docs.conda.io/) and install libraries using the provided environment file.

```bash
conda env create -f packages.yml
```

The folder `1_prepare_data` contains notebooks for the bioinformatic analysis, which is computationally demanding and take about a week to run on a dektop with 24 cores. The `2_analyze_data` folder contains a single notebook which generates the article figures. This notbook runs quickly.





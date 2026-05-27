# kp-evolution-inference

Global comparative study of the evolution of AMR characters in *Klebsiella pneumoniae*. See preprint here! https://www.biorxiv.org/content/10.1101/2025.09.20.677523v2

DOI 10.5281/zenodo.20408311

<img width="1819" height="1246" alt="image" src="https://github.com/user-attachments/assets/3e94e843-342f-4596-b0d2-d158ef858156" />

This repo has three components. The first uses existing (Kleborate) data with HyperTraPS to infer the evolutionary pathways by which Kp acquires different AMR characters across different countries. The second processes newly-sequenced Kp genome data and performs evolutionary analysis on this new data set. The third (much smaller) carries out a synthetic control study examining the performance of HyperTraPS when the "source" data are generated from a reversible evolutionary process (HyperTraPS assumes irreversible changes). 

## Dependencies

Requires R with `hypertrapsct` https://github.com/StochasticBiology/hypertrapsct and a collection of other libraries (automatically generated list): `dplyr`, `ggplot2`, `ggpubr`, `wordcloud`, `lme4`, `phytools`, `countrycode`, `ggbeeswarm`, `ggrepel`, `ggupset`, `tidyr`, `tidyverse`

We also use `LINcoding.py` by Melanie Hennart, which requires Python with `numpy` (and common libraries copy, argparse, itertools, random, re)

The genome assembly needs RagTag https://github.com/malonge/RagTag and Mummer e.g. `brew install mummer`. Processing genomes requires Kleborate https://github.com/klebgenomics/Kleborate .

# General pipeline illustration

| kleborate-analysis/ | new-genomes/ | reversibility-control/ |
|---------------------|--------------|-------------------------|
| **Download and preprocess Kleborate data (hours)**<br>pipeline_preprocessing.R<br>- download.pathogen.watch<br>- preprocess.kleborate<br>- preprocess.metadata<br>- prepare.all.lincoding.files<br>- make_trees.sh (tens of minutes/hours)<br>[uses LINcoding.py, externally produced]<br>- get.tree.metrics (minutes)<br>↓ | **New genome data**<br>[not currently in repo, get from https://osf.io/36r45 ] and existing Kp records [from kleborate-analysis/]<br>↓ | **Synthetic control study on reversibility in HyperTraPS (tens of minutes)**<br>reversibility-control.R (tens of minutes) |
| **Parallelised inference using Kleborate data across countries (days–weeks)**<br>pipeline_run_all.R<br>- run_hypertraps.R<br>- test-phylo.R<br>↓ | **Identify Kp isolates and calculate ANIs (hours)**<br>run-new-genomes.sh<br>- align-batch.sh (minutes)<br>- process-pairs.sh (hours)<br>- process-reports.sh<br>↓ | |
| **Analysis and visualisation (minutes)**<br>- pipeline_plot_all.R [outputs all_models.Rdata, a good save point]<br>- igj.R<br>- igj-followup.R | **Get new transitions and do inference (hours)**<br>process-ANI-all.R<br>↓| |
| | **Plotting and analysis (minutes)**<br>process-ANI-all.R   | |



# Inference using Kleborate data 

In `kleborate-analysis/`.

## Running all models

To reproduce the project first download and preprocess the data with `pipeline_preprocessing.R`, which also invokes `make_trees.sh` to process phylogenetic trees.

The inference is parallelised. Invoke `pipeline_run_all.R` with the clusters you want to use, and run them on separate cores/threads with something like `screen`

```
screen -S cluster.1
Rscript script/pipeline_run_all.R 1
```

When everything is finished, `pipeline_plot_all.R` amalgamates and plots summaries of the data; `igj.R` and plotted with `igj-followup.R` take these outputs forward, forming the first few figures in the manuscript. `test-phylo.R` runs inference with and without phylogenetic information for example countries to assess the impact of phylogenetic correction.

# Analysis of new genome data 

In `new-genomes/`.

❗ _Manual steps required._ First, retrieve `new-kp-genomes.tar.gz` from https://osf.io/36r45 .  Unpack all the outputs into the `number138/` folder. Next, download Klebsiella genomes from Tanzania in FASTA format from here https://pathogen.watch/genomes/all?country=tz&genusId=570 . Unpack into `From_Olav_fixed/`. To pull AMR features in both cases, run

```
kleborate -a *.fasta -o kleborate_results -p kpsc --trim_headers
```

to extract AMR information from these genomes. 

* `align-batch.sh` aligns a collection of new isolates against the Klebsiella reference genome (plus lots of redundant stuff)
* `process-pairs.sh` checks which new isolates look like genuine Klebsiella and does pairwise `dnadiff` comparisons across these to extract ANI-like scores. Then pulls Olav's existing Tanzanian samples and creates (and splits) a script doing all-by-all `dnadiff` comparisons across old and new
* `process-reports.sh` pulls ANI values from all the resulting reports
* `process-ANI-all.R` does the tree building, new transition pulling, inference and prediction assays for this, and produces paper figures
* `run-new-genomes.sh` wraps and illustrates the pipeline

# Reversibility control study 

In `reversibility-control/`.

`reversibility-control.R` explores how much of an issue reversible transitions are in the source data

# Notes on reproducing figures

In the Zotero repository (not on Github, due to filesize limitations), the following files are included (organised by the section of the project that produces them)

| File | Description |
|------|-------------|
| all_models.Rdata | Precomputed output for all HyperTraPS inference across countries |
| igj-all-bubbles.csv | Precomputed "bubble plot" summaries across countries |
| kleborate-dichotomized.csv | Preprocessed KpAMR data across countries |
| samples-by-country.Rdata	| Counts of Kp samples by country |
| upset-details.Rdata	| Details of KpAMR sets by sample (for UpSet plot) |
| Romania.nwk | LINcode-tree linking Romanian samples |
| pca-data-frame.Rdata | Country-wise embedding of summaries in PCA space |
| pca-values-country.Rdata | Particular PCA values by country (for correlation analysis) |
| drug-covariate.Rdata	| PCA statistics with drug use as a covariate |
| prevalence-cor.Rdata | Prevalence/ordering correlation |
| fit-rev-irrev-5.Rdata | Irreversible vs reversible model fits for data subset |
| interaction-graph.Rdata | Inferred interactions between characters |
|--------|---------|
| tree-all.phy | Tree linking old and new Tanzanian Kp genomes |
| genome-feature-summary.Rdata | KpAMR features in old and new Tanzanian genomes |
| fitted-penalised-analysed-138-fit.RData	| Precomputed output for HyperTraPS inference in Tanzanian samples |
| onestep-predictions.Rdata | Predictions for next step from fitted Tanzanian models |
|------|------|
| new-rev-control.Rdata | Outputs from synthetic, reversible model fitting |

The connections to manuscript figures are:

| Figure | Code producing figure | Raw data / preprocessed "save points" |
|------|-------|------|
| Fig 1 | Tasks 7-8 igj-followup.R | samples-by-country.Rdata, upset-details.Rdata, Romania.nwk, igj-all-bubbles.csv, kleborate-dichotomized.csv, all_models.Rdata |
| Fig 2 | Tasks 2,3,7 igj-followup.R | igj-all-bubbles.csv, pca-data-frame.Rdata |
| Fig 3 | Task 4 igj-followup.R | pca-values-country.Rdata |
| Fig 4 | Tasks 5-6 igj-followup.R | drug-covariate.Rdata, samples-by-country.Rdata |
| Fig 5 | Tasks 3,5,6 process-ANI-all.R | tree-all.phy, onestep-predictions.Rdata, genome-feature-summary.Rdata |
| Supp Fig 1 | (A-B) reversibility.R, C Task 11 igj-followup.R | new-rev-control.Rdata, fit-rev-irrev-5.Rdata |
| Supp Fig 2 | Task 9 igj-followup.R | igj-all-bubbles.csv |
| Supp Fig 3 | Task 10 igj-followup.R | prevalence-cor.Rdata |
| Supp Fig 4 | Task 12 igj-followup.R | interaction-graph.Rdata |
| Supp Fig 5 | Task 4 igj-followup.R | pca-values-country.Rdata |
| Supp Fig 6 | Tasks 3,5 process-ANI-all.R | tree-all.phy	|







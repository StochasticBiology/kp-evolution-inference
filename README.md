# kp-evolution-inference

Global comparative study of *Klebsiella pneumoniae* evolution

<img width="800" alt="image" src="https://github.com/user-attachments/assets/c01b7130-6efd-418a-8724-98f657addc06" />

## Dependencies

Requires R with `hypertrapsct` https://github.com/StochasticBiology/hypertrapsct and a collection of other libraries (automatically generated list): `dplyr`, `ggplot2`, `ggpubr`, `wordcloud`, `lme4`, `phytools`, `countrycode`, `ggbeeswarm`, `ggrepel`, `ggupset`, `tidyr`, `tidyverse`

We also use `LINcoding.py` by Melanie Hennart, which requires Python with `numpy` (and common libraries copy, argparse, itertools, random, re)

The genome assembly needs RagTag https://github.com/malonge/RagTag and Mummer e.g. `brew install mummer`

# General pipeline illustration

<img width="800" alt="image" src="https://github.com/user-attachments/assets/6b6840c4-9a5f-4c4c-8a92-751b4e382ded" />

# Inference using Kleborate data 

In `kleborate-analysis/`.

## Running all models

To reproduce the project first download and preprocess the data with `pipeline_preprocessing.R`, which also invokes `make_trees.sh` to process phylogenetic trees.

The inference is parallelised. Invoke `pipeline_run_all.R` with the clusters you want to use, and run them on separate cores/threads with something like `screen`

```
screen -S cluster.1
Rscript script/pipeline_run_all.R 1
```

When everything is finished, `pipeline_plot_all.R` amalgamates and plots summaries of the data; `igj.R` and plotted with `igj-followup.R` take these outputs forward, forming the first few figures in the manuscript.

# Analysis of new genome data 

In `new-genomes/`.

* `align-batch.sh` aligns a collection of new isolates against the Klebsiella reference genome (plus lots of redundant stuff)
* `process-pairs.sh` checks which new isolates look like genuine Klebsiella and does pairwise `dnadiff` comparisons across these to extract ANI-like scores. Then pulls Olav's existing Tanzanian samples and creates (and splits) a script doing all-by-all `dnadiff` comparisons across old and new
* `process-reports.sh` pulls ANI values from all the resulting reports
* `process-ANI-all.R` does the tree building, new transition pulling, inference and prediction assays for this, and produces paper figures
* `run-new-genomes.sh` wraps and illustrates the pipeline

# Reversibility control study 

In `reversibility-control/`.

`reversibility-control.R` explores how much of an issue reversible transitions are in the source data



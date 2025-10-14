# kp-evolution-inference

Global comparative study of the evolution of AMR characters in *Klebsiella pneumoniae*

<img width="1819" height="1246" alt="image" src="https://github.com/user-attachments/assets/3e94e843-342f-4596-b0d2-d158ef858156" />

This repo has three components. The first uses existing (Kleborate) data with HyperTraPS to infer the evolutionary pathways by which Kp acquires different AMR characters across different countries. The second processes newly-sequenced Kp genome data and performs evolutionary analysis on this new data set. The third (much smaller) carries out a synthetic control study examining the performance of HyperTraPS when the "source" data are generated from a reversible evolutionary process (HyperTraPS assumes irreversible changes). 

## Dependencies

Requires R with `hypertrapsct` https://github.com/StochasticBiology/hypertrapsct and a collection of other libraries (automatically generated list): `dplyr`, `ggplot2`, `ggpubr`, `wordcloud`, `lme4`, `phytools`, `countrycode`, `ggbeeswarm`, `ggrepel`, `ggupset`, `tidyr`, `tidyverse`

We also use `LINcoding.py` by Melanie Hennart, which requires Python with `numpy` (and common libraries copy, argparse, itertools, random, re)

The genome assembly needs RagTag https://github.com/malonge/RagTag and Mummer e.g. `brew install mummer`. Processing genomes requires Kleborate https://github.com/klebgenomics/Kleborate .

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



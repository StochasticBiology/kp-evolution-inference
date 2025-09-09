# kp-evolution-inference

Global comparative study of *Klebsiella pneumoniae* evolution

Inference using Kleborate data `kleborate-analysis/`
----

## Dependencies

To run the code you must clone hypertraps-ct and put it in neighbouring folder to this git repository on your computer.

TODO: list all R/python dependencies

## Running all models

To reproduce the project first download and preprocess the data with the following commands:

```
Rscript script/pipeline_preprocessing.R
```

Then to recreate the newick trees run `make_trees.sh` in the exec-folder.

You can then do some statistics on the trees by using `Rscript script/plot_tree_metrics.R`.

To run all the models, set the appropriate number of clusters you want to use in `script/pipeline_run_all.R` and run them on separate cores/threads with something like `screen`

```
screen -S cluster.1
Rscript script/pipeline_run_all.R 1
```

Type ctrl+a then d to detach from the screen and leave it running in the background.

Finally to plot them do:

```
Rscript script/pipeline_plot_all.R
```

The outputs should be pulled with `igj.R` and plotted with `igj-followup.R`, forming the bulk of the figures for the paper. 

Analysis of new genome data `new-genomes/`
-----

Needs RagTag https://github.com/malonge/RagTag and Mummer `brew install mummer` (and Gnuplot). Also needs HyperTraPS-CT.

Currently in Dropbox/Sabrina

* `align-batch.sh` aligns a collection of new isolates against the Klebsiella reference genome (plus lots of redundant stuff)
* `compile-hits-1.sh` checks which new isolates look like genuine Klebsiella and does pairwise `dnadiff` comparisons across these to extract ANI-like scores
* `process-ANI.R` then takes these scores and construct an NJ phylogeny. It also takes Olav's characterisation of AMR features, curates the tree, and runs inference
* `compile-pairs.sh` pulls Olav's existing Tanzanian samples and creates (and splits) a script doing all-by-all `dnadiff` comparisons across old and new
* compile-all-hits and process-ANI-all take a first look at old and new together

Now looking at "number138" -- experiments 1, 3, 8 combined.
Here we have (in set1/)

* `process-pairs.sh` does the first part of `compile-hits-1.sh` then mirrors `compile-pairs.sh`, identifying all "138" Kp isolates then buidling a script to do all-by-all `dnadiff` comparisons
* `process-reports.sh` pulls ANI values from all the resulting reports
* `process-ANI-all.R` does the tree building, new transition pulling, inference and prediction assays for this, and produces paper figures
* `run.sh` illustrating the pipeline

Reversibility control study `reversibility-control/`
-----

`reversibility-control.R` explores how much of an issue reversible transitions are in the source data



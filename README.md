# kp-assembly

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
* `process-ANI-all.R` does the tree building, new transition pulling, inference and prediction assays for this
* `run.sh` illustrating the pipeline

`reversibility-control.R` explores how much of an issue reversible transitions are in the source data

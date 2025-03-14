# kp-assembly

Needs RagTag https://github.com/malonge/RagTag and Mummer `brew install mummer` (and Gnuplot)

Currently in Dropbox/Sabrina

* `align***.sh` XX
* `compile-hits-1.sh` checks which new isolates look like genuine Klebsiella and does pairwise `dnadiff` comparisons across these to extract ANI-like scores
* `process-ANI.R` then takes these scores and construct an NJ phylogeny. It also takes Olav's characterisation of AMR features, curates the tree, and runs inference
* `compile-pairs.sh` pulls Olav's existing Tanzanian samples and creates (and splits) a script doing all-by-all `dnadiff` comparisons across old and new

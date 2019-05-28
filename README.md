# MarginalEffects

This is a R-markdown document about controling bias when ttrying to estimate treatment effects from obervational data.

The document includes

* A brief introduction into Directed Acyclic Graphs (DAGs)
* A DAG-based description of bias from confounding, selective treatement/confounding by indication, self selection into a study, or loss to follow up
* A DAG based description of ways to control bias, in particular adjustment and weighting
* Simulations to compare the performance of adjustment and weigting in different scenarios (nor surpises here, if one is well familiar with DAGs and methods to control bias)
* A brief description of how to use implied conditional independencies to check if a hypothesized DAG is consistent with the oberved data

Disclaimer: This document reflects my current understanding and might contain inaccuracies.

The html document can be viewed [here](http://htmlpreview.github.io/?https://github.com/gbiele/MultipleBiases/blob/master/sim_conf_select_bias.html).

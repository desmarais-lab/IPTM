Functions useful for analyzing the 
Functions jointly analyze the history of interaction (between the sender and receiver) and the contents or topic of text (e.g. email, online messages) using Markov chain Monte Carlo (MCMC).
This package extends the functionality of topic models fit 

The subdirectory `pkg` contains the actual package. The package can be installed with [devtools](https://cran.r-project.org/package=devtools).

```{r}
devtools::install_github("bomin8319/IPTM", subdir = "pkg")
```

Functionality includes:

 - `partial_dependence` which computes the expected prediction made by the random forest if it were marginalized to only depend on a subset of the features. `plot_pd` plots the results.
 - `variable_importance` which computes feature importance for arbitrary loss functions, aggregated across the training data or for individual observations. This may also be used for subsets of the feature space in order to detect interactions.
 - `extract_proximity` and `plot_prox` which computes or extracts proximity matrices and plots them using a biplot given a matrix of principal components of said matrix.

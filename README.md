Functions jointly analyze the history of interaction (between the sender and receiver) and the contents or topic of text (e.g. email, online messages) using Markov chain Monte Carlo (MCMC).

The subdirectory `pkg` contains the actual package. The package can be installed with [devtools](https://cran.r-project.org/package=devtools).

```{r}
devtools::install_github("bomin8319/IPTM", subdir = "pkg")
```

Functionality includes:

 - `MCMC` which runs the Markov chain Monte Carlo algorithm to compute the time-weighted network statistics from the point process model of Perry and Wolfe (2013) and token-topic assignments of the documents over the corpus, given the interaction-pattern (IP) assignment of each document.

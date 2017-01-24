Functions jointly analyze the history of interaction (between the sender and receiver) and the contents or topic of text (e.g. email, online messages) using Markov chain Monte Carlo (MCMC).

The subdirectory `pkg` contains the actual package. The package can be installed with [devtools](https://cran.r-project.org/package=devtools).

```{r}
devtools::install_github("bomin8319/IPTM", subdir = "pkg")
```

Functionality includes:

 - `MCMC` which runs Markov chain Monte Carlo (MCMC) algorithm to compute the time-weighted network statistics from the point process model of Perry and Wolfe (2013) and token-topic assignments of the documents over the corpus, given the interaction-pattern (IP) assignment of each document.

 - `plot_beta` which plots the boxplot of network statistics for each interaction-pattern (IP).

 - `plot_topic` which plots the topic distributions for each interaction-pattern (IP).

- `table_IP` which generates the table that summarizes interaction-pattern (IP) and token-topic assignments. For each interaction-pattern (IP), the topics of highest proportion and their corresponding most likely words are included.
 
The last three functions above are used to make comparison between different interaction-patterns (IP).

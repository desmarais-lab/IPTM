#' @useDynLib IPTM
#' @import stats
#' @import grDevices
#' @import graphics
#' @importFrom Rcpp sourceCpp
#' @importFrom MCMCpack rdirichlet
#' @importFrom mvtnorm rmvnorm dmvnorm
#' @importFrom entropy entropy.empirical
#' @importFrom lda top.topic.words
#' @importFrom reshape melt
#' @importFrom coda mcmc
NULL

#' @title Selected
#' @description Convert matrix output of rmultinom into vector of chosen items
#'
#' @param nSample number of samples to draw
#' @param proportion numeric non-negative vector of length K, specifying the probability for the K classes
#'
#' @return The vector of chosen items from the rmultinom sampling
#'
#' @export
Selected = function(nSample, proportion) {
	# Convert matrix output of rmultinom into vector of chosen items
	#
	# Args 
	#  nSample number of samples to draw
	#  proportion non-negative vector, specifying the probability for the K classes
	#
	# Returns
	#  The vector of chosen items from the rmultinom sampling
	
	multinom.chosen  = rmultinom(nSample, 1, proportion)
	selected.item = vapply(1L:nSample, function(s) {
	               which(multinom.chosen[, s] == TRUE)
	               }, c(1))
	return(selected.item)
}

#' @title Netstats
#' @description Calculate the network statistics (dydadic, triadic, or degree) given the history of interactions
#'
#' @param historyIP list containing the weighted history of interactions for each IP
#' @param node nodelist containing the ID of nodes
#' @param sender specific timepoint that we are calculating the time difference from
#' @param netstat which type of network statistics to use ("dyadic", "triadic", "degree")
#' @param nIP total number of interaction patterns specified by the user
#'
#' @return list of network statistics for specific sender and all possible receivers for each IP
#'
#' @export
Netstats = function(historyIP, node, sender, netstat, nIP) {
  # Calculate the chosen network statistics given the history of interactions
  #
  # Args 
  #  historyIP list containing the weighted history of interactions for each IP
  #  node nodelist containing the ID of nodes
  #  sender specific timepoint that we are calculating the time difference from
  #  netstat which type of network statistics to use ("dyadic", "triadic", "degree")
  #  nIP total number of interaction patterns specified by the user
  #
  # Returns 
  #  list of network statistics for specific sender and all possible receivers for each IP
  intercept = rep(1, length(node) - 1)
  if ("dyadic" %in% netstat) {
    dyadic = Dyadic(historyIP, node, sender, nIP)
  } else {
  	dyadic = NULL
  }
  if ("triadic" %in% netstat) {
    triadic = Triadic(historyIP, node, sender, nIP)
  } else {
  	triadic = NULL
  }
  if ("degree" %in% netstat) {
    degree = Degree(historyIP, node, sender, nIP)
  } else {
  	degree = NULL
  }

  netstatmat = lapply(1:nIP, function(IP) {
  	matrix(cbind(intercept, dyadic[[IP]], triadic[[IP]], degree[[IP]]), nrow = length(node) - 1, byrow = FALSE)
  	})
  return(netstatmat)
}

#' @title TimeWeights
#' @description Calculate the time weights of the network statistics (dydadic, triadic, or degree) given the history of interactions
#'
#' @param historyTW list containing the cumulative count and time of the history of interactions (countmat and timemat)
#' @param node nodelist containing the ID of nodes
#' @param sender specific timepoint that we are calculating the time difference from
#' @param netstat which type of network statistics to use ("dyadic", "triadic", "degree")
#' @param mu parameter of decay
#'
#' @return List of time weights for specific sender and all possible receivers (for each network statistics)
#'
#' @export
TimeWeights = function(historyTW, node, sender, netstat, mu) {
  # Calculate the chosen network statistics given the history of interactions
  #
  # Args 
  #  historyTW list containing the cumulative count and time of the history of interactions (countmat and timemat)
  #  node nodelist containing the ID of nodes
  #  sender specific timepoint that we are calculating the time difference from
  #  netstat which type of network statistics to use ("dyadic", "triadic", "degree")
  #  mu parameter of decay
  #
  # Returns 
  #  List of time weights for specific sender and all possible receivers (for each network statistics)
  intercept = rep(1, length(node) - 1)
  if ("dyadic" %in% netstat) {
    dyadic = DyadicTW(historyTW, node, sender)
  } else {
  	dyadic = NULL
  }
  if ("triadic" %in% netstat) {
    triadic = TriadicTW(historyTW, node, sender)
  } else {
  	triadic = NULL
  }
  if ("degree" %in% netstat) {
    degree = DegreeTW(historyTW, node, sender)
  } else {
  	degree = NULL
  }

  timeweightmat = matrix(cbind(dyadic, triadic, degree), nrow = length(node) - 1, byrow = FALSE)
  timeweightmat[which(timeweightmat == 0)] = -Inf
  TW = cbind(intercept, (1 + exp(timeweightmat))^mu)
  return(TW)
}

#' @title GenerateDocs
#' @description Generate the documents according to the generative process
#'
#' @param nDocs number of documents to be generated
#' @param node nodelist containing the ID of nodes (ID starting from 1)
#' @param vocabulary all vocabularies used over the corpus
#' @param nIP total number of interaction patterns specified by the user
#' @param K total number of topics specified by the user
#' @param xi Poisson parameter for the number of words
#' @param netstat which type of network statistics to use ("dyadic", "triadic", "degree")
#' @param seed an integer value which controls random number generation
#'
#' @return List of edge and text generated according to the generative process
#'
#' @export
GenerateDocs = function(nDocs, node, vocabulary, nIP, K, xi,
						netstat = c("dyadic", "triadic", "degree"), seed = 1) {
  # Generate the documents according to the generative process
  #
  # Args 
  #  nDocs number of documents to be generated
  #  node nodelist containing the ID of nodes (ID starting from 1)
  #  vocabulary all vocabularies used over the corpus
  #  nIP total number of interaction patterns specified by the user
  #  K total number of topics specified by the user
  #  xi Poisson parameter for the number of words
  #  netstat which type of network statistics to use ("dyadic", "triadic", "degree")
  #  seed an integer value which controls random number generation
  #
  # Returns 
  #  List of edge and text generated according to the generative process

	set.seed(seed)
	alpha = 50 / K
  	mvec = rep(1 / K, K)
  	betas = 10
 	W = length(vocabulary)
  	nvec = rep(1, W) / W
  	phi = lapply(1:K, function(k) {
		rdirichlet(1, betas * nvec)
	})
	P = 1 + 2 * ("dyadic" %in% netstat) + 4 * ("triadic" %in% netstat) + 2 * ("degree" %in% netstat)
	b = lapply(1:nIP, function(IP) {
		rmvnorm(1, rep(0, P), diag(P))
	})
	
	mu = rgamma(1, 2, 1)
	delta = rbeta(1, 1, 10)
	
	currentC = sample(1:nIP, K, replace = TRUE)

	t.d = 0
	edge = list()
	text = list()
	p.d = matrix(NA, nrow = nDocs, ncol = nIP)
	
	options(warn = -1)
	for (d in 1:nDocs) {
		N.d = rpois(1, xi)
		theta.d = rdirichlet(1, alpha * mvec)	
		topic.d = Selected(max(1, N.d), theta.d)
		text[[d]] = rep(NA, N.d)
		if (N.d > 0) {
			for (n in 1:N.d){
				text[[d]][n] = vocabulary[Selected(1, phi[[topic.d[n]]])]
				}
			names(text[[d]]) = topic.d
		}

		p.d[d, ] = vapply(1:nIP, function(IP) {
			sum(topic.d %in% which(currentC == IP))
			}, c(1)) / max(1, N.d)	
		history.t = History(edge, p.d, node, t.d, nIP)
		network = lapply(node, function(i) {
				Netstats(history.t[[3]], node, i, netstat, nIP)
				})
		TW = lapply(node, function(i) {
			TimeWeights(history.t[1:2], node, i, netstat, mu)
		})
		X = lapply(1:nIP, function(IP) {
			lapply(node, function(i) {
				network[[i]][[IP]] / TW[[i]]
			})
		})
		XB = lapply(1:nIP, function(IP) {
			MultiplyXBList(X[[IP]], b[[IP]])
			})
		lambda = Reduce('+', lapply(1:nIP, function(IP) {
			p.d[d, IP] * exp(XB[[IP]])
			}))
		
		iJi = lapply(1:nrow(lambda), function(i) {
			vapply(1:ncol(lambda), function(j) {
				rbinom(1, 1, 1 - exp(-delta * lambda[i, j]))	
				}, c(1))
			})
		while (sum(unlist(iJi)) == 0) {
			iJi = lapply(1:nrow(lambda), function(i) {
			vapply(1:ncol(lambda), function(j) {
				rbinom(1, 1, 1 - exp(-delta * lambda[i, j]))
				}, c(1))
			})
		}
		gamma = vapply(1:length(node), function(i) {
				rbeta(1, 1, sum(iJi[[i]]))
				}, c(1))
		LambdaiJi = vapply(1:length(iJi), function(i) {
			sum(vapply(1:nIP, function(IP) {
			p.d[d, IP] * exp(gamma[i] * sum(iJi[[i]] * XB[[IP]][i, ])) 
			}, c(1)))
			}, c(1))		
		LambdaiJi[LambdaiJi == 1] = NA
		Time.inc = sapply(LambdaiJi, function(i) {
			rexp(1, i) 
			})
		i.d = which(Time.inc == min(Time.inc[!is.na(Time.inc)]))
		j.d = which(iJi[[i.d]] == 1)
		j.d[j.d >= i.d] = j.d[j.d >= i.d] + 1
		t.d = t.d + Time.inc[i.d]
		edge[[d]] = list(sender = i.d, receiver = j.d, timestamp = t.d)	
	}
	options(warn = 0)
	return(list(edge = edge, text = text))
}

#' @title AlphamvecOpt
#' @description Optimize the hyperparmeter vector (= alpha * mvec) given the current assignment of C and Z
#'
#' @param K total number of topics specified by the user
#' @param currentZ current state of the assignment of topics 
#' @param alpha Dirichlet concentration prior for topic distribution
#' @param mvec Dirichlet base prior for topic distribution
#'
#' @return The optimized value of the vector (= alpha * mvec)
#'
#' @export
AlphamvecOpt =  function(K, currentZ, alpha, mvec) {
	# Optimize the hyperparmeter vector given the current assignment of Z
	#
	# Args 
	#  K total number of topics specified by the user
	#  currentZ current state of the assignment of topics 
	#  alpha Dirichlet concentration prior for topic distribution
	#  mvec Dirichlet base prior for topic distribution
	#
	# Returns
	#  The optimized value of the vector (= alpha * mvec)
	
	final.vec = list()
	
	iter = 1
	current.vec = alpha * mvec
	n.word = mapply(length, currentZ)
	n.word.table = tabulate(n.word)
	nK.word.list = matrix(NA, nrow = length(currentZ), ncol = K)
	for (d in seq(along = currentZ)) {
		nK.word.list[d, ] = tabulate(currentZ[[d]], K)
	}
		nK.word.table = lapply(1L:K, function(k){
			  tabulate(nK.word.list[,k])
			  })
	while ((abs(alpha - sum(current.vec)) > 0.001) | (iter == 1)) {
	alpha = sum(current.vec)
	S = UpdateDenom(alpha, n.word.table)		
	s = UpdateNum(current.vec, nK.word.table)	
	current.vec = current.vec * s / S
	iter = iter + 1
	}
	final.vec = current.vec	
	return(final.vec)	
}	

#' @title TopicInEqZ
#' @description Calculate topic part of the equation used in multinomial sampling of Z
#'
#' @param K total number of topics specified by the user
#' @param currentZ current state of the assignment of topics 
#' @param alpha Dirichlet concentration prior for document-topic distribution
#' @param mvec Dirichlet base prior for document-topic distribution
#' @param document the id of specific edge of interest
#'
#' @return The vector of constants representing topic part of each K
#'
#' @export
TopicInEqZ = function(K, currentZ, alpha, mvec, document) {
	# Calculate topic-IP part of the equation used in multinomial sampling of Z
	#
	# Args 
 	#  K total number of topics specified by the user
	#  currentZ current state of the assignment of topics 
	#  alpha Dirichlet concentration prior for document-topic distribution
	#  mvec Dirichlet base prior for document-topic distribution
	#  document the id of specific edge of interest
	#
	# Returns 
	#  The vector of constants representing topic part of each K
	topics = currentZ[[document]]
	table.topics = tabulate(topics, K) 

	const = log(table.topics - as.numeric(table.topics > 0) + alpha * mvec)
	return(const)
}

#' @title MCMC
#' @description Iterate Markov Chain Monte Carlo (MCMC) algorithm to sequentially update the assignments of C, Z, and B 
#'
#' @param edge list of document information with 3 elements (element 1 sender, element 2 receiver, element 3 time in unix.time format)
#' @param node nodelist containing the ID of nodes (ID starting from 1)
#' @param textlist list of text (length=number of documents in total) containing the words in each document
#' @param vocabulary all vocabularies used over the corpus
#' @param nIP total number of interaction patterns specified by the user
#' @param K total number of topics specified by the user
#' @param delta.B proposal distribution variance parameter for beta 
#' @param outer size of outer iteration 
#' @param n1 size of first inner iteration for updates of interaction patterns C
#' @param n2 size of second inner iteration for updates of topics K
#' @param n3 size of third inner iteration for updates of beta B
#' @param burn iterations to be discarded at the beginning of beta chain
#' @param thin the thinning interval of beta chain
#' @param netstat which type of network statistics to use ("dyadic", "triadic", "degree")
#' @param seed an integer value which controls random number generation
#' @param plot to plot the convergence diagnostics or not (TRUE/FALSE)
#'
#' @return MMCMC output containing IP assignment, topic assignment, and (beta, mu, delta) chain 
#'
#' @export
MCMC = function(edge, node, textlist, vocabulary, nIP, K, delta.B, 
                outer = 1000, n1 = 3, n2 = 3, n3 = 3300, burn = 300, thin = 10, 
                netstat = c("dyadic", "triadic", "degree"), seed = 1, plot = TRUE) {
  # Iterate MCMC algorithm to sequentially update the assignments of C, Z, and B
  #
  # Args 
  #  edge list of document information with 3 elements (sender, receiver, time)
  #  node nodelist containing the ID of nodes (ID starting from 1)
  #  textlist list of text containing the words in each document
  #  vocabulary all vocabularies used over the corpus
  #  nIP total number of interaction patterns specified by the user
  #  K total number of topics specified by the user
  #  delta.B proposal distribution variance parameter for beta 
  #  outer size of outer iteration
  #  n1 size of first inner iteration for updates of interaction patterns C
  #  n2 size of second inner iteration for updates of topics K
  #  n3 size of third inner iteration for updates of beta B
  #  burn iterations to be discarded at the beginning of beta chain
  #  thin the thinning interval of beta chain
  #  netstat which type of network statistics to use ("dyadic", "triadic", "degree")
  #  seed an integer value which controls random number generation
  #  plot to plot the convergence diagnostics or not (TRUE/FALSE)
  #
  # Returns 
  #  MCMC output containing IP assignment, topic assignment, and (beta, mu, delta) chain
  
  set.seed(seed)

  # initialize alpha, mvec, delta, nvec, eta, lvec, and gammas
	alpha = 50 / K
  	mvec = rep(1 / K, K)
  	betas = 10
 	W = length(vocabulary)
  	nvec = rep(1, W) / W
  	phi = lapply(1:K, function(k) {
		rdirichlet(1, betas * nvec)
	})
		
	mu = rgamma(1, 2, 1)
	delta = rbeta(1, 1, 100)
	
	# initialize C, theta and Z
	currentC = sample(1:nIP, K, replace = TRUE)
	theta = rdirichlet(length(edge), alpha * mvec)
	currentZ = lapply(seq(along = edge), function(d) {
    if (length(textlist[[d]]) > 0) {
    	Selected(length(textlist[[d]]), theta[d, ])
    	} else {
    		Selected(1, theta[d, ])
    	}
    })

    # initialize beta
    sigma = 1
    P = 1 + 2 * ("dyadic" %in% netstat) + 4 * ("triadic" %in% netstat) + 2 * ("degree" %in% netstat)
	bmat = list()
	for (IP in 1L:nIP) {
		bmat[[IP]] = matrix(0, nrow = P, ncol = (n3 - burn) / thin)
  	}

  # to check the convergence  
  if (plot) {
  	logWZ.mat = c()							  
  	alpha.mat = c()
  	entropy.mat = c()
	}
	
  #start outer iteration
  for (o in 1L:outer) {
    cat("outer iteration = ", o, "\n")
    
    # Update the hyperparameter alpha and mvec
    vec = AlphamvecOpt(K, currentZ, alpha, mvec)
    alpha = sum(vec)
    mvec = vec / alpha
    
    # Data augmentation
    p.d = t(vapply(seq(along = edge), function(d) {
		vapply(1:nIP, function(IP) {
	 		sum(currentZ[[d]] %in% which(currentC == IP))
	 		}, c(1)) / length(currentZ[[d]])
    	}, rep(1, 3)))
    	
    iJi = list()
    lambda = list()
    LambdaiJi = list()
    gamma = list()
    for (d in seq(along = edge)) {
    	history.t = History(edge, p.d, node, as.numeric(edge[[d]][3]), nIP)
    	network = lapply(node, function(i) {
				Netstats(history.t[[3]], node, i, netstat, nIP)
				})
		TW = lapply(node, function(i) {
			TimeWeights(history.t[1:2], node, i, netstat, mu)
		})
		X = lapply(1:nIP, function(IP) {
			lapply(node, function(i) {
				network[[i]][[IP]] / TW[[i]]
			})
		})
		XB = lapply(1:nIP, function(IP) {
			MultiplyXBList(X[[IP]], b[[IP]])
			})
		lambda[[d]] = Reduce('+', lapply(1:nIP, function(IP) {
			p.d[d, IP] * exp(XB[[IP]])
			}))
		iJi[[d]] = lapply(1:length(node), function(i) {
			vapply(1:(length(node) - 1), function(j) {
				rbinom(1, 1, 1 - exp(-delta * lambda[[d]][i, j]))	
				}, c(1))
			})
		observedi = as.numeric(edge[[d]][1])
		observedj = as.numeric(unlist(edge[[d]][2]))
		observedj[observedj >= observedi] = observedj[observedj >= observedi] - 1	
		iJi[[d]][[observedi]] = tabulate(observedj, length(node) - 1)
		gamma[[d]] = vapply(1:length(node), function(i) {
			ifelse(sum(iJi[[d]][[i]]) > 0, rbeta(1, 1, sum(iJi[[d]][[i]])^(sum(iJi[[d]][[i]]))), 0)
		}, c(1))
		LambdaiJi[[d]] = vapply(1:length(node), function(i) {
			sum(vapply(1:nIP, function(IP) {
			p.d[d, IP] * exp(sum(iJi[[d]][[i]] * XB[[IP]][i, ])) 
			}, c(1)))
			}, c(1))
		LambdaiJi[[d]][abs(LambdaiJi[[d]] - 1) < 10^(-10)] = NA		
    }	 
    
    cat("inner iteration 1", "\n")
    textlist.raw = unlist(textlist)
    finalZlist.raw = unlist(currentZ)
    table.W = lapply(1L:K, function(k) {
      tabulate(textlist.raw[which(finalZlist.raw == k)], W)
      })
    
    for (i1 in 1L:n1) {
      for (d in seq(along = edge)) {
        textlist.d = textlist[[d]]
        if (length(textlist.d) > 0) {
        	topicpart.d = TopicInEqZ(K, currentZ, alpha, mvec, d)
       		wordpart.d = WordInEqZ(K, textlist.d, table.W, betas, nvec)
        } else {
        	topicpart.d = 0
        	wordpart.d = matrix(0, nrow = length(currentZ[[d]]), ncol = K)
        }
        nonemptyiJi = (LambdaiJi[[d]] * gamma[[d]])[!is.na(LambdaiJi[[d]])]
        sender = as.numeric(edge[[d]][1])
        nonemptyiJi2 = (LambdaiJi[[d]][-sender] * gamma[[d]][-sender])[!is.na(LambdaiJi[[d]][-sender])]
        timeinc = as.numeric(edge[[d]][3]) - ifelse(d > 1 , as.numeric(edge[[d-1]][3]), as.numeric(edge[[d]][3]))
        lambdapart.d = LambdaInEqZ(iJi[[d]], lambda[[d]], nonemptyiJi, nonemptyiJi2, delta, timeinc)
        for (w in 1L:length(currentZ[[d]])) {
          const.Z = topicpart.d + lambdapart.d + wordpart.d[w, ]
          zw.old = currentZ[[d]][w]
          zw.new = Selected(1, exp(const.Z))
          if (zw.new != zw.old) {
            currentZ[[d]][w] = zw.new
            topicpart.d = TopicInEqZ(K, currentZ, alpha, mvec, d)
            if (length(textlist.d) > 0) {	
            wordpart.d = WordInEqZ(K, textlist.d, table.W, betas, nvec)
            }
            table.W = lapply(1L:K, function(k) {
      		tabulate(textlist.raw[which(finalZlist.raw == k)], W)
      		})
      		p.d[d, ] = vapply(1:nIP, function(IP) {
	 		sum(currentZ[[d]] %in% which(currentC == IP))
	 		}, c(1)) / length(currentZ[[d]])
      		LambdaiJi[[d]] = vapply(1:length(node), function(i) {
			sum(vapply(1:nIP, function(IP) {
			p.d[d, IP] * exp(sum(iJi[[d]][[i]] * XB[[IP]][i, ])) 
			}, c(1)))
			}, c(1))
			LambdaiJi[[d]][abs(LambdaiJi[[d]] - 1) < 10^(-10)] = NA	
          }
        }
       }
      }

	
	
	# not touched yet	
    cat("inner iteration 1", "\n")
    lambda.i = list()
    for (i1 in 1L:n1) {
      # C update given Z and B - within each document d
      for (d in seq(along = edge)) {
        currentC.d = currentC[-d]
        edge.byC = SortedC(nIP, currentC.d, edge)
        for (IP in 1L:nIP) {
          history.t = Timediff(edge.byC[[IP]], node, edge[[d]][[3]], lambda)
          X.it = Netstats(history.t, node, edge[[d]][[1]], netstat)
          lambda.i[[IP]] = MultiplyXB(X.it, Beta.mat[[IP]][, (n3 - burn) / thin])
          }
          const.C = log(gammas) + BetaInEqC(nIP, lambda.i, edge[[d]]) +
            TopicInEqC(nIP, K, currentZ, alpha, mvec, d)
          currentC[d] = Selected(1, exp(const.C))
      }
    }
      
    if (plot) {
      entropy.mat = c(entropy.mat, entropy.empirical(currentC))
      alpha.mat = rbind(alpha.mat, alpha)
      logWZ.mat = c(logWZ.mat, 
                    logWZ(nIP, K, currentC[-empty.docs], currentZ[-empty.docs],
                    textlist[-empty.docs], table.W, alpha, mvec, beta, nvec))
      }
    
    cat("inner iteration 3", "\n")
    Beta.old = list()
    Beta.new = list()
    lambda.old = list()
    lambda.new = list()
    X.it.IP = list()
    edge.byC = SortedC(nIP, currentC, edge)
    for (IP in 1L:nIP) {
      Beta.old[[IP]] = Beta.mat[[IP]][, (n3 - burn) / thin]
      X.it.IP[[IP]] = list()
      for (d in seq(along = edge.byC[[IP]])) {
        history.t.IP = Timediff(edge.byC[[IP]], node, edge.byC[[IP]][[d]][[3]], lambda)
        X.it.IP[[IP]][[d]] = Netstats(history.t.IP, node, edge.byC[[IP]][[d]][[1]], netstat)
      }
      lambda.old[[IP]] = MultiplyXBList(X.it.IP[[IP]], Beta.old[[IP]])
      }
    for (i3 in 1L:n3) {
      for (IP in 1L:nIP) {
        Beta.new[[IP]] = rmvnorm(1, Beta.old[[IP]], (delta.B)^2 * diag(P))
        lambda.new[[IP]] = MultiplyXBList(X.it.IP[[IP]], Beta.new[[IP]])
        }
      prior = vapply(1L:nIP, function(IP) {
        dmvnorm(Beta.new[[IP]], rep(0, P), sigma^2 * diag(P), log = TRUE) -
          dmvnorm(Beta.old[[IP]], rep(0, P), sigma^2 * diag(P), log = TRUE)
        }, c(1))
      post = BetaInEqB(nIP, lambda.new, edge.byC) - BetaInEqB(nIP, lambda.old, edge.byC)
      loglike.diff = prior + post
      
      u = log(runif(nIP, 0, 1))
      for (IP in which((u < loglike.diff) == TRUE)) {
        Beta.old[[IP]]  = Beta.new[[IP]]
        lambda.old[[IP]] = lambda.new[[IP]]
        }
      
      if (i3 > burn && i3 %% (thin) == 0) {
        for (IP in 1L:nIP) {
          Beta.mat[[IP]][ , (i3 - burn) / thin] = Beta.old[[IP]]
          }
        }
      }
    }
  
  if (plot) {
    burnin = round(outer / 10)
    par(mfrow = c(2, 2))
  	plot(entropy.mat[-1L:-burnin], type = "l", 
  	     xlab = "(Outer) Iterations", ylab = "Entropy of IP")
  	abline(h = median(entropy.mat[-1L:-burnin]), lty = 1)
  	title("Convergence of Entropy")
  	matplot(alpha.mat[-1L:-burnin,], lty = 1, type = "l", col = 1L:nIP, 
  	        xlab = "(Outer) Iterations", ylab = "alpha")
  	abline(h = apply(alpha.mat[-1L:-burnin,], 2, median), lty = 1, col = 1L:nIP)
	  title("Convergence of Optimized alpha")
	  plot(logWZ.mat[-1L:-burnin], type = "l", 
	       xlab = "(Outer) Iterations", ylab = "logWZ")
	  abline(h = median(logWZ.mat[-1L:-burnin]), lty = 1)
	  title("Convergence of logWZ")
	  matplot(t(Beta.mat[[1]]), lty = 1, col = 1L:P, type = "l", 
	          main = "Traceplot of beta", xlab = "(Inner) Iterations", ylab = "")
	  abline(h = apply(t(Beta.mat[[1]]), 2, median), lty = 1, col = 1L:P)
	  }
  
  chain.final = list(C = currentC, Z = currentZ, B = Beta.mat)

  return(chain.final)
}

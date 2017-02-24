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

#' @title AlphamvecOpt
#' @description Optimize the hyperparmeter vector (= alpha * mvec) given the current assignment of C and Z
#'
#' @param nIP total number of interaction patterns specified by the user
#' @param K total number of topics specified by the user
#' @param currentC current state of the assignment of interaction patterns
#' @param currentZ current state of the assignment of topics 
#' @param alpha Dirichlet concentration prior for topic distribution
#' @param mvec Dirichlet base prior for topic distribution
#'
#' @return The optimized value of the vector (= alpha * mvec)
#'
#' @export
AlphamvecOpt =  function(nIP, K, currentC, currentZ, alpha, mvec) {
	# Optimize the hyperparmeter vector given the current assignment of C and Z
	#
	# Args 
	#  nIP total number of interaction pattern specified by the user
	#  K total number of topics specified by the user
	#  currentC current state of the assignment of interaction patterns
	#  currentZ current state of the assignment of topics 
	#  alpha Dirichlet concentration prior for topic distribution
	#  mvec Dirichlet base prior for topic distribution
	#
	# Returns
	#  The optimized value of the vector (= alpha * mvec)
	
	final.vec = list()
	corpus.byC = SortedZ(nIP, currentC, currentZ)
	
	for (IP in 1L:nIP) {
		iter = 1
		current.vec = alpha[IP] * mvec[,IP]
		n.word = mapply(length, corpus.byC[[IP]])
		n.word.table = tabulate(n.word)
		nK.word.list = matrix(NA, nrow = length(corpus.byC[[IP]]), ncol = K)
		for (d in seq(along = corpus.byC[[IP]])) {
			if (length(corpus.byC[[IP]][[d]]) > 0) {
			nK.word.list[d, ] = tabulateC(corpus.byC[[IP]][[d]], K)
			} else {
				nK.word.list[d, ] = rep(0, K)
			}
		}
		nK.word.table = lapply(1L:K, function(k){
			  tabulate(nK.word.list[,k])
			  })
	while ((abs(alpha[IP] - sum(current.vec)) > 0.001) | (iter == 1)) {
	alpha[IP] = sum(current.vec)
	S = UpdateDenom(alpha[IP], n.word.table)		
	s = UpdateNum(current.vec, nK.word.table)	
	current.vec = current.vec * s / S
	iter = iter + 1
	if (sum(is.na(current.vec)) > 0) {
		stop("ERROR: choose smaller nIP")
		}
	}
	final.vec[[IP]] = current.vec
	}	
	return(final.vec)	
}	

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

#' @title BetaInEqC
#' @description Calculate Beta part of the equation used in multinomial sampling of C
#'
#' @param nIP total number of interaction patterns specified by the user
#' @param lambda.i stochastic intensity of specific sender "i"
#' @param document the id of specific edge of interest
#'
#' @return The vector of constants representing Beta part of each interaction pattern
#'
#' @export
BetaInEqC = function(nIP, lambda.i, document) {
	# Calculate Beta part of the equation used in multinomial sampling of C
	#
	# Args 
	#  nIP total number of interaction patterns specified by the user
	#  lambda.i stochastic intensity of specific sender "i"
	#  document the id of specific document of interest
	#	
	# Returns
	#  The vector of constants representing Beta part of each interaction pattern 
	const = rep(NA, nIP)
	for (IP in 1L:nIP) {
		receiver = document[[2]]
		sender = rep(document[[1]], length(receiver))
		receiver.adjust = which(receiver > sender)
		receiver[receiver.adjust] = receiver[receiver.adjust] - 1 
		const[IP] = sum(lambda.i[[IP]][receiver]) - 
		            length(receiver)*log(sum(exp(lambda.i[[IP]])))
		}
		return(const)
}

#' @title TopicInEqC
#' @description Calculate topic part of the equation used in multinomial sampling of C
#'
#' @param nIP total number of interaction patterns specified by the user
#' @param K total number of topics specified by the user
#' @param currentZ current state of the assignment of topics 
#' @param alpha Dirichlet concentration prior for topic distribution
#' @param mvec Dirichlet base prior for topic distribution
#' @param document the id of specific edge of interest
#'
#' @return The vector of constants representing topic part of each interaction pattern
#'
#' @export
TopicInEqC = function(nIP, K, currentZ, alpha, mvec, document) {
	# Calculate topic part of the equation used in multinomial sampling of C
	#
	# Args 
	#  nIP total number of interaction patterns specified by the user
	#  K total number of topics specified by the user
	#  currentZ current state of the assignment of topics 
	#  alpha Dirichlet concentration prior for topic distribution
	#  mvec Dirichlet base prior for topic distribution
	#  document the id of specific edge of interest
	#
	# Returns
	#  The vector of constants representing topic part of each K
	const = rep(NA, nIP)
	for (IP in 1L:nIP) {
		topics = currentZ[[document]]
		if (length(topics) > 0) { 
			table.topics = tabulateC(topics, K) 
		} else {
			table.topics = rep(0, K)
		}
		num = sum(log(table.topics[topics] - 1 + alpha[IP] * mvec[topics, IP]))
		denom = length(topics) * log(length(topics) - 1 + alpha[IP])
		const[IP] = num - denom
		}
  return(const)
}

#' @title TopicInEqZ
#' @description Calculate topic-IP part of the equation used in multinomial sampling of Z
#'
#' @param K total number of topics specified by the user
#' @param currentC current state of the assignment of interaction patterns 
#' @param currentZ current state of the assignment of topics 
#' @param alpha Dirichlet concentration prior for document-topic distribution
#' @param mvec Dirichlet base prior for document-topic distribution
#' @param document the id of specific edge of interest
#'
#' @return The vector of constants representing topic part of each K
#'
#' @export
TopicInEqZ = function(K, currentC, currentZ, alpha, mvec, document) {
	# Calculate topic-IP part of the equation used in multinomial sampling of Z
	#
	# Args 
  #  K total number of topics specified by the user
  #  currentC current state of the assignment of interaction patterns 
	#  currentZ current state of the assignment of topics 
	#  alpha Dirichlet concentration prior for document-topic distribution
	#  mvec Dirichlet base prior for document-topic distribution
	#  document the id of specific edge of interest
	#
	# Returns 
	#  The vector of constants representing topic part of each K
	topics = currentZ[[document]]
	if (length(topics) > 0) { 
		table.topics = tabulateC(topics, K) 
	} else {
		table.topics = rep(0, K)
	}
	const = log(table.topics - as.numeric(table.topics > 0) + 
	        alpha[currentC[document]] * mvec[, currentC[document]])
	return(const)
}

#' @title Netstats
#' @description Calculate the network statistics (dydadic, triadic, or degree) given the history of interactions
#'
#' @param history list of document information with 3 elements (sender, receiver, time)
#' @param node nodelist containing the ID of nodes
#' @param sender specific timepoint that we are calculating the time difference from
#' @param netstat which type of network statistics to use ("dyadic", "triadic", "degree")
#'
#' @return Matrix of network statistics for specific sender and all possible receivers
#'
#' @export
Netstats = function(history, node, sender, netstat) {
  # Calculate the chosen network statistics given the history of interactions
  #
  # Args 
  #  history list of document information with 3 elements (sender, receiver, time)
  #  node nodelist containing the ID of nodes
  #  sender specific timepoint that we are calculating the time difference from
  #  netstat which type of network statistics to use ("dyadic", "triadic", "degree")
  #
  # Returns 
  #  Matrix of network statistics for specific sender and all possible receivers
  intercept = rep(1, length(node) - 1)
  netstatmat = intercept
  if ("dyadic" %in% netstat) {
    dyadic = Dyadic(history, node, sender)
    netstatmat = cbind(netstatmat, dyadic)
  }
  if ("triadic" %in% netstat) {
    triadic = Triadic(history, node, sender)
    netstatmat = cbind(netstatmat, rowSums(triadic))
  }
  if ("degree" %in% netstat) {
    degree = Degree(history, node, sender)
    netstatmat = cbind(netstatmat, degree)
  }
  netstatmat = matrix(netstatmat, nrow = length(node) - 1, byrow = FALSE)
  return(netstatmat)
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
#' @param lambda parameter of speed at which sender replies, with larger values indicating faster response time
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
#' @return MCMC output containing IP assignment, topic assignment, and beta chain (optional the plots to check the convergence)
#'
#' @export
MCMC = function(edge, node, textlist, vocabulary, nIP, K, delta.B, lambda = 0.05,
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
  #  lambda parameter of response speed with larger values indicating faster response
  #  outer size of outer iteration
  #  n1 size of first inner iteration for updates of interaction patterns C
  #  n2 size of second inner iteration for updates of topics K
  #  n3 size of third inner iteration for updates of beta B
  #  burn iterations to be discarded at the beginning of beta chain
  #  thin the thinning interval of beta chain
  #  seed an integer value which controls random number generation
  #  plot to plot the convergence diagnostics or not (TRUE/FALSE)
  #
  # Returns 
  #  MCMC output containing IP assignment, topic assignment, and beta chain
  
  set.seed(seed)

  # initialize alpha, mvec, delta, nvec, eta, lvec, and gammas
  alpha = rep(50 / K, nIP)
  mvec = matrix(1 / K, nrow = K, ncol = nIP)
  beta = 10
  W = length(vocabulary)
  nvec = rep(1, W) / W
  eta = 50 / nIP
  lvec = rep(1, nIP) / nIP
  gammas = rdirichlet(1, eta * lvec)

  # initialize C, theta and Z
  currentC = Selected(length(edge), gammas)
  theta = lapply(seq(along = edge), function(d) {
    rdirichlet(1, alpha[currentC[d]] * mvec[,currentC[d]])
    })					
  currentZ = lapply(seq(along = edge), function(d) {
    if (length(textlist[[d]]) > 0)
      Selected(length(textlist[[d]]), theta[[d]])
    })
  empty.docs = which(unlist(lapply(currentZ, function(d){ 
    length(d) 
    })) == 0)	

  # initialize beta
  P = 1 + 2 * ("dyadic" %in% netstat) + 1 * ("triadic" %in% netstat) + 2 * ("degree" %in% netstat)
  sigma = 1
  Beta.mat = list()
  for (IP in 1L:nIP) {
    Beta.mat[[IP]] = matrix(0, nrow = P, ncol = (n3 - burn) / thin)
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
    
    #update the hyperparameter alpha and mvec
    vec = AlphamvecOpt(nIP, K, currentC, currentZ, alpha, mvec)
    alpha = unlist(lapply(vec, sum))
    mvec = vapply(1L:nIP, function(IP) {
      vec[[IP]] / alpha[IP]
      }, rep(1, K))
  
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
    
    cat("inner iteration 2", "\n")
    textlist.raw = unlist(textlist)
    finalZlist.raw = unlist(currentZ)
    table.W = lapply(1L:K, function(k) {
      tabulateC(textlist.raw[which(finalZlist.raw == k)], W)
      })
    
    for (i2 in 1L:n2) {
      for (d in (seq(along = edge))[-empty.docs]) { 
        textlist.d = textlist[[d]]
        topicpart.d = TopicInEqZ(K, currentC, currentZ, alpha, mvec, d)
        wordpart.d = WordInEqZ(K, textlist.d, table.W, beta, nvec)
        for (w in 1L:nrow(wordpart.d)) {
          const.Z = topicpart.d + wordpart.d[w, ]
          zw.old = currentZ[[d]][w]
          zw.new = Selected(1, exp(const.Z))
          if (zw.new != zw.old) {
            currentZ[[d]][w] = zw.new
            topicpart.d = TopicInEqZ(K, currentC, currentZ, alpha, mvec, d)	
            wordpart.d = WordInEqZ(K, textlist.d, table.W, beta, nvec)
            table.W[[zw.old]][textlist.d[w]] = table.W[[zw.old]][textlist.d[w]] - 1
            table.W[[zw.new]][textlist.d[w]] = table.W[[zw.new]][textlist.d[w]] + 1
          }
        }
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

#' @title TableBetaIP
#' @description Generate a table summary of the MCMC chain of network statistics coefficients (beta) for each interaction pattern
#'
#' @param MCMCchain a chain obtained using MCMC function
#'
#' @return List of table containing the posterior summary of beta for each interaction pattern
#'
#' @export
TableBetaIP = function(MCMCchain) {
	# Generate a table summary of the MCMC chain of beta for each IP
	#
	# Args 
	#  MCMCchain a chain obtained using MCMC function
	#
	# Returns
	#  List of table containing the posterior summary of beta for each IP

	table.beta = lapply(MCMCchain$B, function(x) {
		summary(mcmc(t(x)))$quantiles[, c(3, 1, 5)]
 	})
 	for (IP in 1L:length(table.beta)) {
 		rownames(table.beta[[IP]]) = c("Intercept", "Send", "Receive", 
 		                               "Triangles", "Outdegree", "Indegree")
 		colnames(table.beta[[IP]]) = c("median", "lower2.5%", "upper97.5%")
 	}
 	return(table.beta)
}

#' @title PlotBetaIP
#' @description Draw a boxplot of the MCMC chain of network statistics (beta) for each interaction pattern
#'
#' @param MCMCchain a chain obtained using MCMC function
#'
#' @return Joint boxplot of the posterior summary of beta (should click for the locator)
#'
#' @export
PlotBetaIP = function(MCMCchain) {
	# Draw a boxplot of the MCMC chain of beta for each IP
	#
	# Args 
	#  MCMCchain a chain obtained using MCMC function
	#
	# Returns
	#  Joint boxplot of the posterior summary of beta (should click for the locator)
	combined.data = list()
	P = nrow(MCMCchain$B[[1]])
	nIP = length(MCMCchain$B)
	for(b in 1L:P){
	combined.data[[b]] = sapply(1L:nIP, function(c){
		cbind(MCMCchain$B[[c]][b,])
		})
		}
	forbox = melt(combined.data)
	boxplot = boxplot(forbox$value ~ forbox$X2 + forbox$L1,
	          at = c(sapply(1L:P, function(x){
	            ((nIP+1) * x - nIP)((nIP+1) * x - 1)
	            })),
	          col = gray.colors(nIP), axes = FALSE, 
	          main = "Comparison of beta coefficients for different IPs")
	abline(h = 0, lty = 1, col = "red")
	axis(2, labels = TRUE)
	box()
	axis(side = 1, line = 0.5, lwd = 0, 
	     at = c(sapply(1L:P, function(x){
	       median(((nIP + 1) * x - nIP)((nIP + 1) * x - 1))
	       })), 
	     labels = c("Intercept", "Send", "Receive", 
	                "Triangles", "Outdegree", "Indegree"))
	legend(locator(1), c(paste("IP", 1L:nIP)), col = gray.colors(nIP), pch = 15)
}

#' @title PlotTopicIP
#' @description Draw a barplot of the topic distributions for each interaction pattern
#'
#' @param MCMCchain MCMCchain a chain obtained using MCMC function
#' @param K total number of topics specified by the user
#'
#' @return Joint barplot of the topic distribution (should click for the locator)
#'
#' @export
PlotTopicIP = function(MCMCchain, K) {
	# Draw a barplot of the topic distributions for each interaction pattern
	#
	# Args 
	#  MCMCchain a chain obtained using MCMC function
	#  K total number of topics specified by the user
	#
	# Returns
	#  Joint barplot of the topic distribution (should click for the locator)
	Zsummary = list()
	nIP = length(MCMCchain$B)
	for (IP in 1L:nIP) {
		Zsummary[[IP]] = list()
		iter = 1
		for (d in which(MCMCchain$C == IP)) {
			Zsummary[[IP]][[iter]] = MCMCchain$Z[[d]]
			iter = iter + 1
			}
		}
	topic.dist = t(sapply(Zsummary, function(x) {
				tabulateC(unlist(x), K) / length(unlist(x))
				}))
	colnames(topic.dist) = c(1L:K)
	barplot(topic.dist, beside = TRUE, xlab = "Topic", ylab = "Proportion", 
	        main = "Topic Distritubitions given IPs")	
	legend(locator(1), c(paste("IP", 1L:nIP)), col = gray.colors(nIP), pch = 15)
	}

#' @title PlotTopic
#' @description Draw a barplot of the topic distributions without considering interaction patterns
#'
#' @param MCMCchain MCMCchain a chain obtained using MCMC function
#' @param K total number of topics specified by the user
#'
#' @return Barplot of the topic distribution
#'
#' @export	
PlotTopic = function(MCMCchain, K) {
	# Draw a barplot of the topic distributions without considering IPs
	#
	# Args 
	#  MCMCchain a chain obtained using MCMC function
	#  K total number of topics specified by the user
	#
	# Returns
	#  Barplot of the topic distribution
	Zsummary = list()
		for (d in seq(along = MCMCchain$C)) {
			Zsummary[[d]] = MCMCchain$Z[[d]]
			}
	topic.dist = t(tabulateC(unlist(Zsummary), K) / length(unlist(Zsummary)))
	colnames(topic.dist) = c(1L:K)
	barplot(topic.dist, beside = TRUE, xlab = "Topic", ylab = "Proportion", 
	        main = "Topic Distritubitions without IPs")
	}

#' @title TableWordIP
#' @description Generate a table that summarizes token-topic assignments with high probabilities for each interaaction pattern
#'
#' @param MCMCchain MCMCchain a chain obtained using MCMC function
#' @param K total number of topics specified by the user
#' @param textlist list of text (length=number of documents in total) containing the words in each document
#' @param vocabulary all vocabularies used over the corpus
#'
#' @return List of table that summarize token-topic assignments with high probabilities (topic proportion > 0.1) for each IP
#'
#' @export
TableWordIP = function(MCMCchain, K, textlist, vocabulary) {
	# Generate a table of token-topic assignments with high probabilities for each IP
	#
	# Args 
	#  MCMCchain a chain obtained using MCMC function
	#  K total number of topics specified by the user
  #  textlist list of text containing the words in each document
  #  vocabulary all vocabularies used over the corpus
	#
	# Returns
	#  List of table that summarize token-topic assignments for each IP
	W = length(vocabulary)
	nIP = length(MCMCchain$B)
	table.word = list()
	for (IP in 1L:nIP) {
		Zsummary = list()
		topic.word = matrix(0, nrow = K, ncol = W)
		colnames(topic.word) = vocabulary
		iter = 1
		for (d in which(MCMCchain$C==IP)) {
			if (length(MCMCchain$Z[[d]])>0){
			Zsummary[[iter]] = MCMCchain$Z[[d]]
			names(Zsummary[[iter]])<- vocabulary[textlist[[d]]]
			iter = iter+1
			}
			}
		topic.dist = t(tabulateC(unlist(Zsummary), K)/length(unlist(Zsummary)))
		colnames(topic.dist) = c(1L:K)
		top.topic = which(topic.dist[, order(topic.dist, decreasing = TRUE)] > 0.1)
		all.word = unlist(Zsummary)
		for (i in seq(along = all.word)){
		matchWZ = which(colnames(topic.word) == names(all.word[i]))
		topic.word[all.word[i], matchWZ] = topic.word[all.word[i], matchWZ] + 1
		}
		table.word[[IP]] = top.topic.words(topic.word, num.words = 10)[, top.topic]
		colnames(table.word[[IP]]) = names(top.topic)
		}
	return(table.word)
}

#' @title GenerateDocs
#' @description Generate the documents according to the generative process
#'
#' @param nDocs number of documents to be generated
#' @param betas network statistics coefficients (beta) for each interaction pattern
#' @param gammas distribution of interaction patterns
#' @param node nodelist containing the ID of nodes (ID starting from 1)
#' @param vocabulary all vocabularies used over the corpus
#' @param alpha Dirichlet concentration prior for topic distribution
#' @param mvec Dirichlet base prior for topic distribution
#' @param beta Dirichlet concentration prior for topic-word distribution
#' @param nvec Dirichlet base prior for topic-word distribution
#' @param lambda parameter of speed at which sender replies, with larger values indicating faster response time
#' @param delta tuning parameter
#' @param netstat which type of network statistics to use ("dyadic", "triadic", "degree")
#' @param seed an integer value which controls random number generation
#'
#' @return edge and text generated according to the generative process
#'
#' @export
GenerateDocs = function(nDocs, betas, gammas, node, vocabulary, alpha, mvec, beta, nvec, lambda, delta = 1,
						netstat = c("dyadic", "triadic", "degree"), seed = 1) {
	
	P = 1 + 2 * ("dyadic" %in% netstat) + 1 * ("triadic" %in% netstat) + 2 * ("degree" %in% netstat)
	nIP = length(betas)
	currentIP = Selected(nDocs, gammas)
	phi = lapply(1:nrow(mvec), function(k) {
		rdirichlet(1, beta * nvec)
	})
	t.d = 0
	edge = list()
	text = list()	
	
	options(warn = -1)
	for (d in 1:nDocs) {
		edge.byC = lapply(1:nIP, function(IP) {
			if (d > 1) edge[which(currentIP[1:(d-1)] == IP)]
		})
		if (d == 1 | length(edge.byC[[currentIP[d]]]) == 0) {
			history.t.IP = matrix(0, nrow = length(node), ncol = length(node))
		} else {
			history.t.IP = Timediff(edge.byC[[currentIP[d]]], node, t.d, lambda)
			}
		X = lapply(node, function(i) {
				Netstats(history.t.IP, node, i, netstat)
				})
		XB = MultiplyXBList(X, betas[[currentIP[d]]])
		iJi = lapply(1:nrow(XB), function(i) {
			sapply(1:ncol(XB), function(j) {
				rbinom(1, 1, 1-exp(-delta*exp(XB[i, j])))	
				})
			})
		LambdaiJi = sapply(1:length(iJi), function(i) {
			exp(sum(iJi[[i]]*XB[i,]))
		})
		LambdaiJi[LambdaiJi == 1] = NA
		Time.inc = sapply(LambdaiJi, function(i) {
			rexp(1, i)
			})
		i.d = which(Time.inc == min(Time.inc[!is.na(Time.inc)]))
		j.d = which(iJi[[i.d]] == 1)
		j.d[j.d >= i.d] = j.d[j.d >= i.d] + 1
		t.d = t.d + Time.inc[i.d]
		edge[[d]] = list(sender = i.d, receiver = j.d, timestamp = t.d)
		
		N.d = rpois(1, 10)
		theta.d = rdirichlet(1, alpha[currentIP[d]] * mvec[,currentIP[d]])	
		text[[d]] = rep(NA, N.d)
		for (n in 1:N.d){
			topic.n = Selected(1, theta.d)
			text[[d]][n] = vocabulary[Selected(1, phi[[topic.n]])]
		}		
	}
	options(warn = 0)
	return(list(edge = edge, text = text))	
}

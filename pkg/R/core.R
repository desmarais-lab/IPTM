#' @useDynLib IPTM
#' @importFrom Rcpp sourceCpp 
NULL

#' @title historymat
#' @description calculate the time difference from previous interactions to specific timepoint
#'
#' @param edge edgelist in the form of matrix with 3 columns (col1:sender, col2:receiver, col3=time in unix.time format)
#' @param node nodelist containing the ID of nodes (ID's starting from 1)
#' @param when specific timepoint that we are calculating the time difference from
#'
#' @return list of time differences from previous interactions to specific timepoint (i.e. 'when') for every combination of nodes
#'
#' @export
historymat = function(edge, node, when) {
  histlist = lapply(node, function(v) {
    lapply(node, function(u) {
      numeric(0)
    })
  })
  edge2 = matrix(edge[edge[, 3] < when, ], byrow = FALSE, ncol = 3)
  if (nrow(edge2) > 0) {
    for (d in 1:nrow(edge2)) {
      histlist[[edge2[d,1]]][[edge2[d,2]]] = c(histlist[[edge2[d, 1]]][[edge2[d, 2]]],
        when - (edge2[d, 3]))
    }
  }
  return(histlist)
}

#' @title netstats
#' @description calculate six network statistiscs, intercept, send, receive, triangle, outdegree and indegree
#'
#' @param histlist list of time differences from previous interactions to specific timepoint (i.e. 'when') for every combination of nodes
#' @param node nodelist containing the ID of nodes (ID's starting from 1)
#' @param sender sender of the specific document
#' @param lambda parameter of speed at which sender replies, with larger values indicating faster response time
#'
#' @return X matrix containing six history-based network statistics 
#'
#' @export
netstats = function(histlist, node, sender, lambda) {
	allxmat  = list()
	for (s in node){
		allxmat[[s]] = list()
		for (r in node){
			allxmat[[s]][[r]] = sum(exp(-lambda * histlist[[s]][[r]]))
		}
	}
	netstatmat = matrix(NA, nrow = length(node) - 1, ncol = 6)
	it = 1
	for (receiver in node[!node == sender]) {
	send = allxmat[[sender]][[receiver]]
	receive = allxmat[[receiver]][[sender]]
	twosend = sum(vapply(node[!node == c(sender, receiver)], function(h) {
				allxmat[[sender]][[h]] * allxmat[[h]][[receiver]]
				}, c(1)))
	tworeceive = sum(vapply(node[!node == c(sender, receiver)], function(h) {
				allxmat[[h]][[sender]] * allxmat[[receiver]][[h]]
				}, c(1)))
	sibling = sum(vapply(node[!node == c(sender, receiver)], function(h) {
				allxmat[[h]][[sender]] * allxmat[[h]][[receiver]]
				}, c(1)))
	cosibling = sum(vapply(node[!node == c(sender, receiver)], function(h) {
				allxmat[[sender]][[h]] * allxmat[[receiver]][[h]]
				}, c(1)))  
	triangle = sum(twosend, tworeceive, sibling, cosibling)
	outdegree = sum(unlist(allxmat[[sender]]))
	indegree = sum(vapply(node, function(i){allxmat[[i]][[receiver]]}, c(1)))	
	netstatmat[it, ] = c(1, send, receive, triangle, outdegree, indegree)
	it = it + 1	
	}
	return(netstatmat)
}

#' @title parupdate
#' @description parameter optimization of alpha and mvec at the same time 
#'
#' @param nIP number of interaction pattern specified by the user
#' @param K number of topics specified by the user
#' @param currentC current state of the assignment of interaction patterns
#' @param currentZ current state of the assignment of topics 
#' @param alpha Dirichlet concentration prior for topic distribution
#' @param mvec Dirichlet base prior for topic distribution
#'
#' @return optimized value of alpha * mvec
#'
#' @export
parupdate =  function(nIP, K, currentC, currentZ, alpha, mvec) {
	finalvec = list()
	corpusCnew = sortedZ(nIP, currentC, currentZ)

	for (IP in 1:nIP) {
		iter=1
		vec = alpha[[IP]] * mvec[[IP]]
		nwords = mapply(length, corpusCnew[[IP]])
		ctable = tabulate(nwords)
		cklist = matrix(NA, nrow = length(corpusCnew[[IP]]), ncol = K)
		for (d in 1:length(corpusCnew[[IP]])) {
			cklist[d,] = tabulate(corpusCnew[[IP]][[d]], nbins = K)
			}
	cktable = lapply(1:K, function(k){
			  tabulate(cklist[,k])
			  })

	while ((abs(alpha[[IP]] - sum(vec)) > 0.001) | (iter == 1)) {
	alpha[[IP]] = sum(vec)
	S = Supdate(alpha[[IP]], ctable)		
	s = Skupdate(vec, cktable)	
	vec = vec * s/S
	iter = iter + 1
	}
	finalvec[[IP]] = vec
	}	
	return(finalvec)	
}	

#' @title selected
#' @description find out the chosen category from the multinomial distribution
#'
#' @param samples number of samples to draw
#' @param proportions numeric non-negative vector of length K, specifying the probability for the K classes
#'
#' @return the chosen category from the rmultinom
#'
#' @export
selected = function(samples, proportions) {
  chosen  = rmultinom(samples, 1, proportions)
  out = vapply(1:samples, function(s) {which(chosen[, s] == TRUE)}, c(1))
  return(out)
}

#' @title betapartC
#' @description calculate the log of beta part used to obtain the constants for multinomial sampling of IP
#'
#' @param nIP number of interaction patterns specified by the user
#' @param lambdai stochastic intensity of specific sender 'i'
#' @param specificedge a row of specific edge of interest (i.e. vector of sender, reciever, and time)
#'
#' @return the log of beta part (first term of Eq. 11)
#'
#' @export
betapartC = function(nIP, lambdai, specificedge) {
  const = rep(NA, nIP)
  for (IP in 1:nIP) {
    edge = matrix(specificedge, ncol = 3)
    receiver = ifelse(edge[, 1] > edge[, 2], edge[, 2],  edge[, 2] - 1)
    const[IP] = lambdai[[IP]][receiver] - log(sum(exp(lambdai[[IP]])))
  }
  return(const)
}

#' @title topicpartC
#' @description calculate the log of topic-IP part used to obtain the constants for multinomial sampling of IP
#'
#' @param nIP number of interaction patterns specified by the user
#' @param K number of topics specified by the user
#' @param currentZ current state of the assignment of topics 
#' @param alpha Dirichlet concentration prior for topic distribution
#' @param mvec Dirichlet base prior for topic distribution
#' @param document specific document of interest 
#'
#' @return the log of topic-IP part (last term of Eq. 11)
#'
#' @export
topicpartC = function(nIP, K, currentZ, alpha, mvec, document) {
  const = rep(NA, nIP)
	for (IP in 1:nIP) {
		topics = currentZ[[document]]
		tabletopics = tabulate(topics, nbins = K)
		num = sum(log(tabletopics[topics] - 1 + alpha[[IP]] * mvec[[IP]][topics]))
		denom = length(topics) * log(length(topics) - 1 + alpha[[IP]])
		const[IP] = num - denom
		}
  return(const)
}

#' @title topicpartZ
#' @description calculate the log of topic-IP part used to obtain the constants for multinomial sampling of K
#'
#' @param currentC current state of the assignment of interaction patterns 
#' @param K number of topics specified by the user
#' @param currentZ current state of the assignment of topics 
#' @param alpha Dirichlet concentration prior for document-topic distribution
#' @param mvec Dirichlet base prior for document-topic distribution
#' @param document specific document of interest 
#'
#' @return the log of topic-IP part (first term of Eq. 13)
#'
#' @export
topicpartZ = function(currentC, K, currentZ, alpha, mvec, document) {
	tabletopics = tabulate(currentZ[[document]], nbins = K)
	const = log(tabletopics - ifelse(tabletopics > 0, 1, 0) + alpha[[currentC[document]]] * mvec[[currentC[document]]])
	return(const)
}

#' @title wordpartZ
#' @description: calculate the log of word-topic part used to obtain the constants for multinomial sampling of K
#'
#' @param K number of topics specified by the user
#' @param textlistd list of text containing the words in specific document d
#' @param tableW summary table of topic-word assignments
#' @param delta Dirichlet concentration prior for topic-word distribution
#' @param nvec Dirichlet base prior for topic-word distribution
#'
#' @return the log of word-topic part (last term of Eq. 13)
#'
#' @export
wordpartZ = function(K, textlistd, tableW, delta, nvec) {
  const = matrix(NA, nrow=length(textlistd), ncol = K)
  for (k in 1:K) {
    num  = log(tableW[[k]][textlistd] -
                 ifelse(tableW[[k]][textlistd] > 0, 1, 0) + delta * nvec[textlistd])
    denom  = log(sum(tableW[[k]]) - sum(ifelse(tableW[[k]] > 0, 1, 0)) + delta)
    const[, k] = num - denom
  }
  return(const)
}

#' @title betapartB
#' @description calculate the log of beta part used to obtain the constants for multinomial sampling of Beta
#'
#' @param nIP number of interaction patterns specified by the user
#' @param lambdai stochastic intensity of specific sender 'i'
#' @param edgeC separated edgelist according to the interaction patterns 
#'
#' @return the log of beta part (used in Eq. 16)
#'
#' @export
betapartB = function(nIP, lambdai, edgeC) {
  const = rep(NA, nIP)
  for (IP in 1:nIP) {
    edge = matrix(edgeC[[IP]], ncol = 3)
    receiver = ifelse(edge[, 1] > edge[, 2], edge[, 2],  edge[, 2] - 1)
    const[IP] = sum(vapply(1:length(receiver), function(r) {
      lambdai[[IP]][r, receiver[r]] - log(sum(exp(lambdai[[IP]][r, ])))
    }, c(1)))
  }
  return(const)
}

#' @title logWZ
#' @description calculate the log of unnormalized constant (corresponding to product of Eq.13 in Section 3.2) to check the convergence
#'
#' @param nIP number of interaction patterns specified by the user
#' @param K number of topics specified by the user
#' @param currentC current state of the assignment of interaction patterns 
#' @param currentZ current state of the assignment of topics 
#' @param textlist list of text (length=number of documents in total) containing the words in each document
#' @param tableW summary table of topic-word assignments
#' @param alpha Dirichlet concentration prior for document-topic distribution
#' @param mvec Dirichlet base prior for document-topic distribution
#' @param delta Dirichlet concentration prior for topic-word distribution
#' @param nvec Dirichlet base prior for topic-word distribution
#'
#' @return the log of unnormalized constant corresponding to product of Eq.13 in Section 3.2
#'
#' @export
logWZ = function(nIP, K, currentC, currentZ, textlist, tableW, alpha, mvec, delta, nvec) {
	finalsum = 0
	for (d in 1:length(currentC)) {
		ktable = tabulate(currentZ[[d]], nbins = K)
		textlistd = textlist[[d]]
		it = 1
		for (k in currentZ[[d]]) {
			part1 = log(tableW[[k]][textlistd[it]] -1 + delta * nvec[textlistd[it]])
			part2 = log(sum(tableW[[k]]) - sum(tableW[[k]]>0) + delta)
			part3 = log(ktable[k] -1 + alpha[[currentC[d]]] * mvec[[currentC[d]]][k])
			part4 = log(sum(ktable) -1 + alpha[[currentC[d]]])
			finalsum = finalsum + part1 - part2 + part3 - part4	
			it = it + 1
	}
	}
	return(finalsum)
}

#' @title MCMC
#' @description
#'
#' @param edge edgelist in the form of matrix with 3 columns (col1:sender, col2:receiver, col3=time in unix.time format)
#' @param node nodelist containing the ID of nodes (ID's starting from 1)
#' @param textlist list of text (length=number of documents in total) containing the words in each document
#' @param vocabulary all vocabularies used over the corpus
#' @param nIP number of interaction patterns specified by the user
#' @param K number of topics pattern specified by the user
#' @param delta_B proposal distribution variance parameter for beta 
#' @param outer size of outer iteration 
#' @param n1 size of first inner iteration for updates of interaction patterns C
#' @param n2 size of second inner iteration for updates of topics K
#' @param n3 size of third inner iteration for updates of betas B
#' @param burn iterations to be discarded at the beginning of the chain
#' @param thin the thinning interval 
#' @param seed an integer value which controls random number generation
#'
#' @return MCMC output containing IP assignment, topic assignment, beta, the log of unnormalized constant, and optimized hyperparameter alpha
#'
#' @export
MCMC = function(edge, node, textlist, vocabulary, nIP, K, delta_B, outer = 200,
  n1 = 3, n2 = 3, n3 = 3300, burn = 300, thin = 3, seed = 1) {

  set.seed(seed)

  # initialize alpha, mvec, delta, nvec, eta, lvec, and gammas
  alpha =  lapply(1:nIP, function(x){50 / K})
  mvec = lapply(1:nIP, function(x){rep(1, K) / K}) 
  delta = 0.01
  W = length(vocabulary)
  nvec = rep(1, W) / W
  eta = 50 / nIP
  lvec = rep(1, nIP) / nIP
  gammas = rdirichlet(1, eta * lvec)

  # initialize C, theta and Z
  currentC = selected(nrow(edge), gammas)
  theta = lapply(1:nrow(edge), function(d) {
  	      rdirichlet(1, alpha[[currentC[d]]] * mvec[[currentC[d]]])
  	      })					
  currentZ = lapply(1:nrow(edge), function(d) {
			 selected(length(textlist[[d]]), theta[[d]])
			 })	  

  # initialize beta
  P = 6
  sigma = 1
  bmat = list()
  for (IP in 1:nIP) {
    bmat[[IP]] = matrix(NA, nrow = P, ncol = (n3 - burn) / thin)
    bmat[[IP]][,] = rmvnorm(1, rep(0, P), sigma^2 * diag(P))
  }

  #to check the convergence  
  logWZmat <- c()							  
  alphamat <- alpha


  #start outer iteration
  for(o in 1:outer) {
    cat("outer iteration = ", o, "\n")
  
  #update the hyperparameter alpha and mvec
   vec = parupdate(nIP, K, currentC, currentZ, alpha, mvec)
   alpha = lapply(1:nIP, function(IP) {sum(vec[[IP]])})
   alphamat = lapply(1:nIP, function(IP) {c(alphamat[[IP]], alpha[[IP]])})
   mvec = lapply(1:nIP, function(IP) {vec[[IP]] / alpha[[IP]]})	

    cat("inner iteration 1", "\n")
    lambdai = list()
    corpusC  = sortedZ(nIP, currentC, currentZ)
   
    for (i1 in 1:n1) {
      # C update given Z and B - within each document d
      for (d in 1:nrow(edge)) {
        currentC2 = currentC[-d]
        edgeC = lapply(1:nIP, function(IP) {
          edge[which(currentC2 == IP), ]
        })
        for (IP in 1:nIP) {
          histlist = historymat(edgeC[[IP]], node, edge[d, 3])
          allxmatlist = netstats(histlist, node, edge[d, 1], 0.05)
          lambdai[[IP]] = multiplyXB(allxmatlist, bmat[[IP]][, (n3 - burn) / thin])
        }
        const = log(gammas)+betapartC(nIP, lambdai, edge[d, ]) +
          topicpartC(nIP, K, currentZ, alpha, mvec, d)
        currentC[d] = selected(1, exp(const))
      }
    }

    cat("inner iteration 2", "\n")
    textlist2 = unlist(textlist)
    finalZlist2 = unlist(currentZ)
    tableW = lapply(1:K, function(k) {
      tabulate(textlist2[which(finalZlist2 == k)], nbins = W)
    })

    for (i2 in 1:n2) {
      for (d in 1:nrow(edge)) {
        textlistd = textlist[[d]]
        topicpartd = topicpartZ(currentC, K, currentZ, alpha, mvec, d)
        wordpartd = wordpartZ(K, textlistd, tableW, delta, nvec)
        for (w in 1:nrow(wordpartd)) {
          const2 = topicpartd + wordpartd[w, ]
          zwold = currentZ[[d]][w]
          zwnew = selected(1, exp(const2))
          if (zwnew != zwold) {
          	currentZ[[d]][w] = zwnew
          	topicpartd = topicpartZ(currentC, K, currentZ, alpha, mvec, d)	
            tableW[[zwold]][textlistd[w]] = tableW[[zwold]][textlistd[w]] - 1
            tableW[[zwnew]][textlistd[w]] = tableW[[zwnew]][textlistd[w]] + 1
            topicpartd = topicpartZ(currentC, K, currentZ, alpha, mvec, d)
            wordpartd = wordpartZ(K, textlistd, tableW, delta, nvec)
          }
        }
      }
    }

    logWZmat = c(logWZmat, logWZ(nIP, K, currentC, currentZ, textlist, tableW,
      alpha, mvec, delta, nvec))

    cat("inner iteration 3", "\n")
    bold = list()
    allxmatlist2 = list()
    for (IP in 1:nIP) {
      bold[[IP]] = bmat[[IP]][, (n3 - burn) / thin]
      edgeC[[IP]] = edge[which(currentC == IP), ]
      allxmatlist2[[IP]] = list()

      for (d in 1:nrow(edgeC[[IP]])) {
        histlist2 = historymat(edgeC[[IP]], node, edgeC[[IP]][d, 3])
        allxmatlist2[[IP]][[d]] = netstats(histlist2, node, edgeC[[IP]][d, 1], 0.05)
      }
    }

    lambdaiold = list()
    lambdainew = list()
    bnew = list()

    for (i3 in 1:n3) {
      for (IP in 1:nIP) {
        bnew[[IP]] = rmvnorm(1, bold[[IP]], (delta_B)^2 * diag(P))
        lambdaiold[[IP]] = t(vapply(1:nrow(edgeC[[IP]]), function(d) {
          multiplyXB(allxmatlist2[[IP]][[d]], bold[[IP]])
        }, rep(1, length(node) - 1)))
        lambdainew[[IP]] = t(vapply(1:nrow(edgeC[[IP]]), function(d) {
          multiplyXB(allxmatlist2[[IP]][[d]], bnew[[IP]])
        }, rep(1, length(node) - 1)))
      }

      prior = vapply(1:nIP, function(IP) {
        dmvnorm(bnew[[IP]], rep(0, P), sigma^2 * diag(P), log = TRUE) -
          dmvnorm(bold[[IP]], rep(0, P), sigma^2 * diag(P), log = TRUE)
      }, c(1))
      post = betapartB(nIP, lambdainew, edgeC) - betapartB(nIP, lambdaiold, edgeC)
      loglikediff = prior + post

      u = log(runif(nIP, 0, 1))
      for (IP in which((u < loglikediff) == TRUE)) {
        bold[[IP]]  = bnew[[IP]]
      }

      if (i3 > burn && i3 %% (thin) == 0) {
        for (IP in 1:nIP) {
          bmat[[IP]][ ,(i3 - burn) / thin] = bold[[IP]]
        }
      }
    }
  }

    out = list(C = currentC, Z = currentZ, B = bmat, L = logWZmat, A = alphamat)

  return(out)
}


#' @title table_betaIP
#' @description summarize the MCMC chain of network statistics (beta) for each interaction pattern, by using a table
#'
#' @param MCMCchain a chain obtained using MCMC function
#'
#' @return list of table containing the posterior summaries of beta for each interaction pattern
#'
#' @export
table_betaIP = function(MCMCchain) {
	require(MCMCpack)
	tables = lapply(MCMCchain$B, function(x) {
		summary(mcmc(t(x)))$quantiles[,c(3,1,5)]
 	})
 	for (IP in 1:length(tables)) {
 		rownames(tables[[IP]]) = c('intercept','send','receive','triangles', 'outdegree', 'indegree')
 		colnames(tables[[IP]]) = c('median', 'lower2.5%', 'upper97.5%')
 	}
 	return(tables)
}

#' @title plot_betaIP
#' @description summarize the MCMC chain of network statistics (beta) for each interaction pattern, by drawing a joint boxplot
#'
#' @param MCMCchain MCMCchain a chain obtained using MCMC function
#'
#' @return plot of the posterior summaries of beta for each interaction pattern (should click for the locator)
#'
#' @export
plot_betaIP = function(MCMCchain) {
	require(reshape)
	data = list()
	P = nrow(MCMCchain$B[[1]])
	nIP = length(MCMCchain$B)
	for(b in 1:P){
	data[[b]] = sapply(1:nIP, function(c){
		cbind(MCMCchain$B[[c]][b,])
		})
		}
	forbox = melt(data)
	boxplot = boxplot(forbox$value ~ forbox$X2 + forbox$L1,
	at = c(sapply(1:P, function(x){((nIP+1) * x - nIP):((nIP+1) * x - 1)})),
	col = gray.colors(nIP), axes = FALSE, main ='Comparison of beta coefficients for different IPs')
	abline(h = 0, lty = 1, col = 'red')
	axis(2, labels = TRUE)
	box()
	axis(side = 1, at = c(sapply(1:P, function(x){median(((nIP + 1) * x - nIP):((nIP + 1) * x - 1))})), line = 0.5, 
	lwd = 0, labels = c('intercept','send','receive','triangles', 'outdegree', 'indegree'))
	legend(locator(1), c(paste('IP', 1:nIP)), col = gray.colors(nIP), pch = 15)
}

#' @title plot_topicIP
#' @description plot the topic distributions for each interaction pattern
#'
#' @param MCMCchain MCMCchain a chain obtained using MCMC function
#' @param K number of topics pattern specified by the user
#'
#' @return joint barplot of the topic distribution (should click for the locator)
#'
#' @export
plot_topicIP = function(MCMCchain, K) {
	Zsummary = list()
	nIP = length(MCMCchain$B)
	for (IP in 1:nIP) {
		Zsummary[[IP]] = list()
		iter = 1
		for (d in which(MCMCchain$C == IP)) {
			Zsummary[[IP]][[iter]] = MCMCchain$Z[[d]]; 
			iter = iter + 1
			}
		}
	topicdist = t(sapply(Zsummary, function(x) {
				tabulate(unlist(x), nbins = K) / length(unlist(x))
				}))
	colnames(topicdist) = c(1:K)
	barplot(topicdist, beside = TRUE, xlab = "topic", ylab = "proportion", main ='Topic Distritubitions given IPs')
	
	legend(locator(1), c(paste('IP', 1:nIP)), col = gray.colors(nIP), pch = 15)
	}

#' @title plot_topic
#' @description  plot the topic distributions without considering interaction patterns
#'
#' @param MCMCchain MCMCchain a chain obtained using MCMC function
#' @param K number of topics pattern specified by the user
#'
#' @return barplot of the topic distribution
#'
#' @export	
plot_topic = function(MCMCchain, K) {
	Zsummary = list()
		for (d in 1:length(MCMCchain$C)) {
			Zsummary[[d]] = MCMCchain$Z[[d]]
			}
	topicdist = t(tabulate(unlist(Zsummary), nbins = K) / length(unlist(Zsummary)))
	colnames(topicdist) = c(1:K)
	barplot(topicdist, beside = TRUE, xlab = "topic", ylab = "proportion", main = 'Topic Distritubitions without IPs')
	}

#' @title table_wordIP
#' @description generate the table that summarizes token-topic assignments of highest probabilities for each interaaction pattern
#'
#' @param MCMCchain MCMCchain a chain obtained using MCMC function
#' @param K number of topics pattern specified by the user
#' @param textlist list of text (length=number of documents in total) containing the words in each document
#' @param vocabulary all vocabularies used over the corpus
#'
#' @return list of tables that summarize token-topic assignments of highest probabilities (topic proportion > 0.1) for each IP
#'
#' @export
table_wordIP = function(MCMCchain, K, textlist, vocabulary) {
	require(lda)
	W = length(vocabulary)
	nIP = length(MCMCchain$B)
	table_word = list()
	for (IP in 1:nIP) {
		Zsummary = list()
		topicword = matrix(0, nrow = K, ncol = W)
		colnames(topicword) = vocabulary
		iter = 1
		for (d in which(MCMCchain$C==IP)) {
			Zsummary[[iter]] = MCMCchain$Z[[d]]
			names(Zsummary[[iter]])<- vocabulary[text[[d]]]
			iter = iter+1
			}
		topicdist = t(tabulate(unlist(Zsummary), nbins=K)/length(unlist(Zsummary)))
		colnames(topicdist)<-c(1:K)
		toptopic = which(topicdist[,order(topicdist, decreasing=TRUE)]>0.1)
		allwords = unlist(Zsummary)
		for (i in 1:length(allwords)){
		matchWZ = which(colnames(topicword)==names(allwords[i]))
		topicword[allwords[i], matchWZ] = topicword[allwords[i], matchWZ]+1
		}
		table_word[[IP]] = top.topic.words(topicword, num.words=10)[, toptopic]
		colnames(table_word[[IP]]) = names(toptopic)
		}
	return(table_word)
}

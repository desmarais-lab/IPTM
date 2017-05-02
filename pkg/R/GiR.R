library(mvtnorm)
library(MCMCpack)
library(entropy)
library(Rcpp)
library(RcppArmadillo)
sourceCpp('~/Desktop/IPTM-master/pkg/src/sampler.cpp')
source('~/Desktop/IPTM-master/pkg/R/core.R')
GenerateDocs = function(nDocs, node, vocabulary, nIP, K, nwords, alpha, mvec, betas, nvec,
						b, delta, currentC, netstat, base.edge, base.text, seed,
						topic_token_assignments = NULL, topic_token_counts = NULL, word_type_topic_counts = NULL, 
						forward = FALSE, backward_init = FALSE, backward = FALSE, base = FALSE) {
	set.seed(seed)
 	W = length(vocabulary)
  	phi = lapply(1L:K, function(k) {
		rdirichlet_cpp(1, betas * nvec)
	})
	L = 3
    P = 1 + L * (2 * ("dyadic" %in% netstat) + 4 * ("triadic" %in% netstat) + 2 *("degree" %in% netstat))
    
    if (base) {
    	t.d = 0
    } else {
    	t.d = 384
    }
	edge = base.edge
	text = base.text
	p.d = matrix(NA, nrow = nDocs, ncol = nIP)
	
	options(warn = -1)
	for (d in 1:nDocs) {
		N.d = nwords
		text[[length(base.text) + d]] = rep(NA, N.d)
		
		if (!backward) {
			theta.d = rdirichlet_cpp(1, alpha * mvec)	
			topic.d = multinom_vec(max(1, N.d), theta.d)
		
			if (N.d > 0) {
				for (n in 1:N.d){
					text[[length(base.text) + d]][n] = multinom_vec(1, phi[[topic.d[n]]])
				}
				names(text[[length(base.text) + d]]) = topic.d
			}
		} else {
			phi.k = rep(NA, K)
			topic.d = topic_token_assignments[[d]]
			word.d = as.numeric(names(topic_token_assignments[[d]]))
			if (N.d > 0) {
				for (n in 1:N.d){
					word_type_topic_counts[word.d [n], topic.d[n]] = word_type_topic_counts[word.d [n], topic.d[n]] - 1
					for (w in 1:W) {
						phi.k[w] = (word_type_topic_counts[w, topic.d[n]] + betas * nvec[w])/(topic_token_counts[topic.d[n]] + betas)
					} 
					text[[length(base.text) + d]][n] = multinom_vec(1, phi.k)
					word_type_topic_counts[text[[d]][n], topic.d[n]] = word_type_topic_counts[text[[d]][n], topic.d[n]] + 1
				}
				names(text[[length(base.text) + d]]) = topic.d
			}
		}

		p.d[d, ] = vapply(1L:nIP, function(IP) {
			sum(topic.d %in% which(currentC == IP))
			}, c(1)) / max(1, N.d)
		if (base & t.d < 384) {
			history.t = lapply(1:nIP, function(IP) {
				lapply(1:3, function(l){
				matrix(0, length(node), length(node))
				})
			})
		} else {	
		history.t = History(edge, p.d, node, t.d + 10^(-10))
		}
		X = lapply(node, function(i) {
  	        Netstats(history.t, node, i, netstat)
            })
   	 	XB = MultiplyXBList(X, b)     
		lambda = lambda_cpp(p.d[d,], XB)
		
		iJi = rbinom_mat((delta * lambda) / (delta * lambda + 1))
		while (sum(iJi) == 0) {
			iJi = rbinom_mat((delta * lambda) / (delta * lambda + 1))
		}
		LambdaiJi = lambdaiJi(p.d[d,], XB, iJi)
		Time.inc = vapply(LambdaiJi, function(lambda) {
			rexp(1, lambda)
			}, c(1))
		i.d = which(Time.inc == min(Time.inc[!is.na(Time.inc)]))
		j.d = which(iJi[i.d,] == 1)
		t.d = t.d + Time.inc[i.d]
		edge[[length(base.edge) + d]] = list(sender = i.d, receiver = j.d, timestamp = t.d)		
		}
	options(warn = 0)
	if (forward) {
		edge = edge[-(1:length(base.edge))]
		text = text[-(1:length(base.text))]
	}
	if (base == TRUE & t.d > 384) {
		cutoff = which_int(384, vapply(1:length(edge), function(d) {edge[[d]][[3]]}, c(1))) - 1
		edge = edge[1:cutoff]
		text = text[1:cutoff]
	}
	return(list(edge = edge, text = text, b = b, base = length(base.edge)))							
} 



GiR_stats = function(GiR_sample, K, currentC, vocabulary, forward = FALSE, backward = FALSE) {
	edge = GiR_sample$edge
	text = GiR_sample$text
	if (backward) {
		edge = edge[-(1:GiR_sample$base)]
		text = text[-(1:GiR_sample$base)]
	}
	
	GiR_stats = c()
	nDocs = length(edge)
	P = length(GiR_sample$b[[1]])
	nIP = length(GiR_sample$b)
	nwords = length(GiR_sample$text[[1]])
	W = length(vocabulary)
	
	GiR_stats[1:P] = Reduce('+', GiR_sample$b) / nIP
	GiR_stats[P + 1] = mean(vapply(1:nDocs, function(d) {
						length(edge[[d]][[2]])
						}, c(1)))
	GiR_stats[P + 2] = mean(vapply(2:nDocs, function(d) {
					 edge[[d]][[3]] - edge[[d-1]][[3]]
					 }, c(1))) 			
	GiR_stats[P + 3] = mean(currentC)
	Tokens_in_Topic = tabulate(vapply(1:nDocs, function(d){
					  as.numeric(names(text[[d]]))
					  }, rep(0, nwords)), K)
	GiR_stats[(P + 4):(P + 3 + nIP)] = vapply(1:nIP, function(IP) {
								 Tokens_in_Topic %*% (currentC == IP)
								 }, c(1))
	GiR_stats[(P + 4 + nIP):(P + 3 + nIP + K)] = Tokens_in_Topic
	Tokens_in_Word = tabulate(vapply(1:nDocs, function(d){
					 text[[d]]
					 }, rep(0, nwords)), W)
	GiR_stats[(P + 4 + nIP + K):(P + 3 + nIP + K + W)] = Tokens_in_Word
	
	return(GiR_stats)
}


Inference = function(edge, node, textlist, vocabulary, nIP, K, sigma_Q, alpha, mvec, betas, nvec, prior.b.mean, prior.b.var, prior.delta, 
					out, n1, n2, n3, burn, thin, netstat, seed, plot = FALSE) {
  
  set.seed(seed)
  # trim the edge so that we only model edges after 16 days
	timeinc = c(as.numeric(edge[[1]][3]), vapply(seq(along = edge)[-1], function(d) {
  	as.numeric(edge[[d]][3]) - as.numeric(edge[[d-1]][3])
 	}, c(1)))

    edge2 = which_int(384, cumsum(timeinc)) : length(edge)
   
  # initialize alpha, mvec, delta, nvec, eta, lvec, and gammas
 	W = length(vocabulary)
  	phi = lapply(1L:K, function(k) {
		rdirichlet_cpp(1, betas * nvec)
	})
	delta = rbeta(1, prior.delta[1], prior.delta[2])
	deltamat = delta
	eta = qnorm(delta)
	
	# initialize C, theta and Z
	currentC = sample(1:nIP, K, replace = TRUE)
	theta = rdirichlet_cpp(length(edge), alpha * mvec)
	currentZ = lapply(seq(along = edge), function(d) {
    if (length(textlist[[d]]) > 0) {
    		multinom_vec(length(textlist[[d]]), theta[d, ])
    	} else {
    		multinom_vec(1, theta[d, ])
    	}
    })
	p.d = t(vapply(seq(along = edge), function(d) {
	  vapply(1L:nIP, function(IP) {
	    sum(currentZ[[d]] %in% which(currentC == IP))
	  }, c(1)) / length(currentZ[[d]])
	}, rep(1, nIP)))

    # initialize beta
    sigma = 1
    L = 3
    P = 1 + L * (2 * ("dyadic" %in% netstat) + 4 * ("triadic" %in% netstat) + 2 *("degree" %in% netstat))
    bmat = list()
	for (IP in 1L:nIP) {
		bmat[[IP]] = matrix(0, nrow = P, ncol = (n3 - burn) / thin)
		bmat[[IP]][, 1:(n3 - burn) / thin] = c(rmvnorm(1, prior.b.mean, prior.b.var))
  	}

    # to check the convergence  
    if (plot) {
     	logWZ.mat = c()							  
     	alpha.mat = c()
     	entropy.mat = c()
    }
  	
    iJi = lapply(seq(along = edge), function(d) {
  	      matrix(0, nrow = length(node), ncol = length(node))
  	      })  
    
    #start outer iteration
    for (o in 1L:out) {
      Beta.old = lapply(bmat, function(b) {
      			rowMeans(b)
         		})  
    
     # Data augmentation
      lambda = list()
      LambdaiJi = list()
      nonemptyiJi = list()
	  observediJi = list()
      for (d in edge2) {
   	 	history.t = History(edge, p.d, node, as.numeric(edge[[d]][3]))
   	 	X = lapply(node, function(i) {
  	        Netstats(history.t, node, i, netstat)
            })
   	 	XB = MultiplyXBList(X, Beta.old)     
		lambda[[d]] = lambda_cpp(p.d[d,], XB)
		
		#calculate the resampling probability
		probij = DataAug_cpp(iJi[[d]], lambda[[d]], delta, timeinc[d])
		iJi[[d]] = rbinom_mat(probij)
		iJi[[d]][as.numeric(edge[[d]][1]),] = tabulateC(as.numeric(unlist(edge[[d]][2])), length(node))
		LambdaiJi[[d]] = lambdaiJi(p.d[d,], XB, iJi[[d]])
		nonemptyiJi[[d]] = LambdaiJi[[d]][!is.na(LambdaiJi[[d]])]
		observediJi[[d]] = LambdaiJi[[d]][as.numeric(edge[[d]][1])]
		}	 
    
    textlist.raw = unlist(textlist)
    table.W = lapply(1L:K, function(k) {
      tabulateC(textlist.raw[which(unlist(currentZ) == k)], W)
      })
    
    for (i1 in 1L:n1) {
      for (d in edge2) { 
        textlist.d = textlist[[d]]
        if (length(textlist.d) > 0) {
        	topicpart.d = TopicInEqZ(K, currentZ[[d]], alpha, mvec, d)
       	wordpart.d = WordInEqZ(K, textlist.d, table.W, betas, nvec)
        } else {
        topicpart.d = 0
        wordpart.d = matrix(0, nrow = length(currentZ[[d]]), ncol = K)
        }
        edgepart.d = EdgeInEqZ(iJi[[d]], lambda[[d]], delta)
        timepart.d = TimeInEqZ(nonemptyiJi[[d]], timeinc[d])
        observed.d = ObservedInEqZ(observediJi[[d]]) 
        fixedpart = topicpart.d + edgepart.d + timepart.d + observed.d 
        for (w in 1L:length(currentZ[[d]])) {
          const.Z = fixedpart + wordpart.d[w, ]
          const.Z = const.Z - max(const.Z)
          zw.old = currentZ[[d]][w]
          zw.new = multinom_vec(1, exp(const.Z))
          if (zw.new != zw.old) {
            currentZ[[d]][w] = zw.new
            topicpart.d = TopicInEqZ(K, currentZ[[d]], alpha, mvec, d)
            if (length(textlist.d) > 0) {	
            wordpart.d = WordInEqZ(K, textlist.d, table.W, betas, nvec)
            }
            table.W = lapply(1L:K, function(k) {
      				  tabulateC(textlist.raw[which(unlist(currentZ) == k)], W)
      				  })
      		p.d[d, ] = vapply(1L:nIP, function(IP) {
	 			sum(currentZ[[d]] %in% which(currentC == IP))
	 			}, c(1)) / length(currentZ[[d]])
      		LambdaiJi[[d]] = lambdaiJi(p.d[d,], XB, iJi[[d]])
			nonemptyiJi[[d]] = LambdaiJi[[d]][!is.na(LambdaiJi[[d]])]
         	observediJi[[d]] = LambdaiJi[[d]][as.numeric(edge[[d]][1])]
          }
        }
       }
      }

    for (i2 in 1L:n2) {
      # C update given Z and B - within each document d
      for (k in unique(unlist(currentZ[edge2]))) { 
        document.k = which(vapply(currentZ, function(d){k %in% d}, c(1)) == 1)
        document.k = document.k[document.k %in% edge2]
        const.C = rep(NA, nIP)
        for (IP in 1L:nIP) {
        	if (!currentC[k] == IP){
          currentC[k] = IP
          for (d in document.k) {
            p.d[d, ] =  vapply(1L:nIP, function(IP) {
              sum(currentZ[[d]] %in% which(currentC == IP))
            }, c(1)) / length(currentZ[[d]])
          }
          for (d in edge2) {
           history.t = History(edge, p.d, node, as.numeric(edge[[d]][3]))
    	   	   X = lapply(node, function(i) {
               Netstats(history.t, node, i, netstat)
               })
    	       XB = MultiplyXBList(X, Beta.old)    
           lambda[[d]] = lambda_cpp(p.d[d,], XB)
		   LambdaiJi[[d]] = lambdaiJi(p.d[d,], XB, iJi[[d]])
		   nonemptyiJi[[d]] = LambdaiJi[[d]][!is.na(LambdaiJi[[d]])]
           observediJi[[d]] = LambdaiJi[[d]][as.numeric(edge[[d]][1])]
          }
          }
          const.C[IP] = sum(vapply(document.k, function(d) {
            EdgeInEqZ(iJi[[d]], lambda[[d]], delta) + TimeInEqZ(nonemptyiJi[[d]], timeinc[d]) + 
    			ObservedInEqZ(observediJi[[d]]) 
          }, c(1))) / length(document.k)
          }
          const.C = const.C - max(const.C)
          currentC[k] = multinom_vec(1, exp(const.C))
      }
    }
    
    if (plot) {
      entropy.mat = c(entropy.mat, entropy.empirical(currentC))
      alpha.mat = rbind(alpha.mat, alpha)
      logWZ.mat = c(logWZ.mat, logWZ(K, currentZ, textlist, table.W, alpha, mvec, betas, nvec))
      }
      
     for (d in edge2) {
     	 p.d[d, ] =  vapply(1L:nIP, function(IP) {
              		 sum(currentZ[[d]] %in% which(currentC == IP))
                     }, c(1)) / length(currentZ[[d]])
     	 history.t = History(edge, p.d, node, as.numeric(edge[[d]][3]))
    	     X = lapply(node, function(i) {
              Netstats(history.t, node, i, netstat)
             })
    	     XB = MultiplyXBList(X, Beta.old)   
    	     lambda[[d]] = lambda_cpp(p.d[d,], XB)
		 LambdaiJi[[d]] = lambdaiJi(p.d[d,], XB, iJi[[d]])
		 nonemptyiJi[[d]] = LambdaiJi[[d]][!is.na(LambdaiJi[[d]])]
         observediJi[[d]] = LambdaiJi[[d]][as.numeric(edge[[d]][1])]
	 }
	 
	 # Beta and delta update
	 prior.old1 = sum(vapply(1L:nIP, function(IP) {
		  		 dmvnorm(Beta.old[[IP]], prior.b.mean, prior.b.var, log = TRUE)
		  		 }, c(1))) 
	 post.old1 = sum(vapply(edge2, function(d) {
	     	     EdgeInEqZ(iJi[[d]], lambda[[d]], delta) + TimeInEqZ(nonemptyiJi[[d]], timeinc[d]) +
    			 	 ObservedInEqZ(observediJi[[d]])
    			 }, c(1))) / length(edge2)
    	 prior.old2 = dbeta(delta, prior.delta[1], prior.delta[2], log = TRUE) 
     post.old2 = sum(vapply(edge2, function(d) {
    			    EdgeInEqZ(iJi[[d]], lambda[[d]], delta)
    			    }, c(1))) / length(edge2)
 
	 options(warn = -1)
	 proposal.var = list()
     for (IP in 1L:nIP) {
     	 proposal.var[[IP]] = cor(t(bmat[[IP]]))
         proposal.var[[IP]][is.na(proposal.var[[IP]])] = 0
         if (sum(eigen(proposal.var[[IP]])$values < 0 ) > 0) {
        	 proposal.var[[IP]] = diag(P)
         }
         }
     proposal.var = lapply(1:nIP, function(IP){diag(P)})
     options(warn = 0)
     for (i3 in 1L:n3) {
     	Beta.new = lapply(1L:nIP, function(IP) {
           rmvnorm(1, Beta.old[[IP]], sigma_Q * proposal.var[[IP]])
         }) 
        for (d in edge2) {
           XB = MultiplyXBList(X, Beta.new)    
           LambdaiJi[[d]] = lambdaiJi(p.d[d,], XB, iJi[[d]])
		   nonemptyiJi[[d]] = LambdaiJi[[d]][!is.na(LambdaiJi[[d]])]
		   observediJi[[d]] = LambdaiJi[[d]][as.numeric(edge[[d]][1])]
         }
         prior.new1 = sum(vapply(1L:nIP, function(IP) {
        		dmvnorm(Beta.new[[IP]], prior.b.mean, prior.b.var, log = TRUE)
         }, c(1))) 
         post.new1 = sum(vapply(edge2, function(d) {
    			    EdgeInEqZ(iJi[[d]], lambda[[d]], delta) + TimeInEqZ(nonemptyiJi[[d]], timeinc[d]) + 
    			    ObservedInEqZ(observediJi[[d]])
    			    }, c(1))) / length(edge2)
    		 loglike.diff = prior.new1 + post.new1 - prior.old1 - post.old1
      
        if (log(runif(1, 0, 1)) < loglike.diff) {
         for (IP in 1L:nIP) {
         Beta.old[[IP]]  = Beta.new[[IP]]
         }
         prior.old1 = prior.new1
         post.old1 = post.new1
         }
        
        	 eta.new = rnorm(1, eta, sigma_Q)
    	 	 delta.new = pnorm(eta.new)
		 prior.new2 = dbeta(delta.new, prior.delta[1], prior.delta[2], log = TRUE)
		 post.new2 = sum(vapply(edge2, function(d) {
    			    EdgeInEqZ(iJi[[d]], lambda[[d]], delta.new)
    			    }, c(1))) / length(edge2)
    	 	 loglike.diff2 = prior.new2 + post.new2 - prior.old2 - post.old2
		  if (log(runif(1, 0, 1)) < loglike.diff2) {
			delta = delta.new
			eta = eta.new
			prior.old2 = prior.new2
			post.old2 = post.new2
			}

         if (i3 > burn & i3 %% (thin) == 0) {
          for (IP in 1L:nIP) {
           bmat[[IP]][ , (i3 - burn) / thin] = Beta.old[[IP]]
           }
           deltamat = c(deltamat, delta)
          }
         }
	 }
    
    if (plot) {
    burnin = round(out / 10)
    par(mfrow = c(2, 2))
  	plot(entropy.mat[-1L:-burnin], type = "l", 
  	     xlab = "(Outer) Iterations", ylab = "Entropy of IP")
  	abline(h = mean(entropy.mat[-1L:-burnin]), lty = 1)
  	title("Convergence of Entropy")
  	matplot(alpha.mat[-1L:-burnin,], lty = 1, type = "l", col = 1L:nIP, 
  	        xlab = "(Outer) Iterations", ylab = "alpha")
  	abline(h = mean(alpha.mat[-1L:-burnin,]), lty = 1, col = 1L:nIP)
	title("Convergence of Optimized alpha")
	plot(logWZ.mat[-1L:-burnin], type = "l", 
	       xlab = "(Outer) Iterations", ylab = "logWZ")
	abline(h = mean(logWZ.mat[-1L:-burnin]), lty = 1)
	title("Convergence of logWZ")
	matplot(bmat[[1]][1,], lty = 1, col = 1L:P, type = "l", 
	          main = "Traceplot of beta", xlab = "(Inner) Iterations", ylab = "")
	abline(h = mean(bmat[[1]][1,]), lty = 1, col = 1L)
	  }
     
  chain.final = list(C = currentC, Z = lapply(edge2, function(d) {currentZ[[d]]}), B = bmat, D = deltamat)

  return(chain.final)
}

#PP_Plots
GiR_PP_Plots = function(Forward_stats, Backward_stats) {
	nms = colnames(Forward_stats)
	
	for (i in 1L:ncol(Forward_stats)) {
		if (nrow(Forward_stats) > 20000) {
			thin = seq(from = floor(nrow(Forward_stats)/10), to = nrow(Forward_stats), length.out = 10000)
			Forward_test = Forward_stats[thin, i]
			Backward_test = Backward_stats[thin, i]
		} else {
			Forward_test = Forward_stats[, i]
			Backward_test = Backward_stats[, i]
		}
		
		all = c(Backward_stats[, i], Forward_stats[, i])
		ylims = c(min(all) - 0.1 * max(abs(all)), max(all) + 0.1 * max(abs(all)))
		xlims = ylims
		
		quantiles = 50
		if (grepl("B_", nms[i]) ) {
			quantiles = 1000
		}
		
		qqplot(x = quantile(Forward_stats[, i], seq(0, 1, length = quantiles)),
			   y = quantile(Backward_stats[, i], seq(0, 1, length = quantiles)),
			   ylim = ylims,
			   xlim = xlims,
			   ylab = "Backward",
			   xlab = "Forward",
			   col = "blue",
			   pch = 19,
			   main = nms[i],
			   cex.lab = 1,
			   cex.axis = 1,
			   cex.main = 1)
		lines(x = xlims, y = ylims, col = "red", lwd = 3)
		text(paste("Backward Mean:", round(mean(Backward_stats[,i]), 4),
				   "\nForward Mean:", round(mean(Forward_stats[,i]), 4),
				   "\nt-test p-value:", round(t.test(Backward_test, Forward_test)$p.value, 4),
				   "\nMann-Whitney p-value:", round(wilcox.test(Backward_test, Forward_test)$p.value,4)),
				   x = xlims[2] - 0.35 * abs(xlims[2] - xlims[1]),
				   y = ylims[1] + 0.15 * abs(ylims[2] - ylims[1]),
				   cex = 0.3)
	}
}      


#PP_Plots
GiR_PP_Plots2 = function(Forward_stats, Backward_stats) {
	nms = colnames(Forward_stats)
	
	for (i in 1L:ncol(Forward_stats)) {
		if (nrow(Forward_stats) > 20000) {
			thin = seq(from = floor(nrow(Forward_stats)/10), to = nrow(Forward_stats), length.out = 10000)
			Forward_test = Forward_stats[thin, i]
			Backward_test = Backward_stats[thin, i]
		} else {
			Forward_test = Forward_stats[, i]
			Backward_test = Backward_stats[, i]
		}
		
		all = c(Backward_stats[, i], Forward_stats[, i])
		ylims = c(min(all) - 0.1 * max(abs(all)), max(all) + 0.1 * max(abs(all)))
		xlims = ylims
		
		quantiles = 50
		if (grepl("B_", nms[i]) ) {
			quantiles = 1000
		}
		
		normalmean = mean(c(quantile(Forward_stats[, i], seq(0, 1, length = quantiles)), 
							quantile(Backward_stats[, i], seq(0, 1, length = quantiles))))
		normalvar = sd(c(quantile(Forward_stats[, i], seq(0, 1, length = quantiles)), 
							quantile(Backward_stats[, i], seq(0, 1, length = quantiles))))

		qqplot(x = pnorm(sort(quantile(Forward_stats[, i], seq(0, 1, length = quantiles))), normalmean, normalvar),
			   y = pnorm(sort(quantile(Backward_stats[, i], seq(0, 1, length = quantiles))), normalmean, normalvar),
			   ylim = pnorm(ylims, normalmean, normalvar),
			   xlim = pnorm(xlims, normalmean, normalvar),
			   ylab = "Backward",
			   xlab = "Forward",
			   col = "blue",
			   pch = 19,
			   main = nms[i],
			   cex.lab = 1,
			   cex.axis = 1,
			   cex.main = 1)
		lines(x = pnorm(xlims, normalmean, normalvar), y = pnorm(ylims, normalmean, normalvar), col = "red", lwd = 3)
		text(paste("Backward Mean:", round(mean(Backward_stats[,i]), 4),
				   "\nForward Mean:", round(mean(Forward_stats[,i]), 4),
				   "\nt-test p-value:", round(t.test(Backward_test, Forward_test)$p.value, 4),
				   "\nMann-Whitney p-value:", round(wilcox.test(Backward_test, Forward_test)$p.value,4)),
				   x = pnorm(xlims[2], normalmean, normalvar) - 0.35 * abs(pnorm(xlims[2], normalmean, normalvar) - pnorm(xlims[1], normalmean, normalvar)),
				   y = pnorm(ylims[1], normalmean, normalvar) + 0.15 * abs(pnorm(ylims[2], normalmean, normalvar) - pnorm(ylims[1], normalmean, normalvar)),
				   cex = 0.3)
	}
}      


# Getting_It_Right
GiR = function(Nsamp = 5000, nDocs = 10, node = 1:4, vocabulary =  c("hi", "hello","bye", "mine", "what"), 
			   nIP = 2, K = 4, nwords = 4, alpha = 2, mvec = rep(1/4, 4), betas = 2, nvec = rep(1/5, 5), 
			   prior.b.mean = c(-3, rep(0, 6)), prior.b.var = diag(7), prior.delta = c(2, 2), sigma_Q = 0.25, 
			   niters = c(1, 1, 1, 50, 0, 1), netstat = "dyadic", generate_PP_plots = TRUE, seed = 1) {
	
	set.seed(seed)
	b = lapply(1L:nIP, function(IP) {
			   c(rmvnorm(1, prior.b.mean, prior.b.var))
			   })
	delta = 	rbeta(1, prior.delta[1], prior.delta[2])
	currentC = sample(1L:nIP, K, replace = TRUE)	 
	base.data = GenerateDocs(nDocs, node, vocabulary, nIP, K, nwords, alpha, mvec, betas, nvec, b, delta, currentC, netstat, seed, 
							 base.edge = list(), base.text = list(), base = TRUE)
	base.edge = base.data$edge	   
	base.text = base.data$text

	#Forward sampling
	Forward_stats = matrix(NA, nrow = Nsamp, ncol = 21)
	colnames(Forward_stats) = c(paste0("B_",1:length(b[[1]])), "Mean_receipients", "Mean_timediff", "Mean_TopicIP", 
						  paste0("Tokens_in_IP_", 1:nIP), paste0("Tokens_in_Topic", 1:K), 
						  paste0("Tokens_in_Word", 1:length(vocabulary)))
	for (i in 1:Nsamp) { 
		if (i %% 5000 == 0) {cat("Forward sampling", i, "\n")}
		set.seed(i)
		b = lapply(1:nIP, function(IP) {
			c(rmvnorm(1, prior.b.mean, prior.b.var))
			})
		delta = rbeta(1, prior.delta[1], prior.delta[2])
		currentC = sample(1L:nIP, K, replace = TRUE)
		Forward_sample = GenerateDocs(nDocs, node, vocabulary, nIP, K, nwords, alpha, mvec, betas, nvec, b, delta, currentC, netstat, 
									  base.edge = base.edge, base.text = base.text, seed, forward = TRUE) 
		Forward_stats[i, ] = GiR_stats(Forward_sample, K, currentC, vocabulary, forward = TRUE, backward = FALSE)
	}
	
	#Backward sampling
	Backward_stats = matrix(NA, nrow = Nsamp, ncol = 21)
	Backward_sample = GenerateDocs(nDocs, node, vocabulary, nIP, K, nwords, alpha, mvec, betas, nvec, b, delta, currentC, netstat, 
								   base.edge = base.edge, base.text = base.text, seed, backward_init = TRUE) 
	for (i in 1:Nsamp) { 
		if (i %% 500 == 0) {cat("Backward sampling", i, "\n")}
		seed = seed + 100
		set.seed(seed)
		Inference_samp = Inference(Backward_sample$edge, node, Backward_sample$text, vocabulary, nIP, K, sigma_Q, 
						  alpha, mvec, betas, nvec, prior.b.mean, prior.b.var, prior.delta,
						  out = niters[1], n1 = niters[2], n2 = niters[3], n3 = niters[4], burn = niters[5], thin = niters[6], 
						  netstat, seed)
		b = lapply(1:nIP, function(IP) {
			rowMeans(Inference_samp$B[[IP]])
		})
		delta = mean(Inference_samp$D)
		currentC = Inference_samp$C
		topic_token_assignments = Inference_samp$Z
		for (d in 1:length(topic_token_assignments)) {
			names(topic_token_assignments[[d]]) = Backward_sample$text[[d + length(base.text)]]
		}
		topic_token_counts = tabulate(unlist(Inference_samp$Z), K)
		word_type_topic_counts = matrix(NA , length(vocabulary), K)
		for (w in 1:length(vocabulary)) {
			word_type_w = which(unlist(Backward_sample$text[-(1:length(base.text))]) == w)
			word_type_topic_counts[w, ] = tabulate(unlist(Inference_samp$Z)[word_type_w], K)
		}
		
		Backward_sample = GenerateDocs(nDocs, node, vocabulary, nIP, K, nwords, alpha, mvec, betas, nvec, b, delta, currentC, netstat, 
								   	   base.edge = base.edge, base.text = base.text, seed, topic_token_assignments = topic_token_assignments, 
								   	   topic_token_counts = topic_token_counts, word_type_topic_counts = word_type_topic_counts, 
								   	   forward = FALSE, backward = TRUE)
		Backward_stats[i, ] = GiR_stats(Backward_sample, K, currentC, vocabulary, forward = FALSE, backward = TRUE)
	}
	
	tstats = rep(0, ncol(Forward_stats))
	wstats = rep(0, ncol(Forward_stats))
	for (j in 1:ncol(Forward_stats)) {
		thin = seq(from = floor(Nsamp / 5), to = Nsamp, length.out = 400)
		Forward_test = Forward_stats[thin, j]
		Backward_test = Backward_stats[thin, j]
		tstats[j] = t.test(Backward_test, Forward_test)$p.value
		wstats[j] = wilcox.test(Backward_test, Forward_test)$p.value
	}
	names(tstats) = names(wstats) = colnames(Forward_stats)						
	if (generate_PP_plots) {
		par(mfrow=c(3,7), oma = c(3,3,3,3), mar = c(5,5,4,1))
		GiR_PP_Plots2(Forward_stats, Backward_stats)
	}			
	return(list(Forward = Forward_stats, Backward = Backward_stats, tstats = tstats, wstats = wstats))	
}

unix.time(TryGiR <- GiR(10^5, seed = 15))
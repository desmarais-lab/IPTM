#' @title GenerateDocs.predict
#' @description Generate a collection of documents according to the generative process of IPTM using Gibbs measure
#'
#' @param nDocs number of documents to be generated
#' @param node nodelist containing the ID of nodes (ID starting from 1)
#' @param vocabulary all vocabularies used over the corpus
#' @param nIP total number of interaction patterns specified by the user
#' @param K total number of topics specified by the user
#' @param nwords number of words in a document (fixed constant for now)
#' @param alpha Dirichlet concentration prior for document-topic distribution
#' @param mvec Dirichlet base prior for document-topic distribution
#' @param betas Dirichlet concentration prior for topic-word distribution
#' @param nvec Dirichlet base prior for topic-word distribution
#' @param b List of coefficients estimating the history effect for each IP 
#' @param delta tuning parameter controlling the number of recipients
#' @param currentC topic-IP assignment
#' @param netstat which type of network statistics to use ("intercept", "dyadic", "triadic", "degree")
#' @param base.edge edges before 384 hours that is used to calculate initial history of interactions
#' @param base.text texts corresponding to base.edge
#' @param topic_token_assignments matrix of topic-token assignments
#' @param forward Logigal indicating whether we are generating forward samples
#' @param backward_init Logigal indicating whether we are generating initial backward samples
#' @param backward Logigal indicating whether we are generating backward samples
#' @param base Logical indicating whether or not we are generating base edges (< 384)
#' @param support support of latent edge vector Ji
#'
#' @return generated edge and text, parameter b used to generate those, and base (if base == TRUE)
#'
#' @export
GenerateDocs.predict = function(nDocs, node, vocabulary, nIP, K, nwords, alpha, mvec, betas, nvec,
                        b, delta, currentC, netstat, base.edge, base.text, topic_token_assignments = NULL,
                        forward = FALSE, backward_init = FALSE, backward = FALSE, base = FALSE, support) {
  W = length(vocabulary)
  phi = lapply(1L:K, function(k) {
    rdirichlet_cpp(1, betas * nvec)
  })
  L = 3
  P = 1 * ("intercept" %in% netstat) + L * (2 * ("dyadic" %in% netstat) + 4 * ("triadic" %in% netstat) + 2 *("degree" %in% netstat))
  
  if (base) {
    t.d = 0
  } else {
    t.d = 384
  }
  edge = base.edge
  text = base.text
  base.length = length(edge)
  p.d = matrix(NA, nrow = nDocs + base.length, ncol = nIP)
  if (!base) {
  for (d in 1:base.length) {
  	 p.d[d, ] = vapply(1L:nIP, function(IP) {
      sum(names(text[[d]]) %in% which(currentC == IP))
    }, c(1)) / max(1, length(text[[d]]))
  }
  }
 history.t = lapply(1:nIP, function(IP) {
        	 	 lapply(1:3, function(l){
         		matrix(0, length(node), length(node))
       			})
     		 }) 
  if (backward) {
   	topic_counts  = tabulateC(unlist(topic_token_assignments), K)  
  }			
  word_type_topic_counts = matrix(0, W, K)
  iJi = list()
  for (d in 1:nDocs) {
    N.d = nwords
    text[[base.length + d]] = rep(NA, N.d)
    iJi[[base.length + d]] = matrix(0, length(node), length(node))
    if (!backward) {
      theta.d = rdirichlet_cpp(1, alpha * mvec)	
      topic.d = multinom_vec(max(1, N.d), theta.d)
        for (n in 1:N.d){
          text[[base.length + d]][n] = multinom_vec(1, phi[[topic.d[n]]])
        }
        names(text[[base.length + d]]) = topic.d
    } else {
      phi.k = rep(NA, K)
      topic.d = topic_token_assignments[[d]]
        for (n in 1:N.d){
          for (w in 1:W) {
            phi.k[w] = (word_type_topic_counts[w, topic.d[n]] + betas * nvec[w]) / (topic_counts[topic.d[n]] + betas)
          } 
          text[[base.length + d]][n] = multinom_vec(1, phi.k)
          word_type_topic_counts[text[[base.length + d]][n], topic.d[n]] = word_type_topic_counts[text[[base.length + d]][n], topic.d[n]] + 1
        }
        names(text[[base.length + d]]) = topic.d
    }
    
    p.d[base.length + d, ] = vapply(1L:nIP, function(IP) {
      sum(names(text[[base.length + d]]) %in% which(currentC == IP))
    }, c(1)) / max(1, N.d)
    if (t.d >= 384) {
      history.t = History(edge, p.d, node, t.d + exp(-745))
    }
    X = lapply(node, function(i) {
      	Netstats(history.t, node, i, netstat)
   		})
    XB = MultiplyXBList(X, b)     
    lambda = lambda_cpp(p.d[base.length + d,], XB)
    
    	for (i in node) {
    		iJi[[base.length + d]][i, -i] = r.gibbs.measure(1, lambda[i, -i], delta, support)
    	}
    LambdaiJi = lambdaiJi(p.d[base.length + d,], XB, iJi[[base.length + d]])
    i.d = multinom_vec(1, LambdaiJi)
    j.d = which(iJi[[base.length + d]][i.d,] == 1)
    t.d = t.d + rexp(1, sum(LambdaiJi))
    edge[[base.length + d]] = list(sender = i.d, receiver = j.d, timestamp = t.d)		
  }
  if (base == TRUE & t.d > 384) {
    cutoff = which_int(384, vapply(1:length(edge), function(d) {edge[[d]][[3]]}, c(1))) - 1
    edge = edge[1:cutoff]
    text = text[1:cutoff]
  }
  return(list(edge = edge[[length(edge)]], text = text[[length(text)]], b = b, d = delta, iJi = iJi[[length(iJi)]]))							
} 

#' @title IPTM_predict.data
#' @description posterior predictive experiments
#'
#' @param D Dth document the user wants to predict
#' @param O number of outer iterations of inference from which to generate predictions
#' @param R the number of iterations to sample predicted data within each outer iteration
#' @param edge list of document information with 3 elements (element 1 sender, element 2 receiver, element 3 time in unix.time format)
#' @param node nodelist containing the ID of nodes (ID starting from 1)
#' @param textlist list of text (length=number of documents in total) containing the words in each document
#' @param vocabulary all vocabularies used over the corpus
#' @param nIP total number of interaction patterns specified by the user
#' @param K total number of topics specified by the user
#' @param sigma_Q proposal distribution variance parameter for beta and delta
#' @param alpha Dirichlet concentration prior for document-topic distribution
#' @param mvec Dirichlet base prior for document-topic distribution
#' @param betas Dirichlet concentration prior for topic-word distribution
#' @param nvec Dirichlet base prior for topic-word distribution
#' @param prior.b.mean mean vector of b in multivariate normal distribution
#' @param prior.b.var covairance matrix of b in multivariate normal distribution
#' @param prior.delta parameter of delta in Normal prior
#' @param out size of outer iterations 
#' @param n_B size of third inner iteration for updates of B
#' @param n_d size of third inner iteration for updates of delta
#' @param burn iterations to be discarded at the beginning of beta and delta chain 
#' @param thinning the thinningning interval of beta and delta chain
#' @param netstat which type of network statistics to use ("intercept", dyadic", "triadic", "degree")
#' @param plot to plot the convergence diagnostics or not (TRUE/FALSE)
#' @param optimize to optimize alpha (Dirichlet concentration prior for document-topic distribution) or not (TRUE/FALSE)
#'
#' @return prediction output 
#'
#' @export
IPTM_predict.data = function(D, O, R, edge, node, textlist, vocabulary, nIP, K, sigma_Q, alpha, mvec, betas, nvec, 
							 prior.b.mean, prior.b.var, prior.delta, 
                                out, n_B, n_d, burn, thinning, netstat, plot = FALSE, optimize = TRUE) {
   
    Inference_samp = IPTM_inference.data(edge[1:(D-1)], node, textlist[1:(D-1)], vocabulary, nIP, K,
    									 sigma_Q, alpha, mvec, betas, nvec, prior.b.mean, prior.b.var, prior.delta,
                               		 out, n_B, n_d, burn, thinning, netstat, plot = FALSE, optimize = TRUE)
    supportD = gibbs.measure.support(length(node) - 1)                          	
	b = lapply(1:nIP, function(IP) {
        rowMeans(Inference_samp$B[[IP]])
    })
    delta = mean(Inference_samp$D)
    currentC = Inference_samp$C
    topic_token_assignments = Inference_samp$Z

	New_sample = GenerateDocs.predict(1, node, vocabulary, nIP, K, nwords = length(textlist[[D]]), alpha, mvec, betas, nvec, b, 
    										 delta, currentC, netstat, base.edge = edge[1:(D-1)], base.text = textlist[1:(D-1)],
    										 topic_token_assignments = topic_token_assignments, 
                             			 forward = TRUE, backward = FALSE, base = FALSE, support = supportD)
 return(New_sample)                                	
}
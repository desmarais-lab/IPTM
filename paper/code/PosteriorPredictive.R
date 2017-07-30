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
IPTM_predict.data = function(D, O, R,
							 edge, node, textlist, vocabulary, nIP, K, sigma_Q, alpha, mvec, betas, nvec, prior.b.mean, prior.b.var, prior.delta, 
                                out, n_B, n_d, burn, thinning, netstat, plot = FALSE, optimize = TRUE) {
   
    Inference_samp = IPTM_inference.data(edge[1:(D-1)], node, textlist[1:(D-1)], vocabulary, nIP, K,
    									 sigma_Q, alpha, mvec, betas, nvec, prior.b.mean, prior.b.var, prior.delta,
                               			 out, n_B, n_d, burn, thinning, netstat, plot = FALSE, optimize = TRUE)
                             	
	b = lapply(1:nIP, function(IP) {
        rowMeans(Inference_samp$B[[IP]])
    })
    delta = mean(Inference_samp$D)
    currentC = Inference_samp$C
    topic_token_assignments = Inference_samp$Z




                                	
}
library(MLmetrics)
library(pROC)
load('Enron.RData')
set.seed(1)
set.seed(1)
missing = matrix(0, nrow = length(Enron$edge), ncol = 3)
missing[sample(556:6613, 1212, replace = FALSE), 1] = 1
missing[sample(556:6613, 1212, replace = FALSE), 2] = 1
missing[sample(556:6613, 1212, replace = FALSE), 3] = 1

setwd('/Users/bomin8319/Desktop/IPTM/experiment/Enron_PPE')

#sender score
AUC = list()
for (nIP in 1:2) {
	AUC[[nIP]] = list()
    	iter = 1
        for (K in c(5, 10, 25, 50, 75, 100)){
        	AUC[[nIP]][[iter]] = matrix(0, nrow = 1212, ncol = 20)
            filename = paste0("Enron_PPE_",nIP,"_",K,".RData")
            load(filename)
            d = 1
for (D in which(missing[,1]==1)) {
    observeds = tabulate(Enron$edge[[D]]$sender, 30)
            for (o in 1:5) {
            for (r in 20:30) {
                    AUC[[nIP]][[iter]][d,r-19] = as.numeric(multiclass.roc(observeds, tabulate(EnronPPE[[o]]$senderpredict[d,r], 30))$auc)
                }
        }
        d = d+1
    }
 iter = iter + 1
 }
}


#receiver score
F1score = list()
for (nIP in 1:2) {
	F1score[[nIP]] = list()
    	iter = 1
        for (K in c(5, 10, 25, 50, 75, 100)){
        	F1score[[nIP]][[iter]] = matrix(0, nrow = 1212, ncol = 20)
            filename = paste0("Enron_PPE_",nIP,"_",K,".RData")
            load(filename)
            d = 1
for (D in which(missing[,2]==1)) {
    observedr = tabulate(Enron$edge[[D]]$receiver, 30)
            for (o in 1:5) {
            for (r in 20:30) {
                	F1 =  F1_Score(observedr, tabulate(EnronPPE[[o]]$receiverpredict[[d]][r,], 30))
					if (is.na(F1)) {F1 = 0}
        			F1score[[nIP]][[iter]][d,r-19] = F1 
                }
        	}
        	d = d+1
    }
    iter = iter + 1
	}
}


#timestamp score
MSE = list()
for (nIP in 1:2) {
	MSE[[nIP]] = list()
    	iter = 1
        for (K in c(5, 10, 25, 50, 75, 100)){
        	MSE[[nIP]][[iter]] = matrix(0, nrow = 1212, ncol = 20)
            filename = paste0("Enron_PPE_",nIP,"_",K,".RData")
            load(filename)
             d = 1
for (D in which(missing[,3]==1)) {
    observedt = Enron$edge[[D]]$unixtime
 for (o in 1:5) {
            for (r in 20:30) {
        			MSE[[nIP]][[iter]][d,r-19] = abs(log(EnronPPE[[o]]$timepredict[d,r]) - log(observedt))/(3600*24)
                }
        	}
        	d = d+1
        }
        iter = iter+1
    }
}


#plot
AUC = AUC / 1212
F1score = F1score / 1212
MSE = MSE / 1212

library(ggplot2)
library(gridExtra)
AUC_new = data.frame(AUC = c(t(AUC))/100, nIP = as.factor(c(sapply(1:3, function(x) rep(x, 6)))), K =as.factor(rep(c(5, 10, 20, 30, 40, 50), 3)))
F1_new = data.frame(F1 = c(t(F1score))/100, nIP = as.factor(c(sapply(1:3, function(x) rep(x, 6)))), K =as.factor(rep(c(5, 10, 20, 30, 40, 50), 3)))
MAD_new = data.frame(MAD = c(t(MSE))/100, nIP = as.factor(c(sapply(1:3, function(x) rep(x, 6)))), K =as.factor(rep(c(5, 10, 20, 30, 40, 50), 3)))

f[[1]] <- ggplot(AUC_new[-13:-18,], aes(K, AUC, col = nIP))+geom_line(aes(group=nIP)) + geom_point(size = 1) 

f[[2]] <- ggplot(F1_new[-13:18,], aes(K, F1, col = nIP))+geom_line(aes(group=nIP)) + geom_point(size = 1) 

f[[3]] <- ggplot(MAD_new[-13:18,], aes(K, MAD, col = nIP))+geom_line(aes(group=nIP)) + geom_point(size = 1) 
 marrangeGrob(f[1:3], nrow = 1, ncol = 3, top = NULL)
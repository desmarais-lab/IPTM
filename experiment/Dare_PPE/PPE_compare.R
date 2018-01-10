library(MLmetrics)
library(pROC)
load('Darenew.RData')
set.seed(1)
missing = matrix(0, nrow = length(Dare$edge), ncol = 3)
missing[sample(391:2247, 371, replace = FALSE), 1] = 1
missing[sample(391:2247, 371, replace = FALSE), 2] = 1
missing[sample(391:2247, 371, replace = FALSE), 3] = 1

setwd('/Users/bomin8319/Desktop/IPTM/experiment/Dare_PPE')

#sender score
AUC = matrix(0, 3, 6)
for (nIP in 1:3) {
    	iter = 1
        for (K in c(5, 10, 20, 30, 40, 50)){
            filename = paste0("Dare_PPE_",nIP,"_",K,".RData")
            load(filename)
            d = 1
for (D in which(missing[,1]==1)) {
    observeds = tabulate(Dare$edge[[D]]$sender, 27)
            for (o in 1:5) {
            for (r in 20:30) {
                    AUC[nIP, iter] = AUC[nIP, iter] + as.numeric(multiclass.roc(observeds, tabulate(DarePPE[[o]]$senderpredict[d,r], 27))$auc)
                }
        }
        d = d+1
    }
 iter = iter + 1
 }
}


#receiver score
F1score = matrix(0, 3, 6)
for (nIP in 1:3) {
    	iter = 1
        for (K in c(5, 10, 20, 30, 40, 50)){
            filename = paste0("Dare_PPE_",nIP,"_",K,".RData")
            load(filename)
            d = 1
for (D in which(missing[,2]==1)) {
    observedr = tabulate(Dare$edge[[D]]$receiver, 27)
            for (o in 1:5) {
            for (r in 20:30) {
                	F1 =  F1_Score(observedr, DarePPE[[o]]$receiverpredict[[d]][r,])
					if (is.na(F1)) {F1 = 0}
        			F1score[nIP, iter]  = F1score[nIP, iter]  + F1 
                }
        	}
        	d = d+1
    }
    iter = iter + 1
	}
}


#timestamp score
MSE = matrix(0, 3, 6)
for (nIP in 1:3) {
    	iter = 1
        for (K in c(5, 10, 20, 30, 40, 50)){
            filename = paste0("Dare_PPE_",nIP,"_",K,".RData")
            load(filename)
             d = 1
for (D in which(missing[,3]==1)) {
    observedt = Dare$edge[[D]]$unixtime
 for (o in 1:5) {
            for (r in 20:30) {
        			MSE[nIP, iter]  = MSE[nIP, iter]  + abs(log(DarePPE[[o]]$timepredict[d,r]) - log(observedt))/3600
                }
        	}
        	d = d+1
        }
        iter = iter+1
    }
}


#plot
AUC = AUC / 371
F1score = F1score / 371
MSE = MSE / 371

library(ggplot2)
library(gridExtra)
AUC_new = data.frame(AUC = c(t(AUC))/100, nIP = as.factor(c(sapply(1:3, function(x) rep(x, 6)))), K =as.factor(rep(c(5, 10, 20, 30, 40, 50), 3)))
F1_new = data.frame(F1 = c(t(F1score))/100, nIP = as.factor(c(sapply(1:3, function(x) rep(x, 6)))), K =as.factor(rep(c(5, 10, 20, 30, 40, 50), 3)))
MAD_new = data.frame(MAD = c(t(MSE))/100, nIP = as.factor(c(sapply(1:3, function(x) rep(x, 6)))), K =as.factor(rep(c(5, 10, 20, 30, 40, 50), 3)))

f[[1]] <- ggplot(AUC_new[-18,], aes(K, AUC, col = nIP))+geom_line(aes(group=nIP)) + geom_point(size = 1) 

f[[2]] <- ggplot(F1_new[-18,], aes(K, F1, col = nIP))+geom_line(aes(group=nIP)) + geom_point(size = 1) 

f[[3]] <- ggplot(MAD_new[-18,], aes(K, MAD, col = nIP))+geom_line(aes(group=nIP)) + geom_point(size = 1) 
 marrangeGrob(f[1:3], nrow = 1, ncol = 3, top = NULL)
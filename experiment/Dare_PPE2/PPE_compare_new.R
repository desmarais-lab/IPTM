library(MLmetrics)
library(pROC)
load('Darenew.RData')
set.seed(1)
missing = matrix(0, nrow = length(Dare$edge), ncol = 3)
missing[sample(391:2247, 371, replace = FALSE), 1] = 1
missing[sample(391:2247, 371, replace = FALSE), 2] = 1
missing[sample(391:2247, 371, replace = FALSE), 3] = 1

setwd('/Users/bomin8319/Desktop/IPTM/experiment/PPE_result')

#sender score
AUC = list()
for (nIP in 1:3) {
	AUC[[nIP]] = matrix(0, 371, 6)
    	iter = 1
        for (K in c(5, 10, 20, 30, 40, 50)){
            filename = paste0("Dare_PPE_",nIP,"_",K,".RData")
            load(filename)
            d = 1
for (D in which(missing[,1]==1)) {
    observeds = tabulate(Dare$edge[[D]]$sender, 27)
      for (o in 1:5) {
    AUC[[nIP]][d,iter] = AUC[[nIP]][d,iter] + as.numeric(multiclass.roc(observeds, tabulate(DarePPE[[o]]$senderpredict[d,10:30],27))$auc)
    }
      d = d+1
    }
 iter = iter + 1
 }
 AUC[[nIP]] = AUC[[nIP]]/5
}

#receiver score
F1score = list()
for (nIP in 1:3) {
	F1score[[nIP]] = matrix(0, 371, 6)
    	iter = 1
        for (K in c(5, 10, 20, 30, 40, 50)){
            filename = paste0("Dare_PPE_",nIP,"_",K,".RData")
            load(filename)
            d = 1
for (D in which(missing[,2]==1)) {
    observedr = tabulate(Dare$edge[[D]]$receiver, 27)
    for (o in 1:5) {
    receivers = colSums(DarePPE[[o]]$receiverpredict[[d]][10:30,])/20
    receivers2 = as.numeric(receivers>0.5)
    if (sum(receivers2)==0) {receivers2[which(receivers == max(receivers))] = 1}
        F1 =  F1_Score(observedr,receivers2)
		if (is.na(F1)) {F1 = 0}
        	F1score[[nIP]][d,iter]  = F1score[[nIP]][d,iter]  + F1 
           	}
           	d = d+1    	
    }
    iter = iter + 1
	}
	F1score[[nIP]] = F1score[[nIP]]/5
}


#timestamp score
MSE = list()
for (nIP in 1:3) {
	MSE[[nIP]] = matrix(0, 371, 6)
    	iter = 1
        for (K in c(5, 10, 20, 30, 40, 50)){
            filename = paste0("Dare_PPE_",nIP,"_",K,".RData")
            load(filename)
             d = 1
for (D in which(missing[,3]==1)) {
    observedt = Dare$edge[[D]]$unixtime - Dare$edge[[D-1]]$unixtime
	 for (o in 1:5) {
        	MSE[[nIP]][d,iter]  = MSE[[nIP]][d,iter]  + median(abs(DarePPE[[o]]$timepredict[d,10:30] - observedt)/3600)
                }
        	d = d+1
        }
        iter = iter+1
    }
    MSE[[nIP]] = MSE[[nIP]]/5
}

AUC_new = NULL
F1_new = NULL
MAD_new = NULL
for (nIP in 1:3){
	iter = 1
	for (K in c(5,10,20,30,40,50)){
			AUC_new = rbind(AUC_new, cbind(mean(AUC[[nIP]][,iter]), nIP, K))
			F1_new = rbind(F1_new, cbind(mean(F1score[[nIP]][,iter]), nIP, K))
			MAD_new = rbind(MAD_new, cbind(mean(MSE[[nIP]][,iter][MSE[[nIP]][,iter]!=Inf]), nIP, K))
			iter = iter + 1
	}	
}
AUC_new = data.frame(AUC = AUC_new[,1], C = as.factor(AUC_new[,2]), K = as.factor(AUC_new[,3]))
F1_new = data.frame(F1 = F1_new[,1], C = as.factor(F1_new[,2]), K = as.factor(F1_new[,3]))
MAD_new = data.frame(MAD = MAD_new[,1], C = as.factor(MAD_new[,2]), K = as.factor(MAD_new[,3]))



#plot
f = list()
library(ggplot2)
library(gridExtra)

f[[1]] <- ggplot(AUC_new, aes(K, AUC, col = C))+geom_line(aes(group = C))+geom_point()
f[[2]] <- ggplot(F1_new, aes(K, F1, col = C))+geom_line(aes(group = C))+geom_point()
f[[3]] <- ggplot(MAD_new, aes(K, MAD, col = C))+geom_line(aes(group = C))+geom_point()

marrangeGrob(f[1:3], nrow = 1, ncol = 3, top = NULL)




#plot
f = list()
library(ggplot2)
library(gridExtra)

f[[1]] <- ggplot(AUC_new, aes(K, AUC, col = C))+geom_boxplot(outlier.size = 0.5, position = position_dodge()) 
f[[2]] <- ggplot(F1_new, aes(K, F1, col = C))+geom_boxplot(outlier.size = 0.5, position = position_dodge())
f[[3]] <- ggplot(MAD_new, aes(K, MAD, col = C))+geom_boxplot(outlier.size = 0.5, position = position_dodge())

marrangeGrob(f[1:3], nrow = 1, ncol = 3, top = NULL)


###########################################################


load('Enron.RData')
set.seed(1)
missing = matrix(0, nrow = length(Enron$edge), ncol = 3)
missing[sample(556:6613, 1212, replace = FALSE), 1] = 1
missing[sample(556:6613, 1212, replace = FALSE), 2] = 1
missing[sample(556:6613, 1212, replace = FALSE), 3] = 1


setwd('/Users/bomin8319/Desktop/IPTM/experiment/PPE_result')

#sender score
AUC = list()
for (nIP in 1:2) {
	AUC[[nIP]] = matrix(0, 1212, 6)
    	iter = 1
        for (K in c(5, 10, 25, 50, 75, 100)){
            filename = paste0("Enron_PPE_",nIP,"_",K,".RData")
            load(filename)
            d = 1
for (D in which(missing[,1]==1)) {
    observeds = tabulate(Enron$edge[[D]]$sender, 30)
      for (o in 1:5) {
    AUC[[nIP]][d,iter] = AUC[[nIP]][d,iter] + as.numeric(multiclass.roc(observeds, tabulate(EnronPPE[[o]]$senderpredict[d,10:30],30))$auc)
    }
      d = d+1
    }
 iter = iter + 1
 }
 AUC[[nIP]] = AUC[[nIP]]/5
}

#receiver score
F1score = list()
for (nIP in 1:2) {
	F1score[[nIP]] = matrix(0, 1212, 6)
    	iter = 1
        for (K in c(5, 10, 25, 50, 75, 100)){
            filename = paste0("Enron_PPE_",nIP,"_",K,".RData")
            load(filename)
            d = 1
for (D in which(missing[,2]==1)) {
    observedr = tabulate(Enron$edge[[D]]$receiver, 30)
    for (o in 1:5) {
    receivers = colSums(EnronPPE[[o]]$receiverpredict[[d]][10:30,])/20
    receivers2 = as.numeric(receivers>0.5)
    if (sum(receivers2)==0) {receivers2[which(receivers == max(receivers))] = 1}
        F1 =  F1_Score(observedr,receivers2)
		if (is.na(F1)) {F1 = 0}
        	F1score[[nIP]][d,iter]  = F1score[[nIP]][d,iter]  + F1 
           	}
           	d = d+1    	
    }
    iter = iter + 1
	}
	F1score[[nIP]] = F1score[[nIP]]/5
}


#timestamp score
MSE = list()
for (nIP in 1:2) {
	MSE[[nIP]] = matrix(0, 1212, 6)
    	iter = 1
        for (K in c(5, 10, 25, 50, 75, 100)){
            filename = paste0("Enron_PPE_",nIP,"_",K,".RData")
            load(filename)
             d = 1
for (D in which(missing[,3]==1)) {
    observedt = Enron$edge[[D]]$timestamp - Enron$edge[[D-1]]$timestamp
	 for (o in 1:5) {
        	MSE[[nIP]][d,iter]  = MSE[[nIP]][d,iter]  + median(abs(EnronPPE[[o]]$timepredict[d,10:30] - observedt)/(3600*24))
                }
        	d = d+1
        }
        iter = iter+1
    }
    MSE[[nIP]] = MSE[[nIP]]/5
}

AUC_new = NULL
F1_new = NULL
MAD_new = NULL
for (nIP in 1:2){
	iter = 1
	for (K in c(5, 10, 25, 50, 75, 100)){
			AUC_new = rbind(AUC_new, cbind(mean(AUC[[nIP]][,iter]), nIP, K))
			F1_new = rbind(F1_new, cbind(mean(F1score[[nIP]][,iter]), nIP, K))
			MAD_new = rbind(MAD_new, cbind(mean(MSE[[nIP]][,iter][MSE[[nIP]][,iter]!=Inf]), nIP, K))
			iter = iter + 1
	}	
}
AUC_new = data.frame(AUC = AUC_new[,1], C = as.factor(AUC_new[,2]), K = as.factor(AUC_new[,3]))
F1_new = data.frame(F1 = F1_new[,1], C = as.factor(F1_new[,2]), K = as.factor(F1_new[,3]))
MAD_new = data.frame(MAD = MAD_new[,1], C = as.factor(MAD_new[,2]), K = as.factor(MAD_new[,3]))



#plot
f = list()
library(ggplot2)
library(gridExtra)

f[[1]] <- ggplot(AUC_new, aes(K, AUC, col = C))+geom_line(aes(group = C))+geom_point()
f[[2]] <- ggplot(F1_new, aes(K, F1, col = C))+geom_line(aes(group = C))+geom_point()
f[[3]] <- ggplot(MAD_new, aes(K, MAD, col = C))+geom_line(aes(group = C))+geom_point()

marrangeGrob(f[1:3], nrow = 1, ncol = 3, top = NULL)


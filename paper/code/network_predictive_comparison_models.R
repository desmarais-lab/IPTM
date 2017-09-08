# read in statistics matrix data
load("netstats.RData")
# array: (documents) X (senders) X (receivers) X (statistics)
# names of these statistics?


# For each model below... The model will be fit up to document d, then
## 1,000 values will be simulated conditional on the input values from document d+1

# Multinomial logit predicting the sender
## One covariate for each network statistic
## Covariate value for statistic k for sender i in document d is the sum over 
## the row i for statistic matrix k for document d

# Logit predicting who is added to the e-mail given the sender
## One covariate value for each network statistic
## Covariate value for statistic k for recipient j in document d, sent by i is 
## the ij element of statistic matrix k for document d

# exponential regression model predicting the ime increment since the last e-mail
## One covariate value for each network statistic
## Covariate value for statistic k for document d is the average statistic 
## value taken over the sender/recipient elements of the statistic matrix k for document d


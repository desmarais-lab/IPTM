library(fields)
library(ergm)
library(anytime)
library(lubridate)
load('Darenew.RData')
sender = sapply(1:2246, function(x) Dare$edge[[x+1]]$sender)
timeinc = sapply(1:2246, function(x) Dare$edge[[x+1]]$unixtime - Dare$edge[[x]]$unixtime)/3600
weeks = sapply(1:2246, function(x) wday(anytime(Dare$edge[[x]]$unixtime,  tz = "America/New_York"), label = TRUE) )
day1 = sapply(1:2246, function(x) hour(ymd_hms(anytime(Dare$edge[[x]]$unixtime, tz = "America/New_York"))))
day = cut(day1, c(0,6,12,18,24), c("0-6", "6-12", "12-18", "18-24"))
day[day1==0] = "0-6"
sender = sender[-which(timeinc==0)]
weeks = weeks[-which(timeinc==0)]
day = day[-which(timeinc==0)]
timeinc = timeinc[-which(timeinc==0)]

par(mfrow = c(1,2))
boxplot(log(timeinc)~weeks, names = c("Sun", "Mon", "Tue", "Wed", "Thur", "Fri", "Sat"))
boxplot(log(timeinc)~day)

Daresurv = data.frame(Y = log(timeinc), sender = as.factor(sender), W = weeks, D = day )

attach(Daresurv)
agam0 <- lm(Daresurv$Y~Daresurv$sender+Daresurv$D+Daresurv$W-1)
norm <- glm(Daresurv$Y~Daresurv$sender+Daresurv$D+Daresurv$W-1, family = gaussian(link = "identity"),control = glm.control(maxit = 100), start = as.numeric(coef(agam0)))
normpred <- simulate(norm, type = "response", se =T)
par(mfrow = c(2,2))
plot(norm) 

#using indicators
Daresurv$W = as.numeric(Daresurv$W %in% c("Sat", "Sun"))
Daresurv$D = as.numeric(Daresurv$D %in% c("12-18", "18-24"))

agam0 <- lm(Daresurv$Y~Daresurv$sender+Daresurv$D+Daresurv$W-1)
norm <- glm(Daresurv$Y~Daresurv$sender+Daresurv$D+Daresurv$W-1, family = gaussian(link = "identity"),control = glm.control(maxit = 100), start = as.numeric(coef(agam0)))
normpred <- simulate(norm, type = "response", se =T)
par(mfrow = c(2,2))
plot(norm) 

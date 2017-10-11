setwd('/Users/bomin8319/Desktop/IPTM/full')
library(fields)
library(ergm)
library(anytime)
load('Darenew.RData')
# 762 - 
Dare$edge = Dare$edge[-which(sapply(Dare$text, function(d){length(d)})==0)]
Dare$text = Dare$text[-which(sapply(Dare$text, function(d){length(d)})==0)]

wday = function(x) {
	as.POSIXlt(x)$wday
}
timeinc = sapply(1:2209, function(x) Dare$edge[[x+1]]$unixtime - Dare$edge[[x]]$unixtime)
weeks = sapply(1:2209, function(x) wday(anytime(Dare$edge[[x]]$unixtime,  tz = "America/New_York")) )
day = sapply(1:2209, function(x) hour(ymd_hms(anytime(Dare$edge[[x]]$unixtime, tz = "America/New_York"))))
day = cut(day, c(0,6,12,18,24), c("0-6", "6-12", "12-18", "18-24"))

weeks = weeks[-which(timeinc==0)]
day = day[-which(timeinc==0)]
timeinc = timeinc[-which(timeinc==0)]


par(mfrow = c(1,2))
boxplot(log(timeinc)~weeks, names = c("Sun", "Mon", "Tue", "Wed", "Thur", "Fri", "Sat"))
boxplot(log(timeinc)~day)


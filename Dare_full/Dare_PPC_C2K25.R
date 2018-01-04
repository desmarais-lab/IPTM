library(IPTM)
load('Darenew.RData')
attach(Dare)
Dare$node = 1:nrow(Dare$node)
Dare$edge = lapply(Dare$edge, function(x){x[1:3]})


nIP = 2
K = 25
#PPC generate
for (i in 1:5) {
	filename = paste0("Dare_full_",2,"_",25,"_ver",i,".RData")
	load(filename)
	PPC = IPTM.PPC(Out = 20, edge = Dare$edge, node = Dare$node, textlist = Dare$text, vocab = Dare$vocab,
                   nIP = nIP, K = K, netstat = c("dyadic", "degree", "triadic"), timestat = c("dayofweek", "timeofday"),
                   inference = Daretest, timeunit = 3600, tz = "America/New_York")
	filename = paste0("PPC", i, "_nIP2_K25_Dare.RData")	
	save(PPC, file = filename)	
}

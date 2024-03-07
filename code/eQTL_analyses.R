source("github/code/functions.R")


#########QTL ANALYSIS########################
run_eqtl("data/BayxSha_SA.yaml","SA")
run_eqtl("data/BayxSha_SW.yaml","SW")
run_eqtl("data/BayxSha_delta.yaml","delta")

#reloaded for analysis
load("data/lodmatrix_delta.RData")
out_delta<-out
load("data/map_delta.RData")
map_delta<-map

load("data/lodmatrix_SA.RData")
out_SA<-out
load("data/map_SA.RData")
map_SA<-map

load("data/lodmatrix_SW.RData")
out_SW<-out
load("data/map_SW.RData")
map_SW<-map

#example figure of what the qtl look like for a given transcript
transcript<-"At1g02100"
par(mar=c(4.1, 4.1, 1.6, 1.1),oma=c(0,0,0,0))
plot(x = out_SA,map = map_SA, lodcolumn = transcript,
     col="forestgreen",xlab="Chromosome",ylab="LOD")
plot(x= out_delta,map=map_delta,lodcolumn =transcript,add = TRUE,col="red" )
plot(x= out_SW,map=map_SW,lodcolumn =transcript,add = TRUE , col = "slateblue")
legend("topright", lwd=2,col=c("forestgreen","slateblue","red","black"),
       c("SA","SW","delta","transcript location"),
       lty=c("solid", "solid", "solid", "dashed"))
title("AT1G02100")
abline(v=2,lwd=3,lty="dashed")


peaks_delta<-fread("data/peaks_delta.txt")
peaks_delta[lodcolumn==transcript]
peaks_sa<-fread("data/peaks_SA.txt")
peaks_sa[lodcolumn==transcript]
peaks_sw<-fread("data/peaks_sw.txt")
peaks_sw[lodcolumn==transcript]

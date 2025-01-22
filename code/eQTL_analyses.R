source("github/code/functions.R")
library(qtl2)

#########QTL ANALYSIS########################
run_eqtl("data/BayxSha_SA.yaml","SA")
run_eqtl("data/BayxSha_SW.yaml","SW")
run_eqtl("data/BayxSha_delta.yaml","delta")

#reloaded for analysis
load("data/eqtl_extra_output/lodmatrix_delta.RData")
out_delta<-out
load("data/eqtl_extra_output/map_delta.RData")
map_delta<-map

load("data/eqtl_extra_output/lodmatrix_SA.RData")
out_SA<-out
load("data/eqtl_extra_output/map_SA.RData")
map_SA<-map

load("data/eqtl_extra_output/lodmatrix_SW.RData")
out_SW<-out
load("data/eqtl_extra_output/map_SW.RData")
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
##########################################################################
transcripts <- c("At3g26830", "At2g14610", "At2g40750", "At2g26020")
titles <- c("PAD3 (AT3G26830)", "PR1 (AT2G14610)", "WRKY54 (AT2G40750)", "PDF1.2b (AT2G26020)")
#gmap[txname=="At2g26020"]
gene_locs<-c(102+71+33,102+19,102+61,102+41)
# Set up the plotting area for a 2x2 grid
par(mfrow=c(2,2), mar=c(2, 2, 1, 1), oma=c(2, 2, 1, 1),cex=0.6)

# Loop through each transcript and create the plots
for (i in seq_along(transcripts)) {
  transcript<-transcripts[i]
  # Calculate the maximum LOD score for the current transcript across all datasets
  max_lod <- max(c(out_SA[, transcript], out_delta[, transcript], out_SW[, transcript]))
  
  # Create the plots with appropriate ylim
  plot(x = out_SA, map = map_SA, lodcolumn = transcript,
       col="#377EB8", xlab="", ylab="", ylim=c(0, max_lod),cex.axis=0.6,lwd=1)
  plot(x = out_delta, map = map_delta, lodcolumn = transcript,
       add = TRUE, col="#E6550D", ylim=c(0, max_lod),cex.axis=0.6,lwd=1)
  plot(x = out_SW, map = map_SW, lodcolumn = transcript,
       add = TRUE, col = "#4DAF4A", ylim=c(0, max_lod),cex.axis=0.6,lwd=1)
  abline(v=gene_locs[i],lwd=2,lty="dashed")
  title(titles[i],cex.main = 1)
}
mtext("Chromosome", side=1, line=-0.5, outer=TRUE,cex=0.6)
mtext("LOD", side=2, line=-0.5, outer=TRUE,cex = 0.6)
##########################################################################
























peaks_delta<-fread("data/peaks_delta.txt")
peaks_delta[lodcolumn==transcript]
peaks_sa<-fread("data/peaks_SA.txt")
peaks_sa[lodcolumn==transcript]
peaks_sw<-fread("data/peaks_sw.txt")
peaks_sw[lodcolumn==transcript]


SAgene_qtls<-qtls[txname %in% c(transcripts)]
ggplot(SAgene_qtls,aes(x=trt,y=rsq,color=type))+
  geom_point()+
  facet_wrap(~txname)+
  scale_fill_grey(start = 0.5,end=0.8)


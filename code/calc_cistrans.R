library(data.table)
library(ggplot2)
library(stringr)

#get peak tables from qtl scripts
deltapeaks<-fread("data/eQTL_output/peaks_delta.txt")
sapeaks<-fread("data/eQTL_output/peaks_SA.txt")
swpeaks<-fread("data/eQTL_output/peaks_SW.txt")
#add treatment group
deltapeaks$trt<-"delta"
sapeaks$trt<-"SA"
swpeaks$trt<-"SW"
peakdf<-rbind(sapeaks,swpeaks,deltapeaks)
setnames(peakdf,old = c("lodcolumn","chr","pos"),new = c("txname","qtl_chr","qtl_pos"))

igmap<-fread("data/imputed_gmap_final.txt")
igmap<-igmap[,.(marker,chr,pos)]
setnames(igmap,names(igmap),c("txname","tx_chr","tx_pos"))
#92,019 qtl across grps, 3,290 dont have useful names, removed, now 88,729 qtls
peaks2<-peakdf[which(peakdf$txname %in% igmap$txname)]
dropped_qtls<-peakdf[!which(peakdf$txname %in% igmap$txname)]
#now i have a table with all the info i need :)
qtls<-merge(peaks,igmap,by="txname",all.x=TRUE)

#label cis and trans, 5cM range for cis
qtls[, type := ifelse(qtl_chr == tx_chr & abs(qtl_pos - tx_pos) < 10, "cis", "trans")]

#hell yeah hell yeah
fwrite(qtls,"data/qtls_notcollapsed_March18_10cM.txt")






both_ids<-collapsed_qtls[,.("cis" %in% type & "trans" %in% type),by= name]
both_ids<-both_ids[V1 == TRUE]
subset<-collapsed_qtls[name %in% both_ids$name]

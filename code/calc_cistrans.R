library(data.table)
library(ggplot2)
library(stringr)

#get peak tables from qtl scripts
deltapeaks<-fread("data/peaks_delta.txt")
sapeaks<-fread("data/peaks_SA.txt")
swpeaks<-fread("data/peaks_SW.txt")
#add treatment group
deltapeaks$trt<-"delta"
sapeaks$trt<-"SA"
swpeaks$trt<-"SW"
peakdf<-rbind(sapeaks,swpeaks,deltapeaks)
setnames(peakdf,old = c("lodcolumn","chr","pos"),new = c("txname","qtl_chr","qtl_pos"))

igmap<-fread("data/imputed_gmap.txt")
igmap<-igmap[,.(marker,chr,pos)]
setnames(igmap,names(igmap),c("txname","tx_chr","tx_pos"))
#92,019 qtl across grps, 3,290 dont have useful names, removed, now 88,729 qtls
peaks<-peakdf[which(peakdf$txname %in% igmap$marker)]

#now i have a table with all the info i need :)
qtls<-merge(peaks,igmap,by="txname",all.x=TRUE)

#label cis and trans, 5cM range for cis
qtls[, type := ifelse(qtl_chr == tx_chr & abs(qtl_pos - tx_pos) < 5, "cis", "trans")]

#hell yeah hell yeah
fwrite(qtls,"data/qtls_notcollapsed_March7.txt")






both_ids<-collapsed_qtls[,.("cis" %in% type & "trans" %in% type),by= name]
both_ids<-both_ids[V1 == TRUE]
subset<-collapsed_qtls[name %in% both_ids$name]

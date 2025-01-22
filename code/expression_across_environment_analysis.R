library(ggplot2)
library(data.table)
library(gridExtra)
library(tidyr)
library(pbapply)
library(ggh4x)
library(forcats)
library(ggforce)

dt<-fread("data/qtls_10cMwindow_fixed_overlap_5.5cm_march18.txt")
check<-dt[rsq>0.5]
highs<-dt[,.SD[rsq==max(rsq)],by = .(trtgroup,type)]
colnames(highs)[3]<-'txname'
#had to pick 2 different transcripts for SA cis and SW cis so they arent the same
#dt[trtgroup == "SW"&type == "cis",.SD, .SDcols=names(dt)][order(-rsq)]
#sa cisAt2g21640
#silwet cis At2g20670
highs
highs<-highs[!((trtgroup %in% c("SA","SW")) & (type == "cis")),]
high1<-dt[type=="cis"&trtgroup == "SA"& txname=="At2g21640",]
high2<-dt[type=="cis"&trtgroup == "SW"& txname=="At2g20670",]
highs<-rbind(highs,high1,high2)
#get transcript abundance data
pheno<-fread("data/rils_pheno.txt")
#read in map of genotypes for the rils, make names compatible 
geno<-fread("data/geno_editednames.txt")
#read in gmap
igmap<-fread("data/imputed_gmap_final.txt")
gmap<-fread("data/original/BayxSha_gmap.csv")
gmap<-gmap[,marker := gsub("\\..*","",marker)]
#if transcript does not have genotype information, 
#just get the genotype information of the nearest marker
highs[, closest_marker_qtl := gmap[.SD, on = .(chr = qtl_chr, pos = qtl_pos), roll = "nearest", x.marker]]

# Define the function
find_closest_marker <- function(highs, igmap, gmap) {
  # Add an empty column for the closest marker based on tx_pos
  highs[, closest_marker_tx := NA_character_]
  for (i in 1:nrow(highs)) {
    # Extract txname for the current row
    current_txname <- highs$txname[i]
    # Find the corresponding tx_pos from igmap
    tx_pos <- igmap[marker == current_txname]$pos
    # Extract tx_chr for the current row to filter gmap
    tx_chr <- igmap[marker == current_txname]$chr
    # Filter gmap for the current chromosome
    gmap_chr <- gmap[chr == tx_chr]
    # If gmap_chr is not empty, proceed
    if (nrow(gmap_chr) > 0) {
      # Calculate the absolute difference between tx_pos and each pos in gmap_chr
      gmap_chr[, diff := abs(tx_pos - pos)]
      # Find the row with the minimum difference
      closest_row <- gmap_chr[which.min(diff)]
      # Update the closestmarker_tx column in highs for the current row
      highs$closest_marker_tx[i] <- closest_row$marker
    }
  }
  
  return(highs)
}

highs_updated <- find_closest_marker(highs, igmap, gmap)

#function to merge the phenotype and genotype information based on the qtl data
merge_pheno_geno<-function(qtldt){
  pbapply(qtldt,1,FUN = function(row){
    tx <- as.character(row[3])
    trtment<-row[4]
    nearestmarker_tx <-as.character(row[10])
    group <- as.character(row[1])
    regtype <- as.character(row[2])
    qtlpos<- as.character(row[7])
    qtlmarker<-as.character(row[9])
    gene_tx <- pheno[, .(ID,treatment, transcript = get(tx))]
    gene_geno <- geno[, .(ID, genotype = get(qtlmarker))]
    genedt <- merge(gene_tx, gene_geno)
    genedt$txname<-rep(tx)
    genedt$closestmarkertotx<-rep(nearestmarker_tx)
    genedt$trtgroup<-rep(group)
    genedt$type<-rep(regtype)
    genedt$qtlmarker<-rep(qtlmarker)
    genedt$qtlpos<-rep(qtlpos)
    genedt$chr<-rep(as.character(row[8]))
    return(genedt)
  })
}
highsgeneinfolist<-merge_pheno_geno(highs_updated)
highsgeneinfo<-rbindlist(highsgeneinfolist)

highsgeneinfo[genotype==-1]<-NA
highsgeneinfo<-na.omit(highsgeneinfo)

##########################################################################################################
#environment single gene plots (FIGURE 6)
##########################################################################################################
unique_combinations<-unique(highsgeneinfo[,.(trtgroup,txname,type)])
highsnodelta<-highsgeneinfo[treatment!='delta']
highsnodelta<-highsnodelta[order(match(trtgroup,c("SA","SW","delta","SA SW","SA delta","SW delta","SA SW delta")))]

ggplot(highsnodelta, aes(x=treatment, y=transcript, color=genotype)) +
  geom_point(aes(group=interaction(ID,txname)),alpha = 0.3) +
  facet_grid2(fct_relevel(trtgroup,"SA","SW","delta","SA delta","SW delta","SA SW","SA SW delta")~type,scales = 'free_y',independent = 'y')+
  theme_bw()+
  scale_color_manual(values = c("#56B4E9","#CC79A7"))+
  theme(
    legend.position = "none",
    text = element_text(size = 6),
    axis.text = element_text(size = 6))+
  ylab("Transcript abundance")+
  xlab("Treatment")+
  stat_summary(
    data = subset(highsnodelta, genotype == "C"),
    aes(group=genotype), 
    fun=mean, geom="line",
    size=0.8, 
    color="#56B4E9") +
  stat_summary(
    data = subset(highsnodelta, genotype == "L"),
    aes(group=genotype), 
    fun=mean, geom="line", 
    size=0.8, 
    color="#CC79A7")


selected_columns <- highsnodelta[, .(txname, type, trtgroup)]
unique_combinations <- unique(selected_columns)


##single gene expression over environment for SA specific genes###
genes<-c("ID","treatment","At3g26830",
         "At2g14610","At2g40750","At2g26020")
rils<-fread("data/rils_pheno.txt")
rils<-rils[,.SD,.SDcols = genes]
parents<-fread("data/parents_pheno_corrected.txt")
parents<-parents[,lapply(.SD,mean),by = .(treatment,ID)]
parents$rep<-NULL
parents<-parents[,.SD,.SDcols = genes]
genesdt<-rbind(rils,parents)
genesdt<-genesdt[treatment != "delta"]

geneslong<-melt(genesdt,id = c("ID","treatment"),variable.name = "transcript")
#FIGURE 3A
ggplot(geneslong,aes(x = treatment, y = value, group = ID))+
  geom_line(alpha = 0.6, linewidth = 0.1)+
  facet_wrap(~transcript, scales = "free_y",
             labeller = as_labeller(c(At3g26830 = "PAD3 (AT3G26830)",
                          At2g14610 = "PR1 (AT2G14610)",
                          At2g40750 = "WRKY54 (AT2G40750)",
                          At2g26020 = "PDF1.2b (AT2G26020)")))+
  theme_bw()+
  theme(panel.background = element_blank(),
        strip.text = element_text(size=7),
        text = element_text(size = 7))+
        ylab("Transcript Abundance")+
        xlab("Treatment")+
  geom_line(data=geneslong[ID=="Bay"],aes(x=treatment,y=value,group=ID),size=1,colour="#56B4E9")+
  geom_line(data=geneslong[ID=="Sha"],aes(x=treatment,y=value,group=ID),size=1,colour="#CC79A7")


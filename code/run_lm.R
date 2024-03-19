library(data.table)
library(parallel)
source("github/plasticity_eQTL/code/functions.R")

#read in map of genotypes for the rils
geno<-fread("data/original/geno.csv")

sapheno<-fread('data/original/reformatted/SApheno_March5.txt')
sapheno$treatment<-"SA"
swpheno<-fread('data/original/reformatted/SWpheno_March5.txt')
swpheno$treatment<-"SW"
deltapheno<-fread('data/original/reformatted/deltapheno_March5.txt')
setnames(deltapheno,"rils","ID")
deltapheno$treatment<-"delta"
pheno<-rbind(sapheno,swpheno,deltapheno)
fwrite(pheno,"rils_pheno.txt")

#adding a column for closest marker in gmap to qtl position
gmap<-fread("data/original/BayxSha_gmap.csv")
setkey(gmap, chr, pos)
# Performing the rolling join
qtls[, closest_marker := gmap[.SD, on = .(chr = qtl_chr, pos = qtl_pos), roll = "nearest", x.marker]]

#Use the lm() function in R to fit a linear model to the data, 
#where the dependent variable is the phenotype of interest and the 
#independent variable is the genotype at the QTL. 
#You can extract the genotype data from the qtl_results data table.

#Once you have fit the linear model, you can extract the effect size 
#of the QTL by examining the coefficient of the qtl_genotypes variable 
#in the model. The coefficient represents the change in the phenotype for 
#each additional copy of the QTL allele.

#edit geno names to match the rest of the input files
colnames(geno)<-gsub("\\.[0-9]+", "",colnames(geno))
gmap$marker<-gsub("\\.[0-9]+", "",gmap$marker)
#making a table for linear model

fwrite(geno,"geno_editednames.txt")
fwrite(gmap,"gmap_editiednames.txt")



clust<-makeCluster(6)
clusterExport(clust,c("qtls","pheno","geno"))
clusterEvalQ(clust,library(data.table)) 
start <- Sys.time()
effectoutput<-parApply(cl = clust, qtls,1, function(x){
  dt<-data.table(transcript=pheno[treatment == x[8],get(x[1])],genotype=geno[,get(x[12])])
  lmout<-(lm(transcript ~ genotype,dt, na.action = "na.exclude"))
  summary(lmout)
})
print( Sys.time() - start )
stopCluster(clust)

saveRDS(effectoutput,"lm_summary")


effectsizes_dt<-readRDS("lm_summary")

rsqs<-unlist(lapply(effectsizes_dt,function(x){
  x[[8]]
}))

dt<-fread("data/qtls_notcollapsed_March18_10cM.txt")
dt<-data.table(dt,rsq=rsqs)
fwrite(dt,"qtls_10cMwindow_notcollapsed_rsq_March18.txt")

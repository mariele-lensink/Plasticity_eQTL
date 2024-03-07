#load packages
library(data.table)
library(readxl)
source("github/code/functions.R")

#read in conversion chart and add NAs in empty spaces
tairnames<-fread("data/transcriptome_conversion.csv",header = T)

#Reformat the replicate files
SApheno<-process_replicates("data/ALL_211_RILs/SA1.xls","data/ALL_211_RILs/SA2.xls")
###had to find which ril is not in carolines files and remove it in mine (RIL 417), not in genotype file either###
SApheno[146,1]
SApheno<-SApheno[-146,]#remove the row for RIL 417
sapheno<-convert_transcript_names(SApheno,tairnames)
setorder(sapheno,ID)
fwrite(sapheno,file = "data/SApheno_March5.txt")#write file

SWpheno<-process_replicates("data/ALL_211_RILs/sw1.xls","data/ALL_211_RILs/sw2.xls")
SWpheno[146,1]
SWpheno<-SWpheno[-146,]
swpheno<-convert_transcript_names(SWpheno,tairnames)
setorder(swpheno,ID)
fwrite(swpheno,file = "data/SWpheno_March5.txt")

#Make the delta phenotype
rils<-sapheno[,1]
sam<-as.matrix(sapheno)
swm<-as.matrix(swpheno)
sam<-sam[,-1]
swm<-swm[,-1]
#calculate difference divided by mean
delta <- as.data.table((sam-swm)/((sam+swm)/2))
delta<-cbind(rils,delta)

fwrite(delta,file = "data/deltapheno_March5.txt")


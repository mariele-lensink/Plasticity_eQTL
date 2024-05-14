#load packages
library(data.table)
library(readxl)
source("github/plasticity_eQTL/code/functions.R")
#read in conversion chart and add NAs in empty spaces
tairnames<-fread("data/original/transcriptome_conversion.csv",header = T)

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
################################################################################################
#without averaging replicates
SA<-as.data.table(process_replicates2("data/original/ALL_211_RILs/SA1.xls","data/original/ALL_211_RILs/SA2.xls"))
SA[,treatment:= "SA"]
setcolorder(SA, c("treatment","rep", "ID", setdiff(names(SA), c("treatment","rep", "ID"))))
SW<-as.data.table(process_replicates2("data/original/ALL_211_RILs/sw1.xls","data/original/ALL_211_RILs/sw2.xls"))
SW[,treatment:="SW"]
setcolorder(SW, c("treatment","rep", "ID", setdiff(names(SW), c("treatment","rep", "ID"))))

pheno<-rbind(SA,SW)
pheno<-convert_transcript_names(pheno,tairnames)
fwrite(pheno,"data/rils_pheno_w_reps.txt")
#######################################################################################################
#####this is WRONG, go back and redo ALL analyses with this data!!!
#formatting parent lines, not averaging replicates
parents<-read_excel("data/original/ALL_211_RILs/Parental CV.xlsx",col_names = FALSE)
parents<-data.table(t(parents[-c(4:6),]))
setnames(parents, as.character(unlist(parents[1,])))
parents <- parents[-1, ]
setnames(parents, old = names(parents)[1:3], new = c("ID","treatment","rep"))
parents[, ID := substr(ID, 1, 3)]
parents<-convert_transcript_names(parents,tairnames)
fwrite(parents,"data/parents_pheno-w_reps.txt")
dt<-fread("data/parents_pheno-w_reps.txt")
dt2<-dt[,lapply(.SD,mean),by = .(ID,treatment), .SDcols = !c("ID", "treatment", "rep")]
fwrite(dt2,"data/parents_pheno.txt")
#######################################################################################################
#CORRECT parental data!
file_names<-c("SA1.xls", "SA2.xls", "sw1.xls", "sw2.xls")
parentdt<-rbindlist(lapply(file_names, process_parents))
parents<-convert_transcript_names(parentdt,tairnames)
setnames(parents,old = "genotype", new = "ID")
fwrite(parents, "data/parents_pheno_corrected.txt")

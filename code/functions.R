#load packages
library(data.table)
library(readxl)
library(qtl2)

get_reformatted_matrix <- function(data) {
  data <- data.table(t(data[-c(1:4), -1]))
  data <- sapply(data, as.numeric)
  colnames(data)[] <- ""
  matrixData <- as.matrix(data)
  return(matrix = matrixData)
}

process_replicates<-function(file1,file2){
  df1 <- read_excel(file1)
  df2 <- read_excel(file2)
  m1<-get_reformatted_matrix(df1)
  m2<-get_reformatted_matrix(df2)
  #get average value between salicylic acid reps
  pheno<-data.table((m1+m2)/2)
  #ad ril ID back to table
  temp<-data.frame(df1)
  txnames<-c(temp[-c(1:4),1])
  #add transcript names back
  colnames(pheno)<-txnames
  #add an ID column with ril numbers
  pheno<-data.frame(ID = as.numeric(c(colnames(df1[-1]))),pheno)
  return(pheno)
}

convert_transcript_names<-function(dt,c_chart){
  colnames(dt)<-sub("^X", "",colnames(dt))
  c_vector<-setNames(as.character(c_chart2[[2]]),as.character(c_chart2[[1]]))
  names(dt) <- ifelse(names(dt) %in% names(c_vector), c_vector[names(dt)], names(dt))
  return(dt)
}

run_eqtl<-function(yaml,tag){
  bayxsha<-read_cross2(yaml,quiet = FALSE)
  #calculate conditional genotype probabilities given the marker data at each putative QTL position
  ##insert pseudomarkers into the genetic map
  map<-insert_pseudomarkers(bayxsha$gmap,step=0)
  #calculate genotype probabilities
  pr<- calc_genoprob(bayxsha,map,error_prob = 0.002)
  ####run the genome scan####
  #scan1() takes as input the genotype probabilities, a matrix of phenotypes, 
  #uses they Haleyy-knott regression 
  out <- scan1(pr, bayxsha$pheno)#output is LOD scores. format is positions x phenotypes
  #find peaks
  peaks<-find_peaks(out, map, threshold=2, drop=1.5)
  fwrite(peaks,file = paste0("data/peaks_",tag,".txt",sep = ''))
  save(map,file =paste0("data/map_",tag,".RData",sep = ''))
  print("finished saving map")
  save(out,file =paste0("data/lodmatrix_",tag,".RData",sep = ''))
  print("finished saving lod (out)")
}


#load packages
library(data.table)
library(readxl)
library(qtl2)
library(car)

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

process_replicates2<-function(file1,file2){
  df1 <- read_excel(file1)
  df2 <- read_excel(file2)
  m1<-get_reformatted_matrix(df1)
  m2<-get_reformatted_matrix(df2)
  dt1 <- data.table(m1)
  dt2 <- data.table(m2)
  ril_numbers <- as.numeric(colnames(df1)[-1])
  temp<-data.frame(df1)
  txnames <- as.character(temp[-c(1:4),1 ])
  dt1 <- data.table(rep = 1, ID = ril_numbers, m1)
  dt2 <- data.table(rep = 2, ID = ril_numbers, m2)
  appended_data <- rbind(dt1, dt2)
  setnames(appended_data, c("rep", "ID", txnames))
  return(appended_data)
}

convert_transcript_names<-function(dt,c_chart){
  colnames(dt)<-sub("^X", "",colnames(dt))
  c_vector<-setNames(as.character(c_chart2[[2]]),as.character(c_chart2[[1]]))
  names(dt) <- ifelse(names(dt) %in% names(c_vector), c_vector[names(dt)], names(dt))
  return(dt)
}

process_parents <- function(file_name) {
  data_path <- paste0("data/original/ALL_211_RILs/", file_name)
  dt <- t(as.data.table(read_excel(data_path, sheet = 2)))
  
  # Process rows and columns
  ids <- sapply(rownames(dt)[-1], substr, 1, 3)
  repnum <- as.character(dt[2, 1])
  trt <- as.character(dt[2, 2])
  dt <- as.data.frame(dt[, -c(1:5)])
  colnames(dt) <- dt[1, ]
  dt <- dt[-1, ]
  
  # Create final data table
  data.table(genotype = ids, treatment = trt, rep = repnum, dt)
}

run_parallel_anova<-function(dt){
  clust<-makeCluster(4)
  dtName <- deparse(substitute(dt))
  print(dtName)
  clusterExport(clust,c(dtName))
  clusterEvalQ(clust, {
    library(data.table)
    library(car)
  })
  start <- Sys.time()
  aov_sumsqs<-rbindlist(parApply(cl=clust,dt[,-c(1:3),with = FALSE], 2, function(x){
    #run and anova on each gene according to batch number
    model1 <- lm(x ~ as.factor(treatment) * as.factor(ID) + as.factor(rep), data = dt)
    #make summary table of anova stats. we are looking at the sum of squares 
    #(sum of the sq of the variation), distance between each ind value and the mean (line of best fit) and then summed
    aov<-Anova(model1,type="III")
    tss<-sum(aov[[1]][2:6])
    #table of sum of squares for each gene and its proportion, factor/residual
    proptable<-data.table(treatment=aov[[1]][2]/tss,
                          replicate=aov[[1]][4]/tss,
                          genotype=aov[[1]][3]/tss,
                          trt.geno.int=aov[[1]][5]/tss,
                          residuals=aov[[1]][6]/tss)
    ##im pretty confident this is wrong but also was 
    return(proptable)
  }))  
  return(aov_sumsqs)
  print( Sys.time() - start )
  stopCluster(clust) 
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

order_split <- function(df) {
  ordered <- df[order(df$qtl_chr, df$qtl_pos), ]
  split(ordered, by = 'qtl_chr')
}

centimorgan_slider<-function(positions,window,stepsize,cms) {
  #set initial window ranges before loop 
  rangemin<-1
  rangemax<-rangemin+window-1
  peakcounts<-c()
  while(rangemax <= cms){
    n<-sum(positions >= rangemin & positions < rangemax)
    peakcounts<-c(peakcounts,n)
    
    rangemin<- rangemin+stepsize
    rangemax<-rangemax+stepsize
    
  }
  return(peakcounts)
}

apply_slider_to_processed_list <- function(qtls_processed, window, stepsize, cms) {
  lapply(qtls_processed, function(lst) {
    lapply(lst, function(i) data.table(centimorgan_slider(i[, qtl_pos], window, stepsize, cms)))
  })
}

reformat_for_ggplot <- function(dt_list) {
  treatment_names <- names(dt_list)
  reformatted_list <- list()
  for(trt in treatment_names) {
    dt <- dt_list[[trt]]
    index <- seq_len(nrow(dt[[1]]))
    for(i in seq_along(dt)) {
      dt[[i]] <- data.table(Count = dt[[i]]$V1, Index = index, Treatment = trt, Chromosome = as.character(i))
    }
    combined_dt <- rbindlist(dt)
    reformatted_list[[trt]] <- combined_dt
  }
  all_combined <- rbindlist(reformatted_list)
  all_combined[, Treatment := factor(Treatment, levels = treatment_names)]
  return(all_combined)
}

generate_random_qtl_data_table <- function(qtl_count, total_cms, permutations = 1000) {
  random_matrix <- matrix(runif(qtl_count * permutations, 0, total_cms), nrow = qtl_count, ncol = permutations)
  random_qtls <- as.data.table(random_matrix)
  return(random_qtls)
}



library(data.table)
library(ggplot2)
library(dbscan)


dt<-fread("qtls_10cMwindow_notcollapsed_rsq_March18.txt")
dt<-dt[,.(txname,trt,type,rsq,lod,qtl_pos,qtl_chr)]

dtlist<-as.list(split(dt,by = c("txname","qtl_chr"), keep.by = T))

condensePairs <- function(dt) {
  # Add a row identifier column
  if (nrow(dt)>1){
    pos_mt<-as.matrix(dt$qtl_pos)
    dt[, group := dbscan(pos_mt, eps = 5.5, minPts = 1)$cluster]
    # Aggregating pairs
    condensed<-dt[,.(
      txname = first(txname),
      trt = list(trt),
      type = first(type),
      rsq = mean(rsq),
      lod = mean(lod),
      qtl_pos = mean(qtl_pos),
      qtl_chr = first(qtl_chr)
    ), by = group]
    condensed$group<-NULL
    
    return(condensed)
  }
  return(dt)
}

qtloverlapfix<-rbindlist(lapply(dtlist,condensePairs))
#FIX treatment group order
desired_order <- c("SA", "SW", "delta")
sort_list <- function(lst, order) {
  lst[order(match(lst, order))]
}
qtloverlapfix[, trt := lapply(trt, sort_list, order = desired_order)]
qtloverlapfix[,trtgroup := sapply(trt, function(x) paste(x,collapse = " "))]
# View the result

fwrite(qtloverlapfix,"qtls_10cMwindow_fixed_overlap_5.5cm_march18.txt")




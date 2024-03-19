library(data.table)

gff<-fread("data/original/TAIR10_GFF3_genes.gff")
gff<-gff[V3=="gene" & !V1 %in% c("ChrC","ChrM"),
         .(V1 =as.numeric(gsub("Chr(\\d)","\\1",V1 )),V4,V5,V9=gsub("ID=([0-9A-Z]+);.+","\\1",V9))]
gff[,V9:=chartr("TG","tg",V9)]
gff[,avg := rowMeans(.SD),.SDcols = c("V4","V5")]
gff<-gff[,.(V1,V9,avg)]
setnames(gff,old= names(gff),new = c("chr","marker","avg"))


gmap<-fread("data/original/BayxSha_gmap.csv")
gmap$marker <- sub("\\.[0-9]+", "",gmap$marker)

imputed_gmap<-merge(gff,gmap,by= c("marker","chr"),all.x=TRUE)
imputed_gmap_list<-split(imputed_gmap,imputed_gmap$chr)

interpolate_pos<- function(dt){
  #now we IMPUTE
  #get anchors for relative calcs
  anchor_indices <- which(!is.na(dt$pos))
  
  #interpolate!!
  for (i in seq_along(anchor_indices)){
    if(i<length(anchor_indices)){
      #range defining rows in between a set of anchors
      start_row <- anchor_indices[i]
      end_row <- anchor_indices[i + 1]
      start_pos<-dt$pos[start_row]#first pos in dt within anchor range
      end_pos<-dt$pos[end_row]#last pos of dt within anchor range
      start_avg<-dt$pos[start_row]#same thing w average column
      end_avg<-dt$avg[end_row]
      
      #calc the new pos values using linear interpolation
      #define rows inside range to interpolate, calculate relative distance and scale by position then add pos back
      dt[(start_row+1):(end_row-1), pos := start_pos + ((avg - start_avg) / (end_avg - start_avg)) * (end_pos - start_pos)]
    }
  }
  return(dt)
}
imputed_by_chr<-lapply(imputed_gmap_list,interpolate_pos)


#finishing the ends of the table that did not fit in a range by using the average interval of pos
interpolate_ends<-function(dt){
  avg_pos_interval<-mean(na.omit(diff(na.omit(dt$pos))))
  
  # get first and last known pos row
  firstrow<- min(which(!is.na(dt$pos)))
  lastrow <- max(which(!is.na(dt$pos)))

  #interpolate head of imputed gmap
  num_rows<-firstrow-1
  dt[1:(firstrow-1), pos := dt$pos[firstrow] - seq(num_rows) * avg_pos_interval]
  #same for tail
  num_rows <- dt[,.N] - lastrow
  dt[(lastrow+1):.N, pos := dt$pos[lastrow] + seq(num_rows) * avg_pos_interval]
  #Interpolate head of imputed gmap

  return(dt)
}
imputed_gmap_final2<-rbindlist(lapply(imputed_by_chr,interpolate_ends))

fwrite(imputed_gmap_final2,"data/imputed_gmap_final.txt")

library(data.table)
library(ggplot2)
library(dplyr)
library(ggupset)


qtls<-fread("data/qtls_10cMwindow_notcollapsed_rsq_March18.txt")
#qtls<-fread("data/qtls_10cMwindow_fixed_overlap_5.5cm_march18.txt")
orthogroups_all<-fread("wrkdir_genespace_w_col0/results/Orthogroups.tsv")
sha_bed<-fread("wrkdir_genespace_w_col0/bed/Sha.bed")
bay_bed<-fread("wrkdir_genespace_w_col0/bed/Bay.bed")
combed<-fread("wrkdir_genespace_w_col0/results/combBed.txt")
combed[20:30,]

qtls$txname<-toupper(qtls$txname)
# Perform the non-equi join
mergedt <- qtls[, {
  # For each txname, check if it is a substring of any Athaliana entries
  matched_entries <- orthogroups_all[grep(txname, Athaliana, fixed = TRUE), .(Athaliana, Bay, Sha)]
  if (nrow(matched_entries) > 0) {
    list(matched_entries$Athaliana[1], matched_entries$Bay[1], matched_entries$Sha[1])
  } else {
    list(NA_character_, NA_character_, NA_character_)
  }
}, by = .(txname)]

# Assign the results back to hightrans
setkey(qtls,txname)
setkey(mergedt,txname)
qtls_w_ortho <- merge(qtls, mergedt, by = "txname", all.x = TRUE)
#turn the paralogs into a list so each individual gene can be accessed
qtls_w_ortho[, c("V1", "V2", "V3") := lapply(.SD, function(x) lapply(x, function(y) strsplit(y, ",\\s*")[[1]])),
             .SDcols = c("V1", "V2", "V3")]
setnames(qtls_w_ortho,old=c("V1","V2","V3"), new = c("Col","Bay","Sha"))
qtls_w_ortho[, Col_count := sapply(Col, function(x) if (length(x) == 1 && is.na(x)) 0 else length(x))]
qtls_w_ortho[, Bay_count := sapply(Bay, function(x) if (length(x) == 1 && is.na(x)) 0 else length(x))]
qtls_w_ortho[, Sha_count := sapply(Sha, function(x) if (length(x) == 1 && is.na(x)) 0 else length(x))]

fwrite(qtls_w_ortho,"notcollapsed_qtls_w_orthos.txt")


qtls_w_ortho[, single_copy_gene := (Bay_count == 1 & Sha_count == 1)]
qtls_w_ortho[, equal_number_paralogs := (Bay_count == Sha_count & Bay_count !=0 & Sha_count !=0)]
qtls_w_ortho[, paralog_number_variation := (Bay_count != Sha_count)]


###########checking to see if one of the colombia paralogs is at the position of the qtl###

#first looking at paralog position in CM on igmap
fetch_positions_and_chromosomes <- function(genes, mapping) {
  pos <- sapply(genes, function(g) mapping[txname == g, tx_pos])
  chr <- sapply(genes, function(g) mapping[txname == g, tx_chr])
  list(pos = pos, chr = chr)
}
qtls_w_ortho[, c("Col_Pos", "Col_Chr") := {
  result <- lapply(Col, fetch_positions_and_chromosomes, mapping = igmap)
  list(lapply(result, `[[`, "pos"), lapply(result, `[[`, "chr"))
}]
check_proximity <- function(pos_list, chr_list, qtl_position, qtl_chr) {
  pos_list <- unlist(pos_list)
  chr_list <- unlist(chr_list)
  any(chr_list == qtl_chr & abs(pos_list - qtl_position) <= 10)
}

# Apply this function to create the new paralog_at_qtl column 
qtls_w_ortho[, paralog_at_qtl := mapply( check_proximity,
                                         pos_list = Col_Pos,
                                         chr_list = Col_Chr,
                                         qtl_position = qtl_pos,
                                         qtl_chr = qtl_chr)]


#now determining if the qtl is in a hotspot
qtls_w_ortho[qtl_chr == 2 & abs(qtl_pos-14.4)>=5 | qtl_chr == 5 & abs(qtl_pos-71)<=5,hotspot := TRUE]
qtls_w_ortho[,hotspot := ifelse(is.na(hotspot),FALSE,hotspot)]
qtls_w_ortho[,no_col_info := ifelse(is.na(Col),TRUE,FALSE)]
###no bay or sha########
qtls_w_ortho[Bay_count==0 & Sha_count == 0, No_Bay_or_Sha :=TRUE]
qtls_w_ortho[,No_Bay_or_Sha := ifelse(is.na(No_Bay_or_Sha),FALSE,No_Bay_or_Sha)]
#changing all cis qtl to be TRUE for paralog_at_qtl
qtls_w_ortho[type=="cis",paralog_at_qtl := TRUE]
qtls_w_ortho[type=='cis',no_col_info := FALSE]
#############################################
qtls_w_ortho$combination <- apply(qtls_w_ortho[, c("equal_number_paralogs",
                                                   "hotspot",
                                                   "no_col_info",
                                                   "paralog_at_qtl",
                                                   "paralog_number_variation",
                                                   "single_copy_gene",
                                                   "No_Bay_or_Sha" )], 1, function(x) {
  paste(names(x)[x], collapse = " ")
})
#test to remove qtls with no bay or sha paralog info, removed 1002
qtls_w_ortho1<-qtls_w_ortho[No_Bay_or_Sha == F]

counts <- qtls_w_ortho1 %>%
  group_by(combination, type) %>%
  summarise(count = n(), .groups = "drop")
#reorder the combination group by frequency
total_counts <- counts %>%
  group_by(combination) %>%
  summarise(total_count = sum(count), .groups = "drop") %>%
  arrange(desc(total_count))

qtls_w_ortho1$combination <- factor(qtls_w_ortho1$combination, 
                                   levels = total_counts$combination)
###########################
ggplot(qtls_w_ortho1, aes(x = as.factor(combination))) +
  axis_combmatrix(sep = " ")+
  geom_boxplot(
    aes(
      group = interaction(combination, type),  # Group by combination and type
      y = rsq,
      x = as.factor(combination),
      color = type),
      alpha = 0.3,
      position = position_dodge(width=0.85),
      width = 0.7)+
  geom_text(
    data = counts,
    aes(
      x = as.factor(combination),
      y = max(qtls_w_ortho$rsq, na.rm = TRUE) * 1.05,  
      label = count,
      group = type,
      color = type,
      angle = 55),
      position = position_dodge(width = 0.9),  
      size = 2,
      show.legend = FALSE)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = c(-0.21,0.82),
        text = element_text(size = 6),
        legend.frame = element_rect(linetype = ),
        axis.title = element_text(size = 7))+
  theme_combmatrix(combmatrix.label.text = element_text(size=6),
                  combmatrix.label.extra_spacing = 2,
                  combmatrix.panel.point.size = 2)+
  xlab("Classification according to Gene and eQTL Characteristics")+
  ylab("Effect Size (R-Squared)")+
  scale_color_manual(values = c("trans" = "navy", "cis" = "cornflowerblue"))
#####################################################
# Function to create a new data.table with transformed list columns
create_transformed_dt <- function(dt) {
  new_dt <- copy(dt)
  list_cols <- names(new_dt)[sapply(new_dt, is.list)]
  new_dt[, (list_cols) := lapply(.SD, function(x) {
    if (is.list(x)) sapply(x, function(y) paste(y, collapse = ","))
    else x
  }), .SDcols = list_cols]
  
  return(new_dt)
}
create_transformed_dt(qtls_w_ortho1)
new_dt_high_trans<-new_dt[rsq > 0.43 & type == 'trans',]
fwrite(new_dt_high_trans,"high_effect_trans_table.txt")



 
##filtered for high effect size
counts <- qtls_w_ortho1[rsq>0.43 & type == 'trans'] %>%
  group_by(combination) %>%
  summarise(count = n(), .groups = "drop")
#reorder the combination group by frequency
total_counts <- counts %>%
  group_by(combination) %>%
  summarise(total_count = sum(count), .groups = "drop") %>%
  arrange(desc(total_count))
##FIGURE 7##
ggplot(qtls_w_ortho1[rsq>0.43 & type =='trans'], aes(x = as.factor(combination))) +
  axis_combmatrix(sep = " ")+
  geom_point(aes(
      y = rsq,
      x = as.factor(combination)))+
  geom_violin(aes(y=rsq,x=as.factor(combination)))+
  geom_text(
    data = counts,
    aes(
      x = as.factor(combination),
      y = max(qtls_w_ortho$rsq, na.rm = TRUE) * 1.05,  
      label = count,
      angle = 55),
    size = 3,
    show.legend = FALSE)+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  xlab("Classification according to Gene and eQTL Characteristics")+
  ylab("Effect Size (R-Squared)")
  #scale_color_manual(values = c("trans" = "navy") #"cis" = "cornflowerblue"))



##TTESTS#
# Perform t-test for cis QTLs
cis_t_test <- t.test(
  rsq ~ equal_number_paralogs == TRUE,
  data = qtls_w_ortho[type == "cis" & No_Bay_or_Sha == FALSE,]
)
# Perform t-test for trans QTLs
trans_t_test <- t.test(
  rsq ~ equal_number_paralogs == TRUE,
  data = qtls_w_ortho[type == "trans" & No_Bay_or_Sha == FALSE,]
)

###CIS DELETION ANALYSIS###
qtls_w_ortho[type=='cis'& equal_number_paralogs==FALSE &  No_Bay_or_Sha == FALSE,.N]/qtls_w_ortho[type=='cis' &  No_Bay_or_Sha == FALSE,.N]

##trans-eqtl explained by paralog at eqtl?
qtls_w_ortho[type=='trans'& paralog_at_qtl==T&No_Bay_or_Sha==F,.N]/qtls_w_ortho[type=='trans'&No_Bay_or_Sha==F,.N]
qtls_w_ortho[type=='trans'& paralog_at_qtl==T&No_Bay_or_Sha==F&paralog_number_variation==T,.N]/qtls_w_ortho[type=='trans'&No_Bay_or_Sha==F,.N]

t.test(
  rsq ~ paralog_at_qtl == TRUE,
  data=qtls_w_ortho[type=='trans' & No_Bay_or_Sha==F,])


#############For single copy genes#########################
t.test(
  rsq ~ paralog_at_qtl == TRUE,
  data=qtls_w_ortho[type=='trans' & No_Bay_or_Sha==F & single_copy_gene==T,])

t.test(
  rsq ~ type,
  data=qtls_w_ortho[No_Bay_or_Sha==F & single_copy_gene==T,]
)






###########this is just to visualize the density of the markers along the genome
ggplot(gmap, aes(x = pos, y = 1, group = chr)) +
  geom_point() +  # Plot points
  geom_line(aes(color = as.factor(chr))) +  # Connect points with lines, color by chromosome
  facet_wrap(~ chr, scales = "free_x", ncol = 1) +  # Facet by chromosome
  labs(x = "Position along chromosome", y = NULL, title = "Marker Positions by Chromosome") +
  theme_minimal() +
  theme(
    axis.title.y = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "none"
  ) +
  scale_x_continuous(breaks = seq(0, max(gmap$pos), by = 10))


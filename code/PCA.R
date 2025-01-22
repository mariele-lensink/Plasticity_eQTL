library(data.table)
library(ggplot2)
library(dplyr)
library(patchwork)
library(cowplot)

#setting up the input tables, rils and parents
rils<-fread("data/rils_pheno.txt")
rils[, type := "ril",]
cols <- c("ID", "treatment", "type")
rest <- setdiff(colnames(rils), cols)
setcolorder(rils, c(cols,rest))

parents<-fread("data/parents_pheno_corrected.txt")
parents<-parents[,lapply(.SD,mean),by = .(treatment,ID)]
parents$rep<-NULL
parents[,type := "parent",]
cols <- c("ID", "treatment", "type")
rest <- setdiff(colnames(parents), cols)
setcolorder(parents, c(cols,rest)) 
colnames(parents)<-gsub("-",".",colnames(parents))

#need to make the delta group for parents (SA+SW)/((SA+SW)/2)
dpbay<-data.table((parents[1,-c(1:3)]-parents[3,-c(1:3)])/((parents[1,-c(1:3)]+parents[3,-c(1:3)])/2))
dpbay<-data.table(ID="Bay",treatment="delta",type="parent",dpbay)
dpsha<-((parents[2,-c(1:3)]-parents[4,-c(1:3)])/((parents[2,-c(1:3)]+parents[4,-c(1:3)])/2))
dpsha<-data.table(ID="Sha",treatment="delta",type = "parent",dpsha)

#combined all rils and parents, all treatment groups into 1 table
dt<-rbind(parents,dpbay,dpsha,rils)

####removing all genes with a cis qtl from the table before doing PCA ######
qtls<-fread("data/qtls_10cMwindow_notcollapsed_rsq_March18.txt")
cistx<-unique(qtls[type=="cis",txname])
dt_cisrm<-dt[,.SD, .SDcols=!cistx]







run_pca_plot <- function(data, treatment_group) {
  # Subset data
  subdt <- data[treatment == treatment_group, -c('treatment', 'type'), with = FALSE]
  
  # Perform PCA
  subdt_pca_scaled <- prcomp(subdt[,-1], scale = TRUE)
  pcs <- as.data.table(data.frame(subdt_pca_scaled$x))
  pcs$ID <- subdt$ID
  
  # Calculate variance percentages
  pca_var <- subdt_pca_scaled$sdev^2
  pca_var_per <- round(pca_var / sum(pca_var) * 100, 2)
  
  #Plot
  p <- ggplot(pcs, aes(x = PC1, y = PC2)) +
    geom_point(size = 0.5) +
    theme_bw() +
    geom_point(data = pcs[ID == "Bay"], aes(col = ID), size = 2.5, colour = "#56B4E9") +
    geom_point(data = pcs[ID == "Sha"], aes(col = ID), size = 2.5, colour = "#CC79A7") +
    xlab(paste0("PC1 - ", pca_var_per[1], "%")) +
    ylab(paste0("PC2 - ", pca_var_per[2], "%")) +
    theme(text = element_text(size = 9), 
          axis.text = element_text(size = 9),
          panel.grid = element_blank())
  
  return(p)
}

# Define datasets and treatment groups
datasets <- list(dataset1 = dt, dataset2 = dt_cisrm)
treatment_groups <- c('SA', 'SW','delta')
plot_list <- list()
#########################################################################################################################
#################################################################################################################################
library(ggplot2)
library(data.table)

run_pca_plot <- function(data, treatment_group) {
  # Subset data
  subdt <- data[treatment == treatment_group, -c('treatment', 'type'), with = FALSE]
  
  # Perform PCA
  subdt_pca_scaled <- prcomp(subdt[, -1, with = FALSE], scale = TRUE)
  pcs <- as.data.table(data.frame(subdt_pca_scaled$x))
  pcs$ID <- subdt$ID
  
  # Calculate variance percentages
  pca_var <- subdt_pca_scaled$sdev^2
  pca_var_per <- round(pca_var / sum(pca_var) * 100, 2)
  
  # Plot
  p <- ggplot(pcs, aes(x = PC1, y = PC2)) +
    geom_point(size = 0.2) +
    theme_bw() +
    geom_point(data = pcs[ID == "Bay"], aes(col = ID), size = 1, colour = "#56B4E9") +
    geom_point(data = pcs[ID == "Sha"], aes(col = ID), size = 1, colour = "#CC79A7") +
    xlab(paste0("PC1 - ", pca_var_per[1], "%")) +
    ylab(paste0("PC2 - ", pca_var_per[2], "%")) +
    theme(text = element_text(size = 6), 
          axis.text = element_text(size = 6),
          axis.text.x = element_text(angle = 45),
          panel.grid = element_blank())
  
  return(p)
}

# Define datasets and treatment groups
datasets <- list(dataset1 = dt, dataset2 = dt_cisrm)
treatment_groups <- c('SA', 'SW', 'delta')

# Apply function across datasets and treatment groups
plots <- lapply(names(datasets), function(ds_name) {
  list_plots <- lapply(treatment_groups, function(tg) {
    run_pca_plot(datasets[[ds_name]], tg)
  })
})
p1 <- plots[[1]][[1]]
p2 <- plots[[1]][[2]]
p3 <- plots[[1]][[3]]
p4 <- plots[[2]][[1]]
p5 <- plots[[2]][[2]]
p6 <- plots[[2]][[3]]
 
grid <- plot_grid(
  p1, p2, p3,
  p4, p5, p6,
  ncol = 3, nrow = 2,
  align = 'hv', 
  axis = 'tb'
)
full_plot <- ggdraw() +
  draw_plot(grid, x = 0.05, y = 0.01, width = 0.95, height = 0.95) +
  draw_label("All Transcripts",angle = 90, x = 0.02, y = 0.67, hjust = 0, vjust = 0, size = 7, fontface = 'bold') +
  draw_label("Transcripts with\ncis-eQTL removed", angle = 90,x = 0.05, y = 0.18, hjust = 0, vjust = 0, size = 7, fontface = 'bold') +
  draw_label("Salicylic Acid", x = 0.26, y = 0.98, hjust = 0.5, vjust = 1, size = 7, fontface = 'bold', angle = 0) +
  draw_label("Silwet", x = 0.58, y = 0.98, hjust = 0.5, vjust = 1, size = 7, fontface = 'bold', angle = 0) +
  draw_label("Delta", x = 0.9, y = 0.98, hjust = 0.5, vjust = 1, size = 7, fontface = 'bold', angle = 0)
full_plot




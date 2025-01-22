source("github/plasticity_eQTL/code/functions.R")
library(ggplot2)
library(parallel)
library(forcats)
library(cowplot)

rils<-fread("data/rils_pheno_w_reps.txt")
parents<-fread("data/parents_pheno_corrected.txt")

aov_sumsqs_rils<-run_parallel_anova(rils)
aov_sumsqs_parents<-run_parallel_anova(parents)

fwrite(aov_sumsqs_rils,"data/propvariance_rils.txt")
fwrite(aov_sumsqs_parents,"data/propvariance_parents_corrected.txt")
##################################################################
#parents
parents<-fread("data/propvariance_parents_corrected.txt")
parents<-parents[order(-residuals,trt.geno.int)]
psub<-parents[,.(treatment,genotypev,trt.geno.int)]
set.seed(123) # For reproducibility
k <- 10 # Adjust this based on your knowledge of the data
clusters <- kmeans(psub[, c(1:3), with = FALSE], centers = k) # Excluding `row_id` from clustering
psub$cluster <- clusters$cluster

total_samples <- 10000
num_clusters <- length(unique(psub$cluster))
samples_per_cluster <- total_samples / num_clusters

# Sample equally from each cluster
set.seed(123)  # Ensure reproducibility
sampled_rows <- psub[, .SD[sample(.N, min(.N, samples_per_cluster))], by = cluster]
sampled_rows[, row_id := .I] 
long_sampled_data <- melt(sampled_rows, id.vars = c("cluster","row_id"), variable.name = "Component", value.name = "Percentage")
long_sampled_data$type<-"Parents"

#rils
rils<-fread("data/propvariance_rils.txt")
rils<-rils[,.(treatment,genotypev,trt.geno.int)]
set.seed(123) # For reproducibility
k <- 10 # Adjust this based on your knowledge of the data
clusters <- kmeans(rils[, c(1:3), with = FALSE], centers = k) # Excluding `row_id` from clustering
rils$cluster <- clusters$cluster

total_samples <- 10000
num_clusters <- length(unique(rils$cluster))
samples_per_cluster <- total_samples / num_clusters

# Sample equally from each cluster
set.seed(123)  # Ensure reproducibility
sampled_rows <- rils[, .SD[sample(.N, min(.N, samples_per_cluster))], by = cluster]
sampled_rows[, row_id := .I] 
long_sampled_data2 <- melt(sampled_rows, id.vars = c("cluster","row_id"), variable.name = "Component", value.name = "Percentage")
long_sampled_data2$type = "RILs"

longdata<-rbind(long_sampled_data,long_sampled_data2)

p1<-ggplot(longdata, aes(x = as.factor(row_id), y = Percentage, fill = Component)) +
  geom_bar(stat = "identity", position = "stack") +
  labs(x = "Transcript", y = "Percent Variance in Transcript\nAbundance (Sum of Squares)", fill = "Component") +
  guides(fill = guide_legend(title = "Component")) +
  facet_wrap(~ type)+
  scale_fill_discrete(labels = c("treatment" = "Treatment",
                                 "genotypev"= "Genotype", 
                                 "trt.geno.int" = "Genotype X\nTreatment"))+
  theme(panel.background = element_blank(),
        legend.position = c(0.2, 0.8),
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank(),
        #axis.title.x =  element_text(margin=margin(t=2, b=1,size = 6)),
        legend.title = element_blank(),
        legend.text = element_text(size = 6),
        legend.key.size = unit(0.2,'cm'),
        legend.key.spacing.y = unit(0.1,'cm'),
        text = element_text(size = 6),
        legend.background = element_blank()
        )


parents$type<-"Parents"
rils$type<-"RILs"
dt<-rbind(parents,rils)
dt_long<-melt(dt,id.vars = "type",variable.name = "Component", value.name = "Value")
means<-dt_long[,.(Average = mean(Value)),by = c("type","Component")]
p2<-ggplot(means,aes(x=as.factor(type),y = Average, fill = Component))+
  geom_bar(stat = "identity",position = "stack")+
  xlab("Genotype")+
  ylab("Average Percent Variance")+
  scale_fill_discrete(labels = c("treatment" = "Treatment",
                                 "genotypev"= "Genotype", 
                                 "trt.geno.int" = "Genotype X\nTreatment",
                                 "replicate" = "Replicate",
                                 "residuals"="Residuals"))+
  theme(panel.background = element_blank(),
        axis.text.x = element_text(angle=50, vjust=1,size = 6,hjust = 1,margin=margin(t=-10,b=1)),  # Adjust top and bottom margins
        legend.text = element_text(size = 6),
        axis.ticks.x = element_blank(),
        legend.title = element_blank(),
        axis.title.x =  element_text(margin=margin(t=1, b=1)),
        axis.title.y = element_text(size=6),
        legend.margin = margin(t=0, b=0, l=-10, r=0),
        legend.key.size = unit(0.2,'cm'),
        legend.key.spacing.y = unit(0.1,'cm'),
        text = element_text(size = 6))

plot_grid(p1, p2, labels = c('A', 'B'), label_size = 12, rel_widths = c(1.9,1))


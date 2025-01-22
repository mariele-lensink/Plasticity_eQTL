library(data.table)
library(ggplot2)
library(rlist)
library(reshape2)
library(patchwork)
library(cowplot)
source("github/plasticity_eQTL/code/functions.R")

qtls<-fread("data/qtls_10cMwindow_notcollapsed_rsq_March18.txt")
qtls<-qtls[type=="cis"]
#reordering peaks by chromosome then position on chromosome (not sure if necessary but it makes me feel better)
qtls_list<-split(qtls,by='trt')

qtls_list_orderded<-lapply(qtls_list, order_split)

windowcounts_split <- apply_slider_to_processed_list(qtls_list_orderded, 10, 1, 105)

sliding_window_counts_ggplot <- reformat_for_ggplot(windowcounts_split)

################## Permutations ##################################################################
#####RANDOMIZE QTLs, 1,000 permutations############################################################################
#get number of qtls
deltaqtlcount<-qtls[trt=='delta',.N]
saqtlcount<-qtls[trt=='SA',.N]
swqtlcount<-qtls[trt=='SW',.N]

#get length of each chromosome, some to total length in centimorgans
totalcms<-max(qtls[qtl_chr==1,qtl_pos])+max(qtls[qtl_chr==2,qtl_pos])+
  max(qtls[qtl_chr==3,qtl_pos])+max(qtls[qtl_chr==4,qtl_pos])+max(qtls[qtl_chr==5,qtl_pos])
#generateing a random uniform distribution of qtls across the genome, same number for each treatment group
deltaperms<-generate_random_qtl_data_table(deltaqtlcount,totalcms,permutations = 1000)
SAperms<-generate_random_qtl_data_table(saqtlcount,totalcms,permutations = 1000)
SWperms<-generate_random_qtl_data_table(swqtlcount,totalcms,permutations = 1000)


#run randomized sliding window on all of the permutations
deltaperm_windowcounts<-as.data.table(apply(deltaperms,2,function(x) centimorgan_slider(x,10,1,430))) 
SAperm_windowcounts<- as.data.table(apply(SAperms,2,function(x) centimorgan_slider(x,10,1,430))) 
SWperm_windowcounts<- as.data.table(apply(SWperms,2,function(x) centimorgan_slider(x,10,1,430))) 
#get 5% cutoffs for values in each group
deltacutoff<-quantile(as.vector(as.matrix(deltaperm_windowcounts)),0.95)
SAcutoff<-quantile(as.vector(as.matrix(SAperm_windowcounts)),0.95)
SWcutoff<-quantile(as.vector(as.matrix(SWperm_windowcounts)),0.95)


ciscounts <- sliding_window_counts_ggplot[!(Count == "0" & Index > 73)]
transcounts <-sliding_window_counts_ggplot[!(Count == "0" & Index > 73)]
custom_colors <- c("SA" = "#377EB8", "delta" = "#E6550D", "SW" = "#4DAF4A")
cis<-ggplot(ciscounts,aes(x=Index,y=Count,color = Treatment))+
  scale_color_manual(values=custom_colors)+
  scale_y_continuous(expand = c(0,0))+
  geom_line()+facet_wrap(~Chromosome,nrow = 1,scales = 'free_x')+
  theme_bw()+
  ylab('# of Transcripts (Cis)')+
  xlab(NULL)+
  theme(text = element_text(size = 6),
        axis.text = element_text(size=6),
        axis.text.x = element_text(angle = 45, vjust = 1.25, hjust=1),
        legend.background = element_rect(fill='transparent'),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.key.size = unit(0.4, "cm"),  # Adjust the key size
        legend.spacing.y = unit(0.05, "cm"),  # Adjust spacing between legend items
        legend.box.margin = margin(1, 1, 1, 1),
        #legend.position=c(0.1,0.75)
        legend.position = 'none')+
  geom_hline(yintercept = deltacutoff, linetype="dashed", color = "#E6550D", linewidth=.4)+
  geom_hline(yintercept = SAcutoff, linetype="dashed", color = "#377EB8", linewidth=.4)+
  geom_hline(yintercept = SWcutoff, linetype="dashed", color = "#4DAF4A", linewidth=.4)
 


trans <- trans + theme(plot.margin = unit(c(0.1, 1, 0.1, 1), "lines"))

# Adjusting margins of the cis plot (reduce top margin)
cis <- cis + theme(plot.margin = unit(c(0.1, 1, 0.11, 1), "lines"))

# Combine the plots with adjusted margins
combined_plot <- plot_grid(trans, cis, labels = c('B', 'C'), label_size = 12, rel_widths = c(1, 3),ncol=1)

plot_grid(fig4a,combined_plot, labels = c('A', ''), label_size = 12, rel_widths = c(1,4))


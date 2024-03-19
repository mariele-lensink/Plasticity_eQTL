library(data.table)
library(ggplot2)
library(tidyr)
library(ggupset)
library(stringr)
library(forcats)
library(dplyr)
library(reshape2)
library(ggpmisc)
library(devtools)
library(plotly)
library(rcartocolor)
library(patchwork)

qtls<-fread("qtls_10cMwindow_notcollapsed_rsq_March18.txt")
x<-data.table(table(qtls[,type,trt]))
colnames(x)[3]<-'txcount'
ggplot(x, aes(fill = type, x = trt, y= txcount))+
  geom_col(position = 'fill')+
  geom_text(aes(label= txcount),position = position_fill(vjust=0.5),size = 3)+
  scale_y_continuous()+
  ylab("Proportion of Cis and Trans eQTLs")+
  xlab("Treatment Group")+
  labs(fill="Regulatory Type")+
  theme_bw()+
  theme(legend.position="none",text = element_text(size = 9),axis.text = element_text(size=9))+
  scale_fill_grey(start = 0.5,end=0.8)


#need to use the collapsed data table for these figures
dt<-fread("qtls_10cMwindow_fixed_overlap_5.5cm_march18.txt")
p1<-ggplot(dt,aes(x=trtgroup,fill=type),scales = 'free',space = 'free')+geom_bar(position = 'dodge')+
  scale_x_discrete(limits=c("SA","SW","delta","SA delta","SW delta","SA SW","SA SW delta"),expand = c(0,0))+
  scale_fill_grey(start = 0.5,end=0.8)+
  scale_y_continuous(limits = c(0,19000), expand = c(0,0))+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        text = element_text(size = 9),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.title.x=element_blank())+
  ylab("# of QTLs")+
  geom_text(stat='count', aes(label=after_stat(count)), vjust=-.2,size = 2,position = position_dodge(width=.9))


p2<-ggplot(dt,aes(x=trtgroup,y = rsq,factor = type))+
  scale_x_discrete(limits=c("SA","SW","delta","SA delta","SW delta","SA SW","SA SW delta"),expand = c(0,0))+
  axis_combmatrix(sep = " ")+
  scale_y_continuous(limits = c(0,1),expand = c(0,0))+
  geom_point(aes(color = type,alpha = 0.10),position=position_jitterdodge(0.2))+
  geom_violin()+
  xlab("Treatment Group Membership(s)")+
  ylab("Effect Size (rsq)")+
  theme_bw()+
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "none",
        text = element_text(size = 9))+
  scale_color_grey(start = 0.5,end=0.8)



p1 + p2 + plot_layout(ncol = 1)


library(reshape2)
source("../../../../../Raquel_Gur/GO/substance/imaging_analysis/summarySE.R")

theme_set(
  theme_bw() +
    theme(legend.position = "top")
  )

library(reshape2)
library(ggplot2)
library(zoo)
library(dplyr)
library(table1)
library(lattice)
library(tidyverse)

dataf<-read.csv("rarecnv_zscores_winsorized8_qced_validated_imputed_apr2024.csv")

dataf<-dataf[,c(5,6,7,8,31:52)]


#long format
dataf_long<-melt(dataf,id=c("rarecnv_id","Group","sex","age"))
names(dataf_long)[5:6]<-c("cnb_test","mean_zscore")

#split column cnb_test and Group
library(stringr)
dataf_long[c('type','test', 'measure')] <- str_split_fixed(dataf_long$cnb_test, '_', 3)

#cnb2$remote<-1
#cnb2$remote[cnb2$platform=="webcnp"]<-0

dataf_long$DelDup<-str_sub(dataf_long$Group, - 3, - 1)
dataf_long$Locus<-str_sub(dataf_long$Group, 1,  3)
###final plots


out<-summarySE(dataf_long,groupvars=c("Locus","DelDup","measure","test"),measurevar="mean_zscore",na.rm=T)


outa<- out[ out$measure=="a" , ]
outs<- out[ out$measure=="s" , ]


#change order
outa$test <- factor(outa$test, levels = c("abf","att","wm","fmem","smem","nvr","spa","eid","edi","adi")) #order needs change for speed
outs$test <- factor(outs$test, levels = c("abf","att","wm","fmem","smem","nvr","spa","eid","edi","adi","sm","mot"))


outa_16p<-outa[outa$Locus=="16p", ]
outa_22q<-outa[outa$Locus=="22q", ]
outs_16p<-outs[outs$Locus=="16p", ]
outs_22q<-outs[outs$Locus=="22q", ]


pd <- position_dodge(0.1) # move them .05 to the left and right


#final plot
library(grid)
grid.newpage()
pushViewport(viewport(layout = grid.layout(nrow = 2, ncol = 2)))
define_region <- function(row, col){
viewport(layout.pos.row = row, layout.pos.col = col)}

#plots
p1<-ggplot(data=outa_22q,aes(x=test, y=mean_zscore, group=DelDup, colour=DelDup)) +
        geom_line() +
        geom_point()+
          geom_errorbar(aes(ymin=mean_zscore-se, ymax=mean_zscore+se), width=.1)+scale_y_continuous(limits=c(-2,0.5))+geom_hline(yintercept=0)+theme(legend.position = "none")

p2<-ggplot(data=outs_22q,aes(x=test, y=mean_zscore, group=DelDup, colour=DelDup)) +
                geom_line() +
                geom_point()+
                  geom_errorbar(aes(ymin=mean_zscore-se, ymax=mean_zscore+se), width=.1)+scale_y_continuous(limits=c(-2,0.5))+geom_hline(yintercept=0)+theme(legend.position = "none")


p3<-ggplot(data=outa_16p,aes(x=test, y=mean_zscore, group=DelDup, colour=DelDup)) +
                        geom_line() +
                        geom_point()+
                          geom_errorbar(aes(ymin=mean_zscore-se, ymax=mean_zscore+se), width=.1)+scale_y_continuous(limits=c(-2,0.5))+geom_hline(yintercept=0)+theme(legend.position = "none")


p4<-ggplot(data=outs_16p,aes(x=test, y=mean_zscore, group=DelDup, colour=DelDup)) +
                                geom_line() +
                                geom_point()+
                                  geom_errorbar(aes(ymin=mean_zscore-se, ymax=mean_zscore+se), width=.1)+scale_y_continuous(limits=c(-2,0.5))+geom_hline(yintercept=0)+theme(legend.position = "none")




  library(grid)
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(nrow = 1, ncol = 2)))
  define_region <- function(row, col){
  viewport(layout.pos.row = row, layout.pos.col = col)}

print(p1, vp = define_region(row = 1, col = 1))
print(p2, vp = define_region(row = 1, col = 2))

grid.newpage()
pushViewport(viewport(layout = grid.layout(nrow = 1, ncol = 2)))
define_region <- function(row, col){
viewport(layout.pos.row = row, layout.pos.col = col)}

print(p3, vp = define_region(row = 1, col = 1))
print(p4, vp = define_region(row = 1, col = 2))

#not used
p5<-ggplot(data=outa,aes(x=cnb_test, y=mean_zscore, group=Group, colour=Group)) +
                                geom_line() +
                                geom_point()+
                                  geom_errorbar(aes(ymin=mean_zscore-se, ymax=mean_zscore+se), width=.1)+scale_y_continuous(limits=c(-4,1))+geom_hline(yintercept=0)+
scale_x_discrete(guide = guide_axis(n.dodge = 2))+ scale_colour_manual(values = c("red", "orange", "blue","violet"))


p6<-ggplot(data=outs,aes(x=cnb_test, y=mean_zscore, group=Group, colour=Group)) +
                                geom_line() +
                                geom_point()+
                                  geom_errorbar(aes(ymin=mean_zscore-se, ymax=mean_zscore+se), width=.1)+scale_y_continuous(limits=c(-4,1))+geom_hline(yintercept=0)+
scale_x_discrete(guide = guide_axis(n.dodge = 2))+ scale_colour_manual(values = c("red", "orange", "blue","violet"))+theme(legend.justification = 'left', legend.position="top")


#barplots
dataf<-read.csv("rarecnv_zscores_winsorized8_qced_validated_imputed_apr2024.csv")

dataf<-dataf[,c(5,6,7,8,53:57)]
 dataf$DelDup<-str_sub(dataf$Group, - 3, - 1)
dataf$Locus<-str_sub(dataf$Group, 1,  3)



out_ciqa<-summarySE(dataf,groupvars=c("Locus","DelDup"),measurevar="imputed_ciq_a",na.rm=T)
out_ciqs<-summarySE(dataf,groupvars=c("Locus","DelDup"),measurevar="imputed_ciq_s",na.rm=T)


#change order
out_ciqa$Locus <- factor(out_ciqa$Locus, levels = c("22q","16p"))
out_ciqa$Group<-paste0(out_ciqa$Locus, out_ciqa$DelDup)
out_ciqa$Group <- factor(out_ciqa$Group, levels = c("22qDel","22qDup","16pDel","16pDup"))

# barplot
p11 <- ggplot(out_ciqa, aes(x=Group, y=imputed_ciq_a,fill="blue")) +
   geom_bar(stat="identity", position=position_dodge()) +
  geom_errorbar(aes(ymin=imputed_ciq_a-se, ymax=imputed_ciq_a+se), width=.2,
                 position=position_dodge(.9))+scale_y_continuous(limits=c(-1,0))+theme(legend.position = "none")+
  theme(axis.text.x = element_text(size=10, face="bold", color = "black"),
        axis.text.y = element_text(size=10, face="bold", color = "black"))+
 font("xlab", size = 12,face="bold")


out_ciqs$Locus <- factor(out_ciqs$Locus, levels = c("22q","16p"))
out_ciqs$Group<-paste0(out_ciqs$Locus, out_ciqs$DelDup)
out_ciqs$Group <- factor(out_ciqs$Group, levels = c("22qDel","22qDup","16pDel","16pDup"))



p12 <- ggplot(out_ciqs, aes(x=Group, y=imputed_ciq_s,fill="blue")) +
                    geom_bar(stat="identity", position=position_dodge()) +
                   geom_errorbar(aes(ymin=imputed_ciq_s-se, ymax=imputed_ciq_s+se), width=.2,
                                  position=position_dodge(.9))+scale_y_continuous(limits=c(-1,0))+theme(legend.position = "none")+
  theme(axis.text.x = element_text(size=10, face="bold", color = "black"),
        axis.text.y = element_text(size=10, face="bold", color = "black"))



                                  grid.newpage()
                                  pushViewport(viewport(layout = grid.layout(nrow = 1, ncol = 2)))
                                  define_region <- function(row, col){
                                  viewport(layout.pos.row = row, layout.pos.col = col)}

                                  print(p11, vp = define_region(row = 1, col = 1))
                                  print(p12, vp = define_region(row = 1, col = 2))



#supplemental figure box and whisker barplots



### box plot with whiskers

dataf<-read.csv("rarecnv_zscores_winsorized8_qced_validated_imputed_apr2024.csv")

dataf<-dataf[,c(5,6,7,8,31:52)]


#long format
dataf_long<-melt(dataf,id=c("rarecnv_id","Group","sex","age"))
names(dataf_long)[5:6]<-c("cnb_test","mean_zscore")

#split column cnb_test and Group
library(stringr)
dataf_long[c('type','test', 'measure')] <- str_split_fixed(dataf_long$cnb_test, '_', 3)

#cnb2$remote<-1
#cnb2$remote[cnb2$platform=="webcnp"]<-0

dataf_long$DelDup<-str_sub(dataf_long$Group, - 3, - 1)
dataf_long$Locus<-str_sub(dataf_long$Group, 1,  3)

data<-dataf_long
ggboxplot(data[which(data$measure=="a"),], x = "test", y = "mean_zscore", color = "Group",palette =c("#ca7dcc", "#1b98e0", "#353436","#02e302"))

ggboxplot(data[which(data$measure=="s"),], x = "test", y = "mean_zscore", color = "Group",palette =c("#ca7dcc", "#1b98e0", "#353436","#02e302"))

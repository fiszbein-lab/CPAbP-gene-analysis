library(tidyverse)
library(reshape2)
library(ggpubr)

setwd("data/")

#load HITindex outcome
HITstat_ctrl_vs_u1amo <- read.table("HITindex_stat_ctrl_vs_u1amo",header = T)
head(HITstat_ctrl_vs_u1amo)

#load annotation file for CPAbP classification
all_bed <- read.table("HeLa_genes_with_atleast_2promoters.txt",header = T)
m_all_bed <- melt(all_bed,measure.vars = c("exon_go1","exon_go2"),value.name = "exon",variable.name = "genomic_order")

#classification
nonCPAbP_go1 <- m_all_bed%>%filter(CPAbP=="N", U1_dependent_CPAbP == "N",genomic_order=="exon_go1")
nonCPAbP_go2 <- m_all_bed%>%filter(CPAbP=="N", U1_dependent_CPAbP == "N",genomic_order=="exon_go2")
CPAbP_go1 <- m_all_bed%>%filter(CPAbP=="Y",genomic_order=="exon_go1")
CPAbP_go2 <- m_all_bed%>%filter(CPAbP=="Y",genomic_order=="exon_go2")


#annotate HITstat outcome and normalize to the mean delta_PSI of genomic_order1
HITstat_nonCPAbP_go1 <- HITstat_ctrl_vs_u1amo%>%subset(exon %in% nonCPAbP_go1$exon)%>%mutate(group="nonCPAbP",genomic_order = "1")%>%filter(bio_significant=="True")
HITstat_nonCPAbP_go2 <- HITstat_ctrl_vs_u1amo%>%subset(exon %in% nonCPAbP_go2$exon)%>%mutate(group="nonCPAbP",genomic_order = "2")%>%filter(bio_significant=="True")
HITstat_nonCPAbP <- rbind(HITstat_nonCPAbP_go1,HITstat_nonCPAbP_go2)
HITstat_nonCPAbP$delta_PSI <- HITstat_nonCPAbP$delta_PSI - median(HITstat_nonCPAbP_go1$delta_PSI)

HITstat_CPAbP_go1 <- HITstat_ctrl_vs_u1amo%>%subset(exon %in% CPAbP_go1$exon)%>%mutate(group="CPAbP",genomic_order = "1")%>%filter(bio_significant=="True")
HITstat_CPAbP_go2 <- HITstat_ctrl_vs_u1amo%>%subset(exon %in% CPAbP_go2$exon)%>%mutate(group="CPAbP",genomic_order = "2")%>%filter(bio_significant=="True")
HITstat_CPAbP <- rbind(HITstat_CPAbP_go1,HITstat_CPAbP_go2)
HITstat_CPAbP$delta_PSI <- HITstat_CPAbP$delta_PSI - median(HITstat_CPAbP_go1$delta_PSI)


#combine normalized results from each group
combined_HITstat <- rbind(HITstat_nonCPAbP,HITstat_CPAbP)
combined_HITstat$group <- factor(combined_HITstat$group, levels=c("nonCPAbP","CPAbP"))

#generate boxplot with t-test
ggplot(combined_HITstat, aes(x=group,y=delta_PSI,fill=genomic_order))+geom_boxplot(notch=T)+theme_minimal()+
  coord_cartesian(ylim=c(-1,1))+
  theme(legend.position = "none")+
  stat_compare_means(method = "t.test",label.y = 1)+
  scale_fill_manual(values = c("lavenderblush4","lightseagreen"))

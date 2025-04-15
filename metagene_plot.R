library(metagene2)
library(GenomicRanges)
library(tidyverse)

setwd("/projectnb/encore/Yun/")


generate_bin_around_TSS <- function(exon){
exon_plus <- exon%>%filter(strand=="+")
exon_plus$lower <- exon_plus$upper+300
exon_plus$upper <- exon_plus$upper-300

exon_minus <- exon%>%filter(strand=="-")
exon_minus$upper <- exon_minus$lower-300
exon_minus$lower <- exon_minus$lower+300

exon <- rbind(exon_plus,exon_minus)


TSS_grange <- makeGRangesFromDataFrame(exon,
                                           keep.extra.columns=F,
                                           ignore.strand=FALSE,
                                           seqinfo=NULL,
                                           seqnames.field="chr",
                                           start.field="upper",
                                           end.field="lower",
                                           strand.field="strand",
                                           starts.in.df.are.0based=FALSE)

return(TSS_grange)
}

antisense_generate_bin_around_TSS <- function(exon){
  exon_plus <- exon%>%filter(strand=="+")
  exon_plus$lower <- exon_plus$upper+300
  exon_plus$upper <- exon_plus$upper-300
  
  exon_minus <- exon%>%filter(strand=="-")
  exon_minus$upper <- exon_minus$lower-300
  exon_minus$lower <- exon_minus$lower+300
  
  exon_plus$strand <- "-"
  exon_minus$strand <- "+"
  
  exon <- rbind(exon_plus,exon_minus)
  
  TSS_grange <- makeGRangesFromDataFrame(exon,
                                         keep.extra.columns=F,
                                         ignore.strand=FALSE,
                                         seqinfo=NULL,
                                         seqnames.field="chr",
                                         start.field="upper",
                                         end.field="lower",
                                         strand.field="strand",
                                         starts.in.df.are.0based=FALSE)
  
  return(TSS_grange)
}


#import CPAbP/nonCPAbP info
all_bed <- read.table("HeLa_genes_with_atleast_2promoters.txt",header = T)
m_all_bed <- melt(all_bed,measure.vars = c("exon_go1","exon_go2"),value.name = "exon",variable.name = "genomic_order")%>%
  separate(exon,c("chr","coordinate"),sep = ":")%>%separate(coordinate,c("upper","lower"),sep = "-",convert = T)

##generate Granges for coordinates of up/down 300nt from TSS - sense
grange_nonCPAbP_go1 <- m_all_bed%>%filter(CPAbP=="N",genomic_order=="exon_go1")%>%generate_bin_around_TSS()
grange_nonCPAbP_go2 <- m_all_bed%>%filter(CPAbP=="N",genomic_order=="exon_go2")%>%generate_bin_around_TSS()
grange_CPAbP_go1 <- m_all_bed%>%filter(CPAbP=="Y",genomic_order=="exon_go1")%>%generate_bin_around_TSS()
grange_CPAbP_go2 <- m_all_bed%>%filter(CPAbP=="Y",genomic_order=="exon_go2")%>%generate_bin_around_TSS()

#combine Granges
TSS_list_CPAbPnonCPAbP <- GRangesList(grange_nonCPAbP_go1,grange_nonCPAbP_go2,grange_CPAbP_go1,grange_CPAbP_go2)
names(TSS_list_CPAbPnonCPAbP) <- c("nonCPAbP_GO1","nonCPAbP_GO2","CPAbP_GO1","CPAbP_GO2")


#import BAM files to look at
bamlist <- c("ctrl1.bam","ctrl2.bam",
             "u1amo1.bam","u1amo2.bam")
names(bamlist) = c("ctrl_rep1", "ctrl_rep2","u1amo_rep1","u1amo_rep2")


#generate metagene object
mg <- metagene2$new(regions=TSS_list_CPAbPnonCPAbP, 
                    bam_files = bamlist,normalization ="RPM",
                    strand_specific =T,bin_count=100)

#antisense
mg_anti <- metagene2$new(regions=anti_TSS_list_CPAbPnonCPAbP, 
                    bam_files = bamlist,normalization ="RPM",
                    strand_specific =T,bin_count=100)

#grouping metadata
design_meta = data.frame(design=mg$get_design_group_names(),
                         Align=c("Ctrl", "Ctrl", "U1AMO", "U1AMO"),
                         Rep=c("rep1", "rep2", "rep1", "rep2"))

#draw metagene plot
mg$produce_metagene(design_metadata=design_meta, facet_by=~region,group_by="Align")+
  scale_color_manual(values=c("#696969","#CC6699"))+
  scale_fill_manual(values=c("#696969","#CC6699"))

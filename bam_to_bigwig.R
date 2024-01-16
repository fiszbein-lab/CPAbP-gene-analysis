library(GenomicAlignments)
library(rtracklayer)


bam2bigwig <- function(input.bam,normalization_factor,output.bw){
  gr <- readGAlignments(input.bam)
  
  #extract reads from chr1-24, X, Y, and MT
  gr.cov <- coverage(gr)[1:25]

    for (i in 1:25) {
    gr.cov@listData[[i]]@values <- gr.cov@listData[[i]]@values/normalization_factor
    
  }
  export.bw(gr.cov,output.bw)
}

bam2bigwig_plus <- function(input.bam,normalization_factor,output.bw){
  gr <- readGAlignments(input.bam)
  gr <- gr[gr@strand=="+"]
  
  gr.cov <- coverage(gr)
  for (i in 1:25) {
    gr.cov@listData[[i]]@values <- gr.cov@listData[[i]]@values/normalization_factor
    
  }
  export.bw(gr.cov,output.bw)
}

bam2bigwig_minus <- function(input.bam,normalization_factor,output.bw){
  gr <- readGAlignments(input.bam)
  gr <- gr[gr@strand=="-"]
  
  gr.cov <- coverage(gr)
  for (i in 1:25) {
    gr.cov@listData[[i]]@values <- gr.cov@listData[[i]]@values/normalization_factor
    
  }
  export.bw(gr.cov,output.bw)
}


#convert bam file to bigwig file
bam2bigwig("sample.bam",4.698,"sample_combined.bw")

#plus/minus strands separately
bam2bigwig_plus("sample.bam",2.426,"sample_plus.bw")
bam2bigwig_minus("sample.bam",2.272,"sample_minus.bw")


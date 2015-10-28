dist_gene <- function(start1, end1, start2, end2){
  dist <- 0
  if( sign(start2 - start1) * sign(start2 - end1)>0  & sign(end2 - start1) * sign(end2 - end1)>0 &sign(start2 - start1) * sign(end2 - start1) >0 & sign(start2 - end1) * sign(end2 - end1)> 0 ){
    dist <- min (abs(start2 - start1),  abs(start2 - end1), abs(end2 - start1),  abs(end2 - end1))
  }
  dist
}

dist_pairwise <- function(x, gene_complete){
  gene_sub <- subset(gene_complete, name2==as.character(x))
  dist <- c()
  if(dim(gene_sub)[1] >1 ){
    for(i in seq_along(gene_sub$bin)){
      for (j in seq_along(gene_sub$bin) ){
      if(j<=i){
        next()
      }
      dist <- c(dist, dist_gene(gene_sub$txStart[i], gene_sub$txEnd[i], gene_sub$txStart[j], gene_sub$txEnd[j] ))
    }
  }
  }else{
    dist <- 0
  }
  
  max(dist)
}

getGene <- function(x, gene_complete){
  gene_sub <- subset(gene_complete, name2==as.character(x))
  len <- with(gene_sub, txEnd-txStart)
  as.character(gene_sub$name)[which.max(len)]
}

gene_complete <- read.table("refGene.txt", header=F, sep="\t")

names(gene_complete) <- c("bin", "name", "chrom", "strand", "txStart", "txEnd", "cdsStart", "cdsEnd", "exonCount", "exonStarts", "exonEnds", "score", "name2", "cdsStartStat", "cdsEndStat", "exonFrames")

with(gene_complete, table(name)) -> temp

gene_complete <- subset(gene_complete, (! name %in% names(temp)[temp>1]) & cdsStartStat=="cmpl" & cdsEndStat=="cmpl" & grepl("NM_", name) & (substr(chrom, 4, 100) %in% c(c(1:22), "X")) )

with(gene_complete, tapply(chrom, as.character(name2), function(x){length(unique(x))})) -> temp
gene_complete <- subset(gene_complete, ! name2 %in% names(temp)[temp>1])

unlist(lapply(unique(gene_complete$name2), dist_pairwise, gene_complete)) -> temp
gene_complete <- subset(gene_complete, ! name2 %in% unique(gene_complete$name2)[temp>50000]) 

lapply(unique(gene_complete$name2), getGene, gene_complete) -> temp
gene_complete <- gene_complete[match(unlist(as.character(temp)), gene_complete$name), ]

write.table(gene_complete[,c(3:13)], "refGene.cleaned.txt", col.names=T, row.names=F, sep="\t", quote=F)

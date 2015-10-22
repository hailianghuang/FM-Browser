source("global.R")
all_gene <- read.csv("dat/refGene.txt", header=T, sep="\t")
names(all_gene) <- c("chr", "strand", "start", "end", "cdsStart", "cdsEnd", "exonCount", "exonStarts", "exonEnds", "score", "name")
temp <-with(all_gene, substr(chr, 4, 100))
temp[temp=="X"] <- 23
all_gene$chr <- as.integer(temp)

all_snp <- read.table("dat/annotation.variant.txt", header=T, sep="\t")
all_snp <- subset(all_snp, tier2=="No")
region <- read.table("dat/HD_region_HG19.txt", header=F)
names(region) <- c("chr", "start", "end", "region")

result <- c()
for (i in unique(all_snp$HD)){
  
  region_start <-region$start[i] - region_flank_5
  region_end <- region$end[i] + region_flank_5
  gene <- subset(all_gene, chr==region$chr[i]& ((start > region_start& start < region_end )|(end > region_start& end < region_end )))
  if(dim(gene)[1]>0){
    result <- rbind(result, data.frame(HD=i, region=paste(sep="", "chr", region[i,]$chr, ":", sprintf("%0.2f", region[i,]$start/1000000), "-", sprintf("%0.2f", region[i,]$end/1000000)), gene= gene$name))
  }
}

write.table(result, "dat/HD.genelist.txt", col.names=T, row.names=T, sep="\t", quote=F)

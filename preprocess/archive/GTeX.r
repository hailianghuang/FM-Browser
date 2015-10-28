

GTeX <- read.table("preprocess/IBD.subset.rank.txt", sep="\t", header=T)

with(GTeX, tapply(rpkm, list(gene, sample_cat), min)) -> temp

write.table(temp, "GTeX.rank.summary.txt", col.names=T, row.names=T, sep="\t", quote=F)

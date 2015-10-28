annotation <- read.table("annotation.variant.txt", header=T, sep="\t")
annotation <- subset(annotation, tier2=="No")
coef <- read.table("all.coef.txt", header=T, sep="\t")

merge(annotation, subset(coef, select=c("add_SNP", "logOR_CD", "logOR_UC")), by.x="variant", by.y="add_SNP") -> variant

variant[with(variant, order(HD, signal)), ] -> variant
write.table(variant, "../dat/variant.txt", col.names=T, row.names=F, sep="\t", quote=F)


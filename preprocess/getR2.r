dat <- read.table("preprocess/all.txt",  header=T)

ret <- c()
for(i in unique(dat$HD)){
  for (j in unique(dat$signal[dat$HD==i])){
    dat_sub <- subset(dat, HD==i&signal==j)
    r <- with(dat_sub, range(position[R^2>0.4]))
    ret <- rbind(ret, data.frame(HD=i, signal=j, chr=dat_sub$chr[1], start = r[1], end = r[2]))
  }
}

write.table(ret, "dat/R2.span.v2.txt", col.names=T, row.names=T, sep="\t", quote=F)

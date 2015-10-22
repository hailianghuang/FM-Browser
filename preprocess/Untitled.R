dat <- read.table("all.txt",  header=T)

ret <- c()
for(i in unique(dat$HD)){
  for (j in unique(dat$signal[dat$HD==i])){
    dat_sub <- subset(dat, HD==i&signal==j)
    ret <- rbind(ret, data.frame(HD=i, signal=j, chr=dat_sub$chr[1], with(dat_sub, range(position[R^2>0.6]))))
  }
}



getRegion <- function(choice){
  
  if(substr(choice, 1, 2)=="HD"){
    choice <- substr(choice, 3, 10)
  }else{
    choice <- HD.genelist$HD[HD.genelist$gene==choice]
  }
  choice
  
}

rectarrows <- function(x0,y0,x1,y1,height,length,...) {
  lwd=par("lwd")
  l0=height*(y1-y0)/sqrt((x1-x0)^2+(y1-y0)^2)
  l1=height*(x1-x0)/sqrt((x1-x0)^2+(y1-y0)^2)
  d0=length*(y1-y0)/sqrt((x1-x0)^2+(y1-y0)^2)
  d1=length*(x1-x0)/sqrt((x1-x0)^2+(y1-y0)^2)
  polygon(x=c(x0+l0,x1+l0-d1,x1,x1-l0-d1,x0-l0),y=c(y0-l1,y1-l1-d0,y1,y1+l1-d0,y0+l1),...)
}

assign_coding <- function(g, gene, snp){
  
  max(grepl(as.character(gene$name)[g], snp$Coding) * snp$P_mean_95)
  
}

assign_TFBS <- function(g, gene, snp, GTeX) {
  ret <- 0;
  snp_sub <- subset(snp, TFBS!="")
  
  GTeX_sub <- GTeX[as.character(gene$name[g]),]
  if( sum(row.names(GTeX)==gene$name[g]) > 0){
    if(     (!GTeX.options()$compareControl | min(GTeX_sub$GUT, GTeX_sub$IMMMUNE)< GTeX_sub$CONTROL) 
            & (!GTeX.options()$topQ | min(GTeX_sub$GUT, GTeX_sub$IMMMUNE) < 5000) 
    ) {
      for(i in seq_along(snp_sub$position)){
        pos <- snp_sub$position[i]
        if(pos > gene$start[g] - flanking()$flanking_5 && pos < gene$end[g] + flanking()$flanking_3 ){
          ret <- max(ret, snp_sub$P_mean_95[i]);
        }
      }
    }
  }
  
  ret
}

assign_roadmap <- function(g, gene, snp, GTeX) {
  ret <- 0;
  snp_sub <- subset(snp, Roadmap!="")
  
  GTeX_sub <- GTeX[as.character(gene$name[g]),]
  if(sum(row.names(GTeX)==gene$name[g]) > 0){
    if(     (!GTeX.options()$compareControl | min(GTeX_sub$GUT, GTeX_sub$IMMMUNE)< GTeX_sub$CONTROL) 
            & (!GTeX.options()$topQ | min(GTeX_sub$GUT, GTeX_sub$IMMMUNE) < 5000) 
    ) {
      for(i in seq_along(snp_sub$position)){
        pos <- snp_sub$position[i]
        if(pos > gene$start[g] - flanking()$flanking_5 && pos < gene$end[g] + flanking()$flanking_3 ){
          ret <- max(ret,  snp_sub$P_mean_95[i]);
        }
      }
    }
  }
  ret;
}

assign_eQTL <- function(g, gene, snp) {
  ret <- 0;
  snp_sub <- subset(snp, eQTL!="")
  eQTL <- snp_sub$eQTL[1]
  if( gene$name[g] %in% unique(gsub('.+\\((.+)\\)', '\\1', strsplit(as.character(eQTL), ",")[[1]])) ){
    ret <- max(ret, max(snp_sub$P_mean_95))
  }
  ret;
}

assign_other <- function(g, gene, snp, GTeX) {
  ret <- 0;
  snp_sub <- snp
  
  GTeX_sub <- GTeX[as.character(gene$name[g]),]
  if( sum(row.names(GTeX)==gene$name[g]) > 0){
    if(     (!GTeX.options()$compareControl | min(GTeX_sub$GUT, GTeX_sub$IMMMUNE)< GTeX_sub$CONTROL) 
            & (!GTeX.options()$topQ | min(GTeX_sub$GUT, GTeX_sub$IMMMUNE) < 5000) 
    ) {
      for(i in seq_along(snp_sub$position)){
        pos <- snp_sub$position[i]
        if(pos > gene$start[g] - flanking()$flanking_5 && pos < gene$end[g] + flanking()$flanking_3 ){
          ret <- max(ret, snp_sub$P_mean_95[i]);
        }
      }
    }
  }
  
  ret
}

getCredible <- function(score){
  
  prob <- score/sum(score)
  if(sum(score)==0){
    prob <- rep(1/length(score), length(score))
  }
  csum <- unlist(lapply(prob, function(x){sum(prob[prob>=x+1E-10])}))
  csum
}

getAssigment <- function(i, all_snp, all_gene, GTeX, region){
  
  region_start <-region$start[i]-50000
  region_end <- region$end[i]+50000
  gene <- subset(all_gene, chr==region$chr[i]& ((start > region_start& start < region_end )|(end > region_start& end < region_end )) )
  
  result <- c()
  if(dim(gene)[1] > 0){
    
    for (j in unique(all_snp$signal[paste(sep="", "HD", all_snp$HD)==as.character(region$region[i])]) ){
      snp <- all_snp[paste(sep="", "HD", all_snp$HD)==as.character(region$region[i]) & all_snp$signal==j, ]
      
      result <- rbind(result, 
                      data.frame(name=gene$name, 
                                 signal=j, 
                                 coding=unlist(lapply(c(1:length(gene$name)), assign_coding, gene, snp)),
                                 tfbs=unlist(lapply(c(1:length(gene$name)), assign_TFBS, gene, snp, GTeX)),
                                 roadmap_score = unlist(lapply(c(1:length(gene$name)), assign_roadmap, gene, snp, GTeX)),
                                 eQTL_score = unlist(lapply(c(1:length(gene$name)), assign_eQTL, gene, snp)), 
                                 other_score = unlist(lapply(c(1:length(gene$name)), assign_other, gene, snp, GTeX)) 
                                 
                      )
      )
    }
    
  }
  list(region=list(start=region_start, end=region_end, chr=region$chr[i]),
       gene=gene, 
       snp= all_snp[paste(sep="", "HD", all_snp$HD)==as.character(region$region[i]), ], 
       result=result)
}

plot_table <- function(choice, all_snp, all_gene, GTeX, region){
  
  i<- getRegion(choice)
  
  if(password()!="ibdgc"){
    plot(0,0, xlim=c(1, 10), ylim=c(1,10), type="n", axes=F, xlab="", ylab="")
    
  }else{
    result <- getAssigment(i, all_snp, all_gene, GTeX, region)
    
    if(is.null(result$result)){
      plot(0,0, xlim=c(1, 10), ylim=c(1,10), type="n", axes=F, xlab="", ylab="")
      text(5,5, "no gene in this region")
      return()
    }
    
    temp <- as.matrix(result$result[,-c(1,2)]) %*% diag(c(weight()$coding, weight()$tfbs, weight()$roadmap, weight()$eQTL, weight()$other))
    
    nsig <- length(unique(result$result$signal))
    ngene <- dim(temp)[1]/nsig
    gene <- result$gene
    gene$score <- tapply(apply(temp, 1, max), list(rep(c(1:ngene), nsig)), sum )
    gene$keep <- getCredible(gene$score) < CI()
    merge(gene, data.frame(GTeX, name=row.names(GTeX)), by="name") -> temp
    subset(temp, select=c("name", "chr", "strand", "start", "end", "score", "keep", "CONTROL", "GUT", "IMMMUNE"))
  }
}

plot_snp <- function(choice, all_snp, all_gene, GTeX, region){
  
  i<- getRegion(choice)
  
  if(password()!="ibdgc"){
    plot(0,0, xlim=c(1, 10), ylim=c(1,10), type="n", axes=F, xlab="", ylab="")
    
  }else{
    result <- getAssigment(i, all_snp, all_gene, GTeX, region)
    
    if(is.null(result$result)){
      plot(0,0, xlim=c(1, 10), ylim=c(1,10), type="n", axes=F, xlab="", ylab="")
      text(5,5, "no gene in this region")
      return()
    }
    
    result$snp
  }
}

plot_region <- function(choice, all_snp, all_gene, GTeX, region){
  
  i<- getRegion(choice)
  
  if(password()!="ibdgc"){
    plot(0,0, xlim=c(1, 10), ylim=c(1,10), type="n", axes=F, xlab="", ylab="")
  }else{
    plot(0,0, xlim=c(1, 10), ylim=c(1,10), type="n", axes=F, xlab="", ylab="")
    text(5,5, "hello")
    
    result <- getAssigment(i, all_snp, all_gene, GTeX, region)
    
    if(is.null(result$result)){
      plot(0,0, xlim=c(1, 10), ylim=c(1,10), type="n", axes=F, xlab="", ylab="")
      text(5,5, "no gene in this region")
      return()
    }
    
    temp <- as.matrix(result$result[,-c(1,2)]) %*% diag(c(weight()$coding, weight()$tfbs, weight()$roadmap, weight()$eQTL, weight()$other))
    nsig <- length(unique(result$result$signal))
    ngene <- dim(temp)[1]/nsig
    gene <- result$gene
    gene$score <- tapply(apply(temp, 1, max), list(rep(c(1:ngene), nsig)), sum )
    gene$keep <- getCredible(gene$score) < CI()
    if(max(gene$score) < 1E-5){
      gene$keep  <- rep(F, length(gene$score))
    }
    
    plot(0,0, xlim=range(result$region$start, result$region$end), ylim=c(-nsig,length(gene$name)), type="n", axes=F,  xlab=paste("chromosome", result$region$chr), ylab="")
    
    #plot gene
    y <- seq_along(gene$name)
    col <- rep("gray", length(gene$name))
    col[gene$keep] <- "black"
    arrows(gene$start[gene$strand=="+"], y[gene$strand=="+"], gene$end[gene$strand=="+"], y[gene$strand=="+"], lwd=3, col=col[gene$strand=="+"], length=0.1)
    arrows(gene$end[gene$strand=="-"],  y[gene$strand=="-"], gene$start[gene$strand=="-"], y[gene$strand=="-"], lwd=3, col=col[gene$strand=="-"], length=0.1)
    for(k in seq_along(gene$start)){
      start <- as.integer(strsplit(as.character(gene$exonStarts[k]), ",")[[1]])
      end <- as.integer(strsplit(as.character(gene$exonEnds[k]), ",")[[1]])
      for(l in seq_along((start))){
        rectarrows(start[l], y[k], end[l], y[k],height=0.5, length=0,  col=col[k], border=col[k])
      }
    }
    
    text(apply(cbind(gene$start, gene$end), 1,mean), y+1, gene$name, col=col)
    
    snp <- result$snp
    col <- rep("darkgray", length(snp$pos))
    col[snp$eQTL !=""] <- "darkgreen"
    col[snp$Roadmap !=""] <- "blue"
    col[snp$TFBS !=""] <- "orange"
    col[snp$Coding!=""] <- "red"
    pch <- rep(1, length(snp$pos))
    pch[snp$functional>0] <- 25
    cex <- exp(snp$P_mean_95)
    
    bg <- rep("white",  length(snp$pos))
    y_snp <- -as.integer(as.factor(result$snp$signal))
    points(snp$position, y_snp, pch=pch, col=col, bg=col, cex=cex)
    axis(1, at=seq(result$region$start, result$region$end, by=(result$region$end-result$region$start)/5),labels = format(round(seq(result$region$start, result$region$end, by=(result$region$end-result$region$start)/5)/1000000, 2), nsmall=2)  )
  }
}

plot_all_region <- function( all_snp, all_gene, GTeX, region, flag=F){
  
  if(password()!="ibdgc"){
    plot(0,0, xlim=c(1, 10), ylim=c(1,10), type="n", axes=F, xlab="", ylab="")
    
  }else{
    
    n_before <- c()
    n_after <- c()
    hd <- c()
    
    withProgress(message = 'Loading', value = 0, {
      
      for(i in unique(all_snp$HD)){
        incProgress(1/length(unique(all_snp$HD)), detail = paste("Region", i, "out of", max(all_snp$HD) ))
        
        result <- getAssigment(i, all_snp, all_gene, GTeX, region)
        if(is.null(result$result)){
          next()
        }
        
        temp <- as.matrix(result$result[,-c(1,2)]) %*% diag(c(weight()$coding, weight()$tfbs, weight()$roadmap, weight()$eQTL, weight()$other))
        nsig <- length(unique(result$result$signal))
        ngene <- dim(temp)[1]/nsig
        gene <- result$gene
        gene$score <- tapply(apply(temp, 1, max), list(rep(c(1:ngene), nsig)), sum )
        gene$keep <- getCredible(gene$score) < CI()
        n_after <- c(n_after, sum(gene$keep))
        n_before <- c(n_before, length(gene$name))
        hd <- c(hd, i)
      }
    })
    
    if(!flag){
      layout(matrix(c(1, 2), 1, 2, byrow = TRUE))
      barplot(t(cbind(GWAS=sum(n_before), FM=sum(n_after))), beside=T, log="y", main="Number of implicated genes")
      barplot(table(n_after), main="Number of implicated genes per region in FM")
    }else{
      data.frame(HD=hd[n_after==1], GWAS=n_before[n_after==1], FM=1)
    }
    
  }
}

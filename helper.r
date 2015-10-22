
HD.genelist <- read.table("dat/HD.genelist.txt", header=T, sep="\t")

getRegion <- function(choice){
  
  if(substr(choice, 1, 2)=="HD"){
    choice <- as.integer(substr(choice, 3, 10))
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

assignProximity <- function(gene, snp, GTeX){

  GTeX_sub <- GTeX[as.character(gene$name),]
  ret <- c()
  for(g in seq_along(gene$name)){
    if( (sum(row.names(GTeX_sub)==gene$name[g]) > 0 ) &  
        (
          (GTeX.options()$compareControl & min(GTeX_sub$GUT[g], GTeX_sub$IMMMUNE[g]) >= GTeX_sub$CONTROL[g]) | 
          (GTeX.options()$topQ & min(GTeX_sub$GUT[g], GTeX_sub$IMMMUNE[g]) > 5000) 
#           (F & min(GTeX_sub$GUT[g], GTeX_sub$IMMMUNE[g]) >= GTeX_sub$CONTROL[g]) | 
#           (F & min(GTeX_sub$GUT[g], GTeX_sub$IMMMUNE[g]) > 5000) 
        )){
        ret <- cbind(ret, rep(F, snp$position))
    }else{
      flanking_5 <- flanking()$flanking_5
      flanking_3 <- flanking()$flanking_3
      if(gene$stran[g] == "-"){
        temp <- flanking_3
        flanking_3 <- flanking_5
        flanking_5 <- temp
      }
      ret <- cbind(ret, sapply(snp$position, function(x){x > gene$start[g] - flanking_5 && x < gene$end[g] +flanking_3 }))
      #ret <- cbind(ret, sapply(snp$position, function(x){x > gene$start[g] -50000 && x < gene$end[g]+50000 }))
    }
  }
  
  colnames(ret) <- gene$name
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
      
      p0 <- weight()$trans
      p1 <- (1 -  weight()$trans) * weight()$proximity
      p2 <- (1-  weight()$trans) *(1- weight()$proximity)
      
      #p0 <- 0
      #p1 <- (1-p0) * 0.1
      #p2 <- (1-p0) * 0.9
      annotation <- subset(snp, select=c("Coding", "TFBS", "Roadmap", "eQTL")) != ""
      prob_mat <- c()
      for (k in seq_along(annotation[,1])){
        anno_k <- annotation[k,]
        if(sum(anno_k)>0){
          prob <- p2 / sum(anno_k)
          prob <- anno_k*prob 
          prob <- c(prob, p1, p0)
        }else{
          prob <- c(rep(0, length(anno_k)), p1+p2, p0)
        }
        prob_mat <- rbind(prob_mat, prob)
      }
      
      coding <- sapply(as.character(gene$name), grepl, snp$Coding)
      if(is.matrix(coding) ){
        if(dim(coding)[2]>1){
          coding <- t(sapply(seq_along(snp$variant), function(x){ ret <- coding[x,]; if(sum(coding[x,])>0){ ret <- coding[x,]/sum(coding[x,])}; ret} ))
        }
      }else{
        if(sum(coding)>0){
          coding <- coding / sum(coding)          
        }
      }
      
      eQTL <- sapply(as.character(gene$name), function(x){x %in% unique(gsub('.+\\((.+)\\)', '\\1', strsplit(as.character(subset(snp, eQTL!="")$eQTL[1]), ",")[[1]])) } ) 
      proximity <- assignProximity(gene, snp, GTeX)
      if(dim(proximity)[2]>1){
        proximity <- t(sapply(seq_along(proximity[,1]), function(x){ ret <- proximity[x,]; if(sum(proximity[x,])>0){ ret <- proximity[x,]/sum(proximity[x,])}; ret} ))
      }    
      score_snp <- coding * prob_mat[,1] +   proximity * prob_mat[,2] +  proximity * prob_mat[,3] +  proximity * prob_mat[,5] +  prob_mat[,4] %o% eQTL
      
      result <- rbind(result, 
                      data.frame(name=gene$name, 
                                 signal=j,
                                 coding =  as.vector(snp$P_mean_95 %*% (coding * prob_mat[,1])),
                                 TFBS =  as.vector(snp$P_mean_95 %*% ( proximity * prob_mat[,2])),
                                 roadmap =  as.vector(snp$P_mean_95 %*% (proximity * prob_mat[,3])),
                                 eQTL =  as.vector(snp$P_mean_95 %*% (prob_mat[,4] %o% eQTL)),
                                 proximity =  as.vector(snp$P_mean_95 %*% ( proximity * prob_mat[,5])),
                                 score=as.vector(snp$P_mean_95 %*% score_snp)
                                 )
                )
    }
  
  }
  list(region=list(start=region_start, end=region_end, chr=region$chr[i]),
       gene=gene, 
       snp= all_snp[paste(sep="", "HD", all_snp$HD)==as.character(region$region[i]), ], 
       result=result)
}

getAssigment_simple <- function(i, all_snp, all_gene, region){
  
  region_start <-region$start[i]-50000
  region_end <- region$end[i]+50000
  gene <- subset(all_gene, chr==region$chr[i]& ((start > region_start& start < region_end )|(end > region_start& end < region_end ) | (start < region_start & end > region_end) ) )
  
  result <- c()
  coding_hit <- c()
  if(dim(gene)[1] > 0){
    for (j in unique(all_snp$signal[paste(sep="", "HD", all_snp$HD)==as.character(region$region[i])]) ){
      snp <- all_snp[paste(sep="", "HD", all_snp$HD)==as.character(region$region[i]) & all_snp$signal==j, ]
      coding <- sapply(as.character(gene$name), grepl, snp$Coding)
      coding_hit <- rbind(coding_hit, snp$P_mean_95 %*% coding)
      
      r <- range(snp$position)
      
      start <- all_gene$start - flanking()$flanking_5
      end <- all_gene$end + flanking()$flanking_3
      start[all_gene$strand=="-"]  <- start[all_gene$strand=="-"] +flanking()$flanking_5 - flanking()$flanking_3
      end[all_gene$strand=="-"]  <- end[all_gene$strand=="-"] -flanking()$flanking_3 +flanking()$flanking_5

      result <- c(result, all_gene$name[ all_gene$chr==region$chr[i] & !(end < r[1] | start>r[2])] )
    }
  }
  
  result <- unique(result)
  if(!is.null(coding_hit)){
    result_coding <- colnames(coding_hit)[apply(coding_hit, 2, max) > 0.1]
    if(length(result)<=5 & length(result_coding)==1){
      result <- result_coding
    }
  }
  
  list(region=list(start=region_start, end=region_end, chr=region$chr[i]),
       gene=gene, 
       snp= all_snp[paste(sep="", "HD", all_snp$HD)==as.character(region$region[i]), ], 
       result=result)
}

plot_table <- function(choice, all_snp, all_gene, GTeX, region){
  
  i<- getRegion(choice)

  result <- getAssigment_simple(i, all_snp, all_gene, region)
  if(is.null(result$gene)){
    plot(0,0, xlim=c(1, 10), ylim=c(1,10), type="n", axes=F, xlab="", ylab="")
    text(5,5, "no gene in this region")
    return()
  }
  gene <- result$gene
  gene$keep <- result$gene$name %in% result$result
  temp <- subset(gene, select=c("name", "strand", "start", "end",  "keep"))
  temp$keep <- as.integer(temp$keep)
  names(temp) <- c("name", "strand", "start", "end", "keep")
  temp
}

plot_snp <- function(choice, all_snp, all_gene, GTeX, region){
  
  i<- getRegion(choice)

    result <- getAssigment_simple(i, all_snp, all_gene, region)
    
    if(is.null(result$result)){
      plot(0,0, xlim=c(1, 10), ylim=c(1,10), type="n", axes=F, xlab="", ylab="")
      text(5,5, "no gene in this region")
      return()
    }
    
  temp <- subset(result$snp, select=c("variant", "signal", "size", "trait.reassigned", "position", "exp_freq_a1", "P_mean_95", "Coding", "TFBS", "Roadmap", "eQTL", "functional"))
  names(temp) <- c("variant", "signal", "size", "trait", "position", "AF", "prob", "coding", "TFBS", "roadmap", "eQTL", "functional")
  temp$prob <- round(temp$prob, digits=3)
  temp
}

plot_region <- function(choice, all_snp, all_gene, GTeX, region){
  
  i<- getRegion(choice)
  
  if(length(i)==0){
    plot(0,0, xlim=c(1, 10), ylim=c(1,10), type="n", axes=F, xlab="", ylab="")
    text(5,5, "Please select a valid region")
    return()
  }
     
  result <- getAssigment_simple(i, all_snp, all_gene, region)
  if(dim(result$gene)[1]==0){
    plot(0,0, xlim=c(1, 10), ylim=c(1,10), type="n", axes=F, xlab="", ylab="")
    text(5,5, "no gene in this region")
    return()
  }
  gene <- result$gene
  gene$keep <- result$gene$name %in% result$result
  plot(0,0, xlim=range(result$region$start, result$region$end), ylim=c(-10,length(gene$name)+1), type="n", axes=F,  xlab=paste("Chromosome", result$region$chr , "(Mb)"), ylab="")
    
  #plot gene
  y <- seq_along(gene$name)
  col <- rep("gray", length(gene$name))
  col[gene$keep] <- "black"
  for(k in seq_along(gene$start)){
    start <- as.integer(strsplit(as.character(gene$exonStarts[k]), ",")[[1]])
    end <- as.integer(strsplit(as.character(gene$exonEnds[k]), ",")[[1]])
    for(l in seq_along((start))){
      rectarrows(start[l], y[k], end[l], y[k],height=0.5, length=0,  col=col[k], border=col[k])
    }
    
    if(gene$strand[k]=="+"){
      intron_start <- c(gene$start[k], end)
      intron_end <- c(start, gene$end[k])
    }else{
      intron_start <- c(gene$end[k], rev(start))
      intron_end <- c(rev(end), gene$start[k])
    }
    for(l in seq_along((intron_start))){
      arrows(intron_start[l],  y[k], intron_end[l] , y[k], lwd=3, col=col[k], length=0.1)
    }
    
  }
    
  text(apply(cbind(pmax(gene$start, result$region$start), pmin(gene$end, result$region$end)), 1,mean), y+1, gene$name, col=col)
    
  snp <- result$snp
  col <- rep("gray20", length(snp$pos))
  col[snp$eQTL !=""] <- "blue"
  col[snp$Roadmap !=""] <- "green"
  col[snp$TFBS !=""] <- "orange"
  col[snp$Coding!=""] <- "red"
  pch <- rep(1, length(snp$pos))
  pch[snp$functional>0] <- 25
  cex <- exp(snp$P_mean_95)
    
  bg <- rep("white",  length(snp$pos))
  y_snp <- -as.integer(as.factor(result$snp$signal))
    
  #points(snp$position, y_snp, pch=pch, col=col, bg=col, cex=cex)
  axis(1, at=seq(result$region$start, result$region$end, by=(result$region$end-result$region$start)/5),labels = format(round(seq(result$region$start, result$region$end, by=(result$region$end-result$region$start)/5)/1000000, 2), nsmall=2)  )
  
#  axis(2, at=y_snp, labels =as.integer(as.factor(result$snp$signal)), tick=F, line=-2)
  
  temp <- paste(sep="", as.integer(as.factor(result$snp$signal)), "-", as.character(result$snp$trait.reassigned), "(", sprintf(fmt="%.1f", result$snp$p_multi) ,")")
  temp <- tapply(temp, y_snp, function(x){x[1]})
  axis(4, at=names(temp), labels = temp, tick=F, las=2, line=-5)
  #abline(h=y_snp, col="gray", lty="dotted")
  temp <- tapply(snp$position, y_snp, range)
  rect(do.call("rbind", temp)[,1], as.integer(names(temp))-0.1, do.call("rbind", temp)[,2], as.integer(names(temp))+0.1, col="gray20", border=NA)
  
  len <- exp(snp$P_mean_95) * 0.15 
  segments(snp$position, y_snp-len, snp$position, y_snp+len, lwd=3, col=col)
  
  #mtext("Signal", 2, at=mean(y_snp), padj=-1)  
  #plot(0,0, xlim=c(1,10), ylim=c(1,10), type="n", axes=F, xlab="", ylab="")
  #legend("bottomright", c("Coding", "TFBS", "Epigenetic", "eQTL", "No function assigned"), col=c("red", "orange", "green", "blue", "darkgray"), pch=c(rep(25, 4), 1), pt.bg=c("red", "orange", "green", "blue", "gray"))
  
}

plotSNP <- function(plot_click, choice, all_snp, all_gene, GTeX, region){
  
  i<- getRegion(choice)
  
  if(length(i)==0){
    return()
  }
  
  result <- getAssigment_simple(i, all_snp, all_gene, region)

  set_temp <- (result$snp)[, -c(2:4, 6)]
  set_temp$signal <- -set_temp$signal
  set_temp$variant <- as.character(set_temp$variant)
  
  temp <-   nearPoints(set_temp, plot_click, xvar = "position", yvar="signal", maxpoints = 1)
  paste("variant: ", temp$variant, 
        "\nsignal: ", -temp$signal,
        "\nposition: ", temp$position,
        "\ntrait: ", temp$trait.reassigned,
        "\nallele frequency: ", temp$exp_freq_a1,
        "\nProbability:", sprintf("%.4g", temp$P_mean_95),
        "\ncoding: ", temp$Coding,
        "\nTFBS:", temp$TFBS,
        "\nEpigenetics:", temp$Roadmap,
        "\neQTL:", temp$eQTL
  )
}

plot_all_region <- function( all_snp, all_gene, GTeX, region, flag=F){
  
    n_before <- c()
    n_after <- c()
    n_before_gene <- c()
    n_after_gene <- c()
    hd <- c()
    
    withProgress(message = 'Loading', value = 0, {

    for(i in unique(all_snp$HD)){
      incProgress(1/length(unique(all_snp$HD)), detail = paste("Region", i, "out of", max(all_snp$HD) ))
      
      if(F){
        result <- getAssigment(i, all_snp, all_gene, GTeX, region)
        
        if(is.null(result$result)){
          next()
        }
        nsig <- length(unique(result$result$signal))
        ngene <- dim(result$result)[1]/nsig
        gene <- result$gene
        gene$score <- tapply(result$result$score, rep(c(1:ngene), nsig), mean )
        gene$keep <- getCredible(gene$score) < CI()
        if(max(gene$score) < 1E-5){
          gene$keep  <- rep(F, length(gene$score))
        }
      }else{
        
        result <- getAssigment_simple(i, all_snp, all_gene, region)
        if(dim(result$gene)[1]==0){
          next()
        }
        gene <- result$gene
        gene$keep <- result$gene$name %in% result$result
        
      }
      
      n_after <- c(n_after, sum(gene$keep))
      n_after_gene <- c(n_after_gene, paste(gene$name[gene$keep], collapse = ","))
      
      start <- all_gene$start - flanking()$flanking_5
      end <- all_gene$end + flanking()$flanking_3
      start[all_gene$strand=="-"]  <- start[all_gene$strand=="-"] +flanking()$flanking_5 - flanking()$flanking_3
      end[all_gene$strand=="-"]  <- end[all_gene$strand=="-"] -flanking()$flanking_3 +flanking()$flanking_5
      
      
      
      n_before <- c(n_before, length(all_gene$name[ all_gene$chr==region$chr[i] & !(end <  R2.span$start[R2.span$HD==i][1] | start> R2.span$end[R2.span$HD==i][1] )] ))
      n_before_gene <- c(n_before_gene, paste(all_gene$name[ all_gene$chr==region$chr[i] & !(end <  R2.span$start[R2.span$HD==i][1] | start> R2.span$end[R2.span$HD==i][1] )], collapse = ","))
      
      hd <- c(hd, i)
    }
    })
  
    if(!flag){
      layout(matrix(c(1, 2,2), 1, 3, byrow = TRUE))
      
      temp <- barplot(cbind(GWAS=sum(n_before), FM=sum(n_after)), ylab="Number of implicated genes", ylim=c(1,max(sum(n_before), sum(n_after))*1.1))
      
      text(temp[1], sum(n_before)*1.05, sum(n_before), cex=1.5)
      text(temp[2],  sum(n_after)*1.05, sum(n_after), cex=1.5)
      
      barplot(table(n_after), ylab="Number of regions", xlab="Number of implicated genes per region (in FM)")
    }else{
      data.frame(HD=hd, GWAS=n_before, FM=n_after, gene_before=n_before_gene, gene_after=n_after_gene)
    }
}

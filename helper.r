
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

getAssigment <- function(i, all_snp, all_gene, region){
  
  region_start <-region$start[i]-50000
  region_end <- region$end[i]+50000
  gene <- subset(all_gene, chr==region$chr[i]& ((start > region_start& start < region_end )|(end > region_start& end < region_end ) | (start < region_start & end > region_end) ) )
  
  result <- c()
  coding_hit <- c()
  if(dim(gene)[1] > 0){
    for (j in unique(all_snp$signal[all_snp$HD==region$HD[i]]) ){
      snp <- all_snp[all_snp$HD==region$HD[i] & all_snp$signal==j, ]
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
    result_coding <- colnames(coding_hit)[apply(coding_hit, 2, max) > 0.5]
    if(length(result)<=999 & length(result_coding)==1){
      result <- result_coding
    }
  }
  
  list(region=list(start=region_start, end=region_end, chr=region$chr[i]),
       gene=gene, 
       snp= all_snp[all_snp$HD==region$HD[i], ], 
       result=result)
}

plotSNP <- function(plot_click, plot_brush, choice, all_snp, all_gene, region){
 
  i<- getRegion(choice)
  
  if(length(i)==0){
    return()
  }
  
  result <- getAssigment(i, all_snp, all_gene, region)
  result <- subset(result$snp, select=c("variant", "signal", "size", "trait", "position", "exp_freq_a1", "P_mean_95","A0", "A1", "logOR_CD", "logOR_UC", "Coding", "TFBS", "Roadmap", "eQTL"))
  
  names(result) <- c("variant", "signal", "size", "trait", "position", "AF", "prob","A0", "A1", "logOR_CD", "logOR_UC", "coding", "TFBS", "roadmap", "eQTL")
  
  result$prob <- round(result$prob, digits=3)
  result$logOR_CD <- round(result$logOR_CD, digits=3)
  result$logOR_UC <- round(result$logOR_UC, digits=3)
  result$variant <- as.character(result$variant)
  result$signal <-  -as.integer(as.factor(result$signal))
  
  selected <- c()
  if(!is.null(plot_brush)){
    selected <- brushedPoints(result, plot_brush, xvar = "position", yvar="signal")
  }else if(!is.null(plot_click)) {
    selected <-  nearPoints(result, plot_click, xvar = "position", yvar="signal")
  }
  
  if(!is.null(selected) && dim(selected)[1]>0){
    selected$signal <- -selected$signal
    selected <- selected[with(selected, order(signal,  -prob)), ]
  }
  selected
}

plot_region <- function(choice, all_snp, all_gene, region){
  
  i<- getRegion(choice)
  
  if(length(i)==0){
    plot(0,0, xlim=c(1, 10), ylim=c(1,10), type="n", axes=F, xlab="", ylab="")
    text(5,5, "Please select a valid region")
    return()
  }
  par(mar=c(5,0,1,5))
  result <- getAssigment(i, all_snp, all_gene, region)
  
  if(dim(result$gene)[1]>0){
  
  
    gene <- result$gene
    gene <- gene[with(gene, order(start)), ]
    gene$keep <- gene$name %in% result$result
    plot(0,0, xlim=range(result$region$start, result$region$end), ylim=c(-10, length(gene$name)+1), type="n", axes=F,  xlab=paste("Chromosome", result$region$chr , "(Mb)"), ylab="", main=paste(sep="", "HD", i))
#    rect(result$region$start, -10, result$region$end, 0 , col="gray90", border=NA)

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
      if(abs(intron_end[l]- intron_start[l]) > 10){
        suppressWarnings(arrows(intron_start[l],  y[k], intron_end[l] , y[k], lwd=3, col=col[k], length=0.1))
      }
    }
    
  }
#  text(apply(cbind(pmax(gene$start, result$region$start), pmin(gene$end, result$region$end)), 1,mean), y+1, gene$name, col=col)
    
    ii_sel <- gene$start>  result$region$start
    if(sum(ii_sel)>0){
    text(pmax(gene$start, result$region$start)[ii_sel], y[ii_sel], adj=1.1, gene$name[ii_sel], col=col[ii_sel])
  }
    if(sum(!ii_sel)>0){
    text(gene$end[!ii_sel], y[!ii_sel], adj=-0.1, gene$name[!ii_sel], col=col[!ii_sel])
  }
  }else{
    plot(0,0, xlim=range(result$region$start, result$region$end), ylim=c(-10, 1), type="n", axes=F,  xlab=paste("Chromosome", result$region$chr , "(Mb)"), ylab="", main=paste(sep="", "HD", i))
#    rect(result$region$start, -10, result$region$end, 0 , col="gray90", border=NA)
  }

  #plot SNPs
  snp <- result$snp
  col <- rep("#808080", length(snp$pos))
  col[snp$eQTL !=""] <- "blue"
  col[snp$Roadmap !=""] <- "green"
  col[snp$TFBS !=""] <- "orange"
  col[snp$Coding!=""] <- "red"

  y_snp <- -as.integer(as.factor(result$snp$signal))
  
  axis(1, at=seq(result$region$start, result$region$end, by=(result$region$end-result$region$start)/5),labels = format(round(seq(result$region$start, result$region$end, by=(result$region$end-result$region$start)/5)/1000000, 2), nsmall=2)  )
  
  #  axis(2, at=y_snp, labels =as.integer(as.factor(result$snp$signal)), tick=F, line=-2)
  
  temp <- paste(sep="", as.integer(as.factor(result$snp$signal)), "-", as.character(result$snp$trait), "(", sprintf(fmt="%.1f", result$snp$p_multi) ,")")
  temp <- tapply(temp, y_snp, function(x){x[1]})
  axis(4, at=names(temp), labels = temp, tick=F, las=2, line=-1)
  abline(h=0, col="gray", lty="dotted")
  temp <- tapply(snp$position, y_snp, range)
  
  rect(do.call("rbind", temp)[,1] -50000, as.integer(names(temp))-0.3, do.call("rbind", temp)[,2]+50000, as.integer(names(temp))+0.3, col="gray95", border="black", lty = "dotted")
  
  len <- exp(snp$P_mean_95) * 0.15 
  segments(snp$position, y_snp-len, snp$position, y_snp+len, lwd=3, col=col)

  mtext("Signal-Trait\n(-log10P)", 4, line=0, at=1, las=1)  
}

plot_table <- function(choice, all_snp, all_gene, region){
  
  i<- getRegion(choice)

  result <- getAssigment(i, all_snp, all_gene, region)
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

plot_all_region <- function( all_snp, all_gene, region){
  
    n_before <- c()
    n_after <- c()
    n_before_gene <- c()
    n_after_gene <- c()
    hd <- c()
    
    withProgress(message = 'Loading', value = 0, {

    for(i in unique(all_snp$HD)){
      incProgress(1/length(unique(all_snp$HD)), detail = paste("Region", i, "out of", max(all_snp$HD) ))
      
      result <- getAssigment(i, all_snp, all_gene, region)
      if(dim(result$gene)[1]==0){

        n_after <- c(n_after, 0)
        n_after_gene <- c(n_after_gene, "")

      }else{
        gene <- result$gene
        gene$keep <- result$gene$name %in% result$result
  
        n_after <- c(n_after, sum(gene$keep))
        n_after_gene <- c(n_after_gene, paste(gene$name[gene$keep], collapse = ", "))
      
        start <- all_gene$start - flanking()$flanking_5
        end <- all_gene$end + flanking()$flanking_3
        start[all_gene$strand=="-"]  <- start[all_gene$strand=="-"] +flanking()$flanking_5 - flanking()$flanking_3
        end[all_gene$strand=="-"]  <- end[all_gene$strand=="-"] -flanking()$flanking_3 +flanking()$flanking_5
      }
      
      n_before <- c(n_before, length(all_gene$name[ all_gene$chr==region$chr[i] & !(end <  region$start[region$HD==i][1] | start> region$end[region$HD==i][1] )] ))
      n_before_gene <- c(n_before_gene, paste(all_gene$name[ all_gene$chr==region$chr[i] & !(end <  region$start[region$HD==i][1] | start> region$end[region$HD==i][1] )], collapse = ", "))
      hd <- c(hd, i)
      
    }
    })
  
    data.frame(HD=hd, GWAS=n_before, FM=n_after, gene_before=n_before_gene, gene_after=n_after_gene)
}

printLegend <- function(choice){
  
  
  i<- getRegion(choice)
  
  ret <- c()
  
  if(length(i)>0){
    ret <- tags$html(
    tags$div(style="width:20%", class="legend", downloadButton("downloadPDF", "Save figure")),
    tags$div(style="width: 1%", class="legend", span(class="bar", style ="background-color:red;")),
    tags$div(style="width:14%", class="legend", "Coding"), 
    tags$div(style="width: 1%", class="legend", span(class="bar", style ="background-color:orange")),
    tags$div(style="width:14%", class="legend", "TFBS"), 
    tags$div(style="width: 1%", class="legend", span(class="bar", style ="background-color:green")),
    tags$div(style="width:14%", class="legend", "Epigenetic"), 
    tags$div(style="width: 1%", class="legend", span(class="bar", style = "background-color:blue")),
    tags$div(style="width:14%", class="legend", "eQTL"), 
    tags$div(style="width: 1%", class="legend", span(class="bar", style = "background-color:#808080;")),
    tags$div(style="width:14%", style="clear:right", class="legend", "No_function")
  )
  }
  
  ret
  
}

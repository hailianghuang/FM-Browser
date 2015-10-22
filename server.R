source("helper.r")

all_gene <- read.csv("dat/refGene.txt", header=T, sep="\t")
names(all_gene) <- c("chr", "strand", "start", "end", "cdsStart", "cdsEnd", "exonCount", "exonStarts", "exonEnds", "score", "name")
temp <-with(all_gene, substr(chr, 4, 100))
temp[temp=="X"] <- 23
all_gene$chr <- as.integer(temp)
all_gene$name <- as.character(all_gene$name)

HD.genelist <- read.table("dat/HD.genelist.txt", header=T, sep="\t")

GTeX <- read.table("dat/GTeX.rank.summary.txt", header=T, row.names=1)
all_snp <- read.table("dat/annotation.variant.txt", header=T, sep="\t")
all_snp <- subset(all_snp, tier2=="No")
all_snp$functional <- apply(all_snp[, c("Coding", "TFBS", "Roadmap", "eQTL")], 1, function(x){sum(as.character(x)!="")>0})
region <- read.table("dat/HD_region_HG19.txt", header=F)
names(region) <- c("chr", "start", "end", "region")

R2.span <- read.table("dat/R2.span.v2.txt", header=T)

shinyServer(
  function(input, output) {

    password <<- reactive({
      as.character(input$password)
    })
    
    flanking <<- reactive({
      ret <- list(
        flanking_5 = input$flanking_5*1000,
        flanking_3 = input$flanking_3*1000
      )
      ret
    })
    
    CI  <<- reactive({
      as.numeric(input$CI)/100
    })
    
    weight <<- reactive({
      ret <- list(
        trans = input$trans/100,
        proximity = input$proximity/100
      )
      ret
    })
    
    GTeX.options <<- reactive({
      ret <- list(
        compareControl = input$compareControl,
        topQ = input$topQ
      )
      ret
    })
    
    output$plot <- renderPlot({
      plot_region(input$HD, all_snp, all_gene, GTeX, region)
    })
    
    output$downloadPDF <- downloadHandler(
      filename = tempfile(fileext=".pdf"),
      content = function(file){

        pdf(file, width=8, height=6)
        plot_region(input$HD, all_snp, all_gene, GTeX, region)
        dev.off()

      },
      contentType = "application/pdf"
    )
    
    output$table <- renderDataTable({
      plot_table(input$HD, all_snp, all_gene, GTeX, region)
    })

    output$snp <- renderDataTable({
      plot_snp(input$HD, all_snp, all_gene, GTeX, region)
    })

    output$plot_all_region <- renderPlot({
      plot_all_region(all_snp, all_gene, GTeX, region)
    })  
    
    output$all_region <- renderDataTable({
      plot_all_region(all_snp, all_gene, GTeX, region, T)
    })
    
    output$info <- renderText({
      plotSNP(input$plot_click, input$HD, all_snp, all_gene, GTeX, region)
    })
    output$legend <- renderPlot({
      plot(0,0, type="n")
      legend("bottomright", c("Coding", "TFBS", "Epigenetic", "eQTL", "No function assigned"), col=c("red", "orange", "green", "blue", "darkgray"), pch=c(rep(25, 4), 1), pt.bg=c("red", "orange", "green", "blue", "gray"))
    })
    
  }
)
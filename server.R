source("helper.r", local=TRUE)

all_gene <- read.csv("dat/refGene.txt", header=T, sep="\t")
names(all_gene) <- c("chr", "strand", "start", "end", "cdsStart", "cdsEnd", "exonCount", "exonStarts", "exonEnds", "score", "name")
temp <-with(all_gene, substr(chr, 4, 100))
temp[temp=="X"] <- 23
all_gene$chr <- as.integer(temp)
all_gene$name <- as.character(all_gene$name)

region <- read.table("dat/HD_region_HG19.txt", header=F)
names(region) <- c("chr", "start", "end", "HD")
region$HD <-  as.integer(gsub("HD", "", as.character(region$HD)))

shinyServer(
  function(input, output) {

    flanking <<- reactive({
      ret <- list(
        flanking_5 = input$flanking_5*1000,
        flanking_3 = input$flanking_3*1000
      )
      ret
    })
    
    output$plot <- renderPlot({
      plot_region(input$HD, all_snp, all_gene, region)
    })
    
    output$downloadPDF <- downloadHandler(
      filename = function() {
        paste('region-', input$HD, "-", Sys.Date(),  '.pdf', sep='')
      },
      content = function(file){
        pdf(file, width=8, height=6)
        plot_region(input$HD, all_snp, all_gene, region)
        dev.off()

      },
      contentType = "application/pdf"
    )
    
    output$gene <- renderDataTable({
      plot_table(input$HD, all_snp, all_gene, region)
    }, options = list(
      pageLength = 100, 
      paging = FALSE)
    )

    output$all_region <- renderDataTable({
      plot_all_region(all_snp, all_gene, region)
    }, options = list(
      pageLength = 100, 
      paging = FALSE)
    )
    
    output$info <-  renderTable({
        plotSNP(input$plot_click, input$plot_brush, input$HD, all_snp, all_gene, region)
    },     caption = "<h4><center> Selected variants </center></h4>",
    caption.placement = getOption("xtable.caption.placement", "top"), 
    caption.color = getOption("xtable.caption.color", "black"), 
    caption.width = getOption("xtable.caption.width", NULL), include.rownames=FALSE)
   
    output$legend <- renderUI({
      printLegend(input$HD)
    })
    
  }
)
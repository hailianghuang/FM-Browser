
HD.genelist <- read.table("dat/HD.genelist.txt", header=T, sep="\t")

choices <- paste(sep="", "HD", unique(HD.genelist$HD))
choices <- c(choices, unique(as.character(HD.genelist$gene)))

shinyUI(fluidPage(
  titlePanel("International IBD Genetics Consortium Fine-mapping project"),
  
  sidebarLayout(
    sidebarPanel(

      selectInput("HD", 
                  label = "Choose a region to display",
                  choices = c("", choices),
                  selected = "", selectize=TRUE),
      sliderInput("flanking_5", 
                  label = "5' flanking (Kbp)",
                  min = 0, max = 100, value =50),
      sliderInput("flanking_3", 
                  label = "3' flanking (Kbp)",
                  min = 0, max = 100, value =50),
      selectInput("CI", 
                  label = "Confidence level (%)",
                  choices = rev(c(50, 80, 90, 95, 99, 99.9)),
                  selected = 95),
      tags$a(href = "http://biorxiv.org/content/early/2015/10/20/028688", "Manuscript and data available at bioRxiv")
#      hr(),
#      checkboxInput(inputId = "weight",
#                    label =strong("Customize weight"),
#                    value = FALSE),
#      conditionalPanel(
#        condition = "input.weight == true ",
#          sliderInput("trans", "Probability for trans-effect (%)", min=0, max=99, value=0),
#          sliderInput("proximity", "Probability for non-assigned function (%)", min=0, max=99, value=10)
#        ),
#      hr(),
#      checkboxInput(inputId = "GTeX",
#                            label = strong("Integrate GTeX"),
#                            value = FALSE),
#      conditionalPanel(
#        condition = "input.GTeX == true ",
#        checkboxInput(inputId = "compareControl",
#                      label = "Better than control",
#                      value = FALSE),
#        checkboxInput(inputId = "topQ",
#                      label = "In top 25% quantile",
#                      value = FALSE)
#      )
      
    ),
    
    mainPanel(
      tabsetPanel(
        tabPanel("Plot",   plotOutput("plot", click = "plot_click"), downloadButton("downloadPDF", "Download figure"), verbatimTextOutput("info"), plotOutput("legend")), 
        tabPanel("Table", dataTableOutput("table")),
        tabPanel("SNP score", dataTableOutput("snp")),
        tabPanel("All region", dataTableOutput("all_region")),
        tabPanel("Plot all region", plotOutput("plot_all_region"))
      )
    )
  )
))

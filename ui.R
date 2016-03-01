
choices <- paste(sep="", "HD", unique(all_snp$HD))
choices <- c(choices, unique(as.character(HD.genelist$gene)))

shinyUI(fluidPage(
  tags$head(
    tags$style(HTML("
table, th, td {
    border: 1px solid black;
                    }
                    .legend {
                      float:left; 
                      height:36px; 
                      line-height: 36px;
                      vertical-align:middle;
                      padding:5px;
                      margin: 0px
                    }
                    span.bar{
                      width:.3rem; 
                      height:36px; 
                      display:block;
                    }

                    "))
    ),
  
  titlePanel("International IBD Genetics Consortium Fine-mapping project"),
  
  tags$div(id="input_container", 
           
    sidebarLayout(
             
    sidebarPanel(

      selectInput("HD", 
                  label = "Choose a region/gene to display",
                  choices = c("", choices),
                  selected = "", selectize=TRUE),
      sliderInput("flanking_5", 
                  label = "5' flanking (Kbp)",
                  min = 0, max = 100, value =50),
      sliderInput("flanking_3", 
                  label = "3' flanking (Kbp)",
                  min = 0, max = 100, value =50),
      tags$a(href = "http://biorxiv.org/content/early/2015/10/20/028688", "Manuscript and results (bioRxiv)")
    , width = 3),
    # load Javascript snippet to parse the query string.
    mainPanel(
      tabsetPanel(
        tabPanel("Regional plot",   plotOutput("plot", click = "plot_click", brush = "plot_brush"), uiOutput("legend"), tags$div(tableOutput("info"))), 
        tabPanel("Gene list", dataTableOutput("gene")),
        tabPanel("Region list", dataTableOutput("all_region"))
      )
    , width=9)
    
    ),
    singleton(tags$script(type="text/javascript", src="js/parse_input.js"))
    
  )
  
))

#Nolan App: Wed 29 January
library(shiny)
library(knitr)

shinyUI(pageWithSidebar(
  headerPanel("Various measures that quantify the contribution of disease risk variants"), 
  
    sidebarPanel(
      
    h3("Input Options"),
    
    wellPanel(
      h4("Input Variables"),
      htmlOutput("inputs"),
      radioButtons("my_method", h4("Choose one:"),
                   list("Sibling relative risk (Lambda_S)" = "risk",
                        "Heritability of liability (h2L)" = "heritability"))
      ),
    
    h3("Upload Data"),
    
    wellPanel(
    helpText("Upload data in the same format as the sample file, downloadable below."),
    a("Sample Data", href="https://www.dropbox.com/sh/bvebwy1fg1mfj86/l6hnK5QYeZ/Sample_Input.csv"),     
    tags$hr(),
    fileInput('file1', 'Choose CSV File from local drive, adjusting parameters if necessary',
              accept=c('text/csv', 'text/comma-separated-values,text/plain')),
    
    checkboxInput('header', 'Header', TRUE),
    radioButtons('sep', 'Separator',
                 c(Comma=',',
                   Semicolon=';',
                   Tab='\t'),
                 'Comma'),
    radioButtons('quote', 'Quote',
                 c(None='',
                   'Double Quote'='"',
                   'Single Quote'="'"),
                 'Double Quote'),
    tags$head(tags$style(type="text/css",
                         "label.radio { display: inline-block; margin:0 10 0 0;  }",
                         ".radio input[type=\"radio\"] { float: none; }"))
        )
    ),
  
  mainPanel(
    tabsetPanel(
      
      tabPanel("Explanation",
               wellPanel(h4("Measures of the Contribution of Independent Risk Variants to Disease"),
                         htmlOutput("explanation")),
               wellPanel(h4("Input Variables"),
               htmlOutput("inputexplanation")),
               wellPanel(h4("Output Variables"),
               htmlOutput("outputexplanation")),
               wellPanel(
                 HTML(knit2html(text="References are from $Witte$ et. al (20xx)")) #Important. Keep this formatting in. For some reason xtable formatting only works if this HTML(knitr) function is used on this page. weird.
               )),
      
      tabPanel("Table",
               h4(textOutput("tabletitle")),
               downloadButton('downloadData', 'Download Table'),
               tags$hr(),
               tableOutput("raw")),
          
      tabPanel("Plot",
               h4(textOutput("plottitle")),
               plotOutput("plot"),
               htmlOutput("plotexplanation")),
                     
      tabPanel("Citation",
               h4("Citation"),
               p(a(href="http://www.complextraitgenomics.com/publications","The Contribution of Genetic Variants to Disease Depends on the Ruler")),
               p("John S. Witte, Peter M. Visscher, Naomi R. Wray"),
               textOutput("cite"))
      
  )
  )))
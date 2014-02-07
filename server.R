#Nolan App: Wed 29 Jan
library(shiny)
library(knitr)
library(xtable)

source('functions.R')
source('plotfunctions.R')

shinyServer(function(input, output) {
  
###Input Data and Calculations###
  Data <- reactive({
    
    inFile <- input$file1
    
    if (is.null(inFile))
      return(NULL)
   
    #Column names:
    name <- c("RiskVariant",
              "RAF",
              "RRBb",
              "RRBB",
              "K",
              "Lambda_S",
              "h2Li",
              "h2Liprop",
              "h2Lipropapprox",
              "lambdasib",
              "proploglambdasib",
              "proplambS2",
              "AUCi",
              "AUCprop",
              "PAR")
    
    df.raw <- read.csv(inFile$datapath, header=input$header, sep=input$sep, quote=input$quote, col.names=name)
   
    #Calculations#
    
    if(input$my_method=="risk"){
      calculation <- get_Vgs_lambdas(df.raw[,1], df.raw[,2], df.raw[,3], df.raw[,4], df.raw[,5],df.raw[,6])
    }
    if(input$my_method=="heritability"){
      calculation <- get_Vgs_h2L(df.raw[,1], df.raw[,2], df.raw[,3], df.raw[,4], df.raw[,5],df.raw[,6])
    }
        
    df.raw$h2Li <- calculation$h2Li
    df.raw$h2Liprop <- calculation$h2Liprop
    df.raw$h2Lipropapprox <- calculation$h2Lipropapprox
    df.raw$lambdasib <- calculation$lambdasib
    df.raw$proploglambdasib <- calculation$proploglambdasib
    df.raw$proplambS2 <- calculation$proplambS2
    df.raw$AUCi <- calculation$AUCi
    df.raw$AUCprop <- calculation$AUCprop
    df.raw$PAR <- calculation$PAR
    
        
    # create a list of data for use in rendering
    info <- list(df.raw=df.raw)
    return(info)
  })
  
### Input Tab ###
  output$inputs <- renderText(
    paste(print(xtable(inputs, include.rownames = FALSE), include.colnames = FALSE, type = "html", html.table.attributes = c("class=table-condensed"), sanitize.text.function=identity, print.results = FALSE),
          tags$script("MathJax.Hub.Queue([\"Typeset\",MathJax.Hub]);"))
  )  
  
### Output: Explanation Tab ###
   output$explanation <-  renderUI( {
    HTML("<p>The contribution of identified independent genetic variants to disease risk 
can be assessed by a number of measures that give values between 0 and 100%. These
different measures have been developed with different goals and from traditionally 
disparate fields, such as epidemiology or quantitative genetics. This calculator allows 
researchers to contrast the most commonly used measures:</p>
         <ul>
         <li> heritability explained, </li>
        <li> sibling relative risk explained, </li>
        <li> familial relative risk explained, </li>
         <li> area under the receiver- operating curve (AUC), </li>
         <li> population attributable risk (PAR). </li>
         </ul>
         <p> These measures can give materially different messages about the impact of 
risk variants on disease. Here we present side-by-side the different models and measures
used to assess the impact of genetic variants on disease.</p>
<p> Key factors that contribute to the differences between measure are the 
importance given to the frequency of a risk variant relative to its effect 
size on disease, and the baseline to which importance is expressed. These factors 
should be explicitly considered when assessing the contribution of a genetic variant to 
complex traits.</p>")
})
  
  
  output$inputexplanation <- renderText(
    paste(print(xtable(parameters[1:6,], include.rownames = FALSE), include.colnames = FALSE, type = "html", html.table.attributes = c("class=table-condensed"), sanitize.text.function=identity, print.results = FALSE),
          tags$script("MathJax.Hub.Queue([\"Typeset\",MathJax.Hub]);"))
  )
  
  output$outputexplanation <- renderText(
    paste(print(xtable(parameters[7:15,], include.rownames = FALSE), include.colnames = FALSE, type = "html", html.table.attributes = c("class=table-condensed"), sanitize.text.function=identity, print.results = FALSE),
          tags$script("MathJax.Hub.Queue([\"Typeset\",MathJax.Hub]);"))
  )
  
### Output: Table Tab ###
  output$raw <- renderTable({
    if (is.null(input$file1)) { return() }
    
    Data()$df.raw       
  })
  
  output$downloadData <- downloadHandler(
    filename = function() { paste("shiny_table", Sys.Date(), '.csv', sep='') },
    content = function(file) {
       write.csv(Data()$df.raw, row.names=F, file)
    })
  
  output$tabletitle <- renderText( {
    if (is.null(input$file1)) { return() }
        "Measures of the Contribution of Independent Risk Variants to Disease"
  })
  
### Output: Plot Tab ###  
  
  PlotData <- reactive({
    inFile <- input$file1
    
    if (is.null(inFile))
      return(NULL)
    
    df.plot <- read.csv(inFile$datapath, header=input$header, sep=input$sep, quote=input$quote)
    
    #Calculations#
    if(input$my_method=="risk"){
      plot_calculation <- get_plot_lambdas(df.plot[,1], df.plot[,2], df.plot[,3], df.plot[,4], df.plot[,5],df.plot[,6])
    }
    if(input$my_method=="heritability"){
      plot_calculation <- get_plot_h2L(df.plot[,1], df.plot[,2], df.plot[,3], df.plot[,4], df.plot[,5],df.plot[,6])
    }
   
  })
  
  output$plottitle <- renderText( {
    if (is.null(input$file1)) { return() }
    "Measures of the Contribution of Independent Risk Variants to Disease"
  })
  
  output$plotexplanation <- renderUI( {
    if (is.null(input$file1)) { return() }
    HTML("<p> Comparison of commonly used measures for assessing the impact of known 
risk variants on disease. The measures are: </p>
<ul>
<li> heritability explained; </li>
<li> approximation of heritability explained; </li>
<li> sibling relative risk explained; </li>
<li> familial relative risk explained (i.e., using log relative risks); and </li>
<li> the proportion based on the area under the curve (AUC) statistic. </li>
</ul>
<p> Each line corresponds to an individual risk variant, indicating the percentage 
of each measure (e.g., total variability) it explains. Lines are different colours 
depending on the relative risk (estimated by the odds ratio, OR) for each variant. 
The total contribution from all inputted variants (assumed to be independent) is given
in parenthesis below the x-axis. Details and formulas for these measures are given in 
the accompanying paper. The y-axis has a squared scale.</p>") })
  
  output$plot <- renderPlot({
    if (is.null(input$file1)) { return() }
    print(PlotData()$plot_calculation)
  })
    
### Output: Citation Tab ###
  output$cite <- renderText({
   HTML("App designed by Cara Nolan and Beben Benyamin")
  })
  
})
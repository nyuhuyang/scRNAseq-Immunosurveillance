library(shiny)

shinyUI(fluidPage(

  titlePanel("Single Cell RNAseq"),
  sidebarLayout(
    sidebarPanel(

      radioButtons(inputId ="Dataset1",
                  label=h3("Data Set"),choices=as.list(samples),                  
                  selected=samples[1]),
      textInput(inputId ="text", label = h3("Input Gene Name"), value = "Enter SINGLE gene only..."),
      actionButton(inputId = "update", "Update Plot !"),
      p("Click the above button to Update Plot After Input Gene Name (Required), e.g. CD8A (Human), Cd8a (Mouse)"),
      textInput(inputId ="threshold", label = h3("low expression threshold"), value = 0),
      selectInput(inputId ="Color1",
                  label=h3("Select a Color for high expression"),choices=list("red", "blue", "green", "purple", "black"),                  
                  selected="red"),
      selectInput(inputId ="dotsize",
                  label=h3("Select a dot size"),choices=list(1,2,3,4,5,6,7),                  
                  selected=2),
      downloadButton('downloadplot', 'Download Plot')
    ),
   # mainPanel(
    #  plotOutput("distPlot",height = "1000px", width = "1000px")
    #)
    mainPanel(
              fluidRow(
                splitLayout(cellWidths = c("100%", "40%"), 
                            plotOutput("distPlot",height = "800px", width = "800px", 
                                       hover = "plot_hover",click = "plot_click",
                                       dblclick = "plot_dblclick")
            #, plotOutput("distPlot2",height = "700px", width = "700px",hover = "plot_hover2",click = "plot_click2",
             #            dblclick = "plot_dblclick2",brush = "plot_brush2")
                            )
        #column(6,plotOutput(outputId="distPlot", width="700px",height="700px")),  
        #column(6,plotOutput(outputId="distPlot2", width="700px",height="700px"))
                )
              ,
              fluidRow(
                column(6,verbatimTextOutput("info3")),
                column(6,verbatimTextOutput("info4"))
              #fluidRow(h4("Plot 1 Values Displayed"),verbatimTextOutput("info")),
              #fluidRow(h4("Plot 2 Values Displayed"),verbatimTextOutput("info2"))
    )
  )
)))

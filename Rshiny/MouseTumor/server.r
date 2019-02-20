library(shiny)
library(knitr)
library(ggplot2)
library(scales)
library(Matrix)
source("./util.R")

shinyServer(function(input, output) {

  observeEvent(input$text,{
          
    output$distPlot <- renderPlot({
        input_list <- ShinnyInput(exp = exp, tsne = tsne,
                                  Dataset1 = input$Dataset1)
        exp_dat = input_list$exp
        coord_dat = input_list$tsne

      ### if there is an update of gene name then do the following below 
    
      if(input$update) {
        ## Printing the Error statement if the gene name is not present 
        validate(
          need(input$text %in% row.names(exp_dat), 'Error: the requested gene was not detected in any cell!')
        )

        gene_exp_coordinates <- PrepareExpTsne(exp_dat=exp_dat, coord_dat=coord_dat, input=input)
      ## Create a table for printing out values in Plot 1 in the below section 
      gene_exp_coordinates_sub=gene_exp_coordinates[,c("x","y","Sample","gene")] 
      colnames(gene_exp_coordinates_sub)=c("tSNE_1","tSNE_2","Sample",input$text)
      rownames(gene_exp_coordinates_sub)=gene_exp_coordinates_sub[,c("Sample")]
      gene_exp_coordinates_sub=gene_exp_coordinates_sub[,-3]

        output$info3 <- renderPrint({
        # With base graphics, need to tell it what the x and y variables are.
        nearPoints(gene_exp_coordinates_sub, input$plot_click,xvar="tSNE_1",yvar="tSNE_2")
        nearPoints(gene_exp_coordinates_sub, input$plot_hover,xvar="tSNE_1",yvar="tSNE_2")
        # nearPoints() also works with hover and dblclick events
      })
    ######################
        g <- Shinnyplot(data = gene_exp_coordinates, input = input)
        ChangeColorScale(g, alpha.use = 1,
                        scaled.expression.threshold = as.numeric(input$threshold),
                        gradient.use = c("grey", input$Color1))
      }
      ### Default view of the plot is there in else statement
      
      else {
        ## Create a table for printing out values in Plot 1 in the below section 
        output$info <- renderText({
          xy_str <- function(e) {
            if(is.null(e)) return("NULL\n")
            paste0("x=", round(e$x, 1), " y=", round(e$y, 1)," Sample=", e$z, "\n")
          }
          paste0(
            "click: ", xy_str(input$plot_click),
            "dblclick: ", xy_str(input$plot_dblclick),
            "hover: ", xy_str(input$plot_hover)
          )
        })

        output$info3 <- renderPrint({
          # With base graphics, need to tell it what the x and y variables are.
          nearPoints(coord_dat, input$plot_click,xvar="tSNE_1",yvar="tSNE_2")
          nearPoints(coord_dat, input$plot_hover,xvar="tSNE_1",yvar="tSNE_2")
          # nearPoints() also works with hover and dblclick events
        })
        gene_exp_coordinates = as.data.frame(tsne[samples[1]])
        colnames(gene_exp_coordinates) = c("x","y")
        gene_exp_coordinates$gene = 1
        g <- Shinnyplot(data = gene_exp_coordinates, input = input)
        ChangeColorScale(g, alpha.use = 1,
                         scaled.expression.threshold = as.numeric(input$threshold),
                         gradient.use = c("grey", input$Color1))
      } ### end of else 
      ##until here 
      
    }) ### end of render Plot 

    ###write the function plot for downloading the figure ###
    plotfunc <- function(){
      if(input$update) {
        input_list <- ShinnyInput(exp = exp, tsne = tsne,
                            Dataset1 = input$Dataset1)
        exp_dat = input_list$exp
        coord_dat = input_list$tsne
        
        gene_exp_coordinates <- PrepareExpTsne(exp_dat=exp_dat, coord_dat=coord_dat, input=input)
        g <- Shinnyplot(data = gene_exp_coordinates, input = input)
        ChangeColorScale(g, alpha.use = 1,
                         scaled.expression.threshold = as.numeric(input$threshold),
                         gradient.use = c("grey", input$Color1))
      }
      ### Default view of the pot is there in else statement
      else {
            input_list <- ShinnyInput(exp = exp, tsne = tsne,
                                    Dataset1 = input$Dataset1)
            coord_dat = as.data.frame(input_list$tsne)
            colnames(coord_dat) =c("x","y")
            coord_dat$gene = 1
            Shinnyplot(data = coord_dat, input = input)
      }
    }
    ### Download plot ###
    
    output$downloadplot <- downloadHandler(
      filename = function() { paste('Plot', '.png') },
      content = function(file) {
        png(file, units="in", width=10, height=7,res=600)
        print(plotfunc())
        dev.off()
      })
  })
})


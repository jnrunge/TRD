# minimally adapted from Andreas Tsouaris's original script


#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#
## The phenotype file must contain 3 columns (row, column and phenotype)
## The phenotype can be a size or ratio or whatever you want
library(shiny)
library(shinydashboard)
library(tidyverse)

# Define UI for application that draws a histogram
ui <- dashboardPage(
    dashboardHeader(title = "Basic dashboard"),
    dashboardSidebar(),
    dashboardBody(
        # Boxes need to be put in a row (or column)
        fluidRow(
            box(width=12,title='PIXL colonies based on size', minimisable=TRUE,
                fileInput("file", label = h3("Phenotype file import :")),
                div('The phenotype file must contain 3 columns (row, column and phenotype) separated by tabs'),
                # sliderInput("sel_range",
                #             "Phenotype range:",
                #             min = min(0),
                #             max = max(6000),
                #             value = c(2000,4000)),
                numericInput("numMin", label = h3("Range minimum"), value = 0),
                numericInput("numMax", label = h3("Range maximum"), value = 500),
                plotOutput("distPlot"),
                verbatimTextOutput('sel_count'),
                fluidRow(column(4,selectInput("select_format", label = h3("Target plate format"), 
                            choices = list("96 colonies" = 96, "384 colonies" = 384, "1536 colonies" = 1536), 
                            selected = 2)),
                         column(4,numericInput("first_col", label = h3("First position column"), value = 1)),
                         column(4,numericInput("first_row", label = h3("First position row"), value = 1))
                
                ),
                fluidRow(
                    column(8, actionButton('process', 'Prepare PIXL program')),
                    column(4, downloadButton('download_PIXL', label='Download'))
                )
                
                
                )
            )
        )
    )

# Define server logic required to draw a histogram
server <- function(input, output) {
   # dat2 <- reactiveVal('NULL')

    reactVals <- reactiveValues(dat2='empty')
    observeEvent(input$file, {
        #print(input$file$name)
        dat <- read.csv(input$file$datapath,sep="\t",header=F,comment.char = "#")
        dat2 <- data.frame(row=dat[,1],
                           col=dat[,2],
                           pheno=dat[,3])
        reactVals$dat2 <- dat2
        output$distPlot <- renderPlot({
            # generate bins based on input$bins from ui.R
            coords <- data.frame(x=c(input$numMin, input$numMax), y=c(0,1))
            # draw the histogram with the specified number of bins
            ggplot()+
                geom_density(data=dat2, aes(x=pheno))+
                geom_rect(data=coords, aes(xmin=input$numMin, xmax=input$numMax), ymin=0, ymax=100, fill='blue', alpha=0.3)+
                xlim(c(min(dat2$pheno), max(dat2$pheno)))+
                ylab('Density')+
                xlab('Phenotype')
        })
        selected <- dat2 %>% filter(between(pheno, input$numMin, input$numMax))
        count <- nrow(selected)
        textOut <- paste(count, 'colonies have been selected')
        output$sel_count <- renderText(textOut)
    })
    
    observeEvent(input$numMin,{
        print(input$numMin)
        if(reactVals$dat2 == 'empty'){
            print('It is empty')
        } else {
            selected <- reactVals$dat2 %>% filter(between(pheno, input$numMin, input$numMax))
            count <- nrow(selected)
            textOut <- paste(count, 'colonies have been selected')
            output$sel_count <- renderText(textOut)
        }
    }
    )
    observeEvent(input$select_format, {
        print(paste('Current format :', input$select_format))
        
    })
    
    observeEvent(input$process, {
        print('Processing')
        first_pos <- c(input$first_col, input$first_row)
        #Getting the input positions
        in_positions <- reactVals$dat2  %>% filter(between(pheno, input$numMin, input$numMax)) %>%
            rename(in_row = row, in_col=col) %>% select(-pheno)
        #Getting the out positions
        max_col <- sqrt(as.numeric(input$select_format)/(2/3))
        max_row <- max_col*(2/3)
        print(max_row)
        target_positions <- data.frame(target_row = rep(1:max_row, max_col),
                                    target_col = rep(1:max_col, each=max_row)) %>%
            filter(target_col >= first_pos[1], !(target_col == first_pos[1] & target_row < first_pos[2]))
        if(nrow(in_positions) <= nrow(target_positions)){
            target_positions <- target_positions[1:nrow(in_positions),]
            movements <- data.frame(source_plate='Source_Plate', in_row=in_positions$in_row, in_col=in_positions$in_col,
                                    target_plate='Target_Plate', target_row=target_positions$target_row, target_col=target_positions$target_col)
           header <- data.frame(source_plate=c('Source_Plate','Target_Plate'), in_row='SBS', in_col=c(nrow(reactVals$dat2), input$select_format),
                                target_plate=c('Source','Target'), target_row='', target_col='')
           movements <- rbind(header, movements)
           colnames(movements) <- header[1,]
           movements <- movements[2:nrow(movements),]
           output$download_PIXL <- downloadHandler(
               filename = 'PIXL_program.csv',
               content=function(file){
                   write.csv(movements, file, row.names=FALSE, col.names=FALSE, sep=',', quote=F)
               }
           )
        } else {
            print('Need more than one plate')
            second_plate=data.frame(target_row = rep(1:max_row, max_col),
                                    target_col = rep(1:max_col, each=max_row)) # assumed empty
            targets_first_plate=nrow(target_positions)
            target_positions=rbind(target_positions, second_plate[1:(nrow(in_positions)-targets_first_plate),])
            target_positions <- target_positions[1:nrow(in_positions),]
            print("1!")
            movements <- data.frame(source_plate='Source_Plate', in_row=in_positions$in_row, in_col=in_positions$in_col,
                                    target_plate=c(rep('Target_Plate',targets_first_plate),
                                                   rep("Target_Plate_2", (nrow(in_positions)-targets_first_plate))), target_row=target_positions$target_row, target_col=target_positions$target_col)
            print("1!")
            header <- data.frame(source_plate=c('Source_Plate','Target_Plate', "Target_Plate_2"), in_row='SBS', in_col=c(nrow(reactVals$dat2), input$select_format, input$select_format),
                                 target_plate=c('Source','Target',"Target"), target_row='', target_col='')
            print("1!")
            movements <- rbind(header, movements)
            print("1!")
            colnames(movements) <- header[1,]
            movements <- movements[2:nrow(movements),]
            output$download_PIXL <- downloadHandler(
                filename = 'PIXL_program.csv',
                content=function(file){
                    write.csv(movements, file, row.names=FALSE, col.names=FALSE, sep=',', quote=F)
                }
            )
        }
        
        
        print(first_pos)
    })
    
}

# Run the application 
shinyApp(ui = ui, server = server)

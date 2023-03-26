library(shiny)
library(ggplot2)
library(devtools)
library(naturalsort)
library(scales)
source_url('https://raw.githubusercontent.com/jnrunge/general/main/functions.R')

setwd("~/Documents/GitHub/TRD/Shiny/TRD")
crosses=list.files("../data/")
crosses=unlist(lapply(crosses,getFirst_v2,split="-"))
getTRD=function(c){
  return(mutate(fread(paste0("../data/",c,"-AF.csv.gz")),AD_A1_rel=AD_A1/sumCount))
}
df_TRD<-getTRD(crosses[1])

chrs=summarise(group_by(df_TRD, chr),maxpos=max(pos))
chrs=chrs[naturalorder(chrs$chr),]



ui <- fluidPage(
  titlePanel("title panel"),
  
  sidebarLayout(
    sidebarPanel("sidebar panel",selectInput("select", label = h3("Select box"), 
                                             choices = crosses, 
                                             selected = 1)),
    mainPanel("main panel",
              fluidRow(
                column(width = 4, class = "well",
                       h4("Brush and double-click to zoom"),
                       plotOutput("plot1", height = 300,
                                  click = "plot1_click",
                                  dblclick = "plot1_dblclick",
                                  brush = brushOpts(
                                    id = "plot1_brush",
                                    resetOnNew = TRUE
                                  )
                       )
                )
                
              ),fluidRow(
                column(width = 6,
                       h4("Points near click"),
                       verbatimTextOutput("click_info")
                ),
                column(width = 6,
                       h4("Brushed points"),
                       verbatimTextOutput("brush_info")
                )
              )),
    
    
  ),
  
  
  hr(),
  fluidRow(column(3, verbatimTextOutput("value")))
)
# Define server logic ----
server <- function(input, output) {
  output$value <- renderPrint({ summary(df_TRD<-getTRD(input$select))})
  ranges <- reactiveValues(x = NULL, y = NULL)
  
  output$plot1 <- renderPlot({
    ggplot(df_TRD<-getTRD(input$select), aes(global_pos, AD_A1_rel))+
      geom_point(alpha=0.1,color="grey")+geom_line(mapping=aes(global_pos, smoothed, color=abs(0.5-smoothed)), inherit.aes=FALSE, linewidth=2)+
      scale_color_viridis_c(option="A", limits = c(0,0.5))+
      ylim(c(0,1))+geom_hline(yintercept = 0.5)+
      geom_vline(xintercept = chrs$global_pos)+theme_bw(16)+ylab("Allele Frequency")+xlab("POS")+theme(legend.position = "none")+
      #geom_hline(yintercept = c(0.4,0.6))+
      ggtitle(input$select)+labs(alpha="Coverage")+scale_x_continuous(labels = comma)+
      coord_cartesian(xlim = ranges$x, ylim = ranges$y, expand = FALSE)
  })
  
  # When a double-click happens, check if there's a brush on the plot.
  # If so, zoom to the brush bounds; if not, reset the zoom.
  observeEvent(input$plot1_dblclick, {
    brush <- input$plot1_brush
    if (!is.null(brush)) {
      ranges$x <- c(brush$xmin, brush$xmax)
      ranges$y <- c(brush$ymin, brush$ymax)
      
    } else {
      ranges$x <- NULL
      ranges$y <- NULL
    }
  })
  
  output$click_info <- renderPrint({
    # Because it's a ggplot2, we don't need to supply xvar or yvar; if this
    # were a base graphics plot, we'd need those.
    nearPoints(getTRD(input$select), input$plot1_click, addDist = TRUE)
  })
  
  output$brush_info <- renderPrint({
    brushedPoints(getTRD(input$select), input$plot1_brush)
  })
}

shinyApp(ui, server)
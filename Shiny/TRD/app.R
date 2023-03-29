library(shiny)
library(ggplot2)
library(devtools)
library(naturalsort)
library(scales)
library(ggtree)
library(ape)
source_url('https://raw.githubusercontent.com/jnrunge/general/main/functions.R')

tree=ape::read.tree("../data/Victor/full2543Matrix.DP10.GQ20.SNPs.99pNonMiss.Biallelic.IDset.vcf.gz.swapped.newick")

setwd("~/Documents/GitHub/TRD/Shiny/TRD")
crosses=list.files("../data/")
crosses=unlist(lapply(crosses,getFirst_v2,split="-"))
crosses=unlist(lapply(crosses,getFirst_v2,split="."))
crosses=unique(crosses)
getTRD=function(c){
  return(mutate(fread(paste0("../data/",c,"-AF.csv.gz")),AD_A1_rel=AD_A1/sumCount))
}
getCov=function(c){
  return(fread(paste0("../data/",c,".bam.depth.gz.summary.gz")))
}
getTRD_regions=function(c){
  if(!file.exists(paste0("../data/",c,"-TRD_regions.csv.gz"))){
    return("none")
  }
  return(fread(paste0("../data/",c,"-TRD_regions.csv.gz")))
}
df_TRD<-getTRD(crosses[1])
df_cov<-getCov(crosses[1])

TRD_regions=getTRD_regions(crosses[1])

chrs=summarise(group_by(df_TRD, chr),maxpos=max(pos))
chrs=chrs[naturalorder(chrs$chr),]

getSelectionForTRDRegions=function(TRD_regions){
  if(!is.data.table(TRD_regions)){
    return("none")
  }
  return(c("none",TRD_regions$ID))
}

ui <- fluidPage(
  titlePanel("title panel"),
  
  sidebarLayout(
    sidebarPanel("sidebar panel",selectInput("select", label = h3("Select box"), 
                                             choices = crosses, 
                                             selected = 1),
                 selectInput("selectTRDregion", label = h3("Select box"), 
                             choices = getSelectionForTRDRegions(TRD_regions), 
                             selected = 1)),
    mainPanel("main panel",
              fluidRow(
                column(width = 8, class = "well",
                       h4("Brush and double-click to zoom"),
                       plotOutput("plot1", height = 300,
                                  click = "plot1_click",
                                  dblclick = "plot1_dblclick",
                                  brush = brushOpts(
                                    id = "plot1_brush",
                                    resetOnNew = TRUE
                                  )
                       ),
                       plotOutput("plot_cov", height = 300,
                                  )
                
                ),
                column(width = 4, class = "well",
                       h4("Phylogeny"),
                       plotOutput("plot_phylo", height = 300,
                       ))
                
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
  rv <- reactiveValues(choices=NULL)
  observeEvent(input$select, {
      
      TRD_regions=getTRD_regions(input$select)
      rv$choices=getSelectionForTRDRegions(TRD_regions)
      updateSelectInput(inputId = "selectTRDregion", choices = rv$choices)
    
  })
  
  output$value <- renderPrint({ summary(df_TRD<-getTRD(input$select))})
  ranges <- reactiveValues(x = NULL, y = NULL)
  
  output$plot1 <- renderPlot({
    
    p<-ggplot(df_TRD<-getTRD(input$select), aes(global_pos, AD_A1_rel))+
      geom_point(alpha=0.1,color="grey")+geom_line(mapping=aes(global_pos, smoothed, color=abs(0.5-smoothed)), inherit.aes=FALSE, linewidth=2)+
      scale_color_viridis_c(option="A", limits = c(0,0.5))+
      ylim(c(0,1))+geom_hline(yintercept = 0.5)+
      geom_vline(xintercept = chrs$global_pos)+theme_bw(16)+ylab("Allele Frequency")+xlab("POS")+theme(legend.position = "none")+
      #geom_hline(yintercept = c(0.4,0.6))+
      ggtitle(input$select)+labs(alpha="Coverage")+scale_x_continuous(labels = comma)+
      coord_cartesian(xlim = ranges$x, ylim = ranges$y, expand = FALSE)
    
    if(is.data.table(TRD_regions<-getTRD_regions(input$select))){
      p=p+geom_label(data = TRD_regions, mapping=aes(y=0.5,x=global_start+(global_end-global_start)/2,label=ID),color="red",size=6,fill="white",
                     inherit.aes = FALSE)
    }
      
    
    if(input$selectTRDregion=="none"){
      ranges$x=NULL
      print(p)
    }else{
      ranges$x=c(TRD_regions$global_start[TRD_regions$ID==input$selectTRDregion],TRD_regions$global_end[TRD_regions$ID==input$selectTRDregion])
      print(p)
    }
  })
  
  output$plot_cov <- renderPlot({
    
    p<-ggplot(df_cov<-getCov(input$select), aes(global_pos,xend=global_pos+999, avg_cov))+
      geom_line(alpha=1,color="darkgrey")+scale_y_log10()+geom_hline(yintercept = mean(df_cov$avg_cov))+
      geom_vline(xintercept = chrs$global_pos)+theme_bw(16)+ylab("Allele Frequency")+xlab("POS")+theme(legend.position = "none")+
      ggtitle("Average coverage")+scale_x_continuous(labels = comma)+
      coord_cartesian(xlim = ranges$x, ylim = ranges$y, expand = FALSE)
    TRD_regions<-getTRD_regions(input$select)
    
    if(input$selectTRDregion=="none"){
      ranges$x=NULL
      print(p)
    }else{
      ranges$x=c(TRD_regions$global_start[TRD_regions$ID==input$selectTRDregion],TRD_regions$global_end[TRD_regions$ID==input$selectTRDregion])
      print(p)
    }
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
  
  output$plot_phylo <- renderPlot({
    # group muzst be sorted by tree$tip.label
    node_df=data.frame(node=1:length(tree$tip.label),
                   group=c("a1",rep("none",length(tree$tip.label)-2),"a2"))
    node_df$group=factor(node_df$group,levels=c("a1","a2","none"))
    tree_groups=dplyr::left_join(tree,node_df, by="node")
    p<-ggtree(tree_groups, 
           branch.length='none', 
           layout='circular',size=0.5)+
      geom_tippoint(aes(color=group,size=group,alpha=group),shape=3)+
      scale_size_manual(values=c(10,10,0))+scale_alpha_manual(values=c(1,1,0))
    
    print(p)
    
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
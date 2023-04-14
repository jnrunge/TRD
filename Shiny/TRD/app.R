library(shiny)
library(ggplot2)
library(devtools)
library(naturalsort)
library(scales)
library(ggtree)
library(ape)
library(DT)
source_url('https://raw.githubusercontent.com/jnrunge/general/main/functions.R')
setwd("~/Documents/GitHub/TRD/Shiny/TRD")
start=TRUE
tree=ape::read.tree("../data/Victor/full2543Matrix.DP10.GQ20.SNPs.99pNonMiss.Biallelic.IDset.vcf.gz.swapped.newick")
df_Strains=fread("../data/Victor/operationalTable_Full2543Sace_Clades.csv")
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
getAlleleSharing=function(c){
  return(fread(paste0("../data/",c,"-AF.csv.gz.allelesharing.csv.gz")))
}
getTRD_regions=function(c){
  if(!file.exists(paste0("../data/",c,"-TRD_regions.csv.gz"))){
    return("none")
  }
  return(fread(paste0("../data/",c,"-TRD_regions.csv.gz")))
}
GFF=ape::read.gff("../data/sace_R64-3-1_20210421.gff.gz")
GFF$attributes=str_replace_all(GFF$attributes, fixed(";"), "; ")
GFF$attributes=str_replace_all(GFF$attributes, fixed("%20"), " ")
GFF$attributes=str_replace_all(GFF$attributes, fixed(","), ", ")
getGFFsubset=function(x, brush_df){
  brush_df=brush_df[x,]

  
  return(subset(GFF, seqid==brush_df$chr & (start<=brush_df$pos & end>=brush_df$pos)))
}

prep_df_AS_filtered=function(brush_df,selectSimilarity){
  df_AS_filtered=df_AS
  
  #print(head(df_AS_filtered))
  
  #print(head(brush_df))
  if(nrow(brush_df)>0){
    df_AS_filtered=subset(df_AS_filtered, paste(`#CHROM`,POS) %in% paste(brush_df$chr,brush_df$pos))
  }
  melted=reshape2::melt(df_AS_filtered, id.vars = c("#CHROM","POS"))
  melted=filter(melted, variable != "chrpos")
  tmp=summarise(group_by(melted, variable), nAll=n())
  vcf_translated_summary=left_join(tmp,summarise(group_by(melted, variable, value), n=n()), by=c("variable"))%>%mutate(p=n/nAll)%>%select(variable,value,p)%>%rename(Strain=variable, Type=value)
  #vcf_translated_summary=vcf_translated_summary[match(tree$tip.label, vcf_translated_summary$Strain),]
  
  #print(head(vcf_translated_summary%>%arrange(-p)))
  
  A1s=vcf_translated_summary$Strain[vcf_translated_summary$Type=="A1_hom" & vcf_translated_summary$p>=selectSimilarity]
  A2s=vcf_translated_summary$Strain[vcf_translated_summary$Type=="A2_hom" & vcf_translated_summary$p>=selectSimilarity]
  return_list=list(df_AS_filtered=df_AS_filtered,A1s=A1s,A2s=A2s)
  
  return_list$A1s=as.character(return_list$A1s)
  return_list$A2s=as.character(return_list$A2s)
  return(return_list)
}

df_TRD=getTRD(crosses[1])

df_cov=getCov(crosses[1])
df_AS=getAlleleSharing(crosses[1])

df_AS_filtered=prep_df_AS_filtered(data.frame(), 0.7)

data_all=list(df_TRD=df_TRD,df_cov=df_cov,df_AS=df_AS,df_AS_filtered=df_AS_filtered)

TRD_regions=getTRD_regions(crosses[1])

chrs=summarise(group_by(df_TRD, chr),maxpos=max(global_pos))
chrs=chrs[naturalorder(chrs$chr),]

getSelectionForTRDRegions=function(TRD_regions){
  if(!is.data.table(TRD_regions)){
    return("none")
  }
  return(c("none",TRD_regions$ID))
}

ui <- fluidPage(
  titlePanel("TRD results exploration"),
  
  sidebarLayout(
    sidebarPanel(width=2,"Selections",selectInput("select", label = h3("Cross"), 
                                             choices = crosses, 
                                             selected = "ChrisC1"),
                 selectInput("selectTRDregion", label = h3("Zoom into"), 
                             choices = getSelectionForTRDRegions(TRD_regions), 
                             selected = "none"),
                 selectInput("selectSimilarity", label = h3("Similarity (phylo"), 
                             choices = c(0.5,0.6,0.7,0.8,0.9), 
                             selected = 0.7)),
    mainPanel(NULL,
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
                       )),
                column(width = 4,
                       h4("Populations with similar alleles"),
                       verbatimTextOutput("pop_info")
                )
                
              ),fluidRow(
                column(width=8, h4("Annotations in brushed region"),
                        DT::dataTableOutput("GFF_subset")),
                
              ),
              ),
    
    
  ),
  
  
)
# Define server logic ----
server <- function(input, output) {
  rv <- reactiveValues(choices=NULL)
  observeEvent(input$select, {
      
      TRD_regions=getTRD_regions(input$select)
      rv$choices=getSelectionForTRDRegions(TRD_regions)
      updateSelectInput(inputId = "selectTRDregion", choices = rv$choices)
      if(start==FALSE){
        data_all=reactiveValues(df_TRD=getTRD(input$select),
                                df_cov=getCov(input$select),
                                df_AS=getAlleleSharing(input$select),
                                df_AS_filtered=prep_df_AS_filtered(brushedPoints(df_TRD, input$plot1_brush), input$selectSimilarity))
      }else{
        start=FALSE
      }
        
      
      
      
      
    
  })
  observeEvent(input$plot1_brush,{
    data_all$df_AS_filtered = prep_df_AS_filtered(brushedPoints(df_TRD, input$plot1_brush),
                                                  input$selectSimilarity)
  })
  if(exists("data_all")){
    df_TRD=data_all$df_TRD
  }
  output$value <- renderPrint({ summary(df_TRD)})
  ranges <- reactiveValues(x = NULL, y = NULL)
  
  output$plot1 <- renderPlot({
    if(exists("data_all")){
      df_TRD=data_all$df_TRD
    }
    
    p<-ggplot(df_TRD, aes(global_pos, AD_A1_rel))+
      geom_point(alpha=0.1,color="grey")+geom_line(mapping=aes(global_pos, smoothed, color=abs(0.5-smoothed)), inherit.aes=FALSE, linewidth=2)+
      scale_color_viridis_c(option="A", limits = c(0,0.5))+
      ylim(c(0,1))+geom_hline(yintercept = 0.5)+
      geom_vline(xintercept = chrs$maxpos)+theme_bw(16)+ylab("Allele Frequency")+xlab("POS")+theme(legend.position = "none")+
      #geom_hline(yintercept = c(0.4,0.6))+
      ggtitle(input$select)+labs(alpha="Coverage")+scale_x_continuous(labels = comma)+
      coord_cartesian(xlim = ranges$x, ylim = ranges$y, expand = FALSE)
    
    if(is.data.table(TRD_regions<-getTRD_regions(input$select))){
      p=p+geom_label(data = TRD_regions, mapping=aes(y=0.5,x=global_start+(global_end-global_start)/2,label=ID),color="red",size=6,fill="white",
                     inherit.aes = FALSE)
    }
      
    
    if(input$selectTRDregion=="none"){
      ranges$x=NULL
      p
    }else{
      ranges$x=c(TRD_regions$global_start[TRD_regions$ID==input$selectTRDregion],TRD_regions$global_end[TRD_regions$ID==input$selectTRDregion])
      p
    }
  })
  
  output$plot_cov <- renderPlot({
    if(exists("data_all")){
      df_cov=data_all$df_cov
    }
    p<-ggplot(df_cov, aes(global_pos,xend=global_pos+999, avg_cov))+
      geom_line(alpha=1,color="darkgrey")+scale_y_log10()+geom_hline(yintercept = mean(df_cov$avg_cov))+
      geom_vline(xintercept = chrs$global_pos)+theme_bw(16)+ylab("Allele Frequency")+xlab("POS")+theme(legend.position = "none")+
      ggtitle("Average coverage")+scale_x_continuous(labels = comma)+
      coord_cartesian(xlim = ranges$x, ylim = ranges$y, expand = FALSE)
    TRD_regions<-getTRD_regions(input$select)
    
    if(input$selectTRDregion=="none"){
      ranges$x=NULL
      p
    }else{
      ranges$x=c(TRD_regions$global_start[TRD_regions$ID==input$selectTRDregion],TRD_regions$global_end[TRD_regions$ID==input$selectTRDregion])
      p
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
    
    data_all$df_AS_filtered = prep_df_AS_filtered(brushedPoints(df_TRD, input$plot1_brush),
                                                  input$selectSimilarity)
    
    df_AS_filtered=data_all$df_AS_filtered$df_AS_filtered
    A1s=data_all$df_AS_filtered$A1s
    A2s=data_all$df_AS_filtered$A2s
    
    #print(head(A1s))
    
    all=rep("none",length(tree$tip.label))
    all[tree$tip.label%in%A1s]="a1"
    all[tree$tip.label%in%A2s]="a2"
    
    node_df=data.frame(node=1:length(tree$tip.label),
                   group=all)
    node_df$group=factor(node_df$group,levels=c("a1","a2","none"))
    tree_groups=dplyr::left_join(tree,node_df, by="node")
    p<-ggtree(tree_groups, 
           branch.length='none', 
           layout='circular',size=0.5)+
      geom_tippoint(aes(color=group,size=group,alpha=group),shape=3)+
      scale_size_manual(values=c(10,10,0))+scale_alpha_manual(values=c(1,1,0))
    
    p
    
  })
  
  output$pop_info<-renderPrint({
    
    data_all$df_AS_filtered = prep_df_AS_filtered(brushedPoints(df_TRD, input$plot1_brush),
                                                  input$selectSimilarity)
    
    df_AS_filtered=data_all$df_AS_filtered$df_AS_filtered
    A1s=data_all$df_AS_filtered$A1s
    A2s=data_all$df_AS_filtered$A2s
    

    
    bind_rows(summarise(group_by(filter(df_Strains, StandardizedName %in% A1s),
                       Clade), n=n()) %>% arrange(-n)%>%mutate(Type="A1_hom"),
              summarise(group_by(filter(df_Strains, StandardizedName %in% A2s),
                                 Clade), n=n()) %>% arrange(-n)%>%mutate(Type="A2_hom"))%>% arrange(-n)
    
  })
  
  output$GFF_subset = DT::renderDataTable({
    if(nrow(brushedPoints(df_TRD, input$plot1_brush))==0){
      GFF
    }else{
      brush_df<-brushedPoints(df_TRD, input$plot1_brush)
      unique(bind_rows(lapply(1:nrow(brush_df),FUN=getGFFsubset,brush_df=brush_df)))
    }
    
  })
  

  

}

shinyApp(ui, server)
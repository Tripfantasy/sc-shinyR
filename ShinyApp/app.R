source("module.R")
options(shiny.maxRequestSize=1000000*1024^2)
ui <- fluidPage(
  theme = shinytheme("darkly"),
  
  tabsetPanel(
    #First tab will act as config for analysis. Setup input, initial assay, etc. 
    tabPanel("Home", "Files and Configuration",
             fileInput("object","Input Object",accept = c(".rds",".Rds")),
             
             selectInput("preprocess_check",
                         label = "Preprocessed?",
                         choices = c(" ","Yes","No"),
                         selected = NULL),
             selectInput("assay",
                         label = "Assay Type",
                         choices = c(" ", "Spatial", "RNA","SCT", "Integrated"),
                         selected = NULL),
             
             actionButton("load_obj",
                          "Load Object")
    ),
    tabPanel("Preprocessing", "Normalization and Filtering",
             textInput("mt",
                       label = "Mitochondrial Filter",
                       value = "5",
                       width = "33%",
                       placeholder = " % MT "
                       ),
             
             textInput("nfL",
                       label = "Minimum Features Per Cell",
                       value = "200",
                       width = "33%",
                       placeholder = "Nfeature per Cell"
                       ),

             textInput("nfU",
                       label = "Maximum Features Per Cell",
                       value = "2500",
                       width = "33%",
                       placeholder = "Nfeature per Cell"),
             
             textInput("sct_features",
                       label = "Desired Variable Features",
                       value = "2000",
                       width = "33%",
                       placeholder = "Variable Features"),
             
             textInput("pcs",
                       label = "Principal Components",
                       value = "10",
                       width = "33%",
                       placeholder = "PCs"),
             
              selectInput("qcplot_check",
                          label = "Add QC plots to report?",
                          choices = c("","Yes","No")),
             
             actionButton("runPreprocess",
                          label = "Start Preprocessing"),
             
             actionButton("qcplots", "Plot Distributions"),
             
             plotOutput("qcplot1"),
             
             tags$div(id = "message", style = "color: red")
             ),

    
    tabPanel("Integration", "Integrate Objects",
            fileInput("custom_input",
                      label = "Files to Integrate",
                      multiple = T),
            actionButton("integrate_data",
                         label = "Integrate Objects"),
            actionButton("plot_integration",
                         label = "Plot Integrated Object"),
            plotOutput("integration_plot")
    ),
    
    tabPanel("Annotation", "Cell Type and Visualization",
             selectInput("sample_type",
                         label = "Sample",
                         choices = c("","Brain","PBMC Test")),
             selectInput("Species",
                         label = "Species",
                         choices = c("","Human","Mouse")),
             actionButton("label_clust",
                          label = "Plot Neighboring Reference"),
             actionButton("transfer_annotations",
                          label = "Annotate Cells"),
             plotOutput("clusplot")
            ),
    
    tabPanel("Differential Expression", "Marker Gene Analysis")
  ))



server <- function(input, output, session) {
  
  ### Load object
  observeEvent(input$load_obj, {
    # Read the object from the fileInput
    seurat <- load_seurat_obj(input$object$datapath)
    
    
    if (is.vector(seurat)){
      showModal(modalDialog(
        title = "Error with file.",
        HTML("<h5>There is an error with the file you uploaded. See below for more details.</h5><br>",
             paste(unlist(seurat), collapse = "<br><br>"))
      ))
    } else {
      #message(class(seurat)) <- For troubleshooting! - DM
      showModal(modalDialog(
        title = "Success!",
        "You've loaded a compatible input object. Please Continue."
      ))
    }
  })




  ### Preprocessing functionality
  
  #Check Fields
  observeEvent(input$runPreprocess, {
        blankfields <- preprocess_fields(mt = input$mt, nfU = input$nfU, nfL = input$nfL, assay = eval(input$assay), 
                                    sct_features = input$sct_features,pcs = input$pcs)
    if(length(blankfields > 0)){
      # Displays message telling user what is missing. (May eventually want error for NaN/Inf cases)
      showModal(modalDialog(
        title = "Missing Input Error",
        print(paste0(blankfields, collapse = " , "))
      ))} else{
        # Populates input variables for preprocessing.
         seurat <- load_seurat_obj(input$object$datapath)
        # message(class(seurat)) <- For troubleshooting! - DM
         mt <<- input$mt
         nfU <<- input$nfU
         nfL <<- input$nfL
         assay <<- input$assay 
         sct_features <<- input$sct_features
         pcs <<- input$pcs 
         
         
         # Preprocess the object (Filter, normalize, cluster, check module for details)
         obj <<- preprocess(obj = seurat,mt, nfU,nfL,assay,sct_features,pcs)
         observeEvent(input$qcplots, {
           
           # Displays distribution of cells given %MT , nFeature, nCount.
           output$qcplot1 <- renderPlot({
             # Note, we group.by orig.ident for 1 plot per metric. Modify this later to allow user to choose what to group.by.  
             VlnPlot(object = obj,features = c("nFeature_RNA","nCount_RNA","percent.mt"), group.by = "orig.ident", cols = "cadetblue") 
           })
         })
      }
  })
### Integration functionality
  observeEvent(input$custom_input,{
    reference <<- load_seurat_obj(input$custom_input$datapath)
  })
  
  observeEvent(input$integrate_data,{
    object.integrated <<- integrate_reference(source = obj,reference = reference)
  })
  
  
  
  observeEvent(input$plot_integration,{
    output$integration_plot <- renderPlot({
      DefaultAssay(object = object.integrated) <- "integrated"
      p1 <- DimPlot(object.integrated,reduction = "umap",group.by = "split_vector", label = T, 
              label.size = 3, repel = T) + ggtitle("Integrated Sources")
      p2 <- DimPlot(object.integrated,reduction = "umap", group.by =  "seurat_clusters",
                    label = T,label.size = 3,repel = T) + ggtitle("Integrated Clusters")
      p1 + p2
    })
  })
  
### Cell Annotation functionality 

  
  observeEvent(input$label_clust,{
    output$clusplot <- renderPlot({
      DimPlot(object.integrated, reduction = "umap", group.by = "seurat_annotations", label = T,
              label.size = 6,repel = T) + ggtitle("Potential Annotations")
    })
  })
  
  observeEvent(input$transfer_annotations,{
    democratic_annotation(object = object.integrated)
  })
 
}
  
  

shinyApp(ui, server)


















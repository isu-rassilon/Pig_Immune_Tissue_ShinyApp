# Loading libraries, make sure to set path and set the virtual environment on R
#setwd("/work/ABG/MUSKAN/ShinyApp/")
#.libPaths("/work/ABG/MUSKAN/ShinyTissues/micromamba/envs/seurat5.3.0/lib/R/library")

#https://muskan-16.shinyapps.io/shinyapp/


#also opton to run via R but needs downloadable files: runGitHub("Pig_Immune_Tissue_ShinyApp", "kapoormuskan")

#https://muskan-16.shinyapps.io/shinyapp/

library(shiny)
library(shinythemes)
library(shinyalert)
library(shinycssloaders)
library(Seurat)
library(ggplot2)
library(patchwork)
library(shinyFeedback)
library(dplyr)
library(scales)
library(presto)
library(DT)
library(memoise)
# Load .rds datasets for 4 tissues

#Load .rds datasets for 4 tissues, for now testing added 2
memo_readRDS <- memoise(function(file_path) {
  readRDS(file_path)
})
# Dataset file mapping (label -> file path)
dataset_choices <- c(
  "Bone Marrow" = "Bone_Marrow_shiny.rds",
  "Spleen" = "Spleen_shiny.rds",
  "Lymph Node" = "Lymph_Node_shiny.rds",
  "Thymus" = "Thymus_shiny.rds"
)

# Reference mapping objects
query_choices <- c(
  "Spleen" = "Spleen_mapped_shiny.rds",
  "Thymus" = "Thymus_mapped_shiny.rds",
  "Lymph Node" = "Lymph_Node_mapped_shiny.rds"
)

#data_list <- list("Bone Marrow" = bone_marrow, "Spleen" = spleen, "Lymph Node" = lymph_node, "Thymus" = thymus)
#gene_list <- unique(c(rownames(bone_marrow), rownames(spleen), rownames(lymph_node), rownames(thymus)))


ui <- fluidPage(
  useShinyFeedback(),
  theme = shinytheme("cerulean"),
  titlePanel(h1("Shiny App to Explore Single Cell RNAseq Tissues")),
  navbarPage(" ",
             tabPanel(icon("home"),
                      fluidRow(
                        column(tags$img(src = "Pig.png", width = "250px", height = "200px"), width = 2),
                        column(
                          br(),
                          p("An open-source application interface available to perform gene expression query, differential gene expression and cell type population comparisons across pig primary and secondary lymphoid organs.", 
                            style = "text-align:justify;color:black;background-color:lavender;padding:15px;border-radius:10px"),
                          br(),
                          
                          tags$ul(
                            p("The app functionalities include:"),
                            tags$li("Home Page: Includes ready-to-use and preprocessed .cloupe files, and cell-type-specific differential gene expression analysis."),
                            tags$li("Gene Expression Page: Visualize expression of selected or typed-in genes across tissues using cell annotations, UMAP, and violin plots."),
                            tags$li("Reference Mapping Page: Compare immune cell types from the reference (Bone Marrow) to predicted cell types in other tissues."),
                            style = "text-align:justify;color:black;background-color:papayawhip;padding:15px;border-radius:10px"),
                          width = 8
                        ),
                        column(
                          br(),
                          tags$img(src = "FAANG_logo_RGBc.png", width = "200px", height = "130px"),
                          br(),
                          p("For more information, visit the",
                            a(href = "https://www.faang.org", "FAANG official website", target = "_blank"),
                            style = "text-align:center;color:black"),
                          width = 2
                        ),
                       hr(),
                        fluidRow(
                          column(12,
                                 h3("Download .cloupe Files")
                          ) ),
                        fluidRow(
                          column(3, downloadButton("download_bm", "Bone Marrow.cloupe")),
                          column(3, downloadButton("download_spleen", "Spleen.cloupe")),
                          column(3, downloadButton("download_thymus", "Thymus.cloupe")),
                          column(3, downloadButton("download_ln", "Lymph Node.cloupe"))
                        ) ),
                      
                      hr(),
                      h3("Differential Gene Expression"),
                      
                      fluidRow(
                        column(2,
                               selectInput("deg_dataset1", "Dataset 1:",
                                           choices = names(dataset_choices), selected = "Bone Marrow")
                        ),
                        column(2,
                               selectInput("deg_dataset2", "Dataset 2:",
                                           choices = names(dataset_choices), selected = "Spleen")
                        ),

                        fluidRow(
                          column(2, numericInput("padj_thresh", "FDR cutoff", value = 0.05, step = 0.01)),
                          column(2, numericInput("logfc_thresh", "logFC cutoff(Dataset1/Dataset2)", value = 0.25, step = 0.01)),
                          
                        ),

                        column(2, uiOutput("deg_group1")),
                        column(2, uiOutput("deg_group2")),
                        column(2, actionButton("run_deg", "Run DEG with Presto", class = "btn-success")),
                        column(2, downloadButton("download_deg", "Download DEG Results"))
                      ),
                      
                      fluidRow(
                        column(12,
                               DT::dataTableOutput("deg_results")
                        )
                      ),
                      useShinyalert(),
                      hr(),
                      p(em("A single-cell immune atlas of primary and secondary lymphoid organs in pigs"), br("--"), style = "text-align:center; font-family: times")
             ),
             
             tabPanel("Gene Expression", icon = icon("signal"),
                      h3(p(em("Visualization of Genes"), icon("signal"), style = "color:blue;text-align:center")),
                      br(),
                      # Select multiple datasets using the loop below
                      column(width = 4,
                             selectInput("datasets", "Select Dataset:", 
                                         choices = names(dataset_choices),
                                         selected = "Bone Marrow",
                                         multiple = TRUE)),
                      # Dropdown selection for genes
                      column(width = 4,
                             selectInput("gene", label = "Select Gene:", 
                                         choices = NULL,
                                         selected = "CD3E"),
                             uiOutput("geneWarning")  # feedback warning inserted here
                      ),
                      column(width = 4,
                             actionButton("update", "Update for dataset", icon = icon("refresh"), class = "btn-primary")),
                      br(), br(),
                      column(width = 12,
                             h4("Gene Expression and Cell Type Annotations Plots"),
                             withSpinner(plotOutput("staticDimPlot", height = 700, width = "1200px")) ,
                             br(),
                             h4("Violin Plot by Cell Type"),

                             withSpinner(plotOutput("violinPlot", height = 1600)),# Changed to plotOutput()

                             withSpinner(plotOutput("violinPlot", height = 500))# Changed to plotOutput()

                      )
                      
             ),
             
             tabPanel("Reference Cell Type Mapping", icon = icon("project-diagram"),
                      fluidPage(
                        fluidRow(
                          column(width = 4,
                                 h5("Reference Dataset: Bone Marrow")
                          ),
                          column(width = 4,
                                 selectInput("query_dataset", "Select Query Dataset:",
                                             choices = names(query_choices), selected = "Thymus")
                          )
                        ),
                        hr(),
                        br(),
                        
                        # Section 1: Mapping Score + Predicted Cell Type IDs
                        fluidRow(
                          column(width = 1),  # left margin
                          column(width = 5,
                                 h4("Mapping Scores"),
                                 withSpinner(plotOutput("mapping_plot", height = 400))
                          ),
                          column(width = 6,
                                 h4("Predicted Cell Type IDs based on Reference Mapping"),
                                 withSpinner(plotOutput("predicted_umap", height = 400))
                          )
                          # right margin
                        ),
                        
                        br(), br(),
                        
                        # Section 2: Cell type dropdown centered
                        fluidRow(
                          column(width = 4, offset = 4,
                                 uiOutput("select_celltype")
                                 
                          ),
                          
                          br(), br(),
                          
                          # Section 3: Highlighted cells + Dot plot
                          fluidRow(
                            column(width = 1),
                            column(width = 5,
                                   h4("Reference and Query Prediction scores"),
                                   withSpinner(plotOutput("highlight_all", height = 400))
                            ),
                            br(), br(), br(),
                            column(width = 6,
                                   h4("Predicted v/s manually Annotated Cell Types"),
                                   withSpinner(plotOutput("dot_plot", height = 450))
                            ),
                            
                          )
                        )
                      )
             )))



server <- function(input, output, session) {
  # ---- Gene Expression Logic ----
  
  
  observe({
    shinyalert(
      title = "Welcome",
      text = "This app allows users to interactively explore scRNAseq data from four healthy pig immune tissues- Bone Marrow, Spleen, Thymus, Lymph Node.",
      type = "info",
      showConfirmButton = TRUE,
      confirmButtonText = "OK"
    )
  })
  #download .cloupe files
  output$download_bm <- downloadHandler(
    filename = function() { "Bone_Marrow.cloupe" },
    content = function(file) {
      file.copy("www/cloupe_files/bm.cloupe", file)
    } )
  output$download_spleen <- downloadHandler(
    filename = function() { "Spleen.cloupe" },
    content = function(file) {
      file.copy("www/cloupe_files/sp.cloupe", file)
    } )
  output$download_thymus <- downloadHandler(
    filename = function() { "Thymus.cloupe" },
    content = function(file) {
      file.copy("www/cloupe_files/th.cloupe", file)
    } )
  output$download_ln <- downloadHandler(
    filename = function() { "Lymph_Node.cloupe" },
    content = function(file) {
      file.copy("www/cloupe_files/ln.cloupe", file)
    })
  
  #reactive function for dataset easy loading for azure team
  datasetInput <- eventReactive(input$update, {
    req(input$datasets)
    selected <- input$datasets
    selected_files <- dataset_choices[selected]
    datasets <- lapply(selected_files, function(file) {
      memo_readRDS(file.path("data", file))
    })
    names(datasets) <- selected
    datasets
  })
  #gene page easy loading
  observeEvent(input$datasets, {
    req(input$datasets)
    files <- dataset_choices[input$datasets]
    genes <- unique(unlist(lapply(files, function(file) {
      obj <- memo_readRDS(file.path("data", file))
      rownames(GetAssayData(obj, layer="data"))
    })))
    updateSelectInput(session, "gene", choices = sort(genes), selected = "CD3E")
  })
  
  output$geneWarning <- renderUI({ NULL })
  
  #DEG##
  
  # Reactive: load selected dataset
 

  # Automatically detect a valid column
  degObj1 <- reactive({
    req(input$deg_dataset1)
    file <- dataset_choices[[input$deg_dataset1]]
    obj <- memo_readRDS(file.path("data", file))
    obj$group_label <- obj$CellType  # assumes column CellType exists
    Idents(obj) <- "group_label"
    obj})
  
  degObj2 <- reactive({
    req(input$deg_dataset2)
    file <- dataset_choices[[input$deg_dataset2]]
    obj <- memo_readRDS(file.path("data", file))
    obj$group_label <- obj$CellType
    Idents(obj) <- "group_label"
    obj})


# Populate group selections
observeEvent(input$deg_dataset1, {
  obj1 <- degObj1()
  output$deg_group1 <- renderUI({
    selectInput("deg_group1_val", "Cell Type in Dataset 1:", choices = unique(Idents(obj1)))
  })})

observeEvent(input$deg_dataset2, {
  obj2 <- degObj2()
  output$deg_group2 <- renderUI({
    selectInput("deg_group2_val", "Cell Type in Dataset 2:", choices = unique(Idents(obj2)))
  })})

degResults <- eventReactive(input$run_deg, {
  obj1 <- degObj1()
  obj2 <- degObj2()
  req(input$deg_group1_val, input$deg_group2_val)
  cells1 <- WhichCells(obj1, idents = input$deg_group1_val)
  cells2 <- WhichCells(obj2, idents = input$deg_group2_val)
  expr1 <- as.matrix(GetAssayData(obj1, layer = "data")[, cells1])
  expr2 <- as.matrix(GetAssayData(obj2, layer = "data")[, cells2])
  expr_combined <- cbind(expr1, expr2)
  meta_combined <- data.frame(group = rep(c("group1", "group2"), c(ncol(expr1), ncol(expr2))))
  raw_results <- presto::wilcoxauc(X = expr_combined, y = meta_combined$group)
# Split into two group-specific results
  df1 <- raw_results[raw_results$group == "group1", ]
  df2 <- raw_results[raw_results$group == "group2", ] 
# Merge by gene (feature)
  merged <- merge(df1, df2, by = "feature", suffixes = c("_g1", "_g2"))
  
  final_df <- data.frame(
    Gene = merged$feature,
    AvgExpr_Dataset1 = signif(merged$avgExpr_g1, 3),
    AvgExpr_Dataset2 = signif(merged$avgExpr_g2, 3),
    logFC = signif(merged$logFC_g1, 3),
    pval = signif(merged$pval_g1, 3),
    padj = signif(merged$padj_g1, 3)
  )
# Filter using UI thresholds
  filtered_df <- final_df[
    final_df$padj < input$padj_thresh &
      abs(final_df$logFC) > input$logfc_thresh,
  ]
  return(filtered_df)
})
# output render table
output$deg_results <- DT::renderDataTable({
  req(degResults())
  DT::datatable(degResults(), options = list(pageLength = 10, scrollX = TRUE), rownames = FALSE)
})

# Download handler
output$download_deg <- downloadHandler(
  filename = function() {
    paste0("DEG_", input$deg_dataset1, "_", input$deg_group1_val, "_vs_",
           input$deg_dataset2, "_", input$deg_group2_val, ".csv") },
  content = function(file) {
    write.csv(degResults(), file, row.names = FALSE) })
#####

#gene expresson page##
output$staticDimPlot <- renderPlot({
  req(input$update)
  selected_datasets <- datasetInput()
  plots <- list()
  warnings <- c()
  
  for (dataset_name in names(selected_datasets)) {
    dataset <- selected_datasets[[dataset_name]]
    if (!(input$gene %in% rownames(dataset))) {
      warnings <- c(warnings, paste("Gene", input$gene, "not found in dataset:", dataset_name))
      next }
    gene_values <- FetchData(dataset, vars = input$gene, layer="data")
    if (all(gene_values == 0, na.rm = TRUE)) {
      warnings <- c(warnings, paste("Gene", input$gene, "has zero expression in dataset:", dataset_name))
      next }
    celltype_col <- if ("celltype" %in% colnames(dataset@meta.data)) {
      "celltype" } else if ("Cell_Type" %in% colnames(dataset@meta.data)) {
      "Cell_Type"} else if ("CellType" %in% colnames(dataset@meta.data)) {
      "CellType" } else {NA }
    if (is.na(celltype_col)) {
      warnings <- c(warnings, paste("No valid cell type found in dataset:", dataset_name))
      next }
    
    gene_plot <- FeaturePlot(dataset, features = input$gene, reduction = "umap") +
      ggtitle(paste("Expression of", input$gene, "in", dataset_name)) +
      theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
    celltype_plot <- DimPlot(dataset, reduction = "umap", group.by = celltype_col, label = TRUE) +
      ggtitle(paste("Cell Type Annotations in", dataset_name)) +
      theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
    plots <- c(plots, list(gene_plot, celltype_plot))
  }
  if (length(warnings) > 0) {
    feedbackWarning("gene", show = TRUE, text = paste(warnings, collapse = " | "))
  } else {
    feedbackWarning("gene", show = FALSE)
  }
  if (length(plots) == 0) {
    plot.new()
    title("This gene is not in any dataset.")
    return()
  }
  wrap_plots(plots, ncol = 2)})

#violin plot for GE
output$violinPlot <- renderPlot({
  req(input$update)
  selected_datasets <- datasetInput()
  plots <- list()
  
  for (dataset_name in names(selected_datasets)) {
    dataset <- selected_datasets[[dataset_name]]
    
    # Skip if gene not in dataset
    if (!(input$gene %in% rownames(dataset))) next
    
    # Check for cell type column
    celltype_col <- if ("celltype" %in% colnames(dataset@meta.data)) {
      "celltype"
    } else if ("Cell_Type" %in% colnames(dataset@meta.data)) {
      "Cell_Type"
    } else if ("CellType" %in% colnames(dataset@meta.data)) {
      "CellType"
    } else {
      next
    }
    
    # Subset data to only nonzero expression for visual clarity
    expr_vals <- FetchData(dataset, vars = input$gene, layer="data")[[1]]
    if (all(expr_vals == 0)) next  # skip if gene not expressed
    
    # Generate violin plot
    vln <- tryCatch({
      VlnPlot(dataset,
              features = input$gene,
              group.by = celltype_col,
              pt.size = 0.1,
              assay = DefaultAssay(dataset),
              layer = "data") +  # use log-normalized data
        theme_minimal() +
        labs(title = paste("Violin Plot of", input$gene, "in", dataset_name)) +
        theme(
          axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
          plot.title = element_text(hjust = 0.5, size = 14)
        )
    }, error = function(e) {
      ggplot() + annotate("text", x = 0.5, y = 0.5, label = paste("Violin failed for", dataset_name))
    })
    
    plots[[dataset_name]] <- vln
  }
  
  if (length(plots) > 0) {
    wrap_plots(plots, ncol = 1)
  } else {
    plot.new()
    title("No valid violin plot to display.")
  }
})
# ---- Reference Mapping Page ----
#this is only to match the predicted_id of bm with other tissue cell types( coz they had jumbled annotaion)
# Reactive function to return selected query object
getQueryObject <- reactive({
  req(input$query_dataset)
  query_choices[[input$query_dataset]]
  memo_readRDS(file.path("data", query_choices[[input$query_dataset]]))
})
# Define dynamic cell type dropdown from predicted.id in the selected query
output$select_celltype <- renderUI({
  obj <- getQueryObject()
  ids <- sort(unique(obj$predicted.id))
  selectInput("celltype", "Select Cell Type Annotation from Reference:",
              choices = ids, selected = ids[1])
})

# Now this works in renderPlot
output$highlight_all <- renderPlot({
  req(input$celltype)
  query_obj <- getQueryObject()
  bm <- memo_readRDS(file.path("data", dataset_choices["Bone Marrow"]))
  
  
  
  
  # b) highlight bone marrow cells
  
  plots <- list( 
    DimPlot(bm, reduction = "umap", 
            cells.highlight = WhichCells(bm, expression = CellType == input$celltype), 
            cols.highlight = "red")+ ggtitle("Bone Marrow (Reference)"))
  
  # c) loop through each query tissue and highlight matching cells
  for (tissue in names(query_choices)) {
    obj <- memo_readRDS(file.path("data", query_choices[[tissue]]))
    
    if (input$query_dataset == tissue) {
      # Sanitize cell type name to match column name
      safe_celltype <- gsub("[^[:alnum:]_]", ".", input$celltype)
      pred_col <- paste0("prediction.score.", safe_celltype)
      
      if (pred_col %in% colnames(obj@meta.data)) {
        plot <- FeaturePlot(obj, features = pred_col, reduction = "umap",
                            min.cutoff = "q5") +
          scale_color_gradientn(colors = c("grey90", "red")) +
          ggtitle(paste(tissue, "-", input$celltype)) +
          theme(legend.position = "right")
      } else {
        plot <- DimPlot(obj, reduction = "umap") +
          ggtitle(paste(tissue, "- No prediction score")) +
          theme(legend.position = "none")
      }
    } else {
      # Greyed-out background plot for non-selected queries
      plot <- DimPlot(obj, reduction = "umap", group.by = "predicted.id") +
        scale_color_manual(values = rep("grey90", length(unique(obj$predicted.id)))) +
        ggtitle(paste(tissue, "- Not selected")) +
        theme(legend.position = "none")
    }
    
    plots <- append(plots, list(plot))
  }
  
  wrap_plots(plots, ncol = 2)
})

#d) mapping score for tissue
output$mapping_plot <- renderPlot({
  obj <- getQueryObject()
  FeaturePlot(obj, features = "MappingScores", reduction = "umap", min.cutoff = 0) +
    scale_color_gradientn(colors = c("white", "yellow", "orange", "red")) +
    ggtitle("Mapping Scores(no cutoff)")
})

#e) predicted cell types for every tssue
output$predicted_umap <- renderPlot({
  obj <- getQueryObject()
  DimPlot(obj, group.by = "predicted.id", reduction = "umap", label = TRUE, repel = TRUE) +
    ggtitle("Predicted Cell Type IDs")
})


#f) dot plot of manual annotaon v/s predicted annottaion
output$dot_plot <- renderPlot({
  obj <- getQueryObject()
  df <- data.frame(CellType = obj$CellType, PredictedID = obj$predicted.id, MappingScore = obj$MappingScores)
  count_df <- df %>% group_by(CellType, PredictedID) %>% summarise(Count = n(), .groups = "drop")
  pct_df <- count_df %>% group_by(CellType) %>% mutate(Percent = 100 * Count / sum(Count)) %>% ungroup()
  score_df <- df %>% group_by(CellType, PredictedID) %>% summarise(MappingScore = mean(MappingScore), .groups = "drop")
  dot_df <- merge(pct_df, score_df, by = c("CellType", "PredictedID"))
  
  ggplot(dot_df, aes(x = PredictedID, y = CellType)) +
    geom_point(aes(size = Percent, fill = MappingScore), shape = 21, color = "black") +
    scale_size_continuous(name = "Percent of Manually Annotated cell types", range = c(1, 7)) +
    scale_fill_gradientn(colors = c("yellow", "orange", "red", "darkred"),
                         limits = c(0.5, 1), oob = squish,
                         name = "Average Mapping Score") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
          axis.text.y = element_text(size = 10),
          axis.title = element_text(size = 12)) +
    labs(title = "Predicted v/s Manually Annotated Cell Types",
         x = "Predicted Cell Type based on Bone Marrow",
         y = "Manual Annotation")
})
}

# Run app
shinyApp(ui, server)

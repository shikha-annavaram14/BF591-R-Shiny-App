library(shiny)
library(bslib)
library(ggplot2)
library(tidyverse)
library(colourpicker)
library(magrittr)
library(DT)
library(pheatmap)
library(DESeq2)
library(fgsea)


ui <- fluidPage(
  titlePanel("Dataset from https://www.science.org/doi/10.1126/sciadv.abd1160#sec-4"),
  tabsetPanel(
    tabPanel(
      "Samples",
      sidebarLayout(
        sidebarPanel(
          "This tab displays the metadata related this to this dataset in a digestable manner. The input for this tab is a csv or tsv sample matrix. ",
          fileInput("sample_info_file", "Upload Sample Information Matrix (CSV)", accept = c(".csv", ".tsv")), width = 3
        ),
        mainPanel(
          tabsetPanel(
            tabPanel(
              "Summary",
              DTOutput("sample_summary"),
              "A summarized version of the metadata"
            ),
            tabPanel(
              "Table",
              DTOutput("sample_table"), 
            ),
            tabPanel(
              "Plots",
                plotOutput("sample_plot", height = 700)
                
              )
            )
          )
        )
      )
    ,
    
    # Counts Matrix Exploration
    tabPanel(
      "Counts",
      sidebarLayout(
        sidebarPanel(
          "This section will take the counts matrix and display it based on user chosen thresholds to show the most significant genes.",
          fileInput("counts_matrix_file", "Upload Normalized Counts Matrix (CSV)", accept = c(".csv", ".tsv")),
          sliderInput("variance_threshold", "Minimum Variance Percentile:", 0, 100, 20),
          sliderInput("nonzero_threshold", "Maximum Zero Samples:", 0, 18, 16), 
          width = 3
        ),
        mainPanel(
          tabsetPanel(
            tabPanel(
              "Summary",
              uiOutput("counts_summary"),
              DTOutput("counts_table")
            ),
            tabPanel(
              "Scatter Plots",
              plotOutput("medvsvar"), 
              plotOutput("medvszeros")
            ),
            tabPanel(
              "Heatmap",
                plotOutput("heatmap", height = 700)
              )
            ,
            tabPanel(
              "PCA",
              sidebarLayout(
                sidebarPanel(
                  selectInput("x_pc", "X-axis Principal Component:", choices = NULL),
                  selectInput("y_pc", "Y-axis Principal Component:", choices = NULL),
                  selectInput("group_by", "Group by Sample Info", choices = NULL), 
                  width = 3
                ),
                mainPanel(
                  plotOutput("pca_plot")
                )
              )
            )
          )
        )
      )
    ),
    
    # Differential Expression
    tabPanel(
      "DE",
        mainPanel(
          tabsetPanel(
            tabPanel(
              "DE Results",
              DTOutput("de_table"),
              downloadButton("download_de_results", "Download DE Results for FGSEA Analysis"),
              "This section takes in the normalized counts from the previous tab and produces the differential expression results. It allows you
              to download this analysis to perform FGSEA analysis."
            ), 
            tabPanel(
              "Volcano Plot",
              plotOutput("volcano_plot")
            )
          )
        )
      )
    ,
    
    # Gene Set Enrichment Analysis
    tabPanel(
      "GSEA",
      sidebarLayout(
        sidebarPanel(
          fileInput("fgsea_file", "Upload FGSEA Analysis from DGE (TSV)", accept = c(".csv", ".tsv")), width = 3,
          "This tab performs gene set enrichment analysis. It takes in a csv or tsv of fgsea results performed on DE data."
        ),
        mainPanel(
        tabsetPanel(
        tabPanel(
          "Top Results",
          sidebarLayout(
            sidebarPanel(
              sliderInput("top_pathways", "Number of Top Pathways:", min = 1, max = 50, value = 10), width = 3
            ),
            mainPanel(
              plotOutput("pathway_barplot", height = 700)
            )
          )
        ),
        tabPanel(
          "Table",
          sidebarLayout(
            sidebarPanel(
              sliderInput("pval_filter", "Adjusted p-value Threshold:", min = 0, max = 0.1, step = 0.0001, value = 0.05),
              radioButtons("nes_filter", "Select Pathway Type:", 
                           choices = list("All" = "all", "Positive NES" = "positive", "Negative NES" = "negative")),
              downloadButton("download_table", "Download Results"), 
              width = 3
            ),
            mainPanel(
              DTOutput("pathway_table")
            )
          )
        ),
        tabPanel(
          "Plots",
          sidebarLayout(
            sidebarPanel(
              sliderInput("scatter_pval_filter", "Adjusted p-value Threshold:", min = 0, max = 0.1, step = 0.0001, value = 0.05), width = 3 
            ),
            mainPanel(
              plotOutput("scatter_plot")
            )
          )
        )
      )
    )
  )
)
)
)


  

server <- function(input, output, session) {
options(shiny.maxRequestSize=30*1024^2)
#SAMPLE INFO EXPLORATION
  #load sample matrix in csv format
  load_sample_data <- reactive({
    req(input$sample_info_file)
    #renaming index 1 column to Sample
    read.csv(input$sample_info_file$datapath) %>% as.data.frame() %>% dplyr::select(Sample = 1, everything())
  })
  
  
  #tab with a summary of the table that includes a summary of the type and values in each column
  output$sample_summary <- renderDT({
    req(load_sample_data())
    table <- load_sample_data()
    #making them factors
    table$Treatment_Time_Point <- as.factor(table$Treatment_Time_Point)
    table$Cell_Type <- as.factor(table$Cell_Type)
    table <- dplyr::select(table, -1)
    #create summary table
    summary_table <- data.frame(
      `Column Name` = colnames(table),
      `Type` = sapply(table, class),
      `Distinct Values` = sapply(table, function(col) {
        if (is.numeric(col)) {
          paste0("Mean: ", round(mean(col, na.rm = TRUE), 2),
                 " (SD: ", round(sd(col, na.rm = TRUE), 2), ")")
        } else {
          paste(unique(col), collapse = ", ")
        }
      }),
      check.names = FALSE
    )
    
    datatable(summary_table, options = list(scrollX = TRUE))
  
  })
  #tab with a data table displaying the sample information, with sortable columns
  output$sample_table <- renderDT({
    req(load_sample_data())
    data.table::as.data.table(load_sample_data())
  })
  #plot to represent data - no numerical categories so plot histogram on how data separates
  output$sample_plot <- renderPlot({
    req(load_sample_data())
    data <- load_sample_data()
    ggplot(data, aes(x = Treatment_Time_Point, fill = Cell_Type)) +
      geom_bar(position = "stack", color = "black") +
      theme_minimal() +
      labs(
        title = "Samples by Treatment Time Point and Cell Type",
        x = "Treatment Time Point",
        y = "Number of Samples",
        fill = "Cell Type"
      )
  })
  
#COUNTS MATRIX EXPLORATION
  #input normalized count matrix
  load_counts_data <- reactive({
    req(input$counts_matrix_file)
    read.csv(input$counts_matrix_file$datapath) %>% as.data.frame() %>% column_to_rownames(var = 'gene') %>% dplyr::select(-1)
  })
  #input contols that filter out genes based on statistical properties
    #slider to include genes with at least X percentile variance
    #slider to include genes with at least X samples that are non zero
  #tab with table or text summarizing effect of filtering
    #number of samples, total number of genes, number and % of genes passing current filter, number and % of genes not passing current filter
  filtered_data <- reactive({
    req(load_counts_data())
    data <- load_counts_data()
    #calculate variance for each row and count number of zeros in each row
    gene_variance <- apply(data, 1, var)
    zero_counts <- rowSums(data == 0)
    #calculate variance threshold
    variance_threshold <- quantile(gene_variance, probs = input$variance_threshold / 100)
    #create filter
    filtered_genes <- which(
      gene_variance >= variance_threshold &
        zero_counts <= input$nonzero_threshold
    )
    
    data[filtered_genes, , drop = FALSE]
    
  })
  #counts text summary
  output$counts_summary <- renderUI({
    req(load_counts_data())
    original <- load_counts_data()
    filtered <- filtered_data()
    
    #calculate all stats needed for display
    total_samples <- ncol(original)
    total_genes <- nrow(original)
    filtered_genes <- nrow(filtered)
    percent_pass <- (filtered_genes / total_genes) * 100
    percent_fail <- 100 - percent_pass
    
  #used html to make a more readable display
   HTML(
    paste(
      "<div style='font-size: 18px;'>",
      "Total Samples: ", total_samples, "<br>",
      "Total Genes: ", total_genes, "<br>",
      "Number of Genes Passing Current Filter: ", filtered_genes, sprintf(" (%.2f%%)", percent_pass), "<br>",
      "Number of Genes Not Passing Current Filter: ", total_genes - filtered_genes, sprintf(" (%.2f%%)", percent_fail),
      "</div>"
    ))
    
  })
  #render filtered table
  output$counts_table <- renderDT({
    req(filtered_data())
    datatable(filtered_data())
  })
    
  #tab with diagnostic scatter plots
    #genes passing filter are marked in darker color, filtered out are lighter color
    #median count vs variance (log scale)
    #median count vs number of zeros
  scatter_filter <- reactive({
    req(load_counts_data())
    data <- load_counts_data()
    #add filter for median
    gene_var <- apply(data, 1, var)
    gene_median <- apply(data, 1, median)
    zero_counts <- rowSums(data == 0)
    
    var_threshold <- quantile(gene_var, input$variance_threshold /100)
    zero_threshold <- input$nonzero_threshold
    #create filter 
    pass_filter <- gene_var >= var_threshold & zero_counts <= (ncol(data) - zero_threshold)
    list(
      variance = gene_var,
      median = gene_median,
      zeros = zero_counts, 
      pass_filter = pass_filter
    )
  })
    output$medvsvar <- renderPlot({
      req(scatter_filter())
      filtered <- scatter_filter()
      #create dataframe that only contains columns needed for plotting
      df <- data.frame(Median = filtered$median, Variance = filtered$variance, PassFilter = filtered$pass_filter) %>% filter(Variance > 0)
      ggplot(df, aes(x = log10(Median), y = Variance, color = PassFilter)) + geom_point() +
        scale_color_manual(values = c("TRUE" = "blue", "FALSE" = "gray")) +
        scale_y_log10() +
        labs(
          title = "Median Count vs Variance",
          x = "log10(Median Count)",
          y = "log10(Variance)"
        )
    })
    
    output$medvszeros <- renderPlot({
      req(scatter_filter())
      filtered <- scatter_filter()
      data <- data.frame(Median = filtered$median, Zeros = filtered$zeros, PassFilter = filtered$pass_filter)
      ggplot(data, aes(x = log10(Median), y = Zeros, color = PassFilter)) + geom_point()+
        scale_color_manual(values = c("TRUE" = "blue", "FALSE" = "gray")) +
        theme_minimal() +
        labs(
          title = "Median Count vs Number of Zeros",
          x = "log10(Median Count)",
          y = "Number of Zeros"
        )
    })

  #tab with clustered heatmap of counts remaining after filtering
    #enable log-transforming counts
    #color bar
    output$heatmap <- renderPlot({
      req(filtered_data())
      data <- filtered_data()
      #log10+1 so there is no log(0)
      pheatmap(log10(data+1), cluster_rows = TRUE, cluster_cols = TRUE, scale = 'row', legend = TRUE)
    })
    
  #tab with pca, include percent variance explained by pc
    #allow used to select which pcs to plot
    
    pca_results <- reactive({
      req(filtered_data())
      data <- filtered_data()
      pca <- prcomp(t(data), scale. = TRUE)
      pca
    })
    #An observer is like a reactive expression in that it can read reactive values and call reactive expressions, 
    #and will automatically re-execute when those dependencies change.
    #in this case the observe is used to update when certain inputs are selected by the user
    #selecting the metadata column to group the pca by 
    observe({
      req(load_sample_data())
      metadata_columns <- colnames(load_sample_data())
      updateSelectInput(session, "group_by", choices = metadata_columns, selected = NULL)
    })
    #selecting which pcs to plot
    observe({
      req(pca_results())
      num_pcs <- length(pca_results()$sdev)
      pc_choices <- paste0("PC", 1:num_pcs)
      updateSelectInput(session, "x_pc", choices = pc_choices, selected = "PC1")
      updateSelectInput(session, "y_pc", choices = pc_choices, selected = "PC2")
    })
    
    output$pca_plot <- renderPlot({
      req(pca_results(), input$x_pc, input$y_pc, input$group_by, load_sample_data())
      
      pca <- pca_results()
      pc_data <- as.data.frame(pca$x)
      pc_data$Sample <- rownames(pc_data)
      meta <- load_sample_data()
      pc_data <- merge(pc_data, meta, by = 'Sample')
      
      
      var_explained <- pca$sdev^2 / sum(pca$sdev^2) * 100
      x_label <- sprintf("%s (%.2f%% variance explained)", input$x_pc, var_explained[as.numeric(gsub("PC", "", input$x_pc))])
      y_label <- sprintf("%s (%.2f%% variance explained)", input$y_pc, var_explained[as.numeric(gsub("PC", "", input$y_pc))])
      
      # group by metadata columns
      ggplot(pc_data, aes(x = !!sym(input$x_pc), y = !!sym(input$y_pc), color = !!sym(input$group_by))) +
        geom_point() +
        labs(x = x_label, y = y_label, title = "PCA Scatter Plot")
    })
#DIFFERENTIAL EXPRESSION
  #input: results of DE, perform in DESeq2
  de_counts <- reactive({
    req(load_counts_data())
    data <- load_counts_data() %>% as.matrix()
  })
  de_results <- reactive({
    req(de_counts(), load_sample_data())
    meta <- load_sample_data() %>% column_to_rownames(var = "Sample")
    dds <- DESeqDataSetFromMatrix(countData = round(de_counts()), colData = meta, design = ~1)
    dds <- DESeq(dds)
    results <- results(dds) %>% as.data.frame()
    results
  })
  #tab with table with sortable feature
  output$de_table <- renderDT({
    req(de_results())
    datatable(de_results())
  })
  output$download_de_results <- downloadHandler(
    filename = function(){
      paste0("differential_expr_results_", Sys.Date(), ".csv")
    },
    content = function(file){
      req(de_results())
      write.csv(de_results(), file, row.names = TRUE)
    }
  )
  #tab with volcano plot
  output$volcano_plot <- renderPlot({
    req(de_results())
    results <- de_results() %>% filter(padj < 0.05 & abs(log2FoldChange) > 1)
    ggplot(results, aes(x = log2FoldChange, y = -log10(padj))) +
      geom_point() +
      labs(title = "DE Volcano Plot", x = "log2FC", y = "-log10(padj)")
  })
  
#FGSEA
  fgsea_data <- reactive({
    req(input$fgsea_file)
    read_tsv(input$fgsea_file$datapath, show_col_types = FALSE) %>% as.data.frame()
  })
  
  output$pathway_barplot <- renderPlot({
    req(fgsea_data())
    res <- fgsea_data() %>% filter(!is.na(NES))
    
    top_pos <- res %>% filter(NES > 0) %>% arrange(desc(NES)) %>% slice_head(n = input$top_pathways)
    top_neg <- res %>% filter(NES < 0) %>% arrange(NES) %>% slice_head(n = input$top_pathways)
    top_paths <- bind_rows(top_pos, top_neg) %>%
        mutate(
          Pathway = str_replace_all(pathway, "_", " ") %>% str_wrap(),
          Direction = ifelse(NES > 0, "Positive", "Negative")
        )

    ggplot(top_paths, aes(x = reorder(Pathway, NES), y = NES, fill = Direction)) +
        geom_col() +
        coord_flip() +
        scale_fill_manual(values = c("Positive" = "red", "Negative" = "blue")) +
        labs(x = 'Pathways', y = "Normalized Enrichment Score (NES)") +
        theme_minimal(base_size = 12) +
        theme(
          axis.text.y = element_text(size = 7),
          plot.title = element_text(size = 8),
          axis.title.x = element_text(size= 6)
        ) + guides(fill="none")
  })
  
  output$pathway_table <- renderDT({
    req(fgsea_data())
    res <- fgsea_data()
    filtered <- res[res$padj < input$pval_filter, ]
    if (input$nes_filter == "positive") {
      filtered <- filtered[filtered$NES > 0, ]
    } else if (input$nes_filter == "negative") {
      filtered <- filtered[filtered$NES < 0, ]
    }
    filtered <- na.omit(filtered)
    datatable(filtered, options = list(pageLength = 10))
  })
  
  output$download_table <- downloadHandler(
    filename = function() { paste0("fgsea_results_", Sys.Date(), ".csv") },
    content = function(file) {
      req(fgsea_data())
      res <- fgsea_data()
      filtered <- res[res$padj < input$pval_filter, ]
      if (input$nes_filter == "positive") {
        filtered <- filtered[filtered$NES > 0, ]
      } else if (input$nes_filter == "negative") {
        filtered <- filtered[filtered$NES < 0, ]
      }
      filtered <- na.omit(filtered)
      write.csv(filtered, file, row.names = FALSE)
    }
  )
  
  output$scatter_plot <- renderPlot({
    req(fgsea_data())
    res <- fgsea_data()
    filtered <- res
    filtered$color <- ifelse(filtered$padj < input$scatter_pval_filter, "blue", "gray")
    
    ggplot(filtered, aes(x = NES, y = -log10(padj), color = color)) +
      geom_point(alpha = 0.6) +
      scale_color_manual(values = c("blue", "gray"), guide = "none") +
      theme_minimal() +
      labs(title = "NES vs Adjusted P-Value", x = "NES", y = "-Log10(padj)")
  })
}


# Run the application
shinyApp(ui = ui, server = server)



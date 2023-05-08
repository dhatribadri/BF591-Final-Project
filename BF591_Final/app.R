#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(plotly)
library(tidyverse)
library(DT)
library(colourpicker) 
library(shinyWidgets)
library(fgsea)
library(gplots)
library(RColorBrewer)


# put all observe events together

gmt <- "~/Desktop/BF591/untitled folder/c2.cp.v7.5.1.symbols.gmt"

# Define UI for app
ui <- fluidPage(
  titlePanel("BF591: RShiny App"),
  h2("Dhatri Badri"),
  
  # Define tabs
  tabsetPanel(
    tabPanel("Samples",
             sidebarLayout(
               sidebarPanel(
                 fileInput(inputId ="file_input",  label= "Upload metadata file"), # should i just do accept=c(".csv", ".tsv")
                 actionButton("submit_button", "Submit")
               ),
               mainPanel(
                 tabsetPanel(
                   tabPanel("Summary",
                            DT::dataTableOutput("summary_output")
                   ),
                   tabPanel("Table",
                            dataTableOutput("table_output")
                   ),
                   tabPanel("Diagnostic Plot",
                            plotlyOutput("plots_output")
                   )
                 )
               )
             )
    ),
    
    tabPanel("Counts",
             sidebarLayout(
               sidebarPanel(
                 fileInput(inputId ="file_input_2",  label= "Upload counts matrix file", accept=c(".csv")), # should i just do accept=c(".csv", ".tsv")
                 sliderInput(inputId="variance_percentile", label="Variance percentile cutoff:",
                             min = 0, max = 100, value = 50),
                 sliderInput(inputId="nonzero_samples", label="Non-zero samples cutoff:",
                             min = 0, max = 60, value = 1),
                 actionButton("submit_counts", "Submit")
               ),
               mainPanel(
                 tabsetPanel(
                   tabPanel("Summary",
                            DT::dataTableOutput("summary_table")
                   ),
                   tabPanel("Scatter plot",
                            plotlyOutput("scatter_plot"),
                            plotlyOutput("scatter_plot_2")
                   ),
                   tabPanel("Heatmap",
                            plotOutput("plot_heatmap")
                   ),
                   tabPanel("PCA",
                            fluidRow(
                              column(4,
                                     selectInput("x_selector", "X-axis Principal Component:", 
                                                 choices = c("PC1", "PC2", "PC3"), selected = "PC1")
                              ),
                              column(4,
                                     selectInput("y_selector", "Y Axis Principal Component:", 
                                                 choices = c("PC1", "PC2", "PC3"), selected = "PC2")
                              )
                            ),
                            plotlyOutput("pca_plot")
                   )
                 )
               )
             )
    ),
    tabPanel("DE",
             sidebarLayout(
               sidebarPanel(
                 fileInput(inputId ="file_input_3",  label= "Upload differential expression file", accept=c(".csv")),
                 p("A volcano plot can be generated with ", tags$b("log2 fold-change"), " on the x-axis and ", tags$b("p-adjusted"), " on the y-axis."),
                 
                 radioButtons("x_axis", "Choose the column for the x-axis", selected="log2FoldChange", choices=c("baseMean", "log2FoldChange", "lfcSE", "stat","pvalue","padj")),
                 radioButtons("y_axis", "Choose the column for the y-axis", selected="padj", choices=c("baseMean", "log2FoldChange", "lfcSE", "stat","pvalue","padj")), 
                 
                 colourInput("bpc", label= "Base point color", value= "#22577A"),
                 colourInput("hpc", label= "Highlight point color", value= "#FFCF56"),
                 
                 sliderInput(inputId = "pval", label = "Select the magnitude of the p adjusted coloring:", min=-300,max=0,value=-150, step=1),
                 actionButton("createplot", label = "Plot", icon = icon("dna"))
               ), 
               mainPanel(
                 tabsetPanel(
                   tabPanel("DE Results",
                            DT::dataTableOutput(outputId = "de_results")
                   ),
                   tabPanel("Plot",
                            plotlyOutput(outputId = "volcano")
                   ),
                   tabPanel("Table",
                            dataTableOutput(outputId = "de_table_output")
                   )
                 )
               )
             )
    ),
    tabPanel("FGSEA",
             sidebarLayout(
               sidebarPanel(
                 fileInput(inputId ="file_input_4",  label= "Upload differential expression file for GSEA"),
                 actionButton("submit_fgsea", "Submit"),
                 tags$div(style = "margin-top: 10px"), # add some space between the submit button and the slider,
                 sliderInput(inputId = "adj_pval_slider", label = " filter table by adjusted p-value :", min = -35, max = 0, value = 10),
                 #p("only for Table Results and Scatter plot")
               ), 
               mainPanel(
                 tabsetPanel(
                   tabPanel("GSEA Results",
                            DT::dataTableOutput(outputId = "gsea_results")
                   ),
                   tabPanel("Barplot of fgsea NES ",
                            sidebarLayout(
                              sidebarPanel(
                                sliderInput(inputId = "num_top_pathways", label = "Number of Top Pathways to Plot:", min = 1, max = 100, value = 10),
                                actionButton("submit_fgsea_barplot", "Submit")
                              ),
                              mainPanel(
                                plotlyOutput(outputId = "fgsea_bar_plot")
                              )
                            )
                   ),
                   tabPanel("Table Results",
                            sidebarLayout(
                              sidebarPanel(
                                radioButtons("nes_direction", "Choose to select all, positive or negative NES pathways", selected="all", choices=c("all", "positive", "negative")),
                                #actionButton("submit_fgsea_table", "Submit"),
                                downloadButton(outputId = "download_fgsea_table", label = "Download Filtered Table")
                              ),
                              mainPanel(
                                DT::dataTableOutput(outputId = "fgsea_table_output")
                              )
                            )
                   ), 
                   tabPanel("Scatter Plot",
                            plotlyOutput(outputId = "fgsea_scatter_plot")
                   )
                 )
               ),
             )
    ),
  )
)

# Define server for app
server <- function(input, output) {
  
  #shinyjs::logjs("hello")
  options(shiny.maxRequestSize = 50*1024^2) # Set maximum upload size to 50 MB
  
  
  # Read in uploaded file as a data frame
  load_data <- reactive({
    req(input$file_input) # require input before proceeding
    
    # Check if uploaded file is either CSV or TSV
    if (!grepl("\\.(csv)$", input$file_input$name, ignore.case = TRUE)) {
      stop("Please upload a CSV file.")
    }
    
    # Read in file as a data frame
    tryCatch(
      {
        if (grepl("\\.csv$", input$file_input$name, ignore.case = TRUE)) {
          read.csv(input$file_input$datapath, header = TRUE)
        } #else {
        #read.table(input$file_input$datapath, header = TRUE, sep = "\t", check.names = FALSE)
        #}
        
      },
      error = function(e) {
        stop("Failed to read the file. Please check if it's well-formatted.")
      }
    )
  })
  
  
  samples_summary <- function(samples_data) {
    # Create a data frame containing column name, type, and mean/distinct values
    names <- colnames(load_data())
    types <- sapply(load_data(), class)
    means_or_values <- sapply(load_data(), function(x) {
      if (class(x) %in% c("integer", "numeric", "double")) { #how do I fix this?
        mean_value <- round(mean(x))
        sd_value <- round(sd(x))
        paste0(mean_value, " (+/-", sd_value, ")")
      } else {
        paste0(paste(unique(x), collapse = ", "))
      }
    })
    
    # Combine into a data frame and return as table
    summary_df <- data.frame("Column name" = names,
                             "Type" = types,
                             "Mean (sd) or Distinct Values" = means_or_values, check.names = FALSE)
    return(summary_df)
  }
  
  
  samples_diagnostic_plot <- function(samples_data) {
    ggplot(samples_data, aes(x = Diagnosis , y = age_of_death, fill= Diagnosis)) +
      geom_violin() +
      labs(title = "Age of death based on diagnosis", x = "Diagnosis", y = "Age of Death")  # add labels to the plot
  }                                                                                                                               
  
  observeEvent(input$submit_button, {
    # Generate summary output
    output$summary_output <- renderDataTable({
      samples_summary(load_data())
    })
    
    # Define output for table component
    output$table_output <- renderDataTable({
      load_data()
    })
    
    output$plots_output <- renderPlotly({
      # Create the plot
      samples_diagnostic_plot(load_data())
    })
    
  })
  
  
  # Read in uploaded file when submit button is clicked
  load_norm_counts <- reactive({
    req(input$file_input_2) # require input before proceeding
    infile <- input$file_input_2
    if (is.null(infile)) {
      return(NULL)
    }
    read.csv(infile$datapath, header = TRUE)
  })
  
  # add variance_filter and non_zero_filter
  counts_filter <- function(counts_data, variance_filter, non_zero_filter) {
    # Gene filtering and summary
    filtered_data <- counts_data %>%
      select(-1) %>%
      filter(
        apply(., 1, var) >= quantile(apply(., 1, var), variance_filter/100),
        rowSums(.[,-1] > 0) >= non_zero_filter
      )
    
    num_samples <- ncol(counts_data) - 1
    total_genes <- nrow(counts_data) - 1
    passed_filter <- nrow(filtered_data)
    failed_filter <- total_genes - passed_filter
    
    counts_summary <- data.frame(
      "Number of samples" = num_samples,
      "Total number of genes" = total_genes,
      "Number and % of genes passing the filter" = paste0(passed_filter, " (", 
                                                          paste0(round(100 * passed_filter/total_genes, 2), "%)")),
      "Number and % of genes not passing the filter" = paste0(failed_filter, " (", 
                                                              paste0(round(100 * failed_filter/total_genes, 2), "%)")), check.names = FALSE
    )
    
    return(counts_summary)
  }
  
  
  observeEvent(input$submit_counts, {
    # Generate summary output
    output$summary_table <- renderDataTable({
      counts_filter(load_norm_counts(),input$variance_percentile,input$nonzero_samples)
    })
  })
  
  
  counts_scatter_plot_1 <- function(counts_data, variance_filter){
    
    counts_data$median_count =  apply(counts_data[-1], 1, median)
    counts_data$variance_count = apply(counts_data[-1], 1, var)
    
    # Check if each gene passes the thresholds
    pass <- counts_data$variance_count >= quantile(counts_data$variance_count, variance_filter/100)
    
    # Create a new column in the data frame indicating whether each gene passed or failed
    counts_data$pass_fail <- ifelse(pass, "pass", "fail")
    
    mean_var_plot <- ggplot(counts_data, mapping=aes(x=log10(variance_count), y=median_count, color=pass_fail)) +
      geom_point(alpha=5/10) + 
      labs(title="median count vs log10(variance)",x="log10(variance)", y="Median Counts") 
    return(mean_var_plot)
    
  }
  
  counts_scatter_plot_2 <- function(counts_data, non_zero_filter) {
    
    counts_data$median_count =  apply(counts_data[-1], 1, median)
    counts_data$non_zero = apply(counts_data[-1] != 0, 1, sum)
    # Check if each gene passes the thresholds
    pass <- (counts_data$non_zero >= non_zero_filter)
    
    # Create a new column in the data frame indicating whether each gene passed or failed
    counts_data$pass_fail <- ifelse(pass, "pass", "fail")
    
    mean_zero_plot <- ggplot(counts_data, mapping=aes(x=non_zero, y=median_count, color=pass_fail)) +
      geom_point(alpha=5/10) + 
      labs(title="median count vs number of zeros",x="Num of zeros", y="Median Counts") 
    return(mean_zero_plot)
  }
  
  plot_heatmap <- function(counts_data, variance_filter, non_zero_filter) {
    
    counts_data$variance_count= apply(counts_data[-1], 1, var)
    counts_data$non_zero = apply(counts_data[-1] != 0, 1, sum)
    # Check if each gene passes the thresholds
    pass <- (counts_data$variance_count >= variance_filter) & (counts_data$non_zero >= non_zero_filter)
    
    # Create a new column in the data frame indicating whether each gene passed or failed
    counts_data$pass_fail <- ifelse(pass, "pass", "fail")
    colnames(counts_data) [1] <- "gene"
    
    # Filter out the rows that passed and select the relevant columns
    filtered_data <- subset(counts_data, pass == TRUE, select = c("gene", names(counts_data)[2:70]))
    
    # Convert the filtered data to a numeric matrix
    filtered_matrix <- data.matrix(filtered_data[, 2:70])
    
    #fd_mat <- as.matrix(fd)
    # Create the heatmap using the heatmap.2 function
    heatmap_color <- brewer.pal(11, 'RdBu')
    #heat_map <- heatmap.2(filtered_matrix, col = heatmap_color, trace = "none", scale = "row", dendrogram = "none", key = TRUE, keysize = 1.5, key.title = NA, symkey = FALSE, density.info = "none", cexCol = 0.8, margins = c(10,10))
    
    heat_map <- heatmap(filtered_matrix, col= heatmap_color,  scale = "row", margins = c(10, 10))
    
    # Add a color bar legend to the plot
    legend("topright", fill = heatmap_color, xpd = TRUE, 
           legend = seq(min(filtered_matrix), max(filtered_matrix), length.out = length(heatmap_color)),
           title = "Color Key", inset=c(0.01, 0.1))
    
    return(heat_map)
  }
  
  
  counts_plot_pca <- function(norm_counts, pc1, pc2) {
    # Run PCA
    pca <- prcomp(t(norm_counts[-1]))
    
    # Calculate variance explained by each PC
    var_explained <- round(100 * pca$sdev^2 / sum(pca$sdev^2), 2)
    
    # Create plot
    pc1_num <- as.numeric(gsub("\\D", "", pc1))
    pc2_num <- as.numeric(gsub("\\D", "", pc2))
    plot_data <- data.frame(
      x = pca$x[, pc1_num],
      y = pca$x[, pc2_num],
      condition = ifelse(grepl("^C", colnames(norm_counts[-1])), "C", "H")
    )
    
    plot <- ggplot(plot_data, aes(x = x, y = y, color = condition)) +
      geom_point() +
      labs(x = paste0(pc1, " (", var_explained[pc1_num], "%)"), 
           y = paste0(pc2, " (", var_explained[pc2_num], "%)")) +
      theme_classic()
    
    return(plot)
  }
  
  
  observeEvent(input$submit_counts, {
    # Define output for plot component
    output$scatter_plot <- renderPlotly({
      # Create the plot
      counts_scatter_plot_1(load_norm_counts(), input$variance_percentile)
    })
    
    output$scatter_plot_2 <- renderPlotly({
      # Create the plot
      counts_scatter_plot_2(load_norm_counts(), input$nonzero_samples)
    })
    
    output$plot_heatmap <- renderPlot({
      # Create the plot
      plot_heatmap(load_norm_counts(), input$variance_percentile, input$nonzero_samples) 
    } ,height = 600, width = 600)
    
    # Define output for plot component
    output$pca_plot <- renderPlotly({
      # Create the plot
      counts_plot_pca(load_norm_counts(), input$x_selector, input$y_selector)
    })
    
  })
  
  load_de <- reactive({
    req(input$file_input_3) # require input before proceeding
    infile <- input$file_input_3
    if (is.null(infile)) {
      return(NULL)
    }
    read.csv(infile$datapath, header = TRUE, stringsAsFactors = TRUE)
  })
  
  volcano_plot <- function(dataf, x_name, y_name, slider, color1, color2) {
    if(is.null(dataf)) {
      return(NULL)
    }     
    dataf<-na.omit(dataf)
    color_label = paste0("padj <", "(1*(10^", slider, "))")
    plot <- ggplot(dataf, aes(x = !!sym(x_name), y=-log10(!!sym(y_name)), color = padj< 1*10^(slider))) +
      geom_point() +
      scale_colour_manual(name = color_label, values = setNames(c(color1, color2, 'grey'),c(F, T, NA))) +
      theme(legend.position = "bottom")
    return(plot)
  }
  
  draw_table <- function(dataf, slider) {
    dataf <- na.omit(dataf)
    new_df <- dataf %>% 
      filter(padj < (1 * 10^slider)) %>% 
      mutate(pvalue = formatC(pvalue, digits = 7),
             padj = formatC(padj, digits = 7))
    return(new_df)
  }
  
  observeEvent(input$createplot, { 
    # Render results table
    output$de_results <- DT::renderDataTable({
      de_data <- load_de()
      if (is.null(de_data)) {
        return(NULL)
      }
      else {
        return(de_data)
      }
    })
    
    # Plot
    output$volcano <- renderPlotly({
      # Create the plot
      volcano_plot(load_de(), input$x_axis, input$y_axis, input$pval, input$bpc, input$hpc)
    })
    
    # Table
    output$de_table_output <- DT::renderDataTable({
      draw_table(load_de(), input$pval)
    })
  })
  
  # Read in uploaded file when submit button is clicked
  load_fgsea <- reactive({
    req(input$file_input_4) # require input before proceeding
    infile <- input$file_input_4
    if (is.null(infile)) {
      return(NULL)
    }
    read.csv(infile$datapath, header = TRUE)
  })
  
  
  run_gsea <- function(fgsea_results, gmt) {
    fgsea_file <- fgsea_results
    colnames(fgsea_file)[1] <- "gene"
    fgsea_results <- drop_na(fgsea_file) %>%
      distinct(gene, log2FoldChange, .keep_all = TRUE)
    fgsea_results_arranged <- arrange(fgsea_results, desc(log2FoldChange)) 
    ranks <- deframe(dplyr::select(fgsea_results_arranged, symbol, log2FoldChange))
    c2pathways <- gmtPathways(gmt)
    fgsea_res <- as_tibble(fgsea(c2pathways, ranks, minSize=15, maxSize=500))
    
    # remove _ in pathway
    fgsea_res <- fgsea_res %>%
      mutate(pathway = factor(pathway)) %>%
      mutate(pathway = str_replace_all(pathway, '_', ' '))
    
    return(fgsea_res)
  }
  
  fgsea_barplot <- function(fgsea_results, num_paths){
    # top positive pathways
    top_x_positive <- fgsea_results %>%
      arrange(desc(NES)) %>%
      filter(row_number() <= num_paths) %>%
      pull(pathway)
    
    # top negative pathways
    top_x_negative <- fgsea_results %>%
      arrange(NES) %>%
      filter(row_number() <= num_paths) %>%
      pull(pathway)
    
    # top  + and - pathways 
    filtered_fgsea_results <- fgsea_results %>%
      filter(pathway %in% c(top_x_positive, top_x_negative)) %>%
      mutate(pathway = factor(pathway)) %>%
      mutate(pathway_name = str_replace_all(pathway, '_', ' '))
    
    fgsea_plt <- filtered_fgsea_results %>%
      mutate(pathway_name = forcats::fct_reorder(pathway_name, NES)) %>%
      ggplot() +
      geom_bar(aes(x=pathway_name, y=NES, fill = NES > 0 ), stat='identity', show.legend = FALSE) +
      scale_fill_manual(values = c('TRUE' = 'blue', 'FALSE' = 'red')) + 
      theme_minimal(base_size = 7) +
      ggtitle('fgsea results for Hallmark MSigDB gene sets') +
      scale_x_discrete(labels = ~ str_wrap(., width = 70)) +
      ylab('Normalized Enrichment Score (NES)') +
      xlab('') +
      coord_flip()
    
    return(fgsea_plt)
  }
  
  adj_p_val <- reactive ({
    input$adj_pval_slider
  })
  
  # Filter data by adjusted p-value and NES direction
  # Define the function for filtering fgsea_results table
  fgsea_filtered_table <- function(fgsea_results, NES_direction, adj_pvalue_slider) {
    
    # Filter by adjusted p-value
    fgsea_results <- fgsea_results %>%
      filter(padj < 10^(adj_pvalue_slider)) %>%
      mutate(pathway = factor(pathway)) %>%
      mutate(pathway = str_replace_all(pathway, '_', ' '))
    
    # Filter by NES direction
    if (NES_direction == "positive") {
      fgsea_results <- fgsea_results %>%
        arrange(desc(NES)) %>%
        filter(NES > 0)# %>%
      #pull(pathway)
      #return(fgsea_results_adjusted)
    } else if (NES_direction == "negative") {
      fgsea_results <- fgsea_results %>%
        arrange(NES) %>%
        filter(NES < 0)# %>%
      #pull(pathway)
      #return(fgsea_results_adjusted)
    } #else if (NES_direction == "all") {
    #return(fgsea_results_adjusted)
    #}
    
    return(fgsea_results)
  }
  
  
  fgsea_scatter_plot_fn <- function(fgsea_results, adj_pvalue_slider) {
    
    # Create scatter plot
    scatter_plt <- ggplot(fgsea_results, aes(x = NES, y = -log10(padj), color = padj < 10^(adj_pvalue_slider))) +
      geom_point() +
      scale_color_manual(values = c("grey", "black")) +
      labs(x = "NES", y = "-log10(padj)", color = "Adjusted p-value < threshold") +
      theme_bw()
    
    return(scatter_plt)
  }
  
  
  #fgsea_results_filtered <- reactive({
  #fgsea_filtered_table(run_gsea(load_fgsea(), gmt), input$nes_direction, adj_p_val())
  #})
  
  
  observeEvent(input$submit_fgsea, { 
    output$gsea_results <- DT::renderDataTable({
      run_gsea(load_fgsea(), gmt)
    })
  })
  
  observeEvent(input$submit_fgsea_barplot, { 
    output$fgsea_bar_plot <- renderPlotly ({
      fgsea_barplot(run_gsea(load_fgsea(), gmt), input$num_top_pathways)
    })
  })
  
  output$fgsea_table_output <- DT::renderDataTable({
    fgsea_filtered_table(fgsea_results= run_gsea(load_fgsea(), gmt), 
                         NES_direction=input$nes_direction, 
                         adj_pvalue_slider=adj_p_val()) 
  })
  
  # Define the file download action
  output$download_fgsea_table <- downloadHandler(
    filename = function() {'fgsea_results.csv'},
    content = function(file) {
      write.csv(fgsea_filtered_table(fgsea_results= run_gsea(load_fgsea(), gmt), 
                                     NES_direction=input$nes_direction, 
                                     adj_pvalue_slider=adj_p_val()), 
                file)
    }
  )
  
  output$fgsea_scatter_plot <- renderPlotly ({
    fgsea_scatter_plot_fn(run_gsea(load_fgsea(), gmt), adj_p_val())
  })
  
  ### THIS IS CRAP
  
  # Define reactive expression for plot data
  plot_data <- reactive({
    pvals_hist <- ggplot(labeled_results) +
      geom_histogram(aes(pvalue), fill="lightblue", color="black", bins=50) +
      labs(y="count") +
      ggtitle("Histogram of raw pvalues obtained from DE analysis (vP0 vs. vAd)") +
      theme_minimal()
    return(pvals_hist)
  })
  
  # Define output for plot component
  output$plots_output <- renderPlotly({
    plot_data() %>%
      ggplot(aes(x = Species, y = value, fill = variable)) +
      geom_col(position = "dodge") +
      labs(x = "Species", y = "Mean value") +
      theme_bw() +
      theme(legend.position = "top")
  })
  
  ### CRAP ENDS
}

# Run the app
shinyApp(ui, server)

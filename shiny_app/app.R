# ConTra Facet Z-Score Shiny App
# This app allows users to select a gene and view its regulators faceted by type

library(shiny)
library(ggplot2)
library(dplyr)
library(readr)
library(tidyr)

# Data loading and processing functions
source("data_functions.R")

# Load data on app startup
message("Loading ConTra datasets...")
data_list <- load_contra_data()
gene_choices <- get_gene_choices(data_list$genes)
message("Data loaded successfully!")

# Define UI
ui <- fluidPage(
  titlePanel("ConTra: Gene Regulator Visualization"),
  
  sidebarLayout(
    sidebarPanel(
      h3("Gene Selection"),
      
      selectInput("gene",
                  "Choose a gene:",
                  choices = gene_choices,
                  selected = gene_choices[1]),
      
      hr(),
      
      h4("Plot Options"),
      checkboxInput("show_gene_overlay", 
                    "Show gene overlay in all panels", 
                    value = TRUE),
      
      numericInput("plot_height", 
                   "Plot height (px):", 
                   value = 400, 
                   min = 200, 
                   max = 800),
      
      hr(),
      
      h4("About"),
      p("This app displays facet z-score plots for gene regulators."),
      p("Z-scores are computed per entity across timepoints (T1-T4)."),
      p("Each panel shows regulators of different types: miRNA, lncRNA, methylation.")
    ),
    
    mainPanel(
      h3("Facet Z-Score Plot"),
      
      div(id = "plot-container",
          uiOutput("plot_output")
      ),
      
      hr(),
      
      h4("Expression Data"),
      tableOutput("expression_table"),
      
      hr(),
      
      h4("Plot Information"),
      verbatimTextOutput("plot_info")
    )
  )
)

# Define server logic
server <- function(input, output, session) {
  
  # Reactive expression for gene expression data
  gene_data <- reactive({
    req(input$gene)
    
    # Get expression data for selected gene and its regulators
    get_gene_expression_data(data_list, input$gene)
  })
  
  # Generate facet z-score plot
  output$facet_plot <- renderPlot({
    req(gene_data())
    
    data <- gene_data()
    
    if (nrow(data) == 0) {
      plot.new()
      text(0.5, 0.5, "No regulator data found for this gene", 
           cex = 1.5, col = "red")
      return()
    }
    
    create_facet_zscore_plot(data, input$gene, input$show_gene_overlay)
    
  }, height = function() input$plot_height)
  
  # Dynamic UI for plot output
  output$plot_output <- renderUI({
    plotOutput("facet_plot", height = paste0(input$plot_height, "px"))
  })
  
  # Expression data table
  output$expression_table <- renderTable({
    req(gene_data())
    
    data <- gene_data()
    if (nrow(data) > 0) {
      # Show a summary of the data
      data %>%
        select(entity, type, timepoint, expression) %>%
        head(20)
    }
  })
  
  # Plot information
  output$plot_info <- renderText({
    req(gene_data())
    
    data <- gene_data()
    n_regulators <- length(unique(data$entity[data$type != "gene"]))
    regulator_types <- unique(data$type[data$type != "gene"])
    
    paste(
      "Selected gene:", input$gene, "\n",
      "Number of regulators found:", n_regulators, "\n",
      "Regulator types:", paste(regulator_types, collapse = ", "), "\n",
      "Total data points:", nrow(data)
    )
  })
}

# Run the application
shinyApp(ui = ui, server = server)
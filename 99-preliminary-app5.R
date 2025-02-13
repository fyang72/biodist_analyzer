# Load necessary libraries
library(here)
library(stringr)
library(InteractiveComplexHeatmap)
library(ComplexHeatmap)
library(circlize)
library(DT)
library(shiny)
library(shinydashboard)
library(dplyr)
library(tidyr)
library(readr)
library(readxl)
library(ggplot2)

# Set working directory and load common scripts
home <- getwd()
common_script_file <- here::here("_common.R")
if (file.exists(common_script_file)) {
  source(common_script_file)
} else {
  warning(sprintf("File %s does not exist.", common_script_file))
}

source("./scripts/01_data_table.R")
source("./scripts/02_viral_analysis.R")
source("./scripts/03_vgc_analysis.R")
source("./scripts/04_rna_analysis.R")

# Define global environment (env) for storing data
env <- new.env()

# Function to extract substring before the underscore
extract_prefix <- function(string) {
  if (str_detect(string, "_")) {
    str_sub(string, 1, str_locate(string, "_")[1] - 1)
  } else {
    string
  }
}

# Define brush action function
brush_action <- function(df, input, output, session, heatmap_id, env) {
  row_index <- unique(unlist(df$row_index))
  col_index <- unique(unlist(df$column_index))
  
  if (heatmap_id == "ht_viral") {
    row_selected <- env$viral$row_index[row_index]
    col_selected <- env$viral$col_index[col_index]
    res <- env$viral$df_filtered_pivoted %>%
      as.data.frame() %>%
      slice(row_selected) %>%
      select(all_of(col_selected))
    output <- outputs_after_brush_action_viral(res, input, output, session, env)    
  } else if (heatmap_id == "ht_vgc") {
    row_selected <- env$vgc$row_index[row_index]
    col_selected <- env$vgc$col_index[col_index]
    res <- env$vgc$df_filtered_pivoted %>%
      as.data.frame() %>%
      slice(row_selected) %>%
      select(all_of(col_selected))
    output <- outputs_after_brush_action_vgc(res, input, output, session, env)    
  } else if (heatmap_id == "ht_rna") {
    row_selected <- env$rna$row_index[row_index]
    col_selected <- env$rna$col_index[col_index]
    res <- env$rna$df_filtered_pivoted %>%
      as.data.frame() %>%
      slice(row_selected) %>%
      select(all_of(col_selected))
    output <- outputs_after_brush_action_rna(res, input, output, session, env)
  } else {
    warning("Unknown heatmap_id")
  }
}

# Define UI components
header <- dashboardHeader(
  title = span("Visualization of BioDistData", style = "font-size: 16px;"),
  dropdownMenuOutput("messageMenu"),
  dropdownMenu(type = "notifications", badgeStatus = "warning",
               notificationItem(text = "5 new users joined today", icon = icon("users")),
               notificationItem(text = "Server load at 90%", icon = icon("exclamation-triangle"))
  ),
  dropdownMenu(type = "tasks", badgeStatus = "danger",
               taskItem(value = 80, color = "aqua", "Data analysis"),
               taskItem(value = 50, color = "green", "Report generation")
  )
)

sidebar <- dashboardSidebar(
  sidebarMenu(
    menuItem("Loading", icon = icon("upload"),
             fileInput("file1", "Choose Excel File", accept = ".xlsx")),
    menuItem("Filtering", icon = icon("filter"),
             selectInput("trt_group", label = "Group:", choices = NULL, multiple = TRUE),
             selectInput("paramcd_f", label = "Analyte:", choices = NULL, multiple = FALSE),
             selectInput("tissue_group", label = "Tissue:", choices = NULL, multiple = TRUE),
             selectInput("biological_matrix", label = "Biological Matrix:", choices = NULL, multiple = TRUE),
             selectInput("atpt_f", label = "Timepoint:", choices = NULL, multiple = TRUE),
             menuItem("More", icon = icon("filter"),
                      checkboxGroupInput("sex_f", label = "Sex:", choices = c("M", "F"), inline = TRUE, width = "100%", selected = c("M", "F")),
                      sliderInput("aval_range", "ValueRange:", min = 0, max = 100, value = c(0, 100))
             )
    ),
    menuItem("Facelifting", icon = icon("paint-brush"),
             radioButtons("log10", label = "Log10:", choices = c("Yes" = "Y", "No" = "N"), inline = TRUE, width = "100%", selected = c("Y")),
             radioButtons("add_aval_txt", "Add cell label:", c("Yes" = TRUE, "No" = FALSE), selected = FALSE, inline = TRUE)
    ),
    actionButton("generate_plot", label = "Generate plots", style = "float:left;color: #fff; background-color: #328332; border-color: #328332")
  )
)

body <- dashboardBody(
  navbarPage(
    id = "bio_dist_data_analysis", 
    title = "", 
    theme = "bootstrap.css",   
    tabPanel(title = "Data", id = "panel_data", uiOutput("ui_data_table")),
    tabPanel(title = "Viral", id = "panel_viral", ui_viral_analysis()), 
    tabPanel(title = "VGC", id = "panel_VGC", ui_vgc_analysis()),  
    tabPanel(title = "RNA", id = "panel_RNA", ui_rna_analysis()), 
    tabPanel(title = "Protein", id = "panel_protein", uiOutput("ui_protein_analysis"))
  )
)

# Define UI
ui <- dashboardPage(header, sidebar, body)

# Define server components
server <- function(input, output, session) {
  
  # Reactive expression to load data from excel file
  data_list <- reactive({
    data_list <- loading_from_excel(input$file1)
    env$data_list <- data_list
    data_list
  })
  
  # Reactive expression to process the loaded data
  input_data <- reactive({
    req(input$file1)
    req(input$log10)
    
    df <- data_list() %>%
      pre_process_input_data() %>% 
      mutate(AVAL = suppressWarnings(as.numeric(AVAL))) %>%
      filter(!is.na(AVAL)) %>%
      mutate(USUBJID = paste0("P", USUBJID))  # in case of integer of ID
     
    if (input$log10 == "Y") {
      df <- df %>% mutate(AVAL = ifelse(AVAL <= 0, 0, log10(AVAL)))
    }
    
    df
  })
 
  # Render UI components
  output$ui_data_table <- renderUI({
    ui_data_table(input_data())
  })
  
  output$ui_viral_shedding <- renderUI({})
  
  output$ui_protein_analysis <- renderUI({})
 
  # Update UI components based on the loaded data
  observe({
    data <- input_data()
 
    updateSelectInput(session, "trt_group", choices = unique(data$GROUP_f), selected = unique(data$GROUP_f) %>% sort())
    updateSelectInput(session, "paramcd_f", choices = unique(data$PARAMCD_f), selected = unique(data$PARAMCD_f)[1])
    updateSliderInput(session, "aval_range", 
                      min = min(data$AVAL, na.rm = TRUE) %>% base::signif(2),
                      max = max(data$AVAL, na.rm = TRUE) %>% base::signif(2),
                      value = range(data$AVAL, na.rm = TRUE) %>% base::signif(2))
  })
  
  observe({
    data <- input_data()  
    shiny::req(input$paramcd_f, input$trt_group)
    atpt_f_choices <- data %>% filter(PARAMCD_f %in% input$paramcd_f, GROUP_f %in% input$trt_group) %>% pull(ATPT_f) %>% unique() %>% sort()
    updateSelectInput(session, "atpt_f", choices = atpt_f_choices, selected = atpt_f_choices)
  })
  
  observe({
    data <- input_data()  
    shiny::req(input$paramcd_f, input$trt_group)
    tissue_choices <- data %>% filter(PARAMCD_f %in% input$paramcd_f, GROUP_f %in% input$trt_group) %>% pull(TISSUECD) %>% unique() %>% sort()
    updateSelectInput(session, "tissue_group", choices = tissue_choices, selected = tissue_choices)
  })
  
  observe({
    data <- input_data()  
    shiny::req(input$paramcd_f, input$trt_group)
    matrix_choices <- data %>% filter(PARAMCD_f %in% input$paramcd_f, GROUP_f %in% input$trt_group, TISSUECD %in% input$tissue_group) %>% pull(MATRIXCD) %>% unique() %>% sort()
    updateSelectInput(session, "biological_matrix", choices = matrix_choices, selected = matrix_choices)
  })

  # Generate plot based on the input data and parameters
  generate_plot <- function(input, output, session, env, param_prefix, heatmap_id) {
    df_filtered <- input_data() %>%
      dplyr::filter(
        GROUP_f %in% input$trt_group,
        PARAMCD_f %in% input$paramcd_f,
        ATPT_f %in% input$atpt_f,
        SEX_f %in% input$sex_f,
        AVAL >= input$aval_range[1] & AVAL <= input$aval_range[2]
      )
    
    if (nrow(df_filtered) == 0) print("no data after filtering")  
    if (nrow(df_filtered) == 0) return(NULL)

    na_threshold = 0.5
    if (param_prefix %in% "viral") {
      #df_filtered = data
      #df_filtered = data
      na_threshold = 0
      tmp = df_filtered %>% distinct(MATRIXCD, ATPT_f)  %>% 
      arrange(MATRIXCD, ATPT_f)  %>%
      mutate(MATRIXCD = paste0(MATRIXCD, "-", ATPT_f))

      df_filtered <- df_filtered %>% 
        mutate(MATRIXCD = paste0(MATRIXCD, "-", ATPT_f))  %>% 
        mutate(MATRIXCD = ordered(MATRIXCD, levels=tmp$MATRIXCD)) 

#     df_pivoted <- df_filtered %>%  
#     dplyr::distinct(USUBJID, MATRIXCD, AVAL) %>% 
#     pivot_wider(
#       id_cols = c("USUBJID"),
#       names_from = c("MATRIXCD"),
#       names_sep = "-",  
#       values_from = "AVAL"
#     )

#     df_pivoted = df_pivoted %>% select(USUBJID, all_of(levels(df_filtered$MATRIXCD)))

#  colnames(df_pivoted)

      # print("LINE 217")
      # print( df_filtered %>% pull(ATPT_f) %>% unique())
      # print(input$atpt_f)
    }

    df_filtered_pivoted <- prepare_dataset_for_heatmap(df_filtered, scale_data = FALSE, na_threshold = na_threshold)

    env[[param_prefix]]$tdata_filtered <- df_filtered
    env[[param_prefix]]$df_filtered_pivoted <- df_filtered_pivoted
    env[[param_prefix]]$row_index <- seq_len(nrow(df_filtered_pivoted))
    env[[param_prefix]]$col_index <- seq_len(ncol(df_filtered_pivoted))

    ht <- create_heatmap2(df_filtered_pivoted, df_filtered, col_scheme = NULL, cluster_rows = FALSE, cluster_columns = FALSE, base_font_size = 8, add_aval_txt = input$add_aval_txt)

    if (!is.null(ht)) {
      makeInteractiveComplexHeatmap(input, output, session, ht_list = ht, heatmap_id = heatmap_id, 
                                    brush_action = function(df, input, output, session) {
                                      brush_action(df, input, output, session, heatmap_id = heatmap_id, env)
                                    })
    } else {
      output[[ht_heatmap]]- renderPlot({
        grid.newpage()
        grid.text("No row exists after filtering.")
      })
    }
  }

  observeEvent(input$generate_plot, {
    req(input$file1)

    param_prefix <- tolower(extract_prefix(input$paramcd_f))
    if (param_prefix == "viral") {
      generate_plot(input, output, session, env, param_prefix, paste0("ht_", param_prefix))
    } else if (param_prefix == "vgc") {
      generate_plot(input, output, session, env, param_prefix, paste0("ht_", param_prefix))
    } else if (param_prefix == "rna") {
      generate_plot(input, output, session, env, param_prefix, paste0("ht_", param_prefix))
    } else {
      warning("Unknown parameter prefix")
    }
  }, ignoreNULL = FALSE)

  # Render message menu
  output$messageMenu <- renderMenu({
    messageData <- data.frame(from = c("loading", "creating heatmap"), message = c("error in loading", "error in creating heatmap"))
    msgs <- apply(messageData, 1, function(row) {
      messageItem(from = row[["from"]], message = row[["message"]])
    })
    dropdownMenu(type = "messages", badgeStatus = "success", .list = msgs)
  })
}

# Run the app
shinyApp(ui, server)

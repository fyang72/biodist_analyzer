###########################################################################
# My note
# https://rstudio.github.io/shinydashboard/behavior.html
###########################################################################
rm(list = ls())

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

# Set working directory and load data
home <- getwd()
common_file <- here::here("_common.R")
if (file.exists(common_file)) {
  source(common_file)
} else {
  warning(paste("File", common_file, "does not exist."))
}

# spec4matrix (an internal library)
data_dir <- here::here(home, "data")
#input_file <- "data_template_based_on_20424913 - 20502226 Biodistribution Tables (ID 5879653)_Review_2025-01-18.xlsx"
#input_file <- "data_template_based_on_20424913 - 20502226 Biodistribution Tables (ID 5879653)_Review_2025-02-06.xlsx"
excel_file <- "data_template_0208.xlsx"

spec4matrix <-
  readxl::read_xlsx(path = here::here(data_dir, excel_file),
                    col_names = FALSE,
                    skip = 4,  # the top 4 lines will be always skipped 
                    sheet = "Biological Matrix") %>%
  janitor::row_to_names(row_number = 1) %>%   suppressWarnings() %>%
  dplyr::rename(
    MATRIX = `Biological Matrix`,
    TISSUE = Tissue,
    MATRIXCD = `Biological Matrix Short`,
    TISSUECD = `Tissue Short`
  ) %>%
  dplyr::rename_all(stringr::str_to_upper)



# Define global environment for storing data
env <- new.env()

# Define function to load data
load_data <- function(file) {
  ext <- tools::file_ext(file$datapath)
  
  
  excel_file <- file$datapath
  if(is.null(excel_file)) {
    return(NULL)
  }

  # Get the names of all sheets in the Excel file
  sheet_names <- excel_sheets(path = file$datapath)

  # Read all sheets into a list of data frames
  data_list <- lapply(sheet_names, function(sheet) {
    read_excel(path = file$datapath,
               col_names = FALSE,
               skip = 4,  # the top 4 lines will be always skipped
               sheet = sheet)  
  })

  # Optionally, name the list elements with the sheet names for easy reference
  names(data_list) <- sheet_names
  
  pre_process_input_data(data_list)
  # s
  # 
  # switch(ext,
  #   #csv = read_csv(file$datapath),
  #   xlsx = read_excel(file$datapath),
  #   stop("Invalid file; Please upload a .xlsx file")
  # )
}

# Define function to create heatmap
# This function creates a heatmap based on the filtered data.
# Parameters:
# - df0: Data frame containing the data to be visualized.
# - trt_group: Character vector of treatment group(s) to filter the data.
# - paramcd: Character vector of analyte(s) to filter the data.
# - atpt_f: Character vector of timepoint(s) to filter the data.
# - sex_f: Character vector of sex(es) to filter the data.
# - aval_range: Numeric vector of length 2 indicating the range of values to filter the data.
# - add_aval_txt: Boolean indicating whether to add text labels to the cells.
# Returns:
# - A heatmap object if the filtered data is not empty, otherwise NULL.

make_heatmap <- function(df0, trt_group, paramcd_f, atpt_f, sex_f, aval_range = c(0, 9999), add_aval_txt = FALSE) {
  # Filter the data based on the provided parameters
  df_filtered <- df0 %>%
    dplyr::filter(
      GROUP_f %in% trt_group,
      PARAMCD_f %in% paramcd_f,
      ATPT_f %in% atpt_f,
      SEX_f %in% sex_f,
      AVAL >= aval_range[1] & AVAL <= aval_range[2]
    )
  
  # print("df_filtered at line 114")
  # print(head(df_filtered))
  # print(head(df0))
  # print(trt_group)
  # print(paramcd_f)
  # print(atpt_f)
  # print(sex_f)
  # 

  if (nrow(df_filtered) == 0) return(NULL)

  df_filtered_pivoted <- prepare_dataset_for_heatmap(
    df_filtered,  
    scale_data = FALSE, 
    na_threshold = 0.5 
  )

  # save into global variables
  env$tdata_filtered <- df_filtered
  env$df_filtered_pivoted <- df_filtered_pivoted
  env$row_index <- seq_len(nrow(df_filtered_pivoted))
  env$col_index <- seq_len(ncol(df_filtered_pivoted))

  create_heatmap2(
    df_filtered_pivoted,
    df0,
    col_scheme = NULL,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    base_font_size = 8,
    add_aval_txt = add_aval_txt
  )
}

# Define function to create boxplot
make_boxplot <- function(res) {
  df0 <- env$tdata_filtered %>%
    filter(USUBJID %in% rownames(res), MATRIXCD %in% colnames(res))

  if (nrow(df0) == 0) return(NULL)

  fig = ggplot(df0 %>% filter(!is.na(GROUP_f)), aes(x = MATRIXCD, y = AVAL, color = GROUP_f, fill = GROUP_f)) +
    geom_boxplot() +
    facet_wrap(~ATPT_f, ncol = 2) +
    xlab("Matrix") +
    ylab("VGC (cp/ug)") +
    scale_fill_manual(values = color.scheme.certara) +
    scale_color_manual(values = color.scheme.certara) +
    theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 1, vjust = 1, face = "plain"))

  fig + theme_prism2(base_size = 20)
}

# Define brush action function
brush_action <- function(df, input, output, session) {
  row_index <- unique(unlist(df$row_index))
  col_index <- unique(unlist(df$column_index))

  row_selected <- env$row_index[row_index]
  col_selected <- env$col_index[col_index]

  res <- env$df_filtered_pivoted %>%
    as.data.frame() %>%
    slice(row_selected) %>%
    select(col_selected)

  output[["boxplot"]] <- renderPlot({
    make_boxplot(res)
  }, height = 700)

  output[["resulting_table"]] <- DT::renderDataTable(
    DT::datatable(res,
                  rownames = TRUE,
                  filter = "top",
                  width = 600,
                  options = list(pageLength = 6, lengthChange = TRUE, width = "100%", scrollX = TRUE))
  )

  output[["note"]] <- renderUI({
    note <- spec4matrix %>%
      filter(MATRIXCD %in% (env$tdata_filtered %>% pull(MATRIXCD) %>% unique())) %>%
      arrange(MATRIXCD) %>%
      unite("ABBR", MATRIXCD, MATRIX, sep = " = ") %>%
      pull(ABBR) %>% paste0(collapse = "; ")
    paste0("Note, ", note, collapse = "")
  })
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
      fluidRow(
      column(width = 12, fileInput("file1", "Choose Excel File", accept = ".xlsx"),
                          # c(".csv", ".xlsx"))
              #downloadLink('downloadData', 'Download Example Dataset')  
       downloadButton("downloadData", "Download Template" #, 
       #style = "padding: 15px 4px 5px 4px;margin: 5px 5px 5px 5px; "
       ))
              #a("Example Dataset", href="./data/data_template_0208.xlsx")
      )
    ),
 
    menuItem("Filtering", icon = icon("filter"),
             selectInput("trt_group", label = "Group:", choices = NULL, multiple = TRUE),
             selectInput("paramcd_f", label = "Analyte:", choices = NULL, multiple = FALSE),   # make it false
             selectInput("tissue_group", label = "Tissue:", choices = NULL, multiple = TRUE),
             selectInput("biological_matrix", label = "Biological Matrix:", choices = NULL, multiple = TRUE) ,
             
              menuItem("More", icon = icon("filter"),
                selectInput("atpt_f", label = "Timepoint:", choices = NULL, multiple = TRUE),             
                checkboxGroupInput("sex_f", label = "Sex:", choices = c("Male", "Female"), inline = TRUE, width = "100%", selected = c("Male", "Female")),
                sliderInput("aval_range", "ValueRange:", min = 0, max = 100, value = c(0, 100))
              )
    ),
    # menuItem("Filtering by Matrix", icon = icon("filter"),
    #          selectInput("tissue_group", label = "Tissue:", choices = NULL, multiple = TRUE),
    #          selectInput("biological_matrix", label = "Biological Matrix:", choices = NULL, multiple = TRUE)    
              
    # ),    
    menuItem("Facelifting", icon = icon("paint-brush"),
             radioButtons("add_aval_txt", "Add cell label:", c("Yes" = TRUE, "No" = FALSE), selected = FALSE, inline = TRUE)
    ),
    actionButton("generate_plot", label = "Generate plots", style = "float:left;color: #fff; background-color: #328332; border-color: #328332")
  )
)

body <- dashboardBody(

  navbarPage(
    id="bio_dist_data_analysis", 
    title="", 
    theme = "bootstrap.css",   
    
    tabPanel(title="Data", id="panel_data", uiOutput("ui_data_table")),
    tabPanel(title="Viral", id="panel_viral", uiOutput("ui_viral_shedding")), 
    #abPanel(title="VGC", id="panel_VGC", uiOutput("ui_vgc_analysis")), 
    

  tabPanel(title="VGC", id="panel_VGC",

  fluidRow(
    column(width = 6, box(title = "Primary-heatmap", width = NULL, solidHeader = TRUE, status = "primary", originalHeatmapOutput("ht", width = 600, height = 350, containment = TRUE))),
    column(width = 6, box(title = "Sub-heatmap", width = NULL, solidHeader = TRUE, status = "primary", subHeatmapOutput("ht", title = NULL, width = 600, height = 350, containment = TRUE))),
    column(width = 12, htmlOutput("note"))
  ),
  tabsetPanel(
    tabPanel("Table",
             fluidRow(
               column(width = 12, box(title = "Result Table in Sub-heatmap", width = NULL, solidHeader = TRUE, status = "primary", DT::dataTableOutput("resulting_table", width = "100%", height = "auto", fill = TRUE)))
             )
    ),
    tabPanel("Boxplot",
             fluidRow(
               column(width = 12, box(title = "Boxplot of Analyte(s) vs Timepoint(s) in Sub-heatmap", width = NULL, solidHeader = TRUE, status = "primary", plotOutput("boxplot", width = "100%", height = 700)))
             )
    ),
    tabPanel("Scatterplot",
             fluidRow(
               column(width = 12, box(title = "Scatterplot ...", width = NULL, solidHeader = TRUE, status = "primary"))
             )
    ),
    tabPanel("Barplot",
             fluidRow(
               column(width = 12, box(title = "Barplot ...", width = NULL, solidHeader = TRUE, status = "primary"))
             )
    ),
    tabPanel("Waterfallplot",
             fluidRow(
               column(width = 12, box(title = "Waterfallplot ...", width = NULL, solidHeader = TRUE, status = "primary"))
             )
    )
  )
), 

  tabPanel(title="RNA", id="panel_RNA", uiOutput("ui_rna_analysis")), 
  tabPanel(title="Protein", id="panel_protein", uiOutput("ui_protein_analysis"))

  ))

# Define UI
ui <- dashboardPage(header, sidebar, body)

# Define server
server <- function(input, output, session) {

  source("./scripts/ui_data_table.R")
  



  output$downloadData <- downloadHandler(
    filename = function() {
      # Use the selected dataset as the suggested file name
      paste0("data_template.xlsx")
    },
    content = function(file) {
      # Write the dataset to the `file` that will be downloaded
      #write.csv(data(), file)       #  library(writexl)    {write_xlsx(data, path = file)}
        library(readxl)
    data_dir <- here::here(HOME, "data")
    excel_file <- "data_template_0208.xlsx"
    # Get the names of all sheets in the Excel file
    sheet_names <- excel_sheets(path = here::here(data_dir, excel_file))
    # Read all sheets into a list of data frames
    data_list <- lapply(sheet_names, function(sheet) {
      read_excel(path = here::here(data_dir, excel_file),
                #col_names = FALSE,
                #skip = 4,
                sheet = sheet)  # the top 4 lines will be always skipped
    })

    # Optionally, name the list elements with the sheet names for easy reference
    names(data_list) <- sheet_names

    
      writexl::write_xlsx(data_list, path = file)
    }
  )






  # Reactive expression to load data
  input_data <- reactive({
    req(input$file1)
    
    load_data(input$file1) %>%
      mutate(AVAL = suppressWarnings(as.numeric(AVAL))) %>%
      filter(!is.na(AVAL)) %>%
      mutate(USUBJID = paste0("P", USUBJID))
    
  })
  

  output$ui_data_table <- renderUI({
    ui_data_table(input_data())
  })
  
  output$ui_viral_shedding <- renderUI({
    #print("in renderUI: titlePlusAgenda")
    
  })
  
  output$ui_vgc_analysis <- renderUI({
    #print("in renderUI: titlePlusAgenda")
    #ui_for_vgc(input_data())
    
  })
  
  output$ui_rna_analysis <- renderUI({
    #print("in renderUI: titlePlusAgenda")
    
  })
    
  output$ui_protein_analysis <- renderUI({
    #print("in renderUI: titlePlusAgenda")
    
  })

  observe({
    data <- input_data()

    # Update selectInput choices based on the loaded data
    updateSelectInput(session, "trt_group", choices = unique(data$GROUP_f), selected = unique(data$GROUP_f))
    updateSelectInput(session, "paramcd_f", choices = unique(data$PARAMCD_f), selected = unique(data$PARAMCD_f)[1])
    updateSliderInput(session, "aval_range", min = min(data$AVAL, na.rm = TRUE), max = max(data$AVAL, na.rm = TRUE), value = range(data$AVAL, na.rm = TRUE))
  })
  
  observe({
    data <- input_data()  
    shiny::req(input$paramcd_f, input$trt_group)
    atpt_f_choices <- data %>% filter(PARAMCD_f %in% input$paramcd_f, GROUP_f %in% input$trt_group) %>% pull(ATPT_f)
     
    updateSelectInput(session, "atpt_f", 
                      choices = atpt_f_choices, 
                      selected = atpt_f_choices[1]
    )
  })
  
  observe({
    data <- input_data()  
    shiny::req(input$paramcd_f, input$trt_group)
    
    tissue_choices <- 
      data %>% filter(
        PARAMCD_f %in% input$paramcd_f, 
        GROUP_f %in% input$trt_group 
      ) %>% pull(TISSUECD)
      
    updateSelectInput(session, "tissue_group", 
                      choices = tissue_choices, 
                      selected = tissue_choices   ) 
  })
  
  observe({
    data <- input_data()  
    shiny::req(input$paramcd_f, input$trt_group)
     
    matrix_choices <- 
      data %>% filter(
        PARAMCD_f %in% input$paramcd_f, 
        GROUP_f %in% input$trt_group, 
        TISSUECD %in% input$tissue_group
      ) %>% pull(MATRIXCD)
      
    updateSelectInput(session, "biological_matrix", 
                      choices = matrix_choices, 
                      selected = matrix_choices   )
  })
  
  
  
  
  observeEvent(input$generate_plot, {
    req(input$file1)
    
    ht <- make_heatmap(
      input_data(),
      trt_group = input$trt_group,
      paramcd_f = input$paramcd_f,
      atpt_f = input$atpt_f,
      sex_f = input$sex_f,
      aval_range = input$aval_range,
      add_aval_txt = input$add_aval_txt
    )
    if (!is.null(ht)) {
      makeInteractiveComplexHeatmap(input, output, session, ht, "ht", brush_action = brush_action)
    } else {
      output$ht_heatmap <- renderPlot({
        grid.newpage()
        grid.text("No row exists after filtering.")
      })
    }
  }, ignoreNULL = FALSE)

  # https://rstudio.github.io/shinydashboard/structure.html#structure-overview
  output$messageMenu <- renderMenu({
    # Create a data frame to hold messages for the dropdown menu.
    # This data frame contains two columns: 'from' and 'message'.
    # 'from' indicates the source of the message, and 'message' contains the message text.
    messageData = 
      data.frame(
        from = c("loading", "creating heatmap"), 
        message = c("error in loading", "error in creating heatmap")
      )
    
    # Generate a list of message items from the messageData data frame.
    # Each row in the data frame is converted to a messageItem.
    msgs <- apply(messageData, 1, function(row) {
      messageItem(from = row[["from"]], message = row[["message"]])
    })

    # Create a dropdown menu with the generated message items.
    dropdownMenu(type = "messages", badgeStatus = "success", .list = msgs)
  })
}

# Run the app
shinyApp(ui, server)

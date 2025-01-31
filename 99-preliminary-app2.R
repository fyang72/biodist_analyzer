
rm(list=ls())
# Load necessary libraries
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

# Set working directory and load data
HOME = getwd()
#data_dir <- here::here(HOME, "data/source/ALXN2350-GLP-MKY/")
#status = "DRAFT"  # Flag for labeling figures as draft
source(here::here(HOME, "_common.R"))
 
#-----------------------------------
# spec4matrix (an internal library)
#-----------------------------------
data_dir <- here::here(HOME, "data")
inFile <- "data_template_based_on_20424913 - 20502226 Biodistribution Tables (ID 5879653)_Review_2025-01-18.xlsx"
spec4matrix = read_xlsx(path = here::here(data_dir, inFile), sheet = "Matrix" ) %>% 
  dplyr::rename(
    #MATRIXID = `CRL ID`, 
    MATRIX = `Biological Matrix`, 
    TISSUE = Tissue, 
    MATRIXCD = `Biological Matrix Short`, 
    TISSUECD = `Tissue Short`
  )
colnames(spec4matrix) <- colnames(spec4matrix) %>% stringr::str_to_upper()    # stringr::str_to_title

# Define global environment for storing data
env = new.env()
 
# Define function to load data
load_data <- function(file) {
  ext <- tools::file_ext(file$datapath)
  switch(ext,
         csv = read_csv(file$datapath),
         xlsx = read_excel(file$datapath),
         stop("Invalid file; Please upload a .csv or .xlsx file")
  )
}

# Define function to create heatmap
make_heatmap <- function(df0, trt_group, paramcd, atpt_f, sex_f, aval_range = c(0, 9999), add_aval_txt = FALSE) {
  df_filtered <- df0 %>%
    filter(
      GROUP_f %in% trt_group,
      PARAMCD %in% paramcd,
      ATPT_f %in% atpt_f,
      SEX_f %in% sex_f,
      AVAL >= aval_range[1] & AVAL <= aval_range[2]
    )  

  if (nrow(df_filtered) == 0) return(NULL)

  df_filtered_pivoted = prepare_dataset_for_heatmap(
    df_filtered,  
    scale_data = FALSE, 
    na_threshold = 0.5 
  )

  # save into global variables
  env$tdata_filtered = df_filtered
  env$res = df_filtered_pivoted
  env$row_index = 1:nrow(df_filtered_pivoted)
  env$col_index = 1:ncol(df_filtered_pivoted)

  create_heatmap2(
    df_filtered_pivoted, 
    df0,
    col_scheme = NULL, # colorRamp2(c(0, 5, 10), c("blue", "white", "red")),                           
    cluster_rows = FALSE, 
    cluster_columns = FALSE, 
    base_font_size = 8, 
    add_aval_txt = add_aval_txt
    )


  # meta4_subj <- tdata_filtered %>%
  #   distinct(USUBJID, GROUP, GROUP_f) %>%
  #   arrange(GROUP, USUBJID)

  # tdata_filtered_pivoted <- tdata_filtered %>%
  #   pivot_wider(id_cols = c("USUBJID"),
  #               names_from = c("MATRIXCD"),
  #               names_sep = "-",
  #               values_from = "AVAL") %>%
  #   unite(ID, USUBJID:USUBJID)

  # janitor_cutoff = 0.5
  # rownames_tdata = tdata_filtered_pivoted$ID
  # tdata_filtered_pivoted = tdata_filtered_pivoted %>% select(-ID) %>% as.matrix()
  # rownames(tdata_filtered_pivoted) = rownames_tdata

  # df1 = tdata_filtered_pivoted
  # saved_colnames <- colnames(df1)
  # df1 <- df1 %>% janitor::remove_empty(which = "cols", cutoff = janitor_cutoff)
  # removed_matrix_lst <- setdiff(saved_colnames, colnames(df1))
  # if (length(removed_matrix_lst) > 0) {
  #   print(paste0("Removed columns:", paste0(removed_matrix_lst, collapse = ", ")))
  # }

  # df1 <- df1 %>% base::scale()
  # tdata_filtered_pivoted = df1

  # tdata_filtered_pivoted_transposed = tdata_filtered_pivoted
  # res <- tdata_filtered_pivoted_transposed

  # m = res
  # l = rep(TRUE, times = nrow(m))

  # # Save data to global environment
  # env$tdata_filtered = tdata_filtered
  # env$tdata_filtered_pivoted_transposed = tdata_filtered_pivoted_transposed
  # env$res = res
  # env$row_index = 1:nrow(res)
  # env$col_index = 1:ncol(res)

  # # Annotations
  # analyte_lst <- colnames(tdata_filtered_pivoted)
  # subject_lst <- rownames(tdata_filtered_pivoted)

  # meta_analyte <- separate(data.frame(analyte_lst), analyte_lst, c("PARAMCD", "ATPT"), sep = "-") %>%
  #   mutate(PARAMCAT = ifelse(PARAMCD %in% c("IL6", "IL10", "IFNY", "MCP1", "TNFA"), "Cytokine", "Complement"))

  # meta_subject <- data.frame(USUBJID = subject_lst) %>%
  #   left_join(df0 %>% distinct(USUBJID, SEX_f, GROUP, GROUP_f))

  # ht = Heatmap(m,
  #              name = "Ratio",
  #              show_row_names = TRUE,
  #              show_column_names = TRUE,
  #              column_names_rot = 45,
  #              cluster_rows = FALSE,
  #              cluster_columns = FALSE,
  #              show_row_dend = FALSE,
  #              cluster_row_slices = FALSE,
  #              column_title = "Matrix",
  #              column_title_side = "bottom",
  #              column_title_gp = gpar(fontsize = 13.2),
  #              column_title_rot = 0,
  #              row_split = factor(meta4_subj$GROUP),
  #              row_names_gp = gpar(fontsize = 12),
  #              cell_fun = function(j, i, x, y, width, height, fill) {
  #                if (add_aval_txt == "yes") {
  #                  grid.text(m[i, j], x, y, gp = gpar(fontsize = 12))
  #                }
  #              })

  # #ht = draw(ht, merge_legend = TRUE)
  # ht
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
  row_index = unique(unlist(df$row_index))
  col_index = unique(unlist(df$column_index))

  row_selected = env$row_index[row_index]
  col_selected = env$col_index[col_index]

  res = env$res %>%
    as.data.frame() %>%
    slice(row_selected) %>%
    select(col_selected)

  output[["boxplot"]] = renderPlot({
    make_boxplot(res)
  }, height = 700)

  output[["resulting_table"]] = DT::renderDataTable(
    DT::datatable(res,
                  rownames = TRUE,
                  filter = "top",
                  width = 600,
                  options = list(pageLength = 6, lengthChange = TRUE, width = "100%", scrollX = TRUE))
  )

  output[["note"]] = renderUI({
    note = spec4matrix %>%
      filter(MATRIXCD %in% (env$tdata_filtered %>% pull(MATRIXCD) %>% unique())) %>%
      arrange(MATRIXCD) %>%
      unite("ABBR", MATRIXCD, MATRIX, sep = " = ") %>%
      pull(ABBR) %>% paste0(collapse = "; ")
    paste0("Note, ", note, collapse = "")
  })
}

# Define UI
ui = dashboardPage(

  dashboardHeader(
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
  ),

  dashboardSidebar(
    sidebarMenu( # ), 
      menuItem("Loading", icon = icon("upload"),
               fileInput("file1", "Choose CSV or Excel File",
                         accept = c(".csv", ".xlsx"))
      ),

      menuItem("Filtering", icon = icon("filter"),
               selectInput("trt_group", label = "Group:", choices = NULL, multiple = TRUE),
               selectInput("paramcd", label = "Analyte:", choices = NULL, multiple = TRUE),
               selectInput("atpt_f", label = "Timepoint:", choices = NULL, multiple = TRUE),
               checkboxGroupInput("sex_f", label = "Sex:", choices = c("Male", "Female"), inline = TRUE, width = "100%", selected = c("Male", "Female")),
               sliderInput("aval_range", "ValueRange:", min = 0, max = 100, value = c(0, 100))
      ), 

    menuItem("Facelifting", icon = icon("paint-brush"),
      radioButtons("add_aval_txt", "Add cell label:", c("Yes" = TRUE, "No" = FALSE), selected = FALSE, inline = TRUE)
    ), 

    actionButton("filter", label = "Generate heatmap", style = "float:left;color: #fff; background-color: #328332; border-color: #328332")
  )
  ), 

  dashboardBody(
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
)
)
 

# Define server
server = function(input, output, session) {

    # Reactive expression to load data
  inputData <- reactive({
    req(input$file1)

    load_data(input$file1) %>%  
    mutate(AVAL = as.numeric(AVAL)) %>%
    mutate(#ARM2 = paste0("Grp", GROUP, "(", as.character(GROUP_f), ")"),
         USUBJID = paste0("P", USUBJID))
  })

  observe({
    data <- inputData()    
    # Update selectInput choices based on the loaded data
    updateSelectInput(session, "trt_group", choices = unique(data$GROUP_f), selected = unique(data$GROUP_f))
    updateSelectInput(session, "paramcd", choices = unique(data$PARAMCD), selected = "VGC_MTH1")
    updateSelectInput(session, "atpt_f", choices = unique(data$ATPT_f), selected = "Terminal Necropsy Weeks 26/27")
    updateSliderInput(session, "aval_range", min = min(data$AVAL, na.rm = TRUE), max = max(data$AVAL, na.rm = TRUE), value = range(data$AVAL, na.rm = TRUE))
  })
  
  observeEvent(input$filter, {
    req(input$file1)
    
    ht = make_heatmap(
      inputData(),
      trt_group = input$trt_group,
      paramcd = input$paramcd,
      atpt_f = input$atpt_f,
      sex_f = input$sex_f,
      aval_range = input$aval_range,
      add_aval_txt = input$add_aval_txt
    )
    if (!is.null(ht)) {
      makeInteractiveComplexHeatmap(input, output, session, ht, "ht", brush_action = brush_action)
    } else {
      output$ht_heatmap = renderPlot({
        grid.newpage()
        grid.text("No row exists after filtering.")
      })
    }
  }, ignoreNULL = FALSE)

  # https://rstudio.github.io/shinydashboard/structure.html#structure-overview
  output$messageMenu <- renderMenu({
    # Code to generate each of the messageItems here, in a list. This assumes
    # that messageData is a data frame with two columns, 'from' and 'message'.
    messageData = 
      data.frame(
        from = c("loading", "creating heatmap"), 
        message = c("error in loading", "error in creating heatmap")
      )
    msgs <- apply(messageData, 1, function(row) {
      messageItem(from = row[["from"]], message = row[["message"]])
    })

    # This is equivalent to calling:
    #   dropdownMenu(type="messages", msgs[[1]], msgs[[2]], ...)
    dropdownMenu(type = "messages", badgeStatus = "success", .list = msgs)
  })
}

# Run the app
shinyApp(ui, server)


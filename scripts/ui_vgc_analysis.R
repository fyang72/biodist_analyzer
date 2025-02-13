
ui_vgc_analysis <- function() { 
  
#tabPanel(title="VGC", id="panel_VGC",   
         #tabPanel(title="Heatmap", id="tab_page1",

tagList(    
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

make_heatmap <- function(df0, env, trt_group, paramcd_f, atpt_f, sex_f, aval_range = c(0, 9999), add_aval_txt = FALSE) {
  # Filter the data based on the provided parameters
  df_filtered <- df0 %>%
    dplyr::filter(
      GROUP_f %in% trt_group,
      PARAMCD_f %in% paramcd_f,
      ATPT_f %in% atpt_f,
      SEX_f %in% sex_f,
      AVAL >= aval_range[1] & AVAL <= aval_range[2]
    )
  
  print("df_filtered at line 114")
  print(head(df_filtered))
  print(head(df0))
  print(trt_group)
  print(paramcd_f)
  print(atpt_f)
  print(sex_f)
  
  
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
  
  ht = create_heatmap2(
    df_filtered_pivoted,
    df0,
    col_scheme = NULL,
    cluster_rows = FALSE,
    cluster_columns = FALSE,
    base_font_size = 8,
    add_aval_txt = add_aval_txt
  )
  
  list(env=env, ht=ht)
  
}




# Define function to create boxplot
make_boxplot <- function(df) {

  if (nrow(df) == 0) return(NULL)
  
  fig = ggplot(df %>% filter(!is.na(GROUP_f)), aes(x = MATRIXCD, y = AVAL, color = GROUP_f, fill = GROUP_f)) +
    geom_boxplot() +
    facet_wrap(~ATPT_f, ncol = 2) +
    xlab("Matrix") +
    ylab("VGC (cp/ug)") +
    scale_fill_manual(values = color.scheme.certara) +
    scale_color_manual(values = color.scheme.certara) +
    theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 1, vjust = 1, face = "plain"))
  
  fig + theme_prism2(base_size = 20)
}



outputs_after_brush_action <- function(res, env, input, output, session) {
  
  output[["boxplot"]] <- renderPlot({
    
   df = env$tdata_filtered %>%
      filter(USUBJID %in% rownames(res), MATRIXCD %in% colnames(res))
   
   make_boxplot(df)
    
  }, height = 700)
  
  output[["resulting_table"]] <- DT::renderDataTable(
    DT::datatable(apply(res, c(1, 2), function(x) signif(x, digits=2)),  ########
                  rownames = TRUE,
                  filter = "top",
                  width = 600,
                  options = list(pageLength = 6, lengthChange = TRUE, width = "100%", scrollX = TRUE))
  )
  
  output[["note"]] <- renderUI({
    
    spec4matrix <- env$data_list$`Biological Matrix` %>% 
      janitor::row_to_names(row_number = 1) %>%   suppressWarnings() %>%
      dplyr::rename(
        MATRIX = `Biological Matrix`,
        TISSUE = Tissue,
        MATRIXCD = `Biological Matrix Short`,
        TISSUECD = `Tissue Short`
      ) %>%
      dplyr::rename_all(stringr::str_to_upper)
    
    note <- spec4matrix %>%
      filter(MATRIXCD %in% (env$tdata_filtered %>% pull(MATRIXCD) %>% unique())) %>%
      arrange(MATRIXCD) %>%
      unite("ABBR", MATRIXCD, MATRIX, sep = " = ") %>%
      pull(ABBR) %>% paste0(collapse = "; ")
    paste0("Note, ", note, collapse = "")
  })
  
  return(output)
}




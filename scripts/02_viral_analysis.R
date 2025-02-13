
ui_viral_analysis <- function() { 
  
#tabPanel(title="viral", id="panel_viral",   
         #tabPanel(title="Heatmap", id="tab_page1",

tagList(    
         fluidRow(
           column(width = 6, box(title = "Primary-heatmap", width = NULL, solidHeader = TRUE, status = "primary", originalHeatmapOutput("ht_viral", width = 600, height = 350, containment = TRUE))),
           column(width = 6, box(title = "Sub-heatmap", width = NULL, solidHeader = TRUE, status = "primary", subHeatmapOutput("ht_viral", title = NULL, width = 600, height = 350, containment = TRUE))),
           column(width = 12, htmlOutput("note_viral"))
         ),
         tabsetPanel(
           tabPanel("Table",
                    fluidRow(
                      column(width = 12, box(title = "Result Table in Sub-heatmap", width = NULL, solidHeader = TRUE, status = "primary", DT::dataTableOutput("resulting_table_viral", width = "100%", height = "auto", fill = TRUE)))
                    )
           ),
           tabPanel("Boxplot",
                    fluidRow(
                      column(width = 12, box(title = "Boxplot of Analyte(s) vs Timepoint(s) in Sub-heatmap", width = NULL, solidHeader = TRUE, status = "primary", plotOutput("boxplot_viral", width = "100%", height = 700)))
                    )
           ),
           tabPanel("Time-profile",
                    fluidRow(
                      column(width = 12, box(title = "Time Profile of Viral Shedding in Multiple Tissues ...", width = NULL, solidHeader = TRUE, status = "primary", plotOutput("timeprofile_viral", width = "100%", height = 700)))
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
# 
# make_heatmap_viral <- function(df0, env, trt_group, paramcd_f, atpt_f, sex_f, aval_range = c(0, 9999), add_aval_txt = FALSE) {
#   # Filter the data based on the provided parameters
#   df_filtered <- df0 %>%
#     dplyr::filter(
#       GROUP_f %in% trt_group,
#       PARAMCD_f %in% paramcd_f,
#       ATPT_f %in% atpt_f,
#       SEX_f %in% sex_f,
#       AVAL >= aval_range[1] & AVAL <= aval_range[2]
#     )
#   
#   print("df_filtered at line 114")
#   print(head(df_filtered))
#   print(head(df0))
#   print(trt_group)
#   print(paramcd_f)
#   print(atpt_f)
#   print(sex_f)
#   
#   
#   if (nrow(df_filtered) == 0) return(NULL)
#   
#   df_filtered_pivoted <- prepare_dataset_for_heatmap(
#     df_filtered,  
#     scale_data = FALSE, 
#     na_threshold = 0.5 
#   )
#   
#   # save into global variables
#   env$viral$tdata_filtered <- df_filtered
#   env$viral$df_filtered_pivoted <- df_filtered_pivoted
#   env$viral$row_index <- seq_len(nrow(df_filtered_pivoted))
#   env$viral$col_index <- seq_len(ncol(df_filtered_pivoted))
#   
#   ht = create_heatmap2(
#     df_filtered_pivoted,
#     df0,
#     col_scheme = NULL,
#     cluster_rows = FALSE,
#     cluster_columns = FALSE,
#     base_font_size = 8,
#     add_aval_txt = add_aval_txt
#   )
#   
#   list(env=env, ht=ht)
#   
# }
# 
# 



outputs_after_brush_action_viral <- function(res, input, output, session, env) {
  
  output[["boxplot_viral"]] <- renderPlot({

   df = env$viral$tdata_filtered %>%
      filter(USUBJID %in% rownames(res), MATRIXCD %in% colnames(res))
   
   make_boxplot_viral(df)
    
  }, height = 700)
  
  
  output[["timeprofile_viral"]] <- renderPlot({

   df = env$viral$tdata_filtered %>%
      filter(USUBJID %in% rownames(res), MATRIXCD %in% colnames(res))
   
   make_timeprofile_for_viral_shedding(df)
    
  }, height = 700)


  
  output[["resulting_table_viral"]] <- DT::renderDataTable(
    DT::datatable(apply(res, c(1, 2), function(x) signif(x, digits=2)),  ########
                  rownames = TRUE,
                  filter = "top",
                  width = 600,
                  options = list(pageLength = 6, lengthChange = TRUE, width = "100%", scrollX = TRUE))
  )
  
  output[["note_viral"]] <- renderUI({
    
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
      filter(MATRIXCD %in% (env$viral$tdata_filtered %>% pull(MATRIXCD) %>% unique())) %>%
      arrange(MATRIXCD) %>%
      unite("ABBR", MATRIXCD, MATRIX, sep = " = ") %>%
      pull(ABBR) %>% paste0(collapse = "; ")
    paste0("Note, ", note, collapse = "")
  })
  
  return(output)
}


# Define function to create boxplot
make_boxplot_viral <- function(df) {
  
  if (nrow(df) == 0) return(NULL)
  
  fig = ggplot(df %>% filter(!is.na(GROUP_f)), aes(x = MATRIXCD, y = AVAL, color = GROUP_f, fill = GROUP_f)) +
    geom_boxplot() +
    facet_wrap(~ATPT_f, ncol = 2) +
    xlab("Matrix") +
    ylab("viral (cp/ug)") +
    scale_fill_manual(values = color.scheme.certara) +
    scale_color_manual(values = color.scheme.certara) +
    theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 1, vjust = 1, face = "plain"))
  
  fig + theme_prism2(base_size = 20)
}







make_timeprofile_for_viral_shedding <- function(
  data, 
   color_by = "TISSUE",
    facet_by = "GROUP_f",
   y_scale_log10 = TRUE) {
  
  # The `df` variable is an internal working data frame used within each function.

  df <- data %>%  
    #filter(PARAMCD_f == PARAMCD_filter) %>% 
    #filter(ATPT_f %in% ATPT_f_filter) %>% 
    mutate(AVAL = as.numeric(AVAL)) %>% 
    drop_na(AVAL)  # Remove rows with NA in AVAL
  
    # Determine the range of the x-axis based on the filtered data
    x_limits <- range(df$ATPTN, na.rm = TRUE)

  p <- ggplot(data = df, aes_string(x = "ATPTN", y = "AVAL", color= color_by)) +
    geom_line(aes(group = USUBJID), color = rgb(0.5, 0.5, 0.5), size = 1, alpha = 0.3) + 
    geom_point() + 
    ggprism::theme_prism() +
    labs(
      x = "Time Point (Hrs)", 
      y = "Concentration (cp/µL)"#, 
      #title = paste("Plot for PARAMCD:", PARAMCD_filter)
    )    # Add labels and title
  
  if (y_scale_log10) {
    p <- p + scale_y_log10()  # Set y-axis to log10 scale
  }
  
  p <- p + 
    #facet_wrap(~GROUP_f) + 
    facet_wrap(as.formula(paste("~", facet_by))) #  ncol=1, scales="free_y")
  
    ggplot2::xlim(x_limits)  # Set x-axis limits based on the filtered data
  
  return(p)
}
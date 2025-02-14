#' @title VGC Analysis Shiny App
#' @description This Shiny application provides a visual and tabular analysis of VGC (Viral Genome Copies) data.
#' @details The application includes multiple tabs for different types of visualizations and a table of results.
#'
#' @section UI Components:
#' - Primary-heatmap: Displays the main heatmap.
#' - Sub-heatmap: Displays a secondary heatmap.
#' - Note: Displays a note related to the VGC data.
#' - Table: Shows a table of results filtered by the sub-heatmap.
#' - Boxplot: Displays a boxplot of analytes vs timepoints.
#' - Scatterplot: Placeholder for scatterplot visualization.
#' - Barplot: Placeholder for barplot visualization.
#' - Waterfallplot: Placeholder for waterfall plot visualization.
#'
#' @section Server Outputs:
#' - boxplot_vgc: Renders a boxplot of VGC data.
#' - resulting_table_vgc: Renders a DataTable of the results.
#' - note_vgc: Renders a note based on the filtered VGC data.
#'
#' @param res The result data to be visualized.
#' @param input Shiny input object.
#' @param output Shiny output object.
#' @param session Shiny session object.
#' @param env Environment containing the VGC data and other necessary data.
#'
#' @return The Shiny output object with the rendered components.
#'
#' @import shiny
#' @import DT
#' @import dplyr
#' @import ggplot2
#' @import janitor
#' @import stringr
#'
library(shiny)
library(DT)
library(dplyr)
library(ggplot2)
library(janitor)
library(stringr)

ui_vgc_analysis <- function() {
  tagList(
    fluidRow(
      column(
        width = 6,
        box(
          title = "Primary-heatmap",
          width = NULL,
          solidHeader = TRUE,
          status = "primary",
          originalHeatmapOutput("ht_vgc", width = 600, height = 350, containment = TRUE)
        )
      ),
      column(
        width = 6,
        box(
          title = "Sub-heatmap",
          width = NULL,
          solidHeader = TRUE, status = "primary",
          subHeatmapOutput("ht_vgc", title = NULL, width = 600, height = 350, containment = TRUE)
        )
      ),
      column(width = 12, htmlOutput("note_vgc"))
    ),
    tabsetPanel(
      tabPanel(
        "Table",
        fluidRow(
          column(
            width = 12,
            box(
              title = "Result Table in Sub-heatmap", width = NULL, solidHeader = TRUE, status = "primary",
              DT::dataTableOutput("resulting_table_vgc", width = "100%", height = "auto", fill = TRUE)
            )
          )
        )
      ),
      tabPanel(
        "Boxplot",
        fluidRow(
          column(
            width = 12,
            box(
              title = "Boxplot of Analyte(s) vs Timepoint(s) in Sub-heatmap", width = NULL, solidHeader = TRUE, status = "primary",
              plotOutput("boxplot_vgc", width = "100%", height = 700)
            )
          )
        )
      ),
      tabPanel(
        "Scatterplot",
        fluidRow(
          column(width = 12, box(title = "Scatterplot of Analyte(s) vs Timepoint(s) in Sub-heatmap", width = NULL, solidHeader = TRUE, status = "primary"))
        )
      ),
      tabPanel(
        "Barplot",
        fluidRow(
          column(width = 12, box(title = "Barplot of Analyte(s) vs Timepoint(s) in Sub-heatmap", width = NULL, solidHeader = TRUE, status = "primary"))
        )
      ),
      tabPanel(
        "Waterfallplot",
        fluidRow(
          column(width = 12, box(title = "Waterfallplot of Analyte(s) vs Timepoint(s) in Sub-heatmap", width = NULL, solidHeader = TRUE, status = "primary"))
        )
      )
    )
  )
}


outputs_after_brush_action_vgc <- function(res, input, output, session, env) {
  output[["boxplot_vgc"]] <- renderPlot(
    {
      df <- env$vgc$tdata_filtered %>%
        filter(USUBJID %in% rownames(res), MATRIXCD %in% colnames(res))

      if (nrow(df) == 0) {
        return(NULL)
      }

      fig <- ggplot(df %>% filter(!is.na(GROUP_f)), aes(x = MATRIXCD, y = AVAL, color = GROUP_f, fill = GROUP_f)) +
        geom_boxplot() +
        facet_wrap(~ATPT_f, ncol = 2) +
        xlab("Matrix") +
        ylab("VGC (cp/ug)") +
        scale_fill_manual(values = color.scheme.certara) +
        scale_color_manual(values = color.scheme.certara) +
        theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 1, vjust = 1, face = "plain"))

      fig + theme_prism2(base_size = 20)
    },
    height = 700
  )

  output[["resulting_table_vgc"]] <- DT::renderDataTable(
    DT::datatable(apply(res, c(1, 2), function(x) signif(x, digits = 2)), ########
      rownames = TRUE,
      filter = "top",
      width = 600,
      options = list(pageLength = 6, lengthChange = TRUE, width = "100%", scrollX = TRUE)
    )
  )
  #
  #   output[["note_vgc"]] <- renderUI({
  #     spec4matrix <- env$data_list$`Biological Matrix` %>%
  #       janitor::row_to_names(row_number = 1) %>%
  #       suppressWarnings() %>%
  #       dplyr::rename(
  #         MATRIX = `Biological Matrix`,
  #         TISSUE = Tissue,
  #         MATRIXCD = `Biological Matrix Short`,
  #         TISSUECD = `Tissue Short`
  #       ) %>%
  #       dplyr::rename_all(stringr::str_to_upper)
  #
  #     note <- spec4matrix %>%
  #       filter(MATRIXCD %in% (env$vgc$tdata_filtered %>% pull(MATRIXCD) %>% unique())) %>%
  #       arrange(MATRIXCD) %>%
  #       unite("ABBR", MATRIXCD, MATRIX, sep = " = ") %>%
  #       pull(ABBR) %>%
  #       paste0(collapse = "; ")
  #     paste0("Note, ", note, collapse = "")
  #   })


  # Define a reactive value to store the cached data
  # spec4matrix_cache <- reactiveVal(NULL)

  # # Update the cache when the data changes
  # observe({
  #   spec4matrix_cache(env$data_list$`Biological Matrix` %>%
  #     janitor::row_to_names(row_number = 1) %>%
  #     suppressWarnings() %>%
  #     dplyr::rename(
  #       MATRIX = `Biological Matrix`,
  #       TISSUE = Tissue,
  #       MATRIXCD = `Biological Matrix Short`,
  #       TISSUECD = `Tissue Short`
  #     ) %>%
  #     dplyr::rename_all(stringr::str_to_upper))
  # })

  # Use the cached data in the renderUI function
  output[["note_vgc"]] <- renderUI({
    #spec4matrix <- spec4matrix_cache()

    note <- env$spec4matrix %>%
      filter(MATRIXCD %in% (env$vgc$tdata_filtered %>% pull(MATRIXCD) %>% unique())) %>%
      arrange(MATRIXCD) %>%
      unite("ABBR", MATRIXCD, MATRIX, sep = " = ") %>%
      pull(ABBR) %>%
      paste0(collapse = "; ")
    paste0("Note, ", note, collapse = "")
  })

  return(output)
}

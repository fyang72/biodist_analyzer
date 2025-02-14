ui_rna_analysis <- function() {
  # tabPanel(title="RNA", id="panel_rna",
  # tabPanel(title="Heatmap", id="tab_page1",

  tagList(
    fluidRow(
      column(width = 6, box(title = "Primary-heatmap", width = NULL, solidHeader = TRUE, status = "primary", originalHeatmapOutput("ht_rna", width = 600, height = 350, containment = TRUE))),
      column(width = 6, box(title = "Sub-heatmap", width = NULL, solidHeader = TRUE, status = "primary", subHeatmapOutput("ht_rna", title = NULL, width = 600, height = 350, containment = TRUE))),
      column(width = 12, htmlOutput("note_rna"))
    ),
    tabsetPanel(
      tabPanel(
        "Table",
        fluidRow(
          column(width = 12, box(title = "Result Table in Sub-heatmap", width = NULL, solidHeader = TRUE, status = "primary", DT::dataTableOutput("resulting_table_rna", width = "100%", height = "auto", fill = TRUE)))
        )
      ),
      tabPanel(
        "Boxplot",
        fluidRow(
          column(width = 12, box(title = "Boxplot of Analyte(s) vs Timepoint(s) in Sub-heatmap", width = NULL, solidHeader = TRUE, status = "primary", plotOutput("boxplot_rna", width = "100%", height = 700)))
        )
      ),
      tabPanel(
        "Scatterplot",
        fluidRow(
          column(width = 12, box(title = "Scatterplot ...", width = NULL, solidHeader = TRUE, status = "primary"))
        )
      ),
      tabPanel(
        "Barplot",
        fluidRow(
          column(width = 12, box(title = "Barplot ...", width = NULL, solidHeader = TRUE, status = "primary"))
        )
      ),
      tabPanel(
        "Waterfallplot",
        fluidRow(
          column(width = 12, box(title = "Waterfallplot ...", width = NULL, solidHeader = TRUE, status = "primary"))
        )
      )
    )
  )
}


outputs_after_brush_action_rna <- function(res, input, output, session, env) {
  output[["boxplot_rna"]] <- renderPlot(
    {
      df <- env$rna$tdata_filtered %>%
        filter(USUBJID %in% rownames(res), MATRIXCD %in% colnames(res))

      if (nrow(df) == 0) {
        return(NULL)
      }

      fig <- ggplot(df %>% filter(!is.na(GROUP_f)), aes(x = MATRIXCD, y = AVAL, color = GROUP_f, fill = GROUP_f)) +
        geom_boxplot() +
        facet_wrap(~ATPT_f, ncol = 2) +
        xlab("Matrix") +
        ylab("RNAxxxxx (cp/ug)") +
        scale_fill_manual(values = color.scheme.certara) +
        scale_color_manual(values = color.scheme.certara) +
        theme(axis.text.x = element_text(size = 10, angle = 45, hjust = 1, vjust = 1, face = "plain"))

      fig + theme_prism2(base_size = 20)
    },
    height = 700
  )

  output[["resulting_table_rna"]] <- DT::renderDataTable(
    DT::datatable(apply(res, c(1, 2), function(x) signif(x, digits = 2)), ########
      rownames = TRUE,
      filter = "top",
      width = 600,
      options = list(pageLength = 6, lengthChange = TRUE, width = "100%", scrollX = TRUE)
    )
  )

  output[["note_rna"]] <- renderUI({
    note <- env$spec4matrix %>%
      filter(MATRIXCD %in% (env$rna$tdata_filtered %>% pull(MATRIXCD) %>% unique())) %>%
      arrange(MATRIXCD) %>%
      unite("ABBR", MATRIXCD, MATRIX, sep = " = ") %>%
      pull(ABBR) %>%
      paste0(collapse = "; ")
    paste0("Note, ", note, collapse = "")
  })

  return(output)
}

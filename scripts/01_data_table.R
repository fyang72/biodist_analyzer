

# output$ui_data_table <- renderUI({
  
  ui_data_table <- function(input_data) {
  
  tabsetPanel(
    tabPanel("Data",
             fluidRow(
               column(width = 12, box(title = "Loaded data", width = NULL, solidHeader = TRUE, status = "primary", 
                                      #DT::dataTableOutput("resulting_table", width = "100%", height = "auto", fill = TRUE)
                                      
                                      DT::renderDataTable(
                                        DT::datatable(input_data,
                                                      rownames = TRUE,
                                                      filter = "top",
                                                      width = 600,
                                                      options = list(pageLength = 6, lengthChange = TRUE, width = "100%", scrollX = TRUE))
                                      )
                                      
               ))
             )
    ),
    tabPanel("QC_Data",
             fluidRow(
               column(width = 6, 
                      box(title = "QC_Group", width = NULL, solidHeader = TRUE, status = "primary", 
                          #plotOutput("boxplot", width = "100%", height = 700)
                          DT::renderDataTable(
                            DT::datatable(input_data %>% distinct(GROUPN, GROUP_f, DOSN, DOSU) %>% dplyr::arrange(GROUPN),
                                          rownames = TRUE,
                                          #filter = "top",
                                          width = 600,
                                          options = list(pageLength = 6, lengthChange = TRUE, width = "100%", scrollX = TRUE))
                          )),
                      
                      box(title = "QC_Analyte", width = NULL, solidHeader = TRUE, status = "primary", 
                          #plotOutput("boxplot", width = "100%", height = 700)
                          DT::renderDataTable(
                            DT::datatable(input_data %>% distinct(PARAMCD_f, PARAM_f) %>% dplyr::arrange(PARAMCD_f),
                                          rownames = TRUE,
                                          #filter = "top",
                                          width = 600,
                                          options = list(pageLength = 6, lengthChange = TRUE, width = "100%", scrollX = TRUE))
                          )), 
                      
                      box(title = "QC_Sex", width = NULL, solidHeader = TRUE, status = "primary", 
                          #plotOutput("boxplot", width = "100%", height = 700)
                          DT::renderDataTable(
                            DT::datatable(input_data %>% distinct(SEXN, SEX_f) %>% dplyr::arrange(SEXN),
                                          rownames = TRUE,
                                          #filter = "top",
                                          width = 600,
                                          options = list(pageLength = 6, lengthChange = TRUE, width = "100%", scrollX = TRUE))
                          ))
               ),
               
               column(width = 6, 
                      box(title = "QC_ATPT ...", width = NULL, solidHeader = TRUE, status = "primary", 
                          DT::renderDataTable(
                            DT::datatable(input_data %>% distinct(ATPTN, ATPT_f) %>% dplyr::arrange(ATPTN),
                                          rownames = TRUE,
                                          #filter = "top",
                                          width = 600,
                                          options = list(pageLength = 6, lengthChange = TRUE, width = "100%", scrollX = TRUE))
                          )), 
                      
                      box(title = "QC_Matrix ...", width = NULL, solidHeader = TRUE, status = "primary", 
                          DT::renderDataTable(
                            DT::datatable(input_data %>% distinct(MATRIXCD, MATRIX, TISSUECD, TISSUE),
                                          rownames = TRUE,
                                          filter = "top",
                                          width = 600,
                                          options = list(pageLength = 6, lengthChange = TRUE, width = "100%", scrollX = TRUE))
                          ))
               )
             )
    )#,
    # tabPanel("QC_Analyte",
    #          fluidRow(
    #            column(width = 12, box(title = "Barplot ...", width = NULL, solidHeader = TRUE, status = "primary"))
    #          )
    # ) 
  )
#})

  }
  
  
  
  
  # 
  # 
  # ui_for_vgc <- function(input_data) {
  # 
  # tabPanel(title="Heatmap", id="tab_page1", 
  #          
  #          fluidRow(
  #            column(width = 6, box(title = "Primary-heatmap", width = NULL, solidHeader = TRUE, status = "primary", originalHeatmapOutput("ht", width = 600, height = 350, containment = TRUE))),
  #            column(width = 6, box(title = "Sub-heatmap", width = NULL, solidHeader = TRUE, status = "primary", subHeatmapOutput("ht", title = NULL, width = 600, height = 350, containment = TRUE))),
  #            column(width = 12, htmlOutput("note"))
  #          ),
  #          tabsetPanel(
  #            tabPanel("Table",
  #                     fluidRow(
  #                       column(width = 12, box(title = "Result Table in Sub-heatmap", width = NULL, solidHeader = TRUE, status = "primary", DT::dataTableOutput("resulting_table", width = "100%", height = "auto", fill = TRUE)))
  #                     )
  #            ),
  #            tabPanel("Boxplot",
  #                     fluidRow(
  #                       column(width = 12, box(title = "Boxplot of Analyte(s) vs Timepoint(s) in Sub-heatmap", width = NULL, solidHeader = TRUE, status = "primary", plotOutput("boxplot", width = "100%", height = 700)))
  #                     )
  #            ),
  #            tabPanel("Scatterplot",
  #                     fluidRow(
  #                       column(width = 12, box(title = "Scatterplot ...", width = NULL, solidHeader = TRUE, status = "primary"))
  #                     )
  #            ),
  #            tabPanel("Barplot",
  #                     fluidRow(
  #                       column(width = 12, box(title = "Barplot ...", width = NULL, solidHeader = TRUE, status = "primary"))
  #                     )
  #            ),
  #            tabPanel("Waterfallplot",
  #                     fluidRow(
  #                       column(width = 12, box(title = "Waterfallplot ...", width = NULL, solidHeader = TRUE, status = "primary"))
  #                     )
  #            )
  #          )
  # )
  # 
  # 
  # }
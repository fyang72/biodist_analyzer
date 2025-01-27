
rm(list=ls())

HOME = getwd()
data_dir <- here::here(HOME, "data/source/ALXN2350-GLP-MKY/")
status = "DRAFT"  #flag for labeling figures as draft
source(here::here(HOME, "_common.R"))

source("./scripts/utils.R")
source("./scripts/my_barplot.R")
source("./scripts/my_heatmap.R")
source("./scripts/my_plot_conc_profile.R")

actionButton_style ="float:left;color: #fff; background-color: #328332; border-color: #328332"
 
# ggplot settings, xgx_theme_set()
options(mrggsave.dir = here("deliv/figure"), mrg.script = "app.qmd")
options(dplyr.summarise.inform = TRUE)

####################################################

dat0 <- readr::read_csv(file=here(HOME, "data", "derived", "alxn2350_biodistribution_data.csv"))

spec_location <- here::here(data_dir, "preclinical2.yml")
spec <- ys_load(spec_location) %>% ys_namespace("plot")

# there is no "P0303" and "P0603" in dat0
dat0 <- dat0 %>%
  yspec_add_factors(spec, STUDYID, ARM, DOSE, SEX, ATPT, BLQ, TISSUE) #%>%  # Refactor table for plots

dat0 <- dat0 %>% 
mutate(ARM2 = paste0("Grp", GROUP, "(", as.character(ARM_f), ")"), 
       ARM2 = ordered(ARM2, levels =  paste0("Grp", dat0$GROUP%>% unique(), "(", levels(dat0$ARM_f), ")")   ), 
       USUBJID = SUBJECT
) 
tdata0 <- dat0 


#######################
tdata0 = tdata # save it
#######################

if (1==2) { 
  tdata <- tdata0   %>% 
    pivot_wider(id_cols=c("USUBJID" ),
                #id_expand=FALSE,
                names_from=c("PARAMCD", "ATPT"),
                names_sep = "-", 
                values_from = "AVAL"
    ) %>%
    tidyr::unite(ID, USUBJID:USUBJID)
  
  rownames_tdata = tdata$ID
  tdata = tdata %>% select(-ID)
  tdata = tdata %>% as.matrix()
  rownames(tdata) = rownames_tdata
  
  ht = ComplexHeatmap::Heatmap(
    #t(scale(t(tdata))), 
    #scale(t(tdata)), 
    t(tdata), 
    name = "value", 
    cluster_rows = FALSE, 
    cluster_columns = FALSE
  )
  
  # ht = draw(ht, merge_legend = TRUE)
  # ht
  
  res <- t(tdata)

}


library(InteractiveComplexHeatmap)
library(ComplexHeatmap)
library(circlize)
library(GetoptLong)
library(DT)
library(shiny) 
library(DT)

env = new.env()

make_heatmap = function(tdata0, trt_group, paramcd, atpt, sex, aval_range=c(0, 9999), 
                        add_aval_txt="no") {
  # 
  # l = res$padj <= fdr & res$baseMean >= base_mean & 
  #   abs(res$log2FoldChange) >= log2fc; l[is.na(l)] = FALSE
  #   
  #   if(sum(l) == 0) return(NULL)
  #   
  #   m = counts(dds, normalized = TRUE)
  #   m = m[l, ]
  
  #l = 1:nrow(m)
  
  #m = res #t(tdata)
  tdata_filtered <-  tdata0 %>% 
    dplyr::filter(
      ARM2 %in% trt_group, 
      PARAMCD %in% paramcd, 
      ATPT %in% atpt, 
      SEX_f %in% sex, 
      AVAL >= aval_range[1] &  AVAL <= aval_range[2]
    )    %>% 
    dplyr::filter(
      ATPT  %in% 92, 
      PARAMCD  %in% "VGC"
    ) %>% 
    
    dplyr::mutate(
      AVAL = log(as.numeric(AVAL))
    )
  
  if (nrow(tdata_filtered)==0) return(NULL)
  
  tdata_filtered_pivoted <- tdata_filtered %>% 
    pivot_wider(id_cols=c("USUBJID" ),
                #id_expand=FALSE,
                names_from=c("MATRIX"),  
                names_sep = "-", 
                values_from = "AVAL"
    ) %>%
    tidyr::unite(ID, USUBJID:USUBJID)
  
  
  ###########################
  
  janitor_cutoff = 0.5
  rownames_tdata = tdata_filtered_pivoted$ID
  tdata_filtered_pivoted = tdata_filtered_pivoted %>% select(-ID)
  tdata_filtered_pivoted = tdata_filtered_pivoted %>% as.matrix()
  rownames(tdata_filtered_pivoted) = rownames_tdata
  
  df1 = tdata_filtered_pivoted
  
  saved_colnames <- colnames(df1) 
  df1 <- df1 %>% janitor::remove_empty(which = "cols", cutoff=janitor_cutoff)
  #df2 <- df2 %>% select(colnames(df1)) #%>% select(-SUBJECT, -ARM_f, -DOSE, -SEX)
  removed_matrix_lst <- setdiff(saved_colnames, colnames(df1))
  if (length(removed_matrix_lst)>0) {
    print(paste0("Removed columns:", paste0(setdiff(saved_colnames, colnames(df1)), collapse=", ")))
  }
  
  df1 <-  df1 %>% base::scale()
  
  tdata_filtered_pivoted = df1
  #############################
  
  

  
  tdata_filtered_pivoted_transposed =  (tdata_filtered_pivoted)    #  t(tdata_filtered_pivoted)  
  res <- tdata_filtered_pivoted_transposed
  
  m = res
  l = rep(TRUE, times = nrow(m))
  
  # save them
  env$tdata_filtered = tdata_filtered
  env$tdata_filtered_pivoted_transposed = tdata_filtered_pivoted_transposed
  
  env$res = res   # save it as a global variable
  env$row_index = 1:nrow(res) # which(l)
  env$col_index = 1:ncol(res)
  

  # for annotations
  analyte_lst <- colnames(tdata_filtered_pivoted)
  subject_lst <- rownames(tdata_filtered_pivoted)
  
  library(tidyr)
  meta_analyte <- 
    tidyr::separate(data.frame(analyte_lst), analyte_lst, c("PARAMCD", "ATPT"), sep = "-") %>% 
    mutate(PARAMCAT = ifelse(PARAMCD %in% c("IL6",  "IL10", "IFNY", "MCP1", "TNFA"), "Cytokine", "Complement"))
  
  meta_subject <- 
    data.frame(USUBJID = subject_lst) %>% 
    left_join(tdata0 %>% distinct(USUBJID, SEX_f, GROUP, ARM_f) ) 
  # 
  # subj_lst <- meta %>% pull(SUBJECT) %>% str_extract("(\\d)+") %>% as.integer()
  # 
  # ht <- Heatmap(
  #   df1, 
  #   cluster_rows = cluster_rows,  
  #   cluster_columns = cluster_columns, 
  #   
  #   #row_split =  factor(meta %>% pull(ARM_f), levels=meta$ARM_f %>% levels()), #factor(paste0(dose_lst, dose_unit), levels=c(paste0(dose_lst %>% unique() %>% sort(), dose_unit))),
  #   row_split =  factor(meta %>% pull(GROUP) ),
  #   
  #   row_names_gp = gpar(fontsize = 12), 
  #   column_names_rot = 45,  #	Rotation of column titles.)       
  #   show_row_dend = FALSE, 
  #   cluster_row_slices = FALSE, 
  #   
  #   heatmap_legend_param = list(
  #     title = paste0("Scaled value"), # , which_paramcd, "(", paste0(which_tissue, collapse=" "), ")"),  
  #     #at = c(-4, 0, 4), 
  #     direction = "horizontal", 
  #     legend_width = unit(4, "cm")
  #   ), 
  #   
  #   top_annotation = HeatmapAnnotation(
  #     Type = data.frame(MATRIX = df1%>% colnames()) %>% 
  #       left_join(df0 %>% distinct(MATRIX, ORDER), by="MATRIX") %>% 
  #       pull(ORDER), #t1[df1%>% colnames()],  
  #     which = c("column"), 
  #     col= list(Type = c("Primary" ="#6d405d", "Secondary"  = "#093b6d", "Other"  = "gray")), #  "#6d405d", "#093b6d",
  #     #name= "sfs",  
  #     #annotation_name_side = "top", 
  #     show_annotation_name = TRUE
  #   ), 
  #   
  #   right_annotation = HeatmapAnnotation(
  #     Dose = meta %>% pull(ARM_f),
  #     
  #     #as.name(paste0("Dose(", dose_unit, ")")) = meta %>% pull(DOSE),
  #     Max.histopath = hist_max, # rep(c(0, 1, 2, 3), 3),
  #     ADA = meta %>% pull(NAB), 
  #     Sex = meta %>% pull(SEX_f),  
  #     Subject = subj_lst, # meta %>% pull(SUBJECT) %>% str_extract("(\\d)+") %>% as.integer(), 
  #     
  #     which = c("row"), 
  #     annotation_legend_param  = list(legend_direction = "horizontal", direction = "horizontal"), # legend_direction = c("vertical", "horizontal"),
  #     col= list(Sex = c("Male" = "#2f71fd", "Female"  = "#f98068" ), 
  #               `Dose(xE12 vg/kg)` = colorRamp2(c(6, 20, 50), c("tan","orange",  "brown")), 
  #               Max.histopath = colorRamp2(c( 1, 2, 3), c("gray90", "gray50", "gray10")), 
  #               Subject = colorRamp2(subj_lst,  gg_color_hue(n=length(subj_lst)))
  #               
  #     )  # c("6" = "gray90", "20" = "gray50", "50" = "gray20"))
  #   ), 
  #   
  #   
  #   cell_fun = function(j, i, x, y, width, height, fill) 
  #   { 
  #     #grid.text(df2[i, j], x, y, gp = gpar(fontsize = 12)), 
  #     grid.points(x, y, size=unit(as.integer(df2[i, j]), "char"),  gp=gpar(fontsize=5))     #  https://stackoverflow.com/questions/68425129/what-does-the-fontsize-represent-in-grid-graphics
  #   }
  #   
  # )
  
  
  
  ht = Heatmap(m, #t(scale(t(m))), 
               name = "Ratio", # "Fold-change",
               
               row_title =  "Subject",   # Y-axis
               row_title_side = "left", # c("left", "right"),
               row_title_gp = gpar(fontsize = 13.2),
               #row_title_rot = switch(row_title_side[1], "left" = 90, "right" = 270),
               
               column_title = "Matrix",    # X-axis
               column_title_side = "bottom", # c("top", "bottom"),
               column_title_gp = gpar(fontsize = 13.2),
               column_title_rot = 0,
               
               
               # right_annotation = HeatmapAnnotation(
               #   type = meta_analyte$PARAMCAT, 
               #   which = "row" 
               # ),
               # top_annotation = HeatmapAnnotation(
               #   which = "column",  
               #   Group = meta_subject$GROUP, 
               #   Sex = meta_subject$SEX_f,
               #   
               #   col= list(
               #     #Group = colorRamp2(seq(1, 7, by=1), c("gray90", "gray80","gray70","gray60","gray50", "gray40", "gray30"))
               #     Group = c(
               #       "1" = "gray90", 
               #       "2" = "gray80",
               #       "3" = "gray70",
               #       "4" = "gray60",
               #       "5" = "gray50", 
               #       "6" = "gray40", 
               #       "7" = "gray30" )
               #     
               #     #Sex = c("Male" = "#2f71fd", "Female"  = "#f98068" ),
               #     #`Dose(xE12 vg/kg)` = colorRamp2(c(6, 20, 50), c("tan","orange",  "brown")),
               #     #Max.histopath = colorRamp2(c( 1, 2, 3), c("gray90", "gray50", "gray10")),
               #     #Subject = colorRamp2(subj_lst,  gg_color_hue(n=length(subj_lst)))
               #     
               #   )
               # ),
               # bottom_annotation = HeatmapAnnotation(
               #   which = "column",  
               #   Sex = meta_subject$SEX_f 
               # ),
               
               show_row_names = FALSE, show_column_names = FALSE, #row_km = 2,
               #column_title = paste0(sum(l), " significant genes with FDR < ", fdr),  
               #column_title = paste0(sum(l), " significant genes with FDR < " ),
               
               cluster_rows = FALSE, 
               cluster_columns = FALSE,
               
               show_row_dend = FALSE,  
               
               cell_fun = function(j, i, x, y, width, height, fill)
               {
                 #grid.text(df2[i, j], x, y, gp = gpar(fontsize = 12)),
                 #grid.points(x, y, size=unit(as.integer(m[i, j]), "char"),  gp=gpar(fontsize=5))     #  https://stackoverflow.com/questions/68425129/what-does-the-fontsize-represent-in-grid-graphics
                 
                 if (add_aval_txt=="yes") {grid.text(m[i, j], x, y, gp = gpar(fontsize = 12))  }
               }
 
  
               )  #+ 
  
  # Heatmap(log10(res$baseMean[l]+1), show_row_names = FALSE, width = unit(5, "mm"),
  #         name = "log10(baseMean+1)", show_column_names = FALSE) +
  # Heatmap(res$log2FoldChange[l], show_row_names = FALSE, width = unit(5, "mm"),
  #         name = "log2FoldChange", show_column_names = FALSE,
  #         col = colorRamp2(c(-2, 0, 2), c("green", "white", "red")))
  
  
  ht = draw(ht, merge_legend = TRUE)
  ht
}

# make the MA-plot with some genes highlighted
make_maplot = function(res#, #highlight = NULL, 
                       #env$res[row_selected, col_selected ], trt_group, paramcd, atpt, sex, aval_range=c(0, 9999)
) {
   
  print(res)
  print(env$tdata_filtered )
  
  
  df0 <- #res[highlight,  ] 
    # res %>%  
    # dplyr::filter(
    #   ARM2 %in% trt_group, 
    #   PARAMCD %in% paramcd, 
    #   ATPT %in% atpt, 
    #   SEX_f %in% sex, 
    #   AVAL >= aval_range[1] &  AVAL <= aval_range[2]
    # )  %>% 
 
    
    env$tdata_filtered  %>% 
    filter(paste0(PARAMCD, "-", ATPT) %in% rownames(res))  %>% 
    filter(USUBJID %in% colnames(res))
    
    #mutate(PARAMCD = as.character(PARAMCD))
  
  if (nrow(df0)==0) return(NULL)
  
  fig = ggplot(df0 %>% filter(!is.na(ARM_f)), aes(x = ATPT, y = AVAL, color = ARM_f, fill= ARM_f)) +
    geom_boxplot( ) + 
    facet_wrap(~PARAMCD, ncol=2) + 
    #scale_y_log10() + 
    #annotation_logticks(sides="l") + 
   
    xlab("Timepoint") + 
    ylab("Ratio from baseline") + 
    scale_fill_manual(values = color.scheme.certara) + 
    scale_color_manual(values = color.scheme.certara)  + 
    theme(axis.text.x = element_text(size=10, angle = 45, hjust = 1, vjust=1, face="plain"))  
  
  fig +  theme_prism2(base_size = 20) 
  
}



# 
# # make the volcano plot with some genes highlited
# make_volcano = function(res, highlight = NULL) {
#   col = rep("#00000020", nrow(res))
#   cex = rep(0.5, nrow(res))
#   names(col) = rownames(res)
#   names(cex) = rownames(res)
#   if(!is.null(highlight)) {
#     col[highlight] = "red"
#     cex[highlight] = 1
#   }
#   x = res$log2FoldChange
#   y = -log10(res$padj)
#   col[col == "red" & x < 0] = "darkgreen"
#   par(mar = c(4, 4, 1, 1))
#   
#   suppressWarnings(
#     plot(x, y, col = col, 
#          pch = 16, 
#          cex = cex,
#          xlab = "log2 fold change", ylab = "-log10(FDR)")
#   )
# }

 

library(shiny)
library(shinydashboard)
body = dashboardBody(
 
  fluidRow(
    column(width = 6,
           box(title = "Primary-heatmap", width = NULL, solidHeader = TRUE, status = "primary",
               originalHeatmapOutput("ht", width = 600, height = 350, containment = TRUE)
           )
    ),
    column(width = 6,
           box(title = "Sub-heatmap", width = NULL, solidHeader = TRUE, status = "primary",
               subHeatmapOutput("ht", title = NULL,  width = 600, height = 350, containment = TRUE)
           )
    )
  ), 
  
  fluidRow(
    column(width = 12,
           box(title = "Boxplot of Analyte(s) vs Timepoint(s) in Sub-heatmap", width = NULL, solidHeader = TRUE, status = "primary",
               plotOutput("ma_plot",   
                          width = "100%",
                          height = 700)
           )   
    )
  ), 
  
  # fluidRow(
  # column(width = 4,
  #        id = "column2",
  #        # box(title = "Sub-heatmap", width = NULL, solidHeader = TRUE, status = "primary",
  #        #     subHeatmapOutput("ht", title = NULL, containment = TRUE)
  #        # ),
  #        box(title = "Output", width = NULL, solidHeader = TRUE, status = "primary",
  #            HeatmapInfoOutput("ht", title = NULL)
  #        ),
  #        box(title = "Note", width = NULL, solidHeader = TRUE, status = "primary",
  #            htmlOutput("note")
  #        ),
  # ),
  # column(width = 4,
  #        box(title = "MA-plot", width = NULL, solidHeader = TRUE, status = "primary",
  #            plotOutput("ma_plot")
  #        ),
  #        box(title = "Volcanno plot", width = NULL, solidHeader = TRUE, status = "primary",
  #            plotOutput("volcanno_plot")
  #        ) # ,
  #        
  # ) ),
  
  fluidRow(
    column(width = 12, 
           box(title = "Result Table in Sub-heatmap", width = NULL, solidHeader = TRUE, status = "primary",
               #DTOutput("res_table", width = "100%", height = "auto", fill = TRUE)
               DT::dataTableOutput("res_table", width = "100%", height = "auto", fill = TRUE)
           )
    ) #, 
    
    # control overflow
    # tags$style("
    #         .content-wrapper, .right-side {
    #             overflow-x: auto;
    #         }
    #         .content {
    #             min-width:1500px;
    #         }
    #     ")
  )
)



library(DT)
library(GetoptLong) # for qq() function


brush_action = function(df, input, output, session) {
  
  #print(df)
  #print(str(df))
  
  row_index = unique(unlist(df$row_index))
  col_index = unique(unlist(df$column_index))
  
  row_selected = env$row_index[row_index]
  col_selected = env$col_index[col_index]
  
    #print("selected")
  # print(env$row_index)
  # print(row_index)
  # print(selected)
  # print(head(df %>% as.data.frame()))
  
  #print(row_selected)
  #print(head(env$res[selected,  ]))
  

  
  res = env$res %>% 
    as.data.frame() %>% 
    dplyr::slice(row_selected) %>% 
    dplyr::select(col_selected)
  
  print("env$res[row_selected, col_selected ]")
  print(res)
  
  output[["ma_plot"]] = renderPlot({
    make_maplot(res #, 
                # trt_group =  input$trt_group, 
                # paramcd = input$paramcd, 
                # atpt = input$atpt, 
                # sex = input$sex, 
                # aval_range = input$aval_range
    )
  }, height = 700 )
  
  output[["volcanno_plot"]] = renderPlot({
    make_volcano(env$res, selected)
  })
  
  output[["res_table"]] = DT::renderDataTable(  # renderDT(
    #formatRound(datatable(res[selected, c("baseMean", "log2FoldChange", "padj")], rownames = TRUE), columns = 1:3, digits = 3)
    
    #formatRound(
      DT::datatable(res, 
                    rownames = TRUE,  
                    filter="top", 
                    width = 600,    # Width/Height in pixels (optional, defaults to automatic sizing)
                    options = list(pageLength = 6, #input$pageLength,
                                   lengthChange = TRUE,
                                   width="100%",
                                   scrollX = TRUE
                    )
      )#, 
      #columns = 1:3, digits = 3
   # )
  )
  
  output[["note"]] = renderUI({
    if(!is.null(df)) {
      HTML(qq("<p>Row indices captured in <b>Output</b> only correspond to the matrix of the differential genes. To get the row indices in the original matrix, you need to perform:</p>
<pre>
l = res$padj <= @{input$fdr} & 
    res$baseMean >= @{input$base_mean} & 
    abs(res$log2FoldChange) >= @{input$log2fc}
l[is.na(l)] = FALSE
which(l)[row_index]
</pre>
<p>where <code>res</code> is the complete data frame from DESeq2 analysis and <code>row_index</code> is the <code>row_index</code> column captured from the code in <b>Output</b>.</p>"))
    }
  })
}

 

ui = dashboardPage(
  dashboardHeader(title = "Heatmap Visualization"),
  dashboardSidebar(
    #selectInput("fdr", label = "Cutoff for FDRs:", c("0.001" = 0.001, "0.01" = 0.01, "0.05" = 0.05)),
    #numericInput("base_mean", label = "Minimal base mean:", value = 0),
    #numericInput("log2fc", label = "Minimal abs(log2 fold change):", value = 0),
    
    selectInput("trt_group", 
                label = "Group:",
                choices=tdata0$ARM2 %>% unique() %>% sort(), 
                selected = tdata0$ARM2 %>% unique() %>% sort(), 
                multiple = TRUE 
    ),
    selectInput("paramcd", 
                label = "Analyte:", 
                choices=tdata0$PARAMCD %>% unique() %>% sort(), 
                selected = tdata0$PARAMCD %>% unique() %>% sort(), 
                multiple = TRUE 
    ),
    
    selectInput("atpt", 
                label = "Timepoint:", 
                choices=tdata0$ATPT %>% unique() %>% sort(), 
                selected = tdata0$ATPT %>% unique() %>% sort(), 
                multiple = TRUE 
    ), 
    
    checkboxGroupInput("sex", 
                       label=paste0("Sex:"), 
                       choices=c("Male", "Female"), 
                       inline=TRUE, 
                       width="100%",
                       selected=c("Male", "Female")), 
    
    sliderInput("aval_range", "ValueRange:",
                min = range(tdata0$AVAL, na.rm=TRUE)[1], 
                max = range(tdata0$AVAL, na.rm=TRUE)[2], 
                value = range(tdata0$AVAL, na.rm=TRUE) # c(0,100)
    ),
    
    radioButtons("add_aval_txt", "Add cell label:", 
                 c("Yes" = "yes",
                   "No" = "no" 
                 ), 
                 selected = "no",
                 inline = TRUE
    ), 
    
    
    actionButton("filter", label = "Generate heatmap", style=actionButton_style)
  ),
  body
)
 

server = function(input, output, session) {
  observeEvent(input$filter, {
    ht = make_heatmap(
      tdata0,
      trt_group =  input$trt_group, 
      paramcd = input$paramcd, 
      atpt = input$atpt, 
      sex = input$sex, 
      aval_range = input$aval_range, 
      add_aval_txt = input$add_aval_txt
    )
    if(!is.null(ht)) {
      makeInteractiveComplexHeatmap(input, output, session, ht, "ht",
                                    brush_action = brush_action)
    } else {
      # The ID for the heatmap plot is encoded as @{heatmap_id}_heatmap, thus, it is ht_heatmap here.
      output$ht_heatmap = renderPlot({
        grid.newpage()
        grid.text("No row exists after filtering.")
      })
    }
  }, ignoreNULL = FALSE)
}


shinyApp(ui, server)







if (1==2) { 
  
  res[selected,  ]
  df0 <- tdata0 %>% 
    filter(PARAMCD %in% c("IP10", "MCP1"))  %>% 
    mutate(PARAMCD = as.character(PARAMCD))
  
  ggplot(df0 %>% filter(!is.na(ARM_f)), aes(x = ATPT, y = AVAL, color = ARM_f, fill= ARM_f)) +
    geom_boxplot( ) + 
    facet_wrap(~PARAMCD) + 
    #scale_y_log10() + 
    #annotation_logticks(sides="l") + 
    theme_prism2( base_size = 20) + 
    xlab("Timepoint") + 
    ylab("Ratio from baseline") + 
    scale_fill_manual(values = color.scheme.certara) + 
    scale_color_manual(values = color.scheme.certara)  + 
    theme(axis.text.x = element_text(size=10, angle = 45, hjust = 1, vjust=1, face="plain"))  
  
  
}

#scale_fill_manual(values = color.scheme.certara[3:length(color.scheme.certara)]) + # gg_color_hue(n=12)) + 
 

# https://www2.cs.sfu.ca/~cameron/Teaching/D-Lib/RegExp.html

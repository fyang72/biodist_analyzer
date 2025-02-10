 
 
 # Have a dataset (dat0) with columns of "STUDYID USUBJID GROUP   SEX MATRIX  ATPT SAMPID PARAMCD AVAL   LLOQ MATRIXCD TISSUECD  DOSE GROUP_f  SEX_f  ATPT_f". 
 # If a column variable ending with "_f", then it means it is a factor. 
 # 
 # Please create a function of heatmap using the following scripts as the starting point, with the features of 
 # 1) using ComplexHeatmap::Heatmap function
 # 2) be able to split the row (subject) using its associated GROUP; and add annotation to each subject with SEX, GROUP, etc.
 # 3) be able to split the column (MATRIXCD) using its associated TISSUE type;  and add annotation to each MATRIXCD with TISSUE. 
 # 4) be able to scale(t(tdata)) or not
 # 5) be able to cluster the row and columns
 #  

   

  #' Create a custom heatmap
  #'
  #' This function generates a heatmap using the ComplexHeatmap package, with options for filtering, scaling, and annotating the data.
  #'
  #' @param data A data frame containing the data to be visualized. The data frame must include the following columns:
  #'   - `ATPT` (numeric or character): Time point for the x-axis.
  #'   - `ATPT_f` (factor or character): Time point factor to filter the data.
  #'   - `AVAL` (numeric): The value to be plotted.
  #'   - `GROUP_f` (factor or character): Group factor for faceting the plot.
  #'   - `MATRIXCD` (factor or character): Biological matrix in biodistribution data
  #'   - `PARAMCD` (character): Parameter code to filter the data. 
  #'   - `USUBJID` (character): Unique subject identifier.
  #' @param ATPT_f_filter A character vector specifying the filter for the ATPT_f column. Default is "Terminal Necropsy Weeks 26/27".
  #' @param PARAMCD_filter A character vector specifying the filter for the PARAMCD column. Default is "VGC_MTH1".
  #' @param scale_data A logical value indicating whether to scale the data. Default is FALSE.
  #' @param col_scheme A color mapping function created by colorRamp2 for the heatmap. Default is colorRamp2(c(0, 5, 10), c("blue", "white", "red")).
  #' @param na_threshold A numeric value specifying the threshold for removing rows and columns with too many NAs. Default is 0.5.
  #' @param cluster_rows A logical value indicating whether to cluster rows. Default is FALSE.
  #' @param cluster_columns A logical value indicating whether to cluster columns. Default is FALSE.
  #'
  #' @return A heatmap object created by the ComplexHeatmap package.
  #'
  #' @examples
  #' \dontrun{
  #' create_heatmap(data, ATPT_f_filter = "Terminal Necropsy Weeks 26/27", PARAMCD_filter = "VGC_MTH1")
  #' }
  #'
  #' @import ComplexHeatmap
  #' @import dplyr
  #' @import tidyr
  #' @import RColorBrewer
  #' @import grDevices
  #' @export
   



# Define the custom heatmap function
prepare_dataset_for_heatmap <- function(
  data,  
  scale_data = FALSE, 
  na_threshold = 0.5 
) { 
  # Prepare the data for heatmap
  #----------------------------------------
  df <- data %>%  
    dplyr::distinct(USUBJID, MATRIXCD, AVAL) %>% 
    dplyr::mutate(
      AVAL = log10(as.numeric(AVAL))     # Log-transform AVAL
    ) %>%
    pivot_wider(
      id_cols = c("USUBJID"),
      names_from = c("MATRIXCD"),
      names_sep = "-",  
      values_from = "AVAL"
    ) %>%
    tidyr::unite(ID, USUBJID:USUBJID)  # Create rownames for heatmap
 
  # Remove rows with too many NAs (e.g., >50% missing)
  df <- df[rowMeans(is.na(df)) <= na_threshold, ]
  
  # Remove columns with too many NAs (e.g., >50% missing)
  df <- df[, colMeans(is.na(df)) <= na_threshold]
  
  # Extract rownames and convert to matrix
  rownames_df <- df$ID
  colnames_df <- setdiff(colnames(df), "ID")
  # 
  df <- df %>% select(-ID) %>% as.matrix()
  rownames(df) <- rownames_df
  
  # Optionally scale the data
  if (scale_data) {
    df <- t(scale(t(df)))
  }

  return(df)
}


# Define the custom heatmap function
create_heatmap2 <- function(df, data,
                           col_scheme = colorRamp2(c(0, 5, 10), c("blue", "white", "red")),                           
                           cluster_rows = FALSE, 
                           cluster_columns = FALSE, 
                           base_font_size = 6, 
                           add_aval_txt = FALSE
    ) {
  

library(ComplexHeatmap)
library(dplyr)
library(tidyr)
library(RColorBrewer)
library(grDevices) # For colorRampPalette()

  # The `df` variable is an internal working data frame used within each function.
 
  
  # Prepare the annotations
  #----------------------------------------
  
  # Prepare row annotations (for subjects)
  row_annotations <- data %>% 
    dplyr::select(USUBJID, GROUPN, SEX_f) %>%  
    dplyr::distinct() %>%
    dplyr::filter(USUBJID %in% rownames(df)) %>%
    dplyr::mutate(USUBJID0 = USUBJID) %>% 
    tidyr::separate(USUBJID, into=c("USUBJID", "PARAMCD"), sep = "-")  %>% 
    column_to_rownames("USUBJID0")
  
  row_anno <- rowAnnotation(
      GROUPN = row_annotations$GROUPN,
      SEX = row_annotations$SEX_f,
      col = list(
        GROUPN = structure(brewer.pal(n = length(unique(row_annotations$GROUPN)), "Set3"), 
                          names = unique(row_annotations$GROUPN)),
        SEX = structure(c("blue", "pink"), names = c("Male", "Female"))  # unique(row_annotations$SEX_f) %>% levels()) # Example for Male/Female
      ), 
      annotation_name_gp = gpar(fontsize = base_font_size),  # Adjust annotation text size
      annotation_legend_param = list(
      GROUPN = list(
        title_gp = gpar(fontsize = base_font_size),  # Adjust legend title text size
        labels_gp = gpar(fontsize = base_font_size)   # Adjust legend labels text size
      ), 
      SEX = list(
        title_gp = gpar(fontsize = base_font_size),  # Adjust legend title text size
        labels_gp = gpar(fontsize = base_font_size)   # Adjust legend labels text size
      )
      
  )
    )
  
  # Prepare column annotations (for MATRIXCD)
  column_annotations <- data %>%
    dplyr::select(MATRIXCD, TISSUECD) %>%
    dplyr::distinct() %>%
    dplyr::filter(MATRIXCD %in% colnames(df)) %>%
    dplyr::mutate(
      MATRIXCD = ordered(MATRIXCD, levels=colnames(df))
    ) %>% 
    dplyr::arrange(MATRIXCD) %>% 
    column_to_rownames("MATRIXCD")
  
# Generate a dynamic color palette for TISSUECD
  unique_tissues <- unique(column_annotations$TISSUECD)
  tissue_colors <- structure( 
    colorRampPalette(brewer.pal(9, "Set1"))(length(unique_tissues)), # Dynamically generate colors
    names = unique_tissues
  )
  
  col_anno <- HeatmapAnnotation(
    TISSUE = column_annotations$TISSUECD,
    col = list(
      TISSUE = tissue_colors  # Map TISSUECD to colors
    ), 
    annotation_name_gp = gpar(fontsize = base_font_size),  # Adjust annotation text size
    annotation_legend_param = list(
    TISSUE = list(
      title_gp = gpar(fontsize = 5),  # Adjust legend title text size
      labels_gp = gpar(fontsize = base_font_size)   # Adjust legend labels text size
    )
  )
  ) 
  
  # Create the heatmap
  # ---------------------
  ht <- ComplexHeatmap::Heatmap(
    df,
    col = col_scheme, 
    column_names_rot = 45,       # Rotate column names
    name = "Value",              # Legend title
    cluster_rows = cluster_rows, # Cluster rows
    cluster_columns = cluster_columns, # Cluster columns
    
    top_annotation = col_anno,   # Add column annotations
    left_annotation = row_anno,  # Add row annotations 
    row_split = row_annotations %>%  dplyr::select(GROUPN),  # Split rows by GROUP
    column_split = column_annotations$TISSUECD, # Split columns by TISSUE
    #cluster_column_slices = FALSE,
    #column_order = colnames(df) #1:ncol(df)
   
    #width = unit(15, "cm"),  # Adjust width
    #height = unit(9, "cm"),  # Adjust height
    row_names_gp = gpar(fontsize = base_font_size),  # Adjust row names text size
    column_names_gp = gpar(fontsize = base_font_size),  # Adjust column names text size

    row_title_gp = gpar(fontsize = base_font_size),  # Adjust row split label text size
    column_title_gp = gpar(fontsize = base_font_size),  # Adjust column split label text size
    heatmap_legend_param = list(
      title = "Value",
      title_gp = gpar(fontsize = base_font_size),
      labels_gp = gpar(fontsize = base_font_size)
    ), 

   row_gap = unit(.25, "mm"),  # Adjust the size of splitting lines
   column_gap = unit(0.25, "mm"), 


  cell_fun = function(j, i, x, y, width, height, fill) {
    if (add_aval_txt) {
    grid.text(df[i, j] %>% base::signif(digits=2), x, y, gp = gpar(fontsize = base_font_size-2))
    }
  }
  )
  
  # Draw the heatmap
  #draw(ht, padding = unit(c(10, 15, 10, 10), "mm"), merge_legend = TRUE) # Adjust padding
  
  return(ht)
  
  # return(list(ht = ht, 
  #             colnames_lst = colnames_df[column_order(ht) %>% unlist() %>% as.integer()], 
  #             rownames_lst = rownames_df[row_order(ht) %>% unlist() %>% as.integer()]
  # )
  # )
}

# 
# 
# # https://github.com/jokergoo/ComplexHeatmap/issues/299
# top_annotation = HeatmapAnnotation(
#             	t1 = anno_block(gp = gpar(fill = 2:5), 
#             		            labels = c("nm1", "nm2", "nm3", "nm4"))),
#             left_annotation = rowAnnotation(
#             	t2 = anno_block(gp = gpar(fill = 2:5), 
#             		            labels = c("nm1", "nm2", "nm3", "nm4"))),
# 







my_heatmap_org <- function(df0, which_tissue, which_paramcd, which_timepoint,  
                       cluster_rows = FALSE,  
                       cluster_columns = TRUE, 
                       log_all_scale = FALSE, 
                       column_wise_scale = TRUE, 
                       row_wise_scale = FALSE,
                       janitor_cutoff = 0.6, 
                       dose_unit = "xE12 vg/kg"
)  {
  #BiocManager::install("ComplexHeatmap")
  
  # https://jokergoo.github.io/ComplexHeatmap-reference/book/a-single-heatmap.html
  # https://jokergoo.github.io/ComplexHeatmap-reference/book/legends.html
  # https://stackoverflow.com/questions/52259684/r-complexheatmap-relocate-and-resize-legend-from-grid-text   # add grid_text
  #ht_global_opt(heatmap_legend_title_gp = gpar(fontsize = 5, fontface = "bold"), heatmap_legend_labels_gp = gpar(fontsize = 5), heatmap_column_names_gp = gpar(fontsize = 5))
 
  
  library(ComplexHeatmap)
  library(janitor) 
  
  
  df <- df0 %>% 
    filter(#TISSUE %in% which_tissue, 
      ATPT %in% c(which_timepoint), 
      PARAMCD %in% c(which_paramcd )
    )  %>% 
    mutate(AVAL = as.numeric(AVAL), 
           HISTOPTH = as.numeric(HISTOPTH)#, 
           #HISTOPTH = ifelse(HISTOPTH %in% c(0), "", as.character(HISTOPTH)) #, 
           #SUBJECT = paste0(SUBJECT, "(", base::substr(SEX, 1, 1), ")")
    )  
  
  if (log_all_scale) {df <- df %>% mutate(AVAL = log(AVAL))  }
  
  # only a single timepoint
  print(paste0(df$ATPT %>% unique(), collapse=","))
  
  #stopifnot(length(df$ATPT %>% unique())==1)
  
  
  #df1 <- df %>% pivot_wider(id_cols=c("SUBJECT", "ARM_f", "DOSE", "SEX"), names_from=MATRIX, values_from = AVAL)  # , "HISTOPTH"
  #df2 <- df %>% pivot_wider(id_cols=c("SUBJECT", "ARM_f", "DOSE", "SEX"), names_from=MATRIX, values_from = HISTOPTH)  # , "HISTOPTH"
  df <- df %>% select(SUBJECT, MATRIX, PARAMCD, AVAL, HISTOPTH)  
  df1 <- df %>% pivot_wider(id_cols=c("SUBJECT"), names_from=MATRIX, values_from = AVAL) %>% as.data.frame() # , "HISTOPTH"
  df2 <- df %>% pivot_wider(id_cols=c("SUBJECT"), names_from=MATRIX, values_from = HISTOPTH) %>% as.data.frame()  # , "HISTOPTH"
  #df2[is.na(df2)] = "x"
  
  saved_colnames <- colnames(df1) 
  df1 <- df1 %>% janitor::remove_empty(which = "cols", cutoff=janitor_cutoff)
  df2 <- df2 %>% select(colnames(df1)) #%>% select(-SUBJECT, -ARM_f, -DOSE, -SEX)
  removed_matrix_lst <- setdiff(saved_colnames, colnames(df1))
  if (length(removed_matrix_lst)>0) {
    print(paste0("Removed columns:", paste0(setdiff(saved_colnames, colnames(df1)), collapse=", ")))
  }
  
  meta <- df1 %>% distinct(SUBJECT) %>% 
    left_join(df0 %>% distinct(SUBJECT, ARM, ARM_f, DOSN, GROUPN, SEX_f, NAB), by="SUBJECT")
  
  rownames(df1) <- df1 %>% pull(SUBJECT)  # subj_lst
  df1 <- df1 %>% select(-SUBJECT)
  
  rownames(df2) <- df2 %>% pull(SUBJECT)  # subj_lst
  df2 <- df2 %>% select(-SUBJECT)  
  
  # subject level information 
  # subj_lst <-  df1 %>% pull(SUBJECT)  # paste0(df$SUBJECT, "(", DOSE, ")")
  # dose_lst <- df1 %>% pull(DOSE)
  # arm_f_lst <- df1 %>% pull(ARM_f)
  # sex_lst <- df1 %>% pull(SEX)
  
  
  
  #meta <- df1 %>% select(SUBJECT, ARM_f, DOSE, SEX)
  # df1 <- df1 %>% #select(-SUBJECT, -ARM_f, -DOSE, -SEX) %>% 
  #   as.data.frame(() 
  
  
  
  # max histogram
  hist_max <-  
    #filter(ATPT %in% c(which_timepoint), 
    # PARAMCD %in% c(which_paramcd ))  %>% 
    #mutate(SUBJECT = SUBJECT %>% str_extract("(\\d)+") %>% as.integer()) %>%  
    data.frame(SUBJECT=rownames(df1)) %>% 
    left_join(df %>% group_by(SUBJECT) %>% 
                summarise(Max_HIST = max(as.numeric(HISTOPTH), na.rm=TRUE)), 
              by="SUBJECT") %>% 
    pull(Max_HIST)
  
  set.seed(123)
  library(circlize)   # for color selection
  #tt  = df1 %>% scale() %>% as.data.frame()
  
  scale_constant = NULL
  if (column_wise_scale | row_wise_scale)  { 
    scaling_direction <- function(df1, column_wise_scale=TRUE, row_wise_scale=FALSE) {
      if (column_wise_scale)  {df1 <-  df1 %>% base::scale()}
      if (row_wise_scale) {df1 <- t(t(df1) %>% base::scale())}
      df1
    }
    
    # tt = bind_rows(
    #   # df1[which(str_detect(rownames(df1), "Whole")) & str_detect(rownames(df1), "vgc"), ] %>% base::scale() %>% as.data.frame(), 
    #   # df1[which(str_detect(rownames(df1), "Whole")) & str_detect(rownames(df1), "mRNA"), ] %>% base::scale() %>% as.data.frame(), 
    #   # df1[which(str_detect(rownames(df1), "Brain")) & str_detect(rownames(df1), "LPA"), ] %>% base::scale() %>% as.data.frame(), 
    #   # df1[which(str_detect(rownames(df1), "Brain")) & str_detect(rownames(df1), "EVV"), ] %>% base::scale() %>% as.data.frame(), 
    #   # df1[which(str_detect(rownames(df1), "Brain")) & str_detect(rownames(df1), "delta4"), ] %>% base::scale() %>% as.data.frame() , 
    #   #  df1[which(str_detect(rownames(df1), "Whole|Brain", negate=TRUE)) & str_detect(rownames(df1), "LPA"), ] %>% base::scale() %>% as.data.frame()  , 
    #   #  df1[which(str_detect(rownames(df1), "Whole|Brain", negate=TRUE)) & str_detect(rownames(df1), "EVV"), ] %>% base::scale() %>% as.data.frame()  , 
    #   #  df1[which(str_detect(rownames(df1), "Whole|Brain", negate=TRUE)) & str_detect(rownames(df1), "delta4"), ] %>% base::scale() %>% as.data.frame()  
    # 
    #   # df1[which(str_detect(rownames(df1), "vgc")), ] %>% scaling_direction(column_wise_scale, row_wise_scale)  %>% as.data.frame(), 
    #   # df1[which(str_detect(rownames(df1), "mRNA")), ] %>% scaling_direction(column_wise_scale, row_wise_scale)  %>% as.data.frame(), 
    #   # df1[which(str_detect(rownames(df1), "LPA")), ] %>% scaling_direction(column_wise_scale, row_wise_scale)  %>% as.data.frame(), 
    #   # df1[which(str_detect(rownames(df1), "EVV")), ] %>% scaling_direction(column_wise_scale, row_wise_scale)  %>% as.data.frame(), 
    #   # df1[which(str_detect(rownames(df1), "delta4")), ] %>% scaling_direction(column_wise_scale, row_wise_scale)  %>% as.data.frame()   
    #   
    #   df1[which(str_detect(rownames(df1), "vgc")), ] %>% base::scale() %>% as.data.frame(), 
    #   df1[which(str_detect(rownames(df1), "mRNA")), ] %>% base::scale() %>% as.data.frame(), 
    #   df1[which(str_detect(rownames(df1), "LPA")), ] %>% base::scale() %>% as.data.frame(), 
    #   df1[which(str_detect(rownames(df1), "EVV")), ] %>% base::scale() %>% as.data.frame(), 
    #   df1[which(str_detect(rownames(df1), "delta4")), ] %>% base::scale() %>% as.data.frame()   
    #   
    # )
    
    tt=NULL
    for (iparam in df0$PARAMCD %>% as.character() %>% unique())  {
      out1 = df1[which(str_detect(rownames(df1), iparam)), ] %>% base::scale()  
      tt = bind_rows(tt, out1 %>% as.data.frame())
      
      out2 = data.frame(PARAMCD=iparam, 
                        MATRIX = colnames(out1), 
                        center = attr(out1,"scaled:center"), 
                        scale = attr(out1,"scaled:scale"))
      rownames(out2) = NULL
      scale_constant = bind_rows(scale_constant, out2)
    }
    
    df1 <- tt[rownames(df1), ]
  }
  
  
  
  
  
  subj_lst <- meta %>% pull(SUBJECT) %>% str_extract("(\\d)+") %>% as.integer()
  
  ht <- Heatmap(
    df1, 
    cluster_rows = cluster_rows,  
    cluster_columns = cluster_columns, 
    
    #row_split =  factor(meta %>% pull(ARM_f), levels=meta$ARM_f %>% levels()), #factor(paste0(dose_lst, dose_unit), levels=c(paste0(dose_lst %>% unique() %>% sort(), dose_unit))),
    row_split =  factor(meta %>% pull(GROUP) ),
    
    row_names_gp = gpar(fontsize = 12), 
    column_names_rot = 45,  #	Rotation of column titles.)       
    show_row_dend = FALSE, 
    cluster_row_slices = FALSE, 
    
    heatmap_legend_param = list(
      title = paste0("Scaled value"), # , which_paramcd, "(", paste0(which_tissue, collapse=" "), ")"),  
      #at = c(-4, 0, 4), 
      direction = "horizontal", 
      legend_width = unit(4, "cm")
    ), 
    
    top_annotation = HeatmapAnnotation(
      Type = data.frame(MATRIX = df1%>% colnames()) %>% 
        left_join(df0 %>% distinct(MATRIX, ORDER), by="MATRIX") %>% 
        pull(ORDER), #t1[df1%>% colnames()],  
      which = c("column"), 
      col= list(Type = c("Primary" ="#6d405d", "Secondary"  = "#093b6d", "Other"  = "gray")), #  "#6d405d", "#093b6d",
      #name= "sfs",  
      #annotation_name_side = "top", 
      show_annotation_name = TRUE
    ), 
    
    right_annotation = HeatmapAnnotation(
      Dose = meta %>% pull(ARM_f),
      
      #as.name(paste0("Dose(", dose_unit, ")")) = meta %>% pull(DOSE),
      Max.histopath = hist_max, # rep(c(0, 1, 2, 3), 3),
      ADA = meta %>% pull(NAB), 
      Sex = meta %>% pull(SEX_f),  
      Subject = subj_lst, # meta %>% pull(SUBJECT) %>% str_extract("(\\d)+") %>% as.integer(), 
      
      which = c("row"), 
      annotation_legend_param  = list(legend_direction = "horizontal", direction = "horizontal"), # legend_direction = c("vertical", "horizontal"),
      col= list(Sex = c("Male" = "#2f71fd", "Female"  = "#f98068" ), 
                `Dose(xE12 vg/kg)` = colorRamp2(c(6, 20, 50), c("tan","orange",  "brown")), 
                Max.histopath = colorRamp2(c( 1, 2, 3), c("gray90", "gray50", "gray10")), 
                Subject = colorRamp2(subj_lst,  gg_color_hue(n=length(subj_lst)))
                
      )  # c("6" = "gray90", "20" = "gray50", "50" = "gray20"))
    ), 
    
    
    cell_fun = function(j, i, x, y, width, height, fill) 
    { 
      #grid.text(df2[i, j], x, y, gp = gpar(fontsize = 12)), 
      grid.points(x, y, size=unit(as.integer(df2[i, j]), "char"),  gp=gpar(fontsize=5))     #  https://stackoverflow.com/questions/68425129/what-does-the-fontsize-represent-in-grid-graphics
    }
    
  )
  
  
  return(list(meta=meta, df1 = df1, scale_constant =scale_constant, ht = ht))
  
}






my_pheatmap <- function(dat) { 
 
  # https://stackoverflow.com/questions/38727324/adding-a-dendrogram-to-a-ggplot2-heatmap
  # https://stackoverflow.com/questions/15505607/diagonal-labels-orientation-on-x-axis-in-heatmaps
  # https://www.datanovia.com/en/lessons/heatmap-in-r-static-and-interactive-visualization/
  # https://stackoverflow.com/questions/36421231/r-markdown-output-size
  #   
  
  df <- dat %>% 
    filter(ATPT %in% c("Week 26/27"), 
           PARAMCD %in% "mRNA")  %>% 
    select(SUBJECT, SEX, DOSE, MATRIX, PARAMCD, AVAL) %>% # , HISTOPTH) %>% 
    
    mutate(AVAL = as.numeric(AVAL), 
           SUBJECT = paste0(SUBJECT, "(", DOSE, ", ", SEX, ")")
    ) %>% 
    select(-DOSE) %>% 
    pivot_wider(id_cols=c("SUBJECT"), names_from=MATRIX, values_from = AVAL)  # , "HISTOPTH"
  
  df_histopath <- dat %>% 
    filter(ATPT %in% c("Week 26/27"), 
           PARAMCD %in% "mRNA")  %>% 
    select(SUBJECT, SEX, DOSE, MATRIX, PARAMCD,   HISTOPTH) %>% # , HISTOPTH) %>% 
    
    mutate(HISTOPTH = as.character(HISTOPTH), 
           HISTOPTH = ifelse(HISTOPTH=="0", "", HISTOPTH), 
           SUBJECT = paste0(SUBJECT, "(", DOSE, ", ", SEX, ")")
    ) %>% 
    select(-DOSE) %>% 
    pivot_wider(id_cols=c("SUBJECT"), names_from=MATRIX, values_from = HISTOPTH)  # , "HISTOPTH"
  
  library(janitor)
  df <- df %>% janitor::remove_empty(which = "cols", cutoff=0.6)
  df_histopath <- df_histopath  %>% select(colnames(df))
  #df_histopath[df_histopath==0] = ""
  # mutate(
  #   across(everything(), ~replace_na(.x, 10^50))
  # )
  
  
  # https://stackoverflow.com/questions/45576805/how-to-replace-all-na-in-a-dataframe-using-tidyrreplace-na
  # select(1:32)   %>%
  #   filter(if_any(everything(), ~!is.na(.)))   # https://stackoverflow.com/questions/41609912/remove-rows-where-all-variables-are-na-using-dplyr
  # filter(!across(everything(), is.na))
  
  # 6503  has 3 histopath
  
  #df <- df %>% select(-4)
  df1 <- df %>%  select(-SUBJECT)  %>% as.matrix()
  rownames(df1) <- df$SUBJECT # paste0(df$SUBJECT, "(", DOSE, ")")
  
  df_histopath1 <- df_histopath %>% select(-SUBJECT) %>% as.matrix()
  rownames(df_histopath1) <- df_histopath$SUBJECT # paste0(df$SUBJECT, "(", DOSE, ")")
  
  
  library(pheatmap)
  pheatmap::pheatmap(df1  %>% scale() , Rowv = TRUE, Colv = TRUE, angle_col = 45, na_col="white" , 
                     display_numbers =  df_histopath1, 
                     number_color = "blue", 
                     fontsize_number = 15
  )
  
}
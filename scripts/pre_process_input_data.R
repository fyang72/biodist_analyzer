

if (1==2) { 
  library(readxl)
  data_dir <- here::here(HOME, "data")
  
  excel_file <- "data_template_0208.xlsx"
  #excel_file <- "CRL31967411 _ 2025Feb10v2 _ data_ template_FY.xlsx"
  
  # Get the names of all sheets in the Excel file
  sheet_names <- excel_sheets(path = here::here(data_dir, excel_file))
  
  # Read all sheets into a list of data frames
  data_list <- lapply(sheet_names, function(sheet) {
    read_excel(path = here::here(data_dir, excel_file),
               col_names = FALSE,
               skip = 4,
               sheet = sheet)  # the top 4 lines will be always skipped
  })
  
  # Optionally, name the list elements with the sheet names for easy reference
  names(data_list) <- sheet_names
  
  
  dat0 = pre_process_input_data(data_list)
  #save(dat0, spec4matrix,file="./data/dat0_0208.RData" )
  
  
  data = dat0 %>% filter(PARAMCD_f %in% "VGC_MTH2")
  df = prepare_dataset_for_heatmap(
    data,  
    scale_data = FALSE, 
    na_threshold = 0.5 
  )
  
  # Define the custom heatmap function
  create_heatmap2(df, data,
                              col_scheme = NULL, #colorRamp2(c(0, 5, 10), c("blue", "white", "red")),                           
                              cluster_rows = FALSE, 
                              cluster_columns = FALSE, 
                              base_font_size = 6, 
                              add_aval_txt = FALSE
  )
    
    
}



  
# load all sheets from excel (See example dataset) 
loading_from_excel <- function(file) {
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
  data_list
}



pre_process_input_data <- function(data_list) {
   
 
# spec4matrix
#spec4matrix = read_xlsx(path = here::here(data_dir, input_file), sheet = "Matrix" )
spec4matrix = data_list$`Biological Matrix` %>% 
  janitor::row_to_names(row_number = 1) %>%   suppressWarnings() %>%
  dplyr::rename( 
    MATRIX = `Biological Matrix`, 
    TISSUE = Tissue, 
    MATRIXCD = `Biological Matrix Short`, 
    TISSUECD = `Tissue Short`
  ) %>% 
  dplyr::rename_with(stringr::str_to_upper)   

# spec4matrix[duplicated(spec4matrix$MATRIX), ]
# spec4matrix[duplicated(spec4matrix$MATRIXCD), ]


# Study Design
meta = NULL
meta$STUDY = parse_meta(data_list$`Study Design`, key="STUDY")
meta$GROUP = parse_meta(data_list$`Study Design`, key="GROUP")
meta$ATPT = parse_meta(data_list$`Study Design`, key="TIMEPOINT")
meta$PARAMCD = parse_meta(data_list$`Study Design`, key="PARAMCD")
meta$SEX = parse_meta(data_list$`Study Design`, key="SEX")



# assembly all
selected_vars = 
  c("STUDYID", "USUBJID", "GROUPN", "GROUP_f", "SEXN", "SEX_f", "DOSN", "DOSU", 
    "MATRIX", "ATPTN", "ATPT_f", "PARAMCD_f", "PARAM_f", "AVAL", "LLOQ")

dat0 <- NULL
sheet_name_lst <-  names(data_list)[which(str_detect(names(data_list), "Viral|VGC|RNA|Protein"))]
for (isheet in sheet_name_lst) {  
  #print(isheet)
  dat0 <- dat0 %>% 
    bind_rows(assembly_input_data(data_list, sheet_names=isheet) %>% add_meta_info(meta, selected_vars))
}

dat0   %>% 
 left_join(spec4matrix %>% select(MATRIX, MATRIXCD, TISSUE, TISSUECD), by=c("MATRIX"))  

}



# meta from "Study Design"
# Define function to parse meta information
parse_meta <- function(meta, key="STUDY") {
  
  library(janitor)
  library(dplyr)
  
  # Define the list of column name vectors
  column_names_list <- list(
    STUDY = c("STUDYID", "TITLE", "NOTE"),
    GROUP = c("GROUPN", "GROUP_f", "DOSN", "DOSU", "MALE", "FEMALE", "NOTE"),
    TIMEPOINT = c("ATPTID", "ATPT_f", "ATPTN", "NOTE"),
    PARAMCD = c("PARAMID", "PARAMCD_f", "PARAM_f", "UNIT", "LLOQ", "NOTE"),
    SEX = c("SEXN", "SEX_f", "NOTE")
  )
  
  # Assign new column names based on the key
  new_colnames <- column_names_list[[key]]
  
  # Handle the case when key is not found
  if (is.null(new_colnames)) {
    new_colnames <- NA
  }
  
  df = meta %>% 
    dplyr::filter(if_any(1, ~ . == key)) %>% 
    janitor::row_to_names(row_number = 1) %>%   suppressWarnings() %>%
    # “snake”, “small_camel”, “big_camel”, “screaming_snake”, “parsed”, “mixed”, 
    # “lower_upper”, “upper_lower”, “swap”, “all_caps”, “lower_camel”, “upper_camel”,
    # “internal_parsing”, “none”, “flip”, “sentence”, “random”, “title”
    clean_names(case="all_caps") %>%   
    dplyr::select(2:NOTE)  %>% 
    stats::setNames(new_colnames)
  
  if (key=="GROUP") {
    df <- df %>% 
      dplyr::distinct(GROUPN, GROUP_f, DOSN, DOSU) %>% 
      dplyr::mutate(GROUPN = as.integer(GROUPN),
                    DOSN = as.numeric(DOSN)) %>% 
      dplyr::arrange(GROUPN)  %>% 
      dplyr::mutate(GROUP_f = ordered(GROUP_f, levels=GROUP_f))
  }
  
  if (key=="TIMEPOINT") {
    df <- df %>% 
      dplyr::distinct(ATPTN, ATPT_f) %>% 
      dplyr::mutate(ATPTN = as.numeric(ATPTN)) %>% 
      dplyr::arrange(ATPTN) %>% 
      dplyr::mutate(ATPT_f = ordered(ATPT_f, levels=ATPT_f))
  }
  
  if (key=="SEX") {
    df <- df %>% 
      dplyr::distinct(SEXN, SEX_f) %>% 
      dplyr::mutate(SEXN = as.integer(SEXN)) %>% 
      dplyr::arrange(SEXN) %>% 
      dplyr::mutate(SEX_f = ordered(SEX_f, levels=SEX_f))
  }
  
  if (key=="PARAMCD") {
    df <- df %>% 
      dplyr::mutate(PARAMID = as.integer(PARAMID)) %>% 
      dplyr::arrange(PARAMID) %>% 
      dplyr::distinct(PARAMCD_f, PARAM_f) %>% 
      dplyr::mutate(PARAMCD_f = ordered(PARAMCD_f, levels=PARAMCD_f), 
                    PARAM_f = ordered(PARAM_f, levels=PARAM_f)) 
  }  
  
  df
}


# assembly_input_data  
assembly_input_data <- function(data_list, sheet_names="Protein (Brain)") {
  
  df = data_list[[sheet_names]]
    
  key_word <- case_when(
    str_detect(toupper(sheet_names), "VIRAL") ~ "VIRAL", 
    str_detect(toupper(sheet_names), "VGC") ~ "VGC", 
    str_detect(toupper(sheet_names), "RNA") ~ "RNA", 
    str_detect(toupper(sheet_names), "PROTEIN") ~ "PROTEIN", 
    TRUE ~ NA
  )
  
  df %>% 
    janitor::row_to_names(row_number = 1) %>%   suppressWarnings() %>% 
    janitor::clean_names(case="all_caps") %>%    
    dplyr::rename(
      USUBJID = SUBJECT,
      MATRIX = BIOLOGICAL_MATRIX, 
      GROUPN = GROUP, 
      SEX_f = SEX,
      ATPT_f = TIME_POINT
    )   %>% 
    dplyr::mutate(
      USUBJID = as.character(USUBJID),
      GROUPN = zoo::na.locf(GROUPN,  fromLast=FALSE) %>% as.integer(),  
      SEX_f = zoo::na.locf(SEX_f,  fromLast=FALSE), 
      MATRIX = zoo::na.locf(MATRIX,  fromLast=FALSE), 
      ATPT_f = zoo::na.locf(ATPT_f,  fromLast=FALSE)
    ) %>% 
    dplyr::select(USUBJID, GROUPN, SEX_f, MATRIX, ATPT_f, starts_with(key_word), starts_with("LLOQ")) %>% 
    tidyr::pivot_longer(cols=starts_with(key_word), names_to = "PARAMCD_f", values_to = "AVAL")  %>%
    tidyr::pivot_longer(cols=starts_with("LLOQ"), names_to = "LLOQ_TYPE", values_to = "LLOQ")  %>% 
    dplyr::select(-LLOQ_TYPE) %>% 
    
    mutate(
      AVAL = suppressWarnings(as.numeric(AVAL)),
      MATRIX = ifelse(stringr::str_sub(MATRIX, -5)=="Right", stringr::str_sub(MATRIX, 1, -8), MATRIX), 
      MATRIX = ifelse(stringr::str_sub(MATRIX, -4)=="Left", stringr::str_sub(MATRIX, 1, -7), MATRIX)  
    )
  
  
  
}


# add_meta_info
add_meta_info <-  function(
    df, 
    meta, 
    selected_vars = 
      c("STUDYID", "USUBJID", "GROUPN", "GROUP_f", "SEXN", "SEX_f", "DOSN", "DOSU",
        "MATRIX", "ATPTN", "ATPT_f", "PARAMCD_f", "PARAM_f", "AVAL", "LLOQ")) {
  
  df <- df %>% 
    
    mutate(STUDYID = meta$STUDY %>% dplyr::pull(STUDYID)%>% unique() %>% first())  %>%  # add study
    
    left_join(meta$GROUP, by="GROUPN")  %>%  # add ordered GROUP_f and numeric DOSE
    
    left_join(meta$ATPT, by="ATPT_f")  %>%  # add numeric ATPT  
    mutate(ATPT_f = ordered(ATPT_f, levels=meta$ATPT %>% pull(ATPT_f) %>% levels())) %>% # ordered ATPT_f
    
    left_join(meta$SEX, by="SEX_f")  %>% # add numeric SEX  
    mutate(SEX_f = ordered(SEX_f, levels=meta$SEX %>% pull(SEX_f) %>% levels()))  %>%      # ordered SEX_f
    
    left_join(meta$PARAMCD, by="PARAMCD_f")  %>% # add PARAM_f
    mutate(
      PARAMCD_f = ordered(PARAMCD_f, levels=meta$PARAMCD %>% pull(PARAMCD_f) %>% levels()), 
      PARAM_f = ordered(PARAM_f, levels=meta$PARAMCD %>% pull(PARAM_f) %>% levels()),       
    )  %>%      # ordered SEX_f
    
    dplyr::select(all_of(selected_vars))  
  
  df
}









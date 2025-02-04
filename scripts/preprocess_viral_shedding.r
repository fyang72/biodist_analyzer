

preprocoess_viral_shedding <- function(df) {
#-----------------------------
# dat_shedding     
#-----------------------------
df = read_xlsx(
  path = here::here(data_dir, inFile), 
  sheet = "Shedding", 
  range = "A5:N592")   

colnames(df) <- colnames(df) %>% stringr::str_to_upper()    # stringr::str_to_title
  
library(zoo)
df = df %>% 
  dplyr::rename(
    USUBJID = SUBJECT,
    MATRIX = `BIOLOGICAL MATRIX`, 
    GROUP = GROUP, 
    SEX = SEX,
    ATPT = `TIME POINT`, 
    LLOQ = `LLOQ (CP/RXN)`, 
    AVAL = `SAMPLE CONC.\r\n(CP/ΜL)`
  )  %>% 
  
  dplyr::mutate(
   STUDYID = as.character("XXXXXX"), 
   USUBJID = as.character(USUBJID),
   GROUP = zoo::na.locf(GROUP,  fromLast=FALSE), 
   GROUP = as.integer(GROUP), 
   SEX = zoo::na.locf(SEX,  fromLast=FALSE), 
   SEX = as.integer(factor(SEX, levels=c("Male", "Female"))),   # Male=1, Female=2
   MATRIX = zoo::na.locf(MATRIX,  fromLast=FALSE), 
   ATPT = zoo::na.locf(ATPT,  fromLast=FALSE), 
   ATPT = case_when(
     ATPT == "Pretest" ~ -3,   # arbitrary
     ATPT == "D1: Predose" ~ 0,  
     ATPT == "D1: 6HR" ~ 6, 
     ATPT == "D1: 24HR" ~ 24,   
     ATPT == "Day 3" ~ 48, 
     ATPT == "Day 8"  ~ 168,  
     ATPT == "Day 14"  ~ 312, 
     ATPT == "Interim Necropsy Week 6"  ~ 1008,
     ATPT == "Terminal Necropsy Weeks 26/27"  ~ 4368, 
     TRUE ~ NA
   ) %>% as.integer(),
   
   SAMPID = NA, 
   PARAMCD = "PK"
   )    
  
# "Feces"  "PBMC"   "Saliva" "Serum"  "Urine" 
df = df %>% 
  mutate(
    MATRIX = case_when(
      #MATRIX == "Serum" ~ "Serum, Protein", 
      TRUE ~ MATRIX
    )
  )


  return(df)

}
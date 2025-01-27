
HOME = getwd()
#-------------------------------
# globally installed on SCP
#-------------------------------
global_installed_package_lst <- 
  c("renv", "azcore",  "dplyr",  "tidyverse", "tidyr", 
    "scales", "glue", "magrittr", "haven", "readxl", "readr"
    )
#library(azcore)
azcore::azcore_library_load(libraries = global_installed_package_lst)

#-------------------------------
# locally installed packages  
#-------------------------------
local_installed_package_lst <-
  c(
  "mrgmisc", "mrggsave", "yspec",    
  #"pmplots", "ggplot2",  #package ‘ggplot2’ 3.4.4 is loaded, but >= 3.5.0 is required by ‘pmplots’
  "DoseFinding", "ComplexHeatmap", "ROCit", "janitor", "here", "ggpubr", "xgxr"
  )
#library(pacman)
pacman::p_load(char=local_installed_package_lst)
 
library(ggplot2, lib.loc = paste0(HOME, "/renv/library/R-4.3/x86_64-pc-linux-gnu"))
library(pmplots, lib.loc = paste0(HOME, "/renv/library/R-4.3/x86_64-pc-linux-gnu"))

library(ggh4x)    # for 
library(ggprism)

# library(survminer)
# library(quarto) 
# library(pmplots)  package ‘ggplot2’ 3.4.4 is loaded, but >= 3.5.0 is required by ‘pmplots’
# library(circlize)   # for color selection
# library(table1)  from Certara 
 
set.seed(1014)

knitr::opts_chunk$set(
  comment = "#>",
  collapse = TRUE,
  # cache = TRUE,
  fig.retina = 2,
  fig.width = 6,
  fig.asp = 2/3,
  fig.show = "hold"
)

options(
  dplyr.print_min = 6,
  dplyr.print_max = 6,
  pillar.max_footer_lines = 2,
  pillar.min_chars = 15,
  stringr.view_n = 6,
  # Temporarily deactivate cli output for quarto 
  cli.num_colors = 0,
  cli.hyperlink = FALSE,
  pillar.bold = TRUE,
  width = 77 # 80 - 3 for #> comment
)

ggplot2::theme_set(ggplot2::theme_gray(12))

# use results: "asis" when setting a status for a chapter
status <- function(type) {
  status <- switch(type,
    polishing = "should be readable but is currently undergoing final polishing",
    restructuring = "is undergoing heavy restructuring and may be confusing or incomplete",
    drafting = "is currently a dumping ground for ideas, and we don't recommend reading it",
    complete = "is largely complete and just needs final proof reading",
    stop("Invalid `type`", call. = FALSE)
  )

  class <- switch(type,
    polishing = "note",
    restructuring = "important",
    drafting = "important",
    complete = "note"
  )

  cat(paste0(
    "\n",
    ":::: status\n",
    "::: callout-", class, " \n",
    "You are reading the work-in-progress second edition of R for Data Science. ",
    "This chapter ", status, ". ",
    "You can find the complete first edition at <https://r4ds.had.co.nz>.\n",
    ":::\n",
    "::::\n"
  ))
}





cement <- function(...) {
  # https://adv-r.hadley.nz/quasiquotation.html
  args <- ensyms(...)
  paste(purrr::map(args, rlang::as_string), collapse = " ")
}

mk_col_title_spec <- function(spec, what, sep="//") { 
  #what = cement({{what}})
  
  title <- paste0(
    spec[[what]][["col"]],
    sep,
    spec[[what]][["short"]]
  )
  
  glue::glue(title, xunit = mrgmisc::parens(spec[[what]][["unit"]]))
}

glue_unit_spec <- function(spec, what) { 
  what = cement({{what}})
  return( glue::glue(spec[[what]][["short"]], xunit = mrgmisc::parens(spec[[what]][["unit"]]))) 
}


# lab <- ys_get_short_unit(spec, parens = TRUE, title_case = TRUE)
#y = pmplots:::col_label(mk_col_title_spec(spec, cement({{.aval}})))
# pmplots:::require_numeric(df, y[1])



gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}


color.scheme.certara = c(
  "#4682ac", "#ee3124", "#fdbb2f", "#6d405d", "#093b6d",
  "#2f71fd", "#336343", "#803333", "#279594", "#ef761b",
  "#29398c", "#32a17e", "#d89a17", "#d64d20", "#9da1bd",
  "#9c8777", "#7059a6", "#e07070", "#475c6b", "#75604D",
  "#067f97", "#b7a148", "#f98068", "#72cbed", "#b8a394",
  "#b35d1b", "#a52f43", "#113df2", "#f2c611", "#52ccbb")



#ggalt  annotation_ticks
theme_prism2 <- function(palette = "black_and_white", base_size = 10,
                         base_family = "sans", base_fontface = "bold",
                         base_line_size = base_size/14, base_rect_size = base_size/14,
                         axis_text_angle = 0, border = FALSE) {
  
  library(ggprism)
  theme_prism(palette = palette,
              base_size = base_size,
              base_family = base_family,
              base_fontface = base_fontface,
              base_line_size = base_line_size,
              base_rect_size = base_rect_size,
              axis_text_angle = axis_text_angle,
              border = border ) +
    
    ggplot2::theme(
      
      #legend.text = ggplot2::element_text(size = legend_text),
      #legend.title = element_blank(),  # ggplot2::element_text(size = legend_title),     # Remove only the legend title by set legend.title = element_blank()
      legend.key =  ggplot2::element_blank() ,
      # left,top, right, bottom.
      legend.position = "bottom",       # legend.position='none', Remove the plot legend   c(0,0) corresponds to the "bottom left" and c(1,1) corresponds to the "top right" position.
      #legend.justification = c(1, 1),     ??
      legend.box = "horizontal" ,  # Horizontal legend box
      legend.background = element_rect(fill="white", size=0.5, linetype="solid", colour ="white"),
      
      panel.background = ggplot2::element_blank(),
      #panel.grid = ggplot2::element_blank(),
      panel.grid.minor = ggplot2::element_blank(),   # element_line(colour = "gray98",size=0.5),   #
      panel.grid.major = element_line(colour = "gray97",size=0.75),  #ggplot2::element_blank()   # element_line(colour = "gray90",size=0.75))
      
      panel.border = element_rect(colour = "black", fill=NA, size=1),
      panel.spacing = unit(0.2, "lines"),
      strip.background = element_rect(colour="transparent", fill="white"),
      
      #plot.margin = unit(c(1,1,1,1), "cm") +
      plot.margin = unit(c(0.1,0.5,0.1,0.1), "cm")
      
    ) #+
  #guides(col=guide_legend(ncol=5, byrow=TRUE))
  
  
}

# theme(
#   panel.grid.major = element_line(colour = "black"),
#   panel.grid.minor = element_line(colour = "grey50")
# )


################################################################################
# Notes
################################################################################


# linux command: 
# orts things in reverse numerical order,  https://stackoverflow.com/questions/1019116/using-ls-to-list-directories-and-their-total-sizes
# du -sh * | sort -hr


# azcore: 
# https://azcollaboration.sharepoint.com/sites/AZ245/R_Governance_ForGxP_Usage/r_package_azcore/SitePages/Getting-Started.aspx
#  get the list of available bundles
# result <-  azcore::azcore_bundle_available() # https://azcollaboration.sharepoint.com/sites/AZ245/R_Governance_ForGxP_Usage/rcore_and_bundles/SitePages/Home.aspx
# library(azcore)
 
# #  load the default version of the community bundle using the :: operator
# azcore::azcore_bundle_load( "community" )
# 
# #  list of loaded bundles
# .bundles
# 
# 
# #  load a single R package
# azcore_library_load( libraries = "one" )
#  
# #  load a sequence of R packages
# azcore_library_load( libraries = c("one", "two") )



if (1==2) { 
  library(dplyr)
  library(magrittr)
  path1 = "/opt/scp/services/azcore/corepackages/3.0.1/libraries/R/4.3.1"
  lib1 = list.dirs(path = path1, full.names = FALSE, recursive = FALSE)
  
  path2 = "/opt/scp/services/azcore/coreutils/0.1.0/libraries/R/4.3.1"
  lib2 = list.dirs(path = path1, full.names = FALSE, recursive = FALSE)
  
  path3 = "/opt/scp/apps/gen-2021a/software/R/4.3.1-foss-2021a-core/lib64/R/site-library"
  lib3 = list.dirs(path = path1, full.names = FALSE, recursive = FALSE)
  
  path4 = "/opt/scp/apps/gen-2021a/software/R/4.3.1-foss-2021a-core/lib64/R/library"
  lib4 = list.dirs(path = path1, full.names = FALSE, recursive = FALSE)

  lib_all <- installed.packages() %>% as.data.frame() %>% select(Package, LibPath) %>% pull(Package)
 
  addl_packages <- c("azcore", "renv")
  
  #from cache:
  from_cache_packages_lst <- c("base",      "compiler",  "datasets",  "graphics",  "grDevices", "grid" ,
   "methods" ,  "parallel",  "splines",   "stats" ,    "stats4" ,   "tcltk"  ,   "tools" ,    "utils")
   
  # library(azcore)
  azcore::azcore_bundle_available()
  path = "/opt/scp/services/azcore/biometrics/0.3.1/libraries/R/4.3.1"
  lib = list.dirs(path = path, full.names = FALSE, recursive = FALSE)   # 313
  
  path =  "/home/kwqz358/01_alxn_projects/cpss-dsu-toxbiodist-analyzer/renv/library/R-4.3/x86_64-pc-linux-gnu"
  lib_installed <- list.dirs(path = path, full.names = FALSE, recursive = FALSE)   # 313  installed.packages() %>% as.data.frame() %>% select(Package, LibPath) %>% pull(Package)
  
  package_list_to_be_removed <- lib_installed[which(lib_installed %in% lib)] %>% unique() %>% sort()
  
  library(fs)
  LIB_PATH <- "/home/kwqz358/01_alxn_projects/cpss-dsu-toxbiodist-analyzer/renv/library/R-4.3/x86_64-pc-linux-gnu"
  
  # not run this line!!!!!!!!!!
  # fs::dir_delete(paste0(LIB_PATH, "/", package_list_to_be_removed ))
  # not run this line!!!!!!!!!!
}


if (1==2) { 
  ##########################################################################################
  # Is it a good practice to push renv folder and renv.lock file to github along with other files
  ##########################################################################################
  # https://forum.posit.co/t/is-it-a-good-practice-to-push-renv-folder-and-renv-lock-file-to-github-along-with-other-files/151208
  
  library(azcore)
    
  # disable automatic snapshots
  auto.snapshot <- getOption("renv.config.auto.snapshot")
  options(renv.config.auto.snapshot = FALSE)
  
  # initialize a new project (with an empty R library)
  renv::init(bare = TRUE)
  
  # install digest 0.6.19
  renv::install("digest@0.6.19")
  
  # save library state to lockfile
  renv::snapshot(exclude= c("palmerpenguins", "pheatmap", "pryr" ))
  
  # remove digest from library
  renv::remove("digest")
  
  # check library status
  renv::status()
  
  # restore lockfile, thereby reinstalling digest 0.6.19
  renv::restore()
  
  # restore automatic snapshots
  options(renv.config.auto.snapshot = auto.snapshot)


}





#####################################



HOME = getwd()
data_dir <- here::here(HOME, "data/source/ALXN2350-GLP-MKY/")
status = "DRAFT"  #flag for labeling figures as draft
#source(here::here(HOME, "_common.R"))

source("./scripts/utils.R")
source("./scripts/my_barplot.R")
source("./scripts/my_heatmap.R")
source("./scripts/my_plot_conc_profile.R")

actionButton_style ="float:left;color: #fff; background-color: #328332; border-color: #328332"

# ggplot settings, xgx_theme_set()
options(mrggsave.dir = here("deliv/figure"), mrg.script = "app.qmd")
options(dplyr.summarise.inform = TRUE)

####################################################

# dat0 <- readr::read_csv(file=here(HOME, "data", "derived", "alxn2340_biodistribution_data.csv"))

if (1==2) { 
  data_dir = paste0(HOME, "/data/source/ALXN2340-GLP-MKY/")
  
  spec_location <- here(data_dir, "preclinical2.yml")
  spec <- ys_load(spec_location) %>% ys_namespace("plot")
  
  library(readr)
  dat0 <- readr::read_csv(paste0(data_dir, "/output/glp_exam.csv"))
  colnames(dat0) <- toupper(colnames(dat0))
  
  dat0 <- dat0 %>% 
    dplyr::rename(
      STUDYID = STUDY_ID, 
      SUBJECT = SUBJECT, 
      SEX = GENDER, 
      
      PARAMCD_GROUP = TEST,
      PARAMCD = TESTCD, 
      MATRIX = BIOLOGICAL_MATRIX, 
      MATRIXCD = TISDETAIL,
      
      ATPT = TIME_POINT, 
      
      AVAL = DV, 
      AVALU =DV_UNIT, 
      BLQ = BLQ, 
      
      ADA = ADA, 
      NAB = NAB_PREDOSE, 
      HISTOPTH = HISTOPATH,   #  max 8 characters
      
      ARM = DOSEGRP, 
      #ARMN = GROUP, 
      DOSE = DOSE.E12.VG.KG, 
      
      PFNO = PFNO, 
      ORDER = ORDER, 
      CAT2 = CAT2  # NOT SURE
    ) %>% 
    
    mutate(
      ARM = ordered(ARM, levels = dat0$DOSEGRP %>% unique()) %>% as.integer(), 
      SEX = case_when(
        SEX == "Male" ~ 0, 
        SEX == "Female" ~ 1
      ), 
      ATPT = case_when(
        ATPT == "Pretest" ~ 0,
        ATPT == "Week 6" ~ 6, 
        ATPT == "Week 26/27" ~ 26, 
        TRUE ~ NA
      ), 
      # ordered(ATPT, levels = c("Pretest", "Week 6" , "Week 26/27"     )) %>% as.integer, 
      
      TIME = ATPT,   # not in numeric
      AVAL = as.numeric(AVAL), 
      PARAMCD_GROUP = ifelse(PARAMCD_GROUP %in% "D4-PGRN", "PGRN", PARAMCD_GROUP),
      TISSUE = ifelse(PARAMCD %in% c("vgc", "mRNA"),  "Whole", TISSUE), 
      MATRIX = ifelse(MATRIX %in% "Lungs" & PARAMCD %in% "mRNA", "Lung", MATRIX)
    )
  
  #dat0 %>% pull(MATRIX) %>% unique()  %>% sort() #
  
  dat0 <- dat0 %>% 
    mutate(MATRIX_GROUP = case_when(
      stringr::str_detect(MATRIX, "Brain") ~ "Brain", 
      stringr::str_detect(MATRIX, "Dorsal Root Ganglion") ~ "Dorsal Root Ganglion", 
      stringr::str_detect(MATRIX, "Liver") ~ "Liver", 
      stringr::str_detect(MATRIX, "Spinal Cord") ~ "Spinal Cord", 
      stringr::str_detect(MATRIX, "Heart") ~ "Heart", 
      stringr::str_detect(MATRIX, "Testis") | stringr::str_detect(MATRIX, "Ovary")  ~ "Testis or Ovary", 
      TRUE ~ MATRIX
    )
    )
  
  #readr::write_csv(dat0, file=here(HOME, "data", "derived", "alxn2340_biodistribution_data2.csv"))
  
  spec_location <- here::here(data_dir, "preclinical2.yml")
  spec <- ys_load(spec_location) %>% ys_namespace("plot")
  
  dat0 <- dat0 %>% 
    yspec_add_factors(spec, STUDYID, ARM, DOSE, PARAMCD, SEX, ATPT, BLQ, TISSUE) #%>%  # Refactor table for plots  
  
  # there is no "P0303" and "P0603" in dat0
  # dat0 <- dat0 %>%
  #   yspec_add_factors(spec, STUDYID, ARM, DOSE, SEX, ATPT, BLQ, TISSUE) #%>%  # Refactor table for plots
  
  dat0 <- dat0 %>% 
    mutate(ARM2 = paste0("Grp", GROUP, "(", as.character(ARM_f), ")"), 
           #ARM2 = ordered(ARM2, levels =  paste0("Grp", dat0$GROUP%>% unique(), "(", levels(dat0$ARM_f), ")")   ), 
           
           #SEX_f = SEX, 
           USUBJID = paste0("P", SUBJECT)
    ) 
  
  
  }

 
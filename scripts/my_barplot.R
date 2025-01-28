


my_barplot = function(
    df,    
     facet_by = "PARAMCD"
) {
  
  # https://stackoverflow.com/questions/52506296/ggplot-geom-point-position-jitterdodge-not-working-when-color-specified
  p <-  
    ggplot(df, aes(x = MATRIXCD, y = AVAL, fill = GROUP_f, group = GROUP_f)) +
    geom_bar(
      position = "dodge", width = 0.8, stat = "summary", fun = "mean",
      color = "black", linewidth = .8
    ) +
    scale_fill_manual(values = color.scheme.certara[3:length(color.scheme.certara)]) + # gg_color_hue(n=12)) + 
    
    # scale_color_discrete( c("gray", "blue", "cyan"), 
    #                      labels = c("Young Adult", "Older Adult", "sfsf")) +
    stat_summary(
      fun.data = mean_sdl, geom = "errorbar", color = "black",
      position = position_dodge(0.8), width = 0.2, linewidth = 0.8
    ) +
    
    #scale_color_manual(values = c("gray", "blue", "cyan", "red", "gray90")) +    # for 0, 1, 2, 3, NA
    #scale_size_manual(values = c(1, 2, 3, 4, 5)) + 
    
    scale_y_log10() + 
    #scale_color_brewer(palette = "Set1") +   # , name = "color_legend_name") +
    #facet_wrap(~TIME, nrow=2)+ theme_bw() +
    
    labs(x="", y= "Expression Levels") + # , fill = fill_legend_txt, color= histopah_legend_txt, size=histopah_legend_txt) + 
    
    theme_prism2() + 
    theme(legend.position="bottom", legend.box="vertical")  +
    theme(axis.text.x = element_text(size=10, angle = 45, hjust = 1, vjust=1, face="plain")) +  # 	Font face ("plain", "italic", "bold", "bold.italic")
    theme(legend.title = element_text())  
  
  #facet_grid(cols = vars(ARM_f), scales = "free") +
  #facet_wrap(as.formula(paste("~", cement({{.param}}))), ncol=1, scales="free_y")
   
   p = p + facet_wrap(as.formula(paste("~", facet_by)), ncol=1, scales="free_y")
   p
 
   
}   




# barplot_for_biodistribution ----------------------------
##' df = dataframe for plotting
##' spec = data specification from yaml file
##' .param = parameter to be used on x axis (default = PARAM or PARAM_GROUP)
##' .aval = analysis value to be used on y axis (default = AVAL)
##' .id = id to identify an individual subject  
##' .time = time point for the sample collection  
##' status = draft or final version of the plot

my_barplot_org = function(
    df,  spec,  
    
    .xmatrix = MATRIX, 
    .yaval = AVAL, 
    .param = PARAMCD, 
    .id = ID,
    .facet = PARAMCD,
    .fill = ARM_f,  
    .histopath = HISTOPTH, 
    
    fill_legend_txt = "Treatement Arm:",
    histopah_legend_txt = "Histopath",
    status = "DRAFT"  #flag for labeling figures as draft
) {
  
  aval_txt <- cement({{.yaval}})
  my_label <- yspec::ys_get_short(spec, short_max = Inf, title_case = TRUE)
  my_unit <- yspec::ys_get_unit(spec, short_max = Inf, title_case = TRUE)
  ylab_txt <- glue::glue(my_label[[aval_txt]], xunit=mrgmisc::parens(my_unit[[aval_txt]]))
  
  
  # https://stackoverflow.com/questions/52506296/ggplot-geom-point-position-jitterdodge-not-working-when-color-specified
  p <-  
    ggplot(df, aes(x = {{.xmatrix}}, y = {{.yaval}}, fill = {{.fill}}, group = {{.fill}})) +
    geom_bar(
      position = "dodge", width = 0.8, stat = "summary", fun = "mean",
      color = "black", linewidth = .8
    ) +
    scale_fill_manual(values = color.scheme.certara[3:length(color.scheme.certara)]) + # gg_color_hue(n=12)) + 
    
    # scale_color_discrete( c("gray", "blue", "cyan"), 
    #                      labels = c("Young Adult", "Older Adult", "sfsf")) +
    stat_summary(
      fun.data = mean_sdl, geom = "errorbar", color = "black",
      position = position_dodge(0.8), width = 0.2, linewidth = 0.8
    ) +
    
    scale_color_manual(values = c("gray", "blue", "cyan", "red", "gray90")) +    # for 0, 1, 2, 3, NA
    scale_size_manual(values = c(1, 2, 3, 4, 5)) + 
    
    scale_y_log10() + 
    #scale_color_brewer(palette = "Set1") +   # , name = "color_legend_name") +
    #facet_wrap(~TIME, nrow=2)+ theme_bw() +
    
    labs(x="", y= "Expression Levels", fill = fill_legend_txt, color= histopah_legend_txt, size=histopah_legend_txt) + 
    
    theme_prism2() + 
    theme(legend.position="bottom", legend.box="vertical")  +
    theme(axis.text.x = element_text(size=10, angle = 45, hjust = 1, vjust=1, face="plain")) +  # 	Font face ("plain", "italic", "bold", "bold.italic")
    theme(legend.title = element_text())  
  
  #facet_grid(cols = vars(ARM_f), scales = "free") +
  #facet_wrap(as.formula(paste("~", cement({{.param}}))), ncol=1, scales="free_y")
  
  facet_txt = cement({{.facet}})
  if (!is.null(facet_txt)) {
    p = p + facet_wrap(as.formula(paste("~", facet_txt)), ncol=1, scales="free_y")
  }
  
  histopath_txt = cement({{.histopath}})
  print(histopath_txt)
  if (histopath_txt %in% colnames(df)) {
    p = p + geom_point( aes(color={{.histopath}}, size={{.histopath}}),   # 
              position = position_jitterdodge(0.3, dodge.width = .8),
              alpha = 0.8#, 
              #show.legend = FALSE
    )  
  }
  
  # if (!is.null(status)) { 
  #   p = p + xgxr::xgx_annotate_status(status) 
  # }
  
  # 
  # library(plotly) 
  # add_aes <- function (mapping, ...) {
  #   new_aes <- structure(append(mapping, as.list(match.call()[-(1:2)])), class = "uneval")
  #   ggplot2:::rename_aes(new_aes)
  # }
  
  #fig$mapping = add_aes(fig$mapping, label=SUBJECT)  
  
  # boxplot for DOSE vs HISTOPTH
  # fig <-  ggplot(df, aes(x = HISTOPTH, y = DOSE, fill=DOSE)) +
  #   geom_boxplot()
  
  
  p
}

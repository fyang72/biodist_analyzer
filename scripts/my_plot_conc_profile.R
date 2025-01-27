
# plot_conc_profile  ----------------------------
##' df = dataframe for plotting
##' spec = data specification from yaml file

##' .x = column to be used on x axis (default = TIME)
##' .y = column to be used on y axis (default = DV)
##' .id = column name to use for grouping data in individual plots
##' .facet = column name of grouping variable for facet wrap
##' .color = column name to use for colouring of data in plots
##' .blq = column name for BLQ samples of data in plots

##' status = draft or final version of the plot
##' logy = put y on a log scale T/F
##' time_units_plot = label for y-axis
##' color_legend_name = name for color_legend_name items
##' shape_legend_name = name for color_legend_name items
my_plot_conc_profile = function(df, 
                             spec,  
                             
                             .x = TIME, 
                             .y = AVAL, 
                             .id = ID,
                             .facet = ARM_f,
                             .color = ARM_f, 
                             .blq = BLQ_f,
                             
                             status = "DRAFT",  #flag for labeling figures as draft
                             logy = TRUE,  
                             time_units_plot = NULL, 
                             
                             color_legend_name = "", 
                             shape_legend_name = "") {
  
  
  x_txt <- cement({{.x}})
  y_txt <- cement({{.y}})
  facet_txt <- cement({{.facet}})
  
  pmplots:::require_numeric(df, x_txt)
  pmplots:::require_numeric(df, y_txt)
  
  print("ok1")
  my_label <- yspec::ys_get_short(spec, short_max = Inf, title_case = TRUE)
  
  print("ok2")
  my_unit <- yspec::ys_get_unit(spec, short_max = Inf, title_case = TRUE)
  xlab <- my_label[[x_txt]]
  ylab <- glue::glue(my_label[[y_txt]], xunit=mrgmisc::parens(my_unit[[y_txt]]))
  
  p = ggplot(data = df, aes(x = {{.x}}, y = {{.y}})) +
    
    geom_line(aes(group = {{.id}}), color = rgb(0.5,0.5,0.5), size = 1, alpha = 0.3) + 
    
    geom_point(data = df, aes(shape = factor({{.blq}}), alpha = 0.3), size = 2, alpha = 0.3) + 
    scale_shape_manual(values=c(1,8), name= shape_legend_name) + 
    scale_color_brewer(palette = "Set1", name = color_legend_name) + 
    
    xgx_stat_ci(aes(x = {{.x}}, color={{.color}}, group={{.color}}), conf_level = 0.95, na.rm = TRUE) + 
    
    theme_prism() + 
    theme(legend.position="bottom", legend.box="vertical")  + 
    theme( legend.title = element_text())  + 
    
    guides(colour = guide_legend(order = 1), # ,  title=legend_name), 
           shape = guide_legend(order = 2))  # , title="ds")
  
  if (!is.null(status)) { 
    p = p + xgx_annotate_status(status) 
  }
  
  if(!is.null(time_units_plot)) {
    p = p + 
      xgx_scale_x_time_units(
        units_dataset = spec[[cement({{.x}})]][["unit"]], 
        units_plot    = time_units_plot)  + 
      labs(y=ylab)
    
  } else {
    p = p + labs(x=xlab, y=ylab)
  }
  
  if (logy) {
    p = p + xgx_scale_y_log10()
  }
  
  p 
  
}

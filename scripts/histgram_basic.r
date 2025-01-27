


#' Creates a Concentration vs Time graph
#'
#' Creates a Concentration vs Time graph using Regeneron data.
#'
#' @param data A data frame with which to create the plots.  The data frame must conform to Regeneron standards and contain the variables: `week`, `ndy`, `concn`, `concni`, `pchg`, `actarmcd`, `ada`
#' @param drug Drug name
#' @param log Should the log axis be used?
#' @param adjust Should an adjustment be applied to the errorbars?
#' @return A ggplot graphic.
#' @import gridExtra
#' @export
#' @examples
#'   inFile <- system.file("extdata/pk_pd_ada_all3.sas7bdat", package = "pkGraph")
#'   library(haven)
#'   theData <- read_sas(inFile)
#'   theData <- dataCheck(theData)
#'   concTime(data = theData)
#'   concTime(data = theData, adjust = TRUE)
base_hist <- function(adpx, xval="WEIGHTBL",
                      group="SEX",
                      binwidth=10,
                      position="dodge",
                      title="Histogram plot of Weight",
                      xlab.txt="Weight(kg)"
                      ) {

  lazy_plot =  lazyeval::interp(~ggplot(adpx, aes(x=WEIGHTBL, color=GROUP, fill=GROUP)) +
  geom_histogram(aes(y=..density..),
                    binwidth=binwidth,
                    colour="black",
                    position=position, #"dodge", #   identity  dodge
                    alpha=0.5,
                    #fill="white",
                    linetype="solid"   # dashed
                    ) +

  geom_density(alpha=0.5) +     # , fill="#FF6666"

  #scale_color_manual(values=c("#999999", "#E69F00", "#56B4E9")) +     # to use custom colors
  #scale_color_grey() +
  #scale_color_brewer(palette="Paired") +    # Continuous colors
  scale_color_brewer(palette="Dark2") +    # Discrete colors
  #scale_color_brewer(palette="Accent") +    # Gradient colors

  #scale_fill_manual(breaks=c("F", "M"), values=c("blue","green","red","orange"))  +
  #scale_fill_brewer() +     # to use color palettes from RColorBrewer package
  #scale_fill_grey() +   #  to use grey color palettes
  scale_fill_brewer(palette="Dark2") +    # Discrete colors

  labs(title=title,x=xlab.txt, y = "Frequency")+     # Density

  theme_classic() +
    #base_theme(font.size=14, legend_position="bottom") +
  #facet_wrap(~SEX, scales = "free") +
  #facet_grid( .~SEX)  +

  geom_vline(data =adpx %>%  group_by(GROUP) %>% dplyr::summarize(meanWT = median(WEIGHTBL, na.rm=TRUE)),
             aes(xintercept = meanWT, color=GROUP), size = 0.8,  linetype="dashed", show.legend=FALSE),  # color = "red",
#  geom_bar() +
#  geom_bar(colour="black", show_guide=FALSE)

  WEIGHTBL = as.name(xval),
  GROUP = as.name(group))
return(lazyeval::lazy_eval(lazy_plot))
  }


#
#  #
##set.seed(1234)
#df <- data.frame(
#  SEX=factor(rep(c("F", "M"), each=200)),
#  WEIGHTBL=round(c(rnorm(200, mean=55, sd=5), rnorm(200, mean=65, sd=5)))
#  )
#
# base_hist(df, xval="WEIGHTBL", group="SEX")
#
##
#
#base_hist(pkpd_final, xval="WEIGHTBL", group="SEX")
#
#
#

#https://stackoverflow.com/questions/70046492/non-linear-regression-line-and-its-computation-failed-in-stat-smooth-argume

# https://www.rdocumentation.org/packages/DoseFinding/versions/1.1-1/topics/DR-Models

# https://www.kristianbrock.com/post/emax-intro/
  
#https://www.rdocumentation.org/packages/DoseFinding/versions/1.1-1/topics/DR-Modelslibrary(DoseFinding)


# STD = c(0, 5, 20, 80, 250, 1000)
# Absorbance = c( 1.268967137,
#  0.736537869,
#  	0.297820701,
#  	0.283878975,
#  	0.228719425,
#  	0.131721101)
# 
# STD = data.frame(STD,Absorbance )

df0 <- tibble(
  Exposure = c(10, 25, 50, 75, 100, 150, 300, 400),
  Response = c(0.03, 0.04, 0.15, 0.12, 0.25, 0.43, 0.54, 0.47)
)
library(DoseFinding)
  fitMod(Exposure, Response, data = df0,  model = "emax")
  fitMod(Exposure, Response, data = df0,  model = "sigEmax")
  
  #e0     0.0594     0.0324
  ## eMax   0.4515     0.0535
  ## ed50 107.0215    11.5670
  ## h      4.1370     1.6398
  
library(tidyverse)
library(tidyr)
library(dplyr)
library(readr)
library(ggplot2)
library(ggpubr)
ggplot(data = df0, aes(x = Exposure, y = Response)) +
  labs(title = "Quantifying PGD2 in cell culture lysates and its enzymatic reactions ",
       caption = "PGD2 ELISA")+
  geom_point(colour = "#69b3a2")+
  stat_smooth(method= "nls", 
              formula = y~e0 + eMax*x^h/(ed50^h+x^h),
              method.args = list( start = c(e0 = 0.0594, eMax = 0.4515, ed50 = 107.0215, h = 4.1370)),
              se=FALSE)+
  xlab(expression(paste("%B/"~B[0])))+
  ylab(expression(paste("Prostaglandin"~ D[2], ~~ " MOX Concentration (pg/ml) ")))  +
  
  theme(plot.background =  element_rect(fill = "transparent"),
        panel.background = element_blank(),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black"))+
  
  theme(legend.spacing.y = unit(0.01, "cm"))+
  theme(legend.position = c(0.77, .91),
        legend.background = element_rect(colour = NA, fill = NA))+
  theme(plot.title = element_text(size = 12, face = "bold.italic"),
        plot.caption = element_text(hjust = 0))
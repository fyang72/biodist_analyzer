# file.path(R.home("bin"), "R")
# "/opt/scp/apps/gen-2021a/software/R/4.3.2-foss-2021a/lib64/R/bin/R"



source("renv/activate.R")
  
#addl.libpath = c( "/opt/scp/services/azcore/biometrics/0.3.1/libraries/R/4.3.1"  )
addl.libpath = c("C:/Users/feng.yang/Documents/FYANG/biodist_analyzer/biodist_analyzer/renv/library/R-4.3/x86_64-w64-mingw32"                 
                 )

.libPaths(unique(c(  addl.libpath, .libPaths() )) )

#library(azcore) 
library(renv)
library(utils) 

#renv::settings$ignored.packages("azcore")

# don't use this: azcore_module_load("CMake/3.20.1-GCCcore-8.2.0")  # for 4.1.0
#azcore::azcore_module_load("CMake/3.20.1-GCCcore-10.3.0")   # for 4.3.1.

# environment restoration example
# if you need to restore environment from renv.lock then use next three lines
#settings$ignored.packages("azcore")
# renv::restore()

#  environment initialization example
#library(renv)
#renv::init(settings = list(ignored.packages=c("azcore")))

# renv::snapshot()

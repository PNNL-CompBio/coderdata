install.packages('BiocManager')
install.packages("remotes")
BiocManager::install(version = "3.16")
remotes::install_github("RGLab/cytolib")
BiocManager::install("cmapR",update=TRUE,ask=FALSE)
install.packages("readr")
install.packages("tidyr")
install.packages("dplyr")


#install.packages("readr")
#install.packages("stringr")
#install.packages("dplyr")
install.packages('reticulate')
                                        #install.packages('tidyr')
install.packages('tidyverse')
install.packages("BiocManager")
BiocManager::install("PharmacoGx",update=TRUE,ask=FALSE)


# Fail the Docker build layer if any required package didn't install.
# install.packages() and BiocManager::install() only WARN on failure, so a
# transient CRAN/Bioconductor outage during `docker build` produced an image
# that built "successfully" but was missing packages -- the failure then
# surfaced hours into a build run instead of at build time.
required_pkgs <- c("reticulate", "tidyverse", "BiocManager", "PharmacoGx")
missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
if (length(missing_pkgs) > 0) stop(paste("Required packages failed to install:", paste(missing_pkgs, collapse = ", ")))

install.packages('curl')
install.packages("rio")
install.packages("readr")
install.packages("curl")
install.packages("stringr")
install.packages("dplyr")
install.packages("XML")
#install.packages('reticulate')
install.packages('tidyr')
install.packages('httr2')
install.packages('data.table')

# Fail the Docker build layer if any required package didn't install.
# install.packages() and BiocManager::install() only WARN on failure, so a
# transient CRAN/Bioconductor outage during `docker build` produced an image
# that built "successfully" but was missing packages -- the failure then
# surfaced hours into a build run instead of at build time.
required_pkgs <- c("curl", "rio", "readr", "stringr", "dplyr", "XML", "tidyr", "httr2", "data.table")
missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
if (length(missing_pkgs) > 0) stop(paste("Required packages failed to install:", paste(missing_pkgs, collapse = ", ")))

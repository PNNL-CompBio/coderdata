install.packages('reticulate', repos='https://cloud.r-project.org')
reticulate::use_virtualenv('/opt/venv', required = TRUE)
install.packages('remotes')
remotes::install_version('rjson', version = '0.2.21', repos = 'https://cloud.r-project.org')
install.packages('synapser', repos = c('http://ran.synapse.org', 'https://cloud.r-project.org'))
install.packages("dplyr")
install.packages("synapser", repos = c("http://ran.synapse.org", "https://cloud.r-project.org"))
install.packages("data.table")
install.packages("R.utils")
install.packages("stringr")
install.packages("tidyr")
install.packages("readr")
install.packages("readxl")

# Fail the Docker build layer if any required package didn't install.
# install.packages() and BiocManager::install() only WARN on failure, so a
# transient CRAN/Bioconductor outage during `docker build` produced an image
# that built "successfully" but was missing packages -- the failure then
# surfaced hours into a build run instead of at build time.
required_pkgs <- c("reticulate", "remotes", "rjson", "synapser", "dplyr", "data.table", "R.utils", "stringr", "tidyr", "readr", "readxl")
missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
if (length(missing_pkgs) > 0) stop(paste("Required packages failed to install:", paste(missing_pkgs, collapse = ", ")))

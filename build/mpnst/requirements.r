# Install the rjson package at version 0.2.21
devtools::install_version("rjson", version = "0.2.21")

# Install other required packages
install.packages("dplyr")
install.packages("data.table")
install.packages("synapser", repos = c("http://ran.synapse.org", "https://cloud.r-project.org"))
install.packages("R.utils")
install.packages("stringr")
install.packages("reticulate")
install.packages("tidyr")
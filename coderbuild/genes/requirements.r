install.packages("BiocManager")
BiocManager::install(c("AnnotationDbi", "org.Hs.eg.db", "biomaRt"), update=FALSE, ask=FALSE)
install.packages('https://cran.r-project.org/src/contrib/Archive/dbplyr/dbplyr_2.3.4.tar.gz', repos = NULL)
install.packages('tidyr')

# Fail the Docker build layer if any required package didn't install
required_pkgs <- c("org.Hs.eg.db", "biomaRt", "tidyr")
missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
if (length(missing_pkgs) > 0) stop(paste("Required packages failed to install:", paste(missing_pkgs, collapse = ", ")))

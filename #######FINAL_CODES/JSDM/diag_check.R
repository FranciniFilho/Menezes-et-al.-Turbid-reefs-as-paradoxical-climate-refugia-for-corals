
# Diagnostic script
cat("Checking packages...\n")
ok <- TRUE
for(pkg in c("readxl", "janitor", "digest", "dplyr")) {
  if(!require(pkg, character.only=TRUE)) {
    cat("Missing:", pkg, "\n")
    ok <- FALSE
  } else {
    cat("Loaded:", pkg, "\n")
    packageVersion(pkg)
  }
}

cat("\nChecking paths...\n")
path <- "C:/Users/rbfra/OneDrive/########PUBLICACOES/############Menezes et al. Mus his distribution and abundance Abrolhos/######22.04.23/DATA/BENTHOS TEMPORAL CLEAN_v3.xlsx"
if(file.exists(path)) {
  cat("Benthos file exists!\n")
} else {
  cat("Benthos file NOT FOUND!\n")
}

dataset_paths <- list(
  CV_02 = "C:/Users/rbfra/OneDrive/########PUBLICACOES/############Menezes et al. Mus his distribution and abundance Abrolhos/######FINAL/#######FINAL_RESULTS/#####output_local_PCA_CV_2_FINAL/dados_abundancia_integrados_long_format.csv",
  CV_30 = "C:/Users/rbfra/OneDrive/########PUBLICACOES/############Menezes et al. Mus his distribution and abundance Abrolhos/######FINAL/#######FINAL_RESULTS/#####output_local_PCA_CV_30_FINAL/dados_abundancia_integrados_long_format.csv",
  CV_ALL = "C:/Users/rbfra/OneDrive/########PUBLICACOES/############Menezes et al. Mus his distribution and abundance Abrolhos/######FINAL/#######FINAL_RESULTS/#####output_local_PCA_CV_all_FINAL/dados_abundancia_integrados_long_format.csv"
)

for(n in names(dataset_paths)) {
  if(file.exists(dataset_paths[[n]])) {
    cat(n, "exists.\n")
  } else {
    cat(n, "NOT FOUND.\n")
  }
}

# analysis_summary.R
# Parses result directories to identify best models.
# Handles both BRT (CSV summary with , or ;) and Bayesian (LOO txt) formats.

# Paths
result_dirs <- c(
    "c:/Users/rbfra/OneDrive/########PUBLICACOES/############Menezes et al. Mus his distribution and abundance Abrolhos/#################FINAL/#######FINAL_RESULTS/BRT_Health_Growth_CV_2",
    "c:/Users/rbfra/OneDrive/########PUBLICACOES/############Menezes et al. Mus his distribution and abundance Abrolhos/#################FINAL/#######FINAL_RESULTS/BRT_Health_Growth_CV_30",
    "c:/Users/rbfra/OneDrive/########PUBLICACOES/############Menezes et al. Mus his distribution and abundance Abrolhos/#################FINAL/#######FINAL_RESULTS/Bayesian_Analyses_Health_Growth_CV_all"
)

report_lines <- c("# Model Analysis Report", "", paste("Generated on", Sys.time()), "")

# Function to process BRT CSV results
process_brt_csv <- function(csv_path) {
    cat("Processing CSV:", csv_path, "\n")
    # Try reading with comma
    df <- tryCatch(read.csv(csv_path, stringsAsFactors = FALSE), error = function(e) NULL)

    # If failed or only 1 column, try semicolon
    if (is.null(df) || ncol(df) <= 1) {
        df <- tryCatch(read.csv(csv_path, sep = ";", stringsAsFactors = FALSE), error = function(e) NULL)
    }

    if (is.null(df)) {
        cat("  Failed to read CSV.\n")
        return(NULL)
    }

    # Check columns
    req_cols <- c("response_variable", "model_name", "brt_cv_deviance")
    if (!all(req_cols %in% names(df))) {
        cat("  Missing columns. Found:", names(df), "\n")
        return(NULL)
    }

    # Filter valid deviance
    df <- df[!is.na(df$brt_cv_deviance), ]
    if (nrow(df) == 0) {
        return(NULL)
    }

    # Split by response variable
    responses <- unique(df$response_variable)
    lines <- c()

    for (resp in responses) {
        sub_df <- df[df$response_variable == resp, ]
        # Find min deviance
        best_idx <- which.min(sub_df$brt_cv_deviance)
        best_row <- sub_df[best_idx, ]

        lines <- c(lines, paste0(
            "- **Response:** ", resp,
            ", **Best Model:** ", best_row$model_name,
            ", **Deviance:** ", best_row$brt_cv_deviance
        ))
    }
    lines
}

# Function to process Bayesian TXT results
process_bayesian_txt <- function(txt_path) {
    # Read table handling row names (header is shifted)
    tbl <- tryCatch(read.table(txt_path, header = TRUE, stringsAsFactors = FALSE), error = function(e) NULL)
    if (is.null(tbl) || nrow(tbl) == 0) {
        return(NULL)
    }

    # The best model is the first row
    best_model <- rownames(tbl)[1]

    # Extract response name from path (parent directory)
    resp_name <- basename(dirname(txt_path))

    paste0("- **Response:** ", resp_name, ", **Best Model:** ", best_model)
}

for (res_dir in result_dirs) {
    cat("Checking directory:", res_dir, "\n")
    report_lines <- c(report_lines, paste0("## Results in `", res_dir, "`"), "")

    # Check for BRT CSV
    csv_file <- file.path(res_dir, "FINAL_MODEL_PERFORMANCE_COMPARISON.csv")
    if (file.exists(csv_file)) {
        lines <- process_brt_csv(csv_file)
        if (!is.null(lines)) report_lines <- c(report_lines, lines)
    } else {
        # Assume Bayesian structure
        outputs_dir <- file.path(res_dir, "Outputs")
        if (dir.exists(outputs_dir)) {
            scenarios <- list.dirs(outputs_dir, recursive = FALSE, full.names = TRUE)
            for (scen in scenarios) {
                scen_name <- basename(scen)
                report_lines <- c(report_lines, paste0("### Scenario: `", scen_name, "`"), "")

                comp_files <- list.files(scen, pattern = "model_comparison_loo_.*\\.txt$", full.names = TRUE, recursive = TRUE)
                for (f in comp_files) {
                    line <- process_bayesian_txt(f)
                    if (!is.null(line)) report_lines <- c(report_lines, line)
                }
                report_lines <- c(report_lines, "")
            }
        }
    }
    report_lines <- c(report_lines, "---", "")
}

report_path <- "c:/Users/rbfra/OneDrive/########PUBLICACOES/############Menezes et al. Mus his distribution and abundance Abrolhos/#################FINAL/#######FINAL_CODES/model_analysis_report.md"
writeLines(report_lines, con = report_path)
cat("Report written to:", report_path, "\n")

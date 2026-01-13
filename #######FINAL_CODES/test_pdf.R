library(ggplot2)
p <- ggplot(mtcars, aes(wt, mpg)) +
    geom_point() +
    labs(title = "Test PDF")
dir <- "c:/Users/rbfra/OneDrive/########PUBLICACOES/############Menezes et al. Mus his distribution and abundance Abrolhos/######FINAL/#######FINAL_RESULTS/PCA_Publication_Figures"
dir.create(dir, showWarnings = FALSE, recursive = TRUE)
ggsave(file.path(dir, "test_output.pdf"), p)
cat("PDF saved successfully\n")

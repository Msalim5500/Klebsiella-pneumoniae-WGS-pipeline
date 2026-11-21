#!/usr/bin/env Rscript
# Generic ggalluvial template. Edit stage1/stage2/stage3 in metadata before use.

suppressPackageStartupMessages({
  library(ggplot2)
  library(ggalluvial)
})

args <- commandArgs(trailingOnly = TRUE)
infile <- if (length(args) >= 1) args[[1]] else "metadata/sankey_metadata_template.tsv"
outfile <- if (length(args) >= 2) args[[2]] else "results/figures/Kp_Sankey_diagram.pdf"

dat <- read.delim(infile, stringsAsFactors = FALSE, check.names = FALSE)
required <- c("isolate_id", "stage1", "stage2", "stage3", "value")
missing <- setdiff(required, names(dat))
if (length(missing)) stop("Missing columns: ", paste(missing, collapse = ", "))
if (any(grepl("REPLACE", as.matrix(dat[, c("stage1", "stage2", "stage3")])))) {
  stop("Replace placeholder categories in the input metadata before plotting.")
}

dat$stage1 <- factor(dat$stage1)
dat$stage2 <- factor(dat$stage2)
dat$stage3 <- factor(dat$stage3)
dat$value <- as.numeric(dat$value)
if (anyNA(dat$value)) stop("The value column must contain numeric values.")

p <- ggplot(dat,
            aes(axis1 = stage1, axis2 = stage2, axis3 = stage3, y = value)) +
  geom_alluvium(aes(fill = stage1), width = 1/12, alpha = 0.8) +
  geom_stratum(width = 1/8, fill = "grey90", colour = "grey35") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3) +
  scale_x_discrete(limits = c("Stage 1", "Stage 2", "Stage 3"), expand = c(.05, .05)) +
  labs(x = NULL, y = "Number of isolates", fill = "Stage 1") +
  theme_minimal(base_size = 11) +
  theme(panel.grid = element_blank())

dir.create(dirname(outfile), recursive = TRUE, showWarnings = FALSE)
ggsave(outfile, p, width = 9, height = 6, device = "pdf")
message("Saved: ", outfile)

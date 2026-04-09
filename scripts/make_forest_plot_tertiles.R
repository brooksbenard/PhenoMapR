#!/usr/bin/env Rscript
#
# Build a forest plot using upper vs lower tertiles of PhenoMapR scores.
# Uses the existing consolidated panel file produced by:
#   scripts/run_tcga_precog_forest_incremental.R
#
# Output:
# - results/tcga_precog_forest_tertiles_results.tsv
# - results/tcga_precog_forest_tertiles.png

suppressPackageStartupMessages({
  library(data.table)
  library(survival)
  library(ggplot2)
})

panel_path <- "results/tcga_precog_forest_panel.tsv"
out_results <- "results/tcga_precog_forest_tertiles_results.tsv"
out_png <- "results/tcga_precog_forest_tertiles.png"

if (!file.exists(panel_path)) {
  stop("Missing panel file: ", panel_path)
}

dt <- fread(panel_path)
req <- c("forest_label", "tcga_code", "stratum", "precog_label_used", "os_time", "os_event", "score")
miss <- setdiff(req, names(dt))
if (length(miss) > 0) {
  stop("Panel file missing required columns: ", paste(miss, collapse = ", "))
}

# Match prior plot behavior: exclude Normal strata from forest plot outputs.
dt <- dt[stratum != "Normal"]

# Patient-level aggregation (avoid counting multiple samples per patient).
dt <- dt[
  is.finite(score) & !is.na(score) & is.finite(os_time) & !is.na(os_event),
  .(
    score = mean(score, na.rm = TRUE),
    os_time = os_time[1],
    os_event = os_event[1],
    tcga_code = tcga_code[1],
    stratum = stratum[1],
    precog_label_used = precog_label_used[1],
    forest_label = forest_label[1]
  ),
  by = .(forest_label, patient)
]

res_list <- list()

labels <- sort(unique(dt$forest_label))
for (lab in labels) {
  sub <- dt[forest_label == lab]
  if (nrow(sub) < 40) next

  q1 <- as.numeric(stats::quantile(sub$score, probs = 1/3, na.rm = TRUE, type = 7))
  q2 <- as.numeric(stats::quantile(sub$score, probs = 2/3, na.rm = TRUE, type = 7))
  if (!is.finite(q1) || !is.finite(q2) || q1 >= q2) next

  sub[, tertile := fifelse(score <= q1, "Low", fifelse(score >= q2, "High", "Mid"))]
  sub <- sub[tertile %in% c("Low", "High")]
  if (nrow(sub) < 20) next

  sub[, tertile := factor(tertile, levels = c("Low", "High"))]

  fit <- survival::coxph(survival::Surv(os_time, os_event) ~ tertile, data = sub)
  hr <- exp(coef(fit))[["tertileHigh"]]
  ci <- exp(confint(fit))[1, ]

  lr <- survival::survdiff(survival::Surv(os_time, os_event) ~ tertile, data = sub)
  p_lr <- 1 - stats::pchisq(lr$chisq, 1)

  res_list[[lab]] <- data.table(
    forest_label = lab,
    tcga_code = unique(sub$tcga_code)[1],
    stratum = unique(sub$stratum)[1],
    precog_label_used = unique(sub$precog_label_used)[1],
    n_samples = nrow(sub),
    n_patients = uniqueN(sub$patient),
    q1 = q1,
    q2 = q2,
    hr = as.numeric(hr),
    ci_low = as.numeric(ci[1]),
    ci_high = as.numeric(ci[2]),
    p_logrank = as.numeric(p_lr)
  )
}

if (length(res_list) == 0) {
  stop("No strata had enough samples for tertile analysis.")
}

res <- rbindlist(res_list, use.names = TRUE, fill = TRUE)
res <- res[is.finite(hr) & is.finite(ci_low) & is.finite(ci_high)]
res <- res[order(hr, decreasing = TRUE)]
res[, sig := ifelse(is.finite(p_logrank) & p_logrank < 0.05, "Significant", "Not significant")]
res$sig <- factor(res$sig, levels = c("Significant", "Not significant"))

fwrite(res, out_results, sep = "\t")

p <- ggplot(res, aes(y = reorder(forest_label, hr), x = hr, color = sig)) +
  geom_vline(xintercept = 1, linetype = "dashed", linewidth = 0.4) +
  geom_errorbar(aes(xmin = ci_low, xmax = ci_high), height = 0.15, linewidth = 0.6) +
  geom_point(size = 2.4) +
  scale_x_log10() +
  scale_color_manual(
    values = c("Significant" = "#B2182B", "Not significant" = "#4D4D4D"),
    name = "Log-rank"
  ) +
  labs(
    title = "TCGA survival by PhenoMapR score (tertile extremes)",
    subtitle = "Upper vs lower tertile (middle tertile excluded); Normal excluded",
    x = "Hazard ratio (High vs Low tertile; log scale)",
    y = "TCGA cancer type / stratum"
  ) +
  theme_minimal(base_size = 12)

ggsave(out_png, p, width = 7.5, height = 4.5, dpi = 200)

cat("Wrote:", out_results, "\n")
cat("Wrote:", out_png, "\n")


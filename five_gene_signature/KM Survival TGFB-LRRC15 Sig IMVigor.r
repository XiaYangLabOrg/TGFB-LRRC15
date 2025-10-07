###################################################################################################################
# LOAD DATA 
###################################################################################################################

library(survival)
library(survminer)
library(ggplot2)
library(dplyr)
# Load files (counts + patient data)
expression_data <- read.csv("./data/IMvigor210_counts_pkg.csv", header = TRUE, row.names = 1, check.names = FALSE)

# are these files different? I am getting different results
# remove NA from binary column. filter to excluded patients in immune phenotype column
# patient_data <- read.csv("./data/excluded.csv", header = TRUE, row.names = 1, check.names = FALSE)
patient_data <- read.csv("./data/IMvigor210_patient.csv", header = TRUE, row.names = 1, check.names = FALSE)
patient_data <- patient_data[!is.na(patient_data$binaryResponse),]
patient_data <- patient_data[patient_data[,"Immune phenotype"] == "excluded",]

# Gene names to uppercase
rownames(expression_data) <- toupper(rownames(expression_data))

###################################################################################################################
# ADD LRRC15 EXPRESSION SORTING (High = top 75%)
###################################################################################################################

lrrc15_expression <- as.numeric(as.character(expression_data["LRRC15", ]))
lrrc15q25 <- unname(quantile(lrrc15_expression, 0.25, na.rm = TRUE))
lrrc15_group <- ifelse(lrrc15_expression >= lrrc15q25, "High", "Low")
patient_data$lrrc15_group <- lrrc15_group[match(rownames(patient_data), colnames(expression_data))]

###################################################################################################################
# SET SIGNATURE
###################################################################################################################

signature_genes <- c("TGFBR2","LRRC15", "WNT5B","MMP2", "SPARC")

###################################################################################################################
# ALIGN (NO extra filtering here; excluded.csv is already filtered)
###################################################################################################################

# Keep all samples present in both tables
common_ids <- intersect(colnames(expression_data), rownames(patient_data))
expression_data_aligned <- expression_data[, common_ids, drop = FALSE]
patient_data_aligned    <- patient_data[common_ids, , drop = FALSE]

expression_data_immunefiltered <- expression_data_aligned
patient_data_immunefiltered    <- patient_data_aligned

###################################################################################################################
# NORMALIZE + WEIGHT GENES
###################################################################################################################

# Keep only signature genes that are present
present <- intersect(toupper(signature_genes), rownames(expression_data_immunefiltered))
if (length(present) == 0) 
  {stop("No genes found")
}

# Filter matrix
expression_mat <- expression_data_immunefiltered[present, , drop = FALSE]

# Transformation
expression_mat <- log1p(expression_mat)  # log(1 + counts)
varmat <- apply(expression_mat, 1, sd, na.rm = TRUE) > 0
expression_mat <- expression_mat[varmat, , drop = FALSE]

expression_mat <- expression_mat[, colSums(!is.na(expression_mat)) > 0, drop = FALSE]
patient_data_immunefiltered <- patient_data_immunefiltered[colnames(expression_mat), , drop = FALSE]

# Z-score per gene across patients
z_scored_data <- t(scale(t(expression_mat)))

# Score (handles single-gene or multi-gene)
if (nrow(z_scored_data) == 1L) {
  gene_weights <- 1
  weighted_scores <- as.numeric(z_scored_data[1, ])
  names(weighted_scores) <- colnames(z_scored_data)
  }  else {
  pca_result <- prcomp(t(z_scored_data))
  gene_weights <- pca_result$rotation[, 1]
  weighted_scores <- colMeans(z_scored_data * gene_weights)
  }

###################################################################################################################
# PLOT KM
###################################################################################################################

# Survival data (expects columns os, censOS)
patient_data_immunefiltered$os <- as.numeric(patient_data_immunefiltered$os)
patient_data_immunefiltered$censOS <- as.integer(patient_data_immunefiltered$censOS)

survival_data <- data.frame(
  time_to_event = patient_data_immunefiltered$os,
  event_status = patient_data_immunefiltered$censOS,
  weighted_scores = weighted_scores
)

survival_data <- survival_data[complete.cases(survival_data), , drop = FALSE]

med <- median(survival_data$weighted_scores, na.rm = TRUE)
survival_data$LRRC15_sig <- ifelse(survival_data$weighted_scores > med, "High", "Low")

print(table(survival_data$LRRC15_sig, useNA = "ifany"))

# Cox model
cox_model <- coxph(Surv(time_to_event, event_status) ~ weighted_scores, data = survival_data)
print(summary(cox_model))

# KM
km_fit <- survfit(Surv(time_to_event, event_status) ~ LRRC15_sig, data = survival_data)

km_plot <- ggsurvplot(
  km_fit, 
  data = survival_data,
  pval = TRUE,
  conf.int = FALSE,
  risk.table = TRUE,
  pval.method = FALSE,
  ggtheme = theme_minimal(base_size = 15),
  risk.table.y.text.col = TRUE,
  risk.table.height = 0.35,
  risk.table.y.text = TRUE,
  palette = c("orangered","#4387E7"),
  title = "Kaplan–Meier Survival Curve",
  xlab = "Time (months)",
  ylab = "Survival Probability",
  pval.size = 6,
  legend.title = "LRRC15 Signature",
  legend.labs = c("High","Low"),
  font.x = c(16,"bold"),
  font.y = c(16,"bold"),
  font.legend = c(16),
  font.tickslab = c(14),
  risk.table.fontsize = 5,
  linetype = "strata",
  break.time.by = 6,
  censor.shape = 124,
  censor.size = 7
)

print(km_plot)

ggsave("./results/KM_IMvigor210.png", plot = km_plot$plot, width = 7, height = 5, dpi = 300)

###################################################################################################################
# PERMUTATION
###################################################################################################################
set.seed(1234)                 
N_random <- 1000              
sig_size <- 5

# KM log-rank p 
km_p_from_scores <- function(scores, surv_df) {
  df <- data.frame(
    time_to_event = surv_df$time_to_event,
    event_status  = surv_df$event_status,
    score         = as.numeric(scores)
  )
  df <- df[complete.cases(df), , drop = FALSE]
  if (nrow(df) < 2) return(NA_real_)
  med <- median(df$score, na.rm = TRUE)
  grp <- ifelse(df$score > med, "High", "Low")
  if (length(unique(grp)) < 2) {
    q25 <- unname(quantile(df$score, 0.25, na.rm = TRUE))
    grp <- ifelse(df$score >= q25, "High", "Low")
  }
  if (length(unique(grp)) < 2) return(NA_real_)
  sdiff <- survdiff(Surv(time_to_event, event_status) ~ grp, data = df)
  pchisq <- 1 - pchisq(sdiff$chisq, length(sdiff$n) - 1)
  as.numeric(pchisq)
}

# Log-rank chi-square and p-value (median split with 25% fallback)
km_logrank_from_scores <- function(scores, surv_df) {
  df <- data.frame(
    time_to_event = surv_df$time_to_event,
    event_status  = surv_df$event_status,
    score         = as.numeric(scores)
  )
  df <- df[complete.cases(df), , drop = FALSE]
  if (nrow(df) < 2) return(list(chisq = NA_real_, p = NA_real_))
  med <- median(df$score, na.rm = TRUE)
  grp <- ifelse(df$score > med, "High", "Low")
  if (length(unique(grp)) < 2) {
    q25 <- unname(quantile(df$score, 0.25, na.rm = TRUE))
    grp <- ifelse(df$score >= q25, "High", "Low")
  }
  if (length(unique(grp)) < 2) return(list(chisq = NA_real_, p = NA_real_))
  sdiff <- survdiff(Surv(time_to_event, event_status) ~ grp, data = df)
  chisq <- as.numeric(sdiff$chisq)
  p <- 1 - pchisq(chisq, length(sdiff$n) - 1)
  list(chisq = chisq, p = as.numeric(p))
}

 

# Scoring Genes
score_gene_set <- function(gset, expr_mat) {
  X <- expr_mat[gset, , drop = FALSE]
  X <- log1p(X)
  keep <- apply(X, 1, sd, na.rm = TRUE) > 0
  X <- X[keep, , drop = FALSE]
  if (nrow(X) == 0) return(rep(NA_real_, ncol(expr_mat)))
  Z <- t(scale(t(X)))
  if (nrow(Z) == 1L) {
    as.numeric(Z[1, ])
  } else {
    pc <- prcomp(t(Z))
    w  <- pc$rotation[, 1]
    as.numeric(colMeans(Z * w))
  }
}

# Genes present in matrix, not in signature
gene_universe <- rownames(expression_data_immunefiltered)
sig_now <- intersect(toupper(signature_genes), gene_universe)
stopifnot(length(sig_now) >= 1)

obs_scores <- score_gene_set(sig_now, expression_data_immunefiltered)
obs_logrank <- km_logrank_from_scores(obs_scores, survival_data)
cat("5-gene signature log-rank: chisq = ", signif(obs_logrank$chisq, 4),
    "; p = ", signif(obs_logrank$p, 4), "\n", sep = "")

gene_pool <- setdiff(gene_universe, sig_now)
gene_pool <- gene_pool[!is.na(gene_pool)]

# Permutations (log-rank chi-square/p-value)
rand_logrank_chisq <- numeric(N_random)
rand_p <- numeric(N_random)
rand_genes <- character(N_random)

pb <- txtProgressBar(min = 0, max = N_random, style = 3)
for (i in seq_len(N_random)) {
  gset <- sample(gene_pool, size = sig_size, replace = FALSE)
  sc   <- score_gene_set(gset, expression_data_immunefiltered)
  lrk <- km_logrank_from_scores(sc, survival_data)
  rand_logrank_chisq[i] <- lrk$chisq
  rand_p[i] <- lrk$p
  rand_genes[i] <- paste(gset, collapse = ";")
  setTxtProgressBar(pb, i)
}
close(pb)

# Empirical p-values (log-rank)
valid_lrk <- is.finite(rand_logrank_chisq)
emp_p_lrk <- mean(rand_logrank_chisq[valid_lrk] >= obs_logrank$chisq)
valid_p <- is.finite(rand_p)
emp_p_by_p <- mean(rand_p[valid_p] <= obs_logrank$p)

cat("Empirical p (log-rank chisq) = ", signif(emp_p_lrk, 4),
    "  using ", sum(valid_lrk), " valid permutations.\n", sep = "")
cat("Empirical p (by p-value) = ", signif(emp_p_by_p, 4),
    "  using ", sum(valid_p), " valid permutations.\n", sep = "")

# Save permutation summary (observed stats + empirical p-values)
perm_summary <- data.frame(
  observed_logrank_chisq = obs_logrank$chisq,
  observed_logrank_p = obs_logrank$p,
  empirical_p_logrank_chisq = emp_p_lrk,
  empirical_p_logrank_p = emp_p_by_p
)
write.csv(perm_summary, "./results/IMvigor_permutation_summary.csv", row.names = FALSE)

# Plot
perm_df <- data.frame(logrank_chisq_random = rand_logrank_chisq,
                      logrank_p_random = rand_p,
                      genes = rand_genes)
write.csv(perm_df, "./results/IMvigor_permutation.csv", row.names = FALSE)

plt2 <- ggplot(perm_df[valid_lrk, ], aes(x = logrank_chisq_random)) +
  geom_histogram(bins = 40) +
  geom_vline(xintercept = obs_logrank$chisq, linetype = 2, color = "red") +
  labs(
    title = "Random 5-gene signature (log-rank chi-square)",
    subtitle = paste0("Observed chi-square = ", signif(obs_logrank$chisq, 3),
                      " | Empirical p = ", signif(emp_p_lrk, 3)),
    x = "Log-rank chi-square", y = "Count"
  ) +
  theme_minimal(base_size = 14)

print(plt2)
ggsave("./results/IMvigor_permutation_hist_logrank.png", plt2, width = 7, height = 5, dpi = 600)

best_idx_p <- order(rand_p)[seq_len(min(10, sum(valid_p)))]
cat("\nTop random sets (smallest log-rank p):\n")
print(data.frame(rank = seq_along(best_idx_p), p = rand_p[best_idx_p], genes = rand_genes[best_idx_p]))

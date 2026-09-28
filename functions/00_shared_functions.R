###############################################################################
# 阶段0·共享管线 Step3: 共享绘图函数库
# 项目: HNSC三方案并行生信分析
# 日期: 2026-07-20
# 说明: 三方案共用的绘图函数，source此文件即可使用
#       包含: KM曲线、timeROC、森林图、热图、风险评分图、火山图等
#
# Positron版: 无硬编码路径，纯函数库
###############################################################################

# =============================================================================
# 1. Kaplan-Meier 生存曲线（支持自定义分组+中位/最佳截断）
# =============================================================================
plot_km <- function(surv_data, risk_data, gene = NULL, title = "", 
                    cutpoint = "median", output_file = NULL, width = 6, height = 5) {
  library(survival)
  library(survminer)
  
  if (!is.null(gene)) {
    expr_values <- risk_data[[gene]]
    if (cutpoint == "median") {
      group <- ifelse(expr_values >= median(expr_values, na.rm = TRUE), "High", "Low")
    } else if (cutpoint == "optimal") {
      df_temp <- data.frame(time = surv_data$OS_time, event = surv_data$OS_status, expr = expr_values)
      res_cut <- surv_cutpoint(df_temp, time = "time", event = "event", variables = "expr")
      group <- surv_categorize(res_cut)$expr
      group <- ifelse(group == "high", "High", "Low")
    }
    surv_data$group <- factor(group, levels = c("Low", "High"))
    title <- ifelse(title == "", paste0(gene, " expression"), title)
  } else {
    surv_data$group <- factor(risk_data$group, levels = c("Low", "High"))
  }
  
  fit <- survfit(Surv(OS_time, OS_status) ~ group, data = surv_data)
  
  p <- ggsurvplot(
    fit, data = surv_data,
    pval = TRUE, pval.method = TRUE,
    conf.int = TRUE,
    risk.table = TRUE,
    legend.title = "Group",
    legend.labs = c("Low", "High"),
    palette = c("#4DBBD5", "#E64B35"),
    title = title,
    xlab = "Time (days)", ylab = "Overall Survival Probability",
    risk.table.height = 0.25,
    ggtheme = theme_bw()
  )
  
  if (!is.null(output_file)) {
    pdf(output_file, width = width, height = height)
    print(p)
    dev.off()
  }
  return(p)
}

# =============================================================================
# 2. Time-dependent ROC曲线
# =============================================================================
plot_time_roc <- function(surv_data, risk_scores, times = c(1, 3, 5), 
                          title = "", output_file = NULL, width = 6, height = 6) {
  library(timeROC)
  library(ggplot2)
  
  # 时间单位转换（天->年）
  times_days <- times * 365
  
  roc_res <- timeROC(
    T = surv_data$OS_time,
    delta = surv_data$OS_status,
    marker = risk_scores,
    cause = 1,
    times = times_days,
    ROC = TRUE
  )
  
  # 整理绘图数据
  plot_df <- data.frame()
  for (i in seq_along(times)) {
    plot_df <- rbind(plot_df, data.frame(
      FPR = 1 - roc_res$TP[, i],
      TPR = roc_res$TP[, i],
      Year = paste0(times[i], "-year"),
      AUC = round(roc_res$AUC[i], 3)
    ))
  }
  plot_df$Year <- factor(plot_df$Year, levels = paste0(times, "-year"))
  
  p <- ggplot(plot_df, aes(FPR, TPR, color = Year)) +
    geom_line(linewidth = 1) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey50") +
    annotate("text", x = 0.6, y = 0.1, 
             label = paste0("AUC:\n", 
                           paste(unique(paste0(plot_df$Year, ": ", plot_df$AUC)), collapse = "\n")),
             size = 3.5, hjust = 0) +
    scale_color_manual(values = c("#E64B35", "#4DBBD5", "#00A087")) +
    labs(title = title, x = "1 - Specificity", y = "Sensitivity") +
    theme_bw() +
    coord_equal()
  
  if (!is.null(output_file)) {
    ggsave(output_file, p, width = width, height = height)
  }
  return(list(plot = p, roc = roc_res))
}

# =============================================================================
# 3. 森林图（单因素/多因素Cox回归）
# =============================================================================
plot_forest <- function(cox_results, title = "", output_file = NULL, width = 8, height = 5) {
  library(ggplot2)
  
  if (!"Variable" %in% colnames(cox_results)) {
    if ("Gene" %in% colnames(cox_results)) {
      cox_results$Variable <- cox_results$Gene
    } else {
      cox_results$Variable <- rownames(cox_results)
    }
  }
  
  cox_results$sig <- ifelse(cox_results$p < 0.05, "sig", "ns")
  
  p <- ggplot(cox_results, aes(x = HR, y = Variable)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "grey50") +
    geom_errorbarh(aes(xmin = lower, xmax = upper, color = sig), height = 0.2) +
    geom_point(aes(color = sig), size = 3) +
    geom_text(aes(label = sprintf("%.2f (%.2f-%.2f), %s", HR, lower, upper, 
                                  ifelse(p < 0.001, "<0.001", round(p, 3)))),
              hjust = -0.1, size = 3) +
    scale_color_manual(values = c("sig" = "#E64B35", "ns" = "grey50")) +
    scale_x_log10() +
    labs(title = title, x = "Hazard Ratio (log scale)", y = "") +
    theme_bw() +
    theme(legend.position = "none",
          axis.text.y = element_text(size = 10))
  
  if (!is.null(output_file)) {
    ggsave(output_file, p, width = width, height = height)
  }
  return(p)
}

# =============================================================================
# 4. 风险评分分布图
# =============================================================================
plot_risk_distribution <- function(surv_data, risk_scores, output_file = NULL, width = 10, height = 8) {
  library(ggplot2)
  library(survminer)
  library(patchwork)
  
  df_temp <- data.frame(time = surv_data$OS_time, event = surv_data$OS_status, risk = risk_scores)
  res_cut <- surv_cutpoint(df_temp, time = "time", event = "event", variables = "risk")
  df_temp$group <- surv_categorize(res_cut)$risk
  df_temp$group <- ifelse(df_temp$group == "high", "High Risk", "Low Risk")
  df_temp <- df_temp[order(df_temp$risk), ]
  df_temp$rank <- seq_len(nrow(df_temp))
  df_status <- df_temp[df_temp$event == 1, ]
  
  p1 <- ggplot(df_temp, aes(rank, risk, color = group)) +
    geom_point(size = 1) +
    scale_color_manual(values = c("Low Risk" = "#4DBBD5", "High Risk" = "#E64B35")) +
    labs(x = "Patient (sorted by risk score)", y = "Risk Score") +
    theme_bw() + theme(legend.position = "none", axis.text.x = element_blank(), axis.ticks.x = element_blank())
  
  p2 <- ggplot(df_status, aes(rank)) +
    geom_segment(aes(xend = rank, y = 0, yend = 1), color = "#E64B35") +
    labs(x = "", y = "OS Event") +
    theme_bw() + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank()) +
    ylim(0, 1)
  
  p <- p1 / p2
  
  if (!is.null(output_file)) {
    ggsave(output_file, p, width = width, height = height)
  }
  return(p)
}

# =============================================================================
# 5. 热图
# =============================================================================
plot_heatmap <- function(expr_matrix, group_info, genes = NULL, 
                         title = "", output_file = NULL, width = 10, height = 8) {
  library(pheatmap)
  
  if (!is.null(genes)) {
    expr_matrix <- expr_matrix[rownames(expr_matrix) %in% genes, ]
  }
  expr_z <- t(scale(t(expr_matrix)))
  annotation_col <- data.frame(Group = group_info$group, row.names = group_info$sample)
  ann_colors <- list(Group = c("Low" = "#4DBBD5", "High" = "#E64B35", "Normal" = "#00A087", "Tumor" = "#E64B35"))
  
  p <- pheatmap(expr_z, annotation_col = annotation_col, annotation_colors = ann_colors,
                show_colnames = FALSE, cluster_cols = TRUE, cluster_rows = TRUE, main = title,
                color = colorRampPalette(c("#3C5488", "white", "#E64B35"))(100))
  
  if (!is.null(output_file)) {
    pdf(output_file, width = width, height = height)
    print(p)
    dev.off()
  }
  return(p)
}

# =============================================================================
# 6. ssGSEA评分
# =============================================================================
calc_ssgsea <- function(expr_matrix, gene_sets) {
  library(GSVA)
  gsva_param <- ssgseaParam(exprData = as.matrix(expr_matrix), geneSets = gene_sets, normalize = TRUE)
  ssgsea_scores <- gsva(gsva_param)
  return(ssgsea_scores)
}

# =============================================================================
# 7. 单基因Cox回归（批量）
# =============================================================================
batch_uni_cox <- function(expr_matrix, surv_data, genes = NULL) {
  library(survival)
  if (is.null(genes)) genes <- rownames(expr_matrix)
  
  common_samples <- intersect(colnames(expr_matrix), surv_data$sample)
  expr_matrix <- expr_matrix[, common_samples]
  surv_data <- surv_data[match(common_samples, surv_data$sample), ]
  
  results <- lapply(genes, function(g) {
    expr <- as.numeric(expr_matrix[g, ])
    if (all(is.na(expr)) || sd(expr, na.rm = TRUE) == 0) return(NULL)
    fit <- coxph(Surv(OS_time, OS_status) ~ expr, data = surv_data)
    s <- summary(fit)
    data.frame(Gene = g, HR = exp(coef(fit)), lower = s$conf.int[, "lower .95"],
               upper = s$conf.int[, "upper .95"], p = s$coefficients[, "Pr(>|z|)"], stringsAsFactors = FALSE)
  })
  results <- do.call(rbind, results)
  results$FDR <- p.adjust(results$p, method = "BH")
  return(results)
}

# =============================================================================
# 8. LASSO-Cox预后模型构建
# =============================================================================
build_lasso_cox <- function(expr_matrix, surv_data, genes, alpha = 1, nfolds = 10,
                            seed = 20260720, foldid = NULL) {
  library(glmnet)
  library(survival)

  common_samples <- intersect(colnames(expr_matrix), surv_data$sample)
  expr_matrix <- expr_matrix[genes, common_samples]
  surv_data <- surv_data[match(common_samples, surv_data$sample), ]

  x <- t(as.matrix(expr_matrix))
  y <- Surv(surv_data$OS_time, surv_data$OS_status)

  # Fix fold assignment for full reproducibility across R/glmnet versions.
  # Without foldid, cv.glmnet uses sample() internally to assign folds,
  # which can produce different lambda.min values across environments even
  # with the same set.seed(), causing HR variation.
  if (is.null(foldid)) {
    set.seed(seed)
    foldid <- sample(rep(1:nfolds, length.out = nrow(x)))
  }
  cat("    Using fixed foldid (seed=", seed, ") for cv.glmnet\n", sep = "")

  cv_fit <- cv.glmnet(x, y, family = "cox", alpha = alpha, nfolds = nfolds,
                      foldid = foldid, cox.ties = "efron")
  cat("    lambda.min:", signif(cv_fit$lambda.min, 5),
      " lambda.1se:", signif(cv_fit$lambda.1se, 5), "\n")

  coef_lasso <- coef(cv_fit, s = "lambda.min")
  selected_genes <- rownames(coef_lasso)[as.numeric(coef_lasso) != 0]

  cat("    LASSO选中基因数:", length(selected_genes), "\n")
  cat("    基因:", paste(selected_genes, collapse = ", "), "\n")

  risk_scores <- as.numeric(x[, selected_genes] %*% coef_lasso[selected_genes, ])
  return(list(cv_fit = cv_fit, coef = coef_lasso, selected_genes = selected_genes,
              risk_scores = risk_scores, samples = common_samples, foldid = foldid))
}

# =============================================================================
# 9. 工具函数：表达矩阵样本对齐
# =============================================================================
align_samples <- function(expr_matrix, surv_data) {
  common <- intersect(colnames(expr_matrix), surv_data$sample)
  expr_matrix <- expr_matrix[, common]
  surv_data <- surv_data[match(common, surv_data$sample), ]
  return(list(expr = expr_matrix, surv = surv_data))
}

cat("共享绘图函数库加载完成!\n")
cat("可用函数: plot_km, plot_time_roc, plot_forest, plot_risk_distribution,\n")
cat("          plot_heatmap, calc_ssgsea, batch_uni_cox, build_lasso_cox, align_samples\n")

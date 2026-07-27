#!/usr/bin/env Rscript
library(metafor)

OUTFILE <- "/gpfs/data/user/shreyags/TE_work/replication_results/new_enrichment/pop_replication_all_nosingletons/TB_vs_pooled_nonTB_metafor_results.txt"
sink(OUTFILE, split = TRUE)

run_TB_vs_pooled_from_OR <- function(df, label) {
    cat(label, "- TB vs. pooled non-TB (from OR/CI)\n")
    if (!"TB" %in% df$Dataset) {
        cat("No 'TB' dataset found skipping.\n")
        return(invisible(NULL))
    }
    # log-OR and variance, calculated from the CI width
    df$yi  <- log(df$OR)
    df$sei <- (log(df$ci_upper) - log(df$ci_lower)) / (2 * 1.96)
    df$vi  <- df$sei^2
    tb_row <- df[df$Dataset == "TB", ]
    nonTB  <- df[df$Dataset != "TB", ]
    # Pool the non-TB populations into one fixed-effect estimate
    nonTB_fit <- rma(yi, vi, data = nonTB, method = "FE")
    pooled_row <- data.frame(
        Dataset = "PooledNonTB",
        yi      = as.numeric(nonTB_fit$b),
        vi      = as.numeric(nonTB_fit$se)^2
    )
    two_strata <- rbind(
        data.frame(Dataset = "TB", yi = tb_row$yi, vi = tb_row$vi),
        pooled_row
    )
    cat("\n--- Log-OR and variance per group ---\n")
    print(two_strata)
    fit <- rma(yi, vi, data = two_strata, method = "FE", slab = two_strata$Dataset)
    cat("\n--- Cochran's Q: TB vs pooled non-TB ---\n")
    print(fit)
    cat("\nOR scale:\n")
    print(data.frame(
        estimate = exp(as.numeric(fit$b)),
        ci.lb    = exp(fit$ci.lb),
        ci.ub    = exp(fit$ci.ub)
    ))
    invisible(fit)
}

#results from "/gpfs/data/user/shreyags/TE_work/replication_results/new_enrichment/pop_replication_all_nosingletons/"
erd_or <- data.frame(
    Dataset  = c("AA", "CAO", "DR", "IE", "TB"),
    OR       = c(0.934947897314429,
                 0.830023400147042,
                 0.968393283682371,
                 0.909102586803235,
                 1.022064306530340),
    ci_lower = c(0.723746985848231,
                 0.731693006552587,
                 0.861304700096378,
                 0.839909108211924,
                 0.825871530992370),
    ci_upper = c(1.200801982340210,
                 0.940045488127009,
                 1.087536617681620,
                 0.983449703555177,
                 1.260350510297230)
)

lrd_or <- data.frame(
    Dataset  = c("AA", "CAO", "DR", "IE", "TB"),
    OR       = c(1.353272976125690,
                 1.749329436100030,
                 1.404265623453000,
                 1.450719147709470,
                 1.220822972630920),
    ci_lower = c(1.066383363513970,
                 1.559361378432190,
                 1.257730341104190,
                 1.347607090712660,
                 0.996360873880348),
    ci_upper = c(1.715603635430380,
                 1.962949344023190,
                 1.567672345599850,
                 1.561658830718050,
                 1.494018592616820)
)

run_TB_vs_pooled_from_OR(erd_or, "ERD")
run_TB_vs_pooled_from_OR(lrd_or, "LRD")

sink()
cat("Results written to:", OUTFILE, "\n")
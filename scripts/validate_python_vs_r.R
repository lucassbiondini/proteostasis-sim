suppressPackageStartupMessages({
  library(survival)
  library(cmprsk)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
csv_path <- if (length(args) >= 1) args[1] else "data/processed/crossval_cohort.csv"
json_path <- if (length(args) >= 2) args[2] else "data/processed/python_reference_metrics.json"
tol <- if (length(args) >= 3) as.numeric(args[3]) else 1e-2

df <- read.csv(csv_path)
py <- fromJSON(json_path)

t_ref <- as.numeric(py$t_ref)

ci <- cuminc(ftime = df$time, fstatus = df$event, group = df$arm, cencode = 0)
ci_t <- timepoints(ci, times = t_ref)$est

r_cif_g0 <- as.numeric(ci_t["0 1", 1])
r_cif_a0 <- as.numeric(ci_t["0 2", 1])
r_cif_g1 <- as.numeric(ci_t["1 1", 1])
r_cif_a1 <- as.numeric(ci_t["1 2", 1])

fit_csh <- coxph(Surv(time, event == 1) ~ arm, data = df, ties = "breslow")
r_hr_csh <- as.numeric(exp(coef(fit_csh)["arm"]))

fit_fg <- crr(ftime = df$time, fstatus = df$event, cov1 = as.matrix(df$arm), failcode = 1, cencode = 0)
r_shr_fg <- as.numeric(exp(fit_fg$coef[1]))

comp <- data.frame(
  Metric = c("CIF_gen_arm0", "CIF_amb_arm0", "CIF_gen_arm1", "CIF_amb_arm1", "HR_CSH_gen", "sHR_FineGray_gen"),
  Python = c(py$cif_gen_arm0, py$cif_amb_arm0, py$cif_gen_arm1, py$cif_amb_arm1, py$hr_csh_gen, py$shr_fg_gen),
  R_GoldStandard = c(r_cif_g0, r_cif_a0, r_cif_g1, r_cif_a1, r_hr_csh, r_shr_fg)
)
comp$AbsDiff <- abs(comp$Python - comp$R_GoldStandard)
comp$RelDiff <- comp$AbsDiff / pmax(abs(comp$R_GoldStandard), 1e-6)
comp$Pass <- comp$AbsDiff < tol | comp$RelDiff < tol

print(comp, row.names = FALSE)

out_report <- list(
  passed = all(comp$Pass),
  tolerance = tol,
  max_abs_diff = max(comp$AbsDiff),
  max_rel_diff = max(comp$RelDiff),
  comparison = comp
)
write_json(out_report, "data/processed/r_crossval_report.json", pretty = TRUE, auto_unbox = TRUE)

if (!all(comp$Pass)) {
  stop("CROSS-VALIDATION FAILED: Python vs R divergence exceeded tolerance!")
}
cat("CROSS-VALIDATION PASSED: Python exact estimators match R (survival + cmprsk) within tolerance.\n")

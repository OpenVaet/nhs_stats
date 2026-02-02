#!/usr/bin/env Rscript
# ------------------------------------------------------------
# NHS staff vaccination model (daily administered dose 1/2/3)
# + Sensitivity analysis on MEDIAN days between doses (via delay kernel)
#
# Core model:
# - N = 1,480,000 staff
# - Dose 2: minimum 21 days after dose 1, with extra-delay ~ Normal(mean, sd)
# - Dose 3: minimum 90 days after dose 2, with extra-delay ~ Normal(mean, sd)
# - Anchors are matched (piecewise scaling). Small mismatches can remain where anchors
#   decrease month-to-month because daily administered doses cannot be negative.
#
# Sensitivity:
# - Runs multiple scenarios varying the delay-kernel median (extra-delay median),
#   while keeping sd fixed (unless you change it).
# - Produces per-scenario monthly dose CSVs + summary tables comparing to baseline.
#
# OUTPUTS:
# - Baseline: data/nhs_staff_vax_model_monthly.csv
# - Sensitivity monthly CSVs: data/sensitivity/nhs_staff_vax_model_monthly__<tag>.csv
# - Sensitivity summary: data/sensitivity/vax_delay_sensitivity_summary.csv
# - Sensitivity diffs vs baseline: data/sensitivity/vax_delay_sensitivity_diff_vs_baseline.csv
# - Tidy long for overlay plotting: data/sensitivity/vax_delay_sensitivity_monthly_long.csv
# ------------------------------------------------------------

suppressPackageStartupMessages({
  library(stats)
  library(utils)
})

N <- 1480000

start_date <- as.Date("2020-12-20")
end_date   <- as.Date("2022-08-31")
dates <- seq(start_date, end_date, by = "day")

# -------------------- Helpers --------------------

# POSIXlt wday: 0=Sun,1=Mon,...,6=Sat
weekday_factor <- function(d) {
  w <- as.POSIXlt(d)$wday
  #        Sun  Mon  Tue  Wed  Thu  Fri  Sat
  fac <- c(0.55,1.15,1.15,1.10,1.10,1.00,0.60)
  fac[w + 1]
}

# Build discrete delay kernel over extra days 0..max_extra with weights ~ Normal(mean, sd).
# Returns weights and the (discrete) median extra-day where CDF crosses 0.5.
make_delay_kernel <- function(mean_delay, sd_delay, max_extra) {
  extra <- 0:max_extra
  w <- dnorm(extra, mean = mean_delay, sd = sd_delay)
  w[!is.finite(w)] <- 0
  w <- pmax(0, w)
  if (sum(w) <= 0) stop("Kernel has zero mass: check mean/sd.")
  w <- w / sum(w)
  cdf <- cumsum(w)
  med_extra <- extra[which(cdf >= 0.5)[1]]
  list(extra = extra, w = w, med_extra = med_extra)
}

# Convolve a daily series with a delayed kernel and a hard minimum gap:
# output[t + min_gap + extra] += daily_in[t] * w[extra]
convolve_with_min_gap <- function(daily_in, min_gap, mean_delay, sd_delay, max_extra = 180) {
  k <- make_delay_kernel(mean_delay, sd_delay, max_extra)
  extra <- k$extra
  w <- k$w

  out <- rep(0, length(daily_in))
  for (t in seq_along(daily_in)) {
    contrib_days <- t + min_gap + extra
    ok <- contrib_days <= length(out)
    if (any(ok)) {
      out[contrib_days[ok]] <- out[contrib_days[ok]] + daily_in[t] * w[ok]
    }
  }
  out
}

# Piecewise scale a daily series so cumulative exactly hits (monotone) anchors
# while keeping realistic weekday seasonality (weekends lower).
scale_daily_to_anchors <- function(dates, daily_raw, anchors, N) {
  anchors <- anchors[order(anchors$date), ]
  anchors$coverage <- pmax(0, pmin(1, anchors$coverage))
  anchors$coverage_mono <- cummax(anchors$coverage)  # cannot decrease if daily>=0

  idx <- match(anchors$date, dates)
  if (any(is.na(idx))) stop("One or more anchor dates are outside modeled date range.")

  daily <- pmax(0, daily_raw)

  # Scale each interval so the total doses in that interval matches the anchor delta
  for (i in seq_len(nrow(anchors) - 1)) {
    i0 <- idx[i]
    i1 <- idx[i + 1]
    if (i1 <= i0) next

    seg <- (i0 + 1):i1

    target_delta <- (anchors$coverage_mono[i + 1] - anchors$coverage_mono[i]) * N
    target_delta <- max(0, target_delta)

    shaped <- daily[seg] * weekday_factor(dates[seg])
    ssum <- sum(shaped)

    if (ssum > 0) {
      daily[seg] <- shaped * (target_delta / ssum)
    } else {
      wfac <- weekday_factor(dates[seg])
      daily[seg] <- (wfac / sum(wfac)) * target_delta
    }
  }

  cum <- cumsum(pmax(0, daily))
  list(daily = pmax(0, daily), cum = cum, cov = cum / N, anchors = anchors)
}

# Build dose 1: monotone spline through anchors -> daily increments -> weekday-shaped scaling
build_dose1_from_anchors <- function(dates, anchors1, N) {
  anchors1 <- anchors1[order(anchors1$date), ]
  anchors1$coverage <- pmax(0, pmin(1, anchors1$coverage))
  anchors1$coverage_mono <- cummax(anchors1$coverage)

  t0 <- min(dates)
  ta <- as.numeric(anchors1$date - t0)
  td <- as.numeric(dates - t0)

  f <- splinefun(x = ta, y = anchors1$coverage_mono, method = "hyman")
  cov_raw <- pmax(0, pmin(1, f(td)))
  cov <- cummax(cov_raw)
  cum <- cov * N

  daily_base <- c(cum[1], diff(cum))
  daily_base <- pmax(0, daily_base)

  scale_daily_to_anchors(dates = dates, daily_raw = daily_base, anchors = anchors1, N = N)
}

verify_anchors <- function(label, anchor_dates, anchor_cov, cum_series, dates, N) {
  idx <- match(anchor_dates, dates)
  if (any(is.na(idx))) stop(paste(label, "verification: anchor date out of range."))

  model_cov <- cum_series[idx] / N
  err_pp <- (model_cov - anchor_cov) * 100
  err_people <- cum_series[idx] - (anchor_cov * N)

  data.frame(
    series = label,
    date = anchor_dates,
    anchor_pct = round(anchor_cov * 100, 3),
    model_pct  = round(model_cov * 100, 3),
    err_pct_pt = round(err_pp, 3),
    err_people = as.integer(round(err_people))
  )
}

# -------------------- Anchors --------------------
anchors_1 <- data.frame(
  date = as.Date(c(
    "2020-12-20",
    "2020-12-31",
    "2021-01-31",
    "2021-02-28",
    "2021-03-07",  # early March
    "2021-05-16",
    "2021-08-31",
    "2021-09-30",
    "2021-10-31",
    "2021-11-30",
    "2021-12-31",
    "2022-01-31",
    "2022-02-28",
    "2022-03-31",
    "2022-04-30",
    "2022-05-31",
    "2022-06-30",
    "2022-07-31",
    "2022-08-31"
  )),
  coverage = c(
    0.000,
    0.050,  # assumption: early ramp by end-Dec 2020
    0.600,
    0.760,
    0.805,
    0.878,
    0.918,
    0.924,
    0.929,
    0.935,
    0.943,
    0.953,
    0.954,
    0.957,
    0.957,
    0.957,
    0.956,
    0.956,
    0.955
  )
)

anchors_2_official <- data.frame(
  date = as.Date(c(
    "2020-12-20",
    "2021-08-31",
    "2021-09-30",
    "2021-10-31",
    "2021-11-30",
    "2021-12-31",
    "2022-01-31",
    "2022-02-28",
    "2022-03-31",
    "2022-04-30",
    "2022-05-31",
    "2022-06-30",
    "2022-07-31",
    "2022-08-31"
  )),
  coverage = c(
    0.000,
    0.878,
    0.890,
    0.899,
    0.907,
    0.915,
    0.920,
    0.923,
    0.930,
    0.931,
    0.932,
    0.931,
    0.932,
    0.929
  )
)

booster_start <- as.Date("2021-09-14")  # mid-Sep launch (assumption)
anchors_3_official <- data.frame(
  date = as.Date(c(
    booster_start,
    "2021-11-30",
    "2021-12-31",
    "2022-01-31",
    "2022-02-28",
    "2022-03-31",
    "2022-04-30",
    "2022-05-31",
    "2022-06-30",
    "2022-07-31",
    "2022-08-31"
  )),
  coverage = c(
    0.000,
    0.609,
    0.744,
    0.773,
    0.781,
    0.793,
    0.797,
    0.801,
    0.802,
    0.804,
    0.802
  )
)

# -------------------- Build dose 1 (fixed across sensitivity) --------------------
m1 <- build_dose1_from_anchors(dates, anchors_1, N)
daily1 <- m1$daily
cum1   <- m1$cum

# Fixed minimum gaps:
min_gap_2 <- 21
min_gap_3 <- 90

# -------------------- Core runner (dose2+3 depend on delay params) --------------------
run_scenario <- function(tag, mean2, sd2, mean3, sd3) {

  # Dose 2 raw from dose 1
  daily2_raw <- convolve_with_min_gap(
    daily_in = daily1,
    min_gap = min_gap_2,
    mean_delay = mean2,
    sd_delay = sd2,
    max_extra = 180
  )

  m2 <- scale_daily_to_anchors(dates, daily2_raw, anchors_2_official, N)

  # Enforce logical constraint: cum2 <= cum1
  cum2 <- pmin(m2$cum, cum1)
  cum2 <- cummax(cum2)
  daily2 <- c(cum2[1], diff(cum2))
  daily2 <- pmax(0, daily2)

  # Dose 3 raw from dose 2
  daily3_raw <- convolve_with_min_gap(
    daily_in = daily2,
    min_gap = min_gap_3,
    mean_delay = mean3,
    sd_delay = sd3,
    max_extra = 240
  )
  daily3_raw[dates < booster_start] <- 0

  m3 <- scale_daily_to_anchors(dates, daily3_raw, anchors_3_official, N)

  # Enforce logical constraint: cum3 <= cum2
  cum3 <- pmin(m3$cum, cum2)
  cum3 <- cummax(cum3)
  daily3 <- c(cum3[1], diff(cum3))
  daily3 <- pmax(0, daily3)

  # After you build constrained daily3 from cum3:
  m3b <- scale_daily_to_anchors(
    dates   = dates,
    daily_raw = daily3,              # use the constrained daily3 as the shape
    anchors = anchors_3_official,
    N = N
  )

  cum3b <- pmin(m3b$cum, cum2)
  cum3b <- cummax(cum3b)
  daily3 <- c(cum3b[1], diff(cum3b))
  daily3 <- pmax(0, daily3)
  cum3 <- cum3b

  # Monthly output
  month_key <- format(dates, "%Y-%m")
  monthly <- aggregate(
    cbind(first_dose = daily1, second_dose = daily2, third_dose = daily3),
    by = list(month = month_key),
    FUN = sum
  )

  monthly$first_dose  <- as.integer(round(monthly$first_dose))
  monthly$second_dose <- as.integer(round(monthly$second_dose))
  monthly$third_dose  <- as.integer(round(monthly$third_dose))

  # Kernel medians (discrete, extra-delay only)
  k2 <- make_delay_kernel(mean2, sd2, 180)
  k3 <- make_delay_kernel(mean3, sd3, 240)

  # Peaks
  p2 <- which.max(monthly$second_dose)
  p3 <- which.max(monthly$third_dose)

  diag <- data.frame(
    tag = tag,
    mean2 = mean2, sd2 = sd2,
    median_extra2 = k2$med_extra,
    median_total2 = min_gap_2 + k2$med_extra,
    mean3 = mean3, sd3 = sd3,
    median_extra3 = k3$med_extra,
    median_total3 = min_gap_3 + k3$med_extra,
    peak2_month = monthly$month[p2],
    peak2_doses = monthly$second_dose[p2],
    peak3_month = monthly$month[p3],
    peak3_doses = monthly$third_dose[p3]
  )

  list(
    tag = tag,
    monthly = monthly,
    cum2 = cum2,
    cum3 = cum3,
    diag = diag
  )
}

# -------------------- Scenarios --------------------
# Interpret these as "median extra-delay" targets (Normal median ≈ mean; discretization/truncation may shift by ~1 day).
# Total median between doses = min_gap + median_extra.
scenarios <- data.frame(
  tag = c("baseline_m60", "m45", "m75"),
  mean2 = c(60, 45, 75),
  sd2   = c(20, 20, 20),
  mean3 = c(60, 45, 75),
  sd3   = c(20, 20, 20),
  stringsAsFactors = FALSE
)

# -------------------- Run all scenarios --------------------
runs <- lapply(seq_len(nrow(scenarios)), function(i) {
  s <- scenarios[i, ]
  run_scenario(s$tag, s$mean2, s$sd2, s$mean3, s$sd3)
})

diag_tbl <- do.call(rbind, lapply(runs, function(x) x$diag))

# Peak month difference vs baseline (in months)
baseline_tag <- "baseline_m60"
baseline_diag <- diag_tbl[diag_tbl$tag == baseline_tag, ]

month_to_num <- function(ym) as.integer(substr(ym,1,4))*12L + as.integer(substr(ym,6,7))
b <- diag_tbl[diag_tbl$tag=="baseline_m60", ]
diag_tbl$peak2_shift_mo <- month_to_num(diag_tbl$peak2_month) - month_to_num(b$peak2_month)
diag_tbl$peak3_shift_mo <- month_to_num(diag_tbl$peak3_month) - month_to_num(b$peak3_month)
print(diag_tbl[,c("tag","median_total2","peak2_month","peak2_shift_mo","median_total3","peak3_month","peak3_shift_mo")], row.names=FALSE)

# -------------------- Write baseline output (for downstream plots) --------------------
baseline_idx <- which(diag_tbl$tag == "baseline_m60")[1]
baseline_monthly <- runs[[baseline_idx]]$monthly

dir.create("data", showWarnings = FALSE, recursive = TRUE)
write.csv(baseline_monthly, "data/nhs_staff_vax_model_monthly.csv", row.names = FALSE)
cat("Wrote baseline: data/nhs_staff_vax_model_monthly.csv\n")

# -------------------- Sensitivity outputs --------------------
dir.create("data/sensitivity", showWarnings = FALSE, recursive = TRUE)

for (r in runs) {
  out_csv <- sprintf("data/sensitivity/nhs_staff_vax_model_monthly__%s.csv", r$tag)
  write.csv(r$monthly, out_csv, row.names = FALSE)
  cat("Wrote: ", out_csv, "\n", sep = "")
}

write.csv(diag_tbl, "data/sensitivity/vax_delay_sensitivity_summary.csv", row.names = FALSE)
cat("Wrote: data/sensitivity/vax_delay_sensitivity_summary.csv\n")

# -------------------- Diff vs baseline (monthly) --------------------
baseline <- baseline_monthly[order(baseline_monthly$month), ]
diff_tbl <- do.call(rbind, lapply(runs, function(r) {
  m <- r$monthly[order(r$monthly$month), ]
  stopifnot(all(m$month == baseline$month))

  data.frame(
    tag = r$tag,
    max_abs_diff_second = max(abs(m$second_dose - baseline$second_dose)),
    max_abs_diff_third  = max(abs(m$third_dose  - baseline$third_dose)),
    sum_abs_diff_second = sum(abs(m$second_dose - baseline$second_dose)),
    sum_abs_diff_third  = sum(abs(m$third_dose  - baseline$third_dose))
  )
}))

write.csv(diff_tbl, "data/sensitivity/vax_delay_sensitivity_diff_vs_baseline.csv", row.names = FALSE)
cat("Wrote: data/sensitivity/vax_delay_sensitivity_diff_vs_baseline.csv\n")

# -------------------- Tidy long output for overlay plotting --------------------
monthly_long <- do.call(rbind, lapply(runs, function(r) {
  m <- r$monthly
  data.frame(
    tag = r$tag,
    month = m$month,
    first_dose  = m$first_dose,
    second_dose = m$second_dose,
    third_dose  = m$third_dose,
    stringsAsFactors = FALSE
  )
}))

write.csv(monthly_long, "data/sensitivity/vax_delay_sensitivity_monthly_long.csv", row.names = FALSE)
cat("Wrote: data/sensitivity/vax_delay_sensitivity_monthly_long.csv\n\n")

# -------------------- Verification vs anchors (baseline only, as before) --------------------
# Dose 1: verify only "given" anchors (exclude the two start assumptions)
official_1_idx <- anchors_1$date >= as.Date("2021-01-31")
v1 <- verify_anchors("dose1",
  anchors_1$date[official_1_idx],
  anchors_1$coverage[official_1_idx],
  cum1, dates, N
)

# Dose 2: verify official anchors from 2021-08-31 onward (exclude the 0 start anchor)
official_2_idx <- anchors_2_official$date >= as.Date("2021-08-31")
v2 <- verify_anchors("dose2",
  anchors_2_official$date[official_2_idx],
  anchors_2_official$coverage[official_2_idx],
  runs[[baseline_idx]]$cum2, dates, N
)

# Dose 3: verify official anchors from 2021-11-30 onward (exclude booster_start 0 anchor)
official_3_idx <- anchors_3_official$date >= as.Date("2021-11-30")
v3 <- verify_anchors("dose3",
  anchors_3_official$date[official_3_idx],
  anchors_3_official$coverage[official_3_idx],
  runs[[baseline_idx]]$cum3, dates, N
)

verif <- rbind(v1, v2, v3)

cat("Verification vs anchors (baseline; percentage-point error and people):\n")
print(verif, row.names = FALSE)

cat("\nMax absolute error (pct-pt) by series (baseline):\n")
print(aggregate(abs(verif$err_pct_pt), by = list(series = verif$series), FUN = max))

cat("\nSensitivity summary (medians are TOTAL days incl. min gaps):\n")
print(diag_tbl, row.names = FALSE)

cat("\nDiff vs baseline (monthly dose counts):\n")
print(diff_tbl, row.names = FALSE)

cat("\nNotes:\n",
    "- Sensitivity varies the *extra-delay* Normal kernel mean (≈ median), and reports the discrete kernel median.\n",
    "- Total median between doses = min_gap + median_extra.\n",
    "- Dose 2 is generated from Dose 1 with a hard minimum gap of 21 days.\n",
    "- Dose 3 is generated from Dose 2 with a hard minimum gap of 90 days and forced to 0 before ",
    format(booster_start), ".\n",
    "- If an anchor decreases month-to-month, the model cannot follow the decrease (administered doses can't be negative),\n",
    "  so it matches a monotone version of those anchors.\n",
    sep = ""
)

cat("\nFirst 10 months of BASELINE modeled administered doses:\n")
print(head(baseline_monthly, 10), row.names = FALSE)

# Verify dose2/dose3 anchors for ALL scenarios (and summarize max abs error)
verif_all <- do.call(rbind, lapply(seq_along(runs), function(i) {
  tag <- runs[[i]]$tag

  # dose2 official anchors
  official_2_idx <- anchors_2_official$date >= as.Date("2021-08-31")
  v2 <- verify_anchors(
    paste0("dose2__", tag),
    anchors_2_official$date[official_2_idx],
    anchors_2_official$coverage[official_2_idx],
    runs[[i]]$cum2, dates, N
  )

  # dose3 official anchors
  official_3_idx <- anchors_3_official$date >= as.Date("2021-11-30")
  v3 <- verify_anchors(
    paste0("dose3__", tag),
    anchors_3_official$date[official_3_idx],
    anchors_3_official$coverage[official_3_idx],
    runs[[i]]$cum3, dates, N
  )

  rbind(v2, v3)
}))

# max abs error per series+scenario
verif_all$scenario <- sub(".*__", "", verif_all$series)
verif_all$series0  <- sub("__.*", "", verif_all$series)

max_err <- aggregate(abs(verif_all$err_pct_pt),
                     by = list(series = verif_all$series0, scenario = verif_all$scenario),
                     FUN = max)
names(max_err)[3] <- "max_abs_err_pctpt"
print(max_err)


# Example: correlate scenario dose totals vs sickness total (replace 'sickness_series' with yours)
# df_sick must have columns: month (YYYY-MM) and sickness_rate_100k (numeric)
# Read it from your exported file if you want:
# df_sick <- read.csv("data/sickness_roll6_prediction_band_2020_2026.csv") |> ...

library(dplyr)

df_sick <- read.csv("data/sickness_roll6_prediction_band_2020_2026.csv", stringsAsFactors = FALSE) |>
  as_tibble() |>
  transmute(
    month = month,
    sickness = as.numeric(sickness_rate_100k)
  )

corr_by_scenario <- do.call(rbind, lapply(runs, function(r) {
  m <- r$monthly |> as_tibble() |> transmute(month = month, doses = second_dose + third_dose)
  d <- inner_join(m, df_sick, by = "month") |> filter(is.finite(doses), is.finite(sickness))
  data.frame(
    tag = r$tag,
    r = if (nrow(d) >= 3) cor(d$doses, d$sickness) else NA_real_,
    n = nrow(d)
  )
}))

print(corr_by_scenario)

# after runs are computed
verif_all <- do.call(rbind, lapply(runs, function(r) {
  official_2_idx <- anchors_2_official$date >= as.Date("2021-08-31")
  v2 <- verify_anchors(paste0("dose2__", r$tag),
                       anchors_2_official$date[official_2_idx],
                       anchors_2_official$coverage[official_2_idx],
                       r$cum2, dates, N)
  official_3_idx <- anchors_3_official$date >= as.Date("2021-11-30")
  v3 <- verify_anchors(paste0("dose3__", r$tag),
                       anchors_3_official$date[official_3_idx],
                       anchors_3_official$coverage[official_3_idx],
                       r$cum3, dates, N)
  rbind(v2, v3)
}))
verif_all$scenario <- sub(".*__", "", verif_all$series)
verif_all$series0  <- sub("__.*", "", verif_all$series)

max_err <- aggregate(abs(verif_all$err_pct_pt),
                     by = list(series = verif_all$series0, scenario = verif_all$scenario),
                     FUN = max)
names(max_err)[3] <- "max_abs_err_pctpt"
print(max_err)

cap_rate <- function(cum_a, cum_b) mean(cum_a >= cum_b - 1e-9)
for (r in runs) {
  cat(r$tag, " cap(cum2==cum1):", round(cap_rate(r$cum2, cum1), 3),
      " cap(cum3==cum2):", round(cap_rate(r$cum3, r$cum2), 3), "\n")
}


corr_window <- function(df, start_ym, end_ym) {
  start <- as.Date(paste0(start_ym, "-01"))
  end   <- as.Date(paste0(end_ym, "-01"))
  df |> filter(month_date >= start, month_date <= end)
}

# build df_sick with month_date to avoid string compare
df_sick2 <- df_sick %>%
  transmute(
    month = month,
    sickness = sickness
  )

corr_by_scenario_window <- do.call(rbind, lapply(runs, function(r) {
  m <- as_tibble(r$monthly) %>%
    transmute(
      month = month,
      month_date = as.Date(paste0(month, "-01")),
      doses = second_dose + third_dose
    )

  d <- inner_join(m, df_sick2, by = "month") %>%
    filter(is.finite(doses), is.finite(sickness)) %>%
    filter(month_date >= as.Date("2020-05-01"), month_date <= as.Date("2022-04-01"))

  data.frame(tag = r$tag, r = cor(d$doses, d$sickness), n = nrow(d))
}))

print(corr_by_scenario_window)
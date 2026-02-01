# ------------------------------------------------------------
# NHS staff vaccination model (daily administered dose 1/2/3)
# - N = 1,480,000 staff
# - Dose 2: minimum 21 days after dose 1, with extra-delay ~ Normal(mean=60, sd=20)
# - Dose 3: minimum 90 days after dose 2, with extra-delay ~ Normal(mean=60, sd=20)
# - Anchors are matched (piecewise scaling). Small mismatches can remain where anchors
#   decrease month-to-month because daily administered doses cannot be negative.
#
# OUTPUT:
#   nhs_staff_vax_model_monthly.csv with:
#     month,first_dose,second_dose,third_dose
#
# VERIFICATION:
#   Prints model vs anchors (pct-pt error and people).
# ------------------------------------------------------------

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

# Convolve a daily series with a delayed kernel:
# output[t + min_gap + extra] += daily_in[t] * w[extra]
convolve_with_min_gap <- function(daily_in, min_gap, mean_delay, sd_delay, max_extra = 180) {
  extra <- 0:max_extra
  w <- dnorm(extra, mean = mean_delay, sd = sd_delay)
  w[extra < 0] <- 0
  w <- w / sum(w)

  out <- rep(0, length(daily_in))
  for (t in seq_along(daily_in)) {
    contrib_days <- t + min_gap + extra
    ok <- contrib_days <= length(out)
    if (any(ok)) out[contrib_days[ok]] <- out[contrib_days[ok]] + daily_in[t] * w[ok]
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

  daily <- daily_raw

  # Ensure no negatives
  daily <- pmax(0, daily)

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
      # if no mass, distribute purely by weekday weights
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

  # Now scale per-anchor-interval to exactly match anchor deltas, adding weekday seasonality
  scaled <- scale_daily_to_anchors(
    dates = dates,
    daily_raw = daily_base,
    anchors = anchors1,
    N = N
  )
  scaled
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
# Interpretation choices (explicit):
# - "in January / February" interpreted as end-of-month.
# - "early March" mapped to 2021-03-07.
# - End-Dec 2020: add a small starter anchor for campaign start.

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

# Dose 2 official anchors (end Aug 2021 onward). Add a start anchor at 0.
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

# Dose 3 (booster) official anchors (end Nov 2021 onward).
# Add booster campaign start at 0 to prevent any boosters before launch.
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

# -------------------- Build dose 1 --------------------

m1 <- build_dose1_from_anchors(dates, anchors_1, N)
daily1 <- m1$daily
cum1   <- m1$cum

# -------------------- Build dose 2 (>=21 days after dose 1) --------------------

min_gap_2  <- 21
mean_delay <- 60
sd_delay   <- 20

daily2_raw <- convolve_with_min_gap(
  daily_in = daily1,
  min_gap = min_gap_2,
  mean_delay = mean_delay,
  sd_delay = sd_delay,
  max_extra = 180
)

m2 <- scale_daily_to_anchors(
  dates = dates,
  daily_raw = daily2_raw,
  anchors = anchors_2_official,
  N = N
)

# Enforce logical constraint: cum2 <= cum1 (cannot have more 2nd doses than 1st doses)
cum2 <- pmin(m2$cum, cum1)
cum2 <- cummax(cum2)
daily2 <- c(cum2[1], diff(cum2))
daily2 <- pmax(0, daily2)

# -------------------- Build dose 3 (>=90 days after dose 2; boosters start mid-Sep 2021) --------------------

min_gap_3 <- 90

daily3_raw <- convolve_with_min_gap(
  daily_in = daily2,
  min_gap = min_gap_3,
  mean_delay = mean_delay,  # reusing 60/20 unless you want different values
  sd_delay = sd_delay,
  max_extra = 240
)

# Zero out any boosters before campaign launch date
daily3_raw[dates < booster_start] <- 0

m3 <- scale_daily_to_anchors(
  dates = dates,
  daily_raw = daily3_raw,
  anchors = anchors_3_official,
  N = N
)

# Enforce logical constraints: cum3 <= cum2
cum3 <- pmin(m3$cum, cum2)
cum3 <- cummax(cum3)
daily3 <- c(cum3[1], diff(cum3))
daily3 <- pmax(0, daily3)

# -------------------- Monthly output --------------------

month <- format(dates, "%Y-%m")

monthly <- aggregate(
  cbind(first_dose = daily1, second_dose = daily2, third_dose = daily3),
  by = list(month = month),
  FUN = sum
)

monthly$first_dose  <- as.integer(round(monthly$first_dose))
monthly$second_dose <- as.integer(round(monthly$second_dose))
monthly$third_dose  <- as.integer(round(monthly$third_dose))

out_csv <- "data/nhs_staff_vax_model_monthly.csv"
write.csv(monthly, out_csv, row.names = FALSE)
cat("Wrote:", out_csv, "\n\n")

# -------------------- Verification vs anchors --------------------
# Dose 1: verify only "given" anchors (exclude the two start assumptions)
official_1_idx <- anchors_1$date >= as.Date("2021-01-31")
v1 <- verify_anchors(
  "dose1",
  anchors_1$date[official_1_idx],
  anchors_1$coverage[official_1_idx],
  cum1, dates, N
)

# Dose 2: verify official anchors from 2021-08-31 onward (exclude the 0 start anchor)
official_2_idx <- anchors_2_official$date >= as.Date("2021-08-31")
v2 <- verify_anchors(
  "dose2",
  anchors_2_official$date[official_2_idx],
  anchors_2_official$coverage[official_2_idx],
  cum2, dates, N
)

# Dose 3: verify official anchors from 2021-11-30 onward (exclude booster_start 0 anchor)
official_3_idx <- anchors_3_official$date >= as.Date("2021-11-30")
v3 <- verify_anchors(
  "dose3",
  anchors_3_official$date[official_3_idx],
  anchors_3_official$coverage[official_3_idx],
  cum3, dates, N
)

verif <- rbind(v1, v2, v3)

cat("Verification vs anchors (percentage-point error and people):\n")
print(verif, row.names = FALSE)

cat("\nMax absolute error (pct-pt) by series:\n")
print(aggregate(abs(verif$err_pct_pt), by = list(series = verif$series), FUN = max))

cat("\nNotes:\n",
    "- Dose 2 is generated from Dose 1 with a hard minimum gap of 21 days.\n",
    "- Dose 3 is generated from Dose 2 with a hard minimum gap of 90 days and forced to 0 before ",
    format(booster_start), ".\n",
    "- If an anchor decreases month-to-month, the model cannot follow the decrease (administered doses can't be negative),\n",
    "  so it will match a monotone version of those anchors and you'll see a small residual error in verification.\n",
    sep = "")

cat("\nFirst 10 months of modeled administered doses:\n")
print(head(monthly, 10), row.names = FALSE)
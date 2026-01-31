# ============================================================================
# NHS Sickness Rate Analysis — FIXED & FUNCTIONAL (v2)
# Breakpoints vs “observable norm” (2010–2018), seasonality-adjusted
# ============================================================================
# The script implements a retrospective time-series analysis of monthly NHS sickness absence rates in England,
# designed to identify departures from an “observable norm” defined by the pre-event period 2010–2018.
# Monthly percentages are converted to a rate per 100,000 staff and a 6-month rolling mean is computed using a symmetric window
# (3 months prior through 2 months after the index month) to provide a smoothed descriptive series.

# To operationalize the “norm,” the script fits an ordinary least squares baseline model on 2010–2018 raw monthly rates
# with a deterministic linear time component and explicit seasonality via month-of-year fixed effects.
# This model is then projected over the full series, producing fitted baseline expectations and both confidence and prediction intervals;
# the baseline expectation is additionally “rolled” to align directly with the 6-month rolling observed series for visual and interpretive comparability.
# Deviations from the norm are computed as both absolute differences (per-100k) and percentage differences (rolling series relative to the rolled baseline expectation),
# which are later used for annotation of salient local highs/lows.

# Structural change is assessed using the Bai–Perron multiple breakpoint framework applied not to the raw sickness rates,
# but to the baseline-adjusted deviation series (raw deviation-from-norm).
# This explicitly targets breakpoints in abnormality relative to the 2010–2018 regime rather than breaks driven by long-run
# secular trend or seasonal structure.

# The script supports two breakpoint specifications: an intercept-only model to detect step changes in deviation (“level”)
# and a time-augmented model to detect changes in deviation slope (“trend”), with the number of breaks selected by BIC under minimum segment-size constraints.

# For each detected breakpoint, approximate 95% confidence intervals for breakpoint timing are computed, and break direction is summarized via local mean deviations
# in windows around the break. Finally, the script estimates simple within-segment linear trends on the smoothed (rolling) series between
# break-defined boundaries for descriptive comparison and produces a publication-style plot combining raw data,
# rolling series, baseline expectation, prediction band, break markers, and deviation labels at local extrema.
# ============================================================================

# ---- Packages ----
# install.packages(c("tidyverse", "lubridate", "slider", "scales", "ggrepel", "strucchange"))
library(tidyverse)
library(lubridate)
library(slider)
library(scales)
library(ggrepel)
library(strucchange)

# ---- User-tunable parameters ----
DATA_PATH <- "data/overall_sickness_rates.csv"

BASELINE_START <- 2010
BASELINE_END   <- 2018
PROJ_START     <- 2019
PROJ_END       <- 2025

# Rolling window: (t-3, t-2, t-1, t, t+1, t+2)  => 6 points
ROLL_BEFORE <- 3
ROLL_AFTER  <- 2
ROLL_COMPLETE <- TRUE

# Breakpoints configuration
MAX_BREAKS <- 10
H_FRACTION <- 0.10    # minimum 10% of data per segment

# Breakpoint model choice:
# - "level": step changes in deviation (dev_100k_raw ~ 1)
# - "trend": changes in deviation trend (dev_100k_raw ~ month_num)
BP_MODEL_TYPE <- "level"

# ---- Helpers ----

safe_month <- function(x) {
  # Tries multiple month formats robustly; returns Date (first of month) or NA
  x <- as.character(x)
  
  d <- suppressWarnings(my(x))  # good for "Jan 2010", "Jan-2010", etc
  
  if (anyNA(d)) {
    d2 <- suppressWarnings(parse_date_time(
      x,
      orders = c(
        "my", "Ym", "Y-m", "Y/m",
        "bY", "Yb", "mY", "m/Y",
        "b-Y", "Y-b", "Ymd", "Y-m-d", "Y/m/d"
      ),
      tz = "UTC"
    ))
    d[is.na(d)] <- d2[is.na(d)]
  }
  
  # normalize to first day of month
  d <- floor_date(as.Date(d), unit = "month")
  as.Date(d)
}

get_optimal_breaks_bic <- function(bp_full_obj) {
  s <- summary(bp_full_obj)
  
  # BIC is typically provided as s$BIC
  if (!is.null(s$BIC)) {
    bic <- s$BIC
  } else if (!is.null(s$RSS) && is.matrix(s$RSS) && "BIC" %in% rownames(s$RSS)) {
    bic <- s$RSS["BIC", ]
  } else {
    stop("Could not extract BIC from breakpoint summary. Check strucchange version/output.")
  }
  
  idx <- which.min(bic)
  
  # Prefer named break counts if present
  if (!is.null(names(bic)) && nzchar(names(bic)[idx])) {
    return(as.integer(names(bic)[idx]))
  }
  
  # Fallback: assume entries correspond to 0,1,2,... breaks
  as.integer(idx - 1)
}

cap_max_breaks <- function(n, h, requested) {
  # maximum feasible breaks given minimum segment fraction h:
  # each segment must have at least floor(h*n) points; so max segments ≈ floor(1/h)
  # => max breaks ≈ floor(1/h) - 1
  max_by_h <- max(0, floor(1 / h) - 1)
  max_by_n <- max(0, n - 2)  # need at least 2 points to define any segmentation sensibly
  min(requested, max_by_h, max_by_n)
}

# ============================================================================
# Read data
# ============================================================================
raw <- read.csv(DATA_PATH, stringsAsFactors = FALSE) |>
  as_tibble()

df <- raw |>
  rename(rate_pct = England) |>
  mutate(
    month    = safe_month(Month),
    rate_pct = as.numeric(rate_pct)
  ) |>
  arrange(month) |>
  mutate(
    rate_100k = rate_pct * 1000,                 # 1% = 1000 per 100k
    month_num = as.numeric(month),
    month_of_year = factor(month(month), levels = 1:12, labels = month.abb)
  )

# Report failed month parses (your warning said 2)
bad_month_rows <- df |> filter(is.na(month))
if (nrow(bad_month_rows) > 0) {
  warning(sprintf("Month parsing failed for %d rows. Showing the Month values:", nrow(bad_month_rows)))
  print(bad_month_rows |> select(Month) |> distinct())
}

# ---- 6-month rolling average ----
df <- df |>
  mutate(
    roll6_100k = slide_dbl(
      rate_100k,
      ~ mean(.x, na.rm = TRUE),
      .before = ROLL_BEFORE,
      .after  = ROLL_AFTER,
      .complete = ROLL_COMPLETE
    )
  )

# ============================================================================
# Baseline model ("observable norm") 2010–2018: seasonality + time
# ============================================================================
df_fit <- df |>
  filter(year(month) >= BASELINE_START, year(month) <= BASELINE_END) |>
  filter(!is.na(rate_100k), !is.na(month_num), !is.na(month_of_year))

fit_2010_2018 <- lm(rate_100k ~ month_num + month_of_year, data = df_fit)

slope_2010_2018 <- unname(coef(fit_2010_2018)["month_num"]) * 365.25

# Predictions: mean + CI + prediction interval
pred_ci <- as_tibble(predict(fit_2010_2018, newdata = df, interval = "confidence", level = 0.95)) |>
  rename(norm_fit_100k = fit, norm_ci_low_100k = lwr, norm_ci_high_100k = upr)

pred_pi <- as_tibble(predict(fit_2010_2018, newdata = df, interval = "prediction", level = 0.95)) |>
  rename(norm_pi_low_100k = lwr, norm_pi_high_100k = upr) |>
  select(norm_pi_low_100k, norm_pi_high_100k)

df <- bind_cols(df, pred_ci, pred_pi) |>
  mutate(
    period = case_when(
      year(month) >= BASELINE_START & year(month) <= BASELINE_END ~ "fit_2010_2018",
      year(month) >= PROJ_START     & year(month) <= PROJ_END     ~ "proj_2019_2025",
      TRUE ~ "other"
    )
  )

# Roll the baseline expectation and bands to match rolling plotted series
df <- df |>
  mutate(
    norm_fit_roll6     = slide_dbl(norm_fit_100k, mean, .before = ROLL_BEFORE, .after = ROLL_AFTER, .complete = ROLL_COMPLETE),
    norm_pi_low_roll6  = slide_dbl(norm_pi_low_100k, mean, .before = ROLL_BEFORE, .after = ROLL_AFTER, .complete = ROLL_COMPLETE),
    norm_pi_high_roll6 = slide_dbl(norm_pi_high_100k, mean, .before = ROLL_BEFORE, .after = ROLL_AFTER, .complete = ROLL_COMPLETE)
  )

# Deviations vs baseline
df <- df |>
  mutate(
    dev_100k_raw  = if_else(!is.na(rate_100k)  & !is.na(norm_fit_100k),  rate_100k  - norm_fit_100k,  NA_real_),
    dev_100k_roll = if_else(!is.na(roll6_100k) & !is.na(norm_fit_roll6), roll6_100k - norm_fit_roll6, NA_real_),
    dev_pct_roll  = if_else(!is.na(roll6_100k) & !is.na(norm_fit_roll6) & norm_fit_roll6 != 0,
                            100 * (roll6_100k - norm_fit_roll6) / norm_fit_roll6,
                            NA_real_),
    dev_label_roll = if_else(!is.na(dev_pct_roll), sprintf("%+.1f%%", dev_pct_roll), NA_character_)
  )

# ============================================================================
# EXPORT: Monthly sickness rolling + prediction band (2020–2025)
#   Output aligns with monthly positivity format: "YYYY-MM,value"
# ============================================================================

OUT_SICKNESS_BANDS <- "data/sickness_roll6_prediction_band_2020_2026.csv"

df_sickness_monthly <- df |>
  mutate(month_key = format(month, "%Y-%m")) |>
  filter(month >= as.Date("2019-01-01"), month <= as.Date("2025-11-01")) |>
  select(
    month_key,
    roll6_100k,
    norm_pi_low_roll6,
    norm_pi_high_roll6
  ) |>
  rename(
    month = month_key,
    sickness_roll6_100k = roll6_100k,
    sickness_pred_low_100k = norm_pi_low_roll6,
    sickness_pred_high_100k = norm_pi_high_roll6
  ) |>
  arrange(month)

write.csv(df_sickness_monthly, OUT_SICKNESS_BANDS, row.names = FALSE)
cat(sprintf("\nExported monthly sickness rolling + prediction band: %s\n", OUT_SICKNESS_BANDS))

print(head(df_sickness_monthly, 24))

# ============================================================================
# STRUCTURAL BREAKS (Bai–Perron) ON DEVIATION-FROM-NORM (RAW)
#   IMPORTANT FIX: refit a FULL breakpoints model with breaks=optimal_n
#   so confint() works (no vcov error).
# ============================================================================
cat("\n========================================\n")
cat("BAI-PERRON STRUCTURAL BREAK ANALYSIS\n")
cat("========================================\n")
cat(sprintf("Baseline norm: %d-%d (seasonality + time)\n", BASELINE_START, BASELINE_END))
cat("Breakpoints run on deviation-from-norm (RAW): dev_100k_raw\n")
cat(sprintf("Model type: %s\n\n", BP_MODEL_TYPE))

df_for_bp <- df |>
  filter(!is.na(dev_100k_raw), !is.na(month), !is.na(month_num)) |>
  arrange(month) |>
  mutate(row_idx = row_number())

bp_formula <- if (BP_MODEL_TYPE == "level") {
  dev_100k_raw ~ 1
} else {
  dev_100k_raw ~ month_num
}

# Cap max breaks to what is feasible (prevents "changed to 9" surprises)
max_breaks_eff <- cap_max_breaks(n = nrow(df_for_bp), h = H_FRACTION, requested = MAX_BREAKS)

bp_full <- breakpoints(
  bp_formula,
  data   = df_for_bp,
  h      = H_FRACTION,
  breaks = max_breaks_eff
)

optimal_n <- get_optimal_breaks_bic(bp_full)
cat(sprintf("BIC-optimal number of breaks: %d\n", optimal_n))

if (optimal_n > 0) {
  
  # REFIT as FULL model with breaks=optimal_n (this is the key fix)
  bp_opt_full <- breakpoints(
    bp_formula,
    data   = df_for_bp,
    h      = H_FRACTION,
    breaks = optimal_n
  )
  
  break_indices <- bp_opt_full$breakpoints
  
  all_breaks <- tibble(
    break_num = seq_along(break_indices),
    row_idx   = break_indices,
    month     = df_for_bp$month[break_indices]
  ) |>
    arrange(month) |>
    mutate(
      mean_before = map_dbl(row_idx, ~ mean(df_for_bp$dev_100k_raw[max(1, .x - 6):max(1, .x - 1)], na.rm = TRUE)),
      mean_after  = map_dbl(row_idx, ~ mean(df_for_bp$dev_100k_raw[min(nrow(df_for_bp), .x + 1):min(nrow(df_for_bp), .x + 6)], na.rm = TRUE)),
      direction = if_else(mean_after > mean_before, "up", "down")
    )
  
  # CI now works (because bp_opt_full is "breakpointsfull")
  bp_conf <- confint(bp_opt_full, level = 0.95)
  
  all_breaks <- all_breaks |>
    mutate(
      ci_lower_month = df_for_bp$month[pmax(1, bp_conf$confint[, 1])],
      ci_upper_month = df_for_bp$month[pmin(nrow(df_for_bp), bp_conf$confint[, 3])]
    )
  
  cat("\n----------------------------------------\n")
  cat("DETECTED STRUCTURAL BREAKS (on dev_100k_raw):\n")
  cat("----------------------------------------\n")
  for (i in 1:nrow(all_breaks)) {
    direction_symbol <- if_else(all_breaks$direction[i] == "up", "↑ UP", "↓ DOWN")
    cat(sprintf("  Break %d: %s  %s\n",
                i,
                format(all_breaks$month[i], "%B %Y"),
                direction_symbol))
    cat(sprintf("           95%% CI: [%s - %s]\n",
                format(all_breaks$ci_lower_month[i], "%b %Y"),
                format(all_breaks$ci_upper_month[i], "%b %Y")))
  }
  
  cat("\nSegment coefficients (optimal model):\n")
  print(coef(bp_opt_full))
  
} else {
  cat("\nNo structural breaks detected by BIC criterion.\n")
  all_breaks <- tibble(
    break_num = integer(),
    row_idx   = integer(),
    month     = as.Date(character()),
    direction = character(),
    ci_lower_month = as.Date(character()),
    ci_upper_month = as.Date(character())
  )
}

cat("========================================\n\n")

# ============================================================================
# LOCAL EXTREMA (for deviation labels)
# ============================================================================
cat("Finding all local highs and lows...\n")

df <- df |>
  mutate(
    is_local_max = roll6_100k > lag(roll6_100k) & roll6_100k > lead(roll6_100k),
    is_local_min = roll6_100k < lag(roll6_100k) & roll6_100k < lead(roll6_100k)
  )

local_extrema <- df |>
  filter(year(month) >= BASELINE_START, !is.na(roll6_100k), !is.na(norm_fit_roll6)) |>
  filter(is_local_max | is_local_min) |>
  mutate(
    extrema_type = if_else(is_local_max, "high", "low"),
    dev_label_extrema = sprintf("%+.1f%%", dev_pct_roll)
  ) |>
  select(month, roll6_100k, norm_fit_roll6, extrema_type, dev_pct_roll, dev_label_extrema)

cat(sprintf("Found %d local extrema from %d onwards\n", nrow(local_extrema), BASELINE_START))

# ============================================================================
# SEGMENT TRENDS (on roll6_100k) with correct boundaries
# ============================================================================
cat("\n========================================\n")
cat("SEGMENT TREND ANALYSIS (on roll6_100k)\n")
cat("========================================\n")

segment_trends <- list()

if (nrow(all_breaks) > 0) {
  break_dates <- sort(all_breaks$month)
  boundaries <- c(min(df$month, na.rm = TRUE), break_dates, max(df$month, na.rm = TRUE))
  
  for (i in seq_len(length(boundaries) - 1)) {
    start_date <- boundaries[i]
    end_date   <- boundaries[i + 1]
    
    seg_name <- if (i == 1) "pre_B1" else sprintf("B%d_to_B%d", i - 1, i)
    
    df_segment <- df |>
      filter(month >= start_date, month < end_date, !is.na(roll6_100k), !is.na(month_num))
    
    if (nrow(df_segment) >= 3) {
      fit_segment <- lm(roll6_100k ~ month_num, data = df_segment)
      slope <- unname(coef(fit_segment)[2]) * 365.25
      segment_trends[[seg_name]] <- list(start = start_date, end = end_date, slope = slope, fit = fit_segment)
      
      cat(sprintf("  %-10s (%s to %s): %+.1f per 100k/year\n",
                  seg_name, format(start_date, "%b %Y"), format(end_date, "%b %Y"), slope))
    }
  }
  
  # Add "Bk_to_end"
  k <- length(break_dates)
  df_k_end <- df |> filter(month >= break_dates[k], !is.na(roll6_100k), !is.na(month_num))
  if (nrow(df_k_end) >= 3) {
    fit_k_end <- lm(roll6_100k ~ month_num, data = df_k_end)
    slope_k_end <- unname(coef(fit_k_end)[2]) * 365.25
    segment_trends[[sprintf("B%d_to_end", k)]] <- list(
      start = break_dates[k],
      end   = max(df$month, na.rm = TRUE),
      slope = slope_k_end,
      fit   = fit_k_end
    )
    cat(sprintf("  %-10s (%s to %s): %+.1f per 100k/year\n",
                sprintf("B%d_to_end", k), format(break_dates[k], "%b %Y"),
                format(max(df$month, na.rm = TRUE), "%b %Y"), slope_k_end))
  }
}

cat("\nBaseline slope (2010-2018):\n")
cat(sprintf("  %+0.1f per 100k/year\n", slope_2010_2018))
cat("========================================\n\n")

# ============================================================================
# MATERIALIZE SEGMENT TRENDS INTO df
#   (creates trend_B1_B2, trend_B2_B3, trend_B3_end, trend_B1_end, trend_B2_end)
# NOTE: Uses base R pipe (|>) compatible code (no { } blocks on RHS).
# ============================================================================

# Always create the columns so ggplot filters won't error
df <- df |>
  mutate(
    trend_B1_B2 = NA_real_,
    trend_B2_B3 = NA_real_,
    trend_B3_end = NA_real_,
    trend_B1_end = NA_real_,
    trend_B2_end = NA_real_
  )

add_segment_trend <- function(df, fit, colname, start_inclusive, end_exclusive = NULL) {
  pred <- predict(fit, newdata = df)
  in_range <- if (is.null(end_exclusive)) {
    df$month >= start_inclusive
  } else {
    df$month >= start_inclusive & df$month < end_exclusive
  }
  df[[colname]] <- ifelse(in_range, pred, NA_real_)
  df
}

# Helper: fit roll6 trend for a date range (base-pipe compatible)
fit_roll_trend <- function(df, start_date, end_date = NULL) {
  d <- df |>
    filter(month >= start_date)

  if (!is.null(end_date)) {
    d <- d |>
      filter(month < end_date)
  }

  d <- d |>
    filter(!is.na(roll6_100k), !is.na(month_num))

  if (nrow(d) < 3) return(NULL)
  lm(roll6_100k ~ month_num, data = d)
}

# Helper: safely fetch a segment fit from segment_trends
get_segment_fit <- function(segment_trends, key) {
  if (is.null(segment_trends[[key]])) return(NULL)
  if (is.null(segment_trends[[key]]$fit)) return(NULL)
  segment_trends[[key]]$fit
}

# Build the trend columns if breaks exist
if (exists("all_breaks") && nrow(all_breaks) > 0) {

  # ---- B1->B2 ----
  if (nrow(all_breaks) >= 2) {
    fit_b1_b2 <- get_segment_fit(segment_trends, "B1_to_B2")
    if (is.null(fit_b1_b2)) {
      fit_b1_b2 <- fit_roll_trend(df, all_breaks$month[1], all_breaks$month[2])
    }
    if (!is.null(fit_b1_b2)) {
      df <- add_segment_trend(df, fit_b1_b2, "trend_B1_B2", all_breaks$month[1], all_breaks$month[2])
    }
  }

  # ---- B2->B3 and B3->end ----
  if (nrow(all_breaks) >= 3) {
    fit_b2_b3 <- get_segment_fit(segment_trends, "B2_to_B3")
    if (is.null(fit_b2_b3)) {
      fit_b2_b3 <- fit_roll_trend(df, all_breaks$month[2], all_breaks$month[3])
    }
    if (!is.null(fit_b2_b3)) {
      df <- add_segment_trend(df, fit_b2_b3, "trend_B2_B3", all_breaks$month[2], all_breaks$month[3])
    }

    # B3->end
    fit_b3_end <- get_segment_fit(segment_trends, "B3_to_end")
    if (is.null(fit_b3_end)) {
      fit_b3_end <- get_segment_fit(segment_trends, sprintf("B%d_to_end", 3))
    }
    if (is.null(fit_b3_end)) {
      fit_b3_end <- fit_roll_trend(df, all_breaks$month[3], NULL)
    }
    if (!is.null(fit_b3_end)) {
      df <- add_segment_trend(df, fit_b3_end, "trend_B3_end", all_breaks$month[3], NULL)
    }
  }

  # ---- B1->end (always) ----
  fit_b1_end <- fit_roll_trend(df, all_breaks$month[1], NULL)
  if (!is.null(fit_b1_end)) {
    df <- add_segment_trend(df, fit_b1_end, "trend_B1_end", all_breaks$month[1], NULL)
    segment_trends[["B1_to_end"]] <- list(
      start = all_breaks$month[1],
      end   = max(df$month, na.rm = TRUE),
      slope = unname(coef(fit_b1_end)[2]) * 365.25,
      fit   = fit_b1_end
    )
  }

  # ---- B2->end (if exists) ----
  if (nrow(all_breaks) >= 2) {
    fit_b2_end <- fit_roll_trend(df, all_breaks$month[2], NULL)
    if (!is.null(fit_b2_end)) {
      df <- add_segment_trend(df, fit_b2_end, "trend_B2_end", all_breaks$month[2], NULL)
      segment_trends[["B2_to_end"]] <- list(
        start = all_breaks$month[2],
        end   = max(df$month, na.rm = TRUE),
        slope = unname(coef(fit_b2_end)[2]) * 365.25,
        fit   = fit_b2_end
      )
    }
  }
}

# ============================================================================
# VISUALIZATION
# ============================================================================

markers <- tibble(
  date = as.Date(c("2020-01-31", "2020-12-08", "2021-11-01", "2022-04-15")),
  short_label = c("COVID", "Vax start", "90% vax", "Full vax")
)

y_range <- range(c(df$rate_100k, df$roll6_100k, df$norm_fit_roll6), na.rm = TRUE)
y_top <- y_range[2]
y_bottom <- y_range[1]

markers <- markers |>
  mutate(y_label = if_else(row_number() %% 2 == 1, y_top * 0.98, y_top * 0.92))

# Colors (kept close to your original)
color_raw <- "grey50"
color_rolling <- "#E69F00"
color_baseline <- "#0072B2"
color_markers <- "#D55E00"
color_breaks <- "#CC79A7"
color_segment_1 <- "#E31A1C"
color_segment_2 <- "#009E73"
color_segment_3 <- "#F0E442"
color_B1_end <- "#7B3294"
color_B2_end <- "#1B9E77"
color_deviation <- "#666666"

# ---- X-axis guides: every 3 months + thicker Jan separators ----
x_min <- floor_date(min(df$month, na.rm = TRUE), unit = "month")
x_max <- ceiling_date(max(df$month, na.rm = TRUE), unit = "month")

x_breaks_3m <- seq(x_min, x_max, by = "3 months")
x_breaks_jan <- seq(floor_date(x_min, "year"), floor_date(x_max, "year"), by = "1 year")

grid_3m <- tibble(x = x_breaks_3m)
grid_jan <- tibble(x = x_breaks_jan)

label_3m <- function(d) {
  d <- as.Date(d)
  ifelse(format(d, "%m") == "01", format(d, "%b\n%Y"), format(d, "%b"))
}

p <- ggplot(df, aes(x = month)) +

  # --- Month/year separators FIRST (so they stay behind everything) ---
  geom_vline(
    data = grid_3m,
    aes(xintercept = x),
    color = "grey92",
    linetype = "dashed",
    linewidth = 0.3
  ) +
  geom_vline(
    data = grid_jan,
    aes(xintercept = x),
    color = "grey80",
    linetype = "solid",
    linewidth = 0.8
  ) +

  # Prediction band for projection period (rolling, approximate)
  geom_ribbon(
    data = df |> filter(period == "proj_2019_2025"),
    aes(ymin = norm_pi_low_roll6, ymax = norm_pi_high_roll6),
    fill = color_baseline, alpha = 0.10
  ) +

  # Raw data (subdued)
  geom_line(aes(y = rate_100k), color = color_raw, linewidth = 0.8, alpha = 0.5) +

  # Rolling average (prominent)
  geom_line(aes(y = roll6_100k), color = color_rolling, linewidth = 1.5, na.rm = TRUE) +

  # Baseline expected (rolling): solid in fit period
  geom_line(
    data = df |> filter(period == "fit_2010_2018"),
    aes(y = norm_fit_roll6),
    color = color_baseline, linewidth = 1.4
  ) +

  # Baseline expected (rolling): dashed in projection
  geom_line(
    data = df |> filter(period == "proj_2019_2025"),
    aes(y = norm_fit_roll6),
    color = color_baseline, linewidth = 1.4, linetype = "dashed"
  ) +

  # Segment trend lines
  geom_line(
    data = df |> filter(!is.na(trend_B1_B2)),
    aes(y = trend_B1_B2),
    color = color_segment_1, linewidth = 1.4
  ) +
  geom_line(
    data = df |> filter(!is.na(trend_B2_B3)),
    aes(y = trend_B2_B3),
    color = color_segment_2, linewidth = 1.4
  ) +
  geom_line(
    data = df |> filter(!is.na(trend_B3_end)),
    aes(y = trend_B3_end),
    color = color_segment_3, linewidth = 1.4
  ) +
  geom_line(
    data = df |> filter(!is.na(trend_B1_end)),
    aes(y = trend_B1_end),
    color = color_B1_end, linewidth = 1.2, linetype = "dashed"
  ) +
  geom_line(
    data = df |> filter(!is.na(trend_B2_end)),
    aes(y = trend_B2_end),
    color = color_B2_end, linewidth = 1.2, linetype = "dotdash"
  ) +

  # Deviation labels at local extrema (NOW above month/year separators)
  geom_label_repel(
    data = local_extrema,
    aes(x = month, y = roll6_100k, label = dev_label_extrema),
    color = color_deviation,
    fill = alpha("white", 0.85),
    size = 2.8,
    fontface = "bold",
    label.padding = unit(0.15, "lines"),
    box.padding = unit(0.3, "lines"),
    point.padding = unit(0.2, "lines"),
    segment.color = color_deviation,
    segment.size = 0.3,
    direction = "y",
    nudge_y = if_else(local_extrema$extrema_type == "high", 100, -100),
    min.segment.length = 0,
    max.overlaps = 50,
    seed = 42
  ) +

  # Event markers
  geom_vline(
    data = markers,
    aes(xintercept = date),
    linetype = "dotted", linewidth = 0.8,
    color = color_markers, alpha = 0.7
  ) +
  geom_label(
    data = markers,
    aes(x = date, y = y_label, label = short_label),
    color = color_markers,
    fill = alpha("white", 0.9),
    size = 2.8,
    fontface = "bold",
    label.padding = unit(0.12, "lines"),
    label.size = 0.2
  ) +

  # Structural breaks (kept above separators; and above labels if overlapping)
  {
    if (nrow(all_breaks) > 0) {
      list(
        geom_vline(
          data = all_breaks,
          aes(xintercept = month),
          linetype = "longdash", linewidth = 1.0,
          color = color_breaks, alpha = 0.8
        ),
        geom_label(
          data = all_breaks |> mutate(
            y_pos = if_else(row_number() %% 2 == 1,
                            y_bottom + (y_top - y_bottom) * 0.15,
                            y_bottom + (y_top - y_bottom) * 0.06),
            label_text = sprintf("B%d %s\n%s",
                                 break_num,
                                 if_else(direction == "up", "↑", "↓"),
                                 format(month, "%b %Y"))
          ),
          aes(x = month, y = y_pos, label = label_text),
          color = color_breaks,
          fill = alpha("white", 0.92),
          size = 2.3,
          fontface = "bold",
          label.padding = unit(0.12, "lines"),
          lineheight = 0.9
        )
      )
    }
  } +

  scale_x_date(
    breaks = x_breaks_3m,
    labels = label_3m,
    expand = expansion(mult = c(0.02, 0.05))
  ) +
  scale_y_continuous(
    labels = label_comma(),
    expand = expansion(mult = c(0.02, 0.12))
  ) +

  labs(
    title = "NHS Sickness Absence Rate (England) - Breaks vs 2010-2018 Norm",
    subtitle = paste(
      "Orange: 6-month rolling average",
      "| Blue: baseline norm (2010–2018, seasonality-adjusted)",
      "| Ribbon: ~95% prediction band (approx., rolled)",
      "| Pink: structural breaks on deviation-from-norm",
      sep = " "
    ),
    x = NULL,
    y = "Sickness rate (per 100,000 staff)",
    caption = {
      cap <- sprintf("%d structural breaks (BIC-optimal) | Baseline slope: %+.0f/yr",
                     nrow(all_breaks), slope_2010_2018)
      if (!is.null(segment_trends[["B1_to_B2"]])) cap <- paste0(cap, sprintf(" | B1→B2: %+.0f/yr", segment_trends[["B1_to_B2"]]$slope))
      if (!is.null(segment_trends[["B2_to_B3"]])) cap <- paste0(cap, sprintf(" | B2→B3: %+.0f/yr", segment_trends[["B2_to_B3"]]$slope))
      if (!is.null(segment_trends[["B1_to_end"]])) cap <- paste0(cap, sprintf(" | B1→end: %+.0f/yr", segment_trends[["B1_to_end"]]$slope))
      cap
    }
  ) +

  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(size = 16, face = "bold", margin = margin(b = 5)),
    plot.subtitle = element_text(size = 9, color = "grey40", margin = margin(b = 10)),
    plot.caption = element_text(size = 8, color = "grey50", hjust = 0, margin = margin(t = 10)),
    axis.text.x = element_text(size = 9, angle = 90, vjust = 0.5, hjust = 1),
    axis.text.y = element_text(size = 10),
    axis.title.y = element_text(size = 11, margin = margin(r = 10)),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "grey85", linewidth = 0.4),
    plot.margin = margin(t = 10, r = 15, b = 25, l = 10),
    legend.position = "none"
  )

print(p)

ggsave(
  "nhs_sickness_analysis_fixed.png",
  plot = p,
  width = 16,
  height = 9,
  dpi = 300,
  bg = "white"
)

cat("\nPlot saved as 'nhs_sickness_analysis_fixed.png'\n")

# ============================================================================
# SUMMARY OUTPUT
# ============================================================================
cat("\n")
cat("========================================================================\n")
cat("                          ANALYSIS SUMMARY                               \n")
cat("========================================================================\n")
cat(sprintf("Structural breaks detected (BIC-optimal): %d\n", nrow(all_breaks)))
if (nrow(all_breaks) > 0) {
  for (i in 1:nrow(all_breaks)) {
    cat(sprintf("  Break %d: %s (%s)\n",
                i,
                format(all_breaks$month[i], "%B %Y"),
                if_else(all_breaks$direction[i] == "up", "upward", "downward")))
  }
}
cat("------------------------------------------------------------------------\n")
cat("TREND SLOPES (per 100k per year):\n")
cat(sprintf("  Baseline norm (2010-2018): %+7.1f\n", slope_2010_2018))

if (length(segment_trends) > 0) {
  for (nm in names(segment_trends)) {
    tr <- segment_trends[[nm]]
    # guard against divide-by-zero if slope baseline is ~0
    multiplier <- if (isTRUE(all.equal(slope_2010_2018, 0))) NA_real_ else tr$slope / slope_2010_2018
    cat(sprintf("  %-24s %+7.1f  (%s)\n",
                paste0(nm, ":"),
                tr$slope,
                ifelse(is.na(multiplier), "n/a vs baseline", sprintf("%.1fx baseline", multiplier))))
  }
}
cat("========================================================================\n")
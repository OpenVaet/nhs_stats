# =============================================================================
# sickness_orgtype_vs_vax_doses.R
# =============================================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(lubridate)
  library(scales)
})

# ---- Force English month names for axis labels (Windows/French OS safe) ----
old_time_locale <- Sys.getlocale("LC_TIME")
on.exit(try(Sys.setlocale("LC_TIME", old_time_locale), silent = TRUE), add = TRUE)
ok <- suppressWarnings(Sys.setlocale("LC_TIME", "C"))
if (is.na(ok) || ok == "") {
  for (loc in c("English_United Kingdom.1252", "English_United States.1252",
                "en_GB.UTF-8", "en_US.UTF-8")) {
    ok <- suppressWarnings(Sys.setlocale("LC_TIME", loc))
    if (!is.na(ok) && ok != "") break
  }
}

SICK_UNIFIED_PATH <- "data/nhs_sickness_rates_unified_2014_2025.csv"
VAX_PATH          <- "data/nhs_staff_vax_model_monthly.csv"
POS_PATH          <- "data/tests_positivity_and_volume.csv"

PLOT_START <- as.Date("2014-01-01")
PLOT_END   <- as.Date("2025-12-01")

# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------
to_num <- function(x) {
  x <- trimws(as.character(x))
  x[x == ""] <- NA_character_
  as.numeric(gsub(",", "", x, fixed = TRUE))
}

parse_month <- function(x) {
  x <- as.character(x)
  as.Date(if_else(
    str_detect(x, "^\\d{4}-\\d{2}$"),
    paste0(x, "-01"),
    x
  ))
}

# -----------------------------------------------------------------------------
# Load unified sickness + aggregate monthly by Org Type
# -----------------------------------------------------------------------------
sick_raw <- readr::read_csv(SICK_UNIFIED_PATH, show_col_types = FALSE)
if (!("Month" %in% names(sick_raw))) {
  stop("Expected column 'Month' (YYYY-MM) in: ", SICK_UNIFIED_PATH)
}

sick_by_type <- sick_raw |>
  transmute(
    month    = as.Date(paste0(Month, "-01")),
    org_type = as.character(`Org Type`),
    sick     = to_num(`FTE Days Sick`),
    avail    = to_num(`FTE Days Available`)
  ) |>
  filter(!is.na(month)) |>
  filter(month >= PLOT_START, month <= PLOT_END) |>
  filter(!is.na(sick), sick != 0) |>
  group_by(month, org_type) |>
  summarise(
    sick_sum  = sum(sick,  na.rm = TRUE),
    avail_sum = sum(avail, na.rm = TRUE),
    .groups = "drop"
  )

# Monthly totals (for total rate)
month_totals <- sick_by_type |>
  group_by(month) |>
  summarise(
    sick_total  = sum(sick_sum,  na.rm = TRUE),
    avail_total = sum(avail_sum, na.rm = TRUE),
    .groups = "drop"
  ) |>
  mutate(
    total_rate_100k = if_else(avail_total > 0, sick_total / avail_total * 1e5, NA_real_)
  )

# Stack contributions by org type (sum to total rate each month)
sick_stack <- sick_by_type |>
  left_join(month_totals |> select(month, avail_total), by = "month") |>
  mutate(contrib_rate_100k = if_else(avail_total > 0, sick_sum / avail_total * 1e5, NA_real_)) |>
  filter(!is.na(contrib_rate_100k))

# -----------------------------------------------------------------------------
# LESS BUSY: collapse small org types into one bucket
# Keep top N by mean contribution; everything else -> "Other org types"
# -----------------------------------------------------------------------------
TOP_N <- 6

top_org <- sick_stack |>
  group_by(org_type) |>
  summarise(m = mean(contrib_rate_100k, na.rm = TRUE), .groups = "drop") |>
  arrange(desc(m)) |>
  slice_head(n = TOP_N) |>
  pull(org_type)

sick_stack <- sick_stack |>
  mutate(org_type_plot = if_else(org_type %in% top_org, org_type, "Other org types")) |>
  group_by(month, org_type_plot) |>
  summarise(contrib_rate_100k = sum(contrib_rate_100k, na.rm = TRUE), .groups = "drop")

org_order <- sick_stack |>
  group_by(org_type_plot) |>
  summarise(m = mean(contrib_rate_100k, na.rm = TRUE), .groups = "drop") |>
  arrange(desc(m)) |>
  pull(org_type_plot)

# Force bucket at the end
org_order <- c(setdiff(org_order, "Other org types"), "Other org types")

sick_stack <- sick_stack |>
  mutate(org_type_plot = factor(org_type_plot, levels = org_order))

# -----------------------------------------------------------------------------
# Load vax monthly doses (robust month parse)
# -----------------------------------------------------------------------------
vax <- readr::read_csv(VAX_PATH, show_col_types = FALSE) |>
  transmute(
    month_raw  = as.character(month),
    month      = parse_month(month_raw),
    first_dose  = to_num(first_dose),
    second_dose = to_num(second_dose),
    third_dose  = to_num(third_dose)
  ) |>
  filter(!is.na(month)) |>
  arrange(month) |>
  mutate(
    first_dose  = replace_na(first_dose, 0),
    second_dose = replace_na(second_dose, 0),
    third_dose  = replace_na(third_dose, 0),
    total_doses = first_dose + second_dose + third_dose
  )

vax_start <- min(vax$month, na.rm = TRUE)
vax_end   <- max(vax$month, na.rm = TRUE)

# -----------------------------------------------------------------------------
# Load positivity (robust month parse)
# -----------------------------------------------------------------------------
pos <- readr::read_csv(POS_PATH, show_col_types = FALSE) |>
  transmute(
    month_raw  = as.character(month),
    month      = parse_month(month_raw),
    pos_monthly = to_num(positivity_mean_7day)
  ) |>
  filter(!is.na(month)) |>
  arrange(month)

# -----------------------------------------------------------------------------
# Join line dataframe (totals + vax + positivity)
# -----------------------------------------------------------------------------
df_line <- month_totals |>
  left_join(vax |> select(month, total_doses), by = "month") |>
  left_join(pos |> select(month, pos_monthly), by = "month") |>
  filter(month >= PLOT_START, month <= PLOT_END) |>
  arrange(month) |>
  mutate(
    # explicit 0 pre-rollout
    total_doses = if_else(month < vax_start, 0, total_doses),
    # post-model: NA (prevents fake flat tail)
    total_doses = if_else(month > vax_end, NA_real_, total_doses)
  )

# -----------------------------------------------------------------------------
# Axis scaling
# -----------------------------------------------------------------------------
left_max_raw <- max(df_line$total_rate_100k, na.rm = TRUE)
left_max <- left_max_raw * 1.30
left_min <- -0.04 * left_max   # small space below zero so the 0-dose line is visible

max_doses <- max(df_line$total_doses, na.rm = TRUE)
k_dose <- if (is.finite(max_doses) && max_doses > 0) left_max / max_doses else NA_real_
df_line <- df_line |>
  mutate(doses_scaled = if (is.finite(k_dose)) total_doses * k_dose else NA_real_)

max_pos <- max(df_line$pos_monthly, na.rm = TRUE)
k_pos <- if (is.finite(max_pos) && max_pos > 0) left_max / max_pos else NA_real_
df_line <- df_line |>
  mutate(pos_scaled = if (is.finite(k_pos)) pos_monthly * k_pos else NA_real_)

# -----------------------------------------------------------------------------
# Markers / breakpoints / background shading
# -----------------------------------------------------------------------------
all_breaks <- tibble(
  break_num = 1:3,
  month = as.Date(c("2019-09-01", "2021-05-01", "2022-12-01"))
) |>
  mutate(label_text = sprintf("B%d  %s", break_num, format(month, "%b %Y")))

markers <- tibble(
  date = as.Date(c("2020-01-31", "2020-12-08", "2021-11-01", "2022-04-15")),
  short_label = c("COVID", "Vax start", "90% vax & mandates", "Mandates dropped")
) |>
  mutate(y_label = if_else(row_number() %% 2 == 1, left_max * 0.985, left_max * 0.94))

# lighter background shading than before
periods_bg <- tibble(
  period = c("no_tests", "corr_high", "corr_none"),
  xmin   = ymd(c(format(PLOT_START, "%Y-%m-%d"), "2020-05-01", "2022-05-01")),
  xmax   = ymd(c("2020-05-01", "2022-05-01", format(PLOT_END %m+% months(1), "%Y-%m-%d"))),
  fill   = c("#F3F3F3", "#EAF6EC", "#F9EAEA")
)

# Year separators only + 6-month labels
x_min <- floor_date(min(df_line$month, na.rm = TRUE), unit = "month")
x_max <- ceiling_date(max(df_line$month, na.rm = TRUE), unit = "month")
x_breaks_jan <- seq(floor_date(x_min, "year"), floor_date(x_max, "year"), by = "1 year")
grid_jan <- tibble(x = x_breaks_jan)

# -----------------------------------------------------------------------------
# Plot (pastel, less busy)
# -----------------------------------------------------------------------------
color_pos     <- "#C77C2E"   # softer orange
color_doses   <- "#4E79A7"   # softer blue
color_breaks  <- "#B07AA1"   # softer purple
color_markers <- "#C77C2E"

p <- ggplot() +

  # Background shading (lighter)
  geom_rect(
    data = periods_bg,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
    fill = I(periods_bg$fill),
    inherit.aes = FALSE,
    alpha = 0.45
  ) +

  # Year separators only
  geom_vline(
    data = grid_jan,
    aes(xintercept = x),
    color = "grey88",
    linewidth = 0.6
  ) +

  # Stacked columns: org-type contributions (pastel palette)
  geom_col(
    data = sick_stack,
    aes(x = month, y = contrib_rate_100k, fill = org_type_plot),
    width = 26,
    alpha = 0.88,
    colour = "white",
    linewidth = 0.12
  ) +
  scale_fill_brewer(palette = "Set3") +

  # Positivity (halo + line)
  geom_line(
    data = df_line,
    aes(x = month, y = pos_scaled),
    color = "white", linewidth = 3.2, na.rm = TRUE
  ) +
  geom_line(
    data = df_line,
    aes(x = month, y = pos_scaled),
    color = color_pos, linewidth = 1.7, na.rm = TRUE
  ) +

  # Vaccine doses (halo + line)
  geom_line(
    data = df_line,
    aes(x = month, y = doses_scaled),
    color = "white", linewidth = 3.8, na.rm = TRUE
  ) +
  geom_line(
    data = df_line,
    aes(x = month, y = doses_scaled),
    color = color_doses, linewidth = 2.0, na.rm = TRUE
  ) +

  # Event markers (smaller)
  geom_vline(
    data = markers,
    aes(xintercept = date),
    linetype = "dashed", linewidth = 0.9,
    color = color_markers, alpha = 0.70
  ) +
  geom_label(
    data = markers,
    aes(x = date, y = y_label, label = short_label),
    color = "grey15",
    fill = alpha("white", 0.92),
    size = 3.2,
    fontface = "bold",
    label.padding = unit(0.15, "lines"),
    label.size = 0.22
  ) +

  # Breakpoints (smaller)
  geom_vline(
    data = all_breaks,
    aes(xintercept = month),
    linetype = "longdash", linewidth = 1.1,
    color = color_breaks, alpha = 0.80
  ) +
  geom_label(
    data = all_breaks |>
      mutate(
        y_pos = if_else(row_number() %% 2 == 1, left_max * 0.38, left_max * 0.28)
      ),
    aes(x = month, y = y_pos, label = label_text),
    color = color_breaks,
    fill = alpha("white", 0.93),
    size = 2.8,
    fontface = "bold",
    label.padding = unit(0.16, "lines"),
    lineheight = 0.95
  ) +

  scale_x_date(
    date_breaks = "6 months",
    date_labels = "%b\n%Y",
    expand = expansion(mult = c(0.01, 0.03))
  ) +

  scale_y_continuous(
    limits = c(left_min, left_max),
    labels = label_comma(),
    name = "Sickness rate (FTE-days sick per 100k FTE-days available)\n(stacked by organisation type contribution)",
    sec.axis = sec_axis(
      ~ if (is.finite(k_dose)) pmax(0, . / k_dose) else NA_real_,
      name = "Vaccine doses administered (monthly)",
      breaks = if (is.finite(max_doses) && max_doses > 0) pretty(c(0, max_doses), n = 6) else waiver(),
      labels = label_number(scale_cut = cut_short_scale())
    )
  ) +

  labs(
    title = "NHS sickness (by organisation type) vs vaccine doses and test positivity",
    subtitle = "Stacked bars: org-type contributions | Blue: doses (right axis) | Orange: positivity (scaled)",
    x = NULL,
    fill = "Org type"
  ) +

  theme_minimal(base_size = 18) +
  theme(
    plot.title = element_text(size = 22, face = "bold", margin = margin(b = 6)),
    plot.subtitle = element_text(size = 13, color = "grey35", margin = margin(b = 12)),
    axis.text.x = element_text(size = 12, angle = 0, vjust = 0.5, hjust = 0.5),
    axis.text.y = element_text(size = 13),
    axis.title.y = element_text(size = 14, margin = margin(r = 10)),
    axis.title.y.right = element_text(size = 14, margin = margin(l = 10)),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "grey90", linewidth = 0.5),
    plot.margin = margin(t = 14, r = 18, b = 18, l = 14),
    legend.position = "bottom",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text  = element_text(size = 11),
    legend.key.height = unit(0.45, "cm"),
    legend.key.width  = unit(0.9, "cm")
  )

print(p)

ggsave(
  "nhs_sickness_orgtype_vs_vax_positivity_substack.png",
  plot = p,
  width = 18,
  height = 10,
  dpi = 320,
  bg = "white"
)

cat("\nPlot saved as 'nhs_sickness_orgtype_vs_vax_positivity_substack.png'\n")

# =============================================================================
# Export merged monthly dataset (wide) - restricted to 2019..2024
# =============================================================================
EXPORT_START <- as.Date("2019-01-01")
EXPORT_END   <- as.Date("2024-12-01")

# Wide by-org-type contributions
stack_wide <- sick_stack |>
  filter(month >= EXPORT_START, month <= EXPORT_END) |>
  mutate(
    month_yyyy_mm = format(month, "%Y-%m"),
    org_type_clean = str_replace_all(as.character(org_type), "[^A-Za-z0-9]+", "_"),
    org_type_clean = str_replace_all(org_type_clean, "_+", "_"),
    org_type_clean = str_replace_all(org_type_clean, "^_|_$", ""),
    col_name = paste0("rate_100k_", org_type_clean)
  ) |>
  select(month_yyyy_mm, col_name, contrib_rate_100k) |>
  pivot_wider(names_from = col_name, values_from = contrib_rate_100k, values_fill = 0)

# Lines data (totals + doses + positivity)
line_out <- df_line |>
  filter(month >= EXPORT_START, month <= EXPORT_END) |>
  transmute(
    month_yyyy_mm = format(month, "%Y-%m"),
    total_rate_100k,
    total_doses,
    pos_monthly
  )

merged_out <- line_out |>
  left_join(stack_wide, by = "month_yyyy_mm") |>
  arrange(month_yyyy_mm)

dir.create("data", showWarnings = FALSE, recursive = TRUE)
out_csv <- "data/nhs_sickness_orgtype_vax_pos_merged_monthly_2019_2024.csv"
readr::write_csv(merged_out, out_csv)

cat("\nMerged monthly dataset saved as: ", out_csv, "\n", sep = "")
cat(sprintf("Rows: %d | Columns: %d\n", nrow(merged_out), ncol(merged_out)))

# =============================================================================
# Stats verification for narrative (B1/B2/B3 periods)
# =============================================================================

B1 <- as.Date("2019-09-01")
B2 <- as.Date("2021-05-01")
B3 <- as.Date("2022-12-01")

# Use the monthly totals + joined series we already built
df_stats <- df_line |>
  arrange(month) |>
  mutate(
    # month index for "per-month" slope
    idx = row_number() - 1L
  )

# Helper: safe correlation
safe_cor <- function(x, y) {
  if (sum(is.finite(x) & is.finite(y)) < 3) return(NA_real_)
  cor(x, y, use = "complete.obs")
}

# Helper: linear slope per month (y ~ idx) inside a period
slope_per_month <- function(df, ycol) {
  y <- df[[ycol]]
  if (sum(is.finite(y)) < 3) return(NA_real_)
  coef(lm(y ~ idx, data = df))[2]
}

# Helper: get peak value + month
peak_info <- function(df, ycol) {
  y <- df[[ycol]]
  if (!any(is.finite(y))) return(tibble(peak_month = as.Date(NA), peak_value = NA_real_))
  i <- which.max(y)
  tibble(peak_month = df$month[i], peak_value = y[i])
}

# Period slicer
slice_period <- function(df, start, end, inclusive_end = TRUE) {
  if (is.na(end)) {
    df |> filter(month > start)
  } else if (inclusive_end) {
    df |> filter(month >= start, month <= end)
  } else {
    df |> filter(month >= start, month < end)
  }
}

# Summarise a period (matches the narrative style stats)
summarise_period <- function(df, label) {
  # Sickness stats
  sick_mean <- mean(df$total_rate_100k, na.rm = TRUE)
  sick_min  <- min(df$total_rate_100k, na.rm = TRUE)
  sick_max  <- max(df$total_rate_100k, na.rm = TRUE)
  sick_sd   <- sd(df$total_rate_100k, na.rm = TRUE)
  sick_slope <- slope_per_month(df, "total_rate_100k")
  sick_peak  <- peak_info(df, "total_rate_100k")

  # Positivity stats (coverage + peak + correlation)
  pos_non_na <- sum(is.finite(df$pos_monthly))
  pos_n      <- nrow(df)
  pos_mean   <- mean(df$pos_monthly, na.rm = TRUE)
  pos_peak   <- peak_info(df, "pos_monthly")
  r_pos_sick <- safe_cor(df$pos_monthly, df$total_rate_100k)

  # Doses stats (including zeros, like your merged series)
  dose_non_zero <- sum(df$total_doses > 0, na.rm = TRUE)
  dose_mean     <- mean(df$total_doses, na.rm = TRUE)
  dose_peak     <- peak_info(df, "total_doses")
  dose_slope    <- slope_per_month(df, "total_doses")
  r_dose_sick   <- safe_cor(df$total_doses, df$total_rate_100k)

  tibble(
    period = label,
    months = nrow(df),

    sick_mean_100k = sick_mean,
    sick_min_100k  = sick_min,
    sick_max_100k  = sick_max,
    sick_sd_100k   = sick_sd,
    sick_slope_per_month = sick_slope,
    sick_peak_month = sick_peak$peak_month,
    sick_peak_100k   = sick_peak$peak_value,

    pos_months_present = pos_non_na,
    pos_months_total   = pos_n,
    pos_mean_pct       = pos_mean,
    pos_peak_month     = pos_peak$peak_month,
    pos_peak_pct       = pos_peak$peak_value,
    r_pos_vs_sick      = r_pos_sick,

    dose_months_nonzero = dose_non_zero,
    dose_mean           = dose_mean,
    dose_peak_month     = dose_peak$peak_month,
    dose_peak           = dose_peak$peak_value,
    dose_slope_per_month = dose_slope,
    r_dose_vs_sick      = r_dose_sick
  )
}

# Build 3 periods (inclusive ends for B1->B2 and B2->B3; post-B3 is strictly after)
df_p1 <- slice_period(df_stats, B1, B2, inclusive_end = TRUE)
df_p2 <- slice_period(df_stats, B2, B3, inclusive_end = TRUE)
df_p3 <- slice_period(df_stats, B3, NA, inclusive_end = TRUE)  # month > B3

stats_tbl <- bind_rows(
  summarise_period(df_p1, "B1_to_B2 (2019-09 .. 2021-05)"),
  summarise_period(df_p2, "B2_to_B3 (2021-05 .. 2022-12)"),
  summarise_period(df_p3, "post_B3 (> 2022-12)")
)

cat("\n=== Narrative stats verification (totals / positivity / doses) ===\n")
print(stats_tbl)

# Also: Org Type shares (average share of total sickness rate, by period)
# share = contrib_rate_100k / total_rate_100k, averaged over months
share_by_period <- function(start, end, label, inclusive_end = TRUE) {
  # monthly total for share denominator
  totals <- month_totals |>
    filter(if (is.na(end)) month > start else if (inclusive_end) month >= start & month <= end else month >= start & month < end) |>
    select(month, total_rate_100k)

  sick_stack |>
    inner_join(totals, by = "month") |>
    mutate(share = if_else(total_rate_100k > 0, contrib_rate_100k / total_rate_100k, NA_real_)) |>
    group_by(org_type) |>
    summarise(
      share_mean = mean(share, na.rm = TRUE),
      contrib_mean_100k = mean(contrib_rate_100k, na.rm = TRUE),
      .groups = "drop"
    ) |>
    arrange(desc(share_mean)) |>
    mutate(period = label) |>
    select(period, org_type, share_mean, contrib_mean_100k)
}

org_shares <- bind_rows(
  share_by_period(B1, B2, "B1_to_B2 (2019-09 .. 2021-05)", inclusive_end = TRUE),
  share_by_period(B2, B3, "B2_to_B3 (2021-05 .. 2022-12)", inclusive_end = TRUE),
  share_by_period(B3, NA, "post_B3 (> 2022-12)")
)

cat("\n=== Top Org Types by average share of sickness rate (per period) ===\n")
org_shares |>
  group_by(period) |>
  slice_head(n = 8) |>
  ungroup() |>
  print(n = Inf)

# Optional: export these verification tables to CSV for audit trail
dir.create("data", showWarnings = FALSE, recursive = TRUE)
readr::write_csv(stats_tbl,  "data/narrative_stats_verification_b1_b2_b3.csv")
readr::write_csv(org_shares, "data/narrative_orgtype_shares_b1_b2_b3.csv")

cat("\nWrote:\n- data/narrative_stats_verification_b1_b2_b3.csv\n- data/narrative_orgtype_shares_b1_b2_b3.csv\n")

# =============================================================================
# Sanity print (matches the narrative windows + rounding)
# =============================================================================

library(dplyr)

# Narrative windows (as used in the verification)
P1_START <- as.Date("2019-09-01"); P1_END <- as.Date("2021-04-01")  # B1 .. just before B2
P2_START <- as.Date("2021-05-01"); P2_END <- as.Date("2022-11-01")  # B2 .. just before B3
P3_START <- as.Date("2022-12-01"); P3_END <- as.Date("2024-12-01")  # from B3 onward (in export)

safe_cor <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  if (sum(ok) < 3) return(NA_real_)
  cor(x[ok], y[ok])
}

peak_row <- function(df, col) {
  v <- df[[col]]
  if (!any(is.finite(v))) return(tibble(month = as.Date(NA), value = NA_real_))
  i <- which.max(v)
  tibble(month = df$month[i], value = v[i])
}

fmt_int <- function(x) ifelse(is.finite(x), format(round(x), big.mark = ","), "NA")
fmt_1   <- function(x) ifelse(is.finite(x), sprintf("%.1f", x), "NA")
fmt_2   <- function(x) ifelse(is.finite(x), sprintf("%.2f", x), "NA")

summ_period <- function(df, label) {
  df <- df %>% arrange(month)
  df <- df %>% mutate(idx = row_number() - 1L)

  # Sickness
  sick_mean <- mean(df$total_rate_100k, na.rm = TRUE)
  sick_min  <- min(df$total_rate_100k, na.rm = TRUE)
  sick_max  <- max(df$total_rate_100k, na.rm = TRUE)
  sick_sd   <- sd(df$total_rate_100k, na.rm = TRUE)
  sick_slope <- if (sum(is.finite(df$total_rate_100k)) >= 3)
    coef(lm(total_rate_100k ~ idx, data = df))[2] else NA_real_
  sick_peak <- peak_row(df, "total_rate_100k")

  # Positivity
  pos_present <- sum(is.finite(df$pos_monthly))
  pos_total   <- nrow(df)
  pos_mean <- mean(df$pos_monthly, na.rm = TRUE)
  pos_peak <- peak_row(df, "pos_monthly")
  r_pos_sick <- safe_cor(df$pos_monthly, df$total_rate_100k)

  # Doses
  dose_nonzero <- sum(df$total_doses > 0, na.rm = TRUE)
  dose_mean <- mean(df$total_doses, na.rm = TRUE)
  dose_peak <- peak_row(df, "total_doses")
  dose_slope <- if (sum(is.finite(df$total_doses)) >= 3)
    coef(lm(total_doses ~ idx, data = df))[2] else NA_real_
  r_dose_sick <- safe_cor(df$total_doses, df$total_rate_100k)

  cat("\n============================================================\n")
  cat(label, "\n")
  cat("Months:", min(df$month), "to", max(df$month), "| n =", nrow(df), "\n\n")

  cat("Sickness rate (per 100k FTE-days available)\n")
  cat("  mean:", fmt_int(sick_mean),
      "| min:", fmt_int(sick_min),
      "| max:", fmt_int(sick_max),
      "| SD:",  fmt_int(sick_sd),
      "| slope/month:", fmt_1(sick_slope), "\n")
  cat("  peak:", format(sick_peak$month, "%Y-%m"), "=", fmt_int(sick_peak$value), "\n\n")

  cat("Test positivity\n")
  cat("  coverage:", pos_present, "/", pos_total,
      "| mean:", fmt_2(pos_mean), "%\n")
  cat("  peak:", format(pos_peak$month, "%Y-%m"), "=", fmt_2(pos_peak$value), "%\n")
  cat("  corr(positivity, sickness) r =", fmt_2(r_pos_sick), "\n\n")

  cat("Vaccination doses (monthly)\n")
  cat("  non-zero months:", dose_nonzero,
      "| mean:", fmt_int(dose_mean),
      "| slope/month:", fmt_int(dose_slope), "\n")
  cat("  peak:", format(dose_peak$month, "%Y-%m"), "=", fmt_int(dose_peak$value), "\n")
  cat("  corr(doses, sickness) r =", fmt_2(r_dose_sick), "\n")
}

# Use df_line as the source (month_totals + total_doses + pos_monthly)
df0 <- df_line %>% filter(month >= as.Date("2019-01-01"), month <= as.Date("2024-12-01"))

p1 <- df0 %>% filter(month >= P1_START, month <= P1_END)
p2 <- df0 %>% filter(month >= P2_START, month <= P2_END)
p3 <- df0 %>% filter(month >= P3_START, month <= P3_END)

summ_period(p1, "B1 -> pre-B2 (2019-09 .. 2021-04)")
summ_period(p2, "B2 -> pre-B3 (2021-05 .. 2022-11)")
summ_period(p3, "Post B3 (2022-12 .. 2024-12)")

# Extra sanity: when do doses stop being non-zero?
last_nonzero <- df0 %>% filter(total_doses > 0) %>% summarise(last = max(month, na.rm = TRUE)) %>% pull(last)
cat("\nLast month with non-zero doses in df_line:", ifelse(is.finite(last_nonzero), format(last_nonzero, "%Y-%m"), "none"), "\n")
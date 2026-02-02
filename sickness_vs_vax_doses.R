# =============================================================================
# sickness_vs_vax_doses.R
# Plot: NHS sickness rate (monthly) vs vaccination doses administered (monthly)
# + adds: COVID test positivity (monthly mean) as a red line (scaled onto left axis)
#
# Inputs:
#   data/sickness_roll6_prediction_band_2020_2026.csv
#     month, sickness_rate_100k, ...
#   data/nhs_staff_vax_model_monthly.csv
#     month, first_dose, second_dose, third_dose
#   data/tests_positivity_and_volume.csv
#     month, positivity_mean_7day, tests_mean_7day_total
#
# Visuals:
#   - Dark grey line: sickness_rate_100k (LEFT axis)
#   - Red line: positivity monthly mean (SCALED to left axis; axis is doses on the right)
#   - Stacked columns: dose1 + dose2 + dose3 (RIGHT axis)
#   - Background shading + breakpoints + event markers kept
#
# Output:
#   nhs_sickness_vs_vax_doses.png
# =============================================================================

library(tidyverse)
library(lubridate)
library(scales)

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

SICK_PATH <- "data/sickness_roll6_prediction_band_2020_2026.csv"
VAX_PATH  <- "data/nhs_staff_vax_model_monthly.csv"
POS_PATH  <- "data/tests_positivity_and_volume.csv"

PLOT_START <- as.Date("2019-01-01")
PLOT_END   <- as.Date("2025-08-01")

# ---- Load sickness (monthly) ----
sick <- read.csv(SICK_PATH, stringsAsFactors = FALSE) |>
  as_tibble() |>
  transmute(
    month = as.Date(paste0(month, "-01")),
    sickness_rate_100k = as.numeric(sickness_rate_100k)
  ) |>
  arrange(month)

# ---- Load vax monthly doses ----
vax <- read.csv(VAX_PATH, stringsAsFactors = FALSE) |>
  as_tibble() |>
  transmute(
    month = as.Date(paste0(month, "-01")),
    first_dose  = as.numeric(first_dose),
    second_dose = as.numeric(second_dose),
    third_dose  = as.numeric(third_dose)
  ) |>
  arrange(month)

# ---- Load monthly positivity mean (no smoothing) ----
pos <- read.csv(POS_PATH, stringsAsFactors = FALSE) |>
  as_tibble() |>
  transmute(
    month = as.Date(paste0(month, "-01")),
    pos_monthly = as.numeric(positivity_mean_7day)
  ) |>
  arrange(month)

# ---- Join + filter window ----
df <- sick |>
  left_join(vax, by = "month") |>
  left_join(pos, by = "month") |>
  filter(month >= PLOT_START, month <= PLOT_END) |>
  arrange(month) |>
  mutate(
    first_dose  = replace_na(first_dose, 0),
    second_dose = replace_na(second_dose, 0),
    third_dose  = replace_na(third_dose, 0),
    total_doses = first_dose + second_dose + third_dose
  )

# =============================================================================
# Axis scaling:
# LEFT axis  = sickness_rate_100k (not forced to 0)
# RIGHT axis = doses (starts at 0)
# Positivity is drawn as red line, scaled onto LEFT axis (since right axis is used).
# =============================================================================

left_min_raw <- min(df$sickness_rate_100k, na.rm = TRUE)
left_max_raw <- max(df$sickness_rate_100k, na.rm = TRUE)

pad <- 0.06 * (left_max_raw - left_min_raw)
left_min <- left_min_raw - pad
left_max <- left_max_raw + pad

max_doses <- max(df$total_doses, na.rm = TRUE)
if (!is.finite(max_doses) || max_doses <= 0) {
  stop("No valid vaccination dose totals found. Check: ", VAX_PATH)
}
k_dose <- (left_max - left_min) / max_doses

df <- df |>
  mutate(
    first_scaled  = left_min + first_dose * k_dose,
    second_scaled = left_min + (first_dose + second_dose) * k_dose,
    third_scaled  = left_min + (first_dose + second_dose + third_dose) * k_dose
  )

# ---- Positivity scaling (onto left axis range) ----
max_pos <- max(df$pos_monthly, na.rm = TRUE)
if (is.finite(max_pos) && max_pos > 0) {
  k_pos <- (left_max - left_min) / max_pos
  df <- df |>
    mutate(pos_scaled = left_min + pos_monthly * k_pos)
} else {
  df <- df |>
    mutate(pos_scaled = NA_real_)
  warning("No valid positivity values found in ", POS_PATH, " (pos line will be omitted).")
}

# =============================================================================
# VISUALIZATION — Substack/print friendly (bigger fonts, cleaner x-axis)
# =============================================================================

# ---- Palette (self-contained) ----
color_sick    <- "grey25"
color_pos     <- "#D55E00"   # positivity line
color_breaks  <- "#CC79A7"
color_markers <- "#D55E00"

# Keep your dose colors (good, colorblind-friendly)
color_d1 <- "#56B4E9"
color_d2 <- "#009E73"
color_d3 <- "#E69F00"

# ---- Breakpoints + labels (single-line labels read better in print) ----
all_breaks <- tibble(
  break_num = 1:3,
  month = as.Date(c("2019-09-01", "2021-05-01", "2022-12-01"))
) |>
  mutate(label_text = sprintf("B%d  %s", break_num, format(month, "%b %Y")))

# ---- Event markers ----
markers <- tibble(
  date = as.Date(c("2020-01-31", "2020-12-08", "2021-11-01", "2022-04-15")),
  short_label = c("COVID", "Vax start", "90% vax & mandates", "Mandates dropped")
)

y_top <- left_max
markers <- markers |>
  mutate(y_label = if_else(row_number() %% 2 == 1, y_top * 0.985, y_top * 0.94))

# ---- Background shading (same as first plot) ----
periods_bg <- tibble(
  period = c("no_tests", "corr_high", "corr_none"),
  xmin   = ymd(c(format(PLOT_START, "%Y-%m-%d"), "2020-05-01", "2022-05-01")),
  xmax   = ymd(c("2020-05-01", "2022-05-01", format(PLOT_END %m+% months(1), "%Y-%m-%d"))),
  fill   = c("#EEEEEE", "#DFF2E1", "#F6DADA")
)
stopifnot(all(!is.na(periods_bg$xmin)), all(!is.na(periods_bg$xmax)))
stopifnot(all(periods_bg$xmax > periods_bg$xmin))

# ---- X-axis: keep only January separators; label every 6 months ----
x_min <- floor_date(min(df$month, na.rm = TRUE), unit = "month")
x_max <- ceiling_date(max(df$month, na.rm = TRUE), unit = "month")
x_breaks_jan <- seq(floor_date(x_min, "year"), floor_date(x_max, "year"), by = "1 year")
grid_jan <- tibble(x = x_breaks_jan)

p <- ggplot(df, aes(x = month)) +

  # --- Background shading for periods (behind everything) ---
  geom_rect(
    data = periods_bg,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = fill),
    inherit.aes = FALSE,
    alpha = 0.70
  ) +
  scale_fill_identity() +

  # --- Year separators only (reduce clutter) ---
  geom_vline(
    data = grid_jan,
    aes(xintercept = x),
    color = "grey85",
    linewidth = 0.7
  ) +

  # --- Stacked vaccination dose columns (scaled to left axis) ---
  # Wider bars for print, but still separated month-to-month
  geom_rect(
    aes(xmin = month - days(14), xmax = month + days(14), ymin = left_min, ymax = first_scaled),
    fill = color_d1, alpha = 0.70
  ) +
  geom_rect(
    aes(xmin = month - days(14), xmax = month + days(14), ymin = first_scaled, ymax = second_scaled),
    fill = color_d2, alpha = 0.70
  ) +
  geom_rect(
    aes(xmin = month - days(14), xmax = month + days(14), ymin = second_scaled, ymax = third_scaled),
    fill = color_d3, alpha = 0.70
  ) +

  # --- Positivity (monthly mean), scaled onto left axis ---
  geom_line(
    aes(y = pos_scaled),
    color = color_pos, linewidth = 1.6, na.rm = TRUE
  ) +

  # --- Sickness line on top ---
  geom_line(
    aes(y = sickness_rate_100k),
    color = color_sick, linewidth = 2.0, na.rm = TRUE
  ) +

  # --- Event markers ---
  geom_vline(
    data = markers,
    aes(xintercept = date),
    linetype = "dashed", linewidth = 1.0,
    color = color_markers, alpha = 0.75
  ) +
  geom_label(
    data = markers,
    aes(x = date, y = y_label, label = short_label),
    color = color_markers,
    fill = alpha("white", 0.92),
    size = 3.6,
    fontface = "bold",
    label.padding = unit(0.18, "lines"),
    label.size = 0.25
  ) +

  # --- Structural breaks + labels ---
  geom_vline(
    data = all_breaks,
    aes(xintercept = month),
    linetype = "longdash", linewidth = 1.2,
    color = color_breaks, alpha = 0.85
  ) +
  geom_label(
    data = all_breaks |>
      mutate(
        y_pos = if_else(row_number() %% 2 == 1,
                        left_min + (left_max - left_min) * 0.36,
                        left_min + (left_max - left_min) * 0.26)
      ),
    aes(x = month, y = y_pos, label = label_text),
    color = color_breaks,
    fill = alpha("white", 0.93),
    size = 3.0,
    fontface = "bold",
    label.padding = unit(0.18, "lines"),
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
    name = "Sickness rate (per 100,000 staff)",
    sec.axis = sec_axis(
      ~ pmax(0, (. - left_min) / k_dose),
      name = "Vaccine doses administered (monthly)",
      breaks = pretty(c(0, max_doses), n = 6),
      labels = label_number(scale_cut = cut_short_scale())
    )
  ) +

  labs(
    title = "NHS staff sickness vs vaccine doses administered",
    subtitle = "Grey: sickness rate | Stacked bars: monthly doses (1st/2nd/3rd) | Red: test positivity (monthly mean, scaled)",
    x = NULL
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
    panel.grid.major.y = element_line(color = "grey88", linewidth = 0.5),
    plot.margin = margin(t = 14, r = 18, b = 18, l = 14),
    legend.position = "none"
  )

print(p)

ggsave(
  "nhs_sickness_vs_vax_doses_substack.png",
  plot = p,
  width = 18,
  height = 10,
  dpi = 320,
  bg = "white"
)

cat("\nPlot saved as 'nhs_sickness_vs_vax_doses_substack.png'\n")


# =============================================================================
# Fiscal-year (Apr-Mar) quarterly aggregation + CSV export
# =============================================================================

# Ensure output folder exists
dir.create("data", showWarnings = FALSE, recursive = TRUE)

df_quarters <- df |>
  mutate(
    # Fiscal year start (April-based)
    fy_start = if_else(month(month) >= 4L, year(month), year(month) - 1L),

    # Label like "2019-20"
    fiscal_year = sprintf("%d-%02d", fy_start, (fy_start + 1L) %% 100L),

    # Fiscal quarters (Apr-Jun = Q1, Jul-Sep = Q2, Oct-Dec = Q3, Jan-Mar = Q4)
    fiscal_quarter = case_when(
      month(month) %in% 4:6   ~ "Q1",
      month(month) %in% 7:9   ~ "Q2",
      month(month) %in% 10:12 ~ "Q3",
      TRUE                    ~ "Q4"
    ),

    q_num = as.integer(str_remove(fiscal_quarter, "Q"))
  ) |>
  group_by(fiscal_year, fiscal_quarter, fy_start, q_num) |>
  summarise(
    positivity_mean_quarter = mean(pos_monthly, na.rm = TRUE),
    doses_sum_quarter       = sum(total_doses, na.rm = TRUE),
    sickness_mean_quarter   = mean(sickness_rate_100k, na.rm = TRUE),
    .groups = "drop"
  ) |>
  arrange(fy_start, q_num) |>
  select(
    fiscal_year,
    fiscal_quarter,
    positivity_mean_quarter,
    doses_sum_quarter,
    sickness_mean_quarter
  )

readr::write_csv(df_quarters, "data/yearly_quarters_positivity_doses_and_sickness.csv")

cat("\nQuarterly fiscal-year CSV saved as 'data/yearly_quarters_positivity_doses_and_sickness.csv'\n")
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
# Markers / breakpoints (kept)
# =============================================================================

all_breaks <- tibble(
  break_num = 1:3,
  month = as.Date(c("2019-09-01", "2021-05-01", "2022-12-01"))
) |>
  mutate(label_text = sprintf("B%d\n%s", break_num, format(month, "%b %Y")))

markers <- tibble(
  date = as.Date(c("2020-01-31", "2020-12-08", "2021-11-01", "2022-04-15")),
  short_label = c("COVID", "Vax start", "90% vax & Mandates", "Drop Mandates")
)

y_top <- left_max
markers <- markers |>
  mutate(y_label = if_else(row_number() %% 2 == 1, y_top * 0.985, y_top * 0.94))

# X-axis: every 3 months + thicker Jan separators
x_min <- floor_date(min(df$month, na.rm = TRUE), unit = "month")
x_max <- ceiling_date(max(df$month, na.rm = TRUE), unit = "month")

x_breaks_3m  <- seq(x_min, x_max, by = "3 months")
x_breaks_jan <- seq(floor_date(x_min, "year"), floor_date(x_max, "year"), by = "1 year")

grid_3m  <- tibble(x = x_breaks_3m)
grid_jan <- tibble(x = x_breaks_jan)

label_3m <- function(d) format(as.Date(d), "%b\n%Y")

# =============================================================================
# Background shading (same as first plot)
# =============================================================================
periods_bg <- tibble(
  period = c("no_tests", "corr_high", "corr_none"),
  xmin   = ymd(c(format(PLOT_START, "%Y-%m-%d"), "2020-05-01", "2022-05-01")),
  xmax   = ymd(c("2020-05-01", "2022-05-01", format(PLOT_END %m+% months(1), "%Y-%m-%d"))),
  fill   = c("#EEEEEE", "#DFF2E1", "#F6DADA")
)

stopifnot(all(!is.na(periods_bg$xmin)), all(!is.na(periods_bg$xmax)))
stopifnot(all(periods_bg$xmax > periods_bg$xmin))

# =============================================================================
# Colors
# =============================================================================
color_sick    <- "grey25"
color_pos     <- "#D55E00"  # red/orange (same as before)
color_breaks  <- "#CC79A7"
color_markers <- "#D55E00"

color_d1 <- "#56B4E9"
color_d2 <- "#009E73"
color_d3 <- "#E69F00"

# =============================================================================
# Plot
# =============================================================================

p <- ggplot(df, aes(x = month)) +
  
  # --- Background shading for periods (behind everything) ---
  geom_rect(
    data = periods_bg,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = fill),
    inherit.aes = FALSE,
    alpha = 0.70
  ) +
  scale_fill_identity() +
  
  # --- Month/year separators FIRST (behind everything else) ---
  geom_vline(
    data = grid_3m,
    aes(xintercept = x),
    color = "grey92", linetype = "dashed", linewidth = 0.3
  ) +
  geom_vline(
    data = grid_jan,
    aes(xintercept = x),
    color = "grey80", linetype = "solid", linewidth = 0.8
  ) +
  
  # --- Stacked vaccination dose columns (scaled to left axis) ---
  geom_rect(
    aes(xmin = month - days(12), xmax = month + days(12), ymin = left_min, ymax = first_scaled),
    fill = color_d1, alpha = 0.75
  ) +
  geom_rect(
    aes(xmin = month - days(12), xmax = month + days(12), ymin = first_scaled, ymax = second_scaled),
    fill = color_d2, alpha = 0.75
  ) +
  geom_rect(
    aes(xmin = month - days(12), xmax = month + days(12), ymin = second_scaled, ymax = third_scaled),
    fill = color_d3, alpha = 0.75
  ) +
  
  # --- Positivity (monthly mean), scaled onto left axis ---
  geom_line(
    aes(y = pos_scaled),
    color = color_pos, linewidth = 1.1, na.rm = TRUE
  ) +
  
  # --- Sickness line on top ---
  geom_line(
    aes(y = sickness_rate_100k),
    color = color_sick, linewidth = 1.4, na.rm = TRUE
  ) +
  
  # --- Event markers (orange dashed vertical lines) ---
  geom_vline(
    data = markers,
    aes(xintercept = date),
    linetype = "dashed", linewidth = 0.9,
    color = color_markers, alpha = 0.75
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
  
  # --- Structural breaks (pink dashed) + labels ---
  geom_vline(
    data = all_breaks,
    aes(xintercept = month),
    linetype = "longdash", linewidth = 1.0,
    color = color_breaks, alpha = 0.85
  ) +
  geom_label(
    data = all_breaks |>
      mutate(
        y_pos = if_else(row_number() %% 2 == 1,
                        left_min + (left_max - left_min) * 0.35,
                        left_min + (left_max - left_min) * 0.25)
      ),
    aes(x = month, y = y_pos, label = label_text),
    color = color_breaks,
    fill = alpha("white", 0.92),
    size = 2.3,
    fontface = "bold",
    label.padding = unit(0.12, "lines"),
    lineheight = 0.9
  ) +
  
  scale_x_date(
    breaks = x_breaks_3m,
    labels = label_3m,
    expand = expansion(mult = c(0.02, 0.05))
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
    title = "NHS staff's sickness vs Vaccine doses administered",
    subtitle = "Dark grey: sickness rate | Stacked bars: monthly doses (1st/2nd/3rd) | Red: test positivity (monthly mean, scaled)",
    x = NULL
  ) +
  
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(size = 16, face = "bold", margin = margin(b = 5)),
    plot.subtitle = element_text(size = 9, color = "grey40", margin = margin(b = 10)),
    axis.text.x = element_text(size = 9, angle = 90, vjust = 0.5, hjust = 1),
    axis.title.y = element_text(size = 11, margin = margin(r = 10)),
    axis.title.y.right = element_text(size = 11, margin = margin(l = 10)),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.major.y = element_line(color = "grey85", linewidth = 0.4),
    legend.position = "none"
  )

print(p)

ggsave(
  "nhs_sickness_vs_vax_doses.png",
  plot = p,
  width = 16,
  height = 9,
  dpi = 300,
  bg = "white"
)

cat("\nPlot saved as 'nhs_sickness_vs_vax_doses.png'\n")


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
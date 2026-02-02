# =============================================================================
# sickness_orgtype_vs_vax_doses.R
# Plot (monthly):
#   1) Sickness rate contributions (FTE Days Sick per 100k FTE Days Available),
#      stacked by Org Type (contributions sum to total monthly rate)
#   2) Vaccine doses administered (monthly) as a line on RIGHT axis
#   3) COVID test positivity (monthly mean) as a red line (scaled onto LEFT axis)
#
# Inputs:
#   data/nhs_sickness_rates_unified_2014_2025.csv
#     Month (YYYY-MM), Org Type, FTE Days Sick, FTE Days Available, ...
#   data/nhs_staff_vax_model_monthly.csv
#     month (YYYY-MM), first_dose, second_dose, third_dose
#   data/tests_positivity_and_volume.csv
#     month (YYYY-MM), positivity_mean_7day, tests_mean_7day_total
#
# Output:
#   nhs_sickness_orgtype_vs_vax_positivity.png
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

SICK_UNIFIED_PATH <- "data/nhs_sickness_rates_unified_2014_2025.csv"
VAX_PATH          <- "data/nhs_staff_vax_model_monthly.csv"
POS_PATH          <- "data/tests_positivity_and_volume.csv"

# Plot window (adjust if needed)
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

# -----------------------------------------------------------------------------
# Load unified sickness + aggregate monthly by Org Type
# -----------------------------------------------------------------------------
sick_raw <- readr::read_csv(SICK_UNIFIED_PATH, show_col_types = FALSE)

# Expect Month in YYYY-MM; fall back to Date if Month missing
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

# Monthly totals for denominator and total rate
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

# Contribution of each Org Type to the total monthly rate:
# contrib_rate_100k = (sick_org / total_available_all) * 100k
sick_stack <- sick_by_type |>
  left_join(month_totals |> select(month, avail_total), by = "month") |>
  mutate(
    contrib_rate_100k = if_else(avail_total > 0, sick_sum / avail_total * 1e5, NA_real_)
  ) |>
  filter(!is.na(contrib_rate_100k))

# Order Org Types by overall mean contribution (stable legend order)
org_order <- sick_stack |>
  group_by(org_type) |>
  summarise(m = mean(contrib_rate_100k, na.rm = TRUE), .groups = "drop") |>
  arrange(desc(m)) |>
  pull(org_type)

sick_stack <- sick_stack |>
  mutate(org_type = factor(org_type, levels = org_order))

# -----------------------------------------------------------------------------
# Load vax monthly doses
# -----------------------------------------------------------------------------
vax <- readr::read_csv(VAX_PATH, show_col_types = FALSE) |>
  transmute(
    month = as.Date(paste0(month, "-01")),
    first_dose  = to_num(first_dose),
    second_dose = to_num(second_dose),
    third_dose  = to_num(third_dose)
  ) |>
  arrange(month) |>
  mutate(
    first_dose  = replace_na(first_dose, 0),
    second_dose = replace_na(second_dose, 0),
    third_dose  = replace_na(third_dose, 0),
    total_doses = first_dose + second_dose + third_dose
  )

# -----------------------------------------------------------------------------
# Load positivity (monthly mean)
# -----------------------------------------------------------------------------
pos <- readr::read_csv(POS_PATH, show_col_types = FALSE) |>
  transmute(
    month = as.Date(paste0(month, "-01")),
    pos_monthly = to_num(positivity_mean_7day)
  ) |>
  arrange(month)

# -----------------------------------------------------------------------------
# Join lines dataframe + filter to plot window
# -----------------------------------------------------------------------------
df_line <- month_totals |>
  left_join(vax |> select(month, total_doses), by = "month") |>
  left_join(pos, by = "month") |>
  filter(month >= PLOT_START, month <= PLOT_END) |>
  arrange(month) |>
  mutate(
    total_doses = replace_na(total_doses, 0)
  )

# -----------------------------------------------------------------------------
# Axis scaling
# LEFT axis  = sickness rates (stacked contributions sum to total_rate_100k)
# RIGHT axis = doses (scaled to left)
# Positivity scaled to LEFT axis
# -----------------------------------------------------------------------------
left_min <- 0

left_max_raw <- max(df_line$total_rate_100k, na.rm = TRUE)
if (!is.finite(left_max_raw) || left_max_raw <= 0) {
  stop("No valid sickness rate found after filtering. Check: ", SICK_UNIFIED_PATH)
}
left_max <- left_max_raw * 1.10  # 10% headroom

max_doses <- max(df_line$total_doses, na.rm = TRUE)
k_dose <- if (is.finite(max_doses) && max_doses > 0) left_max / max_doses else NA_real_

df_line <- df_line |>
  mutate(doses_scaled = if (is.finite(k_dose)) total_doses * k_dose else NA_real_)

max_pos <- max(df_line$pos_monthly, na.rm = TRUE)
k_pos <- if (is.finite(max_pos) && max_pos > 0) left_max / max_pos else NA_real_

df_line <- df_line |>
  mutate(pos_scaled = if (is.finite(k_pos)) pos_monthly * k_pos else NA_real_)

# -----------------------------------------------------------------------------
# Markers / breakpoints (kept from example)
# -----------------------------------------------------------------------------
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

# X-axis: every 3 months + thicker Jan separators, labels YYYY-MM
x_min <- floor_date(min(df_line$month, na.rm = TRUE), unit = "month")
x_max <- ceiling_date(max(df_line$month, na.rm = TRUE), unit = "month")

x_breaks_3m  <- seq(x_min, x_max, by = "3 months")
x_breaks_jan <- seq(floor_date(x_min, "year"), floor_date(x_max, "year"), by = "1 year")

grid_3m  <- tibble(x = x_breaks_3m)
grid_jan <- tibble(x = x_breaks_jan)

label_3m <- function(d) format(as.Date(d), "%Y-%m")

# -----------------------------------------------------------------------------
# Background shading (same logic as example)
# (It’s OK if your plot starts earlier; shading still applies across full range.)
# -----------------------------------------------------------------------------
periods_bg <- tibble(
  period = c("no_tests", "corr_high", "corr_none"),
  xmin   = ymd(c(format(PLOT_START, "%Y-%m-%d"), "2020-05-01", "2022-05-01")),
  xmax   = ymd(c("2020-05-01", "2022-05-01", format(PLOT_END %m+% months(1), "%Y-%m-%d"))),
  fill   = c("#EEEEEE", "#DFF2E1", "#F6DADA")
)

stopifnot(all(!is.na(periods_bg$xmin)), all(!is.na(periods_bg$xmax)))
stopifnot(all(periods_bg$xmax > periods_bg$xmin))

# -----------------------------------------------------------------------------
# Plot
# -----------------------------------------------------------------------------
color_pos     <- "#D55E00"
color_breaks  <- "#CC79A7"
color_markers <- "#D55E00"
color_doses   <- "#0072B2"  # blue line for doses

p <- ggplot() +

  # --- Background shading (behind everything) ---
  geom_rect(
    data = periods_bg,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
    fill = I(periods_bg$fill),
    inherit.aes = FALSE,
    alpha = 0.70
  ) +

  # --- Month/year separators (behind bars/lines) ---
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

  # --- Stacked columns: sickness rate contributions by Org Type ---
  geom_col(
    data = sick_stack,
    aes(x = month, y = contrib_rate_100k, fill = org_type),
    width = 25,
    alpha = 0.85
  ) +
  scale_fill_brewer(palette = "Set3") +

  # --- Positivity (monthly mean), scaled onto left axis ---
  geom_line(
    data = df_line,
    aes(x = month, y = pos_scaled),
    color = color_pos, linewidth = 1.1, na.rm = TRUE
  ) +

  # --- Vaccine doses (monthly), scaled onto left axis (RIGHT axis shows real doses) ---
  geom_line(
    data = df_line,
    aes(x = month, y = doses_scaled),
    color = color_doses, linewidth = 1.2, na.rm = TRUE
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
                        left_max * 0.35,
                        left_max * 0.25)
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
    name = "Sickness rate (FTE Days Sick per 100k FTE Days Available)\n(stacked by Org Type contribution)",
    sec.axis = sec_axis(
      ~ if (is.finite(k_dose)) pmax(0, . / k_dose) else NA_real_,
      name = "Vaccine doses administered (monthly)",
      breaks = if (is.finite(max_doses) && max_doses > 0) pretty(c(0, max_doses), n = 6) else waiver(),
      labels = label_number(scale_cut = cut_short_scale())
    )
  ) +

  labs(
    title = "NHS sickness (by Org Type) vs Vaccine doses and Test positivity",
    subtitle = "Stacked bars: sickness rate contributions by Org Type | Blue: monthly doses (scaled; see right axis) | Red: test positivity (scaled)",
    x = NULL,
    fill = "Org Type"
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
    legend.position = "bottom",
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 9),
    legend.key.height = unit(0.35, "cm")
  )

print(p)

ggsave(
  "nhs_sickness_orgtype_vs_vax_positivity.png",
  plot = p,
  width = 16,
  height = 9,
  dpi = 300,
  bg = "white"
)

cat("\nPlot saved as 'nhs_sickness_orgtype_vs_vax_positivity.png'\n")

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
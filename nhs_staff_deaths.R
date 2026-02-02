# ---- Packages ----
library(readr)
library(dplyr)
library(stringr)
library(ggplot2)
library(lubridate)
library(scales)

# ---- 1) Load files ----
deaths <- read_csv("data/deaths_in_service_long.csv",
                   col_types = cols(
                     year_from_to = col_character(),
                     quarter      = col_character(),  # "Q1".."Q4"
                     deaths       = col_double()
                   ))

denoms <- read_csv("data/april_denominators_2011_2025.csv",
                   col_types = cols(
                     period      = col_character(),   # "201004 to 201104"
                     denominator = col_double()
                   ))

qs <- read_csv("data/yearly_quarters_positivity_doses_and_sickness.csv",
               col_types = cols(
                 fiscal_year = col_character(),
                 fiscal_quarter = col_character(),
                 positivity_mean_quarter = col_double(),
                 doses_sum_quarter = col_double(),
                 sickness_mean_quarter = col_double()
               ))

# Helper: convert "YYYY-YY" -> end year as integer (2012, 2020)
year_from_to_to_end_year <- function(x) {
  end_two  <- as.integer(str_sub(x, 6, 7))
  if_else(end_two < 50L, 2000L + end_two, 1900L + end_two)
}

# Helper: "Q1".."Q4" -> 1..4
quarter_to_num <- function(q) {
  as.integer(str_remove(str_to_upper(q), "Q"))
}

# Helper: get quarter start date for an NHS financial-year quarter
# Q1=Apr (start FY), Q2=Jul, Q3=Oct of (year_end-1), Q4=Jan of (year_end)
quarter_start_date <- function(year_end, q_num) {
  start_year <- if_else(q_num == 4L, year_end, year_end - 1L)
  start_month <- case_when(
    q_num == 1L ~ 4L,
    q_num == 2L ~ 7L,
    q_num == 3L ~ 10L,
    TRUE        ~ 1L
  )
  as.Date(sprintf("%04d-%02d-01", start_year, start_month))
}

# ---- 2) Build quarterly deaths dataset ----
q_dat <- deaths %>%
  mutate(
    year_end = year_from_to_to_end_year(year_from_to),
    q_num    = quarter_to_num(quarter),
    
    # Denominator key: FY April->April
    period   = sprintf("%d04 to %d04", year_end - 1L, year_end),
    
    # Regression time (financial-year-based)
    time_num = (year_end - 1L) + (q_num - 1L) / 4,
    
    # Quarter start date for plotting
    q_start  = quarter_start_date(year_end, q_num)
  )

# ---- 3) Join denominators + quarterly positivity/doses/sickness ----
dat <- q_dat %>%
  left_join(denoms, by = "period") %>%
  mutate(
    rate_per_100k = (deaths / denominator) * 100000
  ) %>%
  left_join(
    qs,
    by = c("year_from_to" = "fiscal_year", "quarter" = "fiscal_quarter")
  ) %>%
  arrange(year_end, q_num) %>%
  mutate(
    period_q = paste0(year_from_to, " ", quarter),
    period_q = factor(period_q, levels = unique(period_q))
  )

# ---- Fit linear trend on 2011–2019 (quarterly points) ----
trend_fit <- lm(
  rate_per_100k ~ time_num,
  data = dat %>% filter(year_end >= 2011, year_end <= 2019)
)

pred_df <- dat %>%
  distinct(period_q, q_start, year_end, q_num, time_num) %>%
  arrange(year_end, q_num) %>%
  mutate(pred_rate = predict(trend_fit, newdata = .))

dat2 <- dat %>%
  left_join(pred_df %>% select(period_q, pred_rate), by = "period_q") %>%
  mutate(excess_pct = (rate_per_100k - pred_rate) / pred_rate * 100)

# Labels (keep as before; you can switch to Q4-only if too busy)
excess_labels <- dat2 %>%
  filter(year_end >= 2021, year_end <= 2025) %>%
  mutate(label = sprintf("%+.1f%%", excess_pct))

pred_df <- pred_df %>%
  mutate(period_q = factor(period_q, levels = levels(dat2$period_q)))

pred_solid  <- pred_df %>% filter(year_end <= 2019)
pred_dashed <- pred_df %>% filter(year_end >= 2020)

# =============================================================================
# Axis scaling (aligned to the monthly script style)
# LEFT axis  = deaths per 100k
# RIGHT axis = doses_sum_quarter (starts at 0), drawn as bars scaled to left axis
# Positivity & sickness are drawn as lines scaled onto left axis
# =============================================================================

left_min_raw <- min(dat2$rate_per_100k, na.rm = TRUE)
left_max_raw <- max(dat2$rate_per_100k, na.rm = TRUE)
pad <- 0.06 * (left_max_raw - left_min_raw)
left_min <- left_min_raw - pad
left_max <- left_max_raw + pad

max_doses_q <- max(dat2$doses_sum_quarter, na.rm = TRUE)
if (!is.finite(max_doses_q) || max_doses_q <= 0) {
  max_doses_q <- NA_real_
}
k_dose <- if (!is.na(max_doses_q)) (left_max - left_min) / max_doses_q else NA_real_

max_pos_q <- max(dat2$positivity_mean_quarter, na.rm = TRUE)
k_pos <- if (is.finite(max_pos_q) && max_pos_q > 0) (left_max - left_min) / max_pos_q else NA_real_

max_sick_q <- max(dat2$sickness_mean_quarter, na.rm = TRUE)
k_sick <- if (is.finite(max_sick_q) && max_sick_q > 0) (left_max - left_min) / max_sick_q else NA_real_

dat2 <- dat2 %>%
  mutate(
    # bars for doses (scaled to left axis)
    doses_scaled = if (!is.na(k_dose)) left_min + doses_sum_quarter * k_dose else NA_real_,
    
    # red positivity line (scaled)
    pos_scaled = if (!is.na(k_pos)) left_min + positivity_mean_quarter * k_pos else NA_real_,
    
    # optional sickness line (scaled)
    sick_scaled = if (!is.na(k_sick)) left_min + sickness_mean_quarter * k_sick else NA_real_
  )

# =============================================================================
# Background shading + vertical separators (quarterly)
# =============================================================================

PLOT_START <- as.Date(min(dat2$q_start, na.rm = TRUE))
PLOT_END   <- as.Date(max(dat2$q_start, na.rm = TRUE))

# periods like the monthly script
periods_bg <- tibble(
  period = c("no_tests", "corr_high", "corr_none"),
  xmin   = ymd(c(format(PLOT_START, "%Y-%m-%d"), "2020-05-01", "2022-05-01")),
  xmax   = ymd(c("2020-05-01", "2022-05-01", format(PLOT_END %m+% months(4), "%Y-%m-%d"))),
  fill   = c("#EEEEEE", "#DFF2E1", "#F6DADA")
)

# vertical grid: every quarter + thicker April separators (FY start)
x_min <- floor_date(PLOT_START, unit = "month")
x_max <- ceiling_date(PLOT_END %m+% months(3), unit = "month")

x_breaks_q  <- seq(x_min, x_max, by = "3 months")     # quarter starts
x_breaks_apr <- seq(ymd(sprintf("%d-04-01", year(x_min) - 1L)),
                    ymd(sprintf("%d-04-01", year(x_max) + 1L)),
                    by = "1 year")

grid_q   <- tibble(x = x_breaks_q)
grid_apr <- tibble(x = x_breaks_apr)

label_q <- function(d) {
  # label each quarter start as "Qx\nYYYY" (FY quarter)
  # simple: show month+year like the other plot
  format(as.Date(d), "%b\n%Y")
}

# =============================================================================
# Colors (same palette logic)
# =============================================================================
color_deaths <- "grey25"
color_pos    <- "#D55E00"  # red/orange
color_sick   <- "grey40"

color_dose_bar <- "#56B4E9"  # single dose bar color

# =============================================================================
# PLOT — Substack/print friendly (larger fonts, cleaner x-axis, less clutter)
# =============================================================================

# ---- Colors (self-contained) ----
color_deaths   <- "grey25"
color_pos      <- "#D55E00"   # positivity (scaled)
color_sick     <- "grey45"    # sickness (scaled, dashed)
color_dose_bar <- "#56B4E9"   # quarterly dose bars
color_trend    <- "#D55E00"   # trend line (keep same family as positivity)

# ---- Excess labels: reduce clutter by labeling Q4 only (edit if you prefer) ----
excess_labels <- dat2 %>%
  filter(year_end >= 2021, year_end <= 2025, q_num == 4L) %>%
  mutate(label = sprintf("%+.1f%%", excess_pct))

# ---- Background shading (same periods) ----
PLOT_START <- as.Date(min(dat2$q_start, na.rm = TRUE))
PLOT_END   <- as.Date(max(dat2$q_start, na.rm = TRUE))

periods_bg <- tibble(
  period = c("no_tests", "corr_high", "corr_none"),
  xmin   = ymd(c(format(PLOT_START, "%Y-%m-%d"), "2020-05-01", "2022-05-01")),
  xmax   = ymd(c("2020-05-01", "2022-05-01", format(PLOT_END %m+% months(4), "%Y-%m-%d"))),
  fill   = c("#EEEEEE", "#DFF2E1", "#F6DADA")
)
stopifnot(all(!is.na(periods_bg$xmin)), all(!is.na(periods_bg$xmax)))
stopifnot(all(periods_bg$xmax > periods_bg$xmin))

# ---- X-axis: keep only FY start separators (April); label every 6 months ----
x_min <- floor_date(PLOT_START, unit = "month")
x_max <- ceiling_date(PLOT_END %m+% months(3), unit = "month")

x_breaks_apr <- seq(
  ymd(sprintf("%d-04-01", year(x_min) - 1L)),
  ymd(sprintf("%d-04-01", year(x_max) + 1L)),
  by = "1 year"
)
grid_apr <- tibble(x = x_breaks_apr)

p <- ggplot(dat2, aes(x = q_start)) +

  # --- Background shading ---
  geom_rect(
    data = periods_bg,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = fill),
    inherit.aes = FALSE,
    alpha = 0.70
  ) +
  scale_fill_identity() +

  # --- FY (April) separators only ---
  geom_vline(
    data = grid_apr,
    aes(xintercept = x),
    color = "grey85",
    linewidth = 0.7
  ) +

  # --- Doses bars (scaled onto left axis) ---
  geom_rect(
    aes(
      xmin = q_start - days(38), xmax = q_start + days(38),
      ymin = left_min, ymax = doses_scaled
    ),
    fill = color_dose_bar, alpha = 0.55, na.rm = TRUE
  ) +

  # --- Positivity scaled (red) ---
  geom_line(
    aes(y = pos_scaled),
    color = color_pos, linewidth = 1.6, na.rm = TRUE
  ) +

  # --- Optional sickness scaled (grey dashed) ---
  geom_line(
    aes(y = sick_scaled),
    color = color_sick, linewidth = 1.2, linetype = "dashed", na.rm = TRUE
  ) +

  # --- Deaths line + points ---
  geom_line(
    aes(y = rate_per_100k),
    color = color_deaths, linewidth = 2.0, na.rm = TRUE
  ) +
  geom_point(
    aes(y = rate_per_100k),
    size = 2.6, color = color_deaths, na.rm = TRUE
  ) +

  # --- Trend: solid 2011–2019, dashed 2020+ (same color family) ---
  geom_line(
    data = pred_solid,
    aes(x = q_start, y = pred_rate, group = 1),
    color = color_trend,
    linewidth = 1.5
  ) +
  geom_line(
    data = pred_dashed,
    aes(x = q_start, y = pred_rate, group = 1),
    color = color_trend,
    linetype = "dashed",
    linewidth = 1.5
  ) +

  # --- Excess labels (Q4 only) ---
  geom_label(
    data = excess_labels,
    aes(x = q_start, y = rate_per_100k, label = label),
    fill = alpha("white", 0.90),
    color = "grey20",
    label.size = 0.25,
    size = 3.2,
    fontface = "bold",
    nudge_y = 0.15
  ) +

  scale_x_date(
    date_breaks = "6 months",
    date_labels = "%b\n%Y",
    expand = expansion(mult = c(0.01, 0.03))
  ) +

  scale_y_continuous(
    limits = c(left_min, left_max),
    labels = label_comma(),
    name = "Deaths per 100,000 (quarterly; normalised by April denominator)",
    sec.axis = sec_axis(
      ~ if (!is.na(k_dose)) pmax(0, (. - left_min) / k_dose) else .,
      name = "Vaccine doses administered (quarterly sum)",
      breaks = if (is.na(max_doses_q)) waiver() else pretty(c(0, max_doses_q), n = 6),
      labels = label_number(scale_cut = cut_short_scale())
    )
  ) +

  labs(
    title = "NHS staff: quarterly deaths vs doses, positivity, and sickness",
    subtitle = "Grey: deaths | Blue bars: doses (quarterly sum) | Orange: positivity (scaled) | Grey dashed: sickness (scaled) | Trend: 2011–2019 (solid), 2020+ (dashed)",
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
  "nhs_staff_quarterly_deaths_substack.png",
  plot = p,
  width = 18,
  height = 10,
  dpi = 320,
  bg = "white"
)

cat("\nPlot saved as 'nhs_staff_quarterly_deaths_substack.png'\n")
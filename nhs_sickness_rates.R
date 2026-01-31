# ============================================================================
# NHS Sickness Rate Analysis - Improved Visualization & Trend Detection
# ============================================================================

# ---- Packages ----
# install.packages(c("tidyverse", "lubridate", "slider", "scales", "ggrepel"))
library(tidyverse)
library(lubridate)
library(slider)
library(scales)
library(ggrepel)  # For non-overlapping labels

# ---- Read data ----
df <- read.csv("data/overall_sickness_rates.csv", stringsAsFactors = FALSE) |>
  as_tibble() |>
  rename(rate_pct = England) |>
  mutate(
    month     = my(Month),
    rate_pct  = as.numeric(rate_pct)
  ) |>
  arrange(month) |>
  mutate(
    # Convert percent -> per 100,000 (multiply by 1000)
    rate_100k = rate_pct * 1000
  )

# ---- 6-month rolling average: (t-3, t-2, t-1, t, t+1, t+2) ----
df <- df |>
  mutate(
    roll6_100k = slide_dbl(rate_100k, mean, .before = 3, .after = 2, .complete = TRUE),
    month_num  = as.numeric(month)
  )

# ---- Linear trend fit on rolling average for 2010-2019 ----
df_fit <- df |>
  filter(year(month) >= 2010, year(month) <= 2019, !is.na(roll6_100k))

fit_2010_2019 <- lm(roll6_100k ~ month_num, data = df_fit)

# ---- Predictions + confidence intervals (2010-2019 trend) ----
pred_all <- as_tibble(predict(fit_2010_2019, newdata = df, interval = "confidence", level = 0.95)) |>
  rename(trend_100k = fit, ci_low_100k = lwr, ci_high_100k = upr)

df <- bind_cols(df, pred_all) |>
  mutate(
    period = case_when(
      year(month) >= 2010 & year(month) <= 2019 ~ "fit_2010_2019",
      year(month) >= 2020 & year(month) <= 2025 ~ "proj_2020_2025",
      TRUE ~ "other"
    ),
    dev_pct = if_else(
      year(month) >= 2020 & year(month) <= 2025 & !is.na(roll6_100k) & !is.na(trend_100k),
      100 * (roll6_100k - trend_100k) / trend_100k,
      NA_real_
    ),
    dev_label = if_else(!is.na(dev_pct), sprintf("%+.1f%%", dev_pct), NA_character_)
  )

# ============================================================================
# REQUIREMENT 2: Detect Exact Trend Inflection Point
# ============================================================================
# 
# Strategy: Find the earliest month (starting from 2020) such that removing
# earlier months does NOT decrease the post-COVID trend slope. We iterate
# backward from the latest data to find when the trend stabilizes.

cat("\n========================================\n")
cat("TREND INFLECTION ANALYSIS\n")
cat("========================================\n\n")

# Get post-2020 data with rolling average
df_post_covid <- df |>
  filter(month >= as.Date("2020-01-01"), !is.na(roll6_100k)) |>
  arrange(month)

# Function to calculate trend slope from a given start month
calc_slope <- function(start_month, data) {
  sub_df <- data |> filter(month >= start_month)
  if (nrow(sub_df) < 6) return(NA_real_)  # Need minimum observations
  fit <- lm(roll6_100k ~ month_num, data = sub_df)
  coef(fit)[2]  # Return slope
}

# Test different start months (from Jan 2020 forward)
candidate_months <- df_post_covid$month
slope_results <- tibble(
  start_month = candidate_months,
  slope = map_dbl(candidate_months, ~calc_slope(.x, df_post_covid))
) |>
  filter(!is.na(slope))

# Find the inflection point: where adding earlier months starts to DECREASE the slope
# (i.e., the earliest month that doesn't drag down the post-COVID upward trend)
slope_results <- slope_results |>
  arrange(start_month) |>
  mutate(
    slope_change = slope - lead(slope),  # Positive if earlier month increases slope
    months_included = n() - row_number() + 1
  )

# The inflection is where including earlier months stops helping the trend
# We look for where the slope stabilizes or starts declining
inflection_analysis <- slope_results |>
  filter(start_month >= as.Date("2020-01-01"), start_month <= as.Date("2022-06-01"))

cat("Slope analysis by starting month:\n")
print(inflection_analysis |> select(start_month, slope, months_included), n = 30)

# Find the optimal inflection point:
# Method: Find where the cumulative trend from that point forward is maximally positive
# while including as many months as possible (balance between signal strength and inclusion)

# Calculate a score: slope * sqrt(months_included) to balance both
inflection_analysis <- inflection_analysis |>
  mutate(
    score = slope * sqrt(months_included),
    slope_per_year = slope * 365.25  # Annualized slope
  )

# Find the earliest month where the slope is at least 95% of the maximum post-COVID slope
max_slope <- max(inflection_analysis$slope, na.rm = TRUE)
threshold_slope <- max_slope * 0.95

inflection_month <- inflection_analysis |>
  filter(slope >= threshold_slope) |>
  arrange(start_month) |>
  slice(1) |>
  pull(start_month)

cat("\n----------------------------------------\n")
cat(sprintf("Maximum post-COVID slope: %.4f per day (%.1f per year)\n", 
            max_slope, max_slope * 365.25))
cat(sprintf("95%% threshold slope: %.4f per day\n", threshold_slope))
cat(sprintf("\n>>> INFLECTION POINT DETECTED: %s <<<\n", format(inflection_month, "%B %Y")))
cat("----------------------------------------\n\n")

# Alternative method: Bai-Perron style - find structural break
# Using a simpler approach: find the month where residuals from 2010-2019 trend
# consistently become positive

df <- df |>
  mutate(
    residual = roll6_100k - trend_100k,
    above_trend = residual > 0
  )

# Find first month of sustained above-trend performance
sustained_above <- df |>
  filter(month >= as.Date("2020-01-01"), !is.na(roll6_100k)) |>
  mutate(
    # Rolling count of consecutive months above trend
    consec_above = slide_dbl(above_trend, sum, .before = 5, .after = 0, .complete = TRUE)
  ) |>
  filter(consec_above >= 6) |>  # At least 6 consecutive months above trend
  slice(1)

if (nrow(sustained_above) > 0) {
  structural_break <- sustained_above$month %m-% months(5)  # Go back 5 months
  cat(sprintf("Structural break (6+ consecutive months above trend): %s\n", 
              format(structural_break, "%B %Y")))
}

# ============================================================================
# REQUIREMENT 3: November 2021+ Trend Line
# ============================================================================

df_nov2021_plus <- df |>
  filter(month >= as.Date("2021-11-01"), !is.na(roll6_100k))

fit_nov2021 <- lm(roll6_100k ~ month_num, data = df_nov2021_plus)

# Add predictions for the Nov 2021+ trend
df <- df |>
  mutate(
    trend_nov2021 = if_else(
      month >= as.Date("2021-11-01"),
      predict(fit_nov2021, newdata = pick(everything())),
      NA_real_
    )
  )

# Report the Nov 2021+ trend
slope_nov2021 <- coef(fit_nov2021)[2] * 365.25  # Annualized
slope_2010_2019 <- coef(fit_2010_2019)[2] * 365.25

cat("\n========================================\n")
cat("TREND COMPARISON\n")
cat("========================================\n")
cat(sprintf("2010-2019 trend: %+.1f per 100k per year\n", slope_2010_2019))
cat(sprintf("Nov 2021+ trend: %+.1f per 100k per year\n", slope_nov2021))
cat(sprintf("Difference: %+.1f per 100k per year (%.1fx steeper)\n", 
            slope_nov2021 - slope_2010_2019,
            slope_nov2021 / slope_2010_2019))
cat("========================================\n\n")

# ============================================================================
# REQUIREMENT 1: Improved Visualization
# ============================================================================

# ---- Prepare label data (less frequent, better positioned) ----
df_labels <- df |>
  filter(period == "proj_2020_2025", !is.na(dev_label), !is.na(roll6_100k)) |>
  mutate(idx = row_number()) |>
  filter((idx - 1) %% 12 == 0)  # Every 12th point (annual) for cleaner look

# ---- Marker dates with improved styling ----
markers <- tibble(
  date = as.Date(c("2020-01-31", "2020-12-08", "2021-11-01", "2022-04-15")),
  label = c(
    "COVID\nJan 2020",
    "First vaccinations\nDec 2020",
    "90% NHS vaccinated\nNov 2021",
    "Fully vaccinated\nApr 2022"
  ),
  short_label = c("COVID", "Vax start", "90% vax", "Full vax")
)

# Calculate y-positions for markers (staggered)
y_range <- range(c(df$rate_100k, df$roll6_100k, df$trend_100k), na.rm = TRUE)
y_top <- y_range[2]
y_bottom <- y_range[1]

markers <- markers |>
  mutate(
    y_line_top = y_top * 1.02,
    y_label = if_else(row_number() %% 2 == 1, y_top * 0.98, y_top * 0.90)
  )

# ---- Add inflection point marker ----
inflection_marker <- tibble(
  date = inflection_month,
  label = sprintf("Trend inflection\n%s", format(inflection_month, "%b %Y"))
)

# ---- Color palette (improved readability) ----
color_raw <- "grey50"
color_rolling <- "#E69F00"  # Orange (colorblind-friendly)
color_trend_2010 <- "#0072B2"  # Blue (colorblind-friendly)
color_trend_nov21 <- "#009E73"  # Teal/Green (colorblind-friendly)
color_markers <- "#D55E00"  # Vermillion (colorblind-friendly)
color_inflection <- "#CC79A7"  # Pink (colorblind-friendly)
color_deviation <- "#882255"  # Dark magenta for deviation labels

# ---- Create the improved plot ----
p <- ggplot(df, aes(x = month)) +
  
  # Confidence band for 2010-2019 projection
  geom_ribbon(
    data = df |> filter(period == "proj_2020_2025"),
    aes(ymin = ci_low_100k, ymax = ci_high_100k),
    fill = color_trend_2010, alpha = 0.15
  ) +
  
  # Raw monthly data (subdued)
  geom_line(aes(y = rate_100k), color = color_raw, linewidth = 0.8, alpha = 0.6) +
  
  # 6-month rolling average (prominent)
  geom_line(aes(y = roll6_100k), color = color_rolling, linewidth = 1.5, na.rm = TRUE) +
  
  # 2010-2019 trend (solid for fit period)
  geom_line(
    data = df |> filter(period == "fit_2010_2019"),
    aes(y = trend_100k),
    color = color_trend_2010, linewidth = 1.4
  ) +
  
  # 2010-2019 projection (dashed for projection period)
  geom_line(
    data = df |> filter(period == "proj_2020_2025"),
    aes(y = trend_100k),
    color = color_trend_2010, linewidth = 1.4, linetype = "dashed"
  ) +
  
  # November 2021+ trend (new - REQUIREMENT 3)
  geom_line(
    data = df |> filter(!is.na(trend_nov2021)),
    aes(y = trend_nov2021),
    color = color_trend_nov21, linewidth = 1.6, linetype = "solid"
  ) +
  
  # Deviation labels (improved - REQUIREMENT 1)
  geom_label_repel(
    data = df_labels,
    aes(y = roll6_100k, label = dev_label),
    color = color_deviation,
    fill = "white",
    size = 3.8,
    fontface = "bold",
    label.padding = unit(0.2, "lines"),
    box.padding = unit(0.4, "lines"),
    point.padding = unit(0.3, "lines"),
    segment.color = color_deviation,
    segment.size = 0.4,
    direction = "y",
    nudge_y = 150,
    min.segment.length = 0,
    max.overlaps = 20
  ) +
  
  # Event markers (improved styling - REQUIREMENT 1)
  geom_vline(
    data = markers,
    aes(xintercept = date),
    linetype = "dotted", linewidth = 0.9,
    color = color_markers, alpha = 0.8
  ) +
  
  # Event labels (positioned at top, cleaner)
  geom_label(
    data = markers,
    aes(x = date, y = y_label, label = short_label),
    color = color_markers,
    fill = alpha("white", 0.9),
    size = 3.2,
    fontface = "bold",
    label.padding = unit(0.15, "lines"),
    label.size = 0.3
  ) +
  
  # Inflection point marker (REQUIREMENT 2)
  geom_vline(
    xintercept = inflection_month,
    linetype = "longdash", linewidth = 1.2,
    color = color_inflection
  ) +
  annotate(
    "label",
    x = inflection_month,
    y = y_bottom + (y_top - y_bottom) * 0.15,
    label = sprintf("Trend inflection\n%s", format(inflection_month, "%b %Y")),
    color = color_inflection,
    fill = alpha("white", 0.9),
    size = 3.5,
    fontface = "bold",
    hjust = 0.5
  ) +
  
  # Scales
  scale_x_date(
    date_breaks = "2 years",
    date_labels = "%Y",
    expand = expansion(mult = c(0.02, 0.05))
  ) +
  scale_y_continuous(
    limits = c(2000, NA),
    labels = label_comma(),
    expand = expansion(mult = c(0.02, 0.08))
  ) +
  
  # Labels
  labs(
    title = "NHS Sickness Absence Rate (England)",
    subtitle = paste0(
      "Monthly rate (grey) | 6-month rolling avg (orange) | ",
      "2010–2019 trend & projection (blue) | Nov 2021+ trend (green)"
    ),
    x = NULL,
    y = "Sickness rate (per 100,000 staff)",
    caption = paste0(
      "Data source: NHS England | ",
      sprintf("2010-2019 trend: %+.0f/yr | Nov 2021+ trend: %+.0f/yr", 
              slope_2010_2019, slope_nov2021)
    )
  ) +
  
  # Theme (clean, professional)
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(size = 18, face = "bold", margin = margin(b = 5)),
    plot.subtitle = element_text(size = 11, color = "grey40", margin = margin(b = 15)),
    plot.caption = element_text(size = 9, color = "grey50", hjust = 0, margin = margin(t = 15)),
    axis.text = element_text(size = 11),
    axis.title.y = element_text(size = 12, margin = margin(r = 10)),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(color = "grey90", linewidth = 0.3),
    panel.grid.major.y = element_line(color = "grey85", linewidth = 0.4),
    plot.margin = margin(t = 15, r = 20, b = 10, l = 10),
    legend.position = "none"
  )

# ---- Display plot ----
print(p)

# ---- Save high-resolution output ----
ggsave(
  "nhs_sickness_analysis_improved.png",
  plot = p,
  width = 14,
  height = 8,
  dpi = 300,
  bg = "white"
)

cat("\nPlot saved as 'nhs_sickness_analysis_improved.png'\n")

# ============================================================================
# SUMMARY OUTPUT
# ============================================================================

cat("\n")
cat("╔════════════════════════════════════════════════════════════════════╗\n")
cat("║                         ANALYSIS SUMMARY                          ║\n")
cat("╠════════════════════════════════════════════════════════════════════╣\n")
cat(sprintf("║  Inflection point detected:  %-37s║\n", format(inflection_month, "%B %Y")))
cat(sprintf("║  2010-2019 trend slope:      %+.1f per 100k per year            ║\n", slope_2010_2019))
cat(sprintf("║  Nov 2021+ trend slope:      %+.1f per 100k per year           ║\n", slope_nov2021))
cat(sprintf("║  Post-COVID acceleration:    %.1fx steeper                      ║\n", slope_nov2021/slope_2010_2019))
cat("╚════════════════════════════════════════════════════════════════════╝\n")
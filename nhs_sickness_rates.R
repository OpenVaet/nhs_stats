# ============================================================================
# NHS Sickness Rate Analysis - Improved Visualization & Trend Detection
# ============================================================================
# 
# This script:
# 1. Improves visual readability (colors, annotations, layout)
# 2. Detects the exact month of post-COVID trend inflection
# 3. Adds a Nov 2021+ trend line in a distinct color
#
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

# ---- Linear trend fit on rolling average for 2010-2018 (unbiased baseline) ----
df_fit <- df |>
  filter(year(month) >= 2010, year(month) <= 2018, !is.na(roll6_100k))

fit_2010_2018 <- lm(roll6_100k ~ month_num, data = df_fit)

# ---- Predictions + confidence intervals (2010-2018 trend) ----
pred_all <- as_tibble(predict(fit_2010_2018, newdata = df, interval = "confidence", level = 0.95)) |>
  rename(trend_100k = fit, ci_low_100k = lwr, ci_high_100k = upr)

df <- bind_cols(df, pred_all) |>
  mutate(
    period = case_when(
      year(month) >= 2010 & year(month) <= 2018 ~ "fit_2010_2018",
      year(month) >= 2019 & year(month) <= 2025 ~ "proj_2019_2025",
      TRUE ~ "other"
    ),
    dev_pct = if_else(
      year(month) >= 2019 & year(month) <= 2025 & !is.na(roll6_100k) & !is.na(trend_100k),
      100 * (roll6_100k - trend_100k) / trend_100k,
      NA_real_
    ),
    dev_label = if_else(!is.na(dev_pct), sprintf("%+.1f%%", dev_pct), NA_character_)
  )


# ============================================================================
# REQUIREMENT 2: Detect ALL Structural Breaks (Residual-Based Heuristic)
# ============================================================================
# 
# Strategy: Find ALL months where the 6-month rolling average changes regime
# relative to the 2010-2018 trend projection. We detect transitions between
# sustained above-trend and sustained below-trend periods.
#
# This is an UNBIASED approach:
# - Baseline: 2010-2018 (excludes any potential late-2019 contamination)
# - Data: 6-month rolling average (smoothed, not raw)
# - No limit on number of breaks detected
#
# ============================================================================

cat("\n========================================\n")
cat("STRUCTURAL BREAK DETECTION (Heuristic)\n")
cat("========================================\n")
cat("Baseline: 2010-2018 trend projection\n")
cat("Data: 6-month rolling average\n")
cat("Detecting ALL regime changes\n\n")

# ---- Calculate residuals from 2010-2018 trend ----
df <- df |>
  mutate(
    residual = roll6_100k - trend_100k,
    above_trend = residual > 0
  )

# ---- Find ALL structural breaks ----
# A break occurs when we transition between sustained above-trend and below-trend

df_breaks <- df |>
  filter(!is.na(roll6_100k), !is.na(above_trend)) |>
  arrange(month) |>
  mutate(
    # Count consecutive months above/below trend (looking backward)
    consec_above = slide_dbl(above_trend, sum, .before = 5, .after = 0, .complete = TRUE),
    consec_below = slide_dbl(!above_trend, sum, .before = 5, .after = 0, .complete = TRUE),
    # Sustained states
    sustained_above = consec_above >= 6,
    sustained_below = consec_below >= 6,
    # Previous states
    prev_sustained_above = lag(sustained_above, default = FALSE),
    prev_sustained_below = lag(sustained_below, default = FALSE),
    # Detect transitions (break = first month of new sustained regime)
    break_up = sustained_above & !prev_sustained_above,
    break_down = sustained_below & !prev_sustained_below
  )

# Extract all breaks
breaks_up <- df_breaks |>
  filter(break_up) |>
  mutate(direction = "up") |>
  select(month, direction, residual)

breaks_down <- df_breaks |>
  filter(break_down) |>
  mutate(direction = "down") |>
  select(month, direction, residual)

all_breaks <- bind_rows(breaks_up, breaks_down) |>
  arrange(month) |>
  mutate(break_num = row_number())

# ---- Report findings ----
cat(sprintf("Found %d structural breaks:\n", nrow(all_breaks)))
cat("----------------------------------------\n")

if (nrow(all_breaks) > 0) {
  for (i in 1:nrow(all_breaks)) {
    direction_symbol <- if_else(all_breaks$direction[i] == "up", "↑ UP", "↓ DOWN")
    cat(sprintf("  Break %d: %s  %s\n", 
                i,
                format(all_breaks$month[i], "%B %Y"),
                direction_symbol))
  }
} else {
  cat("  No sustained breaks detected\n")
}
cat("========================================\n\n")

# For backward compatibility - use first two breaks if they exist
inflection_month_1 <- if (nrow(all_breaks) >= 1) all_breaks$month[1] else NA
inflection_month_2 <- if (nrow(all_breaks) >= 2) all_breaks$month[2] else NA
inflection_ci_lower_1 <- NA
inflection_ci_upper_1 <- NA
inflection_ci_lower_2 <- NA
inflection_ci_upper_2 <- NA
inflection_month <- inflection_month_1
inflection_ci_lower <- NA
inflection_ci_upper <- NA

# ============================================================================
# FIND ALL LOCAL HIGHS AND LOWS (for deviation labels)
# ============================================================================

cat("Finding all local highs and lows...\n")

# Find local extrema in the rolling average
df <- df |>
  mutate(
    # Local max: higher than both neighbors
    is_local_max = roll6_100k > lag(roll6_100k) & roll6_100k > lead(roll6_100k),
    # Local min: lower than both neighbors
    is_local_min = roll6_100k < lag(roll6_100k) & roll6_100k < lead(roll6_100k)
  )

# Get all local maxima and minima from 2010 onwards
local_extrema <- df |>
  filter(year(month) >= 2010, !is.na(roll6_100k)) |>
  filter(is_local_max | is_local_min) |>
  mutate(
    extrema_type = if_else(is_local_max, "high", "low"),
    dev_pct_extrema = 100 * (roll6_100k - trend_100k) / trend_100k,
    dev_label_extrema = sprintf("%+.1f%%", dev_pct_extrema)
  ) |>
  select(month, roll6_100k, trend_100k, extrema_type, dev_pct_extrema, dev_label_extrema)

cat(sprintf("Found %d local extrema (highs and lows) from 2010 onwards\n", nrow(local_extrema)))
cat(sprintf("  - %d highs\n", sum(local_extrema$extrema_type == "high")))
cat(sprintf("  - %d lows\n", sum(local_extrema$extrema_type == "low")))

# ============================================================================
# TREND LINES FOR EACH SEGMENT
# ============================================================================

# --- Trend 1: November 2021+ (post-vaccine plateau) ---
df_nov2021_plus <- df |>
  filter(month >= as.Date("2021-11-01"), !is.na(roll6_100k))

fit_nov2021 <- lm(roll6_100k ~ month_num, data = df_nov2021_plus)

df <- df |>
  mutate(
    trend_nov2021 = if_else(
      month >= as.Date("2021-11-01"),
      predict(fit_nov2021, newdata = pick(everything())),
      NA_real_
    )
  )

# --- Trend 2: COVID surge period (Break 1 to Break 2) ---
if (!is.na(inflection_month_1) && !is.na(inflection_month_2)) {
  df_covid_surge <- df |>
    filter(month >= inflection_month_1, month < inflection_month_2, !is.na(roll6_100k))
  
  if (nrow(df_covid_surge) >= 3) {
    fit_covid_surge <- lm(roll6_100k ~ month_num, data = df_covid_surge)
    
    df <- df |>
      mutate(
        trend_covid_surge = if_else(
          month >= inflection_month_1 & month < inflection_month_2,
          predict(fit_covid_surge, newdata = pick(everything())),
          NA_real_
        )
      )
    slope_covid_surge <- coef(fit_covid_surge)[2] * 365.25
  } else {
    df <- df |> mutate(trend_covid_surge = NA_real_)
    slope_covid_surge <- NA_real_
  }
} else {
  df <- df |> mutate(trend_covid_surge = NA_real_)
  slope_covid_surge <- NA_real_
}

# --- Trend 3: Post-Break 2 to end ---
if (!is.na(inflection_month_2)) {
  df_post_break2 <- df |>
    filter(month >= inflection_month_2, !is.na(roll6_100k))
  
  if (nrow(df_post_break2) >= 3) {
    fit_post_break2 <- lm(roll6_100k ~ month_num, data = df_post_break2)
    
    df <- df |>
      mutate(
        trend_post_break2 = if_else(
          month >= inflection_month_2,
          predict(fit_post_break2, newdata = pick(everything())),
          NA_real_
        )
      )
    slope_post_break2 <- coef(fit_post_break2)[2] * 365.25
  } else {
    df <- df |> mutate(trend_post_break2 = NA_real_)
    slope_post_break2 <- NA_real_
  }
} else {
  df <- df |> mutate(trend_post_break2 = NA_real_)
  slope_post_break2 <- NA_real_
}

# --- Trend 4: Full post-Break1 (Break 1 to end) - for comparison ---
if (!is.na(inflection_month_1)) {
  df_inflection_plus <- df |>
    filter(month >= inflection_month_1, !is.na(roll6_100k))
  
  fit_inflection <- lm(roll6_100k ~ month_num, data = df_inflection_plus)
  
  df <- df |>
    mutate(
      trend_inflection = if_else(
        month >= inflection_month_1,
        predict(fit_inflection, newdata = pick(everything())),
        NA_real_
      )
    )
  slope_inflection <- coef(fit_inflection)[2] * 365.25
} else {
  df <- df |> mutate(trend_inflection = NA_real_)
  slope_inflection <- NA_real_
}

# Report all trends
slope_nov2021 <- coef(fit_nov2021)[2] * 365.25  # Annualized
slope_2010_2018 <- coef(fit_2010_2018)[2] * 365.25

cat("\n========================================\n")
cat("TREND COMPARISON (Annualized Slopes)\n")
cat("========================================\n")
cat(sprintf("Pre-COVID (2010-2018):     %+7.1f per 100k per year\n", slope_2010_2018))
if (!is.na(slope_covid_surge)) {
  cat(sprintf("COVID surge (B1-B2):       %+7.1f per 100k per year\n", slope_covid_surge))
}
if (!is.na(slope_post_break2)) {
  cat(sprintf("Post-Break2:               %+7.1f per 100k per year\n", slope_post_break2))
}
cat(sprintf("Nov 2021+ trend:           %+7.1f per 100k per year\n", slope_nov2021))
if (!is.na(slope_inflection)) {
  cat(sprintf("Full post-Break1:          %+7.1f per 100k per year\n", slope_inflection))
}
cat("========================================\n\n")

# ============================================================================
# VISUALIZATION
# ============================================================================

# ---- Marker dates (known events) ----
markers <- tibble(
  date = as.Date(c("2020-01-31", "2020-12-08", "2021-11-01", "2022-04-15")),
  short_label = c("COVID", "Vax start", "90% vax", "Full vax")
)

# Calculate y-positions for markers (staggered)
y_range <- range(c(df$rate_100k, df$roll6_100k, df$trend_100k), na.rm = TRUE)
y_top <- y_range[2]
y_bottom <- y_range[1]

markers <- markers |>
  mutate(
    y_label = if_else(row_number() %% 2 == 1, y_top * 0.98, y_top * 0.92)
  )

# ---- Color palette ----
color_raw <- "grey50"
color_rolling <- "#E69F00"  # Orange
color_trend_2010 <- "#0072B2"  # Blue
color_trend_nov21 <- "#009E73"  # Teal/Green
color_markers <- "#D55E00"  # Vermillion
color_breaks <- "#CC79A7"  # Pink for structural breaks
color_covid_surge <- "#E31A1C"  # Red
color_trend_inflection <- "#7B3294"  # Purple
color_deviation <- "#666666"  # Grey for deviation labels

# ---- Create the plot ----
p <- ggplot(df, aes(x = month)) +
  
  # Confidence band for 2010-2018 projection
  geom_ribbon(
    data = df |> filter(period == "proj_2019_2025"),
    aes(ymin = ci_low_100k, ymax = ci_high_100k),
    fill = color_trend_2010, alpha = 0.15
  ) +
  
  # Raw monthly data (subdued)
  geom_line(aes(y = rate_100k), color = color_raw, linewidth = 0.8, alpha = 0.5) +
  
  # 6-month rolling average (prominent)
  geom_line(aes(y = roll6_100k), color = color_rolling, linewidth = 1.5, na.rm = TRUE) +
  
  # 2010-2018 trend (solid for fit period)
  geom_line(
    data = df |> filter(period == "fit_2010_2018"),
    aes(y = trend_100k),
    color = color_trend_2010, linewidth = 1.4
  ) +
  
  # 2010-2018 projection (dashed for projection period)
  geom_line(
    data = df |> filter(period == "proj_2019_2025"),
    aes(y = trend_100k),
    color = color_trend_2010, linewidth = 1.4, linetype = "dashed"
  ) +
  
  # November 2021+ trend
  geom_line(
    data = df |> filter(!is.na(trend_nov2021)),
    aes(y = trend_nov2021),
    color = color_trend_nov21, linewidth = 1.4, linetype = "solid"
  ) +
  
  # COVID surge trend (if exists)
  geom_line(
    data = df |> filter(!is.na(trend_covid_surge)),
    aes(y = trend_covid_surge),
    color = color_covid_surge, linewidth = 1.4, linetype = "solid"
  ) +
  
  # Full post-Break1 trend (purple dashed)
  geom_line(
    data = df |> filter(!is.na(trend_inflection)),
    aes(y = trend_inflection),
    color = color_trend_inflection, linewidth = 1.2, linetype = "dashed"
  ) +
  
  # ALL deviation labels at local highs and lows
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
  
  # Event markers (known dates)
  geom_vline(
    data = markers,
    aes(xintercept = date),
    linetype = "dotted", linewidth = 0.8,
    color = color_markers, alpha = 0.7
  ) +
  
  # Event labels
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
  
  # ALL structural breaks (vertical dashed lines)
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
                            y_bottom + (y_top - y_bottom) * 0.12,
                            y_bottom + (y_top - y_bottom) * 0.05),
            label_text = sprintf("B%d %s", break_num, if_else(direction == "up", "↑", "↓"))
          ),
          aes(x = month, y = y_pos, label = label_text),
          color = color_breaks,
          fill = alpha("white", 0.9),
          size = 2.5,
          fontface = "bold",
          label.padding = unit(0.1, "lines")
        )
      )
    }
  } +
  
  # Scales
  scale_x_date(
    date_breaks = "2 years",
    date_labels = "%Y",
    expand = expansion(mult = c(0.02, 0.05))
  ) +
  scale_y_continuous(
    limits = c(2000, NA),
    labels = label_comma(),
    expand = expansion(mult = c(0.02, 0.12))
  ) +
  
  # Labels
  labs(
    title = "NHS Sickness Absence Rate (England) - Structural Break Analysis",
    subtitle = paste0(
      "Rolling avg (orange) | 2010-18 baseline (blue) | COVID surge (red) | ",
      "Nov 2021+ (green) | Breaks (pink dashed)"
    ),
    x = NULL,
    y = "Sickness rate (per 100,000 staff)",
    caption = paste0(
      sprintf("%d structural breaks detected | ", nrow(all_breaks)),
      sprintf("Baseline slope: %+.0f/yr | Nov 2021+: %+.0f/yr", 
              slope_2010_2018, slope_nov2021)
    )
  ) +
  
  # Theme
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(size = 16, face = "bold", margin = margin(b = 5)),
    plot.subtitle = element_text(size = 10, color = "grey40", margin = margin(b = 10)),
    plot.caption = element_text(size = 8, color = "grey50", hjust = 0, margin = margin(t = 10)),
    axis.text = element_text(size = 10),
    axis.title.y = element_text(size = 11, margin = margin(r = 10)),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_line(color = "grey90", linewidth = 0.3),
    panel.grid.major.y = element_line(color = "grey85", linewidth = 0.4),
    plot.margin = margin(t = 10, r = 15, b = 10, l = 10),
    legend.position = "none"
  )

# ---- Display plot ----
print(p)

# ---- Save high-resolution output ----
ggsave(
  "nhs_sickness_analysis_improved.png",
  plot = p,
  width = 16,
  height = 9,
  dpi = 300,
  bg = "white"
)

cat("\nPlot saved as 'nhs_sickness_analysis_improved.png'\n")

# ============================================================================
# SUMMARY OUTPUT
# ============================================================================

cat("\n")
cat("========================================================================\n")
cat("               ANALYSIS SUMMARY                                        \n")
cat("========================================================================\n")
cat(sprintf("Structural breaks detected: %d\n", nrow(all_breaks)))
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
cat(sprintf("  Pre-COVID (2010-2018):    %+7.1f\n", slope_2010_2018))
if (!is.na(slope_covid_surge)) {
  cat(sprintf("  COVID surge (B1-B2):      %+7.1f  (%.1fx baseline)\n", 
              slope_covid_surge, slope_covid_surge/slope_2010_2018))
}
if (!is.na(slope_post_break2)) {
  cat(sprintf("  Post-Break2:              %+7.1f  (%.1fx baseline)\n", 
              slope_post_break2, slope_post_break2/slope_2010_2018))
}
cat(sprintf("  Nov 2021+:                %+7.1f  (%.1fx baseline)\n", 
            slope_nov2021, slope_nov2021/slope_2010_2018))
if (!is.na(slope_inflection)) {
  cat(sprintf("  Full post-Break1:         %+7.1f  (%.1fx baseline)\n", 
              slope_inflection, slope_inflection/slope_2010_2018))
}
cat("========================================================================\n")
# ============================================================================
# positivity_vs_sickness.R
# Plot: NHS sickness (roll6 + prediction interval) vs COVID positivity (roll6)
# + Adds: tests performed (monthly mean of 7-day totals) as COLUMNS with vertical labels
#
# Uses combined file:
#   data/tests_positivity_and_volume.csv
#     month, positivity_mean_7day, tests_mean_7day_total
#
# Restores from first plot:
#   1) Pink dashed breakpoint labels (B1/B2/B3)
#   2) Orange dashed event markers for "COVID" and "Vax start"
#
# Notes:
#   - Sickness roll6 + prediction interval are read from sickness CSV (already computed).
#   - Positivity roll6 (-3..+2) is computed from monthly mean positivity series.
#   - Tests performed are shown as columns in a dedicated "band" at the bottom of the plot
#     (no extra y-axis). Labels are vertical inside each bar.
#   - Plot starts at 2020-05 (so positivity roll6 is defined).
#   - Left axis not forced to 0; right axis starts at 0.
# ============================================================================

library(tidyverse)
library(lubridate)
library(slider)
library(scales)

# ---- Force English month names for axis labels (Windows/French OS safe) ----
old_time_locale <- Sys.getlocale("LC_TIME")
on.exit(try(Sys.setlocale("LC_TIME", old_time_locale), silent = TRUE), add = TRUE)

# Try a few common English locales; fall back to "C" which is usually English.
ok <- suppressWarnings(Sys.setlocale("LC_TIME", "C"))
if (is.na(ok) || ok == "") {
  for (loc in c("English_United Kingdom.1252", "English_United States.1252", "en_GB.UTF-8", "en_US.UTF-8")) {
    ok <- suppressWarnings(Sys.setlocale("LC_TIME", loc))
    if (!is.na(ok) && ok != "") break
  }
}

COMBINED_PATH <- "data/tests_positivity_and_volume.csv"
SICK_PATH     <- "data/sickness_roll6_prediction_band_2020_2026.csv"

PLOT_START <- as.Date("2019-01-01")   # so the "not enough test data" grey period is visible
PLOT_END   <- as.Date("2025-08-01")   # restrict visual to Aug 2025 (included)

# Rolling window for positivity: (t-3, t-2, t-1, t, t+1, t+2)
ROLL_BEFORE <- 3
ROLL_AFTER  <- 2
ROLL_COMPLETE <- TRUE

# ---- Load combined (monthly mean series) and compute positivity roll6 ----
cov <- read.csv(COMBINED_PATH, stringsAsFactors = FALSE) |>
  as_tibble() |>
  transmute(
    month = as.Date(paste0(month, "-01")),
    pos_monthly = as.numeric(positivity_mean_7day),
    tests_7day_mean = as.numeric(tests_mean_7day_total)
  ) |>
  arrange(month) |>
  mutate(
    pos_roll6 = slide_dbl(
      pos_monthly,
      ~ mean(.x, na.rm = TRUE),
      .before = ROLL_BEFORE,
      .after  = ROLL_AFTER,
      .complete = ROLL_COMPLETE
    )
  )

# ---- Load sickness (already roll6 + prediction interval) ----
sick <- read.csv(SICK_PATH, stringsAsFactors = FALSE) |>
  as_tibble() |>
  transmute(
    month = as.Date(paste0(month, "-01")),
    sickness_roll6_100k = as.numeric(sickness_roll6_100k),
    sickness_pred_low_100k = as.numeric(sickness_pred_low_100k),
    sickness_pred_high_100k = as.numeric(sickness_pred_high_100k)
  ) |>
  arrange(month)

# ---- Join + filter window ----
df <- sick |>
  left_join(cov |> select(month, pos_roll6, tests_7day_mean), by = "month") |>
  filter(month >= PLOT_START, month <= PLOT_END) |>
  arrange(month)

# ---- Scale positivity onto sickness axis (shifted) so:
#      - left axis doesn't have to start at 0
#      - right axis can start at 0
left_min_raw <- min(df$sickness_pred_low_100k, df$sickness_roll6_100k, na.rm = TRUE)
left_max_raw <- max(df$sickness_pred_high_100k, df$sickness_roll6_100k, na.rm = TRUE)

pad <- 0.04 * (left_max_raw - left_min_raw)
left_min <- left_min_raw - pad
left_max <- left_max_raw + pad

max_pos <- max(df$pos_roll6, na.rm = TRUE)
if (!is.finite(max_pos) || max_pos <= 0) {
  stop("No valid positivity rolling values after ", START_DATE,
       ". Check tests_positivity_and_volume.csv positivity_mean_7day.")
}

k <- (left_max - left_min) / max_pos

df <- df |>
  mutate(pos_scaled = left_min + pos_roll6 * k)

# ---- Tests bars: allocate a bottom band (no new y-axis) ----
# We'll draw bars in [left_min, left_min + band_height], scaled by tests value.
band_height <- 0.70 * (left_max - left_min) 
band_y0 <- left_min
band_y1 <- left_min + band_height

tests_max <- max(df$tests_7day_mean, na.rm = TRUE)
if (!is.finite(tests_max) || tests_max <= 0) {
  stop("No valid tests values after ", START_DATE,
       ". Check tests_positivity_and_volume.csv tests_mean_7day_total.")
}

df <- df |>
  mutate(
    tests_scaled_height = (tests_7day_mean / tests_max) * band_height,
    tests_bar_top = band_y0 + tests_scaled_height,
    tests_label = if_else(!is.na(tests_7day_mean), comma(round(tests_7day_mean)), NA_character_),
    tests_label_y = band_y0 + pmax(tests_scaled_height * 0.55, band_height * 0.12) # keep label inside bar
  )




# ============================================================================
# CORRELATION: sickness vs positivity (two periods)
# A: 2020-05 .. 2022-04 (inclusive)
# B: 2022-05 .. PLOT_END (inclusive)
# ============================================================================

cor_block <- function(dat, label) {
  d <- dat |>
    filter(!is.na(pos_roll6), !is.na(sickness_roll6_100k))

  n <- nrow(d)
  cat("\n----------------------------------------\n")
  cat("Correlation period:", label, "\n")
  cat("N (complete pairs):", n, "\n")

  if (n < 3) {
    cat("Not enough data points to compute correlation.\n")
    return(invisible(NULL))
  }

  pear <- cor.test(d$pos_roll6, d$sickness_roll6_100k, method = "pearson")
  spear <- cor.test(d$pos_roll6, d$sickness_roll6_100k, method = "spearman", exact = FALSE)

  cat("\nPearson:\n")
  cat(sprintf("  r = %.3f\n", unname(pear$estimate)))
  cat(sprintf("  p-value = %.4g\n", pear$p.value))
  cat(sprintf("  95%% CI = [%.3f, %.3f]\n", pear$conf.int[1], pear$conf.int[2]))

  cat("\nSpearman:\n")
  cat(sprintf("  rho = %.3f\n", unname(spear$estimate)))
  cat(sprintf("  p-value = %.4g\n", spear$p.value))

  invisible(list(pearson = pear, spearman = spear, n = n))
}

# --- Define periods (inclusive bounds) ---
p1_start <- as.Date("2020-05-01")
p1_end   <- as.Date("2022-04-01")
p2_start <- as.Date("2022-05-01")
p2_end   <- PLOT_END

df_p1 <- df |> filter(month >= p1_start, month <= p1_end)
df_p2 <- df |> filter(month >= p2_start, month <= p2_end)

res_p1 <- cor_block(df_p1, "2020-05 to 2022-04")
res_p2 <- cor_block(df_p2, paste0("2022-05 to ", format(p2_end, "%Y-%m")))

# --- Correlation label boxes (robust date handling) ---
corr_labels <- tibble(
  period = c("corr_high", "corr_none"),
  xmin   = ymd(c("2020-05-01", "2022-05-01")),
  xmax   = ymd(c("2022-05-01", format(p2_end %m+% months(1), "%Y-%m-%d")))
) |>
  mutate(
    # robust midpoint (avoids Date + difftime edge cases)
    x = as.Date((as.numeric(xmin) + as.numeric(xmax)) / 2, origin = "1970-01-01"),
    # place the red label a bit lower so it doesn't get lost among event labels
    y = c(
      min(5250, left_max - 0.18 * (left_max - left_min)),
      left_max - 0.16 * (left_max - left_min)
    ),
    label = c(
      sprintf(
        "Correlation (May 2020–Apr 2022)\nPearson r=%.3f (p=%s)\nSpearman ρ=%.3f (p=%s)\nN=%d",
        unname(res_p1$pearson$estimate), fmt_p(res_p1$pearson$p.value),
        unname(res_p1$spearman$estimate), fmt_p(res_p1$spearman$p.value),
        res_p1$n
      ),
      sprintf(
        "Correlation (May 2022–%s)\nPearson r=%.3f (p=%s)\nSpearman ρ=%.3f (p=%s)\nN=%d",
        format(p2_end, "%b %Y"),
        unname(res_p2$pearson$estimate), fmt_p(res_p2$pearson$p.value),
        unname(res_p2$spearman$estimate), fmt_p(res_p2$spearman$p.value),
        res_p2$n
      )
    )
  )

# Sanity checks: if the 2nd row is NA, it explains the missing label immediately
stopifnot(!anyNA(corr_labels$x), !anyNA(corr_labels$y), !anyNA(corr_labels$label))
print(corr_labels)

# ============================================================================
# VISUALIZATION
# ============================================================================

# ---- Breakpoints (hard-coded) + labels (pink dashed + label) ----
all_breaks <- tibble(
  break_num = 1:3,
  month = as.Date(c("2019-09-01", "2021-05-01", "2022-12-01"))
) |>
  mutate(label_text = sprintf("B%d\n%s", break_num, format(month, "%b %Y")))

# ---- Event markers ----
markers <- tibble(
  date = as.Date(c("2020-01-31", "2020-12-08", "2021-11-01", "2022-04-15")),
  short_label = c("COVID", "Vax start", "90% vax & Mandates", "Drop Mandates")
)

y_top <- left_max
markers <- markers |>
  mutate(y_label = if_else(row_number() %% 2 == 1, y_top * 0.985, y_top * 0.94))

# ---- X-axis: every 3 months + thicker Jan separators ----
x_min <- floor_date(min(df$month, na.rm = TRUE), unit = "month")
x_max <- ceiling_date(max(df$month, na.rm = TRUE), unit = "month")

x_breaks_3m <- seq(x_min, x_max, by = "3 months")
x_breaks_jan <- seq(floor_date(x_min, "year"), floor_date(x_max, "year"), by = "1 year")

grid_3m <- tibble(x = x_breaks_3m)
grid_jan <- tibble(x = x_breaks_jan)

label_3m <- function(d) {
  format(as.Date(d), "%b\n%Y")   # e.g., "May\n2021"
}

# ============================================================================
# Background shading for periods + correlation labels
# ============================================================================
periods_bg <- tibble(
  period = c("no_tests", "corr_high", "corr_none"),
  xmin   = ymd(c(format(PLOT_START, "%Y-%m-%d"), "2020-05-01", "2022-05-01")),
  xmax   = ymd(c("2020-05-01", "2022-05-01", format(PLOT_END %m+% months(1), "%Y-%m-%d"))),
  fill   = c("#EEEEEE", "#DFF2E1", "#F6DADA")  # light grey, light green, light red
)

# Sanity checks: if these fail, ggplot will silently drop rows
stopifnot(all(!is.na(periods_bg$xmin)), all(!is.na(periods_bg$xmax)))
stopifnot(all(periods_bg$xmax > periods_bg$xmin))

# Helper to format p-values nicely
fmt_p <- function(p) {
  if (is.na(p)) return("NA")
  if (p < 1e-4) return(sprintf("%.1e", p))
  sprintf("%.4f", p)
}

# ---- Colors (same palette) ----
color_rolling  <- "#E69F00"  # orange
color_baseline <- "#0072B2"  # blue
color_markers  <- "#D55E00"  # orange/red (markers)
color_breaks   <- "#CC79A7"  # pink
color_pos      <- "#D55E00"  # red line for positivity
color_tests    <- "grey70"   # neutral bars

p <- ggplot(df, aes(x = month)) +

  # --- Background shading for periods (behind everything) ---
  geom_rect(
    data = periods_bg,
    aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = fill),
    inherit.aes = FALSE,
    alpha = 0.70
  ) +
  scale_fill_identity() +

  # --- Month/year separators FIRST (behind everything) ---
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

  # --- Tests performed band (bars + vertical labels) ---
  geom_rect(
    aes(xmin = month - days(12), xmax = month + days(12), ymin = band_y0, ymax = tests_bar_top),
    fill = color_tests, alpha = 0.75
  ) +
  geom_text(
    aes(x = month, y = tests_label_y, label = tests_label),
    angle = 90,
    size = 2.5,
    color = "grey20",
    na.rm = TRUE
  ) +

  # Prediction interval (blue ribbon)
  geom_ribbon(
    aes(ymin = sickness_pred_low_100k, ymax = sickness_pred_high_100k),
    fill = color_baseline, alpha = 0.12
  ) +

  # Sickness rolling (orange)
  geom_line(
    aes(y = sickness_roll6_100k),
    color = color_rolling, linewidth = 1.5, na.rm = TRUE
  ) +

  # Positivity rolling (red), scaled
  geom_line(
    aes(y = pos_scaled),
    color = color_pos, linewidth = 1.2, na.rm = TRUE
  ) +

  # Event markers (orange dashed vertical lines)
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

  # Structural breaks (pink dashed) + labels
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

  # --- Correlation labels in background (top area) ---
  geom_label(
    data = corr_labels,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    fill = alpha("white", 0.88),
    color = "grey20",
    size = 2.8,
    label.size = 0.2,
    label.padding = unit(0.15, "lines"),
    lineheight = 0.95
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
      ~ (. - left_min) / k,
      name = "COVID-19 test positivity (%)",
      breaks = sort(unique(c(0, pretty(df$pos_roll6, n = 6))))
    )
  ) +

  labs(
    title = "NHS sickness vs COVID positivity (England)",
    subtitle = paste(
      "Blue band: predicted interval | Orange: sickness roll6 (-3..+2) | Red: positivity roll6 (-3..+2)",
      "| Grey bars: PCR tests (monthly mean of 7-day totals) | Start: 2020-05",
      sep = " "
    ),
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
  "nhs_sickness_vs_positivity.png",
  plot = p,
  width = 16,
  height = 9,
  dpi = 300,
  bg = "white"
)

cat("\nPlot saved as 'nhs_sickness_vs_positivity.png'\n")
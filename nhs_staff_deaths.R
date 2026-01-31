# ---- Packages ----
library(readr)
library(dplyr)
library(stringr)
library(ggplot2)

# ---- 1) Load both files ----
deaths <- read_csv("data/deaths_in_service_long.csv",
                   col_types = cols(
                     year_from_to = col_character(),
                     quarter      = col_character(),
                     deaths       = col_double()
                   ))

denoms <- read_csv("data/april_denominators_2011_2025.csv",
                   col_types = cols(
                     period      = col_character(),
                     denominator = col_double()
                   ))

# Helper: convert "YYYY-YY" (e.g. "2011-12", "2019-20") -> end year as integer (2012, 2020)
year_from_to_to_end_year <- function(x) {
  start_year <- as.integer(str_sub(x, 1, 4))
  end_two    <- as.integer(str_sub(x, 6, 7))
  # assume 2000s for typical data (11..26 etc); adjust if ever needed
  end_year   <- if_else(end_two < 50L, 2000L + end_two, 1900L + end_two)
  
  # sanity: if something like 1999-00 appears, the above still works (end_year=2000)
  end_year
}

# ---- 2) Concatenate yearly deaths (sum Q1..Q4 per year_from_to) ----
yearly_deaths <- deaths %>%
  group_by(year_from_to) %>%
  summarise(deaths_year = sum(deaths, na.rm = TRUE), .groups = "drop") %>%
  mutate(
    year_end = year_from_to_to_end_year(year_from_to),
    # ---- 3) Build the denominator period key per your rule:
    # for 2011 -> "201004 to 201104", for 2012 -> "201104 to 201204", etc.
    period = sprintf("%d04 to %d04", year_end - 1L, year_end)
  )

# ---- 3) Join denominators ----
dat <- yearly_deaths %>%
  left_join(denoms, by = "period")

# ---- 4) Calculate yearly normalized deaths per 100,000 ----
dat <- dat %>%
  mutate(
    rate_per_100k = (deaths_year / denominator) * 100000
  )

# ---- Fit linear trend on 2011–2019 ----
trend_fit <- lm(rate_per_100k ~ year_end, data = dat %>% filter(year_end >= 2011, year_end <= 2019))

# Prediction for all available years
pred_df <- dat %>%
  distinct(year_end, year_from_to) %>%
  arrange(year_end) %>%
  mutate(pred_rate = predict(trend_fit, newdata = .))

# Join predictions and compute excess %
dat2 <- dat %>%
  left_join(pred_df %>% select(year_end, pred_rate), by = "year_end") %>%
  mutate(excess_pct = (rate_per_100k - pred_rate) / pred_rate * 100)

# Labels for excess % (years 2021–2025 only)
excess_labels <- dat2 %>%
  filter(year_end >= 2021, year_end <= 2025) %>%
  mutate(label = sprintf("%+.1f%%", excess_pct))

# Ensure year interval is ordered chronologically on the x-axis
dat2 <- dat2 %>%
  arrange(year_end) %>%
  mutate(year_from_to = factor(year_from_to, levels = unique(year_from_to)))

pred_df <- pred_df %>%
  arrange(year_end) %>%
  mutate(year_from_to = factor(year_from_to, levels = levels(dat2$year_from_to)))

# Split trend into solid (<=2019) and dashed (>=2020)
pred_solid <- pred_df %>% filter(year_end <= 2019)
pred_dashed <- pred_df %>% filter(year_end >= 2020)

# ---- Line chart: actual + trend solid then dashed ----
p <- ggplot(dat2, aes(x = year_from_to)) +
  # Actual series
  geom_line(aes(y = rate_per_100k, group = 1), linewidth = 1.4) +
  geom_point(aes(y = rate_per_100k), size = 2.8) +
  
  # Trend (solid 2011–2019)
  geom_line(
    data = pred_solid,
    aes(y = pred_rate, group = 1),
    color = "red",
    linewidth = 1.4
  ) +
  
  # Trend (dashed 2020+)
  geom_line(
    data = pred_dashed,
    aes(y = pred_rate, group = 1),
    color = "red",
    linetype = "dashed",
    linewidth = 1.4
  ) +
  
  # Excess % labels for 2021–2025
  geom_text(
    data = excess_labels,
    aes(y = rate_per_100k, label = label),
    vjust = -0.9,
    size = 4
  ) +
  
  labs(
    x = "Year interval",
    y = "Deaths per 100,000 (normalized by April year start denominator)",
    title = "NHS Staff, Yearly Deaths per 100,000",
    subtitle = "Red line: solid on linear trend fit 2011–2019, dashed from 2020; labels show excess vs trend for 2021–2025"
  ) +
  scale_y_continuous(limits = c(1, 150), breaks = seq(0, 150, by = 10)) +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 12),
    axis.title  = element_text(size = 14),
    plot.title  = element_text(size = 16, face = "bold"),
    plot.subtitle = element_text(size = 13),
    panel.grid.minor = element_blank()
  )

print(p)

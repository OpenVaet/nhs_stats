# ---- Packages ----
# install.packages(c("tidyverse", "lubridate", "slider", "scales"))
library(tidyverse)
library(lubridate)
library(slider)
library(scales)

# ---- Read data ----
df <- read.csv("data/overall_sickness_rates.csv", stringsAsFactors = FALSE) |>
  as_tibble() |>
  rename(rate_pct = England) |>
  mutate(
    month = my(Month),
    rate_pct = as.numeric(rate_pct)
  ) |>
  arrange(month) |>
  mutate(
    # 1) convert percent -> per 100,000 (i.e., multiply by 1000)
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

fit <- lm(roll6_100k ~ month_num, data = df_fit)

# ---- Predictions + confidence intervals ----
pred_all <- as_tibble(predict(fit, newdata = df, interval = "confidence", level = 0.95)) |>
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

# ---- Label only every 10th projected point (and only where roll6 exists) ----
df_labels <- df |>
  filter(period == "proj_2020_2025", !is.na(dev_label), !is.na(roll6_100k)) |>
  mutate(idx = row_number()) |>
  filter((idx - 1) %% 10 == 0)

# ---- Plot ----
# ---- Marker dates (known metadata / coverage / classification changes) ----
markers <- tibble(
  date = as.Date(c("2020-01-31", "2020-12-08", "2021-11-28", "2022-04-15")),
  label = c(
    "Jan 31, 2020\nCOVID",
    "Dec 2020\nFirst staff vaccinated",
    "Nov 2021\nNHS Staff\n90% vaccinated",
    "Apr 2022\nNHS Staff\nfully vaccinated"
  )
)

y_top <- max(df$rate_100k, df$roll6_100k, df$trend_100k, na.rm = TRUE)

# alternate label heights (high / low / high / low ...)
markers <- markers |>
  mutate(
    y = if_else(row_number() %% 2 == 1, y_top, y_top * 0.93)  # adjust 0.93 to taste
  )

# helper y-position for labels near the top of the plot
y_top <- max(df$rate_100k, df$roll6_100k, df$trend_100k, na.rm = TRUE)

# ---- Plot + markers ----
p <- ggplot(df, aes(x = month)) +
  geom_line(aes(y = rate_100k), color = "black", linewidth = 1.1) +
  geom_line(aes(y = roll6_100k), color = "orange", linewidth = 1.4, na.rm = TRUE) +

  geom_line(
    data = df |> filter(period == "fit_2010_2019"),
    aes(y = trend_100k),
    color = "blue",
    linewidth = 1.3
  ) +
  geom_line(
    data = df |> filter(period == "proj_2020_2025"),
    aes(y = trend_100k),
    color = "blue",
    linewidth = 1.3,
    linetype = "dashed"
  ) +

  geom_text(
    data = df_labels,
    aes(y = roll6_100k, label = dev_label),
    color = "red",
    size = 4.2,
    vjust = -0.6
  ) +

  # ---- Vertical markers + labels ----
  geom_vline(
    data = markers,
    aes(xintercept = date),
    linetype = "dotted",
    linewidth = 1.3,
    color = "red",
    inherit.aes = FALSE
  ) +
  geom_text(
    data = markers,
    aes(x = date, y = y, label = label),
    inherit.aes = FALSE,
    hjust = -0.2,
    vjust = 0.4,
    size = 3.5,
    color = "red"
  ) +

  scale_x_date(
    date_breaks = "1 year",
    date_labels = "%Y",
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  scale_y_continuous(
    limits = c(2000, NA),
    labels = label_number(accuracy = 1)
  ) +
  labs(
    title = "NHS Overall sickness rate (England)",
    subtitle = "Monthly rate (black), past 6-month rolling average (orange), trend 2010–2019 and projection 2020–2025 (blue)",
    x = NULL,
    y = "Sickness rate (per 100,000)"
  ) +
  theme_minimal(base_size = 15) +
  theme(
    plot.title = element_text(size = 18, face = "bold"),
    plot.subtitle = element_text(size = 13),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 13),
    panel.grid.minor = element_blank()
  )

print(p)


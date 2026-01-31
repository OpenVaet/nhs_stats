# sickness_rates_by_categories_plots.R
# Reads:   data/sickness_rates_by_categories.csv
# Writes: outputs/*.png  (one plot per Staff group where Main staff group == "Professionally qualified clinical staff")

# ---- Packages ----
# install.packages(c("tidyverse", "lubridate", "slider", "scales", "readr"))
library(tidyverse)
library(lubridate)
library(slider)
library(scales)
library(readr)

# ---- I/O ----
infile  <- file.path("data", "sickness_rates_by_categories.csv")
out_dir <- "outputs"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Robust read: fix multiline/broken header by collapsing header lines ----
read_fixed_header_csv <- function(path) {
  lines <- readLines(path, warn = FALSE, encoding = "UTF-8")
  lines[1] <- sub("^\ufeff", "", lines[1])  # remove BOM if present

  # find first data row (first two fields then a number)
  first_data <- which(grepl("^[^,]+,[^,]+,[-0-9]", lines))[1]
  if (is.na(first_data) || first_data <= 1) stop("Could not detect header/data boundary in: ", path)

  header <- paste(lines[1:(first_data - 1)], collapse = "")
  new_lines <- c(header, lines[first_data:length(lines)])

  tmp <- tempfile(fileext = ".csv")
  writeLines(new_lines, tmp, useBytes = TRUE)

  # Read everything as character to avoid parsing warnings/issues
  readr::read_csv(
    tmp,
    show_col_types = FALSE,
    col_types = cols(.default = col_character()),
    na = c("", "NA")
  )
}

wide <- read_fixed_header_csv(infile)

# ---- Clean column names ----
names(wide) <- names(wide) |>
  str_replace_all("\\s+", " ") |>
  str_trim() |>
  str_replace_all('^"|"$', "")

# ---- Drop columns with missing/blank names (fixes your error) ----
nm <- names(wide)
keep_named <- !(is.na(nm) | nm == "")
wide <- wide[, keep_named, drop = FALSE]

# ---- Drop columns that are entirely empty/NA (often the ...108 column) ----
wide <- wide |>
  select(where(~ !all(is.na(.x) | .x == "")))

# ---- Wide -> long ----
long <- wide |>
  rename(main_staff_group = `Main staff group`, staff_group = `Staff group`) |>
  pivot_longer(
    cols = -c(main_staff_group, staff_group),
    names_to = "Month",
    values_to = "rate_pct"
  ) |>
  mutate(
    Month = str_replace_all(Month, "\\s+", " ") |> str_trim(),
    month = my(Month),
    rate_pct = as.numeric(rate_pct)
  ) |>
  filter(!is.na(month), !is.na(rate_pct)) |>
  arrange(main_staff_group, staff_group, month) |>
  mutate(
    # convert percent (e.g. 3.97) -> per 100,000 (3.97% = 3.97/100 * 100000 = 3970)
    rate_100k = rate_pct * 1000
  )

# ---- Markers (optional, as in your last version) ----
markers_base <- tibble(
  date = as.Date(c("2020-04-01", "2021-10-01", "2022-02-01", "2023-11-01")),
  label = c(
    "Apr 2020\nRegions remapped",
    "Oct 2021\nUKHSA included",
    "Feb 2022\nDoctor grades revised",
    "Nov 2023\nUKHSA removed"
  )
)

# ---- Plot function ----
make_plot <- function(d, staff_group_name) {
  d <- d |>
    arrange(month) |>
    mutate(
      roll6_100k = slide_dbl(rate_100k, mean, .before = 3, .after = 2, .complete = TRUE),
      month_num  = as.numeric(month)
    )

  # Fit linear trend on 2017–2020 (inclusive) using rolling average
  fit_df <- d |>
    filter(year(month) >= 2017, year(month) <= 2020, !is.na(roll6_100k))

  # if that window is too short (some groups), fall back to all available data up to 2020
  if (nrow(fit_df) < 12) {
    fit_df <- d |>
      filter(year(month) <= 2020, !is.na(roll6_100k))
  }

  fit <- lm(roll6_100k ~ month_num, data = fit_df)

  d <- d |>
    mutate(trend_100k = as.numeric(predict(fit, newdata = d))) |>
    mutate(
      period = case_when(
        year(month) >= 2017 & year(month) <= 2020 ~ "fit_2017_2020",
        year(month) >= 2021 ~ "excess_2021_plus",
        TRUE ~ "other"
      ),
      dev_pct = if_else(
        year(month) >= 2021 & !is.na(roll6_100k) & !is.na(trend_100k),
        100 * (roll6_100k - trend_100k) / trend_100k,
        NA_real_
      ),
      dev_label = if_else(!is.na(dev_pct), sprintf("%+.1f%%", dev_pct), NA_character_)
    )

  # label only every 10th value in 2021+
  df_labels <- d |>
    filter(period == "excess_2021_plus", !is.na(dev_label), !is.na(roll6_100k)) |>
    mutate(idx = row_number()) |>
    filter((idx - 1) %% 10 == 0)

  ggplot(d, aes(x = month)) +
    geom_line(aes(y = rate_100k), color = "black", linewidth = 1.1) +
    geom_line(aes(y = roll6_100k), color = "orange", linewidth = 1.4, na.rm = TRUE) +

    # trend line: solid over 2017–2020, dashed after 2020
    geom_line(
      data = d |> filter(period == "fit_2017_2020"),
      aes(y = trend_100k),
      color = "blue",
      linewidth = 1.3
    ) +
    geom_line(
      data = d |> filter(year(month) >= 2021),
      aes(y = trend_100k),
      color = "blue",
      linewidth = 1.3,
      linetype = "dashed"
    ) +

    # excess labels (2021+ only, every 10 points)
    geom_text(
      data = df_labels,
      aes(y = roll6_100k, label = dev_label),
      color = "red",
      size = 4.2,
      vjust = -0.6
    ) +

    scale_x_date(
      date_breaks = "1 year",
      date_labels = "%Y",
      expand = expansion(mult = c(0.01, 0.02))
    ) +
    scale_y_continuous(
      limits = c(0, NA),
      labels = label_number(accuracy = 1)
    ) +
    labs(
      title = paste0("Overall sickness rate (England) — ", staff_group_name),
      subtitle = "Monthly rate (black), 6-month rolling average (orange), trend 2017–2020 (blue), excess vs trend from 2021+ (labels)",
      x = NULL,
      y = "Sickness rate (per 100,000)"
    ) +

    # white background
    theme_bw(base_size = 15) +
    theme(
      plot.title = element_text(size = 18, face = "bold"),
      plot.subtitle = element_text(size = 13),
      axis.text = element_text(size = 12),
      axis.title = element_text(size = 13),
      panel.grid.minor = element_blank()
    )
}

# ---- Filter requested main staff group and loop staff groups ----
target <- long |> filter(main_staff_group == "Professionally qualified clinical staff")
staff_groups <- sort(unique(target$staff_group))

for (sg in staff_groups) {
  d_sg <- target |> filter(staff_group == sg)

  # make plot
  p <- make_plot(d_sg, sg)

  # safe filename
  safe_sg <- sg |>
    str_replace_all("[^A-Za-z0-9]+", "_") |>
    str_replace_all("_+", "_") |>
    str_replace_all("^_|_$", "")

  outfile <- file.path(out_dir, paste0("sickness_rate_", safe_sg, ".png"))

  ggsave(
    filename = outfile,
    plot = p,
    width = 14,
    height = 6,
    dpi = 150
  )
}

message("Done. Plots written to: ", normalizePath(out_dir))

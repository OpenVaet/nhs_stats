# ============================================================================
# UKHSA COVID-19 (Nation -> England)
# Download:
#   1) PCR positivity (7-day rolling)  -> monthly mean
#   2) PCR tests (7-day total count)   -> monthly mean
#
# Outputs (created only if missing):
#   - data/tests_positivity_rates.csv      (month,value)
#   - data/tests_total_7day_mean.csv       (month,value)
#   - data/tests_positivity_and_volume.csv (month,positivity_mean_7day,tests_mean_7day_total)
# ============================================================================

if (!requireNamespace("httr2", quietly = TRUE)) install.packages("httr2")
if (!requireNamespace("jsonlite", quietly = TRUE)) install.packages("jsonlite")

library(httr2)
library(jsonlite)

dir.create("data", showWarnings = FALSE, recursive = TRUE)

base_api  <- "https://api.ukhsa-dashboard.data.gov.uk"
theme     <- "infectious_disease"
sub_theme <- "respiratory"
topic     <- "COVID-19"
geo_type  <- "Nation"
geography <- "England"

page_size <- 100

metric_url <- function(metric) {
  sprintf(
    "%s/themes/%s/sub_themes/%s/topics/%s/geography_types/%s/geographies/%s/metrics/%s",
    base_api, theme, sub_theme, topic, geo_type, geography, metric
  )
}

fetch_page <- function(url, page) {
  resp <- request(url) |>
    req_url_query(format = "json", page_size = page_size, page = page) |>
    req_perform()

  txt <- resp_body_string(resp)
  fromJSON(txt, simplifyDataFrame = TRUE)
}

fetch_all_results <- function(url) {
  first <- fetch_page(url, 1)

  if (!is.list(first) || is.null(first[["count"]]) || is.null(first[["results"]])) {
    stop("Unexpected response structure from page 1 at: ", url)
  }

  count <- as.integer(first[["count"]])
  total_pages <- ceiling(count / page_size)

  message("\n--- Fetching: ", url)
  message("Total records: ", count, " | Page size: ", page_size, " | Pages: ", total_pages)

  all_results <- vector("list", total_pages)
  all_results[[1]] <- first[["results"]]

  if (total_pages >= 2) {
    for (p in 2:total_pages) {
      dat <- fetch_page(url, p)
      res <- dat[["results"]]

      if (is.null(res) || length(res) == 0) {
        warning("Empty results on page ", p, " (continuing).")
        all_results[[p]] <- NULL
      } else {
        all_results[[p]] <- res
      }

      if (p %% 10 == 0 || p == total_pages) {
        message("Fetched page ", p, " / ", total_pages)
      }
    }
  }

  all_results <- Filter(Negate(is.null), all_results)
  df <- do.call(rbind, all_results)

  if (!is.data.frame(df) || nrow(df) == 0) {
    stop("No rows assembled for: ", url)
  }

  df
}

to_monthly_mean <- function(df) {
  # expects df has year, month, metric_value, date
  df$date <- as.Date(df$date)
  df$month_key <- sprintf("%04d-%02d", df$year, df$month)

  monthly <- aggregate(
    metric_value ~ month_key,
    data = df,
    FUN = function(x) mean(x, na.rm = TRUE)
  )

  monthly <- monthly[order(monthly$month_key), ]
  names(monthly) <- c("month", "value")
  monthly
}

write_if_missing <- function(df, out_file) {
  if (file.exists(out_file)) {
    message("File already exists: ", out_file, " (skipping)")
    return(invisible(FALSE))
  }
  write.csv(df, out_file, row.names = FALSE)
  message("Wrote: ", out_file)
  invisible(TRUE)
}

# ---- Metrics ----
metric_positivity <- "COVID-19_testing_positivity7DayRolling"
metric_tests      <- "COVID-19_testing_PCRcountByDay"

out_pos <- file.path("data", "tests_positivity_rates.csv")
out_tst <- file.path("data", "tests_total_7day_mean.csv")
out_join <- file.path("data", "tests_positivity_and_volume.csv")

# ---- Positivity ----
pos_monthly <- NULL
if (file.exists(out_pos)) {
  pos_monthly <- read.csv(out_pos, stringsAsFactors = FALSE)
} else {
  pos_daily <- fetch_all_results(metric_url(metric_positivity))
  pos_monthly <- to_monthly_mean(pos_daily)
  write_if_missing(pos_monthly, out_pos)
}

# ---- Tests ----
tst_monthly <- NULL
if (file.exists(out_tst)) {
  tst_monthly <- read.csv(out_tst, stringsAsFactors = FALSE)
} else {
  tst_daily <- fetch_all_results(metric_url(metric_tests))
  tst_monthly <- to_monthly_mean(tst_daily)
  write_if_missing(tst_monthly, out_tst)
}

# ---- Combined (month join) ----
combined <- merge(
  transform(pos_monthly, positivity_mean_7day = value)[, c("month", "positivity_mean_7day")],
  transform(tst_monthly, tests_mean_7day_total = value)[, c("month", "tests_mean_7day_total")],
  by = "month",
  all = TRUE
)

combined <- combined[order(combined$month), ]

write_if_missing(combined, out_join)

message("\nPreview combined:")
print(utils::head(combined, 24))
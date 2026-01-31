# Download UKHSA COVID-19 PCR positivity (7-day rolling) for Nation -> England
# and write MONTHLY values to data/tests_positivity_rates.csv (if not already present).
#
# Output columns: month (YYYY-MM), metric_value (monthly mean of daily rolling values)

out_file <- file.path("data", "tests_positivity_rates.csv")

if (file.exists(out_file)) {
  message("File already exists: ", out_file, "\nNothing to do.")
  quit(status = 0)
}

dir.create("data", showWarnings = FALSE, recursive = TRUE)

if (!requireNamespace("httr2", quietly = TRUE)) install.packages("httr2")
if (!requireNamespace("jsonlite", quietly = TRUE)) install.packages("jsonlite")

library(httr2)
library(jsonlite)

base_api  <- "https://api.ukhsa-dashboard.data.gov.uk"
theme     <- "infectious_disease"
sub_theme <- "respiratory"
topic     <- "COVID-19"
geo_type  <- "Nation"
geography <- "England"
metric    <- "COVID-19_testing_positivity7DayRolling"

metric_url <- sprintf(
  "%s/themes/%s/sub_themes/%s/topics/%s/geography_types/%s/geographies/%s/metrics/%s",
  base_api, theme, sub_theme, topic, geo_type, geography, metric
)

# Use a larger page_size to reduce total pages (faster).
page_size <- 100

fetch_page <- function(page) {
  resp <- request(metric_url) |>
    req_url_query(format = "json", page_size = page_size, page = page) |>
    req_perform()

  txt <- resp_body_string(resp)
  fromJSON(txt, simplifyDataFrame = TRUE)
}

# ---- Page 1: determine count + pages ----
first <- fetch_page(1)

if (!is.list(first) || is.null(first[["count"]]) || is.null(first[["results"]])) {
  stop("Unexpected response structure from page 1 at: ", metric_url)
}

count <- as.integer(first[["count"]])
total_pages <- ceiling(count / page_size)

message("Total records: ", count)
message("Page size: ", page_size)
message("Total pages: ", total_pages)

# ---- Fetch all pages ----
all_results <- vector("list", total_pages)
all_results[[1]] <- first[["results"]]

if (total_pages >= 2) {
  for (p in 2:total_pages) {
    dat <- fetch_page(p)
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
  stop("No rows assembled; cannot write CSV.")
}

# ---- Aggregate to monthly ----
# Ensure types are sensible
df$date <- as.Date(df$date)

# Create YYYY-MM key (character, as requested)
df$month_key <- sprintf("%04d-%02d", df$year, df$month)

# Monthly mean of daily rolling positivity
monthly <- aggregate(
  metric_value ~ month_key,
  data = df,
  FUN = function(x) mean(x, na.rm = TRUE)
)

# Tidy + sort
monthly <- monthly[order(monthly$month_key), ]
names(monthly) <- c("month", "value")

write.csv(monthly, out_file, row.names = FALSE)
message("Wrote monthly CSV: ", out_file)

print(head(monthly, 24))
#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package 'data.table' is required. Install with: install.packages('data.table')", call. = FALSE)
  }
  if (!requireNamespace("lubridate", quietly = TRUE)) {
    stop("Package 'lubridate' is required. Install with: install.packages('lubridate')", call. = FALSE)
  }
})

library(data.table)
library(lubridate)

root_dir <- "data/reasons"

canonical_cols <- c(
  "Date","HEE_region_code","HEE_region_name","Org Code","Org Name","Org Type",
  "FTE Days Sick","FTE Days Available","SA Rate (%)","Sort_Date"
)

# Sort_Date is optional
required_cols <- setdiff(canonical_cols, "Sort_Date")

normalize_names <- function(nm) {
  nm <- gsub("^\ufeff", "", nm, useBytes = TRUE) # strip BOM
  nm <- trimws(nm)
  nm <- gsub("\\s+", " ", nm)                    # collapse whitespace
  nm
}

# Aliases per canonical column (add to these as needed)
aliases <- list(
  "Date"               = c("Date", "DATE"),
  "HEE_region_code"    = c("HEE_region_code", "NHSE_Region_Code", "NHSE_REGION_CODE"),
  "HEE_region_name"    = c("HEE_region_name", "NHSE_Region_Name", "NHSE_REGION_NAME"),
  "Org Code"           = c("Org Code", "ORG_CODE"),
  "Org Name"           = c("Org Name", "ORG_NAME"),
  "Org Type"           = c("Org Type", "ORG_TYPE"),
  "FTE Days Sick"      = c("FTE Days Sick", "FTE_DAYS_LOST"),
  "FTE Days Available" = c("FTE Days Available", "FTE_DAYS_AVAILABLE"),
  "SA Rate (%)"        = c("SA Rate (%)", "SICKNESS_ABSENCE_RATE_PERCENT"),
  "Sort_Date"          = c("Sort_Date") # optional
)

# Find CSVs
files <- list.files(root_dir, pattern = "\\.csv$", recursive = TRUE, full.names = TRUE)
if (length(files) == 0L) stop(sprintf("No .csv files found under: %s", root_dir), call. = FALSE)

kept_list <- list()
kept_files <- character()
skipped <- data.table(file = character(), reason = character())

for (f in files) {
  hdr <- tryCatch(fread(f, nrows = 0, showProgress = FALSE), error = function(e) NULL)
  if (is.null(hdr)) {
    skipped <- rbind(skipped, list(f, "cannot read header"))
    next
  }

  cols <- normalize_names(names(hdr))

  # For each canonical col, choose the first alias that exists in this file
  pick_source <- function(canon) {
    als <- normalize_names(aliases[[canon]])
    hit <- als[als %in% cols]
    if (length(hit) == 0) NA_character_ else hit[1]
  }

  source_for_required <- vapply(required_cols, pick_source, character(1))
  missing_required <- required_cols[is.na(source_for_required)]
  if (length(missing_required) > 0L) {
    skipped <- rbind(skipped, list(f, paste0("missing required: ", paste(missing_required, collapse = ", "))))
    next
  }

  # Optional Sort_Date (may be NA)
  source_sort_date <- pick_source("Sort_Date")

  # Build select list: required + optional if present
  select_cols <- c(unname(source_for_required), if (!is.na(source_sort_date)) source_sort_date)

  dt <- tryCatch(fread(f, select = select_cols, showProgress = FALSE), error = function(e) NULL)
  if (is.null(dt)) {
    skipped <- rbind(skipped, list(f, "cannot read data"))
    next
  }

  # Normalize names as read
  setnames(dt, normalize_names(names(dt)))

  # Rename required columns -> canonical
  setnames(dt, old = unname(source_for_required), new = required_cols)

  # Rename optional Sort_Date if present, else create it
  if (!is.na(source_sort_date)) {
    setnames(dt, old = source_sort_date, new = "Sort_Date")
  } else {
    dt[, Sort_Date := NA_character_]
  }

  # Ensure all canonical columns exist & order them
  for (c in setdiff(canonical_cols, names(dt))) dt[, (c) := NA]
  dt <- dt[, ..canonical_cols]

  kept_list[[length(kept_list) + 1L]] <- dt
  kept_files <- c(kept_files, f)
}

cat("\n=== File scan ===\n")
cat(sprintf("Found:   %d\n", length(files)))
cat(sprintf("Kept:    %d\n", length(kept_files)))
cat(sprintf("Skipped: %d\n", nrow(skipped)))
if (nrow(skipped) > 0) {
  cat("\nExamples of skipped files (up to 10):\n")
  print(head(skipped, 10))
}

if (length(kept_list) == 0L) stop("No CSV files matched the required schema (with aliases).", call. = FALSE)

all_dt <- rbindlist(kept_list, use.names = TRUE, fill = TRUE)

# Parse Date values that can be either:
# - "YYYY-MMM" (e.g., "2014-APR")
# - "dd/mm/yyyy" (e.g., "30/06/2022")
parse_mixed_date <- function(x) {
  x <- trimws(x)
  x_up <- toupper(x)

  out <- rep(as.Date(NA), length(x))

  # 1) YYYY-MMM  -> set day = 01
  ok_ym <- grepl("^\\d{4}-[A-Z]{3}$", x_up)
  out[ok_ym] <- as.Date(paste0(x_up[ok_ym], "-01"), format = "%Y-%b-%d")

  # 2) dd/mm/yyyy (UK) -> lubridate::dmy
  still_na <- is.na(out) & grepl("^\\d{1,2}/\\d{1,2}/\\d{4}$", x)
  out[still_na] <- lubridate::dmy(x[still_na], quiet = TRUE)

  # 3) fallback ISO yyyy-mm-dd if it ever appears
  still_na2 <- is.na(out) & grepl("^\\d{4}-\\d{2}-\\d{2}$", x)
  out[still_na2] <- lubridate::ymd(x[still_na2], quiet = TRUE)

  out
}

d <- parse_mixed_date(all_dt[["Date"]])

# Since report is "years & months covered", normalize to month start
d_month <- lubridate::floor_date(d, unit = "month")
valid_d <- d_month[!is.na(d_month)]

cat("\n=== Coverage (years & months) ===\n")
if (length(valid_d) == 0L) {
  cat("Could not parse any dates from Date.\n")
} else {
  cov <- unique(data.table(
    year = lubridate::year(valid_d),
    month_num = lubridate::month(valid_d)
  ))
  setorder(cov, year, month_num)
  cov[, month_lbl := format(as.Date(sprintf("%04d-%02d-01", year, month_num)), "%b")]

  for (yy in unique(cov$year)) {
    cat(sprintf("%d: %s\n", yy, paste(cov[year == yy, month_lbl], collapse = ", ")))
  }
  cat(sprintf("\nEarliest date: %s\n", format(min(valid_d), "%Y-%m-%d")))
  cat(sprintf("Latest date:   %s\n", format(max(valid_d), "%Y-%m-%d")))
}

cat("\n=== Values for 'Org Type' ===\n")
org_types <- sort(unique(na.omit(trimws(all_dt[["Org Type"]]))))
if (length(org_types) == 0L) {
  cat("No non-empty 'Org Type' values found.\n")
} else {
  for (v in org_types) cat(sprintf("- %s\n", v))
}

# --- Uniqueness check on composite key ---
key_cols <- c("Date","HEE_region_code","HEE_region_name","Org Code","Org Name","Org Type")

# (Optional) strip repeated header rows that sometimes appear inside files
all_dt <- all_dt[Date != "Date" & `Org Type` != "Org Type"]

# Count duplicates
dup_summary <- all_dt[, .N, by = key_cols][N > 1]
dup_count_groups <- nrow(dup_summary)
dup_extra_rows <- if (dup_count_groups == 0) 0 else sum(dup_summary$N - 1)

cat("\n=== Uniqueness check (composite key) ===\n")
cat(sprintf("Total rows:                 %d\n", nrow(all_dt)))
cat(sprintf("Duplicate key groups:       %d\n", dup_count_groups))
cat(sprintf("Extra rows beyond unique:   %d\n", dup_extra_rows))

if (dup_count_groups > 0) {
  cat("\nTop duplicate key groups (up to 20):\n")
  print(head(dup_summary[order(-N)], 20))

  # Show full rows for the worst offending key (largest N)
  worst <- dup_summary[which.max(N)]
  cat("\nExample: full rows for most duplicated key:\n")
  print(all_dt[
    Date == worst$Date &
      HEE_region_code == worst$HEE_region_code &
      HEE_region_name == worst$HEE_region_name &
      `Org Code` == worst$`Org Code` &
      `Org Name` == worst$`Org Name` &
      `Org Type` == worst$`Org Type`
  ])
} else {
  cat("OK: all keys are unique.\n")
}

# --- Check duplicates have identical measures ---

key_cols <- c("Date","HEE_region_code","HEE_region_name","Org Code","Org Name","Org Type")

# Drop header rows that sometimes appear inside the data
all_dt <- all_dt[Date != "Date" & `Org Type` != "Org Type"]

# Convert measures to numeric safely (handles chars)
to_num <- function(x) {
  x <- trimws(x)
  x[x == ""] <- NA_character_
  as.numeric(gsub(",", "", x, fixed = TRUE))
}

all_dt[, fte_days_sick_num := to_num(`FTE Days Sick`)]
all_dt[, fte_days_avail_num := to_num(`FTE Days Available`)]

# Duplicate-key groups
dup_groups <- all_dt[, .N, by = key_cols][N > 1]

cat("\n=== Duplicate-key consistency check ===\n")
cat(sprintf("Duplicate key groups: %d\n", nrow(dup_groups)))

if (nrow(dup_groups) == 0) {
  cat("No duplicate keys found.\n")
} else {
  # For each duplicate key group, count distinct values of measures (NA ignored understandingly)
  cons <- all_dt[dup_groups, on = key_cols][
    , .(
      rows = .N,
      distinct_fte_sick = uniqueN(fte_days_sick_num),
      distinct_fte_avail = uniqueN(fte_days_avail_num)
    ),
    by = key_cols
  ]

  # Groups where measures vary
  varying <- cons[distinct_fte_sick > 1 | distinct_fte_avail > 1]

  cat(sprintf("Duplicate groups with varying FTE Days Sick:      %d\n", nrow(cons[distinct_fte_sick > 1])))
  cat(sprintf("Duplicate groups with varying FTE Days Available: %d\n", nrow(cons[distinct_fte_avail > 1])))
  cat(sprintf("Duplicate groups with any variation:              %d\n", nrow(varying)))

  if (nrow(varying) > 0) {
    cat("\nTop varying groups (up to 20):\n")
    print(head(varying[order(-rows)], 20))

    # Show full rows for the worst varying group
    w <- varying[which.max(rows)]
    cat("\nExample: full rows for the most duplicated *varying* key:\n")
    print(all_dt[
      Date == w$Date &
        HEE_region_code == w$HEE_region_code &
        HEE_region_name == w$HEE_region_name &
        `Org Code` == w$`Org Code` &
        `Org Name` == w$`Org Name` &
        `Org Type` == w$`Org Type`,
      .(Date, HEE_region_code, HEE_region_name, `Org Code`, `Org Name`, `Org Type`,
        `FTE Days Sick`, `FTE Days Available`, `SA Rate (%)`, Sort_Date)
    ])
  } else {
    cat("OK: In all duplicate-key groups, both measures are identical.\n")
  }
}

# --- Count duplicate-key groups where FTE Days Sick varies by > 5 ---

key_cols <- c("Date","HEE_region_code","HEE_region_name","Org Code","Org Name","Org Type")

# Drop repeated header rows inside data (if any)
all_dt <- all_dt[Date != "Date" & `Org Type` != "Org Type"]

# Convert to numeric safely
to_num <- function(x) {
  x <- trimws(x)
  x[x == ""] <- NA_character_
  as.numeric(gsub(",", "", x, fixed = TRUE))
}
all_dt[, fte_days_sick_num := to_num(`FTE Days Sick`)]

# Only consider duplicate keys
dup_keys <- all_dt[, .N, by = key_cols][N > 1]

# For each duplicate key group, compute min/max/range ignoring NA
sick_spread <- all_dt[dup_keys, on = key_cols][
  , .(
    rows = .N,
    non_na = sum(!is.na(fte_days_sick_num)),
    min_sick = if (all(is.na(fte_days_sick_num))) NA_real_ else min(fte_days_sick_num, na.rm = TRUE),
    max_sick = if (all(is.na(fte_days_sick_num))) NA_real_ else max(fte_days_sick_num, na.rm = TRUE)
  ),
  by = key_cols
][
  , range_sick := max_sick - min_sick
]

# Threshold
thr <- 5

groups_over_5 <- sick_spread[!is.na(range_sick) & range_sick > thr]
rows_over_5   <- all_dt[groups_over_5, on = key_cols, .N]

cat("\n=== FTE Days Sick variation in duplicate keys ===\n")
cat(sprintf("Duplicate key groups total:                 %d\n", nrow(sick_spread)))
cat(sprintf("Groups with range(FTE Days Sick) > %g:      %d\n", thr, nrow(groups_over_5)))
cat(sprintf("Rows belonging to those groups:             %d\n", rows_over_5))

# Show biggest offenders
cat("\nTop 20 ranges:\n")
print(head(groups_over_5[order(-range_sick)], 20))

cat("\nDone.\n")

# --- Build unified dataset (drop NA/0 sick + dedupe on key) ---

key_cols <- c("Date","HEE_region_code","HEE_region_name","Org Code","Org Name","Org Type")

# Drop repeated header rows that sometimes appear inside the data
all_dt <- all_dt[Date != "Date" & `Org Type` != "Org Type"]

# Mixed date parser: "YYYY-MMM" OR "dd/mm/yyyy" OR "yyyy-mm-dd"
parse_mixed_date <- function(x) {
  x <- trimws(x)
  x_up <- toupper(x)

  out <- rep(as.Date(NA), length(x))

  ok_ym <- grepl("^\\d{4}-[A-Z]{3}$", x_up)
  out[ok_ym] <- as.Date(paste0(x_up[ok_ym], "-01"), format = "%Y-%b-%d")

  still_na <- is.na(out) & grepl("^\\d{1,2}/\\d{1,2}/\\d{4}$", x)
  out[still_na] <- lubridate::dmy(x[still_na], quiet = TRUE)

  still_na2 <- is.na(out) & grepl("^\\d{4}-\\d{2}-\\d{2}$", x)
  out[still_na2] <- lubridate::ymd(x[still_na2], quiet = TRUE)

  out
}

# Parse Date and compute month start for filtering
all_dt[, date_parsed := parse_mixed_date(Date)]
all_dt[, month_start := lubridate::floor_date(date_parsed, unit = "month")]
all_dt[, year := lubridate::year(month_start)]

# Keep 2014..2025 inclusive
all_dt <- all_dt[!is.na(year) & year >= 2014 & year <= 2025]

# Monthly key as YYYY-MM (character)
all_dt[, Month := format(month_start, "%Y-%m")]

# Numeric conversion helper (handles commas, blanks)
to_num <- function(x) {
  x <- trimws(x)
  x[x == ""] <- NA_character_
  as.numeric(gsub(",", "", x, fixed = TRUE))
}

all_dt[, fte_days_sick_num  := to_num(`FTE Days Sick`)]
all_dt[, fte_days_avail_num := to_num(`FTE Days Available`)]
all_dt[, sa_rate_num        := to_num(`SA Rate (%)`)]

# Drop rows with NA or 0 on FTE Days Sick
before_rows <- nrow(all_dt)
all_dt <- all_dt[!is.na(fte_days_sick_num) & fte_days_sick_num != 0]
after_rows <- nrow(all_dt)

# Prefer keeping rows that have Sort_Date (and latest Sort_Date) when deduping
all_dt[, sort_date_num := suppressWarnings(as.numeric(trimws(Sort_Date)))]
all_dt[, has_sort_date := !is.na(sort_date_num)]

data.table::setorder(
  all_dt,
  -has_sort_date,
  -sort_date_num,
  -fte_days_sick_num,
  -fte_days_avail_num
)

# IMPORTANT: dedupe using Month instead of Date (monthly unified)
key_cols_month <- c("Month","HEE_region_code","HEE_region_name","Org Code","Org Name","Org Type")

before_dupes <- nrow(all_dt)
all_dt_unified <- unique(all_dt, by = key_cols_month)
after_dupes <- nrow(all_dt_unified)

cat("\n=== Unification summary ===\n")
cat(sprintf("Rows after year filter:                 %d\n", before_rows))
cat(sprintf("Rows after dropping NA/0 sick:          %d (dropped %d)\n", after_rows, before_rows - after_rows))
cat(sprintf("Rows before dedupe:                     %d\n", before_dupes))
cat(sprintf("Rows after dedupe on monthly key:       %d (removed %d duplicates)\n", after_dupes, before_dupes - after_dupes))

# Replace measures with numeric
all_dt_unified[, `FTE Days Sick` := fte_days_sick_num]
all_dt_unified[, `FTE Days Available` := fte_days_avail_num]
all_dt_unified[, `SA Rate (%)` := sa_rate_num]

# Optional: if you want to *replace* Date with Month, uncomment:
# all_dt_unified[, Date := Month]

# Drop helper columns (keep Month)
all_dt_unified[, c("date_parsed","month_start","year",
                   "fte_days_sick_num","fte_days_avail_num","sa_rate_num",
                   "sort_date_num","has_sort_date") := NULL]

out_path <- "data/nhs_sickness_rates_unified_2014_2025.csv"
dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)

data.table::fwrite(all_dt_unified, out_path)
cat(sprintf("\nWrote unified dataset to: %s\n", out_path))
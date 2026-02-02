#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  if (!requireNamespace("data.table", quietly = TRUE)) {
    stop("Package 'data.table' is required. Install with: install.packages('data.table')", call. = FALSE)
  }
  if (!requireNamespace("lubridate", quietly = TRUE)) {
    stop("Package 'lubridate' is required. Install with: install.packages('lubridate')", call. = FALSE)
  }
  if (!requireNamespace("readxl", quietly = TRUE)) {
    stop("Package 'readxl' is required. Install with: install.packages('readxl')", call. = FALSE)
  }
})

library(data.table)
library(lubridate)
library(readxl)

# -------------------------------------------------------------------
# CONFIG
# -------------------------------------------------------------------
root_dir <- "."  # <- set this to your repository root if needed

out_path <- "data/nhs_sickness_absence_by_reason_pct_unified.csv"

# Exact title line to find (case-insensitive match)
title_text <- "Sickness absence as a percentage of Full Time Equivalent (FTE) days lost due to all reasons:"

# -------------------------------------------------------------------
# HELPERS
# -------------------------------------------------------------------
normalize_str <- function(x) {
  x <- as.character(x)
  x <- gsub("^\ufeff", "", x, useBytes = TRUE)     # strip BOM
  x <- gsub("\u00A0", " ", x, fixed = TRUE)        # NBSP -> space
  x <- gsub("[\u2010\u2011\u2012\u2013\u2014]", "-", x) # fancy dashes -> "-"
  x <- gsub("[\r\n]+", " ", x)                     # collapse newlines
  x <- trimws(x)
  x <- gsub("\\s+", " ", x)                        # collapse whitespace
  x
}

# Parse "April 2019" etc from filename; fallback returns NA
parse_month_from_filename <- function(f) {
  b <- basename(f)
  b <- sub("\\.[Xx][Ll][Ss][Xx]$", "", b)
  b <- normalize_str(b)

  # allow separators between month and year: space, hyphen, underscore
  sep <- "[\\s\\-_]+"

  # Month YYYY (most common)
  pat1 <- paste0(
    "(?i)\\b(",
    "jan(?:uary)?|feb(?:ruary)?|mar(?:ch)?|apr(?:il)?|may|",
    "jun(?:e)?|jul(?:y)?|aug(?:ust)?|sep(?:t(?:ember)?)?|",
    "oct(?:ober)?|nov(?:ember)?|dec(?:ember)?",
    ")\\b", sep, "\\b((?:19|20)\\d{2})\\b"
  )

  # YYYY Month (fallback)
  pat2 <- paste0(
    "(?i)\\b((?:19|20)\\d{2})\\b", sep, "\\b(",
    "jan(?:uary)?|feb(?:ruary)?|mar(?:ch)?|apr(?:il)?|may|",
    "jun(?:e)?|jul(?:y)?|aug(?:ust)?|sep(?:t(?:ember)?)?|",
    "oct(?:ober)?|nov(?:ember)?|dec(?:ember)?",
    ")\\b"
  )

  m <- regexpr(pat1, b, perl = TRUE)
  if (m[1] != -1) {
    s <- regmatches(b, m)
    d <- suppressWarnings(parse_date_time(s, orders = c("B Y", "b Y")))
    if (!is.na(d)) return(as.Date(floor_date(d, "month")))
  }

  m2 <- regexpr(pat2, b, perl = TRUE)
  if (m2[1] != -1) {
    s2 <- regmatches(b, m2)
    d2 <- suppressWarnings(parse_date_time(s2, orders = c("Y B", "Y b")))
    if (!is.na(d2)) return(as.Date(floor_date(d2, "month")))
  }

  as.Date(NA)
}

parse_month_from_text <- function(txt) {
  if (is.null(txt)) return(as.Date(NA))
  txt <- normalize_str(txt)

  m <- regexpr("(?i)\\b(jan(?:uary)?|feb(?:ruary)?|mar(?:ch)?|apr(?:il)?|may|jun(?:e)?|jul(?:y)?|aug(?:ust)?|sep(?:t(?:ember)?)?|oct(?:ober)?|nov(?:ember)?|dec(?:ember)?)\\b\\s+\\b(19|20)\\d{2}\\b",
               txt, perl = TRUE)
  if (m[1] == -1) return(as.Date(NA))
  s <- regmatches(txt, m)

  d <- suppressWarnings(parse_date_time(s, orders = c("B Y", "b Y")))
  if (is.na(d)) return(as.Date(NA))
  as.Date(floor_date(d, "month"))
}

# Find the first occurrence (row,col) of the title text in a small preview block
find_title_cell <- function(preview_df, title_text) {
  mat <- as.matrix(preview_df)
  mat_chr <- matrix(normalize_str(mat), nrow = nrow(mat), ncol = ncol(mat))
  hit <- which(tolower(mat_chr) == tolower(normalize_str(title_text)), arr.ind = TRUE)
  if (nrow(hit) == 0) return(NULL)
  list(row = hit[1, "row"], col = hit[1, "col"])
}

# Extract the % table from one sheet; return a wide data.table or NULL
extract_pct_table_from_sheet <- function(f, sheet) {

  preview <- tryCatch(
    read_excel(
      f, sheet = sheet,
      range = cell_limits(c(1, 1), c(260, 80)),
      col_names = FALSE, col_types = "text",
      .name_repair = "minimal"
    ),
    error = function(e) NULL
  )
  if (is.null(preview)) return(NULL)

  mat <- as.matrix(preview)
  mat_chr <- matrix(normalize_str(mat), nrow = nrow(mat), ncol = ncol(mat))
  mat_low <- tolower(mat_chr)

  # ----------------------------------------------------------------
  # 1) Preferred: find the "percentage..." subtitle row (old format)
  #    then header is next row (this is exactly June 2022 Table 1)
  # ----------------------------------------------------------------
  subtitle_row <- which(apply(mat_low, 1, function(r) {
    any(grepl("sickness absence as a percentage", r, fixed = FALSE), na.rm = TRUE)
  }))

  header_row <- NA_integer_
  if (length(subtitle_row) > 0) {
    header_row <- subtitle_row[1] + 1L
  } else {
    # ----------------------------------------------------------------
    # 2) Fallback: find a row with many Sxx headers (newer formats)
    # ----------------------------------------------------------------
    count_sxx <- function(r) sum(grepl("^\\s*S\\d{2}\\b", r, ignore.case = TRUE), na.rm = TRUE)
    sxx_counts <- apply(mat_chr, 1, count_sxx)
    cand <- which(sxx_counts >= 5)  # needs at least 5 reason columns to be considered header-ish
    if (length(cand) == 0) return(NULL)
    header_row <- cand[1]
  }

  # Read header row
  hdr <- tryCatch(
    read_excel(
      f, sheet = sheet,
      range = cell_limits(c(header_row, 1), c(header_row, 200)),
      col_names = FALSE, col_types = "text",
      .name_repair = "minimal"
    ),
    error = function(e) NULL
  )
  if (is.null(hdr) || ncol(hdr) == 0) return(NULL)

  hdr_vec <- normalize_str(as.character(hdr[1, ]))

  # Reject Table 2 (Count): it includes "FTE days available/lost"
  if (any(grepl("full time equivalent.*days available", hdr_vec, ignore.case = TRUE), na.rm = TRUE) ||
      any(grepl("full time equivalent.*days lost", hdr_vec, ignore.case = TRUE), na.rm = TRUE)) {
    return(NULL)
  }

  non_empty <- which(!is.na(hdr_vec) & hdr_vec != "")
  if (length(non_empty) == 0) return(NULL)
  last_col <- max(non_empty)

  # Read full block (header + data)
  block <- tryCatch(
    read_excel(
      f, sheet = sheet,
      range = cell_limits(c(header_row, 1), c(header_row + 1500, last_col)),
      col_names = FALSE,
      .name_repair = "minimal"
    ),
    error = function(e) NULL
  )
  if (is.null(block) || nrow(block) < 2) return(NULL)

  block_dt <- as.data.table(block)
  setnames(block_dt, paste0("V", seq_len(ncol(block_dt))))

  header <- normalize_str(unlist(block_dt[1], use.names = FALSE))
  data   <- block_dt[-1]
  if (ncol(data) < 2) return(NULL)

  # Find first reason column (S10...)
  reason_start <- which(grepl("^\\s*S\\d{2}\\b", header, ignore.case = TRUE))
  if (length(reason_start) == 0) return(NULL)
  reason_start <- min(reason_start)

  # Staff group = all columns before reason_start (may be blank header)
  staff_cols <- paste0("V", seq_len(reason_start - 1))
  if (length(staff_cols) == 0) return(NULL)

  data[, Staff_group := normalize_str(do.call(fcoalesce, lapply(.SD, as.character))), .SDcols = staff_cols]

  # Keep only real rows
  data <- data[!is.na(Staff_group) & Staff_group != ""]
  
  # Drop spacer rows (if any): rows where all reason cells are NA
  reason_vcols <- paste0("V", reason_start:last_col)

  if (length(reason_vcols) == 1) {
    # rowSums needs >=2D, so handle single column explicitly
    data <- data[!is.na(get(reason_vcols[1]))]
  } else {
    keep <- data[, rowSums(!is.na(.SD)) > 0, .SDcols = reason_vcols]
    data <- data[keep]
  }

  # Rename reason columns
  reason_names <- header[reason_start:last_col]
  reason_names[is.na(reason_names) | reason_names == ""] <- reason_vcols[is.na(reason_names) | reason_names == ""]

  setnames(data, old = reason_vcols, new = reason_names)

  keep_cols <- c("Staff_group", reason_names)
  data[, ..keep_cols]
}

# Wide -> long normalized
wide_to_long <- function(dt_wide) {
  reason_cols <- setdiff(names(dt_wide), "Staff_group")
  m <- melt(
    dt_wide,
    id.vars = "Staff_group",
    measure.vars = reason_cols,
    variable.name = "Reason",
    value.name = "SA_rate_pct",
    variable.factor = FALSE
  )

  m[, Reason := normalize_str(Reason)]
  m[, Staff_group := normalize_str(Staff_group)]

  # split "S10 <label...>" into code + name
  m[, Reason_code := sub("^\\s*([A-Za-z0-9]+)\\b.*$", "\\1", Reason)]
  m[, Reason_name := normalize_str(sub("^\\s*[A-Za-z0-9]+\\s*", "", Reason))]

  # numeric (xlsx may already be numeric, but coerce safely)
  m[, SA_rate_pct := suppressWarnings(as.numeric(SA_rate_pct))]

  # drop blanks
  m <- m[!is.na(SA_rate_pct)]

  m[]
}

# -------------------------------------------------------------------
# MAIN
# -------------------------------------------------------------------
files <- list.files(root_dir, pattern = "\\.[Xx][Ll][Ss][Xx]$", recursive = TRUE, full.names = TRUE)
files <- files[!grepl("[/\\\\]~\\$", files)]  # ignore Excel temp files

if (length(files) == 0L) stop(sprintf("No .xlsx files found under: %s", root_dir), call. = FALSE)

kept <- list()
kept_files <- character()
skipped <- data.table(file = character(), reason = character())

for (f in files) {
  shs <- tryCatch(excel_sheets(f), error = function(e) NULL)
  if (is.null(shs) || length(shs) == 0) {
    skipped <- rbind(skipped, list(f, "cannot read sheets"))
    next
  }

  month_start <- parse_month_from_filename(f)

  if (is.na(month_start)) cat("WARN: month not parsed from filename: ", basename(f), "\n", sep = "")

  found_any <- FALSE

  # Prefer non-Notes sheets when there are multiple tabs
  is_notes <- grepl("^\\s*notes\\s*$", shs, ignore.case = TRUE) |
              grepl("notes", shs, ignore.case = TRUE)

  shs_try <- if (length(shs) > 1 && any(is_notes)) shs[!is_notes] else shs

  # 1st pass: try preferred sheets (non-Notes)
  for (sh in shs_try) {
    dt_wide <- extract_pct_table_from_sheet(f, sh)
    if (is.null(dt_wide)) next

    if (is.na(month_start)) {
      ms2 <- attr(dt_wide, "month_start")
      if (!is.null(ms2) && !is.na(ms2)) month_start <- ms2
    }

    found_any <- TRUE

    dt_long <- wide_to_long(dt_wide)
    dt_long[, Month := if (!is.na(month_start)) format(month_start, "%Y-%m") else NA_character_]
    dt_long[, source_file := f]
    dt_long[, sheet := sh]
    dt_long[, file_mtime := as.POSIXct(file.info(f)$mtime, tz = "UTC")]

    kept[[length(kept) + 1L]] <- dt_long
  }

  # 2nd pass fallback: if nothing found and we skipped Notes, try all sheets
  if (!found_any && !identical(shs_try, shs)) {
    for (sh in shs) {
      dt_wide <- extract_pct_table_from_sheet(f, sh)
      if (is.null(dt_wide)) next

      found_any <- TRUE

      dt_long <- wide_to_long(dt_wide)
      dt_long[, Month := if (!is.na(month_start)) format(month_start, "%Y-%m") else NA_character_]
      dt_long[, source_file := f]
      dt_long[, sheet := sh]
      dt_long[, file_mtime := as.POSIXct(file.info(f)$mtime, tz = "UTC")]

      kept[[length(kept) + 1L]] <- dt_long
    }
  }

  if (found_any) {
    kept_files <- c(kept_files, f)
  } else {
    skipped <- rbind(skipped, list(f, "target % table not found"))
  }
}

cat("\n=== File scan (xlsx) ===\n")
cat(sprintf("Found:   %d\n", length(files)))
cat(sprintf("Kept:    %d\n", length(unique(kept_files))))
cat(sprintf("Skipped: %d\n", nrow(skipped)))
if (nrow(skipped) > 0) {
  cat("\nExamples of skipped files (up to 10):\n")
  print(head(skipped, 10))
}

if (length(kept) == 0L) stop("No tables were extracted from any .xlsx files.", call. = FALSE)

all_dt <- rbindlist(kept, use.names = TRUE, fill = TRUE)

# Basic sanity
cat("\n=== Extract summary ===\n")
cat(sprintf("Rows extracted (long): %d\n", nrow(all_dt)))
cat(sprintf("Distinct staff groups: %d\n", uniqueN(all_dt$Staff_group)))
cat(sprintf("Distinct reason codes: %d\n", uniqueN(all_dt$Reason_code)))
cat(sprintf("Months with Month parsed from filename: %d\n", uniqueN(na.omit(all_dt$Month))))

# -------------------------------------------------------------------
# DEDUPLICATION
# Key: Month + Staff_group + Reason_code
# Keep newest file (mtime). If same mtime, keep the largest SA_rate_pct (arbitrary but deterministic).
# -------------------------------------------------------------------
key_cols <- c("Month", "Staff_group", "Reason_code")

# Order preference: newest file first
setorder(all_dt, -file_mtime, -SA_rate_pct)

before <- nrow(all_dt)
dedup <- unique(all_dt, by = key_cols)
after <- nrow(dedup)

cat("\n=== Dedup summary ===\n")
cat(sprintf("Rows before dedupe: %d\n", before))
cat(sprintf("Rows after dedupe:  %d\n", after))
cat(sprintf("Removed duplicates: %d\n", before - after))

# Report remaining duplicate keys (should be 0)
dups_left <- dedup[, .N, by = key_cols][N > 1]
cat(sprintf("Duplicate key groups after dedupe: %d\n", nrow(dups_left)))

# Final columns (drop provenance if you don’t want them)
# Keep provenance columns because they’re useful to audit; comment these drops if you want provenance kept.
dedup_out <- copy(dedup)
dedup_out[, c("source_file", "sheet", "file_mtime") := NULL]

# Write CSV
dir.create(dirname(out_path), recursive = TRUE, showWarnings = FALSE)
fwrite(dedup_out, out_path)

cat(sprintf("\nWrote normalized+deduped dataset to: %s\n", out_path))
cat("\nDone.\n")

# -------------------------------------------------------------------
# YEAR-MONTH COVERAGE (from Month column)
# -------------------------------------------------------------------
cat("\n=== Coverage (years & months) ===\n")

month_vals <- unique(na.omit(dedup_out$Month))
month_date <- suppressWarnings(as.Date(paste0(month_vals, "-01")))
month_date <- month_date[!is.na(month_date)]

if (length(month_date) == 0L) {
  cat("No parseable Month values found.\n")
} else {
  cov <- unique(data.table(
    year = lubridate::year(month_date),
    month_num = lubridate::month(month_date)
  ))
  setorder(cov, year, month_num)
  cov[, month_lbl := format(as.Date(sprintf("%04d-%02d-01", year, month_num)), "%b")]

  for (yy in unique(cov$year)) {
    cat(sprintf("%d: %s\n", yy, paste(cov[year == yy, month_lbl], collapse = ", ")))
  }

  cat(sprintf("\nEarliest month: %s\n", format(min(month_date), "%Y-%m")))
  cat(sprintf("Latest month:   %s\n", format(max(month_date), "%Y-%m")))
}
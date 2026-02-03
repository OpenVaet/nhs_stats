#!/usr/bin/env Rscript
# =============================================================================
# sickness_reason_vs_vax_doses.R
# Plot (monthly):
#   1) Reason shares of sickness absence (% of FTE days lost due to all reasons),
#      stacked by Reason_code (S10..S99)
#   2) Vaccine doses administered (monthly) as a line on RIGHT axis
#   3) COVID test positivity (monthly mean) as a red line (scaled onto LEFT axis)
#
# Inputs:
#   data/nhs_sickness_absence_by_reason_pct_unified.csv
#     Month (YYYY-MM), Staff_group, Reason_code, Reason_name, SA_rate_pct
#   data/nhs_staff_vax_model_monthly.csv
#     month (YYYY-MM), first_dose, second_dose, third_dose
#   data/tests_positivity_and_volume.csv
#     month (YYYY-MM), positivity_mean_7day, tests_mean_7day_total
#
# Outputs:
#   nhs_sickness_reason_vs_vax_positivity.png
#   data/nhs_sickness_reason_vax_pos_merged_monthly_2019_2024.csv
# =============================================================================

suppressPackageStartupMessages({
	library(tidyverse)
	library(lubridate)
	library(scales)
	library(ggrepel)
})

# ---- Force English month names for axis labels (Windows/French OS safe) ----
old_time_locale <- Sys.getlocale("LC_TIME")
on.exit(try(Sys.setlocale("LC_TIME", old_time_locale), silent = TRUE), add = TRUE)

ok <- suppressWarnings(Sys.setlocale("LC_TIME", "C"))
if (is.na(ok) || ok == "") {
	for (loc in c("English_United Kingdom.1252", "English_United States.1252",
								"en_GB.UTF-8", "en_US.UTF-8")) {
		ok <- suppressWarnings(Sys.setlocale("LC_TIME", loc))
		if (!is.na(ok) && ok != "") break
	}
}

REASON_PATH <- "data/nhs_sickness_absence_by_reason_pct_unified.csv"
VAX_PATH    <- "data/nhs_staff_vax_model_monthly.csv"
POS_PATH    <- "data/tests_positivity_and_volume.csv"

# Plot window (adjust if needed)
PLOT_START <- as.Date("2019-04-01")
PLOT_END   <- as.Date("2025-12-01")

# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------
to_num <- function(x) {
	x <- trimws(as.character(x))
	x[x == ""] <- NA_character_
	as.numeric(gsub(",", "", x, fixed = TRUE))
}

safe_cor <- function(x, y) {
	ok <- is.finite(x) & is.finite(y)
	if (sum(ok) < 3) return(NA_real_)
	cor(x[ok], y[ok])
}

# -----------------------------------------------------------------------------
# Load reason % dataset + aggregate monthly by Reason
# -----------------------------------------------------------------------------
reason_raw <- readr::read_csv(REASON_PATH, show_col_types = FALSE)

req_cols <- c("Month", "Staff_group", "Reason_code", "Reason_name", "SA_rate_pct")
missing <- setdiff(req_cols, names(reason_raw))
if (length(missing) > 0) {
	stop("Missing required columns in ", REASON_PATH, ": ", paste(missing, collapse = ", "))
}

reason_base <- reason_raw %>%
	transmute(
		month = as.Date(paste0(Month, "-01")),
		staff_group = as.character(Staff_group),
		reason_code = as.character(Reason_code),
		reason_name = as.character(Reason_name),
		sa_pct = to_num(SA_rate_pct)
	) %>%
	filter(!is.na(month)) %>%
	filter(month >= PLOT_START, month <= PLOT_END) %>%
	filter(!is.na(sa_pct))

# Prefer "All staff ..." if present, else average across staff groups
all_staff_candidates <- reason_base %>%
	filter(str_detect(tolower(staff_group), "^all\\s+staff")) %>%
	distinct(staff_group) %>%
	pull(staff_group)

use_all_staff <- length(all_staff_candidates) > 0

reason_month <- if (use_all_staff) {
	message("Using staff_group row(s) matching '^All staff' for aggregation.")
	reason_base %>%
		filter(staff_group %in% all_staff_candidates) %>%
		group_by(month, reason_code, reason_name) %>%
		summarise(sa_pct = mean(sa_pct, na.rm = TRUE), .groups = "drop")
} else {
	message("No '^All staff' row found; averaging across staff groups per month+reason.")
	reason_base %>%
		group_by(month, reason_code, reason_name) %>%
		summarise(sa_pct = mean(sa_pct, na.rm = TRUE), .groups = "drop")
}

# Normalize per month so stacks sum to 100 (handles missing reasons)
reason_stack <- reason_month %>%
	group_by(month) %>%
	mutate(
		month_sum = sum(sa_pct, na.rm = TRUE),
		sa_pct_norm = if_else(is.finite(month_sum) & month_sum > 0, sa_pct / month_sum * 100, NA_real_)
	) %>%
	ungroup() %>%
	filter(is.finite(sa_pct_norm))

# Order reasons by overall mean share (stable legend)
reason_order <- reason_stack %>%
	group_by(reason_code, reason_name) %>%
	summarise(m = mean(sa_pct_norm, na.rm = TRUE), .groups = "drop") %>%
	arrange(desc(m)) %>%
	mutate(reason_lbl = paste0(reason_code, " ", reason_name)) %>%
	pull(reason_lbl)

reason_stack <- reason_stack %>%
	mutate(reason_lbl = paste0(reason_code, " ", reason_name),
				 reason_lbl = factor(reason_lbl, levels = reason_order))

# -----------------------------------------------------------------------------
# Keep only reasons that appear in monthly Top 9 (within plot window),
# aggregate all others into "Others"
# -----------------------------------------------------------------------------

# Identify Top 9 reasons per month (by normalized %)
top9_by_month <- reason_stack %>%
	group_by(month) %>%
	arrange(desc(sa_pct_norm)) %>%
	slice_head(n = 9) %>%
	ungroup() %>%
	distinct(reason_lbl)

top_reasons_set <- top9_by_month %>% pull(reason_lbl)

# Collapse everything else to "Others"
reason_stack <- reason_stack %>%
	mutate(reason_grp = if_else(reason_lbl %in% top_reasons_set, as.character(reason_lbl), "Others")) %>%
	group_by(month, reason_grp) %>%
	summarise(sa_pct_norm = sum(sa_pct_norm, na.rm = TRUE), .groups = "drop")

# Stable legend order: by mean share, with Others last
reason_order <- reason_stack %>%
	group_by(reason_grp) %>%
	summarise(m = mean(sa_pct_norm, na.rm = TRUE), .groups = "drop") %>%
	arrange(desc(m)) %>%
	pull(reason_grp)

# Force Others at the end (even if it isn't)
reason_order <- c(setdiff(reason_order, "Others"), "Others")

reason_stack <- reason_stack %>%
	mutate(reason_grp = factor(reason_grp, levels = reason_order))

# -----------------------------------------------------------------------------
# Load vax monthly doses
# -----------------------------------------------------------------------------
vax <- readr::read_csv(VAX_PATH, show_col_types = FALSE) %>%
	transmute(
		month = as.Date(paste0(month, "-01")),
		first_dose  = to_num(first_dose),
		second_dose = to_num(second_dose),
		third_dose  = to_num(third_dose)
	) %>%
	arrange(month) %>%
	mutate(
		first_dose  = replace_na(first_dose, 0),
		second_dose = replace_na(second_dose, 0),
		third_dose  = replace_na(third_dose, 0),
		total_doses = first_dose + second_dose + third_dose
	)

# -----------------------------------------------------------------------------
# Load positivity (monthly mean)
# -----------------------------------------------------------------------------
pos <- readr::read_csv(POS_PATH, show_col_types = FALSE) %>%
	transmute(
		month = as.Date(paste0(month, "-01")),
		pos_monthly = to_num(positivity_mean_7day)
	) %>%
	arrange(month)



# -----------------------------------------------------------------------------
# Build df_line (months present in stack), join vax+pos
# -----------------------------------------------------------------------------
df_line <- reason_stack %>%
	distinct(month) %>%
	arrange(month) %>%
	left_join(vax %>% select(month, total_doses), by = "month") %>%
	left_join(pos, by = "month") %>%
	mutate(total_doses = replace_na(total_doses, 0))

# -----------------------------------------------------------------------------
# Coverage print (year-month)
# -----------------------------------------------------------------------------
cat("\n=== Coverage (years & months) ===\n")
cov <- df_line %>%
	mutate(year = year(month), month_num = month(month)) %>%
	distinct(year, month_num) %>%
	arrange(year, month_num) %>%
	mutate(month_lbl = format(as.Date(sprintf("%04d-%02d-01", year, month_num)), "%b"))

for (yy in unique(cov$year)) {
	cat(sprintf("%d: %s\n", yy, paste(cov %>% filter(year == yy) %>% pull(month_lbl), collapse = ", ")))
}
cat(sprintf("\nEarliest month: %s\n", format(min(df_line$month, na.rm = TRUE), "%Y-%m")))
cat(sprintf("Latest month:   %s\n", format(max(df_line$month, na.rm = TRUE), "%Y-%m")))

# -----------------------------------------------------------------------------
# Axis scaling
# LEFT axis  = percent shares (0..100)
# RIGHT axis = doses (scaled to left)
# Positivity scaled to LEFT axis
# -----------------------------------------------------------------------------
left_max <- 150
left_min <- -10   # small space below 0 so the 0-dose line is visible

max_doses <- max(df_line$total_doses, na.rm = TRUE)
k_dose <- if (is.finite(max_doses) && max_doses > 0) left_max / max_doses else NA_real_

df_line <- df_line %>%
	mutate(doses_scaled = if (is.finite(k_dose)) total_doses * k_dose else NA_real_)

max_pos <- max(df_line$pos_monthly, na.rm = TRUE)
k_pos <- if (is.finite(max_pos) && max_pos > 0) left_max / max_pos else NA_real_

df_line <- df_line %>%
	mutate(pos_scaled = if (is.finite(k_pos)) pos_monthly * k_pos else NA_real_)

# -----------------------------------------------------------------------------
# Markers / breakpoints
# -----------------------------------------------------------------------------
all_breaks <- tibble(
	break_num = 1:3,
	month = as.Date(c("2019-09-01", "2021-05-01", "2022-12-01"))
) %>%
	mutate(label_text = sprintf("B%d\n%s", break_num, format(month, "%b %Y")))

markers <- tibble(
	date = as.Date(c("2020-01-31", "2020-12-08", "2021-11-01", "2022-04-15")),
	short_label = c("COVID", "Vax start", "90% vax & Mandates", "Drop Mandates")
)

# Put marker labels in the top band: ALWAYS > 100 and < left_max
# (alternate slightly to reduce collisions)
markers <- markers %>%
	mutate(
		y_label = if_else(row_number() %% 2 == 1, left_max - 1.4, left_max - 3.2),
		y_label = pmin(left_max - 0.8, pmax(101, y_label))
	)

# X-axis: every 3 months + Jan separators
x_min <- floor_date(min(df_line$month, na.rm = TRUE), unit = "month")
x_max <- ceiling_date(max(df_line$month, na.rm = TRUE), unit = "month")
x_breaks_3m  <- seq(x_min, x_max, by = "3 months")
x_breaks_jan <- seq(floor_date(x_min, "year"), floor_date(x_max, "year"), by = "1 year")
grid_3m  <- tibble(x = x_breaks_3m)
grid_jan <- tibble(x = x_breaks_jan)
label_3m <- function(d) format(as.Date(d), "%Y-%m")

# -----------------------------------------------------------------------------
# Background shading
# -----------------------------------------------------------------------------
periods_bg <- tibble(
	period = c("no_tests", "corr_high", "corr_none"),
	xmin   = ymd(c(format(PLOT_START, "%Y-%m-%d"), "2020-05-01", "2022-05-01")),
	xmax   = ymd(c("2020-05-01", "2022-05-01", format(PLOT_END %m+% months(1), "%Y-%m-%d"))),
	fill   = c("#EEEEEE", "#DFF2E1", "#F6DADA")
)

# =============================================================================
# PLOT — Substack/print friendly (larger fonts, cleaner x-axis, readable labels)
# =============================================================================

# ---- Markers / breakpoints (single-line labels read better) ----
all_breaks <- tibble(
	break_num = 1:3,
	month = as.Date(c("2019-09-01", "2021-05-01", "2022-12-01"))
) %>%
	mutate(label_text = sprintf("B%d  %s", break_num, format(month, "%b %Y")))

markers <- tibble(
	date = as.Date(c("2020-01-31", "2020-12-08", "2021-11-01", "2022-04-15")),
	short_label = c("COVID", "Vax start", "90% vax & mandates", "Mandates dropped")
)

# Put marker labels in the top band: ALWAYS > 100 and < left_max
markers <- tibble(
  date = as.Date(c("2020-01-31", "2020-12-08", "2021-11-01", "2022-04-15")),
  short_label = c("COVID", "Vax start", "90% vax & mandates", "Mandates dropped")
) %>%
  mutate(y = left_max - 1.0)   # all start near the top; ggrepel spreads them

# ---- X-axis: keep only January separators; label every 6 months ----
x_min <- floor_date(min(df_line$month, na.rm = TRUE), unit = "month")
x_max <- ceiling_date(max(df_line$month, na.rm = TRUE), unit = "month")

x_breaks_jan <- seq(floor_date(x_min, "year"), floor_date(x_max, "year"), by = "1 year")
grid_jan <- tibble(x = x_breaks_jan)

# ---- Background shading ----
periods_bg <- tibble(
	period = c("no_tests", "corr_high", "corr_none"),
	xmin   = ymd(c(format(PLOT_START, "%Y-%m-%d"), "2020-05-01", "2022-05-01")),
	xmax   = ymd(c("2020-05-01", "2022-05-01", format(PLOT_END %m+% months(1), "%Y-%m-%d"))),
	fill   = c("#EEEEEE", "#DFF2E1", "#F6DADA")
)

# ---- Softer line colors (stand out but don’t fight the stack) ----
color_pos     <- "#C77C2E"   # softened orange for positivity
color_doses   <- "#4E79A7"   # muted blue for doses
color_breaks  <- "#B07AA1"   # muted purple for breakpoints
color_markers <- "#C77C2E"

# --- Fill palette that scales to any number of categories ---
n_fill <- nlevels(reason_stack$reason_grp)

fill_vals <- scales::hue_pal(l = 70, c = 90)(n_fill)  # pastel-ish, but unlimited
names(fill_vals) <- levels(reason_stack$reason_grp)

# Optional: force "Others" to grey so it reads as "residual"
if ("Others" %in% names(fill_vals)) fill_vals["Others"] <- "grey75"

p <- ggplot() +

	# Background shading
	geom_rect(
		data = periods_bg,
		aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf),
		fill = I(periods_bg$fill),
		inherit.aes = FALSE,
		alpha = 0.70
	) +

	# Year separators only (reduce clutter)
	geom_vline(
		data = grid_jan,
		aes(xintercept = x),
		color = "grey85",
		linewidth = 0.7
	) +

	# Stacked columns: reason shares (sum to 100 each month)
	geom_col(
		data = reason_stack,
		aes(x = month, y = sa_pct_norm, fill = reason_grp),
		width = 26,
		alpha = 0.90,
		colour = "white",
		linewidth = 0.18
	) +

	# Positivity scaled
	geom_line(
		data = df_line,
		aes(x = month, y = pos_scaled),
		color = color_pos, linewidth = 1.6, na.rm = TRUE
	) +

	# Doses scaled (right axis shows actual doses)
	geom_line(
		data = df_line,
		aes(x = month, y = doses_scaled),
		color = color_doses, linewidth = 1.8, na.rm = TRUE
	) +

	# Event markers
	geom_vline(
		data = markers,
		aes(xintercept = date),
		linetype = "dashed", linewidth = 1.0,
		color = color_markers, alpha = 0.70
	) +

	# Structural breaks
	geom_vline(
		data = all_breaks,
		aes(xintercept = month),
		linetype = "longdash", linewidth = 1.2,
		color = color_breaks, alpha = 0.85
	) +
	ggrepel::geom_label_repel(
	  data = markers,
	  aes(x = date, y = y, label = short_label),
	  direction = "y",
	  ylim = c(101, left_max - 0.6),     # keep labels in the top band
	  box.padding = 0.35,
	  point.padding = 0.15,
	  min.segment.length = 0,
	  max.overlaps = Inf,
	  seed = 1,

	  color = "grey15",
	  fill = alpha("white", 0.92),
	  size = 3.6,
	  fontface = "bold",
	  label.padding = unit(0.18, "lines"),
	  label.size = 0.25,
	  segment.color = color_markers,
	  segment.alpha = 0.6
	) +

	scale_x_date(
		date_breaks = "6 months",
		date_labels = "%b\n%Y",
		expand = expansion(mult = c(0.01, 0.03))
	) +

	scale_y_continuous(
		limits = c(left_min, left_max),
		labels = label_number(accuracy = 1),
		name = "Sickness absence share by reason (%)\n(normalised to sum 100 per month)",
		sec.axis = sec_axis(
			~ if (is.finite(k_dose)) pmax(0, . / k_dose) else NA_real_,
			name = "Vaccine doses administered (monthly)",
			breaks = if (is.finite(max_doses) && max_doses > 0) pretty(c(0, max_doses), n = 6) else waiver(),
			labels = label_number(scale_cut = cut_short_scale())
		)
	) +

	# Pastel fill for the stacked columns (Top 9 + Others)
	scale_fill_manual(values = fill_vals) +

	labs(
		title = "NHS sickness reasons (%) vs vaccine doses and test positivity",
		subtitle = "Stacked bars: reason shares (%) | Blue: doses (scaled; right axis) | Orange: positivity (scaled)",
		x = NULL,
		fill = "Reason (Top 9 + Others)"
	) +

	coord_cartesian(clip = "off") +

	theme_minimal(base_size = 18) +
	theme(
		plot.title = element_text(size = 22, face = "bold", margin = margin(b = 6)),
		plot.subtitle = element_text(size = 13, color = "grey35", margin = margin(b = 12)),
		axis.text.x = element_text(size = 12, angle = 0, vjust = 0.5, hjust = 0.5),
		axis.text.y = element_text(size = 13),
		axis.title.y = element_text(size = 14, margin = margin(r = 10)),
		axis.title.y.right = element_text(size = 14, margin = margin(l = 10)),
		panel.grid.minor = element_blank(),
		panel.grid.major.x = element_blank(),
		panel.grid.major.y = element_line(color = "grey88", linewidth = 0.5),
		plot.margin = margin(t = 14, r = 18, b = 18, l = 14),
		legend.position = "bottom",
		legend.title = element_text(size = 12, face = "bold"),
		legend.text  = element_text(size = 10),
		legend.key.height = unit(0.45, "cm"),
		legend.key.width  = unit(0.9, "cm")
	)

print(p)

ggsave(
	"nhs_sickness_reason_vs_vax_positivity_substack.png",
	plot = p,
	width = 18,
	height = 10,
	dpi = 320,
	bg = "white"
)

cat("\nPlot saved as 'nhs_sickness_reason_vs_vax_positivity_substack.png'\n")

# =============================================================================
# Top reasons by period (B1/B2/B3) as quick sanity
# =============================================================================
B1 <- as.Date("2019-09-01")
B2 <- as.Date("2021-05-01")
B3 <- as.Date("2022-12-01")

slice_period <- function(df, start, end, inclusive_end = TRUE) {
	if (is.na(end)) {
		df %>% filter(month > start)
	} else if (inclusive_end) {
		df %>% filter(month >= start, month <= end)
	} else {
		df %>% filter(month >= start, month < end)
	}
}

top_reasons_period <- function(start, end, label, n = 10) {
	dfp <- slice_period(reason_stack, start, end, inclusive_end = TRUE)
	dfp %>%
		group_by(reason_grp) %>%
		summarise(mean_pct = mean(sa_pct_norm, na.rm = TRUE), .groups = "drop") %>%
		arrange(desc(mean_pct)) %>%
		mutate(period = label) %>%
		slice_head(n = n) %>%
		select(period, reason_grp, mean_pct)
}

cat("\n=== Top reasons by average share (per period) ===\n")
bind_rows(
	top_reasons_period(B1, B2, "B1_to_B2 (2019-09 .. 2021-05)"),
	top_reasons_period(B2, B3, "B2_to_B3 (2021-05 .. 2022-12)"),
	top_reasons_period(B3, NA, "post_B3 (> 2022-12)")
) %>% print(n = Inf)

# =============================================================================
# Build monthly TOTAL sickness absence (% of FTE days lost due to all reasons)
# =============================================================================
reason_totals <- reason_month %>%
	group_by(month) %>%
	summarise(
		total_sa_pct = sum(sa_pct, na.rm = TRUE),
		.groups = "drop"
	)

# Add totals into df_line so downstream stats match the old workflow
df_line <- df_line %>%
	left_join(reason_totals, by = "month")

# =============================================================================
# Stats verification for narrative (B1/B2/B3 periods) - REASONS VERSION
# Writes:
#   data/narrative_reason_stats_verification_b1_b2_b3.csv
#   data/narrative_reason_top_reasons_b1_b2_b3.csv
#   data/narrative_reason_stats_verification_b1_b2_b3.txt
# =============================================================================

B1 <- as.Date("2019-09-01")
B2 <- as.Date("2021-05-01")
B3 <- as.Date("2022-12-01")

df_stats <- df_line %>%
	arrange(month) %>%
	mutate(idx = row_number() - 1L)

# slope per month (y ~ idx) inside a period
slope_per_month <- function(df, ycol) {
	y <- df[[ycol]]
	if (sum(is.finite(y)) < 3) return(NA_real_)
	coef(lm(y ~ idx, data = df))[2]
}

# peak value + month
peak_info <- function(df, ycol) {
	y <- df[[ycol]]
	if (!any(is.finite(y))) return(tibble(peak_month = as.Date(NA), peak_value = NA_real_))
	i <- which.max(y)
	tibble(peak_month = df$month[i], peak_value = y[i])
}

slice_period <- function(df, start, end, inclusive_end = TRUE) {
	if (is.na(end)) {
		df %>% filter(month > start)
	} else if (inclusive_end) {
		df %>% filter(month >= start, month <= end)
	} else {
		df %>% filter(month >= start, month < end)
	}
}

summarise_period <- function(df, label) {
	# TOTAL sickness absence (percent of FTE days lost due to all reasons)
	sick_mean <- mean(df$total_sa_pct, na.rm = TRUE)
	sick_min  <- min(df$total_sa_pct, na.rm = TRUE)
	sick_max  <- max(df$total_sa_pct, na.rm = TRUE)
	sick_sd   <- sd(df$total_sa_pct, na.rm = TRUE)
	sick_slope <- slope_per_month(df, "total_sa_pct")
	sick_peak  <- peak_info(df, "total_sa_pct")

	# Positivity
	pos_non_na <- sum(is.finite(df$pos_monthly))
	pos_n      <- nrow(df)
	pos_mean   <- mean(df$pos_monthly, na.rm = TRUE)
	pos_peak   <- peak_info(df, "pos_monthly")
	r_pos_sick <- safe_cor(df$pos_monthly, df$total_sa_pct)

	# Doses
	dose_non_zero <- sum(df$total_doses > 0, na.rm = TRUE)
	dose_mean     <- mean(df$total_doses, na.rm = TRUE)
	dose_peak     <- peak_info(df, "total_doses")
	dose_slope    <- slope_per_month(df, "total_doses")
	r_dose_sick   <- safe_cor(df$total_doses, df$total_sa_pct)

	tibble(
		period = label,
		months = nrow(df),

		sick_mean_pct = sick_mean,
		sick_min_pct  = sick_min,
		sick_max_pct  = sick_max,
		sick_sd_pct   = sick_sd,
		sick_slope_per_month = sick_slope,
		sick_peak_month = sick_peak$peak_month,
		sick_peak_pct   = sick_peak$peak_value,

		pos_months_present = pos_non_na,
		pos_months_total   = pos_n,
		pos_mean_pct       = pos_mean,
		pos_peak_month     = pos_peak$peak_month,
		pos_peak_pct       = pos_peak$peak_value,
		r_pos_vs_sick      = r_pos_sick,

		dose_months_nonzero = dose_non_zero,
		dose_mean           = dose_mean,
		dose_peak_month     = dose_peak$peak_month,
		dose_peak           = dose_peak$peak_value,
		dose_slope_per_month = dose_slope,
		r_dose_vs_sick      = r_dose_sick
	)
}

df_p1 <- slice_period(df_stats, B1, B2, inclusive_end = TRUE)
df_p2 <- slice_period(df_stats, B2, B3, inclusive_end = TRUE)
df_p3 <- slice_period(df_stats, B3, NA, inclusive_end = TRUE)  # month > B3

stats_tbl <- bind_rows(
	summarise_period(df_p1, "B1_to_B2 (2019-09 .. 2021-05)"),
	summarise_period(df_p2, "B2_to_B3 (2021-05 .. 2022-12)"),
	summarise_period(df_p3, "post_B3 (> 2022-12)")
)

# Top reasons by period (uses your reason_stack already collapsed to Top9+Others)
top_reasons_period <- function(start, end, label, n = 10) {
	dfp <- slice_period(reason_stack, start, end, inclusive_end = TRUE)
	dfp %>%
		group_by(reason_grp) %>%
		summarise(mean_pct = mean(sa_pct_norm, na.rm = TRUE), .groups = "drop") %>%
		arrange(desc(mean_pct)) %>%
		mutate(period = label) %>%
		slice_head(n = n) %>%
		select(period, reason_grp, mean_pct)
}

top_reasons_tbl <- bind_rows(
	top_reasons_period(B1, B2, "B1_to_B2 (2019-09 .. 2021-05)"),
	top_reasons_period(B2, B3, "B2_to_B3 (2021-05 .. 2022-12)"),
	top_reasons_period(B3, NA, "post_B3 (> 2022-12)")
)

# ---- Write files ----
dir.create("data", showWarnings = FALSE, recursive = TRUE)

readr::write_csv(stats_tbl,      "data/narrative_reason_stats_verification_b1_b2_b3.csv")
readr::write_csv(top_reasons_tbl,"data/narrative_reason_top_reasons_b1_b2_b3.csv")

# Also write a single text file containing the same printed tables
out_txt <- "data/narrative_reason_stats_verification_b1_b2_b3.txt"
sink(out_txt)
cat("=== Narrative stats verification (TOTAL sickness % / positivity / doses) ===\n\n")
print(stats_tbl, n = Inf)
cat("\n=== Top reasons by average share (per period) ===\n\n")
print(top_reasons_tbl %>% group_by(period) %>% slice_head(n = 10) %>% ungroup(), n = Inf)
sink()

cat("\nWrote:\n",
		"- data/narrative_reason_stats_verification_b1_b2_b3.csv\n",
		"- data/narrative_reason_top_reasons_b1_b2_b3.csv\n",
		"- ", out_txt, "\n", sep = "")
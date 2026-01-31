# reasons_pie_top10_by_staff_group.R
# - Reads all CSVs from data/reasons
# - For each STAFF_GROUP:
#   * 2024–2025 combined (all months)
#   * Top 10 REASON by total FTE_DAYS_LOST_REASON + "Other"
#   * Pie chart with white background
#   * Legend on the right showing Reason + % (so the breakdown is readable)
# - Writes: outputs/reasons/<STAFF_GROUP>.png

# ---- Packages ----
# install.packages(c("tidyverse", "readr", "lubridate", "scales"))
library(tidyverse)
library(readr)
library(lubridate)
library(scales)

# ---- I/O ----
in_dir  <- file.path("data", "reasons")
out_dir <- file.path("outputs", "reasons")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

files <- list.files(in_dir, pattern = "\\.csv$", full.names = TRUE)
if (length(files) == 0) stop("No CSV files found in: ", in_dir)

# ---- Read helper ----
read_one <- function(path) {
  readr::read_csv(
    path,
    show_col_types = FALSE,
    col_types = cols(.default = col_character())
  ) |>
    transmute(
      DATE = dmy(DATE),
      NHSE_CODE = NHSE_CODE,
      ORG_CODE = ORG_CODE,
      STAFF_GROUP = STAFF_GROUP,
      REASON = REASON,
      FTE_DAYS_LOST_REASON = as.numeric(FTE_DAYS_LOST_REASON)
    )
}

df <- map_dfr(files, read_one) |>
  filter(!is.na(DATE), !is.na(FTE_DAYS_LOST_REASON))

# ---- Filter to national totals to avoid double-counting ----
df <- df |>
  filter(
    year(DATE) %in% c(2024, 2025),
    NHSE_CODE == "All NHSE regions",
    ORG_CODE == "All organisations",
    !is.na(STAFF_GROUP), STAFF_GROUP != "",
    !is.na(REASON), REASON != ""
  )

if (nrow(df) == 0) {
  stop("No rows after filtering. Check files contain 2024/2025 and All organisations / All NHSE regions totals.")
}

# ---- Plot function (one staff group) ----
make_pie <- function(d_staff, staff_group_name) {
  by_reason <- d_staff |>
    group_by(REASON) |>
    summarise(total_lost = sum(FTE_DAYS_LOST_REASON, na.rm = TRUE), .groups = "drop") |>
    arrange(desc(total_lost))

  if (nrow(by_reason) == 0 || sum(by_reason$total_lost, na.rm = TRUE) <= 0) return(NULL)

  top10 <- by_reason |> slice_head(n = 10) |> pull(REASON)

  pie_df <- by_reason |>
    mutate(reason_group = if_else(REASON %in% top10, REASON, "Other")) |>
    group_by(reason_group) |>
    summarise(total_lost = sum(total_lost), .groups = "drop") |>
    mutate(
      pct = total_lost / sum(total_lost),
      # Legend label includes percent (clear breakdown)
      legend_label = paste0(reason_group, " — ", percent(pct, accuracy = 0.1))
    ) |>
    arrange(desc(total_lost)) |>
    # keep legend order matching slice order
    mutate(legend_label = factor(legend_label, levels = legend_label),
           reason_group = factor(reason_group, levels = reason_group))

  ggplot(pie_df, aes(x = "", y = total_lost, fill = legend_label)) +
    geom_col(width = 1, color = "white", linewidth = 0.6) +
    coord_polar(theta = "y") +
    labs(
      title = paste0("Top 10 reasons for sickness absence — ", staff_group_name),
      subtitle = "2024–2025 combined (all months, All organisations)",
      fill = NULL
    ) +
    theme_bw(base_size = 16) +
    theme(
      axis.title = element_blank(),
      axis.text  = element_blank(),
      axis.ticks = element_blank(),
      panel.grid = element_blank(),
      plot.title = element_text(size = 18, face = "bold"),
      plot.subtitle = element_text(size = 13),
      legend.position = "right",
      legend.text = element_text(size = 11),
      legend.key.size = unit(0.7, "lines")
    ) +
    guides(
      fill = guide_legend(
        override.aes = list(color = NA),
        byrow = TRUE
      )
    )
}

# ---- Loop staff groups and save ----
staff_groups <- sort(unique(df$STAFF_GROUP))

for (sg in staff_groups) {
  d_sg <- df |> filter(STAFF_GROUP == sg)

  p <- make_pie(d_sg, sg)
  if (is.null(p)) next

  safe_sg <- sg |>
    str_replace_all("[^A-Za-z0-9]+", "_") |>
    str_replace_all("_+", "_") |>
    str_replace_all("^_|_$", "")

  out_file <- file.path(out_dir, paste0("reasons_top10_", safe_sg, ".png"))

  ggsave(out_file, p, width = 14, height = 8, dpi = 150, bg = "white")
}

message("Done. Plots written to: ", normalizePath(out_dir))

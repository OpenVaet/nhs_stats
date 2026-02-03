# =============================================================================
# SIRS Scenario Simulation for NHS Sickness Absence Analysis
# Uses ACTUAL NHS observed data
# =============================================================================
#
# Tests three COVID hypotheses:
#   1. "static"      - COVID adds a constant risk factor
#   2. "exponential" - COVID risk grows (Long COVID / cumulative damage)
#   3. "attenuating" - Normal viral evolution (severe first, then milder)
#
# Combined with four vaccine efficacy scenarios: 0%, 30%, 60%, 95%
#
# =============================================================================

# =============================================================================
# SECTION 1: Load Actual Observed Data
# =============================================================================

cat("Loading observed NHS sickness data...\n")

# Read the actual sickness data
obs_raw <- read.csv("data/sickness_roll6_prediction_band_2020_2026.csv", 
                    stringsAsFactors = FALSE)

observed <- data.frame(
  month = as.Date(paste0(obs_raw$month, "-01")),
  sickness_obs = as.numeric(obs_raw$sickness_rate_100k),
  sickness_roll6 = as.numeric(obs_raw$sickness_roll6_100k),
  pred_low = as.numeric(obs_raw$sickness_pred_low_100k),
  pred_high = as.numeric(obs_raw$sickness_pred_high_100k)
)
observed <- observed[!is.na(observed$month), ]

cat("  Loaded", nrow(observed), "months of observed data\n")
cat("  Range:", format(min(observed$month)), "to", format(max(observed$month)), "\n")
cat("  Sickness range:", round(min(observed$sickness_obs, na.rm=TRUE)), "-", 
    round(max(observed$sickness_obs, na.rm=TRUE)), "per 100k\n\n")

# Also load the merged data with positivity/doses for reference
merged_raw <- read.csv("data/nhs_sickness_orgtype_vax_pos_merged_monthly_2019_2024.csv",
                       stringsAsFactors = FALSE)
merged <- data.frame(
  month = as.Date(paste0(merged_raw$month_yyyy_mm, "-01")),
  total_rate = as.numeric(merged_raw$total_rate_100k),
  doses = as.numeric(merged_raw$total_doses),
  positivity = as.numeric(merged_raw$pos_monthly)
)
merged$doses[is.na(merged$doses)] <- 0

# =============================================================================
# SECTION 2: Calibration from observed baseline
# =============================================================================

# Key dates
COVID_INTRO <- as.Date("2020-01-15")
VAX_START <- as.Date("2020-12-08")
B1_DATE <- as.Date("2019-09-01")
B2_DATE <- as.Date("2021-05-01")
B3_DATE <- as.Date("2022-12-01")

# Calculate baseline parameters from 2019 data (before COVID disruption)
baseline_data <- observed[observed$month < as.Date("2020-01-01"), ]
BASELINE_MEAN <- mean(baseline_data$sickness_obs, na.rm = TRUE)

# Monthly seasonal effects from 2019 data
month_avgs <- aggregate(sickness_obs ~ format(month, "%m"), data = baseline_data, FUN = mean)
names(month_avgs) <- c("month_num", "avg")
month_avgs$month_num <- as.integer(month_avgs$month_num)
month_avgs <- month_avgs[order(month_avgs$month_num), ]

# Center around mean to get seasonal deviations
MONTH_EFFECTS <- month_avgs$avg - BASELINE_MEAN
names(MONTH_EFFECTS) <- month_avgs$month_num

cat("Baseline parameters (from 2019 data):\n")
cat("  Mean:", round(BASELINE_MEAN), "per 100k\n")
cat("  Seasonal amplitude:", round((max(MONTH_EFFECTS) - min(MONTH_EFFECTS))/2), "per 100k\n\n")

# =============================================================================
# SECTION 3: SIRS Model (Euler integration)
# =============================================================================

run_sirs <- function(
    start_date = as.Date("2019-01-01"),
    end_date = as.Date("2025-10-01"),
    covid_scenario = "attenuating",
    vax_eff_infection = 0.60,
    vax_eff_sickness = 0.70,
    vax_final = 0.92,
    # Epi parameters
    beta_base = 0.35,
    infectious_period = 8,
    immunity_duration = 120,
    # Scenario parameters
    growth_rate = 0.4,   # exponential: annual growth
    decay_rate = 0.6,    # attenuating: annual decay
    # Sickness scaling
    severity_acute = 35000,     # acute sickness per I
    severity_persistent = 800,  # persistent elevation per cumulative exposure
    seasonal_beta = 0.25        # seasonal variation in transmission
) {
  
  gamma <- 1 / infectious_period
  omega <- 1 / immunity_duration
  
  # Time setup
  days <- seq(start_date, end_date, by = "day")
  n_days <- length(days)
  
  t_covid <- as.numeric(COVID_INTRO - start_date)
  t_vax <- as.numeric(VAX_START - start_date)
  
  # State vectors
  S <- I <- R <- numeric(n_days)
  S[1] <- 0.999; I[1] <- 0.001; R[1] <- 0
  
  # Cumulative infection for persistent effects
  cumulative_I <- numeric(n_days)
  cumulative_I[1] <- I[1]
  
  # Euler integration
  for (t in 2:n_days) {
    day <- t - 1
    doy <- as.integer(format(days[t], "%j"))
    
    # Seasonal forcing
    seasonal <- 1 + seasonal_beta * cos(2 * pi * (doy - 15) / 365)
    
    # Beta based on scenario
    if (day < t_covid) {
      beta <- 0
    } else {
      days_since <- day - t_covid
      
      beta <- switch(covid_scenario,
        "static" = beta_base * seasonal,
        "exponential" = beta_base * seasonal * (1 + growth_rate * days_since / 365),
        "attenuating" = beta_base * seasonal * exp(-decay_rate * days_since / 365),
        beta_base * seasonal
      )
    }
    
    # Vaccination effect on transmission
    if (day >= t_vax) {
      days_vax <- day - t_vax
      vax_cov <- vax_final * plogis((days_vax - 90) / 30)
      # Vaccine reduces susceptibility
      beta <- beta * (1 - vax_cov * vax_eff_infection * 0.4)
    }
    
    # SIRS equations
    new_inf <- beta * S[t-1] * I[t-1]
    recovery <- gamma * I[t-1]
    waning <- omega * R[t-1]
    
    S[t] <- S[t-1] - new_inf + waning
    I[t] <- I[t-1] + new_inf - recovery
    R[t] <- R[t-1] + recovery - waning
    
    # Clamp and normalize
    S[t] <- max(0, min(1, S[t]))
    I[t] <- max(0, min(1, I[t]))
    R[t] <- max(0, min(1, R[t]))
    total <- S[t] + I[t] + R[t]
    if (total > 0) { S[t] <- S[t]/total; I[t] <- I[t]/total; R[t] <- R[t]/total }
    
    # Cumulative exposure
    cumulative_I[t] <- cumulative_I[t-1] + I[t]
  }
  
  # Build output dataframe
  out <- data.frame(
    date = days,
    month = as.Date(format(days, "%Y-%m-01")),
    S = S, I = I, R = R,
    cumulative_I = cumulative_I
  )
  
  # Add sickness calculation
  out$month_num <- as.integer(format(out$date, "%m"))
  out$baseline <- BASELINE_MEAN + MONTH_EFFECTS[out$month_num]
  
  # Vaccination coverage
  out$days_vax <- as.numeric(out$date - VAX_START)
  out$vax_cov <- ifelse(out$date < VAX_START, 0, 
                        vax_final * plogis((out$days_vax - 90) / 30))
  
  # Severity multiplier by scenario
  out$days_covid <- as.numeric(out$date - COVID_INTRO)
  out$severity_mult <- ifelse(out$date < COVID_INTRO, 0,
    switch(covid_scenario,
      "static" = 1,
      "exponential" = 1 + growth_rate * pmax(0, out$days_covid) / 365,
      "attenuating" = exp(-decay_rate * pmax(0, out$days_covid) / 365),
      1
    )
  )
  
  # COVID sickness components
  # 1. Acute: from current infections
  out$acute_covid <- out$I * severity_acute * out$severity_mult
  
  # 2. Persistent: from cumulative exposure (Long COVID / damage accumulation)
  out$persistent_covid <- switch(covid_scenario,
    "static" = pmin(out$cumulative_I / 365, 3) * severity_persistent,
    "exponential" = pmin(out$cumulative_I / 365, 5) * severity_persistent * 1.5,
    "attenuating" = pmin(out$cumulative_I / 365, 2) * severity_persistent * 0.5,
    0
  )
  
  # Vaccine reduces symptomatic illness
  out$covid_sickness <- (out$acute_covid + out$persistent_covid) * (1 - out$vax_cov * vax_eff_sickness)
  
  # Total sickness
  out$sickness <- out$baseline + out$covid_sickness
  
  # Aggregate to monthly
  monthly <- aggregate(
    cbind(sickness, baseline, covid_sickness, I, vax_cov) ~ month,
    data = out, FUN = mean
  )
  names(monthly) <- c("month", "sickness_rate", "baseline", "covid_contrib", "I_mean", "vax_cov")
  monthly$scenario <- covid_scenario
  monthly$vax_eff <- vax_eff_infection
  
  monthly
}

# =============================================================================
# SECTION 4: Run All Scenarios
# =============================================================================

cat("Running SIRS simulations...\n")

scenarios <- c("static", "exponential", "attenuating")
vax_effs <- c(0, 0.30, 0.60, 0.95)

all_results <- list()
idx <- 1

for (scen in scenarios) {
  for (ve in vax_effs) {
    label <- paste0(scen, "_VE", ve * 100)
    cat("  ", label, "\n")
    
    res <- run_sirs(
      covid_scenario = scen,
      vax_eff_infection = ve,
      vax_eff_sickness = min(ve + 0.10, 1.0)
    )
    res$scenario_label <- label
    all_results[[idx]] <- res
    idx <- idx + 1
  }
}

all_sims <- do.call(rbind, all_results)
cat("  Total scenarios:", length(unique(all_sims$scenario_label)), "\n\n")

# =============================================================================
# SECTION 5: Compute Fit Metrics
# =============================================================================

cat("Computing fit metrics against observed data...\n")

# Merge with observed
combined <- merge(all_sims, observed[, c("month", "sickness_obs")], by = "month")
combined <- combined[!is.na(combined$sickness_rate) & !is.na(combined$sickness_obs), ]

# Metrics by scenario
scenarios_list <- unique(combined$scenario_label)
metrics <- data.frame(
  scenario = character(),
  n = integer(),
  rmse = numeric(),
  mae = numeric(),
  bias = numeric(),
  corr = numeric(),
  stringsAsFactors = FALSE
)

for (s in scenarios_list) {
  d <- combined[combined$scenario_label == s, ]
  if (nrow(d) < 5) next
  
  err <- d$sickness_rate - d$sickness_obs
  metrics <- rbind(metrics, data.frame(
    scenario = s,
    n = nrow(d),
    rmse = round(sqrt(mean(err^2)), 1),
    mae = round(mean(abs(err)), 1),
    bias = round(mean(err), 1),
    corr = round(cor(d$sickness_rate, d$sickness_obs), 3)
  ))
}

metrics <- metrics[order(metrics$rmse), ]

cat("\n")
cat(paste(rep("=", 60), collapse = ""), "\n")
cat("FIT METRICS (sorted by RMSE - lower is better)\n")
cat(paste(rep("=", 60), collapse = ""), "\n\n")
print(metrics, row.names = FALSE)

# Best fit
best <- metrics[1, ]
cat("\n*** BEST FIT:", best$scenario, "***\n")
cat("    RMSE:", best$rmse, "| MAE:", best$mae, "| Correlation:", best$corr, "\n\n")

# Interpretation
cat("INTERPRETATION:\n")
if (grepl("attenuating", best$scenario)) {
  cat("  -> Best fit is ATTENUATING scenario: COVID behaves like a normal virus\n")
  cat("     (severe initial waves, then progressively milder as it adapts)\n")
} else if (grepl("exponential", best$scenario)) {
  cat("  -> Best fit is EXPONENTIAL scenario: cumulative damage / Long COVID\n")
  cat("     (severity grows over time, consistent with persistent harm)\n")
} else {
  cat("  -> Best fit is STATIC scenario: constant added risk\n")
  cat("     (COVID adds a fixed burden that neither improves nor worsens)\n")
}

ve_match <- regmatches(best$scenario, regexpr("VE[0-9]+", best$scenario))
if (length(ve_match) > 0) {
  ve_val <- as.numeric(gsub("VE", "", ve_match))
  if (ve_val >= 60) {
    cat("  -> High vaccine efficacy (", ve_val, "%) provides best fit\n")
  } else if (ve_val >= 30) {
    cat("  -> Moderate vaccine efficacy (", ve_val, "%) provides best fit\n")
  } else {
    cat("  -> LOW/NO vaccine efficacy (", ve_val, "%) provides best fit!\n")
    cat("     This suggests vaccines may not have reduced sickness as expected.\n")
  }
}

# Save metrics
write.csv(metrics, "data/sirs_fit_metrics.csv", row.names = FALSE)

# =============================================================================
# SECTION 6: Generate Plots
# =============================================================================

cat("\nGenerating plots...\n")

# Color scheme
ve_colors <- c("0" = "#E41A1C", "30" = "#FF7F00", "60" = "#4DAF4A", "95" = "#377EB8")

# Plot 1: All scenarios faceted
png("sirs_scenarios_vs_observed.png", width = 1400, height = 1100, res = 100)
par(mfrow = c(3, 1), mar = c(4, 5, 3, 2), oma = c(2, 0, 3, 0))

xlim <- as.Date(c("2019-01-01", "2025-10-01"))
ylim <- c(3500, 7000)

for (scen in c("static", "exponential", "attenuating")) {
  
  title <- switch(scen,
    "static" = "1. Static COVID Risk (constant added burden)",
    "exponential" = "2. Exponential COVID Risk (cumulative damage / Long COVID)",
    "attenuating" = "3. Attenuating COVID Risk (normal viral adaptation)"
  )
  
  plot(NA, xlim = xlim, ylim = ylim, xlab = "", ylab = "Sickness (per 100k)",
       main = title, xaxt = "n", las = 1, cex.main = 1.1)
  
  # X axis
  years <- 2019:2025
  axis(1, at = as.Date(paste0(years, "-01-01")), labels = years)
  
  # Grid
  abline(h = seq(3500, 7000, 500), col = "grey90")
  abline(v = as.Date(paste0(years, "-01-01")), col = "grey90")
  
  # Events
  abline(v = COVID_INTRO, col = "red", lty = 2, lwd = 1.2)
  abline(v = VAX_START, col = "blue", lty = 2, lwd = 1.2)
  abline(v = B1_DATE, col = "#CC79A7", lty = 3)
  abline(v = B2_DATE, col = "#CC79A7", lty = 3)
  abline(v = B3_DATE, col = "#CC79A7", lty = 3)
  
  # Prediction band
  pred_data <- observed[!is.na(observed$pred_low), ]
  polygon(c(pred_data$month, rev(pred_data$month)),
          c(pred_data$pred_low, rev(pred_data$pred_high)),
          col = rgb(0.7, 0.7, 1, 0.3), border = NA)
  
  # Simulated lines
  for (ve in c(0, 30, 60, 95)) {
    label <- paste0(scen, "_VE", ve)
    d <- all_sims[all_sims$scenario_label == label, ]
    if (nrow(d) > 0) {
      lines(d$month, d$sickness_rate, col = ve_colors[as.character(ve)], lwd = 1.8)
    }
  }
  
  # Observed (black, prominent)
  lines(observed$month, observed$sickness_obs, col = "black", lwd = 2.5)
  points(observed$month, observed$sickness_obs, col = "black", pch = 16, cex = 0.4)
  
  # Legend on first panel only
  if (scen == "static") {
    legend("topleft", 
           legend = c("Observed", "VE = 0%", "VE = 30%", "VE = 60%", "VE = 95%"),
           col = c("black", ve_colors),
           lwd = c(2.5, 1.8, 1.8, 1.8, 1.8),
           pch = c(16, NA, NA, NA, NA),
           pt.cex = 0.6,
           cex = 0.85, bg = "white", box.lty = 0)
  }
}

mtext("SIRS Model Scenarios vs Observed NHS Sickness Absence", 
      outer = TRUE, cex = 1.4, font = 2, line = 1)
mtext("Red dashed = COVID intro | Blue dashed = Vax start | Pink dotted = Structural breaks (B1/B2/B3) | Blue band = baseline prediction interval", 
      outer = TRUE, side = 1, cex = 0.85, line = 0.5)

dev.off()
cat("  Saved: sirs_scenarios_vs_observed.png\n")

# Plot 2: Best fits comparison
png("sirs_best_fits.png", width = 1400, height = 700, res = 100)
par(mar = c(5, 5, 4, 2))

top3 <- metrics$scenario[1:min(3, nrow(metrics))]
colors3 <- c("#E41A1C", "#377EB8", "#4DAF4A")

plot(NA, xlim = xlim, ylim = ylim, xlab = "Year", ylab = "Sickness rate (per 100,000 staff)",
     main = "Best-Fitting SIRS Scenarios vs Observed NHS Sickness",
     xaxt = "n", las = 1, cex.main = 1.3)

axis(1, at = as.Date(paste0(years, "-01-01")), labels = years)
abline(h = seq(3500, 7000, 500), col = "grey90")
abline(v = as.Date(paste0(years, "-01-01")), col = "grey90")

# Events
abline(v = COVID_INTRO, col = "grey50", lty = 2)
abline(v = VAX_START, col = "grey50", lty = 2)
text(COVID_INTRO + 20, 6800, "COVID", cex = 0.8, adj = 0)
text(VAX_START + 20, 6800, "Vax", cex = 0.8, adj = 0)

# Prediction band
polygon(c(pred_data$month, rev(pred_data$month)),
        c(pred_data$pred_low, rev(pred_data$pred_high)),
        col = rgb(0.7, 0.7, 1, 0.25), border = NA)

# Observed
lines(observed$month, observed$sickness_obs, col = "black", lwd = 3)
points(observed$month, observed$sickness_obs, col = "black", pch = 16, cex = 0.5)

# Top 3 scenarios
for (i in seq_along(top3)) {
  d <- all_sims[all_sims$scenario_label == top3[i], ]
  lines(d$month, d$sickness_rate, col = colors3[i], lwd = 2.2)
}

legend("topleft",
       legend = c("Observed (NHS data)", top3),
       col = c("black", colors3[1:length(top3)]),
       lwd = c(3, rep(2.2, length(top3))),
       cex = 0.9, bg = "white")

# Add RMSE annotation
rmse_text <- paste(sapply(seq_along(top3), function(i) {
  m <- metrics[metrics$scenario == top3[i], ]
  paste0(top3[i], ": RMSE=", m$rmse)
}), collapse = "  |  ")
mtext(rmse_text, side = 1, line = 3.8, cex = 0.8)

dev.off()
cat("  Saved: sirs_best_fits.png\n")

# Plot 3, 4, 5: Individual full-screen plots for each scenario
for (scen in c("static", "exponential", "attenuating")) {
  
  filename <- paste0("sirs_scenario_", scen, ".png")
  
  title_main <- switch(scen,
    "static" = "Static COVID Risk Scenario",
    "exponential" = "Exponential COVID Risk Scenario", 
    "attenuating" = "Attenuating COVID Risk Scenario"
  )
  
  subtitle <- switch(scen,
    "static" = "COVID adds a constant burden that neither improves nor worsens over time",
    "exponential" = "Cumulative damage / Long COVID: severity grows over time",
    "attenuating" = "Normal viral adaptation: severe initial waves, then progressively milder"
  )
  
  png(filename, width = 1600, height = 900, res = 120)
  par(mar = c(5, 5, 5, 2))
  
  plot(NA, xlim = xlim, ylim = ylim, 
       xlab = "Year", ylab = "Sickness rate (per 100,000 staff)",
       main = "", xaxt = "n", las = 1)
  
  # Title and subtitle
  mtext(title_main, side = 3, line = 2.5, cex = 1.5, font = 2)
  mtext(subtitle, side = 3, line = 1, cex = 0.95, col = "grey30")
  
  # X axis
  axis(1, at = as.Date(paste0(years, "-01-01")), labels = years)
  
  # Grid
  abline(h = seq(3500, 7000, 500), col = "grey90")
  abline(v = as.Date(paste0(years, "-01-01")), col = "grey90")
  
  # Prediction band (baseline expectation)
  polygon(c(pred_data$month, rev(pred_data$month)),
          c(pred_data$pred_low, rev(pred_data$pred_high)),
          col = rgb(0.6, 0.6, 0.9, 0.25), border = NA)
  
  # Event lines with labels
  abline(v = COVID_INTRO, col = "red", lty = 2, lwd = 1.5)
  abline(v = VAX_START, col = "blue", lty = 2, lwd = 1.5)
  abline(v = B1_DATE, col = "#CC79A7", lty = 3, lwd = 1.2)
  abline(v = B2_DATE, col = "#CC79A7", lty = 3, lwd = 1.2)
  abline(v = B3_DATE, col = "#CC79A7", lty = 3, lwd = 1.2)
  
  # Event labels at top
  text(COVID_INTRO, 6900, "COVID\narrives", cex = 0.75, col = "red", pos = 4, font = 2)
  text(VAX_START, 6700, "Vax\nstarts", cex = 0.75, col = "blue", pos = 4, font = 2)
  text(B1_DATE, 6500, "B1", cex = 0.7, col = "#CC79A7", pos = 4)
  text(B2_DATE, 6500, "B2", cex = 0.7, col = "#CC79A7", pos = 4)
  text(B3_DATE, 6500, "B3", cex = 0.7, col = "#CC79A7", pos = 4)
  
  # Simulated lines for this scenario (thicker, more visible)
  lwd_sim <- 2.2
  for (ve in c(0, 30, 60, 95)) {
    label <- paste0(scen, "_VE", ve)
    d <- all_sims[all_sims$scenario_label == label, ]
    if (nrow(d) > 0) {
      lines(d$month, d$sickness_rate, col = ve_colors[as.character(ve)], lwd = lwd_sim)
    }
  }
  
  # Observed data (black, prominent)
  lines(observed$month, observed$sickness_obs, col = "black", lwd = 3)
  points(observed$month, observed$sickness_obs, col = "black", pch = 16, cex = 0.6)
  
  # Legend with RMSE values
  legend_labels <- sapply(c(0, 30, 60, 95), function(ve) {
    label <- paste0(scen, "_VE", ve)
    m <- metrics[metrics$scenario == label, ]
    if (nrow(m) > 0) {
      paste0("VE = ", ve, "% (RMSE: ", m$rmse, ")")
    } else {
      paste0("VE = ", ve, "%")
    }
  })
  
  legend("topleft",
         legend = c("Observed NHS data", legend_labels),
         col = c("black", ve_colors),
         lwd = c(3, rep(lwd_sim, 4)),
         pch = c(16, NA, NA, NA, NA),
         pt.cex = 0.8,
         cex = 0.9, bg = "white", box.lty = 1, box.col = "grey70")
  
  # Add interpretation box
  best_ve <- metrics[grepl(scen, metrics$scenario), ][1, ]
  if (nrow(best_ve) > 0 || !is.na(best_ve$scenario)) {
    interp_text <- paste0("Best fit for this scenario: VE=", 
                          gsub(paste0(scen, "_VE"), "", best_ve$scenario), 
                          "% | RMSE=", best_ve$rmse,
                          " | r=", best_ve$corr)
    mtext(interp_text, side = 1, line = 3.5, cex = 0.85, font = 3)
  }
  
  dev.off()
  cat("  Saved:", filename, "\n")
}

# =============================================================================
# SECTION 7: Summary
# =============================================================================

cat("\n")
cat(paste(rep("=", 60), collapse = ""), "\n")
cat("SIMULATION COMPLETE\n")
cat(paste(rep("=", 60), collapse = ""), "\n")
cat("\nOutput files:\n")
cat("  - sirs_scenarios_vs_observed.png (all scenarios)\n")
cat("  - sirs_best_fits.png (top 3 best-fitting scenarios)\n")
cat("  - data/sirs_fit_metrics.csv (fit statistics)\n")
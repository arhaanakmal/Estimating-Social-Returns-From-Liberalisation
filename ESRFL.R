# libraries
library(WDI)
library(Synth)
library(tidyverse)
library(ggplot2)

# variables and donor pool codes
WDI_INDICATORS <- c(imr = "SP.DYN.IMRT.IN",
                    literacy = "SE.ADT.LITR.ZS",
                    gdppc = "NY.GDP.PCAP.PP.KD",
                    urban = "SP.URB.TOTL.IN.ZS",
                    healthexp = "SH.XPD.CHEX.GD.ZS",
                    eduexp = "SE.XPD.TOTL.GD.ZS",
                    popdens = "EN.POP.DNST",
                    fertility = "SP.DYN.TFRT.IN")


DONOR_COUNTRIES_ISO2 <- c(
  "BGD",  # Bangladesh
  "NPL",  # Nepal
  "PAK",  # Pakistan
  "KHM",  # Cambodia
  "LAO",  # Lao PDR 
  "DZA",  # Algeria
  "MAR",  # Morocco
  "SYR",  # Syria
  "CMR",  # Cameroon
  "CIV",  # Cote d'Ivoire
  "SEN",  # Senegal
  "NER",  # Niger
  "BFA",  # Burkina Faso
  "MLI",  # Mali
  "TZA",  # Tanzania
  "ZMB",  # Zambia
  "MOZ",  # Mozambique
  "ETH"   # Ethiopia
)

INDIA_ISO2 <- "IN"

# time parameters
YEAR_START <- 1970
YEAR_TREATMENT <- 1991
YEAR_END <- 2015
PRE_PERIOD <- YEAR_START:(YEAR_TREATMENT - 1) 
POST_PERIOD <- YEAR_TREATMENT:YEAR_END 

# extract data
all_countries <- c(INDIA_ISO2, DONOR_COUNTRIES_ISO2)

raw <- WDI(country   = all_countries,
           indicator = WDI_INDICATORS,
           start = YEAR_START,
           end = YEAR_END,
           extra = TRUE)

# clean data 
interpolate_series <- function(x) {
  if (all(is.na(x))) return(x)
  approx(seq_along(x), x, xout = seq_along(x), method = "linear", rule = 2)$y}

panel <- raw %>%
  group_by(iso2c) %>%
  mutate(across(c(imr, literacy, gdppc, urban, healthexp,
                  eduexp, popdens, fertility), interpolate_series)) %>%
  ungroup()

# check data coverage
coverage <- panel %>%
  group_by(iso2c, country) %>%
  summarise(n_imr = sum(!is.na(imr)),
            n_literacy = sum(!is.na(literacy)),
            n_gdppc = sum(!is.na(gdppc)),
            n_urban = sum(!is.na(urban)),
            n_healthexp = sum(!is.na(healthexp)),
            n_eduexp = sum(!is.na(eduexp)),
            n_popdens = sum(!is.na(popdens)),
            n_fertility = sum(!is.na(fertility)),
            .groups = "drop")

message("Data coverage by country:")
print(coverage, n = Inf)

# add numeric country IDs (for Synth package)
country_ids <- panel %>%
  distinct(iso2c, country) %>%
  mutate(unit_id = row_number())

panel <- panel %>%
  left_join(country_ids, by = c("iso2c", "country"))

INDIA_ID <- country_ids %>% filter(iso2c == INDIA_ISO2) %>% pull(unit_id)
DONOR_IDS <- country_ids %>% filter(iso2c != INDIA_ISO2) %>% pull(unit_id)

message("India unit_id: ", INDIA_ID)
message("Donor unit_ids: ", paste(DONOR_IDS, collapse = ", "))

# Synth package helper for matrix outputs
synth_predict <- function(synth_out, dp) {
  Y0_raw <- as.matrix(dp$Y0plot)
  Y0     <- matrix(as.numeric(Y0_raw),
                   nrow = nrow(Y0_raw),
                   dimnames = list(rownames(Y0_raw), NULL))
  w      <- matrix(as.numeric(synth_out$solution.w), ncol = 1)
  Y0 %*% w
}

# scm pipeline
run_scm <- function(panel, outcome_var, unit_id_treated, unit_ids_donor,
                    year_start, year_treatment, year_end, label) {
  
  message("\n", strrep("=", 60))
  message("Running SCM for outcome: ", label)
  message(strrep("=", 60))
  
  # define time parameters
  actual_years <- sort(unique(panel$year))
  pre_period <- actual_years[actual_years <  year_treatment]
  post_period <- actual_years[actual_years >= year_treatment]
  all_years <- actual_years
  n_pre <- length(pre_period)
  outcome_match_years <- pre_period[round(seq(1, n_pre, length.out = 5))]
  
  # create dataprep object
  dp <- dataprep(foo = as.data.frame(panel),
                 predictors = c("gdppc", "urban", "healthexp",
                                "eduexp", "popdens", "fertility"),
                 predictors.op = "mean",
                 special.predictors = lapply(outcome_match_years, function(yr){
                   list(outcome_var, yr, "mean")}),
                 dependent = outcome_var,
                 unit.variable = "unit_id",
                 time.variable = "year",
                 treatment.identifier = unit_id_treated,
                 controls.identifier = unit_ids_donor,
                 time.predictors.prior = pre_period,
                 time.optimize.ssr = pre_period,
                 unit.names.variable = "country",
                 time.plot = all_years)
  stopifnot(
    "dataprep returned NULL Y1plot - check that panel years match time.plot" =
      !is.null(dp$Y1plot),
    "dataprep returned NULL Y0plot - check donor IDs and panel coverage" =
      !is.null(dp$Y0plot)
  )
  
  # fit synthetic control
  synth_out <- synth(dp, verbose = FALSE)
  
  # extract weights
  weights_tbl <- as_tibble(synth_out$solution.w, rownames = "unit_id") %>%
    mutate(unit_id = as.integer(unit_id)) %>%
    left_join(country_ids, by = "unit_id") %>%
    rename(weight = w.weight) %>%
    filter(weight > 0.001) %>%
    arrange(desc(weight))
  
  message("Synthetic control weights (", label, "):")
  print(weights_tbl %>% select(country, iso2c, weight))
  
  # extract RMSPE
  synth_vals <- synth_predict(synth_out, dp)
  gaps <- as.numeric(dp$Y1plot) - as.numeric(synth_vals)
  names(gaps) <- rownames(dp$Y1plot)
  pre_rmspe <- sqrt(mean(gaps[as.character(pre_period)]^2, na.rm = TRUE))
  message("Pre-treatment RMSPE: ", round(pre_rmspe, 4))
  
  # create gap series
  gap_series <- tibble(year = as.integer(rownames(dp$Y1plot)),
                       actual = as.numeric(dp$Y1plot),
                       synth = as.numeric(synth_vals),
                       gap = actual - synth)
  
  list(dataprep = dp,
       synth_out = synth_out,
       weights = weights_tbl,
       gap_series = gap_series,
       pre_rmspe = pre_rmspe,
       label = label,
       outcome = outcome_var,
       pre_period = pre_period,
       post_period = post_period)
}

# run scm for imr, literacy
scm_imr <- run_scm(panel, "imr", INDIA_ID, DONOR_IDS, YEAR_START,
                   YEAR_TREATMENT, YEAR_END, "Infant Mortality Rate")

scm_literacy <- run_scm(panel, "literacy", INDIA_ID, DONOR_IDS, YEAR_START,
                        YEAR_TREATMENT, YEAR_END, "Adult Literacy Rate")

# perumatation inference from placebos
run_placebos <- function(scm_result, panel, country_ids, year_treatment) {
  
  message("Running placebo tests for: ", scm_result$label)
  
  india_pre_rmspe <- scm_result$pre_rmspe
  all_years <- c(scm_result$pre_period, scm_result$post_period)
  year_start <- min(scm_result$pre_period)
  year_end <- max(scm_result$post_period)
  outcome_var <- scm_result$outcome
  
  placebo_gaps <- list()
  
  for (placebo_id in DONOR_IDS) {
    
    other_donors <- setdiff(DONOR_IDS, placebo_id)
    controls_for_placebo <- c(INDIA_ID, other_donors)
    
    tryCatch({
      
      actual_years_pl <- sort(unique(panel$year))
      pre_period <- actual_years_pl[actual_years_pl <  year_treatment]
      n_pre <- length(pre_period)
      outcome_match_years <- pre_period[round(seq(1, n_pre, length.out = 5))]
      
      dp_pl <- dataprep(foo = as.data.frame(panel),
                        predictors = c("gdppc", "urban", "healthexp",
                                       "eduexp", "popdens", "fertility"),
                        predictors.op = "mean",
                        special.predictors = lapply(outcome_match_years, 
                                                    function(yr) {
                                                      list(outcome_var, yr, 
                                                           "mean")}),
                        dependent = outcome_var,
                        unit.variable = "unit_id",
                        time.variable = "year",
                        treatment.identifier = placebo_id,
                        controls.identifier = controls_for_placebo,
                        time.predictors.prior = pre_period,
                        time.optimize.ssr = pre_period,
                        unit.names.variable = "country",
                        time.plot = actual_years_pl)
      
      synth_pl <- synth(dp_pl, verbose = FALSE)
      synth_vals_pl <- synth_predict(synth_pl, dp_pl)
      gaps_pl <- as.numeric(dp_pl$Y1plot) - as.numeric(synth_vals_pl)
      names(gaps_pl) <- rownames(dp_pl$Y1plot)
      plot_years_pl <- as.integer(rownames(dp_pl$Y1plot))
      post_years_pl <- plot_years_pl[plot_years_pl >= year_treatment]
      pre_rmspe_pl <- sqrt(mean(gaps_pl[as.character(pre_period)]^2,
                                na.rm = TRUE))
      post_rmspe_pl <- sqrt(mean(gaps_pl[as.character(post_years_pl)]^2,
                                 na.rm = TRUE))
      
      placebo_gaps[[as.character(placebo_id)]] <- tibble(
        unit_id = placebo_id,
        year = plot_years_pl,
        gap = as.numeric(gaps_pl),
        pre_rmspe = pre_rmspe_pl,
        post_rmspe = post_rmspe_pl,
        rmspe_ratio = post_rmspe_pl / pre_rmspe_pl)
      
    }, error = function(e) {
      message("Placebo failed for unit_id", placebo_id, ": ", e$message)
    })
  }
  
  placebos_df <- bind_rows(placebo_gaps) %>%
    left_join(country_ids, by = "unit_id")
  
  # compute india's RMSPE ratio
  india_gaps <- scm_result$gap_series
  india_post_rmspe <- sqrt(mean(india_gaps$gap[india_gaps$year >= year_treatment]^2))
  india_ratio <- india_post_rmspe / india_pre_rmspe
  
  message("India post/pre RMSPE ratio: ", round(india_ratio, 3))
  
  # drop placebos with pre-treatment RMSPE > 2x India's (standard practice)
  threshold <- 2 * india_pre_rmspe
  placebos_clean <- placebos_df %>% filter(pre_rmspe <= threshold)
  n_dropped <- n_distinct(placebos_df$unit_id) - n_distinct(placebos_clean$unit_id)
  message("Placebo units dropped (pre-RMSPE > 2x India's): ", n_dropped)
  
  # p-value: proportion of placebos with RMSPE ratio >= India's
  ratios_clean <- placebos_clean %>%
    distinct(unit_id, rmspe_ratio) %>%
    pull(rmspe_ratio)
  
  p_val <- mean(ratios_clean >= india_ratio)
  message("Permutation p-value: ", round(p_val, 3),
          "  (", sum(ratios_clean >= india_ratio), " of ",
          length(ratios_clean), " placebos >= India)")
  
  list(placebos_df = placebos_df,
       placebos_clean = placebos_clean,
       india_ratio = india_ratio,
       p_value = p_val)
}

placebos_imr <- run_placebos(scm_imr, panel, country_ids, YEAR_TREATMENT)

placebos_literacy <- run_placebos(scm_literacy, panel, country_ids, 
                                  YEAR_TREATMENT)

# leave one out analysis
run_loo <- function(scm_result, panel, country_ids, year_start, year_end) {
  
  message("Leave-one-out robustness for: ", scm_result$label)
  
  # countries with non-trivial weight in the main solution
  loo_units <- scm_result$weights %>%
    filter(weight > 0.05) %>%
    pull(unit_id)
  
  loo_gaps  <- list()
  failed    <- character(0)
  
  for (drop_id in loo_units) {
    drop_name  <- country_ids %>% filter(unit_id == drop_id) %>% pull(country)
    donors_loo <- setdiff(DONOR_IDS, drop_id)
    actual_years_loo    <- sort(unique(panel$year))
    pre_period          <- actual_years_loo[actual_years_loo < YEAR_TREATMENT]
    n_pre               <- length(pre_period)
    outcome_match_years <- pre_period[round(seq(1, n_pre, length.out = 5))]
    
    tryCatch({
      dp_loo <- dataprep(
        foo                   = as.data.frame(panel),
        predictors            = c("gdppc", "urban", "healthexp", 
                                  "eduexp", "popdens", "fertility"),
        predictors.op         = "mean",
        special.predictors    = lapply(outcome_match_years, function(yr) {
          list(scm_result$outcome, yr, "mean")
        }),
        dependent             = scm_result$outcome,
        unit.variable         = "unit_id",
        time.variable         = "year",
        treatment.identifier  = INDIA_ID,
        controls.identifier   = donors_loo,
        time.predictors.prior = pre_period,
        time.optimize.ssr     = pre_period,
        unit.names.variable   = "country",
        time.plot             = actual_years_loo
      )
      synth_loo  <- synth(dp_loo, verbose = FALSE)
      synth_vals <- as.numeric(synth_predict(synth_loo, dp_loo))
      plot_years_loo <- as.integer(rownames(dp_loo$Y1plot))
      
      loo_gaps[[drop_name]] <- tibble(
        year      = plot_years_loo,
        synth_loo = synth_vals,
        dropped   = drop_name
      )
      message("Completed: dropping", drop_name)
      
    }, error = function(e) {
      failed <<- c(failed, drop_name)
      message("Skipped (unfeasible): dropping", drop_name)
    })
  }
  
  if (length(failed) > 0) {
    message("LOO iterations not shown in plot (optimiser failed): ",
            paste(failed, collapse = ", "))
    message("This typically means the synthetic control relies heavily on ",
            "these donors and the remaining pool cannot form a valid counterfactual.")
  }
  
  bind_rows(loo_gaps)
}

loo_imr      <- run_loo(scm_imr,      panel, country_ids, YEAR_START, YEAR_END)
loo_literacy <- run_loo(scm_literacy, panel, country_ids, YEAR_START, YEAR_END)


# plots
theme_paper <- theme_minimal(base_size = 12) +
  theme(
    panel.grid.minor   = element_blank(),
    plot.title         = element_text(face = "bold", size = 13),
    plot.subtitle      = element_text(size = 10, color = "grey40"),
    legend.position    = "bottom",
    axis.line          = element_line(color = "grey60"),
    strip.text         = element_text(face = "bold")
  )

# figure 1: india vs synthetic india
plot_trends <- function(scm_result) {
  df <- scm_result$gap_series %>%
    pivot_longer(c(actual, synth), names_to = "series", values_to = "value") %>%
    mutate(series = recode(series,
                           actual = "India (Actual)",
                           synth  = "Synthetic India"))
  
  ggplot(df, aes(x = year, y = value, linetype = series, color = series)) +
    geom_line(linewidth = 0.9) +
    geom_vline(xintercept = YEAR_TREATMENT - 0.5, linetype = "dashed",
               color = "grey30", linewidth = 0.7) +
    annotate("text", x = YEAR_TREATMENT + 0.5,
             y = max(df$value, na.rm = TRUE) * 0.98,
             label = "1991\nReforms", hjust = 0, size = 3.2, color = "grey30") +
    scale_color_manual(values = c("India (Actual)" = "#c0392b",
                                  "Synthetic India" = "#2980b9")) +
    scale_linetype_manual(values = c("India (Actual)" = "solid",
                                     "Synthetic India" = "dashed")) +
    labs(
      title    = paste("Figure 1:", scm_result$label),
      subtitle = "India vs. Synthetic India, 1980-2015",
      x        = NULL, y        = scm_result$label,
      color    = NULL, linetype = NULL
    ) +
    theme_paper
}

fig1_imr      <- plot_trends(scm_imr)
fig1_literacy <- plot_trends(scm_literacy)

# fig 2: treatment effect gap (India minus synthetic)
plot_gap <- function(scm_result) {
  df <- scm_result$gap_series
  
  ggplot(df, aes(x = year, y = gap)) +
    geom_hline(yintercept = 0, color = "grey50", linewidth = 0.6) +
    geom_vline(xintercept = YEAR_TREATMENT - 0.5, linetype = "dashed",
               color = "grey30", linewidth = 0.7) +
    geom_line(color = "#c0392b", linewidth = 0.9) +
    geom_ribbon(aes(ymin = pmin(gap, 0), ymax = 0),
                fill = "#2980b9", alpha = 0.15) +
    geom_ribbon(aes(ymin = 0, ymax = pmax(gap, 0)),
                fill = "#c0392b", alpha = 0.15) +
    labs(
      title    = paste("Figure 2: Treatment Effect Gap -", scm_result$label),
      subtitle = "India minus Synthetic India (positive = India above synthetic)",
      x        = NULL,
      y        = paste("Gap in", scm_result$label)
    ) +
    theme_paper
}

fig2_imr      <- plot_gap(scm_imr)
fig2_literacy <- plot_gap(scm_literacy)

# fig 3: placebo distribution
plot_placebos <- function(scm_result, placebos_result) {
  
  india_gap_df <- scm_result$gap_series %>%
    select(year, gap) %>%
    mutate(country = "India", is_india = TRUE)
  
  placebo_gap_df <- placebos_result$placebos_clean %>%
    select(year, gap, country) %>%
    mutate(is_india = FALSE)
  
  combined <- bind_rows(india_gap_df, placebo_gap_df)
  
  ggplot() +
    geom_line(data = filter(combined, !is_india),
              aes(x = year, y = gap, group = country),
              color = "grey70", linewidth = 0.4, alpha = 0.7) +
    geom_line(data = filter(combined, is_india),
              aes(x = year, y = gap),
              color = "#c0392b", linewidth = 1.1) +
    geom_hline(yintercept = 0, linetype = "dotted", color = "grey40") +
    geom_vline(xintercept = YEAR_TREATMENT - 0.5, linetype = "dashed",
               color = "grey30", linewidth = 0.7) +
    annotate("text", x = YEAR_START + 1,
             y = max(combined$gap, na.rm = TRUE) * 0.92,
             label = "India", color = "#c0392b",
             fontface = "bold", size = 3.5) +
    labs(
      title    = paste("Figure 3: Placebo Gaps -", scm_result$label),
      subtitle = paste0("Red = India; Grey = donor placebos (pre-RMSPE ≤ 2× India's). ",
                        "p = ", round(placebos_result$p_value, 3)),
      x = NULL, y = paste("Gap in", scm_result$label)
    ) +
    theme_paper
}

fig3_imr      <- plot_placebos(scm_imr,      placebos_imr)
fig3_literacy <- plot_placebos(scm_literacy, placebos_literacy)

# fig 4: leave-one-out robustness
plot_loo <- function(scm_result, loo_df) {
  
  base_synth <- scm_result$gap_series %>% select(year, synth)
  actual_df  <- scm_result$gap_series %>% select(year, actual)
  
  ggplot() +
    geom_line(data = loo_df,
              aes(x = year, y = synth_loo, group = dropped),
              color = "grey70", linewidth = 0.5, alpha = 0.8) +
    geom_line(data = base_synth,
              aes(x = year, y = synth),
              color = "#2980b9", linewidth = 0.9, linetype = "dashed") +
    geom_line(data = actual_df,
              aes(x = year, y = actual),
              color = "#c0392b", linewidth = 1.0) +
    geom_vline(xintercept = YEAR_TREATMENT - 0.5, linetype = "dashed",
               color = "grey30", linewidth = 0.7) +
    labs(
      title    = paste("Figure 4: Leave-One-Out Robustness -", scm_result$label),
      subtitle = "Red = India actual; Blue dashed = main synthetic; Grey = LOO synthetics",
      x = NULL, y = scm_result$label
    ) +
    theme_paper
}

fig4_imr      <- plot_loo(scm_imr,      loo_imr)
fig4_literacy <- plot_loo(scm_literacy, loo_literacy)


# output tables
# table 1: predictor balance
print_balance <- function(scm_result) {
  message("Predictor Balance: ", scm_result$label, " ──")
  tab <- synth.tab(dataprep.res = scm_result$dataprep,
                   synth.res    = scm_result$synth_out)
  print(tab$tab.pred)   
  message("Loss (MSPE)")
  print(tab$tab.loss)
}

print_balance(scm_imr)
print_balance(scm_literacy)

# table 2: country weights
message("SCM Weights: IMR")
print(scm_imr$weights)

message("SCM Weights: Literacy")
print(scm_literacy$weights)

# table 3: RMSPE ratios for inference
rmspe_summary <- tibble(
  Outcome           = c("IMR", "Literacy"),
  India_pre_RMSPE   = c(scm_imr$pre_rmspe,      scm_literacy$pre_rmspe),
  India_RMSPE_ratio = c(placebos_imr$india_ratio, placebos_literacy$india_ratio),
  P_value           = c(placebos_imr$p_value,     placebos_literacy$p_value)
)
message("RMSPE-Based Inference")
print(rmspe_summary)

# exporting plots
output_dir <- " " # add path here
dir.create(output_dir, showWarnings = FALSE)

ggsave(file.path(output_dir, "fig1_imr_trends.pdf"),
       fig1_imr,      width = 8, height = 5)
ggsave(file.path(output_dir, "fig1_literacy_trends.pdf"),
       fig1_literacy, width = 8, height = 5)
ggsave(file.path(output_dir, "fig2_imr_gap.pdf"),
       fig2_imr,      width = 8, height = 5)
ggsave(file.path(output_dir, "fig2_literacy_gap.pdf"),
       fig2_literacy, width = 8, height = 5)
ggsave(file.path(output_dir, "fig3_imr_placebos.pdf"),
       fig3_imr,      width = 8, height = 5)
ggsave(file.path(output_dir, "fig3_literacy_placebos.pdf"),
       fig3_literacy, width = 8, height = 5)
ggsave(file.path(output_dir, "fig4_imr_loo.pdf"),
       fig4_imr,      width = 8, height = 5)
ggsave(file.path(output_dir, "fig4_literacy_loo.pdf"),
       fig4_literacy, width = 8, height = 5)

message("All outputs saved to: ", output_dir, "/")
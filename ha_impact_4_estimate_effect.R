#..........................................................................................
### ++++++++ ESTIMATION OF THE IMPACT OF HUMANITARIAN ASSISTANCE ON MORTALITY +++++++++ ###
#..........................................................................................

#..........................................................................................
## ---------- R CODE TO PERFORM PROPENSITY SCORE MATCHING AND ESTIMATE EFFECT----------- ##
#..........................................................................................

                                          # Written by Francesco Checchi, LSHTM (Jan 2022)
                                          # francesco.checchi@lshtm.ac.uk 



#.........................................................................................
### Reading in required files
#.........................................................................................
    
  #...................................
  ## Predictor data (from earlier code)
   
    # Read values of predictors for model fitting
    x_obs <- read.csv(paste(country, "_x_obs.csv", sep=""), sep="," )


  #...................................
  ## Survey observation datasets (from earlier code)
    
    # Read household-level survey observations
    hh_y_obs <- read.csv(paste(country, "_hh_obs.csv", sep=""), sep="," )

   

#.........................................................................................
### Merge mortality and predictor data
#.........................................................................................
  
  #...................................
  ## Select relevant columns in the data depending on which dependent variable is being modelled
      
    # if CDR is being estimated...
    if (y_hat == "cdr") {
      hh_y_obs <- hh_y_obs[, c("survey_id", "Cluster", "stratum", "n_died", "ptime", "quality_score", "sampling_coverage")]
    }
      
    # if U5DR is being estimated...
    if (y_hat == "cdr_u5") {
      hh_y_obs <- hh_y_obs[, c("survey_id", "Cluster", "stratum", "n_died_u5", "ptime_u5", "quality_score", "sampling_coverage")]
        
    # rename as for CDR to simplify code below
      colnames(hh_y_obs) <- c("survey_id", "Cluster", "stratum", "n_died", "ptime", "quality_score", "sampling_coverage")
    }

  #...................................
  ## Household-level observations: aggregate by survey cluster   
  hh_y_obs[, "n_obs"] <- 1
  hh_y_obs <- aggregate(hh_y_obs[, c("n_died", "ptime", "quality_score", "sampling_coverage", "n_obs")], 
    by = hh_y_obs[, c("survey_id", "stratum", "Cluster")], FUN = sum)
  hh_y_obs[, c("quality_score", "sampling_coverage")] <-   hh_y_obs[, c("quality_score", "sampling_coverage")] /
    hh_y_obs[, "n_obs"]
          
  #...................................
  ## Household-level observations: merge dependent and predictor variables      
  hh_obs <- merge(hh_y_obs, x_obs, by = c("survey_id", "stratum") )
  
    # Eliminate any observations with zero person-time (if U5DR is being analysed, this could occur)
    hh_obs <- subset(hh_obs, ptime > 0)
            
        
#.........................................................................................
### Prepare and explore the dataset
#.........................................................................................
    
  #...................................
  ## Calculate importance weights for survey observations
  hh_obs[, "wt"] <- hh_obs[, "quality_score"] * hh_obs[, "sampling_coverage"]
   
  #...................................
  ## Define, explore and transform exposure if needed
    # Define exposure
    exposure <- "total_actors_rate"
    
    # Explore exposure
    f_hist(exposure, hh_obs, c(NA, NA))
    nrow(hh_obs[which(hh_obs[, exposure] == 0), ])
    
    # Log transform and redefine
    hh_obs[, paste("log_", exposure, sep = "")] <- log(hh_obs[, exposure] + 0.01)
    f_hist(paste("log_", exposure, sep = ""), hh_obs, c(NA, NA))
    exposure <- "log_total_actors_rate"
    
    # # Remove outliers
    # hh_obs[, exposure] <- ifelse(hh_obs[, exposure] < 0, NA, hh_obs[, exposure])
    # f_hist(exposure, hh_obs, c(NA, NA))
    # 
  #...................................
  ## Identify confounders and modify them if needed
    
    # Livelihood to numeric
    hh_obs[, "lhz_cat"] <- as.numeric(factor(hh_obs[, "lhz"]))
    table(hh_obs[, c("lhz", "lhz_cat")])  
  
    # List of confounders to condition on (based on DAG)
    confounders <- c("lhz_cat", "prop_idp_stratum", "dep_rate_stratum", "prop_idp_stratum_lag3",
      "dep_rate_stratum_lag3", "acled_event_rate", "acled_event_rate_lag3",
      "rainfall_rollmean", "ndvi_rollmean","rainfall_rollmean_lag3", "ndvi_rollmean_lag3",
      "tot_wage_cereal_smooth_lag3", "water_price_smooth_rollmean_lag3")
      
    # # Normalise all confounders to 0-1 scale
    # for (i in confounders) {hh_obs[, i] <- hh_obs[, i] / max(hh_obs[, i])}
    
  #...................................
  ## Remove unnecessary variables and missing data

    # Retain only needed variables
    hh_obs <- hh_obs[, c("n_died", "ptime", "wt", "survey_id", "stratum", "Cluster",  exposure, confounders)]
       
    # Remove missing data
    hh_obs <- na.omit(hh_obs)

    # How much data are left?
    length(unique(hh_obs[, "Cluster"]))
    length(unique(hh_obs[, "survey_id"]))    
    
  # #...................................
  # ## Restrict data to common support / overlap
  #   # (see Bia et al. https://journals.sagepub.com/doi/10.1177/1536867X1401400307 )
  # 
  #   # Define formula for exposure as a function of confounders
  #   formula_exp <- formula(paste(paste(exposure)," ~ ", paste(confounders, collapse = " + ")))
  # 
  #   a <- overlap_fun(Y = n_died, treat = log_total_actors_rate, treat_formula = formula_exp, data_set = hh_obs,
  #     n_class = 3, treat_mod = "Normal")
  #      
  
              
#.........................................................................................
### Compute Generalised Propensity Scores 
#.........................................................................................
    # based on Hirano & Imben (2004), https://www.google.com/url?sa=t&rct=j&q=&esrc=s&source=web&cd=&ved=2ahUKEwjozLahmOn1AhXGEcAKHQjOBQsQFnoECAMQAQ&url=https%3A%2F%2Fwww.math.mcgill.ca%2Fdstephens%2FPSMMA%2FArticles%2FHIrano-Imbens-2004.pdf&usg=AOvVaw1fo5Um-exkrcC-rwh2ozZb
    # and Austin (2017), https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5969262/#sim7615-bib-0003 
    
  #...................................
  ## Compute Generalised Propensity Scores
    # Define formula for exposure as a function of confounders
    formula_exp <- formula(paste(paste(exposure)," ~ ", paste(confounders, collapse = " + ")))
    
    # Linear model of exposure as a function of confounders
    fit_exp <- lm(formula = formula_exp, data = hh_obs)
    summary(fit_exp)
    
    # Compute GPS values
    hh_obs[, "gps"] <- dnorm(x = hh_obs[, exposure], mean = fitted.values(fit_exp), 
      sd = summary(fit_exp)$sigma)

  #...................................
  ## Compute GPS stabilised weights
    
    # Numerators
    hh_obs[, "wt_num"] <- dnorm(x = (hh_obs[, exposure] - mean(hh_obs[, exposure]) ) / sd(hh_obs[, exposure]), 
      mean = 0, sd = 1)
    
    # Weights
    hh_obs[, "gps_wt"] <- hh_obs[, "wt_num"] / hh_obs[, "gps"]
      # truncate to avoid extreme weights
      hh_obs[, "gps_wt"] <- ifelse(hh_obs[, "gps_wt"] > quantile(hh_obs[, "gps_wt"], probs = 0.99),
        quantile(hh_obs[, "gps_wt"], probs = 0.99), hh_obs[, "gps_wt"] )
      hh_obs[, "gps_wt"] <- ifelse(hh_obs[, "gps_wt"] < quantile(hh_obs[, "gps_wt"], probs = 0.01),
        quantile(hh_obs[, "gps_wt"], probs = 0.01), hh_obs[, "gps_wt"] )
    
      
#.........................................................................................
### Check balance of confounders with and without adjustment by GPS
#.........................................................................................
    # based on Hirano & Imben (2004), https://www.google.com/url?sa=t&rct=j&q=&esrc=s&source=web&cd=&ved=2ahUKEwjozLahmOn1AhXGEcAKHQjOBQsQFnoECAMQAQ&url=https%3A%2F%2Fwww.math.mcgill.ca%2Fdstephens%2FPSMMA%2FArticles%2FHIrano-Imbens-2004.pdf&usg=AOvVaw1fo5Um-exkrcC-rwh2ozZb
    # with adjustment by GPS, correlations of confounders with outcome should be lower (ideally all < 0.10)
     
  #...................................
  ## Check balance of confounders with adjustment by GPS

    # Divide dataset into quantiles (groups) by exposure, compute medians of each quantile
    exp_q <- quantile(hh_obs[, exposure], probs = seq(0, 1, by = 0.25)) 
      # 4 groups arbitrarily (also avoids data sparsity, which is considerable with 5 groups)
    hh_obs[, "exp_q"] <- cut(hh_obs[, exposure], breaks = exp_q, labels = 1:(length(exp_q) - 1), 
      ordered_result = TRUE, include.lowest = TRUE)
    exp_mids <- aggregate(hh_obs[, exposure], by = list("exp_q" = hh_obs[, "exp_q"]), FUN = median)
    colnames(exp_mids) <- c("exp_q", "median")
    
    # Prepare output for balance-with-adjustment statistics
    out <- expand.grid(as.character(levels(hh_obs[, "exp_q"])), as.character(1:5),  confounders)
    colnames(out) <- c("exp_q", "gps_q", "confounder")
    out[, c("mean_diff", "se_diff", "n_obs", "n_obs_k", "n_obs_nonk")] <- NA
    out[, "exp_q"] <- as.ordered(out[, "exp_q"])
    out[, "gps_q"] <- as.ordered(out[, "gps_q"])
    
    # For each exposure quantile...
    for (k in levels(hh_obs[, "exp_q"])) {
      # compute GPS using exposure model, by fixing exposure to median within quantile
      hh_obs[, "gps_mid"] <- dnorm(x = rep(exp_mids[exp_mids[, "exp_q"] == k, "median"], nrow(hh_obs)), 
        mean = fitted.values(fit_exp), sd = summary(fit_exp)$sigma)
      
      # divide dataset into quintiles by GPS as estimated above
      gps_q <- quantile(hh_obs[, "gps_mid"], probs = seq(0, 1, by = 0.20))
      hh_obs[, "gps_q"] <- cut(hh_obs[, "gps_mid"], breaks = gps_q, labels = 1:(length(gps_q) - 1),
        ordered_result = TRUE, include.lowest = TRUE)
      
      # now, for each GPS quintile...
      for (l in levels(hh_obs[, "gps_q"])) {
        # select all data in this level
        x1 <- subset(hh_obs, gps_q == l)
        
        # split data by whether the exposure quantile == k or not
        x2 <- subset(x1, exp_q == k)
        x3 <- subset(x1, exp_q != k)
        
        # for each confounder...
        for (j in confounders) {
          # find position in output vector
          x4 <- which(out[, "exp_q"] == k & out[, "gps_q"] == l & out[, "confounder"] == j)
          
          # compute difference of means and its standard error
          out[x4, "mean_diff"] <- mean(x2[, j]) - mean(x3[, j])
          out[x4, "se_diff"] <- sqrt( (sd(x2[, j])^2 / nrow(x2)) +  (sd(x3[, j])^2 / nrow(x3))  )
        }  
        
        # capture number of observations within l level
        out[which(out[, "exp_q"] == k & out[, "gps_q"] == l), "n_obs"] <- nrow(x1)
        out[which(out[, "exp_q"] == k & out[, "gps_q"] == l), "n_obs_k"] <- nrow(x2)
        out[which(out[, "exp_q"] == k & out[, "gps_q"] == l), "n_obs_nonk"] <- nrow(x3)
      }
    }

    # Compute t-statistics for association of confounder with exposure, after GPS adjustment
      # remove NaN values from output (#####DUE TO SPARSITY: OK TO DO???? TO REVISIT)  
      out <- na.omit(out)
      
      # calculate weighted mean difference and SE of difference for each exposure quantile and confounder
      out[, "mean_diff_wts"] <- out[, "mean_diff"] * out[, "n_obs"]
      out[, "se_diff_wts"] <- out[, "se_diff"] * out[, "n_obs"]
      out <- aggregate(out[, c("n_obs", "n_obs_k", "n_obs_nonk", "mean_diff_wts", "se_diff_wts")], 
        by = out[, c("exp_q", "confounder")], FUN = sum)
      out[, "wt_mean_diff"] <- out[, "mean_diff_wts"] / out[, "n_obs"]
      out[, "wt_se_diff"] <- out[, "se_diff_wts"] / out[, "n_obs"]
      
      # compute t-statistic for each exposure quantile compared to other quantiles, and for each confounder
      out[, "t_value"] <- out[, "wt_mean_diff"] / out[, "wt_se_diff"]
      out[, "t_value"] <- round(out[, "t_value"], 2)

      # reshape wide and make other changes
      out_adj <- reshape2::dcast(out[, c("exp_q", "confounder", "t_value")], formula = confounder ~ exp_q)
      colnames(out_adj)[colnames(out_adj) != "confounder"] <- paste("adj", levels(hh_obs[, "exp_q"]), sep = "_")
      out_adj[, "confounder"] <- as.character(out_adj[, "confounder"])
      out_adj <- out_adj[order(out_adj[, "confounder"]), ]

  #...................................
  ## Check balance of confounders without adjustment by GPS
            
    # Prepare output for balance-without-adjustment statistics
    out <- expand.grid(as.character(levels(hh_obs[, "exp_q"])), confounders)
    colnames(out) <- c("exp_q", "confounder")
    out[, c("mean_diff", "se_diff", "n_obs", "n_obs_k", "n_obs_nonk")] <- NA
    out[, "exp_q"] <- as.ordered(out[, "exp_q"])

    # For each exposure quantile...
    for (k in levels(hh_obs[, "exp_q"])) {
      
      # split data by whether the exposure quantile == k or not
      x2 <- subset(hh_obs, exp_q == k)
      x3 <- subset(hh_obs, exp_q != k)
      
      # for each confounder...
      for (j in confounders) {
        # find position in output vector
        x4 <- which(out[, "exp_q"] == k & out[, "confounder"] == j)
        
        # compute difference of means and its standard error
        out[x4, "mean_diff"] <- mean(x2[, j]) - mean(x3[, j])
        out[x4, "se_diff"] <- sqrt( (sd(x2[, j])^2 / nrow(x2)) +  (sd(x3[, j])^2 / nrow(x3))  )
      }
      
      # capture number of observations within k quantile
      out[which(out[, "exp_q"] == k), "n_obs"] <- nrow(hh_obs)
        out[which(out[, "exp_q"] == k ), "n_obs_k"] <- nrow(x2)
        out[which(out[, "exp_q"] == k ), "n_obs_nonk"] <- nrow(x3)
    }  
      
    # Compute t-statistics for association of confounder with exposure, after GPS adjustment
      # remove NaN values from output (there shouldn't be any)  
      out <- na.omit(out)
      
      # compute t-statistic for each exposure quantile compared to other quantiles, and for each confounder
      out[, "t_value"] <- out[, "mean_diff"] / out[, "se_diff"]
      out[, "t_value"] <- round(out[, "t_value"], 2)

      # reshape wide and make other changes
      out_unadj <- reshape2::dcast(out[, c("exp_q", "confounder", "t_value")], formula = confounder ~ exp_q)
      colnames(out_unadj)[colnames(out_unadj) != "confounder"] <- paste("unadj", levels(hh_obs[, "exp_q"]), sep = "_")
      out_unadj[, "confounder"] <- as.character(out_unadj[, "confounder"])
      out_unadj <- out_unadj[order(out_unadj[, "confounder"]), ]

    
  #...................................
  ## Compare balance without and with adjustment by GPS
  
    # Merge outputs and save
    out <- merge(out_unadj, out_adj, by = "confounder")
    write.csv(out, file = paste(country, y_hat, "gps_balance.csv", sep = "_"), row.names = FALSE)
    
    # Plot change in t-values for each confounder, before and after adjustment
      # prepare long version of output
      out_unadj[, "method"] <- "unadjusted"
      colnames(out_unadj) <- gsub("unadj_", "", colnames(out_unadj))
      out_adj[, "method"] <- "GPS-adjusted"
      colnames(out_adj) <- gsub("adj_", "", colnames(out_adj))
      x1 <- rbind(out_unadj, out_adj)
      x1 <- reshape2::melt(data = x1, id.vars = c("confounder", "method"), variable.name = "quantile",
        value.name = "t_statistic")
      x1[, "method"] <- factor(x1[, "method"], levels = c("unadjusted", "GPS-adjusted"))
      x1[, "quantile"] <- as.numeric(x1[, "quantile"])
      x1[which(x1[, "quantile"] == 1), "quantile"] <- "1st quartile"
      x1[which(x1[, "quantile"] == 2), "quantile"] <- "2nd quartile"
      x1[which(x1[, "quantile"] == 3), "quantile"] <- "3rd quartile"
      x1[which(x1[, "quantile"] == 4), "quantile"] <- "4th quartile"
      
      # plot
      plot <- ggplot(x1, aes(y = t_statistic, x = method, group = confounder)) +
        scale_y_continuous(name = "t-statistic") +
        scale_x_discrete(expand = c(0.1,0.1)) +
        theme_bw() +
        theme(axis.title = element_text(colour = "grey20")) +
        geom_line(colour = palette_cb[7], alpha = 0.7) +
        annotate(geom = "rect", fill = palette_cb[3], alpha = 0.3, 
          xmin = -Inf, xmax = Inf, ymin = -1.96, ymax = 1.96) +
        facet_wrap(~quantile)
      
      plot
      
      # save plot and estimates
      ggsave(paste(country, y_hat, "gps_balance.png", sep = "_"), 
        dpi = "print", width = 20, height = 10, units = "cm")
      

#.........................................................................................
### Compute dose-response association between exposure and outcome
#.........................................................................................
    # based on Hirano & Imben (2004), https://www.google.com/url?sa=t&rct=j&q=&esrc=s&source=web&cd=&ved=2ahUKEwjozLahmOn1AhXGEcAKHQjOBQsQFnoECAMQAQ&url=https%3A%2F%2Fwww.math.mcgill.ca%2Fdstephens%2FPSMMA%2FArticles%2FHIrano-Imbens-2004.pdf&usg=AOvVaw1fo5Um-exkrcC-rwh2ozZb
    # and Austin (2017), https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5969262/#sim7615-bib-0003 

  #...................................
  ## Plot distributions of exposure (untransformed) and outcome
    
    # Exposure
    plot1 <- ggplot(data = hh_obs) + 
      geom_histogram(aes(x = exp(log_total_actors_rate) ), 
        colour = palette_cb[7], fill = palette_cb[7], alpha = 0.5 ) +
      theme_bw() +
      scale_y_continuous(name = "number of observations") + 
      scale_x_continuous(name = "number of projects per 100,000 people", expand = c(0,0)) +
      theme(axis.title = element_text(colour = "grey20"))
    plot1  
      
    # Outcome
    plot2 <- ggplot(data = hh_obs) + 
      geom_histogram(aes(x = n_died * 10000 / ptime), 
        colour = palette_cb[7], fill = palette_cb[7], alpha = 0.5 ) +
      theme_bw() +
      scale_y_continuous(name = "number of observations") + 
      scale_x_continuous(name = "deaths per 10,000 person-days", expand = c(0,0)) +
      theme(axis.title = element_text(colour = "grey20"))
    plot2  

    # Combine plots and save
    plot <- ggarrange(plot1, plot2, labels = NA, ncol = 1, nrow = 2)
    plot
    ggsave(paste(country, y_hat, "descriptive_plots.png", sep = "_"), 
        dpi = "print", width = 20, height = 13, units = "cm")
    
    
  #...................................
  ## Estimate effect of exposure - Hirano & Imbens method
    # Define formula for outcome as a function of exposure and GPS
    hh_obs[, paste(exposure, "_sq", sep = "")] <- hh_obs[, exposure] ^ 2
    hh_obs[, "gps_sq"] <- hh_obs[, "gps"] ^ 2
    formula_eff <- formula("n_died ~ log_total_actors_rate + log_total_actors_rate_sq +
      gps + gps_sq + log_total_actors_rate : gps")
    
    # Specify survey design
    design_svy <- svydesign(ids = ~Cluster, data = hh_obs)  
    
    # Fit GLM for outcome model
    fit_eff <- svyglm(formula_eff, design = design_svy, family = "poisson", offset = log(ptime))
    summary(fit_eff)

    
  #...................................
  ## Estimate dose-response function - Hirano & Imbens method
    # Specify levels of the exposure that function is evaluated at    
    exp_levels <- quantile(hh_obs[, exposure], probs = seq(0, 0.95, by = 0.05))

    # Prepare prediction dataset
    x1 <- hh_obs[, c("stratum", "Cluster", "ptime")]
    x1[, c(exposure, paste(exposure, "_sq", sep = ""), "gps", "gps_sq", "out_mean", "out_lci", "out_uci")] <- NA

    # Prepare output dose-response dataset
    out <- data.frame("exposure_level" = exp_levels, "outcome_mean" = NA, "outcome_lci" = NA, "outcome_uci" = NA)
    
    # For each exposure level...
    for (i in exp_levels) {
      
      # update prediction dataset
      x1[, exposure] <- i
      x1[, paste(exposure, "_sq", sep = "")] <- i^2
      
      # predict GPS at desired exposure level
      x1[, "gps"] <- dnorm(x = x1[, exposure], mean = fitted.values(fit_exp), sd = summary(fit_exp)$sigma)
      x1[, "gps_sq"] <- x1[, "gps"]^2
      
      # predict outcome at desired exposure level
      x1[, "out_mean"] <- predict(fit_eff, newdata = x1, type = "response")
      x2 <- as.data.frame(predict(fit_eff, newdata = x1, type = "link", se.fit = TRUE))
      x1[, "out_lci"] <- exp(x2[, 1] - 1.96 * x2[, 2])
      x1[, "out_uci"] <- exp(x2[, 1] + 1.96 * x2[, 2])
      
      # mean and 95%CI of outcome for this exposure level (as deaths per 10,000 person-days)
      out[out[, "exposure_level"] == i, "outcome_mean"] <- mean(x1[, "out_mean"]) * 10000
      out[out[, "exposure_level"] == i, "outcome_lci"] <- mean(x1[, "out_lci"]) * 10000
      out[out[, "exposure_level"] == i, "outcome_uci"] <- mean(x1[, "out_uci"]) * 10000
      
    }
 
    # Plot dose-response function
    plot_hi <- ggplot(out, aes(x = exp(exposure_level))) +
      geom_line(aes(y = outcome_mean), colour = palette_cb[7]) +
      geom_ribbon(aes(ymin = outcome_lci, ymax = outcome_uci), alpha = 0.3, fill = palette_cb[7]) +
      theme_bw() +
      theme(axis.title = element_text(colour = "grey20")) +
      scale_y_continuous(name = "deaths per 10,000 person-days") +
      scale_x_continuous(name = "projects per 100,000 population", limits = c(0, NA),
        breaks = seq(0, 100, 10), expand = c(0,0))
    plot_hi
    
    # Save plot and estimates
    ggsave(paste(country, y_hat, "dose_response_hi.png", sep = "_"), 
      dpi = "print", width = 10, height = 10, units = "cm")
    write.csv(out, paste(country, y_hat, "dose_response_hi.csv", sep = "_"), row.names = FALSE)    
    
    
  #...................................
  ## Estimate effect of exposure - propensity weights method
    # Define formula for outcome as a function of exposure
    formula_eff <- formula("n_died ~ log_total_actors_rate")

    # Specify survey design
    design_svy <- svydesign(ids = ~Cluster, weights = ~gps_wt, data = hh_obs)  
    
    # Fit GLM for outcome model
    fit_eff <- svyglm(formula_eff, design = design_svy, family = "poisson", offset = log(ptime))
    summary(fit_eff)

  #...................................
  ## Estimate dose-response function - propensity weights method
    # Specify levels of the exposure that function is evaluated at    
    exp_levels <- quantile(hh_obs[, exposure], probs = seq(0, 0.95, by = 0.05))

    # Prepare prediction dataset
    x1 <- hh_obs[, c("Cluster", "ptime")]
    x1[, c(exposure, "out_mean", "out_lci", "out_uci")] <- NA

    # Prepare output dose-response dataset
    out <- data.frame("exposure_level" = exp_levels, "outcome_mean" = NA, "outcome_lci" = NA, "outcome_uci" = NA)
    
    # For each exposure level...
    for (i in exp_levels) {
      
      # update prediction dataset
      x1[, exposure] <- i
      
      # predict outcome at desired exposure level
      x1[, "out_mean"] <- predict(fit_eff, newdata = x1, type = "response")
      x2 <- as.data.frame(predict(fit_eff, newdata = x1, type = "link", se.fit = TRUE))
      x1[, "out_lci"] <- exp(x2[, 1] - 1.96 * x2[, 2])
      x1[, "out_uci"] <- exp(x2[, 1] + 1.96 * x2[, 2])
      
      # mean and 95%CI of outcome for this exposure level (as deaths per 10,000 person-days)
      out[out[, "exposure_level"] == i, "outcome_mean"] <- mean(x1[, "out_mean"]) * 10000
      out[out[, "exposure_level"] == i, "outcome_lci"] <- mean(x1[, "out_lci"]) * 10000
      out[out[, "exposure_level"] == i, "outcome_uci"] <- mean(x1[, "out_uci"]) * 10000
      
    }
 
    # Plot dose-response function
    plot_wt <- ggplot(out, aes(x = exp(exposure_level))) +
      geom_line(aes(y = outcome_mean), colour = palette_cb[7]) +
      geom_ribbon(aes(ymin = outcome_lci, ymax = outcome_uci), alpha = 0.3, fill = palette_cb[7]) +
      theme_bw() +
      theme(axis.title = element_text(colour = "grey20")) +
      scale_y_continuous(name = "deaths per 10,000 person-days") +
      scale_x_continuous(name = "projects per 100,000 population", limits = c(0, NA),
        breaks = seq(0, 100, 10), expand = c(0,0))
    plot_wt
    
    # Save plot and estimates
    ggsave(paste(country, y_hat, "dose_response_wt.png", sep = "_"), 
      dpi = "print", width = 10, height = 10, units = "cm")
    write.csv(out, paste(country, y_hat, "dose_response_wt.csv", sep = "_"), row.names = FALSE)    
    
    
  #...................................
  ## Produce combined plots
  plot <- ggarrange(plot_hi + 
      scale_y_continuous("deaths per 10,000 person-days", limits = c(0, 1.5), breaks = seq(0, 1.50, 0.25)) +
      theme(plot.margin = margin(t = 1, r = 0, b = 1.5, l = 0) ), 
    plot_wt + 
      scale_y_continuous("deaths per 10,000 person-days", limits = c(0, 1.5), breaks = seq(0, 1.50, 0.25)) +
      theme(axis.title.y = element_blank()) +
      theme(plot.margin = margin(t = 1, r = 1, b = 1.5, l = 0) ), 
    ncol = 2, align = "hv", 
    font.label = list(size = 10, colour = "grey20", face = "bold"), hjust = -0.7, vjust = 4,
    labels = c("GPS adjustment method", "GPS-derived weighting method"))  
  plot
  ggsave(paste(country, y_hat, "dose_response_combi.png", sep = "_"), 
      dpi = "print", width = 25, height = 10, units = "cm")

    
#.........................................................................................
### ENDS
#.........................................................................................
  

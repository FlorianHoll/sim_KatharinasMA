# ------------------------------------------ #
# ----- Power Simulation Katharinas MA ----  #
# ------------------------------------------ #

# === Design der Studie === #

# --- UVs: --- #
# PT (historical vs. psychological) 
# Rotation (40° vs. 160°)
# Anchoring (high vs. low; nur fuer psychological PT)

# --- AVs --- # 
# Visuospatial PT: Reaktionszeiten auf Aufgabe (Erwartung: Hoeher bei größerer Rotation (160°))
# Historical PT: Kompetenz in Historie-Aufgaben (d.h. Parceling durch Summe oder mit WLEs) (Erwartung: Hoehere Kompetenz bei 160°)
# Psychological PT: 
#   Anchoring Effect (Unterschied zu wahrem Wert); Erwartung: Hoeherer Anchoring Effect bei 160°
#   Anchoring Differences (Unterschied zu Schaetzung der Person): kleinere Differenz bei 160°

# === Hier untersucht === #
# Gegeben einer bestimmten Anzahl von Items, Versuchspersonen und Effektgröße, wie groß ist die Wahrscheinlichkeit, dass die
# Schaetzmethode den tatsaechlichen Effekt findet?


library(parallel)  
library(ggplot2)

num_cores = detectCores()
cl <- makeCluster(num_cores)

all_global <- function() {
  lss <- ls(envir = parent.frame())
  for (i in lss) {
    assign(i, get(i, envir = parent.frame()), envir = .GlobalEnv)
  }
}

# push everything to cluster
clusterEvalQ(cl, {
  # ------------------------------------------ #
  # ----- Packages and smaller functions ----- #
  # ------------------------------------------ #
  library(mirt)
  library(lavaan)
  
  z_norm = function(x, y=x) (x - mean(y, na.rm=T)) / sd(y, na.rm=T)
  
  conf_intval = function(x, se) c(x - 1.96 * se, x + 1.96 * se)
  
  is_in_conf_intval = function(x, conf_intval) {
    if (is.na(x) | is.na(conf_intval[1]) | is.na(conf_intval[2])) { 
      return(NA)
    } else {
      if (x >= conf_intval[1] && x <= conf_intval[2]) return(TRUE) else return(FALSE)
    }
  }
  
  meansd = function(x) c(mean(x), sd(x))
  
  return_model_or_NA = function(model_to_fit, df) {
    m = tryCatch(
      sem(model_to_fit, df, estimator="MLR", missing="fiml"),
      error = function(cond) return(NA)
    )
    return(m)
  }
  
  extract_estimates = function(model) {
    if ( !is.na(model)) {
      std_est = subset(standardizedSolution(model, type="std.lv"), op == "~")
      return(as.numeric(std_est[4:5]))
    } else {
      return(NA)
    }
  }
    
  converged = function(model) inspect(model, "converged")
  
  # ------------------------------------------ #
  # --------- IRT MODEL (1PL or 2PL) --------- #
  # ------------------------------------------ #
  irt_model = function(a=1, d, t) 1 / (1 + exp(a * (d - t)))
  
  
  # ------------------------------------------ #
  # -- Function to insert NA at right place -- #
  # ------------------------------------------ #
  insert = function(x, insert_pos, value) {
    
    new = numeric(0)
    
    l = length(x) + length(insert_pos)
    
    for(i in 1:(length(insert_pos)+1)) {
      
      if (i <= length(insert_pos)) {
        if (i == 1) new = append(new, x[1:(insert_pos[i]-1)])
        else new = append(new, x[(insert_pos[i-1]-(i-2)) : (insert_pos[i]-i)])
        new = append(new, value)
      }
      
      else {
        new = append(new, x[(insert_pos[i-1]-(i-2)):length(x)])
      }
    }
    return(new[1:l])
  }
  
  
  # ------------------------------------------ #
  # ------ function to split the dataset ----- #
  # ------------------------------------------ #
  split_test = function(df2, n_splits, intercepts) {
    
    criterion = sqrt(n_splits)/10 
    criterion = max(c(0.002, criterion-.15))
    
    n_items_per_split = ncol(df2) %/% n_splits
    i=0; cat("Iteration \tMSE \tCriterion \n ---------------------------------- \n")
    while(1) {
      shuffled_items = sample(ncol(df2))
      item_split = split(shuffled_items, 1:n_splits) #split items
      
      # if there are more items in the last split, add them to previous
      if (length(item_split) != n_splits) {
        item_split[[length(item_split) - 1]] = append(item_split[[length(item_split) - 1]], 
                                                      item_split[[length(item_split)]])
        item_split[[length(item_split)]] = NULL
      }
      
      split_means = sapply(item_split, function(x) mean(intercepts[x]))
      
      MSE = 1/(length(split_means) - 1) * sum((split_means[-1] - split_means[1]) ^ 2 )
      
      i = i+1
      if (i %% 1000 == 0) criterion = criterion + .02
      
      cat("\r", i, "\t\t", round(MSE, 3), "\t", round(criterion, 3))
      if (MSE < criterion) break
      
    }
    
    return(item_split)
    
  }
  
  
  # ------------------------------------------ #
  # ------ function for simulating data ------ #
  # ------------------------------------------ #
  create_dataset = function(N_participants, N_trials_per_participant, n_splits, hist_effect, ICC) {
    
    N_data_points = N_participants * N_trials_per_participant
    
    # Simulate thetas -- between level (i.e. between participants)
    theta_betw = rnorm(N_participants, 0, sd=sqrt(ICC * 1)) # N ~ (0, ICC*1)
    
    # Simulate thetas -- within level (i.e. within participants)
    theta_within_40 = rnorm(N_participants, 0, sd=sqrt((1-ICC) * 1))
    theta_within_160 = rnorm(N_participants, 0 + hist_effect, sd=sqrt((1-ICC) * 1))
    
    # Add between and within level thetas to form final thetas
    thetas = c(theta_betw + theta_within_40, theta_betw + theta_within_160)
    
    # item balancing
    items = list()
    for(i in 1:N_trials_per_participant) {
      items[[i]] = list(c(i:N_trials_per_participant, 1:i)[1:(N_trials_per_participant/2)], 
                        c(i:N_trials_per_participant, 1:i)[(N_trials_per_participant/2 + 1):(N_trials_per_participant)])
    }
    items = rep(items, 50)
    
    # store everything in data frame
    df = data.frame(ID = rep(1 : N_participants, 2),
                    rot = rep(c(40, 160), each=N_participants),
                    theta = thetas)
    
    # sort by ID
    df = df[order(df$ID),]
    
    # quick check: Did everything work? --> Yep! 
    with(df, aggregate(theta ~ rot, FUN="mean"))
    
    # simulate IRT items
    item_list = paste0("I", 1:N_trials_per_participant)
    intercepts = round(rnorm(N_trials_per_participant, 0, 1), 3)
    if (model_type !="Rasch") alphas = round(runif(N_trials_per_participant, .7, 1.2), 3) else alphas = rep(1, N_trials_per_participant)
    item_params = data.frame(alphas=alphas, difficulty=intercepts)
    
    # simulate responses to items
    for(i in 1:length(item_list)) {
      prob <- irt_model(alphas[i], intercepts[i], df$theta)
      df[item_list[i]] <- rbinom(N_participants*2, 1, prob)
    }
   
    # insert NAs to account for design of experiment (participants only see one version of each item)
    df2 = df[-c(1:3)] # data frame only with columns of interest, i.e. items
    
    for(i in df$ID) {
      for(j in 1:2) {
        no_nas = items[[i]][[j]]
        nas = c(1:N_trials_per_participant)[-no_nas]
        df2[i*2 + j-2, nas] = NA
      }
    }
    
    df[,4:11] = df2
    
    # split data into N parts 
    split_list = split_test(df2, n_splits, intercepts)
    
    dataset = list(df, 
                   list(alphas, intercepts),
                   split_list)
    
    names(dataset) = c("df", "IRT_parameters", "split_list")
    names(dataset$IRT_parameters) = c("alphas", "intercepts")
    
    return(dataset)
    
  }
  
  
  # ------------------------------------------ #
  # --- function for estimating IRT model ---- #
  # ------------------------------------------ #
  estimate_IRT_model = function(dataset_list) {
    
    # --- estimated models and indicators --- #
    #                                         #
    # [1] IRT model with all items            #
    # [2] IRT models with parceling items     #
    # [3] Sum Score with all items            #
    # [4] Sum Score with parceling items      #
    #                                         #
    # --------------------------------------- #
    
    irt_df = dataset_list$df[,4:11]
    
    # --- [1] estimate IRT model with all items --- #
    model_string = paste0("F1 = 1 - ", ncol(irt_df))
    
    model = tryCatch(
      mirt::mirt(irt_df, 
                 model=model_string,
                 itemtype=model_type, 
                 technical = list(removeEmptyRows = TRUE)),
      error = function(cond) return(NA)
    )
    
    # estimate WLEs for single IRT model
    if (!is.na(model)) {
      WLEs_all <- fscores(model, method="WLE", full.scores=T)
      WLEs_rel <- fscores(model, method="WLE", full.scores=T, returnER=T)
    } else WLEs_all = WLEs_rel = NA
    
    # into dataframe
    dataset_list$df["WLEs_all"] = as.numeric(WLEs_all)
    
    
    # --- [2] estimate IRT models with WLEs for each split --- #
    NA_pos = lapply(dataset_list$split, 
                 function(x) {
                    as.numeric(which(rowSums(is.na(irt_df[, x])) == length(x)))
                 }
    )
    
    # estimate IRT model for each split
    models <- lapply(dataset_list$split_list, 
                     function(x) {
                       model_string = paste0("F1 = 1 - ", length(x))
                       m = tryCatch(
                         mirt::mirt(irt_df[,x], 
                                    model=model_string,
                                    itemtype=model_type, 
                                    technical = list(removeEmptyRows = TRUE)),
                         error = function(cond) return(NA)
                       )
                       return(m)
                     }
    )
      
    # estimate WLEs
    WLEs <- lapply(models,
                   function(x) {
                     if (!is.na(x)) {
                       fscores(x, method="WLE", 
                               full.scores = T)
                     } else NA
                   }
    )
    
    WLE_names = paste0("WLE", 1:length(WLEs))
    
    if (length(which(sapply(WLEs, length) == N_participants * n_splits)) != n_splits) {
      for(i in 1:length(WLEs)) {
        if (!is.na(WLEs[[i]])) {
          WLEs[[i]] = insert(as.numeric(WLEs[[i]]), NA_pos[[i]], NA)
        }
      }
    }
      
    names(WLEs) = WLE_names
    
    # into data frame
    dataset_list$df[WLE_names] = sapply(WLEs, function(x) x)
    
    
    # --- [3] sum score with all items --- #
    sum_score_all = as.numeric(rowSums(irt_df, na.rm = T))
    dataset_list$df["sum_score_all"] = sum_score_all
    
    # -- [4] sum scores with parceling items -- # 
    dataset_list$df[paste0("sum_score", 1:n_splits)] = sapply(dataset_list$split_list, 
                                                              function(x) {
                                                                rowSums(irt_df[,x], na.rm=T)
                                                                }
    )
    
    
    # --- output --- #
    dataset_list[[4]] = list(list(WLEs_all, WLEs_rel), WLEs)
    
    names(dataset_list[[4]]) = c("WLE_whole_estimates", "WLE_split_estimates")
    names(dataset_list[[4]]$WLE_whole_estimates) = c("WLEs", "Reliability")
    dataset_list[[5]] = NA_pos
    names(dataset_list)[4:5] = c("WLEs", "NA_pos")
    
    return(dataset_list)
    
  } 
  
  
  # ------------------------------------------ #
  # --- Linear regression and latent models -- #
  # ------------------------------------------ #
  estimate_compare_effects = function(dataset_list) {
    
    # ----------- estimated models ---------- #
    #                                         #
    # [0] linear regression                   #
    # [1] SEM with WLEs                       #
    # [2] SEM with Sum Scores                 # 
    # [3] SEM with single indicator WLEs      #
    # [4] SEM with manifest WLEs              #
    # [5] SEM with manifest Sum Scores        #
    #                                         #
    # --------------------------------------- #
    
    
    # --- data structuring for analyses --- #
    d = dataset_list$df
    
    WLE_norm_names = grep("WLE", names(d), value=T)
    sum_scores_norm_names = grep("sum_score", names(d), value=T)

    WLE_parceling_names = paste0("WLE", 1:n_splits)
    sum_scores_parceling_names = paste0("sum_score", 1:n_splits)
    
    # linear regression 
    d$rot = as.factor(d$rot)
    
    # latent WLE model
    d[paste0(WLE_norm_names, "_norm")] = sapply(d[WLE_norm_names], z_norm)
    d[paste0(sum_scores_norm_names, "_norm")] = sapply(d[sum_scores_norm_names], z_norm)
    
    
    # --- [0] Linear regression model to estimate true effect --- #
    lin_model_norm = lm(scale(theta) ~ rot, d)
    
    # extract effect estimate 
    true_effect = as.numeric(coef(lin_model_norm)[2])
    true_effect_se = summary(lin_model_norm)$coefficients[2,2]
    
    
    # --- [1] Latent model with WLEs as indicators --- #
    lat_model_WLEs_s = paste(
          
          # latent variable 
          "Hist_Comp =~ ", paste0("1*", WLE_parceling_names, collapse = " + "), "\n",
          
          # intercepts
          paste0(WLE_parceling_names, collapse=" + "), "~1", "\n",
          
          # structural model 
          "Hist_Comp ~ rot"
    )
    
    lat_model_WLEs = return_model_or_NA(lat_model_WLEs_s, d)
    

    # --- [2] Latent model with sum scores --- #
    lat_model_sum_scores_s = paste(
      
      # latent variable 
      "Hist_Comp =~ ", paste0("1*", sum_scores_parceling_names, collapse = " + "), "\n",
      
      # intercepts
      paste0(sum_scores_parceling_names, collapse=" + "), "~1", "\n",
      
      # structural model 
      "Hist_Comp ~ rot"
    )
    
    lat_model_sum_scores = return_model_or_NA(lat_model_sum_scores_s, d)
    
    
    # --- [3] Latent model with WLE as single indicator --- #
    WLE_Rel = dataset_list$WLEs$WLE_whole_estimates$Reliability
    lat_model_single_indicator_s = paste(
      
      # latent variable 
      "Hist_Comp =~ ", paste0(sqrt(WLE_Rel), "*", "WLEs_all"), "\n",
      
      # residual error 
      "WLEs_all ~~ ", (1-WLE_Rel) * var(as.numeric(dataset_list$WLEs$WLE_whole_estimates$WLEs)), "*WLEs_all \n",
      
      # structural model 
      "Hist_Comp ~ rot"
    )
    
    lat_model_single_indicator = return_model_or_NA(lat_model_single_indicator_s, d)
    
    
    # --- [4] Latent model with manifest WLEs --- #
    lat_model_manifest_WLEs_s = paste(
      
      # structural model 
      "WLEs_all ~ rot"
    )
    
    lat_model_manifest_WLEs = return_model_or_NA(lat_model_manifest_WLEs_s, d)
    
    
    # --- [5] Latent model with manifest Sum Scores --- #
    lat_model_manifest_Sum_scores_s = paste(
      
      # structural model 
      "sum_score_all ~ rot"
    )
    
    lat_model_manifest_Sum_scores = return_model_or_NA(lat_model_manifest_Sum_scores_s, d)
    
    
    # --- compare to true effect --- #
    latent_models = list(lat_model_WLEs, 
                         lat_model_sum_scores, 
                         lat_model_single_indicator, 
                         lat_model_manifest_WLEs, 
                         lat_model_manifest_Sum_scores)
    
    names_lat_models = c("parceling_WLEs", "parceling_Sum", "SI_WLE", "Manifest_WLE", "Manifest_Sum")
    
    names(latent_models) = names_lat_models
    
    latent_estimates = lapply(latent_models, function(x) extract_estimates(x))
    
    diff = sapply(latent_estimates, function(x) x[1] - as.numeric(true_effect))
    
    in_ci = sapply(latent_models, 
                   function(x) {
                     if (is.na(x)) {
                       NA
                     } else {
                       std_est = subset(standardizedSolution(x, type="std.lv"), op == "~")
                       is_in_conf_intval(true_effect, c(std_est[8], std_est[9])) 
                     }
                   }
    )
    
    simres = data.frame(model = c("LIN_REG", names_lat_models),
                         estimate = c(true_effect, as.numeric(sapply(latent_estimates, function(x) x[1]))),
                         se = c(true_effect_se, as.numeric(sapply(latent_estimates, function(x) x[2]))),
                         in_ci = c(NA, in_ci), 
                         diff = c(0, diff),
                         abs_diff = c(0, abs(diff)))
    
    return(simres)
  }
  
  
})

# ------------------------------------------ #
# ---------- actual simulation  ------------ #
# ------------------------------------------ #
start_simulation = function(N_simulations = 10,
                            N_participants = 100, 
                            N_trials_per_participant = 8, 
                            model_type = "2PL", 
                            n_splits = 2, 
                            hist_effect = .5, 
                            ICC = .1,
                            plot = FALSE
  ) {
  
  # push variables to global environment
  all_global()
  
  # export variables to clusters
  export_to_multiple_clusters = c("N_simulations", 
                                  "N_participants", 
                                  "N_trials_per_participant",
                                  "model_type", 
                                  "n_splits", 
                                  "hist_effect", 
                                  "ICC")
  clusterExport(cl, export_to_multiple_clusters)
  
  # --- SE function --- #
  se = function(x) return ( sd(x, na.rm=T) / sqrt(length(x)) )
  
  # --- create data --- #
  datasets <- parLapply(cl, 1:N_simulations, function(x) {
    create_dataset(N_participants, N_trials_per_participant, n_splits, hist_effect, ICC)})
  
  # --- estimate IRT models --- # 
  datasets <- parLapply(cl, datasets, function(x) estimate_IRT_model(x))
  
  # --- compare effect estimates --- #
  results <- parLapply(cl, datasets, function(x) estimate_compare_effects(x))
  
  # --- convert to one result data frame --- # 
  sim_results = results[[1]]
  if (length(results) > 1) {
    for(i in 2:length(results)) {
      sim_results = rbind(sim_results, results[[i]])
    }
    sim_results$sim_nr = rep(1:N_simulations, each=6)
  }
  
  # --- output --- #
  names_lat_models = c("parceling_WLEs", "parceling_Sum", "SI_WLE", "Manifest_WLE", "Manifest_Sum")
  
  amount_converged = sapply(names_lat_models, 
                            function(x) {
                              N_simulations - length(which(is.na(subset(sim_results, model==x)$estimate)))
                            })
  
  differences = sapply(names_lat_models, 
                       function(x) {
                         mean(with(sim_results, diff[model==x]), na.rm=T)
                       })
  
  differences_se = sapply(names_lat_models, 
                          function(x) {
                            se(with(sim_results, diff[model==x]))
                          })
  
  estimate_in_ci = sapply(names_lat_models,
                          function(x) {
                            sum(as.numeric(as.logical(with(sim_results, in_ci[model==x]))), na.rm = T) 
                          })
  
  names_models_long = c("Parceling-WLE estimates", 
                        "Parceling Sum Scores", 
                        "Single Indicator WLEs", 
                        "Manifest WLE estimates", 
                        "Manifest Sum Scores ")
  
  output_string = paste("\n--------------------------------------------------\n\nCONVERGENCE / ESTIMATABILITY:",
                        paste0("\n  ", names_models_long, ": \t", 
                               round(amount_converged / N_simulations * 100, 3), 
                              "% of the models converged. (", 
                              amount_converged, " out of ", N_simulations, ")", collapse=""),
                        
                        "\n\nMEAN DIFFERENCES FROM TRUE EFFECT (ONLY THOSE THAT CONVERGED):",
                        paste0("\n  ", names_models_long, ": \t", round(differences, 3), 
                               " (SE: ", round(differences_se, 3), ")", collapse=""),
                        
                        "\n\nPROPORTION OF ESTIMATES IN TRUE EFFECT CONFIDENCE INTERVAL (ONLY THOSE THAT CONVERGED):", 
                        paste0("\n  ", names_models_long, ": \t", round(estimate_in_ci / amount_converged * 100, 3), 
                              "% (", estimate_in_ci, " out of ", amount_converged, ")", collapse=""), 
                        
                        "\n\n--------------------------------------------------\n")
  
  cat(output_string)
  
  return(sim_results)
  
}


# ------------------------------------------ #
# --- Functions for finished simulation ---- #
# ------------------------------------------ #

plot_diff_results = function(sim_results) {
  
  p = 
    ggplot(subset(sim_results, model!="LIN_REG"), aes(x=model, y=diff, fill=model)) +
    geom_violin(alpha=.5) + 
    scale_fill_brewer(palette = "Dark2") + 
    geom_boxplot(width=.2, fill="white") +
    theme_minimal() + 
    annotate("text", x=1.5, y=.05, label="True Effect", col="red") +
    geom_hline(yintercept=0, col="red", linetype="dashed") + 
    labs(title= "Differences from True Effect",
         y="Difference", x = "Model") + 
    ylim(-.6, .6)
  
  print(p)
  
}

plot_SEsize = function(sim_results) {
  
  p = 
    ggplot(subset(sim_results, model!="LIN_REG"), aes(x=model, y=se, fill=model)) +
    geom_violin(alpha=.5) + 
    scale_fill_brewer(palette = "Dark2") + 
    geom_boxplot(width=.2, fill="white") +
    theme_minimal() + 
    labs(title= "Size of Standard Error",
         y="Size of SE", x = "Model") 
  print(p)
}

plot_each_result = function(sim_results, N_plots_per_page) {
  
  plot_seq = seq(1, nrow(sim_results), N_plots_per_page*6)
  
  for(i in 1:length(plot_seq)) {
    df = sim_results[plot_seq[i] : (plot_seq[i] + N_plots_per_page * 6 - 1),]
    pl = 
      ggplot(df, aes(x=model, y=estimate, colour=model, ymin=estimate[1]-1.96*se[1], ymax=estimate+1.96*se[1])) + 
      geom_errorbar(aes(ymin=estimate-se, ymax=estimate+se), width=.2, size=1) + 
      geom_line() + 
      geom_point(size=2) + 
      facet_wrap(~sim_nr) +
      ylim(-.5, 1.5) +
      geom_text(data = subset(df, model=="LIN_REG"), 
                aes(y = estimate, label = "True Effect")) +
      geom_hline(data = subset(df, model=="LIN_REG"), 
                 aes(yintercept=estimate), color="black", linetype="dashed", size=.4)
      
    print(pl)
    
    Sys.sleep(15)
  }
  
}

output_results_again = function(sim_results) {
  
  se = function(x) return ( sd(x, na.rm=T) / sqrt(length(x)) )
  
  # --- output --- #
  names_lat_models = c("parceling_WLEs", "parceling_Sum", "SI_WLE", "Manifest_WLE", "Manifest_Sum")
  
  amount_converged = sapply(names_lat_models, 
                            function(x) {
                              N_simulations - length(which(is.na(subset(sim_results, model==x)$estimate)))
                            })
  
  differences = sapply(names_lat_models, 
                       function(x) {
                         mean(with(sim_results, diff[model==x]), na.rm=T)
                       })
  
  differences_se = sapply(names_lat_models, 
                          function(x) {
                            se(with(sim_results, diff[model==x]))
                          })
  
  estimate_in_ci = sapply(names_lat_models,
                          function(x) {
                            sum(as.numeric(as.logical(with(sim_results, in_ci[model==x]))), na.rm = T) 
                          })
  
  names_models_long = c("Parceling-WLE estimates", 
                        "Parceling Sum Scores", 
                        "Single Indicator WLEs", 
                        "Manifest WLE estimates", 
                        "Manifest Sum Scores ")
  
  output_string = paste("\n--------------------------------------------------\n\nCONVERGENCE / ESTIMATABILITY:",
                        paste0("\n  ", names_models_long, ": \t", round(amount_converged / N_simulations * 100, 3), 
                               "% of the models converged. (", amount_converged, " out of ", N_simulations, ")", collapse=""),
                        
                        "\n\nMEAN DIFFERENCES FROM TRUE EFFECT (ONLY THOSE THAT CONVERGED):",
                        paste0("\n  ", names_models_long, ": \t", round(differences, 3), " (SE: ", round(differences_se, 3), ")", collapse=""),
                        
                        "\n\nPROPORTION OF ESTIMATES IN TRUE EFFECT CONFIDENCE INTERVAL (ONLY THOSE THAT CONVERGED):", 
                        paste0("\n  ", names_models_long, ": \t", round(estimate_in_ci / amount_converged * 100, 3), 
                               "% (", estimate_in_ci, " out of ", amount_converged, ")", collapse=""), 
                        
                        "\n\n--------------------------------------------------\n")
  
  cat(output_string)
  
}

write_results_to_csv = function(sim_results, out_path, sim_nr) {
  write.csv2(sim_results, paste0(out_path, "/simulation_", sim_nr, ".csv"), row.names = F)
}

out_path = "C:/Users/Holl/Documents/Promotion/Lehre/WS_20_21_MA/power_sim/finished_csvs"

sim_results1 = start_simulation(N_simulations = 500)
sim_results2 = start_simulation(N_simulations = 500, model_type = "Rasch")



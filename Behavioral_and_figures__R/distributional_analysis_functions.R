# Distributional-analysis helpers for the MULTI behavioral data.
# Functions compute subject-level conditional-accuracy and RT delta plots,
# including previous-congruency contrasts used for delta-delta summaries.
# Author: Davide Gheza

library("lsr") # for quantileCut function


# Summarize RT and accuracy within participant-specific RT quantiles for one
# congruency value of one non-cued dimension.
compute_quantiles <- function(data, congruency_column, n_bins, congruency_value) {
  # Filter by congruency
  filtered_data <- filter(data, !!sym(congruency_column) == congruency_value)
  
  # Compute quantiles for each subject and non-cued dimension congruency
  sbjs_quantile <- filtered_data %>%
    group_by(subid) %>%                        ## quantile is defined within sbj
    mutate(RT_quantile = as.factor(as.numeric(quantileCut(rt, n_bins))))
  
  # Left join quantiles and summarize data
  result <- left_join(filtered_data, sbjs_quantile, by = c(colnames(filtered_data))) %>%
    group_by(RT_quantile, subid) %>%
    summarise(
      accuracy = mean(correct) * 100,
      count = n(),
      mean_rt = mean(rt),
      congruency = congruency_value
    )
  
  return(result)
}


# Compute congruent-minus-incongruent accuracy and RT delta functions for one
# non-cued dimension and each participant.
process_dimension <- function(data, congruency_column, dimension_label, n_bins) {
  # Compute congruent and incongruent data
  congruent_data <- compute_quantiles(data, congruency_column, n_bins, congruency_value = 1)
  incongruent_data <- compute_quantiles(data, congruency_column, n_bins, congruency_value = 0)
  
  # Combine congruent and incongruent data
  combined_data <- bind_rows(congruent_data, incongruent_data)
  
  # Compute average RT for each quantile (and non-cued dimension). Note: average across congruency levels and across subjects! 
  avg_quantile_RT <- combined_data %>%
    group_by(RT_quantile) %>%
    summarise(avg_quantile_RT = weighted.mean(mean_rt, count))
  
  # Merge average quantile RT with combined data
  combined_data <- left_join(combined_data, avg_quantile_RT, by = "RT_quantile")
  
  # Compute delta 
  delta_data <- combined_data %>%
    group_by(RT_quantile, subid) %>%
    reframe(
      avg_quantile_RT = avg_quantile_RT[1],
      acc_con = accuracy[congruency == 1],
      acc_inc = accuracy[congruency == 0],
      delta_acc = accuracy[congruency == 0] - accuracy[congruency == 1],
      delta_rt = mean_rt[congruency==0]-mean_rt[congruency==1],
      non_cued_dimension = dimension_label
    )
  
  return(delta_data)
}


# Apply `process_dimension` to all four non-cued dimensions.
process_all_dimensions <- function(data, n_bins, Labels) {
  # Process each dimension
  delta_data_d1 <- process_dimension(data, paste("congruency",tolower(Labels[1]),sep="_"), Labels[1], n_bins) # e.g., process_dimensions(data, "congruency_color", "Color", n_bins)
  delta_data_d2 <- process_dimension(data, paste("congruency",tolower(Labels[2]),sep="_"), Labels[2], n_bins)
  delta_data_d3 <- process_dimension(data, paste("congruency",tolower(Labels[3]),sep="_"), Labels[3], n_bins)
  delta_data_d4 <- process_dimension(data, paste("congruency",tolower(Labels[4]),sep="_"), Labels[4], n_bins)
  
  # Combine all results
  delta_data_eachNC <- bind_rows(delta_data_d1, delta_data_d2, delta_data_d3, delta_data_d4) %>%
    mutate(non_cued_dimension = factor(non_cued_dimension, levels = Labels, labels = Labels))
  
  return(delta_data_eachNC)
}









# Summarize all four parametric congruency levels and compute participant-level
# contrasts between the fully incongruent condition and levels 1--3.
process_sumcongruency <- function(data, n_bins) {
  # Compute summaries for each parametric congruency level.
  sum_con_0 <- compute_quantiles(data, congruency_column="sum_congruency", n_bins, congruency_value = 0)
  sum_con_1 <- compute_quantiles(data, congruency_column="sum_congruency", n_bins, congruency_value = 1)
  sum_con_2 <- compute_quantiles(data, congruency_column="sum_congruency", n_bins, congruency_value = 2)
  sum_con_3 <- compute_quantiles(data, congruency_column="sum_congruency", n_bins, congruency_value = 3)
  
  # Combine all congruency levels.
  combined_data <- bind_rows(sum_con_0, sum_con_1, sum_con_2, sum_con_3)
  
  # Compute average RT for each quantile. Note: average across congruency levels and across subjects! 
  avg_quantile_RT <- combined_data %>%
    group_by(RT_quantile) %>%
    summarise(avg_quantile_RT = weighted.mean(mean_rt, count))
  
  # Compute average RT by congruency level and quantile. Note: average across subjects only! (needed for Conditional accuracy function plot) 
  avg_cng_quantile_RT <- combined_data %>%
    group_by(RT_quantile, congruency) %>%
    summarise(avg_cng_quantile_RT = weighted.mean(mean_rt, count))
  
  # Merge average quantile RT with combined data
  combined_data <- left_join(combined_data, avg_quantile_RT, by = "RT_quantile")
  combined_data <- left_join(combined_data, avg_cng_quantile_RT, by = c("RT_quantile","congruency"))
  
  # Compute delta
  # contrasting extreme cng levels
  delta_data = group_by(combined_data, RT_quantile, subid) %>%          
    reframe(
      avg_quantile_RT = avg_quantile_RT[1],
      delta_acc_03 = accuracy[congruency==0]-accuracy[congruency==3],
      delta_acc_02 = accuracy[congruency==0]-accuracy[congruency==2],
      delta_acc_01 = accuracy[congruency==0]-accuracy[congruency==1],
      delta_rt_03 = mean_rt[congruency==0]-mean_rt[congruency==3],
      delta_rt_02 = mean_rt[congruency==0]-mean_rt[congruency==2],
      delta_rt_01 = mean_rt[congruency==0]-mean_rt[congruency==1],
    )
  
  return(list("means"=combined_data,"delta"=delta_data))
}








# Summarize one dimension and current congruency value within RT quantiles,
# additionally stratified by that dimension's congruency on the previous trial.
compute_quantiles_delta <- function(data, congruency_column, n_bins, congruency_value, prev_congruency) {
  # Filter by congruency
  filtered_data <- filter(data, !!sym(congruency_column) == congruency_value)
  
  # Compute quantiles for each subject and non-cued dimension congruency
  sbjs_quantile <- filtered_data %>%
    group_by(subid) %>%                        ## quantile is defined within sbj
    mutate(RT_quantile = as.factor(as.numeric(quantileCut(rt, n_bins))))
  
  # Left join quantiles and summarize data
  result <- left_join(filtered_data, sbjs_quantile, by = c(colnames(filtered_data))) %>%
    group_by(RT_quantile, subid, !!sym(prev_congruency)) %>%                              # delta-delta: stratify also by previous congruency level
    summarise(
      accuracy = mean(correct) * 100,
      count = n(),
      mean_rt = mean(rt),
      congruency = congruency_value
    )
  
  return(result)
}

# Compute a participant-level delta-delta function for one dimension: the
# current congruency effect after previous-congruent trials minus the current
# congruency effect after previous-incongruent trials.
process_delta_dimension <- function(data, congruency_column, dimension_label, n_bins, prev_congruency) {
  # Compute congruent and incongruent data
  congruent_data <- compute_quantiles_delta(data, congruency_column, n_bins, congruency_value = 1, prev_congruency)
  incongruent_data <- compute_quantiles_delta(data, congruency_column, n_bins, congruency_value = 0, prev_congruency)
  
  # Combine congruent and incongruent data
  combined_data <- bind_rows(congruent_data, incongruent_data)
  
  # Compute average RT for each quantile (and non-cued dimension). Note: average across congruency levels and across subjects! 
  avg_quantile_RT <- combined_data %>%
    group_by(RT_quantile) %>%
    summarise(avg_quantile_RT = weighted.mean(mean_rt, count))
  
  # Merge average quantile RT with combined data
  combined_data <- left_join(combined_data, avg_quantile_RT, by = "RT_quantile")
  
  # Compute delta 
  delta_data <- combined_data %>%
    group_by(RT_quantile, subid) %>%
    reframe(
      avg_quantile_RT = avg_quantile_RT[1],
      acc_con = accuracy[congruency == 1],
      acc_inc = accuracy[congruency == 0],
      ddelta_acc = (accuracy[congruency == 0 & !!sym(prev_congruency)==1] - accuracy[congruency == 1 & !!sym(prev_congruency)==1]) - 
                  (accuracy[congruency == 0 & !!sym(prev_congruency)==0] - accuracy[congruency == 1 & !!sym(prev_congruency)==0]), # delta-delta: prev.inc. - prev.con. i.e., (CI-CC)-(II-IC)
      ddelta_rt = (mean_rt[congruency==0 & !!sym(prev_congruency)==1] - mean_rt[congruency==1 & !!sym(prev_congruency)==1]) - 
                 (mean_rt[congruency==0 & !!sym(prev_congruency)==0] - mean_rt[congruency==1 & !!sym(prev_congruency)==0]),
      non_cued_dimension = dimension_label
    )
  
  return(delta_data)
}

# Apply `process_delta_dimension` to all four non-cued dimensions.
deltadelta_all_dimensions <- function(data, n_bins, Labels) {
  # Process each dimension
  dd_data_d1 <- process_delta_dimension(data, paste("congruency",tolower(Labels[1]),sep="_"), Labels[1], n_bins, "prev_0_congruent") # e.g., process_dimensions(data, "congruency_color", "Color", n_bins, "prev_0_congruent")
  dd_data_d2 <- process_delta_dimension(data, paste("congruency",tolower(Labels[2]),sep="_"), Labels[2], n_bins, "prev_1_congruent")
  dd_data_d3 <- process_delta_dimension(data, paste("congruency",tolower(Labels[3]),sep="_"), Labels[3], n_bins, "prev_2_congruent")
  dd_data_d4 <- process_delta_dimension(data, paste("congruency",tolower(Labels[4]),sep="_"), Labels[4], n_bins, "prev_3_congruent")
  
  # Combine all results
  deltadelta_data_eachNC <- bind_rows(dd_data_d1, dd_data_d2, dd_data_d3, dd_data_d4) %>%
    mutate(non_cued_dimension = factor(non_cued_dimension, levels = Labels, labels = Labels))
  
  return(deltadelta_data_eachNC)
}

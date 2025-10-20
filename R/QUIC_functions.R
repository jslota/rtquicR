#' Convert time values
#'
#' Converts time values from a funky string of text to proper numerical values in hours.
#' Used by the load_quic_results() function
#'
#' @param QuIC_data formatted QuIC results matrix with funky time values
#' @return Matrix with proper numerical values for time in hours
#' @importFrom stringr str_extract
convert_time_values <- function(QuIC_data) {
  time_str <- rownames(QuIC_data)
  
  # Extract hours and minutes using regex
  hours <- as.numeric(stringr::str_extract(time_str, "\\d+(?=\\s*h)"))
  minutes <- as.numeric(stringr::str_extract(time_str, "\\d+(?=\\s*min)"))
  
  # Replace NA with 0 where needed
  hours[is.na(hours)] <- 0
  minutes[is.na(minutes)] <- 0
  
  # Calculate time in hours
  time_hours <- hours + minutes / 60
  
  # Add to the data and set as rownames
  rownames(QuIC_data) <- time_hours
  
  return(QuIC_data)
}


#' Load and format data frame
#'
#' Takes raw data and formats it for further analysis
#'
#' @param input_file the path to the file containing raw RT-QuIC fluorescence data values in microsoft excel format (.xslx)
#' @param file_type How is the input file formatted? By default, file_type = "clean_table". file_type = "raw_table" for raw data off of the omega in table format.  file_type = "raw_microplate" for raw data off of the omega in microplate format.
#' @param excel_sheet (if more than one sheet) The name of the excel sheet that contains the formatted data
#' @param rows_to_skip (file_type = "raw_table" and file_type = "raw_microplate") The number of header rows to skip in the excel sheet before the data table begins.
#' @param unnecessary_columns (file_type = "raw_table" ONLY) The column number(s) of any unnecessary columns in the excel sheet. The only columns should be the Well ID, and individual reads.
#' @return Properly formatted matrix of QuIC data
#' @examples
#' \dontrun{load_quic_results("raw data/quic_results.xlsx")}
#' \dontrun{load_quic_results("raw data/quic_results.xlsx", filetype="raw_table", excel_sheet="All Cycles", rows_to_skip = 12, unnecessary_columns = c(2,3))}
#' @importFrom dplyr mutate
#' @importFrom readxl read_excel
#' @export
load_quic_results <- function(input_file, file_type = "clean_table", excel_sheet = NULL, rows_to_skip = NULL, unnecessary_columns = NULL) {

  # Check if the file exists
  if (!file.exists(input_file)) {
    stop("The file '", input_file, "' does not exist.", call. = FALSE)
  }

  if (file_type == "clean_table") {
    res <- as.data.frame(readxl::read_excel(input_file))
    if (colnames(res)[1] != "Well") {
      stop("The file '", input_file, "' is not formated as a clean table (First Column should be named 'Well')", call. = FALSE)
    }
  }  else if (file_type == "raw_table") {
    if (is.null(rows_to_skip)) {
      stop("Please specify 'rows_to_skip' when using filetype='raw_table' (Hint, try 'rows_to_skip'=12)", call. = FALSE)
    }
    if (is.null(unnecessary_columns)) {
      stop("Please specify 'unnecessary_columns' when using filetype='raw_table' (Hint, try 'unnecessary_columns'=c(2,3))", call. = FALSE)
    }
    res <- readxl::read_excel(input_file,
                              sheet = excel_sheet,
                              skip = rows_to_skip
                              )[,-unnecessary_columns]
    res <- as.data.frame(dplyr::rename(res, "Well" = `...1`))
    if (colnames(res)[1] != "Well") {
      stop("The file '", input_file, "' is not formated correctly. Check 'rows_to_skip' and 'unnecessary_columns' (First Column should be named 'Well')", call. = FALSE)
    }
  } else if (file_type == "raw_microplate") {
    if (is.null(rows_to_skip)) {
      stop("Please specify 'rows_to_skip' when using filetype='raw_microplate' (Hint, try 'rows_to_skip'=12)", call. = FALSE)
    }
    res <- readxl::read_excel(input_file,
                              sheet = excel_sheet,
                              skip = rows_to_skip)
    out <- data.frame(Well = sprintf("%s%02d", rep(LETTERS[1:8], each = 12), 1:12), # A01 through H12
                      "0 h" = c(as.numeric(res[3,2:13]), as.numeric(res[4,2:13]), as.numeric(res[5,2:13]), as.numeric(res[6,2:13]),
                                as.numeric(res[7,2:13]), as.numeric(res[8,2:13]), as.numeric(res[9,2:13]), as.numeric(res[10,2:13])),
                      check.names = FALSE
    )
    n <- (nrow(res)-10)/12
    new_names <- c("Well", "0 h")
    for (i in 1:n) {
      tmp <- res[(i*12):(i*12+10),]
      timepoint <- gsub("\\).*", "", gsub(".*\\(", "", as.character(tmp[1,1])))
      new_names <- c(new_names, timepoint)
      out <- dplyr::mutate(out, new = c(as.numeric(tmp[4,2:13]), as.numeric(tmp[5,2:13]), as.numeric(tmp[6,2:13]), as.numeric(tmp[7,2:13]),
                                        as.numeric(tmp[8,2:13]), as.numeric(tmp[9,2:13]), as.numeric(tmp[10,2:13]), as.numeric(tmp[11,2:13])))
      colnames(out) <- new_names
    }
    res <- out
    if (colnames(res)[1] != "Well") {
      stop("The file '", input_file, "' is not formated correctly. Check 'rows_to_skip'", call. = FALSE)
    }
  } else {
    stop("Incorrect value for 'file_type'. Valid values are 'clean_table', 'raw_table' and 'raw_microplate'", call. = FALSE)
  }

  rownames(res) <- res[,1]
  res <- as.data.frame(t(res[,-1]))
  res <- convert_time_values(res)
  return(res)
}

#' Format data for signal curve plot
#'
#' Takes matrix of QuIC fluorescence values and formats it with sample info to plot signal against time with ggplot
#' @param plot_data The formatted QuIC fluorescence data that will be used for plotting. Output of 'load_quic_results()'
#' @param plot_samples Matrix that contains matching sample information for each well in the matching QuIC plate
#' @param normalize Normalize method for fluorescence signal. Either normalize="none", normalize="max_RFU_per_plate", or normalize="baseline_RFU_per_well"
#' @param baseline_cycles Must be set only when normalize="baseline_RFU_per_well". Cycles used for baseline calculation. e.g. baseline_cycles=c(13:16)
#' @param smooth Whether to smooth the fluorescence curves. Needed when calculating lag times from slopes.
#' @return Formatted matrix of fluorescence for each well with sample info for signal curve plotting
#' @examples
#' \dontrun{sig_res <- signal_curve(QuIC_res, samples, normalize = "max_RFU_per_plate")}
#' \dontrun{sig_res <- signal_curve(QuIC_res, samples, normalize = "baseline_RFU_per_well", baseline_cycles = c(13:16), smooth = TRUE)}
#' @importFrom reshape2 melt
#' @importFrom signal sgolayfilt
#' @export
signal_curve <- function(plot_data, plot_samples, normalize = "none", baseline_cycles = NULL, smooth = FALSE) {
  
  # Basic checks
  if (!identical(colnames(plot_data), plot_samples$Well)) {
    stop("Well names in plot_data and plot_samples$Well must match.")
  }
  
  if (!identical(colnames(plot_samples)[1:3], c("Well", "Sample", "Dilution"))) {
    stop("plot_samples must have columns: 'Well', 'Sample', and 'Dilution'.")
  }
  
  # Normalization
  if (normalize == "max_RFU_per_plate") {
    plot_data <- 100 * (plot_data - min(plot_data)) / (max(plot_data) - min(plot_data))
    
  } else if (normalize == "baseline_RFU_per_well") {
    if (is.null(baseline_cycles)) {
      stop("normalize = 'baseline_RFU_per_well' requires 'baseline_cycles' (e.g., c(13:16)).")
    }
    
    plot_data[] <- lapply(plot_data, function(col) col / mean(col[baseline_cycles]))
  }
  
  # Smoothing
  if (smooth) {
    if (nrow(plot_data) < 17) {
      warning("Smoothing skipped: fewer than 17 timepoints.")
    } else {
      plot_data[] <- lapply(plot_data, function(col) signal::sgolayfilt(col, p = 3, n = 17))
    }
  }
  
  # Add time column and reshape
  plot_data$Time <- as.numeric(rownames(plot_data))
  plot_long <- reshape2::melt(plot_data, id.vars = "Time", variable.name = "Well", value.name = "Signal")
  
  # Ensure matching types before merging
  plot_long$Well <- as.character(plot_long$Well)
  plot_samples$Well <- as.character(plot_samples$Well)
  
  # Combine with sample metadata
  merged <- merge(plot_long, plot_samples, by = "Well")
  
  return(merged)
}


#' generic signal curve plot
#'
#' makes a generic plot with ggplot from the formatted signal curve plot data
#'
#' @param plot_data output of 'signal_curve()' function
#' @examples
#' \dontrun{plot_signal_curve(plot_data)}
#' @importFrom RColorBrewer brewer.pal
#' @import ggplot2
#' @export
plot_signal_curve <- function(plot_data) {

  if (identical(colnames(plot_data)[1:5],c("Well","Time","Signal","Sample","Dilution"))==FALSE) {
    stop("The file '", deparse(substitute(plot_data)), "' is not formated correctly. First must run 'signal_curve()'", call. = FALSE)
  }

  cols <- RColorBrewer::brewer.pal(n = length(levels(as.factor(plot_data$Dilution))), name = "Set3")

  #make plot
  ggplot2::ggplot(plot_data, aes(x= Time, y = Signal, group = Well, color = as.factor(Dilution))) +
    geom_line() +
    scale_color_manual(values = cols) +
    facet_wrap(~Sample, ncol = 2, dir = "v") +
    theme(panel.background = element_rect(fill = "white"),
          strip.background = element_rect(fill = "white"),
          axis.line = element_line(colour = "black"),
          panel.grid.major = element_line(color = "grey80"),
          panel.grid.minor = element_line(color = "grey92"),
          legend.position = "bottom",
          legend.key = element_rect(fill = "white", colour = "black"))
}

#' calculate max-point fluorescence ratios (MPR)
#'
#' @param data formatted RT-QuIC fluorescence data. Output of 'signal_curve'
#' @examples
#' \dontrun{mpr_res <- calc_mpr(signal_data)}
#' @return Outputs a matrix with the MPR values for each well with matching sample info data
#' @export
calc_mpr <- function(data) {
  if (identical(colnames(data)[1:5],c("Well","Time","Signal","Sample","Dilution"))==FALSE) {
    stop("The file '", deparse(substitute(data)), "' is not formated correctly. First must run 'signal_curve()'", call. = FALSE)
  }
  smpl.dt <- unique(data[,-c(2:3)])
  out <- list()
  for (i in unique(data$Well)) {
    out[[i]]$Well <- i
    out[[i]]$MPR <- max(data[data$Well == i,]$Signal)
  }
  out <- as.data.frame(do.call(rbind, out))
  out$Well <- as.character(out$Well)
  smpl.dt$Well <- as.character(smpl.dt$Well)
  out <- dplyr::left_join(out, smpl.dt, by="Well")
  out$MPR <- as.numeric(out$MPR)
  return(out)
}

#' calculate 1/lag phase
#'
#' calculates the inverse lag phase for each well of the QuIC plate (AKA amyloid formation rate)
#'
#' @param signal_data output of 'signal_curve()' function
#' @param cutoff The time cutoff in hours that will serve as the max time to reach threshold fluorescence
#' @param thresh_method The method for calculating threshold fluorescence. Default is "StdDev" (based on standard deviation). May also use "Mean" (2xmean(fluourescence)), "Max" (10% of max fluorescence), "Manual" (manually specify a number) or "MaxSlope" (based on curve slope instead of signal threshold).
#' @param n_StdDevs Only applies when thresh_method = "StdDev"; The number of standard deviations above baseline for thresholding (10 by default).
#' @param mean_FC Only applies when thresh_method = "Mean"; The fold-change above the baseline for thresholding (2 by default).
#' @param proportion_max Only applies when thresh_method = "Max"; The proportion of Max fluorescence for thresholding (0.1 by default).
#' @param threshold Only applies when thresh_method = "Manual". Numerical value to be used as threshold. 
#' @param thresh_calc_range Range in hours to be used as the basis for calculating threshold. By default thresh_calc_range = c(0,1).
#' @examples
#' \dontrun{lag_data <- calc_lag_phase(signal_data, 40)}
#' \dontrun{lag_data <- calc_lag_phase(signal_data, 40, thresh_method = "2xMean", thresh_calc_range = c(13:16))}
#' \dontrun{lag_data <- calc_lag_phase(signal_data, 40, thresh_method = "Manual", threshold=15000)}
#' \dontrun{lag_data <- calc_lag_phase(signal_data, 40, thresh_method = "MaxSlope")}
#' @return Outputs a matrix with the 1/lag values for each well with matching sample info data
#' @importFrom stats sd
#' @importFrom stats median
#' @importFrom dplyr filter summarise group_by pull reframe mutate left_join
#' @export
calc_lag_phase <- function(signal_data, cutoff, thresh_method = "StdDev", n_StdDevs=10, mean_FC=2, proportion_max=0.1, threshold, thresh_calc_range = c(0,1)) {
  
  # Validate input format
  if (identical(colnames(signal_data)[1:3], c("Well","Time","Signal"))==FALSE) {
    stop("The file '", deparse(substitute(data)), deparse(substitute(sample_info)), "' is not formated correctly. The first the columns must be 'Well', 'Time', and 'Signal'", call. = FALSE)
  }
  
  # Validate cutoff time
  if (is.numeric(cutoff)==FALSE || cutoff < 0 ) {
    stop("Invalid cycle cutoff time. Must be hours in numeric format and > 0.", call. = FALSE)
  }
  signal_data <- dplyr::filter(signal_data, Time < cutoff) # Filter by cutoff time
  
  # Baseline signal() used for threshold methods)
  baseline_signal <- dplyr::filter(signal_data, Time >= thresh_calc_range[1], Time <= thresh_calc_range[2])$Signal
  
  # Compute threshold
  thresh <- switch(thresh_method,
                   "StdDev" = mean(baseline_signal) + n_StdDevs * sd(baseline_signal),
                   "Mean"   = mean_FC * mean(baseline_signal),
                   "Max"    = mean(baseline_signal) + proportion_max * max(signal_data$Signal),
                   "Manual" = {
                     if (!is.numeric(threshold) || threshold < 0 || threshold > 500000) {
                       stop("threshold must be numeric and between 0 and 500000 for 'Manual' method.")
                     }
                     threshold
                   },
                   "MaxSlope" = NA,  # handled below
                   stop("Invalid thresh_method. Must be one of: StdDev, Mean, Max, Manual, MaxSlope.")
  )
  
  #outlier detection
  upper <- stats::median(baseline_signal) + 4*stats::sd(baseline_signal)
  lower <- stats::median(baseline_signal) - 4*stats::sd(baseline_signal)
  
  outlier_wells <- dplyr::filter(signal_data, Time >= thresh_calc_range[1], Time <= thresh_calc_range[2])
  outlier_wells <- dplyr::summarise(dplyr::group_by(outlier_wells, Well), bad = any(Signal < lower | Signal > upper))
  outlier_wells <- dplyr::pull(dplyr::filter(outlier_wells, bad), Well)
  
  if (length(outlier_wells) > 0) {
    message("Potential baseline outliers: ", paste(outlier_wells, collapse = ", "))
  }
  
  # Remove cycles that were used for threshold calculation (in case of abnormally high baseline)
  signal_data <- dplyr::filter(signal_data, Time > thresh_calc_range[2])
  
  lag_data <- dplyr::reframe(dplyr::group_by(signal_data, Well),
                             lag_time = {
                               sig <- Signal
                               time <- Time
                               
                               if (thresh_method == "MaxSlope") {
                                 d_sig <- diff(sig) / diff(time)
                                 max_slope_idx <- which.max(d_sig)
                                 
                                 if (length(max_slope_idx) == 1 && d_sig[max_slope_idx] > 0.25) {
                                   time[max_slope_idx]
                                 } else {
                                   NA_real_
                                 }
                               } else {
                                 first_above <- which(sig > thresh)[1]
                                 if (!is.na(first_above)) {
                                   time[first_above]
                                 } else {
                                   NA_real_
                                 }
                               }
                             }
  )
  
  lag_data <- dplyr::mutate(lag_data,
                            lag_time = ifelse(is.na(lag_time), 0, 1 / lag_time))
  
  # Join back metadata (non-time columns)
  meta_cols <- setdiff(colnames(signal_data), c("Time", "Signal"))
  lag_data <- dplyr::left_join(lag_data, unique(signal_data[meta_cols]), by = "Well")
  
  return(lag_data)
}

#' generic plot of lag phase data
#'
#' Makes a generic plot of lag phase data for each QuIC plate
#'
#' @param lag_data Output of calc_lag_phase() function
#' @examples
#' \dontrun{plot_lag_phase(lag_data)}
#' @importFrom RColorBrewer brewer.pal
#' @import ggplot2
#' @export
plot_lag_phase <- function(lag_data) {

  if (identical(colnames(lag_data)[1:4],c("Well", "lag_time", "Sample","Dilution"))==FALSE) {
    stop("The file '", deparse(substitute(lag_data)), "' is not formated correctly. First must run 'calc_lag_phase()'", call. = FALSE)
  }

  cols <- RColorBrewer::brewer.pal(n = length(levels(as.factor(lag_data$Dilution))), name = "Set3")

  ggplot2::ggplot(rev(lag_data), aes(x=as.factor(Dilution), y=lag_time, color=as.factor(Dilution))) +
    geom_boxplot(color = "black", outlier.color = "white") +
    geom_jitter(size = 2, width = 0.15) +
    scale_color_manual(values = cols) +
    xlab("Dilution") +
    ylab("(lag phase)-1") +
    facet_wrap(~Sample, scales = "free_x") +
    theme(panel.background = element_rect(fill = "white"),
          axis.text.x = element_text(angle = 90),
          axis.line = element_line(color = "black"),
          legend.position = "bottom",
          legend.title = element_blank(),
          legend.key = element_rect(fill = "white"),
          legend.margin = margin(-0.3,0,0,0, "cm"))
}

#' calculate SD50 values
#'
#' Calculates SD50 for each sample from the lag phase data. 3/4 wells for the highest dilution must be
#' positive for the function to return SD50 measurements (This serves as the limit of detection).
#' Otherwise returns the maximum possible value below the limit of detection.
#'
#' @param lag_data Output of `calc_lag_phase()` with columns: Well, lag_time, Sample, Dilution
#' @param starting_dilution Optional numeric. Use only dilutions ≤ this value (e.g., for inhibition at high concentration)
#' @param positivity_threshold Proportion of positive wells required in highest dilution. Default = 0.75
#' @return Data frame with SD50 metrics and result interpretation
#' @examples
#' \dontrun{sd50_res <- calc_SD50(lag_data)}
#' \dontrun{sd50_res <- calc_SD50(lag_data, starting_dilution=1e-04)}
#' @importFrom dplyr case_when distinct mutate left_join
#' @importFrom purrr map_dfr
#' @export
calc_SD50 <- function(lag_data, starting_dilution = NULL, positivity_threshold = 0.75) {
  
  # ---- Input checks ----
  required_cols <- c("Well", "lag_time", "Sample", "Dilution")
  if (!identical(colnames(lag_data)[1:4], required_cols)) {
    stop("Input lag_data must have columns: Well, lag_time, Sample, Dilution")
  }
  
  if (!is.null(starting_dilution)) {
    if (!is.numeric(starting_dilution) || starting_dilution <= 0 || starting_dilution >= 1) {
      stop("starting_dilution must be a numeric value between 0 and 1, e.g. 1e-04")
    }
    lag_data <- lag_data[lag_data$Dilution <= starting_dilution, ]
  }
  
  # ---- Helper for monotonic smoothing ----
  smooth_monotonic <- function(p_vec) {
    while (any(cummin(p_vec) != p_vec)) {
      for (k in 2:length(p_vec)) {
        if (p_vec[k] > p_vec[k - 1]) {
          idx <- which(1:length(p_vec) < k & p_vec < p_vec[k])
          idx <- c(idx, k)
          p_vec[idx] <- mean(p_vec[idx])
        }
      }
    }
    return(p_vec)
  }
  
  # ---- Main function for calculations ----
  get_sd50 <- function(lag_data, sample_id) {
    tmp <- lag_data[lag_data$Sample == sample_id, ]
    tmp <- tmp[order(tmp$Dilution, decreasing = TRUE), ]
    
    dilutions <- unique(tmp$Dilution)
    d_vals <- round(diff(-log10(dilutions)), 3)
    
    if (length(unique(d_vals)) != 1) {
      stop(paste("Dilutions for sample", sample_id, "are not evenly spaced in log10 scale."))
    }
    
    # Filter for positivity in first dilution
    first_dil <- max(tmp$Dilution)
    if (sum(tmp[tmp$Dilution == first_dil, "lag_time"] > 0) < positivity_threshold * sum(tmp$Dilution == first_dil)) {
      tmp$lag_time <- 0  # Remove "positives"
    }
    
    d <- unique(d_vals)
    
    # Find highest dilution with all positive wells
    x0 <- 0
    for (dil in dilutions) {
      pos_count <- sum(tmp$Dilution == dil & tmp$lag_time > 0)
      total_count <- sum(tmp$Dilution == dil)
      if (pos_count == total_count && total_count > 0) {
        x0 <- -log10(dil)
        break
      }
    }
    
    if (x0 == 0) {
      x0 <- -log10(10 * max(dilutions))  # Below LOD
    }
    
    # Compute SD50 using Spearman-Kärber
    included_dilutions <- dilutions[dilutions <= 1.01 * 10^(-x0)]
    p_vector <- numeric()
    n_vector <- numeric()
    
    for (dil in included_dilutions) {
      wells <- tmp[tmp$Dilution == dil, ]
      p_vector <- c(p_vector, mean(wells$lag_time > 0))
      n_vector <- c(n_vector, nrow(wells))
    }
    
    if (length(p_vector) > 1) {
      p_vector <- smooth_monotonic(p_vector)
    }
    
    # If SE = 0, assume 1 well negative at x0 dilution
    p_se <- sqrt(d^2 * sum(p_vector * (1 - p_vector) / (n_vector - 1)))
    if (p_se == 0 && x0 > -log10(10 * max(dilutions))) {
      p_min <- 1 / n_vector[1]
      p_se <- sqrt(d^2 * sum(p_min * (1 - p_min) / (n_vector[1] - 1)))
    }
    
    # Final SD50 calculations
    log10_sd50 <- x0 - d / 2 + d * sum(p_vector)
    
    sd50_res <- data.frame(Sample = sample_id,
                           log10_SD50 = log10_sd50,
                           SE = p_se,
                           CI95_upper = log10_sd50 + 1.96 * p_se,
                           CI95_lower = log10_sd50 - 1.96 * p_se,
                           LOD_min = -log10(10 * max(dilutions)) - d / 2 + d * (1 / n_vector[1]),
                           LOD_max = -log10(min(dilutions)) + d / 2)
    
    # Determine result status
    sd50_res <- dplyr::mutate(sd50_res,
                              Result = dplyr::case_when(log10_SD50 < LOD_min ~ "negative",
                                                        log10_sd50 > LOD_max ~ "saturated",
                                                        .default = "positive"))
    return(sd50_res)
  }
  
  combined_output <- purrr::map_dfr(unique(lag_data$Sample), ~ get_sd50(lag_data, .x))
  
  # Merge with meta-data
  meta_data <- dplyr::distinct(lag_data[,setdiff(names(lag_data), c("Well", "lag_time", "Dilution"))])
  combined_output <- dplyr::left_join(combined_output,
                                      meta_data,
                                      by = "Sample")
  
  return(combined_output)
}

#' calculate AUC values from Lag phase data
#'
#' Calculates lag phase AUC for each sample from the lag phase data. Converts dilution to -log10(dilution) for AUC calculation
#'
#' @param lag_data Output of calc_lag_phase() function
#' @return Outputs AUC lag-phase-curve values for each sample with matching sample info
#' @examples
#' \dontrun{auc_res <- calc_AUC_lag(lag_data)}
#' @importFrom pracma trapz
#' @importFrom stats aggregate
#' @export
calc_AUC_lag <- function(lag_data) {
  if (identical(colnames(lag_data)[1:4],c("Well", "lag_time", "Sample","Dilution"))==FALSE) {
    stop("The file '", deparse(substitute(lag_data)), "' is not formated correctly. First must run 'calc_lag_phase()'", call. = FALSE)
  }
  per_sample_data <- lag_data
  per_sample_data <- per_sample_data[per_sample_data$lag_time != Inf,]
  per_sample_data$Dilution <- -log10(per_sample_data$Dilution)
  per_sample_data <- stats::aggregate(lag_time ~ Dilution+Sample, per_sample_data, mean)

  #Make dataframe for AUC results
  AUC_res <- lag_data[,-c(1,2,4)]
  AUC_res <- unique(AUC_res)
  rownames(AUC_res) <- AUC_res$Sample
  AUC_res$AUC <- 0

  #get AUC values
  for (i in unique(per_sample_data$Sample)) {
    tmp <- per_sample_data[per_sample_data$Sample == i,]
    AUC_res[AUC_res$Sample == i,]$AUC <- pracma::trapz(x=tmp$Dilution, y=tmp$lag_time)
    rm(tmp)
  }
  return(AUC_res)
}

#' calculate AUC values from signal curves
#'
#' Calculates lag phase AUC for each well from the signal curve.
#'
#' @param sig_data Output of signal_curve() function
#' @return Outputs AUC values from signal curves for each well with matching sample info
#' @examples
#' \dontrun{auc_res <- calc_AUC_sig(sig_data)}
#' @importFrom pracma trapz
#' @export
calc_AUC_sig <- function(sig_data) {
  if (identical(colnames(sig_data)[1:5],c("Well","Time","Signal","Sample","Dilution"))==FALSE) {
    stop("The file '", deparse(substitute(sig_data)), "' is not formated correctly. First must run 'signal_curve()'", call. = FALSE)
  }
  #Make dataframe for AUC results
  AUC_res <- sig_data[,-c(2,3)]
  AUC_res <- unique(AUC_res)
  rownames(AUC_res) <- AUC_res$Well
  AUC_res$AUC <- 0

  #Get AUC values
  for (i in unique(sig_data$Well)) {
    tmp <- sig_data[sig_data$Well==i,]
    AUC_res[i,]$AUC <- pracma::trapz(x=tmp$Time, y=tmp$Signal)
  }

  return(AUC_res)
}


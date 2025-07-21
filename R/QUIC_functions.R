#' Convert time values
#'
#' Converts time values from a funky string of text to proper numerical values in hours.
#' Used by the load_quic_results() function
#'
#' @param QuIC_data formatted QuIC results matrix with funky time values
#' @return Matrix with proper numerical values for time in hours
convert_time_values <- function(QuIC_data) {
  QuIC_data$Time <- 0
  for (i in rownames(QuIC_data)) {
    if (length(strsplit(rownames(QuIC_data[i,]), " ")[[1]]) < 3) {
      QuIC_data[i,]$Time <- as.numeric(strsplit(rownames(QuIC_data[i,]), " ")[[1]][1])
    } else {
      QuIC_data[i,]$Time <- as.numeric(strsplit(rownames(QuIC_data[i,]), " ")[[1]][1]) + as.numeric(strsplit(rownames(QuIC_data[i,]), " ")[[1]][3])/60
    }

  }
  rownames(QuIC_data) <- QuIC_data$Time
  rownames(QuIC_data)
  QuIC_data <- QuIC_data[,colnames(QuIC_data) != "Time"]
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
    if (is.null(rows_to_skip)==TRUE) {
      stop("Please specify 'rows_to_skip' when using filetype='raw_table' (Hint, try 'rows_to_skip'=12)", call. = FALSE)
    }
    if (is.null(unnecessary_columns)==TRUE) {
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
    if (is.null(rows_to_skip)==TRUE) {
      stop("Please specify 'rows_to_skip' when using filetype='raw_microplate' (Hint, try 'rows_to_skip'=12)", call. = FALSE)
    }
    res <- readxl::read_excel(input_file,
                              sheet = excel_sheet,
                              skip = rows_to_skip)
    out <- data.frame(Well = c("A01", "A02", "A03", "A04", "A05", "A06", "A07", "A08", "A09", "A10", "A11", "A12",
                               "B01", "B02", "B03", "B04", "B05", "B06", "B07", "B08", "B09", "B10", "B11", "B12",
                               "C01", "C02", "C03", "C04", "C05", "C06", "C07", "C08", "C09", "C10", "C11", "C12",
                               "D01", "D02", "D03", "D04", "D05", "D06", "D07", "D08", "D09", "D10", "D11", "D12",
                               "E01", "E02", "E03", "E04", "E05", "E06", "E07", "E08", "E09", "E10", "E11", "E12",
                               "F01", "F02", "F03", "F04", "F05", "F06", "F07", "F08", "F09", "F10", "F11", "F12",
                               "G01", "G02", "G03", "G04", "G05", "G06", "G07", "G08", "G09", "G10", "G11", "G12",
                               "H01", "H02", "H03", "H04", "H05", "H06", "H07", "H08", "H09", "H10", "H11", "H12"),
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
  
  if ((identical(colnames(plot_data),plot_samples$Well)==FALSE)) {
    stop("The files '", deparse(substitute(plot_data)), "' and '", deparse(substitute(plot_samples)), "' are not formated correctly. Must first run 'load_quic_results()' to generate '",
         deparse(substitute(plot_data)), "'. '", deparse(substitute(plot_samples)), "' must have a column named 'Well' with values that match colnames of '", deparse(substitute(plot_data)), "'.", call. = FALSE)
  }
  if (identical(colnames(plot_samples)[1:3], c("Well","Sample","Dilution"))==FALSE) {
    stop("The file '", deparse(substitute(plot_samples)), "' is not formated correctly. The first the columns must be 'Well', 'Sample', and 'Dilution'", call. = FALSE)
  }
  
  if (normalize == "max_RFU_per_plate") {
    plot_data <- 100*(plot_data - min(plot_data))/(max(plot_data) - min(plot_data))
  } else  if (normalize == "baseline_RFU_per_well") {
    if (is.null(baseline_cycles) == TRUE) {
      stop("You have selected normalize=baseline_RFU_per_plate, but have not provided baseline_cycles. Eg. baseline_cycles=c(13:16)")
    }
    for (i in colnames(plot_data)) {
      plot_data[,i] <- plot_data[,i]/mean(plot_data[baseline_cycles,i])
    }
    
  }
  
  if (smooth == TRUE) {
    for (i in colnames(plot_data)) {
      plot_data[,i] <- signal::sgolayfilt(plot_data[,i], p = 3, n = 17)
    }
    
  }
  
  #format data for plotting
  plot_data$time <- rownames(plot_data)
  plot_data <- reshape2::melt(plot_data, id.vars = "time")
  colnames(plot_data) <- c("Time", "Well", "Signal")
  #convert time to numeric
  plot_data$Time <- as.numeric(plot_data$Time)
  ###add sample info
  plot_data <- merge.data.frame(plot_data, plot_samples, by = "Well")
  return(plot_data)
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
#' @param thresh_method The method for calculating threshold fluorescence. Default is "StdDev" (based on standard deviation). May also use "Mean" (2xmean(fluourescence)), "Max" (10% of max fluorescence) or "Manual" (manually specify a number).
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
#' @importFrom dplyr filter
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
#' @param lag_data Output of calc_lag_phase() function
#' @param starting_dilution The numeric value for the starting dilution used to calculate SD50. Useful if there is inhibition at higher sample concentrations.
#' @param positivity_threshold The proportion of positive wells in the first dilution must be above this threshold for a sample to be considered "positive". 0.75 by default
#' @return Outputs SD50 values, standard error, and 95% confidence interval for each sample with matching sample info. Samples deemed to be "negative" are returned with an SD50 value equivalent to 1 dilution lower than the minimum dilution.
#' @examples
#' \dontrun{sd50_res <- calc_SD50(lag_data)}
#' \dontrun{sd50_res <- calc_SD50(lag_data, starting_dilution=1e-04)}
#' @export
calc_SD50 <- function(lag_data, starting_dilution = NULL, positivity_threshold = 0.75) {

  #Make sure correct input file was provided
  if (identical(colnames(lag_data)[1:4],c("Well", "lag_time", "Sample","Dilution"))==FALSE) {
    stop("The file '", deparse(substitute(lag_data)), "' is not formated correctly. First must run 'calc_lag_phase()'", call. = FALSE)
  }

  #Check for valid starting dilution, then remove less-dilute replicates
  if (is.null(starting_dilution)==FALSE) {
    if (is.numeric(starting_dilution)==FALSE) {
      stop("starting_dilution must be a numerical value. e.g. 1e-04")
    } else if (starting_dilution > 1) {
      stop("Starting dilution must be less than one and greater than 0, e.g. 1e-03, 1e-04, 1e-05...")
    } else if (starting_dilution < 0) {
      stop("Starting dilution must be less than one and greater than 0, e.g. 1e-03, 1e-04, 1e-05...")
    }
    lag_data <- lag_data[lag_data$Dilution <= starting_dilution,]
  }

  #setup output data frame
  keep <- is.element(colnames(lag_data), c("Well", "lag_time", "Dilution")) == FALSE
  sd50_res <- lag_data[,keep]
  sd50_res <- unique(sd50_res)
  rownames(sd50_res) <- sd50_res$Sample
  sd50_res$log10_SD50 <- 0
  sd50_res$SE <- 0
  sd50_res$CI95_upper <- 0
  sd50_res$CI95_lower <- 0
  sd50_res$LOD_min <- 0
  sd50_res$LOD_max <- 0
  sd50_res$Result <- "positive"

  #calculate SD50 for each sample
  for (i in unique(lag_data$Sample)) {
    tmp <- lag_data[lag_data$Sample == i,]

    #Identify d; d=difference between dilutions in log10 scale (e.g. d=1 for 10-fold dilution series)
    d <- unique(round(diff(-log10(unique(tmp[order(tmp$Dilution, decreasing = TRUE),]$Dilution))),3))
    if (length(d) > 1) {
      stop("Please check dilution series... dilutions must be at regular intervals, e.g. 1e-03, 1e-04, 1e-05..")
    }

    #Calculate x0 (lowest dilution to be 4/4 positive)
    dilutions = unique(tmp$Dilution)
    x0 = 0
    for (j in 1:length(dilutions)) {
      if(sum(tmp[tmp$Dilution == dilutions[j],]$lag_time > 0) == nrow(tmp[tmp$Dilution == dilutions[j],])) {
        x0 <- -log10(dilutions[j])
      }
    }
    if (x0 == 0) {
      x0 <- -log10(10*max(dilutions))
    }

    #Remove false positives, need at least 3/4 first-dilution rxns to be positive
    if (sum(tmp[tmp$Dilution == max(tmp$Dilution),]$lag_time > 0) < (positivity_threshold*nrow(tmp[tmp$Dilution == max(dilutions),]))) {
      tmp$lag_time <- 0
    }

    #Calculate SD50; method 1
    #ni <- nrow(tmp)/length(unique(tmp$Dilution)) #number of replicates per dilution
    #Remove dilutions higher than x0
    #tmp <- tmp[tmp$Dilution <= 1.01*10^(-x0),] #the 1.01 factor is just to make sure the x0 dilution is included
    #ri <- sum(tmp$lag_time > 0) #total number of positive wells at X0 dilution and below
    #sd50_res[i,]$log10_SD50 <- (x0 - d/2 + d*(ri/ni)) #Calculate SD50

    #Calculate SD50; method 2
    dilutions.p <- dilutions[dilutions <= 1.01*10^(-x0)] #get all dilutions above x0, including x0
    p.sum <- numeric()
    n <- numeric()
    #p.se <- 0
    for (j in dilutions.p) { #Get p for each dilution >= x0
      p.tmp <- sum(tmp[tmp$Dilution==j,]$lag_time > 0)/nrow(tmp[tmp$Dilution==j,]) #Calculate proportion positive for this dilution
      p.sum <- c(p.sum, p.tmp) #add to sum
      n <- c(n, nrow(tmp[tmp$Dilution==j,]))
      #p.se <- p.se + (p.tmp*(1-p.tmp)/(nrow(tmp[tmp$Dilution==j,])-1))
    }

    #"Smooth" p... only necessary for standard error calculation
    if (length(p.sum) > 1) {
      while (any(cummin(p.sum) != p.sum)) { # while p is non-monotonic
        for (k in 2:length(p.sum)) { # starting from the second dilution and going up,
          if (p.sum[k] > p.sum[k-1]) { # if the current one has higher p than the last one...
            # select all dilutions below the current one that have lower p...
            indices_to_average = c(which(1:length(p.sum) < k & p.sum < p.sum[k]), k)
            # and average them all together with the current value.
            p.sum[indices_to_average] = mean(p.sum[indices_to_average])
          }
        }
      }
    }

    p.se <- sqrt(d^2 * sum(p.sum*(1-p.sum)/(n-1))) # calculate standard error
    if (p.se == 0 && x0 > -log10(10*max(dilutions))) { #If SE=0, set to one replicate negative at the x0 dilution for SE calculation
      p.min <- 1/n[1]
      p.se <- sqrt(d^2 * sum(p.min*(1-p.min)/(n[1]-1)))
    }

    #Calculate SD50
    sd50_res[i,]$log10_SD50 <- (x0 - d/2 + d*sum(p.sum))

    #Calculate standard error and confidence intervals
    sd50_res[i,]$SE <- p.se
    sd50_res[i,]$CI95_upper <- sd50_res[i,]$log10_SD50 + 1.96*sd50_res[i,]$SE
    sd50_res[i,]$CI95_lower <- sd50_res[i,]$log10_SD50 - 1.96*sd50_res[i,]$SE

    #Calculate limits of detection
    sd50_res[i,]$LOD_min <- (-log10((10^d)*max(dilutions)) - d/2 + d*(1/nrow(tmp[tmp$Dilution == max(dilutions),])))
    sd50_res[i,]$LOD_max <- (-log10(min(dilutions)) + d/2)

    #Determine if result is negative, or if RT-QuIC assay was saturated
    if (sd50_res[i,]$log10_SD50 < sd50_res[i,]$LOD_min) {
      sd50_res[i,]$Result <- "negative"
    } else  if (sd50_res[i,]$log10_SD50 == sd50_res[i,]$LOD_max) {
      sd50_res[i,]$Result <- "saturated"
    }

    rm(tmp)
  }
  return(sd50_res)
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

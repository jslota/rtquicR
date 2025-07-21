#' calculate 1/lag phase
#'
#' calculates the inverse lag phase for each well of the QuIC plate (AKA amyloid formation rate)
#'
#' @param data formatted RT-QuIC fluorescence data. Output of 'load_quic_results()'
#' @param sample_info Matrix that contains matching sample information for each well in the matching QuIC plate
#' @param cutoff The time cutoff in hours that will serve as the max time to reach threshold fluorescence
#' @param thresh_method The method for calculating threshold fluorescence. Default is "StdDev" (based on standard deviation). May also use "Mean" (2xmean(fluourescence)), "Max" (10% of max fluorescence) or "Manual" (manually specify a number).
#' @param n_StdDevs Only applies when thresh_method = "StdDev"; The number of standard deviations above baseline for thresholding (10 by default).
#' @param mean_FC Only applies when thresh_method = "Mean"; The fold-change above the baseline for thresholding (2 by default).
#' @param proportion_max Only applies when thresh_method = "Max"; The proportion of Max fluorescence for thresholding (0.1 by default).
#' @param threshold Only applies when thresh_method = "Manual". Numerical value to be used as threshold. 
#' @param thresh_calc_range Range of cycle numbers to be used as the basis for calculating threshold. By default thresh_calc_range = c(1:4).
#' @examples
#' \dontrun{lag_data <- calc_lag_phase(res, samples, 40)}
#' \dontrun{lag_data <- calc_lag_phase(res, samples, 40, thresh_method = "2xMean", thresh_calc_range = c(13:16))}
#' \dontrun{lag_data <- calc_lag_phase(res, samples, 40, thresh_method = "Manual", threshold=15000)}
#' @return Outputs a matrix with the 1/lag values for each well with matching sample info data
#' @importFrom stats sd
#' @importFrom stats median
#' @export
old_calc_lag_phase <- function(data, sample_info, cutoff, thresh_method = "StdDev", n_StdDevs=10, mean_FC=2, proportion_max=0.1, threshold, thresh_calc_range = c(1:4)) {
  
  if ((identical(colnames(data),sample_info$Well)==FALSE)) {
    stop("The files '", deparse(substitute(data)), "' and '", deparse(substitute(sample_info)), "' are not formated correctly. Must first run 'load_quic_results()' to generate '",
         deparse(substitute(data)), "'. '", deparse(substitute(sample_info)), "' must have a column named 'Well' with values that match colnames of '", deparse(substitute(data)), "'.", call. = FALSE)
  }
  if (is.numeric(as.matrix(data))==FALSE) {
    stop("The file '", "' is not formatted correctly. Must only contain numeric values.")
  }
  if (identical(colnames(sample_info)[1:3], c("Well","Sample","Dilution"))==FALSE) {
    stop("The file '", deparse(substitute(data)), deparse(substitute(sample_info)), "' is not formated correctly. The first the columns must be 'Well', 'Sample', and 'Dilution'", call. = FALSE)
  }
  
  #cutoff time
  if (is.numeric(cutoff)==FALSE || cutoff < 0 ) {
    stop("Invalid cycle cutoff. Must be numeric and > 0.", call. = FALSE)
  }
  data <- data[as.numeric(rownames(data)) < as.numeric(cutoff),]
  
  if(thresh_method == "StdDev") {
    #threshold = mean(negative controls) + 10 standard deviation
    thresh <- as.numeric(mean(rowMeans(data)[thresh_calc_range]) + n_StdDevs*mean(apply(data, 1, stats::sd)[thresh_calc_range]))
  } else if(thresh_method == "Mean") {
    thresh <- as.numeric(mean_FC*mean(rowMeans(data)[thresh_calc_range]))
  } else if(thresh_method == "Max") {
    thresh <- as.numeric(mean(rowMeans(data)[thresh_calc_range]) + proportion_max*max(data))
  } else if(thresh_method == "Manual") {
    if (is.numeric(threshold)==FALSE || threshold < 0 || threshold > 500000) {
      stop("Must specify a numeric value for 'thresold' between 0 and 500,000 when thresh_method='Manual'", call. = FALSE)
    }
    thresh <- threshold
  } else  {
    stop("Invalid value for 'thresh_method', valid values are 'StdDev', 'Mean', 'Max', and 'Manual'", call. = FALSE)
  }
  
  #outlier detection
  for (i in thresh_calc_range) {
    upper <- stats::median(as.numeric(data[i,])) + 4*stats::sd(data[i,])
    lower <- stats::median(as.numeric(data[i,])) - 4*stats::sd(data[i,])
    for (j in colnames(data)) {
      if(data[i,j] > upper | data[i,j] < lower) {
        print(paste0("Well ", j, " might be an outlier"))
      }
    }
    
  }
  
  #Remove cycles that were used for threshold calculation (in case of abnormally high baseline)
  data <- data[max(thresh_calc_range):nrow(data),]
  
  #make data frame to collect time-points
  time_data <- data.frame(Well = colnames(data), lag_time = 0)
  
  #get first time-point above threshold
  for (i in 1:ncol(data)) {
    if (is.null(rownames(data[data[,i] > thresh,])) == FALSE) {
      time_data$lag_time[i] <- as.numeric(rownames(data[data[,i] > thresh,])[1])
    } else  {
      time_data$lag_time[i] <- NA
    }
  }
  #1/lag time
  time_data$lag_time <- 1/time_data$lag_time
  #Set negative wells to 1/lag = 0
  time_data[is.na(time_data$lag_time),]$lag_time <- 0
  #get sample info
  time_data <- merge.data.frame(time_data, sample_info, by = "Well")
  
  return(time_data)
}

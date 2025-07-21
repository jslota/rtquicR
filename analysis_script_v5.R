library(rtquicR)
library(readxl)
library(ggplot2)
library(RColorBrewer)
library(dplyr)

#Files to be analyzed
layout_files <- Sys.glob("raw data/*Layout.xlsx")
results_files <- Sys.glob("raw data/*Results.xlsx")

layout_data <- list()
results_data <- list()

#Load and format data for each plate
for (i in 1:length(layout_files)) {
  print(paste0(layout_files[i], " ... ", results_files[i]))
  plate <- strsplit(results_files[i], "_")[[1]][2]
  
  ##################Specify parameters here#######################
  #load and format QuIC data
  res <- load_quic_results(input_file = results_files[i],
                           file_type = "raw_table", #Format of input data; "clean_table", "raw_table", or "raw_microplate"
                           excel_sheet = "All Cycles", #Name of excel sheet with data
                           rows_to_skip = 12, #rows to skip in header of excel sheet
                           unnecessary_columns = c(2,3) #non data columns in excel sheet
  )
  ##################Specify parameters here#######################
  #Add plate number to well IDs
  colnames(res) <- paste0("P", plate, "-", colnames(res))
  #load and format plate layout data
  samples <- read_excel(layout_files[i])
  samples$Well <- paste0("P", plate, "-", samples$Well)
  #Remove empty wells
  if(any(samples$Smpl_Ctrl=="Empty")) {
    empty <- samples %>% filter(Smpl_Ctrl=="Empty") %>% pull(Well)
    samples <- samples %>% filter(!Well %in% empty)
    res <- res[,-which(colnames(res) %in% empty)]
  }
  #add to main dataset
  results_data[[paste("P",plate, sep = "-")]] <- res
  layout_data[[paste("P",plate, sep = "-")]] <- samples
  #remove temp datasets
  rm(res, samples, plate)
}

#Remove poor quality data
# results_data[["P-3"]] <- results_data[["P-3"]][,-which(colnames(results_data[["P-3"]])=="P3-F04")]
# layout_data[["P-3"]] <- layout_data[["P-3"]][-which(layout_data[["P-3"]]$Well=="P3-F04"),]
# results_data[["P-3"]] <- results_data[["P-3"]][,-which(colnames(results_data[["P-3"]])=="P3-G12")]
# layout_data[["P-3"]] <- layout_data[["P-3"]][-which(layout_data[["P-3"]]$Well=="P3-G12"),]

if (dir.exists("signal curve plots/")==FALSE) { dir.create("signal curve plots/") }
###Make signal plots for each QuIC plate
sig_dt <- list()
for (i in names(results_data)) {
  ##################Specify parameters here#######################
  plot_data <- signal_curve(plot_data = results_data[[i]],
                            plot_samples = layout_data[[i]],
                            normalize = "baseline_RFU_per_well", #Normalization method; "none", "max_RFU_per_plate", or "baseline_RFU_per_well"
                            smooth = TRUE, # whether to smooth signal curve
                            baseline_cycles = c(5:9) #Cycles for baseline calculation; only used when normalize="baseline_RFU_per_well"
  )
  ##################Specify parameters here#######################
  sig_dt[[i]] <- plot_data
  print(paste0("Plotting signal curves for plate ", i))
  plot_signal_curve(plot_data)
  ggplot2::ggsave(paste0("signal curve plots/Plate ", i, ".pdf"), width = 6.5, units = "in")
  rm(plot_data)
}


if (dir.exists("lag phase plots/")==FALSE) { dir.create("lag phase plots/") }

###Make lag phase plot for each QuIC plate
full_data <- list()
for (i in names(sig_dt)) {
  #get data
  print(paste0("calculating lag-phase for plate ", i))
  
  ##################Specify parameters here#######################
  #Lag phase calculation
  lag_data <- calc_lag_phase(signal_data = sig_dt[[i]],
                             cutoff = 100, # Max hours cut-off for lag-phase calculation
                             thresh_method = "MaxSlope", #Method for calculating threshold: "StdDev", "Mean", "Max", "Manual", or "MaxSlope"
                             #n_StdDevs = 10, # Specify number of standard deviations above baseline for threshold
                             thresh_calc_range = c(1,3)) #Cycle numbers to use for threshold calculation
  ##################Specify parameters here#######################
  
  full_data[[i]] <- lag_data
  print(paste0("Plotting lag-phase curves for plate ", i))
  plot_lag_phase(lag_data)
  ggplot2::ggsave(paste0("lag phase plots/Plate ", i, " lag.pdf"), width = 9, height = 6.5, units = "in")
  rm(lag_data)
}

sig_dt <- do.call(rbind, sig_dt)

#save results
if (dir.exists("results/")==FALSE) { dir.create("results/") }
write.csv(sig_dt, "results/signal_curve_results.csv", row.names = FALSE)
lag_dat <- do.call(rbind, full_data)
lag_dat <- lag_dat %>% arrange(Sample, desc(Dilution))
write.csv(lag_dat, "results/lag_phase_results.csv", row.names = FALSE)

#AUC of signal curves
auc_sig_res <- sig_dt %>% filter(Smpl_Ctrl=="Smpl")
auc_sig_res <- calc_AUC_sig(sig_data = auc_sig_res)
write.csv(auc_sig_res, "results/auc_signal_curve_results.csv", row.names = FALSE)

#Max-point fluorescent signal data
mpr_res <- sig_dt %>% filter(Smpl_Ctrl=="Smpl")
mpr_res <- calc_mpr(data = mpr_res)
write.csv(mpr_res, "results/mpr_results.csv", row.names = FALSE)

#SD50 analysis
sd50_res <- lag_dat %>% filter(Smpl_Ctrl=="Smpl", is.na(culture_type) == FALSE)
##################Specify parameters here#######################
sd50_res <- calc_SD50(lag_data = sd50_res,
                      starting_dilution = 1e-04, #Starting dilution to be used for sd50 calculation
                      positivity_threshold = 0.75) #The first dilution must have this proportion of positive wells for positive samples
##################Specify parameters here#######################
write.csv(sd50_res, "results/sd50_results.csv", row.names = FALSE)

#AUC of lag phase curves
auc_lag_res <- lag_dat %>% filter(Smpl_Ctrl=="Smpl")
auc_lag_res <- calc_AUC_lag(lag_data = auc_lag_res)
write.csv(auc_lag_res, "results/auc_lag_phase_curve_results.csv", row.names = FALSE)

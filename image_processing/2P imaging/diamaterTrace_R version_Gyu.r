rm(list = ls())

# Full-Width at Half-Maximum (FWHM) of the wave form y(x) and its polarity.
# fwhm results in 'width will be in units of 'x'
# Rev 1.2, April 2006 (Patrick Egan)

fwhm <- function(x, y) {
  # Normalize y
  y <- y / max(y)
  N <- length(y)
  lev50 <- 0.5
  if (y[1] < lev50) {  # find index of center (max or min) of pulse
    centerindex <- which.max(y)
    Pol <- 1
    # print('Pulse Polarity = Positive')
  } 
  else {
    centerindex <- which.min(y)
    Pol <- -1
    # print('Pulse Polarity = Negative')
  }
  
  i <- 2
  while (sign(y[i] - lev50) == sign(y[i - 1] - lev50)) {
    i <- i + 1
  }  # first crossing is between v(i-1) & v(i)
  
  interp <- (lev50 - y[i - 1]) / (y[i] - y[i - 1])
  tlead <- x[i - 1] + interp * (x[i] - x[i - 1])
  
  i <- centerindex + 1  # start search for next crossing at center
  while ((sign(y[i] - lev50) == sign(y[i - 1] - lev50)) && (i <= N - 1)) {
    i <- i + 1
  }
  
  if (i != N) {
    Ptype <- 1
    # print('Pulse is Impulse or Rectangular with 2 edges')
    interp <- (lev50 - y[i - 1]) / (y[i] - y[i - 1])
    ttrail <- x[i - 1] + interp * (x[i] - x[i - 1])
    width <- ttrail - tlead
  } else {
    Ptype <- 2
    # print('Step-Like Pulse, no second edge')
    ttrail <- NA
    width <- NA
  }
  
  return(width)
}

# Input file lists

vessel_paths <- "C:/Users/GB Park/OneDrive - University of Maryland School of Medicine/Desktop/Lab/data/After FIJI/Glutamate/Glutamate 200/test/"
fileArray <- list.files(path = vessel_paths)


# Output file path and name
folder_path <- "C:/Users/GB Park/OneDrive - University of Maryland School of Medicine/Desktop/Lab/data/diamater_numbers/"
output_name<-"250228_anta2_PA1_test"

N_files=length(fileArray)
N_frames=250

allData <- list(
  percentChangeTraces = matrix(NA, nrow = N_frames, ncol = N_files),
  percentChangeRAW=matrix(NA, nrow = N_frames, ncol = N_files),
  DiameterTraces=matrix(NA, nrow = N_frames, ncol = N_files),
  DiameterDeltaRaw=matrix(NA, nrow = N_frames, ncol = N_files),
  FractionalChange=matrix(NA, nrow = N_frames, ncol = N_files),
  zscores=matrix(NA, nrow = N_frames, ncol = N_files),
  sd=matrix(NA,nrow=N_frames,ncol=N_files),
  Do=matrix(NA,nrow=N_frames,ncol=N_files),
  max_change_microns= matrix(NA,nrow=1,ncol=N_files),
  max_baseline_microns=matrix(NA,nrow=1,ncol=N_files)
)

# --------- diameter trace maker -----
for (vessel_iter in seq_along(vessel_paths)) {
  for (file_iter in seq_along(fileArray)) {
    filename <- fileArray[file_iter]
    path <- vessel_paths[vessel_iter]


# ---- Where were the stimulations (which frame(s)?)----
StimFrameArray <- 50

# ---------------------------
filepath <- paste0(path, filename)
DATA <- read.csv(filepath)

# Convert data frame to matrix
DATA <- as.matrix(DATA)

# Column where data starts (3rd column)
dataStart <- 3
DATA <- DATA[, dataStart:ncol(DATA)]  # Scrap the first two columns

# Baseline diameter
StimFrameArray <- 50
Baseline_Diameter <- 30:StimFrameArray[1]

# Frame until you want (first 250 columns)
profile <- DATA[, 1:N_frames]

# End frame (number of columns in profile)
End <- ncol(profile)

# Detection algorithm flag
detection_algorithm <- 0

# Imaging parameters
microns_per_pixel <- 0.21  # microns/pixel
time_per_frame <- 0.33     # seconds per frame

# using peak finder on background subtracted, time-grouped, edge-detected, smoothed images
# requires preprocessing (in ImageJ): 3x3 smoothing (Process->smooth),
# Rolling Ball Background Subtraction (~50pixel radius), Blurring Gaussian
# Kernel, then Edge Finder (alpha = 0.1-0.5) or just fwhm and extract line
# profile for every frame. Use the average intensity projection to pick the
# largest diameter 

if (detection_algorithm == 1) {
  ThisProfilePeakLocs <- list()
  junk <- NULL
  x0 <- profile
  
  for (i in 1:ncol(x0)) {
    # Using findpeaks function from the 'pracma' package (alternative to peakfinder in MATLAB)
    peak_data <- pracma::findpeaks(x0[, i])
    ThisProfilePeakLocs <- peak_data
    
    if (nrow(ThisProfilePeakLocs) == 2) { # if it resolves two peaks indicating the Canny Edge Algorithm resolved a diameter 'edge'
      diameter[i] <- microns_per_pixel * (ThisProfilePeakLocs[nrow(ThisProfilePeakLocs), 1] - ThisProfilePeakLocs[1, 1])
    } else {
      diameter[i] <- NA # if it finds three peaks, don't even try
    }
  }
}

# ----FWHM-----
if (detection_algorithm == 0) {
  x <- seq(0, 1 * nrow(profile), length.out = nrow(profile))
  
  FWHM <- numeric(End)  # Initialize an empty vector for FWHM
  
  for (i in 1:End) {
    xq <- seq(1, 1 * nrow(profile), length.out = 1000)
    y <- spline(x, profile[, i], xout = xq)$y  # Interpolation using spline
    
    FWHM[i] <- fwhm(xq, y)  # Assuming `fwhm` is a custom function or a predefined function in your environment
  }
  
  diameter <- microns_per_pixel * FWHM
}


# ------ outliers ------ 
#but not gonna used
# Calculate the standard deviation of the diameter values
std_FWHM_total <- sd(diameter)

# Identify outlier frames based on 3 standard deviations away from the mean
#mean_diameter <- mean(diameter)
#outlier_threshold <- 3 * std_FWHM_total
#outliers <- abs(diameter - mean_diameter) > outlier_threshold

# Replace outlier frames with NA
#diameter[outliers] <- NA

# --------Basal diameter --------
#(using the 20 frames before the first stimulus)
number_of_baseline_frames <- 20
Do_span <- (StimFrameArray[1] - number_of_baseline_frames)
Do <- mean(diameter[Do_span:StimFrameArray[1]], na.rm = TRUE)  # Baseline mean (ignoring NaN values)
sigma_o <- sd(diameter[Do_span:StimFrameArray[1]], na.rm = TRUE)  # Standard deviation of the baseline (ignoring NaN values)

# Store diameter data
dataOut <- diameter

# ----Calculate z-scores, percent change, and other data----
diameter_trace <- t(diameter) 
z_scored <- (diameter - Do) / sigma_o 
percentChange <- 100 * (diameter - Do) / Do  

# Store results in a list (to mimic a struct in MATLAB)
data <- list(
  DiameterAcrossTime = diameter,
  zscores = z_scored,
  percentChange = percentChange,
  Do = Do,
  sigma_o = sigma_o
)

# Baseline z-scored data for frames up to the first stimulus
data$baseline <- z_scored[1:StimFrameArray[1]]

if (length(StimFrameArray) == 1) {
  data$S1 <- z_scored[(-number_of_baseline_frames + StimFrameArray[1] + 1):length(z_scored)]
}
saveRDS(data, file = paste0(path, '/OUTPUT_', filename, '.RData'))


# ------------------------ FIRST SUBPLOT, FWHM, % change, z-score to baseline mean and std
movmeanweight <- 5

# Create a plot
library(ggplot2)
library(zoo)  # For moving average

windows(title=filename) 
par(mfrow = c(2, 1))  # To create a 3x1 grid for subplots

# Calculate moving average
M <- rollapply(data$DiameterAcrossTime, movmeanweight, mean, fill = NA, align = "center")
# Adjust margins
par(mar = c(5, 4, 4, 2))  # Increase margin space if necessary
options(repr.plot.width = 8, repr.plot.height = 6)
# Create the plot - microns
#plot(data$DiameterAcrossTime, type = 'l', col = 'black', lwd = 1, 
#     xlab = 'frame', ylab = 'microns', xlim = c(0,250), cex.lab = 1.4, cex.axis = 1.4)
#lines(M, col = 'blue', lwd = 2)
#abline(h = Do, col = 'black', lty = 2)
#title('profile FWHM')


# First subplot: % change
# Plot percent change
plot(percentChange, type = 'l', col = 'black', lwd = 1, 
     xlab = 'frame', ylab = '% change', xlim = c(0,250), cex.lab = 1.4, cex.axis = 1.4)
M <- filter(percentChange, rep(1/movmeanweight, movmeanweight), sides = 2)  # Moving average (like movmean)
lines(M, col = 'red', lwd = 2)
abline(h = 0, col = 'black', lty = 2)
title('% change')

# Second subplot: z-score
plot(z_scored, type = 'l', col = 'black', lwd = 1, 
     xlab = 'frame', ylab = 'SD', xlim = c(0,250), cex.lab = 1.4, cex.axis = 1.4)
M <- filter(z_scored, rep(1/movmeanweight, movmeanweight), sides = 2)  # Moving average (like movmean)
lines(M, col = 'magenta', lwd = 2)
abline(h = 0, col = 'black', lty = 2)
title('z-score')


#-------------------


# Perform moving average on percentChange, z_scored
allData$percentChangeTraces[,file_iter] <-filter(percentChange, rep(1/movmeanweight, movmeanweight), sides = 2)  # Moving average (like movmean)
allData$percentChangeRAW[,file_iter] <- percentChange
allData$DiameterTraces[, file_iter] <- diameter
allData$DiameterDeltaRaw[, file_iter] <- diameter - Do
allData$FractionalChange[, file_iter] <- diameter / Do
allData$zscores[, file_iter] <- stats::filter(z_scored, rep(1/movmeanweight, movmeanweight), sides = 2)
allData$sd[, file_iter] <- sigma_o
allData$Do[, file_iter] <- Do

allData$max_change_microns[, file_iter] <- max(diameter[20:length(diameter)])

allData$max_baseline_microns[file_iter] <- max(diameter[1:19])

  }
}

# Store as csv file


full_file_path <- file.path(folder_path, paste0(output_name, ".csv"))

# ✅ Write the matrix to the selected file
write.csv(allData$percentChangeTraces, file = full_file_path, row.names = TRUE)
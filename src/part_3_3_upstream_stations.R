# Load data
library(zoo)
data <- read.table("data/data20102011-10min.txt", header = TRUE)
#  Remove all rows with any missing values
clean_data <- na.omit(data)

# Convert datetime column to proper format
clean_data$datetime <- as.POSIXct(clean_data$datetime, format = "%Y-%m-%d %H:%M:%S")
# Extract cleaned time series
os11 <- clean_data$OS11waterlevel
rpbu <- clean_data$RPBUwaterlevel
vr   <- clean_data$VRwaterlevel

# CCF without error
ccf(os11, rpbu, lag.max = 10, main = "CCF: OS11 → RPBU")
ccf(vr, rpbu, lag.max = 30, main = "CCF: VR → RPBU")
ccf(vr, os11, lag.max = 20, main = "CCF: VR → OS11")

# Safe lagging function
lag_series <- function(series, lag) {
  if (lag >= length(series)) return(rep(NA, length(series)))
  c(rep(NA, lag), head(series, -lag))
}

get_r_squared <- function(predictor, response, lag) {
  if (lag >= length(predictor) || lag >= length(response)) return(NA)
  
  # Align both series
  predictor_lagged <- lag_series(predictor, lag)
  response_aligned <- response[1:(length(response) - lag)]
  predictor_lagged <- predictor_lagged[1:(length(predictor_lagged) - lag)]
  
  # Build data frame
  df <- data.frame(y = response_aligned, x = predictor_lagged)
  df <- na.omit(df)
  if (nrow(df) < 2) return(NA)
  
  model <- lm(y ~ x, data = df)
  return(summary(model)$r.squared)
}


lags <- 0:10  # up to 100 minutes

r2_os11_rpbu <- sapply(lags, function(l) {
  r2 <- get_r_squared(clean_data$OS11waterlevel, clean_data$RPBUwaterlevel, l)
  print(paste("Lag", l, "→ R² =", round(r2, 4)))
  return(r2)
})


r2_vr_rpbu <- sapply(lags, function(l) {
  r2 <- get_r_squared(clean_data$VRwaterlevel, clean_data$RPBUwaterlevel, l)
  print(paste("Lag", l, "→ R² (VR → RPBU) =", round(r2, 4)))
  return(r2)
})

r2_vr_os11 <- sapply(lags, function(l) {
  r2 <- get_r_squared(clean_data$VRwaterlevel, clean_data$OS11waterlevel, l)
  print(paste("Lag", l, "→ R² (VR → OS11) =", round(r2, 4)))
  return(r2)
})


# Plot results
plot(lags * 10, r2_os11_rpbu, type = "b", col = "blue", ylim = c(0, 1),
     xlab = "Lag (minutes)", ylab = "R²", main = "Lagged Regression R²")
lines(lags * 10, r2_vr_rpbu, type = "b", col = "red")
lines(lags * 10, r2_vr_os11, type = "b", col = "green")
legend("bottomright", legend = c("OS11 → RPBU", "VR → RPBU", "VR → OS11"),
       col = c("blue", "red", "green"), lty = 1, pch = 1)    

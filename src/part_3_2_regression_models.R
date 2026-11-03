library(zoo)

data <- read.table("data/RPBU20102011.txt", header = TRUE)

data$date <- as.Date(data$datetime)
daily_max <- aggregate(water ~ date, data = data, FUN = max)
high_tide_days <- daily_max[daily_max$water > 230, ]

selected_days <- head(high_tide_days$date, 2)
filtered_data <- subset(data, date %in% selected_days)
day1_data <- subset(filtered_data, date == selected_days[1])
day1_data$datetime <- as.POSIXct(day1_data$datetime, format = "%Y-%m-%d %H:%M:%S")

fit_model <- function(secs_before_high_tide = 1800) {
  # periods of oscillations
  T1 <- 545; T2 <- 205; T3 <- 85
  omega1 <- 2 * pi / T1
  omega2 <- 2 * pi / T2
  omega3 <- 2 * pi / T3
  
  # select time window before high tide 
  high_time <- day1_data$datetime[which(day1_data$water > 230)[1]]
  fit_window <- subset(day1_data, datetime >= high_time - secs_before_high_tide
                       & datetime <= high_time)
  
  # time in seconds since window start
  t0 <- min(fit_window$datetime)
  fit_window$t_sec <- as.numeric(difftime(fit_window$datetime, t0, units = "secs"))
  
  # add sine and cosine terms
  fit_window$t_sq <- fit_window$t_sec^2
  fit_window$sin1 <- sin(omega1 * fit_window$t_sec)
  fit_window$cos1 <- cos(omega1 * fit_window$t_sec)
  fit_window$sin2 <- sin(omega2 * fit_window$t_sec)
  fit_window$cos2 <- cos(omega2 * fit_window$t_sec)
  fit_window$sin3 <- sin(omega3 * fit_window$t_sec)
  fit_window$cos3 <- cos(omega3 * fit_window$t_sec)
  
  model <- lm(water ~ t_sec + t_sq + sin1 + cos1 + sin2 + cos2 + sin3 + cos3, data = fit_window)
  
  # fit linear model
  fit_window$fitted <- predict(model)
  
  # create plots
  plot(fit_window$datetime, fit_window$water/ 100, 
       pch = 16, col = "black", cex = 0.3,  
       xlab = "Time", ylab = "Water Level (m)", 
       main = paste("Model Fit", secs_before_high_tide, "seconds before High Tide"))
  lines(fit_window$datetime, fit_window$fitted / 100, 
        col = "blue", lwd = 2)
  legend("topleft", legend = c("Observed", "Fitted"), 
         col = c("black", "blue"), lty = c(NA, 1), 
         pch = c(16, NA), cex = 0.8)
  
  # residual plot
  model.stdres <- rstandard(model)
  plot(fit_window$fitted / 100, model.stdres, 
       pch = 16, col = "black", cex = 0.3,  
       xlab = "Predicted Water levels (m)", ylab = "Standardised Residuals", 
       main = paste("Residual Plot for Model with", secs_before_high_tide, "secs Time window"))
  abline(h = 0, col = "blue", lty = 2)
  
  # q-q plot
  plot_title_qq <- paste("Normal Q-Q Plot of Residuals (",secs_before_high_tide, "s window)")
  qqnorm(residuals(model), 
         main = plot_title_qq,
         xlab = "Theoretical Quantiles",
         ylab = "Sample Quantiles")
  qqline(residuals(model), col = "red")
  
  return(list(model = model, fit_window = fit_window))
}

fit_model_15min <- fit_model(900)
fit_model_30min <- fit_model(1800)
fit_model_60min <- fit_model(3600)
fit_model_90min <- fit_model(5400)
fit_model_120min <- fit_model(7200)

durbinWatsonTest(fit_model_15min$model)
durbinWatsonTest(fit_model_30min$model)
durbinWatsonTest(fit_model_60min$model)
durbinWatsonTest(fit_model_90min$model)
durbinWatsonTest(fit_model_120min$model)

do_5min_predictions <- function(model_results_list, alpha = 0.05) {
  # periods of oscillations
  T1 <- 545; T2 <- 205; T3 <- 85
  omega1 <- 2 * pi / T1
  omega2 <- 2 * pi / T2
  omega3 <- 2 * pi / T3
  
  model <- model_results_list$model
  fit_window <- model_results_list$fit_window
  
  # get the prediction time in seconds
  last_fit_t_sec <- max(fit_window$t_sec)
  predict_t_sec <- last_fit_t_sec + 300
  
  # get the prediction time in datetime
  high_time <- day1_data$datetime[which(day1_data$water > 230)[1]]
  predict_time <- high_time + 300
  
  data_for_prediction <- data.frame(
    t_sec = predict_t_sec,
    t_sq = predict_t_sec^2,
    sin1 = sin(omega1 * predict_t_sec),
    cos1 = cos(omega1 * predict_t_sec),
    sin2 = sin(omega2 * predict_t_sec),
    cos2 = cos(omega2 * predict_t_sec),
    sin3 = sin(omega3 * predict_t_sec),
    cos3 = cos(omega3 * predict_t_sec)
  )
  
  y_hat <- predict(model, newdata = data_for_prediction)
  
  # compute confidence interval
  X <- model.matrix(model)
  XtX_inv <- solve(t(X) %*% X)
  n <- nrow(X)
  p <- ncol(X)
  df <- n - p
  SSE <- sum(residuals(model)^2)
  sigma_hat_sq <- SSE / df
  x_tilde <- model.matrix(~ t_sec + t_sq + sin1 + cos1 + sin2 + cos2 + sin3 + cos3, data = data_for_prediction)
  pred_var <- sigma_hat_sq * (1 + x_tilde %*% XtX_inv %*% t(x_tilde))
  t_quant <- qt(1 - alpha / 2, df)
  margin <- t_quant * sqrt(pred_var)
  lower <- y_hat - margin
  upper <- y_hat + margin
  
  # find actual water level
  closest_index <- which.min(abs(day1_data$datetime - predict_time))
  actual_water_level <- day1_data$water[closest_index]
  
  # error calculations
  error <- actual_water_level - y_hat
  abs_error <- abs(error)
  
  return(list(
    predicted = y_hat,
    actual = actual_water_level,
    error = error,
    abs_error = abs_error,
    lower_bound = lower,
    upper_bound = upper,
    margin = margin
  ))
}

prediction_15mins <- do_5min_predictions(fit_model_15min)
prediction_30mins <- do_5min_predictions(fit_model_30min)


data <- read.table("data/RPBU20102011.txt", header = TRUE)
data$water <- data$water / 100 

data$date <- as.Date(data$datetime)
daily_max <- aggregate(water ~ date, data = data, FUN = max)
high_tide_days <- daily_max[daily_max$water > 2.30, ]
selected_day <- high_tide_days$date[1]
day1_data <- subset(data, date == selected_day)
day1_data$datetime <- as.POSIXct(day1_data$datetime, format = "%Y-%m-%d %H:%M:%S")

fit_non_linear_model <- function(time_before_high_tide, model_formula, beta0_init, 
                                 beta1_init) {
  high_time <- day1_data$datetime[which(day1_data$water > 2.3)[1]]
  fit_window <- subset(day1_data, datetime >= high_time - time_before_high_tide & datetime <= high_time)
  t0 <- min(fit_window$datetime)
  fit_window$t_sec <- as.numeric(difftime(fit_window$datetime, t0, units = "secs"))
  
  Y0 <- fit_window$water[1]
  beta2_init <- (beta0_init / Y0) - 1
  
  start_vals <- list(
    beta0 = beta0_init, beta1 = beta1_init, beta2 = beta2_init, 
    alpha1 = 0.2, omega1 = 2 * pi / 545, phi1 = 0,
    alpha2 = 0.2, omega2 = 2 * pi / 205, phi2 = 0,
    alpha3 = 0.2, omega3 = 2 * pi / 85, phi3 = 0
  )
  
  # fit the model using nls
  nls_fit <- nlsLM(
    formula = model_formula,
    data = fit_window,
    start = start_vals,
    control = nls.lm.control(maxiter = 500, ftol = 1e-8)
  )
  
  # plot observed vs. fitted values
  fit_window$fitted_nls <- predict(nls_fit)
  
  plot(fit_window$datetime, fit_window$water,
       pch = 16, col = "gray", cex = 0.4,
       main = paste("Non-Linear Model with", time_before_high_tide, "secs Time Window"),
       xlab = "Time", ylab = "Water Level (m)")
  lines(fit_window$datetime, fit_window$fitted_nls, col = "blue", lwd = 2)
  legend("topleft", legend = c("Observed", "Fitted"), 
         col = c("gray", "blue"), pch = c(16, NA), lty = c(NA, 1), cex = 0.8)
  
  # residual plot
  resid_vals <- residuals(nls_fit)
  fitted_vals <- predict(nls_fit)
  std_resid <- resid_vals / sd(resid_vals)
  
  plot(fitted_vals, std_resid, 
       pch = 16, col = "black", cex = 0.3,  
       xlab = "Predicted Water Levels (m)", 
       ylab = "Standardized Residuals", 
       main = paste("Residual Plot for Model with", time_before_high_tide, "secs Time Window"))
  abline(h = 0, col = "blue", lty = 2)
  
  # q-q plot
  plot_title_qq <- paste("Normal Q-Q Plot of Residuals (",time_before_high_tide, "s window)")
  qqnorm(residuals(nls_fit), 
         main = plot_title_qq,
         xlab = "Theoretical Quantiles",
         ylab = "Sample Quantiles")
  qqline(residuals(nls_fit), col = "red")
  
  return(list(model = nls_fit, fit_window = fit_window))
}

model_formula <- water ~ 
  beta0 / (1 + beta2 * exp(-beta1 * t_sec)) +
  alpha1 * sin(omega1 * t_sec + phi1) +
  alpha2 * sin(omega2 * t_sec + phi2) +
  alpha3 * sin(omega3 * t_sec + phi3)

beta0_init <- 3.5
beta1_init <- 0.0005
nlm_15_mins <- fit_non_linear_model(900, model_formula, beta0_init, beta1_init)
nlm_30_mins <- fit_non_linear_model(1800, model_formula, beta0_init, beta1_init)

beta0_init <- 3
beta1_init <- 0.005
nlm_60_mins <- fit_non_linear_model(3600, model_formula, beta0_init, beta1_init)

beta0_init <- 2
beta1_init <- 0.00000000000005
nlm_90_mins <- fit_non_linear_model(5400, model_formula, beta0_init, beta1_init)

summary(nlm_15_mins$model)
summary(nlm_30_mins$model)
summary(nlm_60_mins$model)
summary(nlm_90_mins$model)

durbinWatsonTest(lm(residuals(nlm_15_mins$model) ~ 1))
durbinWatsonTest(lm(residuals(nlm_30_mins$model) ~ 1))
durbinWatsonTest(lm(residuals(nlm_60_mins$model) ~ 1))
durbinWatsonTest(lm(residuals(nlm_90_mins$model) ~ 1))


do_5min_predictions_nonlinear <- function(model_results_list) {
  
  nls_model <- model_results_list$model
  fit_window <- model_results_list$fit_window
  
  # get coefficients
  coefs <- coef(nls_model)
  
  last_fit_t_sec <- max(fit_window$t_sec)
  predict_t_sec <- last_fit_t_sec + 300  # 5 mins ahead
  
  # compute prediction using the fitted formula
  predicted_water_level <- coefs["beta0"] / (1 + coefs["beta2"] * exp(-coefs["beta1"] * predict_t_sec)) +
    coefs["alpha1"] * sin(coefs["omega1"] * predict_t_sec + coefs["phi1"]) +
    coefs["alpha2"] * sin(coefs["omega2"] * predict_t_sec + coefs["phi2"]) +
    coefs["alpha3"] * sin(coefs["omega3"] * predict_t_sec + coefs["phi3"])
  
  # get actual water level from data
  start_datetime <- min(fit_window$datetime)
  predict_time <- start_datetime + predict_t_sec
  closest_index <- which.min(abs(day1_data$datetime - predict_time))
  actual_water_level <- day1_data$water[closest_index]
  
  # error calculations
  error <- actual_water_level - predicted_water_level
  abs_error <- abs(error)
  percentage_error <- 100 * abs_error / actual_water_level
  
  return(list(
    predicted = predicted_water_level,
    error = error,
    abs_error = abs_error,
    percentage_error = percentage_error
  ))
}

nls_prediction_15mins <- do_5min_predictions_nonlinear(nlm_15_mins)
nls_prediction_30mins <- do_5min_predictions_nonlinear(nlm_30_mins)
nls_prediction_60mins <- do_5min_predictions_nonlinear(nlm_60_mins)
nls_prediction_90mins <- do_5min_predictions_nonlinear(nlm_90_mins)

#Q6

library(minpack.lm)

# Define the high water time as the moment with maximum water level on the selected day
high_time <- day1_data$datetime[which.max(day1_data$water)]



# Suppose high_time and day1_data are already defined from earlier
# Set your time window (e.g., 3600 seconds = 60 minutes)
time_before_high_tide <- 3600

# Subset data in that time window
fit_window <- subset(day1_data, datetime >= high_time - time_before_high_tide & datetime <= high_time)

# Extract Y0 for starting value of beta2
Y0 <- fit_window$water[1]
beta0_init <- 2
beta2_init <- (beta0_init / Y0) - 1


# Define the non-linear model formula
model_formula <- water ~ 
  beta0 / (1 + beta2 * exp(-beta1 * t_sec)) +
  alpha1 * sin(omega1 * t_sec + phi1) +
  alpha2 * sin(omega2 * t_sec + phi2) +
  alpha3 * sin(omega3 * t_sec + phi3)

# Set your initial values (starting values)
start_vals <- list(
  beta0 = 2.0, 
  beta1 = 0.005, 
  beta2 = 3.0,
  alpha1 = 0.2, omega1 = 2 * pi / 545, phi1 = 0,
  alpha2 = 0.2, omega2 = 2 * pi / 205, phi2 = 0,
  alpha3 = 0.2, omega3 = 2 * pi / 85,  phi3 = 0
)

# Prepare time window (e.g., 60 minutes)
fit_window <- subset(day1_data, datetime >= high_time - 3600 & datetime <= high_time)
t0 <- min(fit_window$datetime)
fit_window$t_sec <- as.numeric(difftime(fit_window$datetime, t0, units = "secs"))

# Fit the model using nlsLM (Levenberg-Marquardt)
nls_fit <- nlsLM(
  formula = model_formula,
  data = fit_window,
  start = start_vals,
  control = nls.lm.control(maxiter = 500, ftol = 1e-8)
)

# Summary of the fit
summary(nls_fit)

    

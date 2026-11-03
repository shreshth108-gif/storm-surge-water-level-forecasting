library(zoo)

data <- read.table("data/RPBU20102011.txt", header = TRUE)
head(data, 5)

data$date <- as.Date(data$datetime)

daily_max <- aggregate(water ~ date, data = data, FUN = max)

high_tide_days <- daily_max[daily_max$water > 230, ]

selected_days <- head(high_tide_days$date, 2)

filtered_data <- subset(data, date %in% selected_days)

day1_data <- subset(filtered_data, date == selected_days[1])
day2_data <- subset(filtered_data, date == selected_days[2])

day1_data$datetime <- as.POSIXct(day1_data$datetime, format = "%Y-%m-%d %H:%M:%S")
day2_data$datetime <- as.POSIXct(day2_data$datetime, format = "%Y-%m-%d %H:%M:%S")

# Plot for day 1
plot(day1_data$datetime, day1_data$water / 100,
     type = "p",
     pch = 16, cex = 0.3,
     xlab = "Time", ylab = "Water Level (m)",
     main = paste("Water Levels on", selected_days[1]),
     col = "darkgreen",
     xaxt = "n")
time_ticks <- seq(from = as.POSIXct(paste(selected_days[1], "00:00:00")), 
                  to = as.POSIXct(paste(selected_days[1], "23:59:59")), 
                  by = "2 hours")
axis.POSIXct(1, at = time_ticks, format = "%H:%M")
axis(2)
box()
abline(h = 2.3, col = "red", lty = 2)

# Plot for day 2
plot(day2_data$datetime, day2_data$water / 100,
     type = "p",
     pch = 16, cex = 0.3,
     xlab = "Time", ylab = "Water Level (m)",
     main = paste("Water Levels on", selected_days[2]),
     col = "blue",
     xaxt = "n")
time_ticks <- seq(from = as.POSIXct(paste(selected_days[2], "00:00:00")), 
                  to = as.POSIXct(paste(selected_days[2], "23:59:59")), 
                  by = "2 hours")
axis.POSIXct(1, at = time_ticks, format = "%H:%M")
axis(2)
box()
abline(h = 2.3, col = "red", lty = 2)


# Zoomed plot for day 1
peak_time_day1 <- day1_data$datetime[which.max(day1_data$water)]
zoom_day1 <- subset(day1_data, datetime >= peak_time_day1 - 3600 & datetime <= peak_time_day1 + 3600)

plot(zoom_day1$datetime, zoom_day1$water / 100,
     type = "p", pch = 16, cex = 0.4,
     xlab = "Time", ylab = "Water Level (m)",
     main = paste("Zoomed-In High Tide on", selected_days[1]),
     col = "darkgreen", xaxt = "n")
axis.POSIXct(1, at = seq(min(zoom_day1$datetime), max(zoom_day1$datetime), by = "15 mins"), format = "%H:%M")
axis(2)
box()
abline(h = 2.3, col = "red", lty = 2)


# Zoomed plot for day 2
peak_time_day2 <- day2_data$datetime[which.max(day2_data$water)]
zoom_day2 <- subset(day2_data, datetime >= peak_time_day2 - 3600 & datetime <= peak_time_day2 + 7200)

plot(zoom_day2$datetime, zoom_day2$water / 100,
     type = "p", pch = 16, cex = 0.4,
     xlab = "Time", ylab = "Water Level (m)",
     main = paste("Zoomed-In High Tide on", selected_days[2]),
     col = "blue", xaxt = "n")
axis.POSIXct(1, at = seq(min(zoom_day2$datetime), max(zoom_day2$datetime), by = "15 mins"), format = "%H:%M")
axis(2)
box()
abline(h = 2.3, col = "red", lty = 2)

# Smooth water level using a rolling mean (10 min window)
day1_data$water_smooth <- rollmean(day1_data$water, k = 60, fill = NA, align = "center")
day2_data$water_smooth <- rollmean(day2_data$water, k = 60, fill = NA, align = "center")

# Smoothed data plot for day 1
plot(day1_data$datetime, day1_data$water_smooth / 100,
     type = "l", col = "darkgreen", lwd = 2,
     xlab = "Time", ylab = "Smoothed Water Level (m)",
     main = paste("Smoothed Water Levels on", selected_days[1]),
     xaxt = "n")
time_ticks <- seq(from = as.POSIXct(paste(selected_days[1], "00:00:00")),
                  to   = as.POSIXct(paste(selected_days[1], "23:59:59")),
                  by = "2 hours")
axis.POSIXct(1, at = time_ticks, format = "%H:%M")
axis(2)
box()
abline(h = 2.3, col = "red", lty = 2)

# Smoothed data plot for day 2
plot(day2_data$datetime, day2_data$water_smooth / 100,
     type = "l", col = "blue", lwd = 2,
     xlab = "Time", ylab = "Smoothed Water Level (m)",
     main = paste("Smoothed Water Levels on", selected_days[2]),
     xaxt = "n")
time_ticks <- seq(from = as.POSIXct(paste(selected_days[2], "00:00:00")),
                  to   = as.POSIXct(paste(selected_days[2], "23:59:59")),
                  by = "2 hours")
axis.POSIXct(1, at = time_ticks, format = "%H:%M")
axis(2)
box()
abline(h = 2.3, col = "red", lty = 2)


# Zoomed and smoothed plot for day 1
peak_time_day1 <- day1_data$datetime[which.max(day1_data$water_smooth)]
zoom_day1 <- subset(day1_data, datetime >= peak_time_day1 - 3600 & datetime <= peak_time_day1 + 3600)

plot(zoom_day1$datetime, zoom_day1$water_smooth / 100,
     type = "p", pch = 16, cex = 0.4,
     xlab = "Time", ylab = "Water Level (m)",
     main = paste("Zoomed-In Smoothed High Tide on", selected_days[1]),
     col = "darkgreen", xaxt = "n")
axis.POSIXct(1, at = seq(min(zoom_day1$datetime), max(zoom_day1$datetime), by = "15 mins"), format = "%H:%M")
axis(2)
box()
abline(h = 2.3, col = "red", lty = 2)

# Zoomed and smoothed plot for day 2
peak_time_day2 <- day2_data$datetime[which.max(day2_data$water_smooth)]
zoom_day2 <- subset(day2_data, datetime >= peak_time_day2 - 3600 & datetime <= peak_time_day2 + 7200)

plot(zoom_day2$datetime, zoom_day2$water_smooth / 100,
     type = "p", pch = 16, cex = 0.4,
     xlab = "Time", ylab = "Water Level (m)",
     main = paste("Zoomed-In Smoothed High Tide on", selected_days[2]),
     col = "blue", xaxt = "n")
axis.POSIXct(1, at = seq(min(zoom_day2$datetime), max(zoom_day2$datetime), by = "15 mins"), format = "%H:%M")
axis(2)
box()
abline(h = 2.3, col = "red", lty = 2)


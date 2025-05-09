# Load required libraries
library(tidyverse)
library(lubridate)
library(patchwork)

# Read the uploaded CSV file (now in your working directory or specify full path)
setwd("C:\\Users\\snehi\\Documents\\Time Series Analysis\\Assignment4")
data <- read_csv("transformer_data.csv")

# Rename columns
data <- data %>%
  rename(
    Y_t = Y,
    T_a_t = Ta,
    Phi_s_t = S,
    Phi_l_t = I
  )

# Plot directly using the numeric time index
p1 <- ggplot(data, aes(x = time, y = Y_t)) +
  geom_line(color = "red") +
  labs(title = "Transformer Temperature", y = "Y_t (°C)", x = "")

p2 <- ggplot(data, aes(x = time, y = T_a_t)) +
  geom_line(color = "blue") +
  labs(title = "Outdoor Air Temperature", y = "T_a_t (°C)", x = "")

p3 <- ggplot(data, aes(x = time, y = Phi_s_t)) +
  geom_line(color = "orange") +
  labs(title = "Solar Radiation", y = "Phi_s_t (W/m²)", x = "")

p4 <- ggplot(data, aes(x = time, y = Phi_l_t)) +
  geom_line(color = "green") +
  labs(title = "Transformer Load", y = "Phi_l_t (kA)", x = "Time Step")


# Step 4: Combine plots vertically
print((p1 / p2 / p3 / p4) + plot_layout(guides = "collect"))
glimpse(data)


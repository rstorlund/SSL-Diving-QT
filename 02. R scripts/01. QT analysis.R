#####QT Analysis#####
#This file plots QT intervals over time for
#Steller sea lions performing stationary dives
#Aligns the HR profiles from the start of the dive
#and plots them altogether
#then calculates average change in HR over the dive

#####Load libraries#####
library(dplyr)
library(readr)
library(purrr)
library(zoo)
library(tools)
library(ggplot2)

#####Read in the data#####

#list all CSV files
mydir <- '/Users/rstorlund/Library/CloudStorage/OneDrive-UBC/QT/Analysis/SSL-Diving-QT/01. Raw data'
file_list <- list.files(path=mydir, pattern = "\\.csv$", full.names = TRUE)

#read and combine
all_dives <- file_list %>%
  set_names(nm = basename(.)) %>% 
  map_df(~ read_csv(.x) %>%
           mutate(DiveID = file_path_sans_ext(basename(.x))))

#####Modify for consistency#####

#calculate average QT from 3 measurements
#replace excel calculated QT with R calculated QT
all_dives <- all_dives %>%
  rowwise() %>%
  mutate(QT = mean(c(`1`, `2`, `3`), na.rm = TRUE)) %>%
  ungroup()

#normalize time to start of dive = 0 s
#mark start and end time
starts_ends <- all_dives %>%
  group_by(DiveID) %>%
  summarise(
    # find index of first underwater START
    start_idx = {
      x <- which(Behaviour == "underwater" & Status == "START")
      if (length(x) > 0) x[1] else NA_integer_
    },
    start_time = if (!is.na(start_idx)) Time[start_idx] else NA_real_,
    
    # find first underwater STOP *after* that start_idx
    end_time = if (!is.na(start_idx)) {
      stop_idx <- which(Behaviour == "underwater" & Status == "STOP" & row_number() > start_idx)
      if (length(stop_idx) > 0) Time[stop_idx[1]] else NA_real_
    } else NA_real_
  ) %>%
  select(-start_idx) %>%  # optional cleanup
  ungroup()

#join start/end times back to full table
all_dives <- all_dives %>%
  left_join(starts_ends, by = "DiveID") %>%
  mutate(Time_from_underwater_start = Time - start_time)

# ---- Optional: trim to underwater period only (uncomment if desired) ----
# all_dives <- all_dives %>% filter(!is.na(start_time), !is.na(end_time), Time >= start_time, Time <= end_time)

#build end markers for plotting (end time relative to start, and nearest QT at end)
end_markers <- all_dives %>%
  group_by(DiveID) %>%
  summarise(
    start_time = first(start_time),
    end_time   = first(end_time),
    end_time_rel = if_else(!is.na(start_time) & !is.na(end_time), end_time - start_time, NA_real_),
    # start_QT = if (!all(is.na(start_time))) {
    #   # find QT at the row nearest to start_time
    #   idx <- which.min(abs(Time - start_time[1]))
    #   QT[idx]
    # }else NA_real_,
    end_QT = if (!all(is.na(end_time))) {
      # find QT at the row nearest to end_time
      idx <- which.min(abs(Time - end_time[1]))
      QT[idx]
    } else NA_real_
  ) %>%
  ungroup()

print(end_markers)

#all dives faceted
#0 indicates dive start
#vertical dashed line at each dive's end
ggplot(all_dives, aes(x = Time_from_underwater_start, y = QT, color = DiveID)) +
  geom_point() +
  # add end verticals (use inherit.aes = FALSE because this layer uses different mapping)
  geom_vline(
    data = end_markers %>% filter(!is.na(end_time_rel)),
    aes(xintercept = end_time_rel, color = DiveID),
    linetype = "dashed", inherit.aes = FALSE) +
  geom_vline(
    aes(xintercept = 0, color = DiveID),
    linetype = "dashed", inherit.aes = FALSE) +
  labs(x = "Time since underwater start (s)", y = "QT interval (ms)", color = "Dive ID") +
  theme_minimal() +
  facet_wrap(~DiveID)


###looking at RR.ms vs QT
ggplot(all_dives, aes(x = RR.ms, y = QT, color = DiveID)) +
  geom_point() +
  # add end verticals (use inherit.aes = FALSE because this layer uses different mapping)
  geom_vline(
    data = end_markers %>% filter(!is.na(end_time_rel)),
    aes(xintercept = end_time_rel, color = DiveID),
    linetype = "dashed", inherit.aes = FALSE) +
  geom_vline(
    aes(xintercept = 0, color = DiveID),
    linetype = "dashed", inherit.aes = FALSE) +
  labs(x = "Time since underwater start (s)", y = "QT interval (ms)", color = "Dive ID") +
  theme_minimal() +
  facet_wrap(~DiveID)


#####Corrected QT#####
###basic data
ggplot(all_dives, aes(Time_from_underwater_start, QT)) +
  geom_line() +
  labs(y = "QT interval (s)", x = "Time (s)",
       title = "QT interval over time") +
  theme_minimal() +
  facet_wrap(~DiveID)

#calculate corrected QT using Frederica method (more stable with large changes in HR)
df <- all_dives
df <- df %>%
  mutate(QTcF = QT / RR.ms^(1/3))

ggplot(df, aes(Time_from_underwater_start, QTcF)) +
  geom_line() +
  labs(y = "QTc (Fridericia)", x = "Time (s)",
       title = "QTcF over time") +  
  theme_minimal() +
  facet_wrap(~DiveID)

#QT hysteresis plot
ggplot(df, aes(RR.ms, QT)) +
  geom_path() +  # draws lines in order
  labs(x = "RR interval (s)", y = "QT interval (s)",
       title = "QT–RR hysteresis during dive") +
  theme_minimal() +
  facet_wrap(~DiveID)

#####Smooth QT hysteresis for each dive#####

# ---------------------------
# 0️⃣ Ensure numeric columns
# ---------------------------
df <- all_dives %>%
  mutate(
    Time_from_underwater_start = as.numeric(Time_from_underwater_start),
    RR.ms = as.numeric(RR.ms),
    QT = as.numeric(QT)
  )

# ---------------------------
# 1️⃣ Smooth RR and QT per dive
# ---------------------------
df_smooth <- df %>%
  group_by(DiveID) %>%
  arrange(Time_from_underwater_start) %>%
  group_modify(~{
    
    fit_RR <- loess(RR.ms ~ Time_from_underwater_start, data = .x, span = 0.2)
    fit_QT <- loess(QT ~ Time_from_underwater_start, data = .x, span = 0.3)
    
    .x %>%
      mutate(
        RR_smooth = predict(fit_RR, newdata = .x),
        QT_smooth = predict(fit_QT, newdata = .x)
      )
  }) %>%
  ungroup()

# ---------------------------
# 2️⃣ Interpolate small NAs to avoid breaks
# ---------------------------
df_smooth <- df_smooth %>%
  group_by(DiveID) %>%
  arrange(Time_from_underwater_start) %>%
  mutate(
    RR_smooth = na.approx(RR_smooth, Time_from_underwater_start, na.rm = FALSE),
    QT_smooth = na.approx(QT_smooth, Time_from_underwater_start, na.rm = FALSE)
  ) %>%
  ungroup()

# ---------------------------
# 3️⃣ Optional downsampling for hysteresis plot (1 s bins)
# ---------------------------
df_bins <- df_smooth %>%
  group_by(DiveID) %>%
  mutate(t_bin = floor(Time_from_underwater_start)) %>%
  group_by(DiveID, t_bin) %>%
  summarise(
    RR_bin = median(RR_smooth, na.rm = TRUE),
    QT_bin = median(QT_smooth, na.rm = TRUE),
    time = median(Time_from_underwater_start, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  filter(!is.na(RR_bin) & !is.na(QT_bin)) %>%
  arrange(DiveID, t_bin)

# ---------------------------
# 4️⃣ Plot QT–RR hysteresis loops
# ---------------------------
ggplot(df_bins, aes(x = RR_bin, y = QT_bin, color = time)) +
  geom_path(size = 1) +
  facet_wrap(~DiveID, scales = "free") +
  scale_color_viridis_c(option = "C") +
  labs(
    x = "RR interval (ms)",
    y = "QT interval (ms)",
    color = "Time (s)",
    title = "QT–RR Hysteresis Loops for All Dives"
  ) +
  theme_minimal() +
  theme(strip.text = element_text(face = "bold"))

# ---------------------------
# 5️⃣ QA/QC: Raw vs Smoothed per dive
# ---------------------------
ggplot(df_smooth %>% filter(DiveID == "Hazy 969 125"), aes(Time_from_underwater_start)) +
  geom_line(aes(y = RR.ms), color = "grey70") +
  geom_line(aes(y = RR_smooth), color = "blue") +
  geom_line(aes(y = QT), color = "grey50") +
  geom_line(aes(y = QT_smooth), color = "red") +
  labs(title = "QA/QC: Raw vs Smoothed (Dive1)", y = "ms")

# ---------------------------
# 6️⃣ QA/QC checks
# ---------------------------

# Check for missing smoothed values
na_summary <- df_smooth %>%
  group_by(DiveID) %>%
  summarise(
    n_RR_NA = sum(is.na(RR_smooth)),
    n_QT_NA = sum(is.na(QT_smooth)),
    n_total = n()
  )

print(na_summary)

# Detect extreme values using z-score
df_smooth <- df_smooth %>%
  group_by(DiveID) %>%
  mutate(
    RR_z = (RR_smooth - mean(RR_smooth, na.rm = TRUE))/sd(RR_smooth, na.rm = TRUE),
    QT_z = (QT_smooth - mean(QT_smooth, na.rm = TRUE))/sd(QT_smooth, na.rm = TRUE)
  ) %>%
  ungroup()

outliers <- df_smooth %>%
  filter(abs(RR_z) > 4 | abs(QT_z) > 4)

print(outliers)

# Check slopes (dRR/dt, dQT/dt)
df_smooth <- df_smooth %>%
  group_by(DiveID) %>%
  arrange(Time_from_underwater_start) %>%
  mutate(
    dRR = c(NA, diff(RR_smooth)),
    dQT = c(NA, diff(QT_smooth))
  ) %>%
  ungroup()










#all dives together

#all dives
ggplot(all_dives, aes(x = Time_from_underwater_start, y = QT, color = DiveID)) +
  geom_line() +  # raw QT traces
  geom_smooth(aes(group = DiveID), method = "lm", se = FALSE, linetype = "dashed") +  # linear fit per dive
  geom_vline(
    data = end_markers %>% filter(!is.na(end_time_rel)),
    aes(xintercept = end_time_rel, color = DiveID),
    linetype = "dotted", inherit.aes = FALSE
  ) +
  labs(
    x = "Time since underwater start (s)",
    y = "QT interval (ms)",
    color = "Dive ID"
  ) +
  theme_minimal()
# Notes:
#   geom_smooth(aes(group = DiveID), method = "lm", se = FALSE) fits a separate linear regression for each dive.
# 
# se = FALSE hides confidence intervals; you can set it to TRUE if you want shaded uncertainty.
# 
# linetype = "dashed" makes the fitted lines visually distinct from raw traces.
# 
# End markers (geom_vline) remain to show the underwater stop time.

#faceted plot
ggplot(all_dives, aes(x = Time_from_underwater_start, y = QT)) +
  geom_line(color = "steelblue") +                           # raw QT trace
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed", color = "red") +  # linear trend
  geom_vline(
    data = end_markers %>% filter(!is.na(end_time_rel)),
    aes(xintercept = end_time_rel),
    linetype = "dotted", color = "black"
  ) +
  facet_wrap(~ DiveID, scales = "free_x") +                 # one panel per dive
  labs(
    x = "Time since underwater start (s)",
    y = "QT interval (ms)"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

#with stats
# For each DiveID, fit linear model QT ~ Time_from_underwater_start
fit_stats <- all_dives %>%
  group_by(DiveID) %>%
  summarise(
    fit = list(lm(QT ~ Time_from_underwater_start, data = cur_data())),
    .groups = "drop"
  ) %>%
  rowwise() %>%
  mutate(
    slope = coef(fit)[2],                              # linear slope
    intercept = coef(fit)[1],                          # intercept (optional)
    r_squared = summary(fit)$r.squared,                # R²
    label_text = paste0("slope = ", round(slope, 3),
                        "\nR² = ", round(r_squared, 3))
  ) %>%
  ungroup()

#Step 2: Add labels to the plot
#We can place the labels in the top-right corner of each facet. Using geom_text() with x and y coordinates positioned dynamically (e.g., near the max x and max QT per dive):


# Determine label positions per dive
label_positions <- all_dives %>%
  group_by(DiveID) %>%
  summarise(
    x = max(Time_from_underwater_start, na.rm = TRUE)*0.8,  # 80% along x-axis
    y = max(QT, na.rm = TRUE)*0.95                           # 95% of max QT
  ) %>%
  left_join(fit_stats %>% select(DiveID, label_text), by = "DiveID")

# Plot
ggplot(all_dives, aes(x = Time_from_underwater_start, y = QT)) +
  geom_line(color = "steelblue") +
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed", color = "red") +
  geom_vline(
    data = end_markers %>% filter(!is.na(end_time_rel)),
    aes(xintercept = end_time_rel),
    linetype = "dotted", color = "black"
  ) +
  geom_text(
    data = label_positions,
    aes(x = x, y = y, label = label_text),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = 1,
    size = 3
  ) +
  facet_wrap(~ DiveID, scales = "free_x") +
  labs(
    x = "Time since underwater start (s)",
    y = "QT interval (ms)"
  ) +
  theme_minimal() +
  theme(legend.position = "none")






#####Plot HR over time#####
QT1 <- ggplot(df, aes(x=Time, y=QTav)) +
  geom_point()
QT1
#####Summary stats#####
#####Statistics######









# ---- Optional: trim to underwater period only (uncomment if desired) ----
# all_dives <- all_dives %>% filter(!is.na(start_time), !is.na(end_time), Time >= start_time, Time <= end_time)

# # ---- 6) Build end markers for plotting (end time relative to start, and nearest QT at end) ----
# end_markers <- all_dives %>%
#   group_by(DiveID) %>%
#   summarise(
#     start_time = first(start_time),
#     end_time   = first(end_time),
#     end_time_rel = if_else(!is.na(start_time) & !is.na(end_time), end_time - start_time, NA_real_),
#     # end_QT = if (!all(is.na(end_time))) {
#     #   # find QT at the row nearest to end_time
#     #   idx <- which.min(abs(Time - end_time[1]))
#     #   QT[idx]
#     # } else NA_real_
#     end_QT = if (!all(is.na(end_time))) {
#       # Ensure Time and QT are numeric vectors of the same length
#       if(length(Time) != length(QT)) stop("Time and QT must be the same length")
#       
#       # Find the index of the Time closest to end_time
#       idx <- which.min(abs(Time - end_time[1]))
#       
#       # Return the corresponding QT value
#       QT[idx]
#     } else {
#       NA_real_
#   ) %>%
#   ungroup()


# ---- 6) Build end markers for plotting (end time relative to start, and nearest QT at end) ----
end_markers <- all_dives %>%
  group_by(Individual, DiveID) %>%
  summarise(
    start_time = first(start_time),
    end_time   = first(end_time),
    end_time_rel = if_else(!is.na(start_time) & !is.na(end_time), end_time - start_time, NA_real_),
    end_QT = if (!all(is.na(end_time))) {
      idx <- which.min(abs(Time - end_time))
      QT[idx]
    } else NA_real_,
    .groups = "drop"
  )

print(end_markers)

# ---- 7) Example plot: overlay with a vertical dashed line at each dive's end ----
ggplot(all_dives, aes(x = Time_from_underwater_start, y = QT, color = DiveID)) +
  geom_point() +
  # add end verticals only for the corresponding individual
  geom_vline(
    data = end_markers %>% filter(!is.na(end_time_rel)),
    aes(xintercept = end_time_rel, color = DiveID),
    linetype = "dashed",
    inherit.aes = FALSE
  ) +
  labs(x = "Time since underwater start (s)", y = "QT interval (ms)", color = "Dive ID") +
  theme_minimal() +
  facet_wrap(~Individual)

#Plot linear model by individual
ggplot(all_dives, aes(x = Time_from_underwater_start, y = QT, color = DiveID)) +
  geom_line() +
  # add LM line per individual (using all dives in that facet)
  geom_smooth(
    aes(group = 1),      # group = 1 fits a single LM per facet
    method = "lm",
    color = "black",     # optional: LM line color
    se = FALSE           # remove confidence interval shading
  ) +
  # add end verticals only for the corresponding individual
  geom_vline(
    data = end_markers %>% filter(!is.na(end_time_rel)),
    aes(xintercept = end_time_rel, color = DiveID),
    linetype = "dashed",
    inherit.aes = FALSE
  ) +
  labs(
    x = "Time since underwater start (s)", 
    y = "QT interval (ms)", 
    color = "Dive ID"
  ) +
  theme_minimal() +
  facet_wrap(~Individual)
#linear model does fit the data very well. Sitka shows decrease in QT then increase
#Yasha too
#Hazy okay with LM

#Try fitting a quadratic curve
ggplot(all_dives, aes(x = Time_from_underwater_start, y = QT, color = DiveID)) +
  
  # lines for each dive
  geom_line() +
  
  # quadratic trend per individual (one curve per facet)
  geom_smooth(
    aes(group = 1),            # ignore DiveID, fit per individual/facet
    method = "lm",
    formula = y ~ poly(x, 2),  # quadratic
    color = "black",
    se = FALSE
  ) +
  
  # add end vertical lines for each dive, only for the correct individual
  geom_vline(
    data = end_markers %>% filter(!is.na(end_time_rel)),
    aes(xintercept = end_time_rel, color = DiveID),
    linetype = "dashed",
    inherit.aes = FALSE
  ) +
  
  # labels and theme
  labs(
    x = "Time since underwater start (s)",
    y = "QT interval (ms)",
    color = "Dive ID"
  ) +
  theme_minimal() +
  
  # facet by individual
  facet_wrap(~Individual)

#Check fit of quatratic model
# Fit a quadratic LM per individual
lm_fit <- lm(QT ~ poly(Time_from_underwater_start, 2), data = all_dives %>% filter(Individual == "Ind1"))

# Plot data + fitted values
all_dives %>%
  filter(Individual == "Ind1") %>%
  mutate(fitted_QT = predict(lm_fit, newdata = .)) %>%
  ggplot(aes(x = Time_from_underwater_start)) +
  geom_point(aes(y = QT)) +
  geom_line(aes(y = fitted_QT), color = "red")

all_dives %>%
  group_by(Individual) %>%
  filter(n_distinct(Time_from_underwater_start) > 2) %>%  # only individuals with ≥3 unique points
  ggplot(aes(x = Time_from_underwater_start, y = QT, color = DiveID)) +
  geom_line() +
  geom_smooth(aes(group = 1), method = "lm", formula = y ~ poly(x, 2), color = "black", se = FALSE) +
  facet_wrap(~Individual)

geom_smooth(aes(group = 1),
            method = "lm",
            formula = y ~ ifelse(n_distinct(Time_from_underwater_start) > 2,
                                 poly(x, 2),
                                 x),
            color = "black",
            se = FALSE)

#####TESTING GAM#####
library(ggplot2)
library(dplyr)
library(mgcv)  # needed for GAMs

# ---- 1) Make sure end_markers includes Individual ----
end_markers <- all_dives %>%
  group_by(Individual, DiveID) %>%
  summarise(
    start_time = first(start_time),
    end_time   = first(end_time),
    end_time_rel = if_else(!is.na(start_time) & !is.na(end_time),
                           end_time - start_time,
                           NA_real_),
    end_QT = if (!all(is.na(end_time))) {
      idx <- which.min(abs(Time - end_time))
      QT[idx]
    } else NA_real_,
    .groups = "drop"
  )

# ---- 2) Plot with GAM fit ----
ggplot(all_dives, aes(x = Time_from_underwater_start, y = QT, color = DiveID)) +
  
  # lines for each dive
  geom_line() +
  
  # GAM smooth curve per individual (one per facet)
  geom_smooth(
    aes(group = 1),        # group = 1 ensures one smooth per individual/facet
    method = "gam",
    formula = y ~ s(x, k = 3),  # spline with k=3 basis functions
    color = "black",
    se = FALSE
  ) +
  
  # add end vertical lines per dive, only for the correct individual
  geom_vline(
    data = end_markers %>% filter(!is.na(end_time_rel)),
    aes(xintercept = end_time_rel, color = DiveID),
    linetype = "dashed",
    inherit.aes = FALSE
  ) +
  
  # labels and theme
  labs(
    x = "Time since underwater start (s)",
    y = "QT interval (ms)",
    color = "Dive ID"
  ) +
  theme_minimal() +
  
  # facet by individual
  facet_wrap(~Individual)

#Assess model fit by individual
library(mgcv)

# Example for one individual
ind_data <- all_dives %>% filter(Individual == "Hazy")

gam_fit <- gam(QT ~ s(Time_from_underwater_start, k = 3), data = ind_data)

summary(gam_fit)

# Family: gaussian 
# Link function: identity 
# 
# Formula:
#   QT ~ s(Time_from_underwater_start, k = 3)
# 
# Parametric coefficients:
#   Estimate Std. Error t value Pr(>|t|)    
# (Intercept) 286.8393     0.5639   508.7   <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# Approximate significance of smooth terms:
#   edf Ref.df     F p-value    
# s(Time_from_underwater_start)   1      1 231.4  <2e-16 ***
#   ---
#   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
# 
# R-sq.(adj) =  0.411   Deviance explained = 41.3%
# GCV =  105.9  Scale est. = 105.26    n = 331

#Plot fitted curve vs data
ind_data <- ind_data %>%
  mutate(fitted_QT = predict(gam_fit, newdata = ind_data),
         resid = QT - fitted_QT)

library(ggplot2)

# Fitted vs observed
ggplot(ind_data, aes(x = Time_from_underwater_start)) +
  geom_point(aes(y = QT), color = "blue") +
  geom_line(aes(y = fitted_QT), color = "red", size = 1) +
  labs(y = "QT (ms)", title = "GAM fit for Individual 1")

#Plot residuals
ggplot(ind_data, aes(x = fitted_QT, y = resid)) +
  geom_point() +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(x = "Fitted QT", y = "Residuals") +
  theme_minimal()

#diagnostic plots
par(mfrow = c(2,2))
gam.check(gam_fit)

####GAM FOR ALL####
library(dplyr)
library(ggplot2)
library(mgcv)

# List of individuals
individuals <- c("Hazy", "Sitka", "Yasha")

# Function to fit GAM and produce diagnostics
fit_gam_diagnostics <- function(ind_name, data) {
  
  ind_data <- data %>% filter(Individual == ind_name)
  
  # Skip individuals with too few points
  if(nrow(ind_data) < 4) {
    message("Skipping ", ind_name, ": not enough points to fit GAM")
    return(NULL)
  }
  
  # Fit GAM (spline with 3 basis functions)
  gam_fit <- gam(QT ~ s(Time_from_underwater_start, k = 3), data = ind_data)
  
  # Extract R-squared and deviance explained
  r2 <- summary(gam_fit)$r.sq
  dev <- summary(gam_fit)$dev.expl
  
  cat("\n---", ind_name, "---\n")
  cat("R-squared:", round(r2, 3), "\n")
  cat("Deviance explained:", round(dev, 3), "\n")
  
  # Add fitted values and residuals to data
  ind_data <- ind_data %>%
    mutate(fitted_QT = predict(gam_fit, newdata = ind_data),
           resid = QT - fitted_QT)
  
  # Plot observed vs fitted
  p1 <- ggplot(ind_data, aes(x = Time_from_underwater_start)) +
    geom_point(aes(y = QT), color = "blue") +
    geom_line(aes(y = fitted_QT), color = "red", size = 1) +
    labs(y = "QT (ms)", title = paste("Observed vs Fitted -", ind_name)) +
    theme_minimal()
  
  # Plot residuals vs fitted
  p2 <- ggplot(ind_data, aes(x = fitted_QT, y = resid)) +
    geom_point() +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(x = "Fitted QT", y = "Residuals", title = paste("Residuals -", ind_name)) +
    theme_minimal()
  
  list(gam_fit = gam_fit, obs_vs_fit = p1, resid_plot = p2)
}

# Run for all individuals
gam_results <- lapply(individuals, fit_gam_diagnostics, data = all_dives)

# Example: show plots for Hazy
gam_results[[3]]$obs_vs_fit
gam_results[[3]]$resid_plot

#####QUADRATIC MODEL FOR ALL#####
library(dplyr)
library(ggplot2)

# ---- 1) Build end markers per dive ----
end_markers <- all_dives %>%
  group_by(Individual, DiveID) %>%
  summarise(
    start_time = first(start_time),
    end_time   = first(end_time),
    end_time_rel = if_else(!is.na(start_time) & !is.na(end_time),
                           end_time - start_time,
                           NA_real_),
    end_QT = if (!all(is.na(end_time))) {
      idx <- which.min(abs(Time - end_time))
      QT[idx]
    } else NA_real_,
    .groups = "drop"
  )

# ---- 2) Plot QT over time with quadratic trend ----
p <- ggplot(all_dives, aes(x = Time_from_underwater_start, y = QT, color = DiveID)) +
  
  geom_line() +
  
  # Quadratic smooth per individual
  geom_smooth(
    aes(group = 1),
    method = "lm",
    formula = y ~ poly(x, 2),
    color = "black",
    se = FALSE
  ) +
  
  # Dive end markers
  geom_vline(
    data = end_markers %>% filter(!is.na(end_time_rel)),
    aes(xintercept = end_time_rel, color = DiveID),
    linetype = "dashed",
    inherit.aes = FALSE
  ) +
  
  labs(
    x = "Time since underwater start (s)",
    y = "QT interval (ms)",
    color = "Dive ID"
  ) +
  theme_minimal() +
  facet_wrap(~Individual)

# Print plot
print(p)

# ---- 3) Fit quadratic models per individual and check fit ----
individuals <- unique(all_dives$Individual)

quad_results <- lapply(individuals, function(ind_name) {
  
  ind_data <- all_dives %>% filter(Individual == ind_name)
  
  # Skip if not enough points
  if(nrow(ind_data) < 3) {
    message("Skipping ", ind_name, ": too few points for quadratic")
    return(NULL)
  }
  
  # Fit quadratic model
  quad_fit <- lm(QT ~ poly(Time_from_underwater_start, 2), data = ind_data)
  
  # Extract summary
  fit_summary <- summary(quad_fit)
  r2 <- fit_summary$r.squared
  adj_r2 <- fit_summary$adj.r.squared
  
  cat("\n--- Individual:", ind_name, "---\n")
  cat("R-squared:", round(r2, 3), "\n")
  cat("Adjusted R-squared:", round(adj_r2, 3), "\n")
  
  # Add fitted values and residuals
  ind_data <- ind_data %>%
    mutate(fitted_QT = predict(quad_fit, newdata = ind_data),
           resid = QT - fitted_QT)
  
  # Observed vs fitted plot
  obs_vs_fit <- ggplot(ind_data, aes(x = Time_from_underwater_start)) +
    geom_point(aes(y = QT), color = "blue") +
    geom_line(aes(y = fitted_QT), color = "red", size = 1) +
    labs(y = "QT (ms)", title = paste("Observed vs Fitted -", ind_name)) +
    theme_minimal()
  
  # Residuals vs fitted plot
  resid_plot <- ggplot(ind_data, aes(x = fitted_QT, y = resid)) +
    geom_point() +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(x = "Fitted QT", y = "Residuals", title = paste("Residuals -", ind_name)) +
    theme_minimal()
  
  list(fit = quad_fit, obs_vs_fit = obs_vs_fit, resid_plot = resid_plot)
})

# Show plots for Hazy
quad_results[[which(individuals == "Hazy")]]$obs_vs_fit
quad_results[[which(individuals == "Hazy")]]$resid_plot

# Show plots for Sitka
quad_results[[which(individuals == "Sitka")]]$obs_vs_fit
quad_results[[which(individuals == "Sitka")]]$resid_plot

# Show plots for Yasha
quad_results[[which(individuals == "Yasha")]]$obs_vs_fit
quad_results[[which(individuals == "Yasha")]]$resid_plot

#note left skew of residuals vs fitted
#persists in GAM, therefofe datq driven

#check for cause
ggplot(all_dives, aes(Time_from_underwater_start)) +
  geom_histogram(bins = 30) +
  facet_wrap(~Individual)

#histograms are U shaped indicating high sampling freq early and late in dive
#skew is explained
#don't do anything to try and deal with it, just a result of the way data was
#fewer QT values in the middle because fewer heart beats
#actually makes physiological sense

geom_smooth(
  aes(group = DiveID),
  method = "lm",
  formula = y ~ poly(x, 2),
  se = FALSE
)

#####Try plotting quadratic vs diveID (12 models) instead of Individual#####
dive_ids <- unique(all_dives$DiveID)

quad_dive_results <- lapply(dive_ids, function(dive) {
  
  dive_data <- all_dives %>% filter(DiveID == dive)
  
  # Safety check
  if (n_distinct(dive_data$Time_from_underwater_start) < 3) {
    message("Skipping dive ", dive, ": not enough points for quadratic")
    return(NULL)
  }
  
  # Fit quadratic model
  quad_fit <- lm(
    QT ~ poly(Time_from_underwater_start, 2),
    data = dive_data
  )
  
  fit_sum <- summary(quad_fit)
  
  cat("\n--- Dive:", dive, "---\n")
  cat("Individual:", unique(dive_data$Individual), "\n")
  cat("R-squared:", round(fit_sum$r.squared, 3), "\n")
  cat("Adjusted R-squared:", round(fit_sum$adj.r.squared, 3), "\n")
  
  # Add fitted values + residuals
  dive_data <- dive_data %>%
    mutate(
      fitted_QT = predict(quad_fit, newdata = dive_data),
      resid = QT - fitted_QT
    )
  
  # Observed vs fitted
  p_fit <- ggplot(dive_data, aes(Time_from_underwater_start)) +
    geom_point(aes(y = QT), color = "blue") +
    geom_line(aes(y = fitted_QT), color = "red", linewidth = 1) +
    labs(
      title = paste("Observed vs Fitted – Dive", dive),
      y = "QT (ms)"
    ) +
    theme_minimal()
  
  # Residuals vs fitted
  p_resid <- ggplot(dive_data, aes(fitted_QT, resid)) +
    geom_point() +
    geom_hline(yintercept = 0, linetype = "dashed") +
    labs(
      title = paste("Residuals vs Fitted – Dive", dive),
      x = "Fitted QT",
      y = "Residuals"
    ) +
    theme_minimal()
  
  list(
    fit = quad_fit,
    obs_vs_fit = p_fit,
    resid_plot = p_resid
  )
})

#Plot
library(dplyr)
library(ggplot2)

ggplot(all_dives,
       aes(x = Time_from_underwater_start, y = QT, color = DiveID)) +
  
  geom_line() +
  
  geom_smooth(
    aes(group = DiveID),
    method = "lm",
    formula = y ~ poly(x, 2),
    se = FALSE
  ) +
  
  labs(
    x = "Time since underwater start (s)",
    y = "QT interval (ms)",
    color = "Dive ID"
  ) +
  theme_minimal() +
  facet_wrap(~DiveID, scales = "free_x")

#observed vs fitted
quad_fitted_all <- all_dives %>%
  group_by(DiveID) %>%
  filter(n_distinct(Time_from_underwater_start) >= 3) %>%
  do({
    fit <- lm(QT ~ poly(Time_from_underwater_start, 2), data = .)
    mutate(., fitted_QT = predict(fit, newdata = .))
  }) %>%
  ungroup()

ggplot(quad_fitted_all,
       aes(x = Time_from_underwater_start)) +
  
  geom_point(aes(y = QT), alpha = 0.7) +
  geom_line(aes(y = fitted_QT), color = "red", linewidth = 1) +
  
  labs(
    x = "Time since underwater start (s)",
    y = "QT (ms)"
  ) +
  theme_minimal() +
  facet_wrap(~DiveID, scales = "free_x")

#residuals vs fitted
quad_fitted_all <- quad_fitted_all %>%
  mutate(resid = QT - fitted_QT)

ggplot(quad_fitted_all,
       aes(x = fitted_QT, y = resid)) +
  
  geom_point(alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  
  labs(
    x = "Fitted QT",
    y = "Residuals"
  ) +
  theme_minimal() +
  facet_wrap(~DiveID)

#####CORRECT METHOD EXTRACTING FEATURES FROM EACH DIVE AND MODELLING#####

###Extract metrics
library(dplyr)

dive_summary <- all_dives %>%
  group_by(DiveID, Individual) %>%
  summarise(
    
    Dive_duration = max(Time_from_underwater_start, na.rm = TRUE),
    
    # QT at start: median of first 10% of dive
    QT_start = median(
      QT[Time_from_underwater_start <= 0.1 * Dive_duration],
      na.rm = TRUE
    ),
    
    # Minimum QT
    QT_min = min(QT, na.rm = TRUE),
    
    # Time to minimum QT
    Time_to_QT_min = Time_from_underwater_start[which.min(QT)],
    
    # QT at end: median of last 10% of dive
    QT_end = median(
      QT[Time_from_underwater_start >= 0.9 * Dive_duration],
      na.rm = TRUE
    ),

    .groups = "drop"
  ) %>%
  mutate(
    Delta_QT = QT_start - QT_min,
    QT_recovered = QT_end - QT_min,
    Prop_QT_recovery = (QT_end - QT_min) / (QT_start - QT_min)
  )

###Figure 1
#QT vs time
#one panel per individual
#multiple dives per panel
#NO models emphasized

###Data no models
ggplot(all_dives, aes(x = Time_from_underwater_start, y = QT, color = DiveID)) +
  
  # lines for each dive
  geom_line() +
  
  # add end vertical lines for each dive, only for the correct individual
  geom_vline(
    data = end_markers %>% filter(!is.na(end_time_rel)),
    aes(xintercept = end_time_rel, color = DiveID),
    linetype = "dashed",
    inherit.aes = FALSE
  ) +
  
  # labels and theme
  labs(
    x = "Time since underwater start (s)",
    y = "QT interval (ms)",
    color = "Dive ID"
  ) +
  theme_minimal() +
  
  # facet by individual
  facet_wrap(~Individual)

###Quadratic model
#Try fitting a quatratic curve
ggplot(all_dives, aes(x = Time_from_underwater_start, y = QT, color = DiveID)) +
  
  # lines for each dive
  geom_line() +
  
  # quadratic trend per individual (one curve per facet)
  geom_smooth(
    aes(group = 1),            # ignore DiveID, fit per individual/facet
    method = "lm",
    formula = y ~ poly(x, 2),  # quadratic
    color = "black",
    se = FALSE
  ) +
  
  # add end vertical lines for each dive, only for the correct individual
  geom_vline(
    data = end_markers %>% filter(!is.na(end_time_rel)),
    aes(xintercept = end_time_rel, color = DiveID),
    linetype = "dashed",
    inherit.aes = FALSE
  ) +
  
  # labels and theme
  labs(
    x = "Time since underwater start (s)",
    y = "QT interval (ms)",
    color = "Dive ID"
  ) +
  theme_minimal() +
  
  # facet by individual
  facet_wrap(~Individual)

###Evaluate QT curvature over dive
curvature_summary <- all_dives %>%
  group_by(DiveID, Individual) %>%
  do({
    fit <- lm(QT ~ poly(Time_from_underwater_start, 2), data = .)
    tibble(
      curvature = coef(fit)[3]  # quadratic term
    )
  })

ggplot(curvature_summary,
       aes(x = Individual, y = curvature)) +
  geom_point(size = 3) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(y = "QT curvature (quadratic term)") +
  theme_minimal()

###Timing phenotype
ggplot(dive_summary,
       aes(x = Individual, y = Time_to_QT_min / Dive_duration)) +
  geom_point(size = 3) +
  labs(y = "Relative timing of minimum QT") +
  theme_minimal()

###Look at QT recovery

#Plot QT recovery
library(ggplot2)

ggplot(dive_summary,
       aes(x = Dive_duration, y = Prop_QT_recovery, color = Individual)) +
  geom_point(size = 3) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  geom_smooth(method = "lm", se = FALSE) +
  labs(
    x = "Dive duration (s)",
    y = "Proportional QT recovery"
  ) +
  theme_minimal()
#doesn't seem useful

#Absolute recovery
ggplot(dive_summary,
       aes(x = Dive_duration, y = QT_recovered, color = Individual)) +
  geom_point(size = 3) +
  geom_smooth(method = "lm", se = FALSE) +
  labs(
    x = "Dive duration (s)",
    y = "QT recovered (ms)"
  ) +
  theme_minimal()
#also seems useless






#####comparing QT duration with RR interval#####
ggplot(data = all_dives, aes(x = RR.ms, y = QT, colour = Individual)) +
  geom_point()

###Add linear model
ggplot(all_dives, aes(x = RR.ms, y = QT, colour = Individual)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, colour = "black") +
  labs(x = "RR interval (ms)", y = "QT interval (ms)")

###Try binning data based on common heart rates (RR intervals)
library(dplyr)

all_dives_summary <- all_dives %>%
  mutate(RR_group = case_when(
    RR.ms < 1000 ~ "Fast HR",
    RR.ms < 2000 ~ "Medium HR",
    TRUE ~ "Slow HR"
  )) %>%
  group_by(RR_group) %>%
  summarise(
    mean_QT = mean(QT, na.rm = TRUE),
    sd_QT = sd(QT, na.rm = TRUE),
    n = n()
  )


ggplot(all_dives_summary, aes(x = RR_group, y = mean_QT)) +
  geom_col(fill = "skyblue", width = 0.5) +
  geom_errorbar(aes(ymin = mean_QT - sd_QT, ymax = mean_QT + sd_QT), width = 0.2)

library(ggplot2)
library(dplyr)

# First, summarize by RR group
all_dives_summary <- all_dives %>%
  mutate(RR_group = case_when(
    RR.ms < 1000 ~ "Fast HR",
    RR.ms < 2000 ~ "Medium HR",
    TRUE ~ "Slow HR"
  )) %>%
  group_by(RR_group) %>%
  summarise(
    mean_QT = mean(QT, na.rm = TRUE),
    sd_QT = sd(QT, na.rm = TRUE),
    .groups = "drop"
  )

# Now plot
ggplot(all_dives, aes(x = RR.ms, y = QT, colour = Individual)) +
  geom_point(size = 3, alpha = 0.7) +  # raw points
  # overlay mean ± SD as error bars
  geom_errorbar(
    data = all_dives_summary,
    aes(x = case_when(
      RR_group == "Fast HR" ~ 500,
      RR_group == "Medium HR" ~ 1275,
      RR_group == "Slow HR" ~ 2500
    ),
    ymin = mean_QT - sd_QT,
    ymax = mean_QT + sd_QT
    ),
    inherit.aes = FALSE,
    width = 100,
    colour = "black",
    size = 1
  ) +
  geom_point(
    data = all_dives_summary,
    aes(x = case_when(
      RR_group == "Fast HR" ~ 500,
      RR_group == "Medium HR" ~ 1275,
      RR_group == "Slow HR" ~ 2500
    ),
    y = mean_QT
    ),
    inherit.aes = FALSE,
    colour = "black",
    size = 4
  ) +
  labs(
    x = "RR interval (ms)",
    y = "QT interval (ms)",
    colour = "Individual"
  ) +
  theme_minimal(base_size = 14)

###by individual
library(ggplot2)
library(dplyr)

ggplot(all_dives, aes(x = RR.ms, y = QT, colour = Individual)) +
  geom_point(size = 3, alpha = 0.7) +            # raw points
  geom_smooth(method = "lm", se = FALSE) +       # linear model per individual
  labs(
    x = "RR interval (ms)",
    y = "QT interval (ms)",
    colour = "Individual"
  ) +
  theme_minimal(base_size = 14)

library(ggplot2)
library(dplyr)

# Summarize by RR group and individual
all_dives_summary <- all_dives %>%
  mutate(RR_group = case_when(
    RR.ms < 1000 ~ "Fast HR",
    RR.ms < 2000 ~ "Medium HR",
    TRUE ~ "Slow HR"
  )) %>%
  group_by(Individual, RR_group) %>%
  summarise(
    mean_QT = mean(QT, na.rm = TRUE),
    sd_QT = sd(QT, na.rm = TRUE),
    .groups = "drop"
  )

# Plot
ggplot(all_dives, aes(x = RR.ms, y = QT)) +
  geom_point(size = 3, alpha = 0.7, colour = "steelblue") +  # individual points
  geom_errorbar(
    data = all_dives_summary,
    aes(
      x = case_when(
        RR_group == "Fast HR" ~ 500,
        RR_group == "Medium HR" ~ 1275,
        RR_group == "Slow HR" ~ 2500
      ),
      ymin = mean_QT - sd_QT,
      ymax = mean_QT + sd_QT
    ),
    inherit.aes = FALSE,
    width = 100,
    colour = "black",
    size = 1
  ) +
  geom_point(
    data = all_dives_summary,
    aes(
      x = case_when(
        RR_group == "Fast HR" ~ 500,
        RR_group == "Medium HR" ~ 1275,
        RR_group == "Slow HR" ~ 2500
      ),
      y = mean_QT
    ),
    inherit.aes = FALSE,
    colour = "black",
    size = 4
  ) +
  facet_wrap(~Individual) +  # separate panel per individual
  labs(
    x = "RR interval (ms)",
    y = "QT interval (ms)"
  ) +
  theme_minimal(base_size = 14)

###use density peaks to determine where to calculate means
library(dplyr)
library(stats)

# Find density of RR
dens <- density(all_dives$RR.ms, na.rm = TRUE)

# Find peaks
peaks <- dens$x[which(diff(sign(diff(dens$y))) == -2)]
peaks

###Try k-means clustering
library(dplyr)
library(ggplot2)

all_dives_clean <- all_dives %>%
  filter(!is.na(RR.ms))  # remove rows where RR.ms is NA

set.seed(123)
k_clusters <- kmeans(all_dives_clean$RR.ms, centers = 3)

all_dives_clean <- all_dives_clean %>%
  mutate(RR_group = factor(k_clusters$cluster, labels = c("Cluster 1", "Cluster 2", "Cluster 3")))


# # Step 1: Assign clusters based on RR using k-means
# set.seed(123)  # for reproducibility
# k_clusters <- kmeans(all_dives$RR.ms, centers = 3)
# all_dives <- all_dives %>%
#   mutate(RR_group = factor(k_clusters$cluster, labels = c("Cluster 1", "Cluster 2", "Cluster 3")))

# Step 2: Summarize QT per RR_group and Individual
all_dives_clean_summary <- all_dives_clean %>%
  group_by(Individual, RR_group) %>%
  summarise(
    mean_QT = mean(QT, na.rm = TRUE),
    sd_QT = sd(QT, na.rm = TRUE),
    .groups = "drop"
  )

# Step 3: Compute cluster centers (mean RR) for plotting
cluster_centers <- all_dives_clean %>%
  group_by(Individual, RR_group) %>%
  summarise(RR_center = mean(RR.ms), .groups = "drop")

# Merge centers into summary for plotting
all_dives_clean_summary <- all_dives_clean_summary %>%
  left_join(cluster_centers, by = c("Individual", "RR_group"))

# Step 4: Plot per individual
ggplot(all_dives_clean, aes(x = RR.ms, y = QT)) +
  geom_point(size = 3, alpha = 0.7, colour = "steelblue") +  # raw points
  geom_errorbar(
    data = all_dives_clean_summary,
    aes(
      x = RR_center,
      ymin = mean_QT - sd_QT,
      ymax = mean_QT + sd_QT
    ),
    inherit.aes = FALSE,
    width = 100,
    colour = "black",
    size = 1
  ) +
  geom_point(
    data = all_dives_clean_summary,
    aes(x = RR_center, y = mean_QT),
    inherit.aes = FALSE,
    colour = "black",
    size = 4
  ) +
  facet_wrap(~Individual) +  # one panel per individual
  labs(
    x = "RR interval (ms)",
    y = "QT interval (ms)"
  ) +
  theme_minimal(base_size = 14)

#####APPLY CLUSTERING ANALYSIS TO QT OVER DIVE DURATION#####
library(dplyr)
library(ggplot2)

library(dplyr)
library(ggplot2)

# Remove NAs
all_dives_clean <- all_dives %>% filter(!is.na(RR.ms))

# K-means clustering
set.seed(123)
k_clusters <- kmeans(all_dives_clean$RR.ms, centers = 3)

# Create a data frame of cluster centers to sort
cluster_info <- data.frame(
  cluster = 1:3,
  center = k_clusters$centers
) %>%
  arrange(desc(center)) %>%  # sort so highest RR = Low HR
  mutate(RR_group = c("High HR", "Medium HR", "Low HR"))

# Map cluster number to RR_group
all_dives_clean <- all_dives_clean %>%
  mutate(RR_group = factor(cluster_info$RR_group[k_clusters$cluster], 
                           levels = c("High HR", "Medium HR", "Low HR")))

ggplot(all_dives_clean, aes(x = Time_from_underwater_start, y = QT, colour = RR_group)) +
  geom_point(size = 1, alpha = 0.5) +
  # geom_point(size = 3, alpha = 0.7) +
  # geom_line(aes(group = RR_group), alpha = 0.5) +  # connects points in same RR group
  # Quadratic smooth per individual +
  geom_smooth(
    aes(group = 1),
    method = "lm",
    formula = y ~ poly(x, 2),
    color = "black",
    se = FALSE
  ) +
  facet_wrap(~Individual, scales = "free_x") +    # one panel per individual
  labs(
    x = "Time",
    y = "QT interval (ms)",
    colour = "RR group"
  ) +
  theme_minimal(base_size = 14)

###Looking at the HR groupings
library(dplyr)

all_dives_clean %>%
  group_by(RR_group) %>%
  summarise(
    mean_HR = mean(HR, na.rm = TRUE),
    .groups = "drop"
  )

#by individual
library(dplyr)

all_dives_clean %>%
  group_by(Individual, RR_group) %>%
  summarise(
    mean_HR = mean(HR, na.rm = TRUE),
    sd_HR = sd(HR, na.rm = TRUE),
    n = sum(!is.na(HR)),
    .groups = "drop"
  )






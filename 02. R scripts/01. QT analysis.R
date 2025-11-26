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

# ---- 6) Build end markers for plotting (end time relative to start, and nearest QT at end) ----
end_markers <- all_dives %>%
  group_by(DiveID) %>%
  summarise(
    start_time = first(start_time),
    end_time   = first(end_time),
    end_time_rel = if_else(!is.na(start_time) & !is.na(end_time), end_time - start_time, NA_real_),
    end_QT = if (!all(is.na(end_time))) {
      # find QT at the row nearest to end_time
      idx <- which.min(abs(Time - end_time[1]))
      QT[idx]
    } else NA_real_
  ) %>%
  ungroup()

print(end_markers)

# ---- 7) Example plot: overlay with a vertical dashed line at each dive's end ----
ggplot(all_dives, aes(x = Time_from_underwater_start, y = QT, color = DiveID)) +
  geom_line() +
  # add end verticals (use inherit.aes = FALSE because this layer uses different mapping)
  geom_vline(
    data = end_markers %>% filter(!is.na(end_time_rel)),
    aes(xintercept = end_time_rel, color = DiveID),
    linetype = "dashed", inherit.aes = FALSE
  ) +
  labs(x = "Time since underwater start (s)", y = "QT interval (ms)", color = "Dive ID") +
  theme_minimal()


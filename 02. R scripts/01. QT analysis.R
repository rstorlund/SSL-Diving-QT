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
library(tools)
library(ggplot2)

#####Read in the data#####


#####ChatGPT version#####


# list all CSV files
mydir <- '/Users/rstorlund/Library/CloudStorage/OneDrive-UBC/QT/Analysis/SSL-Diving-QT/01. Raw data'
file_list <- list.files(path=mydir, pattern = "\\.csv$", full.names = TRUE)

# read and combine
all_dives <- file_list %>%
  set_names(nm = basename(.)) %>% 
  map_df(~ read_csv(.x) %>%
           mutate(DiveID = file_path_sans_ext(basename(.x))))

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
    # pick the first Time where Behaviour==underwater AND Status==start
    start_time = if (any(Behaviour == "underwater" & Status == "START", na.rm = TRUE)) {
      Time[ which(Behaviour == "underwater" & Status == "START")[1] ]
    } else NA_real_,
    # pick the first Time where Behaviour==underwater AND Status==stop
    end_time = if (any(Behaviour == "underwater" & Status == "STOP", na.rm = TRUE)) {
      Time[ which(Behaviour == "underwater" & Status == "STOP")[1] ]
    } else NA_real_
  ) %>%
  ungroup()

#Join start/end times back to full table
all_dives <- all_dives %>%
  left_join(starts_ends, by = "DiveID") %>%
  mutate(Time_from_underwater_start = Time - start_time)

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
  geom_point() +
  # add end verticals (use inherit.aes = FALSE because this layer uses different mapping)
  geom_vline(
    data = end_markers %>% filter(!is.na(end_time_rel)),
    aes(xintercept = end_time_rel, color = DiveID),
    linetype = "dashed", inherit.aes = FALSE
  ) +
  labs(x = "Time since underwater start (s)", y = "QT interval (ms)", color = "Dive ID") +
  theme_minimal()

#All dives
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
Notes:
  geom_smooth(aes(group = DiveID), method = "lm", se = FALSE) fits a separate linear regression for each dive.

se = FALSE hides confidence intervals; you can set it to TRUE if you want shaded uncertainty.

linetype = "dashed" makes the fitted lines visually distinct from raw traces.

End markers (geom_vline) remain to show the underwater stop time.

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


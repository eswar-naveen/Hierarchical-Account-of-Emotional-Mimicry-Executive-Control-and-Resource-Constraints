library(tidyverse)
library(zoo)
library(ez)
library(permuco)
library(patchwork)
library(effectsize)

rm(list = ls())

data <- read_csv("rawdata_bacs6.csv", show_col_types = FALSE) %>%
  select(-any_of("Participant")) %>%
  mutate(across(-participant, as.numeric))

data <- data %>%
  group_by(participant) %>%
  slice_head(n = 1400) %>%
  ungroup()

remove_outliers <- function(x) {
  mean_val <- mean(x, na.rm = TRUE)
  sd_val <- sd(x, na.rm = TRUE)
  ifelse(x < (mean_val - 3 * sd_val) | x > (mean_val + 3 * sd_val), NA, x)
}

data_cleaned <- data %>%
  mutate(across(where(is.numeric) & !any_of("participant"), remove_outliers)) %>%
  mutate(across(where(is.numeric) & !any_of("participant"), ~na.approx(.x, na.rm = FALSE, rule = 2)))

nested_data <- data_cleaned %>%
  nest(data = -participant)

duplicates_found <- nested_data$participant[duplicated(nested_data$data)]
data_fixed <- data_cleaned %>% filter(!participant %in% duplicates_found)

write_csv(data_fixed, "Cleaned_Interpolated_Data_Fixed_Exp1.csv")

if ("Time" %in% names(data_fixed)) {
  data_fixed <- data_fixed %>% rename(Timestamp = Time)
}

participant_summaries <- data_fixed %>%
  group_by(participant) %>%
  summarise(across(where(is.numeric), 
                   list(mean = ~mean(.x, na.rm = TRUE), median = ~median(.x, na.rm = TRUE)), 
                   .names = "{.col}_{.fn}"))
write_csv(participant_summaries, "participant_Mean_Median_Exp1.csv")

data_downsized <- data_fixed %>%
  arrange(participant, Timestamp) %>%
  group_by(participant) %>%
  mutate(chunk_id = (row_number() - 1) %/% 100) %>%
  group_by(participant, chunk_id) %>%
  summarise(Timestamp = first(Timestamp),
            across(where(is.numeric) & !any_of("Timestamp"), ~mean(.x, na.rm = TRUE)),
            .groups = "drop") %>%
  select(-chunk_id)

write_csv(data_downsized, "Downsized_Averages_with_Timestamps_Exp1.csv")

averaged_by_timestamp <- data_downsized %>%
  group_by(Timestamp) %>%
  summarise(across(where(is.numeric) & !any_of("participant"), ~mean(.x, na.rm = TRUE)))
write_csv(averaged_by_timestamp, "Averages_By_Timestamp_Exp1.csv")

bins <- seq(from = 0.3, to = 1.0, by = 0.05)
bin_labels <- paste0(seq(300, 950, 50), "-", seq(350, 1000, 50), "ms")

binned_data <- data_fixed %>%
  pivot_longer(cols = -c(participant, Timestamp), names_to = "Condition", values_to = "Value") %>%
  mutate(Bin = cut(Timestamp, breaks = bins, labels = bin_labels, include.lowest = TRUE, right = FALSE)) %>%
  filter(!is.na(Bin), !is.na(Value)) %>%
  group_by(participant, Condition, Bin) %>%
  summarise(Mean_Value = mean(Value, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Bin, values_from = Mean_Value)

binned_data %>%
  group_split(Condition) %>%
  walk(~write_csv(.x %>% select(-Condition), paste0("Binned_", unique(.x$Condition), "_Fixed_Exp1.csv")))

zm_conditions <- c("Happy_Low_ZM", "Happy_High_ZM", "Neutral_Low_ZM", "Neutral_High_ZM")

valid_participants <- binned_data %>%
  filter(Condition %in% zm_conditions) %>%
  drop_na(all_of(bin_labels)) %>%
  count(participant) %>%
  filter(n == 4) %>%
  pull(participant)

n_participants <- length(valid_participants)
cat(sprintf("Retained %d balanced participants for Cluster ANOVA.\n", n_participants))

long_cluster_data <- binned_data %>%
  filter(Condition %in% zm_conditions, participant %in% valid_participants) %>%
  mutate(
    Load = if_else(str_detect(Condition, "High"), "High", "Low"),
    Emotion = if_else(str_detect(Condition, "Happy"), "Happy", "Neutral")
  ) %>%
  arrange(participant, Emotion, Load)

eeg_matrix <- long_cluster_data %>%
  select(all_of(bin_labels)) %>%
  as.matrix()
dimnames(eeg_matrix) <- NULL

design_df <- long_cluster_data %>%
  select(participant, Emotion, Load) %>%
  mutate(
    participant = factor(participant),
    Emotion = factor(Emotion, levels = c("Neutral", "Happy")),
    Load = factor(Load, levels = c("Low", "High"))
  ) %>%
  as.data.frame()

set.seed(42)
cluster_anova <- clusterlm(eeg_matrix ~ Emotion * Load + Error(participant / (Emotion * Load)),
                           data = design_df, np = 10000)

extract_cluster_results <- function(effect_name, effect_obj, bin_labels) {
  uncorrected <- as.data.frame(effect_obj$uncorrected)
  df <- data.frame(
    Effect = effect_name,
    Time_Window = bin_labels,
    F_statistic = uncorrected[, grep("statistic", names(uncorrected), ignore.case = TRUE)[1]],
    P_per_bin = uncorrected[, grep("pvalue$|p_value$", names(uncorrected), ignore.case = TRUE)[1]]
  )
  if (!is.null(effect_obj$clustermass$main)) {
    df$P_clustermass <- effect_obj$clustermass$main[, grep("pvalue|p.value", names(effect_obj$clustermass$main), ignore.case = TRUE)[1]]
  } else {
    df$P_clustermass <- NA
  }
  df %>% mutate(Significant_per_bin = P_per_bin < 0.05, Significant_Clustermass = P_clustermass < 0.05)
}

all_results <- imap_dfr(cluster_anova$multiple_comparison, ~extract_cluster_results(.y, .x, bin_labels)) %>%
  mutate(Time_Bin_Index = rep(1:length(bin_labels), length(unique(Effect))))
write_csv(all_results, "Cluster_ANOVA_Results_ZM_Fixed_Exp1.csv")

long_data_zm <- participant_summaries %>%
  select(participant, ends_with("_ZM_mean")) %>%
  pivot_longer(-participant, names_to = "Condition", values_to = "Value") %>%
  mutate(Load = if_else(str_detect(Condition, "High"), "High", "Low"),
         Emotion = if_else(str_detect(Condition, "Happy"), "Happy", "Neutral")) %>%
  drop_na(Value)

valid_ez_participants <- long_data_zm %>% count(participant) %>% filter(n == 4) %>% pull(participant)
long_data_zm <- long_data_zm %>% filter(participant %in% valid_ez_participants)

anova_result_zm <- ezANOVA(data = long_data_zm, dv = Value, wid = participant, 
                           within = .(Load, Emotion), detailed = TRUE)

anova_results_df <- as.data.frame(anova_result_zm$ANOVA) %>%
  mutate(partial_eta_squared = round(SSn / (SSn + SSd), 4))
write_csv(anova_results_df, "ANOVA_Results_ZM_Fixed_Exp1.csv")

f_crit <- qf(0.95, 1, n_participants - 1)
emotion_bins <- 9:14

fstat_plot_data <- all_results %>%
  mutate(Effect = recode(Effect,
                         "Emotion" = "Emotion (Happy vs Neutral)",
                         "Load" = "Load (High vs Low)",
                         "Emotion:Load" = "Load x Emotion Interaction")) %>%
  mutate(Effect = factor(Effect, levels = c("Emotion (Happy vs Neutral)", "Load (High vs Low)", "Load x Emotion Interaction")))

ribbon_data <- fstat_plot_data %>%
  filter(Effect == "Emotion (Happy vs Neutral)" & Time_Bin_Index %in% emotion_bins) %>%
  group_by(Effect) %>%
  reframe(
    interp_time = approx(Time_Bin_Index, F_statistic, n = 500)$x,
    interp_f = approx(Time_Bin_Index, F_statistic, n = 500)$y
  ) %>%
  rename(Time_Bin_Index = interp_time, F_statistic = interp_f) %>%
  mutate(ymin = f_crit, ymax = pmax(F_statistic, f_crit))

fstat_final <- ggplot(fstat_plot_data, aes(x = Time_Bin_Index, y = F_statistic)) +
  geom_ribbon(data = ribbon_data, aes(ymin = ymin, ymax = ymax), fill = "blue", alpha = 0.2) +
  geom_line(color = "darkblue", linewidth = 1) +
  geom_point(aes(shape = Significant_per_bin, color = Significant_per_bin), size = 2) +
  geom_hline(yintercept = f_crit, color = "red", linetype = "dotted", linewidth = 1) +
  facet_wrap(~Effect, ncol = 1) +
  scale_shape_manual(values = c("TRUE" = 17, "FALSE" = 16)) +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "grey50")) +
  scale_x_continuous(breaks = 1:length(bin_labels), labels = bin_labels) +
  theme_bw() +
  labs(x = "Time Window", y = "F-statistic") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom")

ggsave("Figure4_Fstats_Cluster.png", plot = fstat_final, width = 10, height = 8, dpi = 300)

long_avg <- averaged_by_timestamp %>%
  pivot_longer(-Timestamp, names_to = "Condition", values_to = "Value") %>%
  filter(str_detect(Condition, "ZM")) %>%
  mutate(Emotion = if_else(str_detect(Condition, "Happy"), "Happy", "Neutral"),
         Load = if_else(str_detect(Condition, "High"), "High Load", "Low Load"))

zm_ylim <- range(long_avg$Value, na.rm = TRUE)

plot_emotion <- function(emo_data, title) {
  ggplot(emo_data, aes(x = Timestamp, y = Value, color = Load)) +
    geom_line(linewidth = 1) +
    scale_color_manual(values = c("Low Load" = "darkblue", "High Load" = "forestgreen")) +
    scale_x_continuous(breaks = seq(0.3, 1, by = 0.1)) +
    coord_cartesian(ylim = zm_ylim) +
    theme_bw() +
    labs(title = title, x = "Time (s)", y = "Mean ZM Amplitude (MAV z-score)", color = "Load Condition") +
    theme(legend.position = "bottom")
}

fig_happy <- plot_emotion(long_avg %>% filter(Emotion == "Happy"), "(A) Happy Expressions - ZM Activity")
fig_neutral <- plot_emotion(long_avg %>% filter(Emotion == "Neutral"), "(B) Neutral Expressions - ZM Activity")

figure_line <- fig_happy + fig_neutral + plot_layout(ncol = 2)
ggsave("Figure3_ZM_Lineplots.png", plot = figure_line, width = 12, height = 5, dpi = 300)

zm_means_time <- binned_data %>%
  filter(Condition %in% zm_conditions) %>%
  pivot_longer(-c(participant, Condition), names_to = "Bin", values_to = "Value") %>%
  group_by(Condition, Bin) %>%
  summarise(Mean_EMG_Value = mean(Value, na.rm = TRUE), .groups = "drop") %>%
  mutate(Bin = factor(Bin, levels = bin_labels)) %>%
  mutate(Condition_Label = recode(Condition,
                                  "Happy_Low_ZM" = "Happy Low",
                                  "Happy_High_ZM" = "Happy High",
                                  "Neutral_Low_ZM" = "Neutral Low",
                                  "Neutral_High_ZM" = "Neutral High"))

figure_cond <- ggplot(zm_means_time, aes(x = Bin, y = Mean_EMG_Value, color = Condition_Label, group = Condition_Label)) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = c("Happy High" = "forestgreen", "Happy Low" = "darkblue", "Neutral High" = "orange", "Neutral Low" = "red")) +
  theme_bw() +
  labs(title = "Condition Means Across Time Bins (ZM Muscle) - Exp1", 
       x = "Time Bin", y = "Mean EMG Value", color = "Condition") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom")

ggsave("FigureC3_ConditionMeans.png", plot = figure_cond, width = 10, height = 6, dpi = 300)

cat("--- Finished Analysis for Exp1 ---\n")
library(tidyverse)
library(ggplot2)
library(patchwork)
library(effsize)
library(scales)

setwd("/Users/yengisettyeswarnaveen/Desktop/BACs_Psychopy/Beh/Final")

excluded <- c("P4", "P12", "P13", "P21", "P31")

all_files <- list.files(pattern = "\\.csv$", full.names = TRUE)
valid_files <- all_files[!toupper(tools::file_path_sans_ext(basename(all_files))) %in% excluded]

master_list <- list()

for (f in valid_files) {
  pid <- toupper(tools::file_path_sans_ext(basename(f)))
  df <- read.csv(f, stringsAsFactors = FALSE, na.strings = c("", "NA", "None"))
  
  wm_low_corr <- df$key_resp_3.corr
  wm_low_rt <- df$key_resp_3.rt
  wm_high_corr <- df$key_resp_7.corr
  wm_high_rt <- df$key_resp_7.rt
  vs_low_corr <- df$key_resp_2.corr
  vs_low_rt <- df$key_resp_2.rt
  vs_high_corr <- df$key_resp_5.corr
  vs_high_rt <- df$key_resp_5.rt
  
  master_list[[pid]] <- data.frame(
    Participant = pid,
    WM_Low_Accuracy = mean(wm_low_corr, na.rm = T) * 100,
    WM_High_Accuracy = mean(wm_high_corr, na.rm = T) * 100,
    WM_Low_RT = mean(wm_low_rt[wm_low_corr == 1], na.rm = T),
    WM_High_RT = mean(wm_high_rt[wm_high_corr == 1], na.rm = T),
    VS_Low_Accuracy = mean(vs_low_corr, na.rm = T) * 100,
    VS_High_Accuracy = mean(vs_high_corr, na.rm = T) * 100,
    VS_Low_RT = mean(vs_low_rt[vs_low_corr == 1], na.rm = T),
    VS_High_RT = mean(vs_high_rt[vs_high_corr == 1], na.rm = T)
  )
}

master_df <- bind_rows(master_list) %>% arrange(Participant)

master_df$WM_Low_logRT <- log(master_df$WM_Low_RT)
master_df$WM_High_logRT <- log(master_df$WM_High_RT)
master_df$VS_Low_logRT <- log(master_df$VS_Low_RT)
master_df$VS_High_logRT <- log(master_df$VS_High_RT)

write.csv(master_df, "master_behavioral_data.csv", row.names = FALSE)

# descriptives
mean(master_df$WM_Low_Accuracy, na.rm = T)
sd(master_df$WM_Low_Accuracy, na.rm = T)
mean(master_df$WM_High_Accuracy, na.rm = T)
sd(master_df$WM_High_Accuracy, na.rm = T)

mean(master_df$WM_Low_logRT, na.rm = T)
sd(master_df$WM_Low_logRT, na.rm = T)
mean(master_df$WM_High_logRT, na.rm = T)
sd(master_df$WM_High_logRT, na.rm = T)

# normality
shapiro.test(master_df$WM_Low_Accuracy)
shapiro.test(master_df$WM_High_Accuracy)
shapiro.test(master_df$WM_Low_logRT)
shapiro.test(master_df$WM_High_logRT)

# t-tests
ttest_wm_acc <- t.test(master_df$WM_High_Accuracy, master_df$WM_Low_Accuracy, paired = TRUE)
ttest_wm_rt <- t.test(master_df$WM_High_logRT, master_df$WM_Low_logRT, paired = TRUE)
ttest_vs_acc <- t.test(master_df$VS_High_Accuracy, master_df$VS_Low_Accuracy, paired = TRUE)
ttest_vs_rt <- t.test(master_df$VS_High_logRT, master_df$VS_Low_logRT, paired = TRUE)

ttest_wm_acc
ttest_wm_rt
ttest_vs_acc
ttest_vs_rt

# effect sizes
cohen.d(master_df$WM_High_Accuracy, master_df$WM_Low_Accuracy, paired = TRUE)
cohen.d(master_df$WM_High_logRT, master_df$WM_Low_logRT, paired = TRUE)
cohen.d(master_df$VS_High_Accuracy, master_df$VS_Low_Accuracy, paired = TRUE)
cohen.d(master_df$VS_High_logRT, master_df$VS_Low_logRT, paired = TRUE)

# save t-test results
wm_results <- data.frame(
  Comparison = c("WM Accuracy", "WM logRT"),
  t = c(ttest_wm_acc$statistic, ttest_wm_rt$statistic),
  df = c(ttest_wm_acc$parameter, ttest_wm_rt$parameter),
  p = c(ttest_wm_acc$p.value, ttest_wm_rt$p.value)
)
write.csv(wm_results, "wm_ttest_results.csv", row.names = FALSE)

vs_results <- data.frame(
  Comparison = c("VS Accuracy", "VS logRT"),
  t = c(ttest_vs_acc$statistic, ttest_vs_rt$statistic),
  df = c(ttest_vs_acc$parameter, ttest_vs_rt$parameter),
  p = c(ttest_vs_acc$p.value, ttest_vs_rt$p.value)
)
write.csv(vs_results, "vs_ttest_results.csv", row.names = FALSE)

# plotting
theme_paper <- theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    legend.position = "none",
    plot.margin = margin(12, 12, 12, 12)
  )

load_colors <- c("Low Load" = "#4E9BB9", "High Load" = "#D65F5F")

# reshape for plotting
long_wm_acc <- master_df %>%
  select(Participant, WM_Low_Accuracy, WM_High_Accuracy) %>%
  pivot_longer(-Participant, names_to = "Condition", values_to = "Accuracy") %>%
  mutate(Condition = ifelse(Condition == "WM_Low_Accuracy", "Low Load", "High Load"),
         Condition = factor(Condition, levels = c("Low Load", "High Load")))

long_wm_rt <- master_df %>%
  select(Participant, WM_Low_logRT, WM_High_logRT) %>%
  pivot_longer(-Participant, names_to = "Condition", values_to = "logRT") %>%
  mutate(Condition = ifelse(Condition == "WM_Low_logRT", "Low Load", "High Load"),
         Condition = factor(Condition, levels = c("Low Load", "High Load")))

long_vs_acc <- master_df %>%
  select(Participant, VS_Low_Accuracy, VS_High_Accuracy) %>%
  pivot_longer(-Participant, names_to = "Condition", values_to = "Accuracy") %>%
  mutate(Condition = ifelse(Condition == "VS_Low_Accuracy", "Low Load", "High Load"),
         Condition = factor(Condition, levels = c("Low Load", "High Load")))

long_vs_rt <- master_df %>%
  select(Participant, VS_Low_logRT, VS_High_logRT) %>%
  pivot_longer(-Participant, names_to = "Condition", values_to = "logRT") %>%
  mutate(Condition = ifelse(Condition == "VS_Low_logRT", "Low Load", "High Load"),
         Condition = factor(Condition, levels = c("Low Load", "High Load")))

# WM accuracy plot
summ_wm_acc <- long_wm_acc %>%
  group_by(Condition) %>%
  summarise(M = mean(Accuracy, na.rm = T), SE = sd(Accuracy, na.rm = T) / sqrt(n()))

y_max <- max(long_wm_acc$Accuracy, na.rm = T)
y_range <- y_max - min(long_wm_acc$Accuracy, na.rm = T)
ann_y <- y_max + y_range * 0.08

sig_label_acc <- sprintf("t(%d) = %.2f, p = %.3f", 
                         ttest_wm_acc$parameter, ttest_wm_acc$statistic, ttest_wm_acc$p.value)

p_wm_acc <- ggplot() +
  geom_line(data = long_wm_acc, aes(x = Condition, y = Accuracy, group = Participant),
            color = "grey75", linewidth = 0.4, alpha = 0.8) +
  geom_jitter(data = long_wm_acc, aes(x = Condition, y = Accuracy, color = Condition),
              width = 0.07, size = 2.5, alpha = 0.75) +
  geom_errorbar(data = summ_wm_acc, aes(x = Condition, ymin = M - SE, ymax = M + SE),
                width = 0.10, linewidth = 1.1, color = "black") +
  geom_point(data = summ_wm_acc, aes(x = Condition, y = M),
             size = 5.5, shape = 18, color = "black") +
  annotate("segment", x = 1, xend = 2, y = ann_y, yend = ann_y, linewidth = 0.6) +
  annotate("text", x = 1.5, y = ann_y + y_range * 0.03, label = sig_label_acc, size = 3.8) +
  scale_color_manual(values = load_colors) +
  coord_cartesian(clip = "off") +
  labs(title = "A. Accuracy", x = "Cognitive Load Condition", y = "Accuracy (%)") +
  theme_paper

# WM RT plot
summ_wm_rt <- long_wm_rt %>%
  group_by(Condition) %>%
  summarise(M = mean(logRT, na.rm = T), SE = sd(logRT, na.rm = T) / sqrt(n()))

y_max <- max(long_wm_rt$logRT, na.rm = T)
y_range <- y_max - min(long_wm_rt$logRT, na.rm = T)
ann_y <- y_max + y_range * 0.08

sig_label_rt <- sprintf("t(%d) = %.2f, p = %.3f",
                        ttest_wm_rt$parameter, ttest_wm_rt$statistic, ttest_wm_rt$p.value)

p_wm_rt <- ggplot() +
  geom_line(data = long_wm_rt, aes(x = Condition, y = logRT, group = Participant),
            color = "grey75", linewidth = 0.4, alpha = 0.8) +
  geom_jitter(data = long_wm_rt, aes(x = Condition, y = logRT, color = Condition),
              width = 0.07, size = 2.5, alpha = 0.75) +
  geom_errorbar(data = summ_wm_rt, aes(x = Condition, ymin = M - SE, ymax = M + SE),
                width = 0.10, linewidth = 1.1, color = "black") +
  geom_point(data = summ_wm_rt, aes(x = Condition, y = M),
             size = 5.5, shape = 18, color = "black") +
  annotate("segment", x = 1, xend = 2, y = ann_y, yend = ann_y, linewidth = 0.6) +
  annotate("text", x = 1.5, y = ann_y + y_range * 0.03, label = sig_label_rt, size = 3.8) +
  scale_color_manual(values = load_colors) +
  coord_cartesian(clip = "off") +
  labs(title = "B. Reaction Time", x = "Cognitive Load Condition", y = "log(RT) on Correct Trials") +
  theme_paper

fig_wm <- p_wm_acc | p_wm_rt
ggsave("fig_wm_manipulation_check.png", fig_wm, width = 10, height = 6, dpi = 300, bg = "white")

# VS accuracy plot
summ_vs_acc <- long_vs_acc %>%
  group_by(Condition) %>%
  summarise(M = mean(Accuracy, na.rm = T), SE = sd(Accuracy, na.rm = T) / sqrt(n()))

y_max <- max(long_vs_acc$Accuracy, na.rm = T)
y_range <- y_max - min(long_vs_acc$Accuracy, na.rm = T)
ann_y <- y_max + y_range * 0.08

sig_label_vs_acc <- sprintf("t(%d) = %.2f, p = %.3f",
                            ttest_vs_acc$parameter, ttest_vs_acc$statistic, ttest_vs_acc$p.value)

p_vs_acc <- ggplot() +
  geom_line(data = long_vs_acc, aes(x = Condition, y = Accuracy, group = Participant),
            color = "grey75", linewidth = 0.4, alpha = 0.8) +
  geom_jitter(data = long_vs_acc, aes(x = Condition, y = Accuracy, color = Condition),
              width = 0.07, size = 2.5, alpha = 0.75) +
  geom_errorbar(data = summ_vs_acc, aes(x = Condition, ymin = M - SE, ymax = M + SE),
                width = 0.10, linewidth = 1.1, color = "black") +
  geom_point(data = summ_vs_acc, aes(x = Condition, y = M),
             size = 5.5, shape = 18, color = "black") +
  annotate("segment", x = 1, xend = 2, y = ann_y, yend = ann_y, linewidth = 0.6) +
  annotate("text", x = 1.5, y = ann_y + y_range * 0.03, label = sig_label_vs_acc, size = 3.8) +
  scale_color_manual(values = load_colors) +
  coord_cartesian(clip = "off") +
  labs(title = "A. Accuracy", x = "Cognitive Load Condition", y = "Accuracy (%)") +
  theme_paper

# VS RT plot
summ_vs_rt <- long_vs_rt %>%
  group_by(Condition) %>%
  summarise(M = mean(logRT, na.rm = T), SE = sd(logRT, na.rm = T) / sqrt(n()))

y_max <- max(long_vs_rt$logRT, na.rm = T)
y_range <- y_max - min(long_vs_rt$logRT, na.rm = T)
ann_y <- y_max + y_range * 0.08

sig_label_vs_rt <- sprintf("t(%d) = %.2f, p = %.3f",
                           ttest_vs_rt$parameter, ttest_vs_rt$statistic, ttest_vs_rt$p.value)

p_vs_rt <- ggplot() +
  geom_line(data = long_vs_rt, aes(x = Condition, y = logRT, group = Participant),
            color = "grey75", linewidth = 0.4, alpha = 0.8) +
  geom_jitter(data = long_vs_rt, aes(x = Condition, y = logRT, color = Condition),
              width = 0.07, size = 2.5, alpha = 0.75) +
  geom_errorbar(data = summ_vs_rt, aes(x = Condition, ymin = M - SE, ymax = M + SE),
                width = 0.10, linewidth = 1.1, color = "black") +
  geom_point(data = summ_vs_rt, aes(x = Condition, y = M),
             size = 5.5, shape = 18, color = "black") +
  annotate("segment", x = 1, xend = 2, y = ann_y, yend = ann_y, linewidth = 0.6) +
  annotate("text", x = 1.5, y = ann_y + y_range * 0.03, label = sig_label_vs_rt, size = 3.8) +
  scale_color_manual(values = load_colors) +
  coord_cartesian(clip = "off") +
  labs(title = "B. Reaction Time", x = "Cognitive Load Condition", y = "log(RT) on Correct Trials") +
  theme_paper

fig_vs <- (p_vs_acc | p_vs_rt) +
  plot_annotation(
    title = "Visual Search Task Performance by Load Condition",
    subtitle = paste0("Supplementary | N = ", nrow(master_df)),
    theme = theme(
      plot.title = element_text(face = "bold", size = 15, hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5, color = "grey40")
    )
  )

ggsave("fig_vs_supplementary.png", fig_vs, width = 10, height = 6, dpi = 300, bg = "white")

# heatmap
heatmap_df <- master_df %>%
  select(Participant,
         `WM\nLow Load` = WM_Low_Accuracy,
         `WM\nHigh Load` = WM_High_Accuracy,
         `VS\nLow Load` = VS_Low_Accuracy,
         `VS\nHigh Load` = VS_High_Accuracy) %>%
  mutate(Participant = factor(Participant,
                              levels = Participant[order(as.numeric(gsub("P", "", Participant)))])) %>%
  pivot_longer(-Participant, names_to = "Condition", values_to = "Accuracy") %>%
  mutate(Condition = factor(Condition,
                            levels = c("WM\nLow Load", "WM\nHigh Load", "VS\nLow Load", "VS\nHigh Load")))

fig_heatmap <- ggplot(heatmap_df, aes(x = Condition, y = Participant, fill = Accuracy)) +
  annotate("rect", xmin = 2.5, xmax = 2.52,
           ymin = 0.5, ymax = nlevels(heatmap_df$Participant) + 0.5, fill = "grey30") +
  geom_tile(color = "white", linewidth = 0.6) +
  geom_text(aes(label = sprintf("%.1f", Accuracy),
                color = ifelse(Accuracy < 55 | Accuracy > 90, "white", "black")), size = 3.2) +
  scale_fill_gradient2(low = "#D65F5F", mid = "#FFFBF0", high = "#4E9BB9",
                       midpoint = 70, limits = c(0, 100), name = "Accuracy (%)",
                       breaks = c(0, 25, 50, 75, 100)) +
  scale_color_identity() +
  annotate("text", x = 1.5, y = nlevels(heatmap_df$Participant) + 1.2,
           label = "WM Probe Task", fontface = "bold", size = 4) +
  annotate("text", x = 3.5, y = nlevels(heatmap_df$Participant) + 1.2,
           label = "Visual Search Task", fontface = "bold", size = 4) +
  labs(title = "Appendix: Per-Participant Accuracy",
       subtitle = "Red = low accuracy, Blue = high accuracy",
       x = NULL, y = "Participant") +
  coord_cartesian(clip = "off") +
  theme_classic(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 13, hjust = 0.5),
    plot.subtitle = element_text(size = 10, hjust = 0.5, color = "grey40"),
    axis.text.x = element_text(color = "black", size = 11),
    axis.text.y = element_text(color = "black", size = 9),
    axis.title.y = element_text(face = "bold"),
    axis.ticks.x = element_blank(),
    axis.line = element_blank(),
    legend.position = "right",
    legend.key.height = unit(1.2, "cm"),
    plot.margin = margin(25, 15, 10, 10)
  )

ggsave("fig_appendix_accuracy_heatmap.png", fig_heatmap, width = 8, height = 11, dpi = 300, bg = "white")
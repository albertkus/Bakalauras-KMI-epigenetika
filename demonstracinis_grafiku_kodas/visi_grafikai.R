# ==============================================================================
# BAIGIAMOJO DARBO REZULTATŲ VIZUALIZACIJA
# ==============================================================================

setwd("C:/Users/jūsų_direktorija")

# Paketų tikrinimas
if (!require("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!require("dplyr", quietly = TRUE)) install.packages("dplyr")
if (!require("tidyr", quietly = TRUE)) install.packages("tidyr")

library(ggplot2)
library(dplyr)
library(tidyr)
# ==============================================================================
# 1 GRAFIKAS: Elastic Net hiperparametrų derinimas
# ==============================================================================
if(file.exists("data_viz_elnet.csv")) {
  cv_results <- read.csv("data_viz_elnet.csv")
  
  # Fiksuojame tvarką
  cv_results$Modelis <- factor(cv_results$Modelis, 
                               levels = c("Top_500", "Top_1000", "Top_5000", "Top_10000", "Top_15000", "Top_20000", "Top_25000"))
  
  plot1 <- ggplot(cv_results, aes(x = Modelis, y = MSE_Klaida, group = 1)) +
    geom_line(color = "steelblue", linewidth = 1.2) +
    geom_point(color = "darkblue", size = 4) +
    geom_text(aes(label = paste("Liko aktyvių CpG:\n", Liko_Aktyviu_CpG)), 
              vjust = -0.8, color = "darkred", fontface = "bold", size = 4) +
    labs(
      title = "1 pav. Elastic Net hiperparametrų derinimas",
      subtitle = "MSE paklaidos ir atrinktų CpG santykis",
      x = "Pradinis algoritmui paduotų CpG skaičius (išrūšiuota pagal p-reikšmę)",
      y = "Kryžminės patikros paklaida"
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.1, 0.3))) + 
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 11, color = "gray40"),
      axis.title = element_text(size = 12, face = "bold"),
      axis.text.x = element_text(size = 11, color = "black"),
      axis.text.y = element_text(size = 11),
      panel.grid.minor = element_blank() 
    )
  print(plot1)
}


# ==============================================================================
# 2 GRAFIKAS: Bazinių modelių smuiko (Violin) grafikas
# ==============================================================================
if(file.exists("data_viz_violin.csv")) {
  violin_data <- read.csv("data_viz_violin.csv")
  
  # Fiksuojame tvarką
  violin_data$Modelis <- factor(violin_data$Modelis, 
                                levels = c("McCartney", "EpiScore", "Elastic Net (10-Fold)"))
  
  plot2 <- ggplot(violin_data, aes(x = Modelis, y = Pearson_r, fill = Modelis)) +
    geom_violin(trim = FALSE, alpha = 0.5, color = "black", linewidth = 0.8) +
    geom_jitter(width = 0.15, size = 3, color = "black", alpha = 0.7) +
    stat_summary(fun = mean, geom = "point", shape = 23, size = 4, fill = "white", color = "black", stroke = 1) +
    scale_fill_manual(values = c("McCartney" = "#5DADE2", 
                                 "EpiScore" = "#F5B041", 
                                 "Elastic Net (10-Fold)" = "#AF7AC5")) +
    labs(
      title = "Bazinio Elastic Net ir literatūrinių prediktorių palyginimas",
      subtitle = "Baltais rombais pažymėti modelių vidurkiai, juodais taškais – atskiros iteracijos/rinkiniai",
      x = "KMI prognozavimo algoritmas",
      y = "Pirsono koreliacijos koeficientas (r)"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 11, color = "gray40", face = "italic"),
      axis.text.x = element_text(size = 12, face = "bold"),
      axis.text.y = element_text(size = 11),
      axis.title = element_text(size = 12, face = "bold"),
      panel.grid.major.x = element_blank(),
      legend.position = "none"
    )
  print(plot2)
}


# ==============================================================================
# 3 GRAFIKAS: Elastic Net vs Random Forest metrikų palyginimas (3 dalys)
# ==============================================================================
if(file.exists("data_viz_en_rf.csv")) {
  raw_means <- read.csv("data_viz_en_rf.csv")
  
  mean_data <- raw_means %>%
    pivot_longer(cols = everything(), names_to = "Metric", values_to = "Value") %>%
    mutate(
      Modelis = ifelse(grepl("EN_", Metric), "Elastic Net", "Random Forest"),
      Metrika = case_when(
        grepl("Pearson_r", Metric) ~ "1. Koreliacija (Pearson r)",
        grepl("R2", Metric) ~ "2. Paaiškinta dispersija (R²)",
        grepl("MAE", Metric) ~ "3. Vidutinė paklaida (MAE, kg/m²)"
      )
    )
  
  mean_data$Metrika <- factor(mean_data$Metrika, levels = c("1. Koreliacija (Pearson r)", "2. Paaiškinta dispersija (R²)", "3. Vidutinė paklaida (MAE, kg/m²)"))
  
  plot3 <- ggplot(mean_data, aes(x = Modelis, y = Value, fill = Modelis)) +
    geom_bar(stat = "identity", position = "dodge", width = 0.6, color = "black") +
    geom_text(aes(label = round(Value, 3)), vjust = -0.5, size = 4.5, fontface = "bold") +
    facet_wrap(~Metrika, scales = "free_y") +
    scale_fill_manual(values = c("Elastic Net" = "#66c2a5", "Random Forest" = "#fc8d62")) +
    labs(
      title = "Vidutinių modelių tikslumo metrikų palyginimas (10-Fold)",
      x = NULL,
      y = "Metrikos vertė",
      fill = "Modelis"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
      strip.text = element_text(size = 12, face = "bold"),
      strip.background = element_rect(fill = "gray95", color = "black"),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.text.y = element_text(size = 11, color = "black"),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
      legend.position = "bottom",
      legend.title = element_text(face = "bold"),
      legend.text = element_text(size = 11)
    )
  print(plot3)
}


# ==============================================================================
# 4 GRAFIKAS: Kintamųjų svarba (Lollipop)
# ==============================================================================
if(file.exists("data_viz_rf_importance.csv")) {
  imp_data_raw <- read.csv("data_viz_rf_importance.csv")

  imp_data <- imp_data_raw %>%
    arrange(IncMSE) %>%
    mutate(Feature = factor(Feature, levels = Feature))
  
  levels(imp_data$Feature) <- recode(levels(imp_data$Feature),
                                     "EN_Raw_Score" = "Bazinė Elastic Net prognozė",
                                     "Age" = "Tiriamojo amžius",
                                     "CellType" = "Ląstelių tipas",
                                     "EpiScore_Score" = "EpiScore",
                                     "Gender" = "Tiriamojo lytis",
                                     "McCartney_Score" = "McCartney")
  
  plot4 <- ggplot(imp_data, aes(x = IncMSE, y = Feature)) +
    geom_segment(aes(x = 0, xend = IncMSE, y = Feature, yend = Feature), color = "gray50", linewidth = 1.2) +
    geom_point(color = "#7570b3", size = 6) +
    geom_text(aes(label = round(IncMSE, 1)), hjust = -0.4, size = 4, fontface = "bold") +
    scale_x_continuous(expand = expansion(mult = c(0, 0.15))) +
    labs(
      title = "Metamodelio kintamųjų svarba ansamblyje",
      x = "Vidutinės kvadratinės paklaidos padidėjimas pašalinus kintamąjį (%IncMSE)",
      y = "Kintamasis"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      axis.text.y = element_text(size = 11, face = "bold"),
      axis.title = element_text(size = 12, face = "bold"),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.y = element_blank()
    )
  print(plot4)
}

# ==============================================================================
# 5 GRAFIKAS: Literatūrinių modelių efektyvumas per 7 rinkinius
# ==============================================================================

if(file.exists("data_viz_episcore_full.csv") & file.exists("data_viz_mccartney_full.csv")) {
  
  episc_full <- read.csv("data_viz_episcore_full.csv")
  mccart_full <- read.csv("data_viz_mccartney_full.csv")
  
  episc_full$Model <- "EpiScore"
  mccart_full$Model <- "McCartney"
  lit_comparison <- rbind(episc_full, mccart_full)
  
  plot_data_lit <- pivot_longer(
    lit_comparison,
    cols = c("Pearson_r", "R_Squared", "MAE"),
    names_to = "Metrika",
    values_to = "Reiksme"
  )
  
  plot_data_lit$Metrika <- factor(
    plot_data_lit$Metrika,
    levels = c("Pearson_r", "R_Squared", "MAE"),
    labels = c("Pirsono koreliacija (r)", "Paaiškinta dispersija (R²)", "Vidutinė paklaida (MAE)")
  )
  
  plot5 <- ggplot(plot_data_lit, aes(x = Dataset, y = Reiksme, fill = Model)) +
    geom_bar(stat = "identity", color = "black", alpha = 0.85, position = "dodge") +
    facet_grid(Metrika ~ Model, scales = "free_y") +
    geom_text(aes(label = round(Reiksme, 2)), vjust = -0.5, size = 3, position = position_dodge(width = 0.9)) + 
    labs(
      title = "Literatūrinių KMI prediktorių efektyvumas skirtinguose rinkiniuose",
      x = "Tiriamas duomenų rinkinys",
      y = "Metrikos įvertis"
    ) +
    scale_fill_manual(values = c("EpiScore" = "#F5B041", "McCartney" = "#5DADE2")) +
    theme_minimal() +
    theme(
      legend.position = "bottom", 
      plot.title = element_text(size = 14, face = "bold"),
      axis.text.x = element_text(size = 10, angle = 45, hjust = 1, face = "bold"),
      axis.text.y = element_text(size = 10),
      strip.text = element_text(size = 11, face = "bold", color = "black"), 
      panel.grid.major.x = element_blank(),
      panel.border = element_rect(color = "gray80", fill = NA)
    )
  print(plot5)
} else {
  warning("Nerasti 5 grafiko failai.")
}

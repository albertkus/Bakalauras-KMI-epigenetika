# ==============================================================================
# 1 DALIS: DUOMENŲ ĮKĖLIMAS IR FUNKCIJOS SUKŪRIMAS (7 IMTYS)
# ==============================================================================

setwd("C:/Users/jūsų_direktorija") # svarbu nustatyti savo esamą direktoriją

###############################

# RDS formato duomenų užkrovimas 
# Atvejui jei jie jau buvo užkrauti kitame prediktoriuje - nereikia to daryti pakartotinai
olympic <- readRDS("olympic.rds")
women_training <- readRDS("women_training.rds")

# Dvorkin reikalauja rankinio KMI skaičiavimo pagal esamus duomenis
dvorkin <- readRDS("dvorkin 1.rds")
dvorkin_pheno <- attr(dvorkin, ".annmatrix.cann")
dvorkin_pheno$height_m <- dvorkin_pheno$height * 0.0254
dvorkin_pheno$bmi <- dvorkin_pheno$weight / (dvorkin_pheno$height_m^2)
attr(dvorkin, ".annmatrix.cann") <- dvorkin_pheno
actual_bmi_dv <- attr(dvorkin, ".annmatrix.cann")$bmi
valid_test_dv <- !is.na(actual_bmi_dv)
dvorkin_clean <- dvorkin[, valid_test_dv]
actual_bmi_clean_dv <- actual_bmi_dv[valid_test_dv]

salas6 <- readRDS("salas.rds")
salas12 <- readRDS("salas2.rds")
wang <- readRDS("wang.rds")
zhang <- readRDS("zhang.rds")

###############################


# 1. Įkeliami Elnet EpiScore modelio svoriai
episcore_predictor <- read.csv("BMI_Elnet_EpiScore_weights.csv")
rownames(episcore_predictor) <- episcore_predictor$CpG

# 2. Universali funkcija prediktorių testavimui
evaluate_episcore <- function(dataset_name, dataset_matrix, predictor_df, weight_col = "Coefficient") {
  
  # A. Ištraukiami pacientų klinikinių duomenų lentelę
  pheno_data <- attr(dataset_matrix, ".annmatrix.cann")
  if (is.null(pheno_data)) pheno_data <- dataset_matrix$annmatrix
  
  # B. Automatinis KMI ir amžiaus atradimas
  cols_lower <- tolower(colnames(pheno_data))
  
  bmi_col <- colnames(pheno_data)[grep("bmi", cols_lower)[1]]
  actual_bmi <- as.numeric(pheno_data[[bmi_col]])
  
  age_col <- colnames(pheno_data)[grep("age", cols_lower)[1]]
  age <- as.numeric(pheno_data[[age_col]])
  
  # Išvalomos trūkstamos (NA) reikšmes tiek KMI, tiek amžiaus
  valid_patients <- !is.na(actual_bmi) & !is.na(age)
  data_clean <- dataset_matrix[, valid_patients]
  actual_bmi_clean <- actual_bmi[valid_patients]
  
  # C. Transponuojama matricą (Pacientai = Eilutės, CpG = Stulpeliai)
  data_matrix_t <- t(data_clean)
  
  # D. Dinaminis sankirtos kūrimas
  safe_cpgs <- intersect(predictor_df$CpG, colnames(data_matrix_t))
  
  # E. Surenkami tik reikalingi CpG ir sulygiuojami svoriai
  data_sub <- data_matrix_t[, safe_cpgs]
  weights_ordered <- predictor_df[safe_cpgs, weight_col]
  
  # F. Matricos daugyba
  raw_scores <- as.numeric(data_sub %*% weights_ordered)
  
  # G. Kalibravimas
  calibrator <- lm(actual_bmi_clean ~ raw_scores)
  predicted_bmi <- predict(calibrator)
  
  # H. Suskaičiuojamos metrikos
  pearson_r <- cor(actual_bmi_clean, predicted_bmi)
  r_squared <- summary(calibrator)$r.squared
  mae <- mean(abs(actual_bmi_clean - predicted_bmi))
  
  # Apskaičiuojamas išlaikytas svorio procentas
  total_weight <- sum(abs(predictor_df[[weight_col]]))
  kept_weight <- sum(abs(weights_ordered))
  weight_retained <- (kept_weight / total_weight) * 100
  
  # I. Grąžinami rezultatai
  return(data.frame(
    Dataset = dataset_name,
    N_Patients = length(actual_bmi_clean),
    CpGs_Matched = length(safe_cpgs),
    Weight_Retained_Pct = round(weight_retained, 2),
    Pearson_r = round(pearson_r, 3),
    R_Squared = round(r_squared, 3),
    MAE = round(mae, 3)
  ))
}

# ==============================================================================
# 2 DALIS: MODELIO TESTAVIMAS SU VISAIS DATASETAIS
# ==============================================================================

# Paleidžiama funkcija visoms imtims
results_olympic_elnet <- evaluate_episcore("Olympic", olympic, episcore_predictor, "Coefficient")
results_dvorkin_elnet <- evaluate_episcore("Dvorkin", dvorkin, episcore_predictor, "Coefficient")
results_women_elnet   <- evaluate_episcore("Women Training", women_training, episcore_predictor, "Coefficient")
results_salas6_elnet  <- evaluate_episcore("Salas6", salas6, episcore_predictor, "Coefficient")
results_salas12_elnet  <- evaluate_episcore("Salas12", salas12, episcore_predictor, "Coefficient")
results_wang_elnet    <- evaluate_episcore("Wang", wang, episcore_predictor, "Coefficient")
results_zhang_elnet   <- evaluate_episcore("Zhang", zhang, episcore_predictor, "Coefficient")

# Sujungiami visi rezultatai į vieną lentelę
elnet_comparison_table <- rbind(
  results_olympic_elnet, results_dvorkin_elnet,  results_women_elnet,
  results_salas6_elnet, results_salas12_elnet, results_wang_elnet, results_zhang_elnet
)

print("--- EpiScore KMI prediktoriaus tikslumas (7 imtys) ---")
print(elnet_comparison_table)

# ==============================================================================
# 3 DALIS: VIZUALIZACIJA
# ==============================================================================

if (!require("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!require("tidyr", quietly = TRUE)) install.packages("tidyr")
library(ggplot2)
library(tidyr)


plot_data_elnet <- pivot_longer(
  elnet_comparison_table,
  cols = c("Pearson_r", "R_Squared", "MAE"),
  names_to = "Metrika",
  values_to = "Reiksme"
)

plot_data_elnet$Metrika <- factor(
  plot_data_elnet$Metrika,
  levels = c("Pearson_r", "R_Squared", "MAE"),
  labels = c("Pirsono koreliacija (r)", "Paaiškinta dispersija (R²)", "Vidutinė paklaida (MAE)")
)


ggplot(plot_data_elnet, aes(x = Dataset, y = Reiksme, fill = Dataset)) +
  geom_bar(stat = "identity", color = "black", alpha = 0.85) +
  facet_wrap(~ Metrika, scales = "free_y", ncol = 1) + 
  geom_text(aes(label = round(Reiksme, 2)), vjust = -0.5, size = 3) + 
  labs(
    title = "Elnet EpiScore KMI prediktoriaus efektyvumas (7 duomenų rinkiniai)",
    x = "Tiriamas duomenų rinkinys",
    y = "Metrikos įvertis"
  ) +
  scale_fill_brewer(palette = "Set2") + 
  theme_minimal() +
  theme(
    legend.position = "none", 
    plot.title = element_text(size = 14, face = "bold"),
    plot.subtitle = element_text(size = 11, color = "gray40"),
    axis.text.x = element_text(size = 11, angle = 45, hjust = 1, face = "bold"),
    axis.text.y = element_text(size = 10),
    strip.text = element_text(size = 12, face = "bold", color = "black"), 
    panel.grid.major.x = element_blank()
  )

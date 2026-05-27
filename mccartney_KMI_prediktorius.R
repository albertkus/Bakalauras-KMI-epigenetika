# ==============================================================================
# MCCARTNEY KMI PREDIKTORIAUS ĮVERTINIMAS
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

salas1 <- readRDS("salas.rds")
salas2 <- readRDS("salas2.rds")
wang <- readRDS("wang.rds")
zhang <- readRDS("zhang.rds")

###############################

# 1. Įkeliami McCartney svoriai
mccartney_predictor <- read.csv("bmi_predictor_values_from_mccartney.csv")
rownames(mccartney_predictor) <- mccartney_predictor$CpG

# 2. Universali funkcija
evaluate_mccartney <- function(dataset_name, dataset_matrix) {
  
  # A. Ištraukiama pacientų klinikinių duomenų lentelė
  pheno_data <- attr(dataset_matrix, ".annmatrix.cann")
  if (is.null(pheno_data)) pheno_data <- dataset_matrix$annmatrix
  
  # B. Automatinis KMI ir amžiaus atradimas
  cols_lower <- tolower(colnames(pheno_data))
  
  bmi_col <- colnames(pheno_data)[grep("bmi", cols_lower)[1]]
  actual_bmi <- as.numeric(pheno_data[[bmi_col]])
  
  age_col <- colnames(pheno_data)[grep("age", cols_lower)[1]]
  age <- as.numeric(pheno_data[[age_col]])
  
  # Išvalomi tyriamieji su trūkstamais BMI arba amžiumi
  valid_patients <- !is.na(actual_bmi) & !is.na(age)
  data_clean <- dataset_matrix[, valid_patients]
  actual_bmi_clean <- actual_bmi[valid_patients]
  
  # C. Transponuojame matricą (Pacientai = Eilutės, CpG = Stulpeliai)
  data_matrix_t <- t(data_clean)
  
  # D. Dinaminis sankirtos kūrimas
  safe_cpgs <- intersect(mccartney_predictor$CpG, colnames(data_matrix_t))
  
  # E. Surenkami tik reikalingi CpG ir sulygiuojame svorius
  data_sub <- data_matrix_t[, safe_cpgs]
  weights_ordered <- mccartney_predictor[safe_cpgs, "Weight"]
  
  # F. Matricos daugyba
  raw_scores <- as.numeric(data_sub %*% weights_ordered)
  
  # G. Kalibravimas
  calibrator <- lm(actual_bmi_clean ~ raw_scores)
  predicted_bmi <- predict(calibrator)
  
  # H. Suskaičiuojamos metrikos
  pearson_r <- cor(actual_bmi_clean, predicted_bmi)
  r_squared <- summary(calibrator)$r.squared
  mae <- mean(abs(actual_bmi_clean - predicted_bmi))
  weight_retained <- (sum(abs(weights_ordered)) / sum(abs(mccartney_predictor$Weight))) * 100
  
  # I. Grąžiname rezultatus
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

# 3. Paleidžiama funkcija visiems 7 rinkiniams
results_olympic <- evaluate_mccartney("Olympic", olympic)
results_dvorkin <- evaluate_mccartney("Dvorkin", dvorkin)
results_women   <- evaluate_mccartney("Women Training", women_training)
results_salas1  <- evaluate_mccartney("Salas1", salas1)
results_salas2  <- evaluate_mccartney("Salas2", salas2)
results_wang    <- evaluate_mccartney("Wang", wang)
results_zhang   <- evaluate_mccartney("Zhang", zhang)

# 4. Sujungiami visi rezultatai
final_mccartney_table <- rbind(
  results_olympic, results_dvorkin, results_women,
  results_salas1, results_salas2, results_wang, results_zhang
)

print("--- McCartney BMI Prediktoriaus rezultatai (7 imtys) ---")
print(final_mccartney_table)

# ==============================================================================
# REZULTATŲ VIZUALIZACIJA
# ==============================================================================
if (!require("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!require("tidyr", quietly = TRUE)) install.packages("tidyr")

library(ggplot2)
library(tidyr)

plot_data <- pivot_longer(
  final_mccartney_table,
  cols = c("Pearson_r", "R_Squared", "MAE"),
  names_to = "Metrika",
  values_to = "Reiksme"
)

plot_data$Metrika <- factor(
  plot_data$Metrika,
  levels = c("Pearson_r", "R_Squared", "MAE"),
  labels = c("Pirsono koreliacija (r)", "Paaiškinta dispersija (R²)", "Vidutinė paklaida (MAE)")
)

ggplot(plot_data, aes(x = Dataset, y = Reiksme, fill = Dataset)) +
  geom_bar(stat = "identity", color = "black", alpha = 0.85) +
  facet_wrap(~ Metrika, scales = "free_y", ncol = 1) + 
  geom_text(aes(label = round(Reiksme, 2)), vjust = -0.5, size = 3) + 
  labs(
    title = "McCartney KMI prediktoriaus bazinis efektyvumas (7 duomenų rinkiniai)",
    x = "Tiriamas duomenų rinkinys",
    y = "Vertė"
  ) +
  scale_fill_brewer(palette = "Set3") + 
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

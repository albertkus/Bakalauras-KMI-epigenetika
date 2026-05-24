# ==============================================================================
# 1 DALIS: PAKETŲ ĮDIEGIMAS IR DUOMENŲ PARUOŠIMAS
# ==============================================================================

setwd("C:/Users/jūsų_direktorija")

olympic <- readRDS("olympic.rds")
women_training <- readRDS("women_training.rds")

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


if (!require("glmnet", quietly = TRUE)) install.packages("glmnet")
if (!require("BiocManager", quietly = TRUE)) install.packages("BiocManager")
if (!require("limma", quietly = TRUE)) BiocManager::install("limma")

library(glmnet)
library(limma)

all_datasets_list <- list(
  rownames(olympic), rownames(women_training), rownames(dvorkin),
  rownames(salas1), rownames(salas2), rownames(wang), rownames(zhang)
)
common_cpgs <- Reduce(intersect, all_datasets_list)


extract_pheno <- function(original_matrix, dataset_name) {
  anno <- attr(original_matrix, ".annmatrix.cann")
  if (is.null(anno)) anno <- original_matrix$annmatrix
  if (is.null(anno)) stop(paste("KLAIDA:", dataset_name, "neturi pasiekiamų anotacijų!"))
  
  cols_lower <- tolower(colnames(anno))
  
  # 1. KMI ištraukimas
  bmi_col <- colnames(anno)[grep("bmi", cols_lower)[1]]
  bmi <- as.numeric(anno[[bmi_col]])
  
  # 2. Amžiaus ištraukimas
  if ("age..days.365.25." %in% cols_lower) {
    age_col <- colnames(anno)[cols_lower == "age..days.365.25."]
  } else {
    age_col <- colnames(anno)[grep("age", cols_lower)[1]]
  }
  age <- as.numeric(anno[[age_col]])
  
  # 3. Lyties ištraukimas
  sex_col <- colnames(anno)[grep("sex|gender", cols_lower)[1]]
  if (!is.na(sex_col)) {
    gender <- as.character(anno[[sex_col]])
  } else {
    gender <- rep("Unknown", length(bmi))
  }
  
  # 4. Ląstelių tipo ištraukimas
  cell_col <- colnames(anno)[grep("celltype", cols_lower)[1]]
  if (!is.na(cell_col)) {
    celltype <- as.character(anno[[cell_col]])
  } else {
    celltype <- rep("Unknown", length(bmi))
  }
  
  # 5. QC PARAMETRŲ IŠTRAUKIMAS
  get_qc <- function(name, default_val) {
    idx <- which(cols_lower == name)
    if (length(idx) > 0) return(anno[[idx[1]]]) else return(rep(default_val, length(bmi)))
  }
  
  qc_badsnp <- as.logical(get_qc("qc_badsnp", FALSE))
  qc_badsex <- as.logical(get_qc("qc_badsex", FALSE))
  qc_detection <- as.numeric(get_qc("qc_detection", 1.0)) # 1.0 yra tobula detekcija
  qc_iac <- as.numeric(get_qc("qc_iac", 1.0))             # 1.0 yra tobula koreliacija
  
  df <- data.frame(
    BMI = bmi, 
    Age = age, 
    Gender = as.factor(gender), 
    CellType = as.factor(celltype), 
    Dataset = dataset_name,
    qc_badsnp = qc_badsnp,
    qc_badsex = qc_badsex,
    qc_detection = qc_detection,
    qc_iac = qc_iac,
    stringsAsFactors = FALSE
  )
  return(df)
}

pheno_olympic <- extract_pheno(olympic, "Olympic")
pheno_women <- extract_pheno(women_training, "Women")
pheno_dvorkin <- extract_pheno(dvorkin, "Dvorkin")
pheno_salas1 <- extract_pheno(salas1, "Salas1")
pheno_salas2 <- extract_pheno(salas2, "Salas2")
pheno_wang <- extract_pheno(wang, "Wang")
pheno_zhang <- extract_pheno(zhang, "Zhang")

filter_valid <- function(df) {
  return(!is.na(df$BMI) & !is.na(df$Age))
}

v_olympic <- filter_valid(pheno_olympic)
v_women <- filter_valid(pheno_women)
v_dvorkin <- filter_valid(pheno_dvorkin)
v_salas1 <- filter_valid(pheno_salas1)
v_salas2 <- filter_valid(pheno_salas2)
v_wang <- filter_valid(pheno_wang)
v_zhang <- filter_valid(pheno_zhang)

pheno_combined <- rbind(
  pheno_olympic[v_olympic, ], pheno_women[v_women, ], pheno_dvorkin[v_dvorkin, ],
  pheno_salas1[v_salas1, ], pheno_salas2[v_salas2, ], pheno_wang[v_wang, ], pheno_zhang[v_zhang, ]
)

dna_combined <- cbind(
  olympic[common_cpgs, v_olympic],
  women_training[common_cpgs, v_women],
  dvorkin[common_cpgs, v_dvorkin],
  salas1[common_cpgs, v_salas1],
  salas2[common_cpgs, v_salas2],
  wang[common_cpgs, v_wang],
  zhang[common_cpgs, v_zhang]
)

# ==============================================================================
# 2 DALIS: DUOMENŲ KOKYBĖS KONTROLĖ (QC FILTRAVIMAS)
# ==============================================================================

blogi_meginiai <- with(pheno_combined, {
  (qc_badsnp == TRUE) | 
    (qc_badsex == TRUE) | 
    (qc_detection < 0.95) | 
    (qc_iac < 0.4)
})

blogi_meginiai[is.na(blogi_meginiai)] <- FALSE

pasalinta_meginiu <- sum(blogi_meginiai)

good_qc_idx <- !blogi_meginiai
pheno_combined <- pheno_combined[good_qc_idx, ]
dna_combined <- dna_combined[, good_qc_idx]

pheno_combined$qc_badsnp <- NULL
pheno_combined$qc_badsex <- NULL
pheno_combined$qc_detection <- NULL
pheno_combined$qc_iac <- NULL

# ==============================================================================
# 2.5 DALIS: UNIKALIŲ PACIENTŲ ATRANKA (VISOMS KOHORTOMS)
# ==============================================================================

patient_signatures <- paste(pheno_combined$Dataset, pheno_combined$BMI, pheno_combined$Age, pheno_combined$Gender, pheno_combined$CellType, sep = "_")

is_duplicate <- duplicated(patient_signatures)

isvalyta_dublikatu <- sum(is_duplicate)
print(paste("Visose kohortose bendrai rasta ir pašalinta dublikatų:", isvalyta_dublikatu))

pheno_combined <- pheno_combined[!is_duplicate, ]
dna_combined <- dna_combined[, !is_duplicate]

print(paste("Liko visiškai švarių mėginių tolimesnei analizei:", nrow(pheno_combined)))

# ==============================================================================
# 3 DALIS: BIOLOGINIS IŠANKSTINIS FILTRAVIMAS SU LIMMA
# ==============================================================================

design <- model.matrix(~ BMI + Age + factor(Dataset), data = pheno_combined)

fit <- lmFit(dna_combined, design)
fit <- eBayes(fit)

limma_results <- topTable(fit, coef = "BMI", number = Inf, sort.by = "P")

# ==============================================================================
# 4 DALIS: CPG REITINGAVIMAS PAGAL STIPRUMĄ
# ==============================================================================

total_significant <- sum(limma_results$adj.P.Val <= 0.05, na.rm = TRUE)
print(paste("Iš viso griežtą FDR <= 0.05 ribą peržengė:", total_significant, "CpG taškų."))

top_500_cpgs  <- rownames(limma_results)[1:500]
top_1k_cpgs   <- rownames(limma_results)[1:1000]
top_5k_cpgs   <- rownames(limma_results)[1:5000]
top_10k_cpgs  <- rownames(limma_results)[1:10000]
top_15k_cpgs  <- rownames(limma_results)[1:15000]
top_20k_cpgs  <- rownames(limma_results)[1:20000]
top_25k_cpgs  <- rownames(limma_results)[1:25000]


# ==============================================================================
# 5 DALIS: HIPERPARAMETRŲ DERINIMAS
# ==============================================================================

cpg_lists <- list(
  Top_500 = top_500_cpgs,
  Top_1000 = top_1k_cpgs,
  Top_5000 = top_5k_cpgs,
  Top_10000 = top_10k_cpgs,
  Top_15000 = top_15k_cpgs,
  Top_20000 = top_20k_cpgs,
  Top_25000 = top_25k_cpgs
)

ml_y_target <- pheno_combined$BMI

cv_results <- data.frame(
  Modelis = character(), 
  MSE_Klaida = numeric(), 
  Liko_Aktyviu_CpG = numeric(),
  stringsAsFactors = FALSE
)

for (name in names(cpg_lists)) {
  set.seed(2026)
  print(paste("Treniruojamas modelis su", name, "taškų..."))
  
  selected_cpgs <- cpg_lists[[name]]
  
  ml_x_matrix <- t(dna_combined[selected_cpgs, ])
  
  gc()
  
  fit_cv <- cv.glmnet(
    x = ml_x_matrix,
    y = ml_y_target,
    alpha = 0.5,
    family = "gaussian",
    nfolds = 5
  )

  best_mse <- min(fit_cv$cvm)
  best_coefs <- coef(fit_cv, s = "lambda.min")
  active_cpgs <- sum(best_coefs != 0) - 1 

  cv_results <- rbind(cv_results, data.frame(
    Modelis = name,
    MSE_Klaida = round(best_mse, 3),
    Liko_Aktyviu_CpG = active_cpgs
  ))
}

print(cv_results)

best_idx <- which.min(cv_results$MSE_Klaida)
best_model_str <- as.character(cv_results$Modelis[best_idx])
optimal_n_cpgs <- as.numeric(gsub("Top_", "", best_model_str))

# ==============================================================================
# 5 DALIES REZULTATŲ VIZUALIZACIJA
# ==============================================================================
if (!require("ggplot2", quietly = TRUE)) install.packages("ggplot2")
library(ggplot2)

cv_results$Modelis <- factor(cv_results$Modelis, levels = c("Top_500", "Top_1000", "Top_5000", "Top_10000", "Top_15000", "Top_20000", "Top_25000"))

ggplot(cv_results, aes(x = Modelis, y = MSE_Klaida, group = 1)) +
  geom_line(color = "steelblue", linewidth = 1.2) +
  geom_point(color = "darkblue", size = 4) +
  geom_text(aes(label = paste("Liko aktyvių CpG:\n", Liko_Aktyviu_CpG)), 
            vjust = -0.8, color = "darkred", fontface = "bold", size = 4) +
  labs(
    title = "Elastic Net hiperparametrų derinimas",
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

# ==============================================================================
# 5.5 DALIS: MCCARTNEY IR EPISCORE
# ==============================================================================

dna_combined_t <- t(dna_combined)
# McCartney
mccartney_predictor <- read.csv("bmi_predictor_values_from_mccartney.csv")
rownames(mccartney_predictor) <- mccartney_predictor$CpG
mcc_cpgs <- intersect(mccartney_predictor$CpG, colnames(dna_combined_t))
pheno_combined$McCartney_Score <- as.numeric(dna_combined_t[, mcc_cpgs] %*% mccartney_predictor[mcc_cpgs, "Weight"])

# EpiScore įverčiai
elnet_predictor <- read.csv("BMI_Elnet_EpiScore_weights.csv")
rownames(elnet_predictor) <- elnet_predictor$CpG
epi_cpgs <- intersect(elnet_predictor$CpG, colnames(dna_combined_t))
pheno_combined$EpiScore_Score <- as.numeric(dna_combined_t[, epi_cpgs] %*% elnet_predictor[epi_cpgs, "Coefficient"])


# ==============================================================================
# 6 DALIS: NA REIKŠMIŲ VALYMAS PRIEŠ RANDOM FOREST
# ==============================================================================

pheno_combined$Gender <- as.character(pheno_combined$Gender)
pheno_combined$Gender[is.na(pheno_combined$Gender)] <- "Unknown"
pheno_combined$Gender <- as.factor(pheno_combined$Gender)

pheno_combined$CellType <- as.character(pheno_combined$CellType)
pheno_combined$CellType[is.na(pheno_combined$CellType)] <- "Unknown"
pheno_combined$CellType <- as.factor(pheno_combined$CellType)

pheno_combined$McCartney_Score[is.na(pheno_combined$McCartney_Score)] <- median(pheno_combined$McCartney_Score, na.rm = TRUE)
pheno_combined$EpiScore_Score[is.na(pheno_combined$EpiScore_Score)] <- median(pheno_combined$EpiScore_Score, na.rm = TRUE)

# ==============================================================================
# 7 DALIS: ANSAMBLIS SU 10-FOLD CV
# ==============================================================================
n_cpgs <- optimal_n_cpgs

mega_ensemble_results <- data.frame(
  Foldas = character(),
  EN_Pearson_r = numeric(), EN_R2 = numeric(), EN_MAE = numeric(),
  RF_Pearson_r = numeric(), RF_R2 = numeric(), RF_MAE = numeric(),
  stringsAsFactors = FALSE
)
importance_list <- list()

set.seed(2026)
k_folds <- 10
folds <- sample(rep(1:k_folds, length.out = nrow(pheno_combined)))

for (k in 1:k_folds) {
  print(paste("PASLEPTA DALIS", k, "iš", k_folds))
  
  is_test <- folds == k
  pheno_train <- pheno_combined[!is_test, ]
  pheno_test  <- pheno_combined[is_test, ]
  dna_train <- dna_combined[, !is_test]
  dna_test  <- dna_combined[, is_test]
  
  # LIMMA Atranka
  design_train <- model.matrix(~ BMI + Age + factor(Dataset), data = pheno_train)
  fit <- lmFit(dna_train, design_train)
  fit <- eBayes(fit)
  limma_res_train <- topTable(fit, coef = "BMI", number = n_cpgs, sort.by = "P")
  top_cpgs_train <- rownames(limma_res_train)
  
  ml_x_train <- t(dna_train[top_cpgs_train, ])
  ml_y_train <- pheno_train$BMI
  ml_x_test <- t(dna_test[top_cpgs_train, ])
  gc()
  
  # ==========================================
  # 1-AS LYGIS: ELASTIC NET (EN)
  # ==========================================
  set.seed(2026) 
  fit_cv <- cv.glmnet(x = ml_x_train, y = ml_y_train, alpha = 0.5, family = "gaussian", nfolds = 5, keep = TRUE)
  
  best_lambda_idx <- which(fit_cv$lambda == fit_cv$lambda.min)
  oof_train_preds <- as.numeric(fit_cv$fit.preval[, best_lambda_idx])
  test_raw_preds <- as.numeric(predict(fit_cv, newx = ml_x_test, s = "lambda.min"))
  
  calibrator <- lm(pheno_train$BMI ~ oof_train_preds)
  en_final_test_preds <- predict(calibrator, newdata = data.frame(oof_train_preds = test_raw_preds))
  
  # ==========================================
  # 2-AS LYGIS: RANDOM FOREST METAMODELIS (RF)
  # ==========================================
  rf_train_data <- data.frame(
    BMI = pheno_train$BMI,
    EN_Raw_Score = oof_train_preds,
    McCartney_Score = pheno_train$McCartney_Score,
    EpiScore_Score = pheno_train$EpiScore_Score,
    Age = pheno_train$Age,
    Gender = pheno_train$Gender,
    CellType = pheno_train$CellType
  )
  
  set.seed(2026) 
  rf_model <- randomForest(BMI ~ ., data = rf_train_data, ntree = 3000, importance = TRUE)
  
  importance_list[[paste0("Fold_", k)]] <- importance(rf_model, type = 1) # %IncMSE
  
  rf_test_data <- data.frame(
    EN_Raw_Score = test_raw_preds, 
    McCartney_Score = pheno_test$McCartney_Score,
    EpiScore_Score = pheno_test$EpiScore_Score,
    Age = pheno_test$Age,
    Gender = pheno_test$Gender,
    CellType = pheno_test$CellType
  )
  levels(rf_test_data$Gender) <- levels(rf_train_data$Gender)
  levels(rf_test_data$CellType) <- levels(rf_train_data$CellType)
  
  rf_final_test_preds <- predict(rf_model, newdata = rf_test_data)

  actual_bmi <- pheno_test$BMI
  en_pearson <- cor(actual_bmi, en_final_test_preds)
  en_r2 <- summary(lm(actual_bmi ~ en_final_test_preds))$r.squared
  en_mae <- mean(abs(actual_bmi - en_final_test_preds))
  
  rf_pearson <- cor(actual_bmi, rf_final_test_preds)
  rf_r2 <- summary(lm(actual_bmi ~ rf_final_test_preds))$r.squared
  rf_mae <- mean(abs(actual_bmi - rf_final_test_preds))
  
  mega_ensemble_results <- rbind(mega_ensemble_results, data.frame(
    Foldas = paste("Fold", k),
    EN_Pearson_r = round(en_pearson, 3), EN_R2 = round(en_r2, 3), EN_MAE = round(en_mae, 3),
    RF_Pearson_r = round(rf_pearson, 3), RF_R2 = round(rf_r2, 3), RF_MAE = round(rf_mae, 3)
  ))
}

vidurkiai <- data.frame(
  Foldas = "VIDURKIS",
  EN_Pearson_r = round(mean(mega_ensemble_results$EN_Pearson_r), 3),
  EN_R2 = round(mean(mega_ensemble_results$EN_R2), 3),
  EN_MAE = round(mean(mega_ensemble_results$EN_MAE), 3),
  RF_Pearson_r = round(mean(mega_ensemble_results$RF_Pearson_r), 3),
  RF_R2 = round(mean(mega_ensemble_results$RF_R2), 3),
  RF_MAE = round(mean(mega_ensemble_results$RF_MAE), 3)
)
mega_ensemble_results <- rbind(mega_ensemble_results, vidurkiai)

print(mega_ensemble_results)

# ==============================================================================
# 8 DALIS: KINTAMŲJŲ SVARBA
# ==============================================================================

imp_matrix <- do.call(cbind, importance_list)
avg_importance <- data.frame(
  Kintamasis = rownames(imp_matrix),
  Svarba_IncMSE = rowMeans(imp_matrix)
)

avg_importance <- avg_importance[order(-avg_importance$Svarba_IncMSE), ]
print(avg_importance)
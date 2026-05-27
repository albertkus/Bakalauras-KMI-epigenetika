# ==============================================================================
# DEMONSTRACINIS SKRIPTAS
# Random Forest ansamblio veikimo demonstravimas
# ==============================================================================
# Kadangi pilna DNR metilinimo matricų analizė ir Elastic Net treniravimas 
# reikalauja didelių skaičiavimo resursų, šiame skripte naudojami iš anksto
# išgauti epigenetiniai įverčiai (Elastic Net, McCartney, EpiScore).

if (!require("randomForest", quietly = TRUE)) install.packages("randomForest")
library(randomForest)

rf_data <- read.csv("demo_duomenys.csv", stringsAsFactors = TRUE)

demo_ensemble_results <- data.frame(
  Foldas = character(),
  RF_Pearson_r = numeric(), 
  RF_R2 = numeric(), 
  RF_MAE = numeric(),
  stringsAsFactors = FALSE
)
importance_list <- list()

set.seed(2026) 

# Duomenų padalijimas į 10 atsitiktinių grupių kryžminei patikrai
k_folds <- 10
folds <- sample(rep(1:k_folds, length.out = nrow(rf_data)))


for (k in 1:k_folds) {
  # Mokymo ir testavimo imčių atskyrimas
  is_test <- folds == k
  train_data <- rf_data[!is_test, ]
  test_data  <- rf_data[is_test, ]
  
  train_features <- train_data[, !names(train_data) %in% "Dataset"]
  test_features  <- test_data[, !names(test_data) %in% "Dataset"]
  
  # ==========================================
  # RANDOM FOREST METAMODELIS (RF)
  # ==========================================
  set.seed(2026)
  
  # Apmokomas Random Forest su 3000 medžių
  rf_model <- randomForest(BMI ~ ., data = train_features, ntree = 3000, importance = TRUE)
  
  # Kintamųjų svarbos (%IncMSE) išsaugojimas
  importance_list[[paste0("Fold_", k)]] <- importance(rf_model, type = 1) 
  levels(test_features$Gender) <- levels(train_features$Gender)
  levels(test_features$CellType) <- levels(train_features$CellType)
  
  # Galutinis KMI prognozavimas testinei imčiai
  rf_final_test_preds <- predict(rf_model, newdata = test_features)
  
  # Tikslumo metrikų (Pirsono r, R², MAE) apskaičiavimas
  actual_bmi <- test_features$BMI
  rf_pearson <- cor(actual_bmi, rf_final_test_preds)
  rf_r2 <- summary(lm(actual_bmi ~ rf_final_test_preds))$r.squared
  rf_mae <- mean(abs(actual_bmi - rf_final_test_preds))
  
  # Rezultatų išsaugojimas
  demo_ensemble_results <- rbind(demo_ensemble_results, data.frame(
    Foldas = paste("Fold", k),
    RF_Pearson_r = round(rf_pearson, 3), 
    RF_R2 = round(rf_r2, 3), 
    RF_MAE = round(rf_mae, 3)
  ))
}

# Vidutinių tikslumo metrikų apskaičiavimas
vidurkiai <- data.frame(
  Foldas = "VIDURKIS",
  RF_Pearson_r = round(mean(demo_ensemble_results$RF_Pearson_r), 3),
  RF_R2 = round(mean(demo_ensemble_results$RF_R2), 3),
  RF_MAE = round(mean(demo_ensemble_results$RF_MAE), 3)
)

demo_ensemble_results <- rbind(demo_ensemble_results, vidurkiai)

print(demo_ensemble_results)

# Kintamųjų svarbos apdorojimas ir išvedimas
imp_matrix <- do.call(cbind, importance_list)
avg_importance <- data.frame(
  Kintamasis = rownames(imp_matrix),
  Svarba_IncMSE = round(rowMeans(imp_matrix), 2)
)
avg_importance <- avg_importance[order(-avg_importance$Svarba_IncMSE), ]

print(avg_importance)
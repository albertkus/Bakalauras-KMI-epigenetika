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
  RF_Pearson_r = numeric(), RF_R2 = numeric(), RF_MAE = numeric(),
  stringsAsFactors = FALSE
)
importance_list <- list()

for (k in 1:10) {
  train_set <- rf_data[rf_data$Fold == k & rf_data$Type == "train", ]
  test_set  <- rf_data[rf_data$Fold == k & rf_data$Type == "test", ]
  train_features <- train_set[, !names(train_set) %in% c("Fold", "Type", "Dataset")]
  test_features  <- test_set[, !names(test_set) %in% c("Fold", "Type", "Dataset")]
  
  set.seed(2026)
  
  # Apmokomas Random Forest metamodelis
  rf_model <- randomForest(BMI ~ ., data = train_features, ntree = 3000, importance = TRUE)
  
  # Įrašoma kintamųjų svarba
  importance_list[[paste0("Fold_", k)]] <- importance(rf_model, type = 1)
  
  # Sinchronizuojami kategorinių kintamųjų lygiai
  levels(test_features$Gender) <- levels(train_features$Gender)
  levels(test_features$CellType) <- levels(train_features$CellType)
  
  # Prognozuojama testinė dalis
  rf_final_test_preds <- predict(rf_model, newdata = test_features)
  
  # Apskaičiuojama tiksli dalies statistika
  actual_bmi <- test_features$BMI
  rf_pearson <- cor(actual_bmi, rf_final_test_preds)
  rf_r2 <- summary(lm(actual_bmi ~ rf_final_test_preds))$r.squared
  rf_mae <- mean(abs(actual_bmi - rf_final_test_preds))
  
  demo_ensemble_results <- rbind(demo_ensemble_results, data.frame(
    Foldas = paste("Fold", k),
    RF_Pearson_r = round(rf_pearson, 3), 
    RF_R2 = round(rf_r2, 3), 
    RF_MAE = round(rf_mae, 3)
  ))
}

# Skaičiuojami galutinius kryžminės patikros vidurkiai
vidurkiai <- data.frame(
  Foldas = "VIDURKIS",
  RF_Pearson_r = round(mean(demo_ensemble_results$RF_Pearson_r), 3),
  RF_R2 = round(mean(demo_ensemble_results$RF_R2), 3),
  RF_MAE = round(mean(demo_ensemble_results$RF_MAE), 3)
)
demo_ensemble_results <- rbind(demo_ensemble_results, vidurkiai)

print("==================================================")
print("10-FOLD CV RF REZULTATAI (Sutampa su mega_ensemble_results):")
print("==================================================")
print(demo_ensemble_results)

# Apjungiami kintamųjų svarbos rezultatai
imp_matrix <- do.call(cbind, importance_list)
avg_importance <- data.frame(
  Kintamasis = rownames(imp_matrix),
  Svarba_IncMSE = round(rowMeans(imp_matrix), 3)
)
avg_importance <- avg_importance[order(-avg_importance$Svarba_IncMSE), ]

print(avg_importance)

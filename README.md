# Bakalauras-KMI-epigenetika

Vilniaus universiteto bioinformatikos bakalauro baigiamojo darbo kodas. Projektas skirtas kūno masės indekso (KMI) prognozavimui iš DNR metilinimo mikrogardelių (Illumina EPIC/450k) duomenų, taikant mašininio mokymosi algoritmus.

---

## Failų struktūra


### 1. Pagrindiniai analizės skriptai (Reikalauja pilnų `.rds` DNR matricų)
Šie skriptai atlieka pilną duomenų apjungimą, filtravimą ir modelių apmokymą.
* `ansamble_KMI_prediktorius.R` – Pagrindinis R skriptas, apimantis duomenų QC, *LIMMA* biologinį filtravimą, *Elastic Net* hiperparametrų derinimą, išorinių KMI įverčių integraciją ir *Random Forest* ansamblio 10-fold kryžminę patikrą.
* `elnet_episcore_KMI_prediktorius.R` – Literatūrinio *EpiScore* paremto *Elastic Net* prediktoriaus efektyvumo vertinimas apjungtuose duomenų rinkiniuose.
* `mccartney_KMI_prediktorius.R` – Bazinio *McCartney* epigenetinio KMI prediktoriaus testavimas ir tikslumo metrikų fiksavimas.

### 2. Demonstracinė aplinka nereikalaujanti pilnų DNR matricų
* `/demonstracinis_ansamblio_kodas/demo_ansamblis.R` – Atkuria 100 % identiškus *Random Forest* ansamblio kryžminės patikros rezultatus ($R^2$, Pearson $r$, MAE) ir kintamųjų svarbą. 
* `/demonstracinis_ansamblio_kodas/demo_duomenys.csv` – Iš anksto sugeneruoti ir kalibruoti epigenetiniai įverčiai bei pacientų metaduomenys iš `ansamble_KMI_prediktorius.R`, naudojami `demo_ansamblis.R` veikimui.

### 3. Rezultatų vizualizacija
* `/demonstracinis_grafiku_kodas/visi_grafikai.R` – Skriptas, skirtas sugeneruoti baigiamajame darbe panaudotus grafikus. Veikia nuskaitant iš anksto eksportuotus CSV failus, saugomus tame pačiame aplanke.
* `/demonstracinis_grafiku_kodas/` – Aplankas, kuriame saugomi `data_viz_*.csv` failai, reikalingi grafikų braižymui.

### 4. Literatūrinių modelių svoriai
* `BMI_Elnet_EpiScore_weights.csv` – Originalūs, iš anksto apskaičiuoti *EpiScore* modelio DNR metilinimo žymenų svoriai.
* `bmi_predictor_values_from_mccartney.csv` – Originalūs *McCartney* prediktoriaus biožymenys ir koeficientai.

---

## Naudojimo instrukcija ir paleidimas

### A. Norint paleisti pilną analizę
1. Pasirūpinkite, kad visi pradiniai mikrogardelių duomenys (`.rds` formatu) būtų pagrindinėje darbinėje direktorijoje. (Dėl didelės apimties jie į GitHub nekeliami; nuorodos nurodytos baigiamajame darbe).
2. Atsidarykite `ansamble_KMI_prediktorius.R`.
3. Skripto pradžioje esančioje `setwd(...)` komandoje nurodykite savo darbinės direktorijos kelią.
4. Paleiskite skriptą. Skripto pabaigoje automatiškai gausite visus reikalingus skaičiavimus.

### B. Norint atkartoti ansamblio rezultatus
1. Atsisiųskite aplanką `/demonstracinis_ansamblio_kodas/`.
2. Atsidarykite `demo_ansamblis.R` skriptą.
3. Skripto viršuje nustatykite `setwd(...)` į vietą, kurioje išsisaugojote minėtą aplanką.
4. Paleiskite skriptą. Konsolėje pamatysite iteracijų eigą, 10-fold CV metrikų vidurkius ir kintamųjų svarbos sąrašą.

### C. Norint sugeneruoti darbo grafikus
1. Atsisiųskite aplanką `/demonstracinis_grafiku_kodas/`.
2. Atsidarykite `visi_grafikai.R` skriptą.
3. Nustatykite `setwd(...)` į šio aplanko kelią.
4. Paleiskite skriptą. Programoje *RStudio* vienas po kito bus sugeneruoti baigiamajame darbe naudojami grafikai.

# Bakalauras-KMI-epigenetika
Vilniaus Universiteto bioinformatikos bakalauro baigiamojo darbo kodas. KMI prognozavimas iš DNR metilinimo duomenų.

## Failų struktūra
* `ansamble_KMI_prediktorius.R` - Pagrindinis R skriptas, kuriame realizuotas duomenų apjungimas, Elastic Net modelio treniravimas, Random Forest ansamblio kūrimas ir 10-fold kryžminė patikra.
* `elnet_episcore_KMI_prediktorius.R` - R skriptas, skirtas literatūrinio *EpiScore* (paremto *Elastic Net*) KMI prediktoriaus efektyvumo vertinimui nepriklausomuose duomenų rinkiniuose ir rezultatų vizualizacijai.
* `mccartney_KMI_prediktorius.R` - R skriptas, skirtas literatūrinio *McCartney* KMI prediktoriaus testavimui bei gautų tikslumo metrikų grafiniam atvaizdavimui.
  
### Literatūrinių modelių svorių failai (CSV)
* `BMI_Elnet_EpiScore_weights.csv` - Iš anksto apskaičiuoti *EpiScore* modelio metilinimo žymenų svoriai.
* `bmi_predictor_values_from_mccartney.csv` - Bazinio *McCartney* epigenetinio KMI prediktoriaus biožymenys ir jų koeficientai.

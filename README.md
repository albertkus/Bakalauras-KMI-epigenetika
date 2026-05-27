# Bakalauras-KMI-epigenetika

Vilniaus universiteto bioinformatikos bakalauro baigiamojo darbo kodas. KMI prognozavimas iš DNR metilinimo duomenų.

## Failų struktūra

### R Skriptai
* `ansamble_KMI_prediktorius.R` - Pagrindinis R skriptas, kuriame realizuotas duomenų apjungimas, *Elastic Net* modelio treniravimas, *Random Forest* ansamblio kūrimas ir 10-fold kryžminė patikra.
* `elnet_episcore_KMI_prediktorius.R` - R skriptas, skirtas literatūrinio *EpiScore* (paremto *Elastic Net*) KMI prediktoriaus efektyvumo vertinimui nepriklausomuose duomenų rinkiniuose ir rezultatų vizualizacijai.
* `mccartney_KMI_prediktorius.R` - R skriptas, skirtas bazinio literatūrinio *McCartney* KMI prediktoriaus testavimui bei gautų tikslumo metrikų grafiniam atvaizdavimui.
* `demo_ansamblis.R` - Sutrumpintas demonstracinis skriptas, skirtas greitam *Random Forest* ansamblio veikimo atkartojimui bei kintamųjų svarbos analizei, nereikalaujantis didelių skaičiavimo resursų ar pilnų mikrogardelių matricų.

### Duomenų ir svorių failai (CSV)
* `demo_duomenys.csv` - Paruoštas demonstracinis duomenų rinkinys, naudojamas kartu su `demo_ansamblis.R` skriptu. Jame pateikiami iš anksto apskaičiuoti epigenetiniai įverčiai ir pacientų metaduomenys. Stulpelių reikšmės:
    * `BMI` - Faktinis paciento kūno masės indeksas (kg/m²).
    * `EN_Raw_Score` - Bazinio *Elastic Net* modelio sugeneruotas KMI įvertis (Kryžminės patikros metu gautas pirmojo architektūros lygmens rezultatas).
    * `McCartney_Score` - Literatūrinio *McCartney* prediktoriaus apskaičiuotas KMI įvertis.
    * `EpiScore_Score` - Literatūrinio *EpiScore* prediktoriaus apskaičiuotas KMI įvertis.
    * `Age` - Chronologinis paciento amžius (metais).
    * `Gender` - Paciento lytis.
    * `CellType` - Kraujo ląstelių tipas.
    * `Dataset` - Originalaus duomenų rinkinio pavadinimas.
* `BMI_Elnet_EpiScore_weights.csv` - Iš anksto apskaičiuoti *EpiScore* modelio metilinimo žymenų svoriai.
* `bmi_predictor_values_from_mccartney.csv` - Bazinio *McCartney* epigenetinio KMI prediktoriaus biožymenys ir jų koeficientai.

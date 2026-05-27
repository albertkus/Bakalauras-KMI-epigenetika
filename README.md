# Bakalauras-KMI-epigenetika

Vilniaus universiteto bioinformatikos bakalauro baigiamojo darbo kodas. KMI prognozavimas iš DNR metilinimo duomenų.

## Failų struktūra

### R Skriptai
* `ansamble_KMI_prediktorius.R` - Pagrindinis R skriptas, kuriame realizuotas duomenų apjungimas, *Elastic Net* modelio treniravimas, *Random Forest* ansamblio kūrimas ir 10-fold kryžminė patikra.
* `elnet_episcore_KMI_prediktorius.R` - R skriptas, skirtas literatūrinio *EpiScore* (paremto *Elastic Net*) KMI prediktoriaus efektyvumo vertinimui nepriklausomuose duomenų rinkiniuose ir rezultatų vizualizacijai.
* `mccartney_KMI_prediktorius.R` - R skriptas, skirtas bazinio literatūrinio *McCartney* KMI prediktoriaus testavimui bei gautų tikslumo metrikų grafiniam atvaizdavimui.
* `demo_ansamblis.R` - Sutrumpintas demonstracinis skriptas, skirtas greitam *Random Forest* ansamblio veikimo atkartojimui bei kintamųjų svarbos analizei, nereikalaujantis didelių skaičiavimo resursų ar pilnų mikrogardelių matricų.

### Duomenų ir svorių failai (CSV)
* `demo_duomenys.csv` - Struktūrizuotas duomenų rinkinys, pritaikytas `demo_ansamblis.R` skriptui. Stulpelių reikšmės:
    * `Fold` - Kryžminės patikros iteracijos numeris (1–10).
    * `Type` - Duomenų eilutės paskirtis konkrečiame folde (`train` – modelio mokymui, `test` – nepriklausomam testavimui).
    * `BMI` - Faktinis paciento kūno masės indeksas (kg/m²).
    * `EN_Raw_Score` - Pirmojo lygmens *Elastic Net* modelio išgautas epigenetinis KMI signalas (mokymo imtyje naudojami vidiniai *Out-of-Fold* įverčiai, testavimo – grynosios testinės prognozės).
    * `McCartney_Score` - Istorinio *McCartney* modelio sugeneruotas KMI balas.
    * `EpiScore_Score` - Istorinio *EpiScore* modelio sugeneruotas KMI balas.
    * `Age` - Chronologinis paciento amžius.
    * `Gender` - Paciento lytis.
    * `CellType` - Kraujo ląstelių tipas.
    * `Dataset` - Originalios imties pavadinimas.
* `BMI_Elnet_EpiScore_weights.csv` - Iš anksto apskaičiuoti *EpiScore* modelio metilinimo žymenų svoriai.
* `bmi_predictor_values_from_mccartney.csv` - Bazinio *McCartney* epigenetinio KMI prediktoriaus biožymenys ir jų koeficientai.

library("haven")
library("Hmisc")
library("dplyr")
library("e1071")
library("readxl")
library("afex")
library("readbulk")
library("matrixStats")
library ("writexl")

## load the BARCODE to labcode transcription table into R
BARCODE_to_labcode_complete <- read_excel("./data/BARCODE_to_labcode.xlsx")

BARCODE_to_labcode_NSSI <- BARCODE_to_labcode_complete[1:36,]
BARCODE_to_labcode_NSSI$BARCODE <- as.numeric(BARCODE_to_labcode_NSSI$BARCODE)
BARCODE_to_labcode_NSSI <- arrange(BARCODE_to_labcode_NSSI, BARCODE)

BARCODE_to_labcode_HC <- BARCODE_to_labcode_complete[37:68,]
BARCODE_to_labcode_HC$BARCODE <- as.numeric(BARCODE_to_labcode_HC$BARCODE)
BARCODE_to_labcode_HC <- arrange(BARCODE_to_labcode_HC, BARCODE)

## load general data about the NSSI group into "RawData" 
RawData_long <- read_sav("./data/RawData.sav")
RawData_long <- arrange(RawData_long, BARCODE)
RawData <- RawData_long[RawData_long$BARCODE %in% BARCODE_to_labcode_NSSI$BARCODE,]

## load data about the SITBIG results of the NSSI group into SITBIG
SITBIG_long <- read.csv("./data/Auswertung_SITBIG_final.csv", sep=";")
SITBIG_long <- arrange(SITBIG_long, VPN)
SITBIG_long <- add_row(SITBIG_long, .before=33, VPN = 14031033)#add row for proband 14031033 to SITBIG_long
SITBIG <- SITBIG_long[SITBIG_long$VPN %in% BARCODE_to_labcode_NSSI$BARCODE,]

## load data about the comorbid diagnoses of the NSSI group into Comorbid_Diagnoses
Comorbid_Diagnoses_long <- read.csv("./data/Comorbide_diagnosis.csv", sep=";")
Comorbid_Diagnoses_long <- arrange(Comorbid_Diagnoses_long, Probandencode)
Comorbid_Diagnoses <- Comorbid_Diagnoses_long[Comorbid_Diagnoses_long$Probandencode %in% BARCODE_to_labcode_NSSI$BARCODE,]

## load basic data of the HC group into "basic_data_HC"
basic_data_HC_long <- read_excel("./data/HC_weitere_Daten.xlsx")
basic_data_HC <- basic_data_HC_long[1:32, ]
basic_data_HC <- arrange(basic_data_HC, Code)

## load data about the CTQ results of the control group into "CTQ_HC"
CTQ_HC <- read_sav("./data/CTQ_HC.sav")
CTQ_HC$BARCODE[CTQ_HC$BARCODE == 1113011251] <- 1112008408 
CTQ_HC <- arrange(CTQ_HC, BARCODE)

## load data about the weights of the o-1cm (= first segment) & 1-2cm (= second segment) hair samples into R
first_cm_segment_all_weights <- read_excel("./data/weights_1cm_samples_all_weights.xlsx")
second_cm_segment_weights <- read_excel("./data/weights_2cm_samples.xlsx")

## load raw data of the proteomics analysis into separate tables
SV_proteomics <- first_cm_segment_all_weights <- read_excel("./data/SV1--SV88_20ug_abundance_scores&counts.xlsx")

## check raw data for typos and correct them
RawData[RawData$BARCODE == 14031020, "Groesse"] <- 165 #correct height of no.18 from 65 to 165cm
RawData$Staat <- tolower(RawData$Staat) #label german study participants uniformly as "deutsch"


### CREATE DATA FRAMES FOR THE RESULTS

## create new data frame titled "Results" & copy the sheet- and barcode number of the NSSI group into it
Results <- RawData[c("Bogen", "BARCODE")] #add the NSSI group participants' sheet & barcode numbers to results
label(Results$BARCODE) <- "8- or 10-digit subject code"
Results <- Results %>% rename(Sheet_Nr = Bogen)

## add "Labcode" column from the BARCODE_to_labcode df to "Results"
Results <- merge (Results, BARCODE_to_labcode_NSSI[c("BARCODE", "Lab-Code")], by = "BARCODE")

## add new column 'group' that shows whether participants are in the HC-group (0) or the NSSI-group (1)
Group_membership <- rep("1", length(Results$BARCODE))
Results <- cbind(Results, Group = Group_membership)
label(Results$Group) <- "0 = HC; 1 = NSSI"



### ADD & EVALUATE DATA ABOUT THE NSSI GROUP #####

## add basic personal data of the NSSI-group to "Results"
Results <- cbind(Results, RawData[c("Staat", "Alter", "Geschlecht", "Gewicht", "Groesse")]) #add basic individual data (Staat = nationality, Alter = age, Geschlecht = sex, Gewicht = weigth, Größe = heigth) to "Results"
label(Results$Staat) <- "at time of the study"
label(Results$Alter) <- "years of age"
label(Results$Geschlecht) <- "0 = male; 1 = female"
label(Results$Gewicht) <- "[kg]"
label(Results$Groesse) <- "[cm]"
## translate column names and characters from german to english
Results <- Results %>% 
  rename(
          Nationality = Staat, 
          Age = Alter, 
          Sex = Geschlecht, 
          Weigth = Gewicht, 
          Heigth = Groesse
        )

Results <- Results %>% mutate(Nationality = recode(Nationality, deutsch = 'german', österreich = 'austrian', österreichisch =  'austrian', polnisch = 'polish'))

###calculate BMI and add in new column
BMI_calc <- function(w, h){
  round(w / (h/100)^2, digits = 2)
}                           #add function to calculate BMI(rounded to two decimal places) with given weight in kg and height in cm
Results <- cbind(Results, BMI= mapply(BMI_calc, w = RawData$Gewicht, h = RawData$Groesse)) #add column with BMI
label(Results$BMI) <- "m[kg] / h[m]^2"


## add data about education & employment to Results
Results <- cbind(Results, RawData["Jahre_Schule"]) #add: how many years of school?
label(Results$Jahre_Schule) <- "time spent going to school [years]"
Results <- cbind(Results, hoechster_Schulabschluss = RawData$Abschluss) #add: highest school qualification
label(Results$hoechster_Schulabschluss) <- "0 = none, 1 = 2ndary modern school (Hauptschulabschluss); 2 = General Certificate of Secondary Education (Fachoberschulreife/mittlere Reife); 3 = A-levels (Abitur); 4 = sonstige"
Results <- cbind(Results, hoechster_Berufsschulabschluss_Ausbildung = RawData$Beruf) #add: highest professional qualification
label(Results$hoechster_Berufsschulabschluss_Ausbildung) <- "0 = none or semi-skilled (angelernt); 1 = apprenticeship (Lehre); 2 = masters school (Fach-/Meisterschule); 3 = university (Hochschule)"
Results <- cbind(Results, Angestellt = as.numeric(RawData$Berufl_Sit_A > 0 | RawData$Berufl_Sit_C > 0)) #add a new column with study participants that are employed (Berufl_Sit_A + _C)
label(Results$Angestellt) <- "0 = no; 1 = yes"
Results <- cbind(Results, Selbststaendig = as.numeric(RawData$Berufl_Sit_B > 0)) #add new column with participants that are self-employed (Berufl_Sit_B)
label(Results$Selbstst?ndig) <- "0 = no; 1 = yes"
Results <- cbind(Results, In_Schule_Studium_Ausbildung = as.numeric(RawData$Berufl_Sit_F > 0 | RawData$Berufl_Sit_G > 0 | RawData$Berufl_Sit_H > 0)) #add new column with participants that are students | pupils | in training (Berufl_Sit_F + _G + _H)
label(Results$In_Schule_Studium_Ausbildung) <- "0 = no; 1 = yes"
Results <- cbind(Results, EU_Rente = as.numeric(RawData$Berufl_Sit_J > 0)) #add new column with participants that recieve a disability pension (Berufl_Sit_J)
label(Results$EU_Rente) <- "0 = no; 1 = yes"
Results <- cbind(Results, Arbeitslos = as.numeric(RawData$Berufl_Sit_D > 0 | RawData$Berufl_Sit_E > 0)) #add new column with participants that are unemployed (Berufl_Sit_D + _E)
label(Results$Arbeitslos) <- "0 = no; 1 = yes"

Results <- Results %>% 
  rename(
          School_years = Jahre_Schule, 
          Highest_school_qualification = hoechster_Schulabschluss, 
          Highest_professional_qualification = hoechster_Berufsschulabschluss_Ausbildung, 
          In_employment = Angestellt, 
          Self_employed = Selbststaendig, 
          In_school_college_training = In_Schule_Studium_Ausbildung, 
          Retired = EU_Rente, 
          Unemployed = Arbeitslos
        )

### add data about living situation
Results <- cbind(Results, Aktuelle_Wohnsituation = RawData$akt_Wohnsit)
label(Results$Aktuelle_Wohnsituation) <- "1 = lives alone; 2 = lives with family or partner; 3 = lives with parents; 4 = lives in a flatshare; 5 = lives in therapeutic flatshare or psychiatric family care; 6 = lives in psychiatric transition unit"
Results <- Results %>% rename(Living_situation = Aktuelle_Wohnsituation)

### add data about comorbide diagnoses of the NSSI group to "Results"
Results <- merge (Results, Comorbid_Diagnoses, by.x = "BARCODE", by.y = "Probandencode")

# translate diagnoses
#Results <- Results %>% 
 # rename(
  #  other_kind_of_Bipolar_disorder = andere_Bip, 
   # MD_episode = MD_einzel, 
    #MD_recurring = MD_rez, 
    #Dysthymia = Dysthymie, 
    #Depression_not_further_described = Depression_NNB, 
    #
    #Schizophrenia = Schizophrenie, 
    #Schizophreniform_Disorder = Schizophrenif, 
    #Schizoaffective_Disorder =Schizoaffektiv,
    #Delusional_Disorder = Wahn, 
    #Psychosis_temporary = Kurze_Psych, 
    #Psychosis_lasting = Psych_krank, 
    #Psychosis_substance_induced = Substanz_Psych,
    #Psychosis_not_further_described = Psych_NNB, 
    #Dependency = Abhaengigkeit, 
    #Substance_Abuse = Sub_Missbrauch, 
    #
    #Social_Anxiety_Disorder = Soz_Phobie, 
    #Specific_Phobia = Spez_Phobie, 
    #Schizophrenia = Zwang, 
    #PTSD = PTBS, 
    #Generalized_Anxiety = Gen_Angst,
    #Anxiety_Disorder = Angst_Krank, 
    
    #Anxiety_not_further_described = Angst_NNB 
  
  #)


### add data about medication
Results <- cbind(Results, Long_term_Medication = abs(RawData$MEDI - 1))# which participants regularly take medication?
label(Results$Long_term_Medication) <- "0 = no; 1 = yes" # edit labels: 0 = "no meds", 1 = "takes meds"

## how many different meds do they take?
medcols <- grep("Medikament$", label(RawData), ignore.case = TRUE) # create vector "medcols" that contains the indexes of all the columns about (type of) medication

count_med <- function(data, cols, regexp, select=TRUE){
  output_vec <- rep(NA, nrow(data))
  for (patient in 1:nrow(data)){
    medcount = 0
    for (col in cols){
      if (data[patient,col] != ""){ # only look at the medication columns that actually contain data
        if (label(data[patient,col + 1]) == "Dosis" && grepl("bedarf|abgesetzt", data[patient,col + 1], ignore.case = TRUE) == FALSE){ # don't count meds that are taken "bei Bedarf" or "für die Studie abgesetzt"
          if (grepl(regexp, data[patient,col], ignore.case = TRUE) == select) { 
            medcount = medcount + 1
          }
        }
      }
    }
    output_vec[patient] <- medcount
  }
  return(output_vec)
}

num_med <- count_med(RawData, medcols, "pille|eisen|vitamin|ferro|evakadin", FALSE) # for the number of meds taken: don't count the pill, iron or vitamin supplements

Results <- cbind(Results, n_Meds = num_med) # add new column "n_Meds" to results that shows how many different meds each participant regularly takes
label(Results$n_Meds) <- "Amount of different types of long-term medication"

## which kind of medication?
# 1. psychopharmaceuticals
Results <- cbind(Results, SSRI = count_med(RawData, medcols, "Citalopram|Escitalopram|Fluoxetin|Fluvoxamin|Sertralin|Setralin|Paroxetin")) # SSRI
Results <- cbind(Results, SSNRI = count_med(RawData, medcols, "Venl|Milnaneurax|Duloxetin")) # SSNRI 
Results <- cbind(Results, TCA = count_med(RawData, medcols, "Opipram|Trimipramin")) # tricyclic antidepressants
Results <- cbind(Results, Other_antidepressants = count_med(RawData, medcols, "Valdoxan|Hypnorex|Trazodon|Trazedon")) # other antidepressants
Results <- cbind(Results, Low_potency_neuroleptics = count_med(RawData, medcols, "Chlorprothixen|Truxal|Levomepromazin|Neurocil|Pipamperon|Promethazin|Atosil|Atocil|Dominal")) # low potency typical neuroleptics
Results <- cbind(Results, Medium_potency_neuroleptics = count_med(RawData, medcols, "Perazin")) # medium potency typical neuroleptics
Results <- cbind(Results, Atypical_neuroleptics = count_med(RawData, medcols, "Quet|Quiet|Seroquel|Aripiprazol|Abili|Olanzapin")) # atypical neuroleptics
Results <- cbind(Results, Anxiolytics_Hypnotics = count_med(RawData, medcols, "Lorazepam|Tavor|Zopiclon")) # anxiolytics & hypnotics
Results <- cbind(Results, Antiseizure_drugs = count_med(RawData, medcols, "Lamotrigin|Pregabalin|Gabapentin")) # anticonvulsants
Results <- cbind(Results, Stimulants = count_med(RawData, medcols, "Methylphenidat|Ritalin")) # Ritalin

# 2. non-psychopharmaceuticals
Results <- cbind(Results, Cardiovascular_meds = count_med(RawData, medcols, "Bisoprolol|Metoprolol|Metropolol|Valsartan|Ramipril|Doxazosin")) # meds for the treatment of cardiovascular diseases
Results <- cbind(Results, Respiratory_meds = count_med(RawData, medcols, "Inuvair|Salbutamol")) # meds for the treatment of asthma
Results <- cbind(Results, Thyroid_meds = count_med(RawData, medcols, "Thyroxin|Eferox")) # meds for the treatment of thyroid malfunctions
Results <- cbind(Results, Other_non_psychpharmaceuticals = count_med(RawData, medcols, "Loratadin|Atorvastatin|Triptane")) # other (systemic antihistamines, statins, triptanes)

## add labels to the columns about the types of medication the participants regularly take
label4 <- "Number of drugs of this type that are regularly used by the participant"

label(Results$SSRI) <- label4
label(Results$SSNRI) <- label4
label(Results$TCA) <- label4
label(Results$Other_antidepressants) <- label4
label(Results$Low_potency_neuroleptics) <-label4
label(Results$Medium_potency_neuroleptics) <- label4
label(Results$Atypical_neuroleptics) <- label4
label(Results$Anxiolytics_Hypnotics) <- label4
label(Results$Antiseizure_drugs) <- label4
label(Results$Stimulants) <- label4
label(Results$Cardiovascular_meds) <- label4
label(Results$Respiratory_meds) <- label4
label(Results$Thyroid_meds) <- label4
label(Results$Other_non_psychpharmaceuticals) <- label4



### add data about sport 
Results <- cbind(Results, Sport = as.numeric(RawData$Sport)) # add data about how many minutes per week the participants work out
label(Results$Sport) <- "Estimated minutes per week spent exercising"

### add data about smoking
Results <- cbind(Results, Smoking = abs(RawData$rauchen - 2)) # how many participants are smokers?
label(Results$Smoking) <- "0 = non-smoker; 1 = smoker" # edit label: 0 = smoker, 1 = non-smoker
Results <- cbind(Results, n_Cigs_per_day = RawData$rauchen_anzl) # how many cigarettes a day
label(Results$n_Cigs_per_day) <- "Estimated amount of cigarettes smoked per day"

### add data about hormonal contraception
Results <- cbind(Results, Hormonal_contraception = abs(RawData$Verhuetung -2))
label(Results$Hormonal_contraception) <- "0 = no, 1 = yes"




#### evaluation / scoring of the different tests of the NSSI group

### evaluation of the CTQ 
Umpolung <- function(x) {
       6 - x
   }       #add a function to turn over the inverse items

CTQ_list <- vector(mode = "list", length = 5)
names(CTQ_list) <- c("CTQ_EA", "CTQ_PA", "CTQ_SA", "CTQ_EN", "CTQ_PN")
CTQ_list[[1]] <- c("CTQ_3", "CTQ_8", "CTQ_14", "CTQ_18", "CTQ_25")
CTQ_list[[2]] <- c("CTQ_9", "CTQ_11", "CTQ_12", "CTQ_15", "CTQ_17")
CTQ_list[[3]] <- c("CTQ_20", "CTQ_21", "CTQ_23", "CTQ_24", "CTQ_27")
CTQ_list[[4]] <- c("CTQ_5", "CTQ_7", "CTQ_13", "CTQ_19", "CTQ_28")
CTQ_list[[5]] <- c("CTQ_2", "CTQ_26", "CTQ_1", "CTQ_4", "CTQ_6")

CTQHC_list <- vector(mode = "list", length = 5)
names(CTQHC_list) <- c("CTQ_EA", "CTQ_PA", "CTQ_SA", "CTQ_EN", "CTQ_PN")
CTQHC_list[[1]] <- c("CTQ00300", "CTQ00800", "CTQ01400", "CTQ01800", "CTQ02500")
CTQHC_list[[2]] <- c("CTQ00900", "CTQ01100", "CTQ01200", "CTQ01500", "CTQ01700")
CTQHC_list[[3]] <- c("CTQ02000", "CTQ02100", "CTQ02300", "CTQ02400", "CTQ02700")
CTQHC_list[[4]] <- c("CTQ00500", "CTQ00700", "CTQ01300", "CTQ01900", "CTQ02800")
CTQHC_list[[5]] <- c("CTQ00200", "CTQ02600", "CTQ00100", "CTQ00400", "CTQ00600")


ctq_rowsum <- function(list, data, results_data)
{
  subitems = length(list)
  num_patients = nrow(data)
  for (subtest in 1:subitems)
  {
    temp_RawData <- data[list[[subtest]]]
    if (names(list)[subtest] == "CTQ_EN")
    {
      temp_RawData <- mapply(Umpolung, temp_RawData)
    }
    else if (names(list)[subtest] == "CTQ_PN")
    {
      if (deparse(substitute(list)) == "CTQ_list")
      {
      temp_RawData[c("CTQ_2", "CTQ_26")] <- mapply(Umpolung, temp_RawData[c("CTQ_2", "CTQ_26")])
      }
      else
      {
      temp_RawData[c("CTQ00200", "CTQ02600")] <- mapply(Umpolung, temp_RawData[c("CTQ00200", "CTQ02600")])  
      }
    }
    
    temp_results <- numeric(num_patients)
    for (patient in 1:nrow(data))
    {
      na_count <- sum(is.na(temp_RawData[patient,]))
      
      if (na_count == 0)
      {
        temp_results[patient] <- sum(temp_RawData[patient,])
      }
      else if (na_count == 1)
      {
        temp_RawData[patient,][is.na(temp_RawData[patient,])] <- round(mean(as.numeric(temp_RawData[patient,]), na.rm = TRUE))
        temp_results[patient] <- sum(temp_RawData[patient,])
      }
      else
      {
        temp_results[patient] <- NA
      }
    }
    results_data[, names(list)[subtest]] <- temp_results
  }
  return(results_data)
}


Results <- ctq_rowsum(CTQ_list, RawData, Results)

# label the CTQ subscale columns
label(Results$CTQ_EA) <- "Score subscale emotional abuse (EA); ranging from 5 (no EA whatsoever) to 25 (extreme EA)"
label(Results$CTQ_PA) <- "Score subscale physical abuse (PA); ranging from 5 (no PA whatsoever) to 25 (extreme PA)"
label(Results$CTQ_SA) <- "Score subscale sexual abuse (SA); ranging from 5 (no SA whatsoever) to 25 (extreme SA)"
label(Results$CTQ_EN) <- "Score subscale emotional neglect (EN)); ranging from 5 (no EN whatsoever) to 25 (extreme EN)"
label(Results$CTQ_PN) <- "Score subscale physical neglect (PN)); ranging from 5 (no PN whatsoever) to 25 (extreme PN)"

# add column with total score of the CTQ
Results <-  cbind(Results, CTQ_total_score = rowSums(Results[c("CTQ_EA", "CTQ_PA", "CTQ_SA", "CTQ_EN", "CTQ_PN")])) #overall score of the CTQ (= sum of the subscale scores)
label(Results$CTQ_total_score) <- "Total score of the CTQ; ranging from 25 (no ACEs whatsoever) to 125 (extreme ACEs)"


### evaluation of the BSL-23 
Results <- cbind(Results, BSL_23_BPD_symptoms_1 = round(apply(RawData[grepl("BSL_[1-9]$|BSL_1[0-9]$|BSL_2[0-3]$", names(RawData))], 1, mean), digits = 2)) # calculate the MEAN of each patient's BPD-symptoms items (BSL_1 bis BSL_23)
label(Results$BSL_23_BPD_symptoms_1) <- "Mean score of all BSL-23 items (-> severity of BPD symptoms over the past week); ranging from 0 (no BPD symptoms) to 4 (strongly pronounced BPD symptoms)"
Results <- cbind(Results, BSL_23_BPD_symptoms_2 = rowSums(RawData[grepl("BSL_[1-9]$|BSL_1[0-9]$|BSL_2[0-3]$", names(RawData))])) # calculate the SUM of the BPD-symptom items (BSL_1 bis BSL_23)
label(Results$BSL_23_BPD_symptoms_2) <- "Sum score of all BSL-23 items (-> severity of BPD symptoms over the past week); ranging from 0 (no BPD symptoms) to 4 (strongly pronounced BPD symptoms)"
Results <- cbind(Results, BSL_23_condition = RawData$BSL_Bef) # add column with BSL-23 score for global well-being
label(Results$BSL_23_condition) <- "Score for general condition over the past week, ranging from 0 (very bad) to 10 (very good)"
Results <- cbind(Results, BSL_23_dysfunc_behav_1 = round(apply(RawData[grepl("BSL_2[4-9]$|BSL_3[0-4]$", names(RawData))], 1, mean), digits = 2)) # calculate the MEAN of each patient's items assessing dysregulatory behavior  (BSL_24 bis BSL_34)
label(Results$BSL_23_dysfunc_behav_1) <- "Mean score of all BSL-23 items concerning dysfunctional behavior over the past week; ranging from 0 (no dysfunctional behavior whatsoever) to 4 (strongly pronounced dysfunctional behavior)"
Results <- cbind(Results, BSL_23_dysfunc_behav_2 = rowSums(RawData[grepl("BSL_2[4-9]$|BSL_3[0-4]$", names(RawData))])) # calculate the SUM of the dysregulatory behavior items (BSL_24 bis BSL_34)
label(Results$BSL_23_dysfunc_behav_2) <- "Sum score of all BSL-23 items concerning dysfunctional behavior over the past week; ranging from 0 (no dysfunctional behavior whatsoever) to 4 (strongly pronounced dysfunctional behavior)"


### evaluation of the DSS
Results <- cbind(Results, DSS_Mean = round(apply(RawData[grepl("DSS_[1-9]$|DSS_1[0-9]$|BSL_2[0-1]$", names(RawData))], 1, mean), digits = 2)) # calculate the mean of the 21 DSS-items and add it to Results in a new column titled "DSS_Mean"
label(Results$DSS_Mean) <- "Mean score of the 21 DSS items to assess dissociative symptoms over the past week; ranging from 0 (no dissociative symptoms whatsoever) to 10 (experiences constantly dissociative symptoms)"
Results <- cbind(Results, DSS_Sum = rowSums(RawData[grepl("DSS_[1-9]$|DSS_1[0-9]$|BSL_2[0-1]$", names(RawData))])) # calculate the rowsums of the 21 DSS_items and add it to Results in a new column titled "DSS_Sum"
label(Results$DSS_Sum) <- "Sum score of the 21 DSS items to assess dissociative symptoms over the past week; ranging from 0 (no dissociative symptoms whatsoever) to 210 (experiences constantly dissociative symptoms)"


### evaluation of the SNI
Results <- cbind(Results, SNI_1_partner = as.numeric(RawData$SNI_1 == 1 | RawData$SNI_2 == 1)) # add column that tells us whether the participant has a spouse / is in a serious relationship
label(Results$SNI_1_partner) <- "0 = single; 1 = in a serious relationship" # add information about relationship status
Results <- cbind(Results, SNI_2_children = as.numeric(RawData$SNI_4 > 1)) # add column that tells us with how many of his/her children the participant regularly has contact (meaning at least once in two weeks)
label(Results$SNI_2_children) <- "0 = no regular contact with children; 1 = regular contact with children" ### SAME FOR ALL THE OTHER SNI LABELS!!
Results$SNI_2_children[is.na(Results$SNI_2_children)] <- 0
Results <- cbind(Results, SNI_3_parents = as.numeric(RawData$SNI_6 != 1)) # add column that tells us if the participant has regular contact with AT LEAST ONE of his parents
label(Results$SNI_3_parents) <- "0 = no regular contact with parents; 1 = regular contact with parents"
Results <- cbind(Results, SNI_4_in_laws = as.numeric(RawData$SNI_8 != 1)) # add column that tells us if the participant has regular contact with his parents-in-law
label(Results$SNI_4_in_laws) <- "0 = no regular contact with parents-in-law; 1 = regular contact with parents-in-law"
Results <- cbind(Results, SNI_5_other_relatives = as.numeric(RawData$SNI_10 > 1))  # add column that tells us if the participant has regular contact with other members of his family, apart from his partens and/or children
label(Results$SNI_5_other_relatives) <- "regular contact to other relatives: 0 = no; 1 = yes"
Results <- cbind(Results, SNI_6_close_friends = as.numeric(RawData$SNI_12 > 1)) # add column that tells us if the participant has regular contact with close friends
label(Results$SNI_6_close_friends) <- "regular contact with at least one close friend: 0 = no; 1 = yes"
Results <- cbind(Results, SNI_7_religious_community = as.numeric(RawData$SNI_14 > 1)) # add column that tells us if the participant has regular contact with members of a common religion
label(Results$SNI_7_religious_community) <- "regular contact to member(s) of the same religion: 0 = no; 1 = yes"
Results$SNI_7_religious_community[is.na(Results$SNI_7_religious_community)] <- 0
Results <- cbind(Results, SNI_8_class_mates = as.numeric(RawData$SNI_16 > 1)) # add column that tells us if the participant has regular contact to members of a common class (school or college)
label(Results$SNI_8_class_mates) <- "regular contact with class mate(s): 0 = no; 1 = yes"
Results$SNI_8_class_mates[is.na(Results$SNI_8_class_mates)] <- 0
Results <- cbind(Results, SNI_9_colleagues = as.numeric(RawData$SNI_19 > 1)) # add column that tells us if the participant has regular contact with colleagues 
label(Results$SNI_9_colleagues) <- "regular contact with with collaegue(s): 0 = no; 1 = yes" 
Results$SNI_9_colleagues[is.na(Results$SNI_9_colleagues)] <- 0
Results <- cbind(Results, SNI_10_neighbors = as.numeric(RawData$SNI_20 > 1 )) # add column that tells us if the participant has regular contact to his/her neighbours
label(Results$SNI_10_neighbors) <- "regular contact with neighbor(s): 0 = no; 1 = yes"
Results <- cbind(Results, SNI_11_fellow_volunteers = as.numeric(RawData$SNI_22 > 1)) # add column that tells us if the participant has regular contact to colleagues from a common charity group
label(Results$SNI_11_fellow_volunteers) <- "regular contact with fellow volunteer(s): 0 = no; 1 = yes"
Results$SNI_11_fellow_volunteers[is.na(Results$SNI_11_fellow_volunteers)] <- 0
Results <- cbind(Results, SNI_12_other_groups = as.numeric(RawData$SNI_23 == 2)) # add column that tells us if the participant has regular contact to members from groups other than the ones mentioned above
label(Results$SNI_12_other_groups) <- "regular contact to member(s) of any other kind of group: 0 = no; 1 = yes"
Results$SNI_12_other_groups[is.na(Results$SNI_12_other_groups)] <- 0

Results <- cbind(Results, SNI_Social_network_diversity = rowSums(Results[grepl("SNI", names(Results))])) # calculate the social network diversity of the participants (= the participants' numbers of active social roles) by calculating the rowSums of SNI_1 to SNI_12
label(Results$SNI_Social_network_diversity) <- "Number of different social network domains in which the participant is active; ranging from 0 (no social network) to 12 (very diverse social network)"

temp_SNI_3_parents <- vapply(RawData$SNI_6, function(x) {if (x == 2 | x == 3) 1 else if (x == 4) 2 else 0}, numeric(1)) # calculate the number of regular social interactions with parents (0 - 2) and store it in temp_SNI_3_Eltern
temp_SNI_4_in_laws <- vapply(RawData$SNI_8, function(x) {if (x == 2 | x == 3) 1 else if (x == 4) 2 else 0}, numeric(1)) # calculate the number of regular social interactions with parents-in-law (0 - 2) and store it in temp_SNI_4_Schwiegereltern
temp_SNI_2_and_5_to_12 <- RawData[c("SNI_4", "SNI_10", "SNI_12", "SNI_14", "SNI_16", "SNI_18", "SNI_19", "SNI_20", "SNI_22", "SNI_25")] # create a new data frame with the SNI roles 2 & 5 to 12 and store it in temp_SNI_2_und_5bis12
temp_SNI_2_and_5_to_12[is.na(temp_SNI_2_and_5_to_12)] <- 1 # NA becomes 1 (meaning 0)

Results <- cbind(Results, SNI_Social_network_size = rowSums(temp_SNI_2_and_5_to_12 -1) + temp_SNI_4_in_laws + temp_SNI_3_parents + Results$SNI_1_partner) # calculate the social network size of each participant: add temp_SNI_3_Eltern, temp_SNI_4_Schwiegereltern, Results$SNI_1_Partner and the rowSums of (temp_SNI_2_und_5bis12 - 1) and store it in a new column titled "SNI_Soziales_Netzwerk_Groesse"
label(Results$SNI_Social_network_size) <- "Number of different social interactions the participant has per week"


### evaluation of the HEXACO-60
temp_df_Hex_HH <- data.frame(mapply(Umpolung, RawData[c("Hex_30", "Hex_12", "Hex_60", "Hex_42", "Hex_24", "Hex_48")]), RawData[c("Hex_6", "Hex_54", "Hex_36", "Hex_18")]) # reverse the inversed items of the factor "Honesty-Humility" and put them - together with the non-inversed items of that factor - into temp_df_Hex_HH
Results <- cbind(Results, HEXACO_Honesty_Humility = round(apply(temp_df_Hex_HH, 1, mean), digits = 1)) # calculate the mean of the items in temp_df_Hex_HH and add it in a new column titled "HEXACO_Honesty_Humility" to Results

temp_df_Hex_Emotion <- data.frame(mapply(Umpolung, RawData[c("Hex_53", "Hex_35", "Hex_41", "Hex_59")]), RawData[c("Hex_5", "Hex_29", "Hex_11", "Hex_17", "Hex_23", "Hex_47")]) # create data frame temp_df_Hex_Emotion with all the (reversed) items of the category "Emotionality"
Results <- cbind(Results, HEXACO_Emotionality = round(apply(temp_df_Hex_Emotion, 1, mean), digits = 1)) # calculate the mean of the items in temp_df_Hex_Emotion and add it in a new column titled "HEXACO_Emotionality" to Results

temp_df_Hex_Extr <- data.frame(mapply(Umpolung, RawData[c("Hex_28", "Hex_52", "Hex_10", "Hex_46")]), RawData[c("Hex_4", "Hex_34", "Hex_58", "Hex_16", "Hex_40", "Hex_22")])  # create data frame temp_df_Hex_Extr with all the (reversed) items of the factor "Extraversion"
Results <- cbind(Results, HEXACO_Extraversion = round(apply(temp_df_Hex_Extr, 1, mean), digits = 1)) # calculate the mean of the items in temp_df_Hex_Extr and add it to Results in a new column titled "HEXACO_Extraversion"

temp_df_Hex_Agree <- data.frame(mapply(Umpolung, RawData[c("Hex_9", "Hex_15", "Hex_57", "Hex_21")]), RawData[c("Hex_3", "Hex_27", "Hex_33", "Hex_51", "Hex_39", "Hex_45")]) # reverse the inversed items of the factor "Agreeableness" and store them together with the non-inversed ones in temp_df_Hex_Agree
Results <- cbind(Results, HEXACO_Agreeableness = round(apply(temp_df_Hex_Agree, 1, mean), digits = 1)) # calculate the mean of the "Agreeableness" items and store it in a new column titled "HEXACO_Agreeableness"

temp_df_Hex_Consc <- data.frame(mapply(Umpolung, RawData[c("Hex_26", "Hex_32", "Hex_14", "Hex_20", "Hex_44", "Hex_56")]), RawData[c("Hex_2", "Hex_8", "Hex_38", "Hex_50")]) # create data frame temp_df_Hex_Consc with all the (reversed) items of the factor "Conscientiousness" in it
Results <- cbind(Results, HEXACO_Conscientiousness = round(apply(temp_df_Hex_Consc, 1, mean), digits = 1)) # calculate the mean of the "Conscientiousness" items and store it in new column titled "HEXACO_Conscientiousness"

temp_df_Hex_Openness <- data.frame(mapply(Umpolung, RawData[c("Hex_1", "Hex_31", "Hex_49", "Hex_19", "Hex_55")]), RawData[c("Hex_25", "Hex_7", "Hex_13", "Hex_37", "Hex_43")]) # store the (reversed) items of the factor "Openness to Experience" in temp_df_Hex_Openness
Results <- cbind(Results, HEXACO_Openness = round(apply(temp_df_Hex_Openness, 1 ,mean), digits = 1)) # calculate the mean of the "Openness to Experience" items and store it in new column "HEXACO_Openness"

# label the columns containing the HEXACO results
label(Results$HEXACO_Honesty_Humility) <- "Ranging from 1 (little pronounced personality trait in participant) to 5 (very pronounced personality trait)"
label(Results$HEXACO_Emotionality) <- "Ranging from 1 (little pronounced personality trait in participant) to 5 (very pronounced personality trait)"
label(Results$HEXACO_Extraversion) <- "Ranging from 1 (very introverted: ) to 5 (very pronounced personality trait)"
label(Results$HEXACO_Agreeableness) <- "Ranging from 1 (very little agreeableness: holds grudges, does not forgive others' shortcomings, easily angered) to 5 (very high agreeableness: forgives the wrongs they suffered easily, lenient in judging others, willing to compromise & cooperate, can easily control their temper)"
label(Results$HEXACO_Conscientiousness) <- "Ranging from 1 (little pronounced personality trait in participant) to 5 (very pronounced personality trait)"
label(Results$HEXACO_Openness) <- "Ranging from 1 (little pronounced personality trait in participant) to 5 (very pronounced personality trait)"



#### add data about NSSI

## data about NSSI history
Results <- cbind(Results, age_first_NSSI = SITBIG$age_firsttime)
label(Results$age_first_NSSI) <- "How old was the participant the first time they engaged in NSSI?"

Results <- cbind(Results, age_last_NSSI = SITBIG$age_lasttime)
label(Results$age_last_NSSI) <- "How old was the participant the last time they engaged in NSSI?"

Results <- cbind(Results, lifetime_NSSI = SITBIG$est_NSSI_lifetime)
label(Results$lifetime_NSSI) <- "How often did the participant self-injure in their entire life (estimate)?"

Results <- cbind(Results, NSSI_last_year = as.numeric(SITBIG$est_NSSI_lastyear))
label(Results$NSSI_last_year) <- "How many occasions of NSSI in the previous year?"

Results <- cbind(Results, RawData["SV_1"])
names(Results)[names(Results) == "SV_1"] <- "NSSI_last_3_months_1" #roughly how many NSSI acts in the past 3 months? (assessed with given levels of frequency)
label(Results$NSSI_last_3_months_1) <- "NSSI occured on: 1 = 1-4 days; 2 = 5-20 days; 3 = > 20 days ... over the previous 3 months"

Results <- cbind(Results, RawData["SV_2"]) 
names(Results)[names(Results) == "SV_2"] <- "NSSI_last_3_months_2" #precisely how many NSSI acts in the past 3 months? (estimate of number of days on which participant self-injured)
label(Results$NSSI_last_3_months_2) <- "estimated number of days over the previous 3 months on which NSSI occured"

Results <- cbind(Results, RawData["SV_3"])
names(Results)[names(Results) == "SV_3"] <- "NSSI_regularity_3_months" #how regularly did the participants self-injury in the past 3 months?
label(Results$NSSI_regularity_3_months) <- "1 = < once/month; 2 =  at least once/month but less then once/week; 3 = multiple times per week; 4 = daily"

Results <- cbind(Results, NSSI_last_month = as.numeric(sub(",", ".", SITBIG$est_NSSI_lastmonth, fixed = TRUE)))
label(Results$NSSI_last_month) <- "How many occasions of NSSI in the previous month?"

Results <- cbind(Results, NSSI_last_week = as.numeric(sub(",", ".", SITBIG$est_NSSI_lastweek, fixed = TRUE)))
label(Results$NSSI_last_week) <- "How many occasions of NSSI in the previous week?"

## data about NSSI methods, locations & aims

# change RawData labels to english
label(RawData["SV_4"]) <- as.list(c(SV_4 = "scratching the skin with a sharp object, for example with a compass"))
label(RawData["SV_5"]) <- as.list(c(SV_5 = "scratching the skin without an object, for example using only the nails"))
label(RawData["SV_6"]) <- as.list(c(SV_6 = "cutting the skin with a sharp object, for example a knife or razor blade"))
label(RawData["SV_7"]) <- as.list(c(SV_7 = "hitting, for example hitting yourself with your hands, banging your head on purpose"))
label(RawData["SV_8"]) <- as.list(c(SV_8 = "damaging the the finger- or footnails / the nailbed"))
label(RawData["SV_9"]) <- as.list(c(SV_9 = "pulling out hair"))
label(RawData["SV_10"]) <- as.list(c(SV_10 = "pinching oneself"))
label(RawData["SV_11"]) <- as.list(c(SV_11 = "bitting oneself"))
label(RawData["SV_12"]) <- as.list(c(SV_12 = "burning, for example using a cigarette"))
label(RawData["SV_13"]) <- as.list(c(SV_13 = "scalding, for example using hot water"))
label(RawData["SV_14"]) <- as.list(c(SV_14 = "chafing one's own skin"))
label(RawData["SV_15"]) <- as.list(c(SV_15 = "wound manipulation"))
label(RawData["SV_16"]) <- as.list(c(SV_16 = "tattooing oneself"))
label(RawData["SV_17"]) <- as.list(c(SV_17 = "choking oneself"))
label(RawData["SV_18"]) <- as.list(c(SV_18 = "breaking one's own bones (on purpose)"))
label(RawData["SV_19"]) <- as.list(c(SV_19 = "(oral) intake of corrosive chemicals/toxins"))
label(RawData["SV_20"]) <- as.list(c(SV_20 = "swallowing sharp objects"))
label(RawData["SV_21"]) <- as.list(c(SV_21 = "putting objects under the skin"))
label(RawData["SV_22"]) <- as.list(c(SV_22 = "piercing body parts"))
label(RawData["SV_23"]) <- as.list(c(SV_23 = "destruction of the eye"))
label(RawData["SV_24"]) <- as.list(c(SV_24 = "amputating or castrating oneself"))
label(RawData["SV_25"]) <- as.list(c(SV_25 = "other (please describe"))
label(RawData["SV_26"]) <- as.list(c(SV_26 = "upper arm"))
label(RawData["SV_27"]) <- as.list(c(SV_27 = "lower arm"))
label(RawData["SV_28"]) <- as.list(c(SV_28 = "thigh"))
label(RawData["SV_29"]) <- as.list(c(SV_29 = "calves / ankles"))
label(RawData["SV_30"]) <- as.list(c(SV_30 = "foot / toes"))
label(RawData["SV_31"]) <- as.list(c(SV_31 = "inside of the wrist"))
label(RawData["SV_32"]) <- as.list(c(SV_32 = "hand"))
label(RawData["SV_33"]) <- as.list(c(SV_33 = "fingers"))
label(RawData["SV_34"]) <- as.list(c(SV_34 = "chest / cleavage"))
label(RawData["SV_35"]) <- as.list(c(SV_35 = "stomach"))
label(RawData["SV_36"]) <- as.list(c(SV_36 = "back"))
label(RawData["SV_37"]) <- as.list(c(SV_37 = "genitals"))
label(RawData["SV_38"]) <- as.list(c(SV_38 = "anus"))
label(RawData["SV_39"]) <- as.list(c(SV_39 = "face"))
label(RawData["SV_40"]) <- as.list(c(SV_40 = "head"))
label(RawData["SV_41"]) <- as.list(c(SV_41 = "neck / throat"))
label(RawData["SV_45"]) <- as.list(c(SV_45 = "to avoid school, work or other activities"))
label(RawData["SV_46"]) <- as.list(c(SV_46 = "to get attention or a reaction from one's partner or friends"))
label(RawData["SV_47"]) <- as.list(c(SV_47 = "to avoid doing sth unpleasant that one doesn't want to do"))
label(RawData["SV_48"]) <- as.list(c(SV_48 = "to avoid being among people"))
label(RawData["SV_49"]) <- as.list(c(SV_49 = "to make people act differently or change"))
label(RawData["SV_50"]) <- as.list(c(SV_50 = "to be like someone one admires"))
label(RawData["SV_51"]) <- as.list(c(SV_51 = "to avoid punishment or consequences"))
label(RawData["SV_52"]) <- as.list(c(SV_52 = "to feel something instead of just emptiness, even pain"))
label(RawData["SV_53"]) <- as.list(c(SV_53 = "to punish oneself"))
label(RawData["SV_54"]) <- as.list(c(SV_54 = "to regulate emotions"))
label(RawData["SV_55"]) <- as.list(c(SV_55 = "to reduce stress"))
label(RawData["SV_56"]) <- as.list(c(SV_56 = "out of self-hatred"))
label(RawData["SV_57"]) <- as.list(c(SV_57 = "out of shame or guilt"))

# data about NSSI methods
NSSI_method_sorted_no_cutoff <- sort(colSums(RawData[grep("SV_[4-9]$|SV_1[0-9]|SV_2[0-4]", names(RawData))]), decreasing = TRUE) # sort the NSSI methods (SV_4 to SV_24) by decreasing column sum
NSSI_method_sorted_with_cutoff <- NSSI_method_sorted_no_cutoff[NSSI_method_sorted_no_cutoff > nrow(Results) * 1.5] # don't keep the columns of NSSI methods which participants rather disagree to using, meaning the mean of that column is smaller than 1.5
Results <- cbind(Results, RawData[names(NSSI_method_sorted_with_cutoff)]) # add columns of NSSI methods which participants rather agree to using (mean is greater than 1.5)
Results <- Results %>% rename_at(vars(names(NSSI_method_sorted_with_cutoff)), ~c("NSSI_method_1","NSSI_method_2","NSSI_method_3","NSSI_method_4","NSSI_method_5","NSSI_method_6","NSSI_method_7")) # rename the method columns according to sinking agreement to using that method (method_1 = the method most frequently agreed to using)

# data about NSSI locations
NSSI_loc_sorted_no_cutoff <- sort(colSums(RawData[grep("SV_2[6-9]|SV_3[0-9]|SV_4[01]", names(RawData))]), decreasing = TRUE) # sort the NSSI locations (SV_26 to SV_41) by decreasing column sum
NSSI_loc_sorted_with_cutoff <- NSSI_loc_sorted_no_cutoff[NSSI_loc_sorted_no_cutoff > nrow(Results) * 1.5] # don't keep the columns of NSSI locations which participants rather disagree to injuring at (colsum <= 1.5 * n(participants))
Results <- cbind(Results, RawData[names(NSSI_loc_sorted_with_cutoff)]) # add columns of NSSI locations which participants rather agree to self-injuring at to results
Results <- Results %>% rename_at(vars(names(NSSI_loc_sorted_with_cutoff)), ~c("NSSI_loc_1","NSSI_loc_2","NSSI_loc_3","NSSI_loc_4","NSSI_loc_5","NSSI_loc_6")) # rename the location columns according to decreasing agreement rank (loc_1 = the location most commonly injured at, loc_2 = the second most common, etc.)

# data about NSSI aim
NSSI_aim_sorted_no_cutoff <- sort(colSums(RawData[grep("SV_4[5-9]|SV_5[0-7]", names(RawData))]), decreasing = TRUE) # sort the NSSI aims (=subjective reasons for why participants self-injure) by decreasing column sum
NSSI_aim_sorted_with_cutoff <- NSSI_aim_sorted_no_cutoff[NSSI_aim_sorted_no_cutoff > nrow(Results) * 1.5] # don't keep the columns of NSSI aims which participants rather disagree to 
Results <- cbind(Results, RawData[names(NSSI_aim_sorted_with_cutoff)]) # add columns of NSSI aims which participants rather agree to to Results
Results <- Results %>% rename_at(vars(names(NSSI_aim_sorted_with_cutoff)), ~c("NSSI_aim_1","NSSI_aim_2","NSSI_aim_3","NSSI_aim_4","NSSI_aim_5","NSSI_aim_6","NSSI_aim_7")) # rename the NSSI aim columns according to decreasing agreement rank 

## further description of NSSI
Results <- cbind(Results, most_severe_NSSI_3_months = RawData$SV_42) # add data about how grave the most serious self-injury in the past three months was (SV_42)
label(Results$most_severe_NSSI_3_months) <- "How grave was the most severe self-injury in the past 3 months?: 1 = superficial; 2 = moderate; 3 = severe; e.g. 1 = superficial scratches, scrapes, bruises; 2 = strongly bleeding cuts, 2nd degree burns; 3 = cuts down to the fat tissue, bone fractures, injured tendons, internal bleeding"

Results <- cbind(Results, med_treatm_for_NSSI_3_months = (RawData$SV_43 -1)) # add column that specifies what kind of medical treatment was needed / done for said injury (SV_43)
label(Results$med_treatm_for_NSSI_3_months) <- "What kind of medical treatment was needed/performed for this most severe self-injury: 0 = none; 1 = participant treated the injury themselves; 2 = medical/surgical treatment was performed (by a med. professional); 3 = medical/surgical treatment would've been needed, but was not performed"

Results <- cbind(Results, scars_through_NSSI = (RawData$SV_44 -1)) # add whether participants have scars caused by self-injury (SV-44)
label(Results$scars_through_NSSI) <- "Lasting scars because of NSSI?: 0 = no; 1 = yes"

Results <- cbind(Results, planed_NSSI = (RawData$SV_58 - 1)) # add column that specifies whether participants normally plan their self-injuries...(SV_58)
label(Results$planed_NSSI) <- "Does the participant plan NSSI before engaging in it?: 0 = never; 1 = rarely; 2 = often; 3 = most of the time"

Results <- cbind(Results, impulsive_NSSI = (RawData$SV_59 - 1)) # ...or whether they act impulsively (SV_59)
label(Results$impulsive_NSSI) <- "Does the participant engage in NSSI impulsively?: 0 = never; 1 = rarely; 2 = often; 3 = most of the time"

Results <- cbind(Results, worse_NSSI_than_intended = (RawData$SV_60 - 1)) # did they ever hurt themselves worse than intended? (SV_60)
label(Results$worse_NSSI_than_intended) <- "Was the performed self-injury ever more severe than intended?: 0 = never; 1 = rarely; 2 = often; 3 = most of the time"

Results <- cbind(Results, mean_NSSI_urge = as.numeric(sub(",", ".", SITBIG$Urge_average_intensity, fixed = TRUE)))
label(Results$mean_NSSI_urge) <- "Mean urge to self-injury on a scale from 0 (very weak) to 4 (very strong)"

Results <- cbind(Results, urge_NSSI_3_months = RawData$SV_61) # add column that illuminates the participants' urge to injure themselves over the past three months (SV_61)
label(Results$urge_NSSI_3_months) <- "The urge to engage in NSSI was... in the previous 3 months: 1 = very weak; 2 = weak; 3 = strong; 4 = very strong"

Results <- cbind(Results, NSSI_is_problem = RawData$SV_prob) # add whether participants feel that NSSI is a problem... (SV_prob)
label1 <- "1 = strong disagreement; 2 = disagreement; 3 = so-so; 4 = agreement; 5 = strong agreement"
label(Results$NSSI_is_problem) <- label1

Results <- cbind(Results, wish_to_stop_NSSI = RawData$SV_aufh) #...and whether they want to stop injuring themselves (SV_aufh)
label(Results$wish_to_stop_NSSI) <- label1

Results <- cbind(Results, more_suic._thoughts_w/o_NSSI = RawData$SV_ged) # do the participants feel as if their suicidal thoughts are more severe if they don't self-injure? (SV_ged)
label(Results$more_suic._thoughts_w/o_NSSI) <- label1

## data about pain & NSSI
Results <- cbind(Results, PAIN_in_general = RawData$SVpain_1 - 1) # add column with data about pain perception in general (SVpain_1)
label2 <- "Pain ranging from 0 (no pain at all) to 10 (worst imaginable pain)"
label(Results$PAIN_in_general) <- label2 
Results <- cbind(Results, PAIN_last_NSSI = RawData$SVpain_2 - 1) # add column with data about pain perception during most recent self-injury  (SVpain_2)
label(Results$PAIN_last_NSSI) <- label2
Results <- cbind(Results, PAIN_first_NSSI = RawData$SVpain_3 -1) # add column with data about pain perception during the participant's first self-injury ever (SVpain_3)
label(Results$PAIN_first_NSSI) <- label2
Results <- cbind(Results, PAIN_during_NSSI = RawData$SVpain_4 - 1)  # add column with data about pain perception during self-injury in general (SVpain_4)
label(Results$PAIN_during_NSSI) <- label2
Results <- cbind(Results, PAIN_minutes_after_NSSI = RawData$SVpain_5 - 1) # add column with data about pain a few minutes after having self-injured (SVpain_5)
label(Results$PAIN_minutes_after_NSSI) <- label2
Results <- cbind(Results, time_until_normal_PAIN_perception_after_NSSI = RawData$SVpain_6) # add column with data about how long it takes for the participant's pain perception to return to normal after self-injuring (SVpain_6)
label(Results$time_until_normal_PAIN_perception_after_NSSI) <- "1 = <= 10min; 2 = 10-30min; 3 = 30-60min; 4 = 1-24h; 5 = >24h; 6 = pain perception isn't altered during NSSI"
Results <- cbind(Results, changed_PAIN_perception_since_first_NSSI = RawData$SVpain_7) # add column with data about how the participant's pain perception changed since the first self-injury (SVpain_7)
label(Results$changed_PAIN_perception_since_first_NSSI) <- "1 = pain is much weaker now; 2 = ...is a bit weaker; 3 = ...the same; 4 = ...a bit stronger; 5 = ...much stronger"

## data about suicide ideation / suicide plans / suicide attempts
Results <- cbind(Results, lifetime_SUICIDAL_thoughts_plans_attempts = RawData$Suizid_1) # add a column that specifies whether the participant ever thought about / tried to kill him/herself (Suizid_1)
label(Results$lifetime_SUICIDAL_thoughts_plans_attempts) <- "Did the participant ever think about / plan / attempt to kill themselves?: 1 = never; 2 = thought about it; 3 = planed it; 4 = had a plan and wanted to act on it; 5 = suicide attempt w/o really wanting to die, 6 = suicide attempt with intention to die"
Results <- cbind(Results, n_SUICIDE_attempts = RawData$Suizid_2) # add a column with data about how many suicide attempts the participant made (Suizid_2)
label(Results$n_SUICIDE_attempts) <- "Lifetime count of suicide attempts"
Results <- cbind(Results, last_SUICIDE_attempt = RawData$Suizid_3) # add a column that tells us when the last suicide attempt was
label(Results$last_SUICIDE_attempt) <- "1 = in the previous year; 2 = 1 - 2 years ago; 3 = over 2 years ago"



  

##### ADD & EVALUATE DATA ABOUT THE HEALTHY CONTROL GROUP #####
Results_HC <- Results[0,] # create new dataframe "Results_HC" with the same columns as the Results dataframe with the data about the NSSI group
Results_HC[1:nrow(CTQ_HC),]<- NA # set the number of rows of Results_HC to the number of rows of CTQ_HC and fill in everything with NA by default
Results_HC["Sheet_Nr"] <- CTQ_HC["Bogen"] #add the HC group participants' sheet & barcode numbers to Results_HC
Results_HC["BARCODE"] <- CTQ_HC["BARCODE"]
Results_HC["Age"] <- basic_data_HC["Alter"] #add age
Results_HC["Weigth"] <- basic_data_HC["Gewicht"] #add weight
Results_HC["Heigth"] <- basic_data_HC["Größe"]*100 #add height (in cm)
Results_HC["BMI"] <- basic_data_HC["BMI"] #add BMI
Results_HC["Sex"] <- rep(1, length(Results_HC$Geschlecht)) #add gender
Results_HC["Smoking"] <- basic_data_HC["Raucher"] #add smoking habits
Results_HC$Smoking[grepl("nein", Results_HC$Smoking)] <- 0 # change "no" in the "smoker?" column to 0
Results_HC$Smoking[grepl("^j", Results_HC$Smoking)] <- 1 # change "yes" in the "smoker?" column to 1

## add column to results2 that specifies the group status as "0 = HC
group_status_HC <- rep("0", length(Results_HC$BARCODE)) #add group status
Results_HC["Group"] <- group_status_HC # add the group status of HC as 0 = HC in the column 'group'

## add labcodes of the HC_group
Results_HC[match(BARCODE_to_labcode_HC$BARCODE, Results_HC$BARCODE), "Lab-Code"] <- BARCODE_to_labcode_HC$`Lab-Code`

### evaluation of the CTQ of the HC group

Results_HC <- ctq_rowsum(CTQHC_list, CTQ_HC, Results_HC)

#Results_HC["CTQ_EA"] <- rowSums(CTQ_HC[c("CTQ00300", "CTQ00800", "CTQ01400", "CTQ01800", "CTQ02500")]) #score of the CTQ subscale emotional abuse
#Results_HC["CTQ_PA"] <- rowSums(CTQ_HC[c("CTQ00900", "CTQ01100", "CTQ01200", "CTQ01500", "CTQ01700")]) #score of the CTQ subscale physical abuse
#Results_HC["CTQ_SA"] <- rowSums(CTQ_HC[c("CTQ02000", "CTQ02100", "CTQ02300", "CTQ02400", "CTQ02700")]) #score of the CTQ subscale sexual abuse
#Results_HC["CTQ_EN"] <- rowSums(mapply(Umpolung, CTQ_HC[c("CTQ00500", "CTQ00700", "CTQ01300", "CTQ01900", "CTQ02800")])) #score of the CTQ subscale emotional neglect
#Results_HC["CTQ_PN"] <- (rowSums(mapply(Umpolung, CTQ_HC[c("CTQ00200", "CTQ02600")])) + rowSums(CTQ_HC[c("CTQ00100", "CTQ00400", "CTQ00600")])) #score of the CTQ subscale physical neglect
Results_HC["CTQ_total_score"] <- rowSums(Results_HC[c("CTQ_EA", "CTQ_PA", "CTQ_SA", "CTQ_EN", "CTQ_PN")]) #overall score of the CTQ (=sum of the subscores)



###### MERGE THE NSSI DATAFRAME (=RESULTS) WITH THE HC DATAFRAME (=RESULTS_HC) INTO A NEW DATAFRAME TITLED "COMPLETE_RESULTS_no_bE" #####
Complete_Results_no_Prot <- rbind(Results, zap_labels(Results_HC))

### add the beta Endorphin values of both groups to complete_Results
#Complete_Results <- merge(Complete_Results_no_bE, all_beta_Endorphin[, c("BARCODE", "bE_1.cm", "bE_2.cm")], by="BARCODE")
#label(Complete_Results$bE_1.cm) <- "beta-Endorphin in ng im 0-1cm-Segment"
#label(Complete_Results$bE_2.cm) <- "beta-Endorphin in ng im 1-2cm-Segment"


### PROTEOMICS RESULTS
## create new data frame titled 'proteomics_results' and load barcodes & labcodes into it
# Proteomics_results <- BARCODE_to_labcode_complete[c("BARCODE", "Lab-Code")]

## add columns with all proteins to 'proteomics_results', avoiding doubles
# [,sort(unique(bind_rows(mget(ls(pattern = "^SV.*")))$Description))] <- as.double(NA)

## fill the data sheet with all proteins' abundance scores
# for (row in 1:nrow(Proteomics_results)){
#  proteindf = paste0(Proteomics_results[row,"Lab-Code"], "_proteomics")
#  temp <- tryCatch(
#    {
#      get(proteindf)[, c("Description", grep("Abundance.*", names(get(proteindf)), value = TRUE))]
#    },
#    error=function(e){
#      message(e, appendLF = TRUE)
#      e
#    }
#  )
#  if(inherits(temp, "error")) next
#  for (row2 in 1:nrow(temp)){
#    Proteomics_results[row, as.character(temp[row2,1])] <- as.double(temp[row2,2])
#  }
#}

## calculate for each protein the amount of samples in which it was found and at the results in a new row at the bottom of the data sheet
# nrow_Proteomics_results_sum <- nrow(BARCODE_to_labcode_complete) + 1
# for (col in 3:ncol(Proteomics_results)){
#  temp <- sum(!is.na(Proteomics_results[1:nrow_Proteomics_results_sum - 1,col]))
#  Proteomics_results[nrow_Proteomics_results_sum, col] <- temp
# }

## calculate the percentage of samples that each protein was found in
# all samples, looking at abundance scores
SV_protein_percentage_1 <- SV_proteomics[c(1:1114), c(1:73,150,151,152)]
Amount_1 = rowSums(SV_protein_percentage_1[c(1:1114),4:73] !="", na.rm = TRUE)
Percentage_1 = round(rowSums(SV_protein_percentage_1[c(1:1114),4:73] !="", na.rm = TRUE) / 70 * 100, digits=2)

# all samples, looking at abundance counts (just a test, should be the same result as abundance scores!)
SV_protein_percentage_2 <- SV_proteomics[c(1:1114), c(1,2,3,74:143)]
Amount_2 = rowSums(SV_protein_percentage_2[c(1:1114),4:73] !="", na.rm = TRUE)
Percentage_2 = round(rowSums(SV_protein_percentage_2[c(1:1114),4:73] !="", na.rm = TRUE) / 70 * 100, digits=2)

# only NSSI samples, looking at abundance scores
SV_protein_percentage_NSSI <- SV_proteomics[c(1:1114), c(4:37, 71, 72)]
Amount_NSSI = rowSums(SV_protein_percentage_NSSI[c(1:1114),] !="", na.rm = TRUE)
Percentage_NSSI = round(rowSums(SV_protein_percentage_NSSI[c(1:1114),] !="", na.rm = TRUE) / 36 * 100, digits=2)

# only HC samples, looking at abundance scores
SV_protein_percentage_HC <- SV_proteomics[c(1:1114), c(38:70, 73)]
Amount_HC = rowSums(SV_protein_percentage_HC[c(1:1114),] !="", na.rm = TRUE)
Percentage_HC = round(rowSums(SV_protein_percentage_HC[c(1:1114),] !="", na.rm = TRUE) / 34 * 100, digits=2)

# add columns with the calculated amounts & percentages for both groups as well as each group separately 
SV_protein_percentage_1 <- cbind(SV_protein_percentage_1, Amount_1) #Amount 1 = in how many of the 70 hair samples (both groups) was the protein found? (going by abundance score)
SV_protein_percentage_1 <- cbind(SV_protein_percentage_1, Percentage_1) #Percentage 1 = in which percentage of all hair samples (both groups) was the protein found?(going by abundance score)
SV_protein_percentage_1 <- cbind(SV_protein_percentage_1, Amount_2) #Amount 2 = in how many of the 70 hair samples (both groups) was the protein found? (going by abundance count)
SV_protein_percentage_1 <- cbind(SV_protein_percentage_1, Percentage_2) #Percentage 2 = in which percentage of all hair samples (both groups) was the protein found?(going by abundance count)
SV_protein_percentage_1 <- cbind(SV_protein_percentage_1, Amount_NSSI) #Amount_NSSI = in how many of the 36 NSSI hair samples was the protein found? (going by abundance score)
SV_protein_percentage_1 <- cbind(SV_protein_percentage_1, Percentage_NSSI) #Percentage_NSSI = in which percentage of the 36 NSSI hair samples was the protein found?(going by abundance score)
SV_protein_percentage_1 <- cbind(SV_protein_percentage_1, Amount_HC) #Amount_HC = in how many of the 34 HC hair samples was the protein found? (going by abundance score)
SV_protein_percentage_1 <- cbind(SV_protein_percentage_1, Percentage_HC) #Percentage_HC = in which percentage of the 34 HC hair samples was the protein found?(going by abundance score)

## subset the data frame containing the percentages, removing all proteins that aren't found in at least 50% of each group's samples
temp_percentage_at_least_50 <- SV_protein_percentage_1[SV_protein_percentage_1$Percentage_HC>=50,] #remove all proteins that aren't found in at least 50% of the HC samples, store results in temporary data frame
percentage_at_least_50 <- temp_percentage_at_least_50[temp_percentage_at_least_50$Percentage_NSSI>=50,] #remove all proteins from temporary data frame that aren't found in at least 50% of the NSSI samples
View(percentage_at_least_50) #data frame "percentage_at_least_50" now containts only the proteins that are found in at least 50% of each group's samples

##### STATISTICAL ANALYSIS #####

### have a look at the beta endorphin data
#Histogram_bE <- ggplot(Complete_Results, aes(x = bE_1.cm)) %>%
#  + geom_histogram() %>% # create a histogram of the beta endorphin values    ### FEINTUNING DES AUSSEHENS, SOBALD DIE BETA ENDORPHIN DATEN DA SIND!!
#  + facet_wrap(~ Gruppe) # create separate subgraphs for both groups
#Histogram_bE

### check the NON-LOG beta endorphin for outliers & create a data frame without them   ### VERSION OHNE SKEWNESS DER BETA-ENDORPHIN WERTE
#bE1_0.25 <- quantile(Complete_Results$bE_1.cm, 0.25) # calculate the 0.25-quantile of the beta endorphin values and store it in "bE_0.25"
#bE1_0.75 <- quantile(Complete_Results$bE_1.cm, 0.75) # calculate the 0.75-quantile of the beta endorphin values and store it in "bE_0.75"
#bE1_IQA <- bE1_0.75 - bE1_0.25 # calculate the IQA of the beta endorphin values and store it in "bE_IQA"
#lower_threshold_bE1 <- bE1_0.25 - 1.5 * bE1_IQA # calculate the lower threshold of the beta endorphin values (=0.25-quantile - 1.5*IQA) and store it in "lower_threshold_bE"
#upper_threshold_bE1 <- bE1_0.75 + 1.5 * bE1_IQA # calculate the upper threshold of the beta endorphin values (=0.75-quantile + 1.5*IQA) and store it in "upper_threshold_bE"

#bE1_both_groups_no_outliers <- Complete_Results %>% 
#  filter(bE_1.cm > lower_threshold_bE1 & bE_1.cm < upper_threshold_bE1)



### check the beta endorphin data for skewness and possibly conduct log transformation to correct it
#skewness(Complete_Results$bE_1.cm) # check the beta endorphin data for skewness; if the skewness of the predictor variable is < -0.5 or > +0.5, a log-transformation is probably a good idea

### IF THE DATA IS SKEWED: have a look at the log-transformed beta endorphin values
#Histogram_log_bE1 <- ggplot(Complete_Results, aes(x = log(bE_1.cm))) %>%
#  + geom_histogram() %>% # create a histogram of the log-transformed beta endorphin values
#  + facet_wrap(~ Gruppe) # create separate subgraphs for both groups
#Histogram_log_bE1

### add a new column to Complete_Results that contains the log-transformed beta endorphin values
#Complete_Results <- cbind(Complete_Results, log_bE1 = log(Complete_Results$bE_1.cm))
#label(Complete_Results$log_bE1) <- "log transformierte β-Endorphin-Werte des 0-1cm-Segments"
  
### check log beta endorphin for outliers & create a data frame without them   ### FALLS KEINE SKEWNESS DER ORIGINALEN BETA ENDORPHIN WERTE BESTEHT WEGLASSEN
#log_bE1_0.25 <- quantile(Complete_Results$log_bE1, 0.25) # calculate the 0.25-quantile of the log beta endorphin values and store it in "bE_0.25"
#log_bE1_0.75 <- quantile(Complete_Results$log_bE1, 0.75) # calculate the 0.75-quantile of the log beta endorphin values and store it in "bE_0.75"
#log_bE1_IQA <- log_bE1_0.75 - log_bE1_0.25 # calculate the IQA of the log beta endorphin values and store it in "bE_IQA"
#lower_threshold_log_bE1 <- log_bE1_0.25 - 1.5 * log_bE1_IQA # calculate the lower threshold of the log beta endorphin values (=0.25-quantile - 1.5*IQA) and store it in "lower_threshold_bE"
#upper_threshold_log_bE1 <- log_bE1_0.75 + 1.5 * log_bE1_IQA # calculate the upper threshold of the log beta endorphin values (=0.75-quantile + 1.5*IQA) and store it in "upper_threshold_bE"

#log_bE1_both_groups_no_outliers <- Complete_Results %>% 
#  filter(log_bE1 > lower_threshold_log_bE1 & log_bE1 < upper_threshold_log_bE1) # filter Complete_Results so that a row (meaning a participant's data) is only kept if the log beta endorphin value isn't lower than the lower_threshold_bE or greater than the upper_threshold_bE & store the result in a new data frame titled "log_bE_both_groups_no_outliers" 



#### 1. T-Test

### conduct two-sided unpaired t-test for the beta-endorphin values of the two groups (using bE_both_groups_no_outliers, therefore ignoring beta endorphin outliers)
#t_test_beta_Endorphin <- t.test(bE_1.cm ~ Gruppe, data = bE1_both_groups_no_outliers, null = 0, alternative = "two.sided")
#t_test_beta_Endorphin


#### 2. Correlation

### correlation analysis of the strength of the relationship between beta endorphin and the NSSI frequency in the previous three months before taking the hair sample
#numeric_Hfgkt <- as.numeric(Results$NSSI_Hfgkt_im_letzten_Monat) # store the frequency of NSSI in the past month as numeric values in a vector titled "numeric_Hfgkt"
#Hfgkt_0.25 <- quantile(numeric_Hfgkt, 0.25) # calculate the 0.25-quantile of numeric_Hfgkt and store it in "Hfgkt_0.25"
#Hfgkt_0.75 <- quantile(numeric_Hfgkt, 0.75) # calculate the 0.75-quantile of numeric_Hfgkt and store it in "Hfgkt_0.75"
#Hfgkt_IQA <- Hfgkt_0.75 - Hfgkt_0.25 # calculate the IQA of numeric_Hfgkt and store it in "Hfgkt_IQA"

#lower_threshold_Hfgkt <- Hfgkt_0.25 - 1.5 * Hfgkt_IQA
#upper_threshold_Hfgkt <- Hfgkt_0.75 + 1.5 * Hfgkt_IQA

#numeric_Hfgkt_no_outliers <- log_bE1_both_groups_no_outliers %>%
#  filter(as.numeric(NSSI_Hfgkt_im_letzten_Monat) > lower_threshold_Hfgkt & as.numeric(NSSI_Hfgkt_im_letzten_Monat) < upper_threshold_Hfgkt) # filter the data frame without log beta endorphin outliers to get a new data frame without both beta endorphin and/or  outliers and title it "numeric_Hfgkt_no_Outliers"

#Scatterplot_bE1_Hfgkt <- ggplot(numeric_Hfgkt_no_outliers, aes(x = as.numeric(NSSI_Hfgkt_im_letzten_Monat), y = log_bE1)) %>%
#  + geom_point() %>% # create a scatterplot of the beta endorphin & NSSI frequency pairs   ### FEINTUNING DES AUSSEHENS DES SCATTERPLOTS, SOBALD DIE BETA-ENDORPHIN-DATEN DA SIND!!
#  + geom_smooth(method = "lm", se = FALSE) # add a trendline to the scatterplot
#Scatterplot_bE1_Hfgkt

#cor_bE1_Hfgkt <- cor(numeric_Hfgkt_no_outliers$NSSI_Hfgkt_im_letzten_Monat, numeric_Hfgkt_no_outliers$log_bE1, use = "pairwise.complete.obs", method = "pearson") # calculate the correlation between beta endorphin and NSSI frequency (ignoring outliers)
#cor_bE1_Hfgkt

### correlation analysis of the strength of the relationship between beta endorphin and the severity of NSSI in the previous three months before taking the hair sample
#Scatterplot_bE_Severity <- ggplot(log_bE1_both_groups_no_outliers, aes(x = NSSI_Max_Schweregrad_letzte_drei_Monate, y = log_bE1)) %>%
#  + geom_point() %>% # create a scatterplot of the beta endorphin & NSSI severity pairs   ### FEINTUNING!!!
#  + geom_smooth(method = "lm", se = FALSE) # add a trendline to the scatterplot
#Scatterplot_bE_Severity

#cor_bE_Severity <- cor(log_bE1_both_groups_no_outliers$NSSI_Max_Schweregrad_letzte_drei_Monate, log_bE_both_groups_no_outliers$log_beta_Endorphin, use = "pairwise.complete.obs", method = "pearson") # calculate the correlation between bE and NSSI severity   ### evtl Spearman statt Pearson?!
#cor_bE_Severity

### EVTL WEITERE KORRELATIONSANALYSEN, Z.B. KORRELATION ZWISCHEN BETA ENDORPHIN & SCHMERZEMPFINDEN WÄHREND SV (= Results$SCHMERZ_während_Durchführung_von_SV); BETA ENDORPHIN & SCHMERZ ALLGEMEIN (= Results$SCHMERZ_allgemein); BETA ENDORPHIN & SPORT (confounding!)



#### 3. CTQ & beta endorphin

### have a look at the CTQ data & check for correlation 
# Histogram_CTQ <- ggplot(Complete_Results, aes(x = CTQ_gesamt)) %>%
#   + geom_histogram(binwidth = 5) %>% # create a histogram of the CTQ results
#   + facet_wrap(~ Gruppe) # create separate subgraphs for both groups
# Histogram_CTQ
# 
# Scatterplot_bE_CTQ <- ggplot(bE_both_groups_no_outliers, aes(x = CTQ_gesamt, y = beta_Endorphin, color = Gruppe)) %>%
#   + geom_point() # create a scatterplot of the beta endorphin / CTQ total score pairs (both groups in one scatterplot, but with different coloured dots)
# Scatterplot_bE_CTQ
# 
# cor_bE_CTQ <- cor(bE_both_groups_no_outliers$CTQ_gesamt, bE_both_groups_no_outliers$beta_Endorphin, use = "pairwise.complete.obs", method = "pearson") # calculate the correlation between bE and the total score of the CTQ
# cor_bE_CTQ
# 
# ### test the second hypothesis via regression analysis by predicting the amount of beta-endorphin via the CTQ total score while controlling for NSSI onset, frequency over the past month and amount of exercise
# regression_CTQ_bE <- lm(bE_both_groups_no_outliers ~ CTQ_gesamt + Alter_bei_erstem_NSSI + NSSI_Hfgkt_im_letzten_Monat + Sport, data = Complete_Results)
# summary(regression_CTQ_bE)
# 
# ## test whether one of the CTQ subscales has a particular strong influence on the amount of beta endorphin
# regr_CTQ_EM_bE <- lm(bE_both_groups_no_outliers ~ CTQ_EM + Alter_bei_erstem_NSSI + NSSI_Hfgkt_im_letzten_Monat + Sport, data = Complete_Results)
# summary(regr_CTQ_EM_bE) #CTQ: Emotionale Misshandlung
# 
# regr_CTQ_KM_bE <- lm(bE_both_groups_no_outliers ~ CTQ_KM + Alter_bei_erstem_NSSI + NSSI_Hfgkt_im_letzten_Monat + Sport, data = Complete_Results)
# summary(regr_CTQ_KM_bE) #CTQ: Körperliche Misshandlung
# 
# regr_CTQ_SM_bE <- lm(bE_both_groups_no_outliers ~ CTQ_SM + Alter_bei_erstem_NSSI + NSSI_Hfgkt_im_letzten_Monat + Sport, data = Complete_Results)
# summary(regr_CTQ_SM_bE) #CTQ: Sexuelle Misshandlung
# 
# regr_CTQ_EV_bE <- lm(bE_both_groups_no_outliers ~ CTQ_EV + Alter_bei_erstem_NSSI + NSSI_Hfgkt_im_letzten_Monat + Sport, data = Complete_Results)
# summary(regr_CTQ_EV_bE) #CTQ: Emotionale Vernachlässigung
# 
# regr_CTQ_KV_bE <- lm(bE_both_groups_no_outliers ~ CTQ_KV + Alter_bei_erstem_NSSI + NSSI_Hfgkt_im_letzten_Monat + Sport, data = Complete_Results)
# summary(regr_CTQ_KV_bE) #CTQ: Körperliche Vernachlässigung
# 
# 
# 
# #### 4. pain perception & beta endorphin
# 
# ### have a look at the pain perception DURING the act of NSSI in correlation to the amount of beta endorphin
# Scatterplot_bE_Pain <- ggplot(bE_both_groups_no_outliers, aes(x = SCHMERZ_während_Durchführung_von_SV, y = beta_Endorphin)) %>%
#   + geom_point() # create a scatterplot of the beta endorphin / pain during NSSI pairs 
# Scatterplot_bE_Pain
# 
# cor_bE_Pain <- cor(bE_both_groups_no_outliers$SCHMERZ_während_Durchführung_von_SV, bE_both_groups_no_outliers$beta_Endorphin, use = "pairwise.complete.obs", method = "pearson") # calculate the correlation between bE and the pain perception during NSSI
# cor_bE_Pain
# 
# ### test the third hypothesis via regression analysis by predicting the amount of beta endorphin via the pain perception during the act of NSSI while controlling for frequency and amount of exercise
# regression_Pain_bE <- lm(bE_both_groups_no_outliers ~ SCHMERZ_während_Durchführung_von_SV + NSSI_Hfgkt_im_letzten_Monat + Sport, data = Complete_Results)
# summary(regression_Pain_bE)
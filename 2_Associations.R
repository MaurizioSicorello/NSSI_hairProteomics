library("here")
library("stringr")
library("caret")
library("psych")
library("MASS")
library("missRanger")


source(here("functions", "PLS_helperFuns.R"))

# ANALYSIS PLAN
# 1) permutation test of accuracy and model weights
# 2) test effect of data imputation

# QUESTIONS
# where are most of the clinical/medication data for healthy controls? which dimensions to correlate with the predicted scores?
# What is the best indicator for general psychopathology?
# focus on imputed or unimputed model? (e.g., for correlation of predicted scores with clinical variables. tendency towards complete model)





################################ 
# load & prepare clinical and hair proteomics data


################
# hair data

df_hair <- read.csv2(here("data", "Percentage_at_least_50.csv"))

protNames <- df_hair[, "Accession"]
protFuns <- df_hair[, c("Accession", "Biological_Process", "Cellular_Component", "Molecular_Function")]

abundInd <- which(names(df_hair) == "SV1_Abund."):which(names(df_hair) == "SV88_Abund.")
df_hair <- df_hair[, abundInd]

IDs_hair <- str_split(names(df_hair), "_", simplify = TRUE)[,1]

df_hair <- sapply(df_hair, as.numeric)

df_hair_trans <- data.frame(IDs_hair, t(df_hair))
rownames(df_hair_trans) <- c()
names(df_hair_trans) <- c("ID", protNames)


################
# clinical data

df_clin <- read.csv2(here("data", "Complete_Results.csv"))
names(df_clin)[3] <- "ID" 



################################ 
# Machine Learning models


################
# preparation


# box-cox transformation of proteomics
boxcoxTrans <- function(x){
  
  lambdaTrace <- boxcox(x ~ 1)
  lambda <- lambdaTrace$x[which.max(lambdaTrace$y)]
  
  if(lambda == 0){
    x_trans <- log(x)
  }else{
    x_trans <- (x^lambda - 1) / lambda
  }
  
  unname(x_trans)
  
}

df_hair_trans_noID <- df_hair_trans[,-1]
df_hair_trans[, 2:ncol(df_hair_trans)] <- apply(df_hair_trans_noID, 2, boxcoxTrans)



# check for low variance predictors
sum(describe(df_hair_trans)$sd < 1)


# scale predictors
df_hair_trans[, 2:ncol(df_hair_trans)] <- scale(df_hair_trans[, -1])



################
# PLS

# merge with outcome
df_PLSgroup <- merge(df_clin[, c("ID", "Group")], df_hair_trans, by = "ID") 

# subset complete data
complCaseInd <- describe(df_PLSgroup)$n == nrow(df_PLSgroup)
df_PLSgroup <- df_PLSgroup[, complCaseInd]


PLSmodel_group <- PLSnestedCV(outcome = df_PLSgroup[,2], 
                              predictors = df_PLSgroup[, -(1:2)],
                              nrepeats = 5,
                              nfolds = 5,
                              maxComps = 30,
                              setSeed = 100,
                              classification = TRUE)

PLSmodel_group



# test significance of PLS model
#[takes ~1h to run on my machine, therefore code is commented out]
t1 <- Sys.time()
PLSmodel_group_perm <- permutePLSnestedCV(outcome = df_PLSgroup[,2],
                                      predictors = df_PLSgroup[, -(1:2)],
                                      nrepeats = 1,
                                      nfolds = 5,
                                      nperms = 10,
                                      maxComps = 30,
                                      setSeed = 200,
                                      classification = TRUE)
t2 <- Sys.time()
if(!file.exists(here("results", "PLSperm_group.csv"))){
  write.csv(PLSmodel_group_perm, here("results", "PLSperm_group.csv"), row.names = FALSE)
}else{
  warning("file already exists in folder")
}

PLSperm_group <- read.csv(here("results", "PLSperm_group.csv"))
sum(PLSperm_group >= PLSmodel_group$Accuracy)/length(PLSperm_group$x)


# final PLS model with most frequent number of picked components (=3)
table(PLSmodel_group$outerFolds[,2])

PLSmodel_group_final = plsda(x = PLSmodel_group$dat[,-1], 
                             y = PLSmodel_group$dat[,1],
                             ncomp = 3)
summary(PLSmodel_group_final)

confusionMatrix(predict(PLSmodel_group_final, PLSmodel_group$dat[,-1]), PLSmodel_group$dat[,1])

# compute variable importance and sort
PLS_varImp <- varImp(PLSmodel_group_final)
df_varImp <- data.frame(varImp = PLS_varImp$Overall, Accession = names(df_PLSgroup)[-c(1:2)])
df_varImp <- merge(df_varImp, protFuns, by = "Accession")
df_varImp <- df_varImp[order(df_varImp$varImp, decreasing = TRUE),]
head(df_varImp)

write.csv(df_varImp, here("results", "PLS_varImpAndFunctions.csv"))


################
# correlated predicted scores (cross-validated with Leave-one-out-Cross-validation) to other clinical variables

predictValues <- numeric(nrow(PLSmodel_group$dat))

for(i in 1:nrow(PLSmodel_group$dat)){
  
  PLSmodel_group_LOOCV = plsda(x = PLSmodel_group$dat[-i,-1], 
                               y = PLSmodel_group$dat[-i,1],
                               ncomp = 3)
  predictValues[i] <- predict(PLSmodel_group_LOOCV, PLSmodel_group$dat[i,-1], type = "prob")[,2,1]
  
}

df_predict <- data.frame(ID = df_PLSgroup[, "ID"], predictValues)
confusionMatrix(factor(ifelse(predictValues >= 0.5, 1, 0)), factor(df_PLSgroup$Group))

df_predictClin <- merge(df_predict, df_clin, by = "ID")

# correlation for NSSI participants only
cor(df_predictClin[df_predictClin$Group == 1, c("predictValues", "CTQ_total_score", "Age", "Weigth", "Heigth", "Long_term_Medication", "Smoking", "Sport",
                                                "BSL_23_BPD_symptoms_2")], use = "pairwise.complete.obs")
                   



################
# check random forest as an alternative

fitControl <- trainControl(
  method = "repeatedcv",
  number = 5,
  repeats = 1)

rfGrid <-  expand.grid(mtry = seq(from = 100, to = 200, by = 10))

rfFit <- train(DV ~ ., data = PLSmodel_group$dat, 
                 method = "cforest", 
                 trControl = fitControl,
                 tuneGrid  = rfGrid)
rfFit




################################
# pairwise t-tests (without missings)

nTests <- ncol(df_PLSgroup)-2

pvalues_ttests <- numeric(nTests)

for(i in 3:nTests){
  
  pvalues_ttests[i] <- t.test(df_PLSgroup[df_PLSgroup$Group == 1, i], df_PLSgroup[df_PLSgroup$Group == 0, i])$p.value
  
}

pvalues_ttests_fdr <- p.adjust(pvalues_ttests, method = "fdr")
sum(pvalues_ttests_fdr <= 0.05)

signProtFuns <- protFuns[complCaseInd, ][pvalues_ttests_fdr <= 0.05, ]
test1 <- signProtFuns[,1]





################################
# calculate PLS as a function of % imputed missing values


# merge with DV
df_PLS_withMissings <- merge(df_clin[, c("ID", "Group")], df_hair_trans, by = "ID")


# retain different numbers of variables according to % missings values, impute and perform PLS

allowPercMiss <- seq(from = 0.50, to = 1, by = 0.05)
missingness_ind <- describe(df_PLS_withMissings[,-c(1:2)])$n/nrow(df_PLS_withMissings)

perMiss_out <- data.frame(
  matrix(nrow = length(allowPercMiss), ncol = 3, dimnames = list(NULL, c("percMiss", "nFeatures", "accuracy")))
)


for(i in 1:length(allowPercMiss)){
  
  df_PLS_withMissings_reduced <- df_PLS_withMissings[, missingness_ind >= allowPercMiss[i]]
  
  medImputation <- preProcess(df_PLS_withMissings_reduced, method = "medianImpute")
  df_PLS_withMissings_impute <- predict(medImputation, df_PLS_withMissings_reduced)
  
  PLSmodel_group_imp <- PLSnestedCV(outcome = df_PLS_withMissings_impute[,2], 
                                predictors = df_PLS_withMissings_impute[, -(1:2)],
                                nrepeats = 5,
                                nfolds = 5,
                                maxComps = 30,
                                setSeed = 100,
                                classification = TRUE)
  
  perMiss_out[i, ] <- c(allowPercMiss[i], ncol(df_PLS_withMissings_impute[, -(1:2)]), PLSmodel_group_imp$Accuracy)
  
  
}

ggplot(data = perMiss_out, aes(y = accuracy, x = nFeatures)) +
  geom_point() +
  ylim(0.5, 1) +
  xlab("Number of Predictors")

ggsave(here("figures", "imputationAccuracy.png"), device = "png")





################################
# pairwise t-tests (with missings)

nTests <- ncol(df_PLS_withMissings)-2

pvalues_ttests <- numeric(nTests)

for(i in 3:nTests){
  
  pvalues_ttests[i] <- t.test(df_PLS_withMissings[df_PLS_withMissings$Group == 1, i], df_PLS_withMissings[df_PLS_withMissings$Group == 0, i])$p.value
  
}

pvalues_ttests_fdr <- p.adjust(pvalues_ttests, method = "fdr")
sum(pvalues_ttests_fdr <= 0.05)

signProtFuns <- protFuns[pvalues_ttests_fdr <= 0.05, ]
signProtFuns[,1]


################################
# PLS with stress variables only
df_PLS_stress <- data.frame(df_PLS_withMissings[,1:2], df_PLS_withMissings[, -c(1:2)][, str_detect(protFuns$Biological_Process, "stress")])

stress_medImputation <- preProcess(df_PLS_stress, method = "medianImpute")
df_PLS_stress_imputed <- predict(stress_medImputation, df_PLS_stress)

PLSmodel_group_stress <- PLSnestedCV(outcome = df_PLS_stress_imputed[,2], 
                              predictors = df_PLS_stress_imputed[, -(1:2)],
                              nrepeats = 5,
                              nfolds = 5,
                              maxComps = 30,
                              setSeed = 400,
                              classification = TRUE)

PLSmodel_group_stress

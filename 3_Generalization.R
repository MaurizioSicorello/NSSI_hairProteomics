#### safer to start a new r session before running this script!

library("here")
library("stringr")
library("caret")
library("psych")
library("MASS")
library("missMethods")
library("lsr")
library("Hmisc")

source(here("functions", "PLS_helperFuns.R"))


################################ 
# TO DO

# export and use transformation paremeters from original data



################################ 
# load & prepare hair proteomics data


################
# hair data

df_hair <- read.csv2(here("data", "GT01-GT103_proteins.csv"))
df_primaryFunctions <- read.csv2(here("data", "Final_protein_functions.csv"))[, c("Accession", "Primary_function")]

# subset proteins existent in other dataset
df_hair <- df_hair[df_hair$Accession %in% df_primaryFunctions$Accession, ]
protNames <- df_hair[, "Accession"]

# subset df for abundance scores
abundInd <- which(names(df_hair) == "Abundance..F1..Sample..GT01.H3cm"):which(names(df_hair) == "Abundance..F97..Sample..GT103.H3cm")
df_hair <- df_hair[, abundInd]

# create participant IDs
IDs_hair <- str_split(names(df_hair), "GT", simplify = TRUE)[,2]
IDs_hair <- str_split(IDs_hair, ".H", simplify = TRUE)[,1]

df_hair <- sapply(df_hair, as.numeric)

df_hair_trans <- data.frame(IDs_hair, t(df_hair))
rownames(df_hair_trans) <- c()
names(df_hair_trans) <- c("ID", protNames)

# remove participant with double hair data
df_hair_trans <- df_hair_trans[-which(df_hair_trans$ID == "98B"), ]
df_hair_trans[df_hair_trans$ID == "98A", "ID"] <- 98
df_hair_trans$ID <- as.numeric(df_hair_trans$ID)



################
# clinical data

df_clin <- read.csv2(here("data", "Grouping_GerontoTel_4Maurizio.csv"))
names(df_clin) <- c("ID", "Group", "Coding")
df_clin <- df_clin[!is.na(df_clin$Group), ]


################
# load pls model

PLS_model <- readRDS(here("results", "PLSmodel_generalization.rds"))


################
# load preprocessing parameters of NSSI sample
preproPars <- read.csv(here("results", "preprocessingParameters.csv"))


################################ 
# Preprocessing

# remove variables with more than 50% missings
hist(describe(df_hair_trans)[, "n"])
df_hair_trans <- df_hair_trans[, (describe(df_hair_trans)[, "n"] >= nrow(df_hair_trans)*0.5)]


# box-cox transformation of proteomics
df_hair_trans_noID <- df_hair_trans[,-1]
df_hair_trans_prepro <- df_hair_trans_noID

for(i in 1:ncol(df_hair_trans_prepro)){
  
  if(preproPars$Lambda[i] == 0){
    df_hair_trans_prepro[,i] <- log(df_hair_trans_noID[,i])
  }else{
    df_hair_trans_prepro[,i] <- (df_hair_trans_noID[,i]^preproPars$Lambda[i] - 1) / preproPars$Lambda[i]
  }
  
  df_hair_trans_prepro[,i] <- (df_hair_trans_prepro[,i]-preproPars$Mean[i])/preproPars$SD[i]
}

describe(df_hair_trans_prepro)
mean(describe(df_hair_trans_prepro)$mean)
skew(describe(df_hair_trans_prepro)$mean)
sum(describe(df_hair_trans_prepro)$mean < 0)/ncol(df_hair_trans_prepro)

df_hair_trans <- cbind(df_hair_trans$ID, df_hair_trans_prepro)
names(df_hair_trans)[1] <- "ID"


# median imputation
medImputation <- preProcess(df_hair_trans, method = "medianImpute")
df_hair_trans_imp <- predict(medImputation, df_hair_trans)

# merge with group membership data (imputed)
dfpred_noMiss <- merge(df_clin[, c("ID", "Group")], df_hair_trans_imp, by = "ID")

# merge with group membership data (missings)
dfpred_withMiss <- merge(df_clin[, c("ID", "Group")], df_hair_trans, by = "ID")

# # test what happens when z standardized
# dfpred_noMiss[, 3:ncol(dfpred_noMiss)] <- scale(dfpred_noMiss[, -c(1:2)])
# describe(dfpred_noMiss)

protListPLS <- names(dfpred_noMiss)[-c(1,2)]
if(!file.exists(here("data", "protListPLS_GeneralizationSample.csv"))){
  write.csv(perMiss_out, here("data", "protListPLS_GeneralizationSample.csv"), row.names = FALSE)
}else{
  warning("file already exists in folder")
}


# # add missing variables to repeat PLS model
# missingPreds <- df_primaryFunctions$Accession[!(df_primaryFunctions$Accession %in% names(dfpred_noMiss)[-1])]
# 
# df_missingPreds <- data.frame(
#   matrix(
#     data = 0,
#     nrow = nrow(dfpred_noMiss),
#     ncol = length(missingPreds),
#     dimnames = list(NULL, missingPreds)
#     )
# )
# 
# dfpred_noMiss <- data.frame(dfpred_noMiss, df_missingPreds)


# # order predictors as in model training (important step! PLS model only goes by order, not by name)
# dfpred_noMiss <- dfpred_noMiss[, c("ID", "Group", df_primaryFunctions$Accession)]



################################
# generalization attempt of PLS results

predictCat <- predict(PLS_model, newdata = dfpred_noMiss[, -c(1,2)])
predictProb <- predict(PLS_model, newdata = dfpred_noMiss[, -c(1,2)], type = "prob")[,2,1]
predictScore <- predict(PLS_model, newdata = dfpred_noMiss[, -c(1,2)], type = "raw")[,2,1]

df_results <- data.frame(dfpred_noMiss[,c("ID", "Group")], predictCat, predictProb, predictScore)

describeBy(df_results, group = "Group")

summary(aov(df_results$predictScore ~ as.factor(df_results$Group)))
summary(aov(df_results$predictProb ~ as.factor(df_results$Group)))

etaSquared(aov(df_results$predictProb ~ as.factor(df_results$Group)))

sum(df_results$predictProb >= 0.5)/nrow(df_results)

ggplot(data = df_results, aes(x = predictProb, fill = factor(Group))) +
  
  geom_density(alpha = .4, aes(y=..scaled..)) + 
  geom_rug(aes(x = predictScore, y = 0, color =  factor(Group)), position = position_jitter(height = 0), size = 0.05) +
  
  theme_classic() +
  xlab("Class Probabilities") + ylab("density")

ggsave(here("figures", "generalizationDensityPlot.png"), device = "png")
  

# classification accuracies for different tasks
caret::confusionMatrix(factor(ifelse(dfpred_noMiss$Group == 0, 0, 1)), factor(predictCat)) # HC versus disorder
caret::confusionMatrix(factor(ifelse(dfpred_noMiss$Group == 2 , 1, 0)), factor(predictCat)) # HC & Depressed versus disorder

caret::confusionMatrix(factor(df_results[df_results$Group != 2, "Group"]), factor(df_results[df_results$Group != 2, "predictCat"])) # HC versus depressed

rec_HCvsSuic <- df_results[df_results$Group != 1, "Group"] # HC versus suicide
rec_HCvsSuic <- ifelse(rec_HCvsSuic == 2, 1, 0)
caret::confusionMatrix(factor(rec_HCvsSuic), factor(df_results[df_results$Group != 1, "predictCat"])) 

rec_DepvsSuic <- df_results[df_results$Group != 0, "Group"] # Dep versus suicide [check if correct]
rec_DepvsSuic <- ifelse(rec_HCvsSuic == 1, 0, 1)
caret::confusionMatrix(factor(rec_DepvsSuic), factor(df_results[df_results$Group != 0, "predictCat"]))

suicideVSrest <- factor(ifelse(dfpred_noMiss$Group == 0 | dfpred_noMiss$Group == 1, 0, 1))
pROC_obj <- pROC::roc(suicideVSrest,
                       predictProb,
                smoothed = TRUE,
                # arguments for ci
                ci=TRUE, ci.alpha=0.9, stratified=FALSE,
                # arguments for plot
                plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,
                print.auc=TRUE, show.thres=TRUE)


sens.ci <- pROC::ci.se(pROC_obj)
plot(sens.ci, type="shape", col="lightblue")
## Warning in plot.ci.se(sens.ci, type = "shape", col = "lightblue"): Low
## definition shape.
plot(sens.ci, type="bars")




###############################################
# predicting scores from demographic/clinical characteristics

df_clin2 <- read.csv2(here("data", "Hair_proteomics_cohort_database_4Maurizio.csv"))
df_clin2$ID <- str_sub(df_clin2$ID, 3, -1)

df_inspectProbs <- merge(df_results, df_clin2, by = "ID")

round(cor(df_inspectProbs[df_inspectProbs$Group.x == 2, c("predictProb", "BMI", "Sex", "Age")], use = "pairwise.complete"), 2)
cor.test(df_inspectProbs[df_inspectProbs$Group.x == 2, "predictProb"], df_inspectProbs[df_inspectProbs$Group.x == 2, "BMI"])
cor.test(df_inspectProbs[df_inspectProbs$Group.x == 2, "predictProb"], df_inspectProbs[df_inspectProbs$Group.x == 2, "Sex"])
cor.test(df_inspectProbs[df_inspectProbs$Group.x == 2, "predictProb"], df_inspectProbs[df_inspectProbs$Group.x == 2, "Age"])

df_inspectProbs$GroupNew <- ifelse(df_inspectProbs$Sex == 0, 3, df_inspectProbs$Group.x)

ggplot(data = df_inspectProbs, aes(x = predictProb, fill = factor(GroupNew))) +
  
  geom_density(alpha = .4, aes(y=..scaled..)) + 
  geom_rug(aes(x = predictScore, y = 0, color =  factor(GroupNew)), position = position_jitter(height = 0), size = 0.05) +
  
  theme_classic() +
  xlab("Class Probabilities") + ylab("density")

ggsave(here("figures", "generalizationDensityPlot_gender.png"), device = "png")


###############################################
# bivariate associations

nTests <- ncol(dfpred_withMiss)-2

pvalues_ttests <- numeric(nTests)
cohenD_ttests <- numeric(nTests)
naInd_ttests <- numeric(nTests)
cvAcc_ttests <- numeric(nTests)

# option A
# dfpred_withMiss$Group <- ifelse(dfpred_withMiss$Group == 2 , 1, 0)

# option B
dfpred_withMiss <- dfpred_withMiss[dfpred_withMiss$Group != 1, ]
dfpred_withMiss$Group <- ifelse(dfpred_withMiss$Group == 2, 1, 0)

repeats = 5
folds = 5

for(i in 3:ncol(dfpred_withMiss)){
  
  pvalues_ttests[i-2] <- t.test(dfpred_withMiss[dfpred_withMiss$Group == 1, i], dfpred_withMiss[dfpred_withMiss$Group == 0, i])$p.value
  cohenD_ttests[i-2] <- cohen.d(dfpred_withMiss[, c(2, i)], "Group")$cohen.d[2]
  naInd_ttests[i-2] <- sum(!is.na(dfpred_withMiss[,i]))
  
  cvAcc_ttests[i-2] <- train(
    form =  factor(Group) ~ .,
    data = dfpred_withMiss[,c(2,i)],
    trControl = trainControl(method = "repeatedcv", number = folds, repeats = repeats),
    method = "glm",
    family = "binomial",
    na.action = "na.omit")$results[["Accuracy"]]
  
}

pvalues_ttests_fdr <- p.adjust(pvalues_ttests, method = "fdr")
ttest_signInd_fdr <- ifelse(pvalues_ttests_fdr <= .05, "sign", "")

df_ttestOut <- data.frame("Accession" = preproPars$protList, "N" = naInd_ttests, "cvAcc" = cvAcc_ttests, "cohensD" = cohenD_ttests, "p_uncorr" = pvalues_ttests, 
                          "p_FDR" = pvalues_ttests_fdr, "signInd" = ttest_signInd_fdr)


if(!file.exists(here("results", "tableSX_bivariateAssociations_generalization_verB.csv"))){
  write.csv(df_ttestOut, here("results", "tableSX_bivariateAssociations_generalization_verB.csv"), row.names = FALSE)
}else{
  warning("file already exists in folder")
}


df_NSSI_bivar <- read.csv(here("results", "tableSX_bivariateAssociations.csv"))
df_DEP_bivar <- read.csv(here("results", "tableSX_bivariateAssociations_generalization.csv"))
df_DEP_bivar_B <- read.csv(here("results", "tableSX_bivariateAssociations_generalization_verB.csv"))

df_bivar_bothsample <- merge(df_NSSI_bivar, df_DEP_bivar, by = "Accession")
plot(df_bivar_bothsample$cohensD.x, df_bivar_bothsample$cohensD.y)
cor(df_bivar_bothsample$cohensD.x, df_bivar_bothsample$cohensD.y)

df_bivar_bothsample <- merge(df_NSSI_bivar, df_DEP_bivar_B, by = "Accession")
plot(df_bivar_bothsample$cohensD.x, df_bivar_bothsample$cohensD.y)
cor(df_bivar_bothsample$cohensD.x, df_bivar_bothsample$cohensD.y)

df_DEP_bivar_B[df_DEP_bivar_B$Accession == "E5RH77", ]

which(names(df_hair_trans) == "E5RH77") 

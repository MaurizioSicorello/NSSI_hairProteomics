library("here")
library("stringr")
library("caret")
library("psych")
library("MASS")
library("missRanger")
library("reshape")
library("missMethods")
library("wordcloud")
library("rgl")
library("dplyr")
library("ggExtra")
library("data.table")
library("pROC")
library("plotROC")
library("psychometric")
library("gameofthrones")

source(here("functions", "PLS_helperFuns.R"))


# ANALYSIS PLAN
# 1) permutation test of accuracy and model weights
# 2) test effect of data imputation

# QUESTIONS
# where are most of the clinical/medication data for healthy controls? which dimensions to correlate with the predicted scores?
# What is the best indicator for general psychopathology?
# focus on imputed or unimputed model? (e.g., for correlation of predicted scores with clinical variables. tendency towards complete model)

# NEXT STEPS
# auf zusätzliche HC-Daten warten
# Hauptmodell ist mit allen Variablen und Imputation

# Tors typische Analyse nachmachen: Wie viele randomly sampled proteine braucht man für eine gute vorhersage? (0-200 ausprobieren)

# paper title: Demonstrating the feasibility and applicability of untargeted hair proteomics for human psychopathology



################################ 
# load & prepare clinical and hair proteomics data


################
# hair data

df_hair <- read.csv2(here("data", "Percentage_at_least_50.csv"))
df_primaryFunctions <- read.csv2(here("data", "Final_protein_functions.csv"))[, c("Accession", "Primary_function")]

df_hair <- merge(df_hair, df_primaryFunctions, by = "Accession")

protNames <- df_hair[, "Accession"]
protFuns <- df_hair[, c("Accession", "Primary_function", "Description", "Biological_Process", "Cellular_Component", "Molecular_Function")]

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
# Preprocessing

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

# median imputation
medImputation <- preProcess(df_hair_trans, method = "medianImpute")
df_hair_trans_imp <- predict(medImputation, df_hair_trans)

# merge with group membership data
dfpred_withMiss <- merge(df_clin[, c("ID", "Group")], df_hair_trans, by = "ID")
dfpred_noMiss <- merge(df_clin[, c("ID", "Group")], df_hair_trans_imp, by = "ID")


# hyperparameters for cross-validation
repeats = 5
folds = 5




################################
# naive model based on two cluster solution

# prepare dataframe for model
clustAssign <- read.csv(here("results", "clusterMembership.csv"))

clust1Means <- rowMeans(dfpred_noMiss[,-c(1:2)][, clustAssign$clustAssign == 1])
clust2Means <- rowMeans(dfpred_noMiss[,-c(1:2)][, clustAssign$clustAssign == 2])

df_clustMeans <- data.frame("Group" = dfpred_noMiss$Group, clust1Means, clust2Means)

# compute model
set.seed(sqrt(666))

baselineModel = train(
  form =  factor(Group) ~ .,
  data = df_clustMeans,
  trControl = trainControl(method = "repeatedcv", number = folds, repeats = repeats),
  method = "glm",
  family = "binomial"
)

baselineModel
summary(baselineModel)

# # permutation test
# iterations <- 1000
# 
# accPerm_baseline <- numeric(iterations)
# 
# for(i in 1:iterations){
#   
#   df_clustMeans_perm <- data.frame("Group" = sample(df_clustMeans$Group), df_clustMeans[,-1])
#   
#   accPerm_baseline[i] <- train(
#     form =  Group ~ .,
#     data = df_clustMeans,
#     trControl = trainControl(method = "cv", number = folds),
#     method = "glm",
#     family = "binomial")$results$Accuracy
#   
# }
# 
# write.csv(accPerm_baseline, here("results", "clustPerm_group.csv"), row.names = FALSE)

accPerm_baseline <- read.csv(here("results", "clustPerm_group.csv"))

sum(accPerm_baseline >= baselineModel$results$Accuracy)/nrow(accPerm_baseline)



################
# PLS

PLSmodel_group <- PLSnestedCV(outcome = dfpred_noMiss[,2], 
                              predictors = dfpred_noMiss[, -(1:2)],
                              nrepeats = 5,
                              nfolds = 5,
                              maxComps = 30,
                              setSeed = 100,
                              classification = TRUE)

PLSmodel_group


# # test significance of PLS model
# # [takes ~1h to run on my machine, therefore code is commented out]
# PLSmodel_group_perm <- permutePLSnestedCV(outcome = dfpred_noMiss[,2],
#                                           predictors = dfpred_noMiss[, -(1:2)],
#                                           nrepeats = 1,
#                                           nfolds = 5,
#                                           nperms = 1000,
#                                           maxComps = 30,
#                                           setSeed = 200,
#                                           classification = TRUE)
# if(!file.exists(here("results", "PLSperm_group.csv"))){
#   write.csv(PLSmodel_group_perm, here("results", "PLSperm_group.csv"), row.names = FALSE)
# }else{
#   warning("file already exists in folder")
# }

PLSperm_group <- read.csv(here("results", "PLSperm_group.csv"))
sum(PLSperm_group >= PLSmodel_group$Accuracy)/length(PLSperm_group$x)



################
# final PLS model with most frequent number of picked components (=3)

hist(PLSmodel_group$outerFolds[,2])
table(PLSmodel_group$outerFolds[,2])
14/sum(table(PLSmodel_group$outerFolds[,2]))


PLSmodel_group_final = caret::plsda(x = dfpred_noMiss[,-c(1:2)], 
                             y = factor(dfpred_noMiss[,2]),
                             ncomp = 3)
summary(PLSmodel_group_final)

confusionMatrix(predict(PLSmodel_group_final, dfpred_noMiss[,-c(1:2)]), factor(dfpred_noMiss[, 2]))




################
# ROC curve
ctrl <- trainControl(method="cv",
                     number = 5,
                     summaryFunction=twoClassSummary, 
                     classProbs=T,
                     savePredictions = T)

df_ROC <- dfpred_noMiss[,-1]
df_ROC$Group <- as.factor(ifelse(df_ROC$Group == 1, "NSSI", "HC"))
rocFit <- train(Group ~ ., data=df_ROC, 
               method="pls", preProc=c("center", "scale"), 
               trControl=ctrl,
               tuneGrid = data.frame("ncomp" = 3)
               )

rocFit


# Plot:
png(here("figures", "PLS_ROC.png"))
plot.roc(rocFit$pred$obs, rocFit$pred$NSSI)
dev.off()
# # ggplot alternative
# ggplot(rocFit$pred, 
#        aes(m = NSSI, d = factor(obs, levels = c("HC", "NSSI")))) + 
#   geom_roc(hjust = -0.4, vjust = 1.5) + coord_equal()




################
# check random forest as an alternative

set.seed(20)

fitControl <- trainControl(
  method = "repeatedcv",
  number = 5,
  repeats = 5)

rfGrid <-  expand.grid(mtry = seq(from = 5, to = 200, by = 10))

rfFit <- train(DV ~ ., data = PLSmodel_group$dat, 
               method = "cforest", 
               trControl = fitControl,
               tuneGrid  = rfGrid)
rfFit



################################
# calculate PLS as a function of % imputed missing values


# merge with DV
df_PLS_withMissings <- merge(df_clin[, c("ID", "Group")], df_hair_trans, by = "ID")


# retain different numbers of variables according to % missings values, impute and perform PLS
allowPercMiss <- seq(from = 0.50, to = 1, by = 0.05)
missingness_ind <- describe(df_PLS_withMissings[,-c(1:2)])$n/nrow(df_PLS_withMissings)

perMiss_out <- data.frame(
  matrix(nrow = length(allowPercMiss), ncol = 3+repeats*folds*2, dimnames = list(NULL, c("percMiss", "nFeatures", "accuracy", 
                                                                                         paste0("acc_", c(1:(repeats*folds))),
                                                                                         paste0("nComp_", c(1:(repeats*folds))))))
)


for(i in 1:length(allowPercMiss)){
  
  df_PLS_withMissings_reduced <- df_PLS_withMissings[, c(TRUE, TRUE, missingness_ind >= allowPercMiss[i])]
  
  medImputation <- preProcess(df_PLS_withMissings_reduced, method = "medianImpute")
  df_PLS_withMissings_impute <- predict(medImputation, df_PLS_withMissings_reduced)
  
  PLSmodel_group_imp <- PLSnestedCV(outcome = df_PLS_withMissings_impute[,2], 
                                    predictors = df_PLS_withMissings_impute[, -(1:2)],
                                    nrepeats = repeats,
                                    nfolds = folds,
                                    maxComps = 30,
                                    setSeed = 100,
                                    classification = TRUE)
  
  perMiss_out[i, ] <- c(allowPercMiss[i], 
                        ncol(df_PLS_withMissings_impute[, -(1:2)]), 
                        PLSmodel_group_imp$Accuracy,
                        stack(PLSmodel_group_imp$outerFolds)$values)
  
}

if(!file.exists(here("results", "PLSbyMissings.csv"))){
  write.csv(perMiss_out, here("results", "PLSbyMissings.csv"), row.names = FALSE)
}else{
  warning("file already exists in folder")
}


dfMissPlot <- read.csv(here("results", "PLSbyMissings.csv"))
dfMissPlot$percMiss <- factor(dfMissPlot$percMiss)

dfMissAccPlot <- dfMissPlot[, 1:which(names(dfMissPlot) == "acc_25")]

dfMissAccPlot_long <- melt(data = dfMissAccPlot, id.vars = names(dfMissAccPlot)[1:3])


ggplot(data = dfMissAccPlot_long, aes(y = value, x = nFeatures)) +
  
  geom_point(alpha = 0.1, position = position_jitter(w = 5, h = 0)) +
  
  stat_summary(aes(group = 1), geom = "line", fun = mean, size = 1, colour = "white") +
  stat_summary(aes(group = 1), geom = "point", fun = mean, size = 2.8, colour = "white") +
  
  stat_summary(aes(group = 1), geom = "line", fun = mean, size = 0.7) +
  stat_summary(aes(group = 1), geom = "point", fun = mean, size = 2.5) +
  
  theme_classic() +
  
  ylim(0.5, 1) +
  xlab("Number of Predictors") +
  ylab("Accuracy")

ggsave(here("figures", "imputationAccuracy.png"), device = "png")
ggsave(here("figures", "imputationAccuracy.pdf"), device = "pdf")



################################
# Performance of stress-related proteins

######## PLS with stress variables only
df_PLS_stress <- data.frame(dfpred_noMiss[,1:2], dfpred_noMiss[, -c(1:2)][, str_detect(protFuns$Biological_Process, "stress response")])

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

######## PLS without stress variables
df_PLS_noStress <- data.frame(dfpred_noMiss[,1:2], dfpred_noMiss[, -c(1:2)][, !str_detect(protFuns$Biological_Process, "stress response")])

noStress_medImputation <- preProcess(df_PLS_noStress, method = "medianImpute")
df_PLS_noStress_imputed <- predict(noStress_medImputation, df_PLS_noStress)

PLSmodel_group_noStress <- PLSnestedCV(outcome = df_PLS_noStress_imputed[,2], 
                                       predictors = df_PLS_noStress_imputed[, -(1:2)],
                                       nrepeats = 5,
                                       nfolds = 5,
                                       maxComps = 30,
                                       setSeed = 400,
                                       classification = TRUE)
PLSmodel_group_noStress


######## hair cortisol simulation

# function to simulate data
simCort <- function(effSize, N_HC, N_NSSI){
  
  HC_sim <- rnorm(N_HC)
  NSSI_sim <- rnorm(N_NSSI, mean = effSize)
  
  df_simCort <- data.frame("Group" = c(rep("HC", N_HC), rep("NSSI", N_NSSI)),
                            "value" = c(HC_sim, NSSI_sim))
  
  return(df_simCort)
}

# simulations
iterations = 100
d1_out <- numeric(iterations)
d2_out <- numeric(iterations)


# d = 0.213
for(i in 1:iterations){
  
  set.seed(100+i)
  df_simCort <- simCort(effSize = 0.213, N_HC = 32, N_NSSI = 36)
  
  set.seed(200+i)
  d1_Model = train(
    form =  factor(Group) ~ .,
    data = df_simCort,
    trControl = trainControl(method = "repeatedcv", number = folds, repeats = repeats),
    method = "glm",
    family = "binomial"
  )
  
  d1_out[i] <- d1_Model$results$Accuracy
  
}

mean(d1_out)
sd(d1_out)


# d = 0.478
for(i in 1:iterations){
  
  set.seed(100+i)
  df_simCort <- simCort(effSize = 0.478, N_HC = 32, N_NSSI = 36)
  
  set.seed(200+i)
  d2_Model = train(
    form =  factor(Group) ~ .,
    data = df_simCort,
    trControl = trainControl(method = "repeatedcv", number = folds, repeats = repeats),
    method = "glm",
    family = "binomial"
  )
  
  d2_out[i] <- d2_Model$results$Accuracy
  
}

mean(d2_out)
sd(d2_out)







################
# number of sufficient and necessary features

n_features <- seq(5, 605, by = 5)

######## sufficient number of features
dfIncl_out <- data.frame(
  matrix(nrow = length(n_features),
  ncol = 25,
  dimnames = list(NULL, c(paste0("acc_", c(1:25))))
  )
)

for(i in 1:length(n_features)){
  
  dfIncl_out[i, ] <- PLS_selectPred_rep(data = dfpred_noMiss, n_select = n_features[i], selection = "include")
  cat(paste("\r Iteration", i, "of", length(n_features)))
  
}

dfInclOut <- data.frame(n_features, 
                             "strategy" = rep("inclusion", times = length(n_features)),
                             "meanAcc" = rowMeans(dfIncl_out),
                             dfIncl_out)

if(!file.exists(here("results", "PLSbyInclusion.csv"))){
  write.csv(dfInclOut, here("results", "PLSbyInclusion.csv"), row.names = FALSE)
}else{
  warning("file already exists in folder")
}


######## necessary number of features
dfExcl_out <- data.frame(
  matrix(nrow = length(n_features),
         ncol = 25,
         dimnames = list(NULL, c(paste0("acc_", c(1:25))))
  )
)

for(i in 1:length(n_features)){
  
  dfExcl_out[i, ] <- PLS_selectPred_rep(data = dfpred_noMiss, n_select = n_features[i], selection = "exclude")
  cat(paste("\r Iteration", i, "of", length(n_features)))
  
}

dfExcl_out <- data.frame(n_features, 
                        "strategy" = rep("exclusion", times = length(n_features)),
                        "meanAcc" = rowMeans(dfExcl_out),
                        dfExcl_out)

if(!file.exists(here("results", "PLSbyExclusion.csv"))){
  write.csv(dfExcl_out, here("results", "PLSbyExclusion.csv"), row.names = FALSE)
}else{
  warning("file already exists in folder")
}


######## plot number of necessary and sufficient features
dfInclOut <- read.csv(here("results", "PLSbyInclusion.csv"))
dfExcl_out <- read.csv(here("results", "PLSbyExclusion.csv"))

dfPlotSelection <- rbind(dfInclOut, dfExcl_out)

ggplot(data = dfPlotSelection, aes(y = meanAcc, x = n_features, color = strategy)) +
  
  geom_point() + 
  geom_smooth() +
  
  theme_classic() +
  theme(legend.position = c(0.2, 0.2)) +
  
  scale_x_continuous("Number of best predictors included/excluded", breaks = seq(0, 600, 50)) +
  scale_y_continuous("Accuracy", breaks = seq(0.5, 1, 0.05), limits = c(0.5, 1)) +
  scale_color_discrete(name = "Selection Strategy")

ggsave(here("figures", "selection.png"))



################
# analyze component loadings

library(mixOmics)
PLSmodel_group_final = mixOmics::plsda(X = dfpred_noMiss[,-c(1:2)], 
                                    Y = factor(dfpred_noMiss[,2]),
                                    ncomp = 3)
mixOmics::plotVar(PLSmodel_group_final, comp = c(1,2))
mixOmics::plotVar(PLSmodel_group_final, comp = c(1,3))
mixOmics::plotVar(PLSmodel_group_final, comp = c(2,3))
detach("package:mixOmics", unload = TRUE)




################
# variable importance analyses

# compute variable importance and sort
PLS_varImp <- varImp(PLSmodel_group_final)

cor(PLS_varImp$Overall, abs(df_ttestOut$cohensD))
cor(PLS_varImp$Overall, abs(df_ttestOut$cohensD), method = "spearman")

df_varImp <- data.frame(varImp = PLS_varImp$Overall, Accession = names(dfpred_noMiss)[-c(1:2)])
df_varImp <- merge(df_varImp, protFuns, by = "Accession")
df_varImp <- df_varImp[order(df_varImp$varImp, decreasing = TRUE),]
head(df_varImp)

write.csv(df_varImp, here("results", "PLS_varImpAndFunctions.csv"), row.names = FALSE)


# load variable importance
df_varImp <- read.csv(here("results", "PLS_varImpAndFunctions.csv"))
table(df_varImp$Primary_function)
(78+53)/611

# plot sorted variable importance
ggplot() +
  geom_point(aes(y = sort(df_varImp$varImp), x = c(1:nrow(df_varImp)))) +
  
  ylab("variable importance") + xlab(NULL) +
  
  theme_classic()
ggsave(here("figures", "variableImportance.png"), device = "png")



# molecular functions of top 10 proteins
df_varImp[1:10, "Primary_function"]
df_varImp[1:10, "Biological_Process"]


# mean variable importance by molecular function
df_primFun <- df_varImp[, c("varImp", "Primary_function")]
primFun_unique <- unique(df_primFun$Primary_function)

wordWeight_out <- data.frame(
  matrix(nrow = length(primFun_unique), ncol = 2, dimnames = list(NULL, c("Primary_function", "weight")))
)
wordWeight_out$Primary_function <- primFun_unique 

for(i in 1:nrow(wordWeight_out)){
  
  primFun_weights <- df_primFun[df_varImp$Primary_function == wordWeight_out$Primary_function[i], 1]
  wordWeight_out[i,2] <- sum(primFun_weights)/length(primFun_weights)
  
}





################
# correlated predicted scores (cross-validated with Leave-one-out-Cross-validation) to other clinical variables
# https://www.biorxiv.org/content/10.1101/2020.08.17.255034v1.abstract
predictValues <- numeric(nrow(dfpred_noMiss))

for(i in 1:nrow(dfpred_noMiss)){
  
  PLSmodel_group_LOOCV = plsda(x = dfpred_noMiss[-i,-c(1:2)], 
                               y = factor(dfpred_noMiss[-i,2]),
                               ncomp = 3)
  predictValues[i] <- predict(PLSmodel_group_LOOCV, dfpred_noMiss[i,-c(1:2)], type = "prob")[,2,1]
  
}

df_predict <- data.frame(ID = dfpred_noMiss[, "ID"], predictValues)
confusionMatrix(factor(ifelse(predictValues >= 0.5, 1, 0)), factor(dfpred_noMiss$Group))

df_predictClin <- merge(df_predict, df_clin, by = "ID")



# extract PLS scores from main model
PLSscores <- data.frame("ID" = dfpred_noMiss$ID, 
           "Comp1" = pls::scores(PLSmodel_group_final)[,1], 
           "Comp2" = pls::scores(PLSmodel_group_final)[,2], 
           "Comp3" = pls::scores(PLSmodel_group_final)[,3])

df_predictClin <- merge(df_predictClin, PLSscores, by = "ID")

df_predictClin$BMI <- as.numeric(df_predictClin$BMI)

str(df_predictClin)

# correlate variables with full data BMI Age CTQ
plotCorrs <- rbind(cor(df_predictClin[, c("predictValues", "CTQ_total_score", "BMI", "Age", "BSL_23_BPD_symptoms_2", "BSL_23_dysfunc_behav_2", "BSL_23_condition")], use = "pairwise.complete.obs")[1,2:7],
                    cor(df_predictClin[, c("Comp1", "CTQ_total_score", "BMI", "Age", "BSL_23_BPD_symptoms_2", "BSL_23_dysfunc_behav_2", "BSL_23_condition")], use = "pairwise.complete.obs")[1,2:7],
                    cor(df_predictClin[, c("Comp2", "CTQ_total_score", "BMI", "Age", "BSL_23_BPD_symptoms_2", "BSL_23_dysfunc_behav_2", "BSL_23_condition")], use = "pairwise.complete.obs")[1,2:7],
                    cor(df_predictClin[, c("Comp3", "CTQ_total_score", "BMI", "Age", "BSL_23_BPD_symptoms_2", "BSL_23_dysfunc_behav_2", "BSL_23_condition")], use = "pairwise.complete.obs")[1,2:7])

# plot results
plotCorrs <- as.data.frame(plotCorrs)
names(plotCorrs)[c(1, 4, 5, 6)] <- c("CTQ", "BSL: Symptoms", "BSL: Behavior", "BSL: Personal State")
plotCorrs$`BSL: Personal State` <- abs(plotCorrs$`BSL: Personal State`)
plotCorrs$PLSscore <- factor(c("Overall", "Component 1", "Component 2", "Component 3"), levels = c("Overall", "Component 1", "Component 2", "Component 3"))
plotCorrs_long <- melt(plotCorrs)

plot_corrCI <- CIr(0, n = sum(!is.na(df_predictClin$Age)))

ggplot(data = plotCorrs_long, aes(x = variable, y = value, fill = PLSscore)) +
  geom_bar(stat="identity", position = "dodge", colour="black") +
  
  theme_classic() +
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = plot_corrCI[1], linetype = "dashed") +
  geom_hline(yintercept = plot_corrCI[2], linetype = "dashed") +
  
  scale_y_continuous(name = "Correlation", limits = c(-0.75, 0.75), breaks = seq(-0.8, 0.8, 0.2)) +
  xlab(NULL) +
  theme(legend.title=element_blank()) +
  theme(axis.text.x=element_text(angle = 45, vjust = 0.5)) +
  scale_fill_got_d(option = "Daenerys")

ggsave(here("figures", "dimensionalCorrelations.png"), device = "png")
ggsave(here("figures", "dimensionalCorrelations.pdf"), device = "pdf")


# correlation for NSSI participants only
cor(df_predictClin[df_predictClin$Group == 1, c("predictValues", "Comp1", "Comp2", "Comp3", "CTQ_total_score", "Age", "Weigth", "Heigth", "Long_term_Medication", "Smoking", "Sport",
                                                "BSL_23_BPD_symptoms_2")], use = "pairwise.complete.obs")


# age difference between groups
t.test(df_predictClin[df_predictClin$Group == 0, "Age"], df_predictClin[df_predictClin$Group == 1, "Age"])
ggplot(data = df_predictClin, aes(x = Age, fill = factor(Group))) +
  geom_histogram()

# CTQ difference between groups
t.test(df_predictClin[df_predictClin$Group == 0, "CTQ_total_score"], df_predictClin[df_predictClin$Group == 1, "CTQ_total_score"])
ggplot(data = df_predictClin, aes(x = CTQ_total_score, fill = factor(Group))) +
  geom_histogram()


################################
# age matching

# drop case with missing value
df_ageMatch <- df_predictClin[!is.na(df_predictClin$Age), ]
# drop NSSI patients with lowest ages until equal sample size is reached
for(i in 1:5){
  
  NSSIage <- df_ageMatch[df_ageMatch$Group == 1, c("ID", "Age")]
  NSSIage_order <- NSSIage[order(NSSIage$Age), ]
  NSSI_ID_remove <- NSSIage_order[1:i, "ID"]
  df_ageMatch_small <- df_ageMatch[which(!(df_ageMatch$ID %in% NSSI_ID_remove)), ]
  
  tOut <- t.test(df_ageMatch[df_ageMatch_small$Group == 0, "Age"], df_ageMatch[df_ageMatch_small$Group == 1, "Age"])
  print(tOut)
}

# remove pairs of extreme values in the two groups
HCage <- df_ageMatch[df_ageMatch$Group == 0, c("ID", "Age")]
HCage[order(HCage$Age), ]
ID_remove <- c(NSSI_ID_remove, "SV68", "SV57", "SV75", "SV74", "SV22", "SV32", "SV33", "SV38")
df_ageMatch_small <- df_ageMatch[which(!(df_ageMatch$ID %in% ID_remove)), ]
sum(!(df_ageMatch$ID %in% ID_remove))
t.test(df_ageMatch[df_ageMatch_small$Group == 0, "Age"], df_ageMatch[df_ageMatch_small$Group == 1, "Age"])
25.25-26.629

# rerun main model
df_PLS_ageMatch <- dfpred_noMiss[which(!(dfpred_noMiss$ID %in% ID_remove)), ]

set.seed(777)
PLSmodel_group <- PLSnestedCV(outcome = df_PLS_ageMatch[,2], 
                              predictors = df_PLS_ageMatch[, -(1:2)],
                              nrepeats = 5,
                              nfolds = 5,
                              maxComps = 30,
                              setSeed = 100,
                              classification = TRUE)
PLSmodel_group$Accuracy

# correlate predictions with age
predictValues_matched <- numeric(nrow(df_PLS_ageMatch))

for(i in 1:nrow(df_PLS_ageMatch)){
  
  PLSmodel_group_LOOCV = plsda(x = df_PLS_ageMatch[-i,-c(1:2)], 
                               y = factor(df_PLS_ageMatch[-i,2]),
                               ncomp = 3)
  predictValues_matched[i] <- predict(PLSmodel_group_LOOCV, df_PLS_ageMatch[i,-c(1:2)], type = "prob")[,2,1]
  
}

df_predict_matched <- data.frame(ID = df_PLS_ageMatch[, "ID"], predictValues_matched)
confusionMatrix(factor(ifelse(predictValues >= 0.5, 1, 0)), factor(dfpred_noMiss$Group))

df_predictClin_matched <- merge(df_predict_matched, df_clin, by = "ID")
cor.test(df_predictClin_matched$Age, df_predictClin_matched$predictValues_matched, use = "pairwise.complete.obs")
cor.test(df_predictClin_matched$CTQ_total_score, df_predictClin_matched$predictValues_matched, use = "pairwise.complete.obs")


################################
# predict age

dfpred_noMiss_age <- merge(df_clin[, c("ID", "Age")], dfpred_noMiss[, ])
dfpred_noMiss_age <- dfpred_noMiss_age[!is.na(dfpred_noMiss_age$Age), ]
dfpred_noMiss_age <- dfpred_noMiss_age[which(!(dfpred_noMiss_age$ID %in% ID_remove)), ]
dfpred_noMiss_age_HC <- dfpred_noMiss_age[dfpred_noMiss_age$Group == 0, -c(1,3)]
dfpred_noMiss_age_NSSI <- dfpred_noMiss_age[dfpred_noMiss_age$Group == 1, -c(1,3)]

set.seed(500)
train(Age ~ ., data=dfpred_noMiss_age_HC, 
      method="pls", preProc=c("center", "scale"), 
      trControl=trainControl(method="repeatedcv",
                             repeats = 5,
                             number = 5),
      tuneGrid = data.frame("ncomp" = 3)
)


PLSmodel_age_HC = plsr(data = dfpred_noMiss_age_HC,
                       Age ~ .,
                       ncomp = 3)

# People in the NSSI group are judged as being too old based on the healthy control model (similar to accelerated aging)
predAge_error <- predict(PLSmodel_age_HC, newdata = dfpred_noMiss_age_NSSI[,-1], type = "response", ncomp = 3) - dfpred_noMiss_age_NSSI$Age
t.test(predAge_error)
sd(predict(PLSmodel_age_HC, newdata = dfpred_noMiss_age_NSSI[,-1], type = "response", ncomp = 3) - dfpred_noMiss_age_NSSI$Age)


# repeat with age-matched sample
dfpred_noMiss_age <- merge(df_clin[, c("ID", "Age")], dfpred_noMiss[, ])
dfpred_noMiss_age <- dfpred_noMiss_age[!is.na(dfpred_noMiss_age$Age), ]
dfpred_noMiss_age_HC <- dfpred_noMiss_age[dfpred_noMiss_age$Group == 0, -c(1,3)]
dfpred_noMiss_age_NSSI <- dfpred_noMiss_age[dfpred_noMiss_age$Group == 1, -c(1,3)]

set.seed(1500)
train(Age ~ ., data=dfpred_noMiss_age_HC, 
      method="pls", preProc=c("center", "scale"), 
      trControl=trainControl(method="repeatedcv",
                             repeats = 5,
                             number = 5),
      tuneGrid = data.frame("ncomp" = 3)
)


PLSmodel_age_HC = plsr(data = dfpred_noMiss_age_HC,
                       Age ~ .,
                       ncomp = 3)

# People in the NSSI group are judged as being too old based on the healthy control model (similar to accelerated aging)
predAge_error <- predict(PLSmodel_age_HC, newdata = dfpred_noMiss_age_NSSI[,-1], type = "response", ncomp = 3) - dfpred_noMiss_age_NSSI$Age
t.test(predAge_error)
sd(predict(PLSmodel_age_HC, newdata = dfpred_noMiss_age_NSSI[,-1], type = "response", ncomp = 3) - dfpred_noMiss_age_NSSI$Age)



################################
# t-tests with FDR correction

nTests <- ncol(dfpred_withMiss)-2

pvalues_ttests <- numeric(nTests)
cohenD_ttests <- numeric(nTests)
naInd_ttests <- numeric(nTests)
cvAcc_ttests <- numeric(nTests)

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

df_ttestOut <- data.frame("Accession" = protNames, "N" = naInd_ttests, "cvAcc" = cvAcc_ttests, "cohensD" = cohenD_ttests, "p_uncorr" = pvalues_ttests, 
                          "p_FDR" = pvalues_ttests_fdr, "signInd" = ttest_signInd_fdr,
                          "Primary_function" = factor(protFuns$Primary_function))

if(!file.exists(here("results", "tableSX_bivariateAssociations.csv"))){
  write.csv(df_ttestOut, here("results", "tableSX_bivariateAssociations.csv"), row.names = FALSE)
}else{
  warning("file already exists in folder")
}

df_ttestOut <- read.csv(here("results", "tableSX_bivariateAssociations.csv"))

signVars <- df_ttestOut[df_ttestOut$signInd == "sign", ]
nrow(signVars)
length(unique(signVars$Primary_function))

nrow(signVars[signVars$cvAcc < .70, ])/nrow(signVars)

signVars[signVars$cvAcc > .79, ]

signVars[signVars$Primary_function == "translation", ]

df_ttestOut$signInd <- ifelse(df_ttestOut$signInd  == "sign", "p(FDR) < .05", "N.S.")

ggplot(df_ttestOut, aes(x = cohenD_ttests, y = cvAcc, color = signInd)) +
  
  geom_point() +
  xlim(c(-2.4, 2.4)) +
  
  xlab(expression(paste("Cohen's ", italic("d")))) +  ylab("Accuracy") + 
  
  theme_classic() +
  theme(legend.title=element_blank(), legend.position = c(0.8, 0.15)) 

ggsave(here("figures", "Figure1_singleEffects.png"), device = "png")
ggsave(here("figures", "Figure1_singleEffects.pdf"), device = "pdf")




## failed manhatten plot attempt (basically worked, but currently doesn't look too helpful)

# df_ttestOut <- df_ttestOut[order(df_ttestOut$Primary_function),]
# 
# df_ttestOut <- merge(df_ttestOut, data.frame("Primary_function" = unique(df_ttestOut$Primary_function), 
#            colorInd = rep(c(1:4), length.out = length(unique(df_ttestOut$Primary_function)))))
# 
# 
# ggplot(df_ttestOut, aes(x = c(1:nrow(df_ttestOut)), y = -log10(pvalues_ttests), color = colorInd)) +
#   
#   geom_point() +
#   
#   #scale_color_manual(values = rep(c("#276FBF", "#183059"), unique(length(df_ttestOut$Primary_function)))) +
#   
#   theme_minimal() + 
#   theme(legend.position = "none")





################################


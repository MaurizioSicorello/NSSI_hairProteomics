# paper concept:

# which metabolites make sense to look at?
# how are they distributed?
# are there clusters?
# do they help predict adversity/psychopathology


# questions:

# what percent of missings is tolerable to retain (and impute) variables?
# how about hard to show metabolites? Can we set them to zero and see how important they are?



library("here")
library("psych")
library("corrplot")
library("MASS")
library("cluster")
library("factoextra")
library("dbscan")

df <- read.csv2(here("data", "Percentage_at_least_50.csv"))

#convert abundance scores to numeric
abundInd <- which(names(df) == "SV1_Abund."):which(names(df) == "SV88_Abund.")
df[, abundInd] <- sapply(df[, abundInd], as.numeric)

str(df)
head(df)
describe(df)

# transpose abundance scores
df_trans <- t(df[, abundInd])

# describe data
dfDescribe <- describe(df_trans)

hist(dfDescribe$n/nrow(df_trans))

hist(dfDescribe$mean)
hist(log(dfDescribe$mean))

hist(dfDescribe$sd)
hist(log(dfDescribe$sd))

plot(dfDescribe$sd, dfDescribe$mean)

hist(dfDescribe$skew)

boxplot(df_trans)



numOut <- function(x){
  
  low.quart <- quantile(x, na.rm = T)[2]
  high.quart <- quantile(x, na.rm = T)[4]
  iqr <- high.quart-low.quart
  
  upper.thresh <- high.quart + 3*iqr
  lower.thresh <- low.quart - 3*iqr
  
  num_outliers <- sum((x < lower.thresh | x > upper.thresh), na.rm = T)
  num_outliers
}

hist(apply(df_trans, 2, numOut))


# smaller dataframe of metabolites with few missings

sum(dfDescribe$n/nrow(df_trans) == 1)
sum(dfDescribe$n/nrow(df_trans) > 0.95)
sum(dfDescribe$n/nrow(df_trans) > 0.90)

fullDatInd <- dfDescribe$n/nrow(df_trans) == 1

df_trans_s <- df_trans[, fullDatInd]

sum(complete.cases(df_trans_s))

corrM_noMiss <- cor(df_trans_s, method = "spearman")
corrplot(corrM_noMiss)

mean(corrM_noMiss)
sd(corrM_noMiss)
hist(corrM_noMiss)



# boxcox tranformed data

boxcoxTrans <- function(x){
  
  lambdaTrace <- boxcox(lm(x ~ 1))
  lambda <- lambdaTrace$x[which.max(lambdaTrace$y)]
  
  if(lambda == 0){
    x_trans <- log(x)
  }else{
    x_trans <- (x^lambda - 1) / lambda
  }
  
  unname(x_trans)
  
}


df_trans_s_box <- apply(df_trans_s, 2, boxcoxTrans)

corrM_noMiss_box <- cor(df_trans_s_box)
mean(corrM_noMiss_box)




# z-standardization

df_trans_s_box_z <- apply(df_trans_s_box, 2, scale)
describe(df_trans_s_box_z)


# retranspose for clustering

df_s_box_z <- t(df_trans_s_box_z)


# hierarchical agglomerative clustering

dd <- dist(df_s_box_z, method = "euclidean")
hc <- hclust(dd, method = "ward.D2")

fviz_dend(hc)
fviz_nbclust(df_s_box_z, FUN = hcut, method = "wss")
fviz_nbclust(df_s_box_z, FUN = hcut, method = "silhouette")


# dbscan

kNNdistplot(df_s_box_z, k = 3)
kNNdistplot(df_s_box_z, k = 5)
kNNdistplot(df_s_box_z, k = 8)


res.fpc3 <- dbscan(df_s_box_z, eps = 5.5, minPts = 3)
print(res.fpc3)

res.fpc5 <- dbscan(df_s_box_z, eps = 6, minPts = 5)
print(res.fpc5)

res.fpc8 <- dbscan(df_s_box_z, eps = 6, minPts = 8)
print(res.fpc8)


clustAssign <- res.fpc8$cluster

mean(cor(df_trans_s[, clustAssign == 1]))
sd(cor(df_trans_s[, clustAssign == 1]))

mean(cor(df_trans_s[, clustAssign == 0]))

corrplot(corrM_noMiss_box, order = "hclust", hclust.method = "ward.D2", addrect = 2)

molFunctionClust1 <- df[fullDatInd, ][clustAssign, ][, "Molecular_Function"]
molFunctionClust1

df[, "Molecular_Function"]


# test with PCA how much variance components capture

fa.results <- principal(corrM_noMiss_box)

varExpl <- fa.results$values/201

plot(varExpl)

cumVarExpl <- numeric(length(fa.results$values))
cumVarExpl[1] <- varExpl[1]

for(i in 2:length(cumVarExpl)){
  
  cumVarExpl[i] <- cumVarExpl[i-1] + varExpl[i]
  
}

plot(cumVarExpl)
length(cumVarExpl[cumVarExpl < 0.90])


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



# boxcox tranformed data

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


df_trans_box <- apply(df_trans, 2, boxcoxTrans)

corrM_box <- cor(df_trans_box, use = "pairwise.complete.obs")
corDistr <- corrM_box[upper.tri(corrM_box)]
mean(corDistr)
sd(corDistr)

ggplot() +
  geom_histogram(aes(x = corDistr), color = "black", fill = "grey") +
  
  ylab("Frequency") + xlab("Bivariate Correlations") +
  
  theme_classic()

ggsave(here("figures", "FigS1_correlHist.png"), device = "png")


# z-standardization

df_trans_box_z <- apply(df_trans_box, 2, scale)
describe(df_trans_box_z)


# retranspose for clustering

df_box_z <- t(df_trans_box_z)


# hierarchical agglomerative clustering

dd <- dist(df_box_z, method = "euclidean")
hc <- hclust(dd, method = "ward.D2")

fviz_dend(hc)
fviz_nbclust(df_box_z, FUN = hcut, method = "wss")
fviz_nbclust(df_box_z, FUN = hcut, method = "silhouette")

pdf(here("Figures", "Dendrogram.pdf"))
fviz_dend(
  hc,
  k = 2,
  horiz = TRUE,
  rect = TRUE,
  rect_fill = TRUE,
  rect_border = "jco",
  k_colors = "jco",
  cex = 0.3
)
dev.off()

clustAssign_hc <- cutree(hc, k = 2)

write.csv(data.frame("Accession" = df$Accession, "clustAssign" = clustAssign_hc), here("results", "clusterMembership.csv"), row.names = FALSE)

corM_Clust1 <- cor(df_trans_box_z[, clustAssign_hc == 1], use = "pairwise.complete")
mean(corM_Clust1[upper.tri(corM_Clust1)])
sd(corM_Clust1[upper.tri(corM_Clust1)])

corM_Clust2 <- cor(df_trans_box_z[, clustAssign_hc == 2], use = "pairwise.complete")
mean(corM_Clust2[upper.tri(corM_Clust2)])
sd(corM_Clust2[upper.tri(corM_Clust2)])


# dbscan

df_box_z_imp <- impute_median(df_box_z)

kNNdistplot(df_box_z_imp, k = 3)
kNNdistplot(df_box_z_imp, k = 5)
kNNdistplot(df_box_z_imp, k = 8)


res.fpc3 <- dbscan(df_box_z_imp, eps = 7.8, minPts = 3)
print(res.fpc3)

res.fpc5 <- dbscan(df_box_z_imp, eps = 7.8, minPts = 5)
print(res.fpc5)

res.fpc8 <- dbscan(df_box_z_imp, eps = 7.8, minPts = 8)
print(res.fpc8)


clustAssign <- res.fpc5$cluster

mean(cor(df_trans_box[, clustAssign == 1], use = "pairwise.complete")[upper.tri(cor(df_trans_box[, clustAssign == 1]))])
sd(cor(df_trans_box[, clustAssign == 1], use = "pairwise.complete")[upper.tri(cor(df_trans_box[, clustAssign == 1]))])

mean(cor(df_trans_s[, clustAssign == 0]))

corrplot(corrM_box, order = "hclust", hclust.method = "ward.D2", addrect = 2)

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

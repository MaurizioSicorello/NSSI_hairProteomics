
# nested k-fold CV evaluation of PLS
PLSnestedCV <- function(outcome, predictors, nrepeats, nfolds, maxComps = 30, setSeed = 1000, classification = F){
  
  fitControl <- trainControl(method = "repeatedcv",   
                             number = nfolds, 
                             repeats = nrepeats)
  
  # prepare input df
  df_pred = data.frame(outcome, predictors)
  df_pred = scale(df_pred)
  df_pred = as.data.frame(df_pred[!is.na(df_pred[,1]), ])
  names(df_pred)[1] = "DV"
  if(classification == T){df_pred$DV = as.factor(df_pred$DV)}
  
  # prepare output df
  df_out <- as.data.frame(matrix(nrow = nrepeats*nfolds, ncol = 2))
  df_foldwise <- as.data.frame(matrix(nrow = maxComps, ncol = nrepeats*nfolds))
  plsGrid <- expand.grid(ncomp = seq(1, maxComps))
  count = 0
  
  # main loop
  for(i in 1:nrepeats){
    
    set.seed(setSeed+count)
    plsFolds <- createFolds(df_pred[,1], k = nfolds)
    
    for(j in 1:nfolds){
      
      count = count+1
      
      trainSet <- df_pred[-plsFolds[[j]],]
      testSet <- df_pred[plsFolds[[j]],]
      
      f <- as.formula("DV ~ .")
      
      set.seed(setSeed+count)
      
      plsOptComp <- train(f,
                          data = trainSet,
                          trControl = fitControl,
                          method = "pls",
                          preProcess = c('scale', 'center'),
                          na.action = na.omit,
                          tuneGrid = plsGrid)
      
      if(classification == F){
        Accuracy <- 1 - sum((testSet[,1] - predict(plsOptComp, newdata = testSet))^2)/(var(testSet[,1])*(nrow(testSet)-1))
        df_out[count, 1] <- max(Accuracy, 0)
        df_out[count, 2] <- plsOptComp$bestTune$ncomp
        df_foldwise[, count] <- plsOptComp$results$RMSE
      }else{
        Accuracy <- sum(ifelse(testSet[,1] == predict(plsOptComp, newdata = testSet), 1, 0))/nrow(testSet)
        df_out[count, 1] <- Accuracy
        df_out[count, 2] <- plsOptComp$bestTune$ncomp
        df_foldwise[, count] <- plsOptComp$results$Accuracy
      }
      
    }
  }
  return(list(Accuracy = mean(df_out[,1]), outerFolds = df_out, innerFolds = df_foldwise, dat = df_pred))
}



# permutation test of PLS
permutePLSnestedCV <- function(outcome, predictors, nrepeats, nfolds, nperms, maxComps = 30, setSeed = 1000, classification = F){
  
  savePerm <- numeric(nperms)
  
  for(i in 1:nperms){
    
    outcomePerm <- sample(outcome)
    savePerm[i] <- PLSnestedCV(outcomePerm, predictors, nrepeats, nfolds, maxComps, setSeed = round(runif(1, 1000, 4000)), classification = classification)[[1]]
    
    cat('\r iteration:', i, 'of', nperms)
  }
  
  return(savePerm)
}


# function for cross-validation from DNA methylation mean
meanCrossVal <- function(dv, data, nfolds = 5){
  
  dfOut <- numeric(nfolds)
  
  for(i in 1:nfolds){
    
    dataNoMiss <- data[!is.na(data[, dv]), ]
    folds <- createFolds(dataNoMiss[, dv], k = nfolds)
    trainSet <- dataNoMiss[-folds[[1]], ]
    testSet <- dataNoMiss[folds[[1]], ]
    
    model <- lm(data = trainSet, 
                as.formula(paste(dv, "~", "rowMeans.df_CpG_m.")))
    
    dfOut[i] <- 1 - sum((testSet[,dv] - predict(model, newdata = testSet))^2)/(var(testSet[,dv])*(nrow(testSet)-1))
    
  }
  
  return(dfOut)
  
}


# function that repeats meanCrossVal "nrepeats" times
repeatedMeanCrossVal <- function(dv, data, nrepeats, nfolds = 5){
  
  dfOut <- numeric(nrepeats)
  for(i in 1:nrepeats){
    
    dfOut[i] <- mean(meanCrossVal(dv, meanData))
    
  }
  return(mean(dfOut))
}

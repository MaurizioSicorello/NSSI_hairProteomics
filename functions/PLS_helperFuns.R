
# nested k-fold CV evaluation of PLS
PLSnestedCV <- function(outcome, predictors, nrepeats, nfolds, maxComps = 30, setSeed = 1000, classification = F){
  
  fitControl <- trainControl(method = "cv",   
                             number = nfolds)
  
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


# perform variable selection within cross-validation
PLS_selectPred <- function(data, n_select, nfolds = 5, selection){
  
  # output vector
  accOut <- numeric(nfolds)
  
  # cv split of data
  plsFolds <- createFolds(data[,2], k = nfolds)
  
  # start cv loop
  for(i in 1:nfolds){
    
    # train-test-split
    trainSet <- data[-plsFolds[[i]],-1]
    testSet <- data[plsFolds[[i]],-1]
    
    # perform PLS on training data
    plsda_train <- caret::plsda(x = trainSet[,-1], y = factor(trainSet[,1]), ncomp = 3)
    
    # pick best X predictors based on variable importance
    train_varImp <- data.frame("accession" = row.names(caret::varImp(plsda_train)), caret::varImp(plsda_train))
    varImp_sorted <- train_varImp[order(train_varImp$Overall, decreasing = TRUE), ]
    bestPreds <- varImp_sorted$accession[1:n_select]
    
    # indicator to include or remove best x predictors and group variable
    if(selection == "include"){
      varInd <- c(which(names(trainSet) %in% bestPreds))
    }else if(selection == "exclude"){
      varInd <- which(!(names(trainSet) %in% bestPreds))[-1]
    }else{
      warning("What are you even doing?! You are supposed to write 'exclude' or 'include'!")
    }
    
    # retrain PLS on selected training data
    plsda_retrain <- caret::plsda(x = trainSet[, varInd], y = factor(trainSet[,1]), ncomp = 3)
    
    # evaluate model on test data
    accOut[i] <- sum(ifelse(testSet[,1] == predict(plsda_retrain, newdata = testSet[, varInd]), 1, 0))/nrow(testSet)
    
  }
  
  return(accOut)
  
}


# repeat function "PLS_selectPred" and save results

PLS_selectPred_rep <- function(data, n_select, repeats = 5, nfolds = 5, selection){
  
  dfOut <- data.frame(
    matrix(nrow = nfolds, ncol = repeats)
  )
  
  for(i in 1:repeats){
    
    dfOut[,i] <- PLS_selectPred(data=data, n_select=n_select, nfolds = nfolds, selection=selection)
    
  }
  
  return(stack(dfOut)$values)
  
}

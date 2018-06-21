####################
## Function for visualizing univariate relations
####################

VisualizeRelation <- function(data=deer[["summer"]],model=GLMMs[["summer"]],predvar="dist_to_water",type="RF"){
  len <- 100
  
  isfac <- is.factor(data[[predvar]])
  
  dataclasses <- sapply(data,class)
  
  if(!isfac){
    if(type=="GLMM"){
      standvar <- sprintf("stand_%s",predvar)
    }else{
      standvar <- predvar
    }
    dim <- data[,standvar]
    range <- seq(min(dim),max(dim),length=len)
    
    realmean <- mean(data[,predvar])
    realsd <- sd(data[,predvar])
    
    newdata <- data.frame(temp=range)
    # head(newdata,50)
    names(newdata) <- c(standvar)
    if(type=="GLMM"){ 
      allvars <- names(model@frame)
    }else{
      allvars <- pred.names
    }
    othervars <- allvars[!allvars%in%c(standvar,"used")]
  }else{
    faclevs <- levels(data[[predvar]])
    newdata <- data.frame(temp=factor(faclevs,levels=faclevs))
    names(newdata) <- c(predvar)
    if(type=="GLMM"){ 
      allvars <- names(model@frame)
    }else{
      allvars <- pred.names
    }
    othervars <- allvars[!allvars%in%c(predvar,"used")]
  }
  
  var = othervars[2]
  for(var in othervars){
    thisvar <- data[,var]
    if(is.factor(thisvar)){
      tab <- table(thisvar)
      vals <- names(tab)
      levs <- levels(thisvar)
      mostcom <- vals[which.max(tab)]
      newvec <- factor(rep(mostcom,times=nrow(newdata)),levels=levs)
      newdata[,var] <- newvec
    }else{
      newdata[,var] <- mean(thisvar)
    }
  }
  
  if(type=="GLMM"){
    pred <- plogis(predict(model,newdata))
  }else{
    i=pred.names[3]
    for(i in pred.names){
      if(dataclasses[i]=="integer") newdata[,i] <- as.integer(round(newdata[,i]))
    }
    pred <- numeric(nrow(newdata))
    i=1
    for(i in 1:nrow(newdata)){
      pred[i]<-as.numeric(predict(model,newdata[i,],OOB=TRUE,type="prob")[[1]][,2])
    } 
  }
  
  if(!isfac){
    plot(range,pred,xlab=predictorNames[pred.names==predvar],ylab="Use probability",type="l",lwd=2,xaxt="n")
    ats <- seq(min(range),max(range),length=6)
    if(type=="GLMM"){
      axis(1,ats,labels = round(realmean+ats*realsd))
    }else{
      axis(1,ats,labels = round(ats))
    }
    rug(jitter(data[seq(1,nrow(data),50),standvar]), ticksize = 0.03, side = 1, lwd = 0.5, col = par("fg"))
  }else{
    par(mai=c(1.5,1,1,.2))
    plot(pred~newdata[,1],xlab="",main=predictorNames[pred.names==predvar],ylab="Use probability",lwd=2,las=2)
  }
}


####################
## Function for visualizing interactions
####################


VisualizeInteraction <- function(data=deer[["summer"]],model=GLMMs[["summer"]],var1="dist_to_water",var2="elevation",type="GLMM"){
  len <- 50
  
  dataclasses <- sapply(data,class)
  
  if(type=="GLMM"){
    standvar1 <- sprintf("stand_%s",var1)
    standvar2 <- sprintf("stand_%s",var2)
    realmean1 <- mean(data[,var1])
    realsd1 <- sd(data[,var1])
    realmean2 <- mean(data[,var2])
    realsd2 <- sd(data[,var2])
  }else{
    standvar1 <- var1
    standvar2 <- var2
  }
  
  firstdim <- data[,standvar1]
  seconddim <- data[,standvar2]
  range1 <- seq(min(firstdim),max(firstdim),length=len)
  range2 <- seq(min(seconddim),max(seconddim),length=len)
  newdata <- expand.grid(range1,range2)
  # head(newdata,50)
  names(newdata) <- c(standvar1,standvar2)
  
  if(type=="GLMM"){ 
    allvars <- names(model@frame)
  }else{
    allvars <- pred.names
  }
  othervars <- allvars[!allvars%in%c(standvar1,standvar2,"used")]
  
  var = othervars[2]
  for(var in othervars){
    thisvar <- data[,var]
    if(is.factor(thisvar)){
      tab <- table(thisvar)
      vals <- names(tab)
      levs <- levels(thisvar)
      mostcom <- vals[which.max(tab)]
      newvec <- factor(rep(mostcom,times=nrow(newdata)),levels=levs)
      newdata[,var] <- newvec
    }else{
      newdata[,var] <- mean(thisvar)
    }
  }
  
  if(type=="GLMM"){
    pred <- plogis(predict(model,newdata))
  }else{
    i=pred.names[3]
    for(i in pred.names){
      if(dataclasses[i]=="integer") newdata[,i] <- as.integer(round(newdata[,i]))
    }
    pred <- numeric(nrow(newdata))
    i=1
    for(i in 1:nrow(newdata)){
      pred[i]<-as.numeric(predict(model,newdata[i,],OOB=TRUE,type="prob")[[1]][,2])
    } 
  }
  
  predmat <-  matrix(pred,nrow=len,ncol=len)
  
  par(mai=c(0,0,0,0))
  
  if(type=="GLMM"){
    persp(realmean1+realsd1*range1,realmean2+realsd2*range2,predmat,xlab=var1,ylab=var2,theta = 55, phi = 40, r = sqrt(10), d = 3, 
          ticktype = "detailed", mgp = c(4, 1, 0))
  }else{
    persp(range1,range2,predmat,xlab=var1,ylab=var2,theta = 55, phi = 40, r = sqrt(10), d = 3, 
          ticktype = "detailed", mgp = c(4, 1, 0))
  }
  
}


#####################
## CROSS VALIDATION
#####################

n.folds=3

type= "GLMM" #"RF"
season="summer"
fullmodel=GLMMs[[season]] #RFs[[season]]

CrossValidateByDeer <- function(n.folds,season="summer",type="GLMM",plot=F){
  uniquedeer <- as.character(unique(deer[[season]]$altid))  # list of all unique animals
  ndeer <- length(uniquedeer)  # total number of inds
  folds_df <- data.frame( 
    deer = uniquedeer,
    fold = rep_len(1:n.folds,ndeer)
  )
  foldVector <- folds_df$fold[match(as.character(deer[[season]]$altid),folds_df$deer)]
  
  predictCols <- which(names(deer[[season]])%in%pred.names)
  
  if(type=="RF"){
    fullmodel<-RFs[[season]]
  }else{
    fullmodel <- GLMMs[[season]]
  }
  
  CVresults <- list()
  CVresults$CVpred <- numeric(nrow(deer[[season]]))
  CVresults$realpred <- numeric(nrow(deer[[season]])) 
  CVresults$observed <- numeric(nrow(deer[[season]]))
  
  if(type=="RF"){
    response="used_fac"    #"resp_factor"
  }else{
    response="used"    #"resp_factor"
  }
  
  counter = 1
  
  i=n.folds
  for(i in 1:n.folds){
    if(type=="RF"){
      model <- cforest(formula1, data = deer[[season]][which(foldVector!=i),], controls=cforestControl) 
      predict_CV  <- predict(model,newdata=deer[[season]][which(foldVector==i),],type="prob") 
      predict_real  <-  predict(fullmodel,newdata=deer[[season]][which(foldVector==i),],type="prob")
      REAL <- deer[[season]]$used[which(foldVector==i)]
      j=1
      for(j in 1:length(which(foldVector==i))){
        CVresults$CVpred[counter] <- as.numeric(predict_CV[[j]][,2])
        CVresults$observed[counter] <-  as.numeric(REAL[j])      
        CVresults$realpred[counter] <- as.numeric(predict_real[[j]][,2])   
        counter = counter + 1  
      }
    }else{
      
      model <- glmer(attributes(GLMMs[[season]]@frame)$formula, 
                       family="binomial", data=deer[[season]][which(foldVector!=i),],na.action="na.fail") 
      
      CVresults$CVpred[which(foldVector==i)]  <- plogis(predict(model,newdata=deer[[season]][which(foldVector==i),],allow.new.levels = TRUE)) 
      CVresults$realpred[which(foldVector==i)] <-  predict(fullmodel,newdata=deer[[season]][which(foldVector==i),],allow.new.levels = TRUE)
      CVresults$observed[which(foldVector==i)] <- deer[[season]]$used[which(foldVector==i)]      
    }
    cat(sprintf("fold %s out of %s\n",i,n.folds))
  }
  
  CVresults$CV_RMSE = sqrt(mean((CVresults$observed-CVresults$CVpred)^2))       # root mean squared error for holdout samples in 10-fold cross-validation ...
  CVresults$real_RMSE = sqrt(mean((CVresults$observed-CVresults$realpred)^2))   # root mean squared error for residuals from final model
  
  # realprediction <- predict(rf_model1,newdata=df,type="prob")
  
  if(plot){
    graphics.off()
    par(mfrow=c(2,1))
    par(ask=T)
    pred <- prediction(CVresults$CVpred,CVresults$observed)     # for holdout samples in cross-validation
    perf <- performance(pred,"tpr","fpr")
    auc <- performance(pred,"auc")
    plot(perf, main=sprintf("%s Cross Validation",type))
    
    text(.9,.1,paste("AUC = ",round(auc@y.values[[1]],2),sep=""))
    
    pred <- prediction(CVresults$realpred,CVresults$observed)     # for final model
    perf <- performance(pred,"tpr","fpr")
    auc <- performance(pred,"auc")
    plot(perf, main=sprintf("%s Full Model",type))
    text(.9,.1,paste("AUC = ",round(auc@y.values[[1]],2),sep=""))
  }else{
    pred <- prediction(CVresults$CVpred,CVresults$observed)     # for holdout samples in cross-validation
    perf <- performance(pred,"tpr","fpr")
    CVresults$auc_CV <- performance(pred,"auc")
    pred <- prediction(CVresults$realpred,CVresults$observed)     # for final model
    perf <- performance(pred,"tpr","fpr")
    CVresults$auc_full <- performance(pred,"auc")
  }   
  
  # COHEN KAPPA statistics
  
  thresholds <- seq(0.01,0.99,length=101)   # "artificial" extinction thresholds across which to examine performance
  kappa <- numeric(length(thresholds))
  for(i in 1:length(thresholds)){
    trueLabels <- CVresults$observed
    predLabels <- ifelse(CVresults$CVpred>=thresholds[i],1,0)
    tot <- length(CVresults$observed)
    tp <- length(which((trueLabels==1)&(predLabels==1)))  
    tn <- length(which((trueLabels==0)&(predLabels==0)))
    fp <- length(which((trueLabels==0)&(predLabels==1)))
    fn <- length(which((trueLabels==1)&(predLabels==0)))
    pr_agree <- (tp+tn)/tot    # overall agreement, or accuracy
    pr_agree_rand <- ((tp+fn)/tot)*((tp+fp)/tot)+((fn+tn)/tot)*((fp+tn)/tot)
    kappa[i] <- (pr_agree-pr_agree_rand)/(1-pr_agree_rand)
  }
  #plot(thresholds,kappa,type="l",xlab="Threshold", ylab="Cohen's Kappa", main="Holdout sample performance")
  
  # find threshold value associated with highest Kappa for C-V data
  
  CVresults$cutoff_CV <- thresholds[which.max(kappa)]   # max kappa cutoff
  
  kappa <- numeric(length(thresholds)) 
  for(i in 1:length(thresholds)){
    trueLabels <- CVresults$observed
    predLabels <- ifelse(CVresults$realpred>=thresholds[i],1,0)    
    tot <- length(CVresults$observed)
    tp <- length(which((trueLabels==1)&(predLabels==1)))  
    tn <- length(which((trueLabels==0)&(predLabels==0)))
    fp <- length(which((trueLabels==0)&(predLabels==1)))
    fn <- length(which((trueLabels==1)&(predLabels==0)))
    pr_agree <- (tp+tn)/tot    # overall agreement, or accuracy
    pr_agree_rand <- ((tp+fn)/tot)*((tp+fp)/tot)+((fn+tn)/tot)*((fp+tn)/tot)
    kappa[i] <- (pr_agree-pr_agree_rand)/(1-pr_agree_rand)
  }
  #plot(thresholds,kappa,type="l",xlab="Threshold", ylab="Cohen's Kappa", main="Performance: full model")
  
  CVresults$cutoff_full <- thresholds[which.max(kappa)]   # max kappa cutoff
  
  ### display confusion matrix and kappa for a single threshold
  trueLabels <- CVresults$observed
  predLabels <- ifelse(CVresults$CVpred>=CVresults$cutoff_CV,1,0)    
  tot <- length(CVresults$observed)
  tp <- length(which((trueLabels==1)&(predLabels==1)))  
  tn <- length(which((trueLabels==0)&(predLabels==0)))
  fp <- length(which((trueLabels==0)&(predLabels==1)))
  fn <- length(which((trueLabels==1)&(predLabels==0)))
  pr_agree <- (tp+tn)/tot    # overall agreement, or accuracy
  pr_agree_rand <- ((tp+fn)/tot)*((tp+fp)/tot)+((fn+tn)/tot)*((fp+tn)/tot)
  CVresults$bestkappa_CV <- (pr_agree-pr_agree_rand)/(1-pr_agree_rand)
  
  CVresults$confusionmat <- matrix(c(tn,fn,fp,tp),nrow=2,ncol=2)
  rownames(CVresults$confusionmat) <- c("Actual not used","Actual used")
  colnames(CVresults$confusionmat) <- c("Predicted not used","Predicted used")
  
  CVresults$sensitivity <- tp/(tp+fn)
  CVresults$specificity <- tn/(tn+fp)
  CVresults$toterror <- (fn+fp)/tot
  
  CVresults$CVpred[which(CVresults$CVpred==1)] <- 0.999999
  CVresults$CVpred[which(CVresults$CVpred==0)] <- 0.000001
  CVresults$realpred[which(CVresults$realpred==1)] <- 0.999999
  CVresults$realpred[which(CVresults$realpred==0)] <- 0.000001
  
  realdata = CVresults$observed
  fit_deviance_CV <- mean(-2*(dbinom(CVresults$observed,1,CVresults$CVpred,log=T)-dbinom(realdata,1,realdata,log=T)))
  fit_deviance_real <- mean(-2*(dbinom(CVresults$observed,1,CVresults$realpred,log=T)-dbinom(realdata,1,realdata,log=T)))
  null_deviance <- mean(-2*(dbinom(CVresults$observed,1,mean(CVresults$observed),log=T)-dbinom(realdata,1,realdata,log=T)))
  CVresults$deviance_explained_CV <- (null_deviance-fit_deviance_CV)/null_deviance   # based on holdout samples
  CVresults$deviance_explained_real <- (null_deviance-fit_deviance_real)/null_deviance   # based on full model...
  
  return(CVresults)
}



PlotPerformance <- function(CVresults){
    graphics.off()
    par(mfrow=c(2,1))
    #par(ask=T)
    pred <- prediction(CVresults$CVpred,CVresults$observed)     # for holdout samples in cross-validation
    perf <- performance(pred,"tpr","fpr")
    auc <- performance(pred,"auc")
    plot(perf, main=sprintf("%s Cross Validation",type))
    
    text(.9,.1,paste("AUC = ",round(auc@y.values[[1]],2),sep=""))
    
    pred <- prediction(CVresults$realpred,CVresults$observed)     # for final model
    perf <- performance(pred,"tpr","fpr")
    auc <- performance(pred,"auc")
    plot(perf, main=sprintf("%s Full Model",type))
    text(.9,.1,paste("AUC = ",round(auc@y.values[[1]],2),sep=""))
}








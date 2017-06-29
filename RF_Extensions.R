#########################################
###############  Functions for extending and visualizing results from the Random Forest algorithm in the "party" package
###############      K. Shoemaker 15 Jan 2014


##############################
#######  LOAD PACKAGES

loadPackage <- function(pkg){

  if(pkg %in% rownames(installed.packages()) == FALSE) {suppressMessages(suppressWarnings(install.packages(pkg)))}
  eval(parse(text=sprintf("suppressMessages(suppressWarnings(require(%s)))",pkg)), envir= .GlobalEnv)

}


loadPackage("party")
loadPackage("rms")
loadPackage("ROCR")

# suppressWarnings(require(gbm))
# suppressWarnings(require(dismo))
# suppressWarnings(require(scales))
# suppressWarnings(require(ggplot2))
# suppressWarnings(require(lattice))
# suppressWarnings(require(ipred))
# suppressWarnings(require(party))
# suppressWarnings(require(RColorBrewer))
# suppressWarnings(require(tree))
# suppressWarnings(require(languageR))
# suppressWarnings(require(rms))
# suppressWarnings(require(ROCR))
# suppressWarnings(require(rpart))



################################################################
#####################################  LOAD NEW FUNCTIONS


#################################################
######################  NEW FUNCTION FOR PLOTTING SUBTREES FROM BOOSTED REGRESSION MODELS  ...

## plot_subtree: function to visualize the subtrees for BRT analyses...
##    KTS 7/1/2012  

plot_subtree = function(brt_model, data, pred_cols, lr, scenario_name, response_name, tree_ndx = 1000){
    #  which component tree to visualize??
	sepdist <- 5   
	frame <- pretty.gbm.tree(brt_model,i.tree=tree_ndx)
	pred_names_ndx <- names(data)[pred_cols]
	splitVars <- pred_names_ndx[frame$SplitVar[which(frame$SplitVar!=-1)]+1]
	nsplits <- length(splitVars)                    
	yoffset <- sum(frame$ErrorReduction*1.3)*.02
	graphics.off()
	par(mai=c(.1,.1,.5,.1))
      	# set plot area
	plot(1,1,type="n",xlim=c(-10,20),ylim=c(0,sum(frame$ErrorReduction*1.3)),bty="n",
            xaxt="n",yaxt="n",xlab="",ylab="",main=paste(response_name,scenario_name,"tree",tree_ndx,sep=" "))
	nodecoords <- list(left=list(x=numeric(nsplits),y=numeric(nsplits)),right=list(x=numeric(nsplits),y=numeric(nsplits)))
	nodecoords_base <- c(5,sum(frame$ErrorReduction)*1.2)  ### list(x=c(5,4),y=c(sum(frame$ErrorReduction),sum(frame$ErrorReduction[-1])))  #base...
	nodecoords$left$x[1]=nodecoords_base[1]-sepdist; nodecoords$left$y[1]=nodecoords_base[2]-frame$ErrorReduction[1]
	nodecoords$right$x[1]=nodecoords_base[1]+sepdist; nodecoords$right$y[1]=nodecoords_base[2]-frame$ErrorReduction[1]
	text(nodecoords_base[1],nodecoords_base[2]+yoffset,paste(splitVars[1],"<",round(frame$SplitCodePred[1],2),sep=""))
	points(c(nodecoords_base[1],nodecoords$left$x[1]),c(nodecoords_base[2],nodecoords$left$y[1]),type="l")
	points(c(nodecoords_base[1],nodecoords$right$x[1]),c(nodecoords_base[2],nodecoords$right$y[1]),type="l")
	counter1=1
	splitndx <- numeric(nsplits)
	splitndx[counter1] <- 1      # loop variable at each split...
	for(i in 2:nrow(frame)){
	  if(frame$SplitVar[i]!=-1){   # if not a terminal leaf...
	    counter1=counter1 + 1
	    splitndx[counter1] <- i   # what is the value of the loop variable at this split?
	    sepdist = sepdist*0.8
	    origin_node = which((frame$LeftNode+1)==i|(frame$RightNode+1)==i) 
	    ndx <- which(origin_node==splitndx)  
	    left <- ifelse((frame$LeftNode[origin_node]+1)==i,1,2)
	    text(nodecoords[[left]]$x[ndx],nodecoords[[left]]$y[ndx]+yoffset,paste(splitVars[counter1],"<",round(frame$SplitCodePred[i],2),sep=""))   # split variable...
	    nodecoords$left$x[counter1]=nodecoords[[left]]$x[ndx]-sepdist
	    nodecoords$left$y[counter1]=nodecoords[[left]]$y[ndx]-frame$ErrorReduction[i]
	    nodecoords$right$x[counter1]=nodecoords[[left]]$x[ndx]+sepdist
	    nodecoords$right$y[counter1]=nodecoords[[left]]$y[ndx]-frame$ErrorReduction[i]
	    points(c(nodecoords[[left]]$x[ndx],nodecoords$left$x[counter1]),c(nodecoords[[left]]$y[ndx],nodecoords$left$y[counter1]),type="l")
	    points(c(nodecoords[[left]]$x[ndx],nodecoords$right$x[counter1]),c(nodecoords[[left]]$y[ndx],nodecoords$right$y[counter1]),type="l")   
	  }
	  if(frame$SplitVar[i]==-1){   # if a terminal leaf...
	    origin_node = which((frame$LeftNode+1)==i|(frame$RightNode+1)==i)
	    if(length(origin_node)>0){ 
	      ndx <- which(origin_node==splitndx)  
	      left <- ifelse((frame$LeftNode[origin_node]+1)==i,1,2)
	      text(nodecoords[[left]]$x[ndx],nodecoords[[left]]$y[ndx]-yoffset,paste(round(eval(parse(text=workingModel))$initF+(frame$Prediction[i]/lr),4),sep=""))   # split variable...   
	    }
	  }
	}
}


#############################################
#################### NEW FUNCTION: "find_fraction" for determining the fraction
  #  of the total observations to sample without replacement for each tree in the random forest
  #  or for each term of the boosted model 


find_fraction <- function(dataframe){
  len <- nrow(dataframe)
  nspec <- length(unique(dataframe$spec))  # number of geographic configurations...
  target <- nspec-1
  testarray <- rep(c(1:nspec),each=floor(len/nspec))  
  obj <- optimize(subsample_fn,c(0.005,0.2),target=target,array=testarray)   #(par,target,testarray)
  return(obj$minimum)
}

subsample_fn <- function(par,target,array){
  dif <- numeric(200)
  for(i in 1:200){
    dif[i] <- abs(target-length(unique(sample(array,length(array)*par,replace=F))))
  }
  a <- mean(dif)
  return(a)
}

#################################################
######################  NEW GBM.PERSPEC FUNCTION...

# rewrite the perspective plotting function

gbm.perspec <- function (gbm.object, x = 1, y = 2, pred.means = NULL, x.label = NULL, 
    x.range = NULL, y.label = NULL, z.label = "fitted value", 
    y.range = NULL, z.range = NULL, leg.coords = NULL, ticktype = "detailed", 
    theta = 55, phi = 40, smooth = "none", mask = FALSE, perspective = TRUE, 
    ...) 
{
    if (!require(gbm)) {
        stop("you need to install the gbm package to use this function")
    }
    if (!require(splines)) {
        stop("you need to install the splines package to use this function")
    }
    gbm.call <- gbm.object$gbm.call
    gbm.x <- gbm.call$gbm.x
    n.preds <- length(gbm.x)
    gbm.y <- gbm.call$gbm.y
    pred.names <- gbm.call$predictor.names
    family = gbm.call$family
    have.factor <- FALSE
    x.name <- gbm.call$predictor.names[x]
    if (is.null(x.label)) {
        x.label <- gbm.call$predictor.names[x]
    }
    y.name <- gbm.call$predictor.names[y]
    if (is.null(y.label)) {
        y.label <- gbm.call$predictor.names[y]
    }
    data <- eval(parse(text = gbm.call$dataframe))[, gbm.x]
    n.trees <- gbm.call$best.trees
    if (is.vector(data[, x])) {
        if (is.null(x.range)) {
            x.var <- seq(min(data[, x], na.rm = T), max(data[, 
                x], na.rm = T), length = 50)
        } else {
            x.var <- seq(x.range[1], x.range[2], length = 50)
        }
    } else {
        x.var <- names(table(data[, x]))
        have.factor <- TRUE
    }
    if (is.vector(data[, y])) {
        if (is.null(y.range)) {
            y.var <- seq(min(data[, y], na.rm = T), max(data[, 
                y], na.rm = T), length = 50)
        } else {
            y.var <- seq(y.range[1], y.range[2], length = 50)
        }
    } else {
        y.var <- names(table(data[, y]))
        #if (have.factor) {
        #    stop("at least one marginal predictor must be a vector!")
        #}
        #else {
            have.factor <- TRUE
        #}
    }
    pred.frame <- expand.grid(list(x.var, y.var))
    names(pred.frame) <- c(x.name, y.name)
    pred.rows <- nrow(pred.frame)
    #if (have.factor) {
    #    if (is.factor(pred.frame[, 2])) {
    #        pred.frame <- pred.frame[, c(2, 1)]
    #        x.var <- y.var
    #    }
    #}
    j <- 3
    for (i in 1:n.preds) {
        if (i != x & i != y) {
            if (is.vector(data[, i])) {
                m <- match(pred.names[i], names(pred.means))
                if (is.na(m)) {
                  pred.frame[, j] <- mean(data[, i], na.rm = T)
                }
                else pred.frame[, j] <- pred.means[m]
            }
            if (is.factor(data[, i])) {
                m <- match(pred.names[i], names(pred.means))
                temp.table <- table(data[, i])
                if (is.na(m)) {
                  pred.frame[, j] <- rep(names(temp.table)[2], 
                    pred.rows)
                }
                else {
                  pred.frame[, j] <- pred.means[m]
                }
                pred.frame[, j] <- factor(pred.frame[, j], levels = names(temp.table))
            }
            names(pred.frame)[j] <- pred.names[i]
            j <- j + 1
        }
    }
    prediction <- predict.gbm(gbm.object, pred.frame, n.trees = n.trees, 
        type = "response")
    if (smooth == "model") {
        pred.glm <- glm(prediction ~ ns(pred.frame[, 1], df = 8) * 
            ns(pred.frame[, 2], df = 8), data = pred.frame, family = poisson)
        prediction <- fitted(pred.glm)
    }
    max.pred <- max(prediction)
    min.pred <- min(prediction)
    cat("maximum value = ", round(max.pred, 2), "\n")
    if (is.null(z.range)) {
        if (family == "bernoulli") {
            z.range <- c(0, 1)
        } else if (family == "poisson") {
            z.range <- c(0, max.pred * 1.1)
        } else {
            #z.min <- min(data[, y], na.rm = T)
            #z.max <- max(data[, y], na.rm = T)
            #z.delta <- z.max - z.min
            #z.range <- c(z.min - (1.1 * z.delta), z.max + (1.1 * 
            #    z.delta))
            z.range <- c(min.pred,max.pred)
        }
    }
    if (have.factor == FALSE) {
        pred.matrix <- matrix(prediction, ncol = 50, nrow = 50)
        if (smooth == "average") {
            pred.matrix.smooth <- pred.matrix
            for (i in 2:49) {
                for (j in 2:49) {
                  pred.matrix.smooth[i, j] <- mean(pred.matrix[c((i - 
                    1):(i + 1)), c((j - 1):(j + 1))])
                }
            }
            pred.matrix <- pred.matrix.smooth
        }
        if (mask) {
            mask.trees <- gbm.object$gbm.call$best.trees
            point.prob <- predict.gbm(gbm.object[[1]], pred.frame, 
                n.trees = mask.trees, type = "response")
            point.prob <- matrix(point.prob, ncol = 50, nrow = 50)
            pred.matrix[point.prob < 0.5] <- 0
        }
        if (!perspective) {
            image(x = x.var, y = y.var, z = pred.matrix, zlim = z.range)
        }
        else {
            persp(x = x.var, y = y.var, z = pred.matrix, zlim = z.range, 
                xlab = x.label, ylab = y.label, zlab = z.label, 
                theta = theta, phi = phi, r = sqrt(10), d = 3, 
                ticktype = ticktype, mgp = c(4, 1, 0), ...)
        }
    }
    if (have.factor) {
        factor.list <- names(table(pred.frame[, 1]))
        n <- 1
        if (is.null(z.range)) {
            vert.limits <- c(min.pred * 1.1, max.pred * 1.1)
        } else {
            vert.limits <- z.range
        }
        plot(as.numeric(pred.frame[pred.frame[, 1] == factor.list[1], 2]), 
            prediction[pred.frame[, 1] == factor.list[1]], type = "n", pch=NA,
            ylim = vert.limits, xlab = y.label, ylab = z.label,xaxt="n")
        if(is.factor(pred.frame[, 2])) xaxlab <- names(table(pred.frame[, 2]))
        if(is.vector(pred.frame[, 2])) xaxlab <- seq(min(pred.frame[, 2]),max(pred.frame[, 2]),length=5)
        axis(1,at=c(1:length(xaxlab)),labels=xaxlab)
        for (i in 1:length(factor.list)) {
            factor.level <- factor.list[i]
            #if()
            lines(pred.frame[pred.frame[, 1] == factor.level, 
                2], prediction[pred.frame[, 1] == factor.level], 
                lty = i)
        }
        if (is.null(leg.coords)) {
            #x.max <- max(pred.frame[, 2])
            #x.min <- min(pred.frame[, 2])
            #x.range <- x.max - x.min
            #x.pos <- c(x.min + (0.02 * x.range), x.min + (0.3 * 
            #    x.range))
            #y.max <- max(prediction)
            #y.min <- min(prediction)
            #y.range <- y.max - y.min
            #y.pos <- c(y.min + (0.8 * y.range), y.min + (0.95 * 
            #    y.range))
            legend(locator(1), factor.list, lty = c(1:length(factor.list)),    #x = x.pos, y = y.pos, 
                bty = "n")
        }
        else {
            legend(x = leg.coords[1], y = leg.coords[2], factor.list,     #  
                lty = c(1:length(factor.list)), bty = "n")
        }
    }
}



#######################################################################
########################################################################
#        NEW FUNCTION: makefoldvec: for cross validation in RF and BRT

   # define a fold vector so that the holdout samples are spatially uncorrelated from the training sets
make_foldvec <- function(n.folds,foldvar){
  ##n.folds = 12
  nspecies <- max(as.numeric(as.factor(foldvar)))
  tempfolds1 <- numeric(nspecies)
  temp2 <- c(1:nspecies)
  reps <- floor(nspecies/n.folds)
  counter1 = 1
  counter2 = reps
  foldlink = numeric(nspecies)
  foldlink[1:(n.folds*counter2)] <- rep(c(1:n.folds),each=counter2)
  foldlink[which(foldlink==0)] <- n.folds
  i=1
  for(i in 1:(n.folds)){
    temp3 <- sample(temp2,reps,replace=F)
    tempfolds1[counter1:counter2] <- temp3
    temp2 <- temp2[-which(temp2%in%temp3)]
    counter1 = counter1 + reps
    counter2 = counter2 + reps
  }
  if(length(which(tempfolds1==0))>0) tempfolds1[which(tempfolds1==0)] <- temp2
  foldVector <- foldlink[match(as.numeric(as.factor(foldvar)),tempfolds1)]
  if(n.folds == length(unique(foldvar))) foldVector <- as.numeric(as.factor(foldvar))
  return(foldVector)
}


#############################################################
##########  FUNCTION makeResponseMatrix: for use in Random Forest univariate and interaction plots
#####  makes a new data frame of response variables for prediction, holding all
####   variables constant at mean levels except for the variable(s) of interest...

#attach(oc)

## note factor predictors are ORD, Island, TrophicD, HabMode, ActiveCycle
#pred.names=c("logBodyMass" , "logAllGeogRange" , "P" , "ORD", "Island", "TrophicD" , "HabMode" , "ActivCycle" , "logHomeRange.Indiv.km2" , "logSocGrpSz" , "logPopDen" , "HuPopDen.Change" , "HuPopDen.5p.n.km2" , "logMeanPrecip" , "Temp.Mean.01degC" , "logAET" , "logMean.NPP" , "Mean.HII" , "Max.HII" , "GR.MaxLat.dd" , "GR.MaxLong.dd" , "GR.MinLong.dd" , "GR.MinLat.dd")

makeResponseMatrix <- function(k,RFObject,pred.data,n.preds,predictors,data,rf.x){
  
  factResp = FALSE
  if(is.numeric(pred.data)){
    # initialize the predictor data frame...
    subset2 <- data.frame(nullcol = rep(NA,times=100)) 
    
    # loop through predictor variables  ###change 100 to change bins
    for(v in 1:n.preds){         
      if(v==k){
        subset2 <- cbind(subset2,seq(min(pred.data, na.rm=T),max(pred.data,na.rm=T),length=100))    # if the variable of interest, use the observed values
        names(subset2)[v+1] <- predictors[v]
        #subset2 <- cbind(subset2,)
      } else{
        if(is.numeric(data[,rf.x[v]])){
          subset2 <- cbind(subset2,rep(mean(data[,rf.x[v]],na.rm=T),times=100))   # if NOT the variable of interest, replicate the mean value
        }
        if(is.factor(data[,rf.x[v]])){
          subset2 <- cbind(subset2,factor(rep(names(which.max(table(data[,rf.x[v]]))),times=100),levels=levels(data[,rf.x[v]])))
        }
        names(subset2)[v+1] <- predictors[v]
      }
    }
    
    # compile a data frame (subset2) for prediction. 
    # generate a prediction from random forest model
    
    # check to make sure the predictors are all of the same class
    for(v in 2:(n.preds+1)){
      if(class(data[,rf.x[v-1]])!=class(subset2[,v])) subset2[,v] <- eval(parse(text=paste("as.",class(data[,rf.x[v-1]]),"(subset2[,v])",sep="")))
    }
    
    if(!factResp){
      predictions <- predict(RFObject,newdata=subset2)
    } 
    if(is.factor(predictions)) factResp=TRUE
    if(factResp){
      predictions <- numeric(length(predictions)) #    # 
      for(a in 1:length(predictions)){
        predictions[a] <- as.numeric(predict(RFObject,newdata=subset2[a,],type="prob")[[1]][,2])
      }
    }
    
    
    df <- data.frame(v1=seq(min(pred.data,na.rm=T),max(pred.data,na.rm=T),length=100),v2=predictions)
    names(df) <- c(predictors[k],"y")
  }
  
  if(is.factor(pred.data)){
    # initialize the predictor data frame...
    subset2 <- data.frame(nullcol = rep(NA,times=length(levels(pred.data)))) 
    
    # loop through predictor variables
    for(v in 1:n.preds){       	
      if(v==k){
        #temp <- rep(NA,times=length(levels(pred.data)))
        subset2 <- cbind(subset2,factor(levels(pred.data),levels=levels(pred.data)))    # if the variable of interest, use the observed values
        names(subset2)[v+1] <- predictors[v]
      } else{
        if(is.numeric(data[,rf.x[v]])){
          subset2 <- cbind(subset2,rep(mean(data[,rf.x[v]],na.rm=T),times=length(levels(pred.data))))   # if NOT the variable of interest, replicate the mean value
        }
        if(is.factor(data[,rf.x[v]])){
          subset2 <- cbind(subset2,factor(rep(names(which.max(table(data[,rf.x[v]]))),times=length(levels(pred.data))),levels=levels(data[,rf.x[v]])))
        }
        names(subset2)[v+1] <- predictors[v]
      }
    }
    
    # compile a data frame (subset2) for prediction. 
    # generate a prediction from random forest model
    
    # check to make sure the predictors are all of the same class
    for(v in 2:(n.preds+1)){
      if(class(data[,rf.x[v-1]])!=class(subset2[,v])) subset2[,v] <- eval(parse(text=paste("as.",class(data[,rf.x[v-1]]),"(subset2[,v])",sep="")))
    }
    
    if(!factResp){
      predictions <- predict(RFObject,newdata=subset2)
    } 
    if(is.factor(predictions)) factResp=TRUE
    if(factResp){
      predictions <- numeric(length(predictions)) #    # 
      for(a in 1:length(predictions)){
        predictions[a] <- as.numeric(predict(RFObject,newdata=subset2[a,],type="prob")[[1]][,2])
      }
    }
    
    df <- data.frame(v1=factor(levels(pred.data),levels=levels(pred.data)),v2=predictions)
    names(df) <- c(predictors[k],"y")
  }
  
  return(df) 
}



############################################
##################  NEW FUNCTION: RF_UnivariatePlots
######  For use with random forest object from cforest
######  displays univariate plots, analogous to "gbm.plot" in the "dismo" package  


RF_UnivariatePlots <- function(object, varimp, data,predictors,labels,allpredictors,plot.layout,plot=T){
  
  # gbm.object = eval(parse(text=workingModel))
 # graphics.off()
  if(length(predictors)>1) ask=TRUE else ask=FALSE
  if(plot) par(mfrow=plot.layout,ask=ask)
  smooth = FALSE
  rug = FALSE
  common.scale = TRUE
  write.title = FALSE 
  y.label = "fitted function"
  x.label = NULL
  show.contrib = TRUE
  n.plots = length(predictors)  # min(prod(plot.layout),length(predictors))
  
  returnObj=list()
  
  rf.x <- match(predictors,names(data))   #   which(names(data)%in%predictors) 
  rf.x2 <- match(allpredictors,names(data))   #   which(names(data)%in%allpredictors)    
  #pred.names = names(data[rf.x])  #  predictors
  #n.plots = length(pred.names) 
  #variable.no = 0
  
  original.response <- object@data@get("response")   # OR...   data.cforest@responses@test_trafo
  response.name <- names(object@responses@variables)  # gbm.call$response.name
  if(plot){
    max.plots <- plot.layout[1] * plot.layout[2]
    plot.count <- 0
    n.pages <- 1
    max.vars <- length(predictors)
    if (n.plots > max.vars) {
     n.plots <- max.vars
      warning("reducing no of plotted predictors to maximum available (", 
             max.vars, ")")
    }
  } else{
    max.plots <- length(predictors)
    plot.count <- 0
  }
    
  predictors3 <- list(rep(NA, n.plots))
  responses <- list(rep(NA, n.plots))
  ymin <- numeric(n.plots)
  ymax <- numeric(n.plots)
  RFPredContributions <- names(varimp)[order(varimp,decreasing=T)]
  
  for (j in c(1:n.plots)) { 
    k <- NA     # k is index of "predictors" to use for this plot
    l <- j      # l is index of RFPredContributions to use for this plot
		    # j is index of plot
    while(is.na(k)){
      k <- match(RFPredContributions[l],predictors)     # variable to plot as predictor
      l <- l + 1
    }
    l <- l - 1
    if (is.null(x.label)){ 
      var.name <- RFPredContributions[l]
    } else var.name <- x.label
    pred.data <- data[,which(names(data)==predictors[k])] 

    
    response.matrix <- makeResponseMatrix(k=which(allpredictors==predictors[k]),RFObject=object,
                        pred.data=pred.data,n.preds=length(allpredictors),
                        predictors=allpredictors,data=data,rf.x=rf.x2)   # prediction, holding all other predictors at mean values...
    predictors3[[j]] <- response.matrix[, 1]
    
    #predictors3[[j]] <- factor(predictors3[[j]], levels = levels(pred.data))
    
    if(plot){
      responses[[j]] <- response.matrix[, 2] - mean(response.matrix[,2])
    } else{
      responses[[j]] <- response.matrix[, 2]
    }
    if (j == 1) {
      ymin[j] = min(responses[[j]])
      ymax[j] = max(responses[[j]])
    }else {
      ymin[j] = min(responses[[j]])
      ymax[j] = max(responses[[j]]) #ymin = min(ymin, min(responses[[j]]))
      #ymax = max(ymax, max(responses[[j]]))
    }
  }
  
  for (j in c(1:n.plots)) {
    if (plot.count == max.plots) {
      plot.count = 0
      n.pages <- n.pages + 1
    }
    plot.count <- plot.count + 1
    #if (n.plots == 1) {
    #    k <- match(pred.names[variable.no],RFPredContributions)
    #    if (show.contrib) {
    #        x.label <- paste(var.name, "  (", round(varimp[k], 1), ")", sep = "")
    #    }
    #}else {
    
    k <- NA
    l <- j
    while(is.na(k)){
      k <- match(RFPredContributions[l],predictors)     # variable to plot as predictor
      l <- l + 1
    }
    l <- l - 1

    var.name <- labels[k]   #predictors[k]
    #if (show.contrib) {
    x.label <- paste(var.name, "  (", round(varimp[order(varimp,decreasing=T)][l], 3), ")", sep = "")
    #}else x.label <- var.name
    #}
    #if (common.scale) {
    #inflationFactor <- abs(max(data.cforest.varimp[order(data.cforest.varimp,decreasing=T)])/data.cforest.varimp[order(data.cforest.varimp,decreasing=T)][j])

    if(is.factor(predictors3[[j]])){
      par(las=2) 
    }else{
      par(las=1)
    }

    if(plot){ 
      plot(predictors3[[j]], responses[[j]], ylim = c(ymin[j],
                                                    ymax[j]), type = "l", xlab = x.label, ylab = y.label)
    } else{
      returnObj[[j]] <-  data.frame(pred=predictors3[[j]],resp=responses[[j]]) 
      names(returnObj[[j]])[1] <- RFPredContributions[j]
    }
      
      #}else {
    #    plot(predictors3[[j]], responses[[j]], type = "l", 
    #        xlab = x.label, ylab = y.label)
    #}
    
    #if (smooth & is.vector(predictors3[[j]])) {
    #    temp.lo <- loess(responses[[j]] ~ predictors3[[j]], 
    #        span = 0.3)
    #    lines(predictors3[[j]], fitted(temp.lo), lty = 2, 
    #        col = 2)
    #}
    
    if (plot.count == 1) {
      if (write.title) {
        title(paste(response.name, " - page ", n.pages, 
                    sep = ""))
      }
      #if (rug & is.vector(data[, rf.x[variable.no]])) {
      #    rug(quantile(data[, rf.x[variable.no]], 
      #      probs = seq(0, 1, 0.1), na.rm = TRUE))
      #}
    }else {
      if (write.title & j == 1) {
        title(response.name)
      }
      #if (rug & is.vector(data[, rf.x[k]])) {
      #    rug(quantile(data[, rf.x[k]], probs = seq(0, 
      #      1, 0.1), na.rm = TRUE))
      #}
    }
  }
  # par(op)
  
  if(!plot) return(returnObj)
  
}

#################

#######################################################
###   NEW FUNCTION: "RF_FindInteractions" 
#                            Find interactions from Random forest algorithm


RF_FindInteractions <- function(object,data,predictors,logit=F){

         # gbm.object <- eval(parse(text=workingModel))
         # RFObject <- data.cforest
         # n.trees <- gbm.call$best.trees
         # depth <- gbm.call$interaction.depth
    factResp = FALSE
    rf.x <- which(names(data)%in%predictors)     
    pred.names = names(data[,rf.x])  #  predictors
    original.response <- object@data@get("response")
    n.preds <- length(pred.names)
    cross.tab <- matrix(0, ncol = n.preds, nrow = n.preds)
    cross.tab2 <- matrix(0, ncol = n.preds, nrow = n.preds)
    dimnames(cross.tab) <- list(pred.names, pred.names)
    dimnames(cross.tab2) <- list(pred.names, pred.names)
      #data <- eval(parse(text = gbm.call$dataframe))[, gbm.x]
    for (i in 1:(n.preds - 1)) {
        if (is.vector(data[, rf.x[i]])) {
            x.var <- seq(min(data[, rf.x[i]], na.rm = T), max(data[,rf.x[i]], na.rm = T), length = 20)
        }else {
            x.var <- factor(names(table(data[, rf.x[i]])), levels = levels(data[,rf.x[i]]))
        }
        x.length <- length(x.var)
        cat(i, " ")
        for (j in (i + 1):n.preds) {
            if (is.vector(data[, rf.x[j]])) {
                y.var <- seq(min(data[, rf.x[j]], na.rm = T), max(data[,rf.x[j]], na.rm = T), length = 20)
            }else {
                y.var <- factor(names(table(data[, rf.x[j]])), levels = levels(data[,rf.x[j]]))
            }
            y.length <- length(y.var)
            pred.frame <- expand.grid(list(x.var, y.var))
            names(pred.frame) <- c(pred.names[i], pred.names[j])

            if(class(data[,rf.x[i]])!=class(pred.frame[,1])) pred.frame[,1] <- eval(parse(text=paste("as.",class(data[,rf.x[i]]),"(pred.frame[,1])",sep="")))
            if(class(data[,rf.x[j]])!=class(pred.frame[,2])) pred.frame[,2] <- eval(parse(text=paste("as.",class(data[,rf.x[j]]),"(pred.frame[,2])",sep="")))
            n <- 3
            for (k in 1:n.preds) {
                if (k != i & k != j) {
                  if (is.vector(data[, rf.x[k]])) {
                    pred.frame[, n] <- mean(data[, rf.x[k]], na.rm = T)
                  }else {
                    temp.table <- sort(table(data[, rf.x[k]]), decreasing = TRUE)
                    pred.frame[, n] <- as.factor(rep(names(temp.table)[1], x.length * y.length))
                    pred.frame[, n] <- factor(pred.frame[, n],levels=levels(data[,rf.x[k]]))
                  }
                  names(pred.frame)[n] <- pred.names[k]
                  if(class(data[,rf.x[k]])!=class(pred.frame[,n])) pred.frame[,n] <- eval(parse(text=paste("as.",class(data[,rf.x[k]]),"(pred.frame[,n])",sep="")))

                  n <- n + 1
                }
            }
            if(!factResp){
              prediction <- predict(object,newdata=pred.frame)
            } 
            if(is.factor(prediction)) factResp=TRUE
            if(factResp){
              prediction <- numeric(length(pred.frame[,1])) #    # 
              for(a in 1:length(pred.frame[,1])){
                prediction[a] <- as.numeric(predict(object,newdata=pred.frame[a,],type="prob")[[1]][,2])
              }
            }
            if(logit) prediction <- qlogis(prediction)
            #prediction <- predict(object,newdata=pred.frame) 
            interaction.test.model <- lm(prediction ~ as.factor(pred.frame[,1]) + as.factor(pred.frame[, 2])) 
             #interaction.flag <- round(mean(resid(interaction.test.model)^2) * 5000, 3)
            interaction.flag <- round(sqrt(mean(resid(interaction.test.model)^2)) * 100, 4)
            sumInt <- sum(prediction)
            sumNoInt <- sum(predict(interaction.test.model))
                   # relative entropy(P||Q), is a measure of the information lost when Q is used to approximate P
                   # Q is additive, P is full forest model
            # browser()
		interaction.flag2 <- sum((prediction/sumInt)*log2((prediction/sumInt)/(predict(interaction.test.model)/sumNoInt)))

            #cat(paste(interaction.flag2,"\n"))
            
            cross.tab[i, j] <- interaction.flag
            cross.tab2[i, j] <- interaction.flag2
        }
    }

            # rank in terms of importance
    search.index <- ((n.preds^2) + 1) - rank(cross.tab, ties.method = "first")
    search.index2 <- ((n.preds^2) + 1) - rank(cross.tab2, ties.method = "first")
    n.important <- max(2, round(0.1 * ((n.preds^2)/2), 0))
    var1.names <- rep(" ", n.important)
    var1.index <- rep(0, n.important)
    var2.names <- rep(" ", n.important)
    var2.index <- rep(0, n.important)
    int.size <- rep(0, n.important)
    var1.names2 <- rep(" ", n.important)
    var1.index2 <- rep(0, n.important)
    var2.names2 <- rep(" ", n.important)
    var2.index2 <- rep(0, n.important)
    int.size2 <- rep(0, n.important)
    for (i in 1:n.important) {
        index.match <- match(i, search.index)
        index.match2 <- match(i, search.index2)
        j <- trunc(index.match/n.preds) + 1
	  j2 <- trunc(index.match2/n.preds) + 1
        var1.index[i] <- j
        var1.index2[i] <- j2
        var1.names[i] <- pred.names[j]
        var1.names2[i] <- pred.names[j2]
        k <- index.match%%n.preds
        k2 <- index.match2%%n.preds
        if (k > 0) {
            var2.index[i] <- k
            var2.names[i] <- pred.names[k]
            int.size[i] <- cross.tab[k, j]
        }
        if (k2 > 0) {
            var2.index2[i] <- k2
            var2.names2[i] <- pred.names[k2]
            int.size2[i] <- cross.tab2[k2, j2]
        }

    }
    rank.list <- data.frame(var1.index, var1.names, var2.index, 
        var2.names, int.size)
    rank.list2 <- data.frame(var1.index2, var1.names2, var2.index2, 
        var2.names2, int.size2)
      #cat("\n")
      #return(list(rank.list = rank.list, interactions = cross.tab, gbm.call = gbm.object$gbm.call))

    int_object <- list(rank.list1=rank.list,rank.list2=rank.list2, interactions1 = cross.tab,interactions2 = cross.tab2)
    return(int_object)

}

##############################
##############################

#######################################################
###   NEW FUNCTION: "RF_PlotInteractions2" 
#                        Plots the residuals from "additive" model"


RF_PlotInteractions2 <- function(yvar=1, xvar=6, object,data,predictors){


    theta = 40
    phi = 25

         # gbm.object <- eval(parse(text=workingModel))
         # RFObject <- data.cforest
         # n.trees <- gbm.call$best.trees
         # depth <- gbm.call$interaction.depth
    factResp = FALSE
    rf.x <- which(names(data)%in%predictors)     
    pred.names = names(data[,rf.x])  #  predictors
    original.response <- object@data@get("response")
    n.preds <- length(pred.names)
    #cross.tab <- matrix(0, ncol = n.preds, nrow = n.preds)
    #cross.tab2 <- matrix(0, ncol = n.preds, nrow = n.preds)
    #dimnames(cross.tab) <- list(pred.names, pred.names)
    #dimnames(cross.tab2) <- list(pred.names, pred.names)
      #data <- eval(parse(text = gbm.call$dataframe))[, gbm.x]
    for (i in yvar) {
        if (is.vector(data[, rf.x[i]])) {
            y.var <- seq(min(data[, rf.x[i]], na.rm = T), max(data[,rf.x[i]], na.rm = T), length = 50)
        }else {
            y.var <- factor(names(table(data[, rf.x[i]])), levels = levels(data[,rf.x[i]]))
        }
        y.length <- length(y.var)
        for (j in xvar) {
            if (is.vector(data[, rf.x[j]])) {
                x.var <- seq(min(data[, rf.x[j]], na.rm = T), max(data[,rf.x[j]], na.rm = T), length = 50)
            }else {
                x.var <- factor(names(table(data[, rf.x[j]])), levels = levels(data[,rf.x[j]]))
            }
            x.length <- length(x.var)
            pred.frame <- expand.grid(list(y.var, x.var))
            names(pred.frame) <- c(pred.names[i], pred.names[j])

            if(class(data[,rf.x[i]])!=class(pred.frame[,1])) pred.frame[,1] <- eval(parse(text=paste("as.",class(data[,rf.x[i]]),"(pred.frame[,1])",sep="")))
            if(class(data[,rf.x[j]])!=class(pred.frame[,2])) pred.frame[,2] <- eval(parse(text=paste("as.",class(data[,rf.x[j]]),"(pred.frame[,2])",sep="")))
            n <- 3
            for (k in 1:n.preds) {
                if (k != i & k != j) {
                  if (is.vector(data[, rf.x[k]])) {
                    pred.frame[, n] <- mean(data[, rf.x[k]], na.rm = T)
                  }else {
                    temp.table <- sort(table(data[, rf.x[k]]), decreasing = TRUE)
                    pred.frame[, n] <- as.factor(rep(names(temp.table)[1], x.length * y.length))
                    pred.frame[, n] <- factor(pred.frame[, n],levels=levels(data[,rf.x[k]]))
                  }
                  names(pred.frame)[n] <- pred.names[k]
                  if(class(data[,rf.x[k]])!=class(pred.frame[,n])) pred.frame[,n] <- eval(parse(text=paste("as.",class(data[,rf.x[k]]),"(pred.frame[,n])",sep="")))

                  n <- n + 1
                }
            }
            if(!factResp){
              prediction <- predict(object,newdata=pred.frame)
            } 
            if(is.factor(prediction)) factResp=TRUE
            if(factResp){
              prediction <- numeric(length(pred.frame[,1])) #    # 
              for(a in 1:length(pred.frame[,1])){
                prediction[a] <- as.numeric(predict(object,newdata=pred.frame[a,],type="prob")[[1]][,2])
              }
            }

               #prediction <- predict(object,newdata=pred.frame) 
            interaction.test.model <- lm(prediction ~ as.factor(pred.frame[,1]) + as.factor(pred.frame[, 2]))
               #interaction.flag <- round(mean(resid(interaction.test.model)^2) * 5000, 3)
               #interaction.flag <- round(sqrt(mean(resid(interaction.test.model)^2)) * 100, 4)
            resids <- resid(interaction.test.model)

            resid.matrix <- matrix(resids, ncol = 50, nrow = 50)
                     # matrix(pred.frame)
          
            ylab= pred.names[i]
            xlab= pred.names[j]
            zlab="residual"

            zlim=c(min(resids),max(resids)) 
              
                
            persp(x = x.var, y = y.var, z = resid.matrix, zlim = zlim, 
                xlab = xlab, ylab = ylab, zlab = zlab, 
                theta = theta, phi = phi, r = sqrt(10), d = 3, 
                ticktype = ticktype, mgp = c(4, 1, 0))

            ##cross.tab[i, j] <- interaction.flag
            ##cross.tab2[i, j] <- interaction.flag2
        }
    }
    diff(range(resid.matrix[,1]))
    diff(range(resid.matrix[,50]))
 }

##############################
##############################








#################################
##################### NEW FUNCTION: PLOT IMPORTANCE FOR GBM

summary.gbm <- function (object, cBars = length(object$var.names), n.trees = object$n.trees, 
    plotit = TRUE, order = TRUE, method = relative.influence, 
    normalize = TRUE, ...) 
{
    if (n.trees < 1) {
        stop("n.trees must be greater than 0.")
    }
    if (n.trees > object$n.trees) {
        warning("Exceeded total number of GBM terms. Results use n.trees=", 
            object$n.trees, " terms.\n")
        n.trees <- object$n.trees
    }
    rel.inf <- method(object, n.trees)
    rel.inf[rel.inf < 0] <- 0
    if (order) {
        i <- order(-rel.inf)
    } else {
        i <- 1:length(rel.inf)
    }
    if (cBars == 0) 
        cBars <- min(10, length(object$var.names))
    if (cBars > length(object$var.names)) 
        cBars <- length(object$var.names)
    if (normalize) 
        rel.inf <- 100 * rel.inf/sum(rel.inf)
    if (plotit) {
        barplot(rel.inf[i[cBars:1]], horiz = TRUE, col = rainbow(cBars, 
            start = 3/6, end = 4/6), names = predictorNames[match(object$var.names,predictors)][i[cBars:1]], 
            xlab = "Relative influence", ...)
    }
    return(data.frame(var = object$var.names[i], rel.inf = rel.inf[i]))
}

#################################################
######################  NEW FUNCTION... PLOT INTERACTIONS FOR RANDOM FOREST....

# rewrite the perspective plotting function

RF_InteractionPlots <- function(x=2,y=5,object,data,predictors,family,zlim=NULL){
    factResp=FALSE
    pred.means = NULL
    x.label = NULL 
    x.range = NULL 
    y.label = NULL
    z.label = "fitted value" 
    y.range = NULL
    z.range = NULL
    leg.coords = NULL
    ticktype = "detailed" 
    theta = 55
    phi = 40
    smooth = "none"
    mask = FALSE
    perspective = TRUE 

    rf.x <- which(names(data)%in%predictors)     
    pred.names = names(data[,rf.x])  #  predictors
    original.response <- object@data@get("response")
    n.preds <- length(pred.names)

    rf.y <- which(names(data)==names(original.response))
    have.factor <- FALSE
    x.name <- pred.names[x]
    if (is.null(x.label)) {
        x.label <- x.name
    }
    y.name <- pred.names[y]
    if (is.null(y.label)) {
        y.label <- y.name
    }
    data2 <- data[,rf.x]    #subset out just the predictor variables...
      #n.trees <- gbm.call$best.trees
    if (is.vector(data2[, x])) {
        if (is.null(x.range)) {
            x.var <- seq(min(data2[,x], na.rm = T), max(data2[,x], na.rm = T), length = 50)
        } else {
            x.var <- seq(x.range[1], x.range[2], length = 50)
        }
    } else {
        x.var <- names(table(data2[, x]))
        have.factor <- TRUE
    }
    if (is.vector(data2[, y])) {
        if (is.null(y.range)) {
            y.var <- seq(min(data2[,y], na.rm = T), max(data2[,y], na.rm = T), length = 50)
        } else {
            y.var <- seq(y.range[1], y.range[2], length = 50)
        }
    } else {
        y.var <- names(table(data2[, y]))
        #if (have.factor) {
        #    stop("at least one marginal predictor must be a vector!")
        #}
        #else {
            have.factor <- TRUE
        #}
    }
    pred.frame <- expand.grid(list(x.var, y.var))
    names(pred.frame) <- c(x.name, y.name)
    pred.rows <- nrow(pred.frame)

    #if (have.factor) {
    #    if (is.factor(pred.frame[, 2])) {
    #        pred.frame <- pred.frame[, c(2, 1)]
    #        x.var <- y.var
    #    }
    #}
    j <- 3
    for (i in 1:n.preds) {
        if (i != x & i != y) {
            if (is.vector(data2[, i])) {
                m <- match(pred.names[i], names(pred.means))
                if (is.na(m)) {
                  pred.frame[, j] <- mean(data2[, i], na.rm = T)
                } else pred.frame[, j] <- pred.means[m]
            }
            if (is.factor(data2[, i])) {
                m <- match(pred.names[i], names(pred.means))
                temp.table <- table(data2[, i])
                if (is.na(m)) {
                  pred.frame[, j] <- rep(names(temp.table)[2],pred.rows)
                }else {
                  pred.frame[, j] <- pred.means[m]
                }
                pred.frame[, j] <- factor(pred.frame[, j], levels = names(temp.table))
            }
            names(pred.frame)[j] <- pred.names[i]
            if(class(data2[, i])!=class(pred.frame[,j])) pred.frame[,j] <- eval(parse(text=paste("as.",class(data2[, i]),"(pred.frame[,j])",sep="")))
            j <- j + 1
        } else{
           if(i==x) if(class(data2[, i])!=class(pred.frame[,1])) pred.frame[,1] <- eval(parse(text=paste("as.",class(data2[, i]),"(pred.frame[,1])",sep="")))
           if(i==y) if(class(data2[, i])!=class(pred.frame[,2])) pred.frame[,2] <- eval(parse(text=paste("as.",class(data2[, i]),"(pred.frame[,2])",sep="")))
        }
    }
    if(!factResp){
      prediction <- predict(object,newdata=pred.frame)
    } 
    if(is.factor(prediction)) factResp=TRUE
    if(factResp){
      prediction <- numeric(length(pred.frame[,1])) #    # 
      for(a in 1:length(pred.frame[,1])){
        prediction[a] <- as.numeric(predict(object,newdata=pred.frame[a,],type="prob")[[1]][,2])
      }
    }

     #prediction <- predict(object, newdata=pred.frame)
    #if (smooth == "model") {
    #    pred.glm <- glm(prediction ~ ns(pred.frame[, 1], df = 8) * 
    #        ns(pred.frame[, 2], df = 8), data = pred.frame, family = poisson)
    #    prediction <- fitted(pred.glm)
    #}
    max.pred <- max(prediction)
    min.pred <- min(prediction)
    cat("maximum value = ", round(max.pred, 2), "\n")
    if (is.null(z.range)) {
        if (family == "bernoulli") {
            if(is.null(zlim)){
		 	z.range <- c(0, 1)
		} else{
      		z.range <- zlim
		}
        } else if (family == "poisson") {
            if(is.null(zlim)){
			z.range <- c(0, max.pred * 1.1)
		} else{
			z.range <- zlim
		}
        } else {
            #z.min <- min(data[, y], na.rm = T)
            #z.max <- max(data[, y], na.rm = T)
            #z.delta <- z.max - z.min
            #z.range <- c(z.min - (1.1 * z.delta), z.max + (1.1 * 
            #    z.delta))
            if(is.null(zlim)){
	 		z.range <- c(min.pred,max.pred)
		}else{
			z.range <- zlim
		}
        }
    }
    if (have.factor == FALSE) {
        pred.matrix <- matrix(prediction, ncol = 50, nrow = 50)
        #if (smooth == "average") {
        #    pred.matrix.smooth <- pred.matrix
        #    for (i in 2:49) {
        #        for (j in 2:49) {
        #          pred.matrix.smooth[i, j] <- mean(pred.matrix[c((i - 1):(i + 1)), c((j - 1):(j + 1))])
        #        }
        #    }
        #    pred.matrix <- pred.matrix.smooth
        #}
        #if (mask) {
        #    mask.trees <- gbm.object$gbm.call$best.trees
        #    point.prob <- predict.gbm(gbm.object[[1]], pred.frame, 
        #        n.trees = mask.trees, type = "response")
        #    point.prob <- matrix(point.prob, ncol = 50, nrow = 50)
        #    pred.matrix[point.prob < 0.5] <- 0
        #}
        if (!perspective) {
            image(x = x.var, y = y.var, z = pred.matrix, zlim = z.range)
        } else {
            persp(x = x.var, y = y.var, z = pred.matrix, zlim = z.range, 
                xlab = x.label, ylab = y.label, zlab = z.label, 
                theta = theta, phi = phi, r = sqrt(10), d = 3, 
                ticktype = ticktype, mgp = c(4, 1, 0))
        }
    }
    if (have.factor) {
        factor.list <- names(table(pred.frame[, 1]))
        n <- 1
        if (is.null(z.range)) {
            vert.limits <- c(min.pred * 1.1, max.pred * 1.1)
        } else {
            vert.limits <- z.range
        }
        plot(as.numeric(pred.frame[pred.frame[, 1] == factor.list[1], 2]), 
            prediction[pred.frame[, 1] == factor.list[1]], type = "n", pch=NA,
            ylim = vert.limits, xlab = y.label, ylab = z.label,xaxt="n")
        if(is.factor(pred.frame[, 2])){ 
          xaxlab <- names(table(pred.frame[, 2]))
          horiz.limits <- c(min(as.numeric(pred.frame[, 2])),max(as.numeric(pred.frame[, 2])))
        }
        if(is.vector(pred.frame[, 2])){
          xaxlab <- round(seq(min(pred.frame[, 2]),max(pred.frame[, 2]),length=5),2)
          horiz.limits <- c(min(pred.frame[, 2]),max(pred.frame[, 2]))
        }
        axis(1,at=c(1:length(xaxlab)),labels=xaxlab)
        for (i in 1:length(factor.list)) {
            factor.level <- factor.list[i]
            #if()
            lines(pred.frame[pred.frame[, 1] == factor.level, 
                2], prediction[pred.frame[, 1] == factor.level], 
                lty = i)
        }
        if (is.null(leg.coords)) {
            #x.max <- max(pred.frame[, 2])
            #x.min <- min(pred.frame[, 2])
            #x.range <- x.max - x.min
            #x.pos <- c(x.min + (0.02 * x.range), x.min + (0.3 * 
            #    x.range))
            #y.max <- max(prediction)
            #y.min <- min(prediction)
            #y.range <- y.max - y.min
            #y.pos <- c(y.min + (0.8 * y.range), y.min + (0.95 * 
            #    y.range))
            legend(horiz.limits[1],vert.limits[2], factor.list, lty = c(1:length(factor.list)), adj=1,   #x = x.pos, y = y.pos, 
                bty = "n")
        }
        else {
            legend(x = leg.coords[1], y = leg.coords[2], factor.list,     #  
                lty = c(1:length(factor.list)), bty = "n")
        }
    }

    #### quantify the effect of the interactions (comment this out later!)
    #names(pred.frame)
    #ndx1 <- which(pred.frame$log.tot.patch.area.10-5<0)
    #ndx2 <- which(pred.frame$log.tot.patch.area.10-13>0)
    #a <- prediction[ndx1]
    #b <- prediction[ndx2]
    #diff(range(a))
    #diff(range(b))

}


cforest_crossValidate <- function(data,full.model,foldVector,
                          pred.names, response, formula=formula1,
                          threshold=NULL,binaryresponse=FALSE,fact=FALSE){

	#####################
	# START CV FUNCTION

	results=list()
	if(binaryresponse) par(ask=TRUE)
	counter = 1
	CVprediction <- numeric(nrow(data))
	CVobserved <- numeric(nrow(data))
	realprediction <- numeric(nrow(data))
	realdata <- numeric(nrow(data))

	predictCols <- which(names(data)%in%pred.names)

	responseData <- eval(parse(text=sprintf("data$%s",response)))

	data.controls = cforestControl
	counter=1
	
	if(fact) ndx = 2 else ndx=1

	i=1
	for(i in 1:n.folds){
	  model <- cforest(formula, data = data[which(foldVector!=i),], controls=data.controls) 
	  predict_CV  <- predict(model,newdata=data[which(foldVector==i),],type="prob") 
	  predict_real  <-  predict(full.model,newdata=data[which(foldVector==i),],type="prob")
	  REAL <- eval(parse(text=sprintf("data$%s[which(foldVector==i)]",response)))
	  j=1
	  for(j in 1:length(which(foldVector==i))){
		CVprediction[counter] <- as.numeric(predict_CV[[j]][,ndx])    ### KTS: check this: was [,2 before]
		CVobserved[counter] <-  REAL[j]      
		realprediction[counter] <- as.numeric(predict_real[[j]][,ndx])   
		realdata[counter] <- REAL[j]         
		counter = counter + 1  
	  }
	}

	if(fact){
	  CVobserved = CVobserved-1
	  realdata=realdata-1
	}

	results$CV_RMSE = sqrt(mean((CVobserved-CVprediction)^2))       # root mean squared error for holdout samples in 10-fold cross-validation ...
	results$real_RMSE = sqrt(mean((CVobserved-realprediction)^2))  # root mean squared error for residuals from final model

	# print RMSE statistics
	cat(sprintf("RMSE for C-V (out of bag) data: %s\n",results$CV_RMSE)) 
	cat(sprintf("RMSE for training (in bag) data: %s\n",results$real_RMSE))  


	results$CV_auc <- NULL
	results$real_auc <- NULL

	if(binaryresponse){
	  #graphics.off()
	  par(mfrow=c(2,1))
	  pred <- prediction(CVprediction,CVobserved)     # for holdout samples in cross-validation
	  perf <- performance(pred,"tpr","fpr")
	  auc <- performance(pred,"auc")
	  plot(perf,main="Cross-validation")
	  results$CV_auc <- round(auc@y.values[[1]],2)
	  text(.9,.1,paste("AUC = ",results$CV_auc,sep=""))
	  
	  pred <- prediction(realprediction,CVobserved)     # for final model
	  perf <- performance(pred,"tpr","fpr")
	  auc <- performance(pred,"auc")
	  plot(perf,main="Training")
	  results$real_auc <- round(auc@y.values[[1]],2)
	  text(.9,.1,paste("AUC = ",results$real_auc,sep=""))
	}

	# COHEN KAPPA statistics

	results$CV_maxkappa <- NULL
	results$real_maxkappa <- NULL

	if(binaryresponse){
		#graphics.off()
		par(mfrow=c(2,1))
		thresholds <- seq(0.01,0.99,length=101)   # "artificial" thresholds across which to examine performance
		kappa <- numeric(length(thresholds))
		for(i in 1:length(thresholds)){
		  trueLabels <- CVobserved
		  predLabels <- ifelse(CVprediction>=thresholds[i],1,0)
		  tot <- length(CVobserved)
		  tp <- length(which((trueLabels==1)&(predLabels==1)))      
		  tn <- length(which((trueLabels==0)&(predLabels==0)))
		  fp <- length(which((trueLabels==0)&(predLabels==1)))
		  fn <- length(which((trueLabels==1)&(predLabels==0)))
		  pr_agree <- (tp+tn)/tot    # overall agreement, or accuracy
		  pr_agree_rand <- ((tp+fn)/tot)*((tp+fp)/tot)+((fn+tn)/tot)*((fp+tn)/tot)
		  kappa[i] <- (pr_agree-pr_agree_rand)/(1-pr_agree_rand)
		}
		plot(thresholds,kappa,type="l",xlab="Threshold", ylab="Cohen's Kappa", main="Holdout sample performance")


		# find threshold value associated with highest Kappa for C-V data

		results$CV_maxkappa <- thresholds[which.max(kappa)]  

		cat(sprintf("The most informative threshold (based on Kappa) is: %s\n",results$CV_maxkappa))


		if(binaryresponse){
		  if(is.null(threshold)) threshold=results$CV_maxkappa
		}

		kappa <- numeric(length(thresholds)) 
		for(i in 1:length(thresholds)){
		  trueLabels <- CVobserved
		  predLabels <- ifelse(realprediction>=thresholds[i],1,0)    
		  tot <- length(CVobserved)
		  tp <- length(which((trueLabels==1)&(predLabels==1)))  
		  tn <- length(which((trueLabels==0)&(predLabels==0)))
		  fp <- length(which((trueLabels==0)&(predLabels==1)))
		  fn <- length(which((trueLabels==1)&(predLabels==0)))
		  pr_agree <- (tp+tn)/tot    # overall agreement, or accuracy
		  pr_agree_rand <- ((tp+fn)/tot)*((tp+fp)/tot)+((fn+tn)/tot)*((fp+tn)/tot)
		  kappa[i] <- (pr_agree-pr_agree_rand)/(1-pr_agree_rand)
		}
		plot(thresholds,kappa,type="l",xlab="Threshold", ylab="Cohen's Kappa", main="Performance: full model")

	}
	 
	results$CV_confusion.mat <- NULL 
	results$CV_sensitivity <- NULL
	results$CV_specificity <- NULL
	results$CV_toterror <- NULL
     

	if(binaryresponse){
		cutoff <- results$CV_maxkappa
		### display confusion matrix and kappa for a single threshold
		trueLabels <- CVobserved
		predLabels <- ifelse(CVprediction>=cutoff,1,0)    
		tot <- length(CVobserved)
		tp <- length(which((trueLabels==1)&(predLabels==1)))  
		tn <- length(which((trueLabels==0)&(predLabels==0)))
		fp <- length(which((trueLabels==0)&(predLabels==1)))
		fn <- length(which((trueLabels==1)&(predLabels==0)))
		pr_agree <- (tp+tn)/tot    # overall agreement, or accuracy
		pr_agree_rand <- ((tp+fn)/tot)*((tp+fp)/tot)+((fn+tn)/tot)*((fp+tn)/tot)
		# results$kappa[i] <- (pr_agree-pr_agree_rand)/(1-pr_agree_rand)
		# results$kappa[i]
		results$CV_confusion.mat<-matrix(c(tp,fp,fn,tn),nrow=2,ncol=2)
		colnames(results$CV_confusion.mat) <- c("Model says true","Model says false")
		rownames(results$CV_confusion.mat) <- c("Actually true","Actually false")
		results$CV_sensitivity <- tp/(tp+fn)
		results$CV_specificity <- tn/(tn+fp)
		results$CV_toterror <- (fn+fp)/tot

	}

	if(binaryresponse){
	  CVprediction[which(CVprediction==1)] <- 0.9999
	  CVprediction[which(CVprediction==0)] <- 0.0001
	  realprediction[which(realprediction==1)] <- 0.9999
	  realprediction[which(realprediction==0)] <- 0.0001
	}


	realdata = CVobserved
	fit_deviance_CV <- mean((CVobserved-CVprediction)^2)
	if(binaryresponse) fit_deviance_CV <- mean(-2*(dbinom(CVobserved,1,CVprediction,log=T)-dbinom(realdata,1,realdata,log=T)))
	fit_deviance_real <- mean((CVobserved-realprediction)^2)
	if(binaryresponse) fit_deviance_real <- mean(-2*(dbinom(CVobserved,1,realprediction,log=T)-dbinom(realdata,1,realdata,log=T)))
	null_deviance <- mean((CVobserved-mean(CVobserved))^2)
	if(binaryresponse) null_deviance <- mean(-2*(dbinom(CVobserved,1,mean(CVobserved),log=T)-dbinom(realdata,1,realdata,log=T)))
	results$deviance_explained_CV <- (null_deviance-fit_deviance_CV)/null_deviance   # based on holdout samples
	results$deviance_explained_real <- (null_deviance-fit_deviance_real)/null_deviance   # based on full model...

	cat(sprintf("McFadden's pseudo R-squared for C-V (out of bag) data: %s\n",results$deviance_explained_CV))
	cat(sprintf("McFadden's pseudo R-squared for training (in bag) data: %s\n",results$deviance_explained_real))

	return(results)
}




#############################################################

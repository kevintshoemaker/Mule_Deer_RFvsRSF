
##TODO

# Do I have the most recent version of the veg classes?

###############
## Clear the workspace
###############

rm(list=ls())

###############
## Global vars
###############

SUMMER = TRUE   # FALSE # 

WINTER = !SUMMER


MODELSRUN <- TRUE


##############
# Load packages
###############

library(lme4) #load lme4 package
library(Hmisc)
library(MuMIn)
library(ROCR)
library(rms)
library(RColorBrewer)

#### Read in the RF functions from github

source("RF_Extensions.R")   # change to your script locations


################
# Load functions
################

####################
## Function for visualizing interactions
####################


VisualizeInteraction <- function(data=deer[["summer"]],model=GLMMs[["summer"]],var1="dist_to_water",var2="elevation",type="GLMM"){
  len <- 50
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


################
## Read in Data
################

  #setwd("K:\\public\\Nathan_Jackson\\Methods_Paper\\R Input")

summerdeer=read.csv('summer_main_modeldata_new25May2017.csv')  # read data file into R

summerdeer$altid=factor(summerdeer$altid) # convert catagorical variable to a vector of factor variables
summerdeer$veg_class=factor(summerdeer$veg_class)


winterdeer=read.csv('winter data 6_28_17.csv')  # read data file into R  winter_main_modeldata_new5May2017.csv

winterdeer$altid=factor(winterdeer$altid) # convert catagorical variable to a vector of factor variables
winterdeer$veg_class=factor(winterdeer$veg_class)


##############
# Standardize covariates
##############

summerdeer$stand_dist_to_water=(summerdeer$dist_to_water-mean(summerdeer$dist_to_water))/sd(summerdeer$dist_to_water)
summerdeer$stand_cos_aspect=(summerdeer$cos_aspect-mean(summerdeer$cos_aspect))/sd(summerdeer$cos_aspect)
summerdeer$stand_sin_aspect=(summerdeer$sin_aspect-mean(summerdeer$sin_aspect))/sd(summerdeer$sin_aspect)
summerdeer$stand_elevation=(summerdeer$elevation-mean(summerdeer$elevation))/sd(summerdeer$elevation)
summerdeer$stand_slope=(summerdeer$slope-mean(summerdeer$slope))/sd(summerdeer$slope)


winterdeer$stand_dist_to_water=(winterdeer$dist_to_water-mean(winterdeer$dist_to_water))/sd(winterdeer$dist_to_water)
winterdeer$stand_cos_aspect=(winterdeer$cos_aspect-mean(winterdeer$cos_aspect))/sd(winterdeer$cos_aspect)
winterdeer$stand_sin_aspect=(winterdeer$sin_aspect-mean(winterdeer$sin_aspect))/sd(winterdeer$sin_aspect)
winterdeer$stand_elevation=(winterdeer$elevation-mean(winterdeer$elevation))/sd(winterdeer$elevation)
winterdeer$stand_slope=(winterdeer$slope-mean(winterdeer$slope))/sd(winterdeer$slope)


summerdeer$used_fac=factor(summerdeer$used)

winterdeer$used_fac=factor(winterdeer$used)

deer <- list()

deer[["summer"]] <- summerdeer
deer[["winter"]] <- winterdeer

rm(summerdeer,winterdeer)

###############
# Naming covariables
###############

############### NAMING VARIABLES ############

predictorNames <- c(  "Cos Aspect", # nice readable names
                      "Sin Aspect",
                      "Elevation",
                      "Slope",
                      "Vegetation Class",
                      "Distance to Water"
)

pred.names=c(  "cos_aspect",
               "sin_aspect",
               "elevation",
               "slope",
               "veg_class",
               "dist_to_water"
               
)


#name check
cbind(pred.names,predictorNames)

stand_pred.names <- paste("stand_",pred.names,sep="")



##################
## Run the generalized linear mixed model   [note: kelley's code has * instead of :  -- not sure if it was right!]
##################

#   MODELSRUN <- FALSE
if(MODELSRUN){
  load("summerGLMM.RData")
  load("winterGLMM.RData")
  
  GLMMs <- list()
  GLMMs[["winter"]] <- winterGLMM
  GLMMs[["summer"]] <- summerGLMM
  
  rm(summerGLMM,winterGLMM)
  
}else{
  summerGLMM=glmer(used~stand_dist_to_water + stand_cos_aspect + stand_sin_aspect + 
                              stand_elevation + stand_slope + veg_class + stand_elevation:stand_slope +
                              stand_dist_to_water:stand_slope + stand_dist_to_water:stand_elevation +
                              (1|altid), family="binomial", data=deer[["summer"]],na.action="na.fail") # generalized linear mixed effect
  
  winterGLMM=glmer(used~stand_dist_to_water + stand_cos_aspect + stand_sin_aspect + 
                              stand_elevation + stand_slope + veg_class +  stand_elevation:stand_slope +
                              stand_dist_to_water:stand_slope +
                              (1|altid), family="binomial", data=deer[["winter"]],na.action="na.fail") # generalized linear mixed effect 
                           
    
  GLMMs <- list()
  GLMMs[["winter"]] <- winterGLMM
  GLMMs[["summer"]] <- summerGLMM
  
  # summary(summer_pequop_final)
  save(summerGLMM,file = "summerGLMM.RData")
  save(winterGLMM,file = "winterGLMM.RData")
  
  rm(summerGLMM,winterGLMM)
}



##################
## Run the Random Forest model
##################

#### Define response variable

response="used_fac" 


#### Define our formula (response ~ predictors)

#   MODELSRUN <- FALSE
if(MODELSRUN){
  load("summerRF.RData")
  load("winterRF.RData")
  
  RFs <- list()
  RFs[["winter"]] <- winterRF
  RFs[["summer"]] <- summerRF
  RFs[["winter_importance"]] <- winter_importance
  RFs[["summer_importance"]] <- summer_importance
  
  rm(summerRF,summer_importance,winterRF,winter_importance)
  
}else{
  formula1 <- as.formula(paste(response,"~",paste(pred.names,collapse="+")))
  
  cforestControl <- cforest_unbiased(ntree=500,mtry=3)   # change back to 500!!
  cforestControl@fraction <- 0.03
  
  summerRF <- cforest(formula1, controls=cforestControl, data=deer[["summer"]])
  # get the importance values
  summer_importance<-varimp((summerRF), conditional= FALSE)
  
  
  winterRF <- cforest(formula1, controls=cforestControl, data=deer[["winter"]])
  # get the importance values
  winter_importance<-varimp((winterRF), conditional= FALSE)
  
  RFs <- list()
  RFs[["winter"]] <- winterRF
  RFs[["summer"]] <- summerRF
  RFs[["winter_importance"]] <- winter_importance
  RFs[["summer_importance"]] <- summer_importance
  
  # summary(summer_pequop_final)
  save(summerRF,summer_importance,file = "summerRF.RData")
  save(winterRF,winter_importance,file = "winterRF.RData")
  
  rm(summerRF,summer_importance,winterRF,winter_importance)
}


#   MODELSRUN <- TRUE


#######################
# Visualize importance values
#######################

graphics.off()

svg("importancefig.svg",height=7,width=7)

par(mfcol=c(2,2))
par(mai=c(1,1.5,0.6,0.4))
lengthndx <- length(RFs[[impname]])
col <- rainbow(lengthndx, start = 3/6, end = 4/6)      # rep(brewer.pal(6,"Blues"),each=2)


season <- "summer"
impname <- sprintf("%s_importance",season)

barplot(height=RFs[[impname]][order(RFs[[impname]],decreasing = FALSE)],
        horiz=T,las=1,main="Importance of Predictors, RF",
        xlab="Index of overall importance",col=col,           
        names.arg=predictorNames[match(names(RFs[[impname]]),pred.names)][order(RFs[[impname]],decreasing = FALSE)])


season <- "winter"
impname <- sprintf("%s_importance",season)



barplot(height=RFs[[impname]][order(RFs[[impname]],decreasing = FALSE)],
        horiz=T,las=1,main="Importance of Predictors, RF",
        xlab="Index of overall importance",col=col,           
        names.arg=predictorNames[match(names(RFs[[impname]]),pred.names)][order(RFs[[impname]],decreasing = FALSE)])


##########
# GLM "importance"? (coefficients of standardized variables?)
##########


season="summer"
summ <- summary(GLMMs[[season]])

allvars <- names(summ$coefficients[,1])

# remove interactions
allvars <- allvars[-grep(":",allvars)]

# keep the standardized vars

allvars <- allvars[grep("stand",allvars)]

glm_importance <- summ$coefficients[,1][allvars]

#glm_importance <- glm_importance[order(abs(glm_importance),decreasing = T)]

barplot(height=glm_importance[order(abs(glm_importance),decreasing = FALSE)],
        horiz=T,las=1,main="Standardized coefficients, GLMM",
        xlab="Standardized coefficient",col=col,           
        names.arg=predictorNames[match(names(glm_importance),stand_pred.names)][order(abs(glm_importance),decreasing = FALSE)])




season="winter"
summ <- summary(GLMMs[[season]])

allvars <- names(summ$coefficients[,1])

    # remove interactions
allvars <- allvars[-grep(":",allvars)]

    # keep the standardized vars

allvars <- allvars[grep("stand",allvars)]

glm_importance <- summ$coefficients[,1][allvars]

#glm_importance <- glm_importance[order(abs(glm_importance),decreasing = T)]

barplot(height=glm_importance[order(abs(glm_importance),decreasing = FALSE)],
         horiz=T,las=1,main="Standardized coefficients, GLMM",
         xlab="Standardized coefficient",col=col,           
         names.arg=predictorNames[match(names(glm_importance),stand_pred.names)][order(abs(glm_importance),decreasing = FALSE)])


dev.off()

## note: comparison of RF vs GLMM "importance" values vs regression coefficients of standardized variables for RSF can 
    #  indicate the degree to which nonlinearities and interactions play a role!

## note: RF gives an overall importance of categorical vars, GLMM cannot do that!




#######################
## Visualize univariate relationships
#######################

pred.names

## NOTE: we need to convert to real scale for the axes (unstandardized)

season <- "summer"
VisualizeRelation(data=deer[[season]],model=GLMMs[[season]],predvar="cos_aspect",type="GLMM")
VisualizeRelation(data=deer[[season]],model=GLMMs[[season]],predvar="sin_aspect",type="GLMM")
VisualizeRelation(data=deer[[season]],model=GLMMs[[season]],predvar="elevation",type="GLMM")
VisualizeRelation(data=deer[[season]],model=GLMMs[[season]],predvar="slope",type="GLMM")
VisualizeRelation(data=deer[[season]],model=GLMMs[[season]],predvar="veg_class",type="GLMM")
VisualizeRelation(data=deer[[season]],model=GLMMs[[season]],predvar="dist_to_water",type="GLMM")

season <- "winter"
VisualizeRelation(data=deer[[season]],model=GLMMs[[season]],predvar="cos_aspect",type="GLMM")
VisualizeRelation(data=deer[[season]],model=GLMMs[[season]],predvar="sin_aspect",type="GLMM")
VisualizeRelation(data=deer[[season]],model=GLMMs[[season]],predvar="elevation",type="GLMM")
VisualizeRelation(data=deer[[season]],model=GLMMs[[season]],predvar="slope",type="GLMM")
VisualizeRelation(data=deer[[season]],model=GLMMs[[season]],predvar="veg_class",type="GLMM")
VisualizeRelation(data=deer[[season]],model=GLMMs[[season]],predvar="dist_to_water",type="GLMM")


season="summer"
graphics.off()

svg("summerunivarplots.svg",6,6)
par(mfcol=c(3,2))
par(mai=c(0.8,0.8,0,0))


VisualizeRelation(data=deer[[season]],model=RFs[[season]],predvar="elevation",type="RF")
VisualizeRelation(data=deer[[season]],model=RFs[[season]],predvar="slope",type="RF")
VisualizeRelation(data=deer[[season]],model=RFs[[season]],predvar="dist_to_water",type="RF")


VisualizeRelation(data=deer[[season]],model=GLMMs[[season]],predvar="elevation",type="GLMM")
VisualizeRelation(data=deer[[season]],model=GLMMs[[season]],predvar="slope",type="GLMM")
VisualizeRelation(data=deer[[season]],model=GLMMs[[season]],predvar="dist_to_water",type="GLMM")

dev.off()

graphics.off()
season <- "winter"
svg("winterunivariateplots.svg",6,6)

par(mfcol=c(3,2))
par(mai=c(0.8,0.8,0,0))


VisualizeRelation(data=deer[[season]],model=RFs[[season]],predvar="elevation",type="RF")
VisualizeRelation(data=deer[[season]],model=RFs[[season]],predvar="slope",type="RF")
VisualizeRelation(data=deer[[season]],model=RFs[[season]],predvar="dist_to_water",type="RF")


VisualizeRelation(data=deer[[season]],model=GLMMs[[season]],predvar="elevation",type="GLMM")
VisualizeRelation(data=deer[[season]],model=GLMMs[[season]],predvar="slope",type="GLMM")
VisualizeRelation(data=deer[[season]],model=GLMMs[[season]],predvar="dist_to_water",type="GLMM")

dev.off()


graphics.off()
season="summer"
svg("vegplots.svg",6,6)

par(mfrow=c(2,2))
par(mai=c(1.5,0.8,0,0))

VisualizeRelation(data=deer[[season]],model=RFs[[season]],predvar="veg_class",type="RF")
VisualizeRelation(data=deer[[season]],model=GLMMs[[season]],predvar="veg_class",type="GLMM")



season="winter"
#svg("wintervegplots.svg",6,5)

VisualizeRelation(data=deer[[season]],model=RFs[[season]],predvar="veg_class",type="RF")
VisualizeRelation(data=deer[[season]],model=GLMMs[[season]],predvar="veg_class",type="GLMM")

dev.off()

#######################
## Visualize interactions
########################


graphics.off()

season <- "summer"


VisualizeInteraction(data=deer[[season]],model=GLMMs[[season]],var1="dist_to_water",var2="elevation")

#VisualizeInteraction(data=deer[[season]],model=GLMMs[[season]],var1="stand_elevation",var2="stand_dist_to_water")

#VisualizeInteraction(data=deer[[season]],model=GLMMs[[season]],var1="stand_slope",var2="stand_dist_to_water")

VisualizeInteraction(data=deer[[season]],model=GLMMs[[season]],var2="slope",var1="dist_to_water")

VisualizeInteraction(data=deer[[season]],model=GLMMs[[season]],var1="slope",var2="elevation")



season <- "winter"

VisualizeInteraction(data=deer[[season]],model=GLMMs[[season]],var1="dist_to_water",var2="elevation")

#VisualizeInteraction(data=deer[[season]],model=GLMMs[[season]],var1="stand_elevation",var2="stand_dist_to_water")

#VisualizeInteraction(data=deer[[season]],model=GLMMs[[season]],var1="stand_slope",var2="stand_dist_to_water")

VisualizeInteraction(data=deer[[season]],model=GLMMs[[season]],var2="slope",var1="dist_to_water")

VisualizeInteraction(data=deer[[season]],model=GLMMs[[season]],var1="slope",var2="elevation")

#VisualizeInteraction(data=deer[[season]],model=GLMMs[[season]],var1="stand_slope",var2="stand_elevation")


graphics.off()
svg("disttowater_vs_elevation.svg",6,6)
par(mfrow=c(2,2))

season <- "summer"
VisualizeInteraction(data=deer[[season]],model=RFs[[season]],var1="dist_to_water",var2="elevation",type="RF")
VisualizeInteraction(data=deer[[season]],model=GLMMs[[season]],var1="dist_to_water",var2="elevation",type="GLMM")


season <- "winter"
VisualizeInteraction(data=deer[[season]],model=RFs[[season]],var1="dist_to_water",var2="elevation",type="RF")
VisualizeInteraction(data=deer[[season]],model=GLMMs[[season]],var1="dist_to_water",var2="elevation",type="GLMM")

dev.off()


graphics.off()
svg("slope_vs_elevation.svg",6,6)
par(mfrow=c(2,2))

season <- "summer"
VisualizeInteraction(data=deer[[season]],model=RFs[[season]],var1="slope",var2="elevation",type="RF")
VisualizeInteraction(data=deer[[season]],model=GLMMs[[season]],var1="slope",var2="elevation",type="GLMM")


season <- "winter"
VisualizeInteraction(data=deer[[season]],model=RFs[[season]],var1="slope",var2="elevation",type="RF")
VisualizeInteraction(data=deer[[season]],model=GLMMs[[season]],var1="slope",var2="elevation",type="GLMM")

dev.off()



#################
# MODEL SELECTION: GLMM
#################


if(!MODELSRUN){
  
  #### Summer Global model
  
  summer_pequop_model=glmer(used~stand_dist_to_water +
                              stand_cos_aspect + 
                              stand_sin_aspect +
                              stand_elevation +
                              stand_slope +
                              veg_class +
                              stand_dist_to_water:stand_elevation + 
                              stand_dist_to_water:stand_slope +
                              stand_elevation:stand_slope + 
                              (1|altid),
              family="binomial",data=deer$summer,nAGQ = 1, na.action="na.fail") # generalized linear mixed effect 
  
  
  summer_dred <- dredge(summer_pequop_model, trace = TRUE, rank = "AICc", REML = FALSE)
  summer_dred
  
  model_summary_summer <- summary(summer_pequop_model)
  
  write.csv(model_summary_summer$coefficients, file= "top_summer_main_result19April2017.csv")
  write.csv(summer_dred, file= "summer_main_dredge_result_10April2017.csv")
  
  
  #### Winter Global model
  
  winter_pequop_model=glmer(used~stand_dist_to_water +
                              stand_cos_aspect + 
                              stand_sin_aspect +
                              stand_elevation +
                              stand_slope +
                              veg_class +
                              stand_dist_to_water:stand_elevation + 
                              stand_dist_to_water:stand_slope +
                              stand_elevation:stand_slope + 
                              (1|altid),
                            family="binomial",data=deer$winter,nAGQ = 1, na.action="na.fail") # generalized linear mixed effect 
  
  model_summary_winter <- summary(winter_pequop_model)
  
  winter_dred<- dredge(winter_pequop_model, trace = TRUE, rank = "AICc", REML = FALSE)
  winter_dred
  
  write.csv(model_summary_winter$coefficients, file= "top_winter_main_result19April2017.csv")
  write.csv(winter_dred, file= "winter_main_dredge_result_10April2017.csv")
  

  

}





######################
# RF Interactions
######################

RF_int <- list()

season <- "summer"

# NOTE: this one can take a very long time   ...
RF_int[[season]] <- RF_FindInteractions(object=RFs[[season]],data=deer[[season]],predictors=pred.names)

### plot interaction strength
graphics.off()
lengthndx <- min(9,nrow(RF_int[[season]]$rank.list1))
par(mai=c(0.95,3.1,0.6,0.4))
#ndx <- ndx <- which(predictors%in%pred.names)
barplot(height=(RF_int[[season]]$rank.list1[c(1:min(9,nrow(RF_int[[season]]$rank.list1))),5][c(lengthndx:1)]),
        horiz=T,las=1,main=paste(response, sep=""),
        xlab="Index of interaction strength",col=brewer.pal(lengthndx,"Blues"),           
        names.arg=paste("",predictorNames[match(RF_int[[season]]$rank.list1[,2][c(lengthndx:1)],pred.names)],"\n",predictorNames[match(RF_int[[season]]$rank.list1[,4][c(lengthndx:1)],pred.names)],sep="") )



RF_int[[season]]$rank.list1

season <- "winter"

# NOTE: this one can take a very long time   ...
RF_int[[season]] <- RF_FindInteractions(object=RFs[[season]],data=deer[[season]],predictors=pred.names)

### plot interaction strength
graphics.off()
lengthndx <- min(9,nrow(RF_int[[season]]$rank.list1))
par(mai=c(0.95,3.1,0.6,0.4))
#ndx <- ndx <- which(predictors%in%pred.names)
barplot(height=(RF_int[[season]]$rank.list1[c(1:min(9,nrow(RF_int[[season]]$rank.list1))),5][c(lengthndx:1)]),
        horiz=T,las=1,main=paste(response, sep=""),
        xlab="Index of interaction strength",col=brewer.pal(lengthndx,"Blues"),           
        names.arg=paste("",predictorNames[match(RF_int[[season]]$rank.list1[,2][c(lengthndx:1)],pred.names)],"\n",predictorNames[match(RF_int[[season]]$rank.list1[,4][c(lengthndx:1)],pred.names)],sep="") )



RF_int[[season]]$rank.list1



##### CONDITIONAL INFERENCE TREE  ##################

season="winter"

summer_deer <- ctree(formula=formula1, data=df, controls = ctree_control(mincriterion = 0.99,maxdepth = 4))

plot(summer_deer)

summary(summer_deer)












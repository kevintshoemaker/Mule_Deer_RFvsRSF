
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

spatialdir <- "E:\\Dropbox\\Mule Deer\\Methods paper\\GIS"

##############
# Load packages
###############

library(lme4) #load lme4 package
library(Hmisc)
library(MuMIn)
library(ROCR)
library(rms)
library(RColorBrewer)
library(pROC)    # for running the delong test

library(raster)
library(rgdal)
library(prettymapr)


################
# Load functions
################

#### Read in the RF functions from github

# setwd("E:\\GIT\\Mule_Deer_RFvsRSF")

source("RF_Extensions.R")   # change to your script locations (not necessary if cloned from GITHUB)

source("METHODS_PAPER_ALLFUNCTIONS.R")

################
## Read in Data
################

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
## Run the generalized linear mixed model
##################

MODELSRUN <- TRUE
# setwd("E:\\Dropbox\\Mule Deer\\Methods paper\\CODE_deprecated")
# MODELSRUN <- FALSE
if(MODELSRUN){
  load("summerGLMM.RData")
  load("winterGLMM.RData")
  
  GLMMs <- list()
  GLMMs[["winter"]] <- winterGLMM
  GLMMs[["summer"]] <- summerGLMM
  
  rm(summerGLMM,winterGLMM)
  
}else{
  
  
  ##### GLOBAL MODEL
  
  summer_pequop_model=glmer(used~stand_dist_to_water+(stand_dist_to_water*stand_elevation)+(stand_dist_to_water*stand_slope)
                            +(stand_elevation*stand_slope)+stand_cos_aspect+stand_sin_aspect
                            +stand_elevation+stand_slope+veg_class+(1|altid),family="binomial",data=deer$summer,nAGQ = 1, na.action="na.fail") # generalized linear mixed effect 
  
  summary(summer_pequop_model)
  
       # find the best model using "dredge" function
  
  summer_dred<- dredge(summer_pequop_model, trace = TRUE, rank = "AICc", REML = FALSE)
  summer_dred
  
  write.csv(model_summary$coefficients, file= "top_summer_main_result_cincodemayo.csv")
  write.csv(summer_dred, file= "winter_main_dredge_result_CINCODEMAYO.csv")
  
       # run the best models
  
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

#### Define response variable and formula

response="used_fac" 

formula1 <- as.formula(paste(response,"~",paste(pred.names,collapse="+")))

#### Define the RF settings

cforestControl <- cforest_unbiased(ntree=500,mtry=3)   # change back to 500!!
cforestControl@fraction <- 0.03

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



##### PLOT EXAMPLE CONDITIONAL INFERENCE TREE  ##################

season="winter"

summer_deer <- ctree(formula=formula1, data=df, controls = ctree_control(mincriterion = 0.99,maxdepth = 4))

plot(summer_deer)

summary(summer_deer)



##################
## CROSS VALIDATION
##################

CVsRun <- TRUE

if(CVsRun){
  
  # setwd("E:\\Dropbox\\Mule Deer\\Methods paper\\CODE_deprecated")
  load("CVs.RData")
  
}else{
  CVs <- list()
  
  CVs[["GLMM"]] <- list()
  CVs[["RF"]] <- list()
  
  uniquedeer <- as.character(unique(deer[[season]]$altid))
  n.folds <- length(uniquedeer)
  
  
  type= "GLMM" #"RF"
  season="summer"
  CVs[[type]][[season]] <- CrossValidateByDeer(n.fold=n.folds,season=season,type=type,plot=F)
  
  type= "GLMM" #"RF"
  season="winter"
  CVs[[type]][[season]] <- CrossValidateByDeer(n.fold=n.folds,season=season,type=type,plot=F)
  
  type= "RF" #"RF"
  season="summer"
  CVs[[type]][[season]] <- CrossValidateByDeer(n.fold=n.folds,season=season,type=type,plot=F)
  
  type= "RF" #"RF"
  season="winter"
  CVs[[type]][[season]] <- CrossValidateByDeer(n.fold=n.folds,season=season,type=type,plot=F)
  
  ## Save CV objects
  # setwd("E:\\Dropbox\\Mule Deer\\Methods paper\\CODE_deprecated")
  
  save(CVs,file = "CVs.RData")
  
}

type= "GLMM" #"RF"
season="summer"
PlotPerformance(CVs[[type]][[season]])

type= "GLMM" #"RF"
season="winter"
PlotPerformance(CVs[[type]][[season]])

type= "RF" #"RF"
season="summer"
PlotPerformance(CVs[[type]][[season]])

type= "RF" #"RF"
season="winter"
PlotPerformance(CVs[[type]][[season]])

######################
# run the Delong test
######################

?roc.test

### summer
season="summer"
roc1 <- roc(CVs$GLMM[[season]]$observed,CVs$GLMM[[season]]$CVpred)
roc2 <- roc(CVs$RF[[season]]$observed,CVs$RF[[season]]$CVpred)
roc.test(roc1, roc2)

### winter
season="winter"
roc1 <- roc(CVs$GLMM[[season]]$observed,CVs$GLMM[[season]]$CVpred)
roc2 <- roc(CVs$RF[[season]]$observed,CVs$RF[[season]]$CVpred)
roc.test(roc1, roc2)


# print RMSE statistics
CV_RMSE 
write.csv(CV_RMSE,file="CV_RMSE_10percent.csv")
real_RMSE
write.csv(real_RMSE,file="real_RMSE_10percent.csv")




###################
## DO SPATIAL PROJECTIONS
###################

setwd(spatialdir)

rastercovs <- list()

rastercovs[["winter"]] <- list()
rastercovs[["summer"]] <- list()

##########
# READ IN GIS DATA
##########

## summer cosine
setwd("E:\\Dropbox\\Mule Deer\\Methods paper\\GIS\\Rasters\\summer_cos")
rastercovs[["summer"]][["cos_aspect"]] <- raster("w001001.adf") 
   #rastercovs[["summer"]][["cos_aspect"]] <- trim(rastercovs[["summer"]][["cos_aspect"]]) # note: trim function takes too long to run here..
#plot(rastercovs[["summer"]][["cos_aspect"]])

old_summer_extent <- extent(rastercovs[["summer"]][["cos_aspect"]])
summer_extent <- extent(old_summer_extent@xmin,680000,4580000,4635000)
rastercovs[["summer"]][["cos_aspect"]] <- crop(rastercovs[["summer"]][["cos_aspect"]],summer_extent)
#plot(rastercovs[["summer"]][["cos_aspect"]])

rastercovs[["summer"]][["cos_aspect"]] <- trim(rastercovs[["summer"]][["cos_aspect"]])
summer_extent <- extent(rastercovs[["summer"]][["cos_aspect"]])   # final extent!
rastercovs[["summer"]][["cos_aspect"]] <- crop(rastercovs[["summer"]][["cos_aspect"]],summer_extent)
plot(rastercovs[["summer"]][["cos_aspect"]])

rastercovs[["summer"]][["cos_aspect"]] <- aggregate(rastercovs[["summer"]][["cos_aspect"]], fact=4, fun=modal)
plot(rastercovs[["summer"]][["cos_aspect"]])

summer_crs <- rastercovs[["summer"]][["cos_aspect"]]@crs 

## summer sine
setwd("E:\\Dropbox\\Mule Deer\\Methods paper\\GIS\\Rasters\\summer_sin")
rastercovs[["summer"]][["sin_aspect"]] <- raster("w001001.adf")
rastercovs[["summer"]][["sin_aspect"]] <- crop(rastercovs[["summer"]][["sin_aspect"]],summer_extent)
rastercovs[["summer"]][["sin_aspect"]] <- aggregate(rastercovs[["summer"]][["sin_aspect"]], fact=4, fun=modal)
plot(rastercovs[["summer"]][["sin_aspect"]])

## summer elevation
setwd("E:\\Dropbox\\Mule Deer\\Methods paper\\GIS\\Rasters\\summer_elev")
rastercovs[["summer"]][["elevation"]] <- raster("w001001.adf")
rastercovs[["summer"]][["elevation"]] <- crop(rastercovs[["summer"]][["elevation"]],summer_extent)
rastercovs[["summer"]][["elevation"]] <- aggregate(rastercovs[["summer"]][["elevation"]], fact=4, fun=modal)
plot(rastercovs[["summer"]][["elevation"]])

## summer slope
setwd("E:\\Dropbox\\Mule Deer\\Methods paper\\GIS\\Rasters\\summer_slope")
rastercovs[["summer"]][["slope"]] <- raster("w001001.adf")
rastercovs[["summer"]][["slope"]] <- crop(rastercovs[["summer"]][["slope"]],summer_extent)
rastercovs[["summer"]][["slope"]] <- aggregate(rastercovs[["summer"]][["slope"]], fact=4, fun=modal)
plot(rastercovs[["summer"]][["slope"]])

## summer dist to water
setwd(spatialdir)
rastercovs[["summer"]][["dist_to_water"]] <- raster("Summer_H20.tif")
rastercovs[["summer"]][["dist_to_water"]] <- crop(rastercovs[["summer"]][["dist_to_water"]],summer_extent)
rastercovs[["summer"]][["dist_to_water"]] <- aggregate(rastercovs[["summer"]][["dist_to_water"]], fact=4, fun=modal)
#extent(rastercovs[["summer"]][["dist_to_water"]]) <- summer_extent 
rastercovs[["summer"]][["dist_to_water"]] <- projectRaster(rastercovs[["summer"]][["dist_to_water"]],rastercovs[["summer"]][["slope"]])
plot(rastercovs[["summer"]][["dist_to_water"]])

## summer veg
setwd("E:\\Dropbox\\Mule Deer\\Methods paper\\GIS\\Summer_veg_raster\\Summer_veg_raster")
rastercovs[["summer"]][["veg_class"]] <- raster("summer_veg.tif")
rastercovs[["summer"]][["veg_class"]] <- crop(rastercovs[["summer"]][["veg_class"]],summer_extent)
attrib <- rastercovs[["summer"]][["veg_class"]]@data@attributes
target <- levels(deer$summer$veg_class)
subs <- attrib[[1]][,c("Value","Habitat_Type")]
subs$HabType2 <- c("Riparian","PJ","mtn_mahog","asage","intro_grass","grassland","shb_mead","roads","desert_scrub","Dec_shrub","aspen","NA")
habitats <- rastercovs[["summer"]][["veg_class"]]@data@attributes[[1]]$Habitat_Type
is.factor(rastercovs[["summer"]][["veg_class"]])
rastercovs[["summer"]][["veg_class"]] <- aggregate(rastercovs[["summer"]][["veg_class"]], fact=4, fun=modal)
rastercovs[["summer"]][["veg_class"]] <- projectRaster(rastercovs[["summer"]][["veg_class"]],rastercovs[["summer"]][["slope"]],method="ngb")
rastercovs[["summer"]][["veg_class"]] <- ratify(rastercovs[["summer"]][["veg_class"]])
rat <- levels(rastercovs[["summer"]][["veg_class"]])[[1]]
rat$landcover <- subs$HabType2
rat$code <- subs$Habitat_Type
levels(rastercovs[["summer"]][["veg_class"]]) <- rat

#rastercovs[["summer"]][["veg_class"]]@data@attributes <- attrib
unique(rastercovs$summer$veg_class@data@values)

# rastercovs[["summer"]][["veg_class"]] <- reclassify(rastercovs[["summer"]][["veg_class"]],rcl=c(46,48,NA)) 
# rastercovs[["summer"]][["veg_class"]] <- subs(rastercovs[["summer"]][["veg_class"]],subs[,c("Value","HabType2")],by=1,which=2)
# rastercovs[["summer"]][["veg_class"]] <- as.factor(rastercovs[["summer"]][["veg_class"]])
plot(rastercovs[["summer"]][["veg_class"]])

res(rastercovs[["summer"]][["veg_class"]])

#########
# WINTER!
#########

## winter cosine
setwd("E:\\Dropbox\\Mule Deer\\Methods paper\\GIS\\Rasters")
rastercovs[["winter"]][["cos_aspect"]] <- raster("Winter_Cos2.tif") 
#rastercovs[["winter"]][["cos_aspect"]] <- trim(rastercovs[["winter"]][["cos_aspect"]]) # note: trim function takes too long to run here..
plot(rastercovs[["winter"]][["cos_aspect"]])

# dpath<-"E:\\Dropbox\\Mule Deer\\Methods paper\\GIS\\Rasters\\winter_cos\\w001001.adf" 
# x <- new("GDALReadOnlyDataset", dpath) 
# getDriver(x) 
# getDriverLongName(getDriver(x)) 

rastercovs[["winter"]][["cos_aspect"]] <- trim(rastercovs[["winter"]][["cos_aspect"]])
winter_extent <- extent(rastercovs[["winter"]][["cos_aspect"]])   # final extent!
rastercovs[["winter"]][["cos_aspect"]] <- crop(rastercovs[["winter"]][["cos_aspect"]],winter_extent)
plot(rastercovs[["winter"]][["cos_aspect"]])

rastercovs[["winter"]][["cos_aspect"]] <- aggregate(rastercovs[["winter"]][["cos_aspect"]], fact=4, fun=modal)
plot(rastercovs[["winter"]][["cos_aspect"]])

## winter sine  ## DOESN'T WORK
setwd("E:\\Dropbox\\Mule Deer\\Methods paper\\GIS\\Rasters")
rastercovs[["winter"]][["sin_aspect"]] <- raster("Winter_Sin2.tif")
rastercovs[["winter"]][["sin_aspect"]] <- crop(rastercovs[["winter"]][["sin_aspect"]],winter_extent)
rastercovs[["winter"]][["sin_aspect"]] <- aggregate(rastercovs[["winter"]][["sin_aspect"]], fact=4, fun=modal)
plot(rastercovs[["winter"]][["sin_aspect"]])

## winter elevation
setwd("E:\\Dropbox\\Mule Deer\\Methods paper\\GIS\\Rasters\\winter_elev")
rastercovs[["winter"]][["elevation"]] <- raster("w001001.adf")
rastercovs[["winter"]][["elevation"]] <- crop(rastercovs[["winter"]][["elevation"]],winter_extent)
rastercovs[["winter"]][["elevation"]] <- aggregate(rastercovs[["winter"]][["elevation"]], fact=4, fun=modal)
rastercovs[["winter"]][["elevation"]] <- projectRaster(rastercovs[["winter"]][["elevation"]],rastercovs[["winter"]][["sin_aspect"]])
plot(rastercovs[["winter"]][["elevation"]])

## winter slope 
setwd("E:\\Dropbox\\Mule Deer\\Methods paper\\GIS\\Rasters")
rastercovs[["winter"]][["slope"]] <- raster("Winter_Slope2.tif")
rastercovs[["winter"]][["slope"]] <- crop(rastercovs[["winter"]][["slope"]],winter_extent)
rastercovs[["winter"]][["slope"]] <- aggregate(rastercovs[["winter"]][["slope"]], fact=4, fun=modal)
plot(rastercovs[["winter"]][["slope"]])

## winter dist to water
setwd(spatialdir)
rastercovs[["winter"]][["dist_to_water"]] <- raster("Winter_H2O.tif")
rastercovs[["winter"]][["dist_to_water"]] <- crop(rastercovs[["winter"]][["dist_to_water"]],winter_extent)
rastercovs[["winter"]][["dist_to_water"]] <- aggregate(rastercovs[["winter"]][["dist_to_water"]], fact=4, fun=modal)
plot(rastercovs[["winter"]][["dist_to_water"]])

## winter veg
setwd("E:\\Dropbox\\Mule Deer\\Methods paper\\GIS\\Winter_veg_raster\\Winter_veg_raster")
rastercovs[["winter"]][["veg_class"]] <- raster("winter_veg.tif")
rastercovs[["winter"]][["veg_class"]] <- crop(rastercovs[["winter"]][["veg_class"]],winter_extent)
attrib <- rastercovs[["winter"]][["veg_class"]]@data@attributes
target <- levels(deer$winter$veg_class)
subs <- attrib[[1]][,c("ID","HABITAT_TYPE")]
subs$HabType2 <- c("Riparian","Fir_forest","PJ","mtn_mahog","asage","grssland")
habitats <- rastercovs[["winter"]][["veg_class"]]@data@attributes[[1]]$HABITAT_TYPE
is.factor(rastercovs[["winter"]][["veg_class"]])
rastercovs[["winter"]][["veg_class"]] <- aggregate(rastercovs[["winter"]][["veg_class"]], fact=4, fun=modal)
rastercovs[["winter"]][["veg_class"]] <- projectRaster(rastercovs[["winter"]][["veg_class"]],rastercovs[["winter"]][["slope"]],method="ngb")
rastercovs[["winter"]][["veg_class"]] <- ratify(rastercovs[["winter"]][["veg_class"]])
rat <- levels(rastercovs[["winter"]][["veg_class"]])[[1]]
rat$landcover <- subs$HabType2
rat$code <- subs$Habitat_Type
levels(rastercovs[["winter"]][["veg_class"]]) <- rat

res(rastercovs[["winter"]][["veg_class"]])


###############
# Further GIS data processing
###############

rastercovs$summer <- stack(rastercovs$summer)   # done
rastercovs$winter <- stack(rastercovs$winter)   # done

plot(rastercovs$summer)
plot(rastercovs$winter)


###############
# Make prediction for summer habitat using GLMM 
###############

#GLMMs[["summer"]]

rastercovs$summer$stand_cos_aspect <- (rastercovs$summer$cos_aspect- mean(deer$summer$cos_aspect))/sd(deer$summer$cos_aspect)
plot(rastercovs$summer$stand_cos_aspect)

rastercovs$summer$stand_sin_aspect <- (rastercovs$summer$sin_aspect- mean(deer$summer$sin_aspect))/sd(deer$summer$sin_aspect)

rastercovs$summer$stand_elevation <- (rastercovs$summer$elevation- mean(deer$summer$elevation))/sd(deer$summer$elevation)

rastercovs$summer$stand_slope <- (rastercovs$summer$slope- mean(deer$summer$slope))/sd(deer$summer$slope)

rastercovs$summer$stand_dist_to_water <- (rastercovs$summer$dist_to_water- mean(deer$summer$dist_to_water))/sd(deer$summer$dist_to_water)

is.factor(rastercovs$summer$veg_class)
rat$classes <- rat$landcover
levels(rastercovs$summer$veg_class) <- rat[,c("ID","classes")]

predmaps <- list()
predmaps[["summer"]] <- list()

cat_cov <- c("veg_class")
veg_class <- levels(rastercovs$summer$veg_class)
names(veg_class) <- "veg_class"

factor_list <- list(veg_class)
names(factor_list) <- cat_cov

altid<-factor( 'P003',levels=levels(deer[["summer"]]$altid))
add2<-data.frame(altid)
str(add2)

# predmaps[["summer"]][["GLMM"]]  <- raster::predict(object=rastercovs$summer,model=GLMMs[["summer"]],const=add2,type="prob",factors=levels(rastercovs$summer$veg_class))  # factor_list progress='text',factors=factor_list,
# 
# 
# 
# predmaps[["summer"]][["RF"]]  <- raster::predict(object=rastercovs$summer,model=RFs[["summer"]])  # factors=factor_list

#### can't seem to use predict function- need to use another strategy...

newdata <- data.frame(
  cos_aspect = getValues(rastercovs$summer$cos_aspect),
  sin_aspect = getValues(rastercovs$summer$sin_aspect),
  elevation = as.integer(getValues(rastercovs$summer$elevation)),
  slope = as.integer(getValues(rastercovs$summer$slope)),
  dist_to_water = getValues(rastercovs$summer$dist_to_water),
  veg_class = factor(veg_class$veg_class$classes[match(getValues(rastercovs$summer$veg_class),veg_class$veg_class$ID)],levels=levels(deer$summer$veg_class)),
  stand_cos_aspect = getValues(rastercovs$summer$stand_cos_aspect),
  stand_sin_aspect = getValues(rastercovs$summer$stand_sin_aspect),
  stand_elevation = getValues(rastercovs$summer$stand_elevation),
  stand_slope = getValues(rastercovs$summer$stand_slope),
  stand_dist_to_water = getValues(rastercovs$summer$stand_dist_to_water)
)

newdata$altid = factor( 'P003',levels=levels(deer[["summer"]]$altid))

head(newdata)

newdata$veg_class[1:10]

## deal with NAs?

pred1 <- plogis(predict(GLMMs[["summer"]],newdata))     # predict to the full raster for summer range.  Does not take long!
predmaps[["summer"]][["GLMM"]] <- rastercovs$summer$elevation
predmaps[["summer"]][["GLMM"]] <- setValues(predmaps[["summer"]][["GLMM"]],pred1) 

plot(predmaps$summer$GLMM)

lenth <- ncell(rastercovs$summer$cos_aspect)
pred2 <- numeric(lenth)
jump = 10000
counter=1

while(counter<lenth){
  pred2[counter:min((counter+jump-1),lenth)] <- matrix(unlist(predict(RFs[["summer"]],newdata[counter:min((counter+jump-1),lenth),],type="prob",OOB=TRUE)),ncol=2,byrow=T)[,2]   # this takes a lot longer... OOB=TRUE,
  counter=counter+jump
  cat(sprintf("%d ...",counter-1))
}



 #pred2 <- matrix(unlist(predict(RFs[["summer"]],newdata[1:5,],type="prob",OOB=TRUE)),ncol=2,byrow=T)[,2]   # this takes a lot longer... 
predmaps[["summer"]][["RF"]] <- rastercovs$summer$elevation
predmaps[["summer"]][["RF"]] <- setValues(predmaps[["summer"]][["RF"]],pred2)


###############
# Make prediction for winter habitat using GLMM 
###############

#GLMMs[["winter"]]

rastercovs$winter$stand_cos_aspect <- (rastercovs$winter$cos_aspect- mean(deer$winter$cos_aspect))/sd(deer$winter$cos_aspect)    # something strange here...
plot(rastercovs$winter$stand_cos_aspect)

rastercovs$winter$stand_sin_aspect <- (rastercovs$winter$sin_aspect- mean(deer$winter$sin_aspect))/sd(deer$winter$sin_aspect)

rastercovs$winter$stand_elevation <- (rastercovs$winter$elevation- mean(deer$winter$elevation))/sd(deer$winter$elevation)

rastercovs$winter$stand_slope <- (rastercovs$winter$slope- mean(deer$winter$slope))/sd(deer$winter$slope)

rastercovs$winter$stand_dist_to_water <- (rastercovs$winter$dist_to_water- mean(deer$winter$dist_to_water))/sd(deer$winter$dist_to_water)

is.factor(rastercovs$winter$veg_class)
rat$classes <- rat$landcover
levels(rastercovs$winter$veg_class) <- rat[,c("ID","classes")]

predmaps[["winter"]] <- list()

cat_cov <- c("veg_class")
veg_class <- levels(rastercovs$winter$veg_class)
names(veg_class) <- "veg_class"

factor_list <- list(veg_class)
names(factor_list) <- cat_cov

altid<-factor( 'P003',levels=levels(deer[["winter"]]$altid))
add2<-data.frame(altid)
str(add2)

# predmaps[["winter"]][["GLMM"]]  <- raster::predict(object=rastercovs$winter,model=GLMMs[["winter"]],const=add2,type="prob",factors=levels(rastercovs$winter$veg_class))  # factor_list progress='text',factors=factor_list,
# 
# 
# 
# predmaps[["winter"]][["RF"]]  <- raster::predict(object=rastercovs$winter,model=RFs[["winter"]])  # factors=factor_list


#### can't seem to use predict function- need to use another strategy...

newdata <- data.frame(
  cos_aspect = getValues(rastercovs$winter$cos_aspect),
  sin_aspect = getValues(rastercovs$winter$sin_aspect),
  elevation = as.integer(getValues(rastercovs$winter$elevation)),
  slope = as.integer(getValues(rastercovs$winter$slope)),
  dist_to_water = getValues(rastercovs$winter$dist_to_water),
  veg_class = factor(veg_class$veg_class$classes[match(getValues(rastercovs$winter$veg_class),veg_class$veg_class$ID)],levels=levels(deer$winter$veg_class)),
  stand_cos_aspect = getValues(rastercovs$winter$stand_cos_aspect),
  stand_sin_aspect = getValues(rastercovs$winter$stand_sin_aspect),
  stand_elevation = getValues(rastercovs$winter$stand_elevation),
  stand_slope = getValues(rastercovs$winter$stand_slope),
  stand_dist_to_water = getValues(rastercovs$winter$stand_dist_to_water)
)

newdata$altid = factor( 'P003',levels=levels(deer[["winter"]]$altid))

head(newdata)

newdata$veg_class[1:10]


## deal with NAs?

pred1 <- plogis(predict(GLMMs[["winter"]],newdata))     # predict to the full raster for winter range.  Does not take long!
predmaps[["winter"]][["GLMM"]] <- rastercovs$winter$elevation
predmaps[["winter"]][["GLMM"]] <- setValues(predmaps[["winter"]][["GLMM"]],pred1) 

plot(predmaps[["winter"]][["GLMM"]])

lenth <- ncell(rastercovs$winter$cos_aspect)
pred2 <- numeric(lenth)
jump = 5000
counter=1

while(counter<lenth){
  pred2[counter:min((counter+jump-1),lenth)] <- matrix(unlist(predict(RFs[["winter"]],newdata[counter:min((counter+jump-1),lenth),],type="prob",OOB=TRUE)),ncol=2,byrow=T)[,2]   # this takes a lot longer... OOB=TRUE,
  counter=counter+jump
  cat(sprintf("%d ...",counter-1))
}

  #pred2 <- matrix(unlist(predict(RFs[["winter"]],newdata[1:5,],type="prob",OOB=TRUE)),ncol=2,byrow=T)[,2]   # this takes a lot longer... 
predmaps[["winter"]][["RF"]] <- rastercovs$winter$elevation
predmaps[["winter"]][["RF"]] <- setValues(predmaps[["winter"]][["RF"]],pred2)


## maps of comparison in predictions between the two methods

predmaps[["winter"]][["corLocal"]] <- corLocal(predmaps[["winter"]][["RF"]],predmaps[["winter"]][["GLMM"]],9)
predmaps[["summer"]][["corLocal"]] <- corLocal(predmaps[["summer"]][["RF"]],predmaps[["summer"]][["GLMM"]],9)


## maps of difference in predictions between the two methods

predmaps[["winter"]][["difference"]] <- predmaps[["winter"]][["RF"]]-predmaps[["winter"]][["GLMM"]]
predmaps[["summer"]][["difference"]] <- predmaps[["summer"]][["RF"]]-predmaps[["summer"]][["GLMM"]]

predmaps[["winter"]][["difference2"]] <- focal(predmaps[["winter"]][["difference"]], w=matrix(1, 5, 5), mean, na.rm=T)
predmaps[["summer"]][["difference2"]] <- focal(predmaps[["summer"]][["difference"]], w=matrix(1, 5, 5), mean,na.rm=T)

#focal(x, w=matrix(1, 5, 5), mean)


par(mfrow=c(1,2))
plot(predmaps[["winter"]][["RF"]])
plot(predmaps[["winter"]][["GLMM"]])


par(mfrow=c(1,2))
plot(predmaps[["winter"]][["corLocal"]])
plot(predmaps[["winter"]][["difference2"]])

par(mfrow=c(2,1))
plot(predmaps[["summer"]][["RF"]])
plot(predmaps[["summer"]][["GLMM"]])

par(mfrow=c(2,1))
plot(predmaps[["summer"]][["corLocal"]])
plot(predmaps[["summer"]][["difference2"]])

### save prediction maps as data and gis layers

allPoints <- list()
allPoints[["summer"]] <- readOGR("E:\\Dropbox\\Mule Deer\\Methods paper\\GIS\\Rasters","SummerUsed")
allPoints[["winter"]] <- readOGR("E:\\Dropbox\\Mule Deer\\Methods paper\\GIS\\Rasters","WinterUsed")

plot(rastercovs$summer$elevation)
plot(predmaps$summer$RF)
plot(predmaps$summer$GLMM)
plot(allPoints$summer[seq(1,nrow(allPoints$summer),by=20),],pch=1,cex=0.5,add=T,col=gray(0.1))


### do class check..

sapply(deer$summer,class)
sapply(newdata,class)



## visualize predictions from the two methods ...


## first read in map of points 

?plot.raster

breakpoints <- c(0,0.2,0.4,0.5,0.6,0.7,0.8,0.9,1)
cols <- brewer.pal(length(breakpoints),"Greens")  #   "RdYlGn"
  
#par(mfrow=c(1,2))
graphics.off()
setwd("E:\\Dropbox\\Mule Deer\\Methods paper\\EcolAndEcol_submitted\\figures")
svg(width=5,height=3,file="GLMM_SummerPredictions1.svg")
prettymap(plot(predmaps[["summer"]][["GLMM"]],xaxt="n",yaxt="n",breaks=breakpoints,col=cols),drawarrow = F,mai=c(0,0,0,0.1),bty="n",drawbox = F)
plot(allPoints$summer[seq(1,nrow(allPoints$summer),by=30),],pch=1,cex=0.5,add=T,col=gray(0.1))
dev.off()


breakpoints <- c(0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,1)
cols <- brewer.pal(length(breakpoints),"Greens")  #   "RdYlGn"

graphics.off()
setwd("E:\\Dropbox\\Mule Deer\\Methods paper\\EcolAndEcol_submitted\\figures")
svg(width=5,height=3,file="RF_SummerPredictions1.svg")
prettymap(plot(predmaps[["summer"]][["RF"]],xaxt="n",yaxt="n",breaks=breakpoints,col=cols),drawarrow = F,mai=c(0,0,0,0.1),bty="n",drawbox = F)
plot(allPoints$summer[seq(1,nrow(allPoints$summer),by=30),],pch=1,cex=0.5,add=T,col=gray(0.1))
dev.off()

breakpoints <- c(-1,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1)
cols <- brewer.pal(length(breakpoints),"RdYlGn")  #   ""

graphics.off()
setwd("E:\\Dropbox\\Mule Deer\\Methods paper\\EcolAndEcol_submitted\\figures")
svg(width=5,height=3,file="summer_corLocal.svg")
prettymap(plot(predmaps[["summer"]][["corLocal"]],xaxt="n",yaxt="n",breaks=breakpoints,col=cols),drawarrow = F,mai=c(0,0,0,0.1),bty="n",drawbox = F)
plot(allPoints$summer[seq(1,nrow(allPoints$summer),by=30),],pch=1,cex=0.5,add=T,col=gray(0.1))
dev.off()
graphics.off()

breakpoints <- c(-0.8,-0.5,-0.2,0,0.1,0.2,0.4)
cols <- brewer.pal(length(breakpoints),"PuOr")[length(breakpoints):1]  #   ""

graphics.off()
setwd("E:\\Dropbox\\Mule Deer\\Methods paper\\EcolAndEcol_submitted\\figures")
svg(width=5,height=3,file="summer_difference.svg")
prettymap(plot(predmaps[["summer"]][["difference2"]],xaxt="n",yaxt="n",breaks=breakpoints,col=cols),drawarrow = F,mai=c(0,0,0,0.1),bty="n",drawbox = F)
plot(allPoints$summer[seq(1,nrow(allPoints$summer),by=30),],pch=1,cex=0.5,add=T,col=gray(0.1))
dev.off()
graphics.off()

breakpoints <- c(0,0.2,0.3,0.5,0.6,0.7,0.8,0.9,1)
cols <- brewer.pal(length(breakpoints),"Greens")  #   "RdYlGn"

#par(mfrow=c(1,2))
graphics.off()
setwd("E:\\Dropbox\\Mule Deer\\Methods paper\\EcolAndEcol_submitted\\figures")
svg(width=3,height=7,file="GLMM_WinterPredictions1.svg")
prettymap(plot(predmaps[["winter"]][["GLMM"]],xaxt="n",yaxt="n",breaks=breakpoints,col=cols),drawarrow = F,mai=c(0,0,0,0.1),bty="n",drawbox = F)
plot(allPoints$winter[seq(1,nrow(allPoints$winter),by=20),],pch=1,cex=0.5,add=T,col=gray(0.1))
dev.off()

breakpoints <- c(0,0.2,0.3,0.4,0.5,0.6,0.7,0.8,1)
cols <- brewer.pal(length(breakpoints),"Greens")  #   "RdYlGn"

graphics.off()
setwd("E:\\Dropbox\\Mule Deer\\Methods paper\\EcolAndEcol_submitted\\figures")
svg(width=3,height=7,file="RF_WinterPredictions1.svg")
prettymap(plot(predmaps[["winter"]][["RF"]],xaxt="n",yaxt="n",breaks=breakpoints,col=cols),drawarrow = F,mai=c(0,0,0,0.1),bty="n",drawbox = F)
plot(allPoints$winter[seq(1,nrow(allPoints$winter),by=20),],pch=1,cex=0.5,add=T,col=gray(0.1))
dev.off()


breakpoints <- c(-1,-0.75,-0.5,-0.25,0,0.25,0.5,0.75,1)
cols <- brewer.pal(length(breakpoints),"RdYlGn")  #   ""

graphics.off()
setwd("E:\\Dropbox\\Mule Deer\\Methods paper\\EcolAndEcol_submitted\\figures")
svg(width=3,height=7,file="winter_corLocal.svg.svg")
prettymap(plot(predmaps[["winter"]][["corLocal"]],xaxt="n",yaxt="n",breaks=breakpoints,col=cols),drawarrow = F,mai=c(0,0,0,0.1),bty="n",drawbox = F)
plot(allPoints$winter[seq(1,nrow(allPoints$winter),by=20),],pch=1,cex=0.5,add=T,col=gray(0.1))
dev.off()

breakpoints <- c(-0.6,-0.4,-0.2,0,0.1,0.2,0.6)
cols <- brewer.pal(length(breakpoints),"PuOr")[length(breakpoints):1]  #   ""

graphics.off()
setwd("E:\\Dropbox\\Mule Deer\\Methods paper\\EcolAndEcol_submitted\\figures")
svg(width=3,height=7,file="Winter_difference.svg")
prettymap(plot(predmaps[["winter"]][["difference2"]],xaxt="n",yaxt="n",breaks=breakpoints,col=cols),drawarrow = F,mai=c(0,0,0,0.1),bty="n",drawbox = F)
plot(allPoints$winter[seq(1,nrow(allPoints$winter),by=20),],pch=1,cex=0.5,add=T,col=gray(0.1))
dev.off()
graphics.off()



plot(rastercovs[["summer"]][[1:6]])


#### explore
plot(predmaps[["summer"]][["GLMM"]])
z1 <- drawExtent()
z2 <- crop(rastercovs[["summer"]][[1:6]],y=z1)
plot(z2)


veg_class <- levels(rastercovs$summer$veg_class)
names(veg_class) <- "veg_class"

factor_list <- list(veg_class)
names(factor_list) <- cat_cov

altid<-factor( 'P003',levels=levels(deer[["summer"]]$altid))
add2<-data.frame(altid)
str(add2)
newdata <- data.frame(
  cos_aspect = mean(getValues(rastercovs$summer$cos_aspect)),
  sin_aspect = mean(getValues(rastercovs$summer$sin_aspect)),
  elevation = 1773,
  slope = as.integer(mean(getValues(rastercovs$summer$slope))),
  dist_to_water = mean(getValues(rastercovs$summer$dist_to_water)),
  veg_class = factor(veg_class$veg_class$classes[match(28,veg_class$veg_class$ID)],levels=levels(deer$summer$veg_class)),
  stand_cos_aspect = mean(getValues(rastercovs$summer$stand_cos_aspect)),
  stand_sin_aspect = mean(getValues(rastercovs$summer$stand_sin_aspect)),
  stand_elevation = (0 - mean(deer$summer$elevation))/sd(deer$summer$elevation),
  stand_slope = mean(getValues(rastercovs$summer$stand_slope)),
  stand_dist_to_water = mean(getValues(rastercovs$summer$stand_dist_to_water))
)

plogis(predict(GLMMs[["summer"]],newdata))


# apparently the strong preference for summer region A is due to the interaction effect between slope and elevation (Fig. X)
#  whereby low elevation, high-slope areas are strongly preferred. 


summary(GLMMs[["summer"]])






#### explore winter
plot(predmaps[["winter"]][["difference2"]])
z1 <- drawExtent()
z2 <- crop(rastercovs[["winter"]][[1:6]],y=z1)
plot(z2)

plot(predmaps[["winter"]][["difference2"]])
plot(z1,add=T)

plot(rastercovs$winter$elevation)
plot(z1,add=T)

plot(rastercovs$winter$slope)
plot(z1,add=T)

plot(rastercovs$winter$dist_to_water)
plot(z1,add=T)

table(z2$veg_class@data@values)   # focal area is mostly P-J

veg_class <- levels(rastercovs$winter$veg_class)
names(veg_class) <- "veg_class"

factor_list <- list(veg_class)
names(factor_list) <- cat_cov

altid<-factor( 'P005',levels=levels(deer[["winter"]]$altid))
add2<-data.frame(altid)
str(add2)
newdata <- data.frame(
  cos_aspect = mean(getValues(rastercovs$winter$cos_aspect)),
  sin_aspect = mean(getValues(rastercovs$winter$sin_aspect)),
  elevation = as.integer(1900),
  slope = as.integer(60),
  dist_to_water = 2500,
  veg_class = factor(veg_class$veg_class$classes[match(8,veg_class$veg_class$ID)],levels=levels(deer$winter$veg_class)),
  stand_cos_aspect = mean(getValues(rastercovs$winter$stand_cos_aspect)),
  stand_sin_aspect = mean(getValues(rastercovs$winter$stand_sin_aspect)),
  stand_elevation = (1900 - mean(deer$winter$elevation))/sd(deer$winter$elevation),
  stand_slope = (60 - mean(deer$winter$slope))/sd(deer$winter$slope),
  stand_dist_to_water = (2500 - mean(deer$winter$dist_to_water))/sd(deer$winter$dist_to_water),
  altid = altid
)

plogis(predict(GLMMs[["winter"]],newdata))   ## very high

predict(RFs[["winter"]],newdata,OOB=T,type="prob")   ## 


lapply(deer$winter,class)

lapply(newdata,class)


# apparently the strong preference for winter region B is due to the interaction effect between slope and elevation (Fig. X)
#  whereby low elevation, high-slope areas are strongly preferred. 


getwd()

save(deer,file="deerdata.RData")


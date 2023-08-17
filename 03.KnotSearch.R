##################################################################
# Knot Search 
# - For Lag-Outcome Spline Function Specification in Cross-Basis
##################################################################
#
# "Identifying Sensitive Windows of Exposure to NO2 and 
#             Fetal Growth Trajectories in a Spanish Birth Cohort"
#
#    Kristina W. Whitworth,Alison Rector, Jennifer Ish, 
#        Suneet P. J. Chauhan, Jesús Ibarluzea, Mònica Guxens,
#        Michael D. Swartz, Elaine Symanski, and Carmen Iñiguez
#
#    Epidemiology - May 2022 (Volume 33, Number 3)
#
######################################################################  

#########################
# Load Packages 
#########################
library(INLA)
library(data.table)
library(tidyverse)
library(sf)
library(sp)
library(spdep)
library(dlnm)
library(tsModel)
library(dplyr)
library(foreign)
library(reshape2)
library(janitor)
library(kableExtra)
library(magrittr)
library(summarytools)
library(table1)

##########################
# Read in Data
#######################
NO2_data<-readRDS("CleanData/NO2_global_complete.rds")

#######################################################
# Define Variable lists (Exposure, Outcome, Covariates)  
#########################################################
# Covariates 
covariates_NO2<-c( "m_age", #Maternal age at birth
                   "m_bmi_new", #pre pregnancy BMI
                   "parity3c", #number of previous pregnancies (categorical 3 levels)
                   "marital2c", #Cohabitation - living situation
                   "SC", #combined social class
                   "alc",
                   "smokpreg", 
                   "m_educ", # Maternal Education
                   "p_educ" #paternal education
) 
# NO2 Weekly Exposures
endweek<-34
exposure_NO2<-c(paste0("NO2_w",1:endweek))

#Fetal Growth Outcome Measures 
outcome_names<-c( "ac12",  "ac20_12","ac34_20",
                  "bpd12","bpd20_12","bpd34_20",
                  "efw12","efw20_12", "efw34_20",
                  "fl12","fl20_12","fl34_20")

########################################
# Create Datasets for each outcome
#######################################
#NOTE - this is not strictly necessary, if step is skipped,
# map knot search to full dataset
#create exposure lists specific to weeks for fetal growth trajectories
exposure_12_wk<-c(paste0("NO2_w",1:12))
exposure_20_wk<-c(paste0("NO2_w",1:20))
exposure_34_wk<-c(paste0("NO2_w",1:34))

#make list of exposure weeks for creating core dataframes 
#(12wk, 20wk and 34 wk for each of the 4 fetal growth measurements)
exposures_list<-list(exposure_12_wk,exposure_20_wk,exposure_34_wk,
                     exposure_12_wk,exposure_20_wk,exposure_34_wk,
                     exposure_12_wk,exposure_20_wk,exposure_34_wk,
                     exposure_12_wk,exposure_20_wk,exposure_34_wk)

#create core df function
create_core_df<-function(curr_outcome, curr_exposure){
  core.df<-NO2_data%>%dplyr::select(cohort,
                                    all_of(curr_outcome),
                                    all_of(curr_exposure),
                                    all_of(covariates_NO2),
                                    urban)
  return(core.df)
}
#map over outcomes and exposure lists to make base formulas
core_df_list<-pmap(list(curr_outcome=outcome_names,
                        curr_exposure=exposures_list),
                   create_core_df)
#name formulas in list
names(core_df_list)<-outcome_names

#############################
# Define Knot search function
#############################
library(mvmeta)

knotsearch<-function(myresknots, mydata_knot, n_pers, exposure_per, outcome){
  
  resknots<-myresknots
  
  #create new variable called "test_outcome" same as outcome variable
  mydata_knot$test_outcome=mydata_knot[[outcome]]
  #get the number of rows to loop over from resknots matrix
  tot.knot.vect<-dim(resknots)[[1]]
  #add RMSE column to resknots
  # resknots$RMSE_val<-rep(NA, tot.knot.vect)
  #specify number of lags and weeks
  nlag<-n_pers-1
  nweeks<-n_pers
  #Step 1: Create model for each row of Resknot 
  #set vector to use for increments of exposure
  cenat=c(0,10)
  #specify xlag
  
  core.df<-mydata_knot
  # - store WAIC in Resknot
  #k<-1
  #Step 4: Create Cohort and Meta model for each row of Resknot 
  # - store AIC in Resknot
  
  for(k in 1:tot.knot.vect){
    #select the kth value in resknots$nknots
    curr.nknots<-resknots$nknots[k]
    #set nk to curr.nknots (for use in storring coefficients of cohort model)
    nk<-curr.nknots
    #get current knot vector, begins with column 2 and includes number of col
    # in knot vector, (add curr.nknots-1 to 2)
    curr.knot.vector<-as.numeric(resknots[k,2:(2+curr.nknots-1)])
    lkji<-curr.knot.vector
    
    
    #SETUP for Looping through Cohort models:
    #set vector to use for increments of exposure
    cenat=c(0,10)
    #specify xlag
    xlag<-0:nlag
    coh<-c("Gipuzkoa",  "Sabadell",  "Valencia") 
    ncoh<-length(coh)
    #specify argvarm - the spline type for the predictor dimension of crossbasis
    argvarm<-list(fun="lin")
    #specify arglagm - the spline type for the lag dimension of the crossbasis
    #specify if there is an intercept, here we choose no
    arglagm<-list(fun="ns", knots=lkji, int=F)
    #create matrix of exposure variables for relavent weeks
    Qx<-as.matrix(core.df[,paste("NO2_w", nweeks:1, sep="")])
    #specify the cross basis in the pooled model
    cbxi<-crossbasis(Qx,
                     lag=c(0,nlag),
                     argvar=argvarm,
                     arglag=arglagm)
    
    
    
    coef_curr<- matrix(NA,ncoh,nk+1,dimnames=list(coh))
    # note: dimensions depends of nk and the type of basis, for example if bs, nk+3.
    vcov_curr <- vector("list",ncoh)
    names(vcov_curr ) <- coh
    
    ####################
    # COHORT MODEL LOOP
    ###################
    
    cohort_model_list<-list()
    #m<-3
    for(m in 1:3){
      Qxi<-NULL
      # print(paste0("cohort ", m, ": ", coh[m]))
      cenat=c(0,10)
      #filter to get cohort specific dataframe from cohort data
      my_cohort_data<-core.df%>%filter(cohort==coh[m])
      Qxi<-as.matrix(my_cohort_data[,paste("NO2_w", nweeks:1,sep="")]) # Qxi en  cohorte i
      my_cbi<-crossbasis(Qxi, lag=c(0,nlag), argvar=argvarm, arglag=arglagm)
      
      if(m==2){ 
        modi<-lm(test_outcome ~ 
                   m_age + 
                   m_bmi_new +
                   parity3c +
                   marital2c +
                   SC +
                   alc +
                   smokpreg +
                   m_educ +
                   p_educ +
                   my_cbi,
                 data=my_cohort_data)
      } else {
        modi<-lm(test_outcome ~ 
                   m_age + 
                   m_bmi_new +
                   parity3c +
                   marital2c +
                   SC +
                   alc +
                   smokpreg +
                   m_educ +
                   p_educ +
                   urban +
                   my_cbi,
                 data=my_cohort_data)}
      
      
      #note: this 'if' is necessary because Sabadell don't have rural zone!!! 
      redi<-crossreduce(my_cbi,modi,value=cenat[2],type="var",cen=cenat[1]) #recuding the basis
      rm(my_cbi)
      coef_curr[m,]<-redi$coef    #saving cohort-specific coef and vcov
      vcov_curr[[m]]<-redi$vcov
      
      cohort_model_list[[m]]<-redi
      
    }
    
    names(cohort_model_list)<-coh
    modji<- mvmeta(coef_curr~1,vcov_curr,method="reml",control=list(showiter=F))    
    AIC_curr<-AIC(modji,k=2)   
    resknots$AIC[k]<-AIC_curr
    #Generate RMSE value and store
    # library(MLmetrics)
    # RMSE_curr<-RMSE(modji$summary.fitted.values$mean, mydata_knot[[outcome]])
    # resknots$RMSE_val[m]<-RMSE_curr
  }
  
  
  #get best knot vector
  
  knots.bestAIC <- resknots[which.min(resknots$AIC),]
  #Get a vector of ideal knots
  n<-knots.bestAIC[1,1]
  knot.vector<-NULL
  for(m in 1:n)
  {
    place<-m+1
    knot.vector<-c(knot.vector,knots.bestAIC[place])
  }
  knot.vector<-as.numeric(knot.vector)
  # print("knot search complete")
  return(list(knot.matrix=resknots, best.knots=knot.vector))
}

##############################
# Map Knot Search to Outcomes
###############################

#list of core datasets - alternative is to repeat NO2 dataset 12 times in this list
mydata_list<-list(core_df_list$ac12, core_df_list$ac20_12,
                  core_df_list$ac34_20, core_df_list$bpd12, 
                  core_df_list$bpd20_12, core_df_list$bpd34_20,
                  core_df_list$efw12, core_df_list$efw20_12,
                  core_df_list$efw34_20, core_df_list$fl12, 
                  core_df_list$fl20_12, core_df_list$fl34_20)

# save core data frames list
saveRDS(core_df_list,"CleanData/core_df_list_fetalgrowth_NO2.rds")


#specify exposure variable - week numbers are added in knot search code
exposure_pers=list("NO2_w")

#for each outcome, specify number of exposure periods/lag periods
n_pers_list=as.list(c(12, 20, 34,
                      12, 20, 34,
                      12, 20, 34, 
                      12, 20, 34))

#read in rigid knot search grids - created in '02.KnotSearchGrids.R'
knot_mat_12_rigid<-readRDS("CleanData/knot_mat_12_rigid.rds")
knot_mat_20_rigid<-readRDS("CleanData/knot_mat_20_rigid.rds")
knot_mat_34_rigid<-readRDS("CleanData/knot_mat_34_rigid.rds")

#for each outcome, specify the knot search grid to use
resknot_list<-list(knot_mat_12_rigid, knot_mat_20_rigid, knot_mat_34_rigid,
                   knot_mat_12_rigid, knot_mat_20_rigid, knot_mat_34_rigid,
                   knot_mat_12_rigid, knot_mat_20_rigid, knot_mat_34_rigid,
                   knot_mat_12_rigid, knot_mat_20_rigid, knot_mat_34_rigid)

#create list of outcome variables
outcome_names<-c( "ac12",  "ac20_12","ac34_20",
                  "bpd12","bpd20_12","bpd34_20",
                  "efw12","efw20_12", "efw34_20",
                  "fl12","fl20_12","fl34_20")

#map over outcomes to make base formulas - 'pmap' function requires 'purrr' package
knot_search_results<-pmap(list(myresknots=resknot_list,
                               mydata_knot=mydata_list,
                               n_pers=n_pers_list,
                               exposure_per=c("NO2_w"),
                               outcome=outcome_names),
                          knotsearch)

#add names to list of knot search results
names(knot_search_results)<-outcome_names

#save the full list of knot search results
saveRDS(knot_search_results,"KnotSearch/FetalGrowth_knots.rds")



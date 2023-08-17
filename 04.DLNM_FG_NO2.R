#################################################
# DLNM Estimations
#  Create DLNM fits for each fetal growth outcome
#################################################
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
library(RColorBrewer)
library(geofacet)
library(ggpubr)
library(ggplot2)
library(ggthemes)
library(dplyr)
library(mvmeta)
################
# Read in Data
###############
# read core data frames list
core_df_list<-readRDS("CleanData/core_df_list_fetalgrowth_NO2.rds")

#read in best knot list
knot_search_list<-readRDS("KnotSearch/FetalGrowth_knots.rds")

outcome_names<-c( "ac12",  "ac20_12","ac34_20",
                  "bpd12","bpd20_12","bpd34_20",
                  "efw12","efw20_12", "efw34_20",
                  "fl12","fl20_12","fl34_20")
names(knot_search_list)<-outcome_names

####################################
# Calculate Mean and SD of Outcomes 
## - used for calculating % change 
######################################  

#read in raw data 
load("RawData/Reference_Data_PercChange.RData")
Refs<-Reference_Data_PercChange


#Using predicted values at 12,20 and 34 weeks:

mean_sd_summary <- Refs %>% group_by(OUTCOME) %>% 
  summarise(mean12 = mean(Pred12),
            sd12=sd(Pred12),
            mean20 = mean(Pred20),
            sd20=sd(Pred20),
            mean34 = mean(Pred34),
            sd34=sd(Pred34))

mean_sd_summary_df<-as.data.frame((mean_sd_summary))

mean_sd_summary_df<-mean_sd_summary_df%>%
  mutate(sd_mean_12=sd12/mean12)%>%
  mutate(sd_mean_20=sd20/mean20)%>%
  mutate(sd_mean_34=sd34/mean34)

mean_sd_summary_df$OUTCOME<-ifelse(mean_sd_summary_df$OUTCOME=="whad", "efw", 
                                   ifelse(mean_sd_summary_df$OUTCOME=="dbp", "bpd", 
                                          ifelse(mean_sd_summary_df$OUTCOME=="lf", "fl", "ac")))
#mean_sd_summary_df


#######################################
#########    List of Models   #########
#######################################

outcomes<-c("efw", "bpd", "fl", "ac")
weeks<-c(12,20,34)
# loop over number of outcomes and time points to generate list of models using
# NO2.DLNM

mean_sd_df<-data.frame(outcome=character(), mean_sd_ratio=numeric())
for(i in 1:4){
  for(j in 1:3){
    #chooses the number of weeks from weeks vector
    endweek<-weeks[j]
    outcome_curr<-outcomes[i]
    myoutcome<-paste(outcomes[i], weeks[j], sep="")
    #run mean_sd_complete
    data<-mean_sd_summary_df%>%filter(mean_sd_summary_df$OUTCOME==outcome_curr)%>%select(matches(paste("sd_mean_",endweek,sep="")))
    curr.mean_sd<-data[1,1]
    curr.mean_sd_df<-data.frame(outcome=myoutcome,mean_sd_ratio=curr.mean_sd)
    
    #add current model to the list of models
    mean_sd_df<-bind_rows(mean_sd_df, curr.mean_sd_df)
  }
}

saveRDS(mean_sd_df, "CleanData/mean_sd_vars.rds")

##################################
# FUNCTION TO BUILD META MODEL
#################################
mymodel <- function(mydata_knot, myknot_search, n_pers, outcome){
  # specify best knots
  mybest.knots<-myknot_search$best.knots
  #specify number of lags and weeks
  nlag<-n_pers-1
  nweeks<-n_pers
  #create new variable called "test_outcome" same as outcome variable
  mydata_knot$test_outcome=mydata_knot[[outcome]]
  #specify core data frame
  core.df<-mydata_knot
  #nk is number of knots, length of knot vector
  nk<-length(mybest.knots)
  #########################################
  #SETUP for Looping through Cohort models:
  #########################################
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
  arglagm<-list(fun="ns", knots=mybest.knots, int=F)
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
    redi<-crossreduce(my_cbi,
                      modi,
                      value=cenat[2],
                      #value = 1,
                      type="var",cen=cenat[1]) #recuding the basis
    rm(my_cbi)
    coef_curr[m,]<-redi$coef    #saving cohort-specific coef and vcov
    vcov_curr[[m]]<-redi$vcov
    
    cohort_model_list[[m]]<-redi
    
  }
  
  names(cohort_model_list)<-coh
  curr_meta<- mvmeta(coef_curr~1,vcov_curr,method="reml",control=list(showiter=F)) 
  
  
  #####################
  # Predictions
  #################
  #get the basis with lag attributes
  #blag <- do.call("onebasis",c(list(x=xlag),attr(cbxi,"arglag")))
  
  pred_curr<-crosspred(basis=cbxi,
                       coef = coef(curr_meta),
                       vcov = vcov(curr_meta),
                       bylag=1, #sequence one lag at a time along cross basis
                       at=cenat[2], #calculate at 10 micrograms/m^3 -> cenat=(0,10) -> cenat[2]=10
                       cen=cenat[1]) #center exposure at zero ->cenat=(0,10) -> cenat[1]=0
  
  #save model and pred_cur as .rds files in data folder
  saveRDS(curr_meta, paste0("DLNM_Data/dlnm_meta_model_", outcome, ".rds"))
  saveRDS(pred_curr, paste0("DLNM_Data/pred_data_", outcome, ".rds"))
  
  return(list(dlnm_model=curr_meta, 
              pred_data_forplot=pred_curr,
              pred_data_sums=pred_curr
  ))
}

#######################################
# MAP META BUILD FUNCTION TO OUTCOMES
########################################
n_pers_list=as.list(c(12, 20, 34,
                      12, 20, 34,
                      12, 20, 34, 
                      12, 20, 34))




#create list of outcome variables
outcome_names<-c( "ac12",  "ac20_12","ac34_20",
                  "bpd12","bpd20_12","bpd34_20",
                  "efw12","efw20_12", "efw34_20",
                  "fl12","fl20_12","fl34_20")

#map over outcomes to make base formulas
NO2_models<-pmap(list(mydata_knot=core_df_list,
                      myknot_search=knot_search_list,
                      n_pers=n_pers_list,
                      outcome=outcome_names),
                 mymodel)
names(NO2_models)<-outcome_names
saveRDS(NO2_models,"DLNM_Data/NO2_models.rds")

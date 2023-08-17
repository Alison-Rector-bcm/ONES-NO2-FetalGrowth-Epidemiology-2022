#####################################################################
# Data Preparation for Fetal Growth ~ NO2
##################################################################### 
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

##################
# Load Packages 
##################
library(foreign)
library(dplyr)
library(reshape2)
library(janitor)
library(kableExtra)
library(magrittr)
library(summarytools)
library(table1)
library(flextable)

#############################
# Outcome and Covariate Data
##############################

# define function to seelct variables of interest and rename 
# variables using namekey, returning updated dataset
response_data_prep<-function(raw_response_data, varlist, namekey){
  response.data<-raw_response_data%>%dplyr::select(any_of(varlist))
  response.data<-response.data%>%plyr::rename(replace=namekey, warn_missing=FALSE)
  return(response.data)
}

##########################################
#   Read in exposure and outcome data 
##########################################

#load Fetal Data 2
Fetal2.data<-read.dta("RawData/MatOc_ED_FetalD_NeuroD_2.dta")

#Specify exposure variables
NO2_exp_varlist<-paste0("NO2_w", seq(1:34))

#Specify list of covariates, outcome variables and cohort/id number as keys
NO2_varlist<-c("cohorte", #cohort
               "idnum", #participant ID number
               "birth",
               "edadm", #maternal age
               "imcm4c", #maternal bmi - 4 categories
               "imcm", #continuous maternal bmi
               "paridad3c",  #parity (number of previous pregnancies)
               "sitmarital2c", #cohabitation
               "CSM3", #maternal social class
               "CSP3", #paternal social class
               "alc", # maternal alcohol
               "smokpreg", #maternal smoking during pregnancy
               "estudios3c", #maternal education
               "estudiosp3c", #paternal education
               "tipozonaM", #urbanicity
               "zdbp12", "zdbp20_12", "zdbp34_20", # biparietal diameter -12 wks, 12-20 wks and 20-34 wks gestation
               "zlf12",  "zlf20_12",  "zlf34_20", # Femur Length - 12 wks, 12-20 wks and 20-34 wks gestation
               "zpa12","zpa20_12", "zpa34_20", # Abdominal Circumference -  12 wks, 12-20 wks and 20-34 wks gestation 
               "zwhad12", "zwhad20_12", "zwhad34_20") # Estimated Fetal Weight -  12 wks, 12-20 wks and 20-34 wks gestation

#define namekey to rename variables
namekey<-c(cohorte = "cohort",
           sitmarital2c="marital2c", # marital status
           edadm = "m_age", #maternal age
           estudios3c = "m_educ", #maternal education
           estudiosp3c = "p_educ", # paternal education 
           CSM3 = "SCM", #social class - maternal
           CSP3 = "SCP", #social class - paternal,
           tipozonaM ="urban", # urbanicity - 3 levels
           pd_semejanzas_m_m10sd3 = "IQ_m_m10sd4", #maternal IQ
           sexo = "sex",#child sex
           imcm4c = "m_bmi_4cat", # pre pregnancy BMI - 4 categories
           imcm = "m_bmi", # pre pregnancy BMI - continuous
           paridad3c = "parity3c", # number of previous pregnancies (live birth or fetal death >=22 weeks gestation)
           edMcCarthy_4y = "ch_age", #child age at mccarthy assessment
           zdbp12 = "bpd12",  zdbp20_12 = "bpd20_12", zdbp34_20 = "bpd34_20",
           zlf12 = "fl12",  zlf20_12 = "fl20_12", zlf34_20 = "fl34_20",
           zpa12 = "ac12", zpa20_12 = "ac20_12", zpa34_20 = "ac34_20",
           zwhad12 = "efw12",zwhad20_12 = "efw20_12", zwhad34_20 = "efw34_20"
)

#specify list of variables of interest - exposure, outcome, covariates, observation identifiers
NO2_full_varlist<-c(NO2_varlist, NO2_exp_varlist)

#Use response_data_prep function to create subset of data with renamed variables
NO2.full.data<-response_data_prep(raw_response_data=Fetal2.data, 
                                  varlist=NO2_full_varlist, 
                                  namekey)


#############################
# Refactor Cohort Variable  
#############################

NO2.full.data$cohort<-ifelse(NO2.full.data$cohort=="Asturias", 1,
                             ifelse(NO2.full.data$cohort=="Gipuzkoa", 2,
                                    ifelse(NO2.full.data$cohort=="Sabadell",3,
                                           #comment out 99 for sample, comment out "," for full
                                           ifelse(NO2.full.data$cohort=="Valencia",4,99#,            
                                                  #comment out  comma when using full data
                                                  #ifelse(NO2.full.data$cohort=="Asturias",5,99)   
                                                  #comment out when using full data
                                           ))))

NO2.full.data$cohort<-factor(NO2.full.data$cohort,
                             levels=c(1,2,3,4,5),
                             labels=c("Asturias", "Gipuzkoa", "Sabadell", "Valencia", "Sabadell2"))

#################################
# Exclusion/Inclusion Criteria
################################

#exclude women from Sabadell with ID greater than 657,
# these women were recruited at birth and do not meet inclusion criteria
NO2.full.data<-NO2.full.data%>%
  dplyr::filter(!(cohort=="Sabadell" & idnum>657))

#exclude women who did not get followed to birth (for any reason)
NO2.full.data<-NO2.full.data%>%dplyr::filter(birth!="No")

# excluding observation 631 in Sabadell b/c missing age
NO2.full.data<-NO2.full.data%>%  
  filter(!(cohort=="Sabadell" & idnum==631))  

# Exclude Observations from Asturias if they are in your dataset
NO2.full.data<-NO2.full.data%>%  
  filter(!(cohort=="Astruias"))

##########################################################
# Specifying/Classifying Variables for this analysis
##########################################################

###############################
#Define combined social class
###############################
#create combined social class variable
NO2.full.data$SCM_num<-as.numeric(NO2.full.data$SCM)
NO2.full.data$SCP_num<-as.numeric(NO2.full.data$SCP)
#table variables for comparison
#table(NO2.full.data$SCM,NO2.full.data$SCM, exclude=NULL)
#table(NO2.full.data$SCP_num, NO2.full.data$SCP_num, exclude=NULL)
#set NA to 5
NO2.full.data$SCP_num[is.na(NO2.full.data$SCP_num)==TRUE]<-5
NO2.full.data$SCM_num[is.na(NO2.full.data$SCM_num)==TRUE]<-5
#Define combined SC
NO2.full.data<-NO2.full.data%>%
  mutate(SC=ifelse(NO2.full.data$SCM_num>=NO2.full.data$SCP_num, 
                   NO2.full.data$SCP,NO2.full.data$SCM))

#set SC as factor variable
NO2.full.data$SC<-factor(x=NO2.full.data$SC,
                         levels=c(1,2,3,4),
                         labels=c("High",
                                  "Middle",
                                  "Low", 
                                  "Missing"))

#############################
# Convert Smoking to Factor
#############################
#change smokepreg to factor variables and assign labels
NO2.full.data$smokpreg<-factor(x=NO2.full.data$smokpreg,
                               levels=c(0,1),
                               labels=c("Never Smoked",
                                        "Smoked"))

################################################
# Re-level BMI (Health as Reference Category)
##############################################

#manually change BMI so that normal is set to level 1 (reference group)
NO2.full.data<-NO2.full.data%>%mutate(m_bmi_new=ifelse(m_bmi_4cat=="2.Saludable (18.5<=imc<25)", 1, 
                                                       ifelse(m_bmi_4cat=="1.Bajo peso (imc<18.5)",2, m_bmi_4cat)))
NO2.full.data$m_bmi_new<-factor(x=NO2.full.data$m_bmi_new,
                                levels=c(1,2,3,4),
                                labels=c("Normal (18.5<= bmi <25)",
                                         "Underweight (bmi <18.5)",
                                         "Overweight (25<= bmi <30)",
                                         "Obese (bmi>=30)"))


###########################
# Complete Case Dataset  
###########################
#create list of outcome variables
outcome_names<-c( "ac12",  "ac20_12","ac34_20",
                  "bpd12","bpd20_12","bpd34_20",
                  "efw12","efw20_12", "efw34_20",
                  "fl12","fl20_12","fl34_20")
#create list of exposure names
exposure_NO2<-c(paste0("NO2_w",1:34))
#create list of covariate names
covariates_NO2<-c( "cohort",
                   "m_age", #Maternal age at birth
                   "m_bmi_new", #pre pregnancy BMI
                   "parity3c", #number of previous pregnancies (categorical 3 levels)
                   "marital2c", #Cohabitation - living situation
                   "SC", #combined social class
                   "alc",
                   "smokpreg", 
                   "m_educ", # Maternal Education
                   "p_educ", #paternal education
                   "urban"
) 
FG_NO2_vars<-c(outcome_names, exposure_NO2, covariates_NO2)
#select variables in outcome_names, exposure_NO2 and covariates_NO2

#create dataframe of outcome variables for NO2 fetal growth analysis
NO2.my.vars<-NO2.full.data%>%dplyr::select(all_of(FG_NO2_vars))

NO2.complete<-na.omit(NO2.my.vars)
saveRDS(NO2.complete, "CleanData/NO2_global_complete.rds")

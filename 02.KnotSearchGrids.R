##################################################################
# Create Knot Matrices 
# - For Knot Search for Lag-Outcome Spline Function in Cross-Basis
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

### Create and store knot matrices ###
create_knot_matrix<-function(n_pers){
  #######################################
  #########   FlexKNOT SELECTION    #########
  #######################################
  #Step 1: Determine number of possible knots to consider
  if(n_pers==9){
    nkl<-1:3
    mindistance<-2}
  
  if(n_pers==12){
    nkl<-1:3
    mindistance<-3}
  if(n_pers==20){
    nkl<-1:3
    mindistance<-3}
  if(n_pers==34){
    nkl<-1:3
    mindistance<-3}
  #Step 2: Find possible positions of knots 
  #First possible position is equal to mindistance 
  firstpl<-2
  #Last position is a rounded value no closer than mindistance from last lag
  #this is because for now we are only allowing knots to be located at 
  #multiples of mindistance from starting position
  lastpl<-n_pers-2
  #bydist is the same as mindistance for now - may change later to allow 
  #more flexible knot placements
  bydist<-1
  #create vector of possible knot locations, sequence of values from 
  #first location to last location by bydist
  posl<-seq(firstpl,lastpl,by=bydist) 
  #Step 3: Create grid of all possible knot vectors
  
  #determine the number of possible models to fit by nkl number of knots
  #(nmj<-choose(length(posl),nkl))
  #tot.knot.vect<-sum(nmj)
  
  #create empty data frame to store:
  # -number of knots
  # -knot placement  
  # -model AIC value
  #initiate vector of AIC value
  aicval<-NA
  #initialize empty matrix
  resknots<-data.frame(matrix(NA, nrow=0, ncol=2+length(nkl)))
  #generate knots to put in  
  #loop over possible number of knots to consider
  for(j in 1:length(nkl)){
    #create grid of possible positions for knots
    gridk<-combn(posl,nkl[j])
    #number of possible combinations (should match tot.knot.vect)
    nmax<-ncol(gridk)
    for(n in 1:nmax){
      curr_knots<-as.numeric(gridk[,n])
      resknot_pos<-c(nkl[j],curr_knots,rep(NA,max(nkl)-nkl[j]),NA)
      resknots<-rbind(resknots, resknot_pos)
    }
    #name variables of data.frame
    names(resknots)<-c("nknots",paste("knot",1:length(nkl),sep=""),"AIC")
    
  }
  #change na to 0
  resknots[is.na(resknots)]<-0
  if(length(nkl)==3){
    resknots<-resknots%>%
      mutate(dist_2_1=abs(knot2-knot1))%>%
      mutate(dist_3_2=abs(knot3-knot2))
    
    #remove any row that has a distance less than mindistance
    resknots<-resknots%>%
      filter((dist_2_1==0|dist_2_1>=mindistance)&
               (dist_3_2==0|dist_3_2>=mindistance))
    
    resknots<-resknots%>%dplyr::select(nknots, knot1, knot2, knot3, AIC)
  }
  
  if(length(nkl)==4){
    resknots<-resknots%>%
      mutate(dist_2_1=abs(knot2-knot1))%>%
      mutate(dist_3_2=abs(knot3-knot2))%>%
      mutate(dist_4_3=abs(knot4-knot3))
    
    #remove any row that has a distance less than mindistance
    resknots<-resknots%>%
      filter((dist_2_1==0|dist_2_1>=mindistance)&
               (dist_3_2==0|dist_3_2>=mindistance)&
               (dist_4_3==0|dist_4_3>=mindistance))
    
    resknots<-resknots%>%dplyr::select(nknots, knot1, knot2, knot3, knot4, AIC)
  }
  
  return(resknots)
}
#create flexible knot matrices
#knot_mat_9  <- create_knot_matrix(n_pers=9)
knot_mat_12 <- create_knot_matrix(n_pers=12)
knot_mat_20 <- create_knot_matrix(n_pers=20)
knot_mat_34 <- create_knot_matrix(n_pers=34)

#save flexible knot matrices
#saveRDS(knot_mat_9, "CleanData/knot_mat_9.rds")
saveRDS(knot_mat_12, "CleanData/knot_mat_12.rds")
saveRDS(knot_mat_20, "CleanData/knot_mat_20.rds")
saveRDS(knot_mat_34, "CleanData/knot_mat_34.rds")

#create 'rigid' knot matrices
knot_mat_12_rigid<-knot_mat_12%>%dplyr::filter(knot1%%3==0 & knot2%%3==0 & knot3%%3==0)
knot_mat_20_rigid<-knot_mat_20%>%dplyr::filter(knot1%%3==0 & knot2%%3==0 & knot3%%3==0)
knot_mat_34_rigid<-knot_mat_34%>%dplyr::filter(knot1%%3==0 & knot2%%3==0 & knot3%%3==0)

#save rigid knot matrices
saveRDS(knot_mat_12_rigid, "CleanData/knot_mat_12_rigid.rds")
saveRDS(knot_mat_20_rigid, "CleanData/knot_mat_20_rigid.rds")
saveRDS(knot_mat_34_rigid, "CleanData/knot_mat_34_rigid.rds")
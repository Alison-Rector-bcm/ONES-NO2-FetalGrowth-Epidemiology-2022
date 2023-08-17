#####################################
# Figure Generation
####################################
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
library(grid)
library(gridtext)
library(gridExtra) 
library(gtable)
library(ggplot2)
library(tidyverse)
##################
# Read in Data
#####################
NO2_models<-readRDS("DLNM_Data/NO2_models.rds")
mean_sd_df<-readRDS("CleanData/mean_sd_vars.rds")
#read in best knot list
knot_search_list<-readRDS("KnotSearch/FetalGrowth_knots.rds")

outcome_names<-c( "ac12",  "ac20_12","ac34_20",
                  "bpd12","bpd20_12","bpd34_20",
                  "efw12","efw20_12", "efw34_20",
                  "fl12","fl20_12","fl34_20")
names(knot_search_list)<-outcome_names
###################
# Plot Function
##################
plot_data_function<-function(mydlnm_mod_data, n_pers, outweek_name){
  dlnm_mod_data<-mydlnm_mod_data[[1]]
  sd_div_mean<-mean_sd_df$mean_sd_ratio[mean_sd_df$outcome==outweek_name]
  #meta_data<-meta_data_outcome[[1]]
  library(ggplot2)
  library(viridis)
  # get lower bound of estimated effect, estimated effect and upper bound of
  #estimated effect from the prediction model
  #multiply beta by 10 to get per 10 microgram - 
  #new extraction for plot data does not create estimates per 10
  # lowBw<-round(dlnm_mod_data$matlow[n_pers:1],4)
  # upBw<-round(dlnm_mod_data$mathigh[n_pers:1],4)
  # Bw<-round(dlnm_mod_data$matfit[n_pers:1],4)
  
  lowBw<-signif(dlnm_mod_data$matlow[n_pers:1],2)
  upBw<-signif(dlnm_mod_data$mathigh[n_pers:1],2)
  Bw<-signif(dlnm_mod_data$matfit[n_pers:1],2)
  
  lowBw_percent<-(lowBw*sd_div_mean*100)
  upBw_percent<-(upBw*sd_div_mean*100)
  Bw_percent<-(Bw*sd_div_mean*100)
  
  # create a data frame of the estimated effect and its CI and gestation weeks and
  plot.data<-data.frame(raw_lower=lowBw, raw_Beta=Bw, raw_upper=upBw, 
                        lower=lowBw_percent, Beta=Bw_percent, upper=upBw_percent,
                        weeks=seq(1:n_pers))
  
  # plot.data<-data.frame(lower=lowBw, Beta=Bw, upper=upBw,
  #                       weeks=seq(1:n_pers))
  
  
  #create indicator variable for significant weeks
  plot.data<-plot.data%>%
    #if week has adverse effect, flag as sig.week
    #upper CI is <0 or lower CI is >0 rounded to 2 decimal places
    mutate(sig.weeks=ifelse(upBw<0, .5,ifelse(lowBw>0,.75,0)))
  plot.data$sig.weeks<-as.factor(plot.data$sig.weeks)
  #create end week variable
  last_endweek<-n_pers+1
  plot.data$endweek<-c(2:last_endweek)
  
  print_data<-plot.data
  print_data$outcome=outweek_name
  write.csv(print_data, paste0( "DLNM_Data/raw_beta_", outweek_name, ".csv"))
  return(plot.data)
}

#################
# Map NO2 Plots
library(magrittr)
#get list of pred_meta data
pred_meta_list<-map(NO2_models, magrittr::extract, "pred_data_forplot")

my_n_pers_list<-c(12, 20, 34,
                  12, 20, 34, 
                  12, 20, 34,
                  12, 20, 34)

outcome_names<-c( "ac12",  "ac20_12","ac34_20",
                  "bpd12","bpd20_12","bpd34_20",
                  "efw12","efw20_12", "efw34_20",
                  "fl12","fl20_12","fl34_20")

outcomes<-c("efw", "bpd", "fl", "ac")
weeks<-c(12, 20, 34)
outweek_name_list<-sort(as.vector(outer(outcomes, weeks,paste0)))

plot_data_list<-pmap(list(mydlnm_mod_data=pred_meta_list,
                          n_pers=my_n_pers_list,
                          outweek_name=outweek_name_list),
                     plot_data_function)
names(plot_data_list)<-outcome_names
saveRDS(plot_data_list,"Plot_Data/plot_data_list.rds")

####################
# CUMULATIVE EFFECTS
####################
cumulative_effect_function<-function(plot_data_sum, dlnm_model_sum, outweek_name_sum,
                                     n_pers_sum, best.knots_sum){
  
  nlags<-n_pers_sum-1
  pred_var_forsum<-dlnm_model_sum$pred_data_sums
  #extract sd and mean ratio for outcome week
  sd_div_mean<-mean_sd_df$mean_sd_ratio[mean_sd_df$outcome==outweek_name_sum]
  #get list of significant weeks 
  #(create variable that is a 1 if sig.weeks==0.5, else 0)
  plot_data_sum$num_weeks<-rep(n_pers_sum, n_pers_sum)
  
  sig_week_data<-plot_data_sum%>%
    arrange(plot_data_sum$weeks)%>%
    mutate(J_vector=ifelse(plot_data_sum$sig.weeks==0.5,1,0))
  
  sig_week_data$lag_num<-abs(plot_data_sum$weeks-plot_data_sum$num_weeks)
  
  my_J<-sig_week_data$J_vector
  sig_weeks<-sig_week_data%>%filter(J_vector==1)%>%dplyr::select('weeks')
  #increments of exposure
  cenat=c(0,10)
  #lagfunction
  lagfun<-"ns"
  #lags 'l'
  
  l<-nlags:0
  #create onebaseis C
  best.knots_sum<-as.numeric(unlist(best.knots_sum))
  C<-onebasis(l, 'ns', knots=best.knots_sum, intercept=F)
  
  #extract 'eta' the coefficients from the crosspred 
  eta<-pred_var_forsum$coefficients
  #extract veta
  veta<-pred_var_forsum$vcov
  ################
  # # Check that correct pred has been pulled
  # #cacluated betas
  # beta<-C%*%eta
  # 
  # #create matrix of ones
  # J<-matrix(1, nrow=n_pers_sum, ncol=1)
  # S<-t(J)%*%C%*%eta*5 #note - 5 is the incremental amount of exposures
  # S
  # my_allfit<-pred_var_forsum$allfit
  # my_allfit # should match S
  # 
  # 
  # #calculate cumulative variance
  # varS<-t(J)%*%C%*%veta%*%t(C)%*%J*(5^2) #note - 5 is the incremental amount of exposures
  # seS<-sqrt(varS)
  # seS
  # my_allse<-pred_var_forsum$allse
  # my_allse # should match SeS
  
  #create vector to map betas we want to aggregate
  J_part<-my_J
  
  #calculate sum of consecutive significant betas and associated stand dev
  (beta_part_sums<-t(J_part)%*%C%*%eta*10)  #you can check here, this should match with sum(redi$fit[3:8])
  (var_part_sums<-t(J_part)%*%C%*%veta%*%t(C)%*%J_part*(10^2))
  sd_sums<-sqrt(var_part_sums)
  #get upper and lower bounds for sum betas
  low_beta_part_sums=beta_part_sums-sd_sums
  high_beta_part_sums=beta_part_sums+sd_sums
  
  lowB_sum_percent<-(low_beta_part_sums*sd_div_mean*100)
  upB_sum_percent<-(high_beta_part_sums*sd_div_mean*100)
  B_sum_percent<-(beta_part_sums*sd_div_mean*100)
  
  
  cumulative_effect_df<-data.frame(sig_week_start=min(sig_weeks$weeks),
                                   sig_week_end=max(sig_weeks$weeks),
                                   # B_sum=beta_part_sums,
                                   lowB_sum_percent=signif(lowB_sum_percent,2),
                                   upB_sum_percent=signif(upB_sum_percent,2),
                                   B_sum_percent=signif(B_sum_percent,2))
  print(outweek_name_sum)
  return(cumulative_effect_df)
  
}
#####################
#MAP CUMULATIVE BETAS
########################
rm(weeks)
#specify outweek_name_list
outcomes<-c("efw", "bpd", "fl", "ac")
my_weeks<-c(12, 20, 34)
outweek_name_list<-sort(as.vector(outer(outcomes, my_weeks,paste0)))
#specify best.knots list
library(magrittr)
best.knot_list<-map(knot_search_list, magrittr::extract, "best.knots")


#map over outcomes to make base formulas
sum_sigbeta_list<-pmap(list(plot_data_sum=plot_data_list,
                            dlnm_model_sum=NO2_models,
                            outweek_name_sum=outweek_name_list,
                            n_pers_sum=my_n_pers_list,
                            best.knots_sum=best.knot_list),
                       cumulative_effect_function)

#apply names to sum_sigbeta_list
names(sum_sigbeta_list)<-outcome_names

saveRDS(sum_sigbeta_list, "DLNM_Data/cumulative_sig_beta_NO2_newsig2.rds")
write.csv(sum_sigbeta_list,  "DLNM_Data/cumulative_sig_beta_NO2_newsig2.csv")

####################
# CUMULATIVE ADVERSE
####################
cumulative_adv_effect_function<-function(plot_data_sum, dlnm_model_sum, outweek_name_sum,
                                         n_pers_sum, best.knots_sum){
  
  nlags<-n_pers_sum-1
  pred_var_forsum<-dlnm_model_sum$pred_data_sums
  #extract sd and mean ratio for outcome week
  sd_div_mean<-mean_sd_df$mean_sd_ratio[mean_sd_df$outcome==outweek_name_sum]
  #get list of significant weeks 
  #(create variable that is a 1 if sig.weeks==0.5, else 0)
  plot_data_sum$num_weeks<-rep(n_pers_sum, n_pers_sum)
  
  sig_week_data<-plot_data_sum%>%
    arrange(plot_data_sum$weeks)%>%
    mutate(J_vector=ifelse(plot_data_sum$Beta>0&plot_data_sum$weeks>12,1,0))
  
  sig_week_data$lag_num<-abs(plot_data_sum$weeks-plot_data_sum$num_weeks)
  
  my_J<-sig_week_data$J_vector
  sig_weeks<-sig_week_data%>%filter(J_vector==1)%>%dplyr::select('weeks')
  #increments of exposure
  cenat=c(0,10)
  #lagfunction
  lagfun<-"ns"
  #lags 'l'
  
  l<-nlags:0
  #create onebaseis C
  best.knots_sum<-as.numeric(unlist(best.knots_sum))
  C<-onebasis(l, 'ns', knots=best.knots_sum, intercept=F)
  
  #extract 'eta' the coefficients from the crosspred 
  eta<-pred_var_forsum$coefficients
  #extract veta
  veta<-pred_var_forsum$vcov
  ################
  # # Check that correct pred has been pulled
  # #cacluated betas
  # beta<-C%*%eta
  # 
  # #create matrix of ones
  # J<-matrix(1, nrow=n_pers_sum, ncol=1)
  # S<-t(J)%*%C%*%eta*5 #note - 5 is the incremental amount of exposures
  # S
  # my_allfit<-pred_var_forsum$allfit
  # my_allfit # should match S
  # 
  # 
  # #calculate cumulative variance
  # varS<-t(J)%*%C%*%veta%*%t(C)%*%J*(5^2) #note - 5 is the incremental amount of exposures
  # seS<-sqrt(varS)
  # seS
  # my_allse<-pred_var_forsum$allse
  # my_allse # should match SeS
  
  #create vector to map betas we want to aggregate
  J_part<-my_J
  
  #calculate sum of consecutive significant betas and associated stand dev
  (beta_part_sums<-t(J_part)%*%C%*%eta*10)  #you can check here, this should match with sum(redi$fit[3:8])
  (var_part_sums<-t(J_part)%*%C%*%veta%*%t(C)%*%J_part*(10^2))
  sd_sums<-sqrt(var_part_sums)
  #get upper and lower bounds for sum betas
  low_beta_part_sums=beta_part_sums-sd_sums
  high_beta_part_sums=beta_part_sums+sd_sums
  
  lowB_sum_percent<-(low_beta_part_sums*sd_div_mean*100)
  upB_sum_percent<-(high_beta_part_sums*sd_div_mean*100)
  B_sum_percent<-(beta_part_sums*sd_div_mean*100)
  
  
  cumulative_effect_df<-data.frame(sig_week_start=min(sig_weeks$weeks),
                                   sig_week_end=max(sig_weeks$weeks),
                                   # beta_sum=beta_part_sums,
                                   lowB_sum_percent=round(lowB_sum_percent,4),
                                   upB_sum_percent=round(upB_sum_percent,4),
                                   B_sum_percent=round(B_sum_percent,4))
  print(outweek_name_sum)
  return(cumulative_effect_df)
  
}

###########
# MAP CUMULATIVE ADVERSE
#########################
#rm(weeks)
#specify outweek_name_list
outcomes<-c("efw", "bpd", "fl", "ac")
my_weeks<-c(12, 20, 34)
outweek_name_list<-sort(as.vector(outer(outcomes, my_weeks,paste0)))
#specify best.knots list
library(magrittr)
best.knot_list<-map(knot_search_list, magrittr::extract, "best.knots")


#map over outcomes to make base formulas
sum_advbeta_list<-pmap(list(plot_data_sum=plot_data_list,
                            dlnm_model_sum=NO2_models,
                            outweek_name_sum=outweek_name_list,
                            n_pers_sum=my_n_pers_list,
                            best.knots_sum=best.knot_list),
                       cumulative_adv_effect_function)

#apply names to sum_sigbeta_list
names(sum_advbeta_list)<-outcome_names

saveRDS(sum_advbeta_list, "DLNM_Data/cumulative_adv_beta_NO2_new.rds")
#bind rows to single list adverse
advbeta_df<-bind_rows(sum_advbeta_list, .id="Outcome")
library(kableExtra)
advbeta_df%>%kable()

write.csv(advbeta_df, "DLNM_Data/cumulative_adv_beta_NO2_table_new.csv")

##################################
# BUILD PLOT FUNCTION 
#############################
build_plots<-function(plot_data, plot.title, n_pers){
  
  
  plot.data.temp<-as.data.frame(plot_data)
  #filter data to exclude last week
  plot.data<-plot.data.temp%>%filter(weeks!=n_pers)
  library(ggplot2)
  library(viridis)
  # Specify x-axis specifications
  #specify break values based on n_pers
  mybreak_val=ifelse(n_pers==12, 2, 
                     ifelse(n_pers==20,4,6))
  #create base xseq, from 0 to n_pers by break val
  xseq<-as.vector(seq(0,n_pers, by=mybreak_val))
  xseq_34<-c(xseq,34)
  x_axis_spec=if(n_pers==12|n_pers==20){xseq}else{xseq_34}
  
  # y_axis_min<-round(round(min(plot.data$lower),0)-0.25,2)
  # y_axis_max<-round(round(max(plot.data$upper),0)+0.25,2)
  # my_y_axis_lim=max(abs(y_axis_min), abs(y_axis_max))
  # my_y_axis_lim=min(my_y_axis_lim, 15)
  #my_ylim=c(-my_y_axis_lim, my_y_axis_lim)
  #my_y_break=round(my_y_axis_lim/4,1)
  #specify y axis limits
  ymin=ifelse(n_pers==12,-5, 
              ifelse(n_pers==20,-1.5, -1))
  ymax=ifelse(n_pers==12,5,  
              ifelse(n_pers==20,1.5, 1))
  ybreak=ifelse(n_pers==12,1, 
                ifelse(n_pers==20,0.5,0.25))
  my_ylim=c(ymin, ymax)
  
  
  
  #create base plot with
  #set colors for scale_fill
  colours = c("0"= "gray90", "0.5"="violetred2", "0.75"="slateblue1")
  #get vector of significant weeks to use for vertical line intecepts
  siglines<-plot.data$sig.weeks
  #set base plot with plot.data and x and y variables
  p<-ggplot(data=plot.data, aes(x=weeks, y=Beta))+
    #use geom_rectangle to fill color on significant weeks
    # weeks option fill=sig.weeks identifies factor to use for fill colors,
    # 'alpha' sets transparency level
    # show.legend=FALSE removes legend indicating the alpha level,
    # set ymin/ymax to -/+ infinity to get rectangles to fill entire plot.
    #geom_rect(aes(xmin=weeks-.5,
    #               xmax=endweek-.5,
    #               fill=sig.weeks,
    #               alpha=0.3),
    #           show.legend=FALSE,
    #           ymin=-Inf, ymax=Inf)+
  #   # choose colors from previous setting to fill rectangles,
  #   # remove legend using "guide=FALSE"
  # scale_fill_manual(values=colours, guide=FALSE)+
  #   #set geom line for beta values
  geom_line(aes(x=weeks, y=Beta))+
    #get horizontal dotted line at y =0
    geom_hline(yintercept=0, linetype="dotted")
  #use geom_ribbon to create Confidence interval bands
  p<- p+geom_ribbon(aes(ymin=lower, ymax=upper), linetype=2, alpha=0.4)
  
  p<- p+ 
    #Assign previously set x-axis specifications and y-axis specifications
    scale_x_continuous(breaks=x_axis_spec) +
    scale_y_continuous(breaks=seq(ymin, ymax, ybreak))+
    #make panel.background blank 
    #theme(panel.background=element_blank())+#, 
    #panel.grid.major=element_blank(),
    #panel.grid.minor=element_blank())+
    #plot.background=element_rect(colour="black", fill=NA, size=3))+
    #scale_y_continuous(limits=c(-.6, .6))+
    #scale_y_continuous(limits=c(y_axis_min, y_axis_max))+
    coord_cartesian(ylim=my_ylim, xlim=c(0,n_pers+1))+
    #set x and y lables for plot and title for plot from function input
    xlab("Gestational Week")+
    #ylab("Effect of Exposure as % Change - 10 micrograms NO2")+
    #ylab(expression("%" + " " + change + " " + per + " " + mu +g ))+
    ylab(expression(paste("% Change in fetal growth")))+
    theme(axis.title.y=element_text(size=16))+
    theme(axis.title.x=element_text(size=14))+
    theme(axis.text.x=element_text(size=14))+
    theme(axis.text.y=element_text(size=14))+
    #show x/y axes solid black line
    theme(axis.line=element_line(size = 1, colour = "black", linetype = 1))+
    ggtitle(plot.title)+
    theme(title=element_text(size=12))+
    theme(plot.title = element_text(hjust = 0.5))
  current_plot<-p
  return(current_plot)
}


# plot_data<-plot.data
# plot.title<-"Test Plot - PM2.5 and EFW12"
# n_pers<-12

plot_title_list<-c("Abdominal Circumference (12 wks)", 
                   "Abdominal Circumference (12-20 wks)", 
                   "Abdominal Circumference (20-34 wks)", 
                   "Biparietal Diameter (12 wks)", 
                   "Biparietal Diameter (12-20 wks)", 
                   "Biparietal Diameter (20-34 wks)", 
                   "Estimated Fetal Weight (12 wks)",
                   "Estimated Fetal Weight (12-20 wks)", 
                   "Estimated Fetal Weight (20-34wks",
                   "Femure Length (12 wks)", 
                   "Femure Length (12-20 wks)", 
                   "Femure Length (20-34wks)")
my_n_pers_list<-c(12, 20, 34, 
                  12, 20, 34, 
                  12, 20, 34, 
                  12, 20, 34)

plot_list<-pmap(list(plot_data=plot_data_list,
                     plot.title=plot_title_list,
                     n_pers=my_n_pers_list),
                build_plots)

######################
# CREATE PLOT GRIDS
###################

#Build table of three plots for same outcome at different weeks
#reads in list of data to iput into plot function
#outcome should be character: example "efw"


# myplot_data_list<-list(plot_data_list$efw12,
#                    plot_data_list$bpd12,
#                    plot_data_list$ac20,
#                    plot_data_list$fl12)
#   title_list=c("Estimated Fetal Weight", 
#                "Biparietal Diameter",
#                "Abdominal Circumference",
#                "Femur Length"
#                )
#   week_list=c(12, 12, 12)
#   my_n_pers=12
#   outcome="A. Early Pregnancy Growth"

plot_table_row<-function(myplot_data_list, 
                         title_list, 
                         my_n_pers,
                         outcome){
  
  n<-length(myplot_data_list)
  plot_list<-list()
  for(i in 1:n){
    curr_plot<-build_plots(plot_data=as.data.frame(myplot_data_list[[i]]),
                           plot.title=title_list[i],
                           n_pers=my_n_pers)
    # n_pers=week_list[i])
    
    #remove axis titles from each plot
    curr_plot<-curr_plot+
      theme(axis.title.y=element_blank())+
      theme(plot.title = element_text(hjust = 0.5))+
      #theme(plot.margin = unit(c(1, 1, 1.2, 1), "cm"))+
      theme( plot.background = element_rect(
        fill = NULL,
        colour = "black",
        size = 1.2
      ))
    
    plot_list[[i]]<-curr_plot
  }
  #TODO: double check how to change axes titles
  #TODO: confirm saving
  #Single y -axis title and tick marks on far left plot;
  
  #x-axis each one gets tick mark
  # only one x-axis title across all plots, 
  # title for each plot or caption for each plot indicating outcome/condition
  
  ## Add vertical (90 deg rotate) labels for 'outcome' -> "Early Pregnancy Growth", "Mid Pregnancy Growth"
  #          "Late Pregnancy Growth"
  myoutcome_grob=tableGrob(paste(" \n  "), theme = ttheme_minimal(base_size = 14))
  #my_panel_grobv2=expression(bold(underline("My Sample Header")))
  # 
  # if(outcome=="Early_Pregnancy"){
  #   my_panel_grob=textGrob(paste0(expression(bold("A. Early Pregnancy")),"\n% Change in fetal growth"))}
  # if(outcome=="Mid_Pregnancy"){
  #   my_panel_grob=paste0(expression(bold("B.  PregMidnancy")),"\n% Change in fetal growth")}
  # if(outcome=="Late_Pregnancy"){
  #     my_panel_grob=paste0(expression(bold("C. Late Pregnancy")),"\n% Change in fetal growth")}
  #my_panel_grob=(paste0(expression(bold(outcome)),"\n% Change in fetal growth"))
  
  if(outcome=="Early_Pregnancy"){
    my_grid<-   grid.arrange(grobs=plot_list, 
                             ncol=4,
                             #  top = textGrob(outcome, vjust=1,gp=gpar(fontface="bold", fontsize=14,cex=1.0, vjust=0)))
                             top=richtext_grob("**Early Pregnancy Growth**<br>", vjust=.5),
                             left=richtext_grob("<br>% Change in fetal growth", rot=90, vjust=.5, gp=gpar(fontface="plain", cex=1.0)),
                             # left=textGrob("Test\n Test 2", rot=90, vjust=0.5),
                             bottom=textGrob(" ",#"Gestational Week", 
                                             vjust=0.2, gp=gpar(fontface="plain", cex=1.2)))
    
  }
  
  
  if(outcome=="Mid_Pregnancy"){
    my_grid<-   grid.arrange(grobs=plot_list, 
                             ncol=4,
                             #  top = textGrob(outcome, vjust=1,gp=gpar(fontface="bold", fontsize=14,cex=1.0, vjust=0)))
                             top=richtext_grob("**Mid Pregnancy Growth**<br>", vjust=.5),
                             left=richtext_grob("<br>% Change in fetal growth", rot=90, vjust=.5, gp=gpar(fontface="plain", cex=1.0)),
                             # left=textGrob("Test\n Test 2", rot=90, vjust=0.5),
                             bottom=textGrob(" ",#"Gestational Week", 
                                             vjust=0.2, gp=gpar(fontface="plain", cex=1.2)))
    
  }
  
  
  if(outcome=="Late_Pregnancy"){
    my_grid<-   grid.arrange(grobs=plot_list, 
                             ncol=4,
                             #  top = textGrob(outcome, vjust=1,gp=gpar(fontface="bold", fontsize=14,cex=1.0, vjust=0)))
                             top=richtext_grob("**Late Pregnancy Growth**<br>", vjust=.5),
                             left=richtext_grob("<br>% Change in fetal growth", rot=90, vjust=.5, gp=gpar(fontface="plain", cex=1.0)),
                             # left=textGrob("Test\n Test 2", rot=90, vjust=0.5),
                             bottom=textGrob(" ",#"Gestational Week", 
                                             vjust=0.2, gp=gpar(fontface="plain", cex=1.2)))
    
  }
  my_grid_table<-grid.arrange(my_grid, myoutcome_grob,  ncol=1, heights=c(4,.2))
  # ggsave(my_grid_table,
  #       file=paste0("Figures/plot_table_", outcome, "new_v1.svg"),
  #       width=36,
  #       height=10,
  #       units='cm')
  ggsave(my_grid,
         file=paste0("Figures/plot_table_", outcome, ".svg"),
         width=36,
         height=10,
         units='cm')
  ggsave(my_grid,
         file=paste0("Figures/plot_table_", outcome, ".png"),
         width=36,
         height=10,
         units='cm')
  return(my_grid_table)
}

######################################
# PLOT GRIDS BY GESTATIONAL PERIOD
#################################
early_plot_table<-plot_table_row(
  myplot_data_list=list(plot_data_list$efw12,
                        plot_data_list$bpd12,
                        plot_data_list$ac12,
                        plot_data_list$fl12),
  title_list=c("Estimated Fetal Weight", 
               "Biparietal Diameter",
               "Abdominal Circumference",
               "Femur Length"
  ),
  my_n_pers=12,
  outcome="Early_Pregnancy")



mid_plot_table<-plot_table_row(
  myplot_data_list=list(plot_data_list$efw20_12,
                        plot_data_list$bpd20_12,
                        plot_data_list$ac20_12,
                        plot_data_list$fl20_12),
  title_list=c("Estimated Fetal Weight", 
               "Biparietal Diameter",
               "Abdominal Circumference",
               "Femur Length"
  ),
  my_n_pers=20,
  outcome="Mid_Pregnancy")

late_plot_table<-plot_table_row(
  myplot_data_list=list(plot_data_list$efw34_20,
                        plot_data_list$bpd34_20,
                        plot_data_list$ac34_20,
                        plot_data_list$fl34_20),
  title_list=c("Estimated Fetal Weight", 
               "Biparietal Diameter",
               "Abdominal Circumference",
               "Femur Length"
  ),
  my_n_pers=34,
  outcome="Late_Pregnancy")


my_grid_time<-arrangeGrob(early_plot_table, mid_plot_table, late_plot_table, ncol=1 )

ggsave(my_grid_time,
       file=paste0("Figures/plot_table_allperiods.svg"),
       width=36,
       height=32,
       units='cm')
ggsave(my_grid_time,
       file=paste0("Figures/plot_table_allperiods.png"),
       width=36,
       height=32,
       units='cm')

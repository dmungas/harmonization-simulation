library(tidyverse)
require(readr)
require(MplusAutomation)
require(mirt)
require(tidyr)
require(ggplot2)
require(lubridate)
library(DT)
library(RColorBrewer)
library(ggthemes)
library(officer)
#source("~/Research/Code/runMplus.R")
#source("~/Research/Code/extractVarNames.R")
source("Code/simulation_functions.R")


# *********************** mirt calibration of PITCH TICS ***********************

### mirt calibration of tics

# input PITCH tics dataset
tics <- read.csv("Data/PITCH-response-data-fhl.csv")

# HRS items
vars <- c('UBAK','UDAT','UDAY','UDWR','UIWR','UMON','UNM1','UNM2','UNM5','UNM6',
          'USUB','UYER')
# MHAS items
varsm <- c("UDAY","UFCO2","UFRE1","UMON","UVSC","UWD","UWR1","UWR2","UWR3","UYER")

# select HRS items and persons
hrs <- tics[tics$study_name_short %in% c("HRS_CODA_W6","HRS_W6"),c("newid",
                                                                   "study_wave_number","study_name_short",vars)]
hrs$n_itm <- apply(hrs[,vars],1,function(x) sum(!is.na(x)))
hrs <- hrs[!hrs$n_itm ==0,]

m1hrs <- mirt.model('cog = 1-12') # creates mirt model object
hrs_par <- mirt(hrs[,vars],m1hrs,pars='values') # generates item parameters file that can be edited to guide further analyses

hrscal <- mirt(hrs[,vars],m1hrs,pars=hrs_par) # IRT calibration
# coef(hrscal)
# coef(hrscal, IRTpars = TRUE)

m2hrs <- mirt.model('cog = 1-11')

# HRS calibration excluding delayed word recall
hrs_par_0dr <- mirt(hrs[,vars[!vars %in% "UDWR"]],m2hrs,pars='values')
hrscal_0dr <- mirt(hrs[,vars[!vars %in% "UDWR"]],m2hrs,pars=hrs_par_0dr)



# HRS calibration excluding immediate word recall
hrs_par_0ir <- mirt(hrs[,vars[!vars %in% "UIWR"]],m2hrs,pars='values',
                    technical=list(removeEmptyRows=TRUE))
hrscal_0ir <- mirt(hrs[,vars[!vars %in% "UIWR"]],m2hrs,pars=hrs_par_0ir,
                   technical=list(removeEmptyRows=TRUE))
# coef(hrscal_0ir)

#  mirt calibration of MHAS 

mex <- tics[tics$study_name_short %in% c("MHAS_W1","MHAS_W2"),c("newid",
                                                                "study_wave_number","study_name_short",varsm)]
mex$n_itm <- apply(mex[,varsm],1,function(x) sum(!is.na(x)))
mex <- mex[!mex$n_itm ==0,]

m1mex <- mirt.model('cog = 1-10')
mex_par <- mirt(mex[,varsm],m1mex,pars='values')

mexcal <- mirt(mex[,varsm],m1mex,pars=mex_par)
# coef(mexcal)

#  mirt calibration of combined HRS and MHAS 

hrme <- tics[tics$study_name_short %in% c("HRS_CODA_W6","HRS_W6") |
               tics$study_name_short %in% c("MHAS_W1","MHAS_W2"),c("newid",
                                                                   "study_wave_number","study_name_short",union(vars,varsm))]
hrme$n_itm <- apply(hrme[,union(vars,varsm)],1,function(x) sum(!is.na(x)))
hrme <- hrme[!hrme$n_itm == 0,]

m1hrme <- mirt.model('cog = 1-19')
hrme_par <- mirt(hrme[,union(vars,varsm)],m1hrme,pars='values')


hrmecal <- mirt(data = hrme[,union(vars,varsm)], model = m1hrme, pars = hrme_par)
# hrmecal <- mirt(hrme[,union(vars,varsm)],m1hrme)
# coef(hrmecal)
# parameters(hrmecal)

# coef(hrmecal, IRTpars = TRUE, simplify = TRUE)$items %>%
#   DT::datatable(options = list(pageLength = 19), rownames = TRUE) %>%
#   formatRound(columns = 1:13, digits = 3)

fsc <- fscores(hrmecal,method="EAP")
hrme$ability_est <- fsc

# table(hrme$study_name_short)
hrme$study <- ifelse(grepl("HRS",hrme$study_name_short),"HRS","MHAS")

mean_hrs <- mean(hrme[hrme$study %in% "HRS","ability_est"])
sd_hrs <- sd(hrme[hrme$study %in% "HRS","ability_est"])
mean_mhas <- mean(hrme[hrme$study %in% "MHAS","ability_est"])
sd_mhas <- sd(hrme[hrme$study %in% "MHAS","ability_est"])
mean_sd <- list(mean_hrs,mean_mhas,sd_hrs,sd_mhas)
mean_all <- mean(hrme[hrme$study %in% c("HRS","MHAS"),"ability_est"])
sd_all <- sd(hrme[hrme$study %in% c("HRS","MHAS"),"ability_est"])

mean <- c(mean_hrs,mean_mhas,mean_all)
sd <- c(sd_hrs,sd_mhas,sd_all)

mean_sd <- data.frame(mean,sd)
row.names(mean_sd) <- c("Group 1","Group 2","Combined")

# ----------------------------end mirt calibration -----------------------------

# *************************** Simulation Scenarios *****************************


# Set all scenarios' arguments for consistency

Seed <- 21589
mus <- c(0, -0.24)
sigmas <- c(sd_hrs, sd_mhas)
refGrp <- 1
repTheta <- 1
reps <- 500
fsMeth <- "EAP"
n1 <- 500
n2 <- 500

### scenario 1 - all items shared

items <- union(vars,varsm)

scen_1 <- equateSim(seed = Seed, 
                    grp_mean = mus,
                    grp_sd = sigmas,
                    ref_grp = refGrp,
                    n_rep_theta = repTheta,
                    #n_itm=length(items),
                    pars=hrme_par,
                    n_rep = reps,
                    itms1=items,
                    itms2=items,
                    fsc_method = fsMeth,
                    mod_res_obj=hrmecal,
                    n_samp1 = n1,
                    n_samp2 = n2,
                    save_sims = TRUE,
                    save_log = TRUE,
                    log_file = "Logs/scenario1.txt",
                    verbose = FALSE, 
                    prog_bar = TRUE)

sumstat1 <- summaryStats(scen_1[["summary"]])
longstat1 <- scen_1[["summary"]]
ds1 <- scen_1$datasets


### end scenario 1

### scenario 2 - UMON. UDAY, UYER as linking items, actual items in HRS and MHAS

scen_2 <- equateSim(seed = Seed,
                    grp_mean = mus,
                    grp_sd = sigmas,
                    ref_grp = refGrp,
                    n_rep_theta = repTheta,
                    #n_itm=length(items),
                    pars=hrme_par,
                    n_rep = reps,
                    itms1=vars,
                    itms2=varsm,
                    fsc_method = fsMeth,
                    mod_res_obj=hrmecal,
                    n_samp1 = n1,
                    n_samp2 = n2,
                    save_sims = TRUE,
                    save_log = TRUE,
                    log_file = "Logs/scenario2.txt",
                    verbose = FALSE, 
                    prog_bar = TRUE)

sumstat2 <- summaryStats(scen_2[["summary"]])
longstat2 <- scen_2[["summary"]]
ds2 <- scen_2$datasets

subset(longstat2, group == 0 & type == "raw" & statistic == "est_mean") %>%
  select(avg) %>%
  unlist() %>%
  sd()

### end scenario 2

### scenario 3 - UMON. UDAY, UYER UIWR as linking items

itms2 <- c(varsm,"UIWR")

scen_3 <- equateSim(seed = Seed,
                    grp_mean = mus,
                    grp_sd = sigmas,
                    ref_grp = refGrp,
                    n_rep_theta = repTheta,
                    #n_itm=length(items),
                    pars=hrme_par,
                    n_rep = reps,
                    itms1=vars,
                    itms2=itms2,
                    fsc_method = fsMeth,
                    mod_res_obj=hrmecal,
                    n_samp1 = n1,
                    n_samp2 = n2,
                    save_sims = TRUE,
                    save_log = TRUE,
                    log_file = "Logs/scenario3.txt",
                    verbose = FALSE, 
                    prog_bar = TRUE)

sumstat3 <- summaryStats(scen_3[["summary"]])
longstat3 <- scen_3[["summary"]]
ds3 <- scen_3$datasets

### end scenario 3

### scenario 4 - HRS Items in both studies


scen_4 <- equateSim(seed = Seed,
                    grp_mean = mus,
                    grp_sd = sigmas,
                    ref_grp = refGrp,
                    n_rep_theta = repTheta,
                    #n_itm=length(vars),
                    pars=hrs_par,
                    n_rep = reps,
                    itms1=vars,
                    itms2=vars,
                    fsc_method = fsMeth,
                    mod_res_obj=hrscal,
                    n_samp1 = n1,
                    n_samp2 = n2,
                    save_sims = TRUE,
                    save_log = TRUE,
                    log_file = "Logs/scenario4.txt",
                    verbose = FALSE, 
                    prog_bar = TRUE)

sumstat4 <- summaryStats(scen_4[["summary"]])
longstat4 <- scen_4[["summary"]]
ds4 <- scen_4$datasets

### end scenario 4

### scenario 5 - MHAS Items in both studies


scen_5 <- equateSim(seed = Seed,
                    grp_mean = mus,
                    grp_sd = sigmas,
                    ref_grp = refGrp,
                    n_rep_theta = repTheta,
                    #n_itm=length(varsm),
                    pars=mex_par,
                    n_rep = reps,
                    itms1=varsm,
                    itms2=varsm,
                    fsc_method = fsMeth,
                    mod_res_obj=mexcal,
                    n_samp1 = n1,
                    n_samp2 = n2,
                    save_sims = TRUE,
                    save_log = TRUE,
                    log_file = "Logs/scenario5.txt",
                    verbose = FALSE, 
                    prog_bar = TRUE)

sumstat5 <- summaryStats(scen_5[["summary"]])
longstat5 <- scen_5[["summary"]]
ds5 <- scen_5$datasets

#scen_5_data <- scen_5[["datasets"]]


### end scenario 5

# dat1 <- scen_5_data[scen_5_data$sample == 1,]
# 
# plot(dat1$resid ~ dat1$ability)
# abline(lm(dat1$resid ~ dat1$ability))
# 
# ggplot(data=scen_5_data[scen_5_data$sample == 8,], aes(ability,resid)) + 
#   geom_point() + geom_smooth(method = "loess")


### scenario 6 -  UDAY as linking item


vars6 <- vars[!vars %in% c("UMON",'UYER')]
varsm6 <- varsm[!varsm %in% c("UMON",'UYER')]

scen_6 <- equateSim(seed = Seed,
                    grp_mean = mus,
                    grp_sd = sigmas,
                    ref_grp = refGrp,
                    n_rep_theta = repTheta,
                    #n_itm=length(items),
                    pars=hrme_par,
                    n_rep = reps,
                    itms1=vars,
                    itms2=varsm6,
                    fsc_method = fsMeth,
                    mod_res_obj=hrmecal,
                    n_samp1 = n1,
                    n_samp2 = n2,
                    save_sims = TRUE,
                    save_log = TRUE,
                    log_file = "Logs/scenario6.txt",
                    verbose = FALSE, 
                    prog_bar = TRUE)

sumstat6 <- summaryStats(scen_6[["summary"]])
longstat6 <- scen_6[["summary"]]
ds6 <- scen_6$datasets

### end scenario 6

### scenario 7 - UIWR as linking item, MHAS without UMON UDAY UYER

varsm7 <- c(varsm[!varsm %in% c("UMON",'UYER','UDAY')],'UIWR')


scen_7 <- equateSim(seed = Seed,
                    grp_mean = mus,
                    grp_sd = sigmas,
                    ref_grp = refGrp,
                    n_rep_theta = repTheta,
                    #n_itm=length(items),
                    pars=hrme_par,
                    n_rep = reps,
                    itms1=vars,
                    itms2=varsm7,
                    fsc_method = fsMeth,
                    mod_res_obj=hrmecal,
                    n_samp1 = n1,
                    n_samp2 = n2,
                    save_sims = TRUE,
                    save_log = TRUE,
                    log_file = "Logs/scenario7.txt",
                    verbose = FALSE, 
                    prog_bar = TRUE)

# scen_8 <- equateSim(seed = Seed,grp_mean = mus,
#                     grp_sd = sigmas,ref_grp = refGrp,n_rep_theta = repTheta,n_itm=length(items),
#                     pars=hrme_par,n_rep = reps00,itms1=vars8,itms2=varsm8,fsc_method = fsMeth,
#                     mod_res_obj=hrmecal,n_samp1 = n1,n_samp2 = n2)

sumstat7 <- summaryStats(scen_7[["summary"]])
longstat7 <- scen_7[["summary"]]
ds7 <- scen_7$datasets

### end scenario 7

### scenario 8 - UIWR as linking item, HRS without UMON UDAY UYER

vars8 <- vars[!vars %in% c("UMON",'UYER','UDAY')]
varsm8 <- c(varsm,'UIWR')


scen_8 <- equateSim(seed = Seed,
                    grp_mean = mus,
                    grp_sd = sigmas,
                    ref_grp = refGrp,
                    n_rep_theta = repTheta,
                    #n_itm=length(items),
                    pars=hrme_par,
                    n_rep = reps,
                    itms1=vars8,
                    itms2=varsm8,
                    fsc_method = fsMeth,
                    mod_res_obj=hrmecal,
                    n_samp1 = n1,
                    n_samp2 = n2,
                    save_sims = TRUE,
                    save_log = TRUE,
                    log_file = "Logs/scenario8.txt",
                    verbose = FALSE, 
                    prog_bar = TRUE)

# seed <- 21589
# grp_mean <- c(0,-0.24)
# grp_sd <- c(1,1.1)
# ref_grp <- 1
# n_rep_theta <- 1
# n_itm <- length(items)
# pars <- hrme_par
# n_rep <- 500
# itms1 <- vars8
# itms2 <- varsm8
# fsc_method <- "EAP"
# mod_res_obj <- hrmecal
# n_samp1 <- 500
# n_samp2 <- 500

# scen_8 <- equateSim(seed = Seed,grp_mean=c(mean_hrs,mean_mhas),
#                     grp_sd = sigmas,n_samp=500,n_rep_theta = repTheta,n_itm=length(items),
#                     pars=hrme_par,n_rep = reps00,itms1=vars8,itms2=varsm8,fsc_method = fsMeth,
#                     mod_res_obj=hrmecal)

sumstat8 <- summaryStats(scen_8[["summary"]])
longstat8 <- scen_8[["summary"]]
ds8 <- scen_8$datasets

scen_9 <- equateSim(seed = Seed, 
                    grp_mean = mus,
                    grp_sd = sigmas,
                    ref_grp = refGrp,
                    n_rep_theta = repTheta,
                    #n_itm=length(items),
                    pars=hrme_par,
                    n_rep = reps,
                    itms1=items,
                    itms2=items,
                    fsc_method = fsMeth,
                    mod_res_obj=hrmecal,
                    n_samp1 = n1,
                    n_samp2 = n2,
                    save_sims = TRUE,
                    save_log = TRUE,
                    log_file = "Logs/scenario9.txt",
                    verbose = FALSE, 
                    prog_bar = TRUE,
                    man_override = list(iname = c("UDAY", "UMON"),
                                        parameter = c("b", "b"), 
                                        val = c(0, 1.928)))

sumstat9 <- summaryStats(scen_9[["summary"]])
longstat9 <- scen_9[["summary"]]
ds9 <- scen_9$datasets



### Merge summary statistics from scenarios

sumstat1$scenario <- 1
sumstat2$scenario <- 2
sumstat3$scenario <- 3
sumstat4$scenario <- 4
sumstat5$scenario <- 5
sumstat6$scenario <- 6

# sumstat <- rbind(sumstat1,sumstat2,sumstat3,sumstat4,sumstat5,sumstat6)
# sumstat$scenario_label <- ifelse(sumstat$scenario == 1,"HRS+MHAS_all_shared",
#                                  ifelse(sumstat$scenario == 2,"HRS+MHAS_UMON_UDAY_UYER_shared",
#                                         ifelse(sumstat$scenario == 3,"HRS+MHAS_UMON_UDAY_UYER_UIWR_shared",
#                                                ifelse(sumstat$scenario == 4,"HRS_all_shared",
#                                                       ifelse(sumstat$scenario == 5,"MHAS_all_shared",
#                                                              ifelse(sumstat$scenario == 6,"HRS+HMAS_UDAY_shared",
#                                                                     ifelse(sumstat$scenario == 7,"HRS+MHAS_UIWR_shared_no_MHAS_dates",
#                                                                            ifelse(sumstat$scenario == 7,"HRS+MHAS_UIWR_shared_no_HRS_dates",
#                                                                                   NA))))))))

#load("Results/simulation_results_2020-01-07-18-28.RData")

sumstat7$scenario <- 7
sumstat8$scenario <- 8
sumstat9$scenario <- 9
#sumstat10$scenario <- 10

longstat1$scenario <- 1
longstat2$scenario <- 2
longstat3$scenario <- 3
longstat4$scenario <- 4
longstat5$scenario <- 5
longstat6$scenario <- 6
longstat7$scenario <- 7
longstat8$scenario <- 8
longstat9$scenario <- 9
#longstat10$scenario <- 10

longstat <- bind_rows(longstat1,
                      longstat2,
                      longstat3,
                      longstat4,
                      longstat5,
                      longstat6,
                      longstat7,
                      longstat8,
                      longstat9
                      #,longstat10
                      )

ds1$scenario <- 1
ds2$scenario <- 2
ds3$scenario <- 3
ds4$scenario <- 4
ds5$scenario <- 5
ds6$scenario <- 6
ds7$scenario <- 7
ds8$scenario <- 8
ds9$scenario <- 9
#ds10$scenario <- 10

dsAll <- bind_rows(ds1, ds2, ds3, ds4, ds5, ds6, ds7, ds8, ds9
                   #, ds10
                   )

# sumstat7 <- rbind(sumstat7,sumstat8)
# sumstat7$scenario_label <- ifelse(sumstat7$scenario == 1,"HRS+MHAS_all_shared",
#                                   ifelse(sumstat7$scenario == 2,"HRS+MHAS_UMON_UDAY_UYER_shared",
#                                          ifelse(sumstat7$scenario == 3,"HRS+MHAS_UMON_UDAY_UYER_UIWR_shared",
#                                                 ifelse(sumstat7$scenario == 4,"HRS_all_shared",
#                                                        ifelse(sumstat7$scenario == 5,"MHAS_all_shared",
#                                                               ifelse(sumstat7$scenario == 6,"HRS+HMAS_UDAY_shared",
#                                                                      ifelse(sumstat7$scenario == 7,"HRS+MHAS_UIWR_shared_no_MHAS_dates",
#                                                                             ifelse(sumstat7$scenario == 8,"HRS+MHAS_UIWR_shared_no_HRS_dates",
#                                                                                    NA))))))))

# sumstat <- rbind(sumstat,sumstat7)
sumstat <- rbind(sumstat1,sumstat2,sumstat3,sumstat4,sumstat5,sumstat6, sumstat7, sumstat8, sumstat9
                 #, sumstat10
                 )
sumstat$scenario_label <- factor(sumstat$scenario, 
                                 levels = 1:9, 
                                 labels = c("HRS+MHAS_all_shared",
                                            "HRS+MHAS_UMON_UDAY_UYER_shared",
                                            "HRS+MHAS_UMON_UDAY_UYER_UIWR_shared",
                                            "HRS_all_shared", 
                                            "MHAS_all_shared",
                                            "HRS+HMAS_UDAY_shared",
                                            "HRS+MHAS_UIWR_shared_no_MHAS_dates",
                                            "HRS+MHAS_UIWR_shared_no_HRS_dates",
                                            "HRS+MHAS_UMON(b=1.928)_UDAY(b=0)_UYER(b=-1.928)_shared"))

#save(sumstat,mean_sd,file=paste0("Results/simulation_results_",format(Sys.time(), '%Y-%m-%d-%H-%M'),".RData"))
saveRDS(sumstat,file=paste0("Results/sumstat_",format(Sys.time(), '%Y-%m-%d-%H-%M'),".rds"))
#sumstat <- readRDS("Results/sumstat_2020-07-08-18-13.rds")

save(file=paste0("Results/simulation_results_",format(Sys.time(), '%Y-%m-%d-%H-%M'),".RData"))
#load()


DT::datatable(sumstat, options = list(pageLength = 30), rownames = FALSE) %>%
  formatRound(columns = names(sumstat), digits = 3) %>%
  formatRound(columns = "scenario", digits = 0)

read_docx() %>%  # a new, empty document
  body_add_table(sumstat %>% 
                   mutate_if(is.numeric, round, 3) %>%
                   dplyr::select(scenario, mean, sd, bias, ese, rmse, r_theta_est, mean_st, sd_st, bias_st, ese_st, rmse_st) %>%
                   filter(group == "Combined"), style = "table_template") %>% 
  print(target="Results/Table3.docx")

read_docx() %>%  # a new, empty document
  body_add_table(sumstat %>% 
                   mutate_if(is.numeric, round, 3) %>%
                   dplyr::select(scenario, mean, sd, bias, ese, rmse, r_theta_est, mean_st, sd_st, bias_st, ese_st, rmse_st) %>%
                   filter(group == "Group 1"), style = "table_template") %>% 
  print(target="Results/Table4.docx")

read_docx() %>%  # a new, empty document
  body_add_table(sumstat %>% 
                   mutate_if(is.numeric, round, 3) %>%
                   dplyr::select(scenario, mean, sd, bias, ese, rmse, r_theta_est, mean_st, sd_st, bias_st, ese_st, rmse_st) %>%
                   filter(group == "Group 2"), style = "table_template") %>% 
  print(target="Results/Table5.docx")

# Plot scenarios

sumlong <- sumstat %>%
  dplyr::select(group, mean, mean_st, sd, sd_st, n_rep, scenario) %>%
  pivot_longer(cols = c(mean, mean_st, sd, sd_st),
               names_to = c(".value", "stat"),
               names_sep = c("_")) %>%
  mutate_at(vars(stat), replace_na, "us") %>%
  mutate(low95 = mean - qnorm(.975)*sd/sqrt(n_rep),
         up95 = mean + qnorm(.975)*sd/sqrt(n_rep)) %>%
  mutate(scenF = factor(scenario, levels = 1:9, labels = paste0("Scenario ", 1:9))) %>%
  mutate(Group = factor(group, levels = c("Combined", "Group 1", "Group 2"))) %>%
  mutate(stat = factor(stat, levels = c("us", "st"),
                       labels = c("Unstandardized", "Standardized"))) %>%
  mutate(avg = mean)

thetalines3 <- data.frame(group = rep(c("Group 1", "Group 2", "Combined"), 2),
                         name = rep(c("mean", "mean_st"), each = 3),
                         value = rep(c(0, -.24, -.12), 2)) %>%
  mutate(Group = factor(group, levels = c("Combined", "Group 1", "Group 2")))

thetalines2 <- data.frame(group = rep(c("Group 2", "Combined"), 2),
                         name = rep(c("mean", "mean_st"), each = 2), 
                         value = rep(c(-.24, -.12), 2)) %>%
  mutate(Group = factor(group, levels = c("Combined", "Group 2")))

# Bar plots with 95% CIs (doesn't look that nice)
# sumlong %>% 
#   ggplot(aes(x = scenario, y = mean)) + 
#   geom_col() +
#   facet_grid(Group ~ stat) +
#   geom_hline(data = thetalines, aes(yintercept = value), 
#              lty = 2) +
#   ylab("Mean Ability") +
#   xlab("Scenario") +
#   geom_errorbar(aes(ymin =low95, ymax = up95)) +
#   ylim(-.5, .25)

# Raincloud plots
if(!exists("geom_flat_violin")) {
  source("https://gist.githubusercontent.com/benmarwick/2a1bb0133ff568cbe28d/raw/fb53bd97121f7f9ce947837ef1a4c65a73bffb3f/geom_flat_violin.R")
}

for(s in unique(longstat$type)){
  longstat %>%
    filter(type == s) %>%
    filter(statistic == "est_mean") %>%
    filter(group != 1) %>%
    mutate(scenF = factor(scenario, levels = 1:9, labels = paste0("Scenario ", 1:9))) %>%
    mutate(Group = factor(group, levels = c(0, 2), labels = c("Combined", "Group 2"))) %>%
    ggplot(aes(x = fct_rev(scenF), y = avg, fill = scenF)) +
    geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8, scale = "width") +
    geom_point(aes(y = avg, color = scenF), 
               position = position_jitter(width = .15), 
               size = .5, alpha = 0.8) +
    geom_boxplot(width = .1, outlier.shape = NA, alpha = 0.85, fill = "white", notch = TRUE) +
    ylab("Mean Ability") +
    guides(fill = FALSE) +
    guides(color = FALSE) +
    coord_flip() +
    expand_limits(x = max(longstat$scenario) + 1) +
    xlab("") +
    facet_wrap(~Group, ncol = 2) +
    scale_colour_brewer(name = "Scenario", palette = "Set1") +
    scale_fill_brewer(name = "Scenario", palette = "Set1") +
    geom_hline(data = thetalines2, aes(yintercept = value), 
               lty = 2) +
    theme_few() + 
    theme(panel.background = element_rect(fill = "gray85"),
          axis.text.x = element_text(size = 6))
  
  ggsave(paste0("Plots/Figure_", s, "_2a.tiff"), width = 6.5, height = 4.5) # High quality for submission
  ggsave(paste0("Plots/Figure_", s, "_2a.png"), width = 6.5, height = 4.5) # For Word doc
  
  longstat %>%
    filter(type == s) %>%
    filter(statistic == "est_mean") %>%
    mutate(scenF = factor(scenario, levels = 1:9, labels = paste0("Scenario ", 1:9))) %>%
    mutate(Group = factor(group, levels = 0:2, labels = c("Combined", "Group 1", "Group 2"))) %>%
    ggplot(aes(x = fct_rev(scenF), y = avg, fill = scenF)) +
    geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8, scale = "width") +
    geom_point(aes(y = avg, color = scenF), 
               position = position_jitter(width = .15), 
               size = .5, alpha = 0.8) +
    geom_boxplot(width = .1, outlier.shape = NA, alpha = 0.85, fill = "white", notch = TRUE) +
    ylab("Mean Ability") +
    guides(fill = FALSE) +
    guides(color = FALSE) +
    coord_flip() +
    expand_limits(x = max(longstat$scenario) + 1) +
    xlab("") +
    facet_wrap(~Group, ncol = 3) +
    scale_colour_brewer(name = "Scenario", palette = "Set1") +
    scale_fill_brewer(name = "Scenario", palette = "Set1") +
    geom_hline(data = thetalines3, aes(yintercept = value), 
               lty = 2) +
    theme_few() + 
    theme(panel.background = element_rect(fill = "gray85"),
          axis.text.x = element_text(size = 6))
  
  ggsave(paste0("Plots/Figure_", s, "_3a.tiff"), width = 6.5, height = 4.5) # High quality for submission
  ggsave(paste0("Plots/Figure_", s, "_3a.png"), width = 6.5, height = 4.5) # For Word doc
  
  longstat %>%
    filter(type == s) %>%
    filter(statistic == "est_mean") %>%
    filter(group != 1) %>%
    mutate(scenF = factor(scenario, levels = 1:9, labels = paste0("Scenario ", 1:9))) %>%
    mutate(Group = factor(group, levels = c(0, 2), labels = c("Combined", "Group 2"))) %>%
    ggplot(aes(x = fct_rev(scenF), y = avg, fill = scenF)) +
    geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8, scale = "width") +
    geom_point(aes(y = avg, color = scenF), 
               position = position_jitter(width = .15), 
               size = .5, alpha = 0.8) +
    geom_boxplot(width = .1, outlier.shape = NA, alpha = 0.85, fill = "white", notch = TRUE) +
    ylab("Mean Ability") +
    guides(fill = FALSE) +
    guides(color = FALSE) +
    coord_flip() +
    expand_limits(x = max(longstat$scenario) + 1) +
    xlab("") +
    facet_wrap(~Group, ncol = 2) +
    scale_colour_brewer(name = "Scenario", palette = "Set1") +
    scale_fill_brewer(name = "Scenario", palette = "Set1") +
    theme_few() + 
    theme(panel.background = element_rect(fill = "gray85"),
          axis.text.x = element_text(size = 6))
  
  ggsave(paste0("Plots/Figure_", s, "_2b.tiff"), width = 6.5, height = 4.5) # High quality for submission
  ggsave(paste0("Plots/Figure_", s, "_2b.png"), width = 6.5, height = 4.5) # For Word doc
  
  longstat %>%
    filter(type == s) %>%
    filter(statistic == "est_mean") %>%
    mutate(scenF = factor(scenario, levels = unique(longstat$scenario), labels = paste0("Scenario ", unique(longstat$scenario)))) %>%
    mutate(Group = factor(group, levels = 0:2, labels = c("Combined", "Group 1", "Group 2"))) %>%
    ggplot(aes(x = fct_rev(scenF), y = avg, fill = scenF)) +
    geom_flat_violin(position = position_nudge(x = .2, y = 0), alpha = .8, scale = "width") +
    geom_point(aes(y = avg, color = scenF), 
               position = position_jitter(width = .15), 
               size = .5, alpha = 0.8) +
    geom_boxplot(width = .1, outlier.shape = NA, alpha = 0.85, fill = "white", notch = TRUE) +
    ylab("Mean Ability") +
    guides(fill = FALSE) +
    guides(color = FALSE) +
    coord_flip() +
    expand_limits(x = max(longstat$scenario) + 1) +
    xlab("") +
    facet_wrap(~Group, ncol = 3) +
    scale_colour_brewer(name = "Scenario", palette = "Set1") +
    scale_fill_brewer(name = "Scenario", palette = "Set1") +
    theme_few() + 
    theme(panel.background = element_rect(fill = "gray85"),
          axis.text.x = element_text(size = 6))
  
  ggsave(paste0("Plots/Figure_", s, "_3b.tiff"), width = 6.5, height = 4.5) # High quality for submission
  ggsave(paste0("Plots/Figure_", s, "_3b.png"), width = 6.5, height = 4.5) # For Word doc
}

## Bland-Altman Plots

nPoints <- 500

set.seed(48293)
dsAll %>%
  mutate(Scenario = factor(scenario, levels = 1:9, labels = paste0("Scenario ", 1:9))) %>%
  group_by(scenario) %>%
  sample_n(nPoints) %>%
  ggplot(aes(x = ability, y = resid, colour = factor(group), shape = factor(group))) +
  geom_point() +
  xlab(expression(theta)) +
  ylab("Unstandardized Residual") +
  facet_wrap(~Scenario, ncol = 3) +
  geom_hline(yintercept = c(-.3, .3), lty = 2) +
  scale_colour_manual(values = c("#f1a340d9", "#998ec3d9"), name = "Group") +
  scale_shape_discrete(name = "Group") +
  ylim(-1.65, 1.65) + 
  theme_few() + 
  theme(panel.background = element_rect(fill = "#f7f7f7"))

ggsave("Plots/BlandAltman_us.png", width = 6.5, height = 4.5) # For Word doc
ggsave("Plots/BlandAltman_us.tiff", width = 6.5, height = 4.5) # High quality for submission

set.seed(48293)
dsAll %>%
  mutate(Scenario = factor(scenario, levels = 1:9, labels = paste0("Scenario ", 1:9))) %>%
  group_by(scenario) %>%
  sample_n(nPoints) %>%
  ggplot(aes(x = ability, y = resid_st, colour = factor(group), shape = factor(group))) +
  geom_point() +
  xlab(expression(theta)) +
  ylab("Standardized Residual") +
  facet_wrap(~Scenario, ncol = 3) +
  geom_hline(yintercept = c(-.3, .3), lty = 2) +
  scale_colour_manual(values = c("#f1a340d9", "#998ec3d9"), name = "Group") +
  scale_shape_discrete(name = "Group") +
  ylim(-1.65, 1.65) +
  theme_few() + 
  theme(panel.background = element_rect(fill = "#f7f7f7"))

ggsave("Plots/BlandAltman_st.png", width = 6.5, height = 4.5) # For Word doc
ggsave("Plots/BlandAltman_st.tiff", width = 6.5, height = 4.5) # High quality for submission

badat <- dsAll %>% 
  mutate(w30 = ifelse(resid < .30, 1, 0),
         w30_st = ifelse(resid_st < .30, 1, 0)) %>%
  group_by(group, scenario)

# Proportion of unstandardized residuals outside of .30 by group

table(badat$w30, badat$scenario, badat$group) %>%
  prop.table(margin = c(2, 3))

# Proportion of unstandardized residuals outside of .30 collapsing over group

table(badat$w30, badat$scenario) %>%
  prop.table(margin = 2)

# Proportion of standardized residuals outside of .30 by group

table(badat$w30_st, badat$scenario, badat$group) %>%
  prop.table(margin = c(2, 3))

# Proportion of standardized residuals outside of .30 collapsing over group

table(badat$w30_st, badat$scenario) %>%
  prop.table(margin = 2)



# sumstatt <- readRDS("Results/sumstat_2020-01-07-13-23.rds")
# ---------------------------------End Scenarios -------------------------------





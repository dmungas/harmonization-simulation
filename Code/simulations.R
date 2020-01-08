require(readr)
require(MplusAutomation)
require(mirt)
require(tidyr)
require(ggplot2)
require(lubridate)
source("~/Research/Code/runMplus.R")
source("~/Research/Code/extractVarNames.R")
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

hrmecal <- mirt(hrme[,union(vars,varsm)],m1hrme,pars=hrme_par)
# hrmecal <- mirt(hrme[,union(vars,varsm)],m1hrme)
# coef(hrmecal)
# parameters(hrmecal)

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

### scenario 1 - all items shared

items <- union(vars,varsm)

scen_1 <- equateSim(seed=21589,grp_mean=c(mean_hrs,mean_mhas),
    grp_sd=c(sd_hrs,sd_mhas),n_samp=500,n_rep_theta=1,n_itm=length(items),
    pars=hrme_par,n_rep=100,itms1=items,itms2=items,fsc_method="EAP",
    mod_res_obj=hrmecal)

sumstat1 <- summaryStats(scen_1[["summary"]])

### end scenario 1

### scenario 2 - UMON. UDAY, UYER as linking items, actual items in HRS and MHAS

scen_2 <- equateSim(seed=21589,grp_mean=c(mean_hrs,mean_mhas),
      grp_sd=c(sd_hrs,sd_mhas),n_samp=500,n_rep_theta=1,n_itm=length(items),
      pars=hrme_par,n_rep=500,itms1=vars,itms2=varsm,fsc_method="EAP",
      mod_res_obj=hrmecal)

sumstat2 <- summaryStats(scen_2[["summary"]])

### end scenario 2

### scenario 3 - UMON. UDAY, UYER UIWR as linking items

itms2 <- c(varsm,"UIWR")

scen_3 <- equateSim(seed=21589,grp_mean=c(mean_hrs,mean_mhas),
                    grp_sd=c(sd_hrs,sd_mhas),n_samp=500,n_rep_theta=1,n_itm=length(items),
                    pars=hrme_par,n_rep=500,itms1=vars,itms2=itms2,fsc_method="EAP",
                    mod_res_obj=hrmecal)

sumstat3 <- summaryStats(scen_3[["summary"]])

### end scenario 3

### scenario 4 - HRS Items in both studies


scen_4 <- equateSim(seed=21589,grp_mean=c(mean_hrs,mean_mhas),
                    grp_sd=c(sd_hrs,sd_mhas),n_samp=500,n_rep_theta=1,n_itm=length(vars),
                    pars=hrs_par,n_rep=500,itms1=vars,itms2=vars,fsc_method="EAP",
                    mod_res_obj=hrscal)

sumstat4 <- summaryStats(scen_4[["summary"]])

### end scenario 4

### scenario 5 - MHAS Items in both studies


scen_5 <- equateSim(seed=21589,grp_mean=c(mean_hrs,mean_mhas),
                    grp_sd=c(sd_hrs,sd_mhas),n_samp=500,n_rep_theta=1,n_itm=length(varsm),
                    pars=mex_par,n_rep=500,itms1=varsm,itms2=varsm,fsc_method="EAP",
                    mod_res_obj=mexcal,save_sims=TRUE)

sumstat5 <- summaryStats(scen_5[["summary"]])

scen_5_data <- scen_5[["datasets"]]

### end scenario 5

# dat1 <- scen_5_data[scen_5_data$sample == 1,]
# 
# plot(dat1$resid ~ dat1$ability)
# abline(lm(dat1$resid ~ dat1$ability))
# 
# ggplot(data=scen_5_data[scen_5_data$sample == 8,], aes(ability,resid)) + 
#   geom_point() + geom_smooth(method = "loess")


### scenario 6 - UMON. UDAY, UYER as linking items, actual items in HRS and MHAS


vars6 <- vars[!vars %in% c("UMON",'UYER')]
varsm6 <- varsm[!varsm %in% c("UMON",'UYER')]

scen_6 <- equateSim(seed=21589,grp_mean=c(mean_hrs,mean_mhas),
                    grp_sd=c(sd_hrs,sd_mhas),n_samp=500,n_rep_theta=1,n_itm=length(items),
                    pars=hrme_par,n_rep=500,itms1=vars,itms2=varsm6,fsc_method="EAP",
                    mod_res_obj=hrmecal)

sumstat6 <- summaryStats(scen_6[["summary"]])

### end scenario 6


### Merge summary statistics from scenarios

sumstat1$scenario <- 1
sumstat2$scenario <- 2
sumstat3$scenario <- 3
sumstat4$scenario <- 4
sumstat5$scenario <- 5
sumstat6$scenario <- 6

sumstat <- rbind(sumstat1,sumstat2,sumstat3,sumstat4,sumstat5,sumstat6)
sumstat$scenario_label <- ifelse(sumstat$scenario == 1,"HRS+MHAS_all_shared",
      ifelse(sumstat$scenario == 2,"HRS+MHAS_UMON_UDAY_UYER_shared",
      ifelse(sumstat$scenario == 3,"HRS+MHAS_UMON_UDAY_UYER_UIWR_shared",
      ifelse(sumstat$scenario == 4,"HRS_all_shared",
      ifelse(sumstat$scenario == 5,"MHAS_all_shared",
      ifelse(sumstat$scenario == 6,"HRS+HMAS_UDAY_shared",
      NA))))))

save(sumstat,mean_sd,file=paste0("Results/simulation_results_",format(Sys.time(), '%Y-%m-%d-%H-%M'),".RData"))
# saveRDS(sumstat,file=paste0("Results/sumstat_",format(Sys.time(), '%Y-%m-%d-%H-%M'),".rds"))

# sumstatt <- readRDS("Results/sumstat_2020-01-07-13-23.rds")
# ---------------------------------End Scenarios -------------------------------


require(readr)
require(MplusAutomation)
require(mirt)
require(tidyr)
require(ggplot2)
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
# coef(hrmecal)
parameters(hrmecal)

fsc <- fscores(hrmecal,method="EAP")
hrme$ability_est <- fsc

# table(hrme$study_name_short)
hrme$study <- ifelse(grepl("HRS",hrme$study_name_short),"HRS","MHAS")

mean_hrs <- mean(hrme[hrme$study %in% "HRS","ability_est"])
sd_hrs <- sd(hrme[hrme$study %in% "HRS","ability_est"])
mean_mhas <- mean(hrme[hrme$study %in% "MHAS","ability_est"])
sd_mhas <- sd(hrme[hrme$study %in% "MHAS","ability_est"])

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
      pars=hrme_par,n_rep=100,itms1=vars,itms2=varsm,fsc_method="EAP",
      mod_res_obj=hrmecal)

sumstat2 <- summaryStats(scen_2[["summary"]])

### end scenario 2

### scenario 3 - UMON. UDAY, UYER UIWR as linking items

itms2 <- c(varsm,"UIWR")

scen_3 <- equateSim(seed=21589,grp_mean=c(mean_hrs,mean_mhas),
                    grp_sd=c(sd_hrs,sd_mhas),n_samp=500,n_rep_theta=1,n_itm=length(items),
                    pars=hrme_par,n_rep=100,itms1=vars,itms2=itms2,fsc_method="EAP",
                    mod_res_obj=hrmecal)

sumstat3 <- summaryStats(scen_3[["summary"]])

### end scenario 3

### scenario 4 - HRS Items in both studies


scen_4 <- equateSim(seed=21589,grp_mean=c(mean_hrs,mean_mhas),
                    grp_sd=c(sd_hrs,sd_mhas),n_samp=500,n_rep_theta=1,n_itm=length(vars),
                    pars=hrs_par,n_rep=100,itms1=vars,itms2=vars,fsc_method="EAP",
                    mod_res_obj=hrscal)

sumstat4 <- summaryStats(scen_4[["summary"]])

### end scenario 4

### scenario 5 - MHAS Items in both studies


scen_5 <- equateSim(seed=21589,grp_mean=c(mean_hrs,mean_mhas),
                    grp_sd=c(sd_hrs,sd_mhas),n_samp=500,n_rep_theta=1,n_itm=length(varsm),
                    pars=mex_par,n_rep=100,itms1=varsm,itms2=varsm,fsc_method="EAP",
                    mod_res_obj=mexcal,save_sims=TRUE)

sumstat5 <- summaryStats(scen_5[["summary"]])

scen_5_data <- scen_5[["datasets"]]

### end scenario 5

# ---------------------------------End Scenarios -------------------------------

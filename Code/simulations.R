library(pacman)
p_load(tidyverse, pander, readr, mirt, tidyr, ggplot2, lubridate, 
       DT, RColorBrewer, ggthemes, officer, knitr, progress,
       stringr)
#source("~/Research/Code/runMplus.R")
#source("~/Research/Code/extractVarNames.R")
#source("Code/simulation_functions.R")

if(!dir.exists("Output")) {
  dir.create("Output")
}
outFolderName <- format(Sys.time(), "%Y-%m-%d_%H-%M")
outPath <- file.path("Output", outFolderName)
dir.create(outPath)
logPath <- file.path(outPath, "Logs")
dir.create(logPath)
plotPath <- file.path(outPath, "Plots")
dir.create(plotPath)
resultsPath <- file.path(outPath, "Results")
dir.create(resultsPath)

runOrRead <- "run" # If "read", reads in old files. If "run," generates new simulations.
#readFolder <- "Output/2020-10-05_11-00-AM/Results"

Seed <- 21589
mus <- c(0, -0.24)
#sigmas <- c(sd_hrs, sd_mhas) # set below (after calibration)
refGrp <- 1
repTheta <- 1
reps <- 500
fsMeth <- "EAP"
n1 <- 500
n2 <- 500

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
hrscal_inf <- infoCalc(hrscal)
plot(hrscal_inf$ability,hrscal_inf$information)


m2hrs <- mirt.model('cog = 1-11')

# HRS calibration excluding delayed word recall
hrs_par_0dr <- mirt(hrs[,vars[!vars %in% "UDWR"]],m2hrs,pars='values')
hrscal_0dr <- mirt(hrs[,vars[!vars %in% "UDWR"]],m2hrs,pars=hrs_par_0dr)
# coef(hrscal_0dr)
hrs_0dr_inf <- infoCalc(hrscal_0dr)
plot(hrs_0dr_inf$ability,hrs_0dr_inf$information)

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

sigmas <- c(sd_hrs, sd_mhas)
                                        
# *************************** Simulation Scenarios *****************************

### scenario 1 - all items shared

items <- union(vars,varsm)

if(runOrRead == "run") {
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
                    log_file = file.path(logPath, "scenario1.txt"),
                    verbose = FALSE, 
                    prog_bar = TRUE,
                    scale_lv = "fix_sample_est")
} 

if(runOrRead == "read") {
  scen_1 <- readRDS(file.path(readFolder, "scen_1.Rds"))
}

sumstat1 <- summaryStats(scen_1[["summary"]])
longstat1 <- scen_1[["summary"]]
ds1 <- scen_1$datasets

if(runOrRead == "run") {
saveRDS(scen_1, file=file.path(resultsPath, "scen_1.Rds"))
}

### end scenario 1


### Scenario 2


### scenario 2 - UMON. UDAY, UYER as linking items, actual items in HRS and MHAS

if(runOrRead == "run") {
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
                    log_file = file.path(logPath, "scenario2.txt"),
                    verbose = FALSE, 
                    prog_bar = TRUE,
                    scale_lv = "fix_sample_est")
}

if(runOrRead == "read") {
  scen_2 <- readRDS(file.path(readFolder, "scen_2.Rds"))
}

sumstat2 <- summaryStats(scen_2[["summary"]])
longstat2 <- scen_2[["summary"]]
ds2 <- scen_2$datasets

if(runOrRead == "run") {
saveRDS(scen_2, file=file.path(resultsPath, "scen_2.Rds"))
}
  
### end scenario 2


### Scenario 3


### scenario 3 - UMON. UDAY, UYER UIWR as linking items

itms2 <- c(varsm,"UIWR")

if(runOrRead == "run") {
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
                    log_file = file.path(logPath, "scenario3.txt"),
                    verbose = FALSE, 
                    prog_bar = TRUE,
                    scale_lv = "fix_sample_est")
}

if(runOrRead == "read") {
  scen_3 <- readRDS(file.path(readFolder, "scen_3.Rds"))
}

sumstat3 <- summaryStats(scen_3[["summary"]])
longstat3 <- scen_3[["summary"]]
ds3 <- scen_3$datasets

if(runOrRead == "run") {
saveRDS(scen_3, file=file.path(resultsPath, "scen_3.Rds"))
}
### end scenario 3


### Scenario 4


### scenario 4 - HRS Items in both studies

if(runOrRead == "run") {
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
                    log_file = file.path(logPath, "scenario4.txt"),
                    verbose = FALSE, 
                    prog_bar = TRUE,
                    scale_lv = "fix_sample_est")
}

if(runOrRead == "read") {
  scen_4 <- readRDS(file.path(readFolder, "scen_4.Rds"))
}

sumstat4 <- summaryStats(scen_4[["summary"]])
longstat4 <- scen_4[["summary"]]
ds4 <- scen_4$datasets

if(runOrRead == "run") {
saveRDS(scen_4, file=file.path(resultsPath, "scen_4.Rds"))
}
  
### end scenario 4


### Scenario 5


### scenario 5 - MHAS Items in both studies

if(runOrRead == "run") {
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
                    log_file = file.path(logPath, "scenario5.txt"),
                    verbose = FALSE, 
                    prog_bar = TRUE,
                    scale_lv = "fix_sample_est")
}

if(runOrRead == "read") {
  scen_5 <- readRDS(file.path(readFolder, "scen_5.Rds"))
}

sumstat5 <- summaryStats(scen_5[["summary"]])
longstat5 <- scen_5[["summary"]]
ds5 <- scen_5$datasets

if(runOrRead == "run") {
saveRDS(scen_5, file=file.path(resultsPath, "scen_5.Rds"))
}

### end scenario 5


### Scenario 6


### scenario 6 -  UDAY as linking item

vars6 <- vars[!vars %in% c("UMON",'UYER')]
varsm6 <- varsm[!varsm %in% c("UMON",'UYER')]

if(runOrRead == "run") {
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
                    log_file = file.path(logPath, "scenario6.txt"),
                    verbose = FALSE, 
                    prog_bar = TRUE,
                    scale_lv = "fix_sample_est")
}

if(runOrRead == "read") {
  scen_6 <- readRDS(file.path(readFolder, "scen_6.Rds"))
}

sumstat6 <- summaryStats(scen_6[["summary"]])
longstat6 <- scen_6[["summary"]]
ds6 <- scen_6$datasets

if(runOrRead == "run") {
saveRDS(scen_6, file=file.path(resultsPath, "scen_6.Rds"))
}

### end scenario 6


### Scenario 7


### scenario 7 - UIWR as linking item, MHAS without UMON UDAY UYER

varsm7 <- c(varsm[!varsm %in% c("UMON",'UYER','UDAY')],'UIWR')

if(runOrRead == "run") {
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
                    log_file = file.path(logPath, "scenario7.txt"),
                    verbose = FALSE, 
                    prog_bar = TRUE,
                    scale_lv = "fix_sample_est")
}

if(runOrRead == "read") {
  scen_7 <- readRDS(file.path(readFolder, "scen_7.Rds"))
}

sumstat7 <- summaryStats(scen_7[["summary"]])
longstat7 <- scen_7[["summary"]]
ds7 <- scen_7$datasets

if(runOrRead == "run") {
saveRDS(scen_7, file=file.path(resultsPath, "scen_7.Rds"))
}
### end scenario 7


### Scenario 8


### scenario 8 - UIWR as linking item, HRS without UMON UDAY UYER

vars8 <- vars[!vars %in% c("UMON",'UYER','UDAY')]
varsm8 <- c(varsm,'UIWR')

if(runOrRead == "run") {
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
                    log_file = file.path(logPath, "scenario8.txt"),
                    verbose = FALSE, 
                    prog_bar = TRUE,
                    scale_lv = "fix_sample_est")
}

if(runOrRead == "read") {
  scen_8 <- readRDS(file.path(readFolder, "scen_8.Rds"))
}

sumstat8 <- summaryStats(scen_8[["summary"]])
longstat8 <- scen_8[["summary"]]
ds8 <- scen_8$datasets

if(runOrRead == "run") {
saveRDS(scen_8, file=file.path(resultsPath, "scen_8.Rds"))
}
### end scenario 8


### Scenario 9



if(runOrRead == "run") {
scen_9 <- equateSim(seed = Seed, 
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
                    log_file = file.path(logPath, "scenario9.txt"),
                    verbose = FALSE, 
                    prog_bar = TRUE,
                    scale_lv = "fix_sample_est",
                    man_override = list(iname = c("UDAY", "UMON"),
                                        parameter = c("b", "b"), 
                                        val = c(0, 1.928)))
}

if(runOrRead == "read") {
  scen_9 <- readRDS(file.path(readFolder, "scen_9.Rds"))
}

sumstat9 <- summaryStats(scen_9[["summary"]])
longstat9 <- scen_9[["summary"]]
ds9 <- scen_9$datasets

if(runOrRead == "run") {
saveRDS(scen_9, file=file.path(resultsPath, "scen_9.Rds"))
}
### end scenario 9


### Scenario 10



if(runOrRead == "run") {
scen_10 <- equateSim(seed = Seed, 
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
                    log_file = file.path(logPath, "scenario10.txt"),
                    verbose = FALSE, 
                    prog_bar = TRUE,
                    scale_lv = "fix_sample_est",
                    man_override = list(iname = c("UDAY", "UMON", "UDAY", "UMON", "UYER"),
                                        parameter = c("b", "b", "a1", "a1", "a1"), 
                                        val = c(0, 1.928, 4, 4, 4)))
}

if(runOrRead == "read") {
  scen_10 <- readRDS(file.path(readFolder, "scen_10.Rds"))
}

sumstat10 <- summaryStats(scen_10[["summary"]])
longstat10 <- scen_10[["summary"]]
ds10 <- scen_10$datasets

if(runOrRead == "run") {
saveRDS(scen_10, file=file.path(resultsPath, "scen_10.Rds"))
}
### end scenario 10


nscen <- sum(str_detect(ls(), "^scen_[0-9]"))

### Merge summary statistics from scenarios

if(nscen > 0){
for(ns in 1:nscen){
  #temp <- dynGet(paste0("sumstat", ns)) # dynGet works when knitting
  temp <- get(paste0("sumstat", ns))
  temp$scenario <- ns
  assign(paste0("sumstat", ns), temp)
  
  #temp <- dynGet(paste0("longstat", ns)) # dynGet works when knitting
  temp <- get(paste0("longstat", ns))
  temp$scenario <- ns
  assign(paste0("longstat", ns), temp)
  
  #temp <- dynGet(paste0("ds", ns)) # dynGet works when knitting
  temp <- get(paste0("ds", ns))
  temp$scenario <- ns
  assign(paste0("ds", ns), temp)
}
}

# longstat <- lapply(paste0("longstat", 1:nscen), dynGet) %>%  # dynGet works when knitting
#  do.call(bind_rows, .)
longstat <- lapply(paste0("longstat", 1:nscen), get) %>%
  do.call(bind_rows, .)                    

# dsAll <- lapply(paste0("ds", 1:nscen), dynGet) %>%  # dynGet works when knitting
#  do.call(bind_rows, .)
dsAll <- lapply(paste0("ds", 1:nscen), get) %>%
  do.call(bind_rows, .)
                    
# sumstat <- lapply(paste0("sumstat", 1:nscen), dynGet) %>% # dynGet works when knitting
#  do.call(bind_rows, .)
sumstat <- lapply(paste0("sumstat", 1:nscen), get) %>%
  do.call(bind_rows, .)
                    
scenLabs <- c(
  "HRS+MHAS_all_shared",
  "HRS+MHAS_UMON_UDAY_UYER_shared",
  "HRS+MHAS_UMON_UDAY_UYER_UIWR_shared",
  "HRS_all_shared",
  "MHAS_all_shared",
  "HRS+HMAS_UDAY_shared",
  "HRS+MHAS_UIWR_shared_no_MHAS_dates",
  "HRS+MHAS_UIWR_shared_no_HRS_dates",
  "HRS+MHAS_UMON(b=1.928)_UDAY(b=0)_UYER(b=-1.928)_shared",
  "HRS+MHAS_UMON(b=1.928,a1=4)_UDAY(b=0,a1=4)_UYER(b=-1.928,a1=4)_shared"
)

sumstat$scenario_label <- factor(sumstat$scenario, 
                                 levels = 1:nscen, 
                                 labels = scenLabs[1:nscen])

if(runOrRead == "run") {
saveRDS(sumstat,file=file.path(resultsPath, paste0("sumstat_", outFolderName, ".rds")))
#sumstat <- readRDS("Results/sumstat_2020-07-08-18-13.rds")

save(list = ls(all.names = TRUE),
     file=file.path(resultsPath, paste0("simulation_results_", outFolderName, ".RData")),
     envir = environment())
#load("Output/2020-08-06_09-04-AM/Results/simulation_results_2020-08-06_09-04-AM.RData")

write.csv(sumstat, file=file.path(resultsPath, paste0("sumstat", outFolderName, ".csv")))
}

if(runOrRead == "read") {
  #sumstat <- readRDS(file.path(readFolder, paste0("sumstat_", str_sub(readFolder, 8, -9), ".Rds")))
  write.csv(sumstat, file=file.path(resultsPath, paste0("sumstat", outFolderName, ".csv")))
}
                    
### Compile results                    
                    
DT::datatable(sumstat, options = list(pageLength = 3*nscen), rownames = FALSE) %>%
  formatRound(columns = names(sumstat)[-c(1, ncol(sumstat))], digits = 3) %>%
  formatRound(columns = c("n_rep", "scenario"), digits = 0)


read_docx() %>%  # a new, empty document
  body_add_table(sumstat %>% 
                   mutate_if(is.numeric, round, 3) %>%
                   dplyr::select(scenario, mean, sd, bias, bias_pct, ese, rmse, r_theta_est) %>%
                   filter(group == "Combined"), style = "table_template") %>% 
  print(target=file.path(resultsPath, "Table3.docx"))

read_docx() %>%  # a new, empty document
  body_add_table(sumstat %>% 
                   mutate_if(is.numeric, round, 3) %>%
                   dplyr::select(scenario, mean, sd, bias, bias_pct, ese, rmse, r_theta_est) %>%
                   filter(group == "Group 1"), style = "table_template") %>% 
  print(target=file.path(resultsPath, "Table4.docx"))

read_docx() %>%  # a new, empty document
  body_add_table(sumstat %>% 
                   mutate_if(is.numeric, round, 3) %>%
                   dplyr::select(scenario, mean, sd, bias, bias_pct, ese, rmse, r_theta_est) %>%
                   filter(group == "Group 2"), style = "table_template") %>% 
  print(target=file.path(resultsPath, "Table5.docx"))

# Plot scenarios

sumlong <- sumstat %>%
  dplyr::select(group, mean, mean_st, sd, sd_st, n_rep, scenario) %>%
  pivot_longer(cols = c(mean, mean_st, sd, sd_st),
               names_to = c(".value", "stat"),
               names_sep = c("_")) %>%
  mutate_at(vars(stat), replace_na, "us") %>%
  mutate(low95 = mean - qnorm(.975)*sd/sqrt(n_rep),
         up95 = mean + qnorm(.975)*sd/sqrt(n_rep)) %>%
  mutate(scenF = factor(scenario, levels = 1:nscen, labels = paste0("Scenario ", 1:nscen))) %>%
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

s1cols <- brewer.pal(9, "Set1")
s1cols <- c(s1cols, "#232654")

for(s in unique(longstat$type)){
  rcp <- longstat %>%
    filter(type == s) %>%
    filter(statistic == "est_mean") %>%
    mutate(scenF = factor(scenario, levels = 1:nscen, labels = paste0("Scenario ", 1:nscen))) %>%
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
    scale_colour_manual(name = "Scenario", values = s1cols) + 
    scale_fill_manual(name = "Scenario", values = s1cols) + 
    #scale_colour_brewer(name = "Scenario", palette = "Set3") +
    #scale_fill_brewer(name = "Scenario", palette = "Set3") +
    geom_hline(data = thetalines3, aes(yintercept = value), 
               lty = 2) +
    theme_few() + 
    theme(panel.background = element_rect(fill = "gray85"),
          axis.text.x = element_text(size = 6))
  
  ggsave(file.path(plotPath, paste0("Figure_", s, ".tiff")), width = 6.5, height = 4.5) # High quality for submission
  ggsave(file.path(plotPath, paste0("Figure_", s, ".png")), width = 6.5, height = 4.5) # For Word doc
  

}

## Bland-Altman Plots

nPoints <- 500

set.seed(48293)
ba1 <- dsAll %>%
  mutate(Scenario = factor(scenario, levels = 1:nscen, labels = paste0("Scenario ", 1:nscen))) %>%
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

ggsave(file.path(plotPath, "BlandAltman_us.png"), width = 6.5, height = 4.5) # For Word doc
ggsave(file.path(plotPath, "BlandAltman_us.tiff"), width = 6.5, height = 4.5) # High quality for submission

set.seed(48293)
ba2 <- dsAll %>%
  mutate(Scenario = factor(scenario, levels = 1:nscen, labels = paste0("Scenario ", 1:nscen))) %>%
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

ggsave(file.path(plotPath, "BlandAltman_st.png"), width = 6.5, height = 4.5) # For Word doc
ggsave(file.path(plotPath, "BlandAltman_st.tiff"), width = 6.5, height = 4.5) # High quality for submission

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





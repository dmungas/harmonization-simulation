require(MplusAutomation)
require(mirt)
require(tidyr)
require(ggplot2)
source("~/Research/Code/runMplus.R")
source("~/Research/Code/extractVarNames.R")
source("Code/simulation_functions.R")


tics <- read.csv("Data/PITCH-response-data-fhl.csv")

# extractVarNames(tics)

nms <- c('newid','study_wave_number','study_name_short','UADR','UANI',
    'UBAK','UDAT','UDAY','UDIGB','UDIGF','UDIR','UDRG','UDSB','UDTS','UDWR',
    'UEYE','UFC2','UFCO','UFCO2','UFCO3','UFRE1','UFRE2','UHL1','UHL2','UHL3',
    'UHL4','UIWR','UIWR2','UIWR3','ULET','UMON','UNM1','UNM2','UNM3','UNM4',
    'UNM5','UNM6','UNM7','UNM8','UNMS','UNUM','UNUM2','UPCC','UPCR','UPLC',
    'UPM1','UPM3','UPM4','URE2','URE3','URE4','URE5','URE6','UREG','UREP',
    'UREP2','USEA','USRT','USUB','UTR1','UTR2','UVA1','UVRS','UVSC','UWD',
    'UWLD','UWR1','UWR2','UWR3','UWRT','UYER')

# summary(tics$UADR)




title <- "PITCH IRT Analysis"

analysis <- "estimator = mlr;
"

output <- "sampstat;stdyx;svalues;"

vars <- c('UBAK','UDAT','UDAY','UDWR','UIWR','UMON','UNM1','UNM2','UNM5','UNM6',
     'USUB','UYER')
vars2 <- c('UBAK','UDAT','UDAY','UIWR','UMON','UNM1','UNM2','UNM5','UNM6',
          'USUB','UYER')

varcat <- "
categorical = ubak udat uday udwr uiwr umon unm1 unm2 unm5 unm6 usub,uyer;
"

model_hrs <- "
cog by ubak* udat-uyer;
cog@1;
"

model_hrs_bif <- "
cog by ubak* udat-uyer;
mem by uiwr* udwr (1);
cog@1;
mem@1;
cog with mem@0;
"


rdata <- tics[tics$study_name_short %in% c("HRS_CODA_W6","HRS_W6",]
  # rdata <- tics
  
dirnm <- "/Analysis"
  
  #   ECog Memory AA White - all free (except mrem)
  
  modelout <- "pitch_hrs.inp"
  
  varcat1 <- sub(" udwr","",varcat)
    
  variable <- varcat1
    
res_pitch_hrs <- runMplus(TITLE = title, MODEL = model_cog, usevariables = vars2,
            rdata=rdata, OUTPUT = output, dirnm = dirnm, modelout = modelout,
            VARIABLE = variable, ANALYSIS=analysis)
extract(res_pitch_cog$results, type = "stdyx", summaries=summaries)
rm(res_pitch_hrs)



modelout <- "pitch_hrs_wmslv.inp"

analysis2 <- sub("mlr","wlsmv",analysis)

variable <- varcat1

res_pitch_hrs_wlsmv <- runMplus(TITLE = title, MODEL = model_cog, usevariables = vars2,
        rdata=rdata, OUTPUT = output, dirnm = dirnm, modelout = modelout,
          VARIABLE = variable, ANALYSIS=analysis2)
extract(res_pitch_hrs_wlsmv$results, type = "stdyx", summaries=summaries)
rm(res_pitch_hrs_wlsmv)



modelout <- "pitch_hrs_bif.inp"

variable <- varcat

res_pitch_hrs_bif <- runMplus(TITLE = title, MODEL = model_hrs_bif, usevariables = vars,
    rdata=rdata, OUTPUT = output, dirnm = dirnm, modelout = modelout,
    VARIABLE = variable, ANALYSIS=analysis)
extract(res_pitch_hrs_bif$results, type = "stdyx", summaries=summaries)
rm(res_pitch_hrs_bif)



modelout <- "pitch_hrs_wlsmv_bif.inp"

variable <- varcat

res_pitch_hrs_wlsmv_bif <- runMplus(TITLE = title, MODEL = model_hrs_bif, usevariables = vars,
        rdata=rdata, OUTPUT = output, dirnm = dirnm, modelout = modelout,
        VARIABLE = variable, ANALYSIS=analysis2)
extract(res_pitch_hrs_wlsmv_bif$results, type = "stdyx", summaries=summaries)
rm(res_pitch_hrs_wlsmv_bif)







rdata <- tics[tics$study_name_short %in% c("MHAS_W1","MHAS_W2"),]


varsm <- c("UDAY","UFCO2","UFRE1","UMON","UVSC","UWD","UWR1","UWR2","UWR3","UYER")

varcatm <- "categorical = uday ufco2 ufre1 umon uvsc uwd uwr1 uwr2 uwr3 uyer;"

model_mhas <- "
cog by uday* ufco2-uyer;
cog@1;
"

modelout <- "pitch_mhas_wlsmv.inp"

variable <- varcatm

res_pitch_mhas_wlsmv <- runMplus(TITLE = title, MODEL = model_mhas, usevariables = varsm,
      rdata=rdata, OUTPUT = output, dirnm = dirnm, modelout = modelout,
      VARIABLE = variable, ANALYSIS=analysis2)
extract(res_pitch_mhas_wlsmv$results, type = "stdyx", summaries=summaries)
rm(res_pitch_mhas_wlsmv)

### mirt calibration of tics

hrs <- tics[tics$study_name_short %in% c("HRS_CODA_W6","HRS_W6"),vars]

hrs$n_itm <- apply(hrs,1,function(x) sum(!is.na(x)))

hrs <- hrs[!hrs$n_itm ==0,]

m1hrs <- mirt.model('cog = 1-12')
hrs_par <- mirt(hrs[,1:12],m1hrs,pars='values')

hrscal <- mirt(hrs[,1:12],m1hrs,pars=hrs_par)
coef(hrscal)




require(psych)

discr <- c(2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2)

ds1 <- sim.poly.npl(nvar = 30 ,n = 500,low=-3,high=3,a=discr,
        c=0,z=1,d=NULL, mu=0,sd=1,cat=2)

diff <-  data.frame(ds1[["difficulty"]])
names(diff) <- "difficulty"
disc <-  data.frame(ds1[["discrimination"]])
names(disc) <- "discrimination"


t1 <- data.frame(ds1[["items"]])
t2 <- data.frame(ds1[["theta"]])

names(t2) <- "ability"

t2 <- cbind(t2,t1)




title <- "Simulated IRT Analysis"

analysis <- "estimator = mlr;
"

output <- "sampstat;stdyx;svalues;"

varsim <- c('V1','V2','V3','V4','V5','V6','V7','V8','V9','V10',
          'V11','V12','V13','V14','V15')
varsim2 <- c('V1','V2','V3','V4','V5','V6','V7','V8','V9','V10',
            'V11','V12','V13','V14','V15','V16','V17','V18','V19','V20',
            'V21','V22','V23','V24','V25','V26','V27','V28','V29','V30')

varcatsim <- "
categorical = v1-v15;
"

varcatsim2 <- "
categorical = v1-v30;
"

variable <- varcatsim

model_sim <- "
cog by v1* v2-v15;
cog@1;
"

model_sim_2 <- "
cog by v1*-v15;
cog@1;

v7@1;
v8@1;
v9@1;

[v7@-0.3];
[v8@0.0];
[v9@0.3];
"

model_sim_2 <- "
cog by v1*-v15;
cog@1;

v7@1;
v8@1;
v9@1;

[v7@-0.3];
[v8@0.0];
[v9@0.3];
"

model_sim_3 <- "
cog by v1@1 v2@1 v3@1 v4@1 v5@1 v6@1 v7@1 v8@1 v9@1 v10@1 v11@1 v12@1 v13@1 v14@1 v15@1;
!cog@1;

[v1$1@-2.1];
[v2$1@-1.8];
[v3$1@-1.5];
[v4$1@-1.2];
[v5$1@-0.9];
[v6$1@-0.6];
[v7$1@-0.3];
[v8$1@0.0];
[v9$1@0.3];
[v10$1@0.6];
[v11$1@0.9];
[v12$1@1.2];
[v13$1@1.5];
[v14$1@1.8];
[v15$1@2.1];
"


savedata <- "
FILE = 'irt_sim_fscores_3.dat';
SAVE = FSCORES;
"



rdata <- t2

dirnm <- "/Analysis"

#   ECog Memory AA White - all free (except mrem)

modelout <- "irt_sim.inp"

res_irt_sim <- runMplus(TITLE = title, MODEL = model_sim_2, usevariables = varsim,
      rdata=rdata, OUTPUT = output, dirnm = dirnm, modelout = modelout,
      VARIABLE = varcatsim, ANALYSIS=analysis, delete_data=FALSE,
      SAVEDATA=savedata)
extract(res_irt_sim$results, type = "stdyx", summaries=summaries)
rm(res_irt_sim)


fsc <- readModels("Analysis/irt_sim.out",what="savedata")$savedata

t3 <- cbind(t2,fsc[,c("COG","COG_SE")])

t3$resid <- t3$COG - t3$ability

plot(t3$ability,t3$COG)
plot(t3$ability,t3$resid)


modelout <- "irt_sim_3.inp"

res_irt_sim_3 <- runMplus(TITLE = title, MODEL = model_sim_3, usevariables = varsim,
                        rdata=rdata, OUTPUT = output, dirnm = dirnm, modelout = modelout,
                        VARIABLE = varcatsim, ANALYSIS=analysis, delete_data=FALSE,
                        SAVEDATA=savedata)
extract(res_irt_sim$results, type = "stdyx", summaries=summaries)
rm(res_irt_sim)


fsc <- readModels("Analysis/irt_sim_3.out",what="savedata")$savedata

t4 <- cbind(t2,fsc[,c("COG","COG_SE")])

t4$resid <- t4$COG - t4$ability

plot(t4$ability,t4$COG)
plot(t4$ability,t4$resid)








model_sim_4 <- "
cog by v1-v30@2;
"

for(i in 1:30){
  if (i==1){
    thresh <- paste("[",varsim2[1],"$1","@",diff[1,1]*disc[1,1],"];\n",sep="")
  } else {
    thresh <- paste(thresh,"[",varsim2[i],"$1","@",diff[i,1]*disc[i,1],"];\n",sep="")
  }
}

model_sim_4 <- paste(model_sim_4,thresh,sep="\n")

# Mplus EAP factor scores

modelout <- "irt_sim_mplus_eap.inp"

res_irt_mplus_eap <- runMplus(TITLE = title, MODEL = model_sim_4, usevariables = varsim2,
                          rdata=rdata, OUTPUT = output, dirnm = dirnm, modelout = modelout,
                          VARIABLE = varcatsim2, ANALYSIS=analysis, delete_data=FALSE,
                          SAVEDATA=savedata)
rm(res_irt_mplus_eap)

rm(res_irt_sim,res_irt_sim_3,res_irt_sim_4)

fsc <- readModels("Analysis/irt_sim_mplus_eap.out",what="savedata")$savedata

t4 <- cbind(t2,fsc[,c("COG","COG_SE")])

t4$resid <- t4$COG - t4$ability

plot(t4$ability,t4$COG)
plot(t4$ability,t4$resid)
abline(lm(t4$resid ~ t4$ability))

ggplot(data=t4,aes(x=ability,y=resid)) + geom_point() + geom_smooth(method="lm")





                    
rdata <- t2

savedata2 <- "
FILE = 'irt_sim_fscores_3.dat';
SAVE = FSCORES(50,10);
"


modelout <- "irt_sim_4.inp"

analysis2 <- sub("mlr","bayes",analysis)

res_irt_sim_4 <- runMplus(TITLE = title, MODEL = model_sim_4, usevariables = varsim2,
                          rdata=rdata, OUTPUT = output, dirnm = dirnm, modelout = modelout,
                          VARIABLE = varcatsim2, ANALYSIS=analysis2, delete_data=FALSE,
                          SAVEDATA=savedata2)
extract(res_irt_sim_4$results, type = "stdyx", summaries=summaries)
rm(res_irt_sim_4)


fsc <- readModels("Analysis/irt_sim_4.out",what="savedata")$savedata

t5 <- cbind(t2,fsc[,c("COG.Mean","COG.Median","COG.Standard.Deviation")])

t5$resid <- t5$COG.Median - t5$ability

plot(t5$ability,t5$COG.Median)
plot(t5$ability,t5$resid)
abline(lm(t5$resid ~ t5$ability))

mean(t5$ability)
mean(t5$COG.Median)
sd(t5$ability)
sd(t5$COG.Median)
mean(t5$abil_est_st)
sd(t5$abil_est_st)


t5$abil_est_st <- (t5$COG.Median - mean(t5$COG.Median))*(sd(t5$ability)/sd(t5$COG.Median)) +
    mean(t5$ability) # linear equating to original ability estimate

t5$resid_st <- t5$abil_est_st - t5$ability

plot(t5$ability,t5$abil_est_st)
plot(t5$ability,t5$resid_st)
abline(lm(t5$resid_st ~ t5$ability))



### mirt factor scores
itemnames <- mcal_par[mcal_par$name=="a1","item"]

mdl <- mirt.model('cog = 1-30')


mcal_par <- mirt(t2[,2:31],mdl,pars='values')

mcal_par[mcal_par$name == "a1","value"] <- 2.0
mcal_par[mcal_par$name == "a1","est"] <- FALSE
mcal_par[mcal_par$name == "d","value"] <- diff
mcal_par[mcal_par$name == "d","est"] <- FALSE
mcal_par[mcal_par$item == "GROUP","est"] <- TRUE


mcal <- mirt(t2[,2:31],mdl,itemtype='2PL',pars=mcal_par)
# coef(mcal)

fsc <- fscores(mcal,full.scores=TRUE, plausible.type="MH",plausible.draws=50)
fsc_med <- data.frame(list(fsc[[5]],fsc[[10]],fsc[[15]],fsc[[20]],fsc[[25]],fsc[[30]],fsc[[35]],
    fsc[[40]],fsc[[45]],fsc[[50]]))
names(fsc_med) <- c("D1","D2","D3","D4","D5","D6","D7","D8","D9","D10")
fsc_med$ability_est <- apply(fsc_med,1,median)

t6 <- cbind(t2,fsc_med)

t6$resid <- t6$ability_est - t6$ability

plot(t6$ability,t6$ability_est)
plot(t6$ability,t6$resid)
abline(lm(t6$resid ~ t6$ability))

sd(t6$ability_est)
sd(t6$ability)

# t6[,c("ability","ability_est","resid")]

fsc <- data.frame(fscores(mcal,full.scores=TRUE, method="EAP"))
names(fsc) <- "ability_est"

t7 <- cbind(t2,fsc)
t7$resid <- t7$ability_est - t7$ability


plot(t7$ability,t7$ability_est)
plot(t7$ability,t7$resid)
abline(lm(t7$resid ~ t7$ability))

sd(t7$ability_est)
sd(t7$ability)




### Simulate multiple datasets with same theta

discr <- c(2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2)

#   nrow(diff)

for (i in 1:nrow(diff)) {
    if (i==1){
        diff1 <- diff[i,1]
    } else {
        diff1 <- c(diff1,diff[i,1])
    }
}


theta1 <- rnorm(500,0,1)
#   hist(theta1)


ds <- list()
for (i in 1:100){
    ds2 <- sim.poly.npn(nvar = 30 ,n = 500, a=discr,
        c=0,z=1,low=-2.0,high=2.0, mu=0, sd=1, cat=2,theta=theta1)
    # ds2 <- sim.poly.npn(nvar = 30 ,n = 500, a=discr,
    #     c=0,z=1,d=diff1, mu=0, sd=1, cat=2,theta=theta1)
    
    t1 <- data.frame(ds2[["items"]])
    t2 <- data.frame(ds2[["theta"]])
    
    names(t2) <- "ability"
    
    t2 <- cbind(t2,t1)
    
    ds[[i]] <- t2
}

sim_summ <- list()
for (i in 1 :length(ds)){
    mcal <- mirt(ds[[i]][,2:31],mdl,itemtype='2PL',pars=mcal_par)
    #   coef(mcal)
    fsc <- fscores(mcal,full.scores=TRUE, plausible.type="MH",plausible.draws=50)
    fsc_med <- data.frame(list(fsc[[5]],fsc[[10]],fsc[[15]],fsc[[20]],fsc[[25]],fsc[[30]],fsc[[35]],
                               fsc[[40]],fsc[[45]],fsc[[50]]))
    names(fsc_med) <- c("D1","D2","D3","D4","D5","D6","D7","D8","D9","D10")
    fsc_med$ability_est <- apply(fsc_med,1,mean)
    
    t6 <- cbind(t2,fsc_med)
    
    t6$resid <- t6$ability_est - t6$ability
    
    mean(t6$ability_est)
    sd(t6$ability_est)
    
 
    
    fsc <- data.frame(fscores(mcal,full.scores=TRUE, method="ML"))
    # fsc <- data.frame(fscores(mcal,full.scores=TRUE, method="plausible"))
    names(fsc) <- "ability_est"
    t7 <- cbind(t2,fsc)
    t7$ability_est <- ifelse(t7$ability_est==Inf,3,
        ifelse(t7$ability_est==-Inf,-3,t7$ability_est))
    t7$resid <- t7$ability_est - t7$ability
    
    mean(t7$ability_est)
    sd(t7$ability_est)
    
    
    t6$abil_est_st <- (t6$ability_est - mean(t6$ability_est))*(sd(t6$ability)/sd(t6$ability_est)) +
        mean(t6$ability) # linear equating to true ability metric
    
    t6$resid_st <- t6$abil_est_st - t6$ability
    t6$resid_st <- sqrt((t6$abil_est_st - t6$ability)^2)
    
    sqrt(mean(t6$resid_st^2))
    
    
    
    
    # plot(t6$resid ~ t6$ability)
    
    sim_summ[[i]] <- t6[,c("ability","ability_est","resid","abil_est_st","resid_st")]
    
}

abil <- data.frame(matrix(ncol = length(sim_summ), nrow = 500))
abil_est <- data.frame(matrix(ncol = length(sim_summ), nrow = 500))
res <- data.frame(matrix(ncol = length(sim_summ), nrow = 500))
rmse <- data.frame(matrix(nrow = 1,ncol = length(sim_summ)))
for (i in 1:length(sim_summ)) {
    nm <- paste("dset_",i,sep="")
    abil[,i] <- sim_summ[[i]]$ability
    names(abil)[i] <- nm
    abil_est[,i] <- sim_summ[[i]]$abil_est_st
    names(abil_est)[i] <- nm
    res[,i] <- sim_summ[[i]]$resid_st
    names(res)[i] <- nm
    rmse[1,i] <- sqrt(mean(sim_summ[[i]]$resid_st^2))
    names(rmse)[i] <- nm
}


abil_est$abil_est_m <- apply(abil_est,1,mean)
abil_est$abil_est_sd <- apply(abil_est,1,sd)
res$resid_m <- apply(res,1,mean)
res$resid_sd <- apply(res,1,sd)

summ <- data.frame(abil$dset_1,abil_est$abil_est_m,abil_est$abil_est_sd,
    res$resid_m,res$resid_sd)

plot(summ$res.resid_m ~ summ$abil.dset_1)
plot(summ$abil_est.abil_est_m ~ summ$abil.dset_1)
hist(summ$res.resid_sd)
hist(summ$res.resid_m)


summ <- summ[order(summ$abil.dset_1),]

### example of irtScore from psych
d9 <- sim.irt(30,1000,-2.5,2.5,mod="normal") #dichotomous items
test <- irt.fa(d9$items)
scores <- scoreIrt(test,d9$items)
scores.df <- data.frame(scores,true=d9$theta) #combine the estimates with the true thetas.

scores.df$resid <- scores.df$theta1 - scores.df$true

plot(scores.df$resid ~ scores.df$true)

### mirt based simulation

theta1 <- rnorm(500,0,1)
dataset1 <- simdata(disc, diff, 500, Theta=as.matrix(theta1),  itemtype = '2PL')

# mod <- mirt(dataset1, 1, method = 'MHRM')
# coef(mod)

dataset1 <- data.frame(cbind(theta1,dataset1))
# names(dataset1) <- gsub("Item_","V",names(dataset1))

# dataset1 <- cbind(theta1,dataset1)

mdl <- mirt.model('cog = 1-30')
mcal_par <- mirt(dataset1[,2:31],mdl,pars='values')

mcal_par[mcal_par$name == "a1","value"] <- disc
mcal_par[mcal_par$name == "a1","est"] <- FALSE
mcal_par[mcal_par$name == "d","value"] <- diff
mcal_par[mcal_par$name == "d","est"] <- FALSE
mcal_par[mcal_par$item == "GROUP","est"] <- TRUE


mcal <- mirt(dataset1[,2:31],mdl,itemtype='2PL',pars=mcal_par)
#   coef(mcal)
fsc <- fscores(mcal,full.scores=TRUE, plausible.type="MH",plausible.draws=50)
fsc_med <- data.frame(list(fsc[[5]],fsc[[10]],fsc[[15]],fsc[[20]],fsc[[25]],fsc[[30]],fsc[[35]],
                           fsc[[40]],fsc[[45]],fsc[[50]]))
names(fsc_med) <- c("D1","D2","D3","D4","D5","D6","D7","D8","D9","D10")
fsc_med$ability_est <- apply(fsc_med,1,median)


# fsc <- data.frame(fscores(mcal,full.scores=TRUE, method="MAP"))
# fsc <- data.frame(fscores(mcal,full.scores=TRUE, method="plausible"))
# names(fsc) <- "ability_est"
# t6 <- cbind(t2,fsc)

t6 <- cbind(theta1,fsc_med)

t6$resid <- t6$ability_est - t6$theta1


plot(t6$resid ~ t6$ability)

plot(t6$ability_est ~ t6$theta1)
t6$theta_st <- (t6$theta1 - mean(t6$theta1))/sd(t6$theta1)
t6$abil_est_st <- (t6$ability_est - mean(t6$ability_est))/sd(t6$ability_est)
t6$resid_st <- t6$abil_est_st - t6$theta_st

plot(t6$abil_est_st ~ t6$theta_st)
plot(t6$resid_st ~ t6$theta_st)
abline(lm(t6$resid_st ~ t6$theta_st))

summary(t6$theta1)
summary(t6$ability_est)
mean(t6$theta1)
mean(t6$ability_est)
sd(t6$theta1)
sd(t6$ability_est)

mcal_par[mcal_par$item == "GROUP" & mcal_par$name == "MEAN_1","est"] <- TRUE
mcal_par[mcal_par$item == "GROUP" & mcal_par$name == "COV_11","est"] <- TRUE



################################ Simulate IRT Models ###########################
mdl <- mirt.model('cog = 1-30')

#   mcal_par$item <- gsub("V","Item_",mcal_par$item)

mcal_par <- mirt(ds[[i]][,2:31],mdl,itemtype='2PL',pars='values')
mcal_par[mcal_par$name == "a1","value"] <- disc
mcal_par[mcal_par$name == "a1","est"] <- FALSE
mcal_par[mcal_par$name == "d","value"] <- diff
mcal_par[mcal_par$name == "d","est"] <- FALSE
mcal_par[mcal_par$name == "GROUP","est"] <- TRUE

mcal_par_1 <- mcal_par

mcal_par_1[!mcal_par_1$item %in% c("Item_14","Item_15","Item_16","Item_17"),"est"] <- TRUE

set.seed(21589)
for (j in 1:10) {
    time <- Sys.time()
    theta1 <- rnorm(500,0,1)
    #   hist(theta1)
    
    ds <- list()
    for (i in 1:100){
        ds2 <- data.frame(simdata(disc, diff, 500, Theta=as.matrix(theta1),
            itemtype = '2PL'))
        # ds2 <- sim.poly.npn(nvar = 30 ,n = 500, a=discr,
        #                     c=0,z=1,d=diff, mu=0, sd=1, cat=2,theta=theta1)
        # ds2 <- sim.poly.npn(nvar = 30 ,n = 500, a=discr,
        #     c=0,z=1,d=diff1, mu=0, sd=1, cat=2,theta=theta1)
        
        ds2 <- cbind(theta1,ds2)
        
        ds[[i]] <- ds2
    }
    
    sim_summ <- list()
    for (i in 1:length(ds)){
        mcal <- mirt(ds[[i]][,2:31],mdl,itemtype='2PL',pars=mcal_par_1)
        #   coef(mcal)
        fsc <- fscores(mcal,full.scores=TRUE, plausible.type="MH",plausible.draws=50)
        fsc_med <- data.frame(list(fsc[[5]],fsc[[10]],fsc[[15]],fsc[[20]],fsc[[25]],fsc[[30]],fsc[[35]],
            fsc[[40]],fsc[[45]],fsc[[50]]))
        names(fsc_med) <- c("D1","D2","D3","D4","D5","D6","D7","D8","D9","D10")
        fsc_med$ability_est <- apply(fsc_med,1,median)
        
        
        # fsc <- data.frame(fscores(mcal,full.scores=TRUE, method="MAP"))
        # fsc <- data.frame(fscores(mcal,full.scores=TRUE, method="plausible"))
        # names(fsc) <- "ability_est"
        # t6 <- cbind(t2,fsc)
        
        t6 <- cbind(ds[[i]],fsc_med)
        names(t6) <- sub("theta1","ability",names(t6))
        
        t6$resid <- t6$ability_est - t6$ability
        
        t6$abil_est_st <- (t6$ability_est - mean(t6$ability_est))*(sd(t6$ability)/sd(t6$ability_est)) +
            mean(t6$ability) # linear equating to true ability metric
        
        t6$resid_st <- t6$abil_est_st - t6$ability
        
        # plot(t6$resid ~ t6$ability)
        
        sim_summ[[i]] <- t6[,c("ability","ability_est","resid","abil_est_st","resid_st")]
        
    }
    
    # abil <- data.frame(matrix(ncol = length(sim_summ), nrow = 500))
    # abil_est <- data.frame(matrix(ncol = length(sim_summ), nrow = 500))
    # res <- data.frame(matrix(ncol = length(sim_summ), nrow = 500))
    rmse <- data.frame(matrix(nrow = 1,ncol = length(sim_summ)))
    for (i in 1:length(sim_summ)) {
        nm <- paste("dset_",i,sep="")
        # abil[,i] <- sim_summ[[i]]$ability
        # names(abil)[i] <- nm
        # abil_est[,i] <- sim_summ[[i]]$abil_est_st
        # names(abil_est)[i] <- nm
        # res[,i] <- sim_summ[[i]]$resid_st
        # names(res)[i] <- nm
        rmse[1,i] <- sqrt(mean(sim_summ[[i]]$resid_st^2))
        names(rmse)[i] <- nm
    } # end for i
    if (j==1) {
        rmse_summ <- rmse
    } else {
        rmse_summ <- rbind(rmse_summ,rmse)
    }
    cat(paste("End Iteration - ",j,", Elapsed time: ",Sys.time() - time,"\n",sep=""))
} # end for j

rmse_summ_1 <- rmse_summ
rmse_summ$rmse_avg <- apply(rmse_summ,1,mean)

mean(rmse_summ$rmse_avg)
mean(rmse_summ_1$rmse_avg)

summ1 <- equateSim(seed=21589,n_group=2,grp_mean=c(0.5,-0.5),grp_sd=c(1,1),
    n_samp=500,n_rep_theta=100,n_itm=30,pars=mcal_par_1,n_rep=10,
    itms1=item_list_1,itms2=item_list_2)

seed <- 21589
n_group <- 2
grp_mean <- c(0.5,-0.5)
grp_sd <- c(1,1)
n_samp <- 500
n_rep_theta <- 5
n_itm <- 30
fsc_method <- "EAP"
pars <- mcal_par_1

# Scenario 1

mcal_par_1 <- mcal_par
mcal_par_1[!mcal_par_1$item %in% c("Item_14","Item_15","Item_16","Item_17") & 
    mcal_par_1$name %in% c("a1","d"),"est"] <- TRUE

pars <- mcal_par_1
n_rep <- 2
item_list_1 <- c("Item_1","Item_2","Item_3","Item_4","Item_5","Item_6","Item_7",
    "Item_8","Item_9","Item_10","Item_11","Item_12","Item_13","Item_14",
    "Item_15","Item_16","Item_17","Item_18","Item_19","Item_20","Item_21",
    "Item_22","Item_23","Item_24","Item_25","Item_26",
    "Item_27","Item_28","Item_28","Item_29","Item_30")
item_list_1 <- c("Item_14","Item_15","Item_16","Item_17","Item_18","Item_19","Item_20",
    "Item_21","Item_22","Item_23","Item_24","Item_25","Item_26",
    "Item_27","Item_28","Item_28","Item_29","Item_30")
item_list_2 <- c("Item_1","Item_2","Item_3","Item_4","Item_5","Item_6","Item_7",
    "Item_8","Item_9","Item_10","Item_11","Item_12","Item_13","Item_14",
    "Item_15","Item_16","Item_17")

summ1 <- equateSim(seed=21589,n_group=2,grp_mean=c(0.5,-0.5),grp_sd=c(1,1),
    n_samp=500,n_rep_theta=100,n_itm=30,pars=mcal_par_2,n_rep=10,
    itms1=item_list_1,itms2=item_list_2)

summ1 <- equateSim(seed=21589,n_group=2,grp_mean=c(0.5,-0.5),grp_sd=c(1,1),
    n_samp=500,n_rep_theta=1,n_itm=30,pars=mcal_par_1,n_rep=100,
    itms1=item_list_1,itms2=item_list_2,fsc_method="EAP")

# summ1$rmse_avg <- apply(summ1[,1:(ncol(summ1)-2)],1,mean)
rmse_grp_1_1 <- mean(summ1[summ1$group == 1 & summ1$type == "raw","rmse_avg"])
rmse_grp_2_1 <- mean(summ1[summ1$group == 2 & summ1$type == "raw","rmse_avg"])
rmse_st_grp_1_1 <- mean(summ1[summ1$group == 1 & summ1$type == "standardized","rmse_avg"])
rmse_st_grp_2_1 <- mean(summ1[summ1$group == 2 & summ1$type == "standardized","rmse_avg"])

#   Scenario 2

mcal_par_2 <- mcal_par
mcal_par_2[!mcal_par_2$item %in% c("Item_4","Item_5","Item_6","Item_7") & 
    mcal_par_2$name %in% c("a1","d"),"est"] <- TRUE

item_list_1 <- c("Item_4","Item_5","Item_6","Item_7","Item_18","Item_19","Item_20","Item_21",
                 "Item_22","Item_23","Item_24","Item_25","Item_26",
                 "Item_27","Item_28","Item_28","Item_29","Item_30")
item_list_2 <- c("Item_1","Item_2","Item_3","Item_4","Item_5","Item_6","Item_7",
                 "Item_8","Item_9","Item_10","Item_11","Item_12","Item_13","Item_14",
                 "Item_15","Item_16","Item_17","Item_18")

itms1 <- item_list_1
itms2 <- item_list_2


summ2 <- equateSim(seed=21589,n_group=2,grp_mean=c(0.5,-0.5),grp_sd=c(1,1),
                   n_samp=500,n_rep_theta=100,n_itm=30,pars=mcal_par_2,n_rep=10,
                   itms1=item_list_1,itms2=item_list_2)

summ2 <- equateSim(seed=21589,n_group=2,grp_mean=c(0.5,-0.5),grp_sd=c(1,1),
                   n_samp=500,n_rep_theta=1,n_itm=30,pars=mcal_par_2,n_rep=100,
                   itms1=item_list_1,itms2=item_list_2,fsc_method="EAP")

# summ2$rmse_avg <- apply(summ2[,1:100],1,mean)
rmse_grp_1_2 <- mean(summ2[summ2$group == 1 & summ2$type == "raw","rmse_avg"])
rmse_grp_2_2 <- mean(summ2[summ2$group == 2 & summ2$type == "raw","rmse_avg"])
rmse_st_grp_1_2 <- mean(summ2[summ2$group == 1 & summ2$type == "standardized","rmse_avg"])
rmse_st_grp_2_2 <- mean(summ2[summ2$group == 2 & summ2$type == "standardized","rmse_avg"])

#   Scenario 3

mcal_par_3 <- mcal_par
mcal_par_3[!mcal_par_3$item %in% c("Item_13","Item_14","Item_15","Item_16",
    "Item_17","Item_18") & mcal_par_3$name %in% c("a1","d"),"est"] <- TRUE

item_list_1 <- c("Item_13","Item_14","Item_15","Item_16","Item_17","Item_18","Item_19",
    "Item_20","Item_21","Item_22","Item_23","Item_24","Item_25","Item_26",
    "Item_27","Item_28","Item_29","Item_30")
item_list_2 <- c("Item_1","Item_2","Item_3","Item_4","Item_5","Item_6","Item_7",
    "Item_8","Item_9","Item_10","Item_11","Item_12","Item_13","Item_14",
    "Item_15","Item_16","Item_17","Item_18","Item_19")

summ3 <- equateSim(seed=21589,n_group=2,grp_mean=c(0.5,-0.5),grp_sd=c(1,1),
                   n_samp=500,n_rep_theta=100,n_itm=30,pars=mcal_par_3,n_rep=10,
                   itms1=item_list_1,itms2=item_list_2)
summ3 <- equateSim(seed=21589,n_group=2,grp_mean=c(0.5,-0.5),grp_sd=c(1,1),
                   n_samp=500,n_rep_theta=1,n_itm=30,pars=mcal_par_3,n_rep=100,
                   itms1=item_list_1,itms2=item_list_2,fsc_method="EAP")

# summ3$rmse_avg <- apply(summ3[,1:100],1,mean)
rmse_grp_1_3 <- mean(summ3[summ3$group == 1 & summ3$type == "raw","rmse_avg"])
rmse_grp_2_3 <- mean(summ3[summ3$group == 2 & summ3$type == "raw","rmse_avg"])
rmse_st_grp_1_3 <- mean(summ3[summ3$group == 1 & summ3$type == "standardized","rmse_avg"])
rmse_st_grp_2_3 <- mean(summ3[summ3$group == 2 & summ3$type == "standardized","rmse_avg"])

#   Scenario 4

mcal_par_4 <- mcal_par
mcal_par_4[!mcal_par_4$item %in% c("Item_4","Item_5") & 
    mcal_par_4$name %in% c("a1","d"),"est"] <- TRUE

item_list_1 <- c("Item_4","Item_5","Item_15","Item_16","Item_17","Item_18","Item_19",
                 "Item_20","Item_21","Item_22","Item_23","Item_24","Item_25","Item_26",
                 "Item_27","Item_28","Item_29","Item_30")
item_list_2 <- c("Item_1","Item_2","Item_3","Item_4","Item_5","Item_6","Item_7",
                 "Item_8","Item_9","Item_10","Item_11","Item_12","Item_13","Item_14",
                 "Item_15","Item_16","Item_17","Item_18")

summ4 <- equateSim(seed=21589,n_group=2,grp_mean=c(0.5,-0.5),grp_sd=c(1,1),
                   n_samp=500,n_rep_theta=100,n_itm=30,pars=mcal_par_4,n_rep=10,
                   itms1=item_list_1,itms2=item_list_2)
summ4 <- equateSim(seed=21589,n_group=2,grp_mean=c(0.5,-0.5),grp_sd=c(1,1),
                   n_samp=500,n_rep_theta=1,n_itm=30,pars=mcal_par_4,n_rep=100,
                   itms1=item_list_1,itms2=item_list_2,fsc_method="EAP")

summ4$rmse_avg <- apply(summ4[,1:100],1,mean)
rmse_grp_1_4 <- mean(summ4[summ4$group == 1 & summ4$type == "raw","rmse_avg"])
rmse_grp_2_4 <- mean(summ4[summ4$group == 2 & summ4$type == "raw","rmse_avg"])
rmse_st_grp_1_4 <- mean(summ4[summ4$group == 1 & summ4$type == "standardized","rmse_avg"])
rmse_st_grp_2_4 <- mean(summ4[summ4$group == 2 & summ4$type == "standardized","rmse_avg"])


#   Scenario 5

mcal_par_5 <- mcal_par
mcal_par_5[!mcal_par_5$item %in% c("Item_14","Item_15","Item_16","Item_17") & 
             mcal_par_5$name %in% c("a1","d"),"est"] <- TRUE

item_list_1 <- c("Item_2","Item_4","Item_6","Item_8","Item_10","Item_12","Item_14","Item_15",
                 "Item_16","Item_17","Item_18","Item_20","Item_22","Item_24","Item_26",
                 "Item_28","Item_30")
item_list_2 <- c("Item_1","Item_3","Item_5","Item_7","Item_9","Item_11","Item_13","Item_14",
                 "Item_15","Item_16","Item_17","Item_19","Item_21","Item_23","Item_25","Item_27",
                 "Item_29")

summ5 <- equateSim(seed=21589,n_group=2,grp_mean=c(0.5,-0.5),grp_sd=c(1,1),
                   n_samp=500,n_rep_theta=100,n_itm=30,pars=mcal_par_5,n_rep=10,
                   itms1=item_list_1,itms2=item_list_2)
summ5 <- equateSim(seed=21589,n_group=2,grp_mean=c(0.5,-0.5),grp_sd=c(1,1),
                   n_samp=500,n_rep_theta=1,n_itm=30,pars=mcal_par_5,n_rep=100,
                   itms1=item_list_1,itms2=item_list_2,fsc_method="EAP")

rmse_grp_1_5 <- mean(summ5[summ5$group == 1 & summ5$type == "raw","rmse_avg"])
rmse_grp_2_5 <- mean(summ5[summ5$group == 2 & summ5$type == "raw","rmse_avg"])
rmse_st_grp_1_5 <- mean(summ5[summ5$group == 1 & summ5$type == "standardized","rmse_avg"])
rmse_st_grp_2_5 <- mean(summ5[summ5$group == 2 & summ5$type == "standardized","rmse_avg"])


#   Scenario 6

mcal_par_6 <- mcal_par
mcal_par_6[!mcal_par_6$item %in% c("Item_14","Item_15") & 
             mcal_par_6$name %in% c("a1","d"),"est"] <- TRUE

item_list_1 <- c("Item_2","Item_4","Item_6","Item_8","Item_10","Item_12","Item_14","Item_15",
                 "Item_16","Item_18","Item_20","Item_22","Item_24","Item_26",
                 "Item_28","Item_30")
item_list_2 <- c("Item_1","Item_3","Item_5","Item_7","Item_9","Item_11","Item_13","Item_14",
                 "Item_15","Item_17","Item_19","Item_21","Item_23","Item_25","Item_27",
                 "Item_29")

summ6 <- equateSim(seed=21589,n_group=2,grp_mean=c(0.5,-0.5),grp_sd=c(1,1),
                   n_samp=500,n_rep_theta=100,n_itm=30,pars=mcal_par_6,n_rep=10,
                   itms1=item_list_1,itms2=item_list_2)
summ6 <- equateSim(seed=21589,n_group=2,grp_mean=c(0.5,-0.5),grp_sd=c(1,1),
                   n_samp=500,n_rep_theta=1,n_itm=30,pars=mcal_par_6,n_rep=100,
                   itms1=item_list_1,itms2=item_list_2,fsc_method="EAP")

rmse_grp_1_6 <- mean(summ6[summ6$group == 1 & summ6$type == "raw","rmse_avg"])
rmse_grp_2_6 <- mean(summ6[summ6$group == 2 & summ6$type == "raw","rmse_avg"])
rmse_st_grp_1_6 <- mean(summ6[summ6$group == 1 & summ6$type == "standardized","rmse_avg"])
rmse_st_grp_2_6 <- mean(summ6[summ6$group == 2 & summ6$type == "standardized","rmse_avg"])


#   Scenario 7

mcal_par_7 <- mcal_par
mcal_par_7[!mcal_par_7$item %in% c("Item_12","Item_13","Item_14","Item_15","Item_16","Item_17","Item_18","Item_19") & 
             mcal_par_7$name %in% c("a1","d"),"est"] <- TRUE


item_list_1 <- c("Item_2","Item_4","Item_6","Item_8","Item_10","Item_12","Item_13","Item_14","Item_15",
                 "Item_16","Item_17","Item_18","Item_19","Item_20","Item_22","Item_24","Item_26",
                 "Item_28","Item_30")
item_list_2 <- c("Item_1","Item_3","Item_5","Item_7","Item_9","Item_11","Item_12","Item_13","Item_14",
                 "Item_15","Item_16","Item_17","Item_18","Item_19","Item_21","Item_23","Item_25","Item_27",
                 "Item_29")

summ7 <- equateSim(seed=21589,n_group=2,grp_mean=c(0.5,-0.5),grp_sd=c(1,1),
                   n_samp=500,n_rep_theta=100,n_itm=30,pars=mcal_par_7,n_rep=10,
                   itms1=item_list_1,itms2=item_list_2)
summ7 <- equateSim(seed=21589,n_group=2,grp_mean=c(0.5,-0.5),grp_sd=c(1,1),
                   n_samp=500,n_rep_theta=1,n_itm=30,pars=mcal_par_7,n_rep=100,
                   itms1=item_list_1,itms2=item_list_2,fsc_method="EAP")

rmse_grp_1_7 <- mean(summ7[summ7$group == 1 & summ7$type == "raw","rmse_avg"])
rmse_grp_2_7 <- mean(summ7[summ7$group == 2 & summ7$type == "raw","rmse_avg"])
rmse_st_grp_1_7 <- mean(summ7[summ7$group == 1 & summ7$type == "standardized","rmse_avg"])
rmse_st_grp_2_7 <- mean(summ7[summ7$group == 2 & summ7$type == "standardized","rmse_avg"])


#   Scenario 8

mcal_par_8 <- mcal_par
mcal_par_8[!mcal_par_8$item %in% c("Item_2","Item_4","Item_6","Item_8","Item_10","Item_12","Item_13",
    "Item_14","Item_15","Item_16","Item_17","Item_18","Item_19","Item_20","Item_22","Item_24",
    "Item_26","Item_28","Item_30") & 
    mcal_par_8$name %in% c("a1","d"),"est"] <- TRUE


item_list_1 <- c("Item_2","Item_4","Item_6","Item_8","Item_10","Item_12","Item_13","Item_14","Item_15",
                 "Item_16","Item_17","Item_18","Item_19","Item_20","Item_22","Item_24","Item_26",
                 "Item_28","Item_30")
item_list_2 <- c("Item_2","Item_4","Item_6","Item_8","Item_10","Item_12","Item_13","Item_14","Item_15",
                 "Item_16","Item_17","Item_18","Item_19","Item_20","Item_22","Item_24","Item_26",
                 "Item_28","Item_30")

mcal_par_8 <- mcal_par_8[mcal_par_8$item %in% c("GROUP",item_list_1),]
mcal_par_8$parnum <- 1:nrow(mcal_par_8)

mdl1 <- mirt.model('cog = 1-19')

pars <- mcal_par_8

summ8 <- equateSim(seed=21589,n_group=2,grp_mean=c(0.5,-0.5),grp_sd=c(1,1),
                   n_samp=500,n_rep_theta=1,n_itm=19,pars=mcal_par_8,n_rep=100,
                   itms1=item_list_1,itms2=item_list_2,fsc_method="EAP",model=mdl1)

summ8 <- equateSim(seed=21589,grp_mean=c(0.5,-0.5),grp_sd=c(1,1),
                   n_samp=500,n_rep_theta=1,n_itm=19,pars=mcal_par_8,n_rep=100,
                   itms1=item_list_1,itms2=item_list_2,fsc_method="EAP")


rmse_grp_1_8 <- mean(summ8[summ8$group == 1 & summ8$type == "raw" & 
  summ8$statistic == "rmse","avg"])
rmse_grp_2_8 <- mean(summ8[summ8$group == 2 & summ8$type == "raw" & 
  summ8$statistic == "rmse","avg"])
rmse_st_grp_1_8 <- mean(summ8[summ8$group == 1 & summ8$type == "standardized" & 
  summ8$statistic == "rmse","avg"])
rmse_st_grp_2_8 <- mean(summ8[summ8$group == 2 & summ8$type == "standardized" & 
  summ8$statistic == "rmse","avg"])
mean_grp_1_8 <- mean(summ8[summ8$group == 1 & summ8$type == "raw" & 
  summ8$statistic == "est_mean","avg"])
mean_grp_2_8 <- mean(summ8[summ8$group == 2 & summ8$type == "raw" & 
  summ8$statistic == "est_mean","avg"])
mean_st_grp_1_8 <- mean(summ8[summ8$group == 1 & summ8$type == "standardized" & 
  summ8$statistic == "est_mean","avg"])
mean_st_grp_2_8 <- mean(summ8[summ8$group == 2 & summ8$type == "standardized" & 
  summ8$statistic == "est_mean","avg"])
sd_grp_1_8 <- mean(summ8[summ8$group == 1 & summ8$type == "raw" & 
  summ8$statistic == "est_sd","avg"])
sd_grp_2_8 <- mean(summ8[summ8$group == 2 & summ8$type == "raw" & 
  summ8$statistic == "est_sd","avg"])
sd_st_grp_1_8 <- mean(summ8[summ8$group == 1 & summ8$type == "standardized" & 
  summ8$statistic == "est_sd","avg"])
sd_st_grp_2_8 <- mean(summ8[summ8$group == 2 & summ8$type == "standardized" & 
  summ8$statistic == "est_sd","avg"])


#   Scenario 9

mcal_par_9 <- mcal_par

item_list_1 <- c("Item_1","Item_2","Item_3","Item_4","Item_5","Item_6","Item_7",
                 "Item_8","Item_9","Item_10","Item_11","Item_12","Item_13","Item_14",
                 "Item_15","Item_16","Item_17","Item_18","Item_19","Item_20","Item_21",
                 "Item_22","Item_23","Item_24","Item_25","Item_26",
                 "Item_27","Item_28","Item_28","Item_29","Item_30")
item_list_2 <- c("Item_1","Item_2","Item_3","Item_4","Item_5","Item_6","Item_7",
                 "Item_8","Item_9","Item_10","Item_11","Item_12","Item_13","Item_14",
                 "Item_15","Item_16","Item_17","Item_18","Item_19","Item_20","Item_21",
                 "Item_22","Item_23","Item_24","Item_25","Item_26",
                 "Item_27","Item_28","Item_28","Item_29","Item_30")

mcal_par_9 <- mcal_par_9[mcal_par_9$item %in% c("GROUP",item_list_1),]
mcal_par_9$parnum <- 1:nrow(mcal_par_9)

mdl <- mirt.model('cog = 1-30')

pars <- mcal_par_9

grp_mean <- c(0.5,-0.5)
grp_sd <- c(1,1)
n_samp <- 500
n_rep_theta <- 1
n_itm <- 19
pars <- mcal_par_8
n_rep <- 100
itms1 <- item_list_1
itms2 <- item_list_2
fsc_method <- "EAP"
model <- mdl1

summ9 <- equateSim(seed=21589,grp_mean=c(0,0),grp_sd=c(1,1),
                   n_samp=250,n_rep_theta=1,n_itm=30,pars=mcal_par_9,n_rep=100,
                   itms1=item_list_1,itms2=item_list_2,fsc_method="plausible")

summ9pv <- summ9

summ9eap <- equateSim(seed=21589,grp_mean=c(0,0),grp_sd=c(1,1),
                   n_samp=250,n_rep_theta=1,n_itm=30,pars=mcal_par_9,n_rep=100,
                   itms1=item_list_1,itms2=item_list_2,fsc_method="EAP")

summ9map <- equateSim(seed=21589,grp_mean=c(0,0),grp_sd=c(1,1),
                      n_samp=250,n_rep_theta=1,n_itm=30,pars=mcal_par_9,n_rep=100,
                      itms1=item_list_1,itms2=item_list_2,fsc_method="MAP")

summ9ml <- equateSim(seed=21589,grp_mean=c(0,0),grp_sd=c(1,1),
                      n_samp=250,n_rep_theta=1,n_itm=30,pars=mcal_par_9,n_rep=100,
                      itms1=item_list_1,itms2=item_list_2,fsc_method="ML")


rmse_grp_1_9 <- mean(summ9[summ9$group == 1 & summ9$type == "raw" & 
    summ9$statistic == "rmse","avg"])
rmse_grp_2_9 <- mean(summ9[summ9$group == 2 & summ9$type == "raw" & 
    summ9$statistic == "rmse","avg"])
rmse_st_grp_1_9 <- mean(summ9[summ9$group == 1 & summ9$type == "standardized" & 
    summ9$statistic == "rmse","avg"])
rmse_st_grp_2_9 <- mean(summ9[summ9$group == 2 & summ9$type == "standardized" & 
    summ9$statistic == "rmse","avg"])
mean_grp_1_9 <- mean(summ9[summ9$group == 1 & summ9$type == "raw" & 
    summ9$statistic == "est_mean","avg"])
mean_grp_2_9 <- mean(summ9[summ9$group == 2 & summ9$type == "raw" & 
    summ9$statistic == "est_mean","avg"])
mean_st_grp_1_9 <- mean(summ9[summ9$group == 1 & summ9$type == "standardized" & 
    summ9$statistic == "est_mean","avg"])
mean_st_grp_2_9 <- mean(summ9[summ9$group == 2 & summ9$type == "standardized" & 
    summ9$statistic == "est_mean","avg"])
sd_grp_1_9 <- mean(summ9[summ9$group == 1 & summ9$type == "raw" & 
    summ9$statistic == "est_sd","avg"])
sd_grp_2_9 <- mean(summ9[summ9$group == 2 & summ9$type == "raw" & 
    summ9$statistic == "est_sd","avg"])
sd_st_grp_1_9 <- mean(summ9[summ9$group == 1 & summ9$type == "standardized" & 
    summ9$statistic == "est_sd","avg"])
sd_st_grp_2_9 <- mean(summ9[summ9$group == 2 & summ9$type == "standardized" & 
    summ9$statistic == "est_sd","avg"])

summaryStats <- function(df) {
  rmse_grp_1 <- mean(df[df$group == 1 & df$type == "raw" & 
      df$statistic == "rmse","avg"])
  rmse_grp_2 <- mean(df[df$group == 2 & df$type == "raw" & 
      df$statistic == "rmse","avg"])
  rmse_st_grp_1 <- mean(df[df$group == 1 & df$type == "standardized" & 
      df$statistic == "rmse","avg"])
  rmse_st_grp_2 <- mean(df[df$group == 2 & df$type == "standardized" & 
      df$statistic == "rmse","avg"])
  mean_grp_1 <- mean(df[df$group == 1 & df$type == "raw" & 
      df$statistic == "est_mean","avg"])
  mean_grp_2 <- mean(df[df$group == 2 & df$type == "raw" & 
      df$statistic == "est_mean","avg"])
  mean_st_grp_1 <- mean(df[df$group == 1 & df$type == "standardized" & 
      df$statistic == "est_mean","avg"])
  mean_st_grp_2 <- mean(df[df$group == 2 & df$type == "standardized" & 
      df$statistic == "est_mean","avg"])
  sd_grp_1 <- mean(df[df$group == 1 & df$type == "raw" & 
      df$statistic == "est_sd","avg"])
  sd_grp_2 <- mean(df[df$group == 2 & df$type == "raw" & 
      df$statistic == "est_sd","avg"])
  sd_st_grp_1 <- mean(df[df$group == 1 & df$type == "standardized" & 
      df$statistic == "est_sd","avg"])
  sd_st_grp_2 <- mean(df[df$group == 2 & df$type == "standardized" & 
      df$statistic == "est_sd","avg"])
  rmse <- mean(rmse_grp_1,rmse_grp_2)
  rmse_st <- mean(rmse_st_grp_1,rmse_st_grp_2)
  mean <- mean(mean_grp_1,mean_grp_2)
  mean_st <- mean(mean_st_grp_1,mean_st_grp_2)
  sd <- mean(sd_grp_1,sd_grp_2)
  sd_st <- mean(sd_st_grp_1,sd_st_grp_2)
  
  stat <- data.frame(cbind(rmse,rmse_st,mean,mean_st,sd,sd_st))
  
}

df <- summ9eap

stateap <- summaryStats(summ9eap)
stateap$method <- "EAP"
statmap <- summaryStats(summ9map)
statmap$method <- "MAP"
statpv <- summaryStats(summ9pv)
statpv$method <- "plausible"
statml <- summaryStats(summ9ml)
statml$method <- "ML"

statsumm <- rbind(stateap,statmap,statpv,statml)

equateSim <- function(seed=NULL,n_group=2,grp_mean=c(0.5,-0.5),grp_sd=c(1,1),
    n_samp=500,n_rep_theta=1,n_itm=30,pars=mcal_pars,n_rep=100,
    itms1=item_list_1,itms2=item_list_2,fsc_method="EAP",model=mdl) {
    diff <- pars[pars$name == "d","value"]
    disc <- pars[pars$name == "a1","value"]
    set.seed(seed)
    for (j in 1:n_rep) {
        time <- Sys.time()
        theta1 <- data.frame(rnorm(n_samp,grp_mean[1],grp_sd[1]))
        names(theta1) <- "theta1"
        theta1$group <- 1
        theta2 <- data.frame(rnorm(n_samp,grp_mean[2],grp_sd[2]))
        names(theta2) <- "theta1"
        theta2$group <- 2
        theta1 <- rbind(theta1,theta2)
        
        
        #   hist(theta1)
        
        ds <- list()
        for (i in 1:n_rep_theta){
            ds2 <- data.frame(simdata(disc, diff, n_samp, Theta=as.matrix(theta1$theta1),
                itemtype = '2PL'))
            names(ds2) <- unique(pars[!pars$item == "GROUP","item"])
            # ds2[1:n_samp,"group"] <- 1
            # ds2[(n_samp+1):nrow(ds2),"group"] <- 2
            # ds2 <- sim.poly.npn(nvar = 30 ,n = 500, a=discr,
            #                     c=0,z=1,d=diff, mu=0, sd=1, cat=2,theta=theta1)
            # ds2 <- sim.poly.npn(nvar = 30 ,n = 500, a=discr,
            #     c=0,z=1,d=diff1, mu=0, sd=1, cat=2,theta=theta1)
            
            ds2 <- cbind(theta1,ds2)
            
            for(group in 1:2) {
                if (group == 1) {
                    ds2[ds2$group == 1,!names(ds2) %in% c("theta1","group",itms1)] <- NA
                } else {
                    ds2[ds2$group == 2,!names(ds2) %in% c("theta1","group",itms2)] <- NA
                }
            }
            
            
            ds[[i]] <- ds2
        }
        
        sim_summ <- list()
        for (i in 1:length(ds)){
            mcal <- mirt(ds[[i]][,3:(n_itm+2)],model,itemtype='2PL',pars=pars)
            #   coef(mcal)
            
            # if (fsc_method=="EAP") {
            #   fsc <- fscores(mcal,full.scores=TRUE, method="EAP")
            #   t6 <- cbind(ds[[i]],fsc)
            #   names(t6) <- sub("F1","ability_est",names(t6))
            #   names(t6) <- sub("theta1","ability",names(t6))
            # } else {
            #   fsc <- fscores(mcal,full.scores=TRUE, plausible.type="MH",plausible.draws=50)
            #   fsc_med <- data.frame(list(fsc[[5]],fsc[[10]],fsc[[15]],fsc[[20]],fsc[[25]],fsc[[30]],fsc[[35]],
            #                              fsc[[40]],fsc[[45]],fsc[[50]]))
            #   names(fsc_med) <- c("D1","D2","D3","D4","D5","D6","D7","D8","D9","D10")
            #   fsc_med$ability_est <- apply(fsc_med,1,median)
            #   
            #   t6 <- cbind(ds[[i]],fsc_med)
            #   names(t6) <- sub("theta1","ability",names(t6))
            # }

            fsc <- data.frame(fscores(mcal,full.scores=TRUE, method=fsc_method))
            names(fsc) <- "ability_est"
            t6 <- cbind(ds[[i]],fsc)
            # names(t6) <- sub("F1","ability_est",names(t6))
            names(t6) <- sub("theta1","ability",names(t6))
            
            t6$resid <- t6$ability_est - t6$ability
            
            t6$abil_est_st <- (t6$ability_est - mean(t6$ability_est))*(sd(t6$ability)/sd(t6$ability_est)) +
                mean(t6$ability) # linear equating to true ability metric
            
            t6$resid_st <- t6$abil_est_st - t6$ability
            
            # plot(t6$resid ~ t6$ability)
            
            sim_summ[[i]] <- t6[,c("ability","ability_est","resid","abil_est_st","resid_st","group")]
            
        }
        
        # abil <- data.frame(matrix(ncol = length(sim_summ), nrow = 500))
        # abil_est <- data.frame(matrix(ncol = length(sim_summ), nrow = 500))
        # res <- data.frame(matrix(ncol = length(sim_summ), nrow = 500))
        stat <- data.frame(matrix(nrow = 1,ncol = length(sim_summ)))
        for (i in 1:length(sim_summ)) {
            nm <- paste("dset_",i,sep="")
            # abil[,i] <- sim_summ[[i]]$ability
            # names(abil)[i] <- nm
            # abil_est[,i] <- sim_summ[[i]]$abil_est_st
            # names(abil_est)[i] <- nm
            # res[,i] <- sim_summ[[i]]$resid_st
            # names(res)[i] <- nm
            df <- sim_summ[[i]]
            stat[1,i] <- sqrt(mean(df[df$group==1,"resid"]^2))
            stat[2,i] <- sqrt(mean(df[df$group==2,"resid"]^2))
            stat[3,i] <- sqrt(mean(df[df$group==1,"resid_st"]^2))
            stat[4,i] <- sqrt(mean(df[df$group==2,"resid_st"]^2))
            stat[5,i] <- mean(df[df$group==1,"ability_est"])
            stat[6,i] <- mean(df[df$group==2,"ability_est"])
            stat[7,i] <- mean(df[df$group==1,"abil_est_st"])
            stat[8,i] <- mean(df[df$group==2,"abil_est_st"])
            stat[9,i] <- sd(df[df$group==1,"ability_est"])
            stat[10,i] <- sd(df[df$group==2,"ability_est"])
            stat[11,i] <- sd(df[df$group==1,"abil_est_st"])
            stat[12,i] <- sd(df[df$group==2,"abil_est_st"])
            names(stat)[i] <- nm
            stat[1,"group"] <- 1
            stat[2,"group"] <- 2
            stat[3,"group"] <- 1
            stat[4,"group"] <- 2
            stat[5,"group"] <- 1
            stat[6,"group"] <- 2
            stat[7,"group"] <- 1
            stat[8,"group"] <- 2
            stat[9,"group"] <- 1
            stat[10,"group"] <- 2
            stat[11,"group"] <- 1
            stat[12,"group"] <- 2
            stat[1,"type"] <- "raw"
            stat[2,"type"] <- "raw"
            stat[3,"type"] <- "standardized"
            stat[4,"type"] <- "standardized"
            stat[5,"type"] <- "raw"
            stat[6,"type"] <- "raw"
            stat[7,"type"] <- "standardized"
            stat[8,"type"] <- "standardized"
            stat[9,"type"] <- "raw"
            stat[10,"type"] <- "raw"
            stat[11,"type"] <- "standardized"
            stat[12,"type"] <- "standardized"
            stat[1,"statistic"] <- "rmse"
            stat[2,"statistic"] <- "rmse"
            stat[3,"statistic"] <- "rmse"
            stat[4,"statistic"] <- "rmse"
            stat[5,"statistic"] <- "est_mean"
            stat[6,"statistic"] <- "est_mean"
            stat[7,"statistic"] <- "est_mean"
            stat[8,"statistic"] <- "est_mean"
            stat[9,"statistic"] <- "est_sd"
            stat[10,"statistic"] <- "est_sd"
            stat[11,"statistic"] <- "est_sd"
            stat[12,"statistic"] <- "est_sd"
        } # end for i
        if (j==1) {
            stat_summ <- stat
        } else {
            stat_summ <- rbind(stat_summ,stat)
        }
        cat(paste("Iteration - ",j,", Elapsed time: ",Sys.time() - time,"\n",sep=""))
    } # end for j
    if(n_rep_theta == 1){
      stat_summ$avg <- stat_summ$dset_1
    } else {
      stat_summ$avg <- apply(stat_summ[,1:n_rep],1,mean)
    }
    return(stat_summ)
}

infoSim <- function(seed=NULL,n_group=2,grp_mean=c(0.5,-0.5),grp_sd=c(1,1),
    n_samp=500,pars=mcal_par,itms1=item_list_1,itms2=item_list_2) {
  diff <- pars[pars$name == "d","value"]
  disc <- pars[pars$name == "a1","value"]
  n_itm <- nrow(pars[pars$name == 'a1',])
  links <- intersect(itms1,itms2)
  set.seed(seed)
  
  theta1 <- data.frame(rnorm(n_samp,grp_mean[1],grp_sd[1]))
  names(theta1) <- "theta1"
  theta1$group <- 1
  theta2 <- data.frame(rnorm(n_samp,grp_mean[2],grp_sd[2]))
  names(theta2) <- "theta1"
  theta2$group <- 2
  theta1 <- rbind(theta1,theta2)
  
  ds <- data.frame(simdata(disc, diff, n_samp, Theta=as.matrix(theta1$theta1),
      itemtype = '2PL'))
  names(ds) <- unique(pars[!pars$item == "GROUP","item"])
  ds <- cbind(theta1,ds)
  
  for(group in 1:2) {
    if (group == 1) {
      ds[ds$group == 1,!names(ds) %in% c("theta1","group",itms1)] <- NA
    } else {
      ds[ds$group == 2,!names(ds) %in% c("theta1","group",itms2)] <- NA
    }
  }
  
  pars1 <- pars[pars$item %in% c(itms1,"GROUP"),]
  pars1$parnum <- 1:nrow(pars1)
  pars1[pars1$name=="d" & ! pars1$item %in% links,"est"] <- TRUE
  pars1[pars1$name=="a1" & ! pars2$item %in% links,"est"] <- TRUE
  
  pars2 <- pars[pars$item %in% c(itms2,"GROUP"),]
  pars2$parnum <- 1:nrow(pars2)
  pars2[pars2$name=="d","est"] <- TRUE
  pars2[pars2$name=="a1","est"] <- TRUE
  
  n_itm_1 <- nrow(pars1[pars1$name == 'a1',])
  n_itm_2 <- nrow(pars2[pars2$name == 'a1',])
  
  mdl1 <- mirt.model(paste("cog = 1-",n_itm_1,sep=""))
  mdl2 <- mirt.model(paste("cog = 1-",n_itm_2,sep=""))

  
  Theta <- matrix(seq(-4,4, by = .1))
  
  mcal1 <- mirt(ds[ds$group %in% 1,names(ds) %in% itms1],model=mdl1,
      modelitemtype='2PL',pars=pars1)
  mcal2 <- mirt(ds[ds$group %in% 2,names(ds) %in% itms2],model=mdl2,
      modelitemtype='2PL',pars=pars2)
  tinfo1 <- testinfo(mcal1,Theta)
  tinfo2 <- testinfo(mcal2,Theta)
  tinfo <- data.frame(cbind(Theta,tinfo1,tinfo2))
  names(tinfo) <- c("ability","info_1","info_2")
  
  return(tinfo)
  
  # plot(tinfo$info_1 ~ tinfo$ability, type="l",col="blue")
  # lines(tinfo$info_2 ~ tinfo$ability,col="green")
}


tinfo <- infoSim(seed=21589,n_group=2,grp_mean=c(0.5,-0.5),grp_sd=c(1,1),
    n_samp=500,pars=mcal_par,itms1=item_list_1,itms2=item_list_2)

seed<-21589
n_group<-2
grp_mean<-c(0.5,-0.5)
grp_sd<-c(1,1)
n_samp<-500
pars<-mcal_par
n_rep<-100
itms1<-item_list_1
itms2<-item_list_2
model<-mdl

item_list_1 <- c("Item_2","Item_4","Item_6","Item_8","Item_10","Item_12","Item_14","Item_15",
    "Item_16","Item_17","Item_18","Item_20","Item_22","Item_24","Item_26",
    "Item_28","Item_30")
item_list_2 <- c("Item_1","Item_3","Item_5","Item_7","Item_9","Item_11","Item_13","Item_14",
      "Item_15","Item_16","Item_17","Item_19","Item_21","Item_23","Item_25","Item_27",
      "Item_29")

itms1 <- item_list_1
itms2 <- item_list_2

plot(tinfo$info_1 ~ tinfo$ability, type="l",col="blue")
lines(tinfo$info_2 ~ tinfo$ability,col="green")

tinfo2 <- tinfo %>% gather(key,value,c("info_1","info_2"))
tinfo2$key <- ifelse(tinfo2$key %in% "info_1","Group 1","Group 2")

ggplot(data=tinfo2,aes(x=ability, y=value, color=key)) + geom_line(lwd=1) + xlab("Information")   

#   ----------------------------- End Simulate IRT -----------------------------

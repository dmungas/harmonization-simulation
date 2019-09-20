#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(mirt)
library(dplyr)
library(listviewer)
library(reactR)

# Define UI for application that draws a histogram
ui <- fluidPage(
    
    # Application title
    titlePanel("equateSim"),
    
    
    # Sidebar with a slider input for number of bins 
    sidebarLayout(
        sidebarPanel(
            #h1("If a mirt parameter file is not uploaded, IRT parameters will be simulated."),
            #h2("Select mirt parameter file here:"),
            h2("Optional part"),
            h6("If these are blank, 2PL IRT parameters will be selected randomly."),
            fileInput("parFile", "Upload mirt parameter file in .Rds format", accept = ".Rds"),
            #h2("If not uploading mirt parameter file, enter parameters here (only 2PL available here):"),
            #h3("If left blank, random parameter values will be chosen."),
            
            #h2("Select mirt model results object from mirt analysis here:"),
            #h3("This must be from the same mirt model as parameters file above."),
            fileInput("mmrFile", "Upload mirt model results object in .Rds format", accept = ".Rds"),
            #h2("If not uploading mirt model results file, enter parameters here (only 2PL available here):"),
            #h3("If left blank, random parameter values will be chosen."),
            h3("or (2PL only)"),
            textInput("alphas", "Enter discrimination parameters (separate with commas, spaces, and/or colons):"),
            textInput("betas", "Enter difficulty parameters (separate with commas, spaces, and/or colons):"),
            
            tags$hr(style="border-color: black;"),
            
            h2("Enter simulation rules here:"),
            numericInput("g1m", "True ability mean of group 1:", .5, min = -6, max = 6),
            numericInput("g1s", "True ability standard deviation of group 1:", 1, min = .01, max = 15),
            numericInput("g2m", "True ability mean of group 2:", -.5, min = -6, max = 6),
            numericInput("g2s", "True ability standard deviation of group 2:", 1, min = .01, max = 15),
            numericInput("n_samp", "Sample size of each group in each simulation:", 500, min = 1, max = 1e06),
            numericInput("n_rep_theta", 
                         "Number of repetitions of each simulated sample of true ability values", 
                         1, min = 1, max = 10),
            numericInput("n_itm", "Number of items:", 30, min = 1, max = 1000),
            numericInput("n_rep", "Number of simulated samples of true ability values:", 100, min = 1, max = 1e06),
            textInput("item_list_1", "List of items available for group 1 (separate with commas, spaces, and/or colons):",
                      "1:30"),
            textInput("item_list_2", "List of items available for group 2 (separate with commas, spaces, and/or colons):",
                      "1:30"),
            selectInput("fsc_method", "mirt method for calculating estimated ability values:", 
                        choices = c("EAP", "MAP", "ML", "WLE", "EAPsum", "plausible")),
            numericInput("seed", "Random seed for simulations:", 90210, min = 0, max = 1e10),
            actionButton("runSim", "Run Simulation")),
        
        # Show a table of the results summary
        mainPanel(
            tableOutput("equateSimResults"),
            downloadButton("exportButton1CSV", "Export Simulated Data Summary as CSV file"),
            downloadButton("exportButton1Rds", "Export Simulated Data Summary as Rds file"),
            br(), br(),
            downloadButton("exportButton2CSV", "Export All Simulated Data as CSV file"),
            downloadButton("exportButton2Rds", "Export All Simulated Data as Rds file"),
            h3("Summary of simulation rules:"),
            reactjsonOutput("parSummary", height = '100%')
        )
    ))

# Define server logic required to draw a histogram
server <- function(input, output, session) {
    
    # equateSim is a function the uses R mirt item response theory (IRT) methods
    # to: 1) simulate item level datasets based on item parameters, 2) perform IRT
    # calibrations on the simulated datasets, 3) estimate ability scores for each
    # person in the simulated datasets, and 4) compute metrics to quantify the 
    # correspondence between the estimated ability from the simulated datasets and the 
    # true ability that was used to generate the simulations. Metrics include the
    # root mean square error (RMSE), the mean of the estimated ability scores, and
    # the standard deviation of the estimated ability scores. Multiple simulations 
    # of the same set of simulated true ability values can be generated (n_rep_theta).
    #   This function simulates two groups that can differ in mean and standard
    # deviation of simulated true ability, and different items can be used for each
    # group as long as there is at least one common, linking item.
    # 
    # Parameters:
    #   seed - random seed for simulations (default=NULL)
    #   grp_mean - true ability means of the two groups (default=c(0.5,-0.5))
    #   grp_sd - true ability standard deviations of the two groups (default=c(1,1))
    #   n_rep - number of simulated samples of true ability values (default=100)
    #   n_rep_theta - number of repetitions of each simulated sample of true ability 
    #     values (default=1)
    #   n_samp - sample size of each group in each simulation (default=500)
    #   n_itm - number of items (default=30)
    #   pars - mirt item parameters file (default=mcal_par)
    #     A parameter file can be generated by:
    #       mcal_par <- mirt(df,mdl,pars='values') where df is a dataframe with
    #         item response data and mdl is a mirt model object. The returned file 
    #         (mcal_par in this example) can be edited to change item parameters.
    #   n_rep - number of simulated samples of true ability values (default=100)
    #   itms1 - list of items available for group 1 (default=item_list_1)
    #   itms2 - list of items available for group 2 (default=item_list_2)
    #   fsc_method - mirt method for calculating estimated ability values 
    #     (default="EAP")
    #   mod_res_obj - mirt model results object from mirt analysis. This must
    #     be from the same mirt model as pars.
    
    # equateSim returns a dataframe that has summary statistics (rmse, mean 
    # estimated ability, sd estimated ability) for each group and each simulated 
    # sample of true abilities (n_rep), and ge]nerates these statistics for raw 
    # estimated abilities and for estimated abilities on a metric transformed to 
    # match the true ability metric.
    
    xmpl <- simdata(a = runif(30, 0, 3), d = rnorm(30), N = 30, itemtype = "dich")
    xmod <- mirt.model(paste0("F = 1-", 30))
    mcal_parsD <- mirt(xmpl, xmod, pars = "values")
    
    
    equateSim <- function(seed=100,grp_mean=c(0.5,-0.5),grp_sd=c(1,1),
                          n_samp=500,n_rep_theta=1,n_itm=30,pars=mcal_parsD,n_rep=100,
                          itms1=1:20,itms2=10:30,fsc_method="EAP",
                          mod_res_obj=NULL, save_sims=TRUE) {
        require(dplyr)
        #modelSyntax <- paste0("F = 1-", n_itm)
        extractDiffMatrix <- function(pars){
            diff1 <- pars[grepl("d+",pars$name),c("item","name","value")]
            diff1$item <- as.character(diff1$item)
            diffsum <- diff1 %>% group_by(item) %>% summarise(
                ncat = sum(!is.na(item))
            )
            diff <- matrix(nrow=nrow(diffsum),ncol=max(diffsum$ncat))
            for (j in 1:nrow(diffsum)) {
                vals <- diff1[diff1$item %in% diffsum[j,"item"],"value"]
                if(length(vals) < max(diffsum$ncat)){
                    for (i in (length(vals)+1):max(diffsum$ncat)){
                        vals[i] <- NA
                    }
                }
                diff[j,] <- vals
            }
            return(diff)
        }
        
        diff <- extractDiffMatrix(pars)
        disc <- pars[pars$name == "a1","value"]
        itemtype <- as.character(pars[pars$name =="a1","class"])
        model <- mirt.model(paste("cog = 1-",n_itm,sep=""))
        set.seed(seed)
        j <- 0
        while (j < n_rep) {
            j <- j+1
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
                if (is.null(mod_res_obj)) {
                    ds2 <- data.frame(simdata(a=disc,d=diff,N=n_samp,Theta=as.matrix(theta1$theta1),
                                              itemtype=itemtype))
                } else {
                    ds2 <- data.frame(simdata(model=mod_res_obj,N=n_samp,
                                              Theta=as.matrix(theta1$theta1)))
                }
                
                names(ds2) <- unique(pars[!pars$item == "GROUP","item"])
                # ds2[1:n_samp,"group"] <- 1
                # ds2[(n_samp+1):nrow(ds2),"group"] <- 2
                # ds2 <- sim.poly.npn(nvar = 30 ,n = 500, a=discr,
                #                     c=0,z=1,d=diff, mu=0, sd=1, cat=2,theta=theta1)
                # ds2 <- sim.poly.npn(nvar = 30 ,n = 500, a=discr,
                #     c=0,z=1,d=diff1, mu=0, sd=1, cat=2,theta=theta1)
                
                ds2 <- cbind(theta1,ds2)
                
                i1names <- paste0("Item_", itms1)
                i2names <- paste0("Item_", itms2)
                
                for(group in 1:2) {
                    if (group == 1) {
                        ds2[ds2$group == 1,!names(ds2) %in% c("theta1","group",i1names)] <- NA
                    } else {
                        ds2[ds2$group == 2,!names(ds2) %in% c("theta1","group",i2names)] <- NA
                    }
                }
                ds[[i]] <- ds2
            }
            
            
            
            sim_summ <- list()
            dataset <- list()
            for (i in 1:length(ds)){
                # capture errors due to not all response option occurring in simulated dataset
                pars1 <- tryCatch({
                    mirt(ds[[i]][,3:(n_itm+2)],model=model,pars='values')
                }, warning = function(w){
                    return("warning")
                }, error = function(e){
                    return("error")
                })
                if(!pars1 %in% c("error","warning")) {
                    mcal <- mirt(ds[[i]][,3:(n_itm+2)],model=model,pars=pars1, verbose = FALSE)
                    # mcal <- mirt(ds[[i]][,3:(n_itm+2)],model,itemtype='2PL',pars=pars)
                    #   coef(mcal)
                    # mod_res_obj@Data$model
                    
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
                    
                    t6$ability_est <- ifelse(t6$ability_est %in% c(Inf,-Inf),NA,t6$ability_est)
                    
                    t6$resid <- t6$ability_est - t6$ability
                    
                    t6$abil_est_st <- (t6$ability_est - mean(t6$ability_est,na.rm=TRUE))*
                        (sd(t6$ability)/sd(t6$ability_est,na.rm=TRUE)) +
                        mean(t6$ability) # linear equating to true ability metric
                    
                    t6$resid_st <- t6$abil_est_st - t6$ability
                    t6$sample <- j
                    t6$rep_theta <- i
                    
                    # plot(t6$resid ~ t6$ability)
                    
                    sim_summ[[i]] <- t6[,c("ability","ability_est","resid","abil_est_st","resid_st","group")]
                    
                    if (save_sims == TRUE) {
                        if (j==1 & i == 1) {
                            sim_data <- t6
                        } else {
                            sim_data <- rbind(sim_data,t6)
                        }
                    }
                    
                } else {
                    j <- j-1
                }
            }
            
            ### This code block can be modified to output true ability, estimated
            #   ability, and simulated datasets
            # abil <- data.frame(matrix(ncol = length(sim_summ), nrow = 500))
            # abil_est <- data.frame(matrix(ncol = length(sim_summ), nrow = 500))
            # res <- data.frame(matrix(ncol = length(sim_summ), nrow = 500))
            stat <- data.frame(matrix(nrow = 12,ncol = length(sim_summ)))
            if (length(sim_summ) > 0) {
                for (i in 1:length(sim_summ)) {
                    nm <- paste("dset_",i,sep="")
                    # abil[,i] <- sim_summ[[i]]$ability
                    # names(abil)[i] <- nm
                    # abil_est[,i] <- sim_summ[[i]]$abil_est_st
                    # names(abil_est)[i] <- nm
                    # res[,i] <- sim_summ[[i]]$resid_st
                    # names(res)[i] <- nm
                    df <- sim_summ[[i]]
                    stat[1,i] <- sqrt(mean(df[df$group==1,"resid"]^2,na.rm=TRUE))
                    stat[2,i] <- sqrt(mean(df[df$group==2,"resid"]^2,na.rm=TRUE))
                    stat[3,i] <- sqrt(mean(df[df$group==1,"resid_st"]^2,na.rm=TRUE))
                    stat[4,i] <- sqrt(mean(df[df$group==2,"resid_st"]^2,na.rm=TRUE))
                    stat[5,i] <- mean(df[df$group==1,"ability_est"],na.rm=TRUE)
                    stat[6,i] <- mean(df[df$group==2,"ability_est"],na.rm=TRUE)
                    stat[7,i] <- mean(df[df$group==1,"abil_est_st"],na.rm=TRUE)
                    stat[8,i] <- mean(df[df$group==2,"abil_est_st"],na.rm=TRUE)
                    stat[9,i] <- sd(df[df$group==1,"ability_est"],na.rm=TRUE)
                    stat[10,i] <- sd(df[df$group==2,"ability_est"],na.rm=TRUE)
                    stat[11,i] <- sd(df[df$group==1,"abil_est_st"],na.rm=TRUE)
                    stat[12,i] <- sd(df[df$group==2,"abil_est_st"],na.rm=TRUE)
                    names(stat)[i] <- nm
                    stat[c(1,3,5,7,9,11),"group"] <- 1
                    stat[c(2,4,6,8,10,12),"group"] <- 2
                    stat[c(1,2,5,6,9,10),"type"] <- "raw"
                    stat[c(3,4,7,8,11,12),"type"] <- "standardized"
                    stat[c(1:4),"statistic"] <- "rmse"
                    stat[c(5:8),"statistic"] <- "est_mean"
                    stat[c(9:12),"statistic"] <- "est_sd"
                    stat$samp_num <- j
                } # end for i
            }
            if (j>0) {
                if (j==1) {
                    stat_summ <- stat
                } else {
                    stat_summ <- rbind(stat_summ,stat)
                }
            }
            #cat(paste("Iteration - ",j,", Elapsed time: ",Sys.time() - time,"\n",sep=""))
            incProgress(1/n_rep)
        } # end for j
        if(n_rep_theta == 1){
            stat_summ$avg <- stat_summ$dset_1
        } else {
            stat_summ$avg <- apply(stat_summ[,1:n_rep],1,mean)
        }
        statOut <- list(stat_summ = stat_summ, sim_data = sim_data)
        return(statOut)
    }
    
    # summaryStats is a function that calculates means across simulated datasets 
    # of summary statistics (rmse, rmse_st, mean, mean_st, sd sd_st) returned by 
    # equateSim.
    # 
    # Parameters:
    #   df - data.frame of summary statistics for each simulated dataset returned
    #     by equateSim
    #     
    # summaryStats returns a dataframe withthe means across simulated datasets of 
    # rmse, rmse_st, mean, mean_st, and sd sd_st for each group as well as the 
    # simple average of both groups. It includes a group variable (Group 1,
    # Group 2, Combined) 
    
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
        # rmse12 <- mean(rmse_grp_1,rmse_grp_2)
        # rmse_st12 <- mean(rmse_st_grp_1,rmse_st_grp_2)
        # mean12 <- mean(mean_grp_1,mean_grp_2)
        # mean_st12 <- mean(mean_st_grp_1,mean_st_grp_2)
        # sd12 <- mean(sd_grp_1,sd_grp_2)
        # sd_st12 <- mean(sd_st_grp_1,sd_st_grp_2)
        
        nms <- c("rmse","rmse_st","mean","mean_st","sd","sd_st")
        stat_1 <- data.frame(cbind(rmse_grp_1,rmse_st_grp_1,mean_grp_1,mean_st_grp_1,
                                   sd_grp_1,sd_st_grp_1))
        names(stat_1) <- nms
        stat_1$group <- 1
        stat_2 <- data.frame(cbind(rmse_grp_2,rmse_st_grp_2,mean_grp_2,mean_st_grp_2,
                                   sd_grp_2,sd_st_grp_2))
        names(stat_2) <- nms
        stat_2$group <- 2
        stat_12 <- data.frame(cbind(mean(c(rmse_grp_1,rmse_grp_2)),
                                    mean(c(rmse_st_grp_1,rmse_st_grp_2)),
                                    mean(c(mean_grp_1,mean_grp_2)),
                                    mean(c(mean_st_grp_1,mean_st_grp_2)),
                                    sqrt(0.5*(sd_grp_1^2 + sd_grp_2^2) + (0.5*(mean_grp_1 - mean_grp_2))^2),
                                    sqrt(0.5*(sd_st_grp_1^2 + sd_st_grp_2^2) + (0.5*(mean_st_grp_1 - mean_st_grp_2))^2)))

        names(stat_12) <- nms
        stat_12$group <- 0
        
        stat <- rbind(stat_1,stat_2,stat_12)
        stat$group <- factor(stat$group,levels=c(0,1,2),
                             labels=c("Combined","Group 1","Group 2"))
        return(stat)
    }
    
    alphas <- reactive({
        an <- eval(parse(text=paste0("c(", input$alphas, ")")))
        if(length(an) == 1) an <- rep(an, input$n_itm)
        return(an)
    })
    
    betas <- reactive({
        bn <- eval(parse(text=paste0("c(", input$betas, ")")))
        if(length(bn) == 1) bn <- rep(bn, input$n_itm)
        return(bn)
    })
    
    mcal_pars <- reactive({
        ipf <- input$parFile
        ipfpath <- ipf$datapath
        if(is.null(ipf)) {
            if(input$alphas == "") {
                if(input$betas == "") {
                    aa <- runif(input$n_itm, 0, 3)
                    dd <- rnorm(input$n_itm)
                    xmpl <- simdata(a = aa, d = dd, N = input$n_itm, itemtype = "dich")
                    xmod <- mirt.model(paste0("F = 1-", input$n_itm))
                    out <- mirt(xmpl, xmod, pars = "values", itemtype = "2PL")
                    out$value[out$name == "a1"] <- aa
                    out$value[out$name == "d"] <- dd
                } else {
                    aa <- runif(input$n_itm, 0, 3)
                    dd <- betas()
                    xmpl <- simdata(a = aa, d = dd, N = input$n_itm, itemtype = "dich")
                    xmod <- mirt.model(paste0("F = 1-", input$n_itm))
                    out <- mirt(xmpl, xmod, pars = "values", itemtype = "2PL")
                    out$value[out$name == "a1"] <- aa
                    out$value[out$name == "d"] <- dd
                }
            } else {
                if(input$betas == "") {
                    aa <- alphas()
                    dd <- rnorm(input$n_itm)
                    xmpl <- simdata(a = aa, d = dd, N = input$n_itm, itemtype = "dich")
                    xmod <- mirt.model(paste0("F = 1-", input$n_itm))
                    out <- mirt(xmpl, xmod, pars = "values", itemtype = "2PL")
                    out$value[out$name == "a1"] <- aa
                    out$value[out$name == "d"] <- dd
                } else {
                    aa <- alphas()
                    dd <- betas()
                    xmpl <- simdata(a = aa, d = dd, N = input$n_itm, itemtype = "dich")
                    xmod <- mirt.model(paste0("F = 1-", input$n_itm))
                    out <- mirt(xmpl, xmod, pars = "values", itemtype = "2PL")
                    out$value[out$name == "a1"] <- aa
                    out$value[out$name == "d"] <- dd
                }
            }
        } else out <- readRDS(ipfpath)
        return(out)
    })
    
    items1 <- reactive({
        eval(parse(text=paste0("c(", input$item_list_1, ")")))
    })
    
    items2 <- reactive({
        eval(parse(text=paste0("c(", input$item_list_2, ")")))
    })
    
    simOut <- eventReactive(input$runSim, {
        
        withProgress(message = "Running Simulation", value = 0, {
            
            mro <- input$mmrFile
            
            if(is.null(mro)){
                equateSim(seed = input$seed,
                          grp_mean = c(input$g1m, input$g2m),
                          grp_sd = c(input$g1s, input$g2s),
                          n_samp = input$n_samp,
                          n_rep_theta = input$n_rep_theta,
                          n_itm = input$n_itm,
                          pars = mcal_pars(),
                          n_rep = input$n_rep,
                          itms1 = items1(),
                          itms2 = items2(),
                          fsc_method = input$fsc_method,
                          mod_res_obj = NULL)
            } else {
                mropath <- mro$datapath
                mroRds <- readRDS(mropath)
                equateSim(seed = input$seed,
                          grp_mean = c(input$g1m, input$g2m),
                          grp_sd = c(input$g1s, input$g2s),
                          n_samp = input$n_samp,
                          n_rep_theta = input$n_rep_theta,
                          n_itm = input$n_itm,
                          pars = mcal_pars(),
                          n_rep = input$n_rep,
                          itms1 = items1(),
                          itms2 = items2(),
                          fsc_method = input$fsc_method,
                          mod_res_obj = mroRds)
            }
            
        })
        
        
    })
    
    output$parSummary <- renderReactjson({
        if (input$runSim == 0)
            return()
        
        isolate({
            mro <- input$mmrFile
            if(is.null(mro)){
                ps <- list(seed = input$seed,
                           grp_mean = c(input$g1m, input$g2m),
                           grp_sd = c(input$g1s, input$g2s),
                           n_samp = input$n_samp,
                           n_rep_theta = input$n_rep_theta,
                           n_itm = input$n_itm,
                           pars = mcal_pars(),
                           n_rep = input$n_rep,
                           itms1 = items1(),
                           itms2 = items2(),
                           fsc_method = input$fsc_method,
                           mod_res_obj = NULL)
            } else {
                mropath <- mro$datapath
                mroRds <- readRDS(mropath)
                ps <- list(seed = input$seed,
                           grp_mean = c(input$g1m, input$g2m),
                           grp_sd = c(input$g1s, input$g2s),
                           n_samp = input$n_samp,
                           n_rep_theta = input$n_rep_theta,
                           n_itm = input$n_itm,
                           pars = mcal_pars(),
                           n_rep = input$n_rep,
                           itms1 = items1(),
                           itms2 = items2(),
                           fsc_method = input$fsc_method,
                           mod_res_obj = mroRds)
            }
            reactjson(ps)
        })
    })
    
    output$equateSimResults <- renderTable({
        summaryStats(simOut()$stat_summ)
    })
    
    output$exportButton1CSV <- downloadHandler(
        filename = "equateSim_export_Summary.csv",
        content = function(file) {
            write.csv(simOut()$stat_summ, file)
        }
    )
    output$exportButton1Rds <- downloadHandler(
        filename = "equateSim_export_Summary.Rds",
        content = function(file) {
            saveRDS(simOut()$stat_summ, file)
        }
    )
    
    output$exportButton2CSV <- downloadHandler(
        filename = "equateSim_export_FullData.csv",
        content = function(file) {
            write.csv(simOut()$sim_data, file)
        }
    )
    output$exportButton2Rds <- downloadHandler(
        filename = "equateSim_export_FullData.Rds",
        content = function(file) {
            saveRDS(simOut()$sim_data, file)
        }
    )
}

# Run the application 
shinyApp(ui = ui, server = server)

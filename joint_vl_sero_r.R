rm(list=ls())

library(zoo)
library(data.table)
library(ggplot2)
library(gridExtra)
library(xtable)
library(rjags)
library(MASS)

############################################## 
# data processing and sampler set up
############################################## 

vl_alldat_dt <- readRDS("vl_alldat_dt.rds")

vl_alldat_dt <- vl_alldat_dt[order(pid_new, t)]
vl_alldat_dt_1 <- vl_alldat_dt[, .SD[1], by = c("pid_new")]
vl_alldat_dt_2 <- unique(vl_alldat_dt[, c("pid_new", "L_censored_peak_diag", "R_censored_peak_diag", "day_peak_diag", "day_last_ABNotPos", "day_first_ABPos")])

####

L_cen_index_vt <- which(vl_alldat_dt_2[, L_censored_peak_diag] == 1)
R_cen_index_vt <- which(vl_alldat_dt_2[, R_censored_peak_diag] == 1)
N_cen_index_vt <- which(vl_alldat_dt_2[, L_censored_peak_diag] == 0 & vl_alldat_dt_2[, R_censored_peak_diag] == 0)

####

sero_R_cen_index_vt <- which(is.na(vl_alldat_dt_2[, day_first_ABPos]))
sero_L_0_index_vt <- which(is.na(vl_alldat_dt_2[, day_last_ABNotPos]))
sero_interv_index_vt <- which(!is.na(vl_alldat_dt_2[, day_first_ABPos]) & !is.na(vl_alldat_dt_2[, day_last_ABNotPos]))

vl_alldat_dt_2[, "sero_R_cen" := ifelse(is.na(day_first_ABPos), 1, 0)]
vl_alldat_dt_2[, "sero_L_0" := ifelse(is.na(day_last_ABNotPos), 1, 0)]
vl_alldat_dt_2[, "sero_interv" := ifelse(!is.na(day_first_ABPos) & !is.na(day_last_ABNotPos), 1, 0)]

####

# number of total number of observations for each person
J_vt <- vl_alldat_dt[, .(J = .N), by = "pid_new"]$J

# index of the first observation for each person
i0_vt <- c(1, cumsum(head(J_vt, -1))+1)

n <- length(J_vt)
I_mt <- diag(3)

# sampler set up
nChain <- 2
nAdapt <- 10000
nThin <- 10
nBurn <- nThin*2000
nSave <- nThin*2000

# function to get initial values
getInitsList <- function(chain){
  return(list(".RNG.name" = "base::Super-Duper",  
              ".RNG.seed" = 12345 + 100*chain, 
              "tau_yy" = 1, "delta_Q" = 0.25, 
              "alpha0" = -4, "alpha1" = 8,
              "mu_log" = c(log(5), log(2), log(9)), 
              "Omega_log" = I_mt,
              "beta_vp" = rep(0, dataList$nBaseCov_vp),
              "t_p" = rep(0, dataList$n),
              "log_mt" = cbind(log(vl_alldat_dt_1[, diag_Log10VL_peak]), 
                               log(vl_alldat_dt_1[, day_peak_diag] - vl_alldat_dt_1[, day_firstpos_diag]+0.5),
                               log(vl_alldat_dt_1[, day_lastpos_diag] - vl_alldat_dt_1[, day_peak_diag]+0.5)),
              "tau_yy_sg" = 1, 
              "alpha0_sg" = -4, "alpha1_sg" = 8,
              "t_d_sg" = rep(0, dataList$n),
              "q_sg" = rep(0.75, dataList$n),
              "log_w_d_sg" = log(vl_alldat_dt_1[, day_lastpos_diag] - vl_alldat_dt_1[, day_lastpos_sgRNA] + 0.1),
              "beta_C" = rep(0, 2), 
              "kappa_1" = 12, "kappa_2" = 1/1.2,
              "C" = ifelse((vl_alldat_dt_2[, day_last_ABNotPos] == 28) & is.na(vl_alldat_dt_2[, day_first_ABPos]), 0, 1)
  ))
}

############################################## 
# model fitting
############################################## 

col_name_vt <- c("Covid", "Fever", "Nose", "Throat", "Chest", "GI", "Body", "TasteORSmell")

load.module("glm")

for (k in col_name_vt){
  
  BaseCov_vp <- as.matrix(vl_alldat_dt_1[, c("Age_35", "MaleSex", "HCQ", paste0(k,"AnyGrp")), with = F])
  
  dataList <- list(
    LoD_40_diag = LoD_40_diag,
    LoD_0_diag = LoD_0_diag,
    LoQ_31_diag = LoQ_31_diag,
    LoD_40_sgRNA = LoD_40_sgRNA,
    LoD_0_sgRNA = LoD_0_sgRNA,
    L_cen_index_vt = L_cen_index_vt,
    R_cen_index_vt = R_cen_index_vt,
    N_cen_index_vt = N_cen_index_vt,
    n = length(J_vt),
    J = J_vt, 
    i0 = i0_vt, 
    t = vl_alldat_dt[, t], 
    y = vl_alldat_dt[, y], 
    B = vl_alldat_dt[, B],
    S = vl_alldat_dt[, S],
    Q = vl_alldat_dt[, Q],
    I_mt = I_mt,
    BaseCov_vp = BaseCov_vp,
    nBaseCov_vp = ncol(BaseCov_vp),
    y_sg = vl_alldat_dt[, y_sg], 
    B_sg = vl_alldat_dt[, B_sg],
    S_sg = vl_alldat_dt[, S_sg],
    sero_R_cen_index_vt = sero_R_cen_index_vt,
    sero_L_0_index_vt = sero_L_0_index_vt,
    sero_interv_index_vt = sero_interv_index_vt,
    sero_R_cen_vt = vl_alldat_dt_2[, sero_R_cen],
    sero_L_0_vt = vl_alldat_dt_2[, sero_L_0],
    sero_interv_vt = vl_alldat_dt_2[, sero_interv],
    day_last_ABNotPos_vt = vl_alldat_dt_2[, day_last_ABNotPos],
    day_first_ABPos_vt = vl_alldat_dt_2[, day_first_ABPos],
    day_peak_diag_vt = vl_alldat_dt_2[, day_peak_diag]
  )
  
  #####
  
  modelStart <- Sys.time()
  jags_model <- jags.model("joint_vl_ab_jags.txt",
                           data = dataList,
                           inits = getInitsList,
                           n.chains = nChain, 
                           n.adapt = nAdapt)
  modelEnd <- Sys.time()
  print("Model:")
  print(modelEnd - modelStart)
  
  #####
  
  burnStart <- Sys.time()
  update(jags_model, n.iter = nBurn)
  burnEnd <- Sys.time()
  print("Burn:")
  print(burnEnd - burnStart)
  
  #####
  
  varnames_vt <- c("mu_log", "Sigma_log",
                   "mu_y", "sigma_yy", "delta_Q", "sigma_y",
                   "t_p", "v_p", "w_a", "w_b", "beta_up", "beta_down",
                   "S", "alpha0", "alpha1",
                   "beta_vp",
                   "mu_td", 
                   "gamma_1", "gamma_2", "mu_q", 
                   "mu_lwd", "sigma_lwd", 
                   "mu_y_sg", "sigma_yy_sg", "sigma_y_sg",
                   "t_p_sg", "v_p_sg", "w_a_sg", "w_b_sg", "beta_up_sg", "beta_down_sg",
                   "t_d_sg", "q_sg", "w_d_sg",
                   "S_sg", "alpha0_sg", "alpha1_sg",
                   "C", "beta_C", "p_cens", "kappa_1", "kappa_2", "mu_w_s")
  
  sampStart <- Sys.time()
  coda_samp <- coda.samples(jags_model, 
                            variable.names = varnames_vt, 
                            n.iter = nSave, thin = nThin)
  sampEnd <- Sys.time()
  print("Sample:")
  print(sampEnd - sampStart)
  
  #####
  
  save(col_name_vt, vl_alldat_dt, 
       nChain, nAdapt, nThin, nBurn, nSave, 
       dataList, getInitsList, varnames_vt,
       modelStart, modelEnd, burnStart, burnEnd, sampStart, sampEnd, 
       jags_model, coda_samp, 
       file = paste0("model_results_", k, ".RData"))
  
  print(k)
  
}

##############################
# model diagnostic reports
##############################

load("model_results_Covid.RData")

coda_samp_lt <- list()
for (chain in 1:nChain){
  coda_samp_lt[[chain]] <-  coda_samp[[chain]]
}

coda_samp_lt <- mcmc.list(coda_samp_lt)

modelPar_vt <- c("alpha0", "alpha1",
                 "mu_log[1]", "mu_log[2]", "mu_log[3]", 
                 "Sigma_log[1,1]", "Sigma_log[2,2]", "Sigma_log[3,3]",
                 "Sigma_log[1,2]", "Sigma_log[1,3]", "Sigma_log[2,3]",
                 "sigma_yy", "delta_Q",
                 paste0("beta_vp[", 1:dataList$nBaseCov_vp, "]"),
                 "alpha0_sg", "alpha1_sg",
                 "mu_td", 
                 "gamma_1", "gamma_2", "mu_q",
                 "mu_lwd", "sigma_lwd", 
                 "sigma_yy_sg",
                 paste0("beta_C[", 1:2, "]"),
                 "kappa_1", "kappa_2", "mu_w_s")


ggmcmc(ggs(coda_samp_lt[, modelPar_vt]), file = "model_diag_Covid.pdf")



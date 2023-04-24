rm(list=ls())

library(zoo)
library(data.table)
library(ggplot2)
library(gridExtra)
library(xtable)
library(rjags)
library(MASS)

logit<-function(x){log(x/(1-x))}
expit<-function(x){exp(x)/(1+exp(x))}

############################################## 

dummy_dataset <- readRDS("dummy_dataset.rds")

vl_alldat_dt <- dummy_dataset[order(pid_new, day)]

##############################

LoD_40_diag <- log10((10^(13.59 - 0.279 * 40) + 10^(13.35 - 0.276 * 40))/2)
LoQ_31_diag <- log10((10^(13.59 - 0.279 * 31) + 10^(13.35 - 0.276 * 31))/2)
LoD_0_diag <- log10((10^(13.59 - 0.279 * 0) + 10^(13.35 - 0.276 * 0))/2)

vl_alldat_dt[, "t" := day_since_peak_diag]
vl_alldat_dt[, "y" := diag_Log10VL]
vl_alldat_dt[, "B" := diag_pos]
vl_alldat_dt[, "Q" := ifelse(diag_Log10VL > LoD_40_diag & diag_Log10VL < LoQ_31_diag, 1, 0)]

vl_alldat_dt[, "L_censored_peak_diag" := ifelse(day_firstpos_diag - day_first <= 1 & day_peak_diag - day_firstpos_diag <= 1, 1, 0)]
vl_alldat_dt[, "R_censored_peak_diag" := ifelse(day_last - day_lastpos_diag <= 1 & day_lastpos_diag - day_peak_diag <= 1, 1, 0)]

##############################

LoD_40_sgRNA <- 13.79 - 0.286 * 40
LoD_0_sgRNA <- 13.79 - 0.286 * 0

vl_alldat_dt[, "y_sg" := ifelse(sgRNA_pos == 0, LoD_40_sgRNA, sgRNA_Log10VL)]
vl_alldat_dt[, "B_sg" := sgRNA_pos]

##############################

vl_alldat_dt <- vl_alldat_dt[order(pid_new, t)]

vl_alldat_dt_1 <- vl_alldat_dt[, .SD[1], by = c("pid_new")]
vl_alldat_dt_2 <- unique(vl_alldat_dt[, c("pid_new", "L_censored_peak_diag", "R_censored_peak_diag", "day_peak_diag", "day_last_ABNotPos", "day_first_ABPos")])

####

L_cen_index_vt <- which(vl_alldat_dt_2[, L_censored_peak_diag] == 1)
R_cen_index_vt <- which(vl_alldat_dt_2[, R_censored_peak_diag] == 1)
N_cen_index_vt <- which(vl_alldat_dt_2[, L_censored_peak_diag] == 0 & vl_alldat_dt_2[, R_censored_peak_diag] == 0)

####

sero_R_cen_index_vt <- which(is.na(vl_alldat_dt_2[, day_first_ABPos]) & !is.na(vl_alldat_dt_2[, day_last_ABNotPos]))
sero_L_0_index_vt <- which(!is.na(vl_alldat_dt_2[, day_first_ABPos]) & is.na(vl_alldat_dt_2[, day_last_ABNotPos]))
sero_interv_index_vt <- which(!is.na(vl_alldat_dt_2[, day_first_ABPos]) & !is.na(vl_alldat_dt_2[, day_last_ABNotPos]))

vl_alldat_dt_2[, "sero_R_cen" := ifelse(is.na(day_first_ABPos) & !is.na(day_last_ABNotPos), 1, 0)]
vl_alldat_dt_2[, "sero_L_0" := ifelse(!is.na(day_first_ABPos) & is.na(day_last_ABNotPos), 1, 0)]
vl_alldat_dt_2[, "sero_interv" := ifelse(!is.na(day_first_ABPos) & !is.na(day_last_ABNotPos), 1, 0)]

# number of total number of observations
J_vt <- vl_alldat_dt[, .(J = .N), by = "pid_new"]$J

# location of the first observation for each pp
i0_vt <- c(1, cumsum(head(J_vt, -1))+1)

# vectors and matrices
n <- length(J_vt)
I_mt <- diag(3)

# sample
nChain <- 2
nAdapt <- 10000
nThin <- 10
nBurn <- nThin*2000
nSave <- nThin*2000

# inits list
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


##############################

load.module("glm")

BaseCov_vp <- as.matrix(vl_alldat_dt_1[, c("Age_35", "MaleSex", "HCQ", "CovidAnyGrp"), with = F])

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
  Q = vl_alldat_dt[, Q],
  I_mt = I_mt,
  BaseCov_vp = BaseCov_vp,
  nBaseCov_vp = ncol(BaseCov_vp),
  y_sg = vl_alldat_dt[, y_sg], 
  B_sg = vl_alldat_dt[, B_sg],
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
jags_model <- jags.model("joint_vl_sero_jags.txt",
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



## Antibody Modelling Function for mice
## Create Fig 4-A





## Read parameters
pop <- read.csv("populationParameters.txt",  row.names = 1)

## Set time range
Tmin <- 0
Tmax <- 180 
step_size <- 0.005
times<-c(seq(Tmin,Tmax,step_size))

## ODE Function to model population antibody dynamics
ode_ab_pop <- function(pars) {
  
  k1 <- as.numeric(pars[1]) # Growth rate in phase 1
  k2 <- as.numeric(pars[2]) # Decay rate in phase 2
  A0 <- as.numeric(pars[3]) # Initial condition
  C <- as.numeric(pars[4])  # Carrying capacity for sigmoid growth phase
  
  times <- c(seq(Tmin, Tmax, step_size))
  
  k_time <- data.frame(times = times, k = rep(0, length(times)))
  k_time <- k_time %>% mutate(k = case_when(times <= 35 ~ k1,  # set 35 days as the timing of peak
                                            times > 35 ~ -k2))
  
  k_t <- approxfun(k_time$times, k_time$k, rule = 2)
  
  derivs <- function(times, y, pars, k_t) {
    with(as.list(c(pars, y)), {
      k <- k_t(times)   # Get k value based on time
      if (times < 35) {
        dA <- k * A * (1 - A / C) # Sigmoid growth phase
      } else {
        dA <- k * A # Exponential decay phase
      }
      return(list(c(dA)))
    })
  }
  
  y <- c(A = A0) # Initial observation
  
  out <- ode(y = y, parms = pars, times = times, func = derivs, k_t = k_t)
  as.data.frame(out)
}

## Sampling function to simulate the CI for population dynamics
sample_ab <- function(pop, num, vac){
  
  inds_mean <- which(rownames(pop) %in% c("k1_pop","k2_pop","A0_pop", 'peak_pop')) #row number of parameters
  inds_sd <- which(rownames(pop) %in% c("omega_k1","omega_k2","omega_A0", "omega_peak"))

  pars <- matrix(0, num+1, length(inds_mean))
  
  for (i in 1:length(inds_mean)) {
    mean_par <- pop$value[inds_mean[i]]
    sd_par <- pop$value[inds_sd[i]]
    #sd_par <- pop$se_sa[inds_mean[i]] ####standard error of population parameter
    
    if (i==1) {
      beta_k1_0_1ug <- pop$value[inds_mean[i]+1]
      beta_k1_0_2ug <- pop$value[inds_mean[i]+2]
      beta_k1_0_5ug <- pop$value[inds_mean[i]+3]
      beta_k1_1ug   <- pop$value[inds_mean[i]+4]
      log_k1 = log(mean_par)+(vac=='0.1ug Pfizer')*beta_k1_0_1ug+(vac=='0.2ug Pfizer')*beta_k1_0_2ug+
        (vac=='0.5ug Pfizer')*beta_k1_0_5ug + (vac=='1ug Pfizer')*beta_k1_1ug
      pars[1:num, i] <- exp(rnorm(num, mean = log_k1, sd=sd_par))
      pars[num+1, i] <- exp(log_k1)
    } 
    
    else if (i==2 ) {
      beta_k2_0_1ug <- pop$value[inds_mean[i]+1]
      beta_k2_0_2ug <- pop$value[inds_mean[i]+2]
      beta_k2_0_5ug <- pop$value[inds_mean[i]+3]
      beta_k2_1ug   <- pop$value[inds_mean[i]+4]
      log_k2 = log(mean_par)+(vac=='0.1ug Pfizer')*beta_k2_0_1ug+(vac=='0.2ug Pfizer')*beta_k2_0_2ug+
        (vac=='0.5ug Pfizer')*beta_k2_0_5ug +(vac=='1ug Pfizer')*beta_k2_1ug
      pars[1:num, i] <- exp(rnorm(num, mean = log_k2, sd=sd_par))
      pars[num+1, i] <- exp(log_k2)
    } 
    
    else if (i==3) {
      log_A0 = log(mean_par)
      pars[1:num, i] <- exp(rnorm(num, mean=log_A0, sd=sd_par))
      pars[num+1, i] <- exp(log_A0)
    }
    else if (i==4 ) {
      log_C = log(mean_par)
      pars[1:num, i] <- exp(rnorm(num, mean=log_C, sd=sd_par))
      pars[num+1, i] <- exp(log_C)
    } 
  }
  return(pars)
}

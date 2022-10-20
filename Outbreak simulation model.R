###########################################################################################
###Individual-based negative binomial branching process model adpated from            #####
###https://www.thelancet.com/journals/langlo/article/PIIS2214-109X(20)30074-7/fulltext#####
###########################################################################################
library(purrr)
library(data.table)
library(dplyr)
#paramters:
#num.initial.cases: number of initial cases of each outbreak
#IPfn: assumed incubation period distribution estimated from https://www.sciencedirect.com/science/article/pii/S1755436521000359
#GIfn: assumed generation interval distribution estimated from https://www.sciencedirect.com/science/article/pii/S1755436521000359
#ISOfn: assumed isolation delay distribtuion estimated from the observed data
#prop.asym: asymptomatic ratio estimated from the observed data
#r0community: time-varying reproduction number (from our estimates or from EpiEstim packge) for non-isolated cases 
#rr0isolated: reproduction number for isolated cases which is 0
#disp.iso: dispersion paramter for isolated cases
#disp.community: time-varying or fixed dispersion parameter for non-isolated cases 
#prop.ascertain: effectiveness of contact tracing and is assumed to be 0
#num.add: number of additional source cases added in each generation
#quarantine: irrelevant paramters from the original models
#CC_speed: speed of contact tracing (days), which had no effect on the model output in our siumulation scenarios
dist_setup <- function(dist_shape = NULL, dist_rate = NULL) {
  out <- purrr::partial(rgamma,
                        shape = dist_shape,
                        rate = dist_rate)
  return(out)
}
#--functions to genrate sample of certain distributions distributions--#
GIfn <- dist_setup(dist_shape = 6.7^2/1.8^2,dist_rate =6.7/1.8^2)
IPfn <- dist_setup(dist_shape = 6.8^2/4.1^2,dist_rate =6.8/4.1^2)
ISOfn <- function(n,mean=3.9,sd=3.7){rnorm(n,mean=mean,sd=sd)}

#--initialize the outbreaks--#
outbreak_setup <- function(num.initial.cases, IPfn, GIfn, ISOfn, prop.asym) {
  # sample incubation period for new cases
  IP_samples <- IPfn(num.initial.cases)
  exposure_time <- GIfn(num.initial.cases)
  case_data <- data.table(exposure_time = exposure_time, # Exposure time of 0 for all initial cases
                          asym = purrr::rbernoulli(num.initial.cases, prop.asym),
                          caseid = 1:(num.initial.cases), # set case id
                          infector = 0,
                          missed = TRUE,
                          onset_time = exposure_time+IP_samples,
                          new_cases = NA,
                          generation = 0)
  
  # set isolation time for cluster to minimum time of onset of symptoms + draw from delay distribution
  case_data <- case_data[, iso_time := onset_time + ISOfn(1)][, isolated := FALSE]
  
  case_data$iso_time[case_data$asym] <- Inf
  
  # return
  return(case_data)
}

#--generating the next generation of cases--#
###
outbreak_step <- function(case_data = NULL, disp.iso = NULL, 
                          disp.com = NULL, r0isolated = NULL, 
                          r0community = NULL, prop.asym = NULL, 
                          IPfn = NULL, ISOfn = NULL, 
                          prop.ascertain = NULL, quarantine = NULL, 
                          num.add=NULL,CC_speed=NULL) {
  
  # A vectorised version of isTRUE
  vect_isTRUE <- function(x) {
    purrr::map_lgl(x, isTRUE)
  }
  
  vect_max <- function(x, y) {
    purrr::map2_dbl(x, y, max)
  }
  
  vect_min <- function(x, y) {
    purrr::map2_dbl(x, y, min)
  }
  # For each case in case_data, draw new_cases from a negative binomial distribution
  # with an R0 and dispersion dependent on if isolated=TRUE
  case_data[, new_cases := purrr::map2_dbl(
    ifelse(vect_isTRUE(isolated), 1, disp.com),
    ifelse(vect_isTRUE(isolated),
           0,
           r0community),
    ~ rnbinom(1, size = .x, mu = .y))
  ]
  #
  #add additional source
  case_data.add <- data.table(exposure_time = case_data[generation==max(generation),sample(quantile(exposure_time,c(0.7,0.8,0.9,1)),num.add,replace = T)], # Exposure time of 0 for all initial cases
                              asym = purrr::rbernoulli(num.add, prop.asym),
                              caseid = (nrow(case_data)+1):(num.add+nrow(case_data)), # set case id
                              infector = 999999,
                              missed = TRUE,
                              onset_time = NA,
                              new_cases = NA,
                              generation = case_data[,max(generation)])
  case_data.add[,onset_time:=exposure_time+IPfn(.N)]
  #
  case_data.add <- case_data.add[, iso_time := onset_time + ISOfn(1)][, isolated := FALSE]
  case_data.add$iso_time[case_data.add$asym] <- Inf
  #
  case_data.add[, new_cases := purrr::map2_dbl(
    rep(disp.com,nrow(case_data.add)),
    rep(r0community,nrow(case_data.add)),
    ~ rnbinom(1, size = .x, mu = .y))
  ]
  
  #bind data
  case_data=rbind(case_data,case_data.add)
  # Select cases that have generated any new cases
  new_case_data <- case_data[new_cases > 0]
  # The total new cases generated
  total_new_cases <- case_data[, sum(new_cases,na.rm = T), ]
  
  # If no new cases drawn, outbreak is over so return case_data
  if (total_new_cases == 0) {
    # If everyone is isolated it means that either control has worked or everyone has had a chance to infect but didn't
    case_data$isolated <- TRUE
    
    case_data[caseid %in% case_data.add$caseid , generation:=max(case_data$generation)+1]
    effective_r0 <- 0
    cases_in_gen <- 0
    out <- list(case_data, effective_r0, cases_in_gen)
    names(out) <- c("cases", "effective_r0", "cases_in_gen")
    
    return(out)
  }
  # Compile a data.table for all new cases, new_cases is the amount of people that each infector has infected
  prob_samples <- data.table(
    # time when new cases were exposed, a draw from distribution of generation interval
    exposure_time = unlist(purrr::map2(new_case_data$new_cases, new_case_data$exposure,
                                       function(x, y) {
                                         y+GIfn(x)
                                       })),
    # records the infector of each new person
    infector = unlist(purrr::map2(new_case_data$caseid, new_case_data$new_cases,
                                  function(x, y) {
                                    rep(as.integer(x), as.integer(y))
                                  })),
    # records when infector was isolated
    infector_iso_time = unlist(purrr::map2(new_case_data$iso_time, new_case_data$new_cases,
                                           function(x, y) {
                                             rep(x, as.integer(y))
                                           })),
    # records if infector asymptomatic
    infector_asym = unlist(purrr::map2(new_case_data$asym, new_case_data$new_cases,
                                       function(x, y) {
                                         rep(x, y)
                                       })),
    # draws a sample to see if this person is asymptomatic
    asym = purrr::rbernoulli(n = total_new_cases, p = prop.asym),
    # draws a sample to see if this person is traced
    missed = purrr::rbernoulli(n = total_new_cases, p = 1 - prop.ascertain),
    # sample from the incubation period for each new person
    new_IP_sample = IPfn(total_new_cases),
    isolated = FALSE,
    new_cases = NA,
    generation = case_data[,max(generation)]+1
  )
  
  
  prob_samples <- prob_samples[exposure_time < infector_iso_time][, # filter out new cases prevented by isolation
                                                                  `:=`(# onset of new case is exposure + incubation period sample
                                                                    onset_time = exposure_time + new_IP_sample)]
  
  # If you are not missed, you will be immediately isolated depend on the speed of contact tracing
  prob_samples[, iso_time := ifelse(!vect_isTRUE(missed), infector_iso_time+CC_speed,
                                    # If you are missed and asymptomatic, you will never be isolated
                                    ifelse(vect_isTRUE(asym), Inf,
                                           #if you are missed but symptomatic, you will be isolated with isolation delay
                                           onset_time + ISOfn(1)))]
  
  # Chop out unneeded sample columns
  prob_samples[, c("new_IP_sample", "infector_iso_time", "infector_asym") := NULL]
  # Set new case ids for new people
  prob_samples$caseid <- (nrow(case_data) + 1):(nrow(case_data) + nrow(prob_samples))
  
  ## Number of new cases
  cases_in_gen <- nrow(prob_samples)
  
  ## Estimate the effective r0
  effective_r0 <- nrow(prob_samples) / nrow(case_data[!vect_isTRUE(case_data$isolated)])
  
  # Everyone in case_data so far has had their chance to infect and are therefore considered isolated
  case_data$isolated <- TRUE
  
  # bind original cases + new secondary cases
  case_data <- data.table::rbindlist(list(case_data, prob_samples),
                                     use.names = TRUE)
  
  # Return
  out <- list(case_data, effective_r0, cases_in_gen)
  names(out) <- c("cases", "effective_r0", "cases_in_gen")
  
  return(out)
}

#--incorporating dynamic R and k into the model--#
#R=c()
#k=c()
outbreak_model=function(num.initial.cases=NULL,disp.iso = NULL, disp.com = NULL, r0isolated = NULL, r0community = NULL,
                        prop.asym = NULL, IPfn = NULL, ISOfn = NULL, prop.ascertain = NULL, quarantine = NULL, num.add=NULL,CC_speed=NULL){
  case_data <- outbreak_setup(num.initial.cases = num.initial.cases,
                              IPfn = IPfn,
                              GIfn=GIfn,
                              prop.asym = prop.asym,
                              ISOfn = ISOfn)
  for(i in 1:length(k)){
    disp.com=k[i]
    r0community = R[i]
    out <- outbreak_step(case_data = case_data,
                         disp.iso = disp.iso,
                         disp.com = disp.com,
                         r0isolated = r0isolated,
                         r0community = r0community,
                         IPfn = IPfn,
                         ISOfn = ISOfn,
                         prop.ascertain = prop.ascertain,
                         prop.asym = prop.asym,
                         num.add = num.add,
                         CC_speed = CC_speed)
    case_data <- out[[1]]
  }
  case_data$total_cases=nrow(case_data)
  return(case_data)
}
#--run 50 simulations--#
res=purrr::map(.x = 1:50, ~ outbreak_model(num.initial.cases = 5,
                                            disp.iso = 1,
                                            disp.com = NULL,
                                            r0isolated = 0,
                                            r0community = NULL,
                                            IPfn = IPfn,
                                            ISOfn = ISOfn,
                                            prop.ascertain = 0,
                                            prop.asym = 0.33,
                                            num.add = 5,
                                            CC_speed = sample(0:1,1)))
res.aggregate=rbindlist(res,idcol=TRUE)











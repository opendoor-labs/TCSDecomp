#' State space model
#' Creates a state space model in list form
#' yt = Ht * Bt + e_t
#' Bt = Ft * B_{t=1} + u_t
#'
#' @param par Vector of named parameter values
#' @param yt Univariate time series of data values
#' @param freq Seasonality of the data (1, 4, 12, 53, 365)
#' @param decomp Decomposition model ("tend-cycle-seasonal", "trend-seasonal", "trend-cycle", "trend-noise")
#' @param trend_spec Trend specification (NULL, "rw", "rwd", "2rw"). The default is NULL which will choose the best of all specifications based on the maximum likielhood. 
#' "rw" is the random walk trend. "rwd" is the random walk with random walk drift trend. "2rw" is a 2nd order random walk trend.
#' @param init Initial state values for the Kalman filter
#' @return List of space space matrices
#' @examples
#' SSmodel(par, y, freq, decomp, trend_spec, init)
#' @author Alex Hubbard (hubbard.alex@gmail.com)
#' @export
SSmodel = function(par, yt, freq, decomp, trend_spec, init = NULL){
  yt = matrix(yt, nrow = 1)
  
  #Define the standard deviation of the observation equation
  if(decomp == "trend-noise") {
    sig_e = ifelse(trend_spec == "rw", sd(diff(t(yt))), 
                   ifelse(trend_spec %in% c("rwd", "2rw"), sd(diff(diff(t(yt))))))
  }else if(decomp %in% c("trend-cycle", "trend-seasonal", "trend-cycle-seasonal", "trend-seasonal-cycle")) {
    sig_e = par["sig_e"]
  }else {
    stop("decomp must be one of 'trend-noise', 'trend-cycle', 'trend-seasonal', 'trend-cycle-seasonal', or 'trend-seasonal-cycle'.")
  }
  
  #Define the transition and observation equation matrices based on the trend specification
  if(trend_spec == "rw"){
    #T_t = T_{t-1} + e_t, e_t ~ N(0, sig_t^2)
    #Transition matrix
    Fm = matrix(c(1), nrow = 1, ncol = 1)
    colnames(Fm) = "Tt1"
    rownames(Fm) = "Tt0"
    #Observation matrix
    Hm = matrix(c(1), nrow = 1)
    colnames(Hm) = rownames(Fm)
  }else if(trend_spec == "rwd"){
    #T_t = M_{t-1} + T_{t-1} + e_t, e_t ~ N(0, sig_t^2)
    #M_t = M_{t-1} + n_t, n_t ~ N(0, sig_m^2)
    #Transition matrix
    Fm = rbind(c(1, 1), c(0, 1))
    colnames(Fm) = c("Tt1", "Mt1")
    rownames(Fm) = c("Tt0", "Mt0")
    #Observation matrix
    Hm = matrix(c(1, 0), nrow = 1)
    colnames(Hm) = rownames(Fm)
  }else if(trend_spec == "2rw"){
    #T_t = 2T_{t-1} - T_{t-2} + e_t, e_t ~ N(0, sig_t^2)
    #Transition matrix
    Fm = rbind(c(2, -1), c(1, 0))
    colnames(Fm) = c("Tt1", "Tt2")
    rownames(Fm) = c("Tt0", "Tt1")
    #Observation matrix
    Hm = matrix(c(1, 0), nrow = 1)
    colnames(Hm) = rownames(Fm)
  }
  
  #Define the transition error covariance matrix
  Qm = matrix(0, nrow = nrow(Fm), ncol = ncol(Fm))
  colnames(Qm) = rownames(Qm) = rownames(Fm)
  Qm[rownames(Qm) == "Tt0", colnames(Qm) == "Tt0"] = par["sig_t"]^2
  Qm[rownames(Qm) == "Mt0", colnames(Qm) == "Mt0"] = par["sig_m"]^2
  
  #Define the cycle component
  if(grepl("cycle", decomp)) {
    rho = 1/(1 + exp(-par["rho"])) #Constrain rho (the dampening factor) to be between 0 and 1
    lambda = pi/(1 + exp(-par["lambda"])) #Constrain lambda (the period) to between 0 and pi
    colnames = colnames(Fm)
    rownames = rownames(Fm)
    Cm = rbind(c(rho * cos(lambda), rho * sin(lambda)), c(-rho * sin(lambda), rho * cos(lambda)))
    Fm = as.matrix(Matrix::bdiag(Fm, Cm))
    colnames(Fm) = c(colnames, "Ct1", "Ct2")
    rownames(Fm) = c(rownames, "Ct0", "Ct1")
    Hm = cbind(Hm, matrix(c(1, 0), nrow = 1))
    colnames(Hm) = rownames(Fm)
    Qm = as.matrix(Matrix::bdiag(Qm, diag(2)*par["sig_c"]))
    colnames(Qm) = rownames(Qm) = rownames(Fm)
  }
  
  #Define the seasonal component
  if(grepl("seasonal", decomp)) {
    jiter = as.numeric(unique(gsub("[[:alpha:]]|[[:punct:]]", "", names(par)[grepl("sig_j", names(par))])))
    for(j in jiter) {
      colnames = colnames(Fm)
      rownames = rownames(Fm)
      Sm = rbind(c(cos(2 * pi * j/freq), sin(2 * pi * j/freq)), 
                 c(-sin(2 * pi * j/freq), cos(2 * pi * j/freq)))
      Fm = as.matrix(Matrix::bdiag(Fm, Sm))
      colnames(Fm) = c(colnames, paste0("Stl", j), paste0("Stls", j))
      rownames(Fm) = c(rownames, paste0("St", j), paste0("Sts", j))
    }
    Hm = cbind(Hm, matrix(rep(c(1, 0), length(jiter)), nrow = 1))
    colnames(Hm) = rownames(Fm)
    Qm = as.matrix(Matrix::bdiag(Qm, diag(c(rbind(par[grepl("sig_j", names(par))], par[grepl("sig_j", names(par))])))))
    colnames(Qm) = rownames(Qm) = rownames(Fm)
  }
  
  #Transition equation intercept matrix
  Dm = matrix(0, nrow = nrow(Fm))
  
  #Observaton equation intercept matrix
  Am = matrix(0, nrow = 1, ncol = 1)
  
  #Observation equation error covariance matrix
  Rm = matrix(sig_e^2, nrow = 1, ncol = 1)
  
  #Initial guess for unobserved vector
  if(is.null(init)) {
    B0 = c(rep(mean(c(range(yt[1:ifelse(freq == 1, round(length(y)/4), freq)], na.rm = T)), na.rm = T), length(rownames(Fm)[grepl("Tt", rownames(Fm))])), 
           rep(0, length(rownames(Fm)[grepl("Mt", rownames(Fm))])), 
           rep(0, length(rownames(Fm)[grepl("Ct", rownames(Fm))])), 
           rep(0, length(rownames(Fm)[grepl("St", rownames(Fm))])))
  }else{
    B0 = init[["B0"]]
  }
  names(B0) = c(rownames(Fm)[grepl("Tt", rownames(Fm))], 
                rownames(Fm)[grepl("Mt", rownames(Fm))], 
                rownames(Fm)[grepl("Ct", rownames(Fm))], 
                rownames(Fm)[grepl("St", rownames(Fm))])
  
  #Initial guess for variance of the unobserved vector
  if(is.null(init)) {
    P0 = diag(100, nrow = nrow(Fm), ncol = nrow(Fm))
    rownames(P0) = colnames(P0) = rownames(Fm)
    P0["Tt0", "Tt0"] = var(c(yt[1:ifelse(freq == 1, round(length(y)/4), freq)]), na.rm = T)
    if(trend_spec == "2rw") {
      P0["Tt1", "Tt1"] = P0["Tt0", "Tt0"]
    }
  }else{
    P0 = init[["P0"]]
  }
  return(list(B0 = B0, P0 = P0, At = Am, Dt = Dm, Ht = Hm, Ft = Fm, Rt = Rm, Qt = Qm))
}

#' Trend cycle seasonal decomposition using the Kalman filter
#
#' Estimates a structural time series model using the Kalman filter and maximum likelihood
#' The seasonal and cycle components are assumed to be of a trigonometric form.
#' The function checks 4 trend specifications to decompose a univariate time series
#' into trend, cycle, and/or seasonal components plus noise. The function automatically
#' detects the frequency and checks for a seasonal and cycle component if the user does not specify
#' the frequency or decomposition model. This can be turned off by setting freq or specifying decomp.
#'
#' State space model for decomposition
#' Yt = T_t + C_t + S_t + e_t, e_t ~ N(0, sig_e^2)
#' Y is the data
#' T is the trend component
#' C is the cycle component
#' S is the seasonal component
#' e is the observation error
#' @param y Univariate time series of data values. May also be a 2 column data frame containing a date column.
#' @param freq Seasonality of the data (1, 4, 12, 53, 365)
#' @param decomp Decomposition model ("tend-cycle-seasonal", "trend-seasonal", "trend-cycle", "trend-noise")
#' @param trend_spec Trend specification (NULL, "rw", "rwd", "2rw"). The default is NULL which will choose the best of all specifications based on the maximum likielhood. 
#' "rw" is the random walk trend. "rwd" is the random walk with random walk drift trend. "2rw" is a 2nd order random walk trend.
#' @param optim_methods Vector of 1 to 3 optimization methods in order of preference ("NR", "BFGS", "CG", "BHHH", or "SANN")
#' @param det_obs Set the observation equation error variance to 0 (deterministic observation equation)
#' @param det_trend Set the trend error variance to 0 (deterministic trend)
#' @param det_seasonality Set the seasonality error variances to 0 (deterministic seasonality)
#' @param det_cycle Set thecycle error variance to 0 (deterministic cycle)
#' @param det_drift Set the drift error variance to 0 (deterministic drift)
#' @param maxit Maximum number of iterations for the optimization
#' @param maxtrials Maximum number of optimization trials to get convergence
#' @return List of estimation values including coefficients, convergence code, datea frequency, decomposition used, and trend specification selected.
#' @examples
#' tcs_decomp_estim(y = DT[, c("date", "y")])
#' @author Alex Hubbard (hubbard.alex@gmail.com)
#' @export
tcs_decomp_estim = function (y, freq = NULL, decomp = NULL, trend_spec = NULL, det_obs = F, 
                             det_trend = F, det_seasonality = F, det_cycle = F, det_drift = F, 
                             level = 0.01, optim_methods = c("BFGS", "CG", "NM"), maxit = 1000, maxtrials = 10){
  dates = NULL
  if(is.null(freq)){
    y = data.table::as.data.table(y)
    .SD = data.table::.SD
    datecol = unlist(lapply(colnames(y), function(x){
      if(class(y[, c(x), with = F][[1]]) %in% c("Date")){
        return(x)
      }else{
        return(NULL)
      }
    }))
    if(length(datecol) == 1) {
      #Detect the frequency
      if(length(datecol) > 1){
        stop("Too many date columns. Include only 1 date column or set the frequency manually.")
      }
      datediffs = unique(diff(y[, datecol, with = F][[1]]))
      freq = datediffs[which.max(tabulate(match(diff(y[, datecol, with = F][[1]]), datediffs)))]
      freq = c(365, 52, 12, 4, 1)[which.min(abs(freq -  c(1, 7, 30, 90, 365)))]
      dates = y[, datecol, with = F][[1]]
      y = y[, colnames(y)[colnames(y) != datecol], with = F][[1]]
      rm(datediffs, datecol)
    }else if(length(datecol) > 1){
      stop("Too many date columns detected.")
    }else{
      stop("No date column detected. Include a date column or set the frequency.")
    }
  }
  else{
    if(!is.numeric(freq)){
      stop("Must provide freq as numeric (1 for annual, 4 for quarterly, 12 for monthly, 52 for weekly, 365 for daily).")
    }else if(!freq %in% c(1, 4, 12, 52, 365)){
      stop("Must provide freq as numeric (1 for annual, 4 for quarterly, 12 for monthly, 52 for weekly, 365 for daily).")
    }
  }
  if(level < 0.01 | level > 0.1){
    stop("level must be between 0.01 and 0.1.")
  }
  if(any(!optim_methods %in% c("NR", "BFGS", "BHHH", "SANN", "CG", "NM")) | length(optim_methods) < 1 | length(optim_methods) > 3){
    stop("optim_methods must be a vector of length 1 to 3 containing 'NR', 'BFGS', 'BHHH', 'SANN', 'CG', or 'NM'")
  }
  if(!is.numeric(maxit)){
    stop("maxit must be numeric and greater than 0.")
  }else if(maxit <= 0){
    stop("maxit must be numeric and greater than 0.")
  }
  if(!is.numeric(maxtrials)){
    stop("maxtrials must be numeric and greater than 0.")
  }else if(maxtrials <= 0){
    stop("maxtrials must be numeric and greater than 0.")
  }

  #Set the decomposition
  if(is.null(decomp)){
    #Calculate a periodogram for the data
    pgram = TSA::periodogram(imputeTS::na.kalman(y), plot = F)
    pgram = data.table::data.table(freq = pgram$freq, spec = pgram$spec, period = 1/pgram$freq)[order(-spec), ]
    pgram[, `:=`("d", (spec)/mean(spec, na.rm = T))]
    pgram = pgram[period < length(y), ]

    #Calculate a periodogram for random data
    pgram_base = TSA::periodogram(rnorm(10000), plot = F)
    pgram_base = data.table::data.table(freq = pgram_base$freq, spec = pgram_base$spec, period = 1/pgram_base$spec)[order(-spec), ]
    pgram_base[, `:=`("d", (spec)/mean(spec, na.rm = T))]

    #Find periods that have a significant spectral number
    periods = pgram[which(sapply(1:nrow(pgram), function(x) {
      nrow(pgram_base[d > pgram[x, ]$d, ])/nrow(pgram_base)
    }) <= level), ]$period
    pgram = pgram[, `:=`("d", NULL)][order(period), ]
    pgram[, `:=`("significant", ifelse(period %in% periods,
                                       T, F))]
    decomp = "trend"
    if(length(periods) > 0){
      #Check for seasonality using the periodogram or the unit root test on the frequency differences
      if(freq > 1 & any(abs(periods - freq) < freq/2)){
        decomp = paste0(decomp, "-seasonal")
      }

      #Check for longer term cycle
      if(any(periods > freq * 3)){
        decomp = paste0(decomp, "-cycle")
      }
    }
    if(decomp == "trend"){
      decomp = "trend-noise"
    }
    rm(pgram_base, periods)
  }

  #Define the trend specifications to estiamte
  if(is.null(trend_spec)) {
    iter = c("rw", "rwd", "2rw")
  }else {
    if(!trend_spec %in% c("rw", "rwd", "2rw")){
      stop("trend_spec must be 'rw', 'rwd', or '2rw'.")
    }
    iter = trend_spec
  }

  comb = function(DT1, DT2) {
    return(rbind(DT1, DT2, use.names = T, fill = T))
  }
  cl = parallel::makeCluster(min(c(parallel::detectCores(), length(iter))))
  doSNOW::registerDoSNOW(cl)
  invisible(snow::clusterCall(cl, function(x) .libPaths(x), .libPaths()))
  # `%fun%` = foreach::`%dopar%`
  # fit = foreach::foreach(i = iter, .combine = "comb", .packages = c("data.table", "Matrix", "maxLik", "imputeTS"), .export = c("SSmodel")) %fun% {
  i = iter[1]
    #Set up the initial values
    if(i == "rw"){
      par = c(sig_t = sqrt(1/3 * var(diff(y), na.rm = T)))
    }else if(i %in% c("rwd", "2rw")){
      par = c(sig_t = sqrt(1/7 * var(diff(diff(y)), na.rm = T)))
      if(i == "rwd"){
        par = c(par, sig_m = unname(par["sig_t"]))
      }
    }
    if(grepl("seasonal", decomp)){
      if(freq %in% c(1, 4, 12)){
        seas_freqs = 1:(floor(freq)/2)
      }else if(freq == 52){
        seas_freqs = c(1, 2, 3, 4, 8, 12, 16, 20, 24)
      }else if(freq == 365){
        seas_freqs = c(1, 2, 3, 4, 5, 6, 7, 14, 21, 30, 60, 90, 120, 150, 182)
      }
      par = c(par, sig_j = unname(rep(par["sig_t"]/(2 * length(seas_freqs)), length(seas_freqs))))
      names(par)[grepl("sig_j", names(par))] = paste0("sig_j", seas_freqs)
    }
    if(grepl("cycle", decomp)){
      par = c(par, lambda = log((2 * pi/(freq * 5))/(pi - 2 * pi/(freq * 5))), rho = log((0.5)/(1 - 0.5)), sig_c = unname(par["sig_t"]/2))
    }
    if(grepl("seasonal|cycle", decomp)){
      par = c(par, sig_e = unname(par["sig_t"]))
    }
# 
#     #Set any fixed parameters
#     fixed = NULL
#     if(det_obs == T){
#       par["sig_e"] = 0
#       fixed = c(fixed, "sig_e")
#     }
#     if(det_trend == T){
#       par["sig_t"] = 0
#       fixed = c(fixed, "sig_t")
#     }
#     if(det_seasonality == T){
#       par[grepl("sig_j", names(par))] = 0
#       fixed = c(fixed, names(par)[grepl("sig_j", names(par))])
#     }
#     if(det_cycle == T){
#       par["sig_c"] = 0
#       fixed = c(fixed, "sig_c")
#     }
#     if(det_drift == T){
#       par["sig_m"] = T
#       fixed = c(fixed, "sig_m")
#     }

  #   #Define the objective function
  #   objective = function(par, na_locs, freq, decomp, trend_spec, init = NULL){
  #     yt = matrix(get("y"), nrow = 1)
  #     sp = SSmodel(par, yt, freq, decomp, trend_spec, init)
  #     ans = kalman_filter(B0 = matrix(sp$B0, ncol = 1), P0 = sp$P0, Dt = sp$Dt, At = sp$At, Ft = sp$Ft, Ht = sp$Ht, Qt = sp$Qt, Rt = sp$Rt, yt = yt)
  #     if (!is.null(na_locs)) {
  #       fc = sp$Ht %*% ans$B_tt
  #       yt[, na_locs] = fc[, na_locs]
  #       ans = kalman_filter(matrix(sp$B0, ncol = 1), P0 = sp$P0, Dt = sp$Dt, At = sp$At, Ft = sp$Ft, Ht = sp$Ht, Qt = sp$Qt, Rt = sp$Rt, yt = yt)
  #       assign("y", c(yt), .GlobalEnv)
  #     }
  #     return(ans$loglik)
  #   }
  # 
  #   #Get initial values for Kalman Filter
  #   sp = SSmodel(par, y, freq, decomp, trend_spec = i)
  #   init = kalman_filter(B0 = matrix(sp$B0, ncol = 1), P0 = sp$P0, Dt = sp$Dt, At = sp$At, Ft = sp$Ft, Ht = sp$Ht, Qt = sp$Qt, Rt = sp$Rt, yt = matrix(y, nrow = 1))
  #   init = kalman_smoother(B_tl = init[["B_tl"]], B_tt = init[["B_tt"]], P_tl = init[["P_tl"]], P_tt = init[["P_tt"]], Ft = sp$Ft)
  #   init = list(B0 = init[["B_tt"]][, 1], P0 = init[["P_tt"]][, , 1])
  # 
  #   #Find anhy na values
  #   if(any(is.na(y))){
  #     na_locs = which(is.na(y))
  #   }else{
  #     na_locs = NULL
  #   }
  # 
  #   #Estimate the model
  #   out = tryCatch(maxLik::maxLik(logLik = objective,
  #                                 start = par, method = optim_methods[1], fixed = fixed,
  #                                 finalHessian = F, hess = NULL, control = list(printLevel = 2, iterlim = maxit), init = init, na_locs = na_locs,
  #                                 freq = freq, decomp = decomp, trend_spec = i),
  #                  error = function(err){
  #                    tryCatch(maxLik::maxLik(logLik = objective,
  #                                            start = par, method = optim_methods[min(c(2, length(optim_methods)))], fixed = fixed,
  #                                            finalHessian = F, hess = NULL, control = list(printLevel = 2, iterlim = maxit), init = init, na_locs = na_locs,
  #                                            freq = freq, decomp = decomp, trend_spec = i),
  #                             error = function(err){
  #                               tryCatch(maxLik::maxLik(logLik = objective, start = par, method = optim_methods[min(c(3, length(optim_methods)))], fixed = fixed,
  #                                                       finalHessian = F, hess = NULL, control = list(printLevel = 2, iterlim = maxit), init = init, na_locs = na_locs,
  #                                                       freq = freq, decomp = decomp, trend_spec = i),
  #                                        error = function(err){NULL})
  #                             })
  #                  })
  # 
  #   #Attempt to get convergence if it failed the first time
  #   if(!is.null(out)){
  #     trials = 1
  #     while (out$code != 0 & trials < maxtrials) {
  #       out2 = tryCatch(maxLik::maxLik(logLik = objective, start = coef(out), method = optim_methods[1], fixed = fixed,
  #                                      finalHessian = F, hess = NULL, control = list(printLevel = 2, iterlim = maxit), init = init, na_locs = na_locs, freq = freq,
  #                                      decomp = decomp, trend_spec = i),
  #                       error = function(err){
  #                         tryCatch(maxLik::maxLik(logLik = objective, start = coef(out), method = optim_methods[min(c(2, length(optim_methods)))], fixed = fixed,
  #                                                 finalHessian = F, hess = NULL, control = list(printLevel = 2, iterlim = maxit), init = init, na_locs = na_locs,
  #                                                 freq = freq, decomp = decomp, trend_spec = i),
  #                                  error = function(err){
  #                                    tryCatch(maxLik::maxLik(logLik = objective, start = coef(out), method = optim_methods[min(c(3, length(optim_methods)))], fixed = fixed,
  #                                                            finalHessian = F, hess = NULL, control = list(printLevel = 2, iterlim = maxit), init = init, na_locs = na_locs,
  #                                                            freq = freq, decomp = decomp, trend_spec = i),
  #                                             error = function(err){NULL})
  #                                  })
  #                       })
  #       #End the loop if no parameters changed or estimation failed
  #       if(!is.null(out2) & !all(coef(out) == coef(out2))){
  #         out = out2
  #         trials = trials + 1
  #       }else {
  #         break
  #       }
  #     }
  #   }
  #   rm(init)
  #   gc()
  # 
  #   if(!is.null(out)){
  #     #Retreive the model output
  #     # return(data.table::data.table(model = i, freq = freq, decomp = decomp, convergence = out$code, loglik = out$maximum,
  #     #                               matrix(coef(out), nrow = 1, dimnames = list(NULL, paste0("coef_", names(coef(out)))))))
  #     fit = data.table::data.table(model = i, freq = freq, decomp = decomp, convergence = out$code, loglik = out$maximum,
  #                                  matrix(coef(out), nrow = 1, dimnames = list(NULL, paste0("coef_", names(coef(out))))))
  #   }else{
  #     # return(NULL)
  #     fit = NULL
  #   }
  # # }
  # # snow::stopCluster(cl)
  # 
  # #Select the best model based on the maximum likelihood
  # model_selection = fit[loglik == max(loglik, na.rm = T), ]$model
  # fit = fit[loglik == max(loglik, na.rm = T), ]
  # fit = list(table = fit)
  # if(exists("pgram")){
  #   fit[["periodogram"]] = pgram
  # }else{
  #   fit[["periodogram"]] = NA
  # }
  return(fit)
}

#' Kalman filter an estimated model from tcs_decomp_estim
#' @param y Univariate time series of data values. May also be a 2 column data frame containing a date column.
#' @param model Structural time series model estimated using tcs_decomp_estim.
#' @param plot Logial, whether to plot the output or not.
#' @return List of data tables containing the filtered and smoothed series.
#' @examples
#' tcs_decomp_filter(y = DT[, c("date", "y")], model = model)
#' @author Alex Hubbard (hubbard.alex@gmail.com)
#' @export
tcs_decomp_filter = function(y, model, plot = F, select = NULL){
  y = data.table::as.data.table(y)
  .SD = data.table::.SD
  datecol = unlist(lapply(colnames(y), function(x){
    if(class(y[, c(x), with = F][[1]]) %in% c("Date")){
      return(x)
    }else{
      return(NULL)
    }
  }))
  if(length(datecol) == 1){
    #Detect the frequency
    dates = y[, datecol, with = F][[1]]
    y = y[, colnames(y)[colnames(y) != datecol], with = F][[1]]
    rm(datecol)
  }else if(length(datecol) > 1){
      stop("Too many date columns detected.")
  }else {
    stop("No date column detected. Include a date column or set the frequency.")
  }
  
  #Set the dates if it hasn't been given
  if(is.null(dates)){
    dates = 1:length(y)
  }
  
  #Find any missing values
  if(any(is.na(y))){
    na_locs = which(is.na(y))
  }else{
    na_locs = NULL
  }
  
  #Get model specifications
  freq = model$table$freq
  decomp = model$table$decomp
  trend_spec = model$table$model
  
  #Get the coefficients
  cs = unlist(model$table[, grepl("coef_", colnames(model$table)), with = F])
  cs = cs[!is.na(cs)]
  names(cs) = gsub("coef_", "", names(cs))
  
  #Filter the data
  sp = SSmodel(cs, y, freq, decomp, trend_spec)
  init = kalman_filter(B0 = matrix(sp$B0, ncol = 1), P0 = sp$P0, Dt = sp$Dt, At = sp$At, Ft = sp$Ft, Ht = sp$Ht, Qt = sp$Qt, Rt = sp$Rt, yt = matrix(y, nrow = 1))
  init = kalman_smoother(B_tl = init$B_tl, B_tt = init$B_tt, P_tl = init$P_tl, P_tt = init$P_tt, Ft = sp$Ft)
  init = list(B0 = init[["B_tt"]][, 1], P0 = init[["P_tt"]][, , 1])
  sp = SSmodel(cs, y, freq, decomp, trend_spec, init)
  ans = kalman_filter(B0 = matrix(sp$B0, ncol = 1), P0 = sp$P0, Dt = sp$Dt, At = sp$At, Ft = sp$Ft, Ht = sp$Ht, Qt = sp$Qt, Rt = sp$Rt, yt = matrix(y, nrow = 1))
  if(!is.null(na_locs)){
    fc = sp$Ht %*% ans$B_tt
    y[na_locs] = fc[, na_locs]
    ans = kalman_filter(B0 = matrix(sp$B0, ncol = 1), P0 = sp$P0, Dt = sp$Dt, At = sp$At, Ft = sp$Ft, Ht = sp$Ht, Qt = sp$Qt, Rt = sp$Rt, yt = matrix(y, nrow = 1))
  }
  rownames(ans$B_tt) = rownames(sp$Ft)
  smooth = kalman_smoother(B_tl = ans$B_tl, B_tt = ans$B_tt, P_tl = ans$P_tl, P_tt = ans$P_tt, Ft = sp$Ft)[["B_tt"]]
  rownames(smooth) = rownames(ans$B_tt)
  
  #Retrieve the model output
  iter = c("filter", "smooth")
  cl = parallel::makeCluster(min(c(parallel::detectCores(), length(iter))))
  doSNOW::registerDoSNOW(cl)
  invisible(snow::clusterCall(cl, function(x) .libPaths(x), .libPaths()))
  `%fun%` = foreach::`%dopar%`
  preret = foreach::foreach(i = iter, .packages = c("data.table")) %fun% {
    if(i == "filter"){
      B_tt = ans$B_tt
    }
    else if(i == "smooth"){
      B_tt = smooth
    }
    
    #Get the unobserved series
    series = data.table::data.table(t(B_tt))
    series.l = copy(series)
    series.l[, `:=`(colnames(series.l), lapply(.SD, shift, type = "lag", n = 1)), .SDcols = colnames(series.l)]
    
    #Get the model errors
    errors = t(as.matrix(series)) - sp$Ft %*% t(as.matrix(series.l))
    
    #Get the trend and its errors
    trend = c(sp$Ht[grepl("Tt|Mt", colnames(sp$Ht))] %*% B_tt[which(grepl("Tt|Mt", rownames(B_tt))), ])
    trend_error = c(sp$Ht[grepl("Tt|Mt", colnames(sp$Ht))] %*% errors[grepl("Tt|Mt", rownames(errors)), ])
    
    #Get the cycle and its errors
    if(any(grepl("Ct", rownames(errors)))){
      cycle = c(sp$Ht[grepl("Ct", colnames(sp$Ht))] %*% B_tt[which(grepl("Ct", rownames(B_tt))), ])
      cycle_error = c(sp$Ht[grepl("Ct", colnames(sp$Ht))] %*% errors[grepl("Ct", rownames(errors)), ])
    }else{
      cycle = rep(0, nrow(series))
      cycle_error = rep(0, ncol(errors))
    }
    
    #Get the seasonality and its errors
    if(any(grepl("St", rownames(errors)))){
      seasonal = seasonal = c(sp$Ht[grepl("St", colnames(sp$Ht))] %*% B_tt[which(grepl("St", rownames(B_tt))), ])
      seasonal_error = c(sp$Ht[grepl("St", colnames(sp$Ht))] %*% errors[grepl("St", rownames(errors)), ])
    }else{
      seasonal = rep(0, nrow(series))
      seasonal_error = rep(0, ncol(errors))
    }
    
    #Get the observation and total errors
    observation_error = y - c(sp$Ht %*% t(as.matrix(series)))
    total_error = trend_error + cycle_error + seasonal_error + observation_error
    
    #Combine the filtered series
    toret = data.table::data.table(y = y, trend = trend, cycle = cycle, seasonal = seasonal, 
                                   trend_error = trend_error, cycle_error = cycle_error, seasonal_error = seasonal_error, 
                                   observation_error = observation_error, total_error = total_error)
    
    #Filter out the errors 
    toret[, `:=`("trend_pred", trend - trend_error)]
    toret[, `:=`("cycle_pred", cycle - cycle_error)]
    toret[, `:=`("seasonal_pred", seasonal - seasonal_error)]
    
    #Calculate adusted series
    toret[, `:=`("seasonal_adjusted", y - seasonal)]
    toret[, `:=`("cycle_adjusted", y - cycle)]
    toret[, `:=`("seasonal_cycle_adjusted", y - seasonal - cycle)]
    toret[, `:=`("seasonal_cycle", cycle + seasonal)]
    return(toret)
  }
  names(preret) = iter
  
  #Combine the filtered and smoothed series
  final = rbind(data.table(method = "filter", date = dates, preret[["filter"]]), 
                data.table(method = "smooth", date = dates, preret[["smooth"]]), 
                use.names = T, fill = T)
  
  if(plot == T) {
    for(i in c("filter", "smooth")){
      g1 = ggplot2::ggplot(data.table::melt(final[method == i, ], id.vars = "date", measure.vars = c("y", "trend"))) + 
        ggplot2::ggtitle("Actual vs Trend") + 
        ggplot2::geom_line(ggplot2::aes(x = date, y = value, group = variable, color = variable)) + ggplot2::theme_minimal() + 
        ggplot2::guides(color = ggplot2::guide_legend(title = NULL)) + 
        ggplot2::theme(legend.position = "bottom")
      title = c(ifelse(grepl("cycle", decomp), "Cycle", ""), ifelse(grepl("seasonal", decomp), "Seasonal", ""), "Observation Error")
      title = title[title != ""]
      g2 = ggplot2::ggplot(data.table::melt(final[method == i, ], id.vars = "date", measure.vars = c("cycle", "seasonal", "observation_error"))) + 
        ggplot2::ggtitle(paste(title, collapse = ", ")) + 
        ggplot2::geom_line(ggplot2::aes(x = date, y = value, group = variable, color = variable)) + 
        ggplot2::theme_minimal() + ggplot2::guides(color = ggplot2::guide_legend(title = NULL)) + 
        ggplot2::theme(legend.position = "bottom")
      if(!all(is.na(model[["periodogram"]]))){
        g3 = ggplot2::ggplot(model[["periodogram"]][, .(`Significant Frequencies` = significant, spec = mean(spec)), by = round(period)]) + 
          ggplot2::ggtitle("Periodogram") + ggplot2::scale_x_continuous(name = "Period") + 
          ggplot2::scale_y_continuous(name = "Spec") + 
          ggplot2::geom_col(ggplot2::aes(x = round, y = spec, fill = `Significant Frequencies`, color = `Significant Frequencies`), 
                            position = "dodge", width = 0.75) + 
          ggplot2::theme_minimal() + ggplot2::guides(color = F) + 
          ggplot2::theme(legend.position = "bottom")
        gridExtra::grid.arrange(g1, g2, g3, layout_matrix = rbind(c(1, 2), c(3, 3)), top = grid::textGrob(i, gp = grid::gpar(fontsize = 20, font = 3)))
      }else{
        gridExtra::grid.arrange(g1, g2, nrow = 1, top = grid::textGrob(i, gp = grid::gpar(fontsize = 20, font = 3)))
      }
      Sys.sleep(0.05)
    }
  }
  snow::stopCluster(cl)
  return(final)
}


########
#Call these to build the package
#devtools::document()
#devtools::build_vignettes()
#devtools::install()
#library(projectmap)
#git config remote.origin.url git@github.com:opendoor-labs/TCSDecomp.git

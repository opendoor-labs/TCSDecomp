#' State space model
#' Creates a state space model in list form
#' yt = Ht %*% Bt + e_t
#' Bt = Ft %*% B_{t=1} + u_t
#'
#' @param par Vector of named parameter values
#' @param yt Univariate time series of data values
#' @param freq Seasonality of the data (1, 4, 12, 53, 365)
#' @param decomp Decomposition model ("tend-cycle-seasonal", "trend-seasonal", "trend-cycle", "trend-noise")
#' @param int_order Order of integration (0, 1, 2)
#' @param tend Trend specification ("", "rw", "hp")
#' @param init Initial state values for the Kalman filter
#' @return List of space space matrices
#' @examples
#' SSmodel(par, y, freq, decomp, int_order, trend, init)
#' @author Alex Hubbard (hubbard.alex@gmail.com)
#' @export
SSmodel = function(par, yt, freq, decomp, int_order, trend, init = NULL){
  yt = matrix(yt, nrow = 1)
  #Define the standard deviation of the observation equation
  if(decomp == "trend-noise"){
    sig_e = ifelse(int_order == 0, sd(yt),
                   ifelse(int_order == 1, sd(diff(t(yt))),
                          ifelse(int_order == 2, sd(diff(diff(t(yt)))))))
  }else if(decomp %in% c("trend-cycle", "trend-seasonal", "trend-cycle-seasonal", "trend-seasonal-cycle")){
    sig_e = par["sig_e"]
  }else{
    stop("decomp must be one of 'trend-noise', 'trend-cycle', 'trend-seasonal', 'trend-cycle-seasonal', or 'trend-seasonal-cycle'.")
  }
  
  #Define the transition and observation matrix based on the order of integration
  if(int_order %in% c(0, "0")){
    #T_t = phi1*T_{t-1} + phi2*T_{t-2} + e_t, e_t ~ N(0, sig_t^2)
    
    #Transition matrix
    Fm = rbind(c(par[grepl("phi", names(par))]),
               cbind(diag(1, nrow = length(par[grepl("phi", names(par))]) - 1),
                     matrix(0, nrow = length(par[grepl("phi", names(par))]) - 1)))
    colnames(Fm) = paste0("Tt", 1:ncol(Fm))
    rownames(Fm) = paste0("Tt", 0:(ncol(Fm) - 1))
    
    #Observation matrix
    Hm = matrix(c(1, rep(0, nrow(Fm) - 1)), nrow = 1)
    colnames(Hm) = rownames(Fm)
  }else if(int_order %in% c(1, "1")){
    #T_t = T_{t-1} + e_t, e_t ~ N(0, sig_t^2)
    
    #Transition matrix
    Fm = matrix(c(1), nrow = 1, ncol = 1)
    colnames(Fm) = "Tt1"
    rownames(Fm) = "Tt0"
    
    #Observation matrix
    Hm = matrix(c(1), nrow = 1)
    colnames(Hm) = rownames(Fm)
  }else if(int_order %in% c(2, "2")){
    #trend = "rw": T_t = M_{t-1} + T_{t-1} + e_t, e_t ~ N(0, sig_t^2)
    #              M_t = M_{t-1} + n_t, n_t ~ N(0, sig_m^2)
    #trend = "hp": T_t = 2T_{t-1} - T_{t-2} + e_t, e_t ~ N(0, sig_t^2)
    
    #Transition matrix
    if(trend == "rw"){
      Fm = rbind(c(1, 1), c(0, 1))
      colnames(Fm) = c("Tt1", "Mt1")
      rownames(Fm) = c("Tt0", "Mt0")
    }else if(trend == "hp"){
      Fm = rbind(c(2, -1), c(1, 0))
      colnames(Fm) = c("Tt1", "Tt2")
      rownames(Fm) = c("Tt0", "Tt1")
    }
    
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
  if(grepl("cycle", decomp)){
    rho = 1/(1 + exp(-par["rho"])) #Constrain rho (the dampening) to be between 0 and 1
    lambda = pi/(1 + exp(-par["lambda"])) #Constrain lambda (the period) to between 0 and pi
    colnames = colnames(Fm)
    rownames = rownames(Fm)
    Cm = rbind(c(rho*cos(lambda), rho*sin(lambda)),
               c(-rho*sin(lambda), rho*cos(lambda)))
    Fm = as.matrix(Matrix::bdiag(Fm, Cm))
    colnames(Fm) = c(colnames, "Ct1", "Ct2")
    rownames(Fm) = c(rownames, "Ct0", "Ct1")
    Hm = cbind(Hm, matrix(c(1, 0), nrow = 1))
    colnames(Hm) = rownames(Fm)
    Qm = as.matrix(Matrix::bdiag(Qm, diag(par[grepl("sig_c", names(par))])))
    colnames(Qm) = rownames(Qm) = rownames(Fm)
  }
  
  #Define the seasonal component
  if(grepl("seasonal", decomp)){
    jiter = as.numeric(unique(gsub("[[:alpha:]]|[[:punct:]]", "", names(par)[grepl("sig_j", names(par))])))
    for(j in jiter){
      colnames = colnames(Fm)
      rownames = rownames(Fm)
      Sm = rbind(c(cos(2*pi*j/freq), sin(2*pi*j/freq)),
                 c(-sin(2*pi*j/freq), cos(2*pi*j/freq)))
      Fm = as.matrix(Matrix::bdiag(Fm, Sm))
      colnames(Fm) = c(colnames, paste0("Stl", j), paste0("Stls", j))
      rownames(Fm) = c(rownames, paste0("St", j), paste0("Sts", j))
    }
    Hm = cbind(Hm, matrix(rep(c(1, 0), length(jiter)), nrow = 1))
    colnames(Hm) = rownames(Fm)
    Qm = as.matrix(Matrix::bdiag(Qm, diag(par[grepl("sig_j", names(par))])))
    colnames(Qm) = rownames(Qm) = rownames(Fm)
  }
  
  #Transition equation intercept matrix
  Dm = matrix(0, nrow = nrow(Fm))
  
  #Observaton equation intercept matrix
  Am = matrix(0, nrow = 1, ncol = 1)
  
  #Observation equation error covariance matrix
  Rm = matrix(sig_e^2, nrow = 1, ncol = 1)
  
  #Initial guess for unobserved vector
  if(is.null(init)){
    B0 = c(rep(mean(c(yt[1:ifelse(freq == 1, round(length(y)/4), freq)]), na.rm = T), length(rownames(Fm)[grepl("Tt", rownames(Fm))])),
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
  
  #Initial guess for variance of unobserved vector
  if(is.null(init)){
    P0 = diag(100, nrow = nrow(Fm), ncol = nrow(Fm))
    rownames(P0) = colnames(P0) = rownames(Fm)
    P0["Tt0", "Tt0"] = var(c(yt[1:ifelse(freq == 1, round(length(y)/4), freq)]), na.rm = T)
    if(trend == "hp"){
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
#' @param int_order Order of integration (0, 1, 2)
#' @param blend Blend the model specifications using a minimum mean squared error weighted sum
#' @param level Significance level of statistical tests (0.01, 0.05, 0.1)
#' @param ur_test Unit root test name ("kpss", "adf", "pp")
#' @param ur_type Unit root test deterministi type ("level", "trend")
#' @param optim_methods Vector of 1 to 3 optimization methods in order of preference ("NR", "BFGS", "CG", "BHHH", or "SANN")
#' @param maxit Maximum number of iterations for the optimization
#' @param maxtrials Maximum number of optimization trials to get convergence
#' @return List of estimation values including coefficients, convergence code, datea frequency, decomposition used, and trend specification selected.
#' @examples
#' kf_decomp_estim(y = DT[, c("date", "y")])
#' @author Alex Hubbard (hubbard.alex@gmail.com)
#' @export
kf_decomp_estim = function(y, freq = NULL, decomp = NULL, int_order = NULL,
                           blend = F, level = 0.01, ur_test = "adf", ur_type = "level",
                           optim_methods = c("BFGS", "CG", "NM"), maxit = 1000, maxtrials = 10){
  dates = NULL
  if(is.null(freq)){
    #stop("Must provide freq as numeric (1 for annual, 4 for quarterly, 12 for monthly, 52 for weekly, 365 for daily).")
    y = data.table::as.data.table(y)
    `.SD` = data.table::.SD
    if(any(unlist(y[, lapply(.SD, function(x){class(x) == "Date"})]))){
      #Detect the frequency
      datecol = names(which(unlist(y[, lapply(.SD, function(x){class(x) == "Date"})])))
      if(length(datecol) > 1){
        stop("Too many date columns. Include only 1 date column or set the frequency manually.")
      }
      datediffs = unique(diff(y[, datecol, with = F][[1]]))
      freq = datediffs[which.max(tabulate(match(diff(y[, datecol, with = F][[1]]), datediffs)))]
      freq = c(365, 52, 12, 4, 1)[which.min(abs(freq - c(1, 7, 30, 90, 365)))]
      dates = y[, datecol, with = F][[1]]
      y = y[, colnames(y)[colnames(y) != datecol], with = F][[1]]
      rm(datediffs, datecol)
    }else{
      stop("No date column detected. Include a date column or set the frequency.")
    }
  }else{
    if(!is.numeric(freq)){
      stop("Must provide freq as numeric (1 for annual, 4 for quarterly, 12 for monthly, 52 for weekly, 365 for daily).")
    }else if(!freq %in% c(1, 4, 12, 52, 365)){
      stop("Must provide freq as numeric (1 for annual, 4 for quarterly, 12 for monthly, 52 for weekly, 365 for daily).")
    }
  }
  if(level < 0.01 | level > 0.1){
    stop("level must be between 0.01 and 0.1.")
  }
  if(!ur_test %in% c("kpss", "adf", "pp")){
    stop("ur_test must be one of 'kpss', 'adf', or 'pp'.")
  }
  if(!ur_type %in% c("level", "trend")){
    stop("ur_type must be 'level' or 'trend'.")
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
    pgram = TSA::periodogram(y, plot = F)
    pgram = data.table::data.table(freq = pgram$freq, spec = pgram$spec, period = 1/pgram$freq)[order(-spec), ]
    pgram[, "d" := (spec)/mean(spec, na.rm = T)]
    pgram = pgram[period < length(y), ]
    
    #Calculate a periodogram for random data
    pgram_base = TSA::periodogram(rnorm(10000), plot = F)
    pgram_base = data.table::data.table(freq = pgram_base$freq, spec = pgram_base$spec, period = 1/pgram_base$spec)[order(-spec), ]
    pgram_base[, "d" := (spec)/mean(spec, na.rm = T)]
    
    #Find periods that have a significant spectral number
    periods = pgram[which(sapply(1:nrow(pgram), function(x){
      nrow(pgram_base[d > pgram[x, ]$d, ])/nrow(pgram_base)
    }) <= level), ]$period
    pgram = pgram[, "d" := NULL][order(period), ]
    pgram[, "significant" := ifelse(period %in% periods, T, F)]
    
    decomp = "trend"
    if(length(periods) > 0){
      #Check for seasonality using the periodogram or the unit root test on the frequency differences
      if(freq > 1 & any(abs(periods - freq) < freq/2)){
        decomp = paste0(decomp, "-seasonal")
      }
      #Check for longer term cycle
      if(any(periods > freq*3)){
        decomp = paste0(decomp, "-cycle")
      }
    }
    if(decomp == "trend"){
      decomp = "trend-noise"
    }
    rm(pgram_base, periods)
  }
  
  objective = function(par, yt, freq, decomp, int_order, trend, init = NULL){
    yt = matrix(yt, nrow = 1)
    sp = SSmodel(par, yt, freq, decomp, int_order, trend)
    ans = kf(B0 = matrix(sp$B0, ncol = 1), P0 = sp$P0, Dt = sp$Dt, At = sp$At, Ft = sp$Ft, Ht = sp$Ht, Qt = sp$Qt, Rt = sp$Rt, yt = yt)
    #maxLik maximes the log likelihood, optim minizes the log likelihood
    return(ans$loglik)#return(-ans$logLik)
  }
  
  #Define the trend specifications to estiamte
  if(is.null(int_order)){
    int_order = as.character(forecast::ndiffs(y, alpha = level, test = ur_test, type = ur_type))
    if(int_order %in% c("1", "2")){
      iter = c("1", "2rw", "2hp")
    }else{
      iter = "0"
    }
  }else{
    if(!int_order %in% c(0, 1, 2, "2rw", "2hp")){
      stop("int_order must be NULL, 0, 1, or 2.")
    }
    if(int_order == 2){
      iter == c("2rw", "2hp")
    }else{
      iter = as.character(int_order)
    }
  }
  
  comb = function(DT1, DT2){
    return(rbind(DT1, DT2, use.names = T, fill = T))
  }
  cl = parallel::makeCluster(min(c(parallel::detectCores(), length(iter))))
  doSNOW::registerDoSNOW(cl)
  invisible(snow::clusterCall(cl, function(x) .libPaths(x), .libPaths()))
  `%fun%` = foreach::`%dopar%`
  fit = foreach::foreach(i = iter, .combine = "comb", .packages = c("data.table", "Matrix", "maxLik"), .export = c("SSmodel")) %fun% {
    #Set up the initial values
    if(gsub("rw|hp", "", i) == "0"){
      par = forecast::auto.arima(y, d = 0, max.q = 0, allowdrift = F)
      if(sum(coef(par)[!grepl("sar", names(coef(par)))]) < 1){
        par = c(phi = unname(coef(par)[!grepl("sar", names(coef(par)))]))
      }else{
        par = c(phi = rep(0.5/length(coef(par)[!grepl("sar", names(coef(par)))]), length(coef(par)[!grepl("sar", names(coef(par)))])))
      }
      par = c(par, sig_t = sqrt(var(y)*(1 - sum(par^2))/(2 - sum(par)^2)))
      if(is.na(par["sig_t"])){
        par["sig_t"] = sqrt(1/2*var(y))
      }
    }else if(gsub("rw|hp", "", i) == "1"){
      par = c(sig_t = sqrt(1/3*var(diff(y))))
    }else if(gsub("rw|hp", "", i) == "2"){
      par = c(sig_t = sqrt(1/7*var(diff(diff(y)))))
      if(i == "2rw"){
        par = c(par, sig_m = unname(par["sig_t"]))
      }
    }
    if(grepl("seasonal", decomp)){
      if(freq %in% c(1, 4, 12)){
        seas_freqs = 1:(floor(freq)/2)
      }else if(freq == 52){
        seas_freqs = c(1, 2, 3, 4, 8, 12, 16, 20, 24)
      }else if(freq == 365){
        seas_freqs = c(1, 2, 3, 4, 5, 6, 7, 14, 21, 30, 60, 90, 120, 150, 180)
      }
      par = c(par, sig_j = unname(rep(par["sig_t"]/(2*length(seas_freqs)), 2*length(seas_freqs))))
      names(par)[grepl("sig_j", names(par))] = unlist(lapply(seas_freqs, function(x){paste0(c("sig_j", "sig_js"), x)}))
    }
    if(grepl("cycle", decomp)){
      par = c(par, lambda = log((2*pi/(freq*5))/(pi - 2*pi/(freq*5))), rho = log((0.5)/(1 - 0.5)), sig_c = unname(par["sig_t"]/2), sig_cs = unname(par["sig_t"]/2))
    }
    if(grepl("seasonal|cycle", decomp)){
      par = c(par, sig_e = unname(par["sig_t"]))
    }
    
    #Get initial values for Kalman Filter
    sp = SSmodel(par, y, freq, decomp, int_order = gsub("rw|hp", "", i), trend = gsub("[[:digit:]]", "", i))
    init = kf(B0 = matrix(sp$B0, ncol = 1), P0 = sp$P0, Dt = sp$Dt, At = sp$At, Ft = sp$Ft, Ht = sp$Ht, Qt = sp$Qt, Rt = sp$Rt, yt = matrix(y, nrow = 1))
    init = kf_smoother(B_TL = init[["B_TL"]], B_TT = init[["B_TT"]], P_TL = init[["P_TL"]], P_TT = init[["P_TT"]], Ft = sp$Ft)
    init = list(B0 = init[["B_TT"]][, 1], P0 = init[["P_TT"]][, , 1])
    
    #Estimate the model
    out = tryCatch(maxLik::maxLik(logLik = objective, start = par, method = optim_methods[1],
                                  finalHessian = F, hess = NULL, control = list(printLevel = 2, iterlim = maxit), init = init,
                                  yt = y, freq = freq, decomp = decomp, int_order = gsub("rw|hp", "", i), trend = gsub("[[:digit:]]", "", i)),
                   error = function(err){
                     tryCatch(maxLik::maxLik(logLik = objective, start = par, method = optim_methods[min(c(2, length(optim_methods)))],
                                             finalHessian = F, hess = NULL, control = list(printLevel = 2, iterlim = maxit), init = init,
                                             yt = y, freq = freq, decomp = decomp, int_order = gsub("rw|hp", "", i), trend = gsub("[[:digit:]]", "", i)),
                              error = function(err){
                                tryCatch(maxLik::maxLik(logLik = objective, start = par, method = optim_methods[min(c(3, length(optim_methods)))],
                                                        finalHessian = F, hess = NULL, control = list(printLevel = 2, iterlim = maxit), init = init,
                                                        yt = y, freq = freq, decomp = decomp, int_order = gsub("rw|hp", "", i), trend = gsub("[[:digit:]]", "", i)),
                                         error = function(err){NULL})
                              })
                   })
    
    #Attempt to get convergence if it failed the first time
    if(!is.null(out)){
      trials = 1
      while(out$code != 0 & trials < maxtrials){
        out2 = tryCatch(maxLik::maxLik(logLik = objective, start = coef(out), method = optim_methods[1],
                                       finalHessian = F, hess = NULL, control = list(printLevel = 2, iterlim = maxit), init = init,
                                       yt = y, freq = freq, decomp = decomp, int_order = gsub("rw|hp", "", i), trend = gsub("[[:digit:]]", "", i)),
                        error = function(err){
                          tryCatch(maxLik::maxLik(logLik = objective, start = coef(out), method = optim_methods[min(c(2, length(optim_methods)))],
                                                  finalHessian = F, hess = NULL, control = list(printLevel = 2, iterlim = maxit), init = init,
                                                  yt = y, freq = freq, decomp = decomp, int_order = gsub("rw|hp", "", i), trend = gsub("[[:digit:]]", "", i)),
                                   error = function(err){
                                     tryCatch(maxLik::maxLik(logLik = objective, start = coef(out), method = optim_methods[min(c(3, length(optim_methods)))],
                                                             finalHessian = F, hess = NULL, control = list(printLevel = 2, iterlim = maxit), init = init,
                                                             yt = y, freq = freq, decomp = decomp, int_order = gsub("rw|hp", "", i), trend = gsub("[[:digit:]]", "", i)),
                                              error = function(err){NULL})
                                   })
                        })
        #End the loop if no parameters changed or estimation failed
        if(!is.null(out2) & !all(coef(out) == coef(out2))){
          out = out2
          trials = trials + 1
        }else{
          break
        }
      }
    }
    rm(init)
    gc()
    
    if(!is.null(out)){
      #Retreive the model output
      return(data.table::data.table(model = i, freq = freq, decomp = decomp, convergence = out$code, loglik = out$maximum,
                                    matrix(coef(out), nrow = 1, dimnames = list(NULL, paste0("coef_", names(coef(out))))),
                                    trend = ifelse(i == "0", "autoregressive",
                                                   ifelse(i == "1", "1st order random walk",
                                                          ifelse(i == "2rw", "random walk with random walk drift",
                                                                 ifelse(i == "2hp", "2nd order random walk", NA))))))
    }else{
      return(NULL)
    }
  }
  snow::stopCluster(cl)
  if(blend == T){
    test = do.call("cbind", lapply(fit$model, function(x){
      cs = unlist(fit[model == x, grepl("coef_", colnames(fit)), with = F])
      cs = cs[!is.na(cs)]
      names(cs) = gsub("coef_", "", names(cs))
      sp = SSmodel(cs, y, freq, decomp, int_order = gsub("rw|hp", "", x), trend = gsub("[[:digit:]]", "", x))
      ans = kf(B0 = matrix(sp$B0, ncol = 1), P0 = sp$P0, Dt = sp$Dt, At = sp$At, Ft = sp$Ft, Ht = sp$Ht, Qt = sp$Qt, Rt = sp$Rt, yt = matrix(y, nrow = 1))
      rownames(ans$B_TT) = rownames(sp$Ft)
      return(t(sp$Ht %*% ans$B_TT))
    }))
    colnames(test) = names(fit)
    
    blend_objective = function(par){
      err = y - par[1]*test[, 1] - par[2]*test[, 2] - (1 - par[1] - par[2])*test[, 3]
      return(sum(err)^2)
    }
    blended = optim(par = c(1/3, 1/3), fn = blend_objective, method = "L-BFGS-B",
                    lower = c(0, 0), upper = c(1, 1))
    weights = c(blended$par, 1 - sum(blended$par))
    names(weights) = names(fit)
    rm(blended, blend_objective, test)
  }else{
    model_selection = fit[loglik == max(loglik, na.rm = T), ]$model
    weights = rep(0, length(fit))
    names(weights) = names(fit)
    weights[model_selection] = 1
    fit = fit[loglik == max(loglik, na.rm = T), ]
  }
  fit = merge(fit, data.table::data.table(model = names(weights), weight = weights), by = "model")
  fit = list(table = fit)
  if(exists("pgram")){
    fit[["periodogram"]] = pgram
  }else{
    fit[["periodogram"]] = NA
  }
  return(fit)
}

#' Kalman filter an estimated model from kf_decomp_estim
#' @param y Univariate time series of data values. May also be a 2 column data frame containing a date column.
#' @param model Structural time series model estimated using kf_decomp_estim.
#' @param plot Logial, whether to plot the output or not.
#' @param select If NULL, then use the model selected by kf_dceomp_estim; otherwise, use the one specified with select
#' @return List of data tables containing the filtered and smoothed series.
#' @examples
#' kf_decomp_filter(y = DT[, c("date", "y")], model = model)
#' @author Alex Hubbard (hubbard.alex@gmail.com)
#' @export
kf_decomp_filter = function(y, model, plot = F, select = NULL){
  y = data.table::as.data.table(y)
  `.SD` = data.table::.SD
  if(any(unlist(y[, lapply(.SD, function(x){class(x) == "Date"})]))){
    #Detect the frequency
    datecol = names(which(unlist(y[, lapply(.SD, function(x){class(x) == "Date"})])))
    if(length(datecol) > 1){
      stop("Too many date columns. Include only 1 date column or set the frequency manually.")
    }
    dates = y[, datecol, with = F][[1]]
    y = y[, colnames(y)[colnames(y) != datecol], with = F][[1]]
    rm(datecol)
  }else{
    stop("No date column detected. Include a date column or set the frequency.")
  }
  
  #Set the dates if it hasn't been given
  if(is.null(dates)){
    dates = 1:length(y)
  }
  
  #Get the filtered and smoothed series
  weights = model$table$weight
  names(weights) = model$table$model
  miter = names(weights)[weights != 0]
  if(!is.null(select)){
    miter = select
    weights[names(weights) == select] = 1
    weights[names(weights) != select] = 0
  }
  freq = unique(model$table$freq)
  decomp = unique(model$table$decomp)
  
  cl = parallel::makeCluster(min(c(parallel::detectCores(), length(miter))))
  doSNOW::registerDoSNOW(cl)
  invisible(snow::clusterCall(cl, function(x) .libPaths(x), .libPaths()))
  `%fun%` = foreach::`%do%`
  ret = foreach::foreach(n = miter, .packages = c("data.table")) %fun% {
    cs = unlist(model$table[model == n, grepl("coef_", colnames(model$table)), with = F])
    cs = cs[!is.na(cs)]
    names(cs) = gsub("coef_", "", names(cs))
    
    sp = SSmodel(cs, y, freq, decomp, int_order = gsub("rw|hp", "", n), trend = gsub("[[:digit:]]", "", n))
    init = kf(B0 = matrix(sp$B0, ncol = 1), P0 = sp$P0, Dt = sp$Dt, At = sp$At, Ft = sp$Ft, Ht = sp$Ht, Qt = sp$Qt, Rt = sp$Rt, yt = matrix(y, nrow = 1))
    init = kf_smoother(B_TL = init$B_TL, B_TT = init$B_TT, P_TL = init$P_TL, P_TT = init$P_TT, Ft = sp$Ft)
    init = list(B0 = init[["B_TT"]][, 1], P0 = init[["P_TT"]][, , 1])
    sp = SSmodel(cs, y, freq, decomp, int_order = gsub("rw|hp", "", n), trend = gsub("[[:digit:]]", "", n), init = init)
    ans = kf(B0 = matrix(sp$B0, ncol = 1), P0 = sp$P0, Dt = sp$Dt, At = sp$At, Ft = sp$Ft, Ht = sp$Ht, Qt = sp$Qt, Rt = sp$Rt, yt = matrix(y, nrow = 1))
    rownames(ans$B_TT) = rownames(sp$Ft)
    smooth = kf_smoother(B_TL = ans$B_TL, B_TT = ans$B_TT, P_TL = ans$P_TL, P_TT = ans$P_TT, Ft = sp$Ft)[["B_TT"]]
    rownames(smooth) = rownames(ans$B_TT)
    
    iter = c("filter", "smooth")
    `%fun2%` = foreach::`%dopar%`
    preret = foreach::foreach(i = iter, .packages = c("data.table")) %fun2% {
      if(i == "filter"){
        B_TT = ans$B_TT
      }else if(i == "smooth"){
        B_TT = smooth
      }
      series = data.table::data.table(t(B_TT))
      series.l = copy(series)
      series.l[, colnames(series.l) := lapply(.SD, shift, type = "lag", n = 1), .SDcols = colnames(series.l)]
      
      #Calculate the seriesed and error series
      errors = t(as.matrix(series)) - sp$Ft %*% t(as.matrix(series.l))
      trend = c(sp$Ht[grepl("Tt|Mt", colnames(sp$Ht))] %*% B_TT[which(grepl("Tt|Mt", rownames(B_TT))), ])
      trend_error = c(sp$Ht[grepl("Tt|Mt", colnames(sp$Ht))] %*% errors[grepl("Tt|Mt", rownames(errors)), ])
      if(any(grepl("Ct", rownames(errors)))){
        cycle = c(sp$Ht[grepl("Ct", colnames(sp$Ht))] %*% B_TT[which(grepl("Ct", rownames(B_TT))), ])
        cycle_error = c(sp$Ht[grepl("Ct", colnames(sp$Ht))] %*% errors[grepl("Ct", rownames(errors)), ])
      }else{
        cycle = rep(0, nrow(series))
        cycle_error = rep(0, ncol(errors))
      }
      if(any(grepl("St", rownames(errors)))){
        seasonal = seasonal = c(sp$Ht[grepl("St", colnames(sp$Ht))] %*% B_TT[which(grepl("St", rownames(B_TT))), ])
        seasonal_error = c(sp$Ht[grepl("St", colnames(sp$Ht))] %*% errors[grepl("St", rownames(errors)), ])
      }else{
        seasonal = rep(0, nrow(series))
        seasonal_error = rep(0, ncol(errors))
      }
      observation_error = y - c(sp$Ht %*% t(as.matrix(series)))
      total_error = trend_error + cycle_error + seasonal_error + observation_error
      
      toret = data.table::data.table(y = y, trend = trend, cycle = cycle, seasonal = seasonal,
                                     trend_error = trend_error, cycle_error = cycle_error, seasonal_error = seasonal_error,
                                     observation_error = observation_error, total_error = total_error)
      toret[, "trend_pred" := trend - trend_error]
      toret[, "cycle_pred" := cycle - cycle_error]
      toret[, "seasonal_pred" := seasonal - seasonal_error]
      
      toret[, "seasonal_adjusted" := y - seasonal]
      toret[, "cycle_adjusted" := y - cycle]
      toret[, "seasonal_cycle_adjusted" := y - seasonal - cycle]
      toret[, "seasonal_cycle" := cycle + seasonal]
      return(toret)
    }
    names(preret) = iter
    return(preret)
  }
  names(ret) = miter
  final = rbind(data.table(method = "filter", eval(parse(text = paste(paste0("weights['", miter, "'] * ret[['", miter, "']]$filter"), collapse = " + ")))[, "date" := dates]), 
                data.table(method = "smooth", eval(parse(text = paste(paste0("weights['", miter, "'] * ret[['", miter, "']]$smooth"), collapse = " + ")))[, "date" := dates]),
                use.names = T, fill = T)
  
  if(plot == T){
    for(i in c("filter", "smooth")){
      g1 = ggplot2::ggplot(data.table::melt(final[method == i, ], id.vars = "date", measure.vars = c("y", "trend"))) +
        ggplot2::ggtitle("Actual vs Trend") +
        ggplot2::geom_line(ggplot2::aes(x = date, y = value, group = variable, color = variable)) +
        ggplot2::theme_minimal() + ggplot2::guides(color = ggplot2::guide_legend(title = NULL)) + ggplot2::theme(legend.position = "bottom")
      title = c(ifelse(grepl("cycle", model$decomp), "Cycle", ""), ifelse(grepl("seasonal", model$decomp), "Seasonal", ""), "Observation Error")
      title = title[title != ""]
      g2 = ggplot2::ggplot(data.table::melt(final[method == i, ], id.vars = "date", measure.vars = c("cycle", "seasonal", "observation_error"))) +
        ggplot2::ggtitle(paste(title, collapse = ", ")) +
        ggplot2::geom_line(ggplot2::aes(x = date, y = value, group = variable, color = variable)) +
        ggplot2::theme_minimal() + ggplot2::guides(color = ggplot2::guide_legend(title = NULL)) + ggplot2::theme(legend.position = "bottom")
      if(!all(is.na(model[["periodogram"]]))){
        g3 = ggplot2::ggplot(model[["periodogram"]][, .(`Significant Frequencies` = significant, spec = mean(spec)), by = round(period)]) +
          ggplot2::ggtitle("Periodogram") + ggplot2::scale_x_continuous(name = "Period") + ggplot2::scale_y_continuous(name = "Spec") +
          ggplot2::geom_col(ggplot2::aes(x = round, y = spec, fill = `Significant Frequencies`, color = `Significant Frequencies`), position = "dodge", width = 0.75) +
          ggplot2::theme_minimal() + ggplot2::theme(legend.position = "bottom")
        gridExtra::grid.arrange(g1, g2, g3, layout_matrix = rbind(c(1, 2), c(3, 3)),
                     top = grid::textGrob(i, gp = grid::gpar(fontsize = 20, font = 3)))
      }else{
        gridExtra::grid.arrange(g1, g2, nrow = 1, top = grid::textGrob(i, gp = grid::gpar(fontsize = 20, font = 3)))
      }
      Sys.sleep(0.1)
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

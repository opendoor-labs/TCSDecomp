#' State space model
#' Creates a state space model in list form
#' yt = Ht %*% Bt + e_t
#' Bt = Ft %*% B_{t=1} + u_t
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
#' tcs_ssm(par, y, freq, decomp, trend_spec, init)
#' @author Alex Hubbard (hubbard.alex@gmail.com)
#' @export
tcs_ssm = function(par = NULL, yt = NULL, freq = NULL, decomp = NULL, trend_spec = NULL, 
                   init = NULL, model = NULL, full_seas_freq = F){
  if(!is.null(model)){
    par = unlist(model$table[, grepl("coef_", colnames(mod$table)), with = F])
    names(par) = gsub("coef_", "", names(par))
    yt = model$data
    freq = model$table$freq
    trend_spec = model$table$model
    decomp = model$table$decomp
    full_seas_freq = model$table$full_seas_freq
  }
  
  #Define the seasonal frequencies
  if(full_seas_freq == T | freq %in% c(1, 4, 12)){
    seas_freqs = 1:(floor(freq/2))
  }else if(freq == 52){
    seas_freqs = c(1, 2, 3, 4, 8, 12, 16, 21, 26)
  }else if(freq == 365){
    seas_freqs = c(1, 2, 3, 4, 5, 6, 7, 14, 21, 30, 60, 90, 120, 150, 182)
  }
  
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
  
  #Get the parameters on exogenous data
  if(any(grepl("beta_", names(par)))){
    beta = matrix(par[grepl("beta_", names(par))], nrow = 1, ncol = length(par[grepl("beta_", names(par))]))
    colnames(beta) = gsub("beta_\\.", "", names(par[grepl("beta_", names(par))]))
  }else{
    beta = matrix(0, nrow = 1, ncol = 1)
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
    Qm = as.matrix(Matrix::bdiag(Qm, diag(2)*par["sig_c"]^2))
    colnames(Qm) = rownames(Qm) = rownames(Fm)
  }
  
  #Define the seasonal component
  if(grepl("seasonal", decomp)) {
    jiter = as.numeric(unique(gsub("[[:alpha:]]|[[:punct:]]", "", names(par)[grepl("sig_s", names(par))])))
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
    Qm = as.matrix(Matrix::bdiag(Qm, diag(unlist(lapply(jiter, function(j){rep(par[paste0("sig_s", j)], 2)})))))
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
  }else{
    P0 = init[["P0"]]
  }
  return(list(B0 = B0, P0 = P0, At = Am, Dt = Dm, Ht = Hm, Ft = Fm, Rt = Rm, Qt = Qm, beta = beta))
}

#' Detect cycle and seasonality from the data
#'
#' @param y Univariate time series of data values. May also be a 2 column data frame containing a date column.
#' @param freq Seasonality of the data (1, 4, 12, 52, 365)
#' @return List giving the decomposition and wavelet
#' @param method Method for wavelet analysis comparison ("white.nose", "shuffle", "Fourier.rand", "AR", "ARIMA"). Default is "ARIMA".
#' @param n.sim Number of simulations for wavelet analysis. Default is 100
#' @export
tcs_detect_decomp = function(y, freq, level = 0.01, method = "ARIMA", n.sim = 100){
  #Set the baseline decomposition
  decomp = "trend"
  
  #Test for seasonality and long-term cycle
  capture.output(wave <- WaveletComp::analyze.wavelet(data.frame(y), method = "ARIMA", n.sim = 500, verbose = F))
  wave = data.table::data.table(period = wave$Period, power = wave$Power.avg, pval = wave$Power.avg.pval)
  wave[, "period" := round(period)]
  wave = wave[, .(power = mean(power, na.rm = T), pval = mean(pval, na.rm = T)), by = "period"]
  wave[, "slope" := (shift(power, type = "lead", n = 1) - shift(power, type = "lag"))/(shift(period, type = "lead", n = 1) - shift(period, type = "lag", n = 1))]
  periods = wave[pval <= level, ][(shift(power, type = "lag", n = 1) <= power & shift(power, type = "lead", n = 1) <= power) & 
                                    (shift(sign(slope), type = "lag", n = 1) == 1 & shift(sign(slope), type = "lead", n = 1) == -1), ]$period
  wave[, "Significant Frequencies" := ifelse(period %in% periods, T, F), ]
  
  if(length(periods) > 0){
    #Check for seasonality
    if(freq > 1){
      if(any(abs(periods - freq) < freq/2) | seastests::wo(ts(y, frequency = freq))$stat == T){
        decomp = paste0(decomp, "-seasonal")
      }
    }
    
    #Check for longer term cycle
    if(any(periods > freq * 3)){
      decomp = paste0(decomp, "-cycle")
    }
  }
  
  #If no seasonality or cycle detected, include noise component only
  if(decomp == "trend"){
    decomp = "trend-noise"
  }
  
  return(list(wave = wave, decomp = decomp))
}

#' Detect frequency and dates from the data
#'
#' @param y Univariate time series of data values. May also be a 2 column data frame containing a date column.
#' @param freq Initial setting for the frequency detection
#' @return List giving the dates and frequency of the data
#' @export
tcs_detect_freq = function(y, freq = NULL){
  if(is.ts(y)){
    if(ifelse(is.null(ncol(y)), F, ncol(y) > 1)){
      stop("Data must be a univariate time series.")
    }
    dates = time(y)
    freq = frequency(y)
  }else{
    y = data.table::as.data.table(y)
    datecol = unlist(y[, lapply(.SD, class), .SDcols = colnames(y)])
    datecol = names(datecol[datecol %in% c("Date", "yearmon")])
    if(class(y[, c(datecol), with = F][[1]]) == "yearmon"){
      y[, c(datecol) := as.Date(eval(parse(text = datecol)))]
      is.yearmon = T
    }else{
      is.yearmon = F
    }
    if(length(datecol) == 1) {
      datediffs = unique(diff(unlist(y[, c(datecol), with = F])))
      freq = datediffs[which.max(tabulate(match(diff(y[, c(datecol), with = F][[1]]), datediffs)))]
      freq = c(365, 52, 12, 4, 1)[which.min(abs(freq -  c(1, 7, 30, 90, 365)))]
      dates = y[, c(datecol), with = F][[1]]
      if(is.yearmon == T){
        dates = as.yearmon(dates)
      }
      y = unlist(y[, colnames(y)[colnames(y) != datecol], with = F])
      rm(datediffs, datecol)
    }else if(length(datecol) > 1){
      stop("Too many date columns. Include only 1 date column or set the frequency manually.")
    }else if(length(datecol) == 0 & is.null(freq)){
      stop("No date column detected. Include a date column or set the frequency.")
    }else{
      dates = 1:length(y)
    }
  }
  return(list(data = y, dates = dates, freq = freq))
}

#' Build the date sequence as a Date type
#' 
#' @param y an object created from tcs_detect_freq
#' @author Alex Hubbard (hubbard.alex@gmail.com)
#' @export
tcs_build_dates = function(y){
  dates = tryCatch(as.Date(y$dates), 
                   error = function(err){
                     `%fun%` = lubridate::`%m+%`
                     diff = mean(unique(diff(y$dates)))
                     years = as.numeric(floor(y$dates))
                     parts = as.numeric(y$dates - years)
                     if(y$freq == 1){
                       as.Date(paste0(years, "-01-01")) %fun% years(ceiling((parts/diff)))
                     }else if(y$freq == 4){
                       as.Date(paste0(years, "-01-01")) %fun% months(ceiling((parts/diff))*3)
                     }else if(y$freq == 12){
                       as.Date(paste0(years, "-01-01")) %fun% months(ceiling((parts/diff)))
                     }else if(y$freq == 52){
                       as.Date(paste0(years, "-01-01")) %fun% weeks(ceiling((parts/diff)))
                     }else if(y$freq == 365){
                       as.Date(paste0(years, "-01-01")) %fun% dates(ceiling((parts/diff)))
                     }
                   })
  return(dates)
}

#' Trend cycle seasonal decomposition using the Kalman filter.
#' 
#' Estimates a structural time series model using the Kalman filter and maximum likelihood.
#' The seasonal and cycle components are assumed to be of a trigonometric form.
#' The function checks three trend specifications to decompose a univariate time series
#' into trend, cycle, and/or seasonal components plus noise. The function automatically
#' detects the frequency and checks for a seasonal and cycle component if the user does not specify
#' the frequency or decomposition model. This can be turned off by setting freq or specifying decomp.
#' State space model for decomposition follows
#' Yt = T_t + C_t + S_t + BX_t + e_t, e_t ~ N(0, sig_e^2)
#' Y is the data
#' T is the trend component
#' C is the cycle component
#' S is the seasonal component
#' X is the exogenous data with parameter vector B
#' e is the observation error
#' @param y Univariate time series of data values. May also be a 2 column data frame containing a date column.
#' @param exo Matrix of exogenous variables. Can be used to specify regression effects or other seasonal effects like holidays, etc.
#' @param freq Seasonality of the data (1, 4, 12, 52, 365)
#' #' @param full_seas_freq Whether to use the full range of seasonal frequencies or the default truncated case for higher frequency data
#' @param decomp Decomposition model ("tend-cycle-seasonal", "trend-seasonal", "trend-cycle", "trend-noise")
#' @param trend_spec Trend specification ("rw", "rwd", "2rw"). The default is NULL which will choose the best of all specifications based on the maximum likielhood. 
#' "rw" is the random walk trend. "rwd" is the random walk with random walk drift trend. "2rw" is a 2nd order random walk trend.
#' @param multiplicative If data should be logged to create a multiplicative model
#' @param optim_methods Vector of 1 to 3 optimization methods in order of preference ("NR", "BFGS", "CG", "BHHH", or "SANN")
#' @param det_obs Set the observation equation error variance to 0 (deterministic observation equation)
#' @param det_trend Set the trend error variance to 0 (deterministic trend)
#' @param det_seas Set the seasonality error variances to 0 (deterministic seasonality)
#' @param det_cycle Set thecycle error variance to 0 (deterministic cycle)
#' @param det_drift Set the drift error variance to 0 (deterministic drift)
#' @param maxit Maximum number of iterations for the optimization
#' @param wavelet.method Method for wavelet analysis comparison ("white.nose", "shuffle", "Fourier.rand", "AR", "ARIMA"). Default is "ARIMA".
#' @param wavelet.sim Number of simulations for wavelet analysis. Default is 100
#' @return List of estimation values including coefficients, convergence code, datea frequency, decomposition used, and trend specification selected.
#' @examples
#' tcs_decomp_estim(y = DT[, c("date", "y")])
#' If multiplicative = T, then the data is logged and the original model becomes multiplicative (Y_t = T_t * C_t * S_t * BX_t * e_t)
#' If trend is "rw", the trend model is T_t = T_{t-1} + e_t, e_t ~ N(0, sig_t^2)
#' If trend is "rwd", the trend model is T_t = M_{t-1} + T_{t-1} + e_t, e_t ~ N(0, sig_t^2) with 
#'      M_t = M_{t-1} + n_t, n_t ~ N(0, sig_m^2) where
#' If trend is "2rw", the trend model is T_t = 2T_{t-1} - T_{t-2} + e_t, e_t ~ N(0, sig_t^2)
#' If det_obs = T then sig_e is set to 0 
#' If det_trend = T then the variacne of the trend (sig_t) is set to 0 and is referred to as a smooth trend
#' If det_drift = T then the variance of the drift (sig_m) is set to 0 and is refereed to as a deterministic drift
#' If det_seas = T then the variance all seasonality frequencies j (sig_s) are set to 0 and is referred to as deterministic seasonality
#' If det_cycle = T then the variance of the cyclce (sig_c) is set to 0 and is refreed to as a deterministic cycle
#' @author Alex Hubbard (hubbard.alex@gmail.com)
#' @export
tcs_decomp_estim = function (y, exo = NULL, freq = NULL, full_seas_freq = F, decomp = NULL, trend_spec = NULL,
                             multiplicative = F,
                             det_obs = F, det_trend = F, det_seas = F, det_cycle = F, det_drift = F,
                             wavelet.method = "ARIMA", wavelet.sim = 100, 
                             level = 0.01, optim_methods = c("BFGS", "CG", "NM"), maxit = 10000){
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
  
  #Get the frequency of the data
  y = tcs_detect_freq(y, freq)
  dates = tcs_build_dates(y)
  freq = y$freq
  y = y$data
  
  #Remove leading and trailing NAs
  range = which(!is.na(y))
  y = unname(y[range[1]:range[length(range)]])
  dates = dates[range[1]:range[length(range)]]
  if(!is.null(exo)){
    exo = exo[range[1]:range[length(range)], ]
  }
  
  #Log data if multiplicative is T
  if(multiplicative == T){
    if(!all(y > 0)){
      warning("A multiplicative model was specified but the data is not all positive. Reverting to linear model.")
      multiplicative = F
    }else{
      y = log(y)
    }
  }
  
  #Set the decomposition
  if(is.null(decomp)){
    decomp = tcs_detect_decomp(y, freq, level, wavelet.method, wavelet.sim)
    wave = decomp$wave
    decomp = decomp$decomp
  }else{
    wave = NULL
  }
  
  #Standardize
  mean = mean(y, na.rm = T)
  sd = sd(y, na.rm = T)
  y = (y - mean)/sd
  
  #Define the trend specifications to estiamte
  if(is.null(trend_spec)) {
    iter = c("rw", "rwd", "2rw")
  }else {
    if(!trend_spec %in% c("rw", "rwd", "2rw")){
      stop("trend_spec must be 'rw', 'rwd', or '2rw'.")
    }
    iter = trend_spec
  }
  
  cl = parallel::makeCluster(min(c(parallel::detectCores(), length(iter))))
  doSNOW::registerDoSNOW(cl)
  invisible(snow::clusterCall(cl, function(x) .libPaths(x), .libPaths()))
  `%fun%` = foreach::`%dopar%`
  fit = foreach::foreach(i = iter, .combine = function(DT1, DT2){rbind(DT1, DT2, use.names = T, fill = T)},
                         .packages = c("data.table", "Matrix", "maxLik", "imputeTS"), .export = c("tcs_ssm", "kalman_filter", "kalman_smoother")) %fun% {
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
      if(full_seas_freq == T | freq %in% c(1, 4, 12)){
        seas_freqs = 1:(floor(freq)/2)
      }else if(freq == 52){
        seas_freqs = c(1, 2, 3, 4, 8, 12, 16, 20, 24)
      }else if(freq == 365){
        seas_freqs = c(1, 2, 3, 4, 5, 6, 7, 14, 21, 30, 60, 90, 120, 150, 182)
      }
      par = c(par, sig_s = unname(rep(par["sig_t"]/(2 * length(seas_freqs)), length(seas_freqs))))
      names(par)[grepl("sig_s", names(par))] = paste0("sig_s", seas_freqs)
    }
    if(grepl("cycle", decomp)){
      if(!is.null(wave)){
        period = mean(wave[`Significant Frequencies` == T, ][2:.N]$period)
      }else{
        period = 5*freq
      }
      par = c(par, lambda = log((2 * pi/(period))/(pi - 2 * pi/(period))), rho = log((0.5)/(1 - 0.5)), sig_c = unname(par["sig_t"]/2))
    }
    if(grepl("seasonal|cycle", decomp)){
      par = c(par, sig_e = unname(par["sig_t"]))
    }
    
    #Set any fixed parameters
    fixed = NULL
    if(det_obs == T){
      par["sig_e"] = 0
      fixed = c(fixed, "sig_e")
    }
    if(det_trend == T){
      par["sig_t"] = 0
      fixed = c(fixed, "sig_t")
    }
    if(det_seas == T){
      par[grepl("sig_s", names(par))] = 0
      fixed = c(fixed, names(par)[grepl("sig_s", names(par))])
    }
    if(det_cycle == T){
      par["sig_c"] = 0
      fixed = c(fixed, "sig_c")
    }
    if(det_drift == T){
      par["sig_m"] = 0
      fixed = c(fixed, "sig_m")
    }
    if(is.null(exo)){
      X = t(matrix(0, nrow = length(y), ncol = 1))
      rownames(X) = "X"
      par = c(par, beta_X = 0)
      fixed = c(fixed, "beta_X")
    }else{
      X = t(exo)
      par = c(par, beta_ = coef(lm(y ~ . - 1, data = data.frame(cbind(y, exo)))))
      par[grepl("beta_", names(par))] = 0
    }
    
    #Define the objective function
    objective = function(par, na_locs, freq, decomp, trend_spec, init = NULL){
      yt = matrix(get("y"), nrow = 1)
      sp = tcs_ssm(par = par, yt = yt, freq = freq, decomp = decomp, 
                   trend_spec = trend_spec, init = init, full_seas_freq = full_seas_freq)
      ans = kalman_filter(B0 = matrix(sp$B0, ncol = 1), P0 = sp$P0, Dt = sp$Dt, At = sp$At, Ft = sp$Ft, Ht = sp$Ht, Qt = sp$Qt, Rt = sp$Rt, 
                          yt = yt, X = X, beta = sp$beta)
      if (!is.null(na_locs)) {
        fc = sp$Ht %*% ans$B_tt
        yt[, na_locs] = fc[, na_locs]
        ans = kalman_filter(matrix(sp$B0, ncol = 1), P0 = sp$P0, Dt = sp$Dt, At = sp$At, Ft = sp$Ft, Ht = sp$Ht, Qt = sp$Qt, Rt = sp$Rt, 
                            yt = yt, X = X, beta = sp$beta)
        assign("y", c(yt), .GlobalEnv)
      }
      return(ans$loglik)
    }
    
    #Get initial values for Kalman Filter
    #Get initial values for Kalman Filter
    sp = tcs_ssm(par = par, yt = y, freq = freq, decomp = decomp, trend_spec = i, full_seas_freq = full_seas_freq)
    init = kalman_filter(B0 = matrix(sp$B0, ncol = 1), P0 = sp$P0, Dt = sp$Dt, At = sp$At, Ft = sp$Ft, Ht = sp$Ht, Qt = sp$Qt, Rt = sp$Rt, 
                         yt = matrix(y, nrow = 1), X = X, beta = sp$beta)
    init2 = tryCatch(kalman_smoother(B_tl = init[["B_tl"]], B_tt = init[["B_tt"]], P_tl = init[["P_tl"]], P_tt = init[["P_tt"]], Ft = sp$Ft), 
                     error = function(err){NULL})
    if(is.null(init)){
      init = list(B0 = sp$B0, P0 = init$P0) 
    }else{
      init = list(B0 = init[["B_tt"]][, 1], P0 = init[["P_tt"]][, , 1])
    }
    
    #Find anhy na values
    if(any(is.na(y))){
      na_locs = which(is.na(y))
    }else{
      na_locs = NULL
    }
    
    #Estimate the model
    out = tryCatch(maxLik::maxLik(logLik = objective, 
                                  start = par, method = optim_methods[1], fixed = fixed, 
                                  finalHessian = F, hess = NULL, control = list(printLevel = 2, iterlim = maxit), init = init, na_locs = na_locs, 
                                  freq = freq, decomp = decomp, trend_spec = i), 
                   error = function(err){
                     tryCatch(maxLik::maxLik(logLik = objective, 
                                             start = par, method = optim_methods[min(c(2, length(optim_methods)))], fixed = fixed, 
                                             finalHessian = F, hess = NULL, control = list(printLevel = 2, iterlim = maxit), init = init, na_locs = na_locs, 
                                             freq = freq, decomp = decomp, trend_spec = i),
                              error = function(err){
                                tryCatch(maxLik::maxLik(logLik = objective, start = par, method = optim_methods[min(c(3, length(optim_methods)))], fixed = fixed, 
                                                        finalHessian = F, hess = NULL, control = list(printLevel = 2, iterlim = maxit), init = init, na_locs = na_locs, 
                                                        freq = freq, decomp = decomp, trend_spec = i),
                                         error = function(err){NULL})
                              })
                   })
    rm(init, init2)
    gc()
    
    if(!is.null(out)){
      #Retreive the model output
      ret = data.table::data.table(model = i, freq = freq, full_seas_freq = full_seas_freq, decomp = decomp, multiplicative = multiplicative, 
                                   convergence = out$code, loglik = out$maximum,
                                   matrix(coef(out), nrow = 1, dimnames = list(NULL, paste0("coef_", names(coef(out))))))
      if(is.null(exo)){
        ret[, colnames(ret)[grepl("beta_", colnames(ret))] := NULL]
      }
      return(ret)
    }else{
      return(NULL)
    }
  }
  snow::stopCluster(cl)
  
  #Select the best model based on the maximum likelihood
  fit = list(table = fit[loglik == max(loglik, na.rm = T), ], data = y*sd + mean, dates = dates)
  if(multiplicative == T){
    fit$data = exp(fit$data)
  }
  if(exists("wave")){
    fit[["wavelet"]] = wave
  }else{
    fit[["wavelet"]] = NA
  }
  if(is.null(exo)){
    fit[["exo"]] = NA
  }else{
    fit[["exo"]] = exo
  }
  return(fit)
}

#' Kalman Filter
#' 
#' Kalman filter an estimated model from tcs_decomp_estim object
#' @param y Univariate time series of data values. May also be a 2 column data frame containing a date column.
#' @param model Structural time series model estimated using tcs_decomp_estim.
#' @param plot Logial, whether to plot the output or not.
#' @return List of data tables containing the filtered and smoothed series.
#' @examples
#' tcs_decomp_filter(y = DT[, c("date", "y")], model = model)
#' @author Alex Hubbard (hubbard.alex@gmail.com)
#' @export
tcs_decomp_filter = function(model, y = NULL, exo = NULL, plot = F){
  
  #Get the dates and frequency of the data
  if(is.null(y)){
    dates = model$dates
    freq = model$table$freq
    y = model$data
    X = t(as.matrix(model$exo))
    if(all(is.na(X))){
      X = t(matrix(0, nrow = length(y), ncol = 1))
      rownames(X) = "X"
    }
  }else{
    y = tcs_detect_freq(y, model$freq)
    dates = tcs_build_dates(y)
    freq = y$freq
    y = y$data
    if(is.null(exo)){
      X = t(matrix(0, nrow = length(y), ncol = 1))
      rownames(X) = "X"
    }
  }
  
  #Remove leading and trailing NAs
  range = which(!is.na(y))
  y = unname(y[range[1]:range[length(range)]])
  dates = dates[range[1]:range[length(range)]]
  
  #Apply multiplicative model
  if(model$table$multiplicative == T){
    y = log(y)
  }
  
  #Standardize
  y_mean = mean(y, na.rm = T)
  y_sd = sd(y, na.rm = T)
  y = (y - y_mean)/y_sd
  
  #Find any missing values
  if(any(is.na(y))){
    na_locs = which(is.na(y))
  }else{
    na_locs = NULL
  }
  
  #Get model specifications
  if(freq != model$table$freq){
    warning("Frequency of data does not match that of the model. Using model frequency instead of the detected data frequency.")
    freq = model$table$freq
  }
  decomp = model$table$decomp
  trend_spec = model$table$model
  full_seas_freq = model$table$full_seas_freq
  
  #Get the coefficients
  cs = unlist(model$table[, grepl("coef_", colnames(model$table)), with = F])
  cs = cs[!is.na(cs)]
  names(cs) = gsub("coef_", "", names(cs))
  
  #Filter the data
  sp = tcs_ssm(par = cs, yt = y, freq = freq, decomp = decomp, trend_spec = trend_spec, full_seas_freq = full_seas_freq)
  init = kalman_filter(B0 = matrix(sp$B0, ncol = 1), P0 = sp$P0, Dt = sp$Dt, At = sp$At, Ft = sp$Ft, Ht = sp$Ht, Qt = sp$Qt, Rt = sp$Rt, 
                       yt = matrix(y, nrow = 1), X = X, beta = sp$beta)
  init = kalman_smoother(B_tl = init$B_tl, B_tt = init$B_tt, P_tl = init$P_tl, P_tt = init$P_tt, Ft = sp$Ft)
  init = list(B0 = init[["B_tt"]][, 1], P0 = init[["P_tt"]][, , 1])
  sp = tcs_ssm(par = cs, yt = y, freq = freq, decomp = decomp, trend_spec = trend_spec, init = init, full_seas_freq = full_seas_freq)
  ans = kalman_filter(B0 = matrix(sp$B0, ncol = 1), P0 = sp$P0, Dt = sp$Dt, At = sp$At, Ft = sp$Ft, Ht = sp$Ht, Qt = sp$Qt, Rt = sp$Rt, 
                      yt = matrix(y, nrow = 1), X = X, beta = sp$beta)
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
    
    if(!is.null(exo) | !all(is.na(model$exo))){
      XB = do.call("cbind", lapply(names(cs)[grepl("beta_", names(cs))], function(x){
        cs[x]*X[gsub("beta_\\.", "", x), ]
      }))
      colnames(XB) = rownames(X)
      toret = cbind(toret, XB)
    }
    return(toret)
  }
  names(preret) = iter
  
  #Combine the filtered and smoothed series
  final = rbind(data.table(method = "filter", date = dates, preret[["filter"]]), 
                data.table(method = "smooth", date = dates, preret[["smooth"]]), 
                use.names = T, fill = T)[, "date" := as.Date(date)]

  cols = c("y", "trend", "trend_pred", "seasonal_adjusted", "cycle_adjusted", "seasonal_cycle_adjusted")
  final[, c(cols) := lapply(.SD, function(x){x*y_sd + y_mean}), .SDcols = c(cols), by = "method"]
  cols = c("seasonal", "seasonal_pred", "cycle", "cycle_pred", "seasonal_cycle",
           "trend_error", "cycle_error", "seasonal_error", "observation_error", "total_error", rownames(X))
  final[,  c(cols) := lapply(.SD, function(x){x*y_sd}), .SDcols = c(cols), by = "method"]
  
  if(model$table$multiplicative == T){
    final[, c(colnames(final)[!colnames(final) %in% c("method", "date")]) := lapply(.SD, exp),
          .SDcols = c(colnames(final)[!colnames(final) %in% c("method", "date")])]
  }
  
  if(plot == T) {
    for(i in c("filter", "smooth")){
      g1 = ggplot2::ggplot(data.table::melt(final[method == i, ], id.vars = "date", measure.vars = c("y", "trend"))) + 
        ggplot2::ggtitle("Actual vs Trend") +
        ggplot2::scale_x_date(name = "Date") +
        ggplot2::scale_y_continuous(name = "Value") +
        ggplot2::geom_line(ggplot2::aes(x = date, y = value, group = variable, color = variable)) + ggplot2::theme_minimal() + 
        ggplot2::guides(color = ggplot2::guide_legend(title = NULL)) + 
        ggplot2::theme(legend.position = "bottom")
      title = c(ifelse(grepl("cycle", decomp), "Cycle", ""), ifelse(grepl("seasonal", decomp), "Seasonal", ""), "Observation Error")
      title = title[title != ""]
      g2 = ggplot2::ggplot(data.table::melt(final[method == i, ], id.vars = "date", measure.vars = c("cycle", "seasonal", "observation_error"))) + 
        ggplot2::ggtitle(paste(title, collapse = ", ")) + 
        ggplot2::scale_x_date(name = "Date") +
        ggplot2::scale_y_continuous(name = "Value") +
        ggplot2::geom_line(ggplot2::aes(x = date, y = value, group = variable, color = variable)) + 
        ggplot2::theme_minimal() + ggplot2::guides(color = ggplot2::guide_legend(title = NULL)) + 
        ggplot2::theme(legend.position = "bottom")
      if(!all(is.na(model[["wavelet"]]))){
        g3 = ggplot2::ggplot(model[["wavelet"]]) + 
          ggplot2::ggtitle("Wavelet Power") + ggplot2::scale_x_continuous(name = "Period") + ggplot2::scale_y_continuous(name = "Power") + 
          ggplot2::geom_col(ggplot2::aes(x = period, y = power, fill = `Significant Frequencies`, color = `Significant Frequencies`), 
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
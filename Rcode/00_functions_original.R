###################################
####### Data functions ##########
###################################
# Functions- from (Lindenmeyer and Torrent, 2023)
preparacao = function(X, i) {
  
  if (tipo[i] == 0) {
    
    return(log(X))
  }
  
  if (tipo[i] == 2) {
    
    return(X)
    
  }
  
  if (tipo[i] == 1) {
    
    return(X)
  }
  
  if (tipo[i] == 3) {
    
    return(X)
  }
}

cresc_discreto = function(X) {
  
  Y = c()
  
  for (i in 2:length(X)) {
    
    y = X[i]/X[i-1]-1
    
    Y = append(Y, y)
    
  }
  
  return(Y)
}

# Do the selected IPEA transform to yt
transform_singlestep <- function (yt, transform_id) {
  
  ytlen=length(yt)
  transformed=mat.or.vec(ytlen,1)
  transformed[]=NA
  
  if(transform_id==1) {
    transformed=yt
  }
  if(transform_id==2) {
    transformed[2:ytlen]=yt[2:ytlen]-yt[1:(ytlen-1)]
  }
  if(transform_id==3) {
    transformed[3:ytlen]=(yt[3:ytlen]-yt[2:(ytlen-1)])-(yt[2:(ytlen-1)]-yt[1:(ytlen-2)])
  }
  if(transform_id==4) {
    transformed=log(yt)
  }
  if(transform_id==5) {
    transformed[2:ytlen]=log(yt[2:ytlen])-log(yt[1:(ytlen-1)])
  }
  if(transform_id==6) {
    transformed[3:ytlen]=(log(yt[3:ytlen])-log(yt[2:(ytlen-1)]))-(log(yt[2:(ytlen-1)])-log(yt[1:(ytlen-2)]))
  }
  if(transform_id==7) {
    transformed[2:ytlen]=yt[2:ytlen]/yt[1:(ytlen-1)]-1
  }
  return(transformed)
}

# Data prep function- adapted from Medeiros et al., 2021)
dataprep = function(df,horizon,variable,dataset,lag=12,K=10){
  
  y=df[,variable]
  if(dataset=="B0_AR"){#Data poor
    x=as.matrix(df[,variable])
  }
  if(dataset=="B0_ARDI"){#Data rich
    factors=princomp(scale(df))$scores[,1:K] 
    x = cbind(df[,variable],factors)
  }
  if(dataset=="B1"){
    x = df
  }
  if(dataset=="B2"){
    factors=princomp(scale(df))$scores 
    x = as.matrix(factors)
  }
  if(dataset=="B3"){
    factors=princomp(scale(df))$scores
    x = cbind(df[,variable],factors)
  }
  
  X=embed(as.matrix(x),lag)

  Xin=X[-c((nrow(X)-horizon+1):nrow(X)),,drop=FALSE]
  Xout=X[nrow(X),]
  Xout=t(as.vector(Xout))
  yin=tail(y,nrow(Xin))

  return(list(Xin = Xin, Xout = Xout, yin = yin))

}

###################################
####### DATA-POOR MODELS ##########
###################################

linear_kfold <- function(Xin,yin,horizon,k,lag){
  
  max_lag <- 12
  # Set number of predictions (n-k-h+1)
  n_size <- length(yin)-max_lag-horizon+1
  indices <- sample(1:n_size)  # Shuffle the indices without replacement
  m <-floor(n_size/k) 
  MSE <- rep(0,k)
  kfold_pred <- matrix(NA, nrow =n_size, ncol = (1:lag))
  
  for (j in 1:k){
    start_index <- (j - 1) * m + 1
    end_index <- min(j * m, n_size)
    ind <- indices[start_index:end_index]
    xsub <- rep(TRUE, length(yin))
    xsub[ind] <- FALSE
    
    # Set kfold_test and kfold_train
    kfold_test <- as.matrix(yin[ind])
    sub <- complete.cases(Xin, yin)
    sub <- sub & xsub
    
    kfold_fit = lm(yin[sub] ~ Xin[sub,1:lag,drop=FALSE])
    kfold_coef= coef(kfold_fit)
    kfold_coef[is.na(kfold_coef)] = 0
    
    # Predict results for the whole in-sample set
    for (i in 1:n_size){
      kfold_Xout= Xin[i,]
      kfold_Xout= t(as.vector(kfold_Xout))
      kfold_pred[i] <-  c(1,kfold_Xout[,1:lag])%*%kfold_coef
    }
    kfold_pred = c(rep(NA_real_, max_lag),kfold_pred)
    
    # Choose hyperparameters based on MSE results 
    pred_error <- (kfold_test - kfold_pred[ind])
    MSE[j] <- mean(pred_error^2, na.rm = TRUE)
  }
  
  MSE <- mean(MSE,na.rm = TRUE)
  
  return(MSE)
  
}

func_ar = function (df,horizon,variable,lag_orders,reoptimize_hyperparameters=FALSE,type,best, ...) {
  prep_data = dataprep(df,horizon,variable,dataset="B0_AR")
  Xin = prep_data$Xin
  yin = prep_data$yin
  Xout = prep_data$Xout
  if (reoptimize_hyperparameters==TRUE){
    # Initialize variables to store best model and criteria value
    best = list()
    best_parameter = Inf
  
    for (lag in lag_orders) {
      fit <- lm(yin~Xin[,1:lag])
      if(type=="bic"){
        bic_sel = BIC(fit)
        # Update best model if criteria value is lower
        if (bic_sel < best_parameter) {
          best_parameter <- bic_sel
          best_fit = fit
          best$lag <- lag
        }
      }
      if(type=="aic"){
        aic_sel = AIC(fit)
        # Update best model if criteria value is lower
        if (aic_sel < best_parameter) {
          best_parameter <- aic_sel
          best_fit = fit
          best$lag <- lag
        }
      }
      if(type=="cv"){
        MSE <- linear_kfold(Xin,yin,horizon,k=5,lag)
        if (MSE < best_parameter) {
          best_parameter=MSE
          best_fit = fit
          best$lag <- lag
        }
      }
    }
  }
  if (reoptimize_hyperparameters==FALSE){
    best_fit <- lm(yin~Xin[,1:best$lag])
  }

  coef_opt=coef(best_fit)
  coef_opt[is.na(coef_opt)] = 0
  pred <- c(1,Xout[,1:best$lag])%*%coef_opt
  
  return(list(pred=pred, best=best))
}

func_shrink_poor = function (df,horizon,variable,lag_orders,reoptimize_hyperparameters=FALSE,type,best, ...) {
  prep_data = dataprep(df,horizon,variable,dataset="B0_AR")
  Xin = prep_data$Xin
  yin = prep_data$yin
  Xout = prep_data$Xout
  if (reoptimize_hyperparameters==TRUE){
    # Initialize variables to store best model and criteria value
    # Here loop through different lags is necessary because 
    # ridge regression does not do variable selection on its own
    # lag orders in this case is (2:12) since x should be a matrix with 2 or more columns
    best = list()
    best_parameter = Inf
    for (lag in lag_orders[2:length(lag_orders)]) {
      if(type=="cv"){
        if(regu=="ridge"){
          fit = cv.glmnet(Xin[,1:lag],yin,alpha=0,type.measure = "mse", nfolds = 5)
        }
        if(regu=="lasso"){
          fit = cv.glmnet(Xin[,1:lag],yin,alpha=1,type.measure = "mse", nfolds = 5)
        }
        if(regu=="en"){
          fit  = cv.glmnet(Xin[,1:lag],yin,alpha=0.5,type.measure = "mse", nfolds = 5)
        }
        MSE = fit$cvm[fit$lambda == fit$lambda.min]
        if (MSE < best_parameter) {
          best_parameter = MSE
          best_fit = fit
          best$lag = lag
          best$lambda = fit$cvm[fit$lambda == fit$lambda.min]
        }
      }
    }
  }
  if (reoptimize_hyperparameters==FALSE){
    if(regu=="ridge"){
      best_fit = glmnet(Xin[,1:best$lag],yin,alpha=0)
    }
    if(regu=="lasso"){
      best_fit = glmnet(Xin[,1:best$lag],yin,alpha=1)
    }
    if(regu=="en"){
      best_fit = glmnet(Xin[,1:best$lag],yin,alpha=0.5)
    }
  }
  pred=predict(best_fit,Xout[,1:best$lag],s = best$lambda)

  return(list(pred=pred,best=best))
}

func_rf_poor = function (df,horizon,variable,lag_orders,reoptimize_hyperparameters=FALSE,type,best, ...) {
  best=list()
  prep_data = dataprep(df,horizon,variable,dataset="B0_AR")
  Xin = prep_data$Xin
  yin = prep_data$yin
  Xout = prep_data$Xout
  best_fit = randomForest::randomForest(Xin,yin,importance = TRUE)
  pred=predict(best_fit,Xout)
  
  return(list(pred=pred,best=best))
}

func_bols_poor = function (df,horizon,variable,lag_orders,reoptimize_hyperparameters=FALSE,type,best, ...) {
  # Linear Boosting
  prep_data = dataprep(df,horizon,variable,dataset="B0_AR")
  Xin = prep_data$Xin
  yin = prep_data$yin
  Xout = data.frame(prep_data$Xout)
  concat = data.frame(yin,Xin)
  bctrl <- boost_control(mstop = 300)
  
  if(type=="cv"){
    # Here the variable selection occurs inside the function and from 1:12
    best_fit = mboost::gamboost(yin~.,concat,family=Gaussian(), baselearner="bols", control = bctrl)
    best = list()
    cvm <- cv(model.weights(best_fit), type = "kfold", B = 5)
    cvr <- cvrisk(best_fit, folds = cvm, papply = lapply)
    best$blstop <- mstop(cvr)

  }
  pred=predict(best_fit[best$blstop],Xout)
  
  return(list(pred=pred,best=best))
}  

func_bbs_poor = function (df,horizon,variable,lag_orders,reoptimize_hyperparameters=FALSE,type,best, ...) {
  # Boost with spline (nonlinear)
  prep_data = dataprep(df,horizon,variable,dataset="B0_AR")
  Xin = prep_data$Xin
  yin = prep_data$yin
  Xout = data.frame(prep_data$Xout)
  concat = data.frame(yin,Xin)
  bctrl <- boost_control(mstop = 300)
  df = 4
  
  if(type=="cv"){
    # Here the variable selection occurs inside the function and from 1:12
    best_fit = mboost::gamboost(yin ~ ., concat, family = Gaussian(), baselearner="bbs",dfbase=df,control = bctrl)
    best = list()
    cvm <- cv(model.weights(best_fit), type = "kfold", B = 5)
    cvr <- cvrisk(best_fit, folds = cvm, papply = lapply)
    best$blstop <- mstop(cvr)
  }
  pred=predict(best_fit[best$blstop],Xout)
  
  return(list(pred=pred,best=best))
}  

func_svr_linear_poor = function (df,horizon,variable,lag_orders,reoptimize_hyperparameters=FALSE,type,best, ...) {
  if (reoptimize_hyperparameters==TRUE){
    # Initialize variables to store best model and criteria value
    # Here loop through different lags is necessary because 
    # svr does not do variable selection on its own
    best = list()
    best_parameter = Inf
    for (lag in lag_orders) {
      prep_data = dataprep(df,horizon,variable,lag=lag,dataset="B0_AR")
      Xin = prep_data$Xin
      yin = prep_data$yin
      Xout = prep_data$Xout
      if(type=="cv"){
        obj <- tune(svm,Xin,yin,
                    ranges = list(cost = c(0.1,0.5,1,2,5),
                                  epsilon=seq(0.1, 0.5, by = 0.1)),
                    tunecontrol = tune.control(sampling = "fix")
        )
        fit = svm(Xin,yin, kernel = "linear", type= "eps-regression",cross=5, scale = FALSE,
                  cost = obj$best.parameters$cost,
                  epsilon = obj$best.parameters$epsilon)
        MSE = mean(fit$MSE)
        if (MSE < best_parameter) {
          best_parameter = MSE
          best_fit = fit
          best$lag = lag
          best$Xout = prep_data$Xout
          best$cost = fit$cost
          best$epsilon=fit$epsilon
        }
      }
    }
  }
  if (reoptimize_hyperparameters==FALSE){
    prep_data = dataprep(df,horizon,variable,lag=best$lag,dataset="B0_AR")
    Xin = prep_data$Xin
    yin = prep_data$yin
    best$Xout = prep_data$Xout
    best_fit = svm(Xin,yin, kernel = "linear", type= "eps-regression",cost=best$cost, epsilon=best$epsilon, scale = FALSE)
  }
  
  pred=predict(best_fit,best$Xout)
  
  return(list(pred=pred,best=best))
}  
  
func_svr_rbf_poor = function (df,horizon,variable,lag_orders,reoptimize_hyperparameters=FALSE,type,best, ...) {
  if (reoptimize_hyperparameters==TRUE){
    # Initialize variables to store best model and criteria value
    # Here loop through different lags is necessary because 
    # svr does not do variable selection on its own
    best = list()
    best_parameter = Inf
    for (lag in lag_orders) {
      prep_data = dataprep(df,horizon,variable,lag=lag,dataset="B0_AR")
      Xin = prep_data$Xin
      yin = prep_data$yin
      Xout = prep_data$Xout
      if(type=="cv"){
        obj <- tune(svm,Xin,yin,
                    ranges = list(gamma = c(seq(0.1, 0.9, by = 0.1),(1/ncol(Xin))), 
                                  cost = c(0.1,0.5,1,2,5),
                                  epsilon=seq(0.1, 0.5, by = 0.1)),
                    tunecontrol = tune.control(sampling = "fix")
        )
        fit = svm(Xin,yin, kernel = "radial", type= "eps-regression",cross=5, scale = FALSE,
                  gamma = obj$best.parameters$gamma,
                  cost = obj$best.parameters$cost,
                  epsilon = obj$best.parameters$epsilon)
        MSE = mean(fit$MSE)
        if (MSE < best_parameter) {
          best_parameter = MSE
          best_fit = fit
          best$lag = lag
          best$Xout = prep_data$Xout
          best$cost = fit$cost
          best$epsilon=fit$epsilon
          best$gamma=fit$gamma
        }
      }
    }
  }
  if (reoptimize_hyperparameters==FALSE){
    prep_data = dataprep(df,horizon,variable,lag=best$lag,dataset="B0_AR")
    Xin = prep_data$Xin
    yin = prep_data$yin
    best$Xout = prep_data$Xout
    best_fit = svm(Xin,yin, kernel = "radial", type= "eps-regression",scale = FALSE,
                   gamma=best$gamma,
                   cost=best$cost, 
                   epsilon=best$epsilon)
  }
  
  pred=predict(best_fit,best$Xout)
  
  return(list(pred=pred,best=best))
}    
  
###################################
####### DATA-RICH MODELS ##########
###################################

ardi_kfold <- function(Xin,yin,horizon,k,lag){
  
  # Set number of predictions (n-k-h+1)
  n_size <- length(yin)-lag-horizon+1
  indices <- sample(1:n_size)  # Shuffle the indices without replacement
  m <-floor(n_size/k) 
  MSE <- rep(0,k)
  kfold_pred <- matrix(NA, nrow =n_size, ncol = (1:lag))
  
  for (j in 1:k){
    start_index <- (j - 1) * m + 1
    end_index <- min(j * m, n_size)
    ind <- indices[start_index:end_index]
    xsub <- rep(TRUE, length(yin))
    xsub[ind] <- FALSE
    
    # Set kfold_test and kfold_train
    kfold_test <- as.matrix(yin[ind])
    sub <- complete.cases(Xin, yin)
    sub <- sub & xsub
    
    kfold_fit = lm(yin[sub] ~ Xin[sub,,drop=FALSE])
    kfold_coef= coef(kfold_fit)
    kfold_coef[is.na(kfold_coef)] = 0
    
    # Predict results for the whole in-sample set
    for (i in 1:n_size){
      kfold_Xout= Xin[i,]
      kfold_Xout= t(as.vector(kfold_Xout))
      kfold_pred[i] <-  c(1,kfold_Xout)%*%kfold_coef
    }
    kfold_pred = c(rep(NA_real_, lag),kfold_pred)
    
    # Choose hyperparameters based on MSE results 
    pred_error <- (kfold_test - kfold_pred[ind])
    MSE[j] <- mean(pred_error^2, na.rm = TRUE)
  }
  
  MSE <- mean(MSE,na.rm = TRUE)
  
  return(MSE)
  
}

func_ardi = function (df,horizon,variable,lag_orders,reoptimize_hyperparameters=FALSE,type,best, ...) {
  if (reoptimize_hyperparameters==TRUE){
  # Initialize variables to store best model and criteria value
  best = list()
  best_parameter = Inf
  
    for (K in factor_orders) {
      # Here pf and py are being estimated together!
      for (lag in lag_orders) {
        prep_data = dataprep(df,horizon,variable,lag=lag,K=K,dataset="B0_ARDI")
        Xin = prep_data$Xin
        yin = prep_data$yin
        Xout = prep_data$Xout
        fit <- lm(yin~Xin)
        
        if(type=="bic"){
          bic_sel = BIC(fit)
          # Update best model if criteria value is lower
          if (bic_sel < best_parameter) {
            best_parameter <- bic_sel
            best_fit = fit
            best$lag <- lag
            best$factor <- K
            best$Xout = prep_data$Xout
          }
        }
        if(type=="aic"){
          aic_sel = AIC(fit)
          # Update best model if criteria value is lower
          if (aic_sel < best_parameter) {
            best_parameter <- aic_sel
            best_fit = fit
            best$lag <- lag
            best$factor <- K
            best$Xout = prep_data$Xout
          }
        }
        
        if(type=="cv"){
          MSE <- ardi_kfold(Xin,yin,horizon,k=5,lag)
          if (MSE < best_parameter) {
            best_parameter=MSE
            best_fit = fit
            best$lag <- lag
            best$factor <- K
            best$Xout = prep_data$Xout
          }
        }
      }
    }
  }
  
  if (reoptimize_hyperparameters==FALSE){
    prep_data = dataprep(df,horizon,variable,lag=best$lag,K=best$factor,dataset="B0_ARDI")
    Xin = prep_data$Xin
    yin = prep_data$yin
    best$Xout = prep_data$Xout
    best_fit <- lm(yin~Xin)
  }
  coef_opt=coef(best_fit)
  coef_opt[is.na(coef_opt)] = 0
  pred <- c(1,best$Xout)%*%coef_opt
  
  return(list(pred=pred, best=best))
}

func_shrink_rich = function (df,horizon,variable,lag_orders,reoptimize_hyperparameters=FALSE,type,best, ...) {
  
  if (reoptimize_hyperparameters==TRUE){
  # Initialize variables to store best model and criteria value
  # Here loop through different lags is necessary because 
  # ridge regression does not do variable selection on its own
    best = list()
    best_parameter = Inf
    for (K in factor_orders) {
      for (lag in lag_orders) {
        prep_data = dataprep(df,horizon,variable,lag=lag,K=K,dataset="B0_ARDI")
        Xin = prep_data$Xin
        yin = prep_data$yin
        Xout = prep_data$Xout
        if(type=="cv"){
          if(regu=="ridge"){
            fit = cv.glmnet(Xin,yin,alpha=0,type.measure = "mse", nfolds = 5)
          }
          if(regu=="lasso"){
            fit = cv.glmnet(Xin,yin,alpha=1,type.measure = "mse", nfolds = 5)
          }
          if(regu=="en"){
            fit  = cv.glmnet(Xin,yin,alpha=0.5,type.measure = "mse", nfolds = 5)
          }
          MSE = fit$cvm[fit$lambda == fit$lambda.min]
          if (MSE < best_parameter) {
            best_parameter = MSE
            best_fit = fit
            best$lag = lag
            best$factor <- K
            best$Xout = prep_data$Xout
            best$lambda = fit$cvm[fit$lambda == fit$lambda.min]
          }
        }
      }
    }
  }
  if (reoptimize_hyperparameters==FALSE){
  prep_data = dataprep(df,horizon,variable,lag=best$lag,K=best$factor,dataset="B0_ARDI")
  Xin = prep_data$Xin
  yin = prep_data$yin
  best$Xout = prep_data$Xout
  if(regu=="ridge"){
    best_fit = glmnet(Xin,yin,alpha=0)
    }
  if(regu=="lasso"){
    best_fit = glmnet(Xin,yin,alpha=1)
    }
  if(regu=="en"){
    best_fit = glmnet(Xin,yin,alpha=0.5)
    }
  }
  pred=predict(best_fit,best$Xout,s = best$lambda)
  
  return(list(pred=pred,best=best))
}

func_rf_rich = function (df,horizon,variable,lag_orders,reoptimize_hyperparameters=FALSE,type,best, ...) {
  best=list()
  prep_data = dataprep(df,horizon,variable,dataset="B0_ARDI")
  Xin = prep_data$Xin
  yin = prep_data$yin
  Xout = prep_data$Xout
  best_fit = randomForest::randomForest(Xin,yin,importance = TRUE)
  pred=predict(best_fit,Xout)
  
  return(list(pred=pred,best=best))
}

func_bols_rich = function (df,horizon,variable,lag_orders,reoptimize_hyperparameters=FALSE,type,best, ...) {
  # Linear Boosting
  prep_data = dataprep(df,horizon,variable,dataset="B0_ARDI")
  Xin = prep_data$Xin
  yin = prep_data$yin
  Xout = data.frame(prep_data$Xout)
  concat = data.frame(yin,Xin)
  bctrl <- boost_control(mstop = 300)
  if(type=="cv"){
    # Here the variable selection occurs inside the function and from 1:12
    best_fit = mboost::gamboost(yin~.,concat,family=Gaussian(), baselearner="bols", control = bctrl)
    best = list()
    cvm <- cv(model.weights(best_fit), type = "kfold", B = 5)
    cvr <- cvrisk(best_fit, folds = cvm, papply = lapply)
    best$blstop <- mstop(cvr)
  }
  pred=predict(best_fit[best$blstop],Xout)
  
  return(list(pred=pred,best=best))
}  

func_bbs_rich = function (df,horizon,variable,lag_orders,reoptimize_hyperparameters=FALSE,type,best, ...) {
  # Boost with spline (nonlinear)
  prep_data = dataprep(df,horizon,variable,dataset="B0_ARDI")
  ## smooth P-spline base-learner
  Xin = prep_data$Xin
  yin = prep_data$yin
  Xout = data.frame(prep_data$Xout)
  concat = data.frame(yin,Xin)
  bctrl <- boost_control(mstop = 300)
  df=4
  if(type=="cv"){
    # Here the variable selection occurs inside the function and from 1:12
    best_fit = mboost::gamboost(yin ~ ., concat, family = Gaussian(), baselearner="bbs",dfbase=df,control = bctrl)
    best = list()
    cvm <- cv(model.weights(best_fit), type = "kfold", B = 5)
    cvr <- cvrisk(best_fit, folds = cvm, papply = lapply)
    best$blstop <- mstop(cvr)
  }
  
  pred=predict(best_fit[best$blstop],Xout)
  
  return(list(pred=pred,best=best))
}

func_svr_linear_rich = function (df,horizon,variable,lag_orders,reoptimize_hyperparameters=FALSE,type,best, ...) {
  if (reoptimize_hyperparameters==TRUE){
    # Initialize variables to store best model and criteria value
    # Here loop through different lags is necessary because 
    # svr does not do variable selection on its own
    best = list()
    best_parameter = Inf
    for (K in factor_orders) {
      for (lag in lag_orders) {
        prep_data = dataprep(df,horizon,variable,lag=lag,K=K,dataset="B0_ARDI")
        Xin = prep_data$Xin
        yin = prep_data$yin
        Xout = prep_data$Xout
        if(type=="cv"){
          obj <- tune(svm,Xin,yin,
                      ranges = list(cost = c(0.1,0.5,1,2,5),
                                    epsilon=seq(0.1, 0.5, by = 0.1)),
                      tunecontrol = tune.control(sampling = "fix")
          )
          # Scale =TRUE?
          fit = svm(Xin,yin, kernel = "linear", type= "eps-regression",cross=5, scale = FALSE,
                    cost = obj$best.parameters$cost,
                    epsilon = obj$best.parameters$epsilon)
          MSE = mean(fit$MSE)
          if (MSE < best_parameter) {
            best_parameter = MSE
            best_fit = fit
            best$lag = lag
            best$factor <- K
            best$Xout = prep_data$Xout
            best$cost = fit$cost
            best$epsilon=fit$epsilon
          }
        }
      }
    }
  }
  if (reoptimize_hyperparameters==FALSE){
    prep_data = dataprep(df,horizon,variable,lag=best$lag,K=best$factor,dataset="B0_ARDI")
    Xin = prep_data$Xin
    yin = prep_data$yin
    best$Xout = prep_data$Xout
    best_fit = svm(Xin,yin, kernel = "linear", type= "eps-regression",scale = FALSE,
                   cost=best$cost, 
                   epsilon=best$epsilon)
  }
  
  pred=predict(best_fit,best$Xout)
  
  return(list(pred=pred,best=best))
}    


func_svr_rbf_rich = function (df,horizon,variable,lag_orders,reoptimize_hyperparameters=FALSE,type,best, ...) {
  if (reoptimize_hyperparameters==TRUE){
    # Initialize variables to store best model and criteria value
    # Here loop through different lags is necessary because 
    # svr does not do variable selection on its own
    best = list()
    best_parameter = Inf
    for (K in factor_orders) {
      for (lag in lag_orders) {
        prep_data = dataprep(df,horizon,variable,lag=lag,K=K,dataset="B0_ARDI")
        Xin = prep_data$Xin
        yin = prep_data$yin
        Xout = prep_data$Xout
        if(type=="cv"){
          obj <- tune(svm,Xin,yin,
                      ranges = list(gamma = c(seq(0.1, 0.9, by = 0.1),(1/ncol(Xin))), 
                                    cost = c(0.1,0.5,1,2,5),
                                    epsilon=seq(0.1, 0.5, by = 0.1)),
                      tunecontrol = tune.control(sampling = "fix")
          )
          fit = svm(Xin,yin, kernel = "radial", type= "eps-regression",cross=5, scale = FALSE,
                    gamma = obj$best.parameters$gamma,
                    cost = obj$best.parameters$cost,
                    epsilon = obj$best.parameters$epsilon)
          MSE = mean(fit$MSE)
          if (MSE < best_parameter) {
            best_parameter = MSE
            best_fit = fit
            best$lag = lag
            best$factor <- K
            best$Xout = prep_data$Xout
            best$cost = fit$cost
            best$epsilon=fit$epsilon
            best$gamma=fit$gamma
          }
        }
      }
    }
  }
  if (reoptimize_hyperparameters==FALSE){
    prep_data = dataprep(df,horizon,variable,lag=best$lag,K=best$factor,dataset="B0_ARDI")
    Xin = prep_data$Xin
    yin = prep_data$yin
    best$Xout = prep_data$Xout
    best_fit = svm(Xin,yin, kernel = "radial", type= "eps-regression",scale = FALSE,
                   gamma=best$gamma,
                   cost=best$cost, 
                   epsilon=best$epsilon)
  }
  
  pred=predict(best_fit,best$Xout)
  
  return(list(pred=pred,best=best))
}    
  
  
###################################
### Generating shrinkage schemes###
###################################  
  
func_b1= function (df,horizon,variable,lag_orders,reoptimize_hyperparameters=FALSE,type,best,regu) {
  prep_data = dataprep(df,horizon,variable,dataset="B1")
  Xin = prep_data$Xin
  yin = prep_data$yin
  Xout = prep_data$Xout
  if(regu=="boost"){
    Xout = data.frame(Xout)
    concat = data.frame(yin,Xin)
    bctrl <- boost_control(mstop = 300)
    if(type=="cv"){
      # Here the variable selection occurs inside the function and from 1:12
      best_fit = mboost::gamboost(yin~.,concat,family=Gaussian(), baselearner="bols", control = bctrl)
      best = list()
      cvm <- cv(model.weights(best_fit), type = "kfold", B = 5)
      cvr <- cvrisk(best_fit, folds = cvm, papply = lapply)
      best$blstop <- mstop(cvr)
    }
    pred=predict(best_fit[best$blstop],Xout)
  }
  
  else{
    if(type=="cv"){
      best = list()
      if(regu=="ridge"){
        best_fit = cv.glmnet(Xin,yin,alpha=0,type.measure = "mse", nfolds = 5)
      }
      if(regu=="lasso"){
        best_fit = cv.glmnet(Xin,yin,alpha=1,type.measure = "mse", nfolds = 5)
      }
      if(regu=="en"){
        best_fit  = cv.glmnet(Xin,yin,alpha=0.5,type.measure = "mse", nfolds = 5)
      }
      best$lambda = best_fit$cvm[best_fit$lambda == best_fit$lambda.min]
    }
    pred=predict(best_fit,Xout,s = best$lambda)
  }
  return(list(pred=pred,best=best))
}
  
func_b2= function (df,horizon,variable,lag_orders,reoptimize_hyperparameters=FALSE,type,best,regu) {
  prep_data = dataprep(df,horizon,variable,dataset="B2")
  Xin = prep_data$Xin
  yin = prep_data$yin
  Xout = prep_data$Xout
  if(regu=="boost"){
    Xout = data.frame(Xout)
    concat = data.frame(yin,Xin)
    bctrl <- boost_control(mstop = 300)
    if(type=="cv"){
      # Here the variable selection occurs inside the function and from 1:12
      best_fit = mboost::gamboost(yin~.,concat,family=Gaussian(), baselearner="bols", control = bctrl)
      best = list()
      cvm <- cv(model.weights(best_fit), type = "kfold", B = 5)
      cvr <- cvrisk(best_fit, folds = cvm, papply = lapply)
      best$blstop <- mstop(cvr)
    }
    pred=predict(best_fit[best$blstop],Xout)
  }
  
  else{
    if(type=="cv"){
      best = list()
      if(regu=="ridge"){
        best_fit = cv.glmnet(Xin,yin,alpha=0,type.measure = "mse", nfolds = 5)
      }
      if(regu=="lasso"){
        best_fit = cv.glmnet(Xin,yin,alpha=1,type.measure = "mse", nfolds = 5)
      }
      if(regu=="en"){
        best_fit  = cv.glmnet(Xin,yin,alpha=0.5,type.measure = "mse", nfolds = 5)
      }
      best$lambda = best_fit$cvm[best_fit$lambda == best_fit$lambda.min]
    }
    pred=predict(best_fit,Xout,s = best$lambda)
  }
  return(list(pred=pred,best=best))
}

func_b3= function (df,horizon,variable,lag_orders,reoptimize_hyperparameters=FALSE,type,best,regu) {
  prep_data = dataprep(df,horizon,variable,dataset="B3")
  Xin = prep_data$Xin
  yin = prep_data$yin
  Xout = prep_data$Xout
  if(regu=="boost"){
    Xout = data.frame(Xout)
    concat = data.frame(yin,Xin)
    bctrl <- boost_control(mstop = 300)
    if(type=="cv"){
      # Here the variable selection occurs inside the function and from 1:12
      best_fit = mboost::gamboost(yin~.,concat,family=Gaussian(), baselearner="bols", control = bctrl)
      best = list()
      cvm <- cv(model.weights(best_fit), type = "kfold", B = 5)
      cvr <- cvrisk(best_fit, folds = cvm, papply = lapply)
      best$blstop <- mstop(cvr)
    }
    pred=predict(best_fit[best$blstop],Xout)
  }
  
  else{
    if(type=="cv"){
      best = list()
      if(regu=="ridge"){
        best_fit = cv.glmnet(Xin,yin,alpha=0,type.measure = "mse", nfolds = 5)
      }
      if(regu=="lasso"){
        best_fit = cv.glmnet(Xin,yin,alpha=1,type.measure = "mse", nfolds = 5)
      }
      if(regu=="en"){
        best_fit  = cv.glmnet(Xin,yin,alpha=0.5,type.measure = "mse", nfolds = 5)
      }
      best$lambda = best_fit$cvm[best_fit$lambda == best_fit$lambda.min]
    }
    pred=predict(best_fit,Xout,s = best$lambda)
  }
  return(list(pred=pred,best=best))
}

source("libraries.R")
which_vector_min <- function(...) {
  #takes several vectors of same length as arguments
  #and returns a vector indices_min of same length
  #indices_min[i]= the vector name with the minimum
  #value at index i out of the others
  matrice <- do.call(cbind, list(...))
  indices_min <- apply(matrice, 1, which.min)
  return(indices_min)
}
find_truncnorm_params <- function(params, breaks, probs) {
  mu <- params[1]
  sigma <- params[2]
  p1 <- pnorm(breaks[1], mu, sigma) 
  p2 <- pnorm(breaks[2], mu, sigma)
  p3 <- pnorm(breaks[3], mu, sigma)
  p4 <- pnorm(breaks[4], mu, sigma)
  
  actual_probs <- c(p2 - p1, p3 - p2, p4 - p3)
  actual_probs <- actual_probs / sum(actual_probs) # Normaliser
  
  sum((actual_probs - probs)^2)
}


hill_estimator <- function(x, k) {
  x_sorted <- sort(x, decreasing = TRUE)
  log_ratios <- log(x_sorted[1:k]) - log(x_sorted[k + 1])
  return(mean(log_ratios))
}

cik_logweiss<-function(data,K,p,alpha=.05,int="NCI5"){
  #data: sample of data (vector) 
  #K: a range of k_n
  #p: probability of the quantile to estimate
  #alpha: error rate for the confidence interval
  #names: can be either NCI1, NCI2, NCI3, NCI4, NCI5, GCI1, GCI2
  #default: NCI5
  n<-length(data)
  hills<-Hill(data)$gamma
  cso<-mop(data,c(1),c(1),"RBMOP")
  b<-cso$beta
  rho<-cso$rho
  weiss<-Weissman.q(data,hills,p)$Q[K]
  hills<-hills[K]
  L<-sqrt(K)/log(K/(n*p))
  lambdas<-sqrt(K)*b*hills*((n/K)**rho)
  if(int=="NCI1"){
    c<-0
    sigmas<-hills
  }else if(int=="NCI2"){
    c<-lambdas/(1-rho)
    sigmas<-hills
  }else if(int=="NCI3"){
    c<-lambdas/(1-rho)
    sigmas<-sqrt(hills**2 + 2*lambdas*hills/(sqrt(K)* (1-rho)**2) + lambdas**2/(K*(1-rho)**2 * (1-2*rho)))
  }else if(int=="NCI4"){
    dns<-K/n/p
    c<-lambdas/(1-rho)+lambdas*(1-dns**rho)/(rho*log(dns))
    sigmas<-sqrt(hills**2*(1+1/log(dns)**2))
    }else if(int=="NCI5"){
      dns<-K/n/p
      c<-lambdas/(1-rho)+lambdas*(1-dns**rho)/(rho*log(dns))
      sigmas<- sqrt(hills**2*(1+1/log(dns)**2) + 2*lambdas*hills/(sqrt(K)* (1-rho)**2) + lambdas**2/(K*(1-rho)**2 * (1-2*rho)))
    } else if(int=="GCI1"){
      c<-0
      g1<-qgamma(alpha/2, shape=K, scale = 1/sqrt(K))
      g2<-qgamma(1-alpha/2, shape=K, scale = 1/sqrt(K))
      l_bounds<-log(weiss)-c/L -hills/L * (g2 - sqrt(K))
      u_bounds<- log(weiss) -c/L-hills/L * (g1- sqrt(K))
      return(list(
        "k"=K,
        "error"=alpha,
        "ci"=int,
        "lbounds"=l_bounds,
        "ubounds"=u_bounds,
        "b"=b,
        "rho"=rho
      ))
    } else if(int=="GCI2"){
      c<-lambdas/(1-rho)
      g1<-qgamma(alpha/2, shape=K, scale = 1/sqrt(K))
      g2<-qgamma(1-alpha/2, shape=K, scale = 1/sqrt(K))
      l_bounds<-log(weiss)-c/L -hills/L * (g2 - sqrt(K))
      u_bounds<- log(weiss) -c/L-hills/L * (g1- sqrt(K))
      return(list(
        "k"=K,
        "error"=alpha,
        "ci"=int,
        "lbounds"=l_bounds,
        "ubounds"=u_bounds,
        "b"=b,
        "rho"=rho
      ))
      }
    q<-qnorm(1-alpha/2)
    l_bounds<-log(weiss)-1/L *(sigmas*q +c)
    u_bounds<- log(weiss)+1/L * (sigmas*q-c)
    return(list(
      "k"=K,
      "error"=alpha,
      "ci"=int,
      "lbounds"=l_bounds,
      "ubounds"=u_bounds,
      "b"=b,
      "rho"=rho
    ))
}

cik_beta<-function(data,K,p,a=.85,alpha=.05){
  #data: sample of data (vector) 
  #K: a range of k_n
  #p: probability of the quantile to estimate
  #alpha: error rate for the confidence interval
  #confidence interval by Gardes, Maistre
  n<-length(data)
  hills<-Hill(data)$gamma[K]
  an<-pmax(3,log(K)^a)/n #choice from the simulation section
  dataquant<-sort(data)[n-floor(n*an)]
  tnl<-1/p*qbeta(alpha/2,floor(n*an)+1,n-floor(n*an))
  tnr<- 1/p*qbeta(1-alpha/2,floor(n*an)+1,n-floor(n*an))
  l_bounds<-dataquant*tnl^hills
  u_bounds<-dataquant*tnr^hills
  return(list(
    "k"=K,
    "error"=alpha,
    "lbounds"=l_bounds,
    "ubounds"=u_bounds
  ))
}
algopercentile<-function(n,gamma,rho,p,
                              M=1000,kmax=NULL,distr="burr",alpha=.05,seed=NULL){
      #Algorithm 1 for NCI1, 3 and NCI5
      #n: length of the samples
      #gamma>0,rho<0: parameters of the distribution
      #p: target probability that X>quantile
      #M: number of replications
      #kmax: maximum number of k_n to consider
      #distr: either 'burr', 'frechet', 'invgamma', 'fisher' or 'student'
      #seed: a seed for reproducibility for the replications
      if(is.null(kmax)){
        kmax<-floor(3*n/10)
      }
      set.seed(seed=seed)
      if(distr=="burr"){
      variables<-rburr(n*M,1/gamma,rho)
      quant<-log(qburr(1-p,1/gamma,rho))
      }else if(distr=="frechet"){
        variables<-rfrechet(n*M,1/gamma)
        quant<-log(qfrechet(1-p,1/gamma))
      }else if(distr=="invgamma"){
        variables<-rinvgamma(n*M,-1/rho)
        quant<-log(qinvgamma(1-p,-1/rho))
      }else if(distr=="fisher"){
        df1<-1/gamma #arbitrary value  
        df2<- -2/rho
        variables<-rf(n*M,df1,df2)
        quant<-log(qf(1-p,df1,df2))
      }else if(distr=="student"){
        df<- -2/rho
        variables<-abs(rt(n*M,df))
        quant<-log(qt(1-p/2,df))
      }else{
        stop("The only distributions allowed are 'burr', 'frechet', 'invgamma', 'fisher' and 'student'.")
      }
      #quant=target quantile
      Data<-matrix(variables,nrow = M, ncol= n)
      #Data= matrix of the data
      Errs_n1<-matrix(0,nrow=M,ncol=kmax)
      Errs_n3<-matrix(0,nrow=M,ncol=kmax)
      Errs_n5<-matrix(0,nrow=M,ncol=kmax)
      #matrixes of errors: 
      #M[i,k]=1 means that for replication i and index k
      #the interval M does not contain the target (log) quantile
      for(i in 1:M){
        nci1<-cik_logweiss(Data[i,],1:kmax,p,alpha=alpha,bias = FALSE)
        nci3<-cik_logweiss(Data[i,],1:kmax,p,alpha=alpha,bias=TRUE,enbias = FALSE)
        nci5<-cik_logweiss(Data[i,],1:kmax,p,alpha=alpha,bias=TRUE,enbias=TRUE)
        lbounds_nci1<-nci1$lbounds
        ubounds_nci1<-nci1$ubounds
        lbounds_nci3<-nci3$lbounds
        ubounds_nci3<-nci3$ubounds
        lbounds_nci5<-nci5$lbounds
        ubounds_nci5<-nci5$ubounds
        Errs_n1[i,quant<lbounds_nci1|quant>ubounds_nci1]<-1
        Errs_n3[i,quant<lbounds_nci3|quant>ubounds_nci3]<-1
        Errs_n5[i,quant<lbounds_nci5|quant>ubounds_nci5]<-1
      }
      errs_n1<-colSums(Errs_n1)/M*100
      errs_n3<-colSums(Errs_n3)/M*100
      errs_n5<-colSums(Errs_n5)/M*100
      #errs_n = average vector of percentage fo errors 
      #for all replications and each value of k_n between 1 and kmax
      pn1<-quantile(errs_n1,.75)
      pn3<-quantile(errs_n3,.75)
      pn5<-quantile(errs_n5,.75)
      return(list(
        "Replications"=M,
        "pn1"=unname(pn1),"pn3"=unname(pn3), "pn5"=unname(pn5),
        "quantile"=quant
      ))
      #returns the 75% percentiles for the intervals
      #and the true quantile of the Burr distribution
}


algopercentile_ci<-function(n,gamma,rho,p,
                         M=1000,kmax=NULL,distr="burr",seed=NULL,
                         loginterval=TRUE,
                         interval=cik_logweiss,...){
  #Algorithm 1 for any kind of confidence interval (that returns the same variables as cik_logweiss)
  #n: length of the samples
  #gamma>0,rho<0: parameters of the distribution
  #p: target probability that X>quantile
  #M: number of replications
  #kmax: maximum number of k_n to consider
  #distr: either 'burr', 'frechet', 'invgamma', 'fisher' or 'student'
  #seed: a seed for reproducibility for the replications
  #loginterval: for logarithm or real quantile
  #interval: function computing a confidence interval that should have:
  #data as a vector and a vector of k_n for which the ci is computed
  #the interval function must return a list with value:
  #lbounds for the lower bounds of the ci and ubounds for the upper ones
  #... additional arguments for the considered interval function
  if(is.null(kmax)){
    kmax<-floor(3*n/10)
  }
  set.seed(seed=seed)
  if(distr=="burr"){
    variables<-rburr(n*M,1/gamma,rho)
    quant<-log(qburr(1-p,1/gamma,rho))
  }else if(distr=="frechet"){
    variables<-rfrechet(n*M,1/gamma)
    quant<-log(qfrechet(1-p,1/gamma))
  }else if(distr=="invgamma"){
    variables<-rinvgamma(n*M,-1/rho)
    quant<-log(qinvgamma(1-p,-1/rho))
  }else if(distr=="fisher"){
    df1<-1/gamma #arbitrary value  
    df2<- -2/rho
    variables<-rf(n*M,df1,df2)
    quant<-log(qf(1-p,df1,df2))
  }else if(distr=="student"){
    df<- -2/rho
    variables<-abs(rt(n*M,df))
    quant<-log(qt(1-p/2,df))
  }else{
    stop("The only distributions allowed are 'burr', 'frechet', 'invgamma', 'fisher' and 'student'.")
  }
  if(!loginterval){
    quant<-exp(quant)
  }
  #quant=target quantile
  Data<-matrix(variables,nrow = M, ncol= n)
  #Data= matrix of the data
  Errs_ci<-matrix(0,nrow=M,ncol=kmax)
  #matrixes of errors: 
  #M[i,k]=1 means that for replication i and index k
  #the interval M does not contain the target (log) quantile
  for(i in 1:M){
    ci<-interval(Data[i,],1:kmax,p,...)
    lbounds<-ci$lbounds
    ubounds<-ci$ubounds
    Errs_ci[i,quant<lbounds|quant>ubounds]<-1
  }
  errs_ci<-colSums(Errs_ci)/M*100
  #errs_ci = average vector of percentage fo errors 
  #for all replications and each value of k_n between 1 and kmax
  pci<-quantile(errs_ci,.75)
  return(list(
    "Replications"=M,
    "percentile"=unname(pci),
    "true_quantile"=quant
  ))
  #returns the 75% percentiles for the intervals
  #and the true quantile of the Burr distribution
}

data_for_classif<-function(Data,kmax){
  #transforms the data set into the data D_i for i=1,...kmax
  #in order to apply the classification method
  #data is a vector or a matrix or a data frame whose rows are samples of data
  #kmax is the length of the set to classify
  if (is.vector(Data)) {
    Data <- matrix(Data, nrow = 1)
    vec<-TRUE
  }
  n<-length(Data[1,])
  N<-length(Data[,1])
  Sets<-matrix(0,nrow=N,ncol=kmax)
  for(i in 1:N){
    Datord<-sort(Data[i,])
    Sets[i,]<- c(1:kmax)*(log(Datord[n-c(1:kmax)+1])
        -log(Datord[n-c(1:kmax)]))/hill_estimator(Data[i,],kmax+1)
  }
  return(Sets)
}

training_set<-function(n,gammas,rhos,kmax=NULL,distr="burr",seeds=NULL){
  #generates the training set for the classification method
  #gammas: vector of extreme value indices>0
  #rhos: vector of parameter rho<0 (same length as gammas)
  #kmax: number of k_n to include, default: 3n/10
  #distr: either 'burr', 'frechet', 'invgamma', 'fisher' or 'student'
  #seeds: a vector of seeds for reproducibility of the set
  N<-length(gammas)
  if(is.null(kmax)){
    kmax<-floor(3*n/10)
  }
  if(length(rhos)!=N){
    stop("rhos and gammas should be vectors of the same length")
  }
  if(length(seeds)<N){
    seeds<-rep(seeds[1],N)
  }
  Data_train<-matrix(0,nrow=N,ncol=n)
  for( i in 1:N){
    gamma<-gammas[i]
    rho<-rhos[i]
    set.seed(seeds[i])
    if(distr=="burr"){
      Data_train[i,]<-rburr(n,1/gamma,rho)
    }else if(distr=="frechet"){
      Data_train[i,]<-rfrechet(n,1/gamma)
    }else if(distr=="invgamma"){
      Data_train[i,]<-rinvgamma(n,-1/rho)
    }else if(distr=="fisher"){
      df1<-1/gamma #arbitrary value  
      df2<- -2/rho
      Data_train[i,]<-rf(n,df1,df2)
    }else if(distr=="student"){
      df<- -2/rho
      Data_train[i,]<-abs(rt(n,df))
    }else{
      stop("The only distributions allowed are 'burr', 'frechet', 'invgamma', 'fisher' and 'student'.")
    }
  }
  Sets_train<-data_for_classif(Data_train,kmax)
  return(list("Sets_train"=Sets_train,
              "Data"=Data_train
  ))
}


predict_ci<-function(model,data,K,p,alpha=.05){
  #uses the classification method 'model' to return the confidence interval predicted by the model
  #model: a model of classification of 3 classes (lda, svm, hdda, knn):
  #class 1: NCI1, class 2: NCI3, class 3: NCI5.
  #data: VECTOR for which the confidence interval is predicted
  #kmax: max value of considered k_n
  #p probability of the quantile to predict 
  #alpha: error level
  if(!is.vector(data)){
    stop("data should be a numerical vector.")
  }
  data_matrix <- matrix(data, nrow = 1)
  set_test<-data_for_classif(data_matrix, max(K))
  pred<-predict(model,data.frame(set_test))
  #this version is adapted to LDA, SVM, KNN or HDDA
  #for other versions the "predict" function should have adapted arguments
  pred_label <- if (is.list(pred) && "class" %in% names(pred)) {
    pred$class #for lda method it is a class
  } else {
    pred  #for svm method
  }
  if(pred_label==1){
      return(cik_logweiss(data,K,p,alpha,int="NCI1"))
  }
  if(pred_label==2){
    return(cik_logweiss(data,K,p,alpha,int="NCI3"))
  } else{
    #if pred==3
    return(cik_logweiss(data,K,p,alpha,int="NCI5"))
  }
}


predictnn_ci<-function(keras_model,data,K,p,alpha=.05){
  #uses the classification method 'model' to return the confidence interval predicted by the model
  #it is adapted for the Neural network using tensorflow
  #model: a model of classification of 3 classes that is a keras_model
  #class 1: NCI1, class 2: NCI3, class 3: NCI5.
  #data: VECTOR for which the confidence interval is predicted
  #kmax: max value of considered k_n
  #p probability of the quantile to predict 
  #alpha: error level
  if(!is.vector(data)){
    stop("data should be a numerical vector.")
  }
  data_matrix <- matrix(data, nrow = 1)
  set_test<-data_for_classif(data_matrix, max(K))
  pred<-predict(keras_model,matrix(set_test,nrow=1),verbose=0) #not a data frame
  #this version is adapted to Neural network from tensorflow
  pred_label<-max.col(pred)
    if(pred_label==1){ 
    return(cik_logweiss(data,K,p,alpha,bias=FALSE))
  }
  if(pred_label==2){
    return(cik_logweiss(data,K,p,alpha,bias=TRUE,enbias=FALSE))
  } else{
    #if pred_label==3
    return(cik_logweiss(data,K,p,alpha,bias=TRUE,enbias = TRUE))
  }
}


confidence_logweiss_nor<- function(weiss,alpha,L,c,sigma){
  #function to compute a normal confidence interval for algo_percentile_opt
  #weiss: weissman estimator
  #L normalisation factor
  #c: bias term
  #sigma: sqrt of variance term
  q<-qnorm(1-alpha/2)
  binf<-log(weiss)-1/L *(sigma*q +c)
  bsup<- log(weiss)+1/L * (sigma*q-c)
  return(c(binf,bsup))
}

confidence_logweiss_gam<-function(weiss,k,alpha,L,c,hill){
  #function to compute a gamma confidence interval for algo_percentile_opt
  #weiss: weissman estimator
  #k: subsequence used 
  #L normalisation factor
  #c: bias term
  #hill: hill estimator
  
  g1<-qgamma(alpha/2, shape=k, scale = 1/sqrt(k))
  g2<-qgamma(1-alpha/2, shape=k, scale = 1/sqrt(k))
  binf<-log(weiss)-c/L -hill/L * (g2 - sqrt(k))
  bsup<- log(weiss) -c/L-hill/L * (g1- sqrt(k))
  return(c(binf,bsup))
}

#Changement: percentile_opt -> algopercentile_opt
algopercentiles_opt<-function(model,n,gamma,rho,p,M=1000,kmax=NULL,distr="burr",alpha=.05,seed=NULL){
  #alternative way of computing algo 1 that is faster for some classification methods
  #model: a model of classification of 3 classes (lda, svm, hdda, knn):
  #n: length of the samples
  #gamma>0,rho<0: parameters of the distribution
  #p: target probability that X>quantile
  #M: number of replications
  #kmax: maximum number of k_n to consider
  #distr: either 'burr', 'frechet', 'invgamma', 'fisher' or 'student'
  #seed: a seed for reproducibility for the replications
  if(is.null(kmax)){
    kmax<-floor(3*n/10)
  }
  set.seed(seed=seed)
  if(distr=="burr"){
    Variables<-rburr(n*M,1/gamma,rho)
    quant<-log(qburr(1-p,1/gamma,rho))
  }else if(distr=="frechet"){
    Variables<-rfrechet(n*M,1/gamma)
    quant<-log(qfrechet(1-p,1/gamma))
  }else if(distr=="invgamma"){
    Variables<-rinvgamma(n*M,-1/rho)
    quant<-log(qinvgamma(1-p,-1/rho))
  }else if(distr=="fisher"){
    df1<-1/gamma #arbitrary value  
    df2<- -2/rho
    Variables<-rf(n*M,df1,df2)
    quant<-log(qf(1-p,df1,df2))
  }else if(distr=="student"){
    df<- -2/rho
    Variables<-abs(rt(n*M,df))
    quant<-log(qt(1-p/2,df))
  }else{
    stop("The only distributions allowed are 'burr', 'frechet', 'invgamma', 'fisher' and 'student'.")
  }
  Variables<-matrix(Variables,nrow = M, ncol= n)
  Hills<-matrix(0,nrow=M,ncol=n-1)
  Weiss<-matrix(0,nrow=M,ncol=n-1)
  
  for (i in 1:M){ #Hill and Weissman matrices
    Hills[i,]<-Hill(Variables[i,])$gamma
    Weiss[i,]<-Weissman.q(Variables[i,],Hills[i,],p)$Q
  }
  Variablesord<-matrix(0,nrow = M, ncol= n)
  Sets_test<-matrix(0,nrow=M,ncol = kmax)
  
  for(i in 1:n){
    Variablesord[i,]<-sort(Variables[i,])
    Sets_test[i,]<- c(1:kmax)*(log(Variablesord[i,n-c(1:kmax)+1])
                               -log(Variablesord[i,n-c(1:kmax)]))/Hills[i,kmax+1] 
  }
  preds_model<-predict(model,data.frame(Sets_test))
  #this version is adapted to LDA, SVM, KNN or HDDA
  #for other versions the "predict" function should have adapted arguments
  preds_model <- if (is.list(preds_model) && "class" %in% names(preds_model)) {
    preds_model$class #for lda method it is a class
  } else {
    preds_model  #for svm method
  }
  bhats<-numeric(M)#b estimators
  rhohats<-numeric(M) #rho estimators
  
  for(i in 1:M){ #2nd order parameters vectors bhats and rhohats
    cso<-mop(Variables[i,],c(1),c(1),"RBMOP")
    bhats[i]<-cso$beta
    rhohats[i]<-cso$rho
  }
  
  errs_model<-numeric(kmax)
  for(k in 1:kmax){ 
    err_model<-numeric(M)
    dn<-k/n/p
    L<-sqrt(k)/log(dn) #normalisation factor
    for (x in 1:M){ 
      bhat<-bhats[x] #parameters for the CI
      rhohat<-rhohats[x]
      weissman<-Weiss[x,k]
      hill<-Hills[x,k]
      
      #predictions
      pred_model<-preds_model[x]
      #Confidence intervals to compute
      if(pred_model==1){
        c<-0
        sigmahat<-hill
        ci<-confidence_logweiss_nor(weissman,alpha = alpha,L,c,sigmahat)
        ci_inf<-ci[1]
        ci_sup<-ci[2]
      } else if(pred_model==2){
        lambdahat<-sqrt(k)*bhat*hill*((n/k)**rhohat)
        c<-lambdahat/(1-rhohat)
        sigmahat<-sqrt(hill**2 + 2*lambdahat*hill/(sqrt(k)* (1-rhohat)**2) + lambdahat**2/(k*(1-rhohat)**2 * (1-2*rhohat)))
        ci<-confidence_logweiss_nor(weissman,alpha = alpha,L,c,sigmahat)
        ci_inf<-ci[1]
        ci_sup<-ci[2]
      } else {  #pred_model==3
        lambdahat<-sqrt(k)*bhat*hill*((n/k)**rhohat)
        c<-lambdahat/(1-rhohat)+lambdahat*(1-dn**rhohat)/(rhohat*log(dn))
        sigmahat<- sqrt(hill**2*(1+1/log(dn)**2) + 2*lambdahat*hill/(sqrt(k)* (1-rhohat)**2) + lambdahat**2/(k*(1-rhohat)**2 * (1-2*rhohat)))
        ci<-confidence_logweiss_nor(weissman,alpha = alpha,L,c,sigmahat)
        ci_inf<-ci[1]
        ci_sup<-ci[2]
      }
      if(quant<ci_inf){ #donner une erreur si la borne n'est pas au bon endroit
        err_model[x]<-1
      } else if(quant>ci_sup){
        err_model[x]<-1
      }
    }
    errs_model[k]<-sum(err_model)/M*100
  }
  
  pmodel<-quantile(errs_model,.75)
  
  Li<-list( "percentile"=unname(pmodel),
            "errsmodel"=errs_model,
            "ysmodel"=preds_model
  )
  return(Li)
}

algopercentilesnn_opt<-function(keras_model,n,gamma,rho,p,M=1000,kmax=NULL,distr="burr",alpha=.05,seed=NULL){
  #alternative way of computing algo 1 that is faster for some classification methods and adapted only for
  #Neural network using tensorflow
  #model: a model of classification of 3 classes that is a keras_model
  #n: length of the samples
  #gamma>0,rho<0: parameters of the distribution
  #p: target probability that X>quantile
  #M: number of replications
  #kmax: maximum number of k_n to consider
  #distr: either 'burr', 'frechet', 'invgamma', 'fisher' or 'student'
  #seed: a seed for reproducibility for the replications
  if(is.null(kmax)){
    kmax<-floor(3*n/10)
  }
  set.seed(seed=seed)
  if(distr=="burr"){
    Variables<-rburr(n*M,1/gamma,rho)
    quant<-log(qburr(1-p,1/gamma,rho))
  }else if(distr=="frechet"){
    Variables<-rfrechet(n*M,1/gamma)
    quant<-log(qfrechet(1-p,1/gamma))
  }else if(distr=="invgamma"){
    Variables<-rinvgamma(n*M,-1/rho)
    quant<-log(qinvgamma(1-p,-1/rho))
  }else if(distr=="fisher"){
    df1<-1/gamma #arbitrary value  
    df2<- -2/rho
    Variables<-rf(n*M,df1,df2)
    quant<-log(qf(1-p,df1,df2))
  }else if(distr=="student"){
    df<- -2/rho
    Variables<-abs(rt(n*M,df))
    quant<-log(qt(1-p/2,df))
  }else{
    stop("The only distributions allowed are 'burr', 'frechet', 'invgamma', 'fisher' and 'student'.")
  }
  Variables<-matrix(Variables,nrow = M, ncol= n)
  Hills<-matrix(0,nrow=M,ncol=n-1)
  Weiss<-matrix(0,nrow=M,ncol=n-1)
  
  for (i in 1:M){ #Hill and Weissman matrices
    Hills[i,]<-Hill(Variables[i,])$gamma
    Weiss[i,]<-Weissman.q(Variables[i,],Hills[i,],p)$Q
  }
  Variablesord<-matrix(0,nrow = M, ncol= n)
  Sets_test<-matrix(0,nrow=M,ncol = kmax)
  
  for(i in 1:n){
    Variablesord[i,]<-sort(Variables[i,])
    Sets_test[i,]<- c(1:kmax)*(log(Variablesord[i,n-c(1:kmax)+1])
                               -log(Variablesord[i,n-c(1:kmax)]))/Hills[i,kmax+1] 
  }
  #nn
  preds_keras <- predict(keras_model, Sets_test,verbose=0)
  preds_keras <- max.col(preds_keras)  # return the index of the maximum
  
  bhats<-numeric(M)#b estimators
  rhohats<-numeric(M) #rho estimators
  
  for(i in 1:M){ #2nd order parameters vectors bhats and rhohats
    cso<-mop(Variables[i,],c(1),c(1),"RBMOP")
    bhats[i]<-cso$beta
    rhohats[i]<-cso$rho
  }
  
  errs_keras<-numeric(kmax)
  for(k in 1:kmax){ 
    err_keras<-numeric(M)
    dn<-k/n/p
    L<-sqrt(k)/log(dn) #normalisation factor
    for (x in 1:M){ 
      bhat<-bhats[x] #parameters for the CI
      rhohat<-rhohats[x]
      weissman<-Weiss[x,k]
      hill<-Hills[x,k]
      
      #predictions
      pred_keras<-preds_keras[x]
      #Confidence intervals to compute
      if(pred_keras==1){
        c<-0
        sigmahat<-hill
        ci<-confidence_logweiss_nor(weissman,alpha = alpha,L,c,sigmahat)
        ci_inf<-ci[1]
        ci_sup<-ci[2]
      } else if(pred_keras==2){
        lambdahat<-sqrt(k)*bhat*hill*((n/k)**rhohat)
        c<-lambdahat/(1-rhohat)
        sigmahat<-sqrt(hill**2 + 2*lambdahat*hill/(sqrt(k)* (1-rhohat)**2) + lambdahat**2/(k*(1-rhohat)**2 * (1-2*rhohat)))
        ci<-confidence_logweiss_nor(weissman,alpha = alpha,L,c,sigmahat)
        ci_inf<-ci[1]
        ci_sup<-ci[2]
      } else {  #pred_keras==3
        lambdahat<-sqrt(k)*bhat*hill*((n/k)**rhohat)
        c<-lambdahat/(1-rhohat)+lambdahat*(1-dn**rhohat)/(rhohat*log(dn))
        sigmahat<- sqrt(hill**2*(1+1/log(dn)**2) + 2*lambdahat*hill/(sqrt(k)* (1-rhohat)**2) + lambdahat**2/(k*(1-rhohat)**2 * (1-2*rhohat)))
        ci<-confidence_logweiss_nor(weissman,alpha = alpha,L,c,sigmahat)
        ci_inf<-ci[1]
        ci_sup<-ci[2]
      }
      if(quant<ci_inf){ #donner une erreur si la borne n'est pas au bon endroit
        err_keras[x]<-1
      } else if(quant>ci_sup){
        err_keras[x]<-1
      }
    }
    errs_keras[k]<-sum(err_keras)/M*100
  }
  
  pkeras<-quantile(errs_keras,.75)
  
  Li<-list( "percentile"=unname(pkeras),
           "errskeras"=errs_keras,
           "yskeras"=preds_keras
  )
  return(Li)
}
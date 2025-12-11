
# No change
#use of the intermediate functions
addCommas              <- HDclassif:::addCommas
addCommas_single       <- HDclassif:::addCommas_single
checkTheTypes          <- HDclassif:::checkTheTypes
default_kmeans_control <- HDclassif:::default_kmeans_control
demo_hddc              <- HDclassif:::demo_hddc
demo_hddc_acp          <- HDclassif:::demo_hddc_acp
demo_hddc_crabs        <- HDclassif:::demo_hddc_crabs
getHDclassif.show      <- HDclassif:::getHDclassif.show
hdc_getComplexity      <- HDclassif:::hdc_getComplexity
hdc_getTheModel        <- HDclassif:::hdc_getTheModel
hdc_myEigen            <- HDclassif:::hdc_myEigen
HDC_plot_criteria      <- HDclassif:::HDC_plot_criteria
hdclassif_bic          <- HDclassif:::hdclassif_bic
hdclassif_dim_choice   <- HDclassif:::hdclassif_dim_choice
hdda_control           <- HDclassif:::hdda_control
hdda_prms              <- HDclassif:::hdda_prms
hdda_prms_bis          <- HDclassif:::hdda_prms_bis
hddc_ari               <- HDclassif:::hddc_ari
hddc_control           <- HDclassif:::hddc_control
hddc_e_step            <- HDclassif:::hddc_e_step
hddc_m_step            <- HDclassif:::hddc_m_step
hddc_main              <- HDclassif:::hddc_main
matchTypeAndSetDefault<- HDclassif:::matchTypeAndSetDefault
myAlerts               <- HDclassif:::myAlerts
myCallAlerts           <- HDclassif:::myCallAlerts
plot.hdc               <- HDclassif:::plot.hdc
predict.hdc            <- HDclassif:::predict.hdc
predict.hdmda          <- HDclassif:::predict.hdmda
print.hd               <- HDclassif:::print.hd
print.hdc              <- HDclassif:::print.hdc
setHDclassif.show      <- HDclassif:::setHDclassif.show
simuldata              <- HDclassif:::simuldata
slopeHeuristic         <- HDclassif:::slopeHeuristic

#Changed functions
hdda_prior<-function (data, cls, model = "AkjBkQkDk", graph = FALSE, d_select = "Cattell", 
                 prior=NULL,threshold = 0.2, com_dim = NULL, show = getHDclassif.show(), 
                 scaling = FALSE, cv.dim = 1:10, cv.threshold = c(0.001, 0.005, 
                                                                  0.05, 1:9 * 0.1), cv.vfold = 10, LOO = FALSE, noise.ctrl = 1e-08, 
                 d) 
{
  cls <- as.factor(cls)
  names <- levels(cls)
  if(!is.null(prior) & (sum(prior)!=1) |length(prior)!=length(names) | length(prior[prior<0])!=0){
    print("prior parameter should be a vector of probabilities, prior ignored")
    prior<-as.null(prior)
  }
  if (!missing(d) & missing(d_select)) 
      d_select = d
  call = match.call()
  HDclassif:::hdda_control(call)
  d_select = HDclassif:::myAlerts(d_select, "d_select", "singleCharacterMatch.arg", 
                      "HDDA: ", c("cattell", "bic", "cv"))
  model = HDclassif:::hdc_getTheModel(model)
  ModelNames <- c("AKJBKQKDK", "AKBKQKDK", "ABKQKDK", "AKJBQKDK", 
                  "AKBQKDK", "ABQKDK", "AKJBKQKD", "AKBKQKD", "ABKQKD", 
                  "AKJBQKD", "AKBQKD", "ABQKD", "AJBQD", "ABQD", "ALL")
  Mod2 <- c("AKJBKQKDK", "AKBKQKDK ", "ABKQKDK ", "AKJBQKDK ", 
            "AKBQKDK ", "ABQKDK  ", "AKJBKQKD ", "AKBKQKD ", "ABKQKD  ", 
            "AKJBQKD ", "AKBQKD  ", "ABQKD  ", "AJBQD  ", "ABQD   ")
  names <- levels(cls)
  z <- unclass(cls)
  data <- as.matrix(data)
  K <- max(z)
  if (scaling) {
    data <- scale(data)
    scaling <- list(mu = attr(data, "scaled:center"), sd = attr(data, 
                                                                "scaled:scale"))
  }
  else scaling <- NULL
  N <- nrow(data)
  p <- ncol(data)
  n <- table(z)
  if (LOO) {
    if (model %in% ModelNames[1:14] && !is.null(com_dim)) 
      d_select <- com_dim
    class <- rep(NA, N)
    posterior <- matrix(NA, N, K)
    for (i in 1:N) {
      prms <- NULL
      try(prms <- hdda_prms_prior(data[-i, ], z[-i], model, threshold, 
                            d_select, names, noise.ctrl), silent = TRUE,prior=prior)
      if (!is.null(prms)) {
        res <- NULL
        try(res <- predict.hdc(prms, data[i, ]), silent = TRUE)
        if (!is.null(res)) {
          class[i] <- res$class
          posterior[i, ] <- res$posterior
        }
      }
    }
    class <- factor(class, labels = names, levels = seq_along(names))
    return(list(class = class, posterior = posterior))
  }
  if (d_select == "cv") {
    d_select = "cattell"
    d.max <- if (model %in% ModelNames[7:12]) 
      min(n, p) - 1
    else min(N, p) - 1
    cv.dim <- sort(cv.dim, decreasing = TRUE)
    cv.dim <- cv.dim[cv.dim <= d.max]
    if (length(cv.dim) == 0) 
      stop("cv.dim must be an integer stricly inferior to the dimension.", 
           call. = FALSE)
    cv.threshold <- sort(cv.threshold)
    cv.threshold <- cv.threshold[cv.threshold >= 0 & cv.threshold <= 
                                   1]
    if (length(cv.threshold) == 0) 
      stop("cv.threshold must be a float within 0 and 1.\n", 
           call. = FALSE)
    cv.vfold <- if (cv.vfold < N) 
      cv.vfold
    else N
    u <- sample(1:N)
    ind <- c()
    for (i in 1:cv.vfold) ind[i] <- if (i == 1) 
      floor(N/cv.vfold)
    else floor((N - sum(ind))/(cv.vfold + 1 - i))
    fin <- cumsum(ind)
    debut <- c(1, fin[-cv.vfold] + 1)
    if (model %in% ModelNames[7:14]) {
      n_cv <- length(cv.dim)
      cv.threshold <- rep(0.5, n_cv)
    }
    else {
      n_cv <- length(cv.threshold)
      cv.dim <- rep("cattell", n_cv)
    }
    res <- fails <- rep(0, n_cv)
    N2 <- rep(N, n_cv)
    for (j in 1:cv.vfold) {
      ind <- u[debut[j]:fin[j]]
      prms <- NULL
      i <- 0
      while ((i <- i + 1) <= n_cv && is.null(prms)) {
        try(prms <- hdda_prms_prior(data[-ind, ], z[-ind], 
                              model, cv.threshold[i], cv.dim[i], names, noise.ctrl,prior=prior), 
            silent = TRUE)
        if (!is.null(prms)) {
          try(res[i] <- res[i] + sum(predict.hdc(prms, 
                                                 data[ind, ])$class == cls[ind]), silent = TRUE)
        }
        else {
          N2[i] <- N2[i] - length(ind)
          fails[i] <- fails[i] + 1
        }
      }
      if (i <= n_cv) 
        for (i in i:n_cv) {
          if (model %in% ModelNames[1:6]) {
            d <- HDclassif:::hdclassif_dim_choice(prms$ev, as.vector(table(z[-ind])), 
                                      "cattell", cv.threshold[i], FALSE, noise.ctrl)
          }
          else d <- rep(cv.dim[i], K)
          if (model %in% ModelNames[13:14]) {
            prms$Q <- prms$Q[, 1:d[1]]
          }
          else {
            for (ii in 1:K) {
              if (prms$d[ii] > 1) {
                prms$Q[[ii]] <- prms$Q[[ii]][, 1:d[ii]]
              }
            }
          }
          prms$d <- d
          prms_bis <- hdda_prms_bis(model, prms, p)
          try(res[i] <- res[i] + sum(predict.hdc(prms_bis, 
                                                 data[ind, ])$class == cls[ind]), silent = TRUE)
        }
    }
    if (show) {
      if (model %in% ModelNames[7:14]) 
        cat("\t Model  \t dim\t CV\n")
      else cat("\t Model  \tthreshold\t CV\n")
      for (i in n_cv:1) {
        if (model %in% ModelNames[7:14]) 
          cat("\t", Mod2[model == ModelNames], "\t", 
              cv.dim[i], "\t", res[i]/N2[i] * 100, if (fails[i] > 
                                                       0) 
                paste(" Info: failed", fails, "times"), 
              "\n")
        else cat("\t", Mod2[model == ModelNames], "\t", 
                 cv.threshold[i], "\t\t", res[i]/N2[i] * 100, 
                 if (fails[i] > 0) 
                   paste(" Info: failed", fails, "times"), "\n")
      }
    }
    res <- res/N2 * 100
    res <- res[n_cv:1]
    cv.dim <- cv.dim[n_cv:1]
    cv.threshold <- cv.threshold[n_cv:1]
    if (model %in% ModelNames[7:14]) {
      d <- com_dim <- cv.dim[which.max(res)]
      if (show) 
        cat("Best dimension with respect to the CV results: ", 
            d, ".\n", sep = "")
      if (graph) {
        barplot(res - 100/K, names.arg = cv.dim, offset = 100/K, 
                col = "blue", xlab = "Dimensions", ylab = "Correct classification rate", 
                axes = FALSE, main = paste("Cross-Validation\n(chosen dim=", 
                                           d, ")", sep = ""))
        axis(2, at = floor(100/K + (max(res) - 100/K)/5 * 
                             0:5))
      }
      d_select = d
    }
    else {
      d_select <- "cattell"
      threshold <- cv.threshold[which.max(res)]
      if (show) 
        cat("Best threshold with respect to the CV results: ", 
            threshold, ".\n", sep = "")
      if (graph) {
        barplot(res - 100/K, names.arg = cv.threshold, 
                offset = 100/K, col = "blue", xlab = "Thresholds", 
                ylab = "Correct classification rate", axes = FALSE, 
                main = paste("Cross-Validation\nthreshold=", 
                             threshold, sep = ""))
        axis(2, at = floor(100/K + (max(res) - 100/K)/5 * 
                             0:5))
      }
    }
  }
  if (length(model) > 1) {
    nm = length(model)
    e <- vector(mode = "list", length = nm)
    BIC <- ICL <- c()
    for (i in 1:nm) {
      e[[i]] <- hdda_prms_prior(data, z, model[i], threshold, 
                          d_select, names, noise.ctrl, com_dim,prior=prior)
      BIC[i] <- HDclassif:::hdclassif_bic(e[[i]], p)$bic
      ICL[i] <- HDclassif:::hdclassif_bic(e[[i]], p)$icl
    }
    prms <- e[[which.max(BIC)]]
    prms$BIC <- max(BIC, na.rm = TRUE)
    prms$scaling <- scaling
    prms$threshold <- threshold
    if (show) {
      cat(" # :\t Model \t   BIC\n")
      for (i in 1:nm) {
        if (i < 10) 
          cat(" ")
        wng <- if (any(e[[i]]$b < 1e-05) | any(e[[i]]$a < 
                                               1e-05, na.rm = TRUE)) 
          "info: b < 10e-6"
        else ""
        cat(i, ":\t", Mod2[ModelNames == model[i]], "\t", 
            BIC[i], wng, "\n")
      }
      cat("\nSELECTED: Model ", prms$model, ", BIC=", prms$BIC, 
          ".\n", sep = "")
    }
    if (graph) {
      BIC <- BIC[!is.na(BIC)]
      min_b = min(BIC[BIC != -Inf])
      max_b = max(BIC, na.rm = TRUE)
      BIC[BIC == -Inf] <- min_b
      barplot(BIC - min_b, names.arg = model, offset = min_b, 
              col = "blue", xlab = "models", ylab = "BIC", 
              axes = FALSE, main = paste("BIC for all models\n(chosen model=", 
                                         prms$model, ")", sep = ""))
      axis(2, at = floor(min_b + (max_b - min_b)/5 * 0:5))
    }
    class(prms) <- "hdc"
    return(prms)
  }
  else if (model == "ALL") {
    e <- vector(mode = "list", length = 14)
    BIC <- ICL <- c()
    e[[1]] <- hdda_prms_prior(data, z, ModelNames[1], threshold, 
                        d_select, names, noise.ctrl,prior=prior)
    for (i in 2:6) e[[i]] <- hdda_prms_bis(ModelNames[i], 
                                           e[[1]], p)
    e[[7]] <- hdda_prms_prior(data, z, ModelNames[7], threshold, 
                        d_select, names, noise.ctrl, com_dim,prior=prior)
    for (i in 8:12) e[[i]] <- hdda_prms_bis(ModelNames[i], 
                                            e[[7]], p)
    e[[13]] <- hdda_prms_prior(data, z, ModelNames[13], threshold, 
                         d_select, names, noise.ctrl, com_dim,prior=prior)
    e[[14]] <- hdda_prms_bis(ModelNames[14], e[[13]], p)
    for (i in 1:14) {
      BIC[i] <- HDclassif:::hdclassif_bic(e[[i]], p)$bic
      ICL[i] <- HDclassif:::hdclassif_bic(e[[i]], p)$icl
    }
    prms <- e[[which.max(BIC)]]
    prms$BIC <- max(BIC, na.rm = TRUE)
    prms$scaling <- scaling
    prms$threshold <- threshold
    if (show) {
      cat(" # :\t Model \t   BIC\n")
      for (i in 1:14) {
        if (i < 10) 
          cat(" ")
        wng <- if (any(e[[i]]$b < 1e-05) | any(e[[i]]$a < 
                                               1e-05, na.rm = TRUE)) 
          "info: b < 10e-6"
        else ""
        cat(i, ":\t", Mod2[i], "\t", BIC[i], wng, "\n")
      }
      cat("\nSELECTED: Model ", prms$model, ", BIC=", prms$BIC, 
          ".\n", sep = "")
    }
    if (graph) {
      min_b <- min(BIC[BIC != -Inf])
      max_b <- max(BIC)
      BIC[BIC == -Inf] <- min_b
      barplot(BIC - min_b, names.arg = 1:14, offset = min_b, 
              col = "blue", xlab = "models", ylab = "BIC", 
              axes = FALSE, main = paste("BIC for all models\n(chosen model=", 
                                         prms$model, ")", sep = ""))
      axis(2, at = floor(min_b + (max_b - min_b)/5 * 0:5))
    }
    class(prms) <- "hdc"
    return(prms)
  }
  else {
    prms <- hdda_prms_prior(data, z, model, threshold, d_select, 
                      names, noise.ctrl, com_dim,prior=prior)
    prms$BIC <- HDclassif:::hdclassif_bic(prms, p)$bic
    prms$ICL <- HDclassif:::hdclassif_bic(prms, p)$icl
    prms$scaling <- scaling
    prms$threshold <- threshold
    class(prms) <- "hdc"
    return(prms)
  }
}
hdda_prms_prior<-function (data, cls, model, threshold, method, kname, noise.ctrl, 
          com_dim = NULL,prior=NULL) 
{
  p <- ncol(data)
  N <- nrow(data)
  K <- max(cls)
  com_ev <- NULL
  info <- NULL
  n <- as.vector(table(cls))
  if(is.null(prior)){
  prop <- matrix(n/N, 1, K, dimnames = list(c(""), `Prior probabilities of groups:` = kname))
  }else{
    prop <- matrix(prior, 1, K, dimnames = list(c(""), `Prior probabilities of groups:` = kname))
  }
  mu <- matrix(rowsum(data, cls)/n, K, p, dimnames = list(Class = kname, 
                                                          `Group means:` = paste("V", 1:p, sep = "")))
  if (model %in% c("AKJBKQKD", "AKBKQKD", "ABKQKD", "AKJBQKD", 
                   "AKBQKD", "ABQKD", "AJBQD", "ABQD")) {
    if (N < p) {
      Y <- matrix(0, N, p)
      for (i in 1:K) {
        qui = which(cls == i)
        Y[qui, ] <- (data[qui, ] - matrix(mu[i, ], length(qui), 
                                          p, byrow = TRUE))/sqrt(N)
      }
      if (model %in% c("AJBQD", "ABQD")) {
        donnees <- eigen(tcrossprod(Y), symmetric = TRUE)
      }
      else {
        donnees <- eigen(tcrossprod(Y), symmetric = TRUE, 
                         only.values = TRUE)
      }
    }
    else {
      W <- matrix(0, p, p)
      for (i in 1:K) {
        W <- W + prop[i] * crossprod(data[which(cls == 
                                                  i), ] - matrix(mu[i, ], sum(cls == i), p, byrow = TRUE))/n[i]
      }
      if (model %in% c("AJBQD", "ABQD")) {
        donnees <- eigen(W, symmetric = TRUE)
      }
      else {
        donnees <- eigen(W, symmetric = TRUE, only.values = TRUE)
      }
    }
    ev <- com_ev <- donnees$values
  }
  if (!model %in% c("AJBQD", "ABQD")) {
    if (any(n < p)) {
      Y <- vector(mode = "list", length = K)
    }
    Q <- vector(mode = "list", length = K)
    ev <- matrix(NA, K, min(max(n), p))
    for (i in which(n < p)) {
      Y[[i]] <- (data[which(cls == i), ] - matrix(mu[i, 
      ], sum(cls == i), p, byrow = TRUE))/sqrt(n[i])
      donnees <- eigen(tcrossprod(Y[[i]]), symmetric = TRUE)
      ev[i, 1:n[i]] <- donnees$values
      Q[[i]] <- donnees$vectors
    }
    for (i in which(n >= p)) {
      donnees <- eigen(crossprod(data[which(cls == i), 
      ] - matrix(mu[i, ], sum(cls == i), p, byrow = TRUE))/n[i], 
      symmetric = TRUE)
      ev[i, ] <- donnees$values
      Q[[i]] <- donnees$vectors
    }
  }
  if (model %in% c("AJBQD", "ABQD")) {
    if (!is.null(com_dim)) 
      method <- com_dim
    if (method %in% c("cattell", "bic")) 
      method <- HDclassif:::hdclassif_dim_choice(com_ev, n, method, 
                                     threshold, FALSE, noise.ctrl)
    d <- rep(method, K)
  }
  else if (model %in% c("AKJBKQKD", "AKBKQKD", "ABKQKD", "AKJBQKD", 
                        "AKBQKD", "ABQKD")) {
    if (!is.null(com_dim)) 
      method <- com_dim
    if (method %in% c("cattell", "bic")) 
      method <- HDclassif:::hdclassif_dim_choice(com_ev, n, method, 
                                     threshold, FALSE, noise.ctrl)
    d <- rep(method, K)
    if (d[1] > min(n, p) - 1) {
      d[] <- min(n, p) - 1
      info <- paste("Information: d has been lowered to", 
                    d[1], "because of the class", kname[which.min(n)], 
                    "which has", min(n), "observations.")
    }
    dmax <- if (any(ev < noise.ctrl, na.rm = TRUE)) 
      max(min(unlist(apply(ev < noise.ctrl, 1, which))) - 
            2, 1)
    else Inf
    if (d[1] > dmax) 
      d[] <- dmax
  }
  else {
    d <- HDclassif:::hdclassif_dim_choice(ev, n, method, threshold, FALSE, 
                              noise.ctrl)
  }
  if (model %in% c("AJBQD", "ABQD")) {
    if (N >= p) {
      Q <- matrix(donnees$vectors[, 1:d[1]], p, d[1])
    }
    else {
      Q <- matrix(t(Y) %*% donnees$vectors[, 1:d[1]], p, 
                  d[1])
      normalise <- c()
      for (i in 1:d[1]) normalise[i] <- as.double(crossprod(Q[, 
                                                              i]))
      Q <- Q/matrix(sqrt(normalise), p, d[1], byrow = TRUE)
    }
  }
  else {
    for (i in which(n >= p)) {
      Q[[i]] <- matrix(Q[[i]][, 1:d[i]], p, d[i])
    }
    for (i in which(n < p)) {
      Q[[i]] <- t(Y[[i]]) %*% (Q[[i]][, 1:d[i]])
      normalise <- c()
      for (j in 1:d[i]) normalise[j] <- as.double(crossprod(as.matrix(Q[[i]][, 
                                                                             j])))
      Q[[i]] <- Q[[i]]/matrix(sqrt(normalise), p, d[i], 
                              byrow = TRUE)
    }
  }
  if (model %in% c("AKJBKQKDK", "AKJBQKDK", "AKJBKQKD", "AKJBQKD")) {
    ai <- matrix(NA, K, max(d), dimnames = list(Class = kname, 
                                                `Akj:` = paste("a", 1:max(d), sep = "")))
    for (i in 1:K) ai[i, 1:d[i]] <- ev[i, 1:d[i]]
  }
  else if (model %in% c("AKBKQKDK", "AKBQKDK", "AKBKQKD", "AKBQKD")) {
    ai <- matrix(NA, 1, K, dimnames = list(c("Ak:"), kname))
    for (i in 1:K) ai[i] <- sum(ev[i, 1:d[i]])/d[i]
  }
  else if (model == "AJBQD") {
    ai <- matrix(ev[1:d[1]], 1, d[1], dimnames = list(c("Aj:"), 
                                                      paste("a", 1:d[1], sep = "")))
  }
  else if (model == "ABQD") {
    ai <- matrix(sum(ev[1:d[1]])/d[1], dimnames = list(c("A:"), 
                                                       c("")))
  }
  else {
    a <- 0
    eps <- sum(prop * d)
    for (i in 1:K) a <- a + sum(ev[i, 1:d[i]]) * prop[i]
    ai <- matrix(a/eps, dimnames = list(c("A:"), c("")))
  }
  if (model %in% c("AKJBKQKDK", "AKBKQKDK", "ABKQKDK", "AKJBKQKD", 
                   "AKBKQKD", "ABKQKD")) {
    bi <- matrix(NA, 1, K, dimnames = list(c("Bk:"), kname))
    for (i in which(n >= p)) bi[i] <- sum(ev[i, (d[i] + 1):p])/(p - 
                                                                  d[i])
    for (i in which(n < p)) bi[i] <- sum(ev[i, (d[i] + 1):n[i]])/(p - 
                                                                    d[i])
  }
  else if (model %in% c("ABQD", "AJBQD")) {
    if (N >= p) 
      bi <- matrix(sum(ev[(d[1] + 1):p])/(p - d[1]), dimnames = list(c("B:"), 
                                                                     c("")))
    else bi <- matrix(sum(ev[(d[1] + 1):N])/(N - d[1]), dimnames = list(c("B:"), 
                                                                        c("")))
  }
  else {
    b <- 0
    eps <- sum(prop * d)
    for (i in which(n >= p)) b <- b + sum(ev[i, (d[i] + 1):p]) * 
      prop[i]
    for (i in which(n < p)) b <- b + sum(ev[i, (d[i] + 1):n[i]]) * 
      prop[i]
    bi <- matrix(b/(min(max(n), p) - eps), dimnames = list(c("B:"), 
                                                           c("")))
  }
  d <- matrix(d, 1, K, dimnames = list(c("dim:"), `Intrinsic dimensions of the classes:` = kname))
  class(prop) <- class(mu) <- class(ai) <- class(bi) <- class(d) <- class(ev) <- "hd"
  list(model = model, K = K, d = d, a = ai, b = bi, mu = mu, 
       prop = prop, ev = ev, Q = Q, kname = kname, info = info, 
       N = N, com_ev = com_ev)
}

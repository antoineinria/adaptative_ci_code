load(file = "data/params.RData")
source("functions_ci.R")
load("data/modeles.RData")
alpha<-.05
T<-n*M
nb<-length(rhos_train)
distr<-"student"

#neural network:
library(keras)
library(tensorflow)
library(reticulate)
nn <- keras_model_sequential() %>%
  layer_dense(units = 256, activation = "relu", input_shape = c(ncol(Sets_burr))) %>%
  layer_dropout(rate = 0.3) %>%
  layer_dense(units = 256, activation = "relu") %>%
  layer_dropout(rate = 0.3) %>%
  layer_dense(units = 64, activation = "relu") %>%
  layer_dense(units = 3)  # sortie logits
nn %>% compile(
  optimizer = "adam",
  loss = loss_sparse_categorical_crossentropy(from_logits = TRUE),
  metrics = "accuracy"
)
class_weights <- list("0" = 1, "1" = 1, "2" = 5)

history <- nn %>% fit(
  x = Sets_burr,
  y = as.numeric(ys_burr)-1,#classes: 0: NCI1, 1: NCI3, 2:NCI5
  epochs = 20,
  batch_size = 64,
  class_weight = class_weights
)
prob_nn <- keras_model_sequential() %>%
  nn %>%
  layer_activation("softmax")



# function to compute the algorithm:
quantiles<-function(j){
  rho<-rhos_train[j] 
  G<-gammas_train[j]
  seed<-seeds_res[j]
  plda1<-algopercentile_ci(n,G,rho,p,M,kmax,distr,
                        seed, interval = predict_ci,alpha=alpha,
                        model=lda_model)
  plda2<-algopercentile_ci(n,G,rho,p,M,kmax,distr,
                        seed, interval = predict_ci,alpha=alpha,
                        model=lda_prior)
  psvm<-algopercentiles_opt(svm_prior,n,G,rho,p,M,kmax=kmax,distr=distr,seed=seed)
  phdda<-algopercentile_ci(n,G,rho,p,M,kmax,distr,
                        seed, interval = predict_ci,alpha=alpha,
                        model=hdda_prior)
  pknn<-algopercentiles_opt(knn_model,n,G,rho,p,M,kmax=kmax,distr=distr,alpha=alpha,seed=seed)
  pnn<-algopercentilesnn_opt(prob_nn,n,G,rho,p,M,kmax=kmax,distr=distr,seed=seed)
  return(list(
    "j"=j,
    "plda1"=plda1$percentile,
    "plda2"=plda2$percentile,
    "psvm"=psvm$percentile,
    "phdda"=phdda$percentile,
    "pknn"=pknn$percentile,
    "pnn"=pnn$percentile
    ))
}

quantiles(1)
pslda1_stud<-numeric(nb)
pslda2_stud<-numeric(nb)
pssvm_stud<-numeric(nb)
pshdda_stud<-numeric(nb)
psknn_stud<-numeric(nb)
psnn_stud<-numeric(nb)

for (j in 1:nb) {
  item<-quantiles(j)
  pslda1_stud[item$j] <- item$plda1
  pslda2_stud[item$j] <- item$plda2
  pssvm_stud[item$j] <- item$psvm
  pshdda_stud[item$j] <- item$phdda
  psknn_stud[item$j] <- item$pknn
  psnn_stud[item$j] <- item$pnn
}

# save(pslda1_stud,pslda2_stud,pssvm_stud,pshdda_stud,psknn_stud,psnn_stud,
# file="data/results/res_modeles_stud.RData"
# )


source("functions_ci.R")
source("functions_hdda.R") #for hdda with prior
load("data/train/percentiles_burr.RData")
load("data/params.RData")


#training set based on Algorithm 1 for Burr distributions
ys_burr<-factor(which_vector_min(psn1_burr,psn3_burr,psn5_burr))
table(ys_burr)
Sets_burr<-training_set(n,gammas_train,rhos_train,kmax,"burr",
                        seeds=seeds_trainset)$Sets_train
dim(Sets_burr)

X_burr<-data.frame(Sets_burr)
X_burr$y<-factor(ys_burr)

lda_model<-lda(y~.,data = X_burr)
lda_prior<-lda(y~.,data = X_burr,prior=priors)
svm_prior<-svm(x=X_burr[,-301],y=X_burr$y,kernel="linear",class.weights=svm_weights)
knn_model<-train.kknn(y ~ ., data = X_burr, kmax = 20, kernel = "optimal", scale = TRUE)
#svm_model<-svm(x=X_burr[,-301],y=X_burr$y)
hdda_prior<-hdda_prior(data = X_burr[,-301],d_select="BIC",
                      prior=c(1/4,1/4,1/2),cls = ys_burr)

hdda_prior
#svm_model
svm_prior

#Predictions:

#predictions on training set:
table(ys_burr,predict(lda_model)$class)
table(ys_burr,predict(lda_prior)$class)
#table(ys_burr,predict(svm_model))
table(ys_burr,predict(svm_prior))
table(ys_burr,predict(knn_model,newdata = X_burr))
table(ys_burr,predict(hdda_prior,data = X_burr[,-301])$class)

Sets_invg<-training_set(n,gammas_train,rhos_train,kmax,"invgamma")$Sets_train
load("~data/train/percentiles_invg.RData")
ys_invg<-which_vector_min(psn1_invg,psn3_invg,psn5_invg)
X_invg<-data.frame(Sets_invg)

#prediction on test set with inverse-gamma distribution
table(ys_invg,predict(lda_model,X_invg)$class)
table(ys_invg,predict(lda_prior,X_invg)$class)
#table(ys_invg,predict(svm_model,X_invg))
table(ys_invg,predict(svm_prior,X_invg))
table(ys_invg,predict(kknn_model,newdata = X_invg ))
table(ys_invg,predict(hdda_prior,X_invg)$class)

Sets_fish<-training_set(n,gammas_train,rhos_train,kmax,"fisher")$Sets_train
load("data/train/percentiles_fish.RData")
ys_fish<-which_vector_min(psn1_fish,psn3_fish,psn5_fish)
X_fish<-data.frame(Sets_fish)
table(ys_fish,predict(lda_model,X_fish)$class)
#table(ys_fish,predict(svm_model,X_fish))
table(ys_fish,predict(lda_prior,X_fish)$class)
table(ys_fish,predict(kknn_model,newdata = X_fish ))
table(ys_fish,predict(svm_prior,X_fish))
table(ys_fish,predict(hdda_prior,X_fish)$class)

Sets_stud<-training_set(n,gammas_train,rhos_train,kmax,"student")$Sets_train
load("data/train/percentiles_stud.RData")
ys_stud<-which_vector_min(psn1_stud,psn3_stud,psn5_stud)
X_stud<-data.frame(Sets_stud)
table(ys_stud,predict(lda_model,X_stud)$class)
table(ys_stud,predict(lda_prior,X_stud)$class)
#table(ys_stud,predict(svm_model,X_stud))
table(ys_stud,predict(svm_prior,X_stud))
table(ys_stud,predict(kknn_model,newdata = X_stud ))
table(ys_stud,predict(hdda_prior,X_stud)$class)




save(X_burr,Sets_burr,ys_burr,
    lda_model,lda_prior,
    #svm_model,
    knn_model,svm_prior,hdda_prior,
     file="data/modeles.RData")

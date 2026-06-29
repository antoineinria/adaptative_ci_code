source("functions_ci.R")
n<-1000
kmax<-floor(3*n/10)
p<-1/2/n
N<-10000
M<-1000
priors<-c(1/4,1/4,1/2)
svm_weights <- c("1" = 0.25, "2" = 0.25, "3" = 0.5)


#rho

#Distribution of rho is a truncnorm with approx
#the same proportion of labels NCI1, NCI3, NCI5
#this optimizes the value of mean and variance
#for the truncnorm distribution of rho
inff<- -5
rho1<- -1.25
rho2<- -.6
supp<- -1/10
p1<-1/3
p3<-1/3
p2<-1-p1-p3
breaks <- c(inff, rho1, rho2, supp)
target_probs <- c(p1,p3,p2) #NCI1,NCI5,NCI3


result <- optim(c(-1, 1),
                fn = find_truncnorm_params, 
                breaks = breaks,
                probs = target_probs,
                method = "L-BFGS-B",
                lower = c(-2, .1),
                upper=c(100,2))

mu <- result$par[1]
sigma <- result$par[2]
mu
sigma
nb<-20000
seed_rho<-842
set.seed(seed_rho)
gammas<-runif(nb,1/10,1)
rhos<-rtruncnorm(nb,a=inff,b=supp,mean = mu,sd = sigma)
rhos_train<-rhos[1:10000]
gammas_train<-gammas[1:10000]

seeds_trainset<-seq(40000,49999)
seeds_res<-seq(10001,20003)
seeds_res<-seeds_res[-c(342,4140,8022)] #remove unadapted seeds
seed_nn<-42
save(n,p,N,M,priors,kmax,svm_weights,gammas_train,rhos_train,
     seed_rho,seeds_trainset,seeds_res,seed_nn,
     file="data/params.RData")
#changé; params.RData->data/params.RData
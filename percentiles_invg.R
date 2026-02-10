# library(future)
# library(future.apply)
load(file = "data/params.RData")
source("functions_ci.R")
alpha<-.05
T<-n*M
nb<-length(rhos_train)
distr<-"invgamma"
debb<-1
fin<-nb
# Fonction de traitement
quantiles<-function(j){
  rho<-rhos_train[j] 
  G<-gammas_train[j]
  seed<-seeds_res[j]
  n1<-algopercentile_ci(n,G,rho,p,M,kmax,distr,
                        seed, interval = cik_logweiss,
                        alpha=alpha,loginterval=TRUE,
                        int="NCI1")
  n2<-algopercentile_ci(n,G,rho,p,M,kmax,distr,
                        seed, interval = cik_logweiss,
                        alpha=alpha,loginterval=TRUE,
                        int="NCI2")
  n3<-algopercentile_ci(n,G,rho,p,M,kmax,distr,
                        seed, interval = cik_logweiss,
                        alpha=alpha,loginterval=TRUE,
                        int="NCI3")
  n4<-algopercentile_ci(n,G,rho,p,M,kmax,distr,
                        seed, interval = cik_logweiss,
                        alpha=alpha,loginterval=TRUE,
                        int="NCI4")
  n5<-algopercentile_ci(n,G,rho,p,M,kmax,distr,
                        seed, interval = cik_logweiss,
                        alpha=alpha,loginterval=TRUE,
                              int="NCI5")
  g1<-algopercentile_ci(n,G,rho,p,M,kmax,distr,
                        seed, interval = cik_logweiss,
                        alpha=alpha,loginterval=TRUE,
                        int="GCI1")
  g2<-algopercentile_ci(n,G,rho,p,M,kmax,distr,
                        seed, interval = cik_logweiss,
                        alpha=alpha,loginterval=TRUE,
                        int="GCI2")
  
  be<-algopercentile_ci(n,G,rho,p,M,kmax,distr,
                        seed, interval = cik_beta,
                        alpha=alpha,loginterval=FALSE,
                        a=.85)
  return(list(
    "j"=j,
    "pn1"=n1$percentile,
    "pn2"=n2$percentile,
    "pn3"=n3$percentile,
    "pn4"=n4$percentile,
    "pn5"=n5$percentile,
    "pg1"=g1$percentile,
    "pg2"=g2$percentile,
    "pbe"=be$percentile
  ))
}


psn1_invg<-numeric(nb)
psn2_invg<-numeric(nb)
psn3_invg<-numeric(nb)
psn4_invg<-numeric(nb)
psn5_invg<-numeric(nb)
psg1_invg<-numeric(nb)
psg2_invg<-numeric(nb)

psbe_invg<-numeric(nb)
for (j in 1:nb) {
  item<-quantiles(j)
  psn1_invg[item$j] <- item$pn1
  psn2_invg[item$j] <- item$pn2
  psn3_invg[item$j] <- item$pn3
  psn4_invg[item$j] <- item$pn4
  psn5_invg[item$j] <- item$pn5
  psg1_invg[item$j] <- item$pg1
  psg2_invg[item$j] <- item$pg2
  psbe_invg[item$j] <- item$pbe
}

save(psn1_invg,psn2_invg,psn3_invg,psn4_invg,psn5_invg,
     psg1_invg,psg2_invg,psbe_invg,
file="data/train/percentiles_invg.RData"
)
print("SAVED")


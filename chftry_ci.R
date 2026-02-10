source("functions_ci.R")
load("data/modeles.RData")
library(quantmod)
#library(xts)
#library(zoo)
# data from yahoofinance
getSymbols("CHFTRY=X", src = "yahoo", from = "2022-01-01", to = "2025-08-28")
#save(dats,file = "~/Documents/codes_r_article/codes_finaux/data/real_data/chftry.RData")

data_clean <- na.omit(Ad(`CHFTRY=X`))

log_returns <- diff(log(data_clean))[-1]
names(log_returns) <- "Log_Return"


dats<-as.numeric(coredata(log_returns))
dats<-abs(dats[dats!=0])

n<-length(dats)
n

dord<-sort(dats)

Ks<-1:300
ggplot()+
  geom_line(aes(x=Ks,y=Hill(dats)$gamma[Ks]),color="black",size=1.5)+
  geom_vline(xintercept = 175,linewidth=1.2,linetype="longdash",color="red")+
  labs(
    x = "k",
    y = NULL
  ) +
  theme_linedraw() +
  theme(
    axis.text.x = element_text(size = 30),
    axis.text.y = element_text(size = 30),
    axis.title.x = element_text(size = 30),
    axis.title.y = element_text(size = 30),
    panel.grid.major = element_line(color = "black",linewidth = .3),
    panel.grid.minor = element_line(color = "grey50",linewidth = 0.3)
  )
kmax<-175
quant_1<-(log(n/c(1:kmax)))-(log(n/(kmax+1)))
quant_data<-log(dord[n-c(1:kmax)+1]/dord[n-kmax])
dqq<-data.frame(quant_1,quant_data)
model<-lm(quant_data~quant_1, data=dqq)
coef(model)
ggplot(dqq, aes(x = quant_1, y = quant_data)) +
  geom_point(size=5) +  # Points de données
  geom_smooth(size=3,se=FALSE,method = "lm", color = "grey") + 
  theme_linedraw() +
  labs(x=NULL,y=NULL)+
  theme(
    # #legend.position = "none",
    axis.text.x = element_text(size = 30),
    axis.text.y = element_text(size = 30),
    axis.title.x = element_text(size = 30),
    axis.title.y = element_text(size = 30),
  )
#data visualisation
plot(log_returns$Log_Return,type="l")

classif_dat<-data_for_classif(dats,kmax = 300)
dat_frame<-data.frame(classif_dat)
predict(lda_prior,dat_frame)
predict(lda_model,dat_frame)
predict(knn_model,dat_frame)
predict(hdda_prior,dat_frame)
predict(svm_prior,dat_frame)

#to predict using the Neural network, it needs to be reconstructed
#predict(prob_nn,classif_dat,verbose = 0)

K<-1:300
#intervals
nci1<-cik_logweiss(dats,K=K,p=1/2/n,int="NCI1")
nci2<-cik_logweiss(dats,K=K,p=1/2/n,int="NCI2")
nci3<-cik_logweiss(dats,K=K,p=1/2/n,int="NCI3")
nci4<-cik_logweiss(dats,K=K,p=1/2/n,int="NCI4")
nci5<-cik_logweiss(dats,K=K,p=1/2/n,int="NCI5")
gci1<-cik_logweiss(dats,K=K,p=1/2/n,int="GCI1")
gci2<-cik_logweiss(dats,K=K,p=1/2/n,int="GCI2")

mm<-max(dats)
ylim<-c(1/2*mm,4*mm)
mm
echelle<-seq(0,10,by=.05)
ggplot() +
  geom_ribbon(aes(x=K,ymin = exp(nci5$lbounds), ymax= exp(nci5$ubounds)),fill="navyblue",alpha=.2)+
  geom_line(aes(x = K, y = exp(nci1$lbounds)),linewidth=2, color = "firebrick") +
  geom_line(aes(x = K, y = exp(nci1$ubounds)),linewidth=2, color = "firebrick") +
  # geom_line(aes(x = K, y = exp(nci2$lbounds)), color = "chartreuse") +
  # geom_line(aes(x = K, y = exp(nci2$ubounds)), color = "chartreuse") +
  geom_line(aes(x = K, y = exp(nci3$lbounds)),linewidth=2, color = "forestgreen") +
  geom_line(aes(x = K, y = exp(nci3$ubounds)),linewidth=2, color = "forestgreen") +
  # geom_line(aes(x = K, y = exp(nci4$lbounds)), color = "cyan") +
  # geom_line(aes(x = K, y = exp(nci4$ubounds)), color = "cyan") +
  # geom_line(aes(x = K, y = exp(gci1$lbounds)), color = "red") +
  # geom_line(aes(x = K, y = exp(gci1$ubounds)), color = "red") +
  # geom_line(aes(x = K, y = exp(gci2$lbounds)), color = "springgreen3") +
  # geom_line(aes(x = K, y = exp(gci2$ubounds)), color = "springgreen3") +
  geom_line(aes(x = K, y = exp(nci5$lbounds)),linewidth=3, color = "navyblue") +
  geom_line(aes(x = K, y = exp(nci5$ubounds)),linewidth=3, color = "navyblue") +
  scale_y_continuous(limits = ylim, breaks = echelle,
                     sec.axis = sec_axis(~.,breaks=echelle, name = NULL)
  ) +
  labs(
    x = "k",
    y = "Confidence interval"
  ) +
  theme_linedraw() +
  theme(
    # #legend.position = "none",
    axis.text.x = element_text(size = 30),
    axis.text.y = element_text(size = 30),
    axis.title.x = element_text(size = 30),
    axis.title.y = element_text(size = 30),
    panel.grid.major = element_line(color = "black",linewidth = .3),
    panel.grid.minor = element_line(color = "grey50",linewidth = 0.3)
  )
exp(nci5$lbounds)
exp(nci5$ubounds)

table(exp(nci1$ubounds)-exp(nci1$lbounds)>
        exp(nci5$ubounds)-exp(nci5$lbounds))
nci1$ubounds
nci5$ubounds


exp(nci1$ubounds)-exp(nci1$lbounds)>
  exp(nci5$ubounds)-exp(nci5$lbounds)

exp(nci5$ubounds)-exp(nci5$lbounds)

mean(dats)
max(dats)
cso<-mop(dats,c(1),c(1),"RBMOP")
#second order parameters
cso$rho
cso$beta
rho<-cso$rho
rho
p<-1/2/n
lambda<-sqrt(K)*cso$beta*(n/K)^rho
lambda
d<-K/(n*p)

source("libraries.R")
load("data/params.RData")
load("data/results/res_modeles_burr.RData")

ylim<-c(0,35)
xlim<-c(-3,0)
echx<-seq(-3,0,by=.5)
echelle<-seq(0,100,by=2.5)
sm_vec <- function(x, m=10) {
  stats::filter(x, rep(1/m, m), sides = 2)
}

db=data.frame(x=rhos_train,
              pl1=pslda1_burr,pl2=pslda2_burr,
              phdda=pshdda_burr,
             psvm=pssvm_burr,
             pkn=psknn_burr,
              pnn=psnn_burr
)



dib<-as.data.frame(db[rev(order(db[["x"]])), ])
dbs=data.frame(x=dib$x,
               pl1=sm_vec(dib$pl1),pl2=sm_vec(dib$pl2),
               phdda=sm_vec(dib$phdda),psvm=sm_vec(dib$psvm),
               pkn=sm_vec(dib$pkn),
               pnn=sm_vec(dib$pnn)
)

graph<-ggplot(dbs) +
  geom_line(aes(x = x, y = pl1, color = "LDA1")) +
  geom_line(aes(x = x, y = phdda, color = "HDDA")) +
  geom_line(aes(x = x, y = pnn, color = "NN")) +
  geom_line(aes(x = x, y = pkn, color = "KNN")) +
  geom_line(aes(x = x, y = psvm, color = "SVM")) +
  geom_line(aes(x = x, y = pl2, color = "LDA2")) +

  scale_color_manual(
    name = "",
    values = c(
      "NCI1" = "firebrick",
      "NCI2" = "chartreuse",
      "NCI3" = "darkgreen",
      "NCI4" = "cyan",
      "NCI5" = "navyblue",
      "GCI1" = "red",
      "GCI2" = "springgreen3",
      "LDA1" = "blue",
      "SVM"="purple",
      "LDA2" = "cyan",
      "KNN" = "red",
      "NN"="orange",
      "HDDA" = "forestgreen"
    )
  ) +
  scale_y_continuous(limits = ylim, breaks = echelle,
                     sec.axis = sec_axis(~.,breaks=echelle, name = NULL)
  ) +
  scale_x_continuous(limits=xlim, breaks=echx)+
  labs(
    x = "rho",
    y = NULL
  ) +
  theme_linedraw() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 30),
    axis.text.y = element_text(size = 30),
    axis.title.x = element_text(size = 30),
    axis.title.y = element_text(size = 30),
    panel.grid.major = element_line(color = "black",linewidth = .3),
    panel.grid.minor = element_line(color = "grey50",linewidth = 0.3)
  )
graph



#invgamma
source("libraries.R")
load("data/params.RData")
load("data/results/res_modeles_invg.RData")

ylim<-c(0,35)
xlim<-c(-3,0)
echx<-seq(-3,0,by=.5)
echelle<-seq(0,100,by=2.5)
sm_vec <- function(x, m=10) {
  stats::filter(x, rep(1/m, m), sides = 2)
}

db=data.frame(x=rhos_train,
              pl1=pslda1_invg,pl2=pslda2_invg,
              phdda=pshdda_invg,
              psvm=pssvm_invg,
              pkn=psknn_invg,
              pnn=psnn_invg
)



dib<-as.data.frame(db[rev(order(db[["x"]])), ])
dbs=data.frame(x=dib$x,
               pl1=sm_vec(dib$pl1),pl2=sm_vec(dib$pl2),
               phdda=sm_vec(dib$phdda),psvm=sm_vec(dib$psvm),
               pkn=sm_vec(dib$pkn),
               pnn=sm_vec(dib$pnn)
)

graph<-ggplot(dbs) +
  geom_line(aes(x = x, y = pl1, color = "LDA1")) +
  geom_line(aes(x = x, y = phdda, color = "HDDA")) +
  geom_line(aes(x = x, y = pnn, color = "NN")) +
  geom_line(aes(x = x, y = pkn, color = "KNN")) +
  geom_line(aes(x = x, y = psvm, color = "SVM")) +
  geom_line(aes(x = x, y = pl2, color = "LDA2")) +

  scale_color_manual(
    name = "",
    values = c(
      "NCI1" = "firebrick",
      "NCI2" = "chartreuse",
      "NCI3" = "darkgreen",
      "NCI4" = "cyan",
      "NCI5" = "navyblue",
      "GCI1" = "red",
      "GCI2" = "springgreen3",
      "LDA1" = "blue",
      "SVM"="purple",
      "LDA2" = "cyan",
      "KNN" = "red",
      "NN"="orange",
      "HDDA" = "forestgreen"
    )
  ) +
  scale_y_continuous(limits = ylim, breaks = echelle,
                     sec.axis = sec_axis(~.,breaks=echelle, name = NULL)
  ) +
  scale_x_continuous(limits=xlim, breaks=echx)+
  labs(
    x = "rho",
    y = NULL
  ) +
  theme_linedraw() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 30),
    axis.text.y = element_text(size = 30),
    axis.title.x = element_text(size = 30),
    axis.title.y = element_text(size = 30),
    panel.grid.major = element_line(color = "black",linewidth = .3),
    panel.grid.minor = element_line(color = "grey50",linewidth = 0.3)
  )
graph


#fisher
source("libraries.R")
load("data/params.RData")
load("data/results/res_modeles_fish.RData")

ylim<-c(0,35)
xlim<-c(-3,0)
echx<-seq(-3,0,by=.5)
echelle<-seq(0,100,by=2.5)
sm_vec <- function(x, m=20) {
  stats::filter(x, rep(1/m, m), sides = 2)
}

db=data.frame(x=rhos_train,
              # qn1=qn1s_fish,qn3=qn3s_fish,qn5=qn5s_fish, 
              pl1=pslda1_fish,pl2=pslda2_fish,
              phdda=pshdda_fish,
              psvm=pssvm_fish,
              pkn=psknn_fish,
              pnn=psnn_fish
)



dib<-as.data.frame(db[rev(order(db[["x"]])), ])
dbs=data.frame(x=dib$x,
               pl1=sm_vec(dib$pl1),pl2=sm_vec(dib$pl2),
               phdda=sm_vec(dib$phdda),psvm=sm_vec(dib$psvm),
               pkn=sm_vec(dib$pkn),
               pnn=sm_vec(dib$pnn)
)

graph<-ggplot(dbs) +
  geom_line(aes(x = x, y = pl1, color = "LDA1")) +
  geom_line(aes(x = x, y = phdda, color = "HDDA")) +
  geom_line(aes(x = x, y = pnn, color = "NN")) +
  geom_line(aes(x = x, y = pkn, color = "KNN")) +
  geom_line(aes(x = x, y = psvm, color = "SVM")) +
  geom_line(aes(x = x, y = pl2, color = "LDA2")) +

  scale_color_manual(
    name = "",
    values = c(
      "NCI1" = "firebrick",
      "NCI2" = "chartreuse",
      "NCI3" = "darkgreen",
      "NCI4" = "cyan",
      "NCI5" = "navyblue",
      "GCI1" = "red",
      "GCI2" = "springgreen3",
      "LDA1" = "blue",
      "SVM"="purple",
      "LDA2" = "cyan",
      "KNN" = "red",
      "NN"="orange",
      "HDDA" = "forestgreen"
    )
  ) +
  scale_y_continuous(limits = ylim, breaks = echelle,
                     sec.axis = sec_axis(~.,breaks=echelle, name = NULL)
  ) +
  scale_x_continuous(limits=xlim, breaks=echx)+
  labs(
    x = "rho",
    y = NULL
  ) +
  theme_linedraw() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 30),
    axis.text.y = element_text(size = 30),
    axis.title.x = element_text(size = 30),
    axis.title.y = element_text(size = 30),
    panel.grid.major = element_line(color = "black",linewidth = .3),
    panel.grid.minor = element_line(color = "grey50",linewidth = 0.3)
  )
graph


#student
source("libraries.R")
load("data/params.RData")
load("data/results/res_modeles_stud.RData")

ylim<-c(0,35)
xlim<-c(-3,0)
echx<-seq(-3,0,by=.5)
echelle<-seq(0,100,by=2.5)
sm_vec <- function(x, m=10) {
  stats::filter(x, rep(1/m, m), sides = 2)
}

db=data.frame(x=rhos_train,
              pl1=pslda1_stud,pl2=pslda2_stud,
              phdda=pshdda_stud,
              psvm=pssvm_stud,
              pkn=psknn_stud,
              pnn=psnn_stud
)



dib<-as.data.frame(db[rev(order(db[["x"]])), ])
dbs=data.frame(x=dib$x,
               pl1=sm_vec(dib$pl1),pl2=sm_vec(dib$pl2),
               phdda=sm_vec(dib$phdda),psvm=sm_vec(dib$psvm),
               pkn=sm_vec(dib$pkn),
               pnn=sm_vec(dib$pnn)
)

graph<-ggplot(dbs) +
  geom_line(aes(x = x, y = pl1, color = "LDA1")) +
  geom_line(aes(x = x, y = phdda, color = "HDDA")) +
  geom_line(aes(x = x, y = pnn, color = "NN")) +
  geom_line(aes(x = x, y = pkn, color = "KNN")) +
  geom_line(aes(x = x, y = psvm, color = "SVM")) +
  geom_line(aes(x = x, y = pl2, color = "LDA2")) +

  scale_color_manual(
    name = "",
    values = c(
      "NCI1" = "firebrick",
      "NCI2" = "chartreuse",
      "NCI3" = "darkgreen",
      "NCI4" = "cyan",
      "NCI5" = "navyblue",
      "GCI1" = "red",
      "GCI2" = "springgreen3",
      "LDA1" = "blue",
      "SVM"="purple",
      "LDA2" = "cyan",
      "KNN" = "red",
      "NN"="orange",
      "HDDA" = "forestgreen"
    )
  ) +
  scale_y_continuous(limits = ylim, breaks = echelle,
                     sec.axis = sec_axis(~.,breaks=echelle, name = NULL)
  ) +
  scale_x_continuous(limits=xlim, breaks=echx)+
  labs(
    x = "rho",
    y = NULL
  ) +
  theme_linedraw() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 30),
    axis.text.y = element_text(size = 30),
    axis.title.x = element_text(size = 30),
    axis.title.y = element_text(size = 30),
    panel.grid.major = element_line(color = "black",linewidth = .3),
    panel.grid.minor = element_line(color = "grey50",linewidth = 0.3)
  )
graph


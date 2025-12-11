source("libraries.R")
load("data/params.RData")
load("data/train/percentiles_burr.RData")
load("data/results/res_modeles_burr.RData")

ylim<-c(0,35)
xlim<-c(-3,0)
echx<-seq(-3,0,by=.5)
echelle<-seq(0,100,by=2.5)
sm_vec <- function(x, m=10) {
  stats::filter(x, rep(1/m, m), sides = 2)
}
db=data.frame(x=rhos_train,
              pn1=psn1_burr,pn3=psn3_burr,pn5=psn5_burr,
              plda2=pslda2_burr,
              pbeta=psbe_burr
)
dib<-as.data.frame(db[rev(order(db[["x"]])), ])
dbs=data.frame(x=dib$x,
               pn1=sm_vec(dib$pn1),
               pn3=sm_vec(dib$pn3),
               pn5=sm_vec(dib$pn5),
               plda2=sm_vec(dib$plda2),
               pbeta=sm_vec(dib$pbeta)
)

graph<-ggplot(dbs) +
  geom_line(aes(x = x, y = pn1, color = "NCI1")) +
  geom_line(aes(x = x, y = pn3, color = "NCI3")) +
  geom_line(aes(x = x, y = pn5, color = "NCI5")) +
  geom_line(aes(x = x, y = plda2, color = "LDA2")) +
  geom_line(aes(x = x, y = pbeta, color = "Beta")) +
  
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
      "Beta"="purple",
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
load("data/train/percentiles_invg.RData")
load("data/results/res_modeles_invg.RData")

ylim<-c(0,35)
xlim<-c(-3,0)
echx<-seq(-3,0,by=.5)
echelle<-seq(0,100,by=2.5)
sm_vec <- function(x, m=10) {
  stats::filter(x, rep(1/m, m), sides = 2)
}

db=data.frame(x=rhos_train,
              pn1=psn1_invg,pn3=psn3_invg,pn5=psn5_invg,
              plda2=pslda2_invg, pbeta=psbe_invg
)
dib<-as.data.frame(db[rev(order(db[["x"]])), ])
dbs=data.frame(x=dib$x,
               pn1=sm_vec(dib$pn1),
               pn3=sm_vec(dib$pn3),
               pn5=sm_vec(dib$pn5),
               plda2=sm_vec(dib$plda2),
               pbeta=sm_vec(dib$pbeta)
               
)

graph<-ggplot(dbs) +
  geom_line(aes(x = x, y = pn1, color = "NCI1")) +
  geom_line(aes(x = x, y = pn3, color = "NCI3")) +
  geom_line(aes(x = x, y = pn5, color = "NCI5")) +
  geom_line(aes(x = x, y = plda2, color = "LDA2")) +
  geom_line(aes(x = x, y = pbeta, color = "Beta")) +
  
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
      "Beta"="purple",
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
load("data/train/percentiles_fish.RData")
load("data/results/res_modeles_fish.RData")

ylim<-c(0,35)
xlim<-c(-3,0)
echx<-seq(-3,0,by=.5)
echelle<-seq(0,100,by=2.5)
sm_vec <- function(x, m=20) {
  stats::filter(x, rep(1/m, m), sides = 2)
}

db=data.frame(x=rhos_train,
              pn1=psn1_fish,pn3=psn3_fish,pn5=psn5_fish,
              plda2=pslda2_fish,
              pbeta=psbe_fish
)
dib<-as.data.frame(db[rev(order(db[["x"]])), ])
dbs=data.frame(x=dib$x,
               pn1=sm_vec(dib$pn1),
               pn3=sm_vec(dib$pn3),
               pn5=sm_vec(dib$pn5),
               plda2=sm_vec(dib$plda2),
               pbeta=sm_vec(dib$pbeta)
)

graph<-ggplot(dbs) +
  geom_line(aes(x = x, y = pn1, color = "NCI1")) +
  geom_line(aes(x = x, y = pn3, color = "NCI3")) +
  geom_line(aes(x = x, y = pn5, color = "NCI5")) +
  geom_line(aes(x = x, y = plda2, color = "LDA2")) +
  geom_line(aes(x = x, y = pbeta, color = "Beta")) +
  
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
      "Beta"="purple",
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
load("data/train/percentiles_stud.RData")
load("data/results/res_modeles_stud.RData")

ylim<-c(0,35)
xlim<-c(-3,0)
echx<-seq(-3,0,by=.5)
echelle<-seq(0,100,by=2.5)
sm_vec <- function(x, m=10) {
  stats::filter(x, rep(1/m, m), sides = 2)
}

db=data.frame(x=rhos_train,
              pn1=psn1_stud,pn3=psn3_stud,pn5=psn5_stud,
              plda2=pslda2_stud,
              pbeta=psbe_stud
)
dib<-as.data.frame(db[rev(order(db[["x"]])), ])
dbs=data.frame(x=dib$x,
               pn1=sm_vec(dib$pn1),
               pn3=sm_vec(dib$pn3),
               pn5=sm_vec(dib$pn5),
               plda2=sm_vec(dib$plda2),
               pbeta=sm_vec(dib$pbeta)
)

graph<-ggplot(dbs) +
  geom_line(aes(x = x, y = pn1, color = "NCI1")) +
  geom_line(aes(x = x, y = pn3, color = "NCI3")) +
  geom_line(aes(x = x, y = pn5, color = "NCI5")) +
  geom_line(aes(x = x, y = plda2, color = "LDA2")) +
  geom_line(aes(x = x, y = pbeta, color = "Beta")) +
  
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
      "Beta"="purple",
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

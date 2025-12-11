source("libraries.R")
load("data/params.RData")
load("data/train/percentiles_burr.RData")

ylim<-c(0,35)
xlim<-c(-3,0)
echx<-seq(-3,0,by=.5)
echelle<-seq(0,100,by=2.5)
sm_vec <- function(x, m=10) {
  stats::filter(x, rep(1/m, m), sides = 2)
}

db=data.frame(x=rhos_train,
              pn1=psn1_burr,pn2=psn2_burr,pn3=psn3_burr,pn4=psn4_burr,pn5=psn5_burr,
              pg1=psg1_burr,pg2=psg2_burr
)



dib<-as.data.frame(db[rev(order(db[["x"]])), ])
dbs=data.frame(x=dib$x,
               pn1=sm_vec(dib$pn1),pn2=sm_vec(dib$pn2),
               pn3=sm_vec(dib$pn3),pn4=sm_vec(dib$pn4),
               pn5=sm_vec(dib$pn5),pg1=sm_vec(dib$pg1),pg2=sm_vec(dib$pg2)
)

graph<-ggplot(dbs) +
  geom_line(aes(x = x, y = pn1, color = "NCI1")) +
  geom_line(aes(x = x, y = pn2, color = "NCI2")) +
  geom_line(aes(x = x, y = pn3, color = "NCI3")) +
  geom_line(aes(x = x, y = pn4, color = "NCI4")) +
  geom_line(aes(x = x, y = pn5, color = "NCI5")) +
  geom_line(aes(x = x, y = pg1, color = "GCI1")) +
  geom_line(aes(x = x, y = pg2, color = "GCI2")) +

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
load("data/train/percentiles_invg.RData")

ylim<-c(0,35)
xlim<-c(-3,0)
echx<-seq(-3,0,by=.5)
echelle<-seq(0,100,by=2.5)
sm_vec <- function(x, m=10) {
  stats::filter(x, rep(1/m, m), sides = 2)
}

db=data.frame(x=rhos_train,
              pn1=psn1_invg,pn2=psn2_invg,pn3=psn3_invg,pn4=psn4_invg,pn5=psn5_invg,
              pg1=psg1_invg,pg2=psg2_invg
)



dib<-as.data.frame(db[rev(order(db[["x"]])), ])
dbs=data.frame(x=dib$x,
               pn1=sm_vec(dib$pn1),pn2=sm_vec(dib$pn2),
               pn3=sm_vec(dib$pn3),pn4=sm_vec(dib$pn4),
               pn5=sm_vec(dib$pn5),pg1=sm_vec(dib$pg1),pg2=sm_vec(dib$pg2)
)

graph<-ggplot(dbs) +
  geom_line(aes(x = x, y = pn1, color = "NCI1")) +
  geom_line(aes(x = x, y = pn2, color = "NCI2")) +
  geom_line(aes(x = x, y = pn3, color = "NCI3")) +
  geom_line(aes(x = x, y = pn4, color = "NCI4")) +
  geom_line(aes(x = x, y = pn5, color = "NCI5")) +
  geom_line(aes(x = x, y = pg1, color = "GCI1")) +
  geom_line(aes(x = x, y = pg2, color = "GCI2")) +
  
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
load("data/train/percentiles_fish.RData")

ylim<-c(0,35)
xlim<-c(-3,0)
echx<-seq(-3,0,by=.5)
echelle<-seq(0,100,by=2.5)
sm_vec <- function(x, m=20) {
  stats::filter(x, rep(1/m, m), sides = 2)
}

db=data.frame(x=rhos_train,
              pn1=psn1_fish,pn2=psn2_fish,pn3=psn3_fish,pn4=psn4_fish,pn5=psn5_fish,
              pg1=psg1_fish,pg2=psg2_fish
)



dib<-as.data.frame(db[rev(order(db[["x"]])), ])
dbs=data.frame(x=dib$x,
               pn1=sm_vec(dib$pn1),pn2=sm_vec(dib$pn2),
               pn3=sm_vec(dib$pn3),pn4=sm_vec(dib$pn4),
               pn5=sm_vec(dib$pn5),pg1=sm_vec(dib$pg1),pg2=sm_vec(dib$pg2)
)

graph<-ggplot(dbs) +
  geom_line(aes(x = x, y = pn1, color = "NCI1")) +
  geom_line(aes(x = x, y = pn2, color = "NCI2")) +
  geom_line(aes(x = x, y = pn3, color = "NCI3")) +
  geom_line(aes(x = x, y = pn4, color = "NCI4")) +
  geom_line(aes(x = x, y = pn5, color = "NCI5")) +
  geom_line(aes(x = x, y = pg1, color = "GCI1")) +
  geom_line(aes(x = x, y = pg2, color = "GCI2")) +
  
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
load("data/train/percentiles_stud.RData")

ylim<-c(0,35)
xlim<-c(-3,0)
echx<-seq(-3,0,by=.5)
echelle<-seq(0,100,by=2.5)
sm_vec <- function(x, m=10) {
  stats::filter(x, rep(1/m, m), sides = 2)
}

db=data.frame(x=rhos_train,
              pn1=psn1_stud,pn2=psn2_stud,pn3=psn3_stud,pn4=psn4_stud,pn5=psn5_stud,
              pg1=psg1_stud,pg2=psg2_stud
)



dib<-as.data.frame(db[rev(order(db[["x"]])), ])
dbs=data.frame(x=dib$x,
               pn1=sm_vec(dib$pn1),pn2=sm_vec(dib$pn2),
               pn3=sm_vec(dib$pn3),pn4=sm_vec(dib$pn4),
               pn5=sm_vec(dib$pn5),pg1=sm_vec(dib$pg1),pg2=sm_vec(dib$pg2)
)

graph<-ggplot(dbs) +
  geom_line(aes(x = x, y = pn1, color = "NCI1")) +
  geom_line(aes(x = x, y = pn2, color = "NCI2")) +
  geom_line(aes(x = x, y = pn3, color = "NCI3")) +
  geom_line(aes(x = x, y = pn4, color = "NCI4")) +
  geom_line(aes(x = x, y = pn5, color = "NCI5")) +
  geom_line(aes(x = x, y = pg1, color = "GCI1")) +
  geom_line(aes(x = x, y = pg2, color = "GCI2")) +
  
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


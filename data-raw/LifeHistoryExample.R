

devtools::load_all()
library(LBSPR)

LifeHistoryExample<-new("LifeHistory")
LifeHistoryExample@Linf<-100
LifeHistoryExample@L50<-66
LifeHistoryExample@L95<-67
LifeHistoryExample@MK<-1.5
LifeHistoryExample@LW_A<-0.01
LifeHistoryExample@LW_B<-3

LifeHistoryExample@title<-"Example fish"
LifeHistoryExample@shortDescription<-"Simulated life history of a fish based on B-H invariants"
LifeHistoryExample@speciesName<-"Example fish"
LifeHistoryExample@L_type<-"TL"
LifeHistoryExample@L_units<-"cm"
LifeHistoryExample@Walpha_units<-"g"
LifeHistoryExample@K<-0.2
LifeHistoryExample@M<-0.3
LifeHistoryExample@t0<-0
LifeHistoryExample@Tmax<- floor(-log(0.01)/0.3)
LifeHistoryExample@Steep<-0.99

lbsprSimExample<-lbsprSimWrapper(LifeHistoryExample)

usethis::use_data(LifeHistoryExample, overwrite = TRUE)
usethis::use_data(lbsprSimExample, overwrite = TRUE)


# library(hydromad)
# setwd("U:/")

# data
setwd("C:/Users/kim079/Documents/model_optimisation_framework_v3")
source("packages/sac_with_state_input/R/sacramento.R")


if(!is.loaded("sma_sac")){
  if(Sys.info()["sysname"]=="Linux"){
    dyn.load("packages/sac_with_state_input/src/sac_with_state_input.so")
  } else {
    dyn.load("packages/sac_with_state_input/src/sac_with_state_input.dll")
  }
}


id<-401013
preprocess_dir<-"output/gr4j.calib.param.state.all.sites.preprocess"
input_ts_file<-paste0(preprocess_dir,"/state_error_simulation_data_",id,".csv")
param_file<-paste0(preprocess_dir,"/gr4j_params_",id,".csv")
orig_params<-read.csv(param_file,as.is=T)
area_m2<-orig_params$area_m2

input_ts<-read.csv(input_ts_file,as.is=T)
climate<-data.frame(P=input_ts$P,E=input_ts$E)

sac.run<-function(param,input){
  dt=86400
  climate<-input$climate
  area<-input$area
  # area is in m^2
  # climate is in mm/timestep
  
  xx<-sacramento.sim(climate, uztwm=param[1], uzfwm=param[2], uzk=param[3],
                     pctim=param[4], adimp=param[5], zperc=param[6], rexp=param[7],
                     lztwm=param[8], lzfsm=param[9], lzfpm=param[10], lzsk=param[11],
                     lzpk=param[12], pfree=param[13])/1000*area/dt
  # result is in m^3/sec
}

input<-list(area=area_m2,climate=climate)
param.ranges<-sacramento.ranges() #hydromad::hydromad.getOption("sacramento")
param<-sapply(param.ranges,mean)
out<-sac.run(param = param, input = input)

plot(out,type="l")


# sac.run.hydromad<-function(param,input){
#   dt=86400
#   climate<-input$climate
#   area<-input$area
#   # area is in m^2
#   # climate is in mm/timestep
#   
#   xx<-hydromad::sacramento.sim(climate, uztwm=param[1], uzfwm=param[2], uzk=param[3],
#                      pctim=param[4], adimp=param[5], zperc=param[6], rexp=param[7],
#                      lztwm=param[8], lzfsm=param[9], lzfpm=param[10], lzsk=param[11],
#                      lzpk=param[12], pfree=param[13])/1000*area/dt
#   # result is in m^3/sec
# }
# 
# out.hm<-sac.run.hydromad(param = param, input = input)
# lines(out.hm,col=2,lty=2)
# 
# range(out-out.hm)


# testing new state values
dt=86400
climate<-input$climate
area<-input$area
# area is in m^2
# climate is in mm/timestep

xx<-sacramento.sim(climate, uztwm=param[1], uzfwm=param[2], uzk=param[3],
                   pctim=param[4], adimp=param[5], zperc=param[6], rexp=param[7],
                   lztwm=param[8], lzfsm=param[9], lzfpm=param[10], lzsk=param[11],
                   lzpk=param[12], pfree=param[13])/1000*area/dt

plot(xx,type="l",xlim=c(900,1200))
state_S_ts<-rep(-1e6,length(input$climate[,1])-1)
test_S<-seq(0,700,by=10)
for(i in 1:length(test_S)){
  state_S_ts[1000]<-test_S[i]
  xx<-sacramento.sim(climate, uztwm=param[1], uzfwm=param[2], uzk=param[3],
                     pctim=param[4], adimp=param[5], zperc=param[6], rexp=param[7],
                     lztwm=param[8], lzfsm=param[9], lzfpm=param[10], lzsk=param[11],
                     lzpk=param[12], pfree=param[13],state_S_ts=state_S_ts)/1000*area/dt
  lines(xx,type="l",col=2)
  cat(xx[1000:1010],"\n")
}



# testing return state
xx<-sacramento.sim(climate, uztwm=param[1], uzfwm=param[2], uzk=param[3],
                   pctim=param[4], adimp=param[5], zperc=param[6], rexp=param[7],
                   lztwm=param[8], lzfsm=param[9], lzfpm=param[10], lzsk=param[11],
                   lzpk=param[12], pfree=param[13])/1000*area/dt
plot(xx,type="l")
xx_rs<-sacramento.sim(climate, uztwm=param[1], uzfwm=param[2], uzk=param[3],
                      pctim=param[4], adimp=param[5], zperc=param[6], rexp=param[7],
                      lztwm=param[8], lzfsm=param[9], lzfpm=param[10], lzsk=param[11],
                      lzpk=param[12], pfree=param[13],return_state = T)
xx2<-xx_rs[,1]/1000*area/dt
lines(xx2,col=2,lty=2)

range(xx2-xx)

# testing return state with new state values
state_S_ts<-rep(-1e6,length(input$climate[,1])-1)
state_S_ts[1000]<-400
xx<-sacramento.sim(climate, uztwm=param[1], uzfwm=param[2], uzk=param[3],
                   pctim=param[4], adimp=param[5], zperc=param[6], rexp=param[7],
                   lztwm=param[8], lzfsm=param[9], lzfpm=param[10], lzsk=param[11],
                   lzpk=param[12], pfree=param[13],state_S_ts=state_S_ts)/1000*area/dt
plot(xx,type="l",xlim=c(900,1100))

which(is.infinite(xx))

xx_rs<-sacramento.sim(climate, uztwm=param[1], uzfwm=param[2], uzk=param[3],
                      pctim=param[4], adimp=param[5], zperc=param[6], rexp=param[7],
                      lztwm=param[8], lzfsm=param[9], lzfpm=param[10], lzsk=param[11],
                      lzpk=param[12], pfree=param[13],return_state = T,state_S_ts=state_S_ts)
xx2<-xx_rs[,1]/1000*area/dt
lines(xx2,col=2,lty=2)

range(xx2-xx)

# check adjustment to UZFWC
xx<-sacramento.sim(climate, uztwm=param[1], uzfwm=param[2], uzk=param[3],
                      pctim=param[4], adimp=param[5], zperc=param[6], rexp=param[7],
                      lztwm=param[8], lzfsm=param[9], lzfpm=param[10], lzsk=param[11],
                      lzpk=param[12], pfree=param[13],return_state = F)/1000*area/dt
plot(xx,type="l",xlim=c(900,1100))


xx<-sacramento.sim(climate, uztwm=param[1], uzfwm=param[2], uzk=param[3],
                   pctim=param[4], adimp=param[5], zperc=param[6], rexp=param[7],
                   lztwm=param[8], lzfsm=param[9], lzfpm=param[10], lzsk=param[11],
                   lzpk=param[12], pfree=param[13],return_state = F,state_S_ts=state_S_ts)/1000*area/dt
lines(xx,col=2,lty=2)

xx<-sacramento.sim(climate, uztwm=param[1], uzfwm=param[2], uzk=param[3],
                   pctim=param[4], adimp=param[5], zperc=param[6], rexp=param[7],
                   lztwm=param[8], lzfsm=param[9], lzfpm=param[10], lzsk=param[11],
                   lzpk=param[12], pfree=param[13],return_state = F,state_S2_ts=state_S_ts)/1000*area/dt
lines(xx,col=3,lty=2)

xx<-sacramento.sim(climate, uztwm=param[1], uzfwm=param[2], uzk=param[3],
                   pctim=param[4], adimp=param[5], zperc=param[6], rexp=param[7],
                   lztwm=param[8], lzfsm=param[9], lzfpm=param[10], lzsk=param[11],
                   lzpk=param[12], pfree=param[13],return_state = F,
                   state_S_ts=state_S_ts,state_S2_ts=state_S_ts)/1000*area/dt
lines(xx,col=4,lty=2)

# check return state
xx<-sacramento.sim(climate, uztwm=param[1], uzfwm=param[2], uzk=param[3],
                   pctim=param[4], adimp=param[5], zperc=param[6], rexp=param[7],
                   lztwm=param[8], lzfsm=param[9], lzfpm=param[10], lzsk=param[11],
                   lzpk=param[12], pfree=param[13],return_state = F,state_S2_ts=state_S_ts)/1000*area/dt
plot(xx,type="l",xlim=c(900,1100))
xx_rs<-sacramento.sim(climate, uztwm=param[1], uzfwm=param[2], uzk=param[3],
                   pctim=param[4], adimp=param[5], zperc=param[6], rexp=param[7],
                   lztwm=param[8], lzfsm=param[9], lzfpm=param[10], lzsk=param[11],
                   lzpk=param[12], pfree=param[13],return_state = T,state_S2_ts=state_S_ts)
xx2<-xx_rs[,1]/1000*area/dt
lines(xx2,col=2,lty=2)

range(xx2-xx)


# check lzfpc
xx<-sacramento.sim(climate, uztwm=param[1], uzfwm=param[2], uzk=param[3],
                   pctim=param[4], adimp=param[5], zperc=param[6], rexp=param[7],
                   lztwm=param[8], lzfsm=param[9], lzfpm=param[10], lzsk=param[11],
                   lzpk=param[12], pfree=param[13],return_state = F)/1000*area/dt

xx2<-sacramento.sim(climate, uztwm=param[1], uzfwm=param[2], uzk=param[3],
                   pctim=param[4], adimp=param[5], zperc=param[6], rexp=param[7],
                   lztwm=param[8], lzfsm=param[9], lzfpm=param[10], lzsk=param[11],
                   lzpk=param[12], pfree=param[13],return_state = F,
                   state_S_ts=NA,state_S2_ts=NA,state_S3_ts=state_S_ts)/1000*area/dt
plot(xx,type="l",xlim=c(900,1100))
lines(xx2,col=2,lty=2)

# check return state
xx<-sacramento.sim(climate, uztwm=param[1], uzfwm=param[2], uzk=param[3],
                   pctim=param[4], adimp=param[5], zperc=param[6], rexp=param[7],
                   lztwm=param[8], lzfsm=param[9], lzfpm=param[10], lzsk=param[11],
                   lzpk=param[12], pfree=param[13],return_state = F,state_S3_ts=state_S_ts)/1000*area/dt
plot(xx,type="l",xlim=c(900,1100))
xx_rs<-sacramento.sim(climate, uztwm=param[1], uzfwm=param[2], uzk=param[3],
                      pctim=param[4], adimp=param[5], zperc=param[6], rexp=param[7],
                      lztwm=param[8], lzfsm=param[9], lzfpm=param[10], lzsk=param[11],
                      lzpk=param[12], pfree=param[13],return_state = T,state_S3_ts=state_S_ts)
xx2<-xx_rs[,1]/1000*area/dt
lines(xx2,col=2,lty=2)

range(xx2-xx)

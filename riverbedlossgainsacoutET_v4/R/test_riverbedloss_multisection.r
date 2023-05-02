setwd("C:/Users/kim079/OneDrive - CSIRO/Documents/bradfield/transmission_losses")

source("packages/riverbedlossgainsacoutET_v2/R/riverbedloss.run.r")
source("packages/riverbedlossgainsacoutET_v2/R/sacramento.R")
source("packages/riverbedlossgainsacoutET_v2/R/get.constant.r")


if(!is.loaded("riverbedlossreach_run")){
  if(Sys.info()["sysname"]=="Linux"){
    dyn.load("packages/riverbedlossgainsacoutET_v2/src/riverbedlossgainsacoutET_v2_petrichor.so")
  } else {
    dyn.load("packages/riverbedlossgainsacoutET_v2/src/riverbedlossgainsacoutET_v2.dll")
  }
}

# there are 24 states types (columns) for each subcatment
nstates_subcat<-24


# set.seed(1)

length_ts<-10
sacpar<-sapply(sacramento.ranges(),mean)


parameters<-c()
parameters[1]<-invv<-1 # sec/m
parameters[2]<-alpha<-0.5
parameters[3]<-wetAreaFactor<-1.2
parameters[4]<-fc<-1 #m^3/s/m
parameters[5]<-max_rivbedstore_per_m<-1e5 #m^3/m
parameters[6]<-lossK<-0.6
##################


nsubcat<-c(1,1)

config<-c()
config[1]<-timestep_length<-86400
config[2]<-river_length<-2e5
config[3]<-river_area<-1e5
config[4]<-subcat_area<-2e6

config_l<-c(config,config)
inputs<-cbind(runif(length_ts,0,3),runif(length_ts,0,3),runif(length_ts,0,0))

inputs_2<-cbind(inputs,inputs)

out1<-riverbedloss.run(parameters,config,inputs,nsubcat=1,states=NULL,outputOption=0,sacpar,rep(0,6),initial_riverbedstore=0.01)
out2<-riverbedloss.run(parameters,config_l,inputs_2,nsubcat,states=NULL,outputOption=0,sacpar,rep(0,6),initial_riverbedstore=0.01)

layout(1:2)
plot(out1$outflow,type="l")
plot(out2$outflow,type="l")

out1$outflow-out2$outflow

layout(1:4)
plot(out2$states[,1],type="l")
plot(out2$states[,1+nstates_subcat],type="l")
plot(out2$outflow,type="l")


head(out1$states)
head(out2$states)

head(out2$states[,1])
head(out2$states[,1+nstates_subcat])


##################
nsubcat<-c(1,1,1)

config<-c()
config[1]<-timestep_length<-86400
config[2]<-river_length<-2e5
config[3]<-river_area<-1e5
config[4]<-subcat_area<-2e6

config_l<-c(config,config,config)
inputs<-cbind(runif(length_ts,0,3),runif(length_ts,0,3),runif(length_ts,0,0))

inputs_2<-cbind(inputs,inputs,inputs)

out1<-riverbedloss.run(parameters,config,inputs,nsubcat=1,states=NULL,outputOption=0,sacpar,rep(0,6),initial_riverbedstore=0.01)
out2<-riverbedloss.run(parameters,config_l,inputs_2,nsubcat,states=NULL,outputOption=0,sacpar,rep(0,6),initial_riverbedstore=0.01)

layout(1:2)
plot(out1$outflow,type="l")
plot(out2$outflow,type="l")

out1$outflow-out2$outflow

layout(1:4)
plot(out2$states[,1],type="l")
plot(out2$states[,1+nstates_subcat],type="l")
plot(out2$states[,1+(nstates_subcat*2)],type="l")
plot(out2$outflow,type="l")

out2$states[,1+nstates_subcat]-out2$states[,1+(nstates_subcat*2)]

head(out1$states)
head(out2$states)

head(out2$states[,1])
head(out2$states[,1+nstates_subcat])

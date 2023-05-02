setwd("C:/Users/kim079/OneDrive - CSIRO/Documents/bradfield/transmission_losses")

source("packages/riverbedlossgainsacoutET_v4/R/riverbedloss.run.r")
source("packages/riverbedlossgainsacoutET_v4/R/sacramento.R")
source("packages/riverbedlossgainsacoutET_v4/R/get.constant.r")


if(!is.loaded("riverbedlossreach_run")){
  if(Sys.info()["sysname"]=="Linux"){
    dyn.load("packages/riverbedlossgainsacoutET_v4/src/riverbedlossgainsacoutET_v4.so")
  } else {
    dyn.load("packages/riverbedlossgainsacoutET_v4/src/riverbedlossgainsacoutET_v4.dll")
  }
}

# there are 25 states types (columns) for each subcatment
# 1: discharge
# 2: Smax
# 3: river loss (m^3/timestep)
# 4: river bed/bank store (m^3)
# 5: river bed/bank store drying (m^3/timestep)
# 6: river rain (m^3/timestep)
# 7: river evap (m^3/timestep)
# 8: river drying (m^3/timestep)
# 9: river volume
# 10: runoff (m^3/s)
# 11: river bed/bank store leakage (m^3/timestep)
# 12: rainfall runoff total ET (m^3/s)
# 13: RBS mass balance error (m^3/m)
# remaining 12 states are for routing equations 


nstates_subcat<-25


parameters<-c()
parameters[1]<-invv<-1 # sec/m
parameters[2]<-alpha<-0.5
parameters[3]<-wetAreaFactor<-1.2
parameters[4]<-fc<-1e-6 #m^3/s/m
parameters[5]<-max_rivbedstore_per_m<-1e5 #m^3/m
parameters[6]<-lossK<-0.6


config<-c()
config[1]<-timestep_length<-86400
config[2]<-river_length<-2e5
config[3]<-river_area<-1e5
config[4]<-subcat_area<-2e6

set.seed(1)

length_ts<-1000 #1000

# INPUTS:
# column 1: rain (mm)
# column 2: PET (mm)
# column 3: streamflow (m^3/s)
# inputs<-cbind(runif(length_ts,0,30),runif(length_ts,0,3),runif(length_ts,0,30))
# inputs<-rbind(inputs,cbind(runif(length_ts,0,0),runif(length_ts,0,3),runif(length_ts,0,0)))
inputs<-cbind(runif(length_ts,0,30),runif(length_ts,20,30),runif(length_ts,0,30))


nsubcat<-1

sacpar<-sapply(sacramento.ranges(),mean)
sacinitpar<-rep(0,6) #rep(0.3,6)
# 
# 
out1<-riverbedloss.run(parameters,config,inputs,nsubcat=1,states=NULL,outputOption=0,sacpar,sacinitpar,initial_riverbedstore=0.01)

inputs2<-inputs
inputs2[,1]<-0
inputs2[,2]<-0

out2<-riverbedloss.run(parameters,config,inputs2,nsubcat,states=NULL,outputOption=0,sacpar,rep(0,6),initial_riverbedstore=0.01)

out3<-riverbedloss.run(parameters,config,inputs,nsubcat=1,states=NULL,outputOption=0,sacpar,sacinitpar,initial_riverbedstore=0.01,excl_channel = 1)



# 
# # cat(out1$outflow)
# plot(out1$outflow,type="l")
# 
# head(out1$states.nonrouting)
# head(out1$states.routing)
# 
# plot(out1$states.nonrouting[,10],type="l")
# 
# 
# plot(out1$states.nonrouting[,5],type="l")
# 
plot(out1$outflow,type="l")
lines(out2$outflow,col=2)

range(out1$outflow-out2$outflow)
range(out1$states[,10])
range(out2$states[,10])

plot(out1$states[,10],type="l") # rainfall runoff
plot(out2$states[,10],type="l") # rainfall runoff

# plot(out1$states[,5],type="l") # river bank store AET
# plot(out2$states[,5],type="l") # river bank store AET

plot(out1$outflow,type="l")
lines(out3$outflow,col=2)

mean(out1$outflow)
mean(out3$outflow)

mean(out1$states[,10])
mean(out3$states[,10])


plot(out1$states[,4],type="l")
out1$states[,13]

# test state adjustment
riverbedstore_ts<-rep(-1e6,nrow(inputs)-1)
riverbedstore_ts[400]<-100
out_test_mass_bal<-riverbedloss.run(parameters,config,inputs,nsubcat=1,states=NULL,outputOption=0,sacpar,sacinitpar,initial_riverbedstore=0.01,riverbedstore_ts=riverbedstore_ts)
plot(out_test_mass_bal$states[,4],type="l")
out_test_mass_bal$states[,13]

out_test_mass_bal<-riverbedloss.run(parameters,config,inputs,nsubcat=1,states=NULL,outputOption=0,sacpar,sacinitpar,initial_riverbedstore=0.01,riverbedstore_ts=riverbedstore_ts,max_RBS_mass_balance = 100)
plot(out_test_mass_bal$states[,4],type="l")
out_test_mass_bal$states[,13]

riverbedstore_ts<-rep(-1e6,nrow(inputs)-1)
riverbedstore_ts[400]<-100
riverbedstore_ts[600]<-100
out_test_mass_bal<-riverbedloss.run(parameters,config,inputs,nsubcat=1,states=NULL,outputOption=0,sacpar,sacinitpar,
                                    initial_riverbedstore=0.01,riverbedstore_ts=riverbedstore_ts,max_RBS_mass_balance = NA)
plot(out_test_mass_bal$states[,4],type="l")
out_test_mass_bal$states[401,4]
out_test_mass_bal$states[601,4]
out_test_mass_bal$states[,13]


riverbedstore_ts<-rep(-1e6,nrow(inputs)-1)
riverbedstore_ts[400]<-100
riverbedstore_ts[600]<-100
out_test_mass_bal<-riverbedloss.run(parameters,config,inputs,nsubcat=1,states=NULL,outputOption=0,sacpar,sacinitpar,
                                    initial_riverbedstore=0.01,riverbedstore_ts=riverbedstore_ts,max_RBS_mass_balance = 60)
plot(out_test_mass_bal$states[,4],type="l")
out_test_mass_bal$states[401,4]
out_test_mass_bal$states[601,4]
out_test_mass_bal$states[,13]

riverbedstore_ts<-rep(-1e6,nrow(inputs)-1)
riverbedstore_ts[400]<-100
# riverbedstore_ts[450]<-50
riverbedstore_ts[600]<-100
out_test_mass_bal<-riverbedloss.run(parameters,config,inputs,nsubcat=1,states=NULL,outputOption=0,sacpar,sacinitpar,
                                    initial_riverbedstore=0.01,riverbedstore_ts=riverbedstore_ts,max_RBS_mass_balance = 140)
plot(out_test_mass_bal$states[,4],type="l")
out_test_mass_bal$states[401,4]
out_test_mass_bal$states[601,4]
out_test_mass_bal$states[,13]


riverbedstore_ts<-rep(-1e6,nrow(inputs)-1)
riverbedstore_ts[400]<-0
# riverbedstore_ts[450]<-50
riverbedstore_ts[550]<-0
out_test_mass_bal<-riverbedloss.run(parameters,config,inputs,nsubcat=1,states=NULL,outputOption=0,sacpar,sacinitpar,
                                    initial_riverbedstore=0.01,max_RBS_mass_balance = NA)
plot(out_test_mass_bal$states[,4],type="l",log="y")
out_test_mass_bal<-riverbedloss.run(parameters,config,inputs,nsubcat=1,states=NULL,outputOption=0,sacpar,sacinitpar,
                                    initial_riverbedstore=0.01,riverbedstore_ts=riverbedstore_ts,max_RBS_mass_balance = NA)
lines(out_test_mass_bal$states[,4],col=2,lty=2)
out_test_mass_bal$states[,13]

riverbedstore_ts<-rep(-1e6,nrow(inputs)-1)
riverbedstore_ts[400]<-0
# riverbedstore_ts[450]<-50
riverbedstore_ts[550]<-0
out_test_mass_bal<-riverbedloss.run(parameters,config,inputs,nsubcat=1,states=NULL,outputOption=0,sacpar,sacinitpar,
                                    initial_riverbedstore=0.01,max_RBS_mass_balance = NA)
plot(out_test_mass_bal$states[,4],type="l",log="y")
out_test_mass_bal<-riverbedloss.run(parameters,config,inputs,nsubcat=1,states=NULL,outputOption=0,sacpar,sacinitpar,
                                    initial_riverbedstore=0.01,riverbedstore_ts=riverbedstore_ts,max_RBS_mass_balance = 5)
lines(out_test_mass_bal$states[,4],col=2,lty=2)
out_test_mass_bal$states[,13]

riverbedstore_ts<-rep(-1e6,nrow(inputs)-1)
riverbedstore_ts[400]<-0
# riverbedstore_ts[450]<-50
riverbedstore_ts[550]<-0
out_test_mass_bal<-riverbedloss.run(parameters,config,inputs,nsubcat=1,states=NULL,outputOption=0,sacpar,sacinitpar,
                                    initial_riverbedstore=0.01,max_RBS_mass_balance = NA)
plot(out_test_mass_bal$states[,4],type="l",log="y")
out_test_mass_bal<-riverbedloss.run(parameters,config,inputs,nsubcat=1,states=NULL,outputOption=0,sacpar,sacinitpar,
                                    initial_riverbedstore=0.01,riverbedstore_ts=riverbedstore_ts,max_RBS_mass_balance = 20)
lines(out_test_mass_bal$states[,4],col=2,lty=2)
out_test_mass_bal$states[,13]


riverbedstore_ts<-rep(-1e6,nrow(inputs)-1)
riverbedstore_ts[400]<-0
# riverbedstore_ts[450]<-50
riverbedstore_ts[550]<-0
riverbedstore_ts[600]<-50
riverbedstore_ts[650]<-50
riverbedstore_ts[700]<-10

out_test_mass_bal<-riverbedloss.run(parameters,config,inputs,nsubcat=1,states=NULL,outputOption=0,sacpar,sacinitpar,
                                    initial_riverbedstore=0.01,riverbedstore_ts=riverbedstore_ts,max_RBS_mass_balance = 20)
plot(out_test_mass_bal$states[,4],type="l",log="y")
out_test_mass_bal$states[,13]

stop()

# 
# plot(out3$states[,sac_index],type="l") # rainfall runoff
# 

# test without rainfall runoff

out1<-riverbedloss.run(parameters,config,inputs,nsubcat=1,states=NULL,outputOption=0,sacpar,sacinitpar,initial_riverbedstore=0.01)
out2<-riverbedloss.run(parameters,config,inputs,nsubcat=1,states=NULL,outputOption=0,sacpar,sacinitpar,initial_riverbedstore=0.01,use_RR = 0)

plot(out1$outflow,type="l")
lines(out2$outflow,col=2)

range(out1$outflow-out2$outflow)
range(out1$states[,10])
range(out2$states[,10])

plot(out1$states[,10],type="l") # rainfall runoff
plot(out2$states[,10],type="l") # rainfall runoff



# test 2 subcatchments

nsubcat<-2

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

sac_index<-10
plot(out2$states[,sac_index],type="l") # rainfall runoff
plot(out2$states[,out2$nstates+sac_index],type="l") # rainfall runoff

out2$states[,sac_index]
out2$states[,out2$nstates+sac_index]

plot(out2$outflow,type="l")
lines(out2$states[,sac_index],col=2,lty=2)
lines(out2$states[,out2$nstates+sac_index],col=3,lty=3)

# test multiple sections
nsubcat<-c(1,1)

config<-c()
config[1]<-timestep_length<-86400
config[2]<-river_length<-2e5
config[3]<-river_area<-1e5
config[4]<-subcat_area<-2e6

config_l<-c(config,config)
inputs<-cbind(runif(length_ts,0,0),runif(length_ts,0,0),runif(length_ts,0,30))

inputs_2<-cbind(inputs,inputs)

out1<-riverbedloss.run(parameters,config,inputs,nsubcat=1,states=NULL,outputOption=0,sacpar,rep(0,6),initial_riverbedstore=0.01)
out2<-riverbedloss.run(parameters,config_l,inputs_2,nsubcat,states=NULL,outputOption=0,sacpar,rep(0,6),initial_riverbedstore=0.01)

layout(1:2)
plot(out1$outflow,type="l")
plot(out2$outflow,type="l")

out1$outflow-out2$outflow

layout(1:3)
plot(out2$states[,1],type="l") # outflow for first subcatchment
plot(out2$states[,1+nstates_subcat],type="l") # outflow for second subcatchment
plot(out2$outflow,type="l")


# test 3 subcatchments
nsubcat<-3

config<-c()
config[1]<-timestep_length<-86400
config[2]<-river_length<-2e5
config[3]<-river_area<-1e5
config[4]<-subcat_area<-2e6

config_l<-c(config,config,config)
inputs<-cbind(runif(length_ts,0,3),runif(length_ts,0,3),runif(length_ts,0,30))

inputs_3<-cbind(inputs,inputs,inputs)

out3<-riverbedloss.run(parameters,config_l,inputs_3,nsubcat,states=NULL,outputOption=0,sacpar,sacinitpar,initial_riverbedstore=0.01)

layout(1:3)
plot(out1$outflow,type="l")
plot(out2$outflow,type="l")
plot(out3$outflow,type="l")

mean(out1$outflow)
mean(out2$outflow)
mean(out3$outflow)

head(out3$states.nonrouting)
head(out3$inputs)
out3$states

layout(1:3)
sac_index<-10
plot(out3$states[,sac_index],type="l") # rainfall runoff
plot(out3$states[,out3$nstates+sac_index],type="l") # rainfall runoff
plot(out3$states[,out3$nstates*2+sac_index],type="l") # rainfall runoff

# these should be the same


# test multiple sections

nsubcat<-c(2,3)

config<-c()
config[1]<-timestep_length<-86400
config[2]<-river_length<-2e5
config[3]<-river_area<-1e5
config[4]<-subcat_area<-2e6

config_l<-c(config,config,config,config,config)

inputs<-cbind(runif(length_ts,0,3),runif(length_ts,0,3),runif(length_ts,0,30))
inputs_5<-cbind(inputs,inputs,inputs,inputs,inputs)

out4<-riverbedloss.run(parameters,config_l,inputs_5,nsubcat,states=NULL,outputOption=0,sacpar,sacinitpar,initial_riverbedstore=0.01)

layout(1:2)
plot(out4$outflow,type="l")

out4$states

out4$inputs

# multiple sections again 

nsubcat<-c(2,1,1)

config<-c()
config[1]<-timestep_length<-86400
config[2]<-river_length<-2e5
config[3]<-river_area<-1e5
config[4]<-subcat_area<-2e6

config_l<-c(config,config,config,config)

inputs<-cbind(runif(length_ts,0,3),runif(length_ts,0,3),runif(length_ts,0,30))
inputs_6<-cbind(inputs,inputs,inputs,inputs)

out4<-riverbedloss.run(parameters,config_l,inputs_6,nsubcat,states=NULL,outputOption=0,sacpar,sacinitpar,initial_riverbedstore=0.01)

layout(1:2)
plot(out4$outflow,type="l")

out4$states

out4$inputs


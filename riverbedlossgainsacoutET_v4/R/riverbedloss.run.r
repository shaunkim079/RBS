#' Run river bed/bank store model
#'
#' Caution with the units !!! The following units are used
#'	\itemize{
#'		\item Time is expressed in seconds
#'		\item Flows are expressed in m3/s
#'		\item Volumes are expressed in m3 
#'		\item Climate variables are expressed in mm/timestep
#'	}
#' @export
#' @useDynLib riverbedloss
#' @examples
#'
#'
#'
#'

riverbedloss.run <- function(parameters,config,inputs,nsubcat,states=NULL,outputOption=0,
                             sacpar,sacinitpar,riverbedstore_ts=NA,initial_riverbedstore=0,
                             riverbedstore_fact_ts=NA,initial_riverbedstore_fact=NA,
                             use_RR=1,excl_channel=0,
                             max_RBS_mass_balance=NA){
  if(is.null(dim(inputs))) stop("inputs should have at least 2 dimensions")
	# inputs, remember this is for each subcatchment
	nval     <-as.integer(nrow(inputs))

	# ninputs  <-as.integer(ncol(inputs))
	ninputs<-as.integer(3)
  ninflow  <- ninputs - get.constant("NINPUTSNONROUTING") # number of inflows for each subcatchment
	if(ninflow<0) stop("ninflow <0")
  
	nconfig  <-as.integer(get.constant("NCONFIG")) 

	npar     <-as.integer(ninflow*get.constant("NPARROUTING")+
							get.constant("NPARNONROUTING")) 

	nstates.nonrouting  <- get.constant("NSTATESNONROUTING")
	nstates.routing  	<- get.constant("NSTATESROUTING")
	nstates  <-as.integer(ninflow*nstates.routing+nstates.nonrouting) 

	if(nval==0) stop("nval = 0")
	
	total_nsubcat<-sum(nsubcat)
	if(length(config)/total_nsubcat!=nconfig) stop(sprintf("length(config)/total_nsubcat!=nconfig (%d)\n",
                                          nconfig))

	if(length(parameters)!=npar) stop(sprintf("length(parameters)!=npar (%d)\n",
                                          npar))

	
	if(total_nsubcat != ncol(inputs)/ninputs) stop("problem with the number of inputs")
	
  ierr        <- as.integer(0)  
  nm.config 	<- names(config)
  config_array      <- as.double(unlist(config)) 
  parameters  <- as.double(unlist(parameters)) 
  # inputs      <- as.double(unlist(t(inputs)[1:(nval*ninputs)])) 

  
  # tmp<-matrix(1:36,ncol=6)
  # inputs_array<-c()
  # input_start_cols<-seq(1,ncol(tmp),by=ninputs)
  # for(i in input_start_cols){
  #   inputs_array<-c(inputs_array,c(t(tmp[,i:(i+ninputs-1)])))
  # }
  
  inputs_array<-c()
  # inputs_array<-rep(NA,ncol(inputs)*nval)
  input_start_cols<-seq(1,ncol(inputs),by=ninputs)
  for(i in input_start_cols){
    inputs_array<-c(inputs_array,c(t(inputs[,i:(i+ninputs-1)])))
    # inputs_array[]<-c(t(inputs[,i:(i+ninputs-1)]))
  }
  
  # outputs
  if(is.null(states)){
    states_array<-as.double(rep(0,nval*nstates*total_nsubcat))
  } else {
    states_array<-c()
    states_start_cols<-seq(1,ncol(states),by=nstates)
    for(i in states_start_cols){
      states_array<-c(states_array,c(t(states[,i:(i+nstates-1)])))
    }
    # states<-as.double(states)
    if(length(states_array)!=nval*nstates*total_nsubcat){
      stop("The input states is not the right length (nval*nstates*total_nsubcat = ",nval*nstates*total_nsubcat,")")
    }
    states_array<-as.double(states_array)
  }

  
  if(is.na(riverbedstore_ts[1])){
    riverbedstore_ts<-as.double(rep(-1e6,nval-1))
    use_riverbedstore_ts<-as.integer(0)
  } else {
    stopifnot(length(riverbedstore_ts) == (nval-1))
    if(any(riverbedstore_ts[!(riverbedstore_ts==-1e6)]<0)) stop("shouldn't be riverbedstore < 0")
    riverbedstore_ts<-as.double(riverbedstore_ts)
    use_riverbedstore_ts<-as.integer(1)
  }
  
  if(is.na(riverbedstore_fact_ts[1])){
    riverbedstore_fact_ts<-as.double(rep(-1e6,nval-1))
    use_riverbedstore_fact_ts<-as.integer(0)
  } else {
    stopifnot(length(riverbedstore_fact_ts) == (nval-1))
    if(any(riverbedstore_fact_ts[!(riverbedstore_fact_ts==-1e6)]<0)) stop("shouldn't be riverbedstore < 0")
    riverbedstore_fact_ts<-as.double(riverbedstore_fact_ts)
    use_riverbedstore_fact_ts<-as.integer(1)
    
    # override the abs adjustment
    riverbedstore_ts<-as.double(rep(-1e6,nval-1))
    use_riverbedstore_ts<-as.integer(0)
  }
  
  if(!is.na(initial_riverbedstore_fact)){
    initial_riverbedstore_fact<- as.double(initial_riverbedstore_fact)
  } else {
    initial_riverbedstore_fact<- as.double(-1e6)
  }
  
  if(!is.na(max_RBS_mass_balance)){
    max_RBS_mass_balance<-as.double(max_RBS_mass_balance)
  } else {
    max_RBS_mass_balance<-as.double(-1e6)
  }

	# Run the C code
	#void awrarirrig_run(
	#		int * ierr,
	#		int *nval, int *nconfig, int * ninputs,int * npar,
	#	    int * nstates,
	#		double * config,
	#		double * inputs_array,
	#		double * parameters,
	#		double * states_array,
	#		double * use_irrig,
	#		double * irrigConfig,
	#		double * irrigInputs_array,
	#		double * irrigParameters,
	#		double * irrigStates_array)	
  
	out <- .C("riverbedlossreach_run",
    ierr=ierr,
    nval=nval,nconfig=nconfig,ninputs=ninputs,
	  npar=npar,
    nsections=as.integer(length(nsubcat)),nsubcat=as.integer(nsubcat),
    nstates=nstates,
    config=config_array,
    inputs_array=as.double(inputs_array),
    parameters=parameters,
    states_array=states_array,
    sacpar=as.double(sacpar),
    sacinitpar=as.double(sacinitpar),
    riverbedstore_ts=riverbedstore_ts,
    use_riverbedstore_ts=use_riverbedstore_ts,
    initial_riverbedstore=as.double(initial_riverbedstore),
    outflow_total=as.double(rep(0,nval)),
    riverbedstore_fact_ts=riverbedstore_fact_ts,
    use_riverbedstore_fact_ts=use_riverbedstore_fact_ts,
    initial_riverbedstore_fact=initial_riverbedstore_fact,
    uztwc = double(nval),
    uzfwc = double(nval),
    lztwc = double(nval),
    lzfsc = double(nval),
    lzfpc = double(nval),
    adimc = double(nval),
    sett = double(nval),
    se1 = double(nval),
    se3 = double(nval),
    se4 = double(nval),
    se5 = double(nval),
    roimp = double(nval),
    sdro = double(nval),
    ssur = double(nval),
    sif = double(nval),
    bfp = double(nval),
    bfs = double(nval),
    bfcc = double(nval),
    use_RR=as.integer(use_RR),
    excl_channel=as.integer(excl_channel),
    max_RBS_mass_balance=max_RBS_mass_balance
	)

	if(out$ierr>0){
	  cat("parameters:",parameters,"\n")
	  cat("sacpar:",sacpar,"\n")
	  cat("sacinitpar:",sacinitpar,"\n")
	  # cat("riverbedstore_ts:",riverbedstore_ts,"\n")
	  cat("initial_riverbedstore:",initial_riverbedstore,"\n")
	  # cat("riverbedstore_fact_ts:",riverbedstore_fact_ts,"\n")
	  cat("initial_riverbedstore_fact:",initial_riverbedstore_fact,"\n")
	  cat("riverbedstore indices:",which(riverbedstore_ts!=-1e6),"\n")
	  cat("riverbedstore values:",riverbedstore_ts[which(riverbedstore_ts!=-1e6)],"\n")
	  cat("riverbedstore_fact indices:",which(riverbedstore_fact_ts!=-1e6),"\n")
	  cat("riverbedstore_fact values:",riverbedstore_fact_ts[which(riverbedstore_fact_ts!=-1e6)],"\n")
	  stop(sprintf("\nriverbedloss.run error - C code returned error %d\n",out$ierr))
	}
	  

  if(outputOption==0){
    # browser()
  	# reformat
  	# out$states <- t(matrix(out$states_array,nstates*total_nsubcat,nval))
    # states are written for each subcatchment first, then each timestep, then each state type
    out$states<-matrix(NA,nrow=nval,ncol=nstates*total_nsubcat)
  	for(i in 1:total_nsubcat){
  	  sc_states_array<-out$states_array[((i-1)*nstates*nval+1):(i*nstates*nval)]
  	  out$states[,((i-1)*nstates+1):(i*nstates)]<-t(matrix(sc_states_array,nstates,nval))
  	}
  	states_array<-out$states_array
  	out$states_array <- NULL
  	# out$inputs <- t(matrix(out$inputs_array,ninputs*total_nsubcat,nval))
  	
  	out$inputs<-matrix(NA,nrow=nval,ncol=ninputs*total_nsubcat)
  	for(i in 1:total_nsubcat){
  	  sc_inputs_array<-out$inputs_array[((i-1)*ninputs*nval+1):(i*ninputs*nval)]
  	  out$inputs[,((i-1)*ninputs+1):(i*ninputs)]<-t(matrix(sc_inputs_array,ninputs,nval))
  	}
  	out$inputs_array <- NULL



    out <- list(outflow=out$outflow_total, #out$states[,1],
                inputs = out$inputs,
                parameters = parameters,
                config = config,
                states.nonrouting=out$states[,1:nstates.nonrouting],
                states.routing=out$states[,(nstates.nonrouting+1):nstates],
                states_array=states_array,
                states=out$states,
                nstates=nstates)




    # names(out$parameters)[1:6] <- c("flood.returnflow.coefficient",
    #                                 "flood.ksat","river.monod1","river.monod2","overbank.flow.threshold",
    #                                 "flood.gamma")
    # names(out$config) <- nm.config

    # colnames(out$inputs) <- paste("inp",1:ncol(out$inputs),sep=".")
    # colnames(out$inputs)[1:11] <- c("river.rainfall",
    #                                 "river.evap","floodplain.rainfall","floodplain.evap",
    #                                 "irrigation.diversion","irrigation.returnflow",
    #                                 "urban.diversion","ungauged.inflow.top","ungauged.inflow.bottom",
    #                                 "reservoir.volume","reservoir.area")

    # names.states.nonrouting<-c("outflow","overbank.flow",
    #                            "floodplain.volume","floodplain.area","floodplain.returnflow",
    #                            "river.rainfall.flux","river.evap.flux",
    #                            "floodplain.rainfall.flux","floodplain.evap.flux",
    #                            "floodplain.groundwater.loss","river.groundwater.loss",
    #                            "top.anabranch.loss","bottom.anabranch.loss",
    #                            "previous.reservoir.vol","previous.reservoir.area",
    #                            "reservoir.rainfall.flux","reservoir.evap.flux",
    #                            "reservoir.contribution","reservoir.carryover.flux",
    #                            sprintf("outflow%d",2:7),"river.volume",
    #                            "floodplain.groundwater.max.change.storage",
    #                            "floodplain.groundwater.outflow",
    #                            "floodplain.groundwater.max.infiltration",
    #                            "river.groundwater.max.change.storage",
    #                            "river.groundwater.outflow",
    #                            "river.groundwater.max.infiltration",
    #                            "river.groundwater.max.monod.loss")
    # if(is.null(dim(out$states.nonrouting))){
    #   names(out$states.nonrouting)[1:33] <- names.states.nonrouting
    # } else {
    #   colnames(out$states.nonrouting)[1:33] <- names.states.nonrouting
    # }

    } else if(outputOption==1){
     # out$states <- t(matrix(out$states_array,nstates*total_nsubcat,nval))
     # states_array<-out$states_array
     # out$states_array <- NULL
     # out$inputs <- t(matrix(out$inputs_array,ninputs*total_nsubcat,nval))
     # out$inputs_array <- NULL
     
     # states are written for each subcatchment first, then each timestep, then each state type
     out$states<-matrix(NA,nrow=nval,ncol=nstates*total_nsubcat)
     for(i in 1:total_nsubcat){
       sc_states_array<-out$states_array[((i-1)*nstates*nval+1):(i*nstates*nval)]
       out$states[,((i-1)*nstates+1):(i*nstates)]<-t(matrix(sc_states_array,nstates,nval))
     }
     states_array<-out$states_array
     out$states_array <- NULL
     
     out$inputs<-matrix(NA,nrow=nval,ncol=ninputs*total_nsubcat)
     for(i in 1:total_nsubcat){
       sc_inputs_array<-out$inputs_array[((i-1)*ninputs*nval+1):(i*ninputs*nval)]
       out$inputs[,((i-1)*ninputs+1):(i*ninputs)]<-t(matrix(sc_inputs_array,ninputs,nval))
     }
     out$inputs_array <- NULL

     out <- list(states.nonrouting=out$states[,1:nstates.nonrouting])


     # names.states.nonrouting<-c("outflow","overbank.flow",
     #                            "floodplain.volume","floodplain.area","floodplain.returnflow",
     #                            "river.rainfall.flux","river.evap.flux",
     #                            "floodplain.rainfall.flux","floodplain.evap.flux",
     #                            "floodplain.groundwater.loss","river.groundwater.loss",
     #                            "top.anabranch.loss","bottom.anabranch.loss",
     #                            "previous.reservoir.vol","previous.reservoir.area",
     #                            "reservoir.rainfall.flux","reservoir.evap.flux",
     #                            "reservoir.contribution","reservoir.carryover.flux",
     #                            sprintf("outflow%d",2:7),"river.volume",
     #                            "floodplain.groundwater.max.change.storage",
     #                            "floodplain.groundwater.outflow",
     #                            "floodplain.groundwater.max.infiltration",
     #                            "river.groundwater.max.change.storage",
     #                            "river.groundwater.outflow",
     #                            "river.groundwater.max.infiltration",
     #                            "river.groundwater.max.monod.loss")
     # if(is.null(dim(out$states.nonrouting))){
     #    names(out$states.nonrouting)[1:33] <- names.states.nonrouting
     # } else {
     #    colnames(out$states.nonrouting)[1:33] <- names.states.nonrouting
     # }
     # end of output option 1

  } else if(outputOption==2){

     out <- list(states=out$states, nstates=out$nstates, nval=out$nval)

  }

  class(out) <- "riverbedloss-run"

	return(out)
}


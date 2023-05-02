#include "header.h"
#include "sacramento_state.h"
#include <R.h>

/******************************************************************************
Code to loop through timeseries and implement the timestep code

-- inputs ---
int * ierr,			Error code

int *nval, 			Dimensions of array (nval = number of timesteps)
int *nconfig, 
int * ninputs,
int * npar,
int * nstates, 


-- states ---
double * states_array 	model states

-- Authors --
Shaun Kim, CSIRO



*******************************************************************************/

void riverbedlosssc_runtimestep(
  int * ierr,
  int * nconfig,
  int * ninputs,
  int * npar,
  int * nstates,
  double * config,
  double * inputs,
  double * parameters,
  double * states,
  double * sacpar,
  double * sacinitpar,
  double * riverbedstore,
  double * riverbedstore_fact,
  double * max_RBS_mass_balance);
    				
void riverbedlosssc_run(
		int * ierr,
		int *nval, int *nconfig, int * ninputs,int * npar,
	    int * nstates,
		double * config,
		double * inputs_array,
		double * parameters,
		double * states_array,
		double * sacpar,
		double * sacinitpar,
		double * riverbedstore_ts,
		int * use_riverbedstore_ts,
		double * initial_riverbedstore,
		double * riverbedstore_fact_ts,
		int * use_riverbedstore_fact_ts,
		double * initial_riverbedstore_fact,
		double * uztwc,
		double * uzfwc,
		double * lztwc,
		double * lzfsc,
		double * lzfpc,
		double * adimc,
		double * sett,
		double * se1,
		double * se3,
		double * se4,
		double * se5,
		double * roimp,
		double * sdro,
		double * ssur,
		double * sif,
		double * bfp,
		double * bfs,
		double * bfcc,
		int * use_RR,
		int * excl_channel,
		double * max_RBS_mass_balance)
{
  
  int i,k,debugflag,imax_debug;
	double states[NSTATESNONROUTING+NINFLOWMAX*NSTATESROUTING];
	//double* states=malloc(sizeof(double)*(NSTATESNONROUTING+NINFLOWMAX*NSTATESROUTING));
	
	double subcat_area;
	int n=*nval;
	double etmult=1;
	double dt=1;
	double U[*nval];
	//double uztwc_0,uzfwc_0,lztwc_0,lzfsc_0,lzfpc_0,adimc_0;
	int min_ninc=20;
	double state_S_ts=1e-6, state_S2_ts=1e-6, state_S3_ts=1e-6;
	int use_state_S_ts=0, use_state_S2_ts=0, use_state_S3_ts=0;
	
	double rainfall_river[*nval],evap_river[*nval];
	double sma_sac_vol[*nval];
	double sma_sac_TET[*nval];
	
	double timestep_length=config[0];
	double river_length=config[1];
	double riverbedstore;
	double riverbedstore_fact;
	
	//printf("riverbedlosssc_run\n");
	
	subcat_area = config[3];
	
	debugflag=0;

	if(debugflag>0) 
		printf("\n\t\tnval=%d\n\t\tnconfig=%d\n\t\tninputs=%d\n\t\tnpar=%d\n\t\tnstates=%d\n",
			*nval,*nconfig,*ninputs,*npar,*nstates);
  
	// Initialise states
	for(i=0;i<NSTATESNONROUTING+NINFLOWMAX*NSTATESROUTING;i++) states[i] = states_array[i];
	
	double RBS_mass_balance_error = 0.0;
	
	  // Loop through time series
    for(i=0;i<*nval;i++){
      //printf("\n\t\triverbedlosssc_run : time step %d\n",i);
  		if((debugflag>1) & (i < imax_debug)){ 
  			printf("\n\t\triverbedlosssc_run : time step %d\n",i);
  			printf("\n\t\t-> config = ");
  			for(k=0;k<*nconfig;k++) printf("%0.2f ",config[k]);
  			printf("\n\t\t-> parameters = ");
  			for(k=0;k<*npar;k++) printf("%0.2f ",parameters[k]);
  			printf("\n\t\t-> inputs = ");
  			for(k=0;k<*ninputs;k++) printf("%0.2f ",inputs_array[i* *ninputs+k]);
  		}
  		
  		riverbedstore_fact = -1e6;
  		if (i > 0) {
  		  if (*use_riverbedstore_ts > 0){
  		    if (riverbedstore_ts[i-1] != -1e6) {
  		      
  		      double prev_RBS_mass_balance_error = RBS_mass_balance_error;
  		      RBS_mass_balance_error = prev_RBS_mass_balance_error + ((riverbedstore_ts[i-1] * river_length) - riverbedstore)/river_length;
  		      
  		      if(*max_RBS_mass_balance != -1e6){
  		        if(RBS_mass_balance_error > *max_RBS_mass_balance){
  		          riverbedstore = (*max_RBS_mass_balance - prev_RBS_mass_balance_error) * river_length + riverbedstore;
  		          RBS_mass_balance_error = *max_RBS_mass_balance;
  		        } else if(RBS_mass_balance_error < -*max_RBS_mass_balance){
  		          riverbedstore = (-*max_RBS_mass_balance - prev_RBS_mass_balance_error) * river_length + riverbedstore;
  		          RBS_mass_balance_error = -*max_RBS_mass_balance;
  		        } else {
  		          riverbedstore = riverbedstore_ts[i-1] * river_length;
  		        }
  		      } else {
  		        riverbedstore = riverbedstore_ts[i-1] * river_length;
  		      }
  		      
  		    }
  		    states[12] = RBS_mass_balance_error;
  		    
  		  }
  		  if (*use_riverbedstore_fact_ts > 0){
  		    if (riverbedstore_fact_ts[i-1] != -1e6) {
  		      riverbedstore_fact = riverbedstore_fact_ts[i-1];
  		    }
  		  }
  		} else {
  		  riverbedstore = *initial_riverbedstore * river_length;
  		  if(*initial_riverbedstore_fact != -1e6){
  		    riverbedstore_fact = *initial_riverbedstore_fact;
  		  }
  		  
  		}
  		riverbedlosssc_runtimestep(
  			ierr,
  			nconfig,
  			ninputs,
  			npar,
  			nstates,
  			config,
  			&(inputs_array[i* *ninputs]),
  			parameters,
  			states,
  			sacpar,
  			sacinitpar,
  			&riverbedstore,
  			&riverbedstore_fact,
  			max_RBS_mass_balance);
	
		// store states
		for(k=0;k<*nstates;k++){
		  states_array[i* *nstates+k] = states[k];
		  //printf("states[k]=%0.2f\n ",states[k]);
		}
		//free(states);
			
		if((debugflag>1) & (i < imax_debug))
		{ 
			printf("\n\t\t-> states (non routing) = ");
			for(k=0;k<NSTATESNONROUTING;k++) 
					printf("%0.2f (%d) ",states[k],k);
			
			if(debugflag>2)
			{
				printf("\n\t\t-> states (routing-nonlag) = ");
				for(k=NSTATESNONROUTING;k<4+NSTATESNONROUTING;k++) 
						printf("%0.2f (%d) ",states[k],k);

				printf("\n\t\t-> states (routing-lag) = ");
				for(k=NSTATESNONROUTING+4;k<NSTATESROUTING+NSTATESNONROUTING;k++) 
						printf("%0.2f (%d) ",states[k],k);
			}
			printf("\n\n");
		}

  	if(*ierr>0) return;    
  
	}
  
  
  
  
  // get rainfall and evap
  for(i=0;i<*nval;i++){
    rainfall_river[i] = inputs_array[i* *ninputs];
    evap_river[i] = inputs_array[i* *ninputs + 1];
    U[i] = 0;
  }
  
  if(*use_RR>0)
  {
    sma_sac_state(rainfall_river, evap_river, &n, sacpar, &etmult,
                  &dt, U,
                  uztwc, uzfwc,
                  lztwc, lzfsc, lzfpc, adimc,
                  sett, se1, se3, se4,
                  se5, roimp, sdro, ssur,
                  sif, bfp, bfs, bfcc,
                  &sacinitpar[0],&sacinitpar[1],&sacinitpar[2],&sacinitpar[3],&sacinitpar[4],&sacinitpar[5],
                  &min_ninc,
                  &state_S_ts, &use_state_S_ts,
                  &state_S2_ts, &use_state_S2_ts,
                  &state_S3_ts, &use_state_S3_ts);
                  
    for(i=0;i<*nval;i++){
      if(*excl_channel>0){
        sma_sac_vol[i] = (roimp[i]+sdro[i]+ssur[i]+sif[i]+bfcc[i])/1000*subcat_area/timestep_length;
      } else {
        sma_sac_vol[i] = U[i]/1000*subcat_area/timestep_length;
      }
      
      states_array[i* *nstates] += sma_sac_vol[i]; // add sac to outflow
      states_array[i* *nstates + NSTATESNONROUTING - 4] = sma_sac_vol[i]; // sac output
      sma_sac_TET[i] = sett[i]/1000*subcat_area/timestep_length;
      states_array[i* *nstates + NSTATESNONROUTING - 2] = sma_sac_TET[i]; // sac output
      
      //printf("timestep,sma_sac_vol[i],states_array[i* *nstates + NSTATESNONROUTING - 1] = %d,%0.2f,%0.2f \n",i,sma_sac_vol[i],states_array[i* *nstates + NSTATESNONROUTING - 1]);
    }
                                                                                                        
  }


}

#include "header.h"

/******************************************************************************
Code to loop through subcatchments in section

-- inputs ---
int * ierr,			Error code

int *nval, 			Dimensions of array (nval = number of timesteps)
int *nconfig, 
int * ninputs,
int * npar,
int * nstates, 

double * config, 		configuration data
double * inputs_array,	inputs data
double * parameters,	model parameters
double * states_array,	model states


-- states ---
double * states_array 	model states

-- Authors --
Shaun Kim, CSIRO


*******************************************************************************/


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
    double * max_RBS_mass_balance
    );

void riverbedlosssection_run(
		int * ierr,
		int *nval, int *nconfig, int * ninputs,int * npar, int * nsubcat,
	    int * nstates,
		double * config_array,
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
	//double states[NSTATESNONROUTING+NINFLOWMAX*NSTATESROUTING];
	

	debugflag=0;
	//printf("riverbedlosssection_run\n");

	if(debugflag>0) 
		printf("\n\t\tnval=%d\n\t\tnconfig=%d\n\t\tninputs=%d\n\t\tnpar=%d\n\t\tnstates=%d\n",
			*nval,*nconfig,*ninputs,*npar,*nstates);

	for(i=0;i<*nsubcat;i++){

	  
	  riverbedlosssc_run(
	    ierr,
	    nval, nconfig, ninputs, npar,
	    nstates,
	    &config_array[i* *nconfig],
	    &inputs_array[i* *ninputs* *nval],
	    parameters,
	    &states_array[i* *nstates* *nval],
	    sacpar,
	    sacinitpar,
	    riverbedstore_ts,
	    use_riverbedstore_ts,
	    initial_riverbedstore,
	    riverbedstore_fact_ts,
	    use_riverbedstore_fact_ts,
	    initial_riverbedstore_fact,
	    uztwc,
	    uzfwc,
	    lztwc,
	    lzfsc,
	    lzfpc,
	    adimc,
	    sett,
	    se1,
	    se3,
	    se4,
	    se5,
	    roimp,
	    sdro,
	    ssur,
	    sif,
	    bfp,
	    bfs,
	    bfcc,
	    use_RR,
	    excl_channel,
	    max_RBS_mass_balance);
	  

	}
	
	



}

#include "header.h"

/******************************************************************************
Code to loop through each reach section

-- inputs ---
int * ierr,			Error code

int *nval, 			Dimensions of array (nval = number of timesteps)
int *nconfig, 
int * ninputs,
int * npar,
int * nstates, 

double * config, 		configuration data (see awrarirrig_runtimestep.c)
double * inputs_array,	inputs data (see awrarirrig_runtimestep.c)
double * parameters,	model parameters (see awrarirrig_runtimestep.c)
double * states_array,	model states (see awrarirrig_runtimestep.c)

-- states ---
double * states_array 	model states (see awrarirrig_runtimestep.c)

-- Authors --
Shaun Kim, CSIRO



*******************************************************************************/


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
    double * max_RBS_mass_balance);

void riverbedlossreach_run(
		int * ierr,
		int *nval, int *nconfig, int * ninputs,int * npar,int * nsections,int * nsubcat,
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
		double * outflow_total,
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
		)
{
  
  	int i,k,debugflag,imax_debug;
	//double states[NSTATESNONROUTING+NINFLOWMAX*NSTATESROUTING];
	int config_counter = 0;
	int inputs_counter = 0;
	int states_counter = 0;
	int ii;
	int iii;
	
  //double outflow_total[*nval];
	
	debugflag=0;

	if(debugflag>0) 
		printf("\n\t\tnval=%d\n\t\tnconfig=%d\n\t\tninputs=%d\n\t\tnpar=%d\n\t\tnstates=%d\n",
			*nval,*nconfig,*ninputs,*npar,*nstates);

	//printf("riverbedlossreach_run\n");
	
	for(i=0;i<*nsections;i++){
	  
	  
	  riverbedlosssection_run(
	    ierr,
	    nval, nconfig, ninputs,npar, &nsubcat[i],
	    nstates,
	    &config_array[config_counter],
	    &inputs_array[inputs_counter],
	    parameters,
	    &states_array[states_counter],
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
	  
	  
	  // get discharge for the section
	  for(ii=0;ii<*nval;ii++) outflow_total[ii] = 0;
	  
	  for(iii=0;iii<nsubcat[i];iii++){
	    for(ii=0;ii<*nval;ii++){
	      //outflow_total[ii] += states_array[*nstates* *nval* iii + ii* *nstates];
	      outflow_total[ii] += states_array[states_counter + *nstates* *nval* iii + ii* *nstates];
	      //printf("subcat,timestep,statesarrayindex,outflowtotal = %d,%d,%d,%0.2f \n",iii,ii,*nstates* *nval* iii + ii* *nstates,outflow_total[ii]);
	      
	    }
	  }

	  
	  config_counter += nsubcat[i]* *nconfig;
	  inputs_counter += nsubcat[i]* *ninputs* *nval;
	  states_counter += nsubcat[i]* *nstates* *nval;
	  
	  
    // replacing (overwrites) inflow for the next section (first inflow) downstream
    if(i+1 != *nsections){
      for(ii=0;ii<*nval;ii++){
        inputs_array[inputs_counter + 2 + ii* *ninputs] = outflow_total[ii];
      }
    }

	  

	 
	}
	
	//printf("finished \n");


}

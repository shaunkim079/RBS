//#include "stdafx.h"
#include "header.h"
#include "sacramento.h"
#include <R.h>

/******************************************************************************
Timestep code of river bed/bank store reach model

-- inputs ---
int * ierr,				Error code

int * nconfig, 			dimensions of arrays
int * ninputs,
int * npar,
int * nstates,

double * config       Configuration data

	
double * inputs,      input data. 1D Array ninputs x 1


double * parameters


double * states,
 there are 24 states types (columns) for each subcatment
 1: discharge
 2: Smax
 3: river loss (m^3/timestep)
 4: river bed/bank store (m^3)
 5: river bed/bank store drying (m^3/timestep)
 6: river rain (m^3/timestep)
 7: river evap (m^3/timestep)
 8: river drying (m^3/timestep)
 9: river volume
 10: runoff (m^3/s)
 11: river bed/bank store leakage (m^3/timestep)
 12: rainfall runoff total ET (m^s/s)
 remaining 12 states are for routing equations
	
-- Authors --
Shaun Kim, CSIRO




*******************************************************************************/


/******* Sub routines**********************************************************/
double max(double v1,double v2);
double minmax(double min,double max,double input);
double monod(double p1,double p2,double input);



void muskingum_solve_runtimestep(int * ierr,
                                 double dt,
                                 double inflow,
                                 double K,
                                 double x,
                                 double * routing_volume,
                                 double * routing_outflow);

void lag_runtimestep(int * ierr,
		double dt,
		double inflow,
		double lag,
		double * uh
	);



/******* Main code ************************************************************/
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
		double * max_RBS_mass_balance)
{
  
    // Dimensions
    int i,k,ninflows,debugflag;

    // config data
	double timestep_length,river_length;

    // parameters
	double lag,K,x;

    // input
	double rainfall_river=0, evap_river=0,inflow=0;

   // Other variables
   double river_area,flux;
   
   double subcat_area;
   double Smax,fc,max_rivbedstore,lossK;
   double riverloss;
   
   
   //int n=1;
   double etmult=1;
   //double dt=1;
   double U=0;
   //double uztwc_0,uzfwc_0,lztwc_0,lzfsc_0,lzfpc_0,adimc_0;
   int min_ninc=20;
   double state_S_ts=1e-6, state_S2_ts=1e-6, state_S3_ts=1e-6;
   int use_state_S_ts=0, use_state_S2_ts=0, use_state_S3_ts=0;
   
   double lag_rout_sum, total_rout_sum;
   double rivbedstdrying;
   double wetAreaFactor;
   
   double invv, alpha;
   double riverdrying;
   int lag_max;
   double store_norm;
   double total_rout_deficit;
   double max_rout_store;
   double excess_rout;

    debugflag = 0;
    //printf("riverbedlosssc_runtimestep\n");
    // Check dimensions
    ninflows		= *ninputs-NINPUTSNONROUTING;
    if(ninflows != 1)
    {
        *ierr = 50100101;
        return;
    }
    if(*nstates != NSTATESNONROUTING + NSTATESROUTING*ninflows)
    {
        *ierr = 50100102;
        return;
    }
    //if(*npar != NPARNONROUTING + NPARROUTING*ninflows)
    //{
    //    *ierr = 50100103;
    //    return;
    //}
    if(*nconfig < NCONFIG)
    {
        *ierr = 50100104;
        return;
    }

  // Set config
  timestep_length        = config[0];
  //timestep_length = minmax(1,LARGEVALUE,timestep_length);
  river_length       = config[1];
  //river_length = minmax(1,LARGEVALUE,river_length);
  river_area      = config[2];
  subcat_area = config[3];
   
	// Set parameters
	// IMPORTANT !! assume parameters are per unit river length!!
  
	invv = parameters[0];
	alpha = parameters[1];
	
	lag 	= (1-alpha)*invv*river_length; // Lag in sec
	K 	= alpha*invv*river_length; // Muskingum K
	x 	= 0;
	//printf("lag=%0.2f\n ",lag);
	//printf("K=%0.2f\n ",K);
	
	lag_max = ceil(lag/timestep_length)+1;
	if(lag_max==1) lag_max = 2;
	if(lag_max > (NMAXLAG-1)) lag_max = NMAXLAG-1;
	
	int num_rout_states;
	if(K>0){
	  num_rout_states = lag_max;
	} else {
	  num_rout_states = lag_max-1;
	}
	//lag_max = NMAXLAG;
	//printf("   lag_max=%d\n ",lag_max);
	
	//Smax = parameters[2]*river_length;
	wetAreaFactor = parameters[2];
	fc = parameters[3]*river_length;
	max_rivbedstore = parameters[4]*river_length;
	lossK = parameters[5];
  
  // Set inputs apart from inflows
  rainfall_river        = inputs[0];
  //rainfall_river = minmax(0,5e2,rainfall_river);

  evap_river            = inputs[1];
  //evap_river = minmax(0,5e1,evap_river);
  
  // mass balance:
  double start_rout_storage  = 0;
  for(i=1;i<NMAXLAG;i++){
    start_rout_storage += states[2+NSTATESNONROUTING+i] * timestep_length;
    //printf("i=%d lag_rout_sum=%0.2f\n",i,lag_rout_sum);
    //printf("  states[2+NSTATESNONROUTING+i]=%0.2f\n",states[2+NSTATESNONROUTING+i]);
  }
  start_rout_storage += states[NSTATESNONROUTING];

     // initialise outflow
    states[0] = 0;
		
		
		//lag_rout_sum = 0;
		//for(i=1;i<NMAXLAG;i++){
		//  lag_rout_sum += states[2+NSTATESNONROUTING+i];
		//}
		//total_rout_sum = lag_rout_sum + states[NSTATESNONROUTING];
		//printf("total_rout_sum_0=%0.2f\n ",total_rout_sum);
		
		//for(i=0;i<NMAXLAG;i++){
		//  printf("states[2+NSTATESNONROUTING+i]=%0.2f \n",states[2+NSTATESNONROUTING+i]);
		//}
		//printf("states[NSTATESNONROUTING]=%0.2f \n",states[NSTATESNONROUTING]);
		
    // Loop through tributaries - route inflows
    // (ninflow x current inflow values + ninflow x previous inflow)

		inflow = 0;
		if(inputs[NINPUTSNONROUTING]>0)
				inflow = inputs[NINPUTSNONROUTING];
    //printf("inflow=%0.2f\n ",inflow);

    //for(i=0;i<NMAXLAG;i++){
    //  printf("states[2+NSTATESNONROUTING+i]=%0.2f \n",states[2+NSTATESNONROUTING+i]);
    //}
    //printf("states[NSTATESNONROUTING]=%0.2f \n",states[NSTATESNONROUTING]);
    
			// lag subroutine
			/*
			double check_lag_storage_before  = 0;
			for(i=1;i<NMAXLAG;i++){
			  check_lag_storage_before += states[2+NSTATESNONROUTING+i] * timestep_length;
			  //printf("i=%d lag_rout_sum=%0.2f\n",i,lag_rout_sum);
			  //printf("  states[2+NSTATESNONROUTING+i]=%0.2f\n",states[2+NSTATESNONROUTING+i]);
			}
			double check_lag_inflow = inflow * timestep_length;
			*/
			lag_runtimestep(ierr,timestep_length,
				inflow,lag,
				&(states[2+NSTATESNONROUTING]));
			/*
			double check_lag_storage_after  = 0;
			for(i=1;i<NMAXLAG;i++){
			  check_lag_storage_after += states[2+NSTATESNONROUTING+i] * timestep_length;
			  //printf("i=%d lag_rout_sum=%0.2f\n",i,lag_rout_sum);
			  //printf("  states[2+NSTATESNONROUTING+i]=%0.2f\n",states[2+NSTATESNONROUTING+i]);
			}
			double check_lag_outflow = states[2+NSTATESNONROUTING] * timestep_length;
			*/
			
			//for(i=0;i<NMAXLAG;i++){
			//  printf("states[2+NSTATESNONROUTING+i]=%0.2f \n",states[2+NSTATESNONROUTING+i]);
			//}
			//printf("states[NSTATESNONROUTING]=%0.2f \n",states[NSTATESNONROUTING]);
			
			// Routing subroutine
			/*
			double check_musk_inflow = states[2+NSTATESNONROUTING] * timestep_length;
			double check_musk_storage_before = states[NSTATESNONROUTING];
			*/
			if(K>0.0){
			  muskingum_solve_runtimestep(ierr,
                                 timestep_length,
                                 states[2+NSTATESNONROUTING],
                                       K,
                                       x,
                                       &(states[NSTATESNONROUTING]),
                                       &(states[1+NSTATESNONROUTING]));
			} else {
			  states[NSTATESNONROUTING] = 0.0;
			  states[1+NSTATESNONROUTING] = states[2+NSTATESNONROUTING];
			}
			/*
      double check_musk_storage_after = states[NSTATESNONROUTING];
      double check_musk_outflow = states[1+NSTATESNONROUTING] * timestep_length;
      */
  	  //for(i=0;i<NMAXLAG;i++){
  	  //  printf("states[2+NSTATESNONROUTING+i]=%0.2f \n",states[2+NSTATESNONROUTING+i]);
  	  //}
  	  //printf("states[NSTATESNONROUTING]=%0.2f \n",states[NSTATESNONROUTING]);

			// Sum routing contributions
			states[0]+=states[1+NSTATESNONROUTING];
			//states[25]+=(inflow-states[1+NSTATESNONROUTING])*timestep_length; // river volume
			//if(states[25]<0.0)
			//{
			//	states[0]+=states[25];
			//	states[25] = 0.0;
			//}
			

	// rain evap and loss model need to work dynmically with river volume
	// total river volume is divided up into UH and muskingum routing store
	// UH: &(states[2+NSTATESNONROUTING]) with length NMAXLAG
  // muskingum: &(states[NSTATESNONROUTING])
  // will remove rain evap and loss proportionally
  // Important: UH states are in m3/s so to get volumes we multiply these by timestep length
  lag_rout_sum = 0;
	for(i=1;i<NMAXLAG;i++){
	  lag_rout_sum += states[2+NSTATESNONROUTING+i] * timestep_length;
	  //printf("i=%d lag_rout_sum=%0.2f\n",i,lag_rout_sum);
	  //printf("  states[2+NSTATESNONROUTING+i]=%0.2f\n",states[2+NSTATESNONROUTING+i]);
	}
	total_rout_sum = lag_rout_sum + states[NSTATESNONROUTING];
	//printf("lag_rout_sum=%0.2f states[NSTATESNONROUTING]=%0.2f total_rout_sum=%0.2f\n ",lag_rout_sum,states[NSTATESNONROUTING],lag_rout_sum);
	//printf("states[NSTATESNONROUTING]=%0.2f\n ",states[NSTATESNONROUTING]);
	//printf("total_rout_sum_1=%0.2f\n ",total_rout_sum);
	
	
	// Rainfall and evap fluxes
	states[5] = rainfall_river * 1e-3 *  river_area;
	states[6] = evap_river * 1e-3 * river_area;
	flux = states[6]-states[5];

	// river loss calculations
	Smax = total_rout_sum * wetAreaFactor;
	//printf("total_rout_sum_1=%0.2f wetAreaFactor=%0.2f Smax=%0.2f\n ",total_rout_sum,wetAreaFactor,Smax);
	
	if(*riverbedstore_fact != -1e6){
	  *riverbedstore = *riverbedstore_fact * Smax;
	}
	
	//mass balance:
	double start_riverbedstore = *riverbedstore;
	
	//riverloss = 0;
	//if(Smax>*riverbedstore){
	//  riverloss = (Smax-*riverbedstore) * (1-exp(-lossK));
	//}

	
	
	riverloss = (Smax-*riverbedstore) * (1-exp(-lossK));
	//double orig_riverbedstore = *riverbedstore;
	//printf("Smax=%0.2f *riverbedstore=%0.2f lossK=%0.2f riverloss=%0.2f\n",Smax,*riverbedstore,lossK,riverloss);
	//printf("riverbedstore=%0.2f\n",riverbedstore);
	//if(riverloss<0){
	//  printf("riverloss=%0.2f\n",riverloss);
	//}
	//printf("total_rout_sum=%0.2f\n",total_rout_sum);
	//printf("states[NSTATESNONROUTING]=%0.2f\n",states[NSTATESNONROUTING]);
	double predicted_Smax = (total_rout_sum-riverloss) * wetAreaFactor;
	double predicted_riverbedstore = *riverbedstore+riverloss;
	//printf("riverloss=%0.2f\n",riverloss);
	
	// this solves so that estimated river loss does not increase Smax above RBS ie RBS should instead equal Smax 
	if(riverloss<0){
	  if(predicted_riverbedstore<predicted_Smax){
	    //riverloss = riverloss * 0.5;
	    riverloss = (total_rout_sum * wetAreaFactor - *riverbedstore)/(1+wetAreaFactor);
	    //printf("riverloss=%0.2f\n",riverloss);
	    //predicted_Smax = (total_rout_sum-riverloss) * wetAreaFactor;
	    //predicted_riverbedstore = *riverbedstore+riverloss;
	    //printf("1 predicted_Smax=%0.2f predicted_riverbedstore=%0.2f\n",predicted_Smax,predicted_riverbedstore);
	  }
	} else {
	  if(predicted_riverbedstore>predicted_Smax){
	    //riverloss = riverloss * 0.5;
	    riverloss = (total_rout_sum * wetAreaFactor - *riverbedstore)/(1+wetAreaFactor);
	    //printf("riverloss=%0.2f\n",riverloss);
	    //predicted_Smax = (total_rout_sum-riverloss) * wetAreaFactor;
	    //predicted_riverbedstore = *riverbedstore+riverloss;
	    //printf("2 predicted_Smax=%0.2f predicted_riverbedstore=%0.2f\n",predicted_Smax,predicted_riverbedstore);
	  }
	}
	
	
	if(total_rout_sum>riverloss){
	  if(*riverbedstore+riverloss<0){
	    riverloss = -*riverbedstore;
	    *riverbedstore = 0;
	  } else {
	    *riverbedstore = *riverbedstore+riverloss;
	  }
	  // apportion river volume loss to different routing states
	  if(total_rout_sum>0){
	    /*
	     for(i=1;i<NMAXLAG;i++){
	     states[2+NSTATESNONROUTING+i] -= (riverloss * ((states[2+NSTATESNONROUTING+i]*timestep_length)/total_rout_sum))/timestep_length;
	     }
	     states[NSTATESNONROUTING] -= riverloss * (states[NSTATESNONROUTING]/total_rout_sum);
	     */
	    // this apportions differently for gains
	    if(riverloss>0){
	      double total_rivloss_sum = 0;
	      for(i=1;i<NMAXLAG;i++){
	        total_rivloss_sum = total_rivloss_sum + (riverloss * ((states[2+NSTATESNONROUTING+i]*timestep_length)/total_rout_sum));
	        states[2+NSTATESNONROUTING+i] -= (riverloss * ((states[2+NSTATESNONROUTING+i]*timestep_length)/total_rout_sum))/timestep_length;
	      }
	      total_rivloss_sum = total_rivloss_sum + (riverloss * (states[NSTATESNONROUTING]/total_rout_sum));
	      states[NSTATESNONROUTING] -= riverloss * (states[NSTATESNONROUTING]/total_rout_sum);
	      
	      // printf("check 1 total_rivloss_sum=%0.2f riverloss=%0.2f\n ",total_rivloss_sum,riverloss);
	    } else {
	      // giving higher gains to lower sized river stores
	      max_rout_store = 0;
	      for(i=1;i<lag_max;i++){
	        //printf("i=%d max_rout_store=%0.2f states[2+NSTATESNONROUTING+i]*timestep_length=%0.2f\n ",i,max_rout_store,states[2+NSTATESNONROUTING+i]*timestep_length);
	        max_rout_store = max(states[2+NSTATESNONROUTING+i]*timestep_length,max_rout_store);
	        //printf("max_rout_store=%0.2f\n ",max_rout_store);
	      }
	      max_rout_store = max(states[NSTATESNONROUTING],max_rout_store);
	      //printf("max_rout_store=%0.2f states[NSTATESNONROUTING]=%0.2f\n ",max_rout_store,states[NSTATESNONROUTING]);
	      
	      total_rout_deficit = 0;
	      for(i=1;i<lag_max;i++){
	        total_rout_deficit = total_rout_deficit+(max_rout_store-(states[2+NSTATESNONROUTING+i]*timestep_length));
	      }
	      if(K>0) total_rout_deficit = total_rout_deficit+(max_rout_store-states[NSTATESNONROUTING]);
	      
	      if(-riverloss<=total_rout_deficit){
	        //printf("small gain\n");
	        double total_rivloss_sum = 0;
	        for(i=1;i<lag_max;i++){
	          total_rivloss_sum += (riverloss * ((max_rout_store-(states[2+NSTATESNONROUTING+i]*timestep_length))/total_rout_deficit));
	          states[2+NSTATESNONROUTING+i] -= (riverloss * ((max_rout_store-(states[2+NSTATESNONROUTING+i]*timestep_length))/total_rout_deficit))/timestep_length;
	          //printf("i=%d,states[2+NSTATESNONROUTING+i]=%0.2f\n ",i,states[2+NSTATESNONROUTING+i]);
	          //printf("i=%d ((max_rout_store-(states[2+NSTATESNONROUTING+i]*timestep_length))/total_rout_deficit)=%0.2f\n ",i,((max_rout_store-(states[2+NSTATESNONROUTING+i]*timestep_length))/total_rout_deficit));
	          //printf("total_rout_sum=%0.2f\n ",total_rout_sum);
	          //printf("(states[2+NSTATESNONROUTING+i]*timestep_length)=%0.2f\n ",(states[2+NSTATESNONROUTING+i]*timestep_length));
	        }
	        if(K>0){
	          total_rivloss_sum += riverloss * ((max_rout_store-states[NSTATESNONROUTING])/total_rout_deficit);
	          // printf("check 2 total_rivloss_sum=%0.2f riverloss=%0.2f\n ",total_rivloss_sum,riverloss);
	          
	          states[NSTATESNONROUTING] -= riverloss * ((max_rout_store-states[NSTATESNONROUTING])/total_rout_deficit);
	          //printf("states[NSTATESNONROUTING]=%0.2f\n ",states[NSTATESNONROUTING]);
	        }

	        
	      } else {
	        //printf("large gain\n");
	        //printf("Smax=%0.2f orig_riverbedstore=%0.2f lossK=%0.2f riverloss=%0.2f\n",Smax,orig_riverbedstore,lossK,riverloss);
	        //printf("max_rout_store=%0.2f riverloss=%0.2f total_rout_deficit=%0.2f\n ",max_rout_store,riverloss,total_rout_deficit);
	        excess_rout = -riverloss - total_rout_deficit;
	        double total_rivloss_sum = 0;
	        for(i=1;i<lag_max;i++){
	          //printf("before i=%d states[2+NSTATESNONROUTING+i]=%0.2f\n ",i,states[2+NSTATESNONROUTING+i]);
	          total_rivloss_sum -= ((excess_rout/(num_rout_states) + max_rout_store)/timestep_length - states[2+NSTATESNONROUTING+i])*timestep_length;
	          states[2+NSTATESNONROUTING+i] = (excess_rout/(num_rout_states) + max_rout_store)/timestep_length;
	          //printf("after i=%d states[2+NSTATESNONROUTING+i]=%0.2f\n ",i,states[2+NSTATESNONROUTING+i]);
	        }
	        if(K>0){
	          //printf("before states[NSTATESNONROUTING]=%0.2f\n ",states[NSTATESNONROUTING]);
	          total_rivloss_sum -= (excess_rout/(num_rout_states) + max_rout_store) - states[NSTATESNONROUTING];
	          // printf("check 3 total_rivloss_sum=%0.2f riverloss=%0.2f\n ",total_rivloss_sum,riverloss);
	          
	          states[NSTATESNONROUTING] = excess_rout/(num_rout_states) + max_rout_store;
	          //printf("after states[NSTATESNONROUTING]=%0.2f\n ",states[NSTATESNONROUTING]);
	        }

	      }
	      
	      /*
	       store_norm=0;
	       for(i=1;i<lag_max;i++){
	       store_norm = store_norm + (total_rout_sum-(states[2+NSTATESNONROUTING+i]*timestep_length));
	       }
	       store_norm = store_norm + (total_rout_sum-states[NSTATESNONROUTING]);
	       
	       for(i=1;i<lag_max;i++){
	       //states[2+NSTATESNONROUTING+i] -= (riverloss * (1-((states[2+NSTATESNONROUTING+i]*timestep_length)/total_rout_sum)))/timestep_length;
	       //printf("i=%d (1-((states[2+NSTATESNONROUTING+i]*timestep_length)/total_rout_sum))=%0.2f\n ",i,(1-((states[2+NSTATESNONROUTING+i]*timestep_length)/total_rout_sum)));
	       states[2+NSTATESNONROUTING+i] -= (riverloss * ((total_rout_sum-(states[2+NSTATESNONROUTING+i]*timestep_length))/store_norm))/timestep_length;
	       //printf("i=%d ((total_rout_sum-(states[2+NSTATESNONROUTING+i]*timestep_length))/store_norm)=%0.2f\n ",i,((total_rout_sum-(states[2+NSTATESNONROUTING+i]*timestep_length))/store_norm));
	       //printf("total_rout_sum=%0.2f\n ",total_rout_sum);
	       //printf("(states[2+NSTATESNONROUTING+i]*timestep_length)=%0.2f\n ",(states[2+NSTATESNONROUTING+i]*timestep_length));
	       }
	       //states[NSTATESNONROUTING] -= riverloss * (1-(states[NSTATESNONROUTING]/total_rout_sum));
	       states[NSTATESNONROUTING] -= riverloss * ((total_rout_sum-states[NSTATESNONROUTING])/store_norm);
	       //printf("((total_rout_sum-states[NSTATESNONROUTING])/store_norm)=%0.2f\n ",((total_rout_sum-states[NSTATESNONROUTING])/store_norm));
	       */
	    }
	  } else {
	    // avoid divide by zero
	    double total_rivloss_sum = 0;
	    for(i=1;i<lag_max;i++){
	      total_rivloss_sum += (riverloss/num_rout_states);
	      states[2+NSTATESNONROUTING+i] -= (riverloss/num_rout_states)/timestep_length;
	      //printf("states[2+NSTATESNONROUTING+i]=%0.2f\n ",states[2+NSTATESNONROUTING+i]);
	    }
	    if(K>0){
	      total_rivloss_sum += riverloss/num_rout_states;
	      // printf("check 4 total_rivloss_sum=%0.2f riverloss=%0.2f\n ",total_rivloss_sum,riverloss);
	      
	      states[NSTATESNONROUTING] -= riverloss/num_rout_states;
	      //printf("lag_max=%d\n ",lag_max);
	    }

	    
	  }
	  
	  total_rout_sum -= riverloss;
	} else {
	  riverloss = total_rout_sum;
	  for(i=1;i<NMAXLAG;i++){
	    states[2+NSTATESNONROUTING+i] = 0;
	  }
	  states[NSTATESNONROUTING] = 0;
	  total_rout_sum = 0;
	  
	  *riverbedstore = *riverbedstore+riverloss;

	}
	
	/*
	 if(riverloss<0){
	 if(*riverbedstore<(total_rout_sum * wetAreaFactor)){
	 //printf("riverloss=%0.2f\n",riverloss);
	 //printf("riverbedstore low\n");
	 }
	 } else {
	 if(*riverbedstore>(total_rout_sum * wetAreaFactor)){
	 //printf("riverloss=%0.2f\n",riverloss);
	 //printf("riverbedstore high\n");
	 }
	 }
	 */
	
	
	states[8] = total_rout_sum;
	
  //printf("total_rout_sum=%0.2f\n",total_rout_sum);
  //printf("states[NSTATESNONROUTING]=%0.2f\n",states[NSTATESNONROUTING]);

  states[10] = fc*timestep_length;
  //ensure riverbedstore doesn't go below 0
  if((*riverbedstore - states[10])<0){
    states[10] = *riverbedstore;
    *riverbedstore = 0;
  } else {
    *riverbedstore -= states[10];
  }
  
  //drying out the soil store when river volume is zero
  //if(excessflux>0){
  //  rivbedstdrying = riverbedstore/Smax * excessflux * drycoeff;
  //  riverbedstore -= rivbedstdrying;
  //  //ensure S doesn't go below 0
  //  if(riverbedstore<0){
  //    riverbedstore = 0;
  //  }
  //}

  //printf("riverbedstore=%0.2f\n ",riverbedstore);
  //printf("Smax=%0.2f\n ",Smax);
  // if(flux>0 && *riverbedstore>Smax){
  //   //printf("here\n ");
  //   rivbedstdrying = (1 - Smax/ *riverbedstore) * flux * drycoeff;
  //   riverdrying = flux - rivbedstdrying;
  //   
  // } else {
  //   rivbedstdrying = 0;
  //   riverdrying = flux;
  // }
  
  if(*riverbedstore>=max_rivbedstore){
    rivbedstdrying = evap_river * 1e-3 * subcat_area;
  } else {
    rivbedstdrying = evap_river * 1e-3 * subcat_area * (*riverbedstore/max_rivbedstore);
  }
  
  riverdrying = flux;
  
  //ensure riverbedstore doesn't go below 0
  if(*riverbedstore<rivbedstdrying){
    rivbedstdrying = *riverbedstore;
    *riverbedstore = 0;
  } else {
    *riverbedstore -= rivbedstdrying;
  }
  
  if(riverdrying<total_rout_sum){
    //printf("here \n");

    if(total_rout_sum > 0){
      //printf("here 1 \n");
      // apportion
      for(i=1;i<NMAXLAG;i++){
        states[2+NSTATESNONROUTING+i] -= (riverdrying * ((states[2+NSTATESNONROUTING+i]*timestep_length)/total_rout_sum))/timestep_length;
      }
      states[NSTATESNONROUTING] -= riverdrying * (states[NSTATESNONROUTING]/total_rout_sum);
      total_rout_sum -= riverdrying;
    } else {
      //printf("here 2 \n");

      // apportion rainfall evenly to avoid dividing by zero
      for(i=1;i<lag_max;i++){
        states[2+NSTATESNONROUTING+i] -= (riverdrying/num_rout_states)/timestep_length;
        //printf("states[2+NSTATESNONROUTING+i]=%0.2f\n ",states[2+NSTATESNONROUTING+i]);
      }
      if(K>0) states[NSTATESNONROUTING] -= riverdrying/num_rout_states;
      total_rout_sum -= riverdrying;

    }


  }
  else
  {
    //excessflux = flux - total_rout_sum;
    riverdrying = total_rout_sum;
    //states[5] = rainfall_river * 1e-3 *  river_area;
    //states[6] = evap_river * 1e-3 * river_area;
    //flux = states[6]-states[5];
    states[6] = riverdrying + states[5];
    
    for(i=1;i<NMAXLAG;i++){
      states[2+NSTATESNONROUTING+i] = 0;
    }
    states[NSTATESNONROUTING] = 0;
    total_rout_sum = 0;
  }
  
  //for(i=1;i<NMAXLAG;i++){
  //  printf("states[2+NSTATESNONROUTING+i]=%0.2f \n",states[2+NSTATESNONROUTING+i]);
  //}
  //printf("states[NSTATESNONROUTING]=%0.2f \n",states[NSTATESNONROUTING]);
	
	//Sacramento here, adds to states[0]
	//sma_sac(&rainfall_river, &evap_river, &n, sacpar, &etmult,
  //       &dt, &U,
  //       &sacinitpar[0],&sacinitpar[1],&sacinitpar[2],&sacinitpar[3],&sacinitpar[4],&sacinitpar[5],
  //       &min_ninc,
  //       &state_S_ts, &use_state_S_ts,
  //       &state_S2_ts, &use_state_S2_ts,
  //       &state_S3_ts, &use_state_S3_ts);
  //states[0] += U;
	
	// printf("states[10]=%0.5f monod1=%0.5f monod2=%0.5f states[0]=%0.5f \n",
	// states[10],monod1,monod2,states[0]);
	
	//printf("river_outflow=%0.2f \n",states[0]);
	
	//printf("total_rout_sum=%0.2f\n ",total_rout_sum);
	

	states[1] = Smax;
  states[2] = riverloss;
  states[3] = *riverbedstore;
  states[4] = rivbedstdrying;
  states[7] = riverdrying;
  
  
  
  // mass balance:
  double final_rout_storage = 0;
  for(i=1;i<NMAXLAG;i++){
    final_rout_storage += states[2+NSTATESNONROUTING+i] * timestep_length;
    //printf("i=%d lag_rout_sum=%0.2f\n",i,lag_rout_sum);
    //printf("  states[2+NSTATESNONROUTING+i]=%0.2f\n",states[2+NSTATESNONROUTING+i]);
  }
  final_rout_storage = final_rout_storage + states[NSTATESNONROUTING];
  
  double mass_balance = start_rout_storage + start_riverbedstore - 
    final_rout_storage - *riverbedstore +
    states[5]+(inflow*timestep_length)-
    (states[0]*timestep_length)-
    states[6]-rivbedstdrying-states[10];
  if(mass_balance>1 || mass_balance<-1){
    
    printf("mass_balance=%0.2f\n ",mass_balance);
    printf("  start_rout_storage=%0.2f\n ",start_rout_storage);
    printf("  start_riverbedstore=%0.2f\n ",start_riverbedstore);
    printf("  final_rout_storage=%0.2f\n ",final_rout_storage);
    printf("  *riverbedstore=%0.2f\n ",*riverbedstore);
    printf("  rain=%0.2f\n ",states[5]);
    printf("  inflow=%0.2f\n ",(inflow*timestep_length));
    printf("  outflow=%0.2f\n ",states[0]*timestep_length);
    printf("  evap=%0.2f\n ",states[6]);
    printf("  rivbedstdrying=%0.2f\n ",rivbedstdrying);
    printf("  RBS_leakage=%0.2f\n ",states[10]);
    
    /*
    printf("\n" );
    printf("  check_lag_storage_before=%0.2f\n ",check_lag_storage_before);
    printf("  check_lag_storage_after=%0.2f\n ",check_lag_storage_after);
    printf("  check_lag_inflow=%0.2f\n ",check_lag_inflow);
    printf("  check_lag_outflow=%0.2f\n ",check_lag_outflow);
    printf("  lag mass balance=%0.2f\n ",check_lag_storage_before-check_lag_storage_after+check_lag_inflow-check_lag_outflow);
    */
    /*
    printf("\n" );
    printf("  check_musk_storage_before=%0.2f\n ",check_musk_storage_before);
    printf("  check_musk_storage_after=%0.2f\n ",check_musk_storage_after);
    printf("  check_musk_inflow=%0.2f\n ",check_musk_inflow);
    printf("  check_musk_outflow=%0.2f\n ",check_musk_outflow);
    printf("  musk mass balance=%0.2f\n ",check_musk_storage_before-check_musk_storage_after+check_musk_inflow-check_musk_outflow);
    */
  }
  
  //printf("mass_balance=%0.2f\n ",mass_balance);
  
  //printf("riverloss=%0.2f\n ",riverloss);
	
 	if(debugflag>1)
	{
		printf("\n\t\t-> config = ");
		for(k=0;k<*nconfig;k++) printf("%0.2f ",config[k]);

		printf("\n\t\t-> parameters = ");
		for(k=0;k<*npar;k++) printf("%0.2f ",parameters[k]);
		printf("\n\t\t-> inputs = ");
		for(k=0;k<*ninputs;k++) printf("%0.2f ",inputs[k]);

		printf("\n\t\t-> states (non routing) = ");
		for(k=0;k<NSTATESNONROUTING;k++)
				printf("%0.2f (%d) ",states[k],k);

		printf("\n\t\t-> states (routing-nonlag) = ");
		for(k=NSTATESNONROUTING;k<2+NSTATESNONROUTING;k++)
				printf("%0.2f (%d) ",states[k],k);

		printf("\n\t\t-> states (routing-lag) = ");
		for(k=NSTATESNONROUTING+4;k<NSTATESROUTING+NSTATESNONROUTING;k++)
				printf("%0.2f (%d) ",states[k],k);

		printf("\n\n");
	}
  
  return;
}

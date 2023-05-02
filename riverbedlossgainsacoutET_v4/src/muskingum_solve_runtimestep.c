#include "header.h"

/******************************************************************************
Muskingum routing component of AWRAR. This was modified by Shaun Kim to use a solver.
The model runs multiple inflows (ninflow) toward a downstream point.

dV/dt = I - O
V = k [xI + (1-x) O]
<=>
k [ x dI/dt + (1-x) dO/dt ] = I - O
<=>
dO/dt = I/(1-x)/k - x/(1-x) dI/dt - O/(1-x)/k

if I(t) = I0 (constant inflow), we get
dO/dt = I/(1-x)/k - O/(1-x)/k

if I(t) = (I1-I0)/d*t+I0 (linear inflow), we get
dO/dt = [ (I1-I0)/d*t+I0 ] /(1-x)/k - x/(1-x) (I1-I0)/d - O/(1-x)/k

-- inputs ---
int ierr							Error code
double dt,							time step length (sec)
double inflow,						inflows (m3/s).
double K							First muskingum parameter
double x							Second muskingum parameter

-- states ---
double * previous_inflow,			previous inflow (m3/s)
double * routing_volume,			volumes (m3)
double *outflow,					outflows (m3/s)

-- Authors --
Julien Lerat, CSIRO CLW
Justin Hughes, CSIRO CLW

-- Versions --
2012-10-23 - First translation of Justin's fortran code

-- Original fortran code --

c sample subroutine to include for R
SUBROUTINE rout(nrowdat,ntrib,outRout, musky, lagdat) 
implicit none
INTEGER nrowdat, ntrib, i, n
DOUBLE PRECISION outRout(nrowdat, ntrib+1), musky(ntrib+1,3), 
1 lagdat(nrowdat, ntrib+1)

outRout=0.D0
DO n=1,ntrib+1
outRout(1,n)=lagdat(1,n)
ENDDO

DO i=2,nrowdat
DO n=1,ntrib+1
outRout(i, n)=musky(n,1)*lagdat(i,n)+musky(n,2)*
1    lagdat(i-1,n)+musky(n,3)*outRout(i-1,n)
ENDDO
ENDDO
END SUBROUTINE    

*******************************************************************************/


double calc_root(double K,double x,double inflow,double outflow,double prevV,double dt)
{
  double V,curV,root;
  V = K * (x*inflow + (1-x)*outflow);
  curV = prevV+inflow*dt-outflow*dt;
  
  root = V - curV;
  
  return root;
}

void muskingum_solve_runtimestep(int * ierr,
                           double dt,
                           double inflow,
                           double K,
                           double x,
                           double * routing_volume,
                           double * routing_outflow
)
{
  
  // dV/dt = I - O
  // V = k [xI + (1-x) O]
  
  //printf("previous_inflow=%0.2f \n",*previous_inflow);
  //printf("routing_volume=%0.2f \n",*routing_volume);
  //printf("routing_outflow=%0.2f \n",*routing_outflow);
  
  double prevV;
  prevV = *routing_volume;
  
  double min_outflow=0,max_outflow=prevV/dt+inflow;
  double min_root,max_root;
  
  double m,b;
  double new_outflow,new_root;
  int min_pos=0;
  double tol = 1e-12;
  int counter=0;
  
  min_root = calc_root(K,x,inflow,min_outflow,prevV,dt);
  max_root = calc_root(K,x,inflow,max_outflow,prevV,dt);
  
  if(min_root>0 && max_root>0){
    printf("both routing root bounds are positive \n");
    *ierr=1;
    return;
  }
  if(min_root<0 && max_root<0){
    printf("both routing root bounds are negative \n");
    printf("max_root=%0.2f K=%0.2f x=%0.2f inflow=%0.2f max_outflow=%0.2f prevV=%0.2f dt=%0.2f\n",max_root,K,x,inflow,max_outflow,prevV,dt);
    printf("min_root=%0.2f K=%0.2f x=%0.2f inflow=%0.2f min_outflow=%0.2f prevV=%0.2f dt=%0.2f\n",min_root,K,x,inflow,min_outflow,prevV,dt);
    *ierr=1;
    return;
  }

  if(min_root>0){
    min_pos = 1;
  }
  
  new_outflow = min_outflow;
  new_root = min_root;
  //printf("min_root=%0.2f \n",min_root);
  //printf("max_root=%0.2f \n",max_root);
  
  while(fabs(new_root)>tol && fabs((max_outflow-min_outflow)/max_outflow)>tol){
    counter = counter+1;
    /*
    if(counter>2){
      printf("fabs(new_root)=%0.9f\n",fabs(new_root));
      printf("fabs(max_outflow-min_outflow)=%0.9f\n",fabs(max_outflow-min_outflow));
      printf("fabs((max_outflow-min_outflow)/max_outflow)=%0.9f\n",fabs((max_outflow-min_outflow)/max_outflow));
      printf("max_outflow=%0.9f\n",max_outflow);
      printf("min_outflow=%0.9f\n",min_outflow);
    }
    */
    /*
    if(counter>2){
      printf("prevV=%0.2f\n",prevV);
      printf("inflow=%0.2f\n",inflow*dt);
      printf("new_outflow=%0.2f\n",new_outflow*dt);
    }
    */
    if(counter>100){
      printf("solver is taking too long \n");
      break;
    }
    // get linear parameters
    m = (max_root-min_root)/(max_outflow-min_outflow);
    b = min_root - m*min_outflow;
    
    new_outflow=-b/m;
    if(new_outflow>=max_outflow || new_outflow<=min_outflow){
      //old_outflow = new_outflow;
      // bisect
      new_outflow = (max_outflow+min_outflow)/2;
      /*
      if(counter>2){
        printf("here\n");
        printf("new_outflow=%0.18f\n",new_outflow);
      }
      */
      
    }
    //printf("new_outflow=%0.2f \n",new_outflow);
    new_root = calc_root(K,x,inflow,new_outflow,prevV,dt);
    /*
    if(new_root>0){
      max_outflow = new_outflow;
      max_root = new_root;
    } else {
      min_outflow = new_outflow;
      min_root = new_root;
    }
    */
    
    if((min_pos==1 && new_root>0) || (min_pos==0 && new_root<0)){
      min_outflow = new_outflow;
      //min_root = calc_root(K,x,inflow,min_outflow,prevV,dt);
      min_root = new_root;
    } else {
      max_outflow = new_outflow;
      //max_root = calc_root(K,x,inflow,max_outflow,prevV,dt);
      max_root = new_root;
    }
    /*
    if(counter>2){
      printf("min_outflow=%0.18f\n",min_outflow);
      printf("min_root=%0.18f\n",min_root);
      printf("max_outflow=%0.18f\n",max_outflow);
      printf("max_root=%0.18f\n",max_root);
    }
    */
    
    //printf("new_root=%0.2f \n",new_root);
    
  }
  //printf("counter = %d \n",counter);
  /*
  if(counter>2){
    printf("fabs(new_root)=%0.9f\n",fabs(new_root));
    printf("fabs(max_outflow-min_outflow)=%0.9f\n",fabs(max_outflow-min_outflow));
    printf("fabs((max_outflow-min_outflow)/max_outflow)=%0.9f\n",fabs((max_outflow-min_outflow)/max_outflow));
    printf("max_outflow=%0.9f\n",max_outflow);
    printf("min_outflow=%0.9f\n",min_outflow);
    printf("counter = %d \n",counter);
  }
  */
  
  *routing_outflow = new_outflow;
  *routing_volume = K * (x*inflow + (1-x)*new_outflow);
  /*
  if(counter>100){
    printf("*routing_outflow=%0.2f\n",*routing_outflow);
    printf("*routing_volume=%0.2f\n",*routing_volume);
    printf("inflow=%0.2f\n",inflow);
    printf("new_outflow=%0.2f\n",new_outflow);
    printf("x=%0.2f\n",x);
    printf("K=%0.2f\n",K);
  }
  */

  
}


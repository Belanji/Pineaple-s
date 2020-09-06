#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "cg.h"
#include "pineaple.h"

int main(int argc, char *argv[]){


  double* cg;
  int l1max,l2max;
  int l3min,l3max;
  int l1,l2,l3;


  l1max=150;
  l2max=l1max;

  
  fill_cg(l1max,l2max,cg);
 
}

void fill_cg(int l1max,int l2max, double* cg)
{

  int l3min,l3max;
  int l1,l2,l3;
  int l3step=l1max+l2max;
      
  cg=calloc( l1max*l2max*(l1max+l2max) , sizeof(double) );
  
  //cg_file=fopen("cg.dat","w");

  
  for(l1=0; l1<=l1max; l1++)
    {

      for(l2=0; l2<=l2max; l2++)
	{

	  l3min=abs(l2-l1);
	  l3max=l2+l1;

	  for(l3=l3min; l3<=l3max;l3++)
	    {


	      cg[l2max*l3step*l1+l3step*l2+l3]=clebsch_gordan(l1,l2,l3,0,0);
	

	    
	    };
	};
    };
};

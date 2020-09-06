#include <gsl/gsl_sf.h>
#include <math.h>
#include "cg.h"



double clebsch_gordan( int l1, int l2, int l3,
		       int m1, int m2)
{
  //Calculate the clebsch-gordon coefficients:


  double wigner_3j;

  wigner_3j=gsl_sf_coupling_3j( 2*l1,2*l2,2*l3,2*m1,2*m2,-2*(m1+m2) );

  return pow(-1,m1+m2+l1-l2)*sqrt(2.*l3+1)*wigner_3j;


};



 

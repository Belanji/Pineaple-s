#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <petscsnes.h>
#include "cg.h"
#include "pineaple.h"


char help[]="Drink water!!!";

int main(int argc, char *argv[])
{



  PetscErrorCode ierr;
  Vec  Npm, Vp, Npm_rhs, Vp_rhs;
  Mat Vp_mat, Npm_jac;

   
  
  double *cg;
  double time=3.14;
  
  PetscInt l1max,l2max;
  PetscInt l3min,l3max;
  PetscInt l1,l2,l3;

  struct PNP_constants  constants;

  ierr = PetscInitialize(&argc,& argv,(char*)0,help);
  
  constants.q=1.;
  constants.d=1;
  constants.epsilon=1.;
  constants.V0=2.;
  constants.omega=1.;
  constants.Dc=1;
  constants.N0=1;
  
  constants.l1max=10;
  l1max=constants.l1max;
  l2max=l1max;

  
  
  //Setting up the Matrixes and vectors:

  //Setting up the charge density vectors:
  VecCreate( PETSC_COMM_SELF , &Npm  );
  VecSetSizes(Npm, PETSC_DECIDE, 2*(l1max+1)+4);
  VecSetType(Npm, VECSEQ  );
  VecDuplicate(Npm, &Npm_rhs);
  
  //Important Note: The vector is organized as follows:
  //Npm[0]=sigma_p (lower bound)
  //Npm[1:(l1max+1)]= positive charge density
  //Npm[(l1max+1)+1]=sigma_p (upper bound)

  //Npm[(l1max+1)+2]=sigma_m (lower bound)
  //Npm[(l1max+4):(2*l1max+4)]= negative charge density
  //Npm[l1max+5]=sigma_m (upper bound)


  
  //Setting up Eletric potential Vectors vectors:
  VecCreate( PETSC_COMM_SELF , &Vp  );
  VecSetSizes(Vp, PETSC_DECIDE, l1max+1 );
  VecSetType(Vp, VECSEQ  );
  VecDuplicate( Vp , &Vp_rhs  );




  MatCreateSeqDense(PETSC_COMM_SELF,l1max+1 , l1max+1, NULL, &Vp_mat );
  MatSetUp(Vp_mat);
  
  
  //Setting up rhs vectors:
  


  
  //Sertting up coefficients Matrixes:

  fill_cg(l1max,l2max,cg);

  //Setting up initial consitions:

  
  VecSet(Vp,0);

  
  
  fill_Vp_rhs(Vp_rhs, Vp,Npm, constants, time);
  fill_Vp_Matrix(Vp_mat, Vp, Npm, constants, time);


  KSP  solver;

  KSPCreate(PETSC_COMM_SELF,&solver);
  KSPSetOperators(solver,Vp_mat,Vp_mat);
  KSPSolve(solver, Vp_rhs ,Vp);


  homogeneous_ic(Npm, constants);

  fill_Npm_rhs(Npm_rhs, Vp, Npm,constants,time);
  
};


void homogeneous_ic(Vec Npm,
		    struct PNP_constants constants)
{

  PetscInt l1max=constants.l1max;
  
  PetscInt ip=1;
  PetscInt in=l1max+4;


  PetscScalar N0=constants.N0;
  
  VecSet(Npm,0);

  VecSetValues(Npm,1,&ip,&N0,ADD_VALUES);
  VecSetValues(Npm,1,&in,&N0,ADD_VALUES);

}

void fill_Vp_Matrix(Mat Vp_mat, 
		    const Vec Vp,
		    const Vec Npm,                 
		    struct PNP_constants constants,
		    const double  time)

{

  PetscInt ii,jj,kk;
  PetscInt l1max=constants.l1max;

  PetscScalar rhs;
  const PetscScalar rhs_even=1., rhs_odd=-1.;

  MatZeroEntries(Vp_mat);



  kk=l1max-1;
  for(ii=0; ii<=l1max;ii+=2) MatSetValues(Vp_mat,1,&kk,1,&ii,&rhs_even, ADD_VALUES);
  for(ii=1; ii<=l1max;ii+=2) MatSetValues(Vp_mat,1,&kk,1,&ii,&rhs_odd,  ADD_VALUES);

  
  
  kk=l1max;
  for(ii=0; ii<=l1max;ii++) MatSetValues(Vp_mat,1,&kk,1,&ii,&rhs_even,ADD_VALUES);

  

  for(kk=0; kk<l1max-1;kk++)
    {
      for(ii=0; ii<=l1max;ii++)
	{

	  for(jj=0; jj<=ii-2;jj++)
	    {

	      if ( jj==kk && (jj+ii)%2==0 )
		{

		  rhs=(kk+0.5)*( ii*(ii+1) - kk*(kk+1) )*gammaP(kk);
		  MatSetValues(Vp_mat,1,&kk,1,&ii,&rhs,ADD_VALUES);

		}
	      
	    }
	  
	}

    }

  MatAssemblyBegin(Vp_mat,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(Vp_mat,MAT_FINAL_ASSEMBLY);
  
};


PetscErrorCode fill_Npm_rhs(Vec Npm_rhs, 
                 const Vec Vp_inp,
                 const Vec Npm_inp,
                 const struct PNP_constants constants,
                 const double  time)
{

  PetscInt ii,ip, im;  
  const PetscScalar *Np, *Nm, *Vp;
  PetscScalar * rhs;
  PetscErrorCode ierr;


  const PetscInt l1max=constants.l1max;
  const double omega=constants.omega;
  const double V0=constants.V0;
  const double qd2_epsi=constants.q*constants.d*constants.d/constants.epsilon;

  const double Dc= constants.Dc;



  VecGetArrayRead(Npm_inp, & Np);
  VecGetArrayRead(Vp_inp, & Vp);
  VecGetArray(Npm_rhs, & rhs);


  Nm=&Np[l1max+3];

  printf("Np=%lf\n",Np[1]);
  printf("Nm=%lf\n",Nm[1]);

  VecRestoreArray(Npm_rhs, & rhs);
  VecRestoreArrayRead(Vp_inp, & Vp);
  VecRestoreArrayRead(Npm_inp, & Np);
  
  

}



PetscErrorCode fill_Vp_rhs(Vec Vp_rhs, 
                 const Vec Vp,
                 const Vec Npm,
                 const struct PNP_constants constants,
                 const double  time)

{

  PetscInt ii,ip, im;
  PetscInt l1max=constants.l1max;

  const PetscScalar *Npm_it;
  PetscScalar * k_i;

  PetscErrorCode ierr;
  
  double omega=constants.omega;
  double V0=constants.V0;
  double qd2_epsi=constants.q*constants.d*constants.d/constants.epsilon;

  VecGetArrayRead(Npm, & Npm_it);
  VecGetArray(Vp_rhs, & k_i);
  
  //Filing the rhs:
  for(ii=0; ii<l1max-1; ii++)
    {
      ip=ii+1;
      im=ii+(l1max+3);
	
  
      k_i[ii]=-qd2_epsi * ( Npm_it[ip] - Npm_it[im] );
	 
  
    }

  k_i[l1max-1]=0.5*V0*( cos(omega*time)+I*sin(omega*time) );
  k_i[l1max]= -0.5*V0*( cos(omega*time)+I*sin(omega*time) );


  VecRestoreArray(Vp_rhs, & k_i);
  VecRestoreArrayRead(Npm, & Npm_it);;

  return ierr;
}


 
void fill_cg(PetscInt l1max,PetscInt l2max, double* cg)
{

  PetscInt l3min,l3max;
  PetscInt l1,l2,l3;
  PetscInt l3step=l1max+l2max+1;
      
  cg=(double *)  calloc( (l1max+1)*(l2max+1)*( l1max + l2max + 1 ) , sizeof(double) );
  
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

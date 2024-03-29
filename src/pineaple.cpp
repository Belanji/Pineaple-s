#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <complex.h>
#include <petscsnes.h>
#include <gsl/gsl_sf_legendre.h>
#include "cg.h"
#include "pineaple.h"

const double pi=3.14;
double *cg;


char help[]="Drink water!!!";

int main(int argc, char *argv[])
{



  PetscErrorCode ierr;
  Vec  Npm, Vp, Npm_rhs, Vp_rhs;
  Mat Vp_mat, Npm_jac;


  PetscInt l1max;
  PetscInt l3min,l3max;

  l1max=70;
  //PetscInt ln=1;
  //PetscInt l1=;
  //PetscInt l1[ln]={4};
    
  
  struct PNP_constants  constants;

  ierr = PetscInitialize(&argc,& argv,(char*)0,help);

  constants.V0=0.025;
  constants.d=3e-6;
  constants.T=300;

  constants.omega=60*2*pi;
  constants.Dc=7.5e-11;
  constants.N0=4e16;

  constants.tau1=1e-2;
  constants.tau2=1e-2;

  constants.kappa1=1e-4;
  constants.kappa2=1e-4;

  
  constants.l1max=l1max;
  


  double Nturns=8;  
  double time=0.0;
  double tf=Nturns/constants.omega;
  double timePrint=tf;
  double dt=1/(constants.omega*20000);

  constants.dt=dt;
  constants.charge_adsorbed=fopen("eletric_surfaces.dat","w");
  fprintf(constants.charge_adsorbed,"#t    Surface_1(np   nm)  Surface_2(np   nm)\n");
  
  constants.effective_adsortion=fopen("eletric_time.dat","w");
  fprintf(constants.effective_adsortion,"#t    Sigma1(R+I)  Sigma2(R +i)  E1(R+I) E2(R+I)\n");
  
  //Setting up the Matrixes and vectors:

  //Setting up the charge density vectors:
  VecCreate( PETSC_COMM_SELF , &Npm  );
  VecSetSizes(Npm, PETSC_DECIDE, 2*l1max+6);
  VecSetType(Npm, VECSEQ  );
  VecDuplicate(Npm, &Npm_rhs);

  MatCreateSeqDense(PETSC_COMM_SELF,2*l1max+6 , 2*l1max+6, NULL, &Npm_jac );
  MatSetUp(Npm_jac);

  
  //Important Note: The vector is organized as follows:
  //Npm[0]=sigma_p (lower bound)
  //Npm[1:(l1max+1)]= positive charge density
  //Npm[(l1max+1)+1]=sigma_p (upper bound)

  //Npm[(l1max+1)+2]=sigma_m (lower bound)
  //Npm[(l1max+4):(2*l1max+4)]= negative charge density
  //Npm[2*l1max+5]=sigma_m (upper bound)


  
  //Setting up Eletric potential Vectors vectors:
  VecCreate( PETSC_COMM_SELF , &Vp  );
  VecSetSizes(Vp, PETSC_DECIDE, l1max+1 );
  VecSetType(Vp, VECSEQ  );
  VecDuplicate( Vp , &Vp_rhs  );




  MatCreateSeqDense(PETSC_COMM_SELF,l1max+1 , l1max+1, NULL, &Vp_mat );
  MatSetUp(Vp_mat);
  
  
  //Setting up rhs vectors:
  


  
  //Sertting up coefficients Matrixes:

  //fill_cg(l1max,l1max,cg);

  //Setting up initial consitions:
  VecSet(Vp,0);
  homogeneous_ic(Npm, constants);

  
  print_surface_charge(Npm,Vp, constants, time);
  Print_Charge_Density(Npm,Vp,constants,time,0);

  evolve(Npm, Vp, Npm_rhs, Vp_rhs, Vp_mat, Npm_jac, constants,time,time+timePrint);


  Print_Charge_Density(Npm,Vp,constants,time,1);

  fclose(constants.charge_adsorbed);
  fclose(constants.effective_adsortion);
  return ierr;
};

//double legendreP(double x,
//		 int ii)
//
//{
//
//
//}
PetscErrorCode Print_Charge_Density(const Vec Npm,
				    const Vec Vp_inp,
				    const struct PNP_constants constants,
				    double time,
				    int figNumber)
{

  int ii;
  int jj,jp, jm;
  PetscErrorCode ierr;
  PetscInt l1max=constants.l1max;
  double dx=2./l1max;
  double xx;
  double fp_value,fm_value, fvr_value,fvi_value;
  const PetscScalar *Np, *Vp;
  FILE * snapshot;
  char fig_name[100];

  sprintf(fig_name,"charge_v_%d.dat",figNumber);
  
  snapshot=fopen( fig_name, "w" );
  fprintf(snapshot,"   x     np   nm    real[V]  imag[V]\n");

  
  VecGetArrayRead(Npm, & Np);
  VecGetArrayRead(Vp_inp, & Vp);
  
  for(ii=0; ii<=l1max; ii++)
    {
      fp_value=0;
      fm_value=0;
      fvr_value=0;
      fvi_value=0;
      
      xx=-1.+dx*ii;
      
      for(jj=0; jj<=l1max; jj++)
	{
	  jp=jj+1;
	  jm=jj+l1max+4;
	  
	  fp_value+=real(Np[jp])*gsl_sf_legendre_Pl(jj,xx);
	  fm_value+=real(Np[jm])*gsl_sf_legendre_Pl(jj,xx);
	  fvr_value+=real(Vp[jj])*gsl_sf_legendre_Pl(jj,xx);
	  fvi_value+=imag(Vp[jj])*gsl_sf_legendre_Pl(jj,xx);
	  
	}

      fprintf(snapshot,"%.3lf %.3lf %.3lf %.3lf %.3lf\n",xx,fp_value,fm_value,fvr_value, fvi_value);
      
    }

  VecRestoreArrayRead(Npm, & Np);
  fclose(snapshot);
  
  return ierr;
}


PetscErrorCode homogeneous_ic(Vec Npm,
		    struct PNP_constants constants)
{

  PetscErrorCode ierr;
  PetscInt l1max=constants.l1max;
  
  PetscInt ip=1;
  PetscInt im=l1max+4;


  PetscScalar N0=constants.N0/constants.d;
  
  VecSet(Npm,0);

  VecSetValues(Npm,1,&ip,&N0,ADD_VALUES);
  VecSetValues(Npm,1,&im,&N0,ADD_VALUES);

  return ierr;
}

PetscErrorCode fill_Vp_Matrix(Mat Vp_mat, 
		    const Vec Vp,
		    const Vec Npm,                 
		    struct PNP_constants constants,
		    const double  time)
{
  
  PetscErrorCode ierr;
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

  return ierr;
  
};


PetscErrorCode fill_Npm_Matrix(Mat Npm_mat, 
		    const Vec Vp_inp,
		    const Vec Npm_inp,                 
		    struct PNP_constants constants,
		    const double  time)
{

  PetscInt ii,ip, im;
  PetscInt jj,jp, jm;
  PetscInt kk,kp, km;
  PetscInt pp,qq;
  const PetscScalar *Np, *Nm, *Vp;
  PetscErrorCode ierr;
  PetscScalar rhs;
  //double *cg=

  const PetscInt l1max=constants.l1max;  
  const double Dc= constants.Dc/constants.d;
  const double Dcd2= Dc/constants.d;
  const double qD_kBtd=constants.q*Dcd2/(constants.k_b*constants.T);

  

  const double tau1=constants.tau1;
  const double tau2=constants.tau2;

  const double kappa1=constants.kappa1;
  const double kappa2=constants.kappa2;

  
  VecGetArrayRead(Npm_inp, & Np);
  VecGetArrayRead(Vp_inp, & Vp);

  
  MatZeroEntries(Npm_mat);



  for (ii=0; ii<=l1max;ii++)
    {

      ip=ii+1;
      im=ii+(l1max+4);

      for(kk=0;kk<l1max-1;kk++)
	{
	  
	  kp=kk+1;
	  km=kk+(l1max+4);

	  //Difussion term:
	  for(jj=0; jj<=ii-2; jj++)
	    {
	  
	      if ( jj==kk && (ii+jj)%2==0)
		{
		  
		  rhs=Dcd2*(kk+0.5)*( ii*(ii+1) - kk*(kk+1) )*gammaP(kk);
		  MatSetValues(Npm_mat,1,&kp,1,&ip,&rhs,ADD_VALUES);
		  MatSetValues(Npm_mat,1,&km,1,&im,&rhs,ADD_VALUES);
		}

	    }

	  //First Derivative terms:
	  for(jj=0; jj<=ii-1; jj++)
	    {

	      for(pp=0; pp<=l1max; pp++)
		{
		  for(qq=0; qq<=pp-1; qq++)
		    {

		      if( (ii+jj)%2==1 && (pp+qq)%2==1 )
			{
			  
			  rhs=qD_kBtd*2.*(2*jj+1)*Vp[pp]*clebsch_gordan(jj,kk,qq,0,0)*clebsch_gordan(jj,kk,qq,0,0);

			  MatSetValues(Npm_mat,1,&kp,1,&ip,&rhs,ADD_VALUES);

			  rhs=-rhs;
			  MatSetValues(Npm_mat,1,&km,1,&im,&rhs,ADD_VALUES);

			  
			}

		      
		    }
		    
		}
	      
	    }

	  //Second Derivative terms:
	  for(pp=0; pp<=l1max; pp++)
	    {
	      for(qq=0; qq<=pp-2; qq++)
		{

		  if( (pp+qq)%2==0 )
		    {
			  
		      rhs=qD_kBtd*( pp*(pp+1)-qq*(qq+1) )*Vp[pp]*clebsch_gordan(ii,kk,qq,0,0)*clebsch_gordan(ii,kk,qq,0,0);

		      MatSetValues(Npm_mat,1,&kp,1,&ip,&rhs,ADD_VALUES);

		      rhs=-rhs;
		      MatSetValues(Npm_mat,1,&km,1,&im,&rhs,ADD_VALUES);

			  
		    }

		      
		}
		    
	    }



	  
	}

    }

  //Setting Evolution equations for sigma:

  //Lower boundary sigmas:
  ip=0;
  kp=0;

  km=l1max+3;
  im=l1max+3;

  rhs=-1./tau1;

  
  MatSetValues(Npm_mat,1,&kp,1,&ip,&rhs,ADD_VALUES);
  MatSetValues(Npm_mat,1,&km,1,&im,&rhs,ADD_VALUES);

  for(ii=0; ii<=l1max;ii++)
    {
      ip=ii+1;
      im=ii+(l1max+4);
      
      
      rhs=kappa1*pow(-1,ii);
      MatSetValues(Npm_mat,1,&kp,1,&ip,&rhs,ADD_VALUES);    
      MatSetValues(Npm_mat,1,&km,1,&im,&rhs,ADD_VALUES);
      
    }


  //Filling the upper boundary matrix:
  
  kp=l1max+2;
  ip=l1max+2;

  km=2*l1max+5;
  im=2*l1max+5;

  
  rhs=-1/tau2;
  MatSetValues(Npm_mat,1,&kp,1,&ip,&rhs,ADD_VALUES);
  MatSetValues(Npm_mat,1,&km,1,&im,&rhs,ADD_VALUES);

  
  for(ii=0; ii<=l1max;ii++)
    {
      ip=ii+1;
      im=ii+(l1max+4);

      
      rhs=kappa2;
      MatSetValues(Npm_mat,1,&kp,1,&ip,&rhs,ADD_VALUES);
      MatSetValues(Npm_mat,1,&km,1,&im,&rhs,ADD_VALUES);
      
    }

  
    

  VecRestoreArrayRead(Npm_inp, & Np);
  VecRestoreArrayRead(Vp_inp, & Vp);

  
  MatAssemblyBegin(Npm_mat,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(Npm_mat,MAT_FINAL_ASSEMBLY);

  return ierr;
};



PetscErrorCode fill_Npm_rhs(Vec Npm_rhs, 
                 const Vec Vp_inp,
                 const Vec Npm_inp,
                 const struct PNP_constants constants,
                 const double  time)
{

  PetscInt ii,ip, im;
  PetscInt jj,jp, jm;
  PetscInt kk,kp, km;  
  const PetscScalar *Npm, *Vp;
  PetscScalar * rhs;
  PetscErrorCode ierr;


  const PetscInt l1max=constants.l1max;
  const double omega=constants.omega;
  const double V0=constants.V0;
  const double qd2_epsi=constants.q*constants.d*constants.d/constants.epsilon;

  const double Dcd2= constants.Dc*constants.d*constants.d;


  const double tau1=constants.tau1;
  const double tau2=constants.tau2;

  const double kappa1=constants.kappa1;
  const double kappa2=constants.kappa2;

  
  
  //Set all values to 0.
  VecSet(Npm_rhs,0);

  VecGetArrayRead(Npm_inp, & Npm);
  VecGetArrayRead(Vp_inp, & Vp);
  VecGetArray(Npm_rhs, & rhs);

  //Sigmas_P rhs:
  rhs[0]=Npm[0];
  rhs[l1max+2]=Npm[l1max+2];

  //Sigma_M rhs:
  rhs[l1max+3]=Npm[l1max+3];
  rhs[2*l1max+5]=Npm[2*l1max+5];

  //boundary rhs:
  rhs[l1max]=0;
  rhs[l1max+1]=0;

  rhs[2*l1max+3]=0;
  rhs[2*l1max+4]=0;

  for (ii=0; ii<l1max-1;ii++)
    {
      ip=ii+1;
      im=ii+l1max+4;
      
      rhs[ip]=gammaP(ii)*Npm[ip];
      rhs[im]=gammaP(ii)*Npm[im];

    }
  
  //for (ii=0; ii<=l1max;ii++)
  //  {
  //
  //    ip=ii+1;
  //    in=ii+1;
  //	
  //    
  //    for(jj=0; jj<=ii-2; jj++)
  //	{
  //
  //	  jp=jj+1;
  //	  jn=jj+1;
  //
  //	  for(kk=0;kk<=l1max-2;kk++)
  //	    {
  //
  //	      kp=kk+1;
  //	      kn=kk+1+(l1max+3);
  //
  //
  //	      if ( jj==kk && (ii+jj)%2==0)
  //		{
  //
  //		  rhs[kp]+=(kk+0.5)*( ii*(ii+1) - kk*(kk+1) )*gammaP(kk);
  //		  rhs[kn]+=(kk+0.5)*( ii*(ii+1) - kk*(kk+1) )*gammaP(kk);
  //		}
  //	      
  //	    }
  //
  //	}
  //
  //  }
  //
  ////Filling the sigmas Rhs:
  //
  //rhs[0]=Np[0]/tau1;
  //rhs[l1max+2]=-Np[l1max+2]/tau2;
  //    
  //rhs[l1max+3]=Nm[l1max+3]/tau1;
  //rhs[2*l1max+5]=-Nm[2*l1max+5]/tau2;
  //
  //for(ii=0; ii<=l1max;ii++)
  //  {
  //    ip=ii+1;
  //    in=ii+1;
  //
  //    //sigma_P:
  //    rhs[0]-=kappa1*pow(-1,ii)*Np[ii];
  //    rhs[l1max+2]+=kappa2*Np[ii];
  //    
  //    rhs[l1max+3]-=kappa1*pow(-1,ii)*Nm[ii];
  //    rhs[2*l1max+5]+=kappa2*Np[ii];
  //    
  //  }
  
  VecRestoreArray(Npm_rhs, & rhs);
  VecRestoreArrayRead(Vp_inp, & Vp);
  VecRestoreArrayRead(Npm_inp, & Npm);
  
  
  return ierr;
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
      im=ii+(l1max+4);
	
  
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
	

	    
	    }
	}
    }
}


PetscErrorCode Enforce_Boundary_Conditions(Mat Npm_mat, 
		    const Vec Vp_inp,
		    const Vec Npm_inp,                 
		    struct PNP_constants constants,
		    const double  time)
{

  PetscInt ii,ip, im;
  PetscInt jj,jp, jm;
  PetscInt kk,kp, km;  
  const PetscScalar *Np, *Nm, *Vp;
  PetscErrorCode ierr;
  PetscScalar rhs;
  

  const PetscInt l1max=constants.l1max;

  const double Dc= constants.Dc/constants.d;
  const double qD_kBt=constants.q*Dc/(constants.k_b*constants.T);



  
  const double tau1=constants.tau1;
  const double tau2=constants.tau2;

  const double kappa1=constants.kappa1;
  const double kappa2=constants.kappa2;

  PetscScalar E_lower=0.;
  PetscScalar E_upper=0.;

  
  VecGetArrayRead(Npm_inp, & Np);
  VecGetArrayRead(Vp_inp, & Vp);

  for(ii=0; ii<=l1max;ii++) E_lower+=0.5*pow(-1,ii-1)*ii*(ii+1)*Vp[ii];
  for(ii=0; ii<=l1max;ii++) E_upper+=0.5*ii*(ii+1)*Vp[ii];

  

  kp=l1max;
  km=2*l1max+3;


  ip=0;
  im=l1max+3;
  
  
  rhs=-1./tau1;
  MatSetValues(Npm_mat,1,&kp,1,&ip,&rhs,ADD_VALUES);
  MatSetValues(Npm_mat,1,&km,1,&im,&rhs,ADD_VALUES);
  
  for(ii=0; ii<=l1max;ii++)
    {
      ip=ii+1;
      im=ii+(l1max+4);

      
      rhs=-0.5*Dc*pow(-1,ii-1)*ii*(ii+1)+pow(-1,ii)*(kappa1-qD_kBt*E_lower);
      MatSetValues(Npm_mat,1,&kp,1,&ip,&rhs,ADD_VALUES);


      rhs=-0.5*Dc*pow(-1,ii-1)*ii*(ii+1)+pow(-1,ii)*(kappa1+qD_kBt*E_lower);
      MatSetValues(Npm_mat,1,&km,1,&im,&rhs,ADD_VALUES);
      
    }

  
  kp=l1max+1;
  km=2*l1max+4;

  
  ip=l1max+2;
  im=2*l1max+5;
  
  rhs=1./tau2;
  MatSetValues(Npm_mat,1,&kp,1,&ip,&rhs,ADD_VALUES);
  MatSetValues(Npm_mat,1,&km,1,&im,&rhs,ADD_VALUES);

  
  for(ii=0; ii<=l1max;ii++)
    {
      ip=ii+1;
      im=ii+(l1max+4);

      //sigma_PM Upper:
      rhs=-0.5*Dc*ii*(ii+1)-qD_kBt*E_upper-kappa2;
      MatSetValues(Npm_mat,1,&kp,1,&ip,&rhs,ADD_VALUES);

      rhs=-0.5*Dc*ii*(ii+1)+qD_kBt*E_upper-kappa2;
      MatSetValues(Npm_mat,1,&km,1,&im,&rhs,ADD_VALUES);
           
    }


  VecRestoreArrayRead(Npm_inp, & Np);
  VecRestoreArrayRead(Vp_inp, & Vp);
 
  MatAssemblyBegin(Npm_mat,MAT_FINAL_ASSEMBLY);
  MatAssemblyEnd(Npm_mat,MAT_FINAL_ASSEMBLY);

  return ierr;
}


PetscErrorCode evolve(Vec Npm,
		      Vec Vp,
		      Vec Npm_rhs,
		      Vec Vp_rhs,
		      Mat Vp_mat,
		      Mat Npm_jac,
		      const struct PNP_constants constants,
		      double  time,
		      double tf)
{

  PetscInt l1max=constants.l1max;
  PetscInt ii,ip, im;
  PetscInt jj,jp, jm;
  PetscInt kk,kp, km;  
  PetscErrorCode ierr;
  PetscScalar rhs;
  double dt=constants.dt;
  KSP  solver, vp_solver;


  //Setting up the first eletric potential:
  
  fill_Vp_Matrix(Vp_mat, Vp, Npm, constants, time);
  fill_Vp_rhs(Vp_rhs, Vp,Npm, constants, time);

  
  KSPCreate(PETSC_COMM_SELF,&vp_solver);
  KSPSetOperators(vp_solver,Vp_mat,Vp_mat);
  KSPSolve(vp_solver, Vp_rhs ,Vp);

    
  KSPCreate(PETSC_COMM_SELF,&solver);

  while(time<tf)
    {
      

      fill_Npm_rhs( Npm_rhs, Vp, Npm, constants, time);
      fill_Npm_Matrix(Npm_jac, Vp, Npm, constants, time);
      MatScale(Npm_jac,-dt);
  

      //Filling main DIagonal:
      
      
      kp=0;    
      km=l1max+3;

      rhs=1;
      MatSetValues(Npm_jac,1,&kp,1,&kp,&rhs,ADD_VALUES);
      MatSetValues(Npm_jac,1,&km,1,&km,&rhs,ADD_VALUES);
  
      for(ii=0; ii<l1max-1;ii++)
	{
	  kp=ii+1;
	  km=ii+l1max+4;
	  
	  rhs=gammaP(ii);
	  
	  MatSetValues(Npm_jac,1,&kp,1,&kp,&rhs,ADD_VALUES);
	  MatSetValues(Npm_jac,1,&km,1,&km,&rhs,ADD_VALUES);

	}


      rhs=1;
      kp=l1max+2;
      km=2*l1max+5;
      
      MatSetValues(Npm_jac,1,&kp,1,&kp,&rhs,ADD_VALUES);
      MatSetValues(Npm_jac,1,&km,1,&km,&rhs,ADD_VALUES);

      


      Enforce_Boundary_Conditions(Npm_jac, Vp, Npm, constants, time);
      
      //MatView(Npm_jac,	PETSC_VIEWER_STDOUT_SELF);
      //exit(0);
    
      KSPSetOperators(solver,Npm_jac,Npm_jac);
      KSPSolve(solver, Npm_rhs ,Npm);


      fill_Vp_rhs(Vp_rhs, Vp,Npm, constants, time);
      KSPSolve(vp_solver, Vp_rhs ,Vp);
      
      time+=dt;

      print_surface_charge(Npm, Vp, constants, time);
    }

  
  return ierr;
}

void print_surface_charge(const Vec Npm_inp,
			  const Vec Vp_inp,
			  const struct PNP_constants constants,
			  const double time)
{
  const PetscScalar *Npm, *Vp;
  PetscInt l1max=constants.l1max;
  FILE * charge_adsorbed=constants.charge_adsorbed;
  FILE * effective_adsortion=constants.effective_adsortion;

  VecGetArrayRead(Npm_inp, & Npm);
  VecGetArrayRead(Vp_inp, & Vp);
  
  double bulk_charge_np=2*constants.d*real(Npm[1]);
  double bulk_charge_nm=2*constants.d*real(Npm[l1max+4]);
  
  double total_np=real(Npm[0])+real(Npm[l1max+2])+bulk_charge_np;
  double total_nm=real(Npm[l1max+3])+real(Npm[2*l1max+5])+bulk_charge_nm;
  

  PetscScalar E_lower=0.;
  PetscScalar E_upper=0.;


  for(PetscInt ii=0; ii<=l1max;ii++) E_lower+=0.5*pow(-1,ii-1)*ii*(ii+1)*Vp[ii];
  for(PetscInt ii=0; ii<=l1max;ii++) E_upper+=0.5*ii*(ii+1)*Vp[ii];


  fprintf(charge_adsorbed,"%e %e %e %e %e %e %e %e %e\n",
	  time,
	  real(Npm[0]),
	  imag(Npm[0]),
	  real(Npm[l1max+3]),
	  imag(Npm[l1max+3]),
	  real(Npm[l1max+2]),
	  imag(Npm[l1max+2]),
	  real(Npm[2*l1max+5]),
	  imag(Npm[2*l1max+5])
	  );

  fprintf(effective_adsortion,"%e %e %e %e %e %e %e %e %e\n",
	  time,
	  real(Npm[0]-Npm[l1max+3]),	  
	  imag(Npm[0]-Npm[l1max+3]),	  
	  real(Npm[l1max+2]-Npm[2*l1max+5]),
	  imag(Npm[l1max+2]-Npm[2*l1max+5]),
	  real(-E_lower),
	  imag(-E_lower),
	  real(-E_upper),
	  imag(-E_upper)
	  );


	  
  printf("%e %e %e\n",time,total_np,total_nm);

  VecRestoreArrayRead(Npm_inp, & Npm);
  VecRestoreArrayRead(Vp_inp, & Vp);
}

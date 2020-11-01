

//Usefull functions:

int min (int a ,int b) {

  return a < b ? a : b ;

}

//Delta de Kronecker:
double deltaK (int a ,int b) {

  return a == b ? 1 : 0 ;

}


//Legendre polinomial integeral:
double gammaP(int ii)
{

  
  return 2./(2.*ii+1);

}


struct PNP_constants{

  double q;
  double d;
  double epsilon;
  double V0;
  double omega;
  double Dc;
  double N0;
  double kappa1;
  double kappa2;
  double tau1;
  double tau2;
  double dt;
  
  int l1max;
  
  
};

void fill_cg(int,int,double*);

PetscErrorCode fill_Vp_rhs(Vec Vp_rhs, 
			   Vec Vp,
			   Vec Npm,
			   struct PNP_constants constants,
			   double  time);
  
PetscErrorCode fill_Vp_Matrix(Mat Vp_matrix, 
			      const Vec Vp,
			      const Vec Npm,                 
			      const struct PNP_constants constants,
			      const double  time);

PetscErrorCode fill_Npm_rhs(Vec Npm_rhs, 
			    const Vec Vp_inp,
			    const Vec Npm_inp,
			    const struct PNP_constants constants,
			    const double  time);

PetscErrorCode homogeneous_ic(Vec Npm,
			      struct PNP_constants constants);

PetscErrorCode Enforce_Boundary_Conditions(Mat Npm_mat, 
					   const Vec Vp_inp,
					   const Vec Npm_inp,
					   struct PNP_constants constants,
					   const double  time);

PetscErrorCode fill_Npm_Matrix(Mat Npm_mat, 
			       const Vec Vp_inp,
			       const Vec Npm_inp,                 
			       struct PNP_constants constants,
			       const double  time);
  
PetscErrorCode evolve(Vec Npm,
		      Vec Vp,
		      Vec Npm_rhs,
		      Vec Vp_rhs,
		      Mat Vp_mat,
		      Mat Npm_jac,
		      const struct PNP_constants constants,
		      double  time);

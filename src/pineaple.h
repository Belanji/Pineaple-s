

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

  double d;
  const double q=1.60217662e-19;
  const double epsilon=6.7*8.8541878128e-12;
  const double k_b=1.38064852e-23;
  double T;
  double V0;
  double omega;
  double Dc;
  double N0;
  double kappa1;
  double kappa2;
  double tau1;
  double tau2;
  double dt;
  FILE *charge_adsorbed;
  FILE *effective_adsortion;
  
  int l1max;
  double *cg;
  
};

void fill_cg(int,double*);

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
		      double  time,
		      double tf);


PetscErrorCode Print_Charge_Density(Vec Npm,
				    Vec Vp,
				    const struct PNP_constants constants,
				    double time,
				    int figNumber);

double legendreP(double x,
		 int ii);

void print_surface_charge(const Vec Npm,
			  const Vec Vp,
			  const struct PNP_constants constants,
			  const double time);

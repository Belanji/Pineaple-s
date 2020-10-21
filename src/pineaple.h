

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

  int l1max;
  
  
};

void fill_cg(int,int,double*);

PetscErrorCode  fill_Vp_rhs(Vec Vp_rhs, 
                 Vec Vp,
                 Vec Npm,
                 struct PNP_constants constants,
                 double  time);
  
void fill_Vp_Matrix(Mat Vp_matrix, 
		    const Vec Vp,
		    const Vec Npm,                 
		    const struct PNP_constants constants,
		    const double  time);

PetscErrorCode fill_Npm_rhs(Vec Npm_rhs, 
			    const Vec Vp_inp,
			    const Vec Npm_inp,
			    const struct PNP_constants constants,
			    const double  time);

void homogeneous_ic(Vec Npm,
		    struct PNP_constants constants);

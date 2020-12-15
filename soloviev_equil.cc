#include "constants.h"

double elong, triang, del_zero, R_0, min_rad, inv_asp, solov_A, solov_C, flux_0, solov_alpha;

double theta=0.0;
bool param_is_X=true;
bool stupid_bool=false;

double c0;
double c1;
double c2;
double c3;
double c4;
double c5;
double c6;

double solov_flux(double param)
{
  double X,Y;
  if(param_is_X==true){X=param; Y=( elong*sin(theta) / cos(theta+del_zero*sin(theta)) ) * ( param-1.0 ); }
  else{Y=param; X=1.0+( param*cos(theta+del_zero*sin(theta)) / (elong*sin(theta)) ); }

  if(stupid_bool==true){std::cout << "X : " << X << " Y : " << Y << std::endl;}

  double U_P=0.5*solov_alpha*pow(X,2)*log(X)+0.125*(1.0-solov_alpha)*pow(X,4);
  double U_0=1.0;
  double U_1=pow(X,2);
  double U_2=pow(Y,2)-pow(X,2)*log(X);
  double U_3=pow(X,2)*(pow(X,2)-4.0*pow(Y,2));
  double U_4=2.0*pow(Y,4)-9.0*pow(Y,2)*pow(X,2)-(12.0*pow(Y,2)*pow(X,2)-3.0*pow(X,4))*log(X);
  double U_5=pow(X,6)-12.0*pow(X,4)*pow(Y,2)+8.0*pow(X,2)*pow(Y,4);
  double U_6=8.0*pow(Y,6)-140.0*pow(Y,4)*pow(X,2)+75.0*pow(Y,2)*pow(X,4)-(120.0*pow(Y,4)*pow(X,2)-180.0*pow(Y,2)*pow(X,4)+15.0*pow(X,6))*log(X);

  return U_P+c0*U_0+c1*U_1+c2*U_2+c3*U_3+c4*U_4+c5*U_5+c6*U_6;
}

double U_P(double X,double Y){ return 0.5*solov_alpha*pow(X,2)*log(X)+0.125*(1.0-solov_alpha)*pow(X,4);}
double U_0(double X,double Y){ return 1.0;}
double U_1(double X,double Y){ return pow(X,2);}
double U_2(double X,double Y){ return pow(Y,2)-pow(X,2)*log(X);}
double U_3(double X,double Y){ return pow(X,2)*(pow(X,2)-4.0*pow(Y,2));}
double U_4(double X,double Y){ return 2.0*pow(Y,4)-9.0*pow(Y,2)*pow(X,2)-(12.0*pow(Y,2)*pow(X,2)-3.0*pow(X,4))*log(X);}
double U_5(double X,double Y){ return pow(X,6)-12.0*pow(X,4)*pow(Y,2)+8.0*pow(X,2)*pow(Y,4);}
double U_6(double X,double Y){ return 8.0*pow(Y,6)-140.0*pow(Y,4)*pow(X,2)+75.0*pow(Y,2)*pow(X,4)-(120.0*pow(Y,4)*pow(X,2)-180.0*pow(Y,2)*pow(X,4)+15.0*pow(X,6))*log(X);}

double U_P_X(double X,double Y){ return 0.5*solov_alpha*X*(1.0+2.0*log(X))+0.5*(1.0-solov_alpha)*pow(X,3);}
double U_0_X(double X,double Y){ return 0.0;}
double U_1_X(double X,double Y){ return 2.0*X;}
double U_2_X(double X,double Y){ return -X*(1.0+2.0*log(X));}
double U_3_X(double X,double Y){ return 4.0*pow(X,3)-8.0*X*pow(Y,2);}
double U_4_X(double X,double Y){ return -30.0*X*pow(Y,2)+3.0*pow(X,3)-12.0*X*(2.0*pow(Y,2)-pow(X,2))*log(X);}
double U_5_X(double X,double Y){ return 6.0*pow(X,5)-48.0*pow(X,3)*pow(Y,2)+16.0*X*pow(Y,4);}
double U_6_X(double X,double Y){ return -400.0*pow(Y,4)*X+480.0*pow(Y,2)*pow(X,3)-15.0*pow(X,5)-(240.0*X*pow(Y,4)-720.0*pow(Y,2)*pow(X,3)+90.0*pow(X,5))*log(X);}

double U_P_Y(double X,double Y){ return 0.0;}
double U_0_Y(double X,double Y){ return 0.0;}
double U_1_Y(double X,double Y){ return 0.0;}
double U_2_Y(double X,double Y){ return 2.0*Y;}
double U_3_Y(double X,double Y){ return -8.0*pow(X,2)*Y;}
double U_4_Y(double X,double Y){ return 8.0*pow(Y,3)-18.0*Y*pow(X,2)-24.0*Y*pow(X,2)*log(X);}
double U_5_Y(double X,double Y){ return -24.0*pow(X,4)*Y+32.0*pow(X,2)*pow(Y,3);}
double U_6_Y(double X,double Y){ return 48.0*pow(Y,5)-560.0*pow(Y,3)*pow(X,2)+150.0*Y*pow(X,4)-(480.0*pow(Y,3)*pow(X,2)-360.0*Y*pow(X,4))*log(X);}

double U_P_XX(double X,double Y){ return 0.5*solov_alpha*(3.0+2.0*log(X))+1.5*(1.0-solov_alpha)*pow(X,2);}
double U_0_XX(double X,double Y){ return 0.0;}
double U_1_XX(double X,double Y){ return 2.0;}
double U_2_XX(double X,double Y){ return -3.0-2.0*log(X);}
double U_3_XX(double X,double Y){ return 12.0*pow(X,2)-8.0*pow(Y,2);}
double U_4_XX(double X,double Y){ return -54.0*pow(Y,2)+21.0*pow(X,2)+(36.0*pow(X,2)-24.0*pow(Y,2))*log(X);}
double U_5_XX(double X,double Y){ return 30.0*pow(X,4)-144.0*pow(X,2)*pow(Y,2)+16.0*pow(Y,4);}
double U_6_XX(double X,double Y){ return -640.0*pow(Y,4)+2160.0*pow(X,2)*pow(Y,2)-165.0*pow(X,4)-(240.0*pow(Y,4)-2160.0*pow(X,2)*pow(Y,2)+450.0*pow(X,4))*log(X);}

double U_P_YY(double X,double Y){ return 0.0;}
double U_0_YY(double X,double Y){ return 0.0;}
double U_1_YY(double X,double Y){ return 0.0;}
double U_2_YY(double X,double Y){ return 2.0;}
double U_3_YY(double X,double Y){ return -8.0*pow(X,2);}
double U_4_YY(double X,double Y){ return 24.0*pow(Y,2)-18.0*pow(X,2)-24.0*pow(X,2)*log(X);}
double U_5_YY(double X,double Y){ return -24.0*pow(X,4)+96.0*pow(X,2)*pow(Y,2);}
double U_6_YY(double X,double Y){ return 240.0*pow(Y,4)-1680.0*pow(X,2)*pow(Y,2)+150.0*pow(X,4)-(1440.0*pow(X,2)*pow(Y,2)-360.0*pow(X,4))*log(X);}

int fill_sol(equil_fields *equil,geom_shape geom)
{
  int N_psi=geom.N_psi;
  int N_interp=geom.N_interp;
  int N_theta=geom.N_theta;

  elong=geom.elong;
  triang=geom.triang;
  del_zero=asin(triang);

  R_0=geom.R_0;
  min_rad=geom.min_rad;
  inv_asp=min_rad/R_0;

  solov_A=geom.solov_A;
  solov_C=geom.solov_C;
  flux_0= R_0 * R_0 * ( solov_A + solov_C * R_0 * R_0 ) ;
  solov_alpha=R_0*R_0*solov_A/flux_0;

  //Curvature coefficients
  double N_1=-(1.0+del_zero)*(1.0+del_zero)/(inv_asp*elong*elong);
  double N_2=(1.0-del_zero)*(1.0-del_zero)/(inv_asp*elong*elong);
  double N_3=-elong/(inv_asp*(1.0-triang*triang));

  std::cout << "elong : " << elong << " triang : " << triang << " del_zero : " << del_zero << " R_0 : " << R_0 << " min_rad : " << min_rad << " inv_asp : " << inv_asp << " solov_A : " << solov_A << " solov_C : " << solov_C << " flux_0 : " << flux_0 << " solov_alpha : " << solov_alpha << " N_1 : " << N_1 << " N_2 : " << N_2 << " N_3 : " << N_3 << std::endl;
  
  //Solve for the coefficients of flux function
  Mat A;
  Vec coeffs,sols;
  KSP ksp;
  PetscErrorCode ierr;

  ierr = MatCreate(PETSC_COMM_SELF,&A); CHKERRQ(ierr);
  ierr = MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,7,7);CHKERRQ(ierr);
  ierr = MatSetFromOptions(A);CHKERRQ(ierr);
  ierr = MatSetUp(A);CHKERRQ(ierr);

  ierr = VecCreate(PETSC_COMM_SELF,&sols);CHKERRQ(ierr);
  ierr = VecSetSizes(sols,PETSC_DECIDE,7);
  ierr = VecSetFromOptions(sols);
  ierr = VecCreate(PETSC_COMM_SELF,&coeffs);CHKERRQ(ierr);
  ierr = VecSetSizes(coeffs,PETSC_DECIDE,7);
  ierr = VecSetFromOptions(coeffs);

  double X,Y;
  //Boundary condition 1 (zero outer midplane flux)
  X=1.0+inv_asp; Y=0.0;
  VecSetValue(sols,0,-U_P(X,Y),INSERT_VALUES);
  MatSetValue(A,0,0,U_0(X,Y),INSERT_VALUES);
  MatSetValue(A,0,1,U_1(X,Y),INSERT_VALUES);
  MatSetValue(A,0,2,U_2(X,Y),INSERT_VALUES);
  MatSetValue(A,0,3,U_3(X,Y),INSERT_VALUES);
  MatSetValue(A,0,4,U_4(X,Y),INSERT_VALUES);
  MatSetValue(A,0,5,U_5(X,Y),INSERT_VALUES);
  MatSetValue(A,0,6,U_6(X,Y),INSERT_VALUES);

  //Boundary condition 2 (outer midplane curvature)
  VecSetValue(sols,1,-U_P_YY(X,Y)-N_1*U_P_X(X,Y),INSERT_VALUES);
  MatSetValue(A,1,0,U_0_YY(X,Y)+N_1*U_0_X(X,Y),INSERT_VALUES);
  MatSetValue(A,1,1,U_1_YY(X,Y)+N_1*U_1_X(X,Y),INSERT_VALUES);
  MatSetValue(A,1,2,U_2_YY(X,Y)+N_1*U_2_X(X,Y),INSERT_VALUES);
  MatSetValue(A,1,3,U_3_YY(X,Y)+N_1*U_3_X(X,Y),INSERT_VALUES);
  MatSetValue(A,1,4,U_4_YY(X,Y)+N_1*U_4_X(X,Y),INSERT_VALUES);
  MatSetValue(A,1,5,U_5_YY(X,Y)+N_1*U_5_X(X,Y),INSERT_VALUES);
  MatSetValue(A,1,6,U_6_YY(X,Y)+N_1*U_6_X(X,Y),INSERT_VALUES);

  //Boundary condition 3 (zero inner midplane flux)
  X=1.0-inv_asp;
  VecSetValue(sols,2,-U_P(X,Y),INSERT_VALUES);
  MatSetValue(A,2,0,U_0(X,Y),INSERT_VALUES);
  MatSetValue(A,2,1,U_1(X,Y),INSERT_VALUES);
  MatSetValue(A,2,2,U_2(X,Y),INSERT_VALUES);
  MatSetValue(A,2,3,U_3(X,Y),INSERT_VALUES);
  MatSetValue(A,2,4,U_4(X,Y),INSERT_VALUES);
  MatSetValue(A,2,5,U_5(X,Y),INSERT_VALUES);
  MatSetValue(A,2,6,U_6(X,Y),INSERT_VALUES);

  //Boundary condition 4 (inner midplane curvature)
  VecSetValue(sols,3,-U_P_YY(X,Y)-N_2*U_P_X(X,Y),INSERT_VALUES);
  MatSetValue(A,3,0,U_0_YY(X,Y)+N_2*U_0_X(X,Y),INSERT_VALUES);
  MatSetValue(A,3,1,U_1_YY(X,Y)+N_2*U_1_X(X,Y),INSERT_VALUES);
  MatSetValue(A,3,2,U_2_YY(X,Y)+N_2*U_2_X(X,Y),INSERT_VALUES);
  MatSetValue(A,3,3,U_3_YY(X,Y)+N_2*U_3_X(X,Y),INSERT_VALUES);
  MatSetValue(A,3,4,U_4_YY(X,Y)+N_2*U_4_X(X,Y),INSERT_VALUES);
  MatSetValue(A,3,5,U_5_YY(X,Y)+N_2*U_5_X(X,Y),INSERT_VALUES);
  MatSetValue(A,3,6,U_6_YY(X,Y)+N_2*U_6_X(X,Y),INSERT_VALUES);

  //Boundary condition 5 (high point flux)
  X=1.0-triang*inv_asp; Y=inv_asp*elong;
  VecSetValue(sols,4,-U_P(X,Y),INSERT_VALUES);
  MatSetValue(A,4,0,U_0(X,Y),INSERT_VALUES);
  MatSetValue(A,4,1,U_1(X,Y),INSERT_VALUES);
  MatSetValue(A,4,2,U_2(X,Y),INSERT_VALUES);
  MatSetValue(A,4,3,U_3(X,Y),INSERT_VALUES);
  MatSetValue(A,4,4,U_4(X,Y),INSERT_VALUES);
  MatSetValue(A,4,5,U_5(X,Y),INSERT_VALUES);
  MatSetValue(A,4,6,U_6(X,Y),INSERT_VALUES);

  //Boundary condition 6 (high point slope)
  VecSetValue(sols,5,-U_P_X(X,Y),INSERT_VALUES);
  MatSetValue(A,5,0,U_0_X(X,Y),INSERT_VALUES);
  MatSetValue(A,5,1,U_1_X(X,Y),INSERT_VALUES);
  MatSetValue(A,5,2,U_2_X(X,Y),INSERT_VALUES);
  MatSetValue(A,5,3,U_3_X(X,Y),INSERT_VALUES);
  MatSetValue(A,5,4,U_4_X(X,Y),INSERT_VALUES);
  MatSetValue(A,5,5,U_5_X(X,Y),INSERT_VALUES);
  MatSetValue(A,5,6,U_6_X(X,Y),INSERT_VALUES);

  //Boundary condition 7 (high point curvature)
  VecSetValue(sols,6,-U_P_XX(X,Y)-N_3*U_P_Y(X,Y),INSERT_VALUES);
  MatSetValue(A,6,0,U_0_XX(X,Y)+N_3*U_0_Y(X,Y),INSERT_VALUES);
  MatSetValue(A,6,1,U_1_XX(X,Y)+N_3*U_1_Y(X,Y),INSERT_VALUES);
  MatSetValue(A,6,2,U_2_XX(X,Y)+N_3*U_2_Y(X,Y),INSERT_VALUES);
  MatSetValue(A,6,3,U_3_XX(X,Y)+N_3*U_3_Y(X,Y),INSERT_VALUES);
  MatSetValue(A,6,4,U_4_XX(X,Y)+N_3*U_4_Y(X,Y),INSERT_VALUES);
  MatSetValue(A,6,5,U_5_XX(X,Y)+N_3*U_5_Y(X,Y),INSERT_VALUES);
  MatSetValue(A,6,6,U_6_XX(X,Y)+N_3*U_6_Y(X,Y),INSERT_VALUES);

  ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  ierr = VecAssemblyBegin(sols);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(sols);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(coeffs);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(coeffs);CHKERRQ(ierr);

  KSPCreate(PETSC_COMM_SELF,&ksp);
  KSPSetOperators(ksp,A,A);
  KSPSetFromOptions(ksp);
  KSPSolve(ksp,sols,coeffs);
  KSPDestroy(&ksp);

  int vals_arr[7]={0,1,2,3,4,5,6};
  PetscScalar coeffs_out[7];
  ierr = VecGetValues(coeffs,7,vals_arr,coeffs_out);
  c0=PetscRealPart(coeffs_out[0]);
  c1=PetscRealPart(coeffs_out[1]);
  c2=PetscRealPart(coeffs_out[2]);
  c3=PetscRealPart(coeffs_out[3]);
  c4=PetscRealPart(coeffs_out[4]);
  c5=PetscRealPart(coeffs_out[5]);
  c6=PetscRealPart(coeffs_out[6]);

  for(int iii=0;iii<7;iii++){std::cout << " c_" << iii << "   : " << coeffs_out[iii] << std::endl;}



  double *flux_interp=new double[N_interp];
  double flux_max=std::abs( solov_flux( 1.0 ) );  //This is flux value at midpoint (with param=X). Note it should be <0.

  double rad_max = flux_max * flux_0 ;
  double rad_min = rad_low_rat * flux_max * flux_0 ;
  assert( ( 0.5 * R_0 * R_0 * geom.B_0 * geom.B_0 / solov_A ) > rad_max ) ;
  
  fill_rad( equil, geom, rad_max, rad_min ) ;
  rad_to_psi( equil, geom.N_psi, geom.N_interp );

  

  for(int iii=0;iii<N_interp;iii++){ flux_interp[iii] = ( equil->psi_interp[iii] / flux_0 ) - flux_max ; }

  /*
  //Set up flux vals to get even spacing in rad_var=sqrt(psi) with flux zero at midpoint  
  double flux_spacing=( sqrt( std::abs( flux_max_trunc ) ) - sqrt( flux_min+(flux_max-flux_max_trunc) ) ) / ( N_interp - 1 );
  for(int iii=0;iii<N_interp;iii++){flux_interp[iii]=(sqrt(flux_min+(flux_max-flux_max_trunc))+iii*flux_spacing)*(sqrt(flux_min+(flux_max-flux_max_trunc))+iii*flux_spacing)-flux_max;}
  */
  int half_theta;  if(N_theta%2==0){half_theta=N_theta/2;}else{half_theta=(N_theta-1)/2;}


  double tol = 0.25 * pi ;
  double par_min,par_max;

  equil->maj_rad=new double[geom.N_interp*geom.N_theta];
  equil->height=new double[geom.N_interp*geom.N_theta];
  equil->dens=new double[geom.N_interp*geom.N_theta];
  equil->pres=new double[geom.N_interp*geom.N_theta];
  equil->f_psi=new double[geom.N_interp*geom.N_theta];

  //int digits = std::numeric_limits<double>::digits; // Maximum possible binary digits of doubles
  //std::cout << "Digits are : " << pow(0.5,digits) << std::endl;

  /*for(int iii=0;iii<11*N_interp;iii++)
    {
    theta=geom.theta_grid[5]; param_is_X=false;
    par_min=0.0;
    par_max=inv_asp*elong*sin(theta);
    std::cout << "Flux : " << solov_flux(iii*par_max/(8*N_interp-1)) << std::endl;
    }*/

  /*theta=equil->theta_grid[0];
  stupid_bool=true;
  std::cout << "Flux is : " << solov_flux(1.0+0.5*cos(theta)) << std::endl;
  stupid_bool=false;

  double guess=1.0+0.5*cos(theta);
  double nr_point;
  nr_point=newt_raph(solov_flux,guess,1000);
  std::cout << "nr_point is : " << nr_point << " value : " << solov_flux(nr_point) << "  " << solov_flux(nr_point-1.2e-16)-solov_flux(nr_point) << "  " << solov_flux(nr_point+1.2e-16)-solov_flux(nr_point) << std::endl;*/
  
  double ratio=1.001;

  for(int iii=0;iii<=half_theta;iii++)
    {
      theta=equil->theta_grid[iii];

      std::cout << " theta is : " << theta << std::endl;
      
      if(std::abs(theta+del_zero*sin(theta)-0.5*pi)<tol)
	{ //Y is param
	  param_is_X=false;

	  std::cout << "Y is param" << std::endl;
	  
	  par_min=0.0;
	  par_max=inv_asp*elong*sin(theta);

	  par_max=expand_bracket(solov_flux,par_min,par_max,2.0*par_max,ratio);  //Expand domain slightly


	  for(int jjj=0;jjj<N_interp;jjj++){

	    equil->height[jjj*N_theta+iii]=bisect(solov_flux,par_min,par_max,flux_interp[jjj]);
	    equil->maj_rad[jjj*N_theta+iii]=1.0+(equil->height[jjj*N_theta+iii]*cos(theta+del_zero*sin(theta))/(elong*sin(theta)));
	  }

	}
      else
	{ //X is param
	  param_is_X=true;
	  std::cout << "X is param" << std::endl;
	  
	  if(theta+del_zero*sin(theta)<0.5*pi){
	    par_min=1.0;
	    par_max=1.0+inv_asp*cos(theta+del_zero*sin(theta));

	    par_max=expand_bracket(solov_flux,par_min,par_max,2.0*(par_max-par_min),ratio); //Expand domain slightly
	    //par_max=expand_bracket(solov_flux,par_min,par_max,0.25*(par_max-par_min),2.0-ratio);
	  }
	  else{
	    par_min=1.0+inv_asp*cos(theta+del_zero*sin(theta)); 
	    par_max=1.0;

	    par_min=expand_bracket(solov_flux,par_min,par_max,2.0*(par_max-par_min),2.0-ratio,false); //Expand domain slightly
	  }


	  for(int jjj=0;jjj<N_interp;jjj++){
	    
	    equil->maj_rad[jjj*N_theta+iii]=bisect(solov_flux,par_min,par_max,flux_interp[jjj]);
	    equil->height[jjj*N_theta+iii]=(elong*sin(theta)/cos(theta+del_zero*sin(theta)))*(equil->maj_rad[jjj*N_theta+iii]-1.0);
	  }
	}
    }

  

  //Fill in second half of R,Z matrices, normalise by R_0
  for(int iii=half_theta+1;iii<N_theta;iii++){
    for(int jjj=0;jjj<N_interp;jjj++){
      equil->maj_rad[jjj*N_theta+iii]=equil->maj_rad[jjj*N_theta+(N_theta-iii)];
      equil->height[jjj*N_theta+iii]=-equil->height[jjj*N_theta+(N_theta-iii)];
      }}

  for(int iii=0;iii<N_interp*N_theta;iii++){
    equil->maj_rad[iii] *= R_0 ;
    equil->height[iii] *= R_0 ;
  }


  //Convert flux to psi
  //for(int iii=0;iii<N_interp;iii++){equil->psi_interp[iii]=(flux_interp[iii]+flux_max)*flux_0;}
  //for(int iii=0;iii<N_psi;iii++){equil->psi_grid[iii]=equil->psi_interp[iii*(geom.num_quad+1)];}

  //for(int iii=0;iii<N_interp;iii++){std::cout << geom.psi_interp[iii] << std::endl;}

  //Fill in f_psi,pres,dens
  for( int iii=0 ; iii < N_interp ; iii++ ){ equil->dens[iii*N_theta] = dens( equil->psi_interp[iii] , flux_0 , geom ) ; }
  for( int iii=0 ; iii < N_interp ; iii++ ){ for( int jjj=1 ; jjj < N_theta ; jjj++ ){ equil->dens[iii*N_theta+jjj] = equil->dens[iii*N_theta] ; }}

  double B_0,beta_0;
  B_0=geom.B_0;

  beta_0 = 2.0 * solov_C * rad_max / ( B_0 * B_0 ) ; //beta_0 set by zero-pressure b.c.

  for(int iii=0;iii<N_interp;iii++){
    for(int jjj=0;jjj<N_theta;jjj++){
      equil->f_psi[iii*N_theta+jjj]=sqrt(R_0*B_0*R_0*B_0-2.0*solov_A*equil->psi_interp[iii]);
      equil->pres[iii*N_theta+jjj]=(1.0/mu_0)*(0.5*beta_0*B_0*B_0-solov_C*equil->psi_interp[iii]);
    }}
  

  delete[] flux_interp; flux_interp=NULL;
  
  return 0;
}

#include "constants.h"


int build_matrices(Mat mats[] , geom_shape geom , equil_fields equil , matrix_coeffs *coeffs)
{
  int ierr;
  
  int num_tasks;
  ierr = MPI_Comm_size ( MPI_COMM_WORLD, &num_tasks );CHKERRQ(ierr);
  int remain=(geom.m_range * geom.m_range) % num_tasks;
  int base=(geom.m_range * geom.m_range - remain) / num_tasks;

  int my_rank;
  MPI_Comm_rank(PETSC_COMM_WORLD,&my_rank);

  int my_low;
  int my_high;
  if(my_rank < remain){
    my_low = my_rank * ( base + 1 ) ;
    my_high = ( my_rank + 1 ) * ( base + 1 ) - 1 ; 
  }
  else{
    my_low = remain * ( base + 1 ) + ( my_rank - remain ) * base ;
    my_high = remain * ( base + 1 ) + ( my_rank + 1 - remain ) * base - 1 ;
  }

  //std::cout << "Process " << my_rank << " has a range from " << my_low << " to " << my_high << std::endl;
  
  int row_offset=0, col_offset=0 ;

  gq_mod gq_mod;

  if(geom.quad_type=="gq"){}
  else if(geom.quad_type=="gq_mod"){ build_gq_mod(&gq_mod,geom,equil); }
      
  for(int main_mod=geom.m_min ; main_mod<geom.m_max ; main_mod++)
    {
      col_offset=0;
      
      for(int sec_mod=geom.m_min ; sec_mod<geom.m_max ; sec_mod++)
	{

	  int map_mod = (main_mod - geom.m_min) * geom.m_range + (sec_mod - geom.m_min );
	  
	  if(my_high<my_low){} //Occurs in case my_rank > number of sub matrices
	  else if(my_low <= map_mod && map_mod <= my_high)
	    {
	      build_submatrices(mats,equil,&geom,coeffs,&gq_mod,main_mod,sec_mod,row_offset,col_offset);
	    }
	  else{}

	  if(sec_mod==0){col_offset+=geom.dims_loc[0];}
	  else if(std::abs(sec_mod)==1){col_offset+=geom.dims_loc[1];}
	  else{col_offset+=geom.dims_loc[2];}
	  
	}

      if(main_mod==0){row_offset+=geom.dims_loc[0];}
      else if(std::abs(main_mod)==1){row_offset+=geom.dims_loc[1];}
      else{row_offset+=geom.dims_loc[2];}

      std::cout << "main_mod: " << main_mod << std::endl;
      std::cout << "tor_mod: " << geom.tor_mod << std::endl;
    }
  std::cout << "Finished filling matrix : " << my_rank << std::endl; 
  
  if(geom.quad_type=="gq"){}
  else if(geom.quad_type=="gq_mod"){ delete_gq_mod(&gq_mod,geom.N_psi); }

  
  ierr = MatAssemblyBegin(mats[0],MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(mats[0],MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
 
  ierr = MatAssemblyBegin(mats[1],MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(mats[1],MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  if( geom.Hall_on == 1 ){
    ierr = MatAssemblyBegin(mats[2],MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatAssemblyEnd(mats[2],MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  }
}

void build_gq_mod(gq_mod *gq_mod,geom_shape geom, equil_fields eq)
{
  init_gq_mod(gq_mod,geom,eq);

  for(int iii=0 ; iii < geom.N_psi - 1 ; iii++){
    solve_weights(gq_mod->weights[iii],eq.rad_var[iii+1],eq.rad_var[iii],geom.num_quad,gq_mod->num_div) ;
  }
}

//Lagrange polynomials + derivatives
const double LP_1(double x){return 1.0;}
const double LP_2(double x){return x;}
const double LP_3(double x){return 0.5*(3.0*x*x-1.0);}
const double LP_4(double x){return 0.5*x*(5.0*x*x-3.0);}
const double LP_5(double x){return 0.125*(35.0*x*x*x*x-30.0*x*x+3.0);}
const double LP_6(double x){return 0.125*x*(63.0*x*x*x*x-70.0*x*x+15.0);}


const double dLP_1(double x){return 0.0;}
const double dLP_2(double x){return 1.0;}
const double dLP_3(double x){return 3.0*x;}
const double dLP_4(double x){return 0.5*(15.0*x*x-3.0);}
const double dLP_5(double x){return 0.125*(140.0*x*x*x-60.0*x);}
const double dLP_6(double x){return 0.125*(315.0*x*x*x*x-210.0*x*x+15.0);}

const double (*LP[6])(double) = {LP_1,LP_2,LP_3,LP_4,LP_5,LP_6} ;
const double (*dLP[6])(double) = {dLP_1,dLP_2,dLP_3,dLP_4,dLP_5,dLP_6} ;

int solve_weights(double weights[],double upp, double low, int num_quad, int num_div)
{
  Mat A;
  KSP ksp;
  Vec x,b; // A*x=b
  int ierr;
  const double *gqNeval,*gqNweigh,*gqMeval,*gqMweigh;
  int M = ( num_quad / num_div ) ;

  assign_gq(gqMeval,gqMweigh,M);
  assign_gq(gqNeval,gqNweigh,num_quad);

  double sing_pnt = (upp + low ) / ( low - upp ) ;

  //Create matrix
  ierr = MatCreate(PETSC_COMM_WORLD,&A); CHKERRQ(ierr);
  ierr = MatSetSizes(A,PETSC_DECIDE,PETSC_DECIDE,num_quad,num_quad); CHKERRQ(ierr);
  ierr = MatSetFromOptions(A); CHKERRQ(ierr);
  ierr = MatSetUp(A); CHKERRQ(ierr);

  //Fill matrix with LP values
  for(int iii=0; iii < M ; iii++){
    for(int jjj=0; jjj < num_quad ; jjj++){
      for(int kkk=0; kkk < num_div ; kkk++){
	MatSetValue(A,iii+kkk*M,jjj,LP[iii](gqNeval[jjj]) / pow((sing_pnt - gqNeval[jjj]),kkk) ,INSERT_VALUES);
      }}}

  ierr = MatAssemblyBegin(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(A,MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  //Create solution vector and RHS
  ierr = VecCreate(PETSC_COMM_WORLD,&x);CHKERRQ(ierr);
  ierr = VecSetSizes(x,PETSC_DECIDE,num_quad);CHKERRQ(ierr);
  ierr = VecSetFromOptions(x);CHKERRQ(ierr);
  ierr = VecAssemblyBegin(x);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(x);CHKERRQ(ierr);

  ierr = VecCreate(PETSC_COMM_WORLD,&b);CHKERRQ(ierr);
  ierr = VecSetSizes(b,PETSC_DECIDE,num_quad);CHKERRQ(ierr);
  ierr = VecSetFromOptions(b);CHKERRQ(ierr);

  //Calc LP integrals
  double int_LP[M] ;
  int_LP[0] = 2.0 ;
  for(int iii=1 ; iii < M ; iii++){ int_LP[iii] = 0.0 ; }
  double int_LPd[M],int_LPdd[M];

  for(int iii = 0 ; iii < M ; iii++){
    int_LPd[iii] = LP[iii](sing_pnt) * log( std::abs( ( 1.0 + sing_pnt ) / (1.0 - sing_pnt) ) ) ;
    int_LPdd[iii] = LP[iii](sing_pnt) * ( 2.0 / ( sing_pnt * sing_pnt - 1.0) ) - dLP[iii](sing_pnt) * log( std::abs( ( 1.0 + sing_pnt ) / (1.0 - sing_pnt) ) ) ;
    for(int jjj = 0 ; jjj < M ; jjj++){
      int_LPd[iii] += gqMweigh[jjj] * ( ( LP[iii](gqMeval[jjj]) - LP[iii](sing_pnt) ) / ( sing_pnt - gqMeval[jjj] ) ) ;
      int_LPdd[iii] += gqMweigh[jjj] * ( ( LP[iii](gqMeval[jjj]) - LP[iii](sing_pnt) + dLP[iii](sing_pnt) * ( sing_pnt - gqMeval[jjj] )  ) / ( ( sing_pnt - gqMeval[jjj] ) * ( sing_pnt - gqMeval[jjj] ) ) ) ;
    }}

  //Fill RHS
  for(int iii=0; iii < M ; iii++){
    VecSetValue(b,iii,int_LP[iii],INSERT_VALUES);
    VecSetValue(b,iii+M,int_LPd[iii],INSERT_VALUES);
    VecSetValue(b,iii+2*M,int_LPdd[iii],INSERT_VALUES);
  }

  ierr = VecAssemblyBegin(b);CHKERRQ(ierr);
  ierr = VecAssemblyEnd(b);CHKERRQ(ierr);

  PC pc;

  //Solve linear system
  ierr = KSPCreate(PETSC_COMM_WORLD,&ksp);CHKERRQ(ierr);
  ierr = KSPSetType(ksp,KSPPREONLY);CHKERRQ(ierr);
  ierr = KSPGetPC(ksp,&pc);CHKERRQ(ierr);
  ierr = PCSetType(pc,PCLU);CHKERRQ(ierr);
  ierr = PCFactorSetMatSolverType(pc,MATSOLVERMUMPS);CHKERRQ(ierr);
  ierr = KSPSetOperators(ksp,A,A);CHKERRQ(ierr);
  //ierr = KSPSetFromOptions(ksp); CHKERRQ(ierr); //Uncomment if want to allow users to set options at runtime
  ierr = KSPSolve(ksp,b,x);CHKERRQ(ierr);

  const PetscScalar *weight_tmp ;
  ierr = VecGetArrayRead(x,&weight_tmp);CHKERRQ(ierr);

  for(int iii=0 ; iii<num_quad ; iii++){ weights[iii] = PetscRealPart(weight_tmp[iii]) ; }

  //Clean up
  ierr = KSPDestroy(&ksp);CHKERRQ(ierr);
  ierr = MatDestroy(&A);CHKERRQ(ierr);
  ierr = VecDestroy(&b);CHKERRQ(ierr);
  ierr = VecDestroy(&x);CHKERRQ(ierr);
}

void build_submatrices(Mat mats[],equil_fields eq,geom_shape* geom,matrix_coeffs* coeffs,gq_mod* gq_mod,int main_mod,int sec_mod,int row_offset,int col_offset)
{
  int N_psi=geom->N_psi;
  int N_theta=geom->N_theta;
  int tor_mod=geom->tor_mod;
  
  bool dw_dh[2];

  int dim_sub;
  fill_indices(geom->ind_perp_main,geom->ind_wedge_main,main_mod,N_psi,geom->shape_order);
  fill_indices(geom->ind_perp_sec,geom->ind_wedge_sec,sec_mod,N_psi,geom->shape_order);
  fill_dim(dim_sub,main_mod,N_psi,geom->shape_order);
  
  if(std::abs(main_mod-sec_mod)>geom->m_coup){return;}

  PetscErrorCode ierr;
  PetscInt start,end;

  //Force Matrix 
  build_coeffs(eq,*geom,coeffs,main_mod,sec_mod,"Force");
  
  //perp-perp terms
  for(int i=0;i<geom->num_var_perp;i++){
    for(int j=0;j<geom->num_var_perp;j++){
      
      select_shape(&geom->main_shape,geom->shape_order,i,true);
      select_shape(&geom->sec_shape,geom->shape_order,j,true);
      
      dw_dh[0]=true; dw_dh[1]=true;   
      fill_submatrix(mats[0],geom,gq_mod,eq.rad_interp,coeffs->f_dpdp,geom->ind_perp_main[i],geom->ind_perp_sec[j],dim_sub,row_offset,col_offset,dw_dh,ADD_VALUES);
      
      dw_dh[0]=true; dw_dh[1]=false; 
      fill_submatrix(mats[0],geom,gq_mod,eq.rad_interp,coeffs->f_dpp,geom->ind_perp_main[i],geom->ind_perp_sec[j],dim_sub,row_offset,col_offset,dw_dh,ADD_VALUES);
       
      dw_dh[0]=false; dw_dh[1]=true;
      fill_submatrix(mats[0],geom,gq_mod,eq.rad_interp,coeffs->f_pdp,geom->ind_perp_main[i],geom->ind_perp_sec[j],dim_sub,row_offset,col_offset,dw_dh,ADD_VALUES);
      
      dw_dh[0]=false; dw_dh[1]=false;
      fill_submatrix(mats[0],geom,gq_mod,eq.rad_interp,coeffs->f_pp,geom->ind_perp_main[i],geom->ind_perp_sec[j],dim_sub,row_offset,col_offset,dw_dh,ADD_VALUES);     
    }}

  //perp-wedge
  for(int i=0;i<geom->num_var_perp;i++){
    for(int j=0;j<geom->num_var_wedge;j++){
      
      select_shape(&geom->main_shape,geom->shape_order,i,true);
      select_shape(&geom->sec_shape,geom->shape_order,j,false);

      dw_dh[0]=true; dw_dh[1]=false;
      fill_submatrix(mats[0],geom,gq_mod,eq.rad_interp,coeffs->f_dpw,geom->ind_perp_main[i],geom->ind_wedge_sec[j],dim_sub,row_offset,col_offset,dw_dh,ADD_VALUES);
      
      dw_dh[0]=false; dw_dh[1]=false;
      fill_submatrix(mats[0],geom,gq_mod,eq.rad_interp,coeffs->f_pw,geom->ind_perp_main[i],geom->ind_wedge_sec[j],dim_sub,row_offset,col_offset,dw_dh,ADD_VALUES);
    }}

  //wedge-perp
  for(int i=0;i<geom->num_var_wedge;i++){for(int j=0;j<geom->num_var_perp;j++){
      
      select_shape(&geom->main_shape,geom->shape_order,i,false);
      select_shape(&geom->sec_shape,geom->shape_order,j,true);

      dw_dh[0]=false; dw_dh[1]=true; 
      fill_submatrix(mats[0],geom,gq_mod,eq.rad_interp,coeffs->f_wdp,geom->ind_wedge_main[i],geom->ind_perp_sec[j],dim_sub,row_offset,col_offset,dw_dh,ADD_VALUES);
      
      dw_dh[0]=false; dw_dh[1]=false; 
      fill_submatrix(mats[0],geom,gq_mod,eq.rad_interp,coeffs->f_wp,geom->ind_wedge_main[i],geom->ind_perp_sec[j],dim_sub,row_offset,col_offset,dw_dh,ADD_VALUES);
    }}

  //wedge-wedge
  for(int i=0;i<geom->num_var_wedge;i++){for(int j=0;j<geom->num_var_wedge;j++){
      
      select_shape(&geom->main_shape,geom->shape_order,i,false);
      select_shape(&geom->sec_shape,geom->shape_order,j,false);

      dw_dh[0]=false; dw_dh[1]=false; 
      fill_submatrix(mats[0],geom,gq_mod,eq.rad_interp,coeffs->f_ww,geom->ind_wedge_main[i],geom->ind_wedge_sec[j],dim_sub,row_offset,col_offset,dw_dh,ADD_VALUES);
    }}


  
  int counter;
  if( geom->Hall_on == 1 ){

    //Hall Matrix 
    build_coeffs(eq,*geom,coeffs,main_mod,sec_mod,"Hall");
  
    //perp-perp terms
    for(int i=0;i<geom->num_var_perp;i++){
      for(int j=0;j<geom->num_var_perp;j++){
      
	select_shape(&geom->main_shape,geom->shape_order,i,true);
	select_shape(&geom->sec_shape,geom->shape_order,j,true);

	dw_dh[0]=true; dw_dh[1]=false;
	fill_submatrix(mats[1],geom,gq_mod,eq.rad_interp,coeffs->f_dpp,geom->ind_perp_main[i],geom->ind_perp_sec[j],dim_sub,row_offset,col_offset,dw_dh,ADD_VALUES);
      
	dw_dh[0]=false; dw_dh[1]=false;
	fill_submatrix(mats[1],geom,gq_mod,eq.rad_interp,coeffs->f_pp,geom->ind_perp_main[i],geom->ind_perp_sec[j],dim_sub,row_offset,col_offset,dw_dh,ADD_VALUES);     
      }}

    //perp-wedge
    for(int i=0;i<geom->num_var_perp;i++){
      for(int j=0;j<geom->num_var_wedge;j++){
      
	select_shape(&geom->main_shape,geom->shape_order,i,true);
	select_shape(&geom->sec_shape,geom->shape_order,j,false);

	dw_dh[0]=true; dw_dh[1]=true;
	fill_submatrix(mats[1],geom,gq_mod,eq.rad_interp,coeffs->f_dpdw,geom->ind_perp_main[i],geom->ind_wedge_sec[j],dim_sub,row_offset,col_offset,dw_dh,ADD_VALUES);

	dw_dh[0]=true; dw_dh[1]=false;
	fill_submatrix(mats[1],geom,gq_mod,eq.rad_interp,coeffs->f_dpw,geom->ind_perp_main[i],geom->ind_wedge_sec[j],dim_sub,row_offset,col_offset,dw_dh,ADD_VALUES);

	dw_dh[0]=false; dw_dh[1]=true;
	fill_submatrix(mats[1],geom,gq_mod,eq.rad_interp,coeffs->f_pdw,geom->ind_perp_main[i],geom->ind_wedge_sec[j],dim_sub,row_offset,col_offset,dw_dh,ADD_VALUES);

	dw_dh[0]=false; dw_dh[1]=false;
	fill_submatrix(mats[1],geom,gq_mod,eq.rad_interp,coeffs->f_pw,geom->ind_perp_main[i],geom->ind_wedge_sec[j],dim_sub,row_offset,col_offset,dw_dh,ADD_VALUES);
      }}

    //wedge-perp
    for(int i=0;i<geom->num_var_wedge;i++){
      for(int j=0;j<geom->num_var_perp;j++){
      
	select_shape(&geom->main_shape,geom->shape_order,i,false);
	select_shape(&geom->sec_shape,geom->shape_order,j,true);
      
	dw_dh[0]=false; dw_dh[1]=false; 
	fill_submatrix(mats[1],geom,gq_mod,eq.rad_interp,coeffs->f_wp,geom->ind_wedge_main[i],geom->ind_perp_sec[j],dim_sub,row_offset,col_offset,dw_dh,ADD_VALUES);
      }}

    //wedge-wedge
    for(int i=0;i<geom->num_var_wedge;i++){
      for(int j=0;j<geom->num_var_wedge;j++){
      
	select_shape(&geom->main_shape,geom->shape_order,i,false);
	select_shape(&geom->sec_shape,geom->shape_order,j,false);

	dw_dh[0]=false; dw_dh[1]=true; 
	fill_submatrix(mats[1],geom,gq_mod,eq.rad_interp,coeffs->f_wdw,geom->ind_wedge_main[i],geom->ind_wedge_sec[j],dim_sub,row_offset,col_offset,dw_dh,ADD_VALUES);
	
	dw_dh[0]=false; dw_dh[1]=false; 
	fill_submatrix(mats[1],geom,gq_mod,eq.rad_interp,coeffs->f_ww,geom->ind_wedge_main[i],geom->ind_wedge_sec[j],dim_sub,row_offset,col_offset,dw_dh,ADD_VALUES);
      }}

    
    counter=2;
  }
  else{ counter=1; }



  
  //Inertial Matrix 
  build_coeffs(eq,*geom,coeffs,main_mod,sec_mod,"Inertial");

  //perp-perp terms
  for(int i=0;i<geom->num_var_perp;i++){
    for(int j=0;j<geom->num_var_perp;j++){

      select_shape(&geom->main_shape,geom->shape_order,i,true);
      select_shape(&geom->sec_shape,geom->shape_order,j,true);

      dw_dh[0]=false; dw_dh[1]=false;
      fill_submatrix(mats[counter],geom,gq_mod,eq.rad_interp,coeffs->f_pp,geom->ind_perp_main[i],geom->ind_perp_sec[j],dim_sub,row_offset,col_offset,dw_dh,ADD_VALUES);
    }}

  //wedge-wedge
  for(int i=0;i<geom->num_var_wedge;i++){
    for(int j=0;j<geom->num_var_wedge;j++){

      select_shape(&geom->main_shape,geom->shape_order,i,false);
      select_shape(&geom->sec_shape,geom->shape_order,j,false);

      dw_dh[0]=false; dw_dh[1]=false;
      fill_submatrix(mats[counter],geom,gq_mod,eq.rad_interp,coeffs->f_ww,geom->ind_wedge_main[i],geom->ind_wedge_sec[j],dim_sub,row_offset,col_offset,dw_dh,ADD_VALUES);
    }}

}

void build_coeffs(equil_fields eq,geom_shape geom,matrix_coeffs* coeffs,int main_mod,int sec_mod, std::string which_mat)
{
  int N_interp=geom.N_interp;
  int N_theta=geom.N_theta;

  int half_size=geom.fourier_size_sym-1;
  int fourier_size_sym=geom.fourier_size_sym;
  int fourier_size_full=geom.fourier_size_full;

  std::string deriv_order=geom.deriv_order;

  if( geom.analytical_type == "cylinder_theta" ){
    main_mod = static_cast<int>( geom.tor_mod );
    geom.tor_mod = static_cast<double>( sec_mod ) ;
    sec_mod = main_mod;
  }
  
  int m_diff=main_mod-sec_mod;

  std::complex<double> *result = new std::complex<double>[N_interp];

  double *temp_1 = new double[N_interp*N_theta];
  double *temp_2 = new double[N_interp*N_theta];
  double *temp_3 = new double[N_interp*N_theta];
  

  std::complex<double> *fourier_full_1 = new std::complex<double> [N_interp*fourier_size_full];
  std::complex<double> *fourier_sym_1 = new std::complex<double> [N_interp*fourier_size_sym];
  std::complex<double> *fourier_sym_2 = new std::complex<double> [N_interp*fourier_size_sym];
  std::complex<double> *fourier_sym_3 = new std::complex<double> [N_interp*fourier_size_sym];

  std::complex<double> *input_funcs_2[2];
  std::complex<double> *input_funcs_3[3];
  std::complex<double> *input_funcs_4[4];

  double *input_real_2[2];
  double *input_real_3[3];


  bool fg_real[2]={true,true};

  if( which_mat == "Inertial" )
    {
      // f_pp
      for(int iii=0;iii<N_interp*N_theta;iii++){temp_1[iii]= eq.dens[iii] / (eq.jacob[iii] * eq.g_pp[iii]);}
      calc_fourier_sym(temp_1,fourier_sym_1,N_interp,N_theta);  
      if(m_diff>=0){  for(int iii=0 ; iii<N_interp ; iii++){ coeffs->f_pp[iii] = mu_0 * fourier_sym_1[ iii*fourier_size_sym + m_diff ];}}
      else{ for(int iii=0 ; iii<N_interp ; iii++){ coeffs->f_pp[iii] = mu_0 * std::conj(fourier_sym_1[ iii*fourier_size_sym - m_diff ]);}}


      // f_ww
      for(int iii=0;iii<N_interp*N_theta;iii++){temp_1[iii] = eq.jacob[iii] * eq.dens[iii] ;}
      calc_fourier_sym(temp_1,fourier_sym_1,N_interp,N_theta); 
      if(m_diff>=0){for(int iii=0;iii<N_interp;iii++){ coeffs->f_ww[iii] = mu_0 * fourier_sym_1[iii*fourier_size_sym+m_diff];}}
      else{for(int iii=0;iii<N_interp;iii++){ coeffs->f_ww[iii] = mu_0 * std::conj(fourier_sym_1[iii*fourier_size_sym-m_diff]);}}

      //Change equation variable to s = sqrt(psi)
      double change_var;
      for(int iii=0;iii<N_interp;iii++){
	change_var=1.0/eq.drad_dpsi[iii];
    
	coeffs->f_pp[iii]*=change_var;
	coeffs->f_ww[iii]*=change_var;
      }
    }
  else if( which_mat == "Hall" )
    {
      //Build f_dpp 
      for(int iii=0;iii<N_interp*N_theta;iii++){temp_1[iii] = 1.0 / ( eq.jacob[iii] * eq.g_pp[iii] ) ;}
      input_real_2[0]=eq.g_pp; input_real_2[1]=temp_1;
      dwedge_xi(input_real_2,result,fg_real,geom,eq,main_mod,sec_mod);
      for(int iii=0;iii<N_interp;iii++){coeffs->f_dpp[iii] = result[iii];}

      //Build f_pp
      input_real_3[0]=eq.cc_1; input_real_3[1]=eq.g_pp; input_real_3[2]=temp_1;
      dpol_dwedge_xi(input_real_3,result,geom,eq,main_mod,sec_mod);
      for(int iii=0;iii<N_interp;iii++){coeffs->f_pp[iii] = -result[iii] ;}

      input_real_2[0]=eq.curv_psi; input_real_2[1]=temp_1;
      dwedge_xi(input_real_2,result,fg_real,geom,eq,main_mod,sec_mod);
      for(int iii=0;iii<N_interp;iii++){coeffs->f_pp[iii] += result[iii];}
      
      input_real_2[0]=eq.cc_4; input_real_2[1]=temp_1;
      dpar_xi(input_real_2,result,fg_real,geom,eq,main_mod,sec_mod);
      for(int iii=0;iii<N_interp;iii++){coeffs->f_pp[iii] -= result[iii];}

      //Build f_dpdw   
      for(int iii=0;iii<N_interp*N_theta;iii++){temp_1[iii] = sqrt( eq.g_pp[iii] * eq.mag_sq[iii] ) ;}
      calc_fourier_sym(temp_1,fourier_sym_1,N_interp,N_theta);
      if(m_diff>=0){for(int iii=0;iii<N_interp;iii++){coeffs->f_dpdw[iii] = -fourier_sym_1[iii*fourier_size_sym+m_diff];}}
      else{for(int iii=0;iii<N_interp;iii++){coeffs->f_dpdw[iii] = -std::conj(fourier_sym_1[iii*fourier_size_sym-m_diff]);}}

      //Build f_dpw
      for(int iii=0;iii<N_interp*N_theta;iii++){temp_1[iii] = eq.jacob[iii] * sqrt( eq.g_pp[iii] / eq.mag_sq[iii] ) ;}
        //psi derivative
      deriv_1d(temp_2,temp_1,eq.rad_interp,N_interp,N_theta,true,deriv_order);
      for(int iii=0;iii<N_interp;iii++){for(int jjj=0;jjj<N_theta;jjj++){temp_2[iii*N_theta+jjj] *= eq.drad_dpsi[iii] ;}}
        //
      for(int iii=0;iii<N_interp*N_theta;iii++){temp_2[iii] *= eq.cc_0[iii] ;}
      calc_fourier_sym(temp_2,fourier_sym_1,N_interp,N_theta);
      if(m_diff>=0){for(int iii=0;iii<N_interp;iii++){coeffs->f_dpw[iii] = -fourier_sym_1[iii*fourier_size_sym+m_diff];}}
      else{for(int iii=0;iii<N_interp;iii++){coeffs->f_dpw[iii] = -std::conj(fourier_sym_1[iii*fourier_size_sym-m_diff]);}}

      for(int iii=0;iii<N_interp*N_theta;iii++){temp_1[iii] = eq.cc_1[iii] * eq.jacob[iii] *  sqrt( eq.g_pp[iii] / eq.mag_sq[iii] ) ;}
      calc_fourier_sym(temp_1,fourier_sym_1,N_interp,N_theta);
      calc_fourier_sym(eq.cc_0,fourier_sym_2,N_interp,N_theta);
      input_funcs_2[0]=fourier_sym_2; input_funcs_2[1]=fourier_sym_1; 
      for(int iii=0;iii<N_interp;iii++){coeffs->f_dpw[iii] -= dpol_xi(input_funcs_2,fg_real,main_mod,sec_mod,N_interp,N_theta,iii);}

      //Build f_pdw
      for(int iii=0;iii<N_interp*N_theta;iii++){temp_1[iii] = sqrt( eq.g_pp[iii] * eq.mag_sq[iii] ) ;}
      calc_fourier_sym(temp_1,fourier_sym_1,N_interp,N_theta);
      calc_fourier_sym(eq.cc_1,fourier_sym_2,N_interp,N_theta);
      input_funcs_2[0]=fourier_sym_2; input_funcs_2[1]=fourier_sym_1; 
      for(int iii=0;iii<N_interp;iii++){coeffs->f_pdw[iii] = dpol_xi(input_funcs_2,fg_real,main_mod,sec_mod,N_interp,N_theta,iii);}


      //Build f_pw
      for(int iii=0;iii<N_interp*N_theta;iii++){temp_1[iii] = eq.jacob[iii] * sqrt( eq.g_pp[iii] / eq.mag_sq[iii] ) ;}
        //psi derivative
      deriv_1d(temp_2,temp_1,eq.rad_interp,N_interp,N_theta,true,deriv_order);
      for(int iii=0;iii<N_interp;iii++){for(int jjj=0;jjj<N_theta;jjj++){temp_2[iii*N_theta+jjj] *= eq.drad_dpsi[iii] ;}}
        //
      for(int iii=0;iii<N_interp*N_theta;iii++){temp_2[iii] *= eq.cc_0[iii] ;} 
      calc_fourier_sym(eq.cc_1,fourier_sym_1,N_interp,N_theta);
      calc_fourier_sym(temp_2,fourier_sym_2,N_interp,N_theta);
      input_funcs_2[0]=fourier_sym_1; input_funcs_2[1]=fourier_sym_2;
      for(int iii=0;iii<N_interp;iii++){coeffs->f_pw[iii] = dpol_xi(input_funcs_2,fg_real,main_mod,sec_mod,N_interp,N_theta,iii);}
      
      calc_fourier_sym(eq.cc_0,fourier_sym_2,N_interp,N_theta);
      for(int iii=0;iii<N_interp*N_theta;iii++){temp_1[iii] = eq.cc_1[iii] * eq.jacob[iii] * sqrt( eq.g_pp[iii] / eq.mag_sq[iii] ) ;}  
      calc_fourier_sym(temp_1,fourier_sym_3,N_interp,N_theta);
      input_funcs_3[0]=fourier_sym_1; input_funcs_3[1]=fourier_sym_2; input_funcs_3[2]=fourier_sym_3;
      for(int iii=0;iii<N_interp;iii++){coeffs->f_pw[iii] += dpol_sq_xi(input_funcs_3,main_mod,sec_mod,N_interp,N_theta,iii);}

      for(int iii=0;iii<N_interp*N_theta;iii++){temp_1[iii] = eq.curv_psi[iii] * sqrt( eq.mag_sq[iii] / eq.g_pp[iii] ) ;}
      calc_fourier_sym(temp_1,fourier_sym_3,N_interp,N_theta);
      input_funcs_2[0]=fourier_sym_1; input_funcs_2[1]=fourier_sym_3;
      for(int iii=0;iii<N_interp;iii++){coeffs->f_pw[iii] += dpol_xi(input_funcs_2,fg_real,main_mod,sec_mod,N_interp,N_theta,iii);}

      deriv_1d(temp_2,eq.cc_2,eq.rad_interp,N_interp,N_theta,true,deriv_order);
      for(int iii=0;iii<N_interp;iii++){for(int jjj=0;jjj<N_theta;jjj++){temp_2[iii*N_theta+jjj] *= eq.drad_dpsi[iii] ;}}
      for(int iii=0;iii<N_interp*N_theta;iii++){temp_2[iii] *= eq.jacob[iii] * sqrt( eq.g_pp[iii] / eq.mag_sq[iii] ) ;}
      calc_fourier_sym(temp_2,fourier_sym_1,N_interp,N_theta);
      if(m_diff>=0){for(int iii=0;iii<N_interp;iii++){coeffs->f_pw[iii] += fourier_sym_1[iii*fourier_size_sym+m_diff];}}
      else{for(int iii=0;iii<N_interp;iii++){coeffs->f_pw[iii] += std::conj(fourier_sym_1[iii*fourier_size_sym-m_diff]);}}
      
      calc_fourier_sym(eq.cc_2,fourier_sym_1,N_interp,N_theta);
      for(int iii=0;iii<N_interp*N_theta;iii++){temp_1[iii] = eq.cc_1[iii] * eq.jacob[iii] * sqrt( eq.g_pp[iii] / eq.mag_sq[iii] ) ;}
      calc_fourier_sym(temp_1,fourier_sym_2,N_interp,N_theta);
      input_funcs_2[0]=fourier_sym_1; input_funcs_2[1]=fourier_sym_2;
      for(int iii=0;iii<N_interp;iii++){coeffs->f_pw[iii] -= dpol_xi(input_funcs_2,fg_real,main_mod,sec_mod,N_interp,N_theta,iii);}

      for(int iii=0;iii<N_interp*N_theta;iii++){temp_2[iii] = eq.curv_psi[iii] * eq.curv_psi[iii] * sqrt(  eq.mag_sq[iii] / eq.g_pp[iii] ) / eq.g_pp[iii] ;}
      calc_fourier_sym(temp_2,fourier_sym_2,N_interp,N_theta);
      if(m_diff>=0){for(int iii=0;iii<N_interp;iii++){coeffs->f_pw[iii] -= fourier_sym_2[iii*fourier_size_sym+m_diff];}}
      else{for(int iii=0;iii<N_interp;iii++){coeffs->f_pw[iii] -= std::conj(fourier_sym_2[iii*fourier_size_sym-m_diff]);}}

      if(geom.shear_on){

	std::cout << "Shear on: perp" << std::endl;
    
	for(int iii=0;iii<N_interp*N_theta;iii++){temp_2[iii]=1.0/eq.g_pp[iii];}
	for(int iii=0;iii<N_interp*N_theta;iii++){temp_3[iii]= sqrt( eq.g_pp[iii] / eq.mag_sq[iii] ) ;}
	input_real_3[0]=eq.ones; input_real_3[1]=temp_2; input_real_3[2]=temp_3;
	dpar_sq_xi(input_real_3,result,geom,eq,main_mod,sec_mod);
	for(int iii=0;iii<N_interp;iii++){coeffs->f_pw[iii] += result[iii];}
    
      }
 
      for(int iii=0; iii<N_interp * N_theta ;iii++){ temp_1[iii] = - ( eq.neg_shear[iii] * eq.cc_4[iii] * sqrt( eq.g_pp[iii] / eq.mag_sq[iii] ) ) ;}
      calc_fourier_sym(temp_1,fourier_sym_1,N_interp,N_theta);
      if(m_diff>=0){for(int iii=0;iii<N_interp;iii++){coeffs->f_pw[iii] += fourier_sym_1[iii*fourier_size_sym+m_diff];}}
      else{for(int iii=0;iii<N_interp;iii++){coeffs->f_pw[iii] += std::conj(fourier_sym_1[iii*fourier_size_sym-m_diff]);}}
  
      //Build f_wp
      for(int iii=0;iii<N_interp*N_theta;iii++){temp_1[iii] = eq.jacob[iii] * sqrt( eq.mag_sq[iii] * eq.g_pp[iii] ) ;}
      for(int iii=0;iii<N_interp*N_theta;iii++){temp_2[iii] = eq.g_pp[iii] / eq.mag_sq[iii] ;}
      for(int iii=0;iii<N_interp*N_theta;iii++){temp_3[iii] = 1.0 / ( eq.jacob[iii] * eq.g_pp[iii] ) ;}
      input_real_3[0]=temp_1; input_real_3[1]=temp_2; input_real_3[2]=temp_3;
      dwedge_sq_xi(input_real_3,result,geom,eq,main_mod,sec_mod);
      for(int iii=0;iii<N_interp;iii++){coeffs->f_wp[iii] = -result[iii] ;}

      if(geom.shear_on){

	std::cout << "Shear on: wedge" << std::endl;
    
	for(int iii=0;iii<N_interp*N_theta;iii++){temp_1[iii] = eq.jacob[iii] * sqrt(eq.mag_sq[iii] / eq.g_pp[iii]) ;}
	input_real_3[0]=temp_1; input_real_3[1]=temp_2; input_real_3[2]=temp_3;
	dpar_sq_xi(input_real_3,result,geom,eq,main_mod,sec_mod);
	for(int iii=0;iii<N_interp;iii++){coeffs->f_wp[iii] -= result[iii] ;}

      }

      //Build f_wdw
      for(int iii=0;iii<N_interp*N_theta;iii++){temp_1[iii] = sqrt(eq.g_pp[iii] * eq.mag_sq[iii]) * eq.jacob[iii] ;}
      for(int iii=0;iii<N_interp*N_theta;iii++){temp_2[iii] = sqrt( eq.g_pp[iii] / eq.mag_sq[iii] ) ;}
      input_real_2[0]=temp_1; input_real_2[1]=temp_2;
      dwedge_xi(input_real_2,result,fg_real,geom,eq,main_mod,sec_mod);
      for(int iii=0;iii<N_interp;iii++){coeffs->f_wdw[iii] = result[iii];}


      //Build f_ww
      for(int iii=0;iii<N_interp*N_theta;iii++){temp_2[iii] = eq.jacob[iii] * sqrt( eq.g_pp[iii] / eq.mag_sq[iii] ) ;}
      deriv_1d(temp_3,temp_2,eq.rad_interp,N_interp,N_theta,true,deriv_order);
      for(int iii=0;iii<N_interp;iii++){for(int jjj=0;jjj<N_theta;jjj++){temp_3[iii*N_theta+jjj] *= eq.drad_dpsi[iii] ;}}
      for(int iii=0;iii<N_interp*N_theta;iii++){temp_3[iii] /= eq.jacob[iii] ;}   
      input_real_2[0]=temp_1; input_real_2[1]=temp_3;
      dwedge_xi(input_real_2,result,fg_real,geom,eq,main_mod,sec_mod);
      for(int iii=0;iii<N_interp;iii++){coeffs->f_ww[iii] = result[iii];}
      
      for(int iii=0;iii<N_interp*N_theta;iii++){temp_2[iii] = 1.0 / eq.jacob[iii] ;}
      for(int iii=0;iii<N_interp*N_theta;iii++){temp_3[iii] = eq.cc_1[iii] * eq.jacob[iii] * sqrt( eq.g_pp[iii] / eq.mag_sq[iii] ) ;} 
      input_real_3[0]=temp_1; input_real_3[1]=temp_2; input_real_3[2]=temp_3;
      dwedge_dpol_xi(input_real_3,result,geom,eq,main_mod,sec_mod);
      for(int iii=0;iii<N_interp;iii++){coeffs->f_ww[iii] += result[iii] ;}

      for(int iii=0;iii<N_interp*N_theta;iii++){temp_2[iii] = eq.curv_psi[iii] / sqrt( eq.g_pp[iii] * eq.mag_sq[iii] ) ;}
      input_real_2[0]=temp_1; input_real_2[1]=temp_2;
      dwedge_xi(input_real_2,result,fg_real,geom,eq,main_mod,sec_mod);
      for(int iii=0;iii<N_interp;iii++){coeffs->f_ww[iii] += result[iii];}

      for(int iii=0;iii<N_interp*N_theta;iii++){temp_1[iii] = sqrt( eq.mag_sq[iii] / eq.g_pp[iii] ) * eq.jacob[iii] ;}
      for(int iii=0;iii<N_interp*N_theta;iii++){temp_2[iii] = eq.cc_4[iii] * sqrt( eq.g_pp[iii] / eq.mag_sq[iii] ) ;}
      input_real_2[0]=temp_1; input_real_2[1]=temp_2;
      dpar_xi(input_real_2,result,fg_real,geom,eq,main_mod,sec_mod);
      for(int iii=0;iii<N_interp;iii++){coeffs->f_ww[iii] -= result[iii];}
      
      //Change equation variable to s = sqrt(psi)
      double change_var;
      for(int iii=0;iii<N_interp;iii++){
	change_var=1.0/eq.drad_dpsi[iii];
	
	coeffs->f_pp[iii]*=change_var;
	coeffs->f_dpdw[iii]*=eq.drad_dpsi[iii];
	coeffs->f_pw[iii]*=change_var;
	coeffs->f_wp[iii]*=change_var;
	coeffs->f_ww[iii]*=change_var;
	}
      
      for(int iii=0;iii<N_interp;iii++){
	coeffs->f_dpp[iii]*= imag_unit * hall_const ;
	coeffs->f_pp[iii]*= imag_unit * hall_const ;
	coeffs->f_dpdw[iii]*= imag_unit * hall_const ;
	coeffs->f_dpw[iii]*= imag_unit * hall_const ;
	coeffs->f_pw[iii]*= imag_unit * hall_const ;
	coeffs->f_wp[iii]*= imag_unit * hall_const ;
	coeffs->f_wdw[iii]*= imag_unit * hall_const ;
	coeffs->f_ww[iii]*= imag_unit * hall_const ;
      }
      
    }
  else if( which_mat == "Force" )
    {
      //Build f_dpdp
      calc_fourier_sym(eq.cc_0,fourier_sym_1,N_interp,N_theta);
      if(m_diff>=0){for(int iii=0;iii<N_interp;iii++){coeffs->f_dpdp[iii]=fourier_sym_1[iii*fourier_size_sym+m_diff];}}
      else{for(int iii=0;iii<N_interp;iii++){coeffs->f_dpdp[iii]=std::conj(fourier_sym_1[iii*fourier_size_sym-m_diff]);}}

 
      //Build f_dpp
      calc_fourier_sym(eq.cc_0,fourier_sym_1,N_interp,N_theta);
      calc_fourier_sym(eq.cc_1,fourier_sym_2,N_interp,N_theta);  
      input_funcs_2[0]=fourier_sym_1; input_funcs_2[1]=fourier_sym_2; 
      for(int iii=0;iii<N_interp;iii++){coeffs->f_dpp[iii]=dpol_xi(input_funcs_2,fg_real,main_mod,sec_mod,N_interp,N_theta,iii);}

      calc_fourier_sym(eq.cc_2,fourier_sym_3,N_interp,N_theta);
      if(m_diff>=0){for(int iii=0;iii<N_interp;iii++){coeffs->f_dpp[iii]+=fourier_sym_3[iii*fourier_size_sym+m_diff];}}
      else{for(int iii=0;iii<N_interp;iii++){coeffs->f_dpp[iii]+=std::conj(fourier_sym_3[iii*fourier_size_sym-m_diff]);}}

      //Build f_pdp  
      input_funcs_2[0]=fourier_sym_2; input_funcs_2[1]=fourier_sym_1; 
      for(int iii=0;iii<N_interp;iii++){coeffs->f_pdp[iii]= -dpol_xi(input_funcs_2,fg_real,main_mod,sec_mod,N_interp,N_theta,iii);}

      if(m_diff>=0){for(int iii=0;iii<N_interp;iii++){coeffs->f_pdp[iii]+=fourier_sym_3[iii*fourier_size_sym+m_diff];}}
      else{for(int iii=0;iii<N_interp;iii++){coeffs->f_pdp[iii]+=std::conj(fourier_sym_3[iii*fourier_size_sym-m_diff]);}}
 
      //Build f_pp
      input_funcs_3[0]=fourier_sym_2; input_funcs_3[1]=fourier_sym_1; input_funcs_3[2]=fourier_sym_2;
      for(int iii=0;iii<N_interp;iii++){coeffs->f_pp[iii] = - dpol_sq_xi(input_funcs_3,main_mod,sec_mod,N_interp,N_theta,iii);}

      input_funcs_2[0]=fourier_sym_2; input_funcs_2[1]=fourier_sym_3;
      for(int iii=0;iii<N_interp;iii++){coeffs->f_pp[iii] -= dpol_xi(input_funcs_2,fg_real,main_mod,sec_mod,N_interp,N_theta,iii);}

      input_funcs_2[0]=fourier_sym_3; input_funcs_2[1]=fourier_sym_2;
      for(int iii=0;iii<N_interp;iii++){coeffs->f_pp[iii] += dpol_xi(input_funcs_2,fg_real,main_mod,sec_mod,N_interp,N_theta,iii);}

      //Commented out lines refer to a pressure term.
      //deriv_1d(temp_2,eq.pres,eq.rad_interp,N_interp,N_theta,true,deriv_order);
      //for(int iii=0;iii<N_interp;iii++){for(int jjj=0;jjj<N_theta;jjj++){temp_2[iii*N_theta+jjj] *= eq.drad_dpsi[iii] ;}} 
      //for(int iii=0; iii<N_interp * N_theta ;iii++){ temp_1[iii] = ( eq.cc_2[iii] * eq.curv_psi[iii] / eq.g_pp[iii] ) - ( mu_0 * temp_2[iii] * eq.cc_2[iii] / eq.mag_sq[iii] ) + ( eq.neg_shear[iii] * eq.cc_4[iii] / eq.jacob[iii] ) ;}
      for(int iii=0; iii<N_interp * N_theta ;iii++){ temp_1[iii] = ( eq.cc_2[iii] * eq.curv_psi[iii] / eq.g_pp[iii] ) + ( eq.neg_shear[iii] * eq.cc_4[iii] / eq.jacob[iii] ) ;}
      calc_fourier_sym(temp_1,fourier_sym_1,N_interp,N_theta);
      if(m_diff>=0){for(int iii=0;iii<N_interp;iii++){coeffs->f_pp[iii]+=fourier_sym_1[iii*fourier_size_sym+m_diff];}}
      else{for(int iii=0;iii<N_interp;iii++){coeffs->f_pp[iii]+=std::conj(fourier_sym_1[iii*fourier_size_sym-m_diff]);}}
  
  
      if(geom.shear_on){

	std::cout << "Shear on: perp" << std::endl;
    
	for(int iii=0;iii<N_interp*N_theta;iii++){temp_2[iii]=1.0/eq.g_pp[iii];}
	for(int iii=0;iii<N_interp*N_theta;iii++){temp_3[iii]=1.0/eq.jacob[iii];}
	input_real_3[0]=eq.ones; input_real_3[1]=temp_2; input_real_3[2]=temp_3;
	dpar_sq_xi(input_real_3,result,geom,eq,main_mod,sec_mod);
	for(int iii=0;iii<N_interp;iii++){coeffs->f_pp[iii] -= result[iii];}
    
      }
  

      //Build f_dpw
      for(int iii=0; iii<N_interp * N_theta ;iii++){ temp_1[iii] = sqrt( eq.mag_sq[iii] / eq.g_pp[iii] ) ;}
      input_real_2[0]=eq.g_pp; input_real_2[1]=temp_1;
      dwedge_xi(input_real_2,result,fg_real,geom,eq,main_mod,sec_mod);
      for(int iii=0;iii<N_interp;iii++){coeffs->f_dpw[iii] = result[iii];}
  
  
      //Build f_pw
      input_real_3[0]=eq.cc_1; input_real_3[1]=eq.g_pp; input_real_3[2]=temp_1;
      dpol_dwedge_xi(input_real_3,result,geom,eq,main_mod,sec_mod);
      for(int iii=0;iii<N_interp;iii++){coeffs->f_pw[iii] = - result[iii];}

      input_real_2[0]=eq.curv_psi; input_real_2[1]=temp_1;
      dwedge_xi(input_real_2,result,fg_real,geom,eq,main_mod,sec_mod);
      for(int iii=0;iii<N_interp;iii++){coeffs->f_pw[iii] += result[iii];}

      input_real_2[0]=eq.cc_4; input_real_2[1]=temp_1;
      dpar_xi(input_real_2,result,fg_real,geom,eq,main_mod,sec_mod);
      for(int iii=0;iii<N_interp;iii++){coeffs->f_pw[iii] -= result[iii];}
  
      //Build f_wdp
      for(int iii=0;iii<N_interp*N_theta;iii++){temp_1[iii] = sqrt(eq.g_pp[iii] * eq.mag_sq[iii]) * eq.jacob[iii] ;}
      for(int iii=0;iii<N_interp*N_theta;iii++){temp_2[iii] = 1.0 / eq.jacob[iii] ;}
      input_real_2[0]=temp_1; input_real_2[1]=temp_2;
      dwedge_xi(input_real_2,result,fg_real,geom,eq,main_mod,sec_mod);
      for(int iii=0;iii<N_interp;iii++){coeffs->f_wdp[iii] = - result[iii];}


      //Build f_wp
      input_real_3[0]=temp_1; input_real_3[1]=temp_2; input_real_3[2]=eq.cc_1;
      dwedge_dpol_xi(input_real_3,result,geom,eq,main_mod,sec_mod);
      for(int iii=0;iii<N_interp;iii++){coeffs->f_wp[iii] = - result[iii] ;}

      for(int iii=0;iii<N_interp*N_theta;iii++){temp_2[iii] = eq.cc_2[iii] / eq.mag_sq[iii] ;}
      input_real_2[0]=temp_1; input_real_2[1]=temp_2;
      dwedge_xi(input_real_2,result,fg_real,geom,eq,main_mod,sec_mod);
      for(int iii=0;iii<N_interp;iii++){coeffs->f_wp[iii] -= result[iii];}

      for(int iii=0;iii<N_interp*N_theta;iii++){temp_1[iii] = sqrt(eq.mag_sq[iii] / eq.g_pp[iii]) * eq.jacob[iii] ;}
      for(int iii=0;iii<N_interp*N_theta;iii++){temp_2[iii] = eq.cc_4[iii] / eq.jacob[iii] ;}
      input_real_2[0]=temp_1; input_real_2[1]=temp_2;
      dpar_xi(input_real_2,result,fg_real,geom,eq,main_mod,sec_mod);
      for(int iii=0;iii<N_interp;iii++){coeffs->f_wp[iii] += result[iii];}

      //Build f_ww
      for(int iii=0;iii<N_interp*N_theta;iii++){temp_1[iii] = eq.jacob[iii] * sqrt(eq.mag_sq[iii] * eq.g_pp[iii]) ;}
      for(int iii=0;iii<N_interp*N_theta;iii++){temp_2[iii] = eq.g_pp[iii] / eq.mag_sq[iii] ;}
      for(int iii=0;iii<N_interp*N_theta;iii++){temp_3[iii] = sqrt(eq.mag_sq[iii] / eq.g_pp[iii]) ;}
      input_real_3[0]=temp_1; input_real_3[1]=temp_2; input_real_3[2]=temp_3;
      dwedge_sq_xi(input_real_3,result,geom,eq,main_mod,sec_mod);
      for(int iii=0;iii<N_interp;iii++){coeffs->f_ww[iii] = - result[iii] ;}

      if(geom.shear_on){

	std::cout << "Shear on: wedge" << std::endl;
    
	for(int iii=0;iii<N_interp*N_theta;iii++){temp_1[iii] = eq.jacob[iii] * sqrt(eq.mag_sq[iii] / eq.g_pp[iii]) ;}
	input_real_3[0]=temp_1; input_real_3[1]=temp_2; input_real_3[2]=temp_3;
	dpar_sq_xi(input_real_3,result,geom,eq,main_mod,sec_mod);
	for(int iii=0;iii<N_interp;iii++){coeffs->f_ww[iii] -= result[iii] ;}

      }

      //Change equation variable to s = sqrt(psi)
      double change_var;
      for(int iii=0;iii<N_interp;iii++){
	change_var=1.0/eq.drad_dpsi[iii];
    
	coeffs->f_dpdp[iii]*=eq.drad_dpsi[iii];
	coeffs->f_pp[iii]*=change_var;
	coeffs->f_pw[iii]*=change_var;
	coeffs->f_wp[iii]*=change_var;
	coeffs->f_ww[iii]*=change_var;
      }
    }

  delete[] temp_1;
  delete[] temp_2;
  delete[] temp_3;
  delete[] fourier_full_1;
  delete[] fourier_sym_1;
  delete[] fourier_sym_2;
  delete[] fourier_sym_3;
  delete[] result;

  if( geom.analytical_type == "cylinder_theta" ){
    main_mod = static_cast<int>( geom.tor_mod );
    geom.tor_mod = static_cast<double>( sec_mod ) ;
    sec_mod = main_mod;
  }

}

void fill_submatrix(Mat M,geom_shape *geom,gq_mod* gq_mod,double rad_interp[],std::complex<double> matrix_coeffs[],int main_indices[],int sec_indices[],int dim_sub,int row_offset,int col_offset,bool dw_dh[],InsertMode insert_mode)
{

  int N_psi=geom->N_psi;
  int num_quad=geom->num_quad;
  std::string shape_order=geom->shape_order;
  
  PetscInt start,end;
  MatGetOwnershipRange(M,&start,&end);
  if(end<=row_offset || start>row_offset+dim_sub){return;} //Check that some of the rows to be filled are in the range that this thread owns

  if(geom->quad_type == "gq"){ //Using massive if statement is messy, clean up after proof of concept
    double (* shape_main[2])(double,double,double);
    double (* shape_sec[2])(double,double,double);

    fill_shapes(shape_main,dw_dh[0],geom->main_shape);
    fill_shapes(shape_sec,dw_dh[1],geom->sec_shape);

    //Select shape functions
    double (* shape1[2])(double,double,double);
    double (* shape2[2])(double,double,double);
    double (* shape3[2])(double,double,double);
    double (* shape4[2])(double,double,double);

    shape1[0]=shape_main[0]; shape1[1]=shape_sec[1];
    shape2[0]=shape_main[0]; shape2[1]=shape_sec[0];
    shape3[0]=shape_main[1]; shape3[1]=shape_sec[1];
    shape4[0]=shape_main[1]; shape4[1]=shape_sec[0];
 
    //Fill boundary values (assuming that only boundary values can be '-1' i.e. displacement not known (a priori) to vanish in interior)
    PetscScalar value1,value2,value3,value4;
  
    //Oth gridpoint i.e. first
    //if(start>row_offset+main_indices[0]){}
    //else
    {
      if(!(main_indices[0]==-1))
	{
	  if(!(sec_indices[0]==-1))
	    {
	      value3=integrate1d(matrix_coeffs,rad_interp,shape3,N_psi,num_quad,1,0);
	      MatSetValue(M,main_indices[0]+row_offset,sec_indices[0]+col_offset,value3,insert_mode);	 
	    }

	  if(!(sec_indices[1]==-1))
	    {
	      value4=integrate1d(matrix_coeffs,rad_interp,shape4,N_psi,num_quad,1,0);
	      MatSetValue(M,main_indices[0]+row_offset,sec_indices[1]+col_offset,value4,insert_mode);
	    }
	}
    }

    //if(end>row_offset+main_indices[1] && start<=row_offset+main_indices[1])
    {
      if(!(main_indices[1]==-1))
	{
	  //1st gridpoint
	  if(!(sec_indices[0]==-1))
	    {
	      value1=integrate1d(matrix_coeffs,rad_interp,shape1,N_psi,num_quad,1,0);
	      MatSetValue(M,main_indices[1]+row_offset,sec_indices[0]+col_offset,value1,insert_mode);
	    }

	  if(!(sec_indices[1]==-1))
	    {
	      value2=integrate1d(matrix_coeffs,rad_interp,shape2,N_psi,num_quad,1,0);
	      value3=integrate1d(matrix_coeffs,rad_interp,shape3,N_psi,num_quad,2,1);
	      MatSetValue(M,main_indices[1]+row_offset,sec_indices[1]+col_offset,value2+value3,insert_mode);
	    }

	  if(!(sec_indices[2]==-1))
	    {
	      value4=integrate1d(matrix_coeffs,rad_interp,shape4,N_psi,num_quad,2,1);
	      MatSetValue(M,main_indices[1]+row_offset,sec_indices[2]+col_offset,value4,insert_mode);
	    }
	}
    }

    //Fill interior values; sort out "start" value (may be greater than 2)
    int iii=2;
    while(iii<N_psi-2 ) //&& main_indices[iii]<end)
      {
	if(!(main_indices[iii]==-1))
	  {
	    if(!(sec_indices[iii-1]==-1))
	      {
		value1=integrate1d(matrix_coeffs,rad_interp,shape1,N_psi,num_quad,iii,iii-1);
		MatSetValue(M,main_indices[iii]+row_offset,sec_indices[iii-1]+col_offset,value1,insert_mode);
	      }

	    if(!(sec_indices[iii]==-1))
	      {
		value2=integrate1d(matrix_coeffs,rad_interp,shape2,N_psi,num_quad,iii,iii-1);
		value3=integrate1d(matrix_coeffs,rad_interp,shape3,N_psi,num_quad,iii+1,iii);
		MatSetValue(M,main_indices[iii]+row_offset,sec_indices[iii]+col_offset,value2+value3,insert_mode);
	      }

	    if(!(sec_indices[iii+1]==-1))
	      {
		value4=integrate1d(matrix_coeffs,rad_interp,shape4,N_psi,num_quad,iii+1,iii);
		MatSetValue(M,main_indices[iii]+row_offset,sec_indices[iii+1]+col_offset,value4,insert_mode);
	      }
	  }
      
	iii++;
      }
  

    //if(end>row_offset+main_indices[N_psi-2] && start<=row_offset+main_indices[N_psi-2])
    {
      if(!(main_indices[N_psi-2]==-1))
	{
	  //(N-2)th
	  if(!(sec_indices[N_psi-3]==-1))
	    {
	      value1=integrate1d(matrix_coeffs,rad_interp,shape1,N_psi,num_quad,N_psi-2,N_psi-3);
	      MatSetValue(M,main_indices[N_psi-2]+row_offset,sec_indices[N_psi-3]+col_offset,value1,insert_mode);
	    }

	  if(!(sec_indices[N_psi-2]==-1))
	    {
	      value2=integrate1d(matrix_coeffs,rad_interp,shape2,N_psi,num_quad,N_psi-2,N_psi-3);
	      value3=integrate1d(matrix_coeffs,rad_interp,shape3,N_psi,num_quad,N_psi-1,N_psi-2);
	      MatSetValue(M,main_indices[N_psi-2]+row_offset,sec_indices[N_psi-2]+col_offset,value2+value3,insert_mode);
	    }
  
	  if(!(sec_indices[N_psi-1]==-1))
	    {
	      value4=integrate1d(matrix_coeffs,rad_interp,shape4,N_psi,num_quad,N_psi-1,N_psi-2);
	      MatSetValue(M,main_indices[N_psi-2]+row_offset,sec_indices[N_psi-1]+col_offset,value4,insert_mode);
	    }
	}
    }
	  
    //if(end>row_offset+main_indices[N_psi-1] && start<=row_offset+main_indices[N_psi-1])
    {
      //(N-1)th i.e. final
      if(!(main_indices[N_psi-1]==-1))
	{
	  if(!(sec_indices[N_psi-2]==-1))
	    {
	      value1=integrate1d(matrix_coeffs,rad_interp,shape1,N_psi,num_quad,N_psi-1,N_psi-2);
	      MatSetValue(M,main_indices[N_psi-1]+row_offset,sec_indices[N_psi-2]+col_offset,value1,insert_mode);
	    }
      
	  if(!(sec_indices[N_psi-1]==-1))
	    {
	      value2=integrate1d(matrix_coeffs,rad_interp,shape2,N_psi,num_quad,N_psi-1,N_psi-2);
	      MatSetValue(M,main_indices[N_psi-1]+row_offset,sec_indices[N_psi-1]+col_offset,value2,insert_mode);
	    }
	}
    }
  }
  else if(geom->quad_type == "gq_mod"){
    
    void (*shape_main[2])(std::vector<double>*,double,double) ;
    void (*shape_sec[2])(std::vector<double>*,double,double) ;

    fill_shapes(shape_main,dw_dh[0],geom->main_shape);
    fill_shapes(shape_sec,dw_dh[1],geom->sec_shape);

    //Select shape functions
    void (*shape1[2])(std::vector<double>*,double,double) ;
    void (*shape2[2])(std::vector<double>*,double,double) ;
    void (*shape3[2])(std::vector<double>*,double,double) ;
    void (*shape4[2])(std::vector<double>*,double,double) ;

    shape1[0]=shape_main[0]; shape1[1]=shape_sec[1];
    shape2[0]=shape_main[0]; shape2[1]=shape_sec[0];
    shape3[0]=shape_main[1]; shape3[1]=shape_sec[1];
    shape4[0]=shape_main[1]; shape4[1]=shape_sec[0];

    //Set wt_gq flag to zero so that gq_mod is used initially
    gq_mod->wt_gq = 0 ;
 
    //Fill boundary values (assuming that only boundary values can be '-1' i.e. displacement not known (a priori) to vanish in interior)
    PetscScalar value1,value2,value3,value4;
  
    //Oth gridpoint i.e. first
    //if(start>row_offset+main_indices[0]){}
    //else
    {
      if(!(main_indices[0]==-1))
	{
	  if(!(sec_indices[0]==-1))
	    {
	      value3=integrate1d(matrix_coeffs,gq_mod,rad_interp,shape3,N_psi,num_quad,1,0);
	      MatSetValue(M,main_indices[0]+row_offset,sec_indices[0]+col_offset,value3,insert_mode);	 
	    }

	  if(!(sec_indices[1]==-1))
	    {
	      value4=integrate1d(matrix_coeffs,gq_mod,rad_interp,shape4,N_psi,num_quad,1,0);
	      MatSetValue(M,main_indices[0]+row_offset,sec_indices[1]+col_offset,value4,insert_mode);
	    }
	}
    }

    //if(end>row_offset+main_indices[1] && start<=row_offset+main_indices[1])
    {
      if(!(main_indices[1]==-1))
	{
	  //1st gridpoint
	  if(!(sec_indices[0]==-1))
	    {
	      value1=integrate1d(matrix_coeffs,gq_mod,rad_interp,shape1,N_psi,num_quad,1,0);
	      MatSetValue(M,main_indices[1]+row_offset,sec_indices[0]+col_offset,value1,insert_mode);
	    }

	  if(!(sec_indices[1]==-1))
	    {
	      value2=integrate1d(matrix_coeffs,gq_mod,rad_interp,shape2,N_psi,num_quad,1,0);
	      value3=integrate1d(matrix_coeffs,gq_mod,rad_interp,shape3,N_psi,num_quad,2,1);
	      MatSetValue(M,main_indices[1]+row_offset,sec_indices[1]+col_offset,value2+value3,insert_mode);
	    }

	  if(!(sec_indices[2]==-1))
	    {
	      value4=integrate1d(matrix_coeffs,gq_mod,rad_interp,shape4,N_psi,num_quad,2,1);
	      MatSetValue(M,main_indices[1]+row_offset,sec_indices[2]+col_offset,value4,insert_mode);
	    }
	}
    }

    //Fill interior values; sort out "start" value (may be greater than 2)
    int iii=2;
    while(iii<N_psi-2 ) //&& main_indices[iii]<end)
      {
	if(!(main_indices[iii]==-1))
	  {
	    if(!(sec_indices[iii-1]==-1))
	      {
		value1=integrate1d(matrix_coeffs,gq_mod,rad_interp,shape1,N_psi,num_quad,iii,iii-1);
		MatSetValue(M,main_indices[iii]+row_offset,sec_indices[iii-1]+col_offset,value1,insert_mode);
	      }

	    if(!(sec_indices[iii]==-1))
	      {
		value2=integrate1d(matrix_coeffs,gq_mod,rad_interp,shape2,N_psi,num_quad,iii,iii-1);
		value3=integrate1d(matrix_coeffs,gq_mod,rad_interp,shape3,N_psi,num_quad,iii+1,iii);
		MatSetValue(M,main_indices[iii]+row_offset,sec_indices[iii]+col_offset,value2+value3,insert_mode);
	      }

	    if(!(sec_indices[iii+1]==-1))
	      {
		value4=integrate1d(matrix_coeffs,gq_mod,rad_interp,shape4,N_psi,num_quad,iii+1,iii);
		MatSetValue(M,main_indices[iii]+row_offset,sec_indices[iii+1]+col_offset,value4,insert_mode);
	      }
	  }
      
	iii++;
      }
  

    //if(end>row_offset+main_indices[N_psi-2] && start<=row_offset+main_indices[N_psi-2])
    {
      if(!(main_indices[N_psi-2]==-1))
	{
	  //(N-2)th
	  if(!(sec_indices[N_psi-3]==-1))
	    {
	      value1=integrate1d(matrix_coeffs,gq_mod,rad_interp,shape1,N_psi,num_quad,N_psi-2,N_psi-3);
	      MatSetValue(M,main_indices[N_psi-2]+row_offset,sec_indices[N_psi-3]+col_offset,value1,insert_mode);
	    }

	  if(!(sec_indices[N_psi-2]==-1))
	    {
	      value2=integrate1d(matrix_coeffs,gq_mod,rad_interp,shape2,N_psi,num_quad,N_psi-2,N_psi-3);
	      value3=integrate1d(matrix_coeffs,gq_mod,rad_interp,shape3,N_psi,num_quad,N_psi-1,N_psi-2);
	      MatSetValue(M,main_indices[N_psi-2]+row_offset,sec_indices[N_psi-2]+col_offset,value2+value3,insert_mode);
	    }
  
	  if(!(sec_indices[N_psi-1]==-1))
	    {
	      value4=integrate1d(matrix_coeffs,gq_mod,rad_interp,shape4,N_psi,num_quad,N_psi-1,N_psi-2);
	      MatSetValue(M,main_indices[N_psi-2]+row_offset,sec_indices[N_psi-1]+col_offset,value4,insert_mode);
	    }
	}
    }
	  
    //if(end>row_offset+main_indices[N_psi-1] && start<=row_offset+main_indices[N_psi-1])
    {
      //(N-1)th i.e. final
      if(!(main_indices[N_psi-1]==-1))
	{
	  if(!(sec_indices[N_psi-2]==-1))
	    {
	      value1=integrate1d(matrix_coeffs,gq_mod,rad_interp,shape1,N_psi,num_quad,N_psi-1,N_psi-2);
	      MatSetValue(M,main_indices[N_psi-1]+row_offset,sec_indices[N_psi-2]+col_offset,value1,insert_mode);
	    }
      
	  if(!(sec_indices[N_psi-1]==-1))
	    {
	      value2=integrate1d(matrix_coeffs,gq_mod,rad_interp,shape2,N_psi,num_quad,N_psi-1,N_psi-2);
	      MatSetValue(M,main_indices[N_psi-1]+row_offset,sec_indices[N_psi-1]+col_offset,value2,insert_mode);
	    }
	}
    }
  }

}



std::complex<double> integrate1d(std::complex<double> func_grid[],double rad_interp[],double (*shape[])(double,double,double),int N_psi,int num_quad,int upper_point,int lower_point)
{
  assert( upper_point >= lower_point );

  if(shape[0] == N_NULL || shape[1] == N_NULL){ return 0.0; }

  //Guassian quadrature
  std::complex<double> result = 0.0;
  double upper = rad_interp[upper_point * (num_quad+1)] ;
  double lower = rad_interp[lower_point * (num_quad+1)] ;

  const double* gqNeval;
  const double* gqNweigh;

  assign_gq(gqNeval,gqNweigh,num_quad);

  for(int iii=0 ; iii < num_quad; iii++){
    result += gqNweigh[iii] * func_grid[lower_point * (num_quad+1) + iii + 1] * shape[0](rad_interp[lower_point * (num_quad+1) + iii + 1],upper,lower) * shape[1](rad_interp[lower_point * (num_quad+1) + iii + 1],upper,lower) ;
  }

  result *= 0.5 * ( upper - lower );

  return result;
}

std::complex<double> integrate1d(std::complex<double> func_grid[],gq_mod* gq_mod,double rad_interp[],void (*shape[])(std::vector<double>*,double,double),int N_psi,int num_quad,int upper_point,int lower_point)
{
  assert( upper_point >= lower_point );
  
  std::complex<double> result = 0.0;
  double upper = rad_interp[upper_point * (num_quad+1)] ;
  double lower = rad_interp[lower_point * (num_quad+1)] ;

  if(shape[0] == N_NULL_mod || shape[1] == N_NULL_mod){ return 0.0; } //Does this comparison always work? Is it based on addresses?

  std::vector<double> shape1,shape2;
  shape[0](&shape1,upper,lower);
  shape[1](&shape2,upper,lower);

  std::vector<double> final_pol = combine_pols(shape1,shape2) ;

  const double* gqNeval;
  const double* gqNweigh;

  assign_gq(gqNeval,gqNweigh,num_quad);

  double shape_out;
  if(gq_mod->wt_gq == 1){ //Just do regular gq integration
      
    for(int iii=0 ; iii < num_quad; iii++){
      
      shape_out = 0.0 ;
      for( int jjj = 0 ; jjj < final_pol.size() ; jjj++ ){
	shape_out += final_pol[jjj] * pow(rad_interp[lower_point * (num_quad+1) + iii + 1],jjj) ;
      }
      
      result += gqNweigh[iii] * func_grid[lower_point * (num_quad+1) + iii + 1] * shape_out ;
    }

    result *= 0.5 * ( upper - lower );

    return result;
  }
  else{
    std::complex<double> result_compare = 0.0 ;
    if(final_pol.size() == 1 || final_pol.size() == 2){ //Case of polynomial < O(2)
      for(int iii=0 ; iii < num_quad; iii++){
	//Calculate quadrature via gq_mod
	shape_out = 0.0 ;
	for( int jjj = 0 ; jjj < final_pol.size() ; jjj++ ){
	  shape_out += final_pol[jjj] * pow(rad_interp[lower_point * (num_quad+1) + iii + 1],jjj) ;
	}

	result += gq_mod->weights[lower_point][iii] * func_grid[lower_point * (num_quad+1) + iii + 1] * shape_out ;
	
	//Calc normal gq to compare
	result_compare += ( gqNweigh[iii] - gq_mod->weights[lower_point][iii] ) * func_grid[lower_point * (num_quad+1) + iii + 1] * shape_out ; 
      }

      result *= 0.5 * ( upper - lower );
      result_compare *= 0.5 * ( upper - lower );

      if( std::abs( result_compare ) < gq_mod->tol ){ gq_mod->wt_gq = 1 ;}

      return result;
    }
    else{
      for(int iii=0 ; iii < num_quad; iii++){
	//Calculate quadrature via gq_mod

	//Ln/cn terms
	shape_out = 0.0 ;
	for( int jjj = 0 ; jjj < 2 ; jjj++ ){
	  shape_out += final_pol[jjj] * pow(rad_interp[lower_point * (num_quad+1) + iii + 1],jjj) ;
	}
	
	//Constant/linear parts of quadrature done via gq_mod
	result += gq_mod->weights[lower_point][iii] * func_grid[lower_point * (num_quad+1) + iii + 1] * shape_out ;

	//Take diff with normal gq on ln/cn parts to compare
	result_compare += ( gqNweigh[iii] - gq_mod->weights[lower_point][iii] ) * func_grid[lower_point * (num_quad+1) + iii + 1] * shape_out ;

	//Higher order terms
	shape_out = 0.0 ;
	for( int jjj = 2 ; jjj < final_pol.size() ; jjj++ ){
	  shape_out += final_pol[jjj] * pow(rad_interp[lower_point * (num_quad+1) + iii + 1],jjj) ;
	}
	
	//Higher order terms done via regular gq
	result += gqNweigh[iii] * func_grid[lower_point * (num_quad+1) + iii + 1] * shape_out ;
      }

      result *= 0.5 * ( upper - lower );
      result_compare *= 0.5 * ( upper - lower );

      if( std::abs( result_compare ) < gq_mod->tol ){ gq_mod->wt_gq = 1 ; }

      return result; 
    }
  }
  
  return 0.0;
}

std::vector<double> combine_pols(std::vector<double> shape1,std::vector<double> shape2)
{
  int len1 = shape1.size() ;
  int len2 = shape2.size() ;
  std::vector<double> out_vec (len1 + len2 - 1,0.0) ;

  for(int iii = 0; iii < len1 ; iii++){
    for(int jjj = 0 ; jjj < len2 ; jjj++){
      out_vec[iii+jjj] += shape1[iii] * shape2[jjj] ;
    }}

  return out_vec ;
}


void select_shape(std::string* out_shape,std::string shape_order,int shape_num,bool is_perp)
{
  if(shape_order=="NHLC"){
    assert(shape_num==0);

    if(is_perp){*out_shape="NHLN";}
    else{*out_shape="NHCN";}
  }
  else if(shape_order=="HLC"){
    assert(shape_num==0);

    if(is_perp){*out_shape="HLN";}
    else{*out_shape="NHCN";}
  }
  else if(shape_order=="LN"){
    assert(shape_num==0);

    if(is_perp){*out_shape="NHLN";}
    else{*out_shape="NHLN";}
  }
  else if(shape_order=="HLN"){
    assert(shape_num==0);

    *out_shape="HLN";
  }
  else if(shape_order=="NHQL"){
    if(is_perp){assert(shape_num<=1);}
    else{assert(shape_num==0);}

    if(is_perp){
      if(shape_num==0){*out_shape="NHQD_1";}
      else{*out_shape="NHQD_2";}}
    else{*out_shape="NHLN";}
  }
  else if(shape_order=="HQL"){
    if(is_perp){assert(shape_num<=1);}
    else{assert(shape_num==0);}

    if(is_perp){
      if(shape_num==0){*out_shape="HQD_1";}
      else{*out_shape="HQD_2";}}
    else{*out_shape="NHLN";}
  }
  else if(shape_order=="HQD"){
    assert(shape_num<=1);

    if(shape_num==0){*out_shape="HQD_1";}
    else{*out_shape="HQD_2";}
  }
  else if(shape_order=="NHCQ"){
    if(is_perp){assert(shape_num<=1);}
    else{assert(shape_num<=1);}

    if(is_perp){
      if(shape_num==0){*out_shape="NHCB_1";}
      else{*out_shape="NHCB_2";}}
    else{
      if(shape_num==0){*out_shape="NHQD_1";}
      else{*out_shape="NHQD_2";}}
  }
  else if(shape_order=="CB"){
    if(is_perp){assert(shape_num<=1);}
    else{assert(shape_num<=1);}

    if(shape_num==0){*out_shape="NHCB_1";}
    else{*out_shape="NHCB_2";}
  }
  else if(shape_order=="NHQC"){
    if(is_perp){assert(shape_num<=2);}
    else{assert(shape_num<=1);}

    if(is_perp){
      if(shape_num==0){*out_shape="NHQR_1";}
      else if(shape_num==1){*out_shape="NHQR_2";}
      else{*out_shape="NHQR_3";}}
    else{
      if(shape_num==0){*out_shape="NHCB_1";}
      else{*out_shape="NHCB_2";}}
  }
  else{std::cout << "That shape order is not recognised (select_shape)" << std::endl;}
}

void fill_shapes(double (*shape[])(double,double,double),bool deriv,std::string shape_order)
{
  /****************************************************************************************************************************************************************************/
  //LINEAR

  if(shape_order=="NHCN"){
    if(deriv){
      shape[0]=N_NULL;
      shape[1]=N_NULL;
    }
    else{
      shape[0]=N_NULL;
      shape[1]=N_A_nhcn;
    }  
  }
  else if(shape_order=="NHLN"){
    if(deriv){
      shape[0]=dN_A_nhln;
      shape[1]=dN_B_nhln;
    }
    else{
      shape[0]=N_A_nhln;
      shape[1]=N_B_nhln;
    } 
  }
  else if(shape_order=="HLN"){
    if(deriv){
      shape[0]=dN_A_hln;
      shape[1]=dN_B_hln;
    }
    else{
      shape[0]=N_A_hln;
      shape[1]=N_B_hln;
    }
  }
  else if(shape_order=="NHQD_1"){
    if(deriv){
      shape[0]=dN_A_nhqd;
      shape[1]=dN_B_nhqd;
    }
    else{
      shape[0]=N_A_nhqd;
      shape[1]=N_B_nhqd;
    } 
  }
  else if(shape_order=="NHQD_2"){
    if(deriv){
      shape[0]=N_NULL;
      shape[1]=dN_C_nhqd;
    }
    else{
      shape[0]=N_NULL;
      shape[1]=N_C_nhqd;
    } 
  }
  else if(shape_order=="HQD_1"){
    if(deriv){
      shape[0]=dN_A_hqd;
      shape[1]=dN_B_hqd;
    }
    else{
      shape[0]=N_A_hqd;
      shape[1]=N_B_hqd;
    } 
  }
  else if(shape_order=="HQD_2"){
    if(deriv){
      shape[0]=N_NULL;
      shape[1]=dN_C_hqd;
    }
    else{
      shape[0]=N_NULL;
      shape[1]=N_C_hqd;
    } 
  }
  else if(shape_order=="NHCB_1"){
    if(deriv){
      shape[0]=dN_A_nhcb;
      shape[1]=dN_B_nhcb;
    }
    else{
      shape[0]=N_A_nhcb;
      shape[1]=N_B_nhcb;
    } 
  }
  else if(shape_order=="NHCB_2"){
    if(deriv){
      shape[0]=dN_C_nhcb;
      shape[1]=dN_D_nhcb;
    }
    else{
      shape[0]=N_C_nhcb;
      shape[1]=N_D_nhcb;
    } 
  }
  else if(shape_order=="NHQR_1"){
    if(deriv){
      shape[0]=dN_A_nhqr;
      shape[1]=dN_B_nhqr;
    }
    else{
      shape[0]=N_A_nhqr;
      shape[1]=N_B_nhqr;
    } 
  }
  else if(shape_order=="NHQR_2"){
    if(deriv){
      shape[0]=N_NULL;
      shape[1]=dN_C_nhqr;
    }
    else{
      shape[0]=N_NULL;
      shape[1]=N_C_nhqr;
    } 
  }
  else if(shape_order=="NHQR_3"){
    if(deriv){
      shape[0]=dN_D_nhqr;
      shape[1]=dN_E_nhqr;
    }
    else{
      shape[0]=N_D_nhqr;
      shape[1]=N_E_nhqr;
    } 
  }
  else{std::cout << "That shape order is not recognised (fill_shapes)" << std::endl;}

}


void fill_shapes(void (*shape[])(std::vector<double>*,double,double),bool deriv,std::string shape_order)
{
  /****************************************************************************************************************************************************************************/
  //LINEAR

  if(shape_order=="NHCN"){
    if(deriv){
      shape[0]=N_NULL_mod;
      shape[1]=N_NULL_mod;
    }
    else{
      shape[0]=N_NULL_mod;
      shape[1]=N_A_nhcn_mod;
    }  
  }
  else if(shape_order=="NHLN"){
    if(deriv){
      shape[0]=dN_A_nhln_mod;
      shape[1]=dN_B_nhln_mod;
    }
    else{
      shape[0]=N_A_nhln_mod;
      shape[1]=N_B_nhln_mod;
    } 
  }
  else if(shape_order=="HLN"){
    if(deriv){
      shape[0]=dN_A_hln_mod;
      shape[1]=dN_B_hln_mod;
    }
    else{
      shape[0]=N_A_hln_mod;
      shape[1]=N_B_hln_mod;
    }
  }
  else if(shape_order=="NHQD_1"){
    if(deriv){
      shape[0]=dN_A_nhqd_mod;
      shape[1]=dN_B_nhqd_mod;
    }
    else{
      shape[0]=N_A_nhqd_mod;
      shape[1]=N_B_nhqd_mod;
    } 
  }
  else if(shape_order=="NHQD_2"){
    if(deriv){
      shape[0]=N_NULL_mod;
      shape[1]=dN_C_nhqd_mod;
    }
    else{
      shape[0]=N_NULL_mod;
      shape[1]=N_C_nhqd_mod;
    } 
  }
  else if(shape_order=="HQD_1"){
    if(deriv){
      shape[0]=dN_A_hqd_mod;
      shape[1]=dN_B_hqd_mod;
    }
    else{
      shape[0]=N_A_hqd_mod;
      shape[1]=N_B_hqd_mod;
    } 
  }
  else if(shape_order=="HQD_2"){
    if(deriv){
      shape[0]=N_NULL_mod;
      shape[1]=dN_C_hqd_mod;
    }
    else{
      shape[0]=N_NULL_mod;
      shape[1]=N_C_hqd_mod;
    } 
  }
  else if(shape_order=="NHCB_1"){
    if(deriv){
      shape[0]=dN_A_nhcb_mod;
      shape[1]=dN_B_nhcb_mod;
    }
    else{
      shape[0]=N_A_nhcb_mod;
      shape[1]=N_B_nhcb_mod;
    } 
  }
  else if(shape_order=="NHCB_2"){
    if(deriv){
      shape[0]=dN_C_nhcb_mod;
      shape[1]=dN_D_nhcb_mod;
    }
    else{
      shape[0]=N_C_nhcb_mod;
      shape[1]=N_D_nhcb_mod;
    } 
  }
  else if(shape_order=="NHQR_1"){
    if(deriv){
      shape[0]=dN_A_nhqr_mod;
      shape[1]=dN_B_nhqr_mod;
    }
    else{
      shape[0]=N_A_nhqr_mod;
      shape[1]=N_B_nhqr_mod;
    } 
  }
  else if(shape_order=="NHQR_2"){
    if(deriv){
      shape[0]=N_NULL_mod;
      shape[1]=dN_C_nhqr_mod;
    }
    else{
      shape[0]=N_NULL_mod;
      shape[1]=N_C_nhqr_mod;
    } 
  }
  else if(shape_order=="NHQR_3"){
    if(deriv){
      shape[0]=dN_D_nhqr_mod;
      shape[1]=dN_E_nhqr_mod;
    }
    else{
      shape[0]=N_D_nhqr_mod;
      shape[1]=N_E_nhqr_mod;
    } 
  }
  else{std::cout << "That shape order is not recognised (fill_shapes)" << std::endl;}

}



void calc_fourier_sym(double grid_func[],std::complex<double> fourier_func[],int N_psi,int N_theta)
{
  fftw_complex *out;
  fftw_plan plan;

  double *in = new double[N_theta];
  int size=(N_theta/2)+1; if(!(N_theta%2==0) /*If N_theta is odd*/){size=((N_theta-1)/2)+1;}
  int fourier_size=(N_theta/2); if(!(N_theta%2==0) /*If N_theta is odd*/){fourier_size=((N_theta-1)/2)+1;}

  out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * size);
  for(int iii=0;iii<size;iii++){out[iii][0]=0.0; out[iii][1]=0.0;}
  
  plan = fftw_plan_dft_r2c_1d(N_theta, in, out, FFTW_ESTIMATE); //ESTIMATE , MEASURE , PATIENT , EXHAUSTIVE - only changes how quickly the plan runs, not how accurate the results are (results don't change)
  
  for(int iii=0;iii<size;iii++){out[iii][0]=0.0; out[iii][1]=0.0;}
  
  for(int jjj=0;jjj<N_psi;jjj++)
    {

      for(int iii=0;iii<N_theta;iii++){in[iii]=grid_func[jjj*N_theta+iii];}

      fftw_execute(plan);

      

      for(int iii=0;iii<fourier_size;iii++){fourier_func[jjj*fourier_size+iii]=(out[iii][0]+imag_unit*out[iii][1])/static_cast<double>(N_theta);}
 
    } 

  delete[] in;
  fftw_destroy_plan(plan);
  fftw_free(out);
}

void calc_fourier_full(std::complex<double> grid_func[],std::complex<double> fourier_func[],int N_psi,int N_theta)
{ 
  fftw_complex *out;
  fftw_plan plan;

  fftw_complex *in = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N_theta);

  out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * N_theta);
  for(int iii=0;iii<N_theta;iii++){out[iii][0]=0.0; out[iii][1]=0.0;}
  
  plan = fftw_plan_dft_1d(N_theta,in,out,FFTW_FORWARD,FFTW_ESTIMATE); //ESTIMATE , MEASURE , PATIENT , EXHAUSTIVE - only changes how quickly the plan runs, not how accurate the results are (results don't change)

  for(int iii=0;iii<N_theta;iii++){out[iii][0]=0.0; out[iii][1]=0.0;}
  
  for(int jjj=0;jjj<N_psi;jjj++)
    {

      for(int iii=0;iii<N_theta;iii++){in[iii][0]=std::real(grid_func[jjj*N_theta+iii]); in[iii][1]=std::imag(grid_func[jjj*N_theta+iii]);}

      fftw_execute(plan);

      //Since the N_theta/2 and -N_theta/2 term are ambiguous, just go to |N_theta/2 - 1|. Also order from -(N_theta/2 - 1) to (N_theta/2 - 1)
      if(N_theta%2==0)
	{
	  for(int iii=0;iii<(N_theta/2)-1;iii++){fourier_func[jjj*(N_theta-1)+iii]=(out[(N_theta/2)+1+iii][0]+imag_unit*out[(N_theta/2)+1+iii][1])/static_cast<double>(N_theta);} //Negative modes
	  for(int iii=(N_theta/2)-1;iii<N_theta-1;iii++){fourier_func[jjj*(N_theta-1)+iii]=(out[iii-(N_theta/2)+1][0]+imag_unit*out[iii-(N_theta/2)+1][1])/static_cast<double>(N_theta);} //Non-negative modes
	}
      else
	{
	  for(int iii=0;iii<(N_theta-1)/2;iii++){fourier_func[jjj*(N_theta)+iii]=(out[((N_theta+1)/2)+iii][0]+imag_unit*out[((N_theta+1)/2)+iii][1])/static_cast<double>(N_theta);} //Negative modes
	  for(int iii=(N_theta-1)/2;iii<N_theta;iii++){fourier_func[jjj*(N_theta)+iii]=(out[iii-((N_theta-1)/2)][0]+imag_unit*out[iii-((N_theta-1)/2)][1])/static_cast<double>(N_theta);} //Non-negative modes
	}
    }
  

  fftw_destroy_plan(plan);
  fftw_free(out);
  fftw_free(in);
}


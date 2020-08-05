#include "constants.h"

void init_equil(equil_fields *equil,int N_psi,int N_interp,int N_theta)
{
  equil->psi_grid=new double[N_psi];
  equil->psi_interp=new double[N_interp];
  equil->theta_grid=new double[N_theta];
  equil->rad_var=new double[N_psi];
  equil->rad_interp=new double[N_interp];
  equil->drad_dpsi=new double[N_interp];

  equil->maj_rad=new double[N_interp*N_theta];
  equil->height=new double[N_interp*N_theta];
  
  equil->dens=new double[N_interp*N_theta];
  equil->pres=new double[N_interp*N_theta];
  equil->f_psi=new double[N_interp*N_theta];

  //////////////////////////////////////// Dubious necessity
  equil->dpsi_dR=new double[N_interp*N_theta];
  equil->dpsi_dZ=new double[N_interp*N_theta];
  equil->dth_dR=new double[N_interp*N_theta];
  equil->dth_dZ=new double[N_interp*N_theta];
  //////////////////////////////////////////
  
  equil->g_pp=new double[N_interp*N_theta];
  equil->g_pt=new double[N_interp*N_theta];
  equil->g_tt=new double[N_interp*N_theta];
  equil->g_phph=new double[N_interp*N_theta];	
  equil->jacob=new double[N_interp*N_theta];
  equil->mag_sq=new double[N_interp*N_theta];
  
  equil->curv_psi=new double[N_interp*N_theta];
  equil->neg_shear=new double[N_interp*N_theta];
  equil->j_dot_b=new double[N_interp*N_theta];
  equil->ones=new double[N_interp*N_theta];

  equil->cc_0=new double[N_interp*N_theta];
  equil->cc_1=new double[N_interp*N_theta];
  equil->cc_2=new double[N_interp*N_theta];
  equil->cc_3=new double[N_interp*N_theta];
  equil->cc_4=new double[N_interp*N_theta];
}

void init_coeffs(matrix_coeffs* coeffs,geom_shape geom)
{
  int N_interp=geom.N_interp;
  int N_theta=geom.N_theta;
  
  coeffs->f_dpdp=new std::complex<double>[N_interp*N_theta]; clean_grid(coeffs->f_dpdp,N_interp*N_theta);
  coeffs->f_dpp=new std::complex<double>[N_interp*N_theta]; clean_grid(coeffs->f_dpp,N_interp*N_theta);
  coeffs->f_pdp=new std::complex<double>[N_interp*N_theta]; clean_grid(coeffs->f_pdp,N_interp*N_theta);
  coeffs->f_pp=new std::complex<double>[N_interp*N_theta]; clean_grid(coeffs->f_pp,N_interp*N_theta);
  coeffs->f_dpw=new std::complex<double>[N_interp*N_theta]; clean_grid(coeffs->f_dpw,N_interp*N_theta);
  coeffs->f_pw=new std::complex<double>[N_interp*N_theta]; clean_grid(coeffs->f_pw,N_interp*N_theta);
  coeffs->f_dwdp=new std::complex<double>[N_interp*N_theta]; clean_grid(coeffs->f_dwdp,N_interp*N_theta);
  coeffs->f_dwp=new std::complex<double>[N_interp*N_theta]; clean_grid(coeffs->f_dwp,N_interp*N_theta);
  coeffs->f_wdp=new std::complex<double>[N_interp*N_theta]; clean_grid(coeffs->f_wdp,N_interp*N_theta);
  coeffs->f_wp=new std::complex<double>[N_interp*N_theta]; clean_grid(coeffs->f_wp,N_interp*N_theta);
  coeffs->f_dww=new std::complex<double>[N_interp*N_theta]; clean_grid(coeffs->f_dww,N_interp*N_theta);
  coeffs->f_ww=new std::complex<double>[N_interp*N_theta]; clean_grid(coeffs->f_ww,N_interp*N_theta);
}

void init_gq_mod(gq_mod* gq_mod, geom_shape geom, equil_fields eq)
{
  gq_mod->weights = new double*[geom.N_psi-1] ;
  for(int i=0 ; i < geom.N_psi - 1 ; i++){
    gq_mod->weights[i] = new double[geom.num_quad] ;
  }

  std::string value="";
  read_in("num_div",value);
  if( !(value == "") ) { gq_mod->num_div = std::stoi( value ) ; }
  else{ gq_mod->num_div = 3 ; }

  gq_mod->wt_gq = 0 ;
  gq_mod->tol = 1.0e-10 ;

  //Since comparing only linear/constant terms need to normalise the tolerance 
  if(geom.shape_order=="NHLC"){ gq_mod->tol /= ( eq.rad_var[1] - eq.rad_var[0] ) ; }
  else if(geom.shape_order=="NHQL"){ gq_mod->tol /= pow(( eq.rad_var[1] - eq.rad_var[0] ),3) ; }
  else if(geom.shape_order=="NHCQ"){ gq_mod->tol /= pow(( eq.rad_var[1] - eq.rad_var[0] ),5) ; } 
}

void delete_gq_mod(gq_mod* gq_mod, int N_psi)
{
  for(int i=0 ; i < N_psi - 1 ; i++){
    delete[] gq_mod->weights[i]; gq_mod->weights[i] = NULL ;
  }

  delete[] gq_mod->weights; gq_mod->weights = NULL ;
}

void delete_full_geom_shape(geom_shape* geom)
{
  delete[] geom->pol_pos; geom->pol_pos=NULL;
  
  for(int i=0;i<geom->num_var_perp;i++){
    delete[] geom->ind_perp_main[i]; geom->ind_perp_main[i]=NULL;
    delete[] geom->ind_perp_sec[i]; geom->ind_perp_sec[i]=NULL;
  }
  for(int i=0;i<geom->num_var_wedge;i++){
    delete[] geom->ind_wedge_main[i]; geom->ind_wedge_main[i]=NULL;
    delete[] geom->ind_wedge_sec[i]; geom->ind_wedge_sec[i]=NULL;
  }
  
  delete[] geom->ind_perp_main; geom->ind_perp_main=NULL;
  delete[] geom->ind_perp_sec; geom->ind_perp_sec=NULL;
  delete[] geom->ind_wedge_main; geom->ind_wedge_main=NULL;
  delete[] geom->ind_wedge_sec; geom->ind_wedge_sec=NULL;
}

void delete_full_equil_fields(equil_fields* eq)
{
  delete[] eq->psi_grid; eq->psi_grid=NULL;
  delete[] eq->theta_grid; eq->theta_grid=NULL;
  delete[] eq->rad_var; eq->rad_var=NULL;
  delete[] eq->rad_interp; eq->rad_interp=NULL;
  delete[] eq->psi_interp; eq->psi_interp=NULL;
  delete[] eq->drad_dpsi; eq->drad_dpsi=NULL;
  delete[] eq->maj_rad; eq->maj_rad=NULL;
  delete[] eq->height; eq->height=NULL;
  
  delete[] eq->f_psi; eq->f_psi=NULL;
  delete[] eq->dens; eq->dens=NULL;
  delete[] eq->pres; eq->pres=NULL;
  
  delete[] eq->dpsi_dR; eq->dpsi_dR=NULL;
  delete[] eq->dpsi_dZ; eq->dpsi_dZ=NULL;
  delete[] eq->dth_dR; eq->dth_dR=NULL;
  delete[] eq->dth_dZ; eq->dth_dZ=NULL;
  
  delete[] eq->g_pp; eq->g_pp=NULL;
  delete[] eq->g_pt; eq->g_pt=NULL;
  delete[] eq->g_tt; eq->g_tt=NULL;
  delete[] eq->g_phph; eq->g_phph=NULL;
  delete[] eq->jacob; eq->jacob=NULL;
  delete[] eq->mag_sq; eq->mag_sq=NULL;
  delete[] eq->ones; eq->ones=NULL;
  delete[] eq->j_dot_b; eq->j_dot_b=NULL;
  delete[] eq->curv_psi; eq->curv_psi=NULL;
  delete[] eq->neg_shear; eq->neg_shear=NULL;
  
  delete[] eq->cc_0; eq->cc_0=NULL;
  delete[] eq->cc_1; eq->cc_1=NULL;
  delete[] eq->cc_2; eq->cc_2=NULL;
  delete[] eq->cc_3; eq->cc_3=NULL;
  delete[] eq->cc_4; eq->cc_4=NULL;
}

void delete_full_matrix_coeffs(matrix_coeffs* coeffs)
{
  delete[] coeffs->f_dpdp; coeffs->f_dpdp=NULL;
  delete[] coeffs->f_dpp; coeffs->f_dpp=NULL;
  delete[] coeffs->f_pdp; coeffs->f_pdp=NULL;
  delete[] coeffs->f_pp; coeffs->f_pp=NULL;
  delete[] coeffs->f_dpw; coeffs->f_dpw=NULL;
  delete[] coeffs->f_pw; coeffs->f_pw=NULL;
  delete[] coeffs->f_dwdp; coeffs->f_dwdp=NULL;
  delete[] coeffs->f_dwp; coeffs->f_dwp=NULL;
  delete[] coeffs->f_wdp; coeffs->f_wdp=NULL;
  delete[] coeffs->f_wp; coeffs->f_wp=NULL;
  delete[] coeffs->f_dww; coeffs->f_dww=NULL;
  delete[] coeffs->f_ww; coeffs->f_ww=NULL;
}

void assign_gq( const double*& gq_eval , const double*& gq_weigh , int num_quad){
  
  if(num_quad == 4){ gq_eval = gq4eval ; gq_weigh = gq4weigh ; }
  else if(num_quad == 6){ gq_eval = gq6eval ; gq_weigh = gq6weigh ; }
  else if(num_quad == 12){ gq_eval = gq12eval ; gq_weigh = gq12weigh ; }
  else if(num_quad == 18){ gq_eval = gq18eval ; gq_weigh = gq18weigh ; }
}



void udsym(equil_fields* eq, geom_shape geom)
{
  int half_theta;
  if(geom.N_theta%2==0){half_theta=(geom.N_theta-2)/2;}
  else{half_theta=(geom.N_theta-1)/2;}
  int N_theta=geom.N_theta;

  //Set up-down values to their average
  for(int iii=0;iii<geom.N_interp;iii++){
    for(int jjj=1;jjj<=half_theta;jjj++){
      eq->maj_rad[iii*N_theta+jjj]=0.5*(eq->maj_rad[iii*N_theta+jjj]+eq->maj_rad[(iii+1)*N_theta-jjj]); eq->maj_rad[(iii+1)*N_theta-jjj]=eq->maj_rad[iii*N_theta+jjj];
      eq->height[iii*N_theta+jjj]=0.5*(eq->height[iii*N_theta+jjj]-eq->height[(iii+1)*N_theta-jjj]); eq->height[(iii+1)*N_theta-jjj]=-eq->height[iii*N_theta+jjj]; //Note height should be up-down antisymmetric
      eq->f_psi[iii*N_theta+jjj]=0.5*(eq->f_psi[iii*N_theta+jjj]+eq->f_psi[(iii+1)*N_theta-jjj]); eq->f_psi[(iii+1)*N_theta-jjj]=eq->f_psi[iii*N_theta+jjj];
      eq->dens[iii*N_theta+jjj]=0.5*(eq->dens[iii*N_theta+jjj]+eq->dens[(iii+1)*N_theta-jjj]); eq->dens[(iii+1)*N_theta-jjj]=eq->dens[iii*N_theta+jjj];
      eq->pres[iii*N_theta+jjj]=0.5*(eq->pres[iii*N_theta+jjj]+eq->pres[(iii+1)*N_theta-jjj]); eq->pres[(iii+1)*N_theta-jjj]=eq->pres[iii*N_theta+jjj];
      eq->g_pp[iii*N_theta+jjj]=0.5*(eq->g_pp[iii*N_theta+jjj]+eq->g_pp[(iii+1)*N_theta-jjj]); eq->g_pp[(iii+1)*N_theta-jjj]=eq->g_pp[iii*N_theta+jjj];
      eq->g_pt[iii*N_theta+jjj]=0.5*(eq->g_pt[iii*N_theta+jjj]+eq->g_pt[(iii+1)*N_theta-jjj]); eq->g_pt[(iii+1)*N_theta-jjj]=eq->g_pt[iii*N_theta+jjj];
      eq->g_tt[iii*N_theta+jjj]=0.5*(eq->g_tt[iii*N_theta+jjj]+eq->g_tt[(iii+1)*N_theta-jjj]); eq->g_tt[(iii+1)*N_theta-jjj]=eq->g_tt[iii*N_theta+jjj];
      eq->g_phph[iii*N_theta+jjj]=0.5*(eq->g_phph[iii*N_theta+jjj]+eq->g_phph[(iii+1)*N_theta-jjj]); eq->g_phph[(iii+1)*N_theta-jjj]=eq->g_phph[iii*N_theta+jjj];
      eq->jacob[iii*N_theta+jjj]=0.5*(eq->jacob[iii*N_theta+jjj]+eq->jacob[(iii+1)*N_theta-jjj]); eq->jacob[(iii+1)*N_theta-jjj]=eq->jacob[iii*N_theta+jjj];
      eq->mag_sq[iii*N_theta+jjj]=0.5*(eq->mag_sq[iii*N_theta+jjj]+eq->mag_sq[(iii+1)*N_theta-jjj]); eq->mag_sq[(iii+1)*N_theta-jjj]=eq->mag_sq[iii*N_theta+jjj];
      eq->j_dot_b[iii*N_theta+jjj]=0.5*(eq->j_dot_b[iii*N_theta+jjj]+eq->j_dot_b[(iii+1)*N_theta-jjj]); eq->j_dot_b[(iii+1)*N_theta-jjj]=eq->j_dot_b[iii*N_theta+jjj];
      eq->curv_psi[iii*N_theta+jjj]=0.5*(eq->curv_psi[iii*N_theta+jjj]+eq->curv_psi[(iii+1)*N_theta-jjj]); eq->curv_psi[(iii+1)*N_theta-jjj]=eq->curv_psi[iii*N_theta+jjj];
      eq->neg_shear[iii*N_theta+jjj]=0.5*(eq->neg_shear[iii*N_theta+jjj]+eq->neg_shear[(iii+1)*N_theta-jjj]); eq->neg_shear[(iii+1)*N_theta-jjj]=eq->neg_shear[iii*N_theta+jjj];

      eq->cc_0[iii*N_theta+jjj]=0.5*(eq->cc_0[iii*N_theta+jjj]+eq->cc_0[(iii+1)*N_theta-jjj]); eq->cc_0[(iii+1)*N_theta-jjj]=eq->cc_0[iii*N_theta+jjj];
      eq->cc_1[iii*N_theta+jjj]=0.5*(eq->cc_1[iii*N_theta+jjj]+eq->cc_1[(iii+1)*N_theta-jjj]); eq->cc_1[(iii+1)*N_theta-jjj]=eq->cc_1[iii*N_theta+jjj];
      eq->cc_2[iii*N_theta+jjj]=0.5*(eq->cc_2[iii*N_theta+jjj]+eq->cc_2[(iii+1)*N_theta-jjj]); eq->cc_2[(iii+1)*N_theta-jjj]=eq->cc_2[iii*N_theta+jjj];
      eq->cc_3[iii*N_theta+jjj]=0.5*(eq->cc_3[iii*N_theta+jjj]+eq->cc_3[(iii+1)*N_theta-jjj]); eq->cc_3[(iii+1)*N_theta-jjj]=eq->cc_3[iii*N_theta+jjj];
      eq->cc_4[iii*N_theta+jjj]=0.5*(eq->cc_4[iii*N_theta+jjj]+eq->cc_4[(iii+1)*N_theta-jjj]); eq->cc_4[(iii+1)*N_theta-jjj]=eq->cc_4[iii*N_theta+jjj];
    }
  }
}


int calc_hermitian_diff(Mat M,int row_dim,int col_dim)
{
  Mat M_tran,M_fin;
  Vec M_max;
  PetscErrorCode ierr;

  double *workspace = new double[row_dim*col_dim];
  PetscScalar *vals = new PetscScalar[row_dim*col_dim]; clean_grid(vals,row_dim*col_dim);
  PetscInt *row_ind,*col_ind; row_ind = new PetscInt[row_dim]; col_ind = new PetscInt[col_dim];

  for(int iii=0;iii<row_dim;iii++){row_ind[iii]=iii;}
  for(int iii=0;iii<col_dim;iii++){col_ind[iii]=iii;}
  MatGetValues(M,row_dim,row_ind,col_dim,col_ind,vals);


  for(int iii=0;iii<row_dim;iii++){
    for(int jjj=0;jjj<col_dim;jjj++){
      workspace[iii*col_dim+jjj]=std::abs(vals[iii*col_dim+jjj]-std::conj(vals[jjj*col_dim+iii]))/(std::abs(vals[iii*col_dim+jjj])+std::abs(vals[jjj*col_dim+iii]));
    }}

  double max=0.0;
  int row_count=0,col_count=0;
  for(int iii=0;iii<row_dim;iii++){
    for(int jjj=0;jjj<col_dim;jjj++){
      if(workspace[iii*col_dim+jjj]>max){max=workspace[iii*col_dim+jjj];row_count=iii;col_count=jjj;}
    }}

  std::cout << "Max herm difference is (normalised): " << max << " on (row,column): (" << row_count << "," << col_count << ")" << std::endl;

  MatInfo info;
  MatGetInfo(M,MAT_GLOBAL_MAX,&info);

  std::cout << "There are : "  << info.nz_allocated << " non-zeros " << std::endl;
  std::cout << "There are : "  << info.nz_used << " non-zeros " << std::endl;
  std::cout << "There are : "  << info.nz_unneeded << " non-zeros " << std::endl;
  
  delete[] workspace; delete[] vals; delete[] row_ind; delete[] col_ind;

  return 0;
}

void normalise(PetscScalar perp_values[],PetscScalar wedge_values[],int N_r,int m_range,int num_eig,bool rotate)
{
  if(rotate==true)
    {
      PetscScalar max_val_perp,max_val_wedge;
      for(int jjj=0;jjj<num_eig;jjj++)
	{
	  max_val_perp=0.0; max_val_wedge=0.0;
      
	  for(int iii=0;iii<m_range*N_r;iii++)
	    {
	      if(std::abs(perp_values[jjj*m_range*N_r+iii])>std::abs(max_val_perp)){max_val_perp=perp_values[jjj*m_range*N_r+iii];};
	      if(std::abs(wedge_values[jjj*m_range*N_r+iii])>std::abs(max_val_wedge)){max_val_wedge=wedge_values[jjj*m_range*N_r+iii];};
	    }

	  if(std::abs(max_val_perp)<std::abs(max_val_wedge) && !(std::abs(max_val_perp)==0.0)){max_val_perp*=std::abs(max_val_wedge)/std::abs(max_val_perp);}
	  else if(std::abs(max_val_perp)==0.0){max_val_perp=max_val_wedge;}
	    
	  for(int iii=0;iii<m_range*N_r;iii++)
	    {
	      perp_values[jjj*m_range*N_r+iii]/=max_val_perp;
	      wedge_values[jjj*m_range*N_r+iii]/=max_val_perp;
	    }
	}
    }
  else if(rotate==false)
    {
      double max_val;
      for(int jjj=0;jjj<num_eig;jjj++)
	{
	  max_val=0.0;
      
	  for(int iii=0;iii<m_range*N_r;iii++)
	    {
	      if(std::abs(perp_values[jjj*m_range*N_r+iii])>max_val){max_val=std::abs(perp_values[jjj*m_range*N_r+iii]);};
	      if(std::abs(wedge_values[jjj*m_range*N_r+iii])>max_val){max_val=std::abs(wedge_values[jjj*m_range*N_r+iii]);};
	    }

	  for(int iii=0;iii<m_range*N_r;iii++)
	    {
	      perp_values[jjj*m_range*N_r+iii]/=max_val;
	      wedge_values[jjj*m_range*N_r+iii]/=max_val;
	    }
	}
    }
}

void sort_eigenvecs(Vec eigenvecs[],PetscScalar xi_perp[],PetscScalar xi_wedge[],geom_shape geom,int num_vecs,int keep_indices[])
{
 const PetscScalar *tmp_arr;


  for(int iii=0;iii<num_vecs;iii++)
    { 
      VecGetArrayRead(eigenvecs[keep_indices[iii]],&tmp_arr);

      for(int jjj=0;jjj<geom.m_range;jjj++)
	{
	  fill_indices(geom.ind_perp_main,geom.ind_wedge_main,jjj+geom.m_min,geom.N_psi,geom.shape_order);

	  for(int kkk=0;kkk<geom.N_psi;kkk++)
	    {
	      if(geom.ind_perp_main[0][kkk]==-1){xi_perp[iii*geom.m_range*geom.N_psi+jjj*geom.N_psi+kkk]=0.0;}
	      else{xi_perp[iii*geom.m_range*geom.N_psi+jjj*geom.N_psi+kkk]=tmp_arr[geom.pol_pos[jjj]+geom.ind_perp_main[0][kkk]];}

	      if(geom.ind_wedge_main[0][kkk]==-1){xi_wedge[iii*geom.m_range*geom.N_psi+jjj*geom.N_psi+kkk]=0.0;}
	      else{xi_wedge[iii*geom.m_range*geom.N_psi+jjj*geom.N_psi+kkk]=tmp_arr[geom.pol_pos[jjj]+geom.ind_wedge_main[0][kkk]];}
	    }


	}

      VecRestoreArrayRead(eigenvecs[keep_indices[iii]],&tmp_arr);
    }
}



double newt_raph(double (*func)(double),double guess,int N_max,double tol)
{
  int counter=0;
  double deriv_val=0.0;
  double param_diff=1.2e-16;

  while(counter<N_max)
    {
      if( std::abs(func(guess))<tol ){ return guess; }

      deriv_val=(func(guess+param_diff)-func(guess-param_diff))/(2.0*param_diff);

      if(deriv_val==0){std::cout << "Newt-Raph reached turning point" << std::endl; return guess;}

      guess-=func(guess) / deriv_val ;
      
      counter++;
    }

  std::cout << "Newt-Raph failed to converge within tol in max iterations" << std::endl;
  return guess;
}

double expand_bracket(double (*func)(double),double a,double b,double brac_size,double ratio,bool a_fixed)
{
  assert(a<=b);
  assert(ratio != 1.0);

  if(ratio>1.0 && a_fixed==true || ratio<1.0 && a_fixed==false)
    {
      while((b-a)<brac_size)
	{
	  if( func(a)*func(b)<0.0 ){
	    if( a_fixed ){ return b; }
	    else{ return a; }
	  }

	  if( a_fixed ){ b*=ratio; }
	  else{ a*=ratio; }
	}
    }
  else
    {
      while((b-a)>brac_size)
	{
	  if( func(a)*func(b)<0.0 ){
	    if( a_fixed ){ return b; }
	    else{ return a; }
	  }

	  if( a_fixed ){ b*=ratio; }
	  else{ a*=ratio; }
	}
    }

  std::cout << "Max/min bracket size reached" << std::endl;
  if( a_fixed ){ return b; }
  else{ return a; }
}


double bisect(double (*func)(double),double a,double b,double target,double tol,int N_max)
{
  assert(a<b);

  //if(!((func(a)-target)*(func(b)-target)<0.0))
  //{std::cout << func(a) << "  " << a << "  " << func(b) << "  " << b << "  " << target << std::endl;}
  
  if(func(a)-target == 0.0){return a;}
  else if(func(b)-target == 0.0){return b;}
  //else{assert((func(a)-target)*(func(b)-target)<0.0);}


  int counter=0;
  double mid;
  while(counter<N_max)
    {
      mid=(a+b)*0.5;
      if( (b-a)*0.5 < tol || func(mid)-target==0.0 ){return mid;}

      if( (func(a)-target)*(func(mid)-target) > 0.0 ){a=mid;}
      else{b=mid;}

      counter++;
    }

  std::cout << "Max number of iterations exceeded" << std::endl;
  return mid;
}

void trim(std::string & s)
{
  size_t p = s.find_first_not_of(" \t");
  s.erase(0, p);

  p = s.find_last_not_of(" \t");
  if (std::string::npos != p)
    s.erase(p+1);
}

void read_in(std::string token_name,std::string & value)
{
  std::string line, token;
  std::string file_name="./input.txt"; //path is from where program is run

  std::ifstream myfile(file_name.c_str());

  if (myfile)  // same as: if (myfile.good())
    {
      while (getline( myfile, line ))  // same as: while (getline( myfile, line ).good())
	{
	  size_t str_pos = line.find("//",0);
	  if (str_pos!=std::string::npos) line =  line.substr(0,str_pos);
	  str_pos = line.find("%",0);
	  if (str_pos!=std::string::npos) line =  line.substr(0,str_pos);

	  str_pos = line.find("=",0);

	  if (str_pos!=std::string::npos) token = line.substr(0,str_pos);
	  trim(token);


	  if (token.compare(token_name)==0)
	    {
	      if (str_pos!=std::string::npos) value = line.substr(str_pos+1,line.size());
	      trim(value);
	      if(value==""){std::cout << "No value entered for " << token << std::endl;}
	      return;
	    }
	}

      std::cout << "Unable to find pair to key " << token_name << std::endl;
      value="";
      myfile.close();
    }
  else{std::cout << "Non-existant file" << std::endl;}	
}



void fill_dim(int &dim_loc,int pol_mode,int N_psi,std::string shape_order)
{
  if(shape_order=="NHLC" || shape_order=="HLC"){
    dim_loc=2*N_psi-3;
  }
  else if(shape_order=="LN"){
    if(std::abs(pol_mode)==1){dim_loc=2*N_psi-2;}
    else{dim_loc=2*N_psi-3;}
  }
  else if(shape_order=="HLN"){
    if(std::abs(pol_mode)==1){dim_loc=2*N_psi-2;}
    else{dim_loc=2*N_psi-3;}
  }
  else if(shape_order=="NHQL" || shape_order=="HQL"){
    if(std::abs(pol_mode)==1){dim_loc=3*N_psi-3;}
    else{dim_loc=3*N_psi-4;}
  }
  else if(shape_order=="HQD"){
    if(std::abs(pol_mode)==1){dim_loc=4*N_psi-4;}
    else{dim_loc=4*N_psi-5;}
  }
  else if(shape_order=="NHCQ" || shape_order=="HCQ"){
    if(std::abs(pol_mode)==1){dim_loc=4*N_psi-3;}
    else{dim_loc=4*N_psi-5;}
  }
  else if(shape_order=="CB")
    {
      if(std::abs(pol_mode)==1){dim_loc=4*N_psi-3;}
      else{dim_loc=4*N_psi-4;}
    }
  else if(shape_order=="NHQC")
    {
      if(std::abs(pol_mode)==1){dim_loc=5*N_psi-4;}
      else{dim_loc=5*N_psi-4;}
    }
  else{std::cout << "That shape order is not recognised (fill_dim)" << std::endl;}

}

void fill_indices(int *perp_indices[],int *wedge_indices[],int pol_mode,int N_psi,std::string shape_order)
{
  
  if(shape_order=="NHLC" || shape_order=="HLC")
    {	  
      perp_indices[0][0]=-1; perp_indices[0][N_psi-1]=-1; //Normal vanishes at centre and outer wall
      wedge_indices[0][0]=0; wedge_indices[0][N_psi-1]=-1; //Wedge defined on half-grid. N_psi-1 term is 0 as outside the domain.

      for(int iii=1;iii<N_psi-1;iii++)
	{
	  perp_indices[0][iii]=2*iii-1;
	  wedge_indices[0][iii]=2*iii;
	}
    }
  else if(shape_order=="LN" || shape_order=="HLN")
    {
      if(std::abs(pol_mode)==1){
      perp_indices[0][0]=-1; perp_indices[0][N_psi-1]=-1; //Normal vanishes at centre and outer wall
      wedge_indices[0][0]=0; wedge_indices[0][N_psi-1]=2*N_psi - 3;
      
      for(int iii=1;iii<N_psi-1;iii++)
	{
	  perp_indices[0][iii]=2*iii-1;
	  wedge_indices[0][iii]=2*iii;
	}
      }
      else{
      perp_indices[0][0]=-1; perp_indices[0][N_psi-1]=-1; //Normal vanishes at centre and outer wall
      wedge_indices[0][0]=-1; wedge_indices[0][N_psi-1]=2*N_psi - 4; //Wedge vanishes at centre.
      
      for(int iii=1;iii<N_psi-1;iii++)
	{
	  perp_indices[0][iii]=2*iii-2;
	  wedge_indices[0][iii]=2*iii-1;
	}
      }
      
    }
  else if(shape_order=="NHQL" || shape_order=="HQL")
    {	  
      if(std::abs(pol_mode)==1){
	perp_indices[0][0]=-1; perp_indices[0][N_psi-1]=-1; //Normal vanishes at centre and outer wall
	perp_indices[1][0]=1; perp_indices[1][N_psi-1]=-1; 
	wedge_indices[0][0]=0; wedge_indices[0][N_psi-1]=3*N_psi-4; //Wedge non-zero at centre and non-zero at outer wall
      
	for(int iii=1;iii<N_psi-1;iii++)
	  {
	    perp_indices[0][iii]=3*iii-1;
	    perp_indices[1][iii]=3*iii+1;
	    wedge_indices[0][iii]=3*iii;
	  }
      }
      else{
	perp_indices[0][0]=-1; perp_indices[0][N_psi-1]=-1; //Normal vanishes at centre and outer wall
	perp_indices[1][0]=0; perp_indices[1][N_psi-1]=-1; 
	wedge_indices[0][0]=-1; wedge_indices[0][N_psi-1]=3*N_psi-5; //Wedge vanishes at centre and non-zero at outer wall
      
	for(int iii=1;iii<N_psi-1;iii++)
	  {
	    perp_indices[0][iii]=3*iii-2;
	    perp_indices[1][iii]=3*iii;
	    wedge_indices[0][iii]=3*iii-1;
	  }
      }
    }
  else if(shape_order=="HQD")
    {	  
      if(std::abs(pol_mode)==1){
	perp_indices[0][0]=-1; perp_indices[0][N_psi-1]=-1; //Normal vanishes at centre and outer wall
	perp_indices[1][0]=1; perp_indices[1][N_psi-1]=-1; 
	wedge_indices[0][0]=0; wedge_indices[0][N_psi-1]=4*N_psi-5; //Wedge non-zero at centre and non-zero at outer wall
	wedge_indices[1][0]=2; wedge_indices[1][N_psi-1]=-1;
	
	for(int iii=1;iii<N_psi-1;iii++)
	  {
	    perp_indices[0][iii]=4*iii-1;
	    perp_indices[1][iii]=4*iii+1;
	    wedge_indices[0][iii]=4*iii;
	    wedge_indices[1][iii]=4*iii+2;
	  }
      }
      else{
	perp_indices[0][0]=-1; perp_indices[0][N_psi-1]=-1; //Normal vanishes at centre and outer wall
	perp_indices[1][0]=0; perp_indices[1][N_psi-1]=-1; 
	wedge_indices[0][0]=-1; wedge_indices[0][N_psi-1]=4*N_psi-6; //Wedge non-zero at centre and non-zero at outer wall
	wedge_indices[1][0]=1; wedge_indices[1][N_psi-1]=-1;
	
	for(int iii=1;iii<N_psi-1;iii++)
	  {
	    perp_indices[0][iii]=4*iii-2;
	    perp_indices[1][iii]=4*iii;
	    wedge_indices[0][iii]=4*iii-1;
	    wedge_indices[1][iii]=4*iii+1;
	  }
      }
    }
  else if(shape_order=="NHCQ" || shape_order=="HCQ")
    {	  
      if(std::abs(pol_mode)==1){
	perp_indices[0][0]=-1; perp_indices[0][N_psi-1]=-1; //Normal vanishes at centre and outer wall
	perp_indices[1][0]=1; perp_indices[1][N_psi-1]=4*N_psi-4; 
	wedge_indices[0][0]=0; wedge_indices[0][N_psi-1]=4*N_psi-5; //Wedge non-zero at centre and non-zero at outer wall
	wedge_indices[1][0]=2; wedge_indices[1][N_psi-1]=-1; 
      
	for(int iii=1;iii<N_psi-1;iii++)
	  {
	    perp_indices[0][iii]=4*iii-1;
	    perp_indices[1][iii]=4*iii+1;
	    wedge_indices[0][iii]=4*iii;
	    wedge_indices[1][iii]=4*iii+2;
	  }
      }
      else{
	perp_indices[0][0]=-1; perp_indices[0][N_psi-1]=-1; //Normal vanishes at centre and outer wall
	perp_indices[1][0]=-1; perp_indices[1][N_psi-1]=4*N_psi-6; 
	wedge_indices[0][0]=-1; wedge_indices[0][N_psi-1]=4*N_psi-7; //Wedge non-zero at centre and non-zero at outer wall
	wedge_indices[1][0]=0; wedge_indices[1][N_psi-1]=-1; 
      
	for(int iii=1;iii<N_psi-1;iii++)
	  {
	    perp_indices[0][iii]=4*iii-3;
	    perp_indices[1][iii]=4*iii-1;
	    wedge_indices[0][iii]=4*iii-2;
	    wedge_indices[1][iii]=4*iii;
	  }
      }
    }
  else if(shape_order=="CB")
    {	  
      if(std::abs(pol_mode)==1){
	perp_indices[0][0]=-1; perp_indices[0][N_psi-1]=-1; //Normal vanishes at centre and outer wall
	perp_indices[1][0]=0; perp_indices[1][N_psi-1]=4*N_psi-6; 
	wedge_indices[0][0]=1; wedge_indices[0][N_psi-1]=4*N_psi-5; //Wedge non-zero at centre and non-zero at outer wall
	wedge_indices[1][0]=-1; wedge_indices[1][N_psi-1]=4*N_psi-4; 
      
	for(int iii=1;iii<N_psi-1;iii++)
	  {
	    perp_indices[0][iii]=4*iii-2;
	    perp_indices[1][iii]=4*iii-1;
	    wedge_indices[0][iii]=4*iii;
	    wedge_indices[1][iii]=4*iii+1;
	  }
      }
      else{
	perp_indices[0][0]=-1; perp_indices[0][N_psi-1]=-1; //Normal vanishes at centre and outer wall
	perp_indices[1][0]=-1; perp_indices[1][N_psi-1]=4*N_psi-7; 
	wedge_indices[0][0]=-1; wedge_indices[0][N_psi-1]=4*N_psi-6; //Wedge non-zero at centre and non-zero at outer wall
	wedge_indices[1][0]=0; wedge_indices[1][N_psi-1]=4*N_psi-5; 
      
	for(int iii=1;iii<N_psi-1;iii++)
	  {
	    perp_indices[0][iii]=4*iii-3;
	    perp_indices[1][iii]=4*iii-2;
	    wedge_indices[0][iii]=4*iii-1;
	    wedge_indices[1][iii]=4*iii;
	  }
      }
    }
  else if(shape_order=="NHQC")
    {	  
      if(std::abs(pol_mode)==1){
	perp_indices[0][0]=-1; perp_indices[0][N_psi-1]=-1; //Normal vanishes at centre and outer wall
	perp_indices[1][0]=2; perp_indices[1][N_psi-1]=-1; //Normal midpoints
	perp_indices[2][0]=0; perp_indices[2][N_psi-1]=5*N_psi-7;  //Normal derivative
	wedge_indices[0][0]=1; wedge_indices[0][N_psi-1]=5*N_psi-6; //Wedge non-zero at centre and non-zero at outer wall
	wedge_indices[1][0]=-1; wedge_indices[1][N_psi-1]=5*N_psi-5; //Wedge derivative
      
	for(int iii=1;iii<N_psi-1;iii++)
	  {
	    perp_indices[0][iii]=5*iii-2;
	    perp_indices[1][iii]=5*iii+2;
	    perp_indices[2][iii]=5*iii-1;
	    wedge_indices[0][iii]=5*iii;
	    wedge_indices[1][iii]=5*iii+1;
	  }
      }
      else{
	perp_indices[0][0]=-1; perp_indices[0][N_psi-1]=-1; //Normal vanishes at centre and outer wall
	perp_indices[1][0]=2; perp_indices[1][N_psi-1]=-1;  //Normal midpoints
	perp_indices[2][0]=0; perp_indices[2][N_psi-1]=5*N_psi-7;  //Normal derivative
	wedge_indices[0][0]=-1; wedge_indices[0][N_psi-1]=5*N_psi-6; //Wedge non-zero at centre and non-zero at outer wall
	wedge_indices[1][0]=1; wedge_indices[1][N_psi-1]=5*N_psi-5; //Wedge derivative
      
	for(int iii=1;iii<N_psi-1;iii++)
	  {
	    perp_indices[0][iii]=5*iii-2;
	    perp_indices[1][iii]=5*iii+2;
	    perp_indices[2][iii]=5*iii-1;
	    wedge_indices[0][iii]=5*iii;
	    wedge_indices[1][iii]=5*iii+1;
	  }
      }
    }
  else{std::cout << "That shape order is not recognised (fill_indices)" << std::endl;}
  
}

void fill_2dgrid(double (&foo)(double,double),double grid[],double psi_vals[],int N_psi,double theta_vals[],int N_theta)
{
  for(int iii=0;iii<N_psi;iii++)
    {
      for(int jjj=0;jjj<N_theta;jjj++)
	{
	  grid[iii*N_theta+jjj]=foo(psi_vals[iii],theta_vals[jjj]);
	}
    }
}

void fill_2dgrid(std::complex<double> (&foo)(double,double),std::complex<double> grid[],double psi_vals[],int N_psi,double theta_vals[],int N_theta)
{
  for(int iii=0;iii<N_psi;iii++)
    {
      for(int jjj=0;jjj<N_theta;jjj++)
	{
	  grid[iii*N_theta+jjj]=foo(psi_vals[iii],theta_vals[jjj]);
	}
    }
}

void clean_grid(std::complex<double> grid[],int size)
{for(int iii=0;iii<size;iii++){grid[iii]=0.0;}}

void clean_grid(double grid[],int size)
{for(int iii=0;iii<size;iii++){grid[iii]=0.0;}}



void print_grid(double grid[],int N_psi,int N_theta)
{
  for(int i=0;i<N_psi;i++)
    {
      for(int j=0;j<N_theta;j++)
	{
	  std::cout << grid[i*N_theta+j] << "   ";
	}
      std::cout << std::endl;
    }
    std::cout << std::endl;
}

void print_grid(std::complex<double> grid[],int N_psi,int N_theta)
{
  double max=0.0;
   for(int iii=0;iii<N_psi*N_theta;iii++)
    {
      if(std::abs(grid[iii])>max){max=std::abs(grid[iii]);}
    }
  
  for(int i=0;i<N_psi;i++)
    {
      for(int j=0;j<N_theta;j++)
	{
	  std::cout << "(" << std::real(grid[i*N_theta+j]) << "," << std::imag(grid[i*N_theta+j]) << ")" << "  ";
	}
      std::cout << std::endl;
    }
    std::cout << std::endl;
}




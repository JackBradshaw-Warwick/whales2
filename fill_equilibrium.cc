#include "constants.h"

void fill_theta(equil_fields *equil,geom_shape geom)
{
  int N_interp = geom.N_interp ;
  int N_psi = geom.N_psi ;
  int N_theta = geom.N_theta ;

  double flux_max = pol_flux_theta( geom.min_rad ) ;
  double flux_min = pol_flux_theta( rad_low_rat * geom.min_rad ) ;


  fill_rad( equil, geom, flux_max, flux_min ) ;
  rad_to_psi( equil, geom.N_psi, geom.N_interp );


  double *radii=new double[N_interp];
  for(int iii=0;iii<N_interp;iii++){
    radii[iii] = sqrt( 2.0 * equil->psi_interp[iii] / geom.B_0 ); //bisect( pol_flux_theta, 0.0, geom.min_rad * 1.001, equil->psi_interp[iii] ) ; //Slightly larger than [rad_low_rat * geom.min_rad , geom.min_rad]  to ensure root is inside domain
  }


  for(int iii=0;iii<N_interp;iii++){for(int jjj=0;jjj<N_theta;jjj++){
      equil->maj_rad[iii*N_theta+jjj] = radii[iii] * cos(equil->theta_grid[jjj]) ;
      equil->height[iii*N_theta+jjj] = radii[iii] * sin(equil->theta_grid[jjj]) ;
    }}

  fill_2dgrid(dens,equil->dens,equil->psi_interp,N_interp,equil->theta_grid,N_theta);

  double *grad_psi=new double[N_interp];
  deriv_1d(grad_psi,equil->psi_interp,radii,N_interp,geom.deriv_order);

  for( int iii=0 ; iii < N_interp ; iii++ ){ equil->pres[iii*N_theta] = ( 0.5 / mu_0 ) * ( (geom.beta_0 + 1.0) * geom.B_0 * geom.B_0 - (grad_psi[iii] / radii[iii]) * (grad_psi[iii] / radii[iii]) ) ; }
  for( int iii=0 ; iii < N_interp ; iii++ ){ for( int jjj=1 ; jjj < N_theta ; jjj++ ){ equil->pres[iii*N_theta+jjj] = equil->pres[iii*N_theta] ; }}
    
  for(int iii=0 ; iii< N_interp * N_theta ; iii++){ equil->f_psi[iii] = 0.0 ; }

  for(int iii=0;iii<N_interp;iii++){for(int jjj=0;jjj<N_theta;jjj++){ 
      equil->g_pp[iii*N_theta+jjj] = grad_psi[iii] * grad_psi[iii];
      equil->g_pt[iii*N_theta+jjj] = 0.0;
      equil->g_tt[iii*N_theta+jjj] = 1.0;
      equil->g_phph[iii*N_theta+jjj] = 1.0 / ( radii[iii] * radii[iii] );
      equil->jacob[iii*N_theta+jjj] = radii[iii] / grad_psi[iii];
      equil->mag_sq[iii*N_theta+jjj] = (grad_psi[iii] / radii[iii]) * (grad_psi[iii] / radii[iii]);
      equil->curv_psi[iii*N_theta+jjj] = 0.0;
      equil->neg_shear[iii*N_theta+jjj] = 0.0;
      equil->j_dot_b[iii*N_theta+jjj] = 0.0;
      
      equil->ones[iii*N_theta+jjj] = 1.0;
    }}

  fill_cc(equil,N_interp,N_theta);
  
  equil->dpsi_dR=new double[0];
  equil->dpsi_dZ=new double[0];
  equil->dth_dR=new double[0];
  equil->dth_dZ=new double[0];
  

  delete[] radii; radii=NULL;
  delete[] grad_psi; grad_psi=NULL;
}

void fill_screw(equil_fields *equil,geom_shape geom)
{

  double flux_max = pol_flux_screw( geom.min_rad ) ;
  double flux_min = pol_flux_screw( rad_low_rat * geom.min_rad ) ;


  fill_rad( equil, geom, flux_max, flux_min ) ;
  rad_to_psi( equil, geom.N_psi, geom.N_interp );

    
  double *radii=new double[geom.N_interp];
  for(int iii=0;iii<geom.N_interp;iii++){
    radii[iii]=bisect(pol_flux_screw,0.0,geom.min_rad*1.001,equil->psi_interp[iii]); //Slightly larger than [rad_low_rat * geom.min_rad , geom.min_rad]  to ensure root is inside domain
  }

  for(int iii=0;iii<geom.N_interp;iii++){for(int jjj=0;jjj<geom.N_theta;jjj++){
      equil->maj_rad[iii*geom.N_theta+jjj]=radii[iii]*cos(equil->theta_grid[jjj]);
      equil->height[iii*geom.N_theta+jjj]=radii[iii]*sin(equil->theta_grid[jjj]);
    }}

  fill_2dgrid(dens,equil->dens,equil->psi_interp,geom.N_interp,equil->theta_grid,geom.N_theta);


  
  double *grad_psi=new double[geom.N_interp]; //Equal to B_\theta
  deriv_1d(grad_psi,equil->psi_interp,radii,geom.N_interp,geom.deriv_order);

  fill_2dgrid(f_psi,equil->f_psi,radii,geom.N_interp,equil->theta_grid,geom.N_theta);

  /*for(int iii=0;iii<geom.N_interp;iii++){for(int jjj=0;jjj<geom.N_theta;jjj++){
    equil->f_psi[iii*geom.N_theta+jjj]=-grad_psi[iii]*safety(radii[iii])/radii[iii];
    }}*/

  //for(int iii=0;iii<geom.N_interp;iii++){std::cout << "radius : " << (equil->f_psi[iii*geom.N_theta]) << " " << -2.0*radii[iii] << std::endl;}

  fill_2dgrid(pres,equil->pres,radii,geom.N_interp,equil->theta_grid,geom.N_theta);

  /*double *grad_psi=new double[geom.N_interp];
    deriv_1d(grad_psi,equil->psi_interp,radii,geom.N_interp,geom.deriv_order);

    double beta_0,B_0;
    read_in("beta_0",value);  beta_0 = std::stod( value );
    read_in("B_0",value);  B_0 = std::stod( value );

    for(int iii=0;iii<geom.N_interp;iii++){equil->pres[iii*geom.N_theta]=(0.5/mu_0)*((beta_0+1.0)*B_0-(grad_psi[iii]/radii[iii]));}
    for(int iii=0;iii<geom.N_interp;iii++){for(int jjj=1;jjj<geom.N_theta;jjj++){equil->pres[iii*geom.N_theta+jjj]=equil->pres[iii*geom.N_theta];}}*/

  //for(int iii=0;iii<geom.N_psi;iii++){std::cout << radii[iii*(geom.num_quad+1)] << std::endl;}

  delete[] radii; radii=NULL;
  delete[] grad_psi; grad_psi=NULL;
}

void fill_full(equil_fields *equil,geom_shape geom)
{/*
   std::string value;

   double flux_max,flux_min;
   read_in("flux_max",value);  flux_max = std::stod( value );

   flux_min=1.0e-8*flux_max;
  
   double sqrt_diff=(sqrt(flux_max)-sqrt(flux_min))/(geom->N_interp-1);
   for(int iii=0;iii<geom->N_interp;iii++){geom->psi_interp[iii]=(sqrt(flux_min)+iii*sqrt_diff)*(sqrt(flux_min)+iii*sqrt_diff);}

   for(int iii=0;iii<geom->N_psi;iii++){geom->psi_grid[iii]=geom->psi_interp[iii*(geom->num_quad+1)];}

   fill_2dgrid(maj_rad,equil->maj_rad,geom->psi_interp,geom->N_interp,geom->theta_grid,geom->N_theta);
   fill_2dgrid(height,equil->height,geom->psi_interp,geom->N_interp,geom->theta_grid,geom->N_theta);
   fill_2dgrid(dens,equil->dens,geom->psi_interp,geom->N_interp,geom->theta_grid,geom->N_theta);
   fill_2dgrid(pres,equil->pres,geom->psi_interp,geom->N_interp,geom->theta_grid,geom->N_theta);
   fill_2dgrid(f_psi,equil->f_psi,geom->psi_interp,geom->N_interp,geom->theta_grid,geom->N_theta);*/
}

int fill_whales(equil_fields *equil,geom_shape *geom)
{
  //Open HDF file for reading data from
  
  hid_t file_id;
  herr_t status;

  std::stringstream filepath;
  filepath << geom->input_filepath ;
      
  //Open HDF5 datafile
  file_id = H5Fopen(filepath.str().c_str(), H5F_ACC_RDONLY, H5P_DEFAULT);
  
  hid_t eqgroup_id, N_psi_id, N_theta_id, num_quad_id, m_min_id, m_max_id, m_coup_id, tor_mod_id, shape_order_id, deriv_order_id, interp_order_id;
  
  //Open group
  eqgroup_id = H5Gopen2(file_id, "/equil", H5P_DEFAULT);

  N_psi_id = H5Aopen( eqgroup_id, "N_psi", H5P_DEFAULT );
  status = H5Aread( N_psi_id, H5T_NATIVE_INT, &geom->N_psi );
  status = H5Aclose( N_psi_id );
  
  N_theta_id = H5Aopen( eqgroup_id, "N_theta", H5P_DEFAULT );
  status = H5Aread( N_theta_id, H5T_NATIVE_INT, &geom->N_theta );
  status = H5Aclose( N_theta_id );

  num_quad_id = H5Aopen( eqgroup_id, "num_quad", H5P_DEFAULT );
  status = H5Aread( num_quad_id, H5T_NATIVE_INT, &geom->num_quad );
  status = H5Aclose( num_quad_id );

  m_min_id = H5Aopen( eqgroup_id, "m_min", H5P_DEFAULT );
  status = H5Aread( m_min_id, H5T_NATIVE_INT, &geom->m_min );
  status = H5Aclose( m_min_id );

  m_max_id = H5Aopen( eqgroup_id, "m_max", H5P_DEFAULT );
  status = H5Aread( m_max_id, H5T_NATIVE_INT, &geom->m_max );
  status = H5Aclose( m_max_id );
  
  m_coup_id = H5Aopen( eqgroup_id, "m_coup", H5P_DEFAULT );
  status = H5Aread( m_coup_id, H5T_NATIVE_INT, &geom->m_coup );
  status = H5Aclose( m_coup_id );
  
  tor_mod_id = H5Aopen( eqgroup_id, "tor_mod", H5P_DEFAULT );
  status = H5Aread( tor_mod_id, H5T_NATIVE_DOUBLE, &geom->tor_mod );
  status = H5Aclose( tor_mod_id );

  hid_t string_id ;
  string_id = H5Tcopy(H5T_C_S1);
  status = H5Tset_size(string_id, H5T_VARIABLE);
  
  shape_order_id = H5Aopen( eqgroup_id, "shape_order", H5P_DEFAULT );
  status = H5Aread( shape_order_id, string_id , &geom->shape_order );
  status = H5Aclose( shape_order_id );
  
  deriv_order_id = H5Aopen( eqgroup_id, "deriv_order", H5P_DEFAULT );
  status = H5Aread( deriv_order_id, string_id, &geom->deriv_order );
  status = H5Aclose( deriv_order_id );
  
  interp_order_id = H5Aopen( eqgroup_id, "interp_order", H5P_DEFAULT );
  status = H5Aread( interp_order_id, string_id, &geom->interp_order );
  status = H5Aclose( interp_order_id );

  status = H5Tclose(string_id);

  //RECALCULATE GEOM VALUES
  calc_geom(geom);

  //Read in equil
  hid_t f_psi_id, dens_id, pres_id, g_pp_id, g_pt_id, g_tt_id, g_phph_id, jacob_id, mag_sq_id, j_dot_b_id, curv_psi_id, neg_shear_id;
  
  f_psi_id = H5Dopen2( eqgroup_id , "f_psi" , H5P_DEFAULT );  
  status = H5Dread( f_psi_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, equil->f_psi );
  status = H5Dclose( f_psi_id );

  dens_id = H5Dopen2( eqgroup_id , "dens" , H5P_DEFAULT );  
  status = H5Dread( dens_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, equil->dens );
  status = H5Dclose( dens_id );
  
  pres_id = H5Dopen2( eqgroup_id , "pres" , H5P_DEFAULT );  
  status = H5Dread( pres_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, equil->pres );
  status = H5Dclose( pres_id );
  
  g_pp_id = H5Dopen2( eqgroup_id , "g_pp" , H5P_DEFAULT );  
  status = H5Dread( g_pp_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, equil->g_pp );
  status = H5Dclose( g_pp_id );
  
  g_pt_id = H5Dopen2( eqgroup_id , "g_pt" , H5P_DEFAULT );  
  status = H5Dread( g_pt_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, equil->g_pt );
  status = H5Dclose( g_pt_id );
  
  g_tt_id = H5Dopen2( eqgroup_id , "g_tt" , H5P_DEFAULT );  
  status = H5Dread( g_tt_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, equil->g_tt );
  status = H5Dclose( g_tt_id );
  
  g_phph_id = H5Dopen2( eqgroup_id , "g_phph" , H5P_DEFAULT );  
  status = H5Dread( g_phph_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, equil->g_phph );
  status = H5Dclose( g_phph_id );
  
  jacob_id = H5Dopen2( eqgroup_id , "jacob" , H5P_DEFAULT );  
  status = H5Dread( jacob_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, equil->jacob );
  status = H5Dclose( jacob_id );
  
  mag_sq_id = H5Dopen2( eqgroup_id , "mag_sq" , H5P_DEFAULT );  
  status = H5Dread( mag_sq_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, equil->mag_sq );
  status = H5Dclose( mag_sq_id );
  
  j_dot_b_id = H5Dopen2( eqgroup_id , "j_dot_b" , H5P_DEFAULT );  
  status = H5Dread( j_dot_b_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, equil->j_dot_b );
  status = H5Dclose( j_dot_b_id );

  curv_psi_id = H5Dopen2( eqgroup_id , "curv_psi" , H5P_DEFAULT );  
  status = H5Dread( curv_psi_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, equil->curv_psi );
  status = H5Dclose( curv_psi_id );
  
  neg_shear_id = H5Dopen2( eqgroup_id , "neg_shear" , H5P_DEFAULT );  
  status = H5Dread( neg_shear_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, equil->neg_shear );
  status = H5Dclose( neg_shear_id );

  //Read in equil grids
  hid_t rad_interp_id, psi_interp_id,  drad_dpsi_id, maj_rad_id, height_id;

  rad_interp_id = H5Dopen2( eqgroup_id , "rad_interp" , H5P_DEFAULT );  
  status = H5Dread( rad_interp_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, equil->rad_interp );
  status = H5Dclose( rad_interp_id );

  psi_interp_id = H5Dopen2( eqgroup_id , "psi_interp" , H5P_DEFAULT );  
  status = H5Dread( psi_interp_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, equil->psi_interp );
  status = H5Dclose( psi_interp_id );

  drad_dpsi_id = H5Dopen2( eqgroup_id , "drad_dpsi" , H5P_DEFAULT );  
  status = H5Dread( drad_dpsi_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, equil->drad_dpsi );
  status = H5Dclose( drad_dpsi_id );

  maj_rad_id = H5Dopen2( eqgroup_id , "maj_rad" , H5P_DEFAULT );  
  status = H5Dread( maj_rad_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, equil->maj_rad );
  status = H5Dclose( maj_rad_id );

  height_id = H5Dopen2( eqgroup_id , "height" , H5P_DEFAULT );  
  status = H5Dread( height_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, equil->height );
  status = H5Dclose( height_id );

  status = H5Gclose(eqgroup_id);


  //Read in grid
  hid_t gridgroup_id, psi_grid_id, theta_grid_id, rad_var_id;

  //Open group
  gridgroup_id = H5Gopen2(file_id, "/grid", H5P_DEFAULT);
  
  psi_grid_id = H5Dopen2( gridgroup_id , "psi" , H5P_DEFAULT );  
  status = H5Dread( psi_grid_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, equil->psi_grid );
  status = H5Dclose( psi_grid_id );
  
  theta_grid_id = H5Dopen2( gridgroup_id , "theta" , H5P_DEFAULT );  
  status = H5Dread( theta_grid_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, equil->theta_grid );
  status = H5Dclose( theta_grid_id );
  
  rad_var_id = H5Dopen2( gridgroup_id , "rad_var" , H5P_DEFAULT );  
  status = H5Dread( rad_var_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, equil->rad_var );
  status = H5Dclose( rad_var_id );

  status = H5Gclose(gridgroup_id);
   
  status = H5Fclose(file_id);

  return 0;
}
  

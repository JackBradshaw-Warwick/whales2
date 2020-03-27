#include "constants.h"

int save_equil_full(geom_shape geom,equil_fields equil)
{
  //Creates HDF file and saves basic data
  save_equil_part( geom, equil);

  int N_psi=geom.N_psi;
  int N_theta=geom.N_theta;
  int N_interp=geom.N_interp;
  int num_quad=geom.num_quad;

  //Open HDF file for writing data into
  
  hid_t file_id;
  herr_t status;

  std::stringstream filepath;
  filepath << geom.output_dir << geom.equil_filename ;
      
  //Create HDF5 datafile
  file_id = H5Fopen(filepath.str().c_str(), H5F_ACC_RDWR, H5P_DEFAULT);

  hid_t eqgroup_id, eqspace_id, gridspace_id, rad_interp_id, psi_interp_id, drad_dpsi_id, maj_rad_id, height_id, f_psi_id, dens_id, pres_id, g_pp_id, g_pt_id, g_tt_id, g_phph_id, jacob_id, mag_sq_id, j_dot_b_id, curv_psi_id, neg_shear_id;
  hsize_t eq_dims[2]={N_interp,N_theta};
  hsize_t grid_dims[1]={N_interp};
  
  //Create group for data input
  eqgroup_id = H5Gopen2(file_id, "/equil", H5P_DEFAULT);

  //Create the dataspace for the dataset.
  gridspace_id = H5Screate_simple(1,grid_dims,NULL);

  //Create, write and close the dataset.
  rad_interp_id = H5Dcreate2 (eqgroup_id,"rad_interp",H5T_NATIVE_DOUBLE, gridspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite (rad_interp_id, H5T_NATIVE_DOUBLE, gridspace_id, gridspace_id, H5P_DEFAULT, equil.rad_interp);
  status = H5Dclose( rad_interp_id );
  
  psi_interp_id = H5Dcreate2 (eqgroup_id,"psi_interp",H5T_NATIVE_DOUBLE, gridspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite (psi_interp_id, H5T_NATIVE_DOUBLE, gridspace_id, gridspace_id, H5P_DEFAULT, equil.psi_interp);
  status = H5Dclose( psi_interp_id );
  
  drad_dpsi_id = H5Dcreate2 (eqgroup_id,"drad_dpsi",H5T_NATIVE_DOUBLE, gridspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite (drad_dpsi_id, H5T_NATIVE_DOUBLE, gridspace_id, gridspace_id, H5P_DEFAULT, equil.drad_dpsi);
  status = H5Dclose( drad_dpsi_id );

  status = H5Sclose(gridspace_id);

  //Create the dataspace for the dataset.
  eqspace_id = H5Screate_simple(2,eq_dims,NULL);
  
  maj_rad_id = H5Dcreate2 (eqgroup_id,"maj_rad",H5T_NATIVE_DOUBLE, eqspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite (maj_rad_id, H5T_NATIVE_DOUBLE, eqspace_id, eqspace_id, H5P_DEFAULT, equil.maj_rad);
  status = H5Dclose( maj_rad_id );
  
  height_id = H5Dcreate2 (eqgroup_id,"height",H5T_NATIVE_DOUBLE, eqspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite (height_id, H5T_NATIVE_DOUBLE, eqspace_id, eqspace_id, H5P_DEFAULT, equil.height);
  status = H5Dclose( height_id );
  
  f_psi_id = H5Dcreate2 (eqgroup_id,"f_psi",H5T_NATIVE_DOUBLE, eqspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite (f_psi_id, H5T_NATIVE_DOUBLE, eqspace_id, eqspace_id, H5P_DEFAULT, equil.f_psi);
  status = H5Dclose( f_psi_id );

  dens_id = H5Dcreate2 (eqgroup_id,"dens",H5T_NATIVE_DOUBLE, eqspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite (dens_id, H5T_NATIVE_DOUBLE, eqspace_id, eqspace_id, H5P_DEFAULT, equil.dens);
  status = H5Dclose( dens_id );
  
  pres_id = H5Dcreate2 (eqgroup_id,"pres",H5T_NATIVE_DOUBLE, eqspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite (pres_id, H5T_NATIVE_DOUBLE, eqspace_id, eqspace_id, H5P_DEFAULT, equil.pres);
  status = H5Dclose( pres_id );
  
  g_pp_id = H5Dcreate2 (eqgroup_id,"g_pp",H5T_NATIVE_DOUBLE, eqspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite (g_pp_id, H5T_NATIVE_DOUBLE, eqspace_id, eqspace_id, H5P_DEFAULT, equil.g_pp);
  status = H5Dclose( g_pp_id );
  
  g_pt_id = H5Dcreate2 (eqgroup_id,"g_pt",H5T_NATIVE_DOUBLE, eqspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite (g_pt_id, H5T_NATIVE_DOUBLE, eqspace_id, eqspace_id, H5P_DEFAULT, equil.g_pt);
  status = H5Dclose( g_pt_id );
  
  g_tt_id = H5Dcreate2 (eqgroup_id,"g_tt",H5T_NATIVE_DOUBLE, eqspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite (g_tt_id, H5T_NATIVE_DOUBLE, eqspace_id, eqspace_id, H5P_DEFAULT, equil.g_tt);
  status = H5Dclose( g_tt_id );
  
  g_phph_id = H5Dcreate2 (eqgroup_id,"g_phph",H5T_NATIVE_DOUBLE, eqspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite (g_phph_id, H5T_NATIVE_DOUBLE, eqspace_id, eqspace_id, H5P_DEFAULT, equil.g_phph);
  status = H5Dclose( g_phph_id );
  
  jacob_id = H5Dcreate2 (eqgroup_id,"jacob",H5T_NATIVE_DOUBLE, eqspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite (jacob_id, H5T_NATIVE_DOUBLE, eqspace_id, eqspace_id, H5P_DEFAULT, equil.jacob);
  status = H5Dclose( jacob_id );
  
  mag_sq_id = H5Dcreate2 (eqgroup_id,"mag_sq",H5T_NATIVE_DOUBLE, eqspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite (mag_sq_id, H5T_NATIVE_DOUBLE, eqspace_id, eqspace_id, H5P_DEFAULT, equil.mag_sq);
  status = H5Dclose( mag_sq_id );
  
  j_dot_b_id = H5Dcreate2 (eqgroup_id,"j_dot_b",H5T_NATIVE_DOUBLE, eqspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite (j_dot_b_id, H5T_NATIVE_DOUBLE, eqspace_id, eqspace_id, H5P_DEFAULT, equil.j_dot_b);
  status = H5Dclose( j_dot_b_id );
  
  curv_psi_id = H5Dcreate2 (eqgroup_id,"curv_psi",H5T_NATIVE_DOUBLE, eqspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite (curv_psi_id, H5T_NATIVE_DOUBLE, eqspace_id, eqspace_id, H5P_DEFAULT, equil.curv_psi);
  status = H5Dclose( curv_psi_id );
  
  neg_shear_id = H5Dcreate2 (eqgroup_id,"neg_shear",H5T_NATIVE_DOUBLE, eqspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite (neg_shear_id, H5T_NATIVE_DOUBLE, eqspace_id, eqspace_id, H5P_DEFAULT, equil.neg_shear);
  status = H5Dclose( neg_shear_id );


  status = H5Sclose(eqspace_id);
  status = H5Gclose(eqgroup_id);
  
  /*****************************************************************************************************************************************************************************************/
  
  status = H5Fclose(file_id);

  return 0;
}





int save_equil_part(geom_shape geom,equil_fields equil)
{
  int N_psi=geom.N_psi;
  int N_theta=geom.N_theta;
  int N_interp=geom.N_interp;
  int num_quad=geom.num_quad;


  double *k_par=new double[geom.m_range*N_psi*(N_theta+1)];
  double *alf_vel=new double[geom.m_range*N_psi*(N_theta+1)];			   
  double *alf_freq=new double[geom.m_range*N_psi*(N_theta+1)];

  for(int m_count=0;m_count<geom.m_range;m_count++){
    for(int iii=0;iii<N_psi;iii++){
      for(int jjj=0;jjj<N_theta;jjj++){
	if( geom.analytical_type == "cylinder_theta" ){
	  k_par[m_count*N_psi*(N_theta+1)+iii*(N_theta+1)+jjj]=(equil.f_psi[iii*N_theta+jjj]*equil.g_phph[iii*N_theta+jjj]*(geom.m_min+m_count)+(geom.tor_mod/equil.jacob[iii*N_theta+jjj]))/sqrt(equil.mag_sq[iii*N_theta+jjj]);
	  k_par[m_count*N_psi*(N_theta+1)+iii*(N_theta+1)+N_theta]=k_par[m_count*N_psi*(N_theta+1)+iii*(N_theta+1)];
	}
	else{
	  k_par[m_count*N_psi*(N_theta+1)+iii*(N_theta+1)+jjj]=(equil.f_psi[iii*N_theta+jjj]*equil.g_phph[iii*N_theta+jjj]*geom.tor_mod+((geom.m_min+m_count)/equil.jacob[iii*N_theta+jjj]))/sqrt(equil.mag_sq[iii*N_theta+jjj]);
	  k_par[m_count*N_psi*(N_theta+1)+iii*(N_theta+1)+N_theta]=k_par[m_count*N_psi*(N_theta+1)+iii*(N_theta+1)];
	}
	
	alf_vel[m_count*N_psi*(N_theta+1)+iii*(N_theta+1)+jjj]=sqrt(equil.mag_sq[iii*N_theta+jjj]/(mu_0*equil.dens[iii*N_theta+jjj]));
	alf_vel[m_count*N_psi*(N_theta+1)+iii*(N_theta+1)+N_theta]=alf_vel[m_count*N_psi*(N_theta+1)+iii*(N_theta+1)];
	
	alf_freq[m_count*N_psi*(N_theta+1)+iii*(N_theta+1)+jjj]=k_par[m_count*N_psi*(N_theta+1)+iii*(N_theta+1)+jjj]*alf_vel[m_count*N_psi*(N_theta+1)+iii*(N_theta+1)+jjj];
	alf_freq[m_count*N_psi*(N_theta+1)+iii*(N_theta+1)+N_theta]=alf_freq[m_count*N_psi*(N_theta+1)+iii*(N_theta+1)];
      }}}


     
    
  /***************************************************************************************************************************************************************************************/
  /***************************************************************************************************************************************************************************************/
  //Create HDF file for writing data into
  
  hid_t file_id;
  herr_t status;

  std::stringstream filepath;
  filepath << geom.output_dir << geom.equil_filename ; 
      
  //Create HDF5 datafile
  file_id = H5Fcreate(filepath.str().c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT); //Overwrite file if it already exists

  /***************************************************************************************************************************************************************************************/
  //Write grid data
  
  //GRIDPOINTS
  hid_t gridgroup_id,gridspace_id,Rgrid_id,Zgrid_id,rad_var_id,psi_id,theta_id;
  hsize_t grid_dims[2]={N_psi,N_theta};
  double *grid_in=new double [N_psi*N_theta];

  //Create group for grid data input
  gridgroup_id = H5Gcreate2(file_id, "/grid", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  //Create dataspace for outputting grid data
  gridspace_id = H5Screate_simple(2,grid_dims,NULL);

  //Create the dataset, write xdata and close.
  Rgrid_id = H5Dcreate2(gridgroup_id, "R", H5T_NATIVE_DOUBLE, gridspace_id,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  for(int iii=0;iii<N_psi;iii++){for(int jjj=0;jjj<N_theta;jjj++){grid_in[iii*N_theta+jjj]=equil.maj_rad[iii*(num_quad+1)*N_theta+jjj];}}
  status = H5Dwrite (Rgrid_id, H5T_NATIVE_DOUBLE, gridspace_id, gridspace_id, H5P_DEFAULT, grid_in);
  status = H5Dclose(Rgrid_id);

  //Create the dataset, write ydata and close.
  Zgrid_id = H5Dcreate2(gridgroup_id, "Z", H5T_NATIVE_DOUBLE, gridspace_id,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  for(int iii=0;iii<N_psi;iii++){for(int jjj=0;jjj<N_theta;jjj++){grid_in[iii*N_theta+jjj]=equil.height[iii*(num_quad+1)*N_theta+jjj];}}
  status = H5Dwrite (Zgrid_id, H5T_NATIVE_DOUBLE, gridspace_id, gridspace_id, H5P_DEFAULT, grid_in);
  status = H5Dclose(Zgrid_id);

  delete[] grid_in;
  status = H5Sclose(gridspace_id);

  hsize_t grid_dim[1]={N_psi};
  gridspace_id = H5Screate_simple(1,grid_dim,NULL);

  psi_id = H5Dcreate2(gridgroup_id, "psi", H5T_NATIVE_DOUBLE, gridspace_id,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite (psi_id, H5T_NATIVE_DOUBLE, gridspace_id, gridspace_id, H5P_DEFAULT, equil.psi_grid);
  status = H5Dclose(psi_id);

  rad_var_id = H5Dcreate2(gridgroup_id, "rad_var", H5T_NATIVE_DOUBLE, gridspace_id,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite (rad_var_id, H5T_NATIVE_DOUBLE, gridspace_id, gridspace_id, H5P_DEFAULT, equil.rad_var);
  status = H5Dclose(rad_var_id);

  status = H5Sclose(gridspace_id);

  grid_dim[0]=N_theta;
  gridspace_id = H5Screate_simple(1,grid_dim,NULL);

  theta_id = H5Dcreate2(gridgroup_id, "theta", H5T_NATIVE_DOUBLE, gridspace_id,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite (theta_id, H5T_NATIVE_DOUBLE, gridspace_id, gridspace_id, H5P_DEFAULT, equil.theta_grid);
  status = H5Dclose(theta_id);

  status = H5Sclose(gridspace_id);
  status = H5Gclose(gridgroup_id);


  /***************************************************************************************************************************************************************************************/
  //Write alfven data
  
  //ALFVEN DATA
  hid_t alfgroup_id,alfspace_id,k_id,vel_id,freq_id;
  hsize_t alf_dims[3]={geom.m_range,N_psi,N_theta+1};

  //Create group for data input
  alfgroup_id = H5Gcreate2(file_id, "/alfven", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      
  //Create the data space for the dataset.
  alfspace_id = H5Screate_simple(3,alf_dims,NULL);

  k_id = H5Dcreate2(alfgroup_id, "k_par", H5T_NATIVE_DOUBLE, alfspace_id,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite (k_id, H5T_NATIVE_DOUBLE, alfspace_id, alfspace_id, H5P_DEFAULT, k_par);
  status = H5Dclose(k_id);

  vel_id = H5Dcreate2(alfgroup_id, "vel", H5T_NATIVE_DOUBLE, alfspace_id,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite (vel_id, H5T_NATIVE_DOUBLE, alfspace_id, alfspace_id, H5P_DEFAULT, alf_vel);
  status = H5Dclose(vel_id);

  freq_id = H5Dcreate2(alfgroup_id, "freq", H5T_NATIVE_DOUBLE, alfspace_id,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Dwrite (freq_id, H5T_NATIVE_DOUBLE, alfspace_id, alfspace_id, H5P_DEFAULT, alf_freq);
  status = H5Dclose(freq_id);
      
      
  status = H5Sclose(alfspace_id);
  status = H5Gclose(alfgroup_id);
  
  delete[] k_par;
  delete[] alf_vel;
  delete[] alf_freq;

  /*****************************************************************************************************************************************************************************************/

  hid_t eqgroup_id, eqspace_id, N_psi_id, N_theta_id, num_quad_id, m_min_id, m_max_id, m_coup_id, tor_mod_id, shape_order_id, deriv_order_id, interp_order_id;
  hsize_t eq_dims[1]={1};
  
  //Create group for data input
  eqgroup_id = H5Gcreate2(file_id, "/equil", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

  //Create the data space for the attribute.
  eqspace_id = H5Screate_simple(1,eq_dims,NULL);

  //Create, write and close the dataset attribute.
  N_psi_id = H5Acreate2 (eqgroup_id,"N_psi",H5T_NATIVE_INT, eqspace_id, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Awrite( N_psi_id, H5T_NATIVE_INT, &geom.N_psi);
  status = H5Aclose( N_psi_id);

  N_theta_id = H5Acreate2 (eqgroup_id,"N_theta",H5T_NATIVE_INT, eqspace_id, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Awrite( N_theta_id, H5T_NATIVE_INT, &geom.N_theta);
  status = H5Aclose( N_theta_id);

  num_quad_id = H5Acreate2 (eqgroup_id,"num_quad",H5T_NATIVE_INT, eqspace_id, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Awrite( num_quad_id, H5T_NATIVE_INT, &geom.num_quad);
  status = H5Aclose( num_quad_id);

  m_min_id = H5Acreate2 (eqgroup_id,"m_min",H5T_NATIVE_INT, eqspace_id, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Awrite( m_min_id, H5T_NATIVE_INT, &geom.m_min);
  status = H5Aclose( m_min_id);
  
  m_max_id = H5Acreate2 (eqgroup_id,"m_max",H5T_NATIVE_INT, eqspace_id, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Awrite( m_max_id, H5T_NATIVE_INT, &geom.m_max);
  status = H5Aclose( m_max_id);

  m_coup_id = H5Acreate2 (eqgroup_id,"m_coup",H5T_NATIVE_INT, eqspace_id, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Awrite( m_coup_id, H5T_NATIVE_INT, &geom.m_coup);
  status = H5Aclose( m_coup_id);
 
  tor_mod_id = H5Acreate2 (eqgroup_id,"tor_mod",H5T_NATIVE_DOUBLE, eqspace_id, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Awrite( tor_mod_id, H5T_NATIVE_DOUBLE, &geom.tor_mod);
  status = H5Aclose( tor_mod_id);

  hid_t string_id ;
  string_id = H5Tcopy(H5T_C_S1);
  status = H5Tset_size(string_id, H5T_VARIABLE);

  shape_order_id = H5Acreate2 (eqgroup_id,"shape_order", string_id, eqspace_id, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Awrite( shape_order_id, string_id, &geom.shape_order );
  status = H5Aclose( shape_order_id);

  
  deriv_order_id = H5Acreate2 (eqgroup_id,"deriv_order",string_id, eqspace_id, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Awrite( deriv_order_id, string_id, &geom.deriv_order);
  status = H5Aclose( deriv_order_id);
  
  interp_order_id = H5Acreate2 (eqgroup_id,"interp_order",string_id, eqspace_id, H5P_DEFAULT, H5P_DEFAULT);
  status = H5Awrite( interp_order_id, string_id, &geom.interp_order);
  status = H5Aclose( interp_order_id);

  status = H5Tclose(string_id);

  status = H5Sclose(eqspace_id);
  status = H5Gclose(eqgroup_id);
  
  /*****************************************************************************************************************************************************************************************/
  
  status = H5Fclose(file_id);

  return 0;
}


int load_mats(geom_shape geom , Mat mats[] )
{
  int ierr;
  
  std::stringstream tmp_filename_E;
  tmp_filename_E <<  geom.mats_dir << geom.mats_filename << "_E.dat"; 
	
  PetscViewer mat_in_E;
  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,tmp_filename_E.str().c_str(),FILE_MODE_READ,&mat_in_E); CHKERRQ(ierr);
  ierr = MatLoad(mats[0],mat_in_E);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&mat_in_E); CHKERRQ(ierr);

  ierr = MatAssemblyBegin(mats[0],MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(mats[0],MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  int counter;
  if(geom.Hall_on==1){
    std::stringstream tmp_filename_H;
    tmp_filename_H <<  geom.mats_dir << geom.mats_filename << "_H.dat"; 
	
    PetscViewer mat_in_H;
    ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,tmp_filename_H.str().c_str(),FILE_MODE_READ,&mat_in_H); CHKERRQ(ierr);
    ierr = MatLoad(mats[1],mat_in_H);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&mat_in_H); CHKERRQ(ierr);

    ierr = MatAssemblyBegin(mats[1],MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
    ierr = MatAssemblyEnd(mats[1],MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

    counter=2;
  }
  else{
    counter=1;
  }


  std::stringstream tmp_filename_F;
  tmp_filename_F << geom.mats_dir << geom.mats_filename << "_F.dat"; 
	
  PetscViewer mat_in_F;
  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,tmp_filename_F.str().c_str(),FILE_MODE_READ,&mat_in_F); CHKERRQ(ierr);
  ierr = MatLoad(mats[counter],mat_in_F);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&mat_in_F); CHKERRQ(ierr);
 
  ierr = MatAssemblyBegin(mats[counter],MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);
  ierr = MatAssemblyEnd(mats[counter],MAT_FINAL_ASSEMBLY);CHKERRQ(ierr);

  std::cout << "Mats read" << std::endl;

  return ierr;
}

int save_mats(geom_shape geom , Mat mats[] )
{
  int ierr;

  std::stringstream tmp_filename_E;
  tmp_filename_E <<  geom.mats_dir << geom.mats_filename << "_E.dat"; 
	
  PetscViewer mat_out_E;
  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,tmp_filename_E.str().c_str(),FILE_MODE_WRITE,&mat_out_E); CHKERRQ(ierr);
  ierr = MatView(mats[0],mat_out_E);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&mat_out_E); CHKERRQ(ierr);

  int counter;
  if(geom.Hall_on==1){
    std::stringstream tmp_filename_H;
    tmp_filename_H <<  geom.mats_dir << geom.mats_filename << "_H.dat"; 
	
    PetscViewer mat_out_H;
    ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,tmp_filename_H.str().c_str(),FILE_MODE_WRITE,&mat_out_H); CHKERRQ(ierr);
    ierr = MatView(mats[1],mat_out_H);CHKERRQ(ierr);
    ierr = PetscViewerDestroy(&mat_out_H); CHKERRQ(ierr);

    counter=2;
  }
  else{
    counter=1;
  }

  std::stringstream tmp_filename_F;
  tmp_filename_F << geom.mats_dir << geom.mats_filename << "_F.dat"; 
	
  PetscViewer mat_out_F;
  ierr = PetscViewerBinaryOpen(PETSC_COMM_WORLD,tmp_filename_F.str().c_str(),FILE_MODE_WRITE,&mat_out_F); CHKERRQ(ierr);
  ierr = MatView(mats[counter],mat_out_F);CHKERRQ(ierr);
  ierr = PetscViewerDestroy(&mat_out_F); CHKERRQ(ierr);

  std::cout << "Mats read" << std::endl;

  return ierr;
}

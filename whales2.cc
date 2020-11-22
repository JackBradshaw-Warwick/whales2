#include "constants.h"

double m;


int main(int argc,char* argv[])
{
  clock_t start_time,end_time;
  start_time=clock();
  
  /***************************************************************************************************************************************************************************************/
  /***************************************************************************************************************************************************************************************/


  Mat E,H,F; /* E \omega^2 + H \omega + F = 0 */
  EPS eps; /* eigenproblem solver context */
  PEP pep; /* Polynomial solver */
  Mat *mats;
  
  equil_fields equil;
  geom_shape geom;
  matrix_coeffs coeffs;

  PetscInt nconv,my_rank;
  PetscErrorCode ierr;
  

  ierr = SlepcInitialize(&argc,&argv,(char*)0,(char*)0); if ( ierr ){ return ierr; }
  MPI_Comm_rank(PETSC_COMM_WORLD,&my_rank);

  build_geom(&geom);

  build_equil(&equil,&geom);

  if(my_rank==0){
    if ( geom.write_equil == 1 ){ save_equil_full(geom,equil); }
    else { save_equil_part(geom,equil); }
  }

  /***************************************************************************************************************************************************************************************/
  /***************************************************************************************************************************************************************************************/
 
  ierr = MatCreate(PETSC_COMM_WORLD,&E); CHKERRQ(ierr);
  if(geom.Hall_on == 1){ ierr = MatCreate(PETSC_COMM_WORLD,&H); CHKERRQ(ierr); }
  ierr = MatCreate(PETSC_COMM_WORLD,&F); CHKERRQ(ierr);

  if( geom.load_mats == 1 )
    {
      
      if(geom.Hall_on == 1){
	mats = new Mat[3];
	mats[0] = F ; mats[1] = H ; mats[2] = E ;
      }
      else{
	mats = new Mat[2];
	mats[0] = F ; mats[1] = E ;
      }

      load_mats( geom , mats ) ;

    }
  else
    {
      init_coeffs(&coeffs,geom);
      
      ierr = MatSetSizes(E,PETSC_DECIDE,PETSC_DECIDE,geom.dim,geom.dim);CHKERRQ(ierr);
      ierr = MatSetFromOptions(E);CHKERRQ(ierr);
      ierr = MatSetUp(E);CHKERRQ(ierr);

      if(geom.Hall_on == 1){
	ierr = MatSetSizes(H,PETSC_DECIDE,PETSC_DECIDE,geom.dim,geom.dim);CHKERRQ(ierr);
	ierr = MatSetFromOptions(H);CHKERRQ(ierr);
	ierr = MatSetUp(H);CHKERRQ(ierr);
      }
      
      ierr = MatSetSizes(F,PETSC_DECIDE,PETSC_DECIDE,geom.dim,geom.dim);CHKERRQ(ierr);
      ierr = MatSetFromOptions(F);CHKERRQ(ierr);
      ierr = MatSetUp(F);CHKERRQ(ierr);

      if(geom.Hall_on == 1){
	mats = new Mat[3];
	mats[0] = F ; mats[1] = H ; mats[2] = E ;
      }
      else{
	mats = new Mat[2];
	mats[0] = F ; mats[1] = E ;
      }

      build_matrices( mats , geom , equil , &coeffs) ;

      //delete_full_matrix_coeffs(&coeffs);

      /*
	make_hermitian(&E,dim);
	if(geom.Hall_on == 1){
	make_hermitian(&H,dim);
	}
	make_hermitian(&F,dim);
      */

      /*
	calc_hermitian_diff(E,dim,dim);
	if(geom.Hall_on == 1){
	calc_hermitian_diff(H,dim,dim);
	}
	calc_hermitian_diff(F,dim,dim);
      */

      if(!(geom.write_mats==0)){ save_mats( geom , mats ); }
    }

  /*
    ierr = MatView(E,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
    if(geom.Hall_on == 1){
    ierr = MatView(H,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
    }
    ierr = MatView(F,PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);
  */

  clock_t time1;
  if( my_rank == 0 ){ time1 = clock() ;}
  
  if( geom.Hall_on == 1 ){
    ierr = MatScale(mats[0],-1.0);CHKERRQ(ierr);

    ierr = PEPCreate(PETSC_COMM_WORLD,&pep);CHKERRQ(ierr);
    ierr = PEPSetOperators(pep,3,mats);CHKERRQ(ierr);
    ierr = PEPSetProblemType(pep,PEP_GENERAL);CHKERRQ(ierr);
  
    ierr = PEPSetFromOptions(pep);CHKERRQ(ierr);
    //ierr = PEPSetUp(pep);CHKERRQ(ierr);

    ierr = PEPSolve(pep);CHKERRQ(ierr);
    ierr = PEPGetConverged(pep,&nconv);CHKERRQ(ierr);
  }
  else{
    ierr = EPSCreate(PETSC_COMM_WORLD,&eps);CHKERRQ(ierr);
    ierr = EPSSetOperators(eps,F,E);CHKERRQ(ierr);
    ierr = EPSSetProblemType(eps,EPS_GHIEP);CHKERRQ(ierr);
  
    ierr = EPSSetFromOptions(eps);CHKERRQ(ierr);
    //ierr = EPSSetUp(eps);CHKERRQ(ierr);

    ierr = EPSSolve(eps);CHKERRQ(ierr);
    ierr = EPSGetConverged(eps,&nconv);CHKERRQ(ierr);
  }

  time1=clock()-time1;
  if(my_rank==0){std::cout << "Solver took " << static_cast<double>(time1)/static_cast<double>(CLOCKS_PER_SEC) << " seconds" << std::endl;}
  if(my_rank==0){std::cout << "There are: " << nconv << " successfully converged eigenvalues" << std::endl;}

  ierr=MatDestroy(&E);CHKERRQ(ierr);
  if( geom.Hall_on == 1){ ierr=MatDestroy(&H);CHKERRQ(ierr); }
  ierr=MatDestroy(&F);CHKERRQ(ierr);

  //Move all output data to root process
  Vec *eigenvec_tmp = new Vec[nconv];
  Vec *eigenvec = new Vec[nconv];
  PetscScalar *eigenval = new PetscScalar[nconv];

  for(PetscInt iii=0;iii<nconv;iii++)
    {
      ierr = VecCreate(PETSC_COMM_WORLD,&eigenvec_tmp[iii]);CHKERRQ(ierr);
      ierr = VecSetSizes(eigenvec_tmp[iii],PETSC_DECIDE,geom.dim);
      ierr = VecSetFromOptions(eigenvec_tmp[iii]);
      ierr = VecAssemblyBegin(eigenvec_tmp[iii]);CHKERRQ(ierr);
      ierr = VecAssemblyEnd(eigenvec_tmp[iii]);CHKERRQ(ierr);

      if( geom.Hall_on == 1 ){
	ierr = PEPGetEigenpair(pep,iii,&eigenval[iii],NULL,eigenvec_tmp[iii],NULL);CHKERRQ(ierr);
      }
      else{
	ierr = EPSGetEigenpair(eps,iii,&eigenval[iii],NULL,eigenvec_tmp[iii],NULL);CHKERRQ(ierr);
      }
    }

  IS index;
  ierr = ISCreateStride(PETSC_COMM_WORLD,geom.dim,0,1,&index);

  for(PetscInt iii=0;iii<nconv;iii++)
    {
      //Initialise locally held vector (sequential on each thread)
      ierr = VecCreate(PETSC_COMM_SELF,&eigenvec[iii]);CHKERRQ(ierr);
      ierr = VecSetType(eigenvec[iii],VECSEQ);CHKERRQ(ierr);
      ierr = VecSetSizes(eigenvec[iii],PETSC_DECIDE,geom.dim);
      ierr = VecSetFromOptions(eigenvec[iii]);
      ierr = VecAssemblyBegin(eigenvec[iii]);CHKERRQ(ierr);
      ierr = VecAssemblyEnd(eigenvec[iii]);CHKERRQ(ierr);

      //Scatter values from parallel to sequential vectors
      VecScatter ctx_r;
      ierr = VecScatterCreate(eigenvec_tmp[iii],index,eigenvec[iii],index,&ctx_r);
      ierr = VecScatterBegin(ctx_r,eigenvec_tmp[iii],eigenvec[iii],INSERT_VALUES,SCATTER_FORWARD);
      ierr = VecScatterEnd(ctx_r,eigenvec_tmp[iii],eigenvec[iii],INSERT_VALUES,SCATTER_FORWARD);
      ierr = VecScatterDestroy(&ctx_r);

    }
  ierr = ISDestroy(&index);

  //Delete parallel-stored vectors
  for(int iii=0;iii<nconv;iii++){ierr=VecDestroy(&eigenvec_tmp[iii]);CHKERRQ(ierr);}
  delete[] eigenvec_tmp;


  /***************************************************************************************************************************************************************************************/
  /***************************************************************************************************************************************************************************************/

  if( geom.Hall_on == 1 ){
    
  }
  else{
    //Output EPS information
    std::stringstream EPS_dat;
    EPS_dat << geom.output_dir << "EPS.dat";
    PetscViewer view_EPS_dat;
    ierr = PetscViewerASCIIOpen(PETSC_COMM_WORLD,EPS_dat.str().c_str(),&view_EPS_dat);CHKERRQ(ierr);
    ierr = EPSView(eps,view_EPS_dat);CHKERRQ(ierr);
				     
    ierr = PetscViewerDestroy(&view_EPS_dat);
    ierr = EPSDestroy(&eps);CHKERRQ(ierr);
  }
  
  /************************************************************************************************************************************************************************************/
  /************************************************************************************************************************************************************************************/
  //Only master thread participates in post-processing
  
  if(!(my_rank==0))
    {
      for(int iii=0;iii<nconv;iii++){ierr=VecDestroy(&eigenvec[iii]);CHKERRQ(ierr);}
      delete[] eigenvec;
      delete[] eigenval;
    }
  else
    {
      int keep_indices[nconv];

      PetscInt nconv_mod=0;
      PetscReal polmode_norm[geom.m_range];
      double *rel_pow_vals = new double[nconv*geom.m_range];
      int *all_pol_mods= new int[nconv];

      const PetscScalar *tmp_arr;
      
      for(PetscInt iii=0;iii<nconv;iii++)
	{
	  //Find "main" pol_mod by calculating relative Fourier powers
	  int pol_mod = find_pol_mode(polmode_norm,geom,eigenvec,iii) ;
	    
	  keep_indices[nconv_mod]=iii;
	  all_pol_mods[iii]=pol_mod;

	  for(int jjj=0;jjj<geom.m_range;jjj++)
	    {
	      rel_pow_vals[nconv_mod*geom.m_range+jjj]=polmode_norm[jjj];
	    }
		      
	  nconv_mod++;
	}

      //ierr = VecView(eigenvec[keep_indices[0]],PETSC_VIEWER_STDOUT_WORLD);CHKERRQ(ierr);

      std::cout << "Number of relevant modes is: " << nconv_mod << std::endl;

      /************************************************************************************************************************************************************************************/

      hid_t file_id;
      herr_t status;

      std::stringstream FILE;
      FILE << geom.output_dir << geom.results_filename;
      
      if(geom.results_filename == geom.equil_filename){
	//Open HDF5 datafile
	file_id = H5Fopen(FILE.str().c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
      }
      else{
	//Create HDF5 datafile
	file_id = H5Fcreate(FILE.str().c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
      }
      
      hid_t datagroup_id;
      
      //EIGENVALUE DATA
      hid_t ev_id_r, ev_id_i,evspace_id,pol_mod_id;
      hsize_t ev_dims[1]={nconv_mod};
      double *ev_data_in = new double[ev_dims[0]]; //Create data array for data input
  
      //Create group for data input
      datagroup_id = H5Gcreate2(file_id, "/data", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

      //Create the data space for the attribute.
      evspace_id = H5Screate_simple(1,ev_dims,NULL);

      //AUTHOR NOTE: Maybe shouldn't take up all attribute space with these?
      //Create, write and close the dataset attribute.
      ev_id_r = H5Acreate2 (datagroup_id,"eigenvalue_r",H5T_NATIVE_DOUBLE, evspace_id, H5P_DEFAULT, H5P_DEFAULT);
      for(int iii=0;iii<ev_dims[0];iii++){ev_data_in[iii]=PetscRealPart(eigenval[keep_indices[iii]]);}
      status = H5Awrite(ev_id_r,H5T_NATIVE_DOUBLE,ev_data_in);
      status = H5Aclose(ev_id_r);
      
      ev_id_i = H5Acreate2 (datagroup_id,"eigenvalue_i",H5T_NATIVE_DOUBLE, evspace_id, H5P_DEFAULT, H5P_DEFAULT);
      for(int iii=0;iii<ev_dims[0];iii++){ev_data_in[iii]=PetscImaginaryPart(eigenval[keep_indices[iii]]);}
      status = H5Awrite(ev_id_i,H5T_NATIVE_DOUBLE,ev_data_in);
      status = H5Aclose(ev_id_i);

      pol_mod_id = H5Acreate2 (datagroup_id,"pol_mod",H5T_NATIVE_DOUBLE, evspace_id, H5P_DEFAULT, H5P_DEFAULT);
      for(int iii=0;iii<ev_dims[0];iii++){ev_data_in[iii]=static_cast<double>(all_pol_mods[keep_indices[iii]]);}
      status = H5Awrite(pol_mod_id,H5T_NATIVE_DOUBLE,ev_data_in);
      status = H5Aclose(pol_mod_id);

      status = H5Sclose(evspace_id);
      delete[] eigenval;
      delete[] ev_data_in;
      delete[] all_pol_mods;


      //EIGENVECTOR (INTERPOLATED) DATA
      hid_t dataspace_id,xidata_id_r,nudata_id_r,xidata_id_i,nudata_id_i,b_par_id_r,b_par_id_i,b_perp_id_r,b_perp_id_i,b_wedge_id_r,b_wedge_id_i;
      hsize_t data_dims[3]={nconv_mod,geom.N_psi,geom.N_theta};
      
      //Create the data space for the dataset.
      dataspace_id = H5Screate_simple(3,data_dims,NULL);

      xidata_id_r = H5Dcreate2(datagroup_id, "xi_r", H5T_NATIVE_DOUBLE, dataspace_id,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      status = H5Dclose(xidata_id_r);

      xidata_id_i = H5Dcreate2(datagroup_id, "xi_i", H5T_NATIVE_DOUBLE, dataspace_id,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      status = H5Dclose(xidata_id_i);

      nudata_id_r = H5Dcreate2(datagroup_id, "nu_r", H5T_NATIVE_DOUBLE, dataspace_id,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      status = H5Dclose(nudata_id_r);

      nudata_id_i = H5Dcreate2(datagroup_id, "nu_i", H5T_NATIVE_DOUBLE, dataspace_id,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      status = H5Dclose(nudata_id_i);

      b_par_id_r = H5Dcreate2(datagroup_id, "b_par_r", H5T_NATIVE_DOUBLE, dataspace_id,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      status = H5Dclose(b_par_id_r);

      b_par_id_i = H5Dcreate2(datagroup_id, "b_par_i", H5T_NATIVE_DOUBLE, dataspace_id,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      status = H5Dclose(b_par_id_i);

      b_perp_id_r = H5Dcreate2(datagroup_id, "b_perp_r", H5T_NATIVE_DOUBLE, dataspace_id,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      status = H5Dclose(b_perp_id_r);

      b_perp_id_i = H5Dcreate2(datagroup_id, "b_perp_i", H5T_NATIVE_DOUBLE, dataspace_id,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      status = H5Dclose(b_perp_id_i);

      b_wedge_id_r = H5Dcreate2(datagroup_id, "b_wedge_r", H5T_NATIVE_DOUBLE, dataspace_id,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      status = H5Dclose(b_wedge_id_r);

      b_wedge_id_i = H5Dcreate2(datagroup_id, "b_wedge_i", H5T_NATIVE_DOUBLE, dataspace_id,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      status = H5Dclose(b_wedge_id_i);

      //MODE DATA
      hid_t modspace_id,mod_id;
      hsize_t mod_dim[1]={geom.m_range};
      double *mod_in = new double[mod_dim[0]]; //Create data array for data input

      //Create the data space for the attribute.
      modspace_id = H5Screate_simple(1,mod_dim,NULL);

      //Create, write and close the dataset attribute.
      mod_id = H5Acreate2 (datagroup_id,"m_val",H5T_NATIVE_DOUBLE, modspace_id, H5P_DEFAULT, H5P_DEFAULT);
      for(int iii=0;iii<mod_dim[0];iii++){mod_in[iii]=geom.m_min+iii;}
      status = H5Awrite(mod_id,H5T_NATIVE_DOUBLE,mod_in);
      status = H5Aclose(mod_id);

      status = H5Sclose(modspace_id);
      delete[] mod_in;
      
      //SEPERATE POLOIDAL MODE DATA
      hid_t polspace_id,xipol_id_r,nupol_id_r,xipol_id_i,nupol_id_i;
      hsize_t pol_dims[3]={nconv_mod,geom.m_range,geom.N_psi};

      //Create the data space for the dataset.
      polspace_id = H5Screate_simple(3,pol_dims,NULL);

      xipol_id_r = H5Dcreate2(datagroup_id, "xi_r_polmodes", H5T_NATIVE_DOUBLE, polspace_id,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      status = H5Dclose(xipol_id_r);
      
      xipol_id_i = H5Dcreate2(datagroup_id, "xi_i_polmodes", H5T_NATIVE_DOUBLE, polspace_id,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      status = H5Dclose(xipol_id_i);
      
      nupol_id_r = H5Dcreate2(datagroup_id, "nu_r_polmodes", H5T_NATIVE_DOUBLE, polspace_id,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      status = H5Dclose(nupol_id_r);
      
      nupol_id_i = H5Dcreate2(datagroup_id, "nu_i_polmodes", H5T_NATIVE_DOUBLE, polspace_id,H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      status = H5Dclose(nupol_id_i);
      

      //RELATIVE MODE POWER DATA
      hid_t pow_id,powspace_id;
      hsize_t pow_dims[2]={nconv_mod,geom.m_range};
      //Create the data space for the attribute.
      powspace_id = H5Screate_simple(2,pow_dims,NULL);

      //Create, write and close the dataset attribute.
      pow_id = H5Dcreate2 (datagroup_id,"mode_power",H5T_NATIVE_DOUBLE, powspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
      status = H5Dwrite(pow_id,H5T_NATIVE_DOUBLE, powspace_id, powspace_id, H5P_DEFAULT, rel_pow_vals);
      status = H5Dclose(pow_id);

      status = H5Sclose(powspace_id);
      delete[] rel_pow_vals;

      /************************************************************************************************************************************************************************************/

      //EIGENVECTOR (INTERPOLATED) DATA
      double *data_in = new double[geom.N_psi*geom.N_theta]; //Create data array for data input
      double *pol_in = new double[geom.m_range*geom.N_psi]; //Create data array for data input
      
      PetscScalar *xi_j = new PetscScalar[geom.N_psi*geom.m_range];
      PetscScalar *xi_full = new PetscScalar[geom.N_psi*geom.N_theta];

      PetscScalar *nu_j = new PetscScalar[geom.N_psi*geom.m_range];
      PetscScalar *nu_full = new PetscScalar[geom.N_psi*geom.N_theta];

      PetscScalar *b_par = new PetscScalar[geom.N_psi*geom.N_theta];
      PetscScalar *b_perp = new PetscScalar[geom.N_psi*geom.N_theta];
      PetscScalar *b_wedge = new PetscScalar[geom.N_psi*geom.N_theta];

      /***************************************************************************************************************************************************************************************/


      
      //Create the data space for the dataset.
      dataspace_id = H5Screate_simple(3,data_dims,NULL);

      hsize_t countf[3] = { 1 , 1 , 1 } ;
      hsize_t blockf[3] = { 1 , geom.N_psi , geom.N_theta } ;
      hsize_t blockfp[3] = { 1 , geom.m_range , geom.N_psi } ;

      hsize_t countd[1] = { 1 } ;
      hsize_t blockd[1] = { geom.N_psi * geom.N_theta } ;
      hsize_t offsetd[1] = { 0 } ;
      hid_t memdata_id;
      hsize_t memdata_dims[1]={geom.N_psi*geom.N_theta};
      memdata_id = H5Screate_simple(1,memdata_dims,NULL);
      status = H5Sselect_hyperslab (memdata_id, H5S_SELECT_SET, offsetd, NULL, countd, blockd);

      hsize_t blockp[1] = { geom.m_range * geom.N_psi } ;
      hid_t mempol_id;
      hsize_t mempol_dims[1]={ geom.m_range * geom.N_psi };
      mempol_id = H5Screate_simple(1,mempol_dims,NULL);
      status = H5Sselect_hyperslab (mempol_id, H5S_SELECT_SET, offsetd, NULL, countd, blockp);
      
      for(int iii = 0 ; iii < nconv_mod ; iii++){


	//Select eigenval and extract
	sort_eigenvecs(eigenvec,xi_j,nu_j,geom,keep_indices,iii);

	//Destroy vec. Needs rethinking a little if actually filter modes (with keep indices)
	ierr = VecDestroy(&eigenvec[iii]); CHKERRQ(ierr);	

	//Rotate result so that Re(\xi_\perp) =/= 0 along line \theta=0
	rotate(xi_j,nu_j,geom.N_psi,geom.m_range);

	//Set eigenvector storage to zero
	for(int iii=0;iii<geom.N_psi*geom.N_theta;iii++){xi_full[iii]=0.0;nu_full[iii]=0.0;}
      
	for(int jjj=0;jjj<geom.N_psi;jjj++){
	  for(int m_count=0;m_count<geom.m_range;m_count++){
	    for(int kkk=0;kkk<geom.N_theta;kkk++){
		xi_full[jjj*geom.N_theta+kkk]+=xi_j[m_count*(geom.N_psi)+jjj]*(cos((geom.m_min+m_count)*equil.theta_grid[kkk])+PETSC_i*sin((geom.m_min+m_count)*equil.theta_grid[kkk]));
		nu_full[jjj*geom.N_theta+kkk]+=nu_j[m_count*(geom.N_psi)+jjj]*(cos((geom.m_min+m_count)*equil.theta_grid[kkk])+PETSC_i*sin((geom.m_min+m_count)*equil.theta_grid[kkk]));
	      }
	  }
	}

	convert_to_mag(b_par,b_perp,b_wedge,xi_full,nu_full,equil,geom,1,geom.shape_order);

	//Define space in file into which data is written
	hsize_t offsetf[3] = { iii , 0 , 0 } ;
	status = H5Sselect_hyperslab (dataspace_id, H5S_SELECT_SET, offsetf, NULL, countf, blockf);
	  
	//Create and write the displacement data.
	xidata_id_r = H5Dopen2(datagroup_id, "xi_r", H5P_DEFAULT);	
	for(int jjj=0;jjj<data_dims[1]*data_dims[2];jjj++){data_in[jjj]=PetscRealPart(xi_full[jjj]);}
	status = H5Dwrite(xidata_id_r, H5T_NATIVE_DOUBLE, memdata_id, dataspace_id, H5P_DEFAULT, data_in);
	status = H5Dclose(xidata_id_r);

	xidata_id_i = H5Dopen2(datagroup_id, "xi_i", H5P_DEFAULT);
	for(int jjj=0;jjj<data_dims[1]*data_dims[2];jjj++){data_in[jjj]=PetscImaginaryPart(xi_full[jjj]);}
	status = H5Dwrite(xidata_id_i, H5T_NATIVE_DOUBLE, memdata_id, dataspace_id, H5P_DEFAULT, data_in);
	status = H5Dclose(xidata_id_i);

	nudata_id_r = H5Dopen2(datagroup_id, "nu_r", H5P_DEFAULT);
	for(int jjj=0;jjj<data_dims[1]*data_dims[2];jjj++){data_in[jjj]=PetscRealPart(nu_full[jjj]);}
	status = H5Dwrite(nudata_id_r, H5T_NATIVE_DOUBLE, memdata_id, dataspace_id, H5P_DEFAULT, data_in);
	status = H5Dclose(nudata_id_r);

	nudata_id_i = H5Dopen2(datagroup_id, "nu_i", H5P_DEFAULT);
	for(int jjj=0;jjj<data_dims[1]*data_dims[2];jjj++){data_in[jjj]=PetscImaginaryPart(nu_full[jjj]);}
	status = H5Dwrite(nudata_id_i, H5T_NATIVE_DOUBLE, memdata_id, dataspace_id, H5P_DEFAULT, data_in);
	status = H5Dclose(nudata_id_i);

	b_par_id_r = H5Dopen2(datagroup_id, "b_par_r", H5P_DEFAULT);
	for(int jjj=0;jjj<data_dims[1]*data_dims[2];jjj++){data_in[jjj]=PetscRealPart(b_par[jjj]);}
	status = H5Dwrite (b_par_id_r, H5T_NATIVE_DOUBLE, memdata_id, dataspace_id, H5P_DEFAULT, data_in);
	status = H5Dclose(b_par_id_r);

	b_par_id_i = H5Dopen2(datagroup_id, "b_par_i", H5P_DEFAULT);
	for(int jjj=0;jjj<data_dims[1]*data_dims[2];jjj++){data_in[jjj]=PetscImaginaryPart(b_par[jjj]);}
	status = H5Dwrite (b_par_id_i, H5T_NATIVE_DOUBLE, memdata_id, dataspace_id, H5P_DEFAULT, data_in);
	status = H5Dclose(b_par_id_i);

	b_perp_id_r = H5Dopen2(datagroup_id, "b_perp_r", H5P_DEFAULT);
	for(int jjj=0;jjj<data_dims[1]*data_dims[2];jjj++){data_in[jjj]=PetscRealPart(b_perp[jjj]);}
	status = H5Dwrite (b_perp_id_r, H5T_NATIVE_DOUBLE, memdata_id, dataspace_id, H5P_DEFAULT, data_in);
	status = H5Dclose(b_perp_id_r);

	b_perp_id_i = H5Dopen2(datagroup_id, "b_perp_i", H5P_DEFAULT);
	for(int jjj=0;jjj<data_dims[1]*data_dims[2];jjj++){data_in[jjj]=PetscImaginaryPart(b_perp[jjj]);}
	status = H5Dwrite (b_perp_id_i, H5T_NATIVE_DOUBLE, memdata_id, dataspace_id, H5P_DEFAULT, data_in);
	status = H5Dclose(b_perp_id_i);

	b_wedge_id_r = H5Dopen2(datagroup_id, "b_wedge_r", H5P_DEFAULT);
	for(int jjj=0;jjj<data_dims[1]*data_dims[2];jjj++){data_in[jjj]=PetscRealPart(b_wedge[jjj]);}
	status = H5Dwrite (b_wedge_id_r, H5T_NATIVE_DOUBLE, memdata_id, dataspace_id, H5P_DEFAULT, data_in);
	status = H5Dclose(b_wedge_id_r);

	b_wedge_id_i = H5Dopen2(datagroup_id, "b_wedge_i", H5P_DEFAULT);
	for(int jjj=0;jjj<data_dims[1]*data_dims[2];jjj++){data_in[jjj]=PetscImaginaryPart(b_wedge[jjj]);}
	status = H5Dwrite (b_wedge_id_i, H5T_NATIVE_DOUBLE, memdata_id, dataspace_id, H5P_DEFAULT, data_in);
	status = H5Dclose(b_wedge_id_i);

	//Define space in file into which data is written
	status = H5Sselect_hyperslab (polspace_id, H5S_SELECT_SET, offsetf, NULL, countf, blockfp);

	//Write polspace data
	xipol_id_r = H5Dopen2(datagroup_id, "xi_r_polmodes", H5P_DEFAULT);
	for(int jjj=0;jjj<pol_dims[1]*pol_dims[2];jjj++){pol_in[jjj]=PetscRealPart(xi_j[jjj]);}
	status = H5Dwrite (xipol_id_r, H5T_NATIVE_DOUBLE, mempol_id, polspace_id, H5P_DEFAULT, pol_in);
	status = H5Dclose(xipol_id_r);

	xipol_id_i = H5Dopen2(datagroup_id, "xi_i_polmodes", H5P_DEFAULT);
	for(int jjj=0;jjj<pol_dims[1]*pol_dims[2];jjj++){pol_in[jjj]=PetscImaginaryPart(xi_j[jjj]);}
	status = H5Dwrite (xipol_id_i, H5T_NATIVE_DOUBLE, mempol_id, polspace_id, H5P_DEFAULT, pol_in);
	status = H5Dclose(xipol_id_i);

	nupol_id_r = H5Dopen2(datagroup_id, "nu_r_polmodes", H5P_DEFAULT);
	for(int jjj=0;jjj<pol_dims[1]*pol_dims[2];jjj++){pol_in[jjj]=PetscRealPart(nu_j[jjj]);}
	status = H5Dwrite (nupol_id_r, H5T_NATIVE_DOUBLE, mempol_id, polspace_id, H5P_DEFAULT, pol_in);
	status = H5Dclose(nupol_id_r);

	nupol_id_i = H5Dopen2(datagroup_id, "nu_i_polmodes", H5P_DEFAULT);
	for(int jjj=0;jjj<pol_dims[1]*pol_dims[2];jjj++){pol_in[jjj]=PetscImaginaryPart(nu_j[jjj]);}
	status = H5Dwrite (nupol_id_i, H5T_NATIVE_DOUBLE, mempol_id, polspace_id, H5P_DEFAULT, pol_in);
	status = H5Dclose(nupol_id_i);	
      }
 
      status = H5Sclose(dataspace_id);
      status = H5Sclose(polspace_id);
      
      delete[] eigenvec;
      delete[] data_in;
      delete[] xi_full;
      delete[] nu_full;

      delete[] b_par;
      delete[] b_perp;
      delete[] b_wedge;

      delete[] pol_in;
      delete[] xi_j;
      delete[] nu_j;   

      //Close datagroup
      status = H5Gclose(datagroup_id);

      //Close file
      status = H5Fclose(file_id);
    }

  delete_full_geom_shape(&geom);
  delete_full_equil_fields(&equil);
 
  
  /***************************************************************************************************************************************************************************************/
  /***************************************************************************************************************************************************************************************/
  
  ierr = SlepcFinalize();


  if(my_rank==0){end_time=clock();std::cout << "Whole program took " << (end_time-start_time)/static_cast<double>(CLOCKS_PER_SEC) << " seconds" << std::endl;}
  
  return ierr;
}


#include "constants.h"

void read_geom(geom_shape *geom)
{
  //Read from input deck
  std::string value;

  read_in("fill_type",geom->fill_type);
  if( geom->fill_type == "" ) { geom->fill_type = "analytic" ; }
    
  if( geom->fill_type == "analytic" ){  
    read_in("analytical_type",geom->analytical_type);
    if( geom->analytical_type == "" ) { geom->analytical_type = "cylinder_theta" ; }  
  }
  else if( geom->fill_type == "numerical" ){
    read_in("numerical_type",geom->numerical_type);
    if( geom->numerical_type == "" ) { geom->numerical_type = "whales" ; }

    read_in("input_filepath",geom->input_filepath);
    if( geom->input_filepath == "" ) { geom->input_filepath = "./Results/data.h5" ; }  
  }
  else { std::cout << "Fill_type is not recognised as it has been entered" << std::endl; }

  if( geom->numerical_type == "whales"){
    hid_t file_id;
    herr_t status;

    std::stringstream filepath;
    filepath << geom->input_filepath ;
      
    //Open HDF5 datafile
    file_id = H5Fopen(filepath.str().c_str(), H5F_ACC_RDONLY, H5P_DEFAULT); 
  
    hid_t eqgroup_id;
    //ints and doubles
    hid_t N_psi_id, N_theta_id, num_quad_id, m_min_id, m_max_id, m_coup_id, tor_mod_id ,  dens_form_id , A_dens_id , B_dens_id , mu_id , nu_id , min_rad_id , beta_0_id , B_0_id , R_0_id , elong_id , triang_id , alpha_sol_id , solov_A_id , solov_C_id , flux_max_id , shear_on_id , Hall_on_id ; 
    //strings
    hid_t quad_type_id , shape_order_id, deriv_order_id, interp_order_id ;

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

    dens_form_id = H5Aopen( eqgroup_id, "dens_form", H5P_DEFAULT );
    status = H5Aread( dens_form_id, H5T_NATIVE_INT, &geom->dens_form );
    status = H5Aclose( dens_form_id );

    A_dens_id = H5Aopen( eqgroup_id, "A_dens", H5P_DEFAULT );
    status = H5Aread( A_dens_id, H5T_NATIVE_DOUBLE, &geom->A_dens );
    status = H5Aclose( A_dens_id );

    B_dens_id = H5Aopen( eqgroup_id, "B_dens", H5P_DEFAULT );
    status = H5Aread( B_dens_id, H5T_NATIVE_DOUBLE, &geom->B_dens );
    status = H5Aclose( B_dens_id );

    mu_id = H5Aopen( eqgroup_id, "mu", H5P_DEFAULT );
    status = H5Aread( mu_id, H5T_NATIVE_DOUBLE, &geom->mu );
    status = H5Aclose( mu_id );

    nu_id = H5Aopen( eqgroup_id, "nu", H5P_DEFAULT );
    status = H5Aread( nu_id, H5T_NATIVE_DOUBLE, &geom->nu );
    status = H5Aclose( nu_id );

    min_rad_id = H5Aopen( eqgroup_id, "min_rad", H5P_DEFAULT );
    status = H5Aread( min_rad_id, H5T_NATIVE_DOUBLE, &geom->min_rad );
    status = H5Aclose( min_rad_id );

    beta_0_id = H5Aopen( eqgroup_id, "beta_0", H5P_DEFAULT );
    status = H5Aread( beta_0_id, H5T_NATIVE_DOUBLE, &geom->beta_0 );
    status = H5Aclose( beta_0_id );

    B_0_id = H5Aopen( eqgroup_id, "B_0", H5P_DEFAULT );
    status = H5Aread( B_0_id, H5T_NATIVE_DOUBLE, &geom->B_0 );
    status = H5Aclose( B_0_id );

    R_0_id = H5Aopen( eqgroup_id, "R_0", H5P_DEFAULT );
    status = H5Aread( R_0_id, H5T_NATIVE_DOUBLE, &geom->R_0 );
    status = H5Aclose( R_0_id );

    elong_id = H5Aopen( eqgroup_id, "elong", H5P_DEFAULT );
    status = H5Aread( elong_id, H5T_NATIVE_DOUBLE, &geom->elong );
    status = H5Aclose( elong_id );

    triang_id = H5Aopen( eqgroup_id, "triang", H5P_DEFAULT );
    status = H5Aread( triang_id, H5T_NATIVE_DOUBLE, &geom->triang );
    status = H5Aclose( triang_id );

    alpha_sol_id = H5Aopen( eqgroup_id, "alpha_sol", H5P_DEFAULT );
    status = H5Aread( alpha_sol_id, H5T_NATIVE_DOUBLE, &geom->alpha_sol );
    status = H5Aclose( alpha_sol_id );

    solov_A_id = H5Aopen( eqgroup_id, "solov_A", H5P_DEFAULT );
    status = H5Aread( solov_A_id, H5T_NATIVE_DOUBLE, &geom->solov_A );
    status = H5Aclose( solov_A_id );

    solov_C_id = H5Aopen( eqgroup_id, "solov_C", H5P_DEFAULT );
    status = H5Aread( solov_C_id, H5T_NATIVE_DOUBLE, &geom->solov_C );
    status = H5Aclose( solov_C_id );

    flux_max_id = H5Aopen( eqgroup_id, "flux_max", H5P_DEFAULT );
    status = H5Aread( flux_max_id, H5T_NATIVE_DOUBLE, &geom->flux_max );
    status = H5Aclose( flux_max_id );
  
    shear_on_id = H5Aopen( eqgroup_id, "shear_on", H5P_DEFAULT );
    status = H5Aread( shear_on_id, H5T_NATIVE_INT, &geom->shear_on );
    status = H5Aclose( shear_on_id );

    Hall_on_id = H5Aopen( eqgroup_id, "Hall_on", H5P_DEFAULT );
    status = H5Aread( Hall_on_id, H5T_NATIVE_INT, &geom->Hall_on );
    status = H5Aclose( Hall_on_id );

    hid_t string_id ;
    string_id = H5Tcopy(H5T_C_S1);
    status = H5Tset_size(string_id, H5T_VARIABLE);

    char* temp;
    
    quad_type_id = H5Aopen( eqgroup_id, "quad_type", H5P_DEFAULT );
    status = H5Aread( quad_type_id, string_id , &temp );
    status = H5Aclose( quad_type_id );
    geom->quad_type = temp ; 
  
    shape_order_id = H5Aopen( eqgroup_id, "shape_order", H5P_DEFAULT );
    status = H5Aread( shape_order_id, string_id , &temp );
    status = H5Aclose( shape_order_id );
    geom->shape_order = temp ;
  
    deriv_order_id = H5Aopen( eqgroup_id, "deriv_order", H5P_DEFAULT );
    status = H5Aread( deriv_order_id, string_id, &temp );
    status = H5Aclose( deriv_order_id );
    geom->deriv_order = temp ;
  
    interp_order_id = H5Aopen( eqgroup_id, "interp_order", H5P_DEFAULT );
    status = H5Aread( interp_order_id, string_id, &temp );
    status = H5Aclose( interp_order_id );
    geom->interp_order = temp ;

    status = H5free_memory(temp);
    status = H5Tclose(string_id);

    status = H5Gclose(eqgroup_id);

    status = H5Fclose(file_id);

    read_in("output_dir",geom->output_dir);
    if( geom->output_dir == "" ) { geom->output_dir = "./Results/" ; }

    read_in("results_filename",geom->results_filename);
    if( geom->results_filename == "" ) { geom->results_filename = "data.h5" ; }

    read_in("write_equil",value);
    if( !(value == "") ) { geom->write_equil = std::stoi( value ) ; }
    else { geom->write_equil = 0 ; }

    read_in("equil_filename",geom->equil_filename);
    if( geom->equil_filename == "" ) { geom->equil_filename = geom->results_filename ; }

    read_in("load_mats",value);
    if( !(value == "") ) { geom->load_mats = std::stoi( value ) ; }
    else { geom->load_mats = 0 ; }

    //If loading matrices then need to load whales equilibrium
    if( geom->load_mats == 1 ){
      geom->fill_type = "numerical";
      geom->numerical_type = "whales" ;

      read_in("input_filepath",geom->input_filepath);
      if( geom->input_filepath == "" ) { geom->input_filepath = "./Results/data.h5" ; }  
    }

    geom->write_mats = 0 ; //Default value
    if( geom->load_mats == 0  ){ // If mats file exists then no need to write 
      read_in("write_mats",value);
      if( !(value == "") ) { geom->write_mats = std::stoi( value ) ; }

      if( geom->write_mats == 1 ){ //If writing mats then need to write whales equil
	geom->write_equil = 1 ;
      }   
    }

    if( geom->load_mats == 1 || geom->write_mats == 1){    
      read_in("mats_dir",geom->mats_dir);
      if( geom->mats_dir == "" ) { geom->mats_dir = geom->output_dir ; }
  
      read_in("mats_filename",geom->mats_filename);
      if( geom->mats_filename == "" ) { geom->mats_filename = "mat" ; }
    }
  }
  else{ //If not loading equilibrium from file, read from input.txt
    
    read_in("N_psi",value);
    if( !(value == "") ) { geom->N_psi = std::stoi( value ) ; }
    else { geom->N_psi = 4 ; }
  
    read_in("N_theta",value);
    if( !(value == "") ) { geom->N_theta = std::stoi( value ) ; }
    else { geom->N_theta = 4 ; }

    read_in("num_quad",value);
    if( !(value == "") ) { geom->num_quad = std::stoi( value ) ; }
    else { geom->num_quad = 4 ; }

    if(!(geom->num_quad == 4 || geom->num_quad == 6 || geom->num_quad == 12 || geom->num_quad == 18)){
      std::cout << "That value for num_quad is invalid. Setting num_quad = 4." << std::endl;
      geom->num_quad = 4 ;
    }

    read_in("quad_type",value);
    if( !(value == "") ) { geom->quad_type = value ; }
    else { geom->quad_type = "gq" ; }
  
    read_in("m_min",value);
    if( !(value == "") ) { geom->m_min = std::stoi( value ) ; }
    else { geom->m_min = 0 ; }
  
    read_in("m_max",value);
    if( !(value == "") ) { geom->m_max = std::stoi( value ) ; }
    else { geom->m_max = 1 ; }
  
    //Calculated here for convenience
    geom->m_range = geom->m_max - geom->m_min ;
  
    read_in("m_coup",value);
    if( !(value == "") ) { geom->m_coup = std::stoi( value ) ; }
    else { geom->m_coup = geom->m_range ; }

    read_in("tor_mod",value);
    if( !(value == "") ) { geom->tor_mod = std::stod( value ) ; }
    else { geom->tor_mod = 0.0 ; }

    read_in("shape_order",geom->shape_order); //Can read directly in as string
    if( geom->shape_order == "" ) { geom->shape_order = "NHLC" ; } //Default
    
    read_in("deriv_order",geom->deriv_order);
    if( geom->deriv_order == "" ) { geom->deriv_order = "Quadratic" ; }
    
    read_in("interp_order",geom->interp_order);
    if( geom->interp_order == "" ) { geom->interp_order = "Cubic_pol" ; }

    if( geom->analytical_type == "cylinder_theta" || geom->analytical_type == "cylinder_screw" ){ geom->m_coup = 0 ; }

    read_in("dens_form",value);
    if( !(value == "") ) { geom->dens_form = std::stoi( value ) ; }
    else { geom->dens_form = 0 ; }

    read_in("A_dens",value);
    if( !(value == "") ) { geom->A_dens = std::stod( value ) ; }
    else { geom->A_dens = 1.0 ; }

    assert( geom->A_dens > 0.0 );

    read_in("B_dens",value);
    if( !(value == "") ) { geom->B_dens = std::stod( value ) ; }
    else { geom->B_dens = 0.0 ; }

    assert( geom->B_dens < 1.0 );

    read_in("mu",value);
    if( !(value == "") ) { geom->mu = std::stod( value ) ; }
    else { geom->mu = 0.0 ; }

    assert( geom->mu >= 0.0 );

    read_in("nu",value);
    if( !(value == "") ) { geom->nu = std::stod( value ) ; }
    else { geom->nu = 1.0 ; }

    read_in("min_rad",value);
    if( !(value == "") ) { geom->min_rad = std::stod( value ) ; }
    else { geom->min_rad = 1.0 ; }

    read_in("beta_0",value);
    if( !(value == "") ) { geom->beta_0 = std::stod( value ) ; }
    else { geom->beta_0 = 0.0 ; }

    read_in("B_0",value);
    if( !(value == "") ) { geom->B_0 = std::stod( value ) ; }
    else { geom->B_0 = 1.0 ; }

    read_in("R_0",value);
    if( !(value == "") ) { geom->R_0 = std::stod( value ) ; }
    else { geom->R_0 = 3.0 ; }

    read_in("elong",value);
    if( !(value == "") ) { geom->elong = std::stod( value ) ; }
    else { geom->elong = 1.0 ; }

    read_in("triang",value);
    if( !(value == "") ) { geom->triang = std::stod( value ) ; }
    else { geom->triang = 0.0 ; }
    
    read_in("alpha_sol",value);
    if( !(value == "") ) { geom->alpha_sol = std::stod( value ) ; }
    else { geom->alpha_sol = 1.0 ; }
    
    read_in("solov_A",value);
    if( !(value == "") ) { geom->solov_A = std::stod( value ) ; }
    else { geom->solov_A = 1.0 ; }

    read_in("solov_C",value);
    if( !(value == "") ) { geom->solov_C = std::stod( value ) ; }
    else { geom->solov_C = 0.0 ; } 

    read_in("flux_max",value);
    if( !(value == "") ) { geom->flux_max = std::stod( value ) ; }
    else { geom->flux_max = 1.0 ; }

    read_in("output_dir",geom->output_dir);
    if( geom->output_dir == "" ) { geom->output_dir = "./Results/" ; }

    read_in("results_filename",geom->results_filename);
    if( geom->results_filename == "" ) { geom->results_filename = "data.h5" ; }

    read_in("write_equil",value);
    if( !(value == "") ) { geom->write_equil = std::stoi( value ) ; }
    else { geom->write_equil = 0 ; }

    read_in("equil_filename",geom->equil_filename);
    if( geom->equil_filename == "" ) { geom->equil_filename = geom->results_filename ; }

    read_in("load_mats",value);
    if( !(value == "") ) { geom->load_mats = std::stoi( value ) ; }
    else { geom->load_mats = 0 ; }

    //If loading matrices then need to load whales equilibrium
    if( geom->load_mats == 1 ){
      geom->fill_type = "numerical";
      geom->numerical_type = "whales" ;

      read_in("input_filepath",geom->input_filepath);
      if( geom->input_filepath == "" ) { geom->input_filepath = "./Results/data.h5" ; }  
    }

    geom->write_mats = 0 ; //Default value
    if( geom->load_mats == 0  ){ // If mats file exists then no need to write 
      read_in("write_mats",value);
      if( !(value == "") ) { geom->write_mats = std::stoi( value ) ; }

      if( geom->write_mats == 1 ){ //If writing mats then need to write whales equil
	geom->write_equil = 1 ;
      }   
    }

    if( geom->load_mats == 1 || geom->write_mats == 1){    
      read_in("mats_dir",geom->mats_dir);
      if( geom->mats_dir == "" ) { geom->mats_dir = geom->output_dir ; }
  
      read_in("mats_filename",geom->mats_filename);
      if( geom->mats_filename == "" ) { geom->mats_filename = "mat" ; }
    }

    read_in("shear_on",value);
    if( !(value == "") ) { geom->shear_on = std::stoi( value ) ; }
    else { geom->shear_on = 1 ; }

    read_in("Hall_on",value);
    if( !(value == "") ) { geom->Hall_on = std::stoi( value ) ; }
    else { geom->Hall_on = 1 ; }
  }
}

void calc_geom(geom_shape *geom)
{
  //Calculate from values read in
  
  geom->m_range = geom->m_max - geom->m_min ;

  int half_size=(geom->N_theta/2)-1; if(!(geom->N_theta%2==0)){ half_size=(geom->N_theta-1)/2; }
  geom->fourier_size_sym=half_size+1;
  geom->fourier_size_full=2*half_size+1;

  geom->N_interp=geom->N_psi+(geom->N_psi-1)*geom->num_quad;

  //Calculate dimensions and matrix boundaries of m-modes
  geom->dim = 0 ;
  geom->pol_pos = new int[geom->m_range] ; geom->pol_pos[0] = 0 ;

  for(int m_count=0;m_count<3;m_count++){ fill_dim(geom->dims_loc[m_count],m_count,geom->N_psi,geom->shape_order); }
  
  int dim_temp=0;
  geom->pol_pos[0]=0;
  for(int m_count=0;m_count<geom->m_range;m_count++){
    fill_dim(dim_temp,geom->m_min+m_count,geom->N_psi,geom->shape_order);
    geom->dim+=dim_temp;
    if(!(m_count==geom->m_range-1)){geom->pol_pos[m_count+1]=geom->pol_pos[m_count]+dim_temp;}
  }

  //For higher order elements add more variables 
  if(geom->shape_order=="NHLC" || geom->shape_order=="HLC" || geom->shape_order=="LN" || geom->shape_order=="HLN"){geom->num_var_perp=1; geom->num_var_wedge=1;}
  else if(geom->shape_order=="NHQL" || geom->shape_order=="HQL"){geom->num_var_perp=2; geom->num_var_wedge=1;}
  else if(geom->shape_order=="NHCQ" || geom->shape_order=="HCQ" || geom->shape_order=="CB" || geom->shape_order=="HQD"){geom->num_var_perp=2; geom->num_var_wedge=2;}
  else if(geom->shape_order=="NHQC" ){geom->num_var_perp=3; geom->num_var_wedge=2;}
  else{std::cout << "That shape order is not recognised (calc_geom)" << std::endl;}

  
  geom->ind_perp_main=new int*[geom->num_var_perp]; geom->ind_perp_sec=new int*[geom->num_var_perp];
  geom->ind_wedge_main=new int*[geom->num_var_wedge]; geom->ind_wedge_sec=new int*[geom->num_var_wedge];

  for(int i=0;i<geom->num_var_perp;i++){geom->ind_perp_main[i]=new int[geom->N_psi]; geom->ind_perp_sec[i]=new int[geom->N_psi];}
  for(int i=0;i<geom->num_var_wedge;i++){geom->ind_wedge_main[i]=new int[geom->N_psi]; geom->ind_wedge_sec[i]=new int[geom->N_psi]; }
}

void build_geom(geom_shape *geom)
{
  read_geom(geom);
  calc_geom(geom);
}


void build_equil(equil_fields *equil,geom_shape *geom)
{
  init_equil(equil,geom->N_psi,geom->N_interp,geom->N_theta);

  //Theta is some angle, periodic in 2*pi. Grid is equally spaced in theta but not necessarily poloidal arclength
  for(int iii=0 ; iii < geom->N_theta ; iii++){ equil->theta_grid[iii] = iii * 2.0 * pi / geom->N_theta ;} 
  
  if( geom->fill_type == "analytic" ){
    
    if( geom->analytical_type == "cylinder_theta" ){ fill_theta(equil,*geom); }
    else if( geom->analytical_type == "cylinder_screw" ){ fill_screw(equil,*geom); }
    else if( geom->analytical_type == "soloviev" ){ fill_sol(equil,*geom); }
    else if( geom->analytical_type == "solov_simp" ){ fill_sol_simp(equil,*geom); }
    else if( geom->analytical_type == "full" ){ fill_full(equil,*geom); }
    else{ std::cout << "Analytical type not recognised" << std::endl; }

    
    if(!( geom->analytical_type =="cylinder_theta" )){

      //Convert psi to rad.
      //psi_to_rad(equil,geom->N_psi,geom->N_interp);
  
      build_jacobian(equil,*geom);
      build_equil_funcs(equil,*geom);
    }
    
    //Is this necessary?
    udsym(equil, *geom);
    
  }
  else if( geom->fill_type == "numerical" ){

    if( geom->numerical_type == "whales" ){ fill_whales(equil,geom); }
    else{ std::cout << "Numerical type not recognised" << std::endl; }

  }
  else{ std::cout << "Fill value not recognised" << std::endl; }
 
}


void build_jacobian(equil_fields *equil,geom_shape geom)
{
  std::string deriv_order=geom.deriv_order;
    
  int N_interp=geom.N_interp;
  int N_theta=geom.N_theta;
  
  double *dR_dpsi=new double[N_interp*N_theta];
  double *dR_dth=new double[N_interp*N_theta];
  double *dZ_dpsi=new double[N_interp*N_theta];
  double *dZ_dth=new double[N_interp*N_theta];

  deriv_1d(dR_dpsi,equil->maj_rad,equil->rad_interp,N_interp,N_theta,true,deriv_order); //Do derivative in terms of sqrt_psi as this will give more accurate result (since R~sqrt(psi))
  for(int iii=0;iii<N_interp;iii++){for(int jjj=0;jjj<N_theta;jjj++){ dR_dpsi[iii*N_theta+jjj] *= equil->drad_dpsi[iii] ;}}
  deriv_ang(dR_dth,equil->maj_rad,equil->theta_grid,N_theta,N_interp,false);
  deriv_1d(dZ_dpsi,equil->height,equil->rad_interp,N_interp,N_theta,true,deriv_order);
  for(int iii=0;iii<N_interp;iii++){for(int jjj=0;jjj<N_theta;jjj++){ dZ_dpsi[iii*N_theta+jjj] *= equil->drad_dpsi[iii] ;}}
  deriv_ang( dZ_dth, equil->height, equil->theta_grid, N_theta, N_interp, false);
  

  double *determ=new double[N_interp*N_theta];
  for(int iii=0;iii<N_interp*N_theta;iii++){determ[iii]=dR_dpsi[iii]*dZ_dth[iii]-dR_dth[iii]*dZ_dpsi[iii];}

  for(int iii=0;iii<N_interp*N_theta;iii++)
    {
      if(determ[iii]==0.0 && dZ_dth[iii]==0.0){equil->dpsi_dR[iii]=0.0;}
      else if(determ[iii]==0.0 && !(dZ_dth[iii]==0.0)){std::cout << "Error: determinant element = 0" << std::endl; return;}
      else{equil->dpsi_dR[iii]=dZ_dth[iii]/determ[iii];}
    }

  for(int iii=0;iii<N_interp*N_theta;iii++)
    {
      if(determ[iii]==0.0 && dR_dth[iii]==0.0){equil->dpsi_dZ[iii]=0.0;}
      else if(determ[iii]==0.0 && !(dR_dth[iii]==0.0)){std::cout << "Error: determinant element = 0" << std::endl; return;}
      else{equil->dpsi_dZ[iii]=-dR_dth[iii]/determ[iii];}
    }

  for(int iii=0;iii<N_interp*N_theta;iii++)
    {
      if(determ[iii]==0.0 && dZ_dpsi[iii]==0.0){equil->dth_dR[iii]=0.0;}
      else if(determ[iii]==0.0 && !(dZ_dpsi[iii]==0.0)){std::cout << "Error: determinant element = 0" << std::endl; return;}
      else{equil->dth_dR[iii]=-dZ_dpsi[iii]/determ[iii];}
    }

  for(int iii=0;iii<N_interp*N_theta;iii++)
    {
      if(determ[iii]==0.0 && dR_dpsi[iii]==0.0){equil->dth_dZ[iii]=0.0;}
      else if(determ[iii]==0.0 && !(dR_dpsi[iii]==0.0)){std::cout << "Error: determinant element = 0" << std::endl; return;}
      else{equil->dth_dZ[iii]=dR_dpsi[iii]/determ[iii];}
    }

  delete[] dR_dpsi;
  delete[] dR_dth;
  delete[] dZ_dpsi;
  delete[] dZ_dth;
  delete[] determ;
}

void build_equil_funcs(equil_fields *equil,geom_shape geom)
{
  std::string deriv_order=geom.deriv_order;
    
  int N_interp=geom.N_interp;
  int N_theta=geom.N_theta;
  
  if( geom.analytical_type == "cylinder_screw" ){
    for(int iii=0;iii<N_interp*N_theta;iii++){
      equil->g_pp[iii]= equil->dpsi_dR[iii] * equil->dpsi_dR[iii] + equil->dpsi_dZ[iii] * equil->dpsi_dZ[iii] ;
      equil->g_pt[iii]=0.0;
      equil->g_tt[iii]=1.0/(equil->maj_rad[iii]*equil->maj_rad[iii]+equil->height[iii]*equil->height[iii]);
      equil->g_phph[iii]=1.0;
    }
  }
  else{
    for(int iii=0;iii<N_interp*N_theta;iii++){
      equil->g_pp[iii]=equil->dpsi_dR[iii]*equil->dpsi_dR[iii]+equil->dpsi_dZ[iii]*equil->dpsi_dZ[iii];
      equil->g_pt[iii]=equil->dpsi_dR[iii]*equil->dth_dR[iii]+equil->dpsi_dZ[iii]*equil->dth_dZ[iii];
      equil->g_tt[iii]=equil->dth_dR[iii]*equil->dth_dR[iii]+equil->dth_dZ[iii]*equil->dth_dZ[iii];
      equil->g_phph[iii]=1.0/(equil->maj_rad[iii]*equil->maj_rad[iii]);
    }
  }

  for(int iii=0;iii<N_interp*N_theta;iii++){
    equil->jacob[iii]=1.0/sqrt((equil->g_pp[iii]*equil->g_tt[iii]-equil->g_pt[iii]*equil->g_pt[iii])*equil->g_phph[iii]);
    equil->mag_sq[iii]=(equil->f_psi[iii]*equil->f_psi[iii]+equil->g_pp[iii])*equil->g_phph[iii];
    equil->ones[iii]=1.0;
  }


  delete[] equil->dpsi_dR; equil->dpsi_dR=NULL;
  delete[] equil->dpsi_dZ; equil->dpsi_dZ=NULL;
  delete[] equil->dth_dR; equil->dth_dR=NULL;
  delete[] equil->dth_dZ; equil->dth_dZ=NULL;
  double *temp_1=new double[N_interp*N_theta];
  double *temp_2=new double[N_interp*N_theta];

  //Fill neg_shear
  for(int iii=0;iii<N_interp*N_theta;iii++){
    temp_1[iii]=equil->f_psi[iii]*equil->jacob[iii]*equil->g_phph[iii];
  }

  deriv_1d(temp_2,temp_1,equil->rad_interp,N_interp,N_theta,true,deriv_order);
  for(int iii=0;iii<N_interp;iii++){for(int jjj=0;jjj<N_theta;jjj++){temp_2[iii*N_theta+jjj]*=equil->drad_dpsi[iii];}}

  for(int iii=0;iii<N_interp*N_theta;iii++){
    equil->neg_shear[iii]=temp_2[iii]/equil->jacob[iii];
    temp_1[iii]*=(equil->g_pt[iii]/equil->g_pp[iii]);
  }

  deriv_ang(temp_2,temp_1,equil->theta_grid,N_theta,N_interp,false);
  for(int iii=0;iii<N_interp*N_theta;iii++){
    equil->neg_shear[iii]+=temp_2[iii]/equil->jacob[iii];
  }

  //Fill curv_psi
  for(int iii=0;iii<N_interp*N_theta;iii++){
    temp_1[iii]=equil->pres[iii];
  }
  
  deriv_1d(temp_2,temp_1,equil->rad_interp,N_interp,N_theta,true,deriv_order);
  for(int iii=0;iii<N_interp;iii++){for(int jjj=0;jjj<N_theta;jjj++){temp_2[iii*N_theta+jjj]*=equil->drad_dpsi[iii];}}
  
  for(int iii=0;iii<N_interp*N_theta;iii++){
    equil->curv_psi[iii]=2.0*mu_0*temp_2[iii];
    temp_1[iii]=equil->mag_sq[iii];
  }
  
  deriv_1d(temp_2,temp_1,equil->rad_interp,N_interp,N_theta,true,deriv_order);
  for(int iii=0;iii<N_interp;iii++){for(int jjj=0;jjj<N_theta;jjj++){temp_2[iii*N_theta+jjj]*=equil->drad_dpsi[iii];}}
  
  for(int iii=0;iii<N_interp*N_theta;iii++){
    equil->curv_psi[iii]+=temp_2[iii];
  }
  deriv_ang(temp_2,temp_1,equil->theta_grid,N_theta,N_interp,false);
  for(int iii=0;iii<N_interp*N_theta;iii++){
    equil->curv_psi[iii]+=equil->g_pt[iii]*temp_2[iii]/equil->g_pp[iii];
  }

  for(int iii=0;iii<N_interp*N_theta;iii++){
    equil->curv_psi[iii]*=equil->g_pp[iii]/equil->mag_sq[iii];
  }


  //Fill j_dot_b
  for(int iii=0;iii<N_interp*N_theta;iii++){
    temp_1[iii]=equil->f_psi[iii];
  }
  deriv_1d(temp_2,temp_1,equil->rad_interp,N_interp,N_theta,true,deriv_order);
  for(int iii=0;iii<N_interp;iii++){for(int jjj=0;jjj<N_theta;jjj++){temp_2[iii*N_theta+jjj]*=equil->drad_dpsi[iii];}}

  for(int iii=0;iii<N_interp*N_theta;iii++){
    equil->j_dot_b[iii]=temp_2[iii]*equil->g_pp[iii]*equil->g_phph[iii];
    temp_1[iii]=equil->jacob[iii]*equil->g_pt[iii]*equil->g_phph[iii];
  }

  deriv_ang(temp_2,temp_1,equil->theta_grid,N_theta,N_interp,false);

  for(int iii=0;iii<N_interp*N_theta;iii++){
    equil->j_dot_b[iii]-=temp_2[iii]*equil->f_psi[iii]/equil->jacob[iii];
    temp_1[iii]=equil->jacob[iii]*equil->g_pp[iii]*equil->g_phph[iii];
  }

  deriv_1d(temp_2,temp_1,equil->rad_interp,N_interp,N_theta,true,deriv_order);
  for(int iii=0;iii<N_interp;iii++){for(int jjj=0;jjj<N_theta;jjj++){temp_2[iii*N_theta+jjj]*=equil->drad_dpsi[iii];}}

  for(int iii=0;iii<N_interp*N_theta;iii++){
    equil->j_dot_b[iii]-=temp_2[iii]*equil->f_psi[iii]/equil->jacob[iii];
  }

  delete[] temp_1; delete[] temp_2;

  fill_cc(equil,N_interp,N_theta);

  if(geom.Hall_on == 1){
    for(int iii=0;iii<N_interp*N_theta;iii++){ equil->pres[iii] = 0.0 ; }
  }

}


void rad_to_psi(equil_fields *equil,int N_psi,int N_interp)
{
  //rad_var=sqrt(psi) 
  for( int iii=0 ; iii < N_psi ; iii++ ){ equil->psi_grid[iii] = equil->rad_var[iii] * equil->rad_var[iii] ;}

  for( int iii=0 ; iii < N_interp ; iii++ ){ equil->psi_interp[iii] = equil->rad_interp[iii] * equil->rad_interp[iii] ;}

  for( int iii=0 ; iii < N_interp ; iii++ ){ equil->drad_dpsi[iii] = 1.0 / ( 2.0 * equil->rad_interp[iii] ) ;}
}

void fill_cc(equil_fields *equil,int N_interp,int N_theta)
{
  for( int iii=0 ; iii < N_interp * N_theta ; iii++ ){  
    equil->cc_0[iii] = equil->mag_sq[iii] / equil->jacob[iii] ;
    
    equil->cc_1[iii] = equil->g_pt[iii] / equil->g_pp[iii] ;
    
    equil->cc_2[iii] = equil->mag_sq[iii] * equil->curv_psi[iii] / (equil->g_pp[iii] * equil->jacob[iii]) ;
    
    equil->cc_3[iii] = 1.0 ;
    
    equil->cc_4[iii] = ( equil->g_pp[iii] * equil->neg_shear[iii] - equil->j_dot_b[iii] ) / equil->mag_sq[iii] ;    
  }
}

void fill_rad(equil_fields *equil,geom_shape geom, double flux_max, double flux_min)
{
  int N_interp = geom.N_interp ;
  int N_psi = geom.N_psi ;

  //Even spacing
  /*double diff = ( sqrt(flux_max) - sqrt(flux_min) ) / ( N_interp - 1 ) ;
    for(int iii=0 ; iii < N_interp ; iii++){ equil->rad_interp[iii] = sqrt(flux_min) + iii * diff ; }
    for(int iii=0 ; iii < N_psi ; iii++){ equil->rad_var[iii] = equil->psi_interp[ iii * (geom.num_quad + 1) ] ; }*/

  //Equally spaced grid
  double diff = ( sqrt(flux_max) - sqrt(flux_min) ) / ( N_psi - 1 ) ;
  for(int iii=0 ; iii < N_psi ; iii++){
    equil->rad_var[iii] = sqrt(flux_min) + iii * diff ;
    equil->rad_interp[ iii * (geom.num_quad + 1) ] = equil->rad_var[iii] ;
  }

  const double* gqNeval;
  const double* gqNdummy;

  assign_gq(gqNeval,gqNdummy,geom.num_quad);

  //Gaussian quadrature spacing 
  for(int iii=0 ; iii < N_psi - 1 ; iii++){
    for(int jjj=0 ; jjj < geom.num_quad ; jjj++){
      equil->rad_interp[ iii * (geom.num_quad + 1) + jjj + 1 ] = 0.5 * ( equil->rad_var[iii+1] + equil->rad_var[iii] + ( equil->rad_var[iii+1] - equil->rad_var[iii] ) * gqNeval[jjj] )    ;
    }}
}

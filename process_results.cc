#include "constants.h"

void convert_to_mag(std::complex<double> b_par[],std::complex<double> b_perp[],std::complex<double> b_wedge[],std::complex<double> xi_perp[],std::complex<double> xi_wedge[],equil_fields eq,geom_shape geom,int num_sols,std::string shape_order)
{
  std::string deriv_order="Quadratic";
  int N_psi=geom.N_psi;
  int N_theta=geom.N_theta;
  std::complex<double> imag_unit=1.0i;
  double tor_mod=static_cast<double>(geom.tor_mod);

  std::string value;
  read_in("analytical_type",value);
	
  if(value=="cylinder_theta" || value=="cylinder_screw"){tor_mod=static_cast<double>(geom.m_min);}

  
  //if(shape_order=="NHLC")
    {
      std::complex<double> *perp_temp=new std::complex<double>[N_psi*N_theta];
      std::complex<double> *wedge_temp=new std::complex<double>[N_psi*N_theta];
      std::complex<double> *workspace_1=new std::complex<double>[N_psi*N_theta];
      std::complex<double> *workspace_2=new std::complex<double>[N_psi*N_theta];
      double *workspace_doub=new double[N_psi*N_theta];
      
      for(int iii=0;iii<num_sols;iii++){
	
	for(int jjj=0;jjj<N_psi;jjj++){
	  for(int kkk=0;kkk<N_theta;kkk++){ //Normalise to regular \xi_\perp and \xi_\wedge. 
	    perp_temp[jjj*N_theta+kkk]=xi_perp[iii*N_psi*N_theta+jjj*N_theta+kkk] / eq.jacob[jjj*N_theta+kkk] ;
	    wedge_temp[jjj*N_theta+kkk]=xi_wedge[iii*N_psi*N_theta+jjj*N_theta+kkk] / sqrt( eq.g_pp[jjj*N_theta+kkk] / eq.mag_sq[jjj*N_theta+kkk] ) ;
	  }}

	/******************************************************************************************************************************************************************************************/
	//b-perp
	  
	deriv_ang(workspace_1,perp_temp,eq.theta_grid,N_theta,N_psi,false);

	for(int jjj=0;jjj<N_psi;jjj++){
	  for(int kkk=0;kkk<N_theta;kkk++){
	    b_perp[iii*N_psi*N_theta+jjj*N_theta+kkk]=eq.f_psi[jjj*N_theta+kkk]*eq.g_phph[jjj*N_theta+kkk]*imag_unit*tor_mod*perp_temp[jjj*N_theta+kkk]+(1.0/eq.jacob[jjj*N_theta+kkk])*workspace_1[jjj*N_theta+kkk];
	  }}

	/******************************************************************************************************************************************************************************************/
	//b-wedge
	
	deriv_ang(workspace_1,wedge_temp,eq.theta_grid,N_theta,N_psi,false);

	for(int jjj=0;jjj<N_psi;jjj++){
	  for(int kkk=0;kkk<N_theta;kkk++){
	    b_wedge[iii*N_psi*N_theta+jjj*N_theta+kkk]=(eq.g_pp[jjj*N_theta+kkk]/eq.mag_sq[jjj*N_theta+kkk])*(eq.f_psi[jjj*N_theta+kkk]*eq.g_phph[jjj*N_theta+kkk]*imag_unit*tor_mod*wedge_temp[jjj*N_theta+kkk]+(1.0/eq.jacob[jjj*N_theta+kkk])*workspace_1[jjj*N_theta+kkk]-eq.neg_shear[jjj*N_theta+kkk]*perp_temp[jjj*N_theta+kkk]);
	  }}

	/******************************************************************************************************************************************************************************************/
	//b-parallel

	//-(div B^2) dot \vec{xi_perp}
	deriv_1d(workspace_doub,eq.mag_sq,eq.psi_grid,N_psi,N_theta,true,deriv_order);
	for(int jjj=0;jjj<N_psi;jjj++){for(int kkk=0;kkk<N_theta;kkk++){b_par[iii*N_psi*N_theta+jjj*N_theta+kkk]=-workspace_doub[jjj*N_theta+kkk];}}
	deriv_ang(workspace_doub,eq.mag_sq,eq.theta_grid,N_theta,N_psi,false);
	for(int jjj=0;jjj<N_psi;jjj++){for(int kkk=0;kkk<N_theta;kkk++){b_par[iii*N_psi*N_theta+jjj*N_theta+kkk]-=(eq.g_pt[jjj*N_theta+kkk]/eq.g_pp[jjj*N_theta+kkk])*workspace_doub[jjj*N_theta+kkk];}}
	for(int jjj=0;jjj<N_psi;jjj++){for(int kkk=0;kkk<N_theta;kkk++){b_par[iii*N_psi*N_theta+jjj*N_theta+kkk]*=perp_temp[jjj*N_theta+kkk];}}

	//-B^2 (div dot \vec{xi_perp})
	for(int jjj=0;jjj<N_psi;jjj++){for(int kkk=0;kkk<N_theta;kkk++){workspace_1[jjj*N_theta+kkk]=eq.jacob[jjj*N_theta+kkk]*perp_temp[jjj*N_theta+kkk];}}
	deriv_1d(workspace_2,workspace_1,eq.psi_grid,N_psi,N_theta,true,deriv_order);
	for(int jjj=0;jjj<N_psi;jjj++){for(int kkk=0;kkk<N_theta;kkk++){b_par[iii*N_psi*N_theta+jjj*N_theta+kkk]-=(eq.mag_sq[jjj*N_theta+kkk]/eq.jacob[jjj*N_theta+kkk])*workspace_2[jjj*N_theta+kkk];}}

	for(int jjj=0;jjj<N_psi;jjj++){for(int kkk=0;kkk<N_theta;kkk++){workspace_1[jjj*N_theta+kkk]=(eq.g_pt[jjj*N_theta+kkk]/eq.g_pp[jjj*N_theta+kkk])*eq.jacob[jjj*N_theta+kkk]*perp_temp[jjj*N_theta+kkk];}}
	deriv_ang(workspace_2,workspace_1,eq.theta_grid,N_theta,N_psi,false);
	for(int jjj=0;jjj<N_psi;jjj++){for(int kkk=0;kkk<N_theta;kkk++){b_par[iii*N_psi*N_theta+jjj*N_theta+kkk]-=(eq.mag_sq[jjj*N_theta+kkk]/eq.jacob[jjj*N_theta+kkk])*workspace_2[jjj*N_theta+kkk];}}

	//-div dot (B^2 \vec{xi_wedge})
	for(int jjj=0;jjj<N_psi;jjj++){for(int kkk=0;kkk<N_theta;kkk++){b_par[iii*N_psi*N_theta+jjj*N_theta+kkk]-=eq.g_pp[jjj*N_theta+kkk]*eq.g_phph[jjj*N_theta+kkk]*imag_unit*tor_mod*wedge_temp[jjj*N_theta+kkk];}}

	deriv_ang(workspace_1,wedge_temp,eq.theta_grid,N_theta,N_psi,false);
	for(int jjj=0;jjj<N_psi;jjj++){for(int kkk=0;kkk<N_theta;kkk++){b_par[iii*N_psi*N_theta+jjj*N_theta+kkk]+=(eq.f_psi[jjj*N_theta+kkk]/eq.jacob[jjj*N_theta+kkk])*workspace_1[jjj*N_theta+kkk];}}

	//-mu j dot T xi_perp	  
	if(value=="cylinder_theta" || value=="cylinder_screw"){}
	else{
	  for(int jjj=0;jjj<N_psi;jjj++){for(int kkk=0;kkk<N_theta;kkk++){workspace_1[jjj*N_theta+kkk]=(1.0-eq.jacob[jjj*N_theta+kkk]*eq.jacob[jjj*N_theta+kkk]*eq.g_pp[jjj*N_theta+kkk]*eq.g_phph[jjj*N_theta+kkk]*eq.g_tt[jjj*N_theta+kkk])/(eq.jacob[jjj*N_theta+kkk]*eq.g_pt[jjj*N_theta+kkk]);}}
	  deriv_ang(workspace_2,workspace_1,eq.theta_grid,N_theta,N_psi,false);
	  for(int jjj=0;jjj<N_psi;jjj++){for(int kkk=0;kkk<N_theta;kkk++){b_par[iii*N_psi*N_theta+jjj*N_theta+kkk]-=perp_temp[jjj*N_theta+kkk]*workspace_2[jjj*N_theta+kkk]/eq.jacob[jjj*N_theta+kkk];}}
	}

	for(int jjj=0;jjj<N_psi;jjj++){for(int kkk=0;kkk<N_theta;kkk++){workspace_1[jjj*N_theta+kkk]=eq.jacob[jjj*N_theta+kkk]*eq.g_pp[jjj*N_theta+kkk]*eq.g_phph[jjj*N_theta+kkk];}}
	deriv_1d(workspace_2,workspace_1,eq.psi_grid,N_psi,N_theta,true,deriv_order);
	for(int jjj=0;jjj<N_psi;jjj++){for(int kkk=0;kkk<N_theta;kkk++){b_par[iii*N_psi*N_theta+jjj*N_theta+kkk]+=perp_temp[jjj*N_theta+kkk]*workspace_2[jjj*N_theta+kkk]/eq.jacob[jjj*N_theta+kkk];}}

	deriv_1d(workspace_doub,eq.f_psi,eq.psi_grid,N_psi,N_theta,true,deriv_order);
	for(int jjj=0;jjj<N_psi;jjj++){for(int kkk=0;kkk<N_theta;kkk++){b_par[iii*N_psi*N_theta+jjj*N_theta+kkk]+=eq.f_psi[jjj*N_theta+kkk]*eq.g_phph[jjj*N_theta+kkk]*workspace_doub[jjj*N_theta+kkk]*perp_temp[jjj*N_theta+kkk];}}

      }
      
      delete[] perp_temp;delete[] wedge_temp;delete[] workspace_1;delete[] workspace_2;delete[] workspace_doub;
    }
}



int find_pol_mode(PetscReal polmode_norm[],int m_range,int m_min,Vec eigenvectors[],int particular_val,int dim_locs[])
{
  PetscErrorCode ierr;
  const PetscScalar *tmp_arr;
  
  ierr = VecGetArrayRead(eigenvectors[particular_val],&tmp_arr);CHKERRQ(ierr);
  
  for(int jjj=0;jjj<m_range;jjj++)
    {
      polmode_norm[jjj]=0.0;
    }
	  
  int m_count=m_min;
  int k_offset=0;
  int k_range=0;
  
  for(int jjj=0;jjj<m_range;jjj++)
    {
      //Iterate through one "m-value" of the eigenvalue
      if(m_count==0){k_range=dim_locs[0];}
      else if(std::abs(m_count)==1 || std::abs(m_count)==2){k_range=dim_locs[1];}
      else{k_range=dim_locs[2];}

      for(int kkk=0;kkk<k_range;kkk++)
	{
	  polmode_norm[jjj]+=PetscRealPart(tmp_arr[k_offset+kkk])*PetscRealPart(tmp_arr[k_offset+kkk])+PetscImaginaryPart(tmp_arr[k_offset+kkk])*PetscImaginaryPart(tmp_arr[k_offset+kkk]);
	}
      k_offset+=k_range;
      m_count+=1; //Check this properly as added without testing, but makes sense!!
    }

  int max_index=m_min;
  for(int jjj=0;jjj<m_range;jjj++)
    {
      if(polmode_norm[jjj]>polmode_norm[(max_index-m_min)]){max_index=m_min+jjj;}
    }

  return max_index;
}


int find_rad_mode(double xi_perp[],int length,int eigenmode_number)
{
  int output=0;
  /*for(int iii=0;iii<length-1;iii++)
    {
      if((xi_perp[eigenmode_number*length+iii]/xi_perp[eigenmode_number*length+iii+1])<0.0){output+=1;}
      }*/
  
  double max=0.0;
  for(int iii=0;iii<length;iii++)
    {
      if(std::abs(xi_perp[eigenmode_number*length+iii])>max){max=std::abs(xi_perp[eigenmode_number*length+iii]);}
    }

  for(int iii=0;iii<length-1;iii++)
    {
      if((xi_perp[eigenmode_number*length+iii]/xi_perp[eigenmode_number*length+iii+1])<0.0){output+=1;}
    }
  
  int local_zeros[output+2];
  int counter=1;
  for(int iii=0;iii<length-1;iii++)
    {
      if((xi_perp[eigenmode_number*length+iii]/xi_perp[eigenmode_number*length+iii+1])<0.0){local_zeros[counter]=iii; counter++;}
    }
  local_zeros[0]=0;
  local_zeros[output+1]=length;

  int minus_counter=0;
  for(int iii=0;iii<output+1;iii++)
    {
      double local_max=0.0;
      for(int jjj=local_zeros[iii];jjj<local_zeros[iii+1];jjj++)
	{
	  if(std::abs(xi_perp[eigenmode_number*length+jjj])>local_max){local_max=std::abs(xi_perp[eigenmode_number*length+jjj]);}
	}
      if(local_max<0.001*max){minus_counter+=1;}
    }
  
    output-=minus_counter;

    return (output+1);
  /*


  double max=0.0;
  for(int iii=0;iii<length;iii++)
    {
      if(std::abs(xi_perp[eigenmode_number*length+iii])>max){max=std::abs(xi_perp[eigenmode_number*length+iii]);}
    }

  int counter=0;
  for(int iii=0;iii<length;iii++)
    {
      if(std::abs(xi_perp[eigenmode_number*length+iii])>std::abs(xi_perp[eigenmode_number*length+iii+1]) && std::abs(xi_perp[eigenmode_number*length+iii])>std::abs(xi_perp[eigenmode_number*length+iii-1])){counter++;}
    }
  double local_maxes[counter]={};
  counter=0;
  for(int iii=0;iii<length;iii++)
    {
      if(std::abs(xi_perp[eigenmode_number*length+iii])>std::abs(xi_perp[eigenmode_number*length+iii+1]) && std::abs(xi_perp[eigenmode_number*length+iii])>std::abs(xi_perp[eigenmode_number*length+iii-1])){local_maxes[counter]=std::abs(xi_perp[eigenmode_number*length+iii]);counter++;}
    }

  for(int iii=0;iii<counter;iii++)
    {
      if(local_maxes[iii]>0.01*max){output++;}
    }

  if(output==1)
    {
      std::cout << std::endl;
      for(int iii=0;iii<counter;iii++)
	{
	  std::cout << local_maxes[iii] << std::endl;
	}
      std::cout << std::endl;
    }
  
  
    return (output);*/
}




int find_mode_number(double xi_perp[],int length,int eigenmode_number)
{
  int output=0;
  /*for(int iii=0;iii<length-1;iii++)
    {
      if((xi_perp[eigenmode_number*length+iii]/xi_perp[eigenmode_number*length+iii+1])<0.0){output+=1;}
      }*/
  
  double max=0.0;
  for(int iii=0;iii<length;iii++)
    {
      if(std::abs(xi_perp[eigenmode_number*length+iii])>max){max=std::abs(xi_perp[eigenmode_number*length+iii]);}
    }

  for(int iii=0;iii<length-1;iii++)
    {
      if((xi_perp[eigenmode_number*length+iii]/xi_perp[eigenmode_number*length+iii+1])<0.0){output+=1;}
    }
  
  int local_zeros[output+2];
  int counter=1;
  for(int iii=0;iii<length-1;iii++)
    {
      if((xi_perp[eigenmode_number*length+iii]/xi_perp[eigenmode_number*length+iii+1])<0.0){local_zeros[counter]=iii; counter++;}
    }
  local_zeros[0]=0;
  local_zeros[output+1]=length;

  int minus_counter=0;
  for(int iii=0;iii<output+1;iii++)
    {
      double local_max=0.0;
      for(int jjj=local_zeros[iii];jjj<local_zeros[iii+1];jjj++)
	{
	  if(std::abs(xi_perp[eigenmode_number*length+jjj])>local_max){local_max=std::abs(xi_perp[eigenmode_number*length+jjj]);}
	}
      if(local_max<0.001*max){minus_counter+=1;}
    }
  
    output-=minus_counter;

    return (output+1);
  /*


  double max=0.0;
  for(int iii=0;iii<length;iii++)
    {
      if(std::abs(xi_perp[eigenmode_number*length+iii])>max){max=std::abs(xi_perp[eigenmode_number*length+iii]);}
    }

  int counter=0;
  for(int iii=0;iii<length;iii++)
    {
      if(std::abs(xi_perp[eigenmode_number*length+iii])>std::abs(xi_perp[eigenmode_number*length+iii+1]) && std::abs(xi_perp[eigenmode_number*length+iii])>std::abs(xi_perp[eigenmode_number*length+iii-1])){counter++;}
    }
  double local_maxes[counter]={};
  counter=0;
  for(int iii=0;iii<length;iii++)
    {
      if(std::abs(xi_perp[eigenmode_number*length+iii])>std::abs(xi_perp[eigenmode_number*length+iii+1]) && std::abs(xi_perp[eigenmode_number*length+iii])>std::abs(xi_perp[eigenmode_number*length+iii-1])){local_maxes[counter]=std::abs(xi_perp[eigenmode_number*length+iii]);counter++;}
    }

  for(int iii=0;iii<counter;iii++)
    {
      if(local_maxes[iii]>0.01*max){output++;}
    }

  if(output==1)
    {
      std::cout << std::endl;
      for(int iii=0;iii<counter;iii++)
	{
	  std::cout << local_maxes[iii] << std::endl;
	}
      std::cout << std::endl;
    }
  
  
    return (output);*/
}

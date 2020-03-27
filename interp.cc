#include "constants.h"

void interpolate_radial(PetscScalar grid_values[],PetscScalar halfway_values[],int gridpoints,int nconv,std::string order)
{
  if(order=="linear")
    {
      for(int jjj=0;jjj<nconv;jjj++)
	{
	  for(int iii=0;iii<gridpoints-1;iii++)
	    {
	      halfway_values[jjj*(gridpoints-1)+iii]=0.5*(grid_values[jjj*gridpoints+iii]+grid_values[jjj*gridpoints+iii+1]);
	    }
	}
    }
}

void interpolate_bilinear(double r_grid[],double r_grid_polate[],int N_r,double theta_grid[],double theta_grid_polate[],int N_theta,PetscScalar orig_vals[],PetscScalar polate_vals[],int nconv)
{
  for(int iii=0;iii<nconv;iii++)
    {
      for(int jjj=0;jjj<N_r-1;jjj++)
	{
	  for(int kkk=0;kkk<N_theta-1;kkk++)
	    {
	      PetscScalar val1=(r_grid[jjj+1]-r_grid_polate[jjj])*(theta_grid[kkk+1]-theta_grid_polate[kkk])*orig_vals[iii*N_r*N_theta+jjj*N_theta+kkk];
	      PetscScalar val2=(r_grid_polate[jjj]-r_grid[jjj])*(theta_grid[kkk+1]-theta_grid_polate[kkk])*orig_vals[iii*N_r*N_theta+(jjj+1)*N_theta+kkk];
	      PetscScalar val3=(r_grid[jjj+1]-r_grid_polate[jjj])*(theta_grid_polate[kkk]-theta_grid[kkk])*orig_vals[iii*N_r*N_theta+jjj*N_theta+kkk+1];
	      PetscScalar val4=(r_grid_polate[jjj]-r_grid[jjj])*(theta_grid_polate[kkk]-theta_grid[kkk])*orig_vals[iii*N_r*N_theta+(jjj+1)*N_theta+kkk+1];
	      polate_vals[iii*(N_r-1)*(N_theta-1)+jjj*(N_theta-1)+kkk]=(1.0/((theta_grid[kkk+1]-theta_grid[kkk])*(r_grid[jjj+1]-r_grid[jjj])))*(val1+val2+val3+val4);
	    }
	}
    }
}

std::complex<double> interp_lagrange(std::complex<double> grid_func[],double grid[],double int_point,int count_lower,int num_order)
{
  std::complex<double> result=0.0;
  std::complex<double> result_temp=1.0;
  
  //Num_order -> 1 = linear, 2 = quadratic etc
  if(num_order%2==0)
    {
      int size_either_side=num_order/2;
      for(int iii=count_lower-size_either_side;iii<=count_lower+size_either_side;iii++)
	{
	  for(int jjj=count_lower-size_either_side;jjj<=count_lower+size_either_side;jjj++)
	    {
	      if(jjj==iii){}
	      else{result_temp*=(int_point-grid[jjj])/(grid[iii]-grid[jjj]);} //
	    }
	  result_temp*=grid_func[iii];
	  result+=result_temp;
	  result_temp=1.0;
	}
    }
  else
    {
      int size_either_side=(num_order-1)/2;
      for(int iii=count_lower-size_either_side;iii<=count_lower+size_either_side+1;iii++)
	{
	  for(int jjj=count_lower-size_either_side;jjj<=count_lower+size_either_side+1;jjj++)
	    {
	      if(jjj==iii){}
	      else{result_temp*=(int_point-grid[jjj])/(grid[iii]-grid[jjj]);}
	    }
	  result_temp*=grid_func[iii];
	  result+=result_temp;
	  result_temp=1.0;
	}
    }

  return result;
}

std::complex<double> interpolate_1d(std::complex<double> grid_func[],double grid[],int N_psi,double int_point,std::string order)
{
  std::complex<double> result;
  int count_lower;
  if(int_point==grid[0]){return grid_func[0];}
  for(int iii=0;iii<N_psi-1;iii++)
    {
      if(int_point==grid[iii+1]){return grid_func[iii+1];}
      else if(grid[iii]<int_point && int_point<grid[iii+1]){count_lower=iii; iii=N_psi;} //Replace with more efficient search at some point - say binary
    }



  if(order=="Linear_pol")
    {
      assert(N_psi>=2);
      result=interp_lagrange(grid_func,grid,int_point,count_lower,1);
    }
  else if(order=="Quadratic_pol")
    {
      assert(N_psi>=3);
      if(count_lower==0){count_lower=1;}
      result=interp_lagrange(grid_func,grid,int_point,count_lower,2);
    }
  else if(order=="Cubic_pol")
    {
      assert(N_psi>=4);
      if(count_lower==0){count_lower=1;}
      if(count_lower==N_psi-2){count_lower=N_psi-3;}

      result=interp_lagrange(grid_func,grid,int_point,count_lower,3);
    }
  else if(order=="Quartic_pol")
    {
      assert(N_psi>=5);
      if(count_lower==0 || count_lower==1){count_lower=2;}
      if(count_lower==N_psi-2){count_lower=N_psi-3;}

      result=interp_lagrange(grid_func,grid,int_point,count_lower,4);
    }
  else if(order=="All")
  {
    if(N_psi%2==0){count_lower=(N_psi/2)-1;}
    else{count_lower=(N_psi-1)/2;}

    result=interp_lagrange(grid_func,grid,int_point,count_lower,N_psi-1);
  }
  else{std::cout << "Interpolation order not valid" << std::endl;}



  return result;
}

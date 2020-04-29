#include "constants.h"

void deriv_1d(double deriv_grid[],double dep_var[],double ind_var[],int num_points,std::string order)
{
  clean_grid(deriv_grid,num_points);
  
  double *diff_grid = new double[num_points-1];
  for(int iii=0;iii<num_points-1;iii++){diff_grid[iii]=ind_var[iii+1]-ind_var[iii];}

  if(order=="Linear")
    {
      assert(num_points>=2);

      //Boundaries
      deriv_grid[0]=(dep_var[1]-dep_var[0])/diff_grid[0];
      deriv_grid[(num_points-1)]=(dep_var[num_points-1]-dep_var[num_points-2])/diff_grid[num_points-2];

      //Interior
      for(int iii=1;iii<num_points-1;iii++)
	{
	  deriv_grid[iii]=(dep_var[iii+1]-dep_var[iii-1])/(diff_grid[iii]+diff_grid[iii-1]);
	}

    }  
  else if(order=="Quadratic")
    {
      assert(num_points>=3);

      //Boundaries
      deriv_grid[0]=(-diff_grid[1]*(2.0*diff_grid[0]+diff_grid[1])*dep_var[0]+(diff_grid[0]+diff_grid[1])*(diff_grid[0]+diff_grid[1])*dep_var[1]-diff_grid[0]*diff_grid[0]*dep_var[2])/(diff_grid[0]*diff_grid[1]*(diff_grid[0]+diff_grid[1]));
      deriv_grid[num_points-1]=((2.0*diff_grid[num_points-2]+diff_grid[num_points-3])*diff_grid[num_points-3]*dep_var[num_points-1]-(diff_grid[num_points-2]+diff_grid[num_points-3])*(diff_grid[num_points-2]+diff_grid[num_points-3])*dep_var[num_points-2]+diff_grid[num_points-2]*diff_grid[num_points-2]*dep_var[num_points-3])/(diff_grid[num_points-2]*diff_grid[num_points-3]*(diff_grid[num_points-2]+diff_grid[num_points-3]));

      //Interior
      for(int iii=1;iii<num_points-1;iii++)
	{
	  deriv_grid[iii]=(diff_grid[iii-1]*diff_grid[iii-1]*dep_var[iii+1]+(diff_grid[iii]+diff_grid[iii-1])*(diff_grid[iii]-diff_grid[iii-1])*dep_var[iii]-diff_grid[iii]*diff_grid[iii]*dep_var[iii-1])/(diff_grid[iii]*diff_grid[iii-1]*(diff_grid[iii]+diff_grid[iii-1]));
	}
    }
  else if(order=="Cubic")
    {
      assert(num_points>=4);

      //Boundaries
      deriv_grid[0]=ffd_cb(dep_var,diff_grid,0);
      deriv_grid[1]=ffd_cb(dep_var,diff_grid,1);
      deriv_grid[2]=fd_cb(dep_var,diff_grid,2);
      deriv_grid[num_points-2]=bd_cb(dep_var,diff_grid,num_points-2);
      deriv_grid[num_points-1]=bbd_cb(dep_var,diff_grid,num_points-1);

      //Interior
      for(int iii=3;iii<num_points-2;iii++)
	{
	  deriv_grid[iii]=cd_cb(dep_var,diff_grid,iii);
	}
    }
  else{std::cout << "That is not a valid order" << std::endl;}

  delete[] diff_grid;
}


void deriv_1d(double deriv_grid[],double dep_var[],double ind_var[],int num_ind,int num_ignor,bool ind_first,std::string order)
{
  clean_grid(deriv_grid,num_ind*num_ignor);
  
  double *dep_temp = new double[num_ind];
  double *deriv_temp = new double[num_ind];
  
  if(ind_first==true)
    {
      for(int iii=0;iii<num_ignor;iii++){
	for(int jjj=0;jjj<num_ind;jjj++){dep_temp[jjj]=dep_var[jjj*num_ignor+iii];}
	deriv_1d(deriv_temp,dep_temp,ind_var,num_ind,order);
	for(int jjj=0;jjj<num_ind;jjj++){deriv_grid[jjj*num_ignor+iii]=deriv_temp[jjj];}
      }	
    }
  else
    {
      for(int iii=0;iii<num_ignor;iii++){
	for(int jjj=0;jjj<num_ind;jjj++){dep_temp[jjj]=dep_var[iii*num_ind+jjj];}
	deriv_1d(deriv_temp,dep_temp,ind_var,num_ind,order);
	for(int jjj=0;jjj<num_ind;jjj++){deriv_grid[iii*num_ind+jjj]=deriv_temp[jjj];}
      }	
    }

  delete[] dep_temp;
  delete[] deriv_temp;
}



void deriv_1d(std::complex<double> deriv_grid[],std::complex<double> dep_var[],double ind_var[],int num_points,std::string order)
{
  clean_grid(deriv_grid,num_points);
  
  double *diff_grid = new double[num_points-1];
  for(int iii=0;iii<num_points-1;iii++){diff_grid[iii]=ind_var[iii+1]-ind_var[iii];}

  if(order=="Linear")
    {
      assert(num_points>=2);

      //Boundaries
      deriv_grid[0]=(dep_var[1]-dep_var[0])/diff_grid[0];
      deriv_grid[(num_points-1)]=(dep_var[num_points-1]-dep_var[num_points-2])/diff_grid[num_points-2];

      //Interior
      for(int iii=1;iii<num_points-1;iii++)
	{
	  deriv_grid[iii]=(dep_var[iii+1]-dep_var[iii-1])/(diff_grid[iii]+diff_grid[iii-1]);
	}

    }  
  else if(order=="Quadratic")
    {
      assert(num_points>=3);

      //Boundaries
      deriv_grid[0]=(-diff_grid[1]*(2.0*diff_grid[0]+diff_grid[1])*dep_var[0]+(diff_grid[0]+diff_grid[1])*(diff_grid[0]+diff_grid[1])*dep_var[1]-diff_grid[0]*diff_grid[0]*dep_var[2])/(diff_grid[0]*diff_grid[1]*(diff_grid[0]+diff_grid[1]));
      deriv_grid[num_points-1]=((2.0*diff_grid[num_points-2]+diff_grid[num_points-3])*diff_grid[num_points-3]*dep_var[num_points-1]-(diff_grid[num_points-2]+diff_grid[num_points-3])*(diff_grid[num_points-2]+diff_grid[num_points-3])*dep_var[num_points-2]+diff_grid[num_points-2]*diff_grid[num_points-2]*dep_var[num_points-3])/(diff_grid[num_points-2]*diff_grid[num_points-3]*(diff_grid[num_points-2]+diff_grid[num_points-3]));

      //Interior
      for(int iii=1;iii<num_points-1;iii++)
	{
	  deriv_grid[iii]=(diff_grid[iii-1]*diff_grid[iii-1]*dep_var[iii+1]+(diff_grid[iii]+diff_grid[iii-1])*(diff_grid[iii]-diff_grid[iii-1])*dep_var[iii]-diff_grid[iii]*diff_grid[iii]*dep_var[iii-1])/(diff_grid[iii]*diff_grid[iii-1]*(diff_grid[iii]+diff_grid[iii-1]));
	}
    }
  else{std::cout << "That is not a valid order" << std::endl;}

  delete[] diff_grid;
}



void deriv_1d(std::complex<double> deriv_grid[],std::complex<double> dep_var[],double ind_var[],int num_ind,int num_ignor,bool ind_first,std::string order)
{
  clean_grid(deriv_grid,num_ind*num_ignor);
  
  std::complex<double> *dep_temp = new std::complex<double>[num_ind];
  std::complex<double> *deriv_temp = new std::complex<double>[num_ind];
  
  if(ind_first==true)
    {
      for(int iii=0;iii<num_ignor;iii++){
	for(int jjj=0;jjj<num_ind;jjj++){dep_temp[jjj]=dep_var[jjj*num_ignor+iii];}
	deriv_1d(deriv_temp,dep_temp,ind_var,num_ind,order);
	for(int jjj=0;jjj<num_ind;jjj++){deriv_grid[jjj*num_ignor+iii]=deriv_temp[jjj];}
      }	
    }
  else
    {
      for(int iii=0;iii<num_ignor;iii++){
	for(int jjj=0;jjj<num_ind;jjj++){dep_temp[jjj]=dep_var[iii*num_ind+jjj];}
	deriv_1d(deriv_temp,dep_temp,ind_var,num_ind,order);
	for(int jjj=0;jjj<num_ind;jjj++){deriv_grid[iii*num_ind+jjj]=deriv_temp[jjj];}
      }	
    }

  delete[] dep_temp;
  delete[] deriv_temp;
}


double ffd_cb(double dep_var[],double diff[],int deriv)
{
  double x=diff[deriv];
  double y=diff[deriv+1];
  double z=diff[deriv+2];

  double a=-(3.0*x*x+4.0*x*y+y*y+2.0*x*z+y*z)*y*(y+z)*z;
  double b=(x+y)*(x+y+z)*(x+y)*(x+y+z)*z;
  double c=-x*(x+y+z)*x*(x+y+z)*(y+z);
  double d=x*(x+y)*x*(x+y)*y;
  
  double quot=x*y*z*(x+y)*(y+z)*(x+y+z);
  
  return (1.0/quot)*(a*dep_var[deriv]+b*dep_var[deriv+1]+c*dep_var[deriv+2]+d*dep_var[deriv+3]);
}

double fd_cb(double dep_var[],double diff[],int deriv)
{
  double x=diff[deriv-1];
  double y=diff[deriv];
  double z=diff[deriv+1];

  double a=-y*(y+z)*y*z*(y+z);
  double b=(-2.0*x*y+y*y-x*z+y*z)*z*(x+y)*(x+y+z);
  double c=x*(y+z)*x*(y+z)*(x+y+z);
  double d=-x*y*x*y*(x+y);
  
  double quot=x*y*z*(x+y)*(y+z)*(x+y+z);
  
  return (1.0/quot)*(a*dep_var[deriv-1]+b*dep_var[deriv]+c*dep_var[deriv+1]+d*dep_var[deriv+2]);
}

double cd_cb(double dep_var[],double diff[],int deriv)
{
  double w=diff[deriv-2];
  double x=diff[deriv-1];
  double y=diff[deriv];
  double z=diff[deriv+1];

  double a=(2.0*x*y-y*y+x*z-y*z)*z*(x+y)*(x+y+z);
  double b=(-2.0*w*y-2.0*x*y+y*y-w*z-x*z+y*z)*z*(w+x+y)*(w+x+y+z);
  double c=(-w*x-x*x+w*y+2.0*x*y+w*z+2.0*x*z)*w*(x+y+z)*(w+x+y+z);
  double d=(w*x+x*x-w*y-2.0*x*y)*w*(x+y)*(w+x+y);
  
  double quot=w*z*(x+y)*(w+x+y)*(x+y+z)*(w+x+y+z);
  
  return (1.0/quot)*(a*dep_var[deriv-2]+b*dep_var[deriv-1]+c*dep_var[deriv+1]+d*dep_var[deriv+2]);
}

double bd_cb(double dep_var[],double diff[],int deriv)
{
  double x=diff[deriv-2];
  double y=diff[deriv-1];
  double z=diff[deriv];

  double a=y*z*y*z*(y+z);
  double b=-(x+y)*z*z*(x+y)*(x+y+z);
  double c=(-x*y-y*y+x*z+2.0*y*z)*x*(y+z)*(x+y+z);
  double d=y*(x+y)*x*y*(x+y);
  
  double quot=x*y*z*(x+y)*(y+z)*(x+y+z);
  
  return (1.0/quot)*(a*dep_var[deriv-2]+b*dep_var[deriv-1]+c*dep_var[deriv]+d*dep_var[deriv+1]);
}

double bbd_cb(double dep_var[],double diff[],int deriv)
{
  double x=diff[deriv-3];
  double y=diff[deriv-2];
  double z=diff[deriv-1];

  double a=-z*(y+z)*y*z*(y+z);
  double b=z*(x+y+z)*z*(x+y)*(x+y+z);
  double c=-(y+z)*(x+y+z)*x*(y+z)*(x+y+z);
  double d=(x*y+y*y+2.0*x*z+4.0*y*z+3.0*z*z)*x*y*(x+y);
  
  double quot=x*y*z*(x+y)*(y+z)*(x+y+z);
  
  return (1.0/quot)*(a*dep_var[deriv-3]+b*dep_var[deriv-2]+c*dep_var[deriv-1]+d*dep_var[deriv]);
}

void deriv_ang(double deriv_grid[],double dep_var[],double ind_var[],int num_points)
{
  clean_grid(deriv_grid,num_points);
  
  int size=(num_points/2)+1; if(!(num_points%2==0)){size=((num_points-1)/2)+1;}
    
  fftw_complex *out_r2c;  out_r2c = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * size);
  fftw_plan plan_r2c;
  double *in_r2c = new double[num_points];
 
  
  plan_r2c = fftw_plan_dft_r2c_1d(num_points, in_r2c, out_r2c, FFTW_ESTIMATE);

  for(int iii=0;iii<size;iii++){out_r2c[iii][0]=0.0; out_r2c[iii][1]=0.0;}
  for(int iii=0;iii<num_points;iii++){in_r2c[iii]=dep_var[iii];}
  fftw_execute(plan_r2c);

  fftw_destroy_plan(plan_r2c); 

  //for(int iii=0;iii<size;iii++){std::cout << out[iii][0] << " + i " << out[iii][1] << std::endl;}

  fftw_complex *in_c2r; in_c2r = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * size);
  fftw_plan plan_c2r;
  double *out_c2r = new double[num_points];

  plan_c2r = fftw_plan_dft_c2r_1d(num_points, in_c2r, out_c2r, FFTW_ESTIMATE);

  for(int iii=0;iii<num_points;iii++){out_c2r[iii]=0.0;}
  for(int iii=0;iii<size;iii++){in_c2r[iii][0]=-iii*out_r2c[iii][1]; in_c2r[iii][1]=iii*out_r2c[iii][0];}

  fftw_free(out_r2c);
  
  fftw_execute(plan_c2r);

  for(int iii=0;iii<num_points;iii++){deriv_grid[iii]=out_c2r[iii]/static_cast<double>(num_points);}

  fftw_destroy_plan(plan_c2r);
  fftw_free(in_c2r);
  delete[] in_r2c;
  delete[] out_c2r;
}

void deriv_ang(double deriv_grid[],double dep_var[],double ind_var[],int num_ind,int num_ignor,bool ind_first)
{
  clean_grid(deriv_grid,num_ind*num_ignor);
  
  double *dep_temp = new double[num_ind];
  double *deriv_temp = new double[num_ind];

  if(ind_first==true)
    {
      for(int iii=0;iii<num_ignor;iii++){
	for(int jjj=0;jjj<num_ind;jjj++){dep_temp[jjj]=dep_var[jjj*num_ignor+iii];}
	deriv_ang(deriv_temp,dep_temp,ind_var,num_ind);
	for(int jjj=0;jjj<num_ind;jjj++){deriv_grid[jjj*num_ignor+iii]=deriv_temp[jjj];}
      }	
    }
  else
    {
      for(int iii=0;iii<num_ignor;iii++){
	for(int jjj=0;jjj<num_ind;jjj++){dep_temp[jjj]=dep_var[iii*num_ind+jjj];}
	deriv_ang(deriv_temp,dep_temp,ind_var,num_ind);
	for(int jjj=0;jjj<num_ind;jjj++){deriv_grid[iii*num_ind+jjj]=deriv_temp[jjj];}
      }	
    }

  delete[] dep_temp;
  delete[] deriv_temp;
}


void deriv_ang(std::complex<double> deriv_grid[],std::complex<double> dep_var[],double ind_var[],int num_points)
{
  clean_grid(deriv_grid,num_points);
    
  fftw_complex *out_r2c;  out_r2c = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * num_points);
  fftw_complex *in_r2c;  in_r2c = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * num_points);
  fftw_plan plan_r2c;

  plan_r2c = fftw_plan_dft_1d(num_points, in_r2c, out_r2c,FFTW_FORWARD, FFTW_ESTIMATE);

  for(int iii=0;iii<num_points;iii++){out_r2c[iii][0]=0.0; out_r2c[iii][1]=0.0;}
  for(int iii=0;iii<num_points;iii++){in_r2c[iii][0]=std::real(dep_var[iii]); in_r2c[iii][1]=std::imag(dep_var[iii]);}

  fftw_execute(plan_r2c);
  fftw_destroy_plan(plan_r2c);
  fftw_free(in_r2c);

  fftw_complex *in_c2r; in_c2r = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * num_points);
  fftw_complex *out_c2r; out_c2r = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * num_points);
  fftw_plan plan_c2r;

  plan_c2r = fftw_plan_dft_1d(num_points, in_c2r, out_c2r,FFTW_BACKWARD,FFTW_ESTIMATE);

  for(int iii=0;iii<num_points;iii++){out_c2r[iii][0]=0.0; out_c2r[iii][1]=0.0;}
  if(num_points%2==0){
  for(int iii=0;iii<(num_points/2);iii++){in_c2r[iii][0]=-iii*out_r2c[iii][1]; in_c2r[iii][1]=iii*out_r2c[iii][0];}
  for(int iii=(num_points/2)+1;iii<num_points;iii++){in_c2r[iii][0]=-(iii-num_points)*out_r2c[iii][1]; in_c2r[iii][1]=(iii-num_points)*out_r2c[iii][0];}} //Ignore the ambiguous point
  else{
    for(int iii=0;iii<((num_points-1)/2)+1;iii++){in_c2r[iii][0]=-iii*out_r2c[iii][1]; in_c2r[iii][1]=iii*out_r2c[iii][0];}
    for(int iii=((num_points-1)/2)+1;iii<num_points;iii++){in_c2r[iii][0]=-(iii-num_points)*out_r2c[iii][1]; in_c2r[iii][1]=(iii-num_points)*out_r2c[iii][0];}
  }

  fftw_free(out_r2c);
  
  fftw_execute(plan_c2r);

  for(int iii=0;iii<num_points;iii++){deriv_grid[iii]=(out_c2r[iii][0]+imag_unit*out_c2r[iii][1])/static_cast<double>(num_points);}
  
  fftw_destroy_plan(plan_c2r);
  fftw_free(in_c2r);
  fftw_free(out_c2r);
}

void deriv_ang(std::complex<double> deriv_grid[],std::complex<double> dep_var[],double ind_var[],int num_ind,int num_ignor,bool ind_first)
{
  clean_grid(deriv_grid,num_ind*num_ignor);
  
  std::complex<double> *dep_temp = new std::complex<double>[num_ind];
  std::complex<double> *deriv_temp = new std::complex<double>[num_ind];

  if(ind_first==true)
    {
      for(int iii=0;iii<num_ignor;iii++){
	for(int jjj=0;jjj<num_ind;jjj++){dep_temp[jjj]=dep_var[jjj*num_ignor+iii];}
	deriv_ang(deriv_temp,dep_temp,ind_var,num_ind);
	for(int jjj=0;jjj<num_ind;jjj++){deriv_grid[jjj*num_ignor+iii]=deriv_temp[jjj];}
      }	
    }
  else
    {
      for(int iii=0;iii<num_ignor;iii++){
	for(int jjj=0;jjj<num_ind;jjj++){dep_temp[jjj]=dep_var[iii*num_ind+jjj];}
	deriv_ang(deriv_temp,dep_temp,ind_var,num_ind);
	for(int jjj=0;jjj<num_ind;jjj++){deriv_grid[iii*num_ind+jjj]=deriv_temp[jjj];}
      }	
    }

  delete[] dep_temp;
  delete[] deriv_temp;
}


std::complex<double> dpol_xi(std::complex<double> *funcs_fourier[],bool fg_sym[],int main_mod,int sec_mod,int N_psi,int N_theta,int psi_index)
{
  //Calculates F \partial_\theta ( G \xi) as \sum_{f,g s.t. f+g+sec_mod=main_mod} F_f G_g i(g+sec_mod)  where m_2=sec_mod m_diff=main_mod-sec_mod
  assert(psi_index<N_psi);

  int half_size=(N_theta/2)-1; if(!(N_theta%2==0)){half_size=(N_theta-1)/2;}
  int fourier_size_sym=half_size+1;
  int fourier_size_full=2*half_size+1;
  std::complex<double> result=0.0;
  int m_diff=main_mod-sec_mod;

  int f_index;
  int g_index;

  if(fg_sym[0]==true && fg_sym[1]==true)
    {
      if(m_diff>0)
	{
	  //f modes are positive, g are negative
	  for(int m_count=0;m_count<=half_size-m_diff-1;m_count++)
	    {
	      f_index=psi_index*fourier_size_sym+half_size-m_count;
	      g_index=psi_index*fourier_size_sym+std::abs(m_count-half_size+m_diff);
	  
	      result+=static_cast<double>(sec_mod+m_count-half_size+m_diff)*funcs_fourier[0][f_index]*std::conj(funcs_fourier[1][g_index]);
	    }
	  //f modes are pos, g are pos
	  for(int m_count=half_size-m_diff;m_count<=half_size;m_count++)
	    {
	      if(m_count<0){m_count=0;} //Checks needed for cases when |m_diff|>half_size (not needed for other loops as they cannot happen)
	      if(m_count-half_size+m_diff>(fourier_size_sym-1)){break;}
	      
	      f_index=psi_index*fourier_size_sym+half_size-m_count;
	      g_index=psi_index*fourier_size_sym+m_count-half_size+m_diff;

	      result+=static_cast<double>(sec_mod+m_count-half_size+m_diff)*funcs_fourier[0][f_index]*funcs_fourier[1][g_index];	  
	    }
	  //f modes are neg, g are pos
	  for(int m_count=half_size+1;m_count<=2*half_size-m_diff;m_count++)
	    {
	      f_index=psi_index*fourier_size_sym+std::abs(half_size-m_count);
	      g_index=psi_index*fourier_size_sym+m_count-half_size+m_diff;
	  
	      result+=static_cast<double>(sec_mod+m_count-half_size+m_diff)*std::conj(funcs_fourier[0][f_index])*funcs_fourier[1][g_index];
	    }
	}
      else
	{
	  //f modes are positive, g are negative
	  for(int m_count=0;m_count<=half_size+m_diff;m_count++)
	    {
	      f_index=psi_index*fourier_size_sym+half_size+m_diff-m_count;
	      g_index=psi_index*fourier_size_sym+std::abs(m_count-half_size);

	      result+=static_cast<double>(sec_mod+m_count-half_size)*funcs_fourier[0][f_index]*std::conj(funcs_fourier[1][g_index]);
	    }
	  //f modes are neg, g are neg
	  for(int m_count=half_size+m_diff+1;m_count<=half_size;m_count++)
	    {
	      if(m_count<0){m_count=0;} //Checks needed for cases when |m_diff|>half_size (not needed for other loops as they cannot happen)
	      if(std::abs(half_size+m_diff-m_count)>(fourier_size_sym-1)){break;}
	      
	      f_index=psi_index*fourier_size_sym+std::abs(half_size+m_diff-m_count);
	      g_index=psi_index*fourier_size_sym+std::abs(m_count-half_size);
	  
	      result+=static_cast<double>(sec_mod+m_count-half_size)*std::conj(funcs_fourier[0][f_index])*std::conj(funcs_fourier[1][g_index]);
	    }
	  //f modes are neg, g are pos
	  for(int m_count=half_size+1;m_count<=2*half_size+m_diff;m_count++)
	    {
	      f_index=psi_index*fourier_size_sym+std::abs(half_size+m_diff-m_count);
	      g_index=psi_index*fourier_size_sym+m_count-half_size;
	  
	      result+=static_cast<double>(sec_mod+m_count-half_size)*std::conj(funcs_fourier[0][f_index])*funcs_fourier[1][g_index];
	    }
	}
    }
  else if(fg_sym[0]==true && fg_sym[1]==false)
    {
      if(m_diff>0)
	{
	  //f modes are positive, g are negative
	  for(int m_count=0;m_count<=half_size-m_diff-1;m_count++)
	    {
	      f_index=psi_index*fourier_size_sym+half_size-m_count;
	      g_index=psi_index*fourier_size_full+m_count+m_diff;
	  
	      result+=static_cast<double>(sec_mod+m_count-half_size+m_diff)*funcs_fourier[0][f_index]*funcs_fourier[1][g_index];      
	    }
	  //f modes are pos, g are pos
	  for(int m_count=half_size-m_diff;m_count<=half_size;m_count++)
	    {
	      if(m_count<0){m_count=0;} //Checks needed for cases when |m_diff|>half_size (not needed for other loops as they cannot happen)
	      if(m_count+m_diff>(fourier_size_full-1)){break;}
	      
	      f_index=psi_index*fourier_size_sym+half_size-m_count;
	      g_index=psi_index*fourier_size_full+m_count+m_diff;

	      result+=static_cast<double>(sec_mod+m_count-half_size+m_diff)*funcs_fourier[0][f_index]*funcs_fourier[1][g_index];     
	    }
	  //f modes are neg, g are pos
	  for(int m_count=half_size+1;m_count<=2*half_size-m_diff;m_count++)
	    {
	      f_index=psi_index*fourier_size_sym+std::abs(half_size-m_count);
	      g_index=psi_index*fourier_size_full+m_count+m_diff;
	  
	      result+=static_cast<double>(sec_mod+m_count-half_size+m_diff)*std::conj(funcs_fourier[0][f_index])*funcs_fourier[1][g_index];  
	    }
	}
      else
	{
	  //f modes are positive, g are negative
	  for(int m_count=0;m_count<=half_size+m_diff;m_count++)
	    {
	      f_index=psi_index*fourier_size_sym+half_size+m_diff-m_count;
	      g_index=psi_index*fourier_size_full+m_count;

	      result+=static_cast<double>(sec_mod+m_count-half_size)*funcs_fourier[0][f_index]*funcs_fourier[1][g_index];	      
	    }
	  //f modes are neg, g are neg
	  for(int m_count=half_size+m_diff+1;m_count<=half_size;m_count++)
	    {
	      if(m_count<0){m_count=0;} //Checks needed for cases when |m_diff|>half_size (not needed for other loops as they cannot happen)
	      if(std::abs(half_size+m_diff-m_count)>(fourier_size_sym-1)){break;}
	      
	      f_index=psi_index*fourier_size_sym+std::abs(half_size+m_diff-m_count);
	      g_index=psi_index*fourier_size_full+m_count;
	  
	      result+=static_cast<double>(sec_mod+m_count-half_size)*std::conj(funcs_fourier[0][f_index])*funcs_fourier[1][g_index];      
	    }
	  //f modes are neg, g are pos
	  for(int m_count=half_size+1;m_count<=2*half_size+m_diff;m_count++)
	    {
	      f_index=psi_index*fourier_size_sym+std::abs(half_size+m_diff-m_count);
	      g_index=psi_index*fourier_size_full+m_count;
	  
	      result+=static_cast<double>(sec_mod+m_count-half_size)*std::conj(funcs_fourier[0][f_index])*funcs_fourier[1][g_index];     
	    }
	}
    }
  else if(fg_sym[0]==false && fg_sym[1]==true)
    {
      if(m_diff>0)
	{
	  //f modes are positive, g are negative
	  for(int m_count=0;m_count<=half_size-m_diff-1;m_count++)
	    {
	      f_index=psi_index*fourier_size_full+2*half_size-m_count;
	      g_index=psi_index*fourier_size_sym+std::abs(m_count-half_size+m_diff);
	  
	      result+=static_cast<double>(sec_mod+m_count-half_size+m_diff)*funcs_fourier[0][f_index]*std::conj(funcs_fourier[1][g_index]);  
	    }
	  //f modes are pos, g are pos
	  for(int m_count=half_size-m_diff;m_count<=half_size;m_count++)
	    {
	      if(m_count<0){m_count=0;} //Checks needed for cases when |m_diff|>half_size (not needed for other loops as they cannot happen)
	      if(m_count-half_size+m_diff>(fourier_size_sym-1)){break;}
	      
	      f_index=psi_index*fourier_size_full+2*half_size-m_count;
	      g_index=psi_index*fourier_size_sym+m_count-half_size+m_diff;

	      result+=static_cast<double>(sec_mod+m_count-half_size+m_diff)*funcs_fourier[0][f_index]*funcs_fourier[1][g_index];      
	    }
	  //f modes are neg, g are pos
	  for(int m_count=half_size+1;m_count<=2*half_size-m_diff;m_count++)
	    {
	      f_index=psi_index*fourier_size_full+2*half_size-m_count;
	      g_index=psi_index*fourier_size_sym+m_count-half_size+m_diff;
	  
	      result+=static_cast<double>(sec_mod+m_count-half_size+m_diff)*funcs_fourier[0][f_index]*funcs_fourier[1][g_index];      
	    }
	}
      else
	{
	  //f modes are positive, g are negative
	  for(int m_count=0;m_count<=half_size+m_diff;m_count++)
	    {
	      f_index=psi_index*fourier_size_full+2*half_size+m_diff-m_count;
	      g_index=psi_index*fourier_size_sym+std::abs(m_count-half_size);

	      result+=static_cast<double>(sec_mod+m_count-half_size)*funcs_fourier[0][f_index]*std::conj(funcs_fourier[1][g_index]);      
	    }
	  //f modes are neg, g are neg
	  for(int m_count=half_size+m_diff+1;m_count<=half_size;m_count++)
	    {
	      if(m_count<0){m_count=0;} //Checks needed for cases when |m_diff|>half_size (not needed for other loops as they cannot happen)
	      if(2*half_size+m_diff-m_count<0){break;}
	      
	      f_index=psi_index*fourier_size_full+2*half_size+m_diff-m_count;
	      g_index=psi_index*fourier_size_sym+std::abs(m_count-half_size);
	  
	      result+=static_cast<double>(sec_mod+m_count-half_size)*funcs_fourier[0][f_index]*std::conj(funcs_fourier[1][g_index]);      
	    }
	  //f modes are neg, g are pos
	  for(int m_count=half_size+1;m_count<=2*half_size+m_diff;m_count++)
	    {
	      f_index=psi_index*fourier_size_full+2*half_size+m_diff-m_count;
	      g_index=psi_index*fourier_size_sym+m_count-half_size;
	  
	      result+=static_cast<double>(sec_mod+m_count-half_size)*funcs_fourier[0][f_index]*funcs_fourier[1][g_index];
	    }
	}
    }
  else if(fg_sym[0]==false && fg_sym[1]==false)
    {
      if(m_diff>0)
	{
	  //f modes are positive, g are negative
	  for(int m_count=0;m_count<=half_size-m_diff-1;m_count++)
	    {
	      f_index=psi_index*fourier_size_full+2*half_size-m_count;
	      g_index=psi_index*fourier_size_full+m_count+m_diff;
	  
	      result+=static_cast<double>(sec_mod+m_count-half_size+m_diff)*funcs_fourier[0][f_index]*funcs_fourier[1][g_index];
	    }
	  //f modes are pos, g are pos
	  for(int m_count=half_size-m_diff;m_count<=half_size;m_count++)
	    {
	      if(m_count<0){m_count=0;} //Checks needed for cases when |m_diff|>half_size (not needed for other loops as they cannot happen)
	      if(m_diff+m_count>(fourier_size_full-1)){break;}
	      
	      f_index=psi_index*fourier_size_full+2*half_size-m_count;
	      g_index=psi_index*fourier_size_full+m_count+m_diff;

	      result+=static_cast<double>(sec_mod+m_count-half_size+m_diff)*funcs_fourier[0][f_index]*funcs_fourier[1][g_index];	  
	    }
	  //f modes are neg, g are pos
	  for(int m_count=half_size+1;m_count<=2*half_size-m_diff;m_count++)
	    {
	      f_index=psi_index*fourier_size_full+2*half_size-m_count;
	      g_index=psi_index*fourier_size_full+m_count+m_diff;
	  
	      result+=static_cast<double>(sec_mod+m_count-half_size+m_diff)*funcs_fourier[0][f_index]*funcs_fourier[1][g_index];
	    }
	}
      else
	{
	  //f modes are positive, g are negative
	  for(int m_count=0;m_count<=half_size+m_diff;m_count++)
	    {
	      f_index=psi_index*fourier_size_full+2*half_size+m_diff-m_count;
	      g_index=psi_index*fourier_size_full+m_count;

	      result+=static_cast<double>(sec_mod+m_count-half_size)*funcs_fourier[0][f_index]*funcs_fourier[1][g_index];
	    }
	  //f modes are neg, g are neg
	  for(int m_count=half_size+m_diff+1;m_count<=half_size;m_count++)
	    {
	      if(m_count<0){m_count=0;} //Checks needed for cases when |m_diff|>half_size (not needed for other loops as they cannot happen)
	      if(2*half_size+m_diff-m_count<0){break;}
	      
	      f_index=psi_index*fourier_size_full+2*half_size+m_diff-m_count;
	      g_index=psi_index*fourier_size_full+m_count;
	  
	      result+=static_cast<double>(sec_mod+m_count-half_size)*funcs_fourier[0][f_index]*funcs_fourier[1][g_index];
	    }
	  //f modes are neg, g are pos
	  for(int m_count=half_size+1;m_count<=2*half_size+m_diff;m_count++)
	    {
	      f_index=psi_index*fourier_size_full+2*half_size+m_diff-m_count;
	      g_index=psi_index*fourier_size_full+m_count;
	  
	      result+=static_cast<double>(sec_mod+m_count-half_size)*funcs_fourier[0][f_index]*funcs_fourier[1][g_index];
	    }
	}
    }

  //Multiply by i
  result*=imag_unit;

  return result;
}


std::complex<double> dpol_sq_xi(std::complex<double> *funcs_fourier[], int main_mod,int sec_mod,int N_psi,int N_theta,int psi_index)
{
  //Calculates f(\theta) \partial_\theta ( g(\theta) \partial_\theta (h(\theta) \xi))
  int half_size=(N_theta/2)-1; if(!(N_theta%2==0)){half_size=(N_theta-1)/2;}
  int fourier_size_sym=half_size+1;
  int fourier_size_full=2*half_size+1;
  
  std::complex<double> *inner_funcs[2]; //Create array for funcs g and h
  inner_funcs[0]=funcs_fourier[1];
  inner_funcs[1]=funcs_fourier[2];

  bool fg_sym[2]={true,true};

  std::complex<double> *first_deriv = new std::complex<double>[N_psi*fourier_size_full]; //This is wasteful as mostly matrix will be zeroes - this is to facilitate use with fdgxi, but should be amended to something more efficient in future versions
  for(int iii=0;iii<N_psi*fourier_size_full;iii++){first_deriv[iii]=0.0;}
  
  for(int m_prime=0;m_prime<fourier_size_full;m_prime++)
    {
      first_deriv[psi_index*fourier_size_full+m_prime]=dpol_xi(inner_funcs,fg_sym,m_prime-half_size+sec_mod,sec_mod,N_psi,N_theta,psi_index);
    }
  
  std::complex<double> *final_funcs[2];
  final_funcs[0]=funcs_fourier[0];
  final_funcs[1]=first_deriv;

  fg_sym[1]=false;
  
  std::complex<double> result=dpol_xi(final_funcs,fg_sym,main_mod,sec_mod,N_psi,N_theta,psi_index); //Here sec_mod is kept separate, so need fdgxi

  delete[] first_deriv;
  
  return result;    
}




void dpar_xi(double *input_funcs[],std::complex<double> result[],bool fg_sym[],geom_shape geom,equil_fields eq,int main_mod,int sec_mod)
{
  //Calculates
  clean_grid(result,geom.N_interp);

  int m_diff=main_mod-sec_mod;
 
  double *temp_1 = new double[geom.N_interp*geom.N_theta];
  for(int iii=0;iii<geom.N_interp*geom.N_theta; iii++){
    temp_1[iii]=input_funcs[0][iii]/eq.jacob[iii];}

  std::complex<double> *fourier_funcs[2];
  std::complex<double> *fourier_1,*fourier_2;
  
  if(fg_sym[0]==false || fg_sym[1]==false){
    //fourier_1 = new std::complex<double> [geom.N_interp*geom.fourier_size_full];
    //calc_fourier_full(temp_1,fourier_1,geom.N_interp,geom.N_theta); fourier_funcs[0] = fourier_1;
    //calc_fourier_full(input_funcs[1],fourier_1,geom.N_interp,geom.N_theta); fourier_funcs[1] = fourier_1;
  }
  else{
    fourier_1 = new std::complex<double> [geom.N_interp*geom.fourier_size_sym];
    fourier_2 = new std::complex<double> [geom.N_interp*geom.fourier_size_sym];
    
    calc_fourier_sym(temp_1,fourier_1,geom.N_interp,geom.N_theta);
    fourier_funcs[0]=fourier_1;
    calc_fourier_sym(input_funcs[1],fourier_2,geom.N_interp,geom.N_theta);
    fourier_funcs[1]=fourier_2;
  }


  for(int iii=0;iii<geom.N_interp;iii++){ result[iii] = dpol_xi(fourier_funcs,fg_sym,main_mod,sec_mod,geom.N_interp,geom.N_theta,iii); }
  

  for(int iii=0;iii<geom.N_interp*geom.N_theta; iii++){
    temp_1[iii]=input_funcs[0][iii]*input_funcs[1][iii]*eq.f_psi[iii]*eq.g_phph[iii];}

  if(fg_sym[0]==false || fg_sym[1]==false){
    //calc_fourier_full(temp_1,fourier_1,geom.N_interp,geom.N_theta);
    //int half_size=geom.fourier_size_sym-1;
    //for(int iii=0;iii<geom.N_interp;iii++){ result[iii] += imag_unit*geom.tor_mod*fourier_1[iii*geom.fourier_size_full+half_size+m_diff]; }
  }
  else{
    calc_fourier_sym(temp_1,fourier_1,geom.N_interp,geom.N_theta);
    if(m_diff>=0){ for(int iii=0;iii<geom.N_interp;iii++){ result[iii] += imag_unit*geom.tor_mod*fourier_1[iii*geom.fourier_size_sym+m_diff]; } }
    else{ for(int iii=0;iii<geom.N_interp;iii++){ result[iii] += imag_unit*geom.tor_mod*std::conj(fourier_1[iii*geom.fourier_size_sym-m_diff]); } }
  }

  delete[] temp_1;
  delete[] fourier_1;
  delete[] fourier_2;
}



void dpar_sq_xi(double *input_funcs[],std::complex<double> result[],geom_shape geom,equil_fields eq,int main_mod,int sec_mod)
{
  clean_grid(result,geom.N_psi);
  
  //Calculates 

  int fourier_size_sym=geom.fourier_size_sym;
  int fourier_size_full=geom.fourier_size_full;
  bool fg_sym[2]={true,true};
  int m_diff=main_mod-sec_mod;

  double *temp_1 = new double[geom.N_interp*geom.N_theta];

  std::complex<double> *fourier_funcs_2[2];
  std::complex<double> *fourier_funcs_3[3];
  std::complex<double> *fourier_1,*fourier_2,*fourier_3;

  fourier_1 = new std::complex<double> [geom.N_interp*geom.fourier_size_sym];
  fourier_2 = new std::complex<double> [geom.N_interp*geom.fourier_size_sym];
  fourier_3 = new std::complex<double> [geom.N_interp*geom.fourier_size_sym];

  /////////
  for( int iii = 0 ; iii < geom.N_interp * geom.N_theta ; iii++ ){temp_1[iii] = input_funcs[0][iii] / eq.jacob[iii];}
  calc_fourier_sym(temp_1,fourier_1,geom.N_interp,geom.N_theta);
  for( int iii = 0 ; iii < geom.N_interp * geom.N_theta ; iii++ ){temp_1[iii] = input_funcs[1][iii] / eq.jacob[iii];}
  calc_fourier_sym(temp_1,fourier_2,geom.N_interp,geom.N_theta);
  calc_fourier_sym(input_funcs[2],fourier_3,geom.N_interp,geom.N_theta);
  fourier_funcs_3[0] = fourier_1; fourier_funcs_3[1] = fourier_2; fourier_funcs_3[2] = fourier_3;
  for(int iii=0;iii<geom.N_interp;iii++){ result[iii] = dpol_sq_xi(fourier_funcs_3,main_mod,sec_mod,geom.N_interp,geom.N_theta,iii); }
  
  ////////
  for( int iii = 0 ; iii < geom.N_interp * geom.N_theta ; iii++ ){temp_1[iii] = input_funcs[1][iii] * input_funcs[2][iii] * eq.f_psi[iii] * eq.g_phph[iii] ;}
  calc_fourier_sym(temp_1,fourier_2,geom.N_interp,geom.N_theta);
  fourier_funcs_2[0] = fourier_1; fourier_funcs_2[1] = fourier_2;
  for(int iii=0;iii<geom.N_interp;iii++){ result[iii] += geom.tor_mod * imag_unit * dpol_xi(fourier_funcs_2,fg_sym,main_mod,sec_mod,geom.N_interp,geom.N_theta,iii); }
    
  ////////
  for( int iii = 0 ; iii < geom.N_interp * geom.N_theta ; iii++ ){temp_1[iii] = input_funcs[0][iii] * input_funcs[1][iii] * eq.f_psi[iii] * eq.g_phph[iii] / eq.jacob[iii] ;}
  calc_fourier_sym(temp_1,fourier_1,geom.N_interp,geom.N_theta);
  fourier_funcs_2[0] = fourier_1; fourier_funcs_2[1] = fourier_3;
  for(int iii=0;iii<geom.N_interp;iii++){ result[iii] += geom.tor_mod * imag_unit * dpol_xi(fourier_funcs_2,fg_sym,main_mod,sec_mod,geom.N_interp,geom.N_theta,iii); }
  
  delete [] fourier_3;
  delete [] fourier_2;
  
  ////////
  for( int iii = 0 ; iii < geom.N_interp * geom.N_theta ; iii++ ){temp_1[iii] = input_funcs[0][iii] * input_funcs[1][iii] * input_funcs[2][iii] * eq.f_psi[iii] * eq.f_psi[iii] * eq.g_phph[iii] * eq.g_phph[iii] ;}
  calc_fourier_sym(temp_1,fourier_1,geom.N_interp,geom.N_theta);
  if(m_diff>=0){ for(int iii=0;iii<geom.N_interp;iii++){ result[iii] -= geom.tor_mod * geom.tor_mod * fourier_1[iii*geom.fourier_size_sym+m_diff]; } }
  else{ for(int iii=0;iii<geom.N_interp;iii++){ result[iii] -= geom.tor_mod * geom.tor_mod * std::conj(fourier_1[iii*geom.fourier_size_sym-m_diff]); } }
  
  delete[] temp_1;
  delete[] fourier_1;
}



void dwedge_xi(double *input_funcs[],std::complex<double> result[],bool fg_sym[],geom_shape geom,equil_fields eq,int main_mod,int sec_mod)
{
  //ASSUMES ALL INPUT FUNCS ARE REAL!
  //Calculates derivative in wedge direction
  clean_grid(result,geom.N_interp);

  int m_diff=main_mod-sec_mod;

  //Poloidal derivative
  double *temp_1 = new double[geom.N_interp*geom.N_theta];
  for(int iii=0;iii<geom.N_interp*geom.N_theta; iii++){
    temp_1[iii]= - input_funcs[0][iii] * eq.f_psi[iii] / (eq.jacob[iii] * eq.g_pp[iii]) ;}

  std::complex<double> *fourier_funcs[2];
  std::complex<double> *fourier_1,*fourier_2;
  
  fourier_1 = new std::complex<double> [geom.N_interp*geom.fourier_size_sym];
  fourier_2 = new std::complex<double> [geom.N_interp*geom.fourier_size_sym];
    
  calc_fourier_sym(temp_1,fourier_1,geom.N_interp,geom.N_theta);
  fourier_funcs[0]=fourier_1;
  calc_fourier_sym(input_funcs[1],fourier_2,geom.N_interp,geom.N_theta);
  fourier_funcs[1]=fourier_2;


  for(int iii=0;iii<geom.N_interp;iii++){ result[iii] = dpol_xi(fourier_funcs,fg_sym,main_mod,sec_mod,geom.N_interp,geom.N_theta,iii) ;}
  

  for(int iii=0;iii<geom.N_interp*geom.N_theta; iii++){
    temp_1[iii]= input_funcs[0][iii] * input_funcs[1][iii] * eq.g_phph[iii] ;}

  calc_fourier_sym(temp_1,fourier_1,geom.N_interp,geom.N_theta);
  if(m_diff>=0){ for(int iii=0;iii<geom.N_interp;iii++){ result[iii] += imag_unit*geom.tor_mod*fourier_1[iii*geom.fourier_size_sym+m_diff]; } }
  else{ for(int iii=0;iii<geom.N_interp;iii++){ result[iii] += imag_unit*geom.tor_mod*std::conj(fourier_1[iii*geom.fourier_size_sym-m_diff]); } }

  delete[] temp_1;
  delete[] fourier_1;
  delete[] fourier_2;
}


void dwedge_sq_xi(double *input_funcs[],std::complex<double> result[],geom_shape geom,equil_fields eq,int main_mod,int sec_mod)
{
  clean_grid(result,geom.N_psi);
  
  //Calculates 

  int fourier_size_sym=geom.fourier_size_sym;
  int fourier_size_full=geom.fourier_size_full;
  bool fg_sym[2]={true,true};
  int m_diff=main_mod-sec_mod;

  double *temp_1 = new double[geom.N_interp*geom.N_theta];

  std::complex<double> *fourier_funcs_2[2];
  std::complex<double> *fourier_funcs_3[3];
  std::complex<double> *fourier_1,*fourier_2,*fourier_3;

  fourier_1 = new std::complex<double> [geom.N_interp*geom.fourier_size_sym];
  fourier_2 = new std::complex<double> [geom.N_interp*geom.fourier_size_sym];
  fourier_3 = new std::complex<double> [geom.N_interp*geom.fourier_size_sym];

  /////////dpol_sq
  for( int iii = 0 ; iii < geom.N_interp * geom.N_theta ; iii++ ){temp_1[iii] = input_funcs[0][iii] / (eq.jacob[iii] * eq.g_pp[iii]) ;}
  calc_fourier_sym(temp_1,fourier_1,geom.N_interp,geom.N_theta);
  for( int iii = 0 ; iii < geom.N_interp * geom.N_theta ; iii++ ){temp_1[iii] = input_funcs[1][iii] * eq.f_psi[iii] * eq.f_psi[iii] / (eq.jacob[iii] * eq.g_pp[iii]) ;}
  calc_fourier_sym(temp_1,fourier_2,geom.N_interp,geom.N_theta);
  calc_fourier_sym(input_funcs[2],fourier_3,geom.N_interp,geom.N_theta);
  fourier_funcs_3[0] = fourier_1; fourier_funcs_3[1] = fourier_2; fourier_funcs_3[2] = fourier_3;
  for(int iii=0;iii<geom.N_interp;iii++){ result[iii] = dpol_sq_xi(fourier_funcs_3,main_mod,sec_mod,geom.N_interp,geom.N_theta,iii) ;}

  delete [] fourier_3;
  
  ////////dpol_dphi
  for( int iii = 0 ; iii < geom.N_interp * geom.N_theta ; iii++ ){ temp_1[iii] = input_funcs[0][iii] * eq.f_psi[iii] / (eq.jacob[iii] * eq.g_pp[iii]) ;}
  calc_fourier_sym(temp_1,fourier_1,geom.N_interp,geom.N_theta);
  for( int iii = 0 ; iii < geom.N_interp * geom.N_theta ; iii++ ){ temp_1[iii] = input_funcs[1][iii] * input_funcs[2][iii] * eq.g_phph[iii] ;}
  calc_fourier_sym(temp_1,fourier_2,geom.N_interp,geom.N_theta);
  fourier_funcs_2[0] = fourier_1; fourier_funcs_2[1] = fourier_2;
  for(int iii=0;iii<geom.N_interp;iii++){ result[iii] -= geom.tor_mod * imag_unit * dpol_xi(fourier_funcs_2,fg_sym,main_mod,sec_mod,geom.N_interp,geom.N_theta,iii); }
    
  ////////dphi_dpol
  for( int iii = 0 ; iii < geom.N_interp * geom.N_theta ; iii++ ){ temp_1[iii] = input_funcs[0][iii] * input_funcs[1][iii] * eq.g_phph[iii] / (eq.jacob[iii] * eq.g_pp[iii]) ;}
  calc_fourier_sym(temp_1,fourier_1,geom.N_interp,geom.N_theta);
  for( int iii = 0 ; iii < geom.N_interp * geom.N_theta ; iii++ ){ temp_1[iii] = eq.f_psi[iii] * input_funcs[2][iii] ;}
  calc_fourier_sym(temp_1,fourier_2,geom.N_interp,geom.N_theta);
  fourier_funcs_2[0] = fourier_1; fourier_funcs_2[1] = fourier_2;
  for(int iii=0;iii<geom.N_interp;iii++){ result[iii] -= geom.tor_mod * imag_unit * dpol_xi(fourier_funcs_2,fg_sym,main_mod,sec_mod,geom.N_interp,geom.N_theta,iii); }
  
  delete [] fourier_2;
  
  ////////dphi_sq
  for( int iii = 0 ; iii < geom.N_interp * geom.N_theta ; iii++ ){temp_1[iii] = input_funcs[0][iii] * input_funcs[1][iii] * input_funcs[2][iii] * eq.g_phph[iii] * eq.g_phph[iii] ;}
  calc_fourier_sym(temp_1,fourier_1,geom.N_interp,geom.N_theta);
  if(m_diff>=0){ for(int iii=0;iii<geom.N_interp;iii++){ result[iii] -= geom.tor_mod * geom.tor_mod * fourier_1[iii*geom.fourier_size_sym+m_diff]; } }
  else{ for(int iii=0;iii<geom.N_interp;iii++){ result[iii] -= geom.tor_mod * geom.tor_mod * std::conj(fourier_1[iii*geom.fourier_size_sym-m_diff]); } }
  
  delete[] temp_1;
  delete[] fourier_1;
}


void dwedge_dpol_xi(double *input_funcs[],std::complex<double> result[],geom_shape geom,equil_fields eq,int main_mod,int sec_mod)
{
  clean_grid(result,geom.N_psi);
  
  //Calculates 

  int fourier_size_sym=geom.fourier_size_sym;
  int fourier_size_full=geom.fourier_size_full;
  bool fg_sym[2]={true,true};
  int m_diff=main_mod-sec_mod;

  double *temp_1 = new double[geom.N_interp*geom.N_theta];

  std::complex<double> *fourier_funcs_2[2];
  std::complex<double> *fourier_funcs_3[3];
  std::complex<double> *fourier_1,*fourier_2,*fourier_3;

  fourier_1 = new std::complex<double> [geom.N_interp*geom.fourier_size_sym];
  fourier_2 = new std::complex<double> [geom.N_interp*geom.fourier_size_sym];
  fourier_3 = new std::complex<double> [geom.N_interp*geom.fourier_size_sym];

  /////////dpol_dphi
  for( int iii = 0 ; iii < geom.N_interp * geom.N_theta ; iii++ ){ temp_1[iii] = input_funcs[0][iii] * input_funcs[1][iii] * eq.g_phph[iii] ;}
  calc_fourier_sym(temp_1,fourier_1,geom.N_interp,geom.N_theta);
  calc_fourier_sym(input_funcs[2],fourier_2,geom.N_interp,geom.N_theta);
  fourier_funcs_2[0] = fourier_1; fourier_funcs_2[1] = fourier_2;
  for(int iii=0;iii<geom.N_interp;iii++){ result[iii] = geom.tor_mod * imag_unit * dpol_xi(fourier_funcs_2,fg_sym,main_mod,sec_mod,geom.N_interp,geom.N_theta,iii) ;}

  /////////dpol_sq
  for( int iii = 0 ; iii < geom.N_interp * geom.N_theta ; iii++ ){temp_1[iii] = input_funcs[0][iii] / (eq.jacob[iii] * eq.g_pp[iii]) ;}
  calc_fourier_sym(temp_1,fourier_1,geom.N_interp,geom.N_theta);
  for( int iii = 0 ; iii < geom.N_interp * geom.N_theta ; iii++ ){temp_1[iii] = input_funcs[1][iii] * eq.f_psi[iii] ;}
  calc_fourier_sym(temp_1,fourier_2,geom.N_interp,geom.N_theta);
  calc_fourier_sym(input_funcs[2],fourier_3,geom.N_interp,geom.N_theta);
  fourier_funcs_3[0] = fourier_1; fourier_funcs_3[1] = fourier_2; fourier_funcs_3[2] = fourier_3;
  for(int iii=0;iii<geom.N_interp;iii++){ result[iii] -= dpol_sq_xi(fourier_funcs_3,main_mod,sec_mod,geom.N_interp,geom.N_theta,iii) ;}


  delete [] temp_1;
  delete [] fourier_1;
  delete [] fourier_2;
  delete [] fourier_3;
}

void dpol_dwedge_xi(double *input_funcs[],std::complex<double> result[],geom_shape geom,equil_fields eq,int main_mod,int sec_mod)
{
  clean_grid(result,geom.N_psi);
  
  //Calculates 

  int fourier_size_sym=geom.fourier_size_sym;
  int fourier_size_full=geom.fourier_size_full;
  bool fg_sym[2]={true,true};
  int m_diff=main_mod-sec_mod;

  double *temp_1 = new double[geom.N_interp*geom.N_theta];

  std::complex<double> *fourier_funcs_2[2];
  std::complex<double> *fourier_funcs_3[3];
  std::complex<double> *fourier_1,*fourier_2,*fourier_3;

  fourier_1 = new std::complex<double> [geom.N_interp*geom.fourier_size_sym];
  fourier_2 = new std::complex<double> [geom.N_interp*geom.fourier_size_sym];
  fourier_3 = new std::complex<double> [geom.N_interp*geom.fourier_size_sym];

  /////////dpol_dphi
  calc_fourier_sym(input_funcs[0],fourier_1,geom.N_interp,geom.N_theta);
  for( int iii = 0 ; iii < geom.N_interp * geom.N_theta ; iii++ ){ temp_1[iii] = input_funcs[1][iii] * input_funcs[2][iii] * eq.g_phph[iii] ;}
  calc_fourier_sym(temp_1,fourier_2,geom.N_interp,geom.N_theta);
  fourier_funcs_2[0] = fourier_1; fourier_funcs_2[1] = fourier_2;
  for(int iii=0;iii<geom.N_interp;iii++){ result[iii] = geom.tor_mod * imag_unit * dpol_xi(fourier_funcs_2,fg_sym,main_mod,sec_mod,geom.N_interp,geom.N_theta,iii) ;}

  /////////dpol_sq
  calc_fourier_sym(input_funcs[0],fourier_1,geom.N_interp,geom.N_theta);
  for( int iii = 0 ; iii < geom.N_interp * geom.N_theta ; iii++ ){temp_1[iii] = input_funcs[1][iii] * eq.f_psi[iii] / (eq.jacob[iii] * eq.g_pp[iii]) ;}
  calc_fourier_sym(temp_1,fourier_2,geom.N_interp,geom.N_theta);
  calc_fourier_sym(input_funcs[2],fourier_3,geom.N_interp,geom.N_theta);
  fourier_funcs_3[0] = fourier_1; fourier_funcs_3[1] = fourier_2; fourier_funcs_3[2] = fourier_3;
  for(int iii=0;iii<geom.N_interp;iii++){ result[iii] -= dpol_sq_xi(fourier_funcs_3,main_mod,sec_mod,geom.N_interp,geom.N_theta,iii) ;}


  delete [] temp_1;
  delete [] fourier_1;
  delete [] fourier_2;
  delete [] fourier_3;
}


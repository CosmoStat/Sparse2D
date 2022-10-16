/******************************************************************************
**                   Copyright (C) 2014 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Fred Ngole
**
**    Date:  24/09/2014
**    
**    File:  sr_util.cc
**
*******************************************************************************
**
**    DESCRIPTION   
**    ----------- 
**    Tools used in sprite.cc and useful beyond the scope of super-resolution             
**
**
******************************************************************************/
#include "sr_util.h"

int compute_centroid(Ifloat *img,double sig,double *cent,int niter)
{
  int Nx = (*img).nx();
  int Ny = (*img).ny();
  cent[0]=0;
  cent[1]=0;

  // Loop to compute the centroid 
  int i,j,k;
  Ifloat w(Nx,Ny);
  Ifloat wimg(Nx,Ny);
  for (i=0;i<Nx;i++) 
    for (j=0;j<Ny;j++) w(i,j)=1;

  for (k=0; k < niter; k++)
    {
      double tot=0; 
      for (i=0;i<Nx;i++) 
	for (j=0;j<Ny;j++) 
	  {
	    if (k>0)
	      {
		w(i,j)=exp(-(pow(i+1-cent[0],2)+pow(j+1-cent[1],2))/(2*pow(sig,2)));
	      }
	    wimg(i,j)=(*img)(i,j)*w(i,j);
	   
	    tot+=wimg(i,j);
	  }
      cent[0]=0;
      cent[1]=0;
      for (i=0;i<Nx;i++)
	{
	  for (j=0;j<Ny;j++)
	    {
	      cent[0]+=wimg(i,j)*(i+1);
	      cent[1]+=wimg(i,j)*(j+1);
	    }
	} 
      cent[0]/=tot;
      cent[1]/=tot;
      
    }
  
  // Memory free
  w.free();
  wimg.free();
  if (isnan(cent[0]) || isinf(cent[0]) || isnan(cent[1]) || isinf(cent[1]))
    {
      return 0;
    }
  else return 1; 
}

double gauss2d(double amp,double x,double y,double theta,double sigx,double sigy)
{
  double X = x*cos(theta) + y*sin(theta);
  double Y = -x*sin(theta) + y*cos(theta);
  double output = amp*exp(-(pow(X/sigx,2)+pow(Y/sigy,2))/2);
  return output;
}

void gauss2darray(double amp,double theta,double sigx,double sigy, double xcen,double ycen,Ifloat &im_out)
{
  int Nx = im_out.nx(), Ny = im_out.ny();
  int i,j;
  for (i=0;i<Nx;i++)
    for (j=0;j<Ny;j++)
      im_out(i,j)=gauss2d(amp,i-xcen,j-ycen,theta,sigx,sigy);
}

void dist_map(Ifloat &map,double xcen,double ycen)
{
  int Nx = map.nx(), Ny = map.ny();
  int i,j;
  for (i=0;i<Nx;i++)
    for (j=0;j<Ny;j++)
      map(i,j)=pow(i-xcen,2)+pow(j-ycen,2);
}

double gauss2d_der1sig(Ifloat res,Ifloat gauss_approx,Ifloat distances_map,double sig) // 1st derivative of a quadratic data attachment in a 2d gaussian fitting with respect to sigma^2
{
  int Nx = res.nx(), Ny = res.ny();
  double gamma = pow(sig,2);
  Ifloat temp(Nx,Ny);
  hadamard(res,gauss_approx,temp);
  double der = dot(temp,distances_map);
  der = -der/(2*pow(gamma,2));
  temp.free();
  return der;
}

double gauss2d_der1sigx(Ifloat res, Ifloat gauss_approx, double sigx, double xcen) // 1st derivative of a quadratic data attachment in a 2d gaussian fitting with respect to sigmax^2
{
  int Nx = res.nx(), Ny = res.ny();
  double gamma = pow(sigx,2);
  double der=0;
  int i,j;
  for (i=0;i<Nx;i++)
    for (j=0;j<Ny;j++)
      der = der - pow(i-xcen,2)*gauss_approx(i,j)*res(i,j);
  der = der/(2*pow(gamma,2));
  return der;
}

double gauss2d_der1sigy(Ifloat res, Ifloat gauss_approx, double sigy, double ycen) // 1st derivative of a quadratic data attachment in a 2d gaussian fitting with respect to sigmay^2
{
  int Nx = res.nx(), Ny = res.ny();
  double gamma = pow(sigy,2);
  double der=0;
  int i,j;
  for (i=0;i<Nx;i++)
    for (j=0;j<Ny;j++)
      der = der - pow(j-ycen,2)*gauss_approx(i,j)*res(i,j);
  der = der/(2*pow(gamma,2));
  return der;
}

double gauss2d_der2sig(Ifloat res,Ifloat gauss_approx,Ifloat distances_map,double sig) // 2nd derivative of a quadratic data attachment in a 2d gaussian fitting with respect to sigma^2
{
  int Nx = res.nx(), Ny = res.ny();
  double gamma = pow(sig,2);
  int i,j;
  Ifloat temp(Nx,Ny);
  
  for (i=0;i<Nx;i++)
    for (j=0;j<Ny;j++)
      temp(i,j) = res(i,j)-distances_map(i,j)*(res(i,j)-gauss_approx(i,j))/(4*gamma);
  hadamard(temp,gauss_approx,temp);
  double der = dot(temp,distances_map);
  der = der/pow(gamma,3);
  return der;
}

double gauss2d_der2sigx(Ifloat res, Ifloat gauss_approx, double sigx, double xcen) // 2nd derivative of a quadratic data attachment in a 2d gaussian fitting with respect to sigmax^2
{
  int Nx = res.nx(), Ny = res.ny();
  double gamma = pow(sigx,2);
  double der=0;
  int i,j;
  for (i=0;i<Nx;i++)
    for (j=0;j<Ny;j++)
      der = der + pow(i-xcen,2)*gauss_approx(i,j)*(res(i,j)-((pow(i-xcen,2)/(4*gamma))*(res(i,j)-gauss_approx(i,j))))/pow(gamma,3);
  
  
  return der;
}

double gauss2d_der2sigy(Ifloat res, Ifloat gauss_approx, double sigy, double ycen) // 2nd derivative of a quadratic data attachment in a 2d gaussian fitting with respect to sigmay^2
{
  int Nx = res.nx(), Ny = res.ny();
  double gamma = pow(sigy,2);
  double der=0;
  int i,j;
  for (i=0;i<Nx;i++)
    for (j=0;j<Ny;j++)
      der = der + pow(j-ycen,2)*gauss_approx(i,j)*(res(i,j)-((pow(j-ycen,2)/(4*gamma))*(res(i,j)-gauss_approx(i,j))))/pow(gamma,3);
  return der;
}

double gauss2d_der1xcen(Ifloat res, Ifloat gauss_approx, double sig, double xcen)
{
  int Nx = res.nx(), Ny = res.ny();
  double der=0;
  int i,j;
  for (i=0;i<Nx;i++)
    for (j=0;j<Ny;j++)
      der = der-(i-xcen)*gauss_approx(i,j)*res(i,j);
  der = der/pow(sig,2);
  return der;
}

double gauss2d_der1ycen(Ifloat res, Ifloat gauss_approx, double sig, double ycen)
{
  int Nx = res.nx(), Ny = res.ny();
  double der=0;
  int i,j;
  for (i=0;i<Nx;i++)
    for (j=0;j<Ny;j++)
      der = der-(j-ycen)*gauss_approx(i,j)*res(i,j);
  der = der/pow(sig,2);
  return der;
}

double gauss2d_der2xcen(Ifloat res, Ifloat gauss_approx, double sig,double xcen)
{
  int Nx = res.nx(), Ny = res.ny();
  double der=0;
  int i,j;
  
  for (i=0;i<Nx;i++)
    for (j=0;j<Ny;j++)
      der = der+(gauss_approx(i,j)/pow(sig,2))*((1-pow((i-xcen)/sig,2))*res(i,j) + gauss_approx(i,j)*pow((i-xcen)/sig,2));
  return der;
}

double gauss2d_der2ycen(Ifloat res, Ifloat gauss_approx, double sig,double ycen)
{
  int Nx = res.nx(), Ny = res.ny();
  double der=0;
  int i,j;  
  for (i=0;i<Nx;i++)
    for (j=0;j<Ny;j++)
      der = der+(gauss_approx(i,j)/pow(sig,2))*((1-pow((j-ycen)/sig,2))*res(i,j) + gauss_approx(i,j)*pow((j-ycen)/sig,2));
  return der;
}

double gauss2d_der2xcenycen(Ifloat res, Ifloat gauss_approx, double sig, double xcen, double ycen)
{
  int Nx = res.nx(), Ny = res.ny();
  double der=0;
  int i,j; 
  for (i=0;i<Nx;i++)
    for (j=0;j<Ny;j++)
      der = der + (i-xcen)*(j-ycen)*(gauss_approx(i,j)-res(i,j))*gauss_approx(i,j);
  der = der/pow(sig,4);
  return der;
}


void gauss2D_fitting_first_guess(Ifloat im,double*sigx,double*sigy,double*xcen,double*ycen,double*amp,int smooth_ker)
{
  int Nx = im.nx(), Ny = im.ny();
  Ifloat kern(smooth_ker,smooth_ker);
  int i,j;
  for (i=0;i<smooth_ker;i++)
    for (j=0;j<smooth_ker;j++)
      {
	kern(i,j) = 1.0/pow((double) smooth_ker,(double) 2);	
      }
  Ifloat im_smooth(Nx,Ny);
  psf_convol(im, kern, im_smooth,True,False);
  int ind;
  *amp = im_smooth.max(ind);
  *xcen = floor(ind/Nx);
  *ycen = ind%Nx;
  double del = im(*xcen,*ycen)/exp(1.0);              //1/e value
  i=0;
  int i0 = *ycen;
  while (((i0+i) < Ny) && ((i0-i+1) > 0) && (abs(im(*xcen,i0+i)) > abs(del)) && (abs(im(*xcen,i0-i)) > abs(del))) 
    i+=1;
  *sigx = i;
  i=0;
   i0 = *xcen;
  while (((i0+i) < Nx) && ((i0-i+1) > 0) && (abs(im(i0+i,*ycen)) > abs(del)) && (abs(im(i0-i,*ycen)) > abs(del))) 
    i+=1;
  *sigy = i;
  kern.free();
  im_smooth.free();
}

void gauss2d_fit(Ifloat im,double *xcen,double *ycen,double *amp,double sig0,double *sig,int nb_iter_alter,int nb_subiter,double *mse,int max_test) // Isotropic 2D gaussian fitting grdient descent/Newton method 
{
  int Nx = im.nx(), Ny = im.ny();
  int k,l,i,j;
  //*amp = amp0;
  *sig = sig0;
  Ifloat res(Nx,Ny);
  Ifloat gauss_approx(Nx,Ny);
  Ifloat distance_map(Nx,Ny);
  dist_map(distance_map,*xcen,*ycen);
  for (l=0;l<nb_iter_alter;l++)
    {
      for (k=0;k<nb_subiter;k++)
	{
      	  gauss2darray(*amp,0,*sig,*sig,*xcen,*ycen,gauss_approx);
	  for (i=0;i<Nx;i++)
	    for (j=0;j<Ny;j++)
	      res(i,j) = im(i,j) - gauss_approx(i,j);
	  /*if (mse!=NULL)
	  // mse[k*nb_sub_iter+l] = res.energy();*/
	  double der1 = gauss2d_der1sig(res,gauss_approx,distance_map,*sig);
	  if (der1==0)
	    {
	      cout << "Warning: 1st derivative = 0 in 2D gaussian fitting; exit the current loop..."<< endl;
	      break;
	    }
	  double der2 =  abs(gauss2d_der2sig(res,gauss_approx,distance_map,*sig));
	  double temp_sig =  *sig;
	  int count=0;
	  double mu;
	  if (abs(der2)>0)
	    mu = 1.0/abs(der2);
	  else
	    {
	      while ((der2==0) && (count < max_test))
		{
		  cout << "Warning: 2nd derivative = 0 in 2D gaussian fitting"<< endl;
		  count++;
		  temp_sig = sqrt(pow(temp_sig,2)+pow(*sig,2)/10);
		  der2 = abs(gauss2d_der2sig(res,gauss_approx,distance_map,temp_sig));
		}
	      if (der2==0)
		{
		  cout << "Warning: 2nd derivative = 0 in 2D gaussian fitting; exit the current loop..."<< endl;
		  break;
		}
	    }
	  double temp = pow(*sig,2) - mu*der1;
	  count = 0;
	  while ((temp<=0) && (count <max_test))
	    {
	      cout << "Warning: sigma^2 <= 0 in 2D gaussian fitting; reducing the step..."<< endl;
	      count++;
	      mu=mu/2;
	      temp = pow(*sig,2) - mu*der1;
	    }
	  if (temp>0)
	    *sig = sqrt(temp);
	}
      gauss2darray(1,0,*sig,*sig,*xcen,*ycen,gauss_approx);
      *amp = dot(gauss_approx,im)/gauss_approx.energy();
      for (k=0;k<nb_subiter;k++)
	{
	  gauss2darray(*amp,0,*sig,*sig,*xcen,*ycen,gauss_approx);
	  for (i=0;i<Nx;i++)
	    for (j=0;j<Ny;j++)
	      res(i,j) = im(i,j) - gauss_approx(i,j);
	  double gradx = gauss2d_der1xcen(res,gauss_approx,*sig,*xcen);
	  double hessx = gauss2d_der2xcen(res,gauss_approx,*sig,*xcen);
	  double tol=1.0;
	  double mu = MIN(1.0/abs(hessx),tol/abs(gradx));
	  *xcen = *xcen -mu*gradx;
	}
      for (k=0;k<nb_subiter;k++)
	{
	  gauss2darray(*amp,0,*sig,*sig,*xcen,*ycen,gauss_approx);
	  for (i=0;i<Nx;i++)
	    for (j=0;j<Ny;j++)
	      res(i,j) = im(i,j) - gauss_approx(i,j);
	  double grady = gauss2d_der1ycen(res,gauss_approx,*sig,*ycen);
	  double hessy = gauss2d_der2ycen(res,gauss_approx,*sig,*ycen);
	  double tol=1.0;
	  double mu = MIN(1.0/abs(hessy),tol/abs(grady));
	  *ycen = *ycen -mu*grady;
	}
    } 
  gauss_approx.free();
  res.free();
  distance_map.free();
}
 
void gauss2d_fit_2(Ifloat im,double *xcen,double *ycen,double *amp,double *sigx,double *sigy,int nb_iter_alter,int nb_subiter,double *mse,int max_test) // 2D gaussian fitting grdient descent/Newton method 
// *sigx and *sigy have to be initialized with reasonable first guess by the calling routine
{
  int Nx = im.nx(), Ny = im.ny();
  int k,l,i,j;
  Ifloat res(Nx,Ny);
  Ifloat gauss_approx(Nx,Ny);
  Ifloat distance_map(Nx,Ny);
  dist_map(distance_map,*xcen,*ycen);
  for (l=0;l<nb_iter_alter;l++)
    {
      for (k=0;k<nb_subiter;k++)
	{
      	  gauss2darray(*amp,0,*sigx,*sigy,*xcen,*ycen,gauss_approx);
	  for (i=0;i<Nx;i++)
	    for (j=0;j<Ny;j++)
	      res(i,j) = im(i,j) - gauss_approx(i,j);
	  /*if (mse!=NULL)
	  // mse[k*nb_sub_iter+l] = res.energy();*/
	    double der1 = gauss2d_der1sigx(res,gauss_approx,*sigx,*xcen);
	  
	  if (der1==0)
	    {
	      cout << "Warning: 1st derivative = 0 in 2D gaussian fitting; exit the current loop..."<< endl;
	      break;
	    }
	  double der2 =  abs(gauss2d_der2sigx(res,gauss_approx,*sigx,*xcen));
	  double temp_sig =  *sigx;
	  int count=0;
	  double mu;
	  if (abs(der2)>0)
	    mu = 1.0/abs(der2);
	  else
	    {
	      while ((der2==0) && (count < max_test))
		{
		  cout << "Warning: 2nd derivative = 0 in 2D gaussian fitting"<< endl;
		  count++;
		  temp_sig = sqrt(pow(temp_sig,2)+pow(*sigx,2)/10);
		  der2 = abs(gauss2d_der2sigx(res,gauss_approx,temp_sig,*xcen));
		}
	      if (der2==0)
		{
		  cout << "Warning: 2nd derivative = 0 in 2D gaussian fitting; exit the current loop..."<< endl;
		  break;
		}
	    }
	  double temp = pow(*sigx,2) - mu*der1;
	  count = 0;
	  while ((temp<=0) && (count <max_test))
	    {
	      cout << "Warning: sigmax^2 <= 0 in 2D gaussian fitting; reducing the step..."<< endl;
	      count++;
	      mu=mu/2;
	      temp = pow(*sigx,2) - mu*der1;
	    }
	  if (temp>0)
	    *sigx = sqrt(temp);
	}
      for (k=0;k<nb_subiter;k++)
	{
      	  gauss2darray(*amp,0,*sigx,*sigy,*xcen,*ycen,gauss_approx);
	  for (i=0;i<Nx;i++)
	    for (j=0;j<Ny;j++)
	      res(i,j) = im(i,j) - gauss_approx(i,j);
	  /*if (mse!=NULL)
	  // mse[k*nb_sub_iter+l] = res.energy();*/
	    double der1 = gauss2d_der1sigy(res,gauss_approx,*sigy,*ycen);
	  
	  if (der1==0)
	    {
	      cout << "Warning: 1st derivative = 0 in 2D gaussian fitting; exit the current loop..."<< endl;
	      break;
	    }
	  double der2 =  abs(gauss2d_der2sigy(res,gauss_approx,*sigy,*ycen));
	  double temp_sig =  *sigy;
	  int count=0;
	  double mu;
	  if (der2>0)
	    mu = 1.0/der2;
	  else
	    {
	      while ((der2==0) && (count < max_test))
		{
		  cout << "Warning: 2nd derivative = 0 in 2D gaussian fitting"<< endl;
		  count++;
		  temp_sig = sqrt(pow(temp_sig,2)+pow(*sigy,2)/10);
		  der2 = abs(gauss2d_der2sigy(res,gauss_approx,temp_sig,*ycen));
		}
	      if (der2==0)
		{
		  cout << "Warning: 2nd derivative = 0 in 2D gaussian fitting; exit the current loop..."<< endl;
		  break;
		}
	    }
	  double temp = pow(*sigy,2) - mu*der1;
	  count = 0;
	  while ((temp<=0) && (count <max_test))
	    {
	      cout << "Warning: sigmay^2 <= 0 in 2D gaussian fitting; reducing the step..."<< endl;
	      count++;
	      mu=mu/2;
	      temp = pow(*sigy,2) - mu*der1;
	    }
	  if (temp>0)
	    *sigy = sqrt(temp);
	}
      gauss2darray(1,0,*sigx,*sigy,*xcen,*ycen,gauss_approx);
      *amp = dot(gauss_approx,im)/gauss_approx.energy();
      for (k=0;k<nb_subiter;k++)
	{
	  gauss2darray(*amp,0,*sigx,*sigy,*xcen,*ycen,gauss_approx);
	  for (i=0;i<Nx;i++)
	    for (j=0;j<Ny;j++)
	      res(i,j) = im(i,j) - gauss_approx(i,j);
	  double gradx = gauss2d_der1xcen(res,gauss_approx,*sigx,*xcen);
	  double hessx = gauss2d_der2xcen(res,gauss_approx,*sigx,*xcen);
	  double tol=1.0;
	  double mu = MIN(1.0/abs(hessx),tol/abs(gradx));
	  *xcen = *xcen -mu*gradx;
	}
      for (k=0;k<nb_subiter;k++)
	{
	  gauss2darray(*amp,0,*sigx,*sigy,*xcen,*ycen,gauss_approx);
	  for (i=0;i<Nx;i++)
	    for (j=0;j<Ny;j++)
	      res(i,j) = im(i,j) - gauss_approx(i,j);
	  double grady = gauss2d_der1ycen(res,gauss_approx,*sigy,*ycen);
	  double hessy = gauss2d_der2ycen(res,gauss_approx,*sigy,*ycen);
	  double tol=1.0;
	  double mu = MIN(1.0/abs(hessy),tol/abs(grady));
	  *ycen = *ycen -mu*grady;
	}
	} 
  gauss_approx.free();
  res.free();
  distance_map.free();
}

int compute_centroid_arr(fltarray *data,double sig,double **cent,int niter)
{
  int Nz = (*data).nz(), Nx = (*data).nx(), Ny = (*data).ny();
  int k,i,j;
  Ifloat temp_img(Nx,Ny);
  double *cent_temp = (double*)malloc(2*sizeof(double));
  int flag = 1;
  for (k=0;k<Nz;k++)
    {
      for (i=0;i<Nx;i++)
	for(j=0;j<Ny;j++) temp_img(i,j)=(*data)(i,j,k);
      
      int flag_temp = compute_centroid(&temp_img,sig,cent_temp,niter);
      if (flag_temp==0) flag=0;
      cent[0][k]=cent_temp[0];
      cent[1][k]=cent_temp[1];
    }
  // Memory free 
  temp_img.free();
  free(cent_temp);
  return flag;
}

void thresholding(Ifloat *data_in,Ifloat *data_out,Ifloat *thresh,bool thresh_type)
{
  // If thresh_type = 0 => hard thresholding ; otherwise => soft thresholding ; default is hard
  int Nx = (*data_in).nx(),Ny = (*data_in).ny();
  int i,j;
  for (i=0;i<Nx;i++)
    for(j=0;j<Ny;j++)
      {
	if (fabs((*data_in)(i,j)) < (*thresh)(i,j)) (*data_out)(i,j) = 0;
	else (*data_out)(i,j) = (*data_in)(i,j) - (*thresh)(i,j)*thresh_type*sign_num((*data_in)(i,j));
      }
}

void thresholding(Ifloat *data_in,Ifloat *data_out,double thresh,bool thresh_type)
{
  // If thresh_type = 0 => hard thresholding ; otherwise => soft thresholding ; default is hard
  int Nx = (*data_in).nx(),Ny = (*data_in).ny();
  int i,j;
  for (i=0;i<Nx;i++)
    for(j=0;j<Ny;j++)
      {
	if (fabs((*data_in)(i,j)) < thresh) (*data_out)(i,j) = 0;
	else (*data_out)(i,j) = (*data_in)(i,j) - thresh*thresh_type*sign_num((*data_in)(i,j));
      }
}

void thresholding(fltarray *data_in,fltarray *data_out,fltarray* thresh,bool thresh_type)
{
  // If thresh_type = 0 => hard thresholding ; otherwise => soft thresholding ; default is hard
  int Nx = (*data_in).nx(),Ny = (*data_in).ny(),Nz=(*data_in).nz();
  int i,j,k;
  for (i=0;i<Nx;i++)
    for(j=0;j<Ny;j++)
      for(k=0;k<Nz;k++)
	{
	 
	  if (fabs((*data_in)(i,j,k)) < (*thresh)(i,j,k)) (*data_out)(i,j,k) = 0;
	  else (*data_out)(i,j,k) = (*data_in)(i,j,k) - (*thresh)(i,j,k)*thresh_type*sign_num((*data_in)(i,j,k));
	}
}

void circ_thresh(Ifloat *data_in,Ifloat *data_out,double r,double *cent)
{
  int Nx=data_in->nx(),Ny=data_in->ny(),i,j;
  for (i=0;i<Nx;i++)
    for(j=0;j<Ny;j++)
      {
	(*data_out)(i,j)=(*data_in)(i,j);
	if ((pow(i-cent[0],2)+pow(j-cent[1],2)) > pow(r,2)) (*data_out)(i,j)=0;
      }
}

void circ_thresh(fltarray *data_in,fltarray *data_out,double r,double **cent)
{
  int Nx=(*data_in).nx(),Ny=(*data_in).ny(),Nz=(*data_in).nz(),i,j,k;
  for (k=0;k<Nz;k++)
    for (i=0;i<Nx;i++)
      for(j=0;j<Ny;j++)
	{
	  (*data_out)(i,j,k)=(*data_in)(i,j,k);
	  if ((pow(i-cent[0][k],2)+pow(j-cent[1][k],2)) > pow(r,2)) (*data_out)(i,j,k)=0;
	}
}

void flux_est(Ifloat *data, double r,double *cent,double*flux)
{
  int Nx=data->nx(),Ny=data->ny();
  Ifloat data_temp(Nx,Ny);
  circ_thresh(data,&data_temp,r,cent);
  *flux = data_temp.total();
  data_temp.free();
}

void flux_est(fltarray *data, double r,double **cent,double*flux,int ref_im_ind)
{
  int Nx=data->nx(),Ny=data->ny(),Nz=data->nz(),i,j,k;
  fltarray data_temp(Nx,Ny,Nz);
  Ifloat im_temp(Nx,Ny);
  circ_thresh(data,&data_temp,r,cent);
  for (i=0;i<Nx;i++)
    for(j=0;j<Ny;j++) im_temp(i,j)=data_temp(i,j,ref_im_ind);
  flux[ref_im_ind] = im_temp.total();
  for (k=0;k<Nz;k++) 
    {
      for (i=0;i<Nx;i++)
	for(j=0;j<Ny;j++) im_temp(i,j)=data_temp(i,j,k);
      if (k!=ref_im_ind)
	{
	  flux[k] = im_temp.total();
	  flux[k] = flux[k]/flux[ref_im_ind];
	}  
    }
  flux[ref_im_ind] = 1;
  im_temp.free();
  data_temp.free();
}

int wl_trans(Ifloat *Dat,mr_opt opt,MultiResol* MR_Data)
{
  type_transform transf = opt.transf;
  int nb_sc = opt.nb_sc;
  int nb_usc = opt.nb_usc;
  sb_type_norm Norm = opt.Norm;
  type_undec_filter U_Filter = opt.U_Filter;
  type_sb_filter SB_Filter = opt.SB_Filter;
  type_border Bord = opt.Bord;
  Bool Verbose = False;
  FilterAnaSynt *FAS = opt.FAS;
  MR_Data->alloc ((*Dat).nl(), (*Dat).nc(), nb_sc, transf, FAS, Norm, nb_usc,U_Filter);
  MR_Data->Border = Bord;
  MR_Data->Verbose = Verbose;
  MR_Data->transform ((*Dat));
  
  if (MR_Data != NULL) return 1;
  else return 0;
}
int sign_num (double a)
{
  if (a==0) return 0;
  else return a/fabs(a);
}

mr_opt mr_opt_init(int Nx,type_transform transf,int nb_sc,int nb_usc,sb_type_norm Norm,type_undec_filter U_Filter,type_sb_filter SB_Filter,type_border Bord,Bool Verbose )
{
  FilterAnaSynt *FAS ;
  mr_opt opt ;
  opt.transf=transf ;
  if (nb_sc<0)  
    opt.nb_sc = floor(log(Nx)/log(2))-1 ;
  else 
    opt.nb_sc = nb_sc;
  if (nb_usc<0 || nb_usc<opt.nb_sc)  
    opt.nb_usc = opt.nb_sc ;
  else 
    opt.nb_usc = nb_usc;
  
  opt.Norm = Norm ;
  opt.U_Filter = U_Filter ;
  opt.SB_Filter = SB_Filter ; 
  opt.Bord = Bord ;
  opt.FAS = new FilterAnaSynt();
  if ((transf == TO_MALLAT) || (transf == TO_UNDECIMATED_MALLAT))
    {
      (opt.FAS)->Verbose = Verbose;
      (opt.FAS)->alloc(SB_Filter);
    }
  return opt;
}
void wl_thresholding(MultiResol*wav_coeff,fltarray *thresh,bool thresh_type)
{
  int Nx = (*wav_coeff).size_band_nc(0),Ny=(*wav_coeff).size_band_nl(0),Nz=(*wav_coeff).nbr_band()-1,i,j,k; // The coarse scale is not thresholded
  fltarray temp_coeff(Nx,Ny,Nz), coeff_thresh(Nx,Ny,Nz);
  for (i=0;i<Nx;i++)
    for (j=0;j<Ny;j++)
      for(k=0;k<Nz;k++) temp_coeff(i,j,k) = (*wav_coeff)(k,j,i);
  thresholding(&temp_coeff,&coeff_thresh,thresh,thresh_type);
  for (i=0;i<Nx;i++)
    for (j=0;j<Ny;j++)
      for(k=0;k<Nz;k++)  (*wav_coeff)(k,j,i)=coeff_thresh(i,j,k);
  temp_coeff.free();
  coeff_thresh.free();
}

double noise_est(Ifloat *data)
{
  int Nx = data->nx(),Ny=data->ny(),i,j; 
  Ifloat temp(Nx,Ny);
  Ifloat wl_temp(Nx,Ny);
  mr_opt den_opt = mr_opt_init(Nx,TO_UNDECIMATED_MALLAT,2,2,NORM_L2);
  MultiResol MR_Data;
  int flag = wl_trans(data,den_opt,&MR_Data);
  for (i=0;i<Nx;i++)
    for (j=0;j<Ny;j++) wl_temp(i,j) = MR_Data(0,j,i);
  
  double sig = get_noise(wl_temp,METHOD_MAD);
  MR_Data.free();
  
  temp.free();
  wl_temp.free();
  return sig;
}

int noise_est(fltarray* data,fltarray* noise_arr)
{
  int Nx = data->nx(),Ny=data->ny(),Nz=data->nz(),i,j,k; 
  Ifloat temp(Nx,Ny);
  Ifloat wl_temp(Nx,Ny);
  mr_opt den_opt = mr_opt_init(Nx,TO_UNDECIMATED_MALLAT,2,2,NORM_L2);
  for (k=0;k<Nz;k++)
    {
      for (i=0;i<Nx;i++)
	for (j=0;j<Ny;j++) temp(i,j) = (*data)(i,j,k);
      MultiResol MR_Data;
      int flag = wl_trans(&temp,den_opt,&MR_Data);
      for (i=0;i<Nx;i++)
	for (j=0;j<Ny;j++) wl_temp(i,j) = MR_Data(0,j,i);
      double sig = get_noise(wl_temp,METHOD_MAD);
      for (i=0;i<Nx;i++)
	for (j=0;j<Ny;j++) (*noise_arr)(i,j,k)=sig;
      
      MR_Data.free();
    }
  
  temp.free();
  wl_temp.free();
  if (noise_arr != NULL) return 1;
  else return 0;
}



int wl_gnoise_est(MultiResol* wav_coeff,fltarray* noise_arr) 
// Wavelet coefficient of an array dominated with noise 
{
  int Nx = (*wav_coeff).size_band_nc(0),Ny=(*wav_coeff).size_band_nl(0),Nz=(*wav_coeff).nbr_band()-1,i,j,k; 
  Ifloat temp(Nx,Ny);
  for (k=0;k<Nz;k++)
    {
      for (i=0;i<Nx;i++)
	for (j=0;j<Ny;j++) temp(i,j) = (*wav_coeff)(k,j,i);
      double sig = get_noise (temp,METHOD_MAD);
      
      if(sig < 0)
	cout << "Noise estimation failed"<< endl;
      for (i=0;i<Nx;i++)
	for (j=0;j<Ny;j++) (*noise_arr)(i,j,k)=sig;
    }
  temp.free();
  return 1;
}

void wl_gnoise_est_dat(Ifloat *data,mr_opt opt,fltarray* noise_arr) // This function assumes uncorrelated nearly gaussian noise
{
  MultiResol wav_coeff;
  int Nx = data->nx(),Ny = data->ny(),k,i,j;
  int nb_band = wav_coeff.computeNbBand(Nx,Ny,opt.nb_sc,opt.transf,opt.nb_usc);
  double* wl_nrm = new double[nb_band];
  mr_norm(opt,wl_nrm,nb_band,Nx,Ny);
  double sig = noise_est(data);
  
  for (k=0;k<nb_band-1;k++)
    {
      for (i=0;i<Nx;i++)
	for (j=0;j<Ny;j++) (*noise_arr)(i,j,k)=wl_nrm[k]*sig;
    }
  delete [] wl_nrm;
  
}
void wl_filter(Ifloat *img,Ifloat*img_filt,mr_opt opt,double nsig,bool thresh_type,fltarray* noise_map,fltarray* coeff_thresh,Ifloat *coarse_scale)
{
  MultiResol wav_coeff;

  int flag_trans = wl_trans(img,opt,&wav_coeff);
  int Nx = wav_coeff.size_band_nc(0),Ny=wav_coeff.size_band_nl(0),Nz=wav_coeff.nbr_band()-1,i,j,k; // The coarse scale is not thresholded
  fltarray noise_arr(Nx,Ny,Nz);
  if (noise_map!=NULL)
    {
      for(k=0;k<Nz;k++)
	{ 
	  for (i=0;i<Nx;i++)
	    for (j=0;j<Ny;j++)
	      {
		noise_arr(i,j,k)=(*noise_map)(i,j,k)*nsig;
	      }
	  
	}
    }
  else
    {
      wl_gnoise_est_dat(img,opt,&noise_arr);
      for (i=0;i<Nx;i++)
	for (j=0;j<Ny;j++)
	  for(k=0;k<Nz;k++)  noise_arr(i,j,k)=noise_arr(i,j,k)*nsig;
    }

  wl_thresholding(&wav_coeff,&noise_arr,thresh_type);
   
  if (coarse_scale!=NULL)
    for (i=0;i<Nx;i++)
      for (j=0;j<Ny;j++)
	  (*coarse_scale)(i,j) = wav_coeff(Nz,j,i);
  wav_coeff.recons(*img_filt,opt.Bord);
  
  if (coeff_thresh!=NULL)
     for (i=0;i<Nx;i++)
       for (j=0;j<Ny;j++) 
	 for (k=0;k<Nz;k++) (*coeff_thresh)(i,j,k) = wav_coeff(k,j,i);
  
  wav_coeff.free();
  noise_arr.free();
  
}

void wl_filter_analysis(Ifloat img,Ifloat&img_filt,mr_opt opt,double nsig,bool thresh_type,fltarray std_map,MultiResol &coeff_init,int nb_iter,double mu,double *mse)
{
  int Nx = coeff_init.size_band_nc(0),Ny=coeff_init.size_band_nl(0),Nz=coeff_init.nbr_band(),i,j,k,l;
  MultiResol xold;
  MultiResol x;
  MultiResol z;
  MultiResol xtemp;
  MultiResol grad_step;
  fltarray thresholds(Nx,Ny,Nz-1);
    
  for (k=0;k<Nz-1;k++)
    for (j=0;j<Ny;j++)
      for (i=0;i<Nx;i++)
	thresholds(i,j,k)=nsig*std_map(i,j,k);
  z.alloc(Ny,Nx,opt.nb_sc,opt.transf);
  for (k=0;k<Nz;k++)
    for (j=0;j<Ny;j++)
      for (i=0;i<Nx;i++)
	z(k,j,i) = coeff_init(k,j,i);
  double t=1;
    
  xold.alloc(Ny,Nx,opt.nb_sc,opt.transf);
  //grad_step.alloc(Ny,Nx,opt.nb_sc,opt.transf);
  x.alloc(Ny,Nx,opt.nb_sc,opt.transf);
  xtemp.alloc(Ny,Nx,opt.nb_sc,opt.transf);
 
  Ifloat temp_recons(Nx,Ny);
    for (l=0;l<nb_iter;l++)
    {
        z.recons(temp_recons,opt.Bord);
        for (j=0;j<Ny;j++)
            for (i=0;i<Nx;i++)
                img_filt(i,j) = img(i,j)-temp_recons(i,j);
         if (mse!=NULL)
             mse[l]=img_filt.energy();
        
      
        int flag_trans = wl_trans(&img_filt,opt,&grad_step);
        
      // ----------- Gradient Step ------------- //
      for (k=0;k<Nz;k++)
          for (j=0;j<Ny;j++)
              for (i=0;i<Nx;i++)
              {
                  x(k,j,i) = z(k,j,i)+mu*grad_step(k,j,i);
                  xtemp(k,j,i)=x(k,j,i);
              }
        grad_step.free(); // Allocated in wl_trans
      // ----------- Proximal step ------------ //
      wl_thresholding(&xtemp,&thresholds,thresh_type);
     
      for (k=0;k<Nz;k++)
	for (j=0;j<Ny;j++)
	  for (i=0;i<Nx;i++)
	      x(k,j,i) = x(k,j,i)-xtemp(k,j,i);
      double t_temp = (1+sqrt(4*pow(t,2)+1))/2;
      double lambda = 1 + (t-1)/t_temp;
      t = t_temp;
     for (k=0;k<Nz;k++)
	for (j=0;j<Ny;j++)
	  for (i=0;i<Nx;i++)
	    {
	      z(k,j,i) = xold(k,j,i)+lambda*(x(k,j,i)-xold(k,j,i));
	      xold(k,j,i) = x(k,j,i);
	    }
    }
    
  z.recons(temp_recons,opt.Bord);
  for (j=0;j<Ny;j++)
    for (i=0;i<Nx;i++)
      img_filt(i,j) = img(i,j)-temp_recons(i,j);
  for (k=0;k<Nz;k++)
    for (j=0;j<Ny;j++)
      for (i=0;i<Nx;i++)
	coeff_init(k,j,i) = z(k,j,i);
  temp_recons.free();
  x.free();
  xtemp.free();
  xold.free();
  z.free();
  thresholds.free();
}


void wl_arr_filter(fltarray *data,fltarray*data_filt,mr_opt opt,double nsig,bool thresh_type)
{
  int Nz = data->nz(), Nx = data->nx(), Ny = data->ny(),i,j,k;
  for (k=0;k<Nz;k++)
    {
      Ifloat temp_im(Nx,Ny);
      Ifloat temp_im_filt(Nx,Ny);
      for (i=0;i<Nx;i++)
	for(j=0;j<Ny;j++) temp_im(i,j)=(*data)(i,j,k);
     
      wl_filter(&temp_im,&temp_im_filt,opt,nsig,thresh_type);
       
      for (i=0;i<Nx;i++)
	for(j=0;j<Ny;j++) (*data_filt)(i,j,k)=temp_im_filt(i,j);
      temp_im.free();
      temp_im_filt.free();	  
    }  	 
}


void coadd_frames(fltarray Data, Ifloat &Ima,Bool GetMedian)
{
  int Nx, Ny, Nz;
  Nx = Data.nx();
  Ny = Data.ny();
  Nz = Data.nz();
  int x,y,z;
  double Val;
  if (GetMedian == False)
    {
    
      for (x=0; x < Nx; x++)
	for (y=0; y < Ny; y++)
	  {	
	    Val = 0.;
	    for (z=0; z < Nz; z++) Val += Data(x,y,z);
	    if (Nz > 0) Ima(x,y) = (float) (Val / (double) Nz);
	    else Ima(x,y) = 0;
	  }
    }
  else
    {
      
      float *V = new float [Nz];
      for (x=0; x < Nx; x++)
	for (y=0; y < Ny; y++)
	  {
	    if (Nz > 0) 
	      {
		for (z=0; z < Nz; z++)  V[z] = Data(x,y,z);
	   
		Ima(x,y) = get_median(V, Nz);
	      }
	    else Ima(x,y) = 0.;
	 
	  }

      delete [] V;

    }
   
   
}

/*********************************************************************/

void coadd_frames(fltarray Dat, Ifloat &Ima, float Zoom,
		      Bool GetMedian,int MaxDist, 
		      int Surface, type_interp_corr TypeInterp,
		      Bool OptNoOffset,Bool MeanSub,Bool Verbose)
{

  int PosX = -1;
  int PosY = -1;
  int ValRet, i, nz = Dat.nz();
  int Nx1 = (int)(Zoom * Dat.nx());
  int Ny1 = (int)(Zoom * Dat.ny());
  int k,l;
  fltarray DatOut(Nx1, Ny1, Dat.nz());   
  xcorr_def Xcor;
   Xcor.xc_np = nz;
   Xcor.xc_x = new double [nz];
   Xcor.xc_y = new double [nz];
   Xcor.xc_level = new double [nz];
   for (i=0; i < nz; i++) Xcor.xc_x[i] = Xcor.xc_y[i] = Xcor.xc_level[i] = 0.;
   Xcor.xc_dx = MaxDist;
   Xcor.xc_dy = MaxDist;
   if ((Surface < 1) || (Surface > Dat.nx()/2))
   {
      Xcor.xc_hx = MAX(0, Dat.nx()/2 -1 - MaxDist);
      Xcor.xc_hy = MAX(0, Dat.ny()/2 -1 - MaxDist);
   }
   else
   {
      Xcor.xc_hx = Surface;
      Xcor.xc_hy = Surface;
   }
   if (PosX < 0) PosX =  Dat.nx() / 2;
   if (PosY < 0) PosY =  Dat.ny() / 2;
   int D1 = MIN(PosX-MaxDist,PosY-MaxDist);
   int D2 = MIN(Dat.nx()-PosX-1-MaxDist,Dat.ny()-PosY-1-MaxDist);
   int DM = MIN(D1,D2);
   if (DM < 2)
   {
      cout << "Error: search position cannot be at the border of the image ... " << endl;
      exit(-1);
   } 
   if (Surface >= DM) 
   {
      Surface = DM-1;
      Xcor.xc_hx = Surface;
      Xcor.xc_hy = Surface;
      cout << "Warning: surface must be decreased: new value = " << Surface << endl;
   } 
       
   //Xcor.xc_hx = 10;
   //Xcor.xc_hy = 10;
   Xcor.xc_method = 0;
   char *interp_method= (char*)StringCorrInterp(TypeInterp);
   
   float **TabFrame = new float * [nz];
   float * TabMean = new float [nz];
   for (i = 0; i < nz; i++) TabFrame[i] = Dat.buffer() + i*Dat.nx()*Dat.ny();
   
   
   if (OptNoOffset == False) // otherwise we calculte the offset frame by correlation
   {  
      Ifloat ImaFrame;
      Xcor.xc_init = 0;
      
      // Test position and surface
      
      for (i = 0; i < nz; i++) 
      {
          Xcor.xc_x[i] = Xcor.xc_y[i] = 0.;
	  if (MeanSub == True)
	  {
	     ImaFrame.alloc(TabFrame[i], Dat.ny(), Dat.nx());
	     TabMean[i] = (float) average(ImaFrame);
 	     for (k=0; k < ImaFrame.nl(); k++)
	     for (l=0; l < ImaFrame.nc(); l++)  ImaFrame(k,l) -= TabMean[i];
	  }
      }
      float *Pattern = TabFrame[0];
      
      ValRet = cube_get_offset(TabFrame, Pattern, 
                               Dat.nx(), Dat.ny(), Dat.nz(), 
			       PosX, PosY, &Xcor);
      for (i = 0; i < nz; i++) 
      {
 	  ImaFrame.alloc(TabFrame[i], Dat.ny(), Dat.nx());
	  if (MeanSub == True)
	  {
 	     for (k=0; k < ImaFrame.nl(); k++)
	     for (l=0; l < ImaFrame.nc(); l++)  ImaFrame(k,l) += TabMean[i];
	  }
      }
      // cout <<"ValRet = " << ValRet << endl; 
      if (Verbose == True)
        for (i = 0; i < nz; i++)
           cout << "Offset Frame: " << i+1 << ", dx = " << Xcor.xc_x[i] << " dy = " << Xcor.xc_y[i] << endl;     

      
   }
   
  
   
   if ((OptNoOffset == False) || (Zoom != 1))
   {
       if (Verbose == True) cout << "Resample cube ... " << endl;
       
      
       int ValRet, i;
    
    
    float **TabFrameOut;
    
    TabFrameOut = new float * [nz];   
    for (i = 0; i < nz; i++) 
        TabFrameOut[i] = DatOut.buffer() + i*DatOut.nx()*DatOut.ny();
    
    ValRet = cube_resample(TabFrame, TabFrameOut, &Xcor, Dat.nx(), Dat.ny(), Dat.nz(),
                           DatOut.nx(), DatOut.ny(), interp_method); 
    
    coadd_frames(DatOut, Ima, GetMedian);
    free(TabFrameOut);
       
   }
else 
    {
      coadd_frames(Dat, Ima);
    }
  DatOut.free();
  free(Xcor.xc_x);
  free(Xcor.xc_y);
  free(Xcor.xc_level);
  free(TabFrame);
  free(TabMean);
   
}

void shift_est(fltarray *data, double sig_gfit, double *sig_vec,double **cent,double **shift,int ref_im_ind,double nsig)
{
  int Nx  = data->nx();
  int Ny  = data->ny();
  int Nz  = data->nz();
  fltarray thresh(Nx,Ny,Nz);
  fltarray data_filt(Nx,Ny,Nz);
  int i,j,k;
  for(k=0;k<Nz;k++)
    {
      Ifloat temp(Nx,Ny);
      for (i=0;i<Nx;i++)
	for (j=0;j<Ny;j++) temp(i,j)=(*data)(i,j,k);
      int psnr=floor(temp.max()/sig_vec[k]);
      if (psnr<nsig)
	{
	  for (i=0;i<Nx;i++)
	    for (j=0;j<Ny;j++) thresh(i,j,k)=(psnr-1)*sig_vec[k]; 
	}
      else 
	{
	  for (i=0;i<Nx;i++)
	    for (j=0;j<Ny;j++) thresh(i,j,k)=nsig*sig_vec[k]; 
	}
      temp.free();
    }
  thresholding(data,&data_filt,&thresh);
  int flag_cen = compute_centroid_arr(&data_filt,sig_gfit,cent);
  thresh.free();    
  for(k=0;k<Nz;k++)
    {
      shift[0][k] = cent[0][k]-cent[0][ref_im_ind];
      shift[1][k] = cent[1][k]-cent[1][ref_im_ind];
    }
  data_filt.free();
}

double sinc(double x)
{
  double out;
  if (x==0) out = 1;
  else out = sin(x)/x;
  return out;
}

void lanczos(double u,double*mask,int mask_rad)
{
  int i;
  
  for (i=0;i<2*mask_rad+1;i++)  mask[i] = sinc(Pi*(u-(i-mask_rad)))*sinc(Pi*(u-(i-mask_rad))/mask_rad);
	      
}

void lanczos(double u[],double**mask,int mask_rad)
{
  int i,j;
  for (i=0;i<2*mask_rad+1;i++) 
    for (j=0;j<2*mask_rad+1;j++) mask[i][j] = sinc(Pi*(u[0]-(i-mask_rad)))*sinc(Pi*(u[0]-(i-mask_rad))/mask_rad)* sinc(Pi*(u[1]-(j-mask_rad)))*sinc(Pi*(u[1]-(j-mask_rad))/mask_rad);
}

void lanczos(double u[],Ifloat &mask,int mask_rad)
{
  int i,j;
  for (i=0;i<2*mask_rad+1;i++) 
    for (j=0;j<2*mask_rad+1;j++) mask(i,j) = sinc(Pi*(u[0]-(i-mask_rad)))*sinc(Pi*(u[0]-(i-mask_rad))/mask_rad)* sinc(Pi*(u[1]-(j-mask_rad)))*sinc(Pi*(u[1]-(j-mask_rad))/mask_rad);
}


void lanczos_upsamp(Ifloat &im,int upfact,Ifloat &im_out, double *u,int kern_rad)
{
  int i,j,k;
    Bool flag=False;
    if (u == NULL)
    {
      flag = True;
      double * u = new double[2];
      u[0]=0;
      u[1]=0;
    }
    int Nx  = im.nx();
    int Ny  = im.ny();
    Ifloat im_temp(Nx,Ny*upfact);
    double** hlin = new double*[upfact];
    double** hcol = new double*[upfact];
    double* buff_in_lines = new double[Ny]; // For lines upsampling
    double* buff_out_lines = new double[Ny];
    double* buff_in_cols = new double[Nx]; // For columns upsampling
    double* buff_out_cols = new double[Nx];
    for(i = 0; i < upfact; i++)
    {
        hlin[i] = new double[2*kern_rad+1];
        hcol[i] = new double[2*kern_rad+1];
    }
   
    for (i=0;i<upfact;i++)
    {
        lanczos(u[0]+(double)(i)/upfact,hlin[i],kern_rad);
        lanczos(u[1]+(double)(i)/upfact,hcol[i],kern_rad);
    }
    // Lines upasmpling
    for (i=0;i<Nx;i++)
      {
	for (j=0;j<Ny;j++)
	  buff_in_lines[j] = im(i,j);
	for (k=0;k<upfact;k++)
	  {
	    convol1D_cent(buff_in_lines,buff_out_lines,hlin[k],Ny,2*kern_rad+1);
	    for (j=0;j<Ny;j++)
	      im_temp(i,j*upfact+k) = buff_out_lines[j];
	  }	
      }
    
    // Columns upasmpling
    for (j=0;j<Ny*upfact;j++)
      {
	for (i=0;i<Nx;i++)
	  buff_in_cols[i] = im_temp(i,j);
	for (k=0;k<upfact;k++)
	  {
	    convol1D_cent(buff_in_cols,buff_out_cols,hcol[k],Nx,2*kern_rad+1);
	    for (i=0;i<Nx;i++)
	      im_out(i*upfact+k,j) = buff_out_cols[i];
	  }	
	  }
    
    for(i = 0; i < upfact-1; i++)
    {
        delete [] hlin[i];
        delete [] hcol[i];
    }
    delete [] hlin;
    delete [] hcol;
    delete [] buff_in_lines;
    delete [] buff_out_lines;
    delete [] buff_in_cols;
    delete [] buff_out_cols;
    im_temp.free();
    if (flag==True)
      delete [] u;
      
}

void lanczos_stacking(fltarray &data,Ifloat &im_upsamp,int upfact,double **u,int kern_rad,Bool mean_en)
{
  int Nx  = data.nx();
  int Ny  = data.ny();
  int Nz  = data.nz();
  Ifloat im_lr_temp(Nx,Ny);
  Ifloat im_hr_temp(Nx*upfact,Ny*upfact);
  double * uk = new double[2];
  int i,j,k;
  fltarray data_upsamp(Nx*upfact,Ny*upfact,Nz);
 
  for (k=0;k<Nz;k++)
    {
      for (i=0;i<Nx;i++)
	for (j=0;j<Ny;j++)
	  im_lr_temp(i,j) = data(i,j,k);
       
      uk[0] = -u[0][k];
      uk[1] = -u[1][k];
      
      lanczos_upsamp(im_lr_temp,upfact,im_hr_temp,uk,kern_rad);
      for (i=0;i<Nx*upfact;i++)
	for (j=0;j<Ny*upfact;j++)
	  data_upsamp(i,j,k)=im_hr_temp(i,j);
    }
  
  im_stacking(data_upsamp,im_upsamp,mean_en);
  im_lr_temp.free();
  im_hr_temp.free();
  data_upsamp.free();
  delete [] uk;
}

void im_stacking(fltarray &data,Ifloat &im_stack,Bool mean_en)
{
  int Nx  = data.nx();
  int Ny  = data.ny();
  int Nz  = data.nz();
  int i,j,k;
  float * stack = new float[Nz];
  for (i=0;i<Nx;i++)
    for (j=0;j<Ny;j++)
      {
	 im_stack(i,j)=0;
	if (mean_en==True)
	  {
	    for (k=0;k<Nz;k++)
	      im_stack(i,j) =  im_stack(i,j)+data(i,j,k);
	    im_stack(i,j) =  im_stack(i,j)/Nz;
	  }
	else
	  {
	    for (k=0;k<Nz;k++)
	      stack[k]=data(i,j,k);
	    im_stack(i,j) = get_median(stack,Nz);
	    
	  }
      }
  delete [] stack;
}

void convol_decim_mat(Ifloat mask,int D,Ifloat mat,int Nx, int Ny) // Nx Ny Size of input image
{
    int i,j,i_max = mask.nx(),j_mask=mask.ny();
    
}

void decim(Ifloat *img,Ifloat *img_filt,Ifloat *img_dec,int D,Ifloat* mask,Bool filt_en,Bool fft_en)
{
  int Nx = img->nx(),Ny = img->ny(),i,j;
  
  if ((filt_en==True))
    if (fft_en==True) psf_convol((*img), (*mask), (*img_filt),True,False);
  if (D==1) (*img_dec)=(*img);
  else
    {
      if (filt_en==True)
	{
	  for(i=0;i<floor(Nx/D);i++)
	    for(j=0;j<floor(Ny/D);j++) (*img_dec)(i,j)=(*img_filt)(D*i,D*j);
	}
      else
	{
	  for(i=0;i<floor(Nx/D);i++)
	    for(j=0;j<floor(Ny/D);j++) (*img_dec)(i,j)=(*img)(D*i,D*j);	  
	}
    }	
}

void rotate(Ifloat *input,Ifloat *output,int dir) // Same as in IDL
{
  dir = dir%8;
  int Nx = input->nx(), Ny = input->ny(),x,y,x1,y1;
  switch(dir){
  case 0 : 
    for (x=0;x<Nx;x++)
      for (y=0;y<Ny;y++)
	{
	  x1=x;
	  y1=y;
	  (*output)(x1,y1)=(*input)(x,y);
	}
    break;
  case 1 : 
    for (x=0;x<Nx;x++)
      for (y=0;y<Ny;y++)
	{
	  x1=Ny-y-1;
	  y1=x;
	  (*output)(x1,y1)=(*input)(x,y);
	}
    break;
  case 2 : 
    for (x=0;x<Nx;x++)
      for (y=0;y<Ny;y++)
	{
	  x1=Nx-x-1;
	  y1=Ny-y-1;
	  (*output)(x1,y1)=(*input)(x,y);
	}
    break;
  case 3 : 
    for (x=0;x<Nx;x++)
      for (y=0;y<Ny;y++)
	{
	  x1=y;
	  y1=Nx-x-1;
	  (*output)(x1,y1)=(*input)(x,y);
	}
    break;
  case 4 : 
    for (x=0;x<Nx;x++)
      for (y=0;y<Ny;y++)
	{
	  x1=y;
	  y1=x;
	  (*output)(x1,y1)=(*input)(x,y);
	}
    break;
  case 5 : 
    for (x=0;x<Nx;x++)
      for (y=0;y<Ny;y++)
	{
	  x1=Nx-x-1;
	  y1=y;
	  (*output)(x1,y1)=(*input)(x,y);
	}
    break;
  case 6 : 
    for (x=0;x<Nx;x++)
      for (y=0;y<Ny;y++)
	{
	  x1=Ny-y-1;
	  y1=Nx-x-1;
	  (*output)(x1,y1)=(*input)(x,y);
	}
    break;
  case 7 : 
    for (x=0;x<Nx;x++)
      for (y=0;y<Ny;y++)
	{
	  x1=x;
	  y1=Ny-y-1;
	  (*output)(x1,y1)=(*input)(x,y);
	}
    break;
  }
    
}

void transpose_decim(Ifloat*im_in,Ifloat*im_out,int D)
{
  int Nxlr = im_in->nx(),Nylr = im_in->ny(),Nxhr = im_out->nx(),Nyhr = im_out->ny(),i,j;
  im_in->reform (Nxlr*Nylr);
  im_out->reform (Nxhr*Nyhr);
  for (i=0;i<Nxlr*Nylr;i++)
    {
      int Li = i/Nxlr;
      int Ci = i - Li*Nxlr;
      (*im_out)(Li*Nxlr*pow((double) D,(double)2)+Ci*D)=(*im_in)(i);
    }
  im_in->reform (Nxlr,Nylr);
  im_out->reform(Nxhr,Nyhr);
}



void power_meth(void(*opname)(double*,double*),int size_input,int *nb_iter,double *spec_rad,int nb_iter_max,double tol)
{
  double L_old=0;
  double L=1;
  *nb_iter=0;
  double * input = (double*)malloc(size_input*sizeof(double));
  //double * output = (double*)malloc(size_output*sizeof(double));
  int i;
  
  for (i=0;i<size_input;i++)
    {
      input[i]=rand();
    }
  double l0 = norm2(input,size_input);
  scale (input,size_input,1.0/l0);
    while ((fabs(L_old-L)/L <tol) & (*nb_iter <nb_iter_max))
      {
	scale(input,size_input,1.0/L);
	(*opname)(input,input);
	L_old = L;
	L = norm2(input,size_input);
	(*nb_iter)++;
      }
  *spec_rad = L*(1+tol); // Better have a larger value
  if (*nb_iter==nb_iter_max) cout << "Warning: max number of iteration reached in spectral radius estimation"<< endl;
  free(input);
}

double norm2 (double * x,int siz)
{
  int i;
  double out=0;
  for (i=0;i<siz;i++) out+=pow(x[i],2);
  out = sqrt(out);
  return out;
			    
}

double dot(Ifloat Im1,Ifloat Im2)
{
  int i,j,Nx=Im1.nx(),Ny=Im1.ny();
  double dot_prod=0;
  for (i=0;i<Nx;i++)
    for (j=0;j<Ny;j++) 
      dot_prod = dot_prod+Im1(i,j)*Im2(i,j);
  return dot_prod;
}

void hadamard(Ifloat Im1,Ifloat Im2,Ifloat &Im_out)
{
  int i,j,Nx=Im1.nx(),Ny=Im1.ny();
  for (i=0;i<Nx;i++)
    for (j=0;j<Ny;j++) 
      Im_out(i,j) = Im1(i,j)*Im2(i,j);
}

void scale (double * x,double siz,double a)
{
  int i;
  for (i=0;i<siz;i++) x[i]*=a;
}

void ineq_cons(Ifloat &A, Ifloat B,double b) // projects A onto  E = {X/X(i,j)>=b*B(i,j)}
{
  int i,j,Nx=A.nx(),Ny=B.ny();
  for (i=0;i<Nx;i++)
    for (j=0;j<Ny;j++) 
      if (A(i,j)<b*B(i,j)) A(i,j) = b*B(i,j);
}


void reverse(double*u,int length,double*ur)
{
  int i;
  for(i=0;i<length;i++)
    ur[i]=u[length-1-i];
}

void holes(double*h,int sc_nb,int length,double*hout)
{
  int i;
  for (i=0;i<length;i++)
    hout[i*(int)pow((double) 2,(double) sc_nb)]=h[i];
}

void convol1D(double*in,double*out,double*h,int length,int filt_length)
// Edges zero convention
{
  int i,k;
  for (i=0;i<length;i++)
    {
      out[i]=0;
      for (k=0;k<MIN(i,filt_length);k++)
	out[i] = out[i]+in[i-k]*h[filt_length-1-k];
    }
}

void correl2D(Ifloat im1,Ifloat im2, Ifloat &im_correl)
// Edge zero convention
{
    
    int i,j,k,l,Nx = im1.nx(),Ny=im1.ny(),icen=(Nx-1)/2,jcen=(Ny-1)/2;
    for (i=0;i<Nx; i++) {
        for (j=0; j<Ny; j++) {
            im_correl(i,j)=0;
            for (k=MAX(i-Nx,icen-Nx)+1; k<MIN(i,icen)+1; k++) {
                for (l=MAX(j-Ny,jcen-Ny)+1; l<MIN(j,jcen)+1; l++) {
                    im_correl(i,j) = im_correl(i,j)+im1(i-k,j-l)*im2(icen-k,jcen-l);
                }
            }
        }
    }
    
}

void convol1D_cent(double*in,double*out,double*h,int length,int filt_length)
// Edges zero convention
{
  int i,k;
  int rad = (filt_length-1)/2; // The filter length is supposed to be odd
  for (i=0;i<length;i++)
    {
      out[i]=0;
      for (k=MAX(-rad,i-(length-1));k<MIN(i,rad);k++)
	out[i] = out[i]+in[i-k]*h[rad-k];
    }
}


void sep_filt2d(Ifloat* im_in,Ifloat*im_out,double *h,double *g,int lh,int lg)
{
  int i,j,k,Nx,Ny;
  Nx = im_in->nx();
  Ny = im_in->ny();
  Ifloat im_temp(Nx,Ny);
  
  // Columns filtering
  for(i=0;i<Nx;i++)      
      for(j=0;j<Ny;j++)
	for (k=0;k<MIN(j,lg);k++)
	  im_temp(i,j) = im_temp(i,j) + (*im_in)(i,k)*g[lg-k-1];
   
  // Lines filtering
  for(j=0;j<Ny;j++)      
      for(i=0;i<Nx;i++)
	for (k=0;k<MIN(i,lh);k++)
	  (*im_out)(i,j) = (*im_out)(i,j) + im_temp(k,j)*h[lh-k-1];
   im_temp.free();
}

void grad(Ifloat im,Ifloat &im_gradx,Ifloat &im_grady,double* hx,double* hy,int lhx,int lhy)
{
  double *filt_id = new double[3];
  filt_id[0]=0;
  filt_id[1]=1;
  filt_id[2]=0;
  sep_filt2d(&im,&im_gradx,hx,filt_id,lhx,3);
  sep_filt2d(&im,&im_gradx,filt_id,hy,3,lhy);
  delete [] filt_id; 
}

void transp_sep_filt2d(Ifloat* im_in,Ifloat*im_out,double *h,double *g,int lh,int lg)
{
  double *hr = new double[lh];
  double *gr = new double[lg];
  reverse(h,lh,hr);
  reverse(g,lg,gr);
  sep_filt2d(im_in,im_out,hr,gr,lh,lg);
  delete [] hr;
  delete [] gr;
}

void grad_cons_noise_est(Ifloat im_in,Ifloat& im_out,Ifloat &res,double * sig,double gamma,int nb_iter,double * mse_iter,double mu, double *hx,double *hy,int lhx,int lhy)
{
  int i,j,k,Nx,Ny;
  Nx = im_in.nx();
  Ny = im_in.ny();
  Ifloat im_gradx(Nx,Ny);
  Ifloat im_grady(Nx,Ny);
  Ifloat im_tgradx(Nx,Ny);
  Ifloat im_tgrady(Nx,Ny);
  if (hx==NULL)
    {
      lhx = 3;
      hx = new double[lhx];
      hx[0]=-1/sqrt(2);
      hx[1]=0;
      hx[2]=1/sqrt(2);
    }
  if (hy==NULL)
    {
      lhy = 3;
      hy = new double[lhy];
      hy[0]=-1/sqrt(2);
      hy[1]=0;
      hy[2]=1/sqrt(2);
    }
  double *filt_id = new double[3];
  filt_id[0]=0;
  filt_id[1]=1;
  filt_id[2]=0;
  for (k=0;k<nb_iter;k++)
   {
     
     for (i=0;i<Nx;i++)
       for (j=0;j<Ny;j++)
	 res(i,j) = im_in(i,j)-im_out(i,j);
     if(mse_iter!=NULL)
       mse_iter[k] = res.energy()/(Nx*Ny);
     
     grad(im_out,im_gradx,im_grady,hx,hy,lhx,lhy);
     
     transp_sep_filt2d(&im_gradx,&im_tgradx,hx,filt_id,lhx,lhx);
     transp_sep_filt2d(&im_grady,&im_tgrady,filt_id,hy,lhy,lhy);    
     for (i=0;i<Nx;i++)
       for (j=0;j<Ny;j++)
	 im_out(i,j) = im_out(i,j) + mu*(res(i,j)-gamma*(im_tgradx(i,j)+im_tgrady(i,j)));
   }
  *sig = get_noise(res,METHOD_MAD);
  im_gradx.free();
  im_grady.free();
  im_tgradx.free();
  im_tgrady.free();
  delete [] hx;
  delete [] filt_id; 
}


void mad_res_noise_est(Ifloat im_in,Ifloat& im_out,Ifloat &res,double * sig,int nb_iter,double * mse_iter,double *sig_iter,double mu, int nsig,Bool sig_up,int thresh_type)
{
  int i,j,k,Nx,Ny;
  Nx = im_in.nx();
  Ny = im_in.ny();
  *sig = get_noise(im_in,METHOD_MAD);
  
  for (k=0;k<nb_iter;k++)
   {
     
     for (i=0;i<Nx;i++)
       for (j=0;j<Ny;j++)
	 res(i,j) = im_in(i,j)-im_out(i,j);
     if(mse_iter!=NULL)
       mse_iter[k] = res.energy()/(Nx*Ny); 
     for (i=0;i<Nx;i++)
       for (j=0;j<Ny;j++)
	 im_out(i,j) = im_out(i,j) + mu*res(i,j);
     if (sig_up==True)
       *sig = get_noise(res,METHOD_MAD);
     if((sig_iter!=NULL)&&(sig_up==True))
       sig_iter[k] = *sig; 
     thresholding(&im_out,&im_out,(nsig)*(*sig),thresh_type);
   }
  *sig = get_noise(res,METHOD_MAD);
 
}

/*void randomn(double *x,double sig,int length,double mean)

{
  int i;  
  random_device rd;
  default_random_engine generator(rd());
  normal_distribution<double> distribution(mean,sig);
  for(i=0;i<length;i++)
    {  
      x[i] = distribution(generator);     
    }
    }*/

void randomngsl(double *x,double sig,int length,double mean)
{
  int i;  
  gsl_rng *rng;
  const gsl_rng_type * T;
  T = gsl_rng_default;
  rng = gsl_rng_alloc (T);
  int s = renewSeed();
  gsl_rng_set(rng, s);
  
  for(i=0;i<length;i++)
    {  
      x[i] =gsl_ran_gaussian(rng,sig)+mean;     
    }
}

/*void randomn(Ifloat *mat, double sig,double mean)
{
  random_device rd;
  default_random_engine generator(rd());
  normal_distribution<double> distribution(mean,sig);
  int i,j,Nx,Ny;
  Nx = mat->nx();
  Ny = mat->ny();
  for (i=0;i<Nx;i++)
    for (j=0;j<Ny;j++)
      {
      	double x;
	(*mat)(i,j) = distribution(generator);
      }
      }*/

void randomngsl(Ifloat *mat, double sig,double mean)
{
   gsl_rng *rng;
  const gsl_rng_type * T;
  T = gsl_rng_default;
  rng = gsl_rng_alloc (T);
  int s = renewSeed();
  gsl_rng_set(rng, s);
  int i,j,Nx,Ny;
  Nx = mat->nx();
  Ny = mat->ny();
  for (i=0;i<Nx;i++)
    for (j=0;j<Ny;j++)
      {
      	double x;
	(*mat)(i,j) = gsl_ran_gaussian(rng,sig)+mean;
      }
}

/*void randomn(fltarray *mat,double *sig,double mean)
{
  random_device rd;
  default_random_engine generator(rd());
  int i,j,k,Nx,Ny,Nz;
  Nx = mat->nx();
  Ny = mat->ny();
  Nz = mat->nz();
  for (k=0;k<Nz;k++)
    {
    normal_distribution<double> distribution(mean,sig[k]);
    for (i=0;i<Nx;i++)
      for (j=0;j<Ny;j++)
	{
	  
	  (*mat)(i,j,k) = distribution(generator);
	  
	}
    }
    }*/

void randomngsl(fltarray *mat,double *sig,double mean)
{
   gsl_rng *rng;
  const gsl_rng_type * T;
  T = gsl_rng_default;
  rng = gsl_rng_alloc (T);
  int s = renewSeed();
  gsl_rng_set(rng, s);
  int i,j,k,Nx,Ny,Nz;
  Nx = mat->nx();
  Ny = mat->ny();
  Nz = mat->nz();
  for (k=0;k<Nz;k++)
    {
      for (i=0;i<Nx;i++)
	for (j=0;j<Ny;j++)
	  {
	    
	    (*mat)(i,j,k) = gsl_ran_gaussian(rng,sig[k])+mean;
	    
	    
	  }
    }
}

/*void randomn(fltarray *mat,double *sig,double*mean)
{
  random_device rd;
  default_random_engine generator(rd());
  int i,j,k,Nx,Ny,Nz;
  Nx = mat->nx();
  Ny = mat->ny();
  Nz = mat->nz();
  for (k=0;k<Nz;k++)
    {
    normal_distribution<double> distribution(mean[k],sig[k]);
    for (i=0;i<Nx;i++)
      for (j=0;j<Ny;j++)
	{
	  
	  (*mat)(i,j,k) = distribution(generator);
	}
    }
    }*/

void randomngsl(fltarray *mat,double *sig,double*mean)
{
   gsl_rng *rng;
  const gsl_rng_type * T;
  T = gsl_rng_default;
  rng = gsl_rng_alloc (T);
  int s = renewSeed();
  gsl_rng_set(rng, s);
  int i,j,k,Nx,Ny,Nz;
  Nx = mat->nx();
  Ny = mat->ny();
  Nz = mat->nz();
  for (k=0;k<Nz;k++)
    {
    for (i=0;i<Nx;i++)
      for (j=0;j<Ny;j++)
	{
	  
	  (*mat)(i,j,k) = gsl_ran_gaussian(rng,sig[k])+mean[k];
	}
    }
}

double var_iter(double x,int n,double *mean,double*M)
{
  double delta = x - *mean;
  *mean = *mean+delta/(n+1);
  *M = *M+delta*(x-*mean);
  double var = 0;
  if (n>0)
    var = *M/n;
  return var;
}

void mr_support_filt(Ifloat img,Ifloat &img_filt, mr_opt opt,double nsig,int nb_iter,double*mse,Bool Pos_coeff,double lambda,Bool coarse_cons,Bool pos_cons, Bool drop_coarse,Bool iso_cons, double sig)
{
  int Nx = img.nx(), Ny = img.ny(),i,j,k,n;
  if (!coarse_cons) 
    lambda =1;
  if (drop_coarse) 
    lambda = 0;
  Ifloat resi(Nx,Ny);
  Ifloat zeros_mat(Nx,Ny);
  Ifloat mask(Nx,Ny);
  double mu = 1.0;
  if(iso_cons== True)
    {
      double *cent = new double[2];
      int flag = compute_centroid(&img,sig,cent,10);
      for (i=0;i<Nx;i++) 
	for (j=0;j<Ny;j++) 
	  {
	    mask(i,j)=exp((pow(i+1-cent[0],2)+pow(j+1-cent[1],2))/(2*pow(sig/6,2)));
	  }
      free(cent);
      mu = mu/pow(mask.max(),2);
      
      
    }
  for (i=0;i<Nx;i++)
      for (j=0;j<Ny;j++)
	img_filt(i,j) = 0; // The algorithm is initialized at 0
  // Support setting
  
  MultiResol mr_data;
  MultiResol mr_data_filt;
  MultiResol mr_res;
  int flag = wl_trans(&img,opt,&mr_data);
  flag = wl_trans(&img,opt,&mr_res);
  int nb_band = mr_data.nbr_band();
  fltarray supp(Nx,Ny,nb_band);
  for (i=0;i<Nx;i++)
    for (j=0;j<Ny;j++)
      for (k=0;k<nb_band;k++)
	supp(i,j,k) = 1;
  Ifloat suppk(Nx,Ny);
  
  for (k=0;k<nb_band-1;k++)
    {
      Ifloat bandk = mr_data.extract_band (k);
      Ifloat bandk_out(Nx,Ny);
      Ifloat res(Nx,Ny);
      double sig=0;
      mad_res_noise_est(bandk,bandk_out,res,&sig);
      Bool abs_en = (Bool)(!Pos_coeff);
      
      check_ineq(bandk,nsig*sig,abs_en,suppk);
      for (i=0;i<Nx;i++)
	for (j=0;j<Ny;j++)
	  supp(i,j,k) = suppk(j,i);
      bandk.free();
      bandk_out.free();
      res.free();
    }
  
  if (drop_coarse==True or coarse_cons==True)
    {
      for (i=0;i<Nx;i++)
	for (j=0;j<Ny;j++)
	  supp(i,j,nb_band-1) =0;
    }
  
  fltarray supp_cons(Nx,Ny,nb_band);
  for (k=0;k<nb_band-1;k++)
    for (i=0;i<Nx;i++)
      for (j=0;j<Ny;j++)
	supp_cons(i,j,k) = supp(i,j,k);   
  for (i=0;i<Nx;i++)
    for (j=0;j<Ny;j++)
      supp_cons(i,j,nb_band-1) =lambda;

  for (k=0;k<nb_band;k++)
    for (i=0;i<Nx;i++)
      for (j=0;j<Ny;j++)
	mr_data(k,j,i) = mr_data(k,j,i)*supp(i,j,k);  
  
  for (n=0;n<nb_iter;n++)
    {
      if (mse!=NULL)
	 mse[n] = 0;
      wl_trans(&img_filt,opt,&mr_data_filt);
      for (k=0;k<nb_band;k++)
	for (i=0;i<Nx;i++)
	  for (j=0;j<Ny;j++)
	    {
	      mr_res(k,j,i) = mr_data(k,j,i) - mr_data_filt(k,j,i)*supp_cons(i,j,k);
	      if (mse!=NULL)
		mse[n] = mse[n]+ pow(mr_res(k,j,i),2);
	    }
      mr_res.recons(resi,opt.Bord);
      for (i=0;i<Nx;i++)
	for (j=0;j<Ny;j++)
	  img_filt(i,j) = img_filt(i,j)+mu*(resi(i,j)-iso_cons*img_filt(i,j)*pow(mask(i,j),2));
      if (pos_cons==True)
	ineq_cons(img_filt, zeros_mat); 
    }
  resi.free();
  zeros_mat.free();
  supp.free();
  suppk.free();  
  supp_cons.free();
  mr_data.free();
  mr_data_filt.free();
  mr_res.free();
  mask.free();
}

void check_ineq(Ifloat img,double thresh,Bool abs_en,Ifloat &flag)
{
   int Nx = img.nx(), Ny = img.ny(),i,j;
   for (i=0;i<Nx;i++)
    for (j=0;j<Ny;j++)
      {
	flag(i,j)=0;
	if(abs_en==True)
	  {
	    if (fabs(img(i,j))>thresh)
	      flag(i,j)=1;
	  }
	else
	  if (img(i,j)>thresh)
	    flag(i,j)=1;
      }
}

void decim_conv_mat(Ifloat mask,int im_siz[],int decim_fact,Ifloat &output,double flux,double sig) 
{
  int i,j,l,m,Nx=mask.nx(),Ny = mask.ny();
  int mask_rad = (Nx-1)/2; // The kernel mask is supposed to have odd sides lengths
  Ifloat mask_extend(im_siz[0],im_siz[1]);
  Ifloat mask_rot(Nx,Ny);
  rotate(&mask,&mask_rot,2);
  int shift[2];
  int D = decim_fact;
  for (i=0;i<im_siz[0];i++)
    for (j=0;j<im_siz[1];j++)
    {
      shift[0] = i;
      shift[1] = j;
      convmask_shift_extend(mask_rot,im_siz,shift,mask_extend);
      for(l=0;l<floor(im_siz[0]/decim_fact);l++)
	    for(m=0;m<floor(im_siz[1]/decim_fact);m++) 
	      if ((D*l<shift[0]+mask_rad+1) or (D*l>shift[0]-mask_rad-1) or (D*m<shift[1]+mask_rad+1) or (D*m>shift[1]-mask_rad-1))
		{
		  output(j*im_siz[0]+i,m*floor(im_siz[0]/decim_fact)+l) = output(j*im_siz[0]+i,m*floor(im_siz[0]/decim_fact)+l)+mask_extend(D*l,D*m)*flux/sig; // ! output has to be correctly initiliazed 
		}
      
      }
  mask_extend.free();
  mask_rot.free();
}

void convmask_shift_extend(Ifloat mask,int im_siz[],int shift[],Ifloat &output)
{
  
  int Nx = mask.nx(),Ny = mask.ny(),i,j;
  int mask_rad = (Nx-1)/2; // The kernel mask is supposed to have odd sides lengths
  for (i = 0;i<Nx;i++)
    for (j = 0;j<Ny;j++)
      {
	if((i+shift[0]-mask_rad>-1) or (i+shift[0]-mask_rad<im_siz[0]) or (j+shift[1]-mask_rad>-1) or (i+shift[1]-mask_rad<im_siz[1]))
	  output(i+shift[0]-mask_rad,j+shift[1]-mask_rad) = mask(Nx-i,Ny-j);	
      }
    
}

void noise_map_comp(double **shifts, double *sig, double *flux,int nb_im,int lancrad, int im_siz[],mr_opt opt, int decim_fact,fltarray &output)
// Computes a noise standard deviation map
{
  int k,i,j,l,m;
  MultiResol mr_data; 
  int nb_band = mr_data.computeNbBand(im_siz[0],im_siz[1],opt.nb_sc,opt.transf,opt.nb_usc);
  Ifloat decim_conv(im_siz[0]*im_siz[1],im_siz[0]*im_siz[1]/pow((double) decim_fact,(double) 2));
  Ifloat im_temp(im_siz[0],im_siz[1]);
  Ifloat mask(2*lancrad+1,2*lancrad+1);
  double utemp[2];
 
  for(k=0;k<nb_im;k++)
    {
      utemp[0] = shifts[0][k];
      utemp[1] = shifts[1][k];
      lanczos(utemp,mask,lancrad);
      decim_conv_mat(mask,im_siz,decim_fact,decim_conv,flux[k],sig[k]);  
    }
   
  for(l=0;l<im_siz[0]*im_siz[1]/pow((double)decim_fact,(double)2);l++)
    {
      for(i=0;i<im_siz[0];i++)
	for(j=0;j<im_siz[1];j++)
	  {
	    im_temp(i,j) = decim_conv(j*im_siz[0]+i,l);
	  }
      int flag = wl_trans(&im_temp,opt,&mr_data);
      for (m=0;m<nb_band;m++)
	for(i=0;i<im_siz[0];i++)
	  for(j=0;j<im_siz[1];j++)
	    output(i,j,m) = sqrt(pow(output(i,j,m),2)+pow(mr_data(m,j,i),2));
      }
  decim_conv.free();
  mask.free();
  im_temp.free();
  mr_data.free();
}

int renewSeed()
{
  FILE *file = fopen("/dev/urandom","r");   
  int seed;
  fread(&seed, sizeof(int), 1, file);
  fclose(file);
  //printf("Random generator seed = %d\n", seed);
  return seed;
}

void mr_norm(mr_opt opt,double*nrm,int nb_band,int Nx,int Ny)
{
  if (Ny==-1)
    Ny = Nx;
  Ifloat dirac(Nx,Ny);
  dirac(floor(Nx/2),floor(Ny/2))=1;
  MultiResol MR_Data;
  int flag = wl_trans(&dirac,opt,&MR_Data);
  int k;
  for (k=0;k<nb_band;k++)
    {
      Ifloat band = MR_Data.extract_band (k);
      nrm[k] = sqrt(band.energy());
    }
  MR_Data.free();
}

void shell(unsigned long n,double* a)
/* 
   Sort an array a[1..n] into ascending order by Shell's method (diminishing increment sort). n is input; a s replaced on output by its sorted rarrangement.
 */
{
  unsigned long i,j,inc;
  double v;
  inc=1; // Determine the starting increment.
  do {
    inc *=3;
    inc++;
  } while (inc<=n);
  do {   // Loop over the partial sorts.
    inc/=3;
    for (i=inc;i<n;i++) {  // Outer loop of straight insertion.
      v = a[i];
      j=i;
      while (a[j-inc]>v) {    // Inner loop of straight insertion.
	a[j] = a[j-inc];
	j -= inc;
	if (j< inc) break;
	}
      a[j]=v;
    }
  } while (inc>1);
}

double median(double *arr,int n)
{
  shell(n,arr);
  double med = arr[n/2];
  return med;
}


/****************************************************************************/

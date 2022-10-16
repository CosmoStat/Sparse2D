/******************************************************************************
**                   Copyright (C) 1997 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Yanling Fang && Yves Bobichon && Jean-Luc Starck
**
**    Date:  97/12/04
**    
**    File:  MR_RayFilter.cc
**
*******************************************************************************
**
**    DESCRIPTION  filter an image with Rayleigh noise
**    ----------- 
**                 
**    PARAMETRES    
**    ----------    
** 
**    USAGE: mr_rfilter option image output
**    
**  
**
**
**
******************************************************************************/
 
#include "MR_Obj.h"
#include "MR_Rayleigh.h"
#include "MR_Filter.h"

static void mr_sthreshold1 (MultiResol &MR_Data, fltarray &Threshold_N,
                            fltarray &Threshold_P);
static void mr_sthreshold2 (MultiResol &MR_Data, fltarray &Threshold_N1,
                            fltarray &Threshold_N2,
		            fltarray &Threshold_P1 , fltarray &Threshold_P2);
static int test_index(int position, int shift, int dimension);


		     
void mr_sfilter(Ifloat &Imag, Ifloat &Result, StatRayleigh & NoiseModel, 
                fltarray & Tab_Epsilon_Threshold,
		float Epsilon_Ratio, int MaxIter, float Sigma_Converge)

{
   int i,j,iteration;
   int Nl = Imag.nl();
   int Nc = Imag.nc();
   int Nbr_Plan = NoiseModel.nbr_scale();
   int NbrBand = NoiseModel.nbr_band();
   type_noise_speckle  Stat_Noise=NoiseModel.which_noise();
   type_transform Transform = NoiseModel.type_trans();
   Ifloat  Residual(Nl, Nc, "residual"); 
   fltarray TN1,TN2,TP1,TP2;
   float Sigma,Delta_Sigma,Ref_Sigma;
   float Old_Sigma=0;
   MultiResol MR_Ima (Nl, Nc, Nbr_Plan, Transform, "MR_Ima");
   TN1.alloc(NbrBand-1);  
   TP1.alloc(NbrBand-1); 
         
 
   if (NoiseModel.Verbose == True) 
     cout  << "compute thresholds min max for each scale " << endl;
   NoiseModel.find_threshold(NbrBand, Tab_Epsilon_Threshold, TN1, TP1);
   if (Epsilon_Ratio!=1.)
   {
       TN2.alloc(NbrBand-1);       
       TP2.alloc(NbrBand-1);  
       for (int s = 0; s < NbrBand-1; s++) Tab_Epsilon_Threshold(s)/=Epsilon_Ratio;
       NoiseModel.find_threshold(NbrBand, Tab_Epsilon_Threshold, TN2, TP2);
   }

   if (NoiseModel.Verbose == True)
   {
      for (int s = 0; s < NbrBand-1; s++)
      {
	if (Epsilon_Ratio==1.) cout << "band " << s << " : " << "Epsilon= "<< Tab_Epsilon_Threshold(s) << " n= " << TN1(s) << "... p = " << TP1(s) << endl;
	else cout << "band " << s << ": " << "[Eps.1;Eps.2]=["<< Tab_Epsilon_Threshold(s)*Epsilon_Ratio << ","<< Tab_Epsilon_Threshold(s) << "]->[n1,n2]=[" << TN1(s) << "," << TN2(s) << "]  [p1,p2]=[" << TP1(s) << "," << TP2(s) << "]" << endl;
      }
   }

	INFO_X(Imag, "IMAG000");

   switch (Stat_Noise)
   {
      case NOISE_RAYLEIGH:
   	  // compute initial normalisation image -> last smooth image 
   	  MR_Ima.transform(Imag);
  	  Result = MR_Ima.band(NbrBand-1);
  	  if (Transform !=  TO_PAVE_BSPLINE)
   	  {
   	    float F1 = flux(Imag);
    	    float F2 = F1 / flux(Result);
    	    for (i = 0; i<Nl; i++)
    	    for (j = 0; j<Nc; j++)  Result(i,j) *= F2;
  	  }
  	  iteration=0;
  	  Delta_Sigma=10.;
  	  while ((ABS(Delta_Sigma) >  Sigma_Converge)&&(iteration<MaxIter))
  	  {
	     // compute normalised image -> residual
	     for (i = 0; i<Nl; i++)
	     for (j = 0; j<Nc; j++) 
	     {
	        if (Result(i,j) > FLOAT_EPSILON) Residual (i,j) = Imag(i,j) / Result(i,j) ;
	        else 
		{
		  int shift;
		  float Mean=0.;
		  for (shift = 1; Mean == 0.; shift++)
		    Mean = 0.5*(Result(i, test_index (j, -shift, Nc)) + Result(i, test_index (j, shift, Nc)));
		  Residual(i,j)= Imag(i,j) / Mean;
		}
	     }
	     Sigma=sigma(Residual);

	    // compute wavelet transform of Residual
	    MR_Ima.transform(Residual);

	    // threshold the wavelet transform of Residual
            if (Epsilon_Ratio==1.) mr_sthreshold1 (MR_Ima, TN1, TP1);
	    else mr_sthreshold2 (MR_Ima, TN1, TN2, TP1,TP2);

	    // residual reconstruction
	    MR_Ima.recons(Residual);

	    // compute significant residual
	    Residual *= Result; 

	   // filtered residual becomes the new image of normalisation 
           Result = Residual;

	   // next iteration
	   Delta_Sigma=Sigma-Old_Sigma;
	   Old_Sigma=Sigma;
	   iteration++;
	   if (NoiseModel.Verbose == True)
	   {
	    if (iteration==1)
	        cout << "iteration " << iteration << " -> Sigma=" << Sigma  <<  endl; 
	    else
	        cout << "iteration " << iteration << " -> Sigma=" << Sigma << " Delta Sigma = "<< Delta_Sigma <<  endl; 
	   }
        }
        break;
    case NOISE_LOG_RAYLEIGH:
        // compute log of input image
       for (i = 0; i<Nl; i++)
       for (j = 0; j<Nc; j++)
       {
	  if ( Imag(i,j) > FLOAT_EPSILON) Imag(i,j) = log(Imag(i,j));
	  else  
	  {
	      int shift;
	      float Mean=-1.;
	      for (shift = 1; Mean <= 0; shift++)
		  Mean = 0.5*(Imag(i,test_index (j, -shift, Nc) ) + Imag(i,test_index (j, shift, Nc) ));
	      Imag(i,j)=log (Mean);
	  }
	}
    
       // save original image in residual
       // !!! if NOISE_LOG_RAYLEIGH residual = log(Imag) !!!
       Residual = Imag;
       Result.init();
       iteration=0;
       Delta_Sigma=10.;
       Ref_Sigma=(float)PI*sqrt(1./24.);

      while ((ABS(Delta_Sigma) > Sigma_Converge)&&(iteration<MaxIter))
      {
	// compute wavelet transform
 	MR_Ima.transform(Residual);
	
	// threshold the wavelet transform
	if (Epsilon_Ratio==1.) mr_sthreshold1 (MR_Ima, TN1, TP1);
	else mr_sthreshold2 (MR_Ima, TN1, TN2, TP1,TP2);
	
	// image reconstruction
	MR_Ima.recons(Residual);
 	// compute new solution
	Result += Residual;
 	for (i = 0; i<Nl; i++) 
	for (j = 0; j<Nc; j++) Residual(i,j) = Imag(i,j) - Result(i,j);
	Sigma = Residual.sigma();
	Delta_Sigma = Sigma - Ref_Sigma;
	iteration++;
	if (NoiseModel.Verbose == True)  
	 cout << "iteration " << iteration << " -> Sigma=" << Sigma << " Sigma-Sigma_Ref = "<< Delta_Sigma <<  endl; 
    }
    
    // compute the final image -> exp 
    for (i = 0; i<Nl; i++) 
    for (j = 0; j<Nc; j++) Result(i,j) = exp(Result(i,j));
    break;
  case NOISE_LAPLACE:
    Imag *=Imag ;
    
    // compute initial normalisation image -> last smooth image 
    MR_Ima.transform(Imag);
    Result = MR_Ima.band(NbrBand-1);
    
    iteration=0;
    Delta_Sigma=10.;
    Ref_Sigma=1.;
    while ((ABS(Delta_Sigma) >  Sigma_Converge)&&(iteration<MaxIter))
      {
	// compute normalised image -> residual
	for (i = 0; i<Nl; i++)
	for (j = 0; j<Nc; j++) 
	{
	      if (Result(i,j) > FLOAT_EPSILON) Residual (i,j) = Imag(i,j) / Result(i,j);
	      else 
              {
		  int shift;
		  float Mean=0.;
		  for (shift = 1; Mean == 0.; shift++)
		    Mean = 0.5*(Result(i, test_index (j, -shift, Nc)) + Result(i, test_index (j, shift, Nc)));
		  Residual(i,j)= Imag(i,j) / Mean;
	      }
	 }
	Sigma=sigma(Residual);
		
	// compute wavelet transform of Residual
	MR_Ima.transform(Residual);
	
	// threshold the wavelet transform of Residual
	if (Epsilon_Ratio==1.) mr_sthreshold1 (MR_Ima, TN1, TP1);
	else mr_sthreshold2 (MR_Ima, TN1, TN2, TP1,TP2);
	
	// image reconstruction
	MR_Ima.recons(Residual);
	
	// compute  significant residual
	Residual *= Result; 
	
	// filtered residual becomes the new image of normalisation 
        Result = Residual;
	
	// next iteration
	Delta_Sigma=Sigma-Ref_Sigma;
	iteration++;
	if (NoiseModel.Verbose == True) cout << "iteration " << iteration << " -> Sigma=" << Sigma << " Sigma-Sigma_Ref = "<< Delta_Sigma <<  endl; 
      }
    
    // compute the final image -> sqrt
    for (i = 0; i<Nl; i++)
    for (j = 0; j<Nc; j++) 
    {
       if (Result(i,j)>0.) Result (i,j) = sqrt (Result(i,j)) ;
       else 
       {
	  int shift;
	  float Mean=-1.;
	  for (shift = 1; Mean < 0; shift++)
		Mean = 0.5*(Result (i,test_index (j, -shift, Nc) ) + Result (i,test_index (j, shift, Nc) ));
	  Result (i,j) = sqrt(Mean);
       }
    }
    break;
  default:
    cerr << "unknown Noise type" << endl;
  }
}



/*********************************************************************/

static void mr_sthreshold1 (MultiResol &MR_Data, fltarray &Threshold_N, fltarray &Threshold_P)
{
  int s,i,j;
  int Nbr_Plan = MR_Data.nbr_band();
  int Nl=MR_Data.size_ima_nl();
  int Nc=MR_Data.size_ima_nc();
  
  for (s = 0; s < (Nbr_Plan-1); s++)
    {      
      for (i = 0; i < Nl; i++)
	for (j = 0; j < Nc; j++)
	  if ( ( MR_Data(s,i,j) > Threshold_N(s)) && ( MR_Data(s,i,j)< Threshold_P(s)) ) MR_Data(s,i,j) = 0.;
    }
  
}

/***********************************************************************/


static void mr_sthreshold2 (MultiResol &MR_Data, fltarray &Threshold_N1, fltarray &Threshold_N2,
		     fltarray &Threshold_P1 , fltarray &Threshold_P2)
{
  int s,i,j;
  int FS=0;
  int NS=1;
  int Nbr_Plan = MR_Data.nbr_band();
  type_transform Transform = MR_Data.Type_Transform;
  if (Transform == TO_UNDECIMATED_MALLAT) 
  {
     FS = 2;
     NS = 3;
  }
  int Nl=MR_Data.size_ima_nl();
  int Nc=MR_Data.size_ima_nc();
  for (s = 0; s < (Nbr_Plan-1); s++)
    {    
      if (s <= FS)
	{
	  for (i = 0; i < Nl; i++)
	    for (j = 0; j < Nc; j++)
	      {
		// if w(0,i,j) is not significant : w(0,i,j) -> 0
		if ((MR_Data(s,i,j)>Threshold_N1(s)) && (MR_Data(s,i,j)<Threshold_P1(s))) MR_Data(s,i,j) = 0.;
		else 
		  {
		    // else if w(0.i,j) is signicant : look if the corresponding coeff w(1,i,j) is significant
		    if ( ( MR_Data(s+NS,i,j) < Threshold_N1(s+NS)) || ( MR_Data(s+NS,i,j) > Threshold_P1(s+NS))) 
		      {
			// if w(1,i,j) is also significant : apply soft thresholding for w(0,i,j)
			if ((MR_Data(s,i,j) > Threshold_P1(s)) && ( MR_Data(s,i,j)< Threshold_P2(s)))
			  // if w(0,i,j) is positive
			  MR_Data(s,i,j)*=(MR_Data(s,i,j) - Threshold_P1(s))/(Threshold_P2(s)-Threshold_P1(s));
			else
			  {
			    // if w(0,i,j) is negative
			    if ((MR_Data(s,i,j)<Threshold_N1(s)) && (MR_Data(s,i,j)>Threshold_N2(s)))
			      MR_Data(s,i,j)*=(MR_Data(s,i,j)-Threshold_N1(s))/(Threshold_N2(s)-Threshold_N1(s));
			  }
		      }
		    else 
		      {
			//if w(1,i,j) is not significant we don't keep w(0,i,j) -> 0
			MR_Data(s,i,j) = 0.;
		      }
		  }
	      }
	}
      else 
	{
	  // for other scales (s>0) apply standard soft thresholding
	  for (i = 0; i < Nl; i++)
	    for (j = 0; j < Nc; j++)
	      {
		if ( ( MR_Data(s,i,j) > Threshold_N1(s)) && ( MR_Data(s,i,j)< Threshold_P1(s)) )
		  MR_Data(s,i,j) = 0.;
		else 
		  {
		    if ( ( MR_Data(s,i,j) > Threshold_P1(s)) && ( MR_Data(s,i,j)< Threshold_P2(s)) )
		      MR_Data(s,i,j) *=  (MR_Data(s,i,j) - Threshold_P1(s)) / (Threshold_P2(s) - Threshold_P1(s));
		    else
		      {
			if ( ( MR_Data(s,i,j) < Threshold_N1(s)) && ( MR_Data(s,i,j) > Threshold_N2(s)) )
			  MR_Data(s,i,j) *=  (MR_Data(s,i,j) - Threshold_N1(s)) / (Threshold_N2(s) - Threshold_N1(s));
		      }
		  }
	      }
	}
    }
}
/***********************************************************************/

static int test_index(int position, int shift, int dimension)
{
  if ( (position+shift)<0 ) return (position-shift);
  
  else
  {
    if ((position+shift)>=dimension)
    {
      return (position-shift);
    }
    else
    {
      return position+shift;
    }
  }
}

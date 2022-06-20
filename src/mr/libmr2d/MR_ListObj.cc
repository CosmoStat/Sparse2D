
#include "MR_ListObj.h"

#define DEBUG_OBJ 0
 
/****************************************************************************/

ListObj2D::ListObj2D ()
{
      Nb_TotalObj = 0;
      Nb_TotalSubObj = 0;
      NameImagIn = NULL;
      Nl = Nc = 0;
      KeepIma = False;
      ReconsMethod = DEF_RECOBJ_METHOD;
      Nb_iter_rec = DEF_MAX_ITER_RECOBJ;
      ErrorRec = DEF_ERROR_RECOBJ;
      TabObj = NULL;
      FirstUseScale=1;
      LastUseScale=-1;
      Verbose = False;
      FluxMult = 1.;
      Deconv = False;
      Fwhm = -1;
      ZoomPSF = 1.;
      VisionModel=MVM_RUEBIJ;
      F_obj = NULL;
      MVM = NULL;
      IsotropSeparation = False;
      AnitropicSeparationParam=1.5;
}

/****************************************************************************/

ListObj2D::~ListObj2D()
{
      // cout << "Begin deallocation ListObj2D" << endl;
      if (TabObj != NULL) delete [] TabObj;
      if (NameImagIn != NULL) delete NameImagIn;
      Nb_TotalObj = 0;
      NameImagIn = NULL;
      Nl = Nc = 0;
      KeepIma = False;
      ReconsMethod = DEF_RECOBJ_METHOD;
      Nb_iter_rec = DEF_MAX_ITER_RECOBJ;
      ErrorRec = DEF_ERROR_RECOBJ;
      TabObj = NULL;
      FirstUseScale=1;
      LastUseScale=-1;
      Verbose = False;
      // cout << "End deallocation ListObj2D" << endl;
}

/****************************************************************************/

static int ind_ok(int ind, int N)
{
  int Val = 1;
  if ((ind < 0) || (ind >= N)) Val = 0;
  return Val;
}
/****************************************************************************/

static void box_obj_sig_and_noise(Ifloat &Im_rec, MRNoiseModel & ModelData,
                       t_Gauss_2D &Gauss, int Nlo, int Nco, int Nl, int Nc,
                       int dx, int dy, float Flux, int & Sizebox,
                       float & SigBox, float & SigNoise, float &FluxBox)
{
    double Val;
    int k, b=0;
    int Xmax = (int) (Gauss.mx + 0.5);
    int Ymax = (int) (Gauss.my + 0.5);
     
    FluxBox = Im_rec(Ymax, Xmax); 
    SigBox = FluxBox*FluxBox;
    SigNoise=0.;
    if (ModelData.UseRmsMap == True) 
       SigNoise += ModelData.RmsMap(Ymax+dy,Xmax+dx, I_MIRROR)
                          * ModelData.RmsMap(Ymax+dy,Xmax+dx, I_MIRROR);

#if DEBUG_OBJ
   cout << "box_obj_sig_and_noise" << endl;
    cout << "   Nlo = " << Nlo << " Nco = " << Nco << endl;
    cout << "   Ymax  = " << Ymax  << " Xmax   = " << Xmax  << endl;
    cout << "   Flux = " << Flux << endl;
#endif
   
   do  
   {
      b++;
      for (k = -b; k <= b; k++)
      {
         if (ind_ok(Ymax-b,Nlo) && ind_ok(Xmax+k, Nco))
	 { 
            Val =  Im_rec(Ymax-b, Xmax+k);
            FluxBox += Val;
            SigBox += Val*Val;

            if (ModelData.UseRmsMap == True)
            {
                Val = ModelData.RmsMap(Ymax-b+dy,Xmax+dx+k, I_MIRROR);
                SigNoise +=Val*Val;
            }
	 }
	 if (ind_ok(Ymax+b,Nlo) && ind_ok(Xmax+k, Nco))
 	 {
            Val = Im_rec(Ymax+b, Xmax+k);
            FluxBox += Val;
            SigBox += Val*Val;
            if (ModelData.UseRmsMap == True)
  	    {
               Val = ModelData.RmsMap(Ymax+b+dy,Xmax+dx+k, I_MIRROR);
               SigNoise +=Val*Val;
	    }
	 }
      }
      for (k = -b + 1; k <= b - 1; k++)
      {
         if (ind_ok(Xmax-b,Nco) && ind_ok(Ymax+k,Nlo))
 	 { 
            Val = Im_rec(Ymax+k, Xmax-b);
            FluxBox += Val;
            SigBox += Val*Val;
            if (ModelData.UseRmsMap == True) 
 	    {
                Val = ModelData.RmsMap(Ymax+k+dy,Xmax+dx-b, I_MIRROR);
                SigNoise += Val*Val;
	    }
	 }
	 if (ind_ok(Xmax+b,Nco) && ind_ok(Ymax+k,Nlo))
 	 { 
            Val = Im_rec(Ymax+k, Xmax+b);
            FluxBox += Val;
            SigBox += Val*Val;
            if (ModelData.UseRmsMap == True) 
 	    {
               Val = ModelData.RmsMap(Ymax+k+dy,Xmax+dx+b, I_MIRROR);
               SigNoise +=Val*Val;
            }
	 }
      } // end for
   } while (FluxBox < 0.9 * Flux);

#if DEBUG_OBJ
    cout << "b = " << b  << endl;
#endif

   Sizebox = 2*b+1;
   float S2 =  Sizebox*Sizebox;
   SigBox = SigBox / S2 -  (FluxBox/S2*FluxBox/S2);
   if (SigBox >= 0.) SigBox = sqrt(SigBox);
   else SigBox = 0.;
   
   if (ModelData.UseRmsMap == True) SigNoise = sqrt(SigNoise)/(float) S2;
}
       
/****************************************************************************/

double ListObj2D::get_aperture_flux_obj(Ifloat &Data,  
      double  PosX,  double PosY, double SigmaX, double SigmaY)
{
   int Nl = Data.nl();
   int Nc = Data.nc();
   float MaxSigma = MAX(SigmaX,SigmaY);
   // cout << "MaxSigma " << MaxSigma <<  endl;
   int i,j;
   int Depx = (int) (PosX-KSigmaAper*MaxSigma);
   int Depy = (int) (PosY-KSigmaAper*MaxSigma);
   int Endx = (int) (PosX+KSigmaAper*MaxSigma);
   int Endy = (int) (PosY+KSigmaAper*MaxSigma);
   double  Flux=0.;
   if (Depx < 0) Depx = 0;
   if (Depy < 0) Depy = 0;
   if (Endx >= Nc) Endx = Nc-1;
   if (Endy >= Nl) Endy = Nl-1;

#if DEBUG_OBJ
 cout << "Depx = " << Depx << " Depy = " << Depy << " Endx = " << Endx  << " Endy = " << Endy << endl; 
#endif

   if (BgrMethod == BGR_VALUE)
   {
      for (i=Depy; i <= Endy; i++)
      for (j=Depx; j <= Endx; j++) Flux += Data(i,j) -  BackgroundValue;
   }
   else
   {
      for (i=Depy; i <= Endy; i++)
      for (j=Depx; j <= Endx; j++) Flux += Data(i,j) - BgrModel(i,j);
   }
#if DEBUG_OBJ
 cout << " Flux aperture = " << Flux << endl;
#endif
   return Flux;
}

/****************************************************************************/

void ListObj2D::init_background(Ifloat &Data)
{
   switch (BgrMethod)
   {
      case BGR_VALUE:
            if (BackgroundValue < 0) BackgroundValue = average(Data);
	    break;
      case BGR_IMA:
            if ((BgrModel.nl() != Data.nl()) || (BgrModel.nc() != Data.nc()))
	    {
	        cerr << "Error: background image must have the same size as the input data ... " << endl;
		exit(-1);
	    }
	    break;
      case BGR_MR:
           {
	      // Background Estimation from the  wavelet transform
              int Nl = Data.nl();
              int Nc = Data.nc();              
	      int Nbr_Plan=1;
              int Ns = MIN(Nl,Nc);
              while (BGRSize < Ns) 
              {
                 Nbr_Plan ++;
                 Ns /= 2;
              }
              if (Nbr_Plan < 3) Nbr_Plan = 3;
              if (Verbose == True)
                cout << "Background estimation: number of scales = " << Nbr_Plan << endl;
              MultiResol MR_Data (Nl, Nc, Nbr_Plan, TM_TO_PYR, "MR_Transform");
              MR_Data.Border=I_MIRROR;
              MR_Data.transform (Data);
              for (int b=0; b < Nbr_Plan-1; b++) MR_Data.band(b).init();
	      BgrModel.alloc(Nl, Nc,"model");
              MR_Data.recons(BgrModel);
              if (Verbose == True) io_write_ima_float("xx_bgr.fits", BgrModel);
	   }
	   break;
      default: 
           cerr << "Error: unknown background estimation method ... " << endl;
	   exit(-1);
	   break;   
   }
}

/****************************************************************************/

void ListObj2D::create_list(Ifloat &Data,  MRNoiseModel & ModelData, 
                            MultiResol & MR_Data, Bool WriteAllObj, 
                            Bool InfoSubObj, Bool WriteObjFullSize)
{
    int i,j,k,l,Nb_obj,dx,dy;
    MultiResol W0;
    double Flux;
    t_Gauss_2D Gauss;
    double Eps_ErrorRec;
    int Num_inTotal=0;
    char NameObj[80];
    int Nb_ech = ModelData.nbr_band()-1;
    int AddBorderX=0;
    int AddBorderY=0;
    Bool ApplyDeconv = False;
    Ifloat ImaObjFullSize;
    
    if ((VisionModel == MVM_RUEBIJ) && (F_obj == NULL))
    {
       cout << "Error: F_obj is not allocated in create_list ... " << endl;
       exit(-1);
    }
    else if ((VisionModel != MVM_RUEBIJ) && (MVM == NULL))
    {
       cout << "Error:  MVM is not allocated in create_list ... " << endl;
       exit(-1);
    }

    // When a deconvolution is performed, the image must be bigger
    // in order to take into account PSF effect
    if (ReconsMethod == GRAD_PSF)
    {
       AddBorderY = Psf.nl();
       AddBorderX = Psf.nc();
       ApplyDeconv = True;
    }
    else if (ReconsMethod == GRAD_PSF_XMM)
    {
       AddBorderY = (int) (XMM_PSF_NL / ZoomPSF);
       AddBorderX = (int) (XMM_PSF_NC / ZoomPSF);
       ApplyDeconv = True;
    }
    
    if ((Fwhm > 0) && (ReconsMethod == GRAD_PSF))
    {
       IMGauss = im_gaussian(Psf.nl(), Psf.nc(), Fwhm*ZoomPSF);
       norm_flux(IMGauss);
       psf_convol (Psf, IMGauss);
    }
    if (AperPhot == True) init_background(Data);

   // verify LastUseScale and FirstUseScale
   if (LastUseScale < 0) LastUseScale = Nb_ech-1;
   else if (LastUseScale >= Nb_ech) LastUseScale = Nb_ech-1;
   if (FirstUseScale < 0) FirstUseScale = 0;
   if (FirstUseScale > LastUseScale) FirstUseScale = LastUseScale;

   // count the total number of objects
   Nb_TotalObj = 0;
   Nl = Data.nl();
   Nc = Data.nc();
   if (WriteObjFullSize == True) ImaObjFullSize.alloc(Nl,Nc,"ImaObjFullSize");
   if (IsotropSeparation == True) Ima_SumObjAnisotopic.alloc(Nl,Nc,"ImaAnisotropObjFullSize");
   
   if (VisionModel == MVM_RUEBIJ) 
   {
      for (i= LastUseScale ; i>= FirstUseScale ; i--) 
                        Nb_TotalObj += (F_obj->Nb_arbre)[i];
   }
   else Nb_TotalObj = MVM->nbr_obj();

   // initialization
   TabObj = new Object_2D[Nb_TotalObj];
   Ima_SumObj.alloc(Nl, Nc, "Objects sum");
   Ifloat Im_rec(Nl, Nc, "Reconstructed_object");

if (Verbose == True)
{
   cout << "Total number of objects: " << Nb_TotalObj << endl;
   cout << "FirstUseScale = " << FirstUseScale + 1<< endl;
   cout << "LastUseScale  = " << LastUseScale + 1 << endl;
}
    for (i= LastUseScale ; i>= FirstUseScale ; i--)
    {
      info_obj2d *CurrObj=NULL;
      Nb_obj = (VisionModel == MVM_RUEBIJ) ? 
             (F_obj->Nb_arbre)[i] :  (MVM->nbr_obj_in_scale)(i);

      for(j=0; j < Nb_obj; j++)
      {
	  Arbre *TreeObj=NULL;
          float PosMaxCoefX=0,PosMaxCoefY=0;
          Bool RecObj = False;
          int Nlo,Nco;
          if (VisionModel == MVM_RUEBIJ)
          {
  	      TreeObj = F_obj->get_obj (i, j);
              if ((!TreeObj->S_obj) || (InfoSubObj == True)) 
              {
                 RecObj = True;  
 	         F_obj->creat_wavelet (W0, i,j, dx, dy, AddBorderX, AddBorderY);
              }
 	      PosMaxCoefX = TreeObj->Xmax+TreeObj->Xpos;
	      PosMaxCoefY = TreeObj->Ymax+TreeObj->Ypos;
          }
          else
          {
               CurrObj = MVM->get_obj_info(i,j+1);
	       dx = CurrObj->DepImaX;
	       dy = CurrObj->DepImaY;
  	       if ((CurrObj->SubObj == False) || (InfoSubObj == True)) 
               {
                   RecObj = True;
                   MVM->get_mrobj(i, j+1, MR_Data, W0);
               }
	  }

          if (RecObj == True)
	  {
              Nlo = W0.size_ima_nl();
	      Nco = W0.size_ima_nc();
          // cout << "create obj: " << i << " " <<  j << " "<< dx << " " << dy << endl;
          //cout << "echelle ancetre = " << TreeObj->ancetre()->Num_ech << endl;

	  Im_rec.resize(Nlo, Nco);
          // Im_rec.init();
	  
	  // In case of deconvolution, we need the input data corresponding
	  // to the object aera
	  if ((ApplyDeconv == True) || (ReconsMethod == GRAD_FIX_STEP))
	  {
  	     for (k = 0; k < Nlo; k++)
	     for (l = 0; l < Nco; l++)
	     {
	        if ((k+dy < Data.nl()) &&  (l+dx < Data.nc())
		     && (k+dy >= 0) &&  (l+dx >= 0))
 	             Im_rec(k,l) = Data(k+dy,l+dx);
	        else Im_rec(k,l) = 0.;
	     }
	  }
	  
 	 // Reconstruction

#if DEBUG_OBJ	 
	 cout << "create obj: " << i << " " <<  j << " "<< dx << " " << dy << endl;
	 cout << "PosMaxCoefX = " << PosMaxCoefX << "  PosMaxCoefY = " << PosMaxCoefY << endl;
#endif
	 if (VisionModel == MVM_RUEBIJ)
               recons_obj(W0, Im_rec, Eps_ErrorRec, Flux, Gauss, 
	                  PosMaxCoefX, PosMaxCoefY);
 	 else recons_obj(W0, Im_rec, Eps_ErrorRec, Flux, Gauss);

	 if (Flux < FLOAT_EPSILON)
	 {
	    cerr << "Error: reconstructed object flux = 0 ... " << endl;
	    W0.write("xx_err.mr");
	    exit(-1);
	 }
  	 W0.free();
#if DEBUG_OBJ		    
	 cout << "Flux rec = " << Flux << endl;	    
#endif
  
	Eps_ErrorRec = Eps_ErrorRec*Flux;
 	TabObj[Num_inTotal].NumObj = Num_inTotal+1;
        if (VisionModel == MVM_RUEBIJ)
        {
	   TabObj[Num_inTotal].NumScale = i+1;
	   TabObj[Num_inTotal].NumObjScale = j+1;
	   TabObj[Num_inTotal].PosMaxCoef_X = TreeObj->Xmax+TreeObj->Xpos;
	   TabObj[Num_inTotal].PosMaxCoef_Y = TreeObj->Ymax+TreeObj->Ypos;
	   TabObj[Num_inTotal].ValMaxCoef = TreeObj->Val_max;
	}
        else
        {
 	   TabObj[Num_inTotal].NumScale = CurrObj->ScaleMax+1;
	   TabObj[Num_inTotal].NumObjScale =  CurrObj->ObjNumberInScale;
	   TabObj[Num_inTotal].PosMaxCoef_X = CurrObj->PosMaxImaX;
	   TabObj[Num_inTotal].PosMaxCoef_Y = CurrObj->PosMaxImaY;
	   TabObj[Num_inTotal].ValMaxCoef = CurrObj->Max;
        }

  	if (ModelData.which_noise() != NOISE_EVENT_POISSON)
	   TabObj[Num_inTotal].SNR_ValMaxCoef = 
              ABS(TabObj[Num_inTotal].ValMaxCoef) 
                /    ModelData.sigma(i,TabObj[Num_inTotal].PosMaxCoef_Y,
                                         TabObj[Num_inTotal].PosMaxCoef_X);
        else 
        {
           TabObj[Num_inTotal].SNR_ValMaxCoef = 
                 ModelData.prob_signal(TabObj[Num_inTotal].ValMaxCoef, i,
                                         TabObj[Num_inTotal].PosMaxCoef_Y,
                                         TabObj[Num_inTotal].PosMaxCoef_X);
        }
                                               
	TabObj[Num_inTotal].Nlo = Nlo;
	TabObj[Num_inTotal].Nco = Nco;
	TabObj[Num_inTotal].Nli = Nl;
	TabObj[Num_inTotal].Nci = Nc;
	TabObj[Num_inTotal].SigmaX = Gauss.sX;
	TabObj[Num_inTotal].SigmaY = Gauss.sY;
	TabObj[Num_inTotal].PosX = Gauss.mx + dx;
	TabObj[Num_inTotal].PosY = Gauss.my + dy;
	TabObj[Num_inTotal].DepX = dx;
	TabObj[Num_inTotal].DepY = dy;
	TabObj[Num_inTotal].Angle = Gauss.teta*RAD_DEG;
	if ( (Gauss.sX > AnitropicSeparationParam * Gauss.sY) 
	   || (Gauss.sY > AnitropicSeparationParam * Gauss.sX)) TabObj[Num_inTotal].Isotrop = False;
	
	TabObj[Num_inTotal].FluxObjRec = Flux;
	if (AperPhot == True)
	    Flux = get_aperture_flux_obj(Data, 
	                TabObj[Num_inTotal].PosX,TabObj[Num_inTotal].PosY,
			TabObj[Num_inTotal].SigmaX,TabObj[Num_inTotal].SigmaY);
	TabObj[Num_inTotal].Flux = Flux;
       
	TabObj[Num_inTotal].Magnitude = 
                   (Flux <= 0) ? 0. :-2.5 * (float) log10((double) Flux);
	TabObj[Num_inTotal].ValPixMax = Gauss.amp;
	TabObj[Num_inTotal].ErrorRec = Eps_ErrorRec;

	if (VisionModel == MVM_RUEBIJ)
        {
           if (TreeObj->S_obj) 
	   {
	      TabObj[Num_inTotal].SubObj = True;
	      Nb_TotalSubObj ++;
	   }
        }
	else if (CurrObj->SubObj == True) 
	{
           info_obj2d *FatherObj = MVM->get_obj_info (CurrObj->NumFatherObj);
	   TabObj[Num_inTotal].SubObj = True;
	   TabObj[Num_inTotal].FatherScale = FatherObj->ScaleMax+1;
	   TabObj[Num_inTotal].FatherNumberInScale  = FatherObj->ObjNumberInScale;
	   Nb_TotalSubObj ++;
	}

#if DEBUG_OBJ	 
   cout << "Xmax =  " <<  Gauss.mx+dx << " Ymax = " <<  Gauss.my+dy << endl;
#endif
	// Flux Error calculation
	// search the 90% flux box (Sizebox), the signal standard deviation
        // in this box and the noise standard deviation in this box
	// in the Poisson case, 
	// SNR = 90% Flux Source / sqrt(flux data in the box which contains  
	//                              the 90% of the flux) 
	// (and no SNR = (Nev_box - Nev_bgr) / sqrt(sigma^2_source +sigma^2_noise)
	// in the Gaussian case, SNR = sigma_signal / sigma_noise
	int Sizebox;
	float SigBox, SigNoise, FluxBox;

        // OUTPUT: Sizebox= size of the box containing 90% of the flux
        //         SigBox=standard deviation in this box
        //         SigNoise=standard deviation of the noise in this box
        //                  calculated only if ModelData.UseRmsMap==True
        // FluxBox = flux inside the box
        box_obj_sig_and_noise(Im_rec, ModelData, Gauss,  Nlo, Nco,  Nl,  Nc,
                        dx,  dy,  TabObj[Num_inTotal].FluxObjRec, 
			Sizebox,  SigBox, SigNoise, FluxBox);

#if DEBUG_OBJ
   cout << "Xmax =  " <<  Gauss.mx+dx << " Ymax = " <<  Gauss.my+dy << endl;
#endif

       switch (ModelData.which_noise())
       {
          case NOISE_GAUSSIAN: 
          case NOISE_POISSON:
          case NOISE_GAUSS_POISSON:
          case NOISE_MULTI: 
              // the error is proportional to the number of pixels
	      // of the box which contains 90% of the flux
	      // We consider that 90% of the flux calculation implies
	      // 90% of the error    
              TabObj[Num_inTotal].ErrorFlux = 
                    sqrt((double) Sizebox*Sizebox) * ModelData.SigmaNoise / 0.9;
              
              TabObj[Num_inTotal].ErrorFlux =
               sqrt(TabObj[Num_inTotal].ErrorFlux*TabObj[Num_inTotal].ErrorFlux 
                     + Eps_ErrorRec*Eps_ErrorRec);
              TabObj[Num_inTotal].SNR_Obj = SigBox / ModelData.SigmaNoise;
              break;
          case NOISE_NON_UNI_ADD: 
	   if (ModelData.UseRmsMap == True) 
	   {
 	      if (SigNoise > FLOAT_EPSILON)
                  TabObj[Num_inTotal].SNR_Obj = SigBox / SigNoise;
	      else  TabObj[Num_inTotal].SNR_Obj = 10000.;
	      
	      // previous SigNoise was estimated only in the box 
	      // which contain 90% of the flux. We need here the 
	      // error in the box which contains the object
 	      TabObj[Num_inTotal].ErrorFlux = Sizebox*Sizebox*SigNoise  / 0.9;
	   }
	   else 
	   {
	      TabObj[Num_inTotal].ErrorFlux = 0.;
	      TabObj[Num_inTotal].SNR_Obj = 0.;
	   }
	   break;
        case NOISE_EVENT_POISSON:
             { 
              // the error is proportional to root square of the flux 
              int ii,jj;
              
	      SigNoise = 0.;
	      for (k = Sizebox/2-1; k <= Sizebox/2+1; k++)
	      for (l = Sizebox/2-1; l <= Sizebox/2+1; l++)
	      {   
	         ii = (int) (TabObj[Num_inTotal].PosY + k);
	         jj = (int) (TabObj[Num_inTotal].PosX + l);
	         if ( ind_ok(ii,Nl) && ind_ok(jj,Nc) )
	         {
	            SigNoise += Data(ii, jj);
	         }
	      }
	      SigNoise = sqrt(SigNoise);
	      TabObj[Num_inTotal].ErrorFlux = SigNoise / 0.9;
	      if (SigNoise > 0) 
  	         TabObj[Num_inTotal].SNR_Obj = FluxBox / SigNoise;
  	      else TabObj[Num_inTotal].SNR_Obj = 10000;
 	    }
 	   break;
        case NOISE_UNDEFINED:
        case NOISE_UNI_UNDEFINED:
        case NOISE_NON_UNI_MULT:     
          default:
	      TabObj[Num_inTotal].ErrorFlux = 0.;
	      TabObj[Num_inTotal].SNR_Obj = 0.;
              break;
       }


	// if a transform has been applied to an image, the inverse
	// transform must be performed on the reconstructed object,
	// This inverse transform can not be done in the case of
	// of Poisson noise with few events.
	// Several parameters must be reestimated:
	//   the flux, the error flux,  the moment, ...
	// SNR estimations don't have to be recalculated.
	if ((ModelData.TransImag == True) &&
            (ModelData.which_noise() != NOISE_EVENT_POISSON))
	{
	   float ValIma;
	   float ValImaTrans;

	   double FluxInvTrans=0.;
	   for (k=0;k<Nlo;k++)
	   for (l=0;l<Nco;l++)
	   {
	      ValIma = Data(k+dy,l+dx,I_CONT);
	      ValImaTrans = ModelData.val_transform(ValIma);
	      Im_rec(k,l) = ValIma - 
                     ModelData.val_invtransform( ValImaTrans - Im_rec(k,l));
	      FluxInvTrans += Im_rec(k,l);
	   }
	   Estime_param_gauss_2D(Gauss, Im_rec);
	   TabObj[Num_inTotal].SigmaX = Gauss.sX;
	   TabObj[Num_inTotal].SigmaY = Gauss.sY;
	   TabObj[Num_inTotal].PosX = Gauss.mx + dx;
	   TabObj[Num_inTotal].PosY = Gauss.my + dy;
	   TabObj[Num_inTotal].Angle = Gauss.teta*RAD_DEG;
	   TabObj[Num_inTotal].ErrorFlux = FluxInvTrans * 
                 TabObj[Num_inTotal].ErrorFlux / TabObj[Num_inTotal].Flux;
	   TabObj[Num_inTotal].Flux = FluxInvTrans;
	   TabObj[Num_inTotal].Magnitude = 
                   (FluxInvTrans <= 0) ? 0. :-2.5 * (float) log10((double) FluxInvTrans);
	   TabObj[Num_inTotal].ValPixMax = Gauss.amp;
	   TabObj[Num_inTotal].Surface = Sizebox*Sizebox;
	}
        

if (Verbose == True)
{
cout<< "Object number: " << Num_inTotal+1 << " object:" << i+1 <<","<< j+1 << endl;
if (VisionModel == MVM_RUEBIJ) 
       cout << "   Last Scale = " << TreeObj->ancetre()->Num_ech + 1 << endl;
else cout << "   Last Scale = " <<  CurrObj->RootScale+1 << endl;
cout<< "    Size = "<<"  "<< Nlo << "x" << Nco << endl;
cout<< "    Depx = " << dx << " Depy = " << dy << endl;
cout<< "    PosX = " << Gauss.mx + dx << " PosY = " << Gauss.my + dy << endl;
cout<< "    SigmaX = " << Gauss.sX << " SigmaY = " << Gauss.sY << endl;
cout<< "    Angle = " << Gauss.teta*RAD_DEG << endl;
if (AperPhot == True) 
cout<< "    Aper. Phot. Flux = " << TabObj[Num_inTotal].Flux << endl;	     
cout<< "    Reconstructed object Flux = " << TabObj[Num_inTotal].FluxObjRec << endl;
cout<< "    ErrorFlux = " << TabObj[Num_inTotal].ErrorFlux << endl;
cout<< "    Reconstruction Error " <<  TabObj[Num_inTotal].ErrorRec << endl;
cout<< "    SNR Max Coef= " <<	TabObj[Num_inTotal].SNR_ValMaxCoef << endl; 
cout<< "    Pos Max CoefX= " << TabObj[Num_inTotal].PosMaxCoef_X ;  
cout<< "    Pos Max CoefY= " <<  TabObj[Num_inTotal].PosMaxCoef_Y << endl;
cout<< "    SNR Obj= " <<	TabObj[Num_inTotal].SNR_Obj << endl; 
cout<< "    SizeBox= " << Sizebox << endl;

if (TabObj[Num_inTotal].SubObj == True) cout<< "    Sub-object = True" << endl;
else cout<< "    Sub-object = False" << endl;
} 

	if (WriteAllObj == True)
	{
	    sprintf(NameObj, "ima_obj_%d_%d.fits", i+1, j+1);
            if (WriteObjFullSize == True) 
            {
                ImaObjFullSize.init();
                add_image (ImaObjFullSize, Im_rec, dx, dy);
	        io_write_ima_float(NameObj, ImaObjFullSize);
            }
            else io_write_ima_float(NameObj, Im_rec);
 	}

	// Additionne a l'image reconstruite finale si non sous-objet
	if (TabObj[Num_inTotal].SubObj == False) 
	{
             if ((IsotropSeparation == False) || (TabObj[Num_inTotal].Isotrop == True))
	           add_image (Ima_SumObj, Im_rec, dx, dy);
	     else  add_image (Ima_SumObjAnisotopic, Im_rec, dx, dy);
         }
	
	if (KeepIma == True)
	{
	    TabObj[Num_inTotal].Image->alloc(Nlo, Nco, "TabObj");
	    *(TabObj[Num_inTotal].Image) = Im_rec;
	}
	Num_inTotal++;
     }
     else Nb_TotalSubObj ++;
    }
  }
  Nb_TotalObj=Num_inTotal;
}

/****************************************************************************/
// 
// void ListObj2D::create_list(MR2D_Seg & F_obj,Ifloat &Data,MultiResol & MR_Data, 
// 		   MRNoiseModel & ModelData, Bool WriteAllObj, Bool InfoSubObj)
// {
//     int i,j,k,l,Nb_obj,dx,dy;
//     MultiResol W0;
//     double Flux;
//     t_Gauss_2D Gauss;
//     double Eps_ErrorRec;
//     int Num_inTotal=0;
//     char NameObj[80];
//     int Nb_ech = ModelData.nbr_band()-1;
//     int AddBorderX=0;
//     int AddBorderY=0;
//     Bool ApplyDeconv = False;
// 
//     // When a deconvolution is performed, the image must be bigger
//     // in order to take into account PSF effect
//     if (ReconsMethod == GRAD_PSF)
//     {
//        AddBorderY = Psf.nl();
//        AddBorderX = Psf.nc();
//        ApplyDeconv = True;
//     }
//     else if (ReconsMethod == GRAD_PSF_XMM)
//     {
//        AddBorderY = (int) (XMM_PSF_NL / ZoomPSF);
//        AddBorderX = (int) (XMM_PSF_NC / ZoomPSF);
//        ApplyDeconv = True;
//     }
//     
//     if (Fwhm > 0)
//     {
//        IMGauss = im_gaussian(Psf.nl(), Psf.nc(), Fwhm*ZoomPSF);
//        norm_flux(IMGauss);
//        psf_convol (Psf, IMGauss);
//     }
//     if (AperPhot == True) init_background(Data);
// 
//    // verify LastUseScale and FirstUseScale
//    if (LastUseScale < 0) LastUseScale = Nb_ech-1;
//    else if (LastUseScale >= Nb_ech) LastUseScale = Nb_ech-1;
//    if (FirstUseScale < 0) FirstUseScale = 0;
//    if (FirstUseScale > LastUseScale) FirstUseScale = LastUseScale;
// 
//    // count the total number of objects
//    Nb_TotalObj = F_obj.nbr_obj();
//  
//    // initialization
//    Nl = Data.nl();
//    Nc = Data.nc();
//    TabObj = new Object_2D[Nb_TotalObj];
//    Ima_SumObj.alloc(Nl, Nc, "Objects sum");
//    Ifloat Im_rec(Nl, Nc, "im_rec");
// 
// if (Verbose == True)
// {
//    cout << "Total number of objects: " << Nb_TotalObj << endl;
//    cout << "Number of sub-objects: " << F_obj.nbr_sub_obj() << endl;
//    cout << "FirstUseScale = " << FirstUseScale + 1<< endl;
//    cout << "LastUseScale  = " << LastUseScale + 1 << endl;
// }
//     for (i= LastUseScale ; i>= FirstUseScale ; i--)
//     {
//       info_obj2d *CurrObj;
//       Nb_obj = F_obj.nbr_obj_in_scale(i);
//       for(j=0; j < Nb_obj; j++)
//       {
//           CurrObj = F_obj.get_obj_info(i,j+1);
// 	  dx = CurrObj->DepImaX;
// 	  dy = CurrObj->DepImaY;
//   	  if ((CurrObj->SubObj == False) || (InfoSubObj == True))
// 	  {
// 	  // Construit l'ondelette de l'objet
//           F_obj.get_mrobj(i, j+1, MR_Data, W0);
//      // cout << "create obj: " << i << " " <<  j << " "<< dx << " " << dy << endl;
//     // cout << "echelle ancetre = " << TreeObj->ancetre()->Num_ech << endl;
// 
// 	  int Nlo = W0.size_ima_nl();
// 	  int Nco = W0.size_ima_nc();
// 	  Im_rec.resize(Nlo, Nco);
// 
//           		
// 	// Reconstruction
//  	recons_obj(W0, Im_rec, Eps_ErrorRec, Flux, Gauss);
// 	if (Flux < FLOAT_EPSILON)
// 	{
// 	   cerr << "Error: reconstructed object flux = 0 ... " << endl;
// 	   W0.write("xx_err.mr");
// 	   exit(-1);
// 	}
//   	W0.free();
//   
// 	Eps_ErrorRec = Eps_ErrorRec*Flux;
//  
//  	TabObj[Num_inTotal].NumObj = Num_inTotal+1;
// 	TabObj[Num_inTotal].NumScale = CurrObj->ScaleMax+1;
// 	TabObj[Num_inTotal].NumObjScale =  CurrObj->ObjNumberInScale;
// 
// 	TabObj[Num_inTotal].PosMaxCoef_X = CurrObj->PosMaxImaX;
// 	TabObj[Num_inTotal].PosMaxCoef_Y = CurrObj->PosMaxImaY;
// 	TabObj[Num_inTotal].ValMaxCoef = CurrObj->Max;
//  	
//   	if (ModelData.which_noise() != NOISE_EVENT_POISSON)
// 	   TabObj[Num_inTotal].SNR_ValMaxCoef = 
//               ABS(TabObj[Num_inTotal].ValMaxCoef) 
//                 /    ModelData.sigma(i,TabObj[Num_inTotal].PosMaxCoef_Y,
//                                          TabObj[Num_inTotal].PosMaxCoef_X);
//         else 
//         {
//            TabObj[Num_inTotal].SNR_ValMaxCoef = 
//                  ModelData.prob_signal(TabObj[Num_inTotal].ValMaxCoef, i,
//                                          TabObj[Num_inTotal].PosMaxCoef_Y,
//                                          TabObj[Num_inTotal].PosMaxCoef_X);
//         }
//                                                
// 	TabObj[Num_inTotal].Nlo = Nlo;
// 	TabObj[Num_inTotal].Nco = Nco;
// 	TabObj[Num_inTotal].Nli = Nl;
// 	TabObj[Num_inTotal].Nci = Nc;
// 	TabObj[Num_inTotal].SigmaX = Gauss.sX;
// 	TabObj[Num_inTotal].SigmaY = Gauss.sY;
// 	TabObj[Num_inTotal].PosX = Gauss.mx + dx;
// 	TabObj[Num_inTotal].PosY = Gauss.my + dy;
// 	TabObj[Num_inTotal].DepX = dx;
// 	TabObj[Num_inTotal].DepY = dy;
// 	TabObj[Num_inTotal].Angle = Gauss.teta*RAD_DEG;
//  	TabObj[Num_inTotal].FluxObjRec = Flux;
// 	if (AperPhot == True)
// 	    Flux = get_aperture_flux_obj(Data, 
// 	                TabObj[Num_inTotal].PosX,TabObj[Num_inTotal].PosY,
// 			TabObj[Num_inTotal].SigmaX,TabObj[Num_inTotal].SigmaY);
// 	TabObj[Num_inTotal].Flux = Flux;
// 
// 	TabObj[Num_inTotal].Magnitude = 
//                    (Flux <= 0) ? 0. :-2.5 * (float) log10((double) Flux);
// 	TabObj[Num_inTotal].ValPixMax = Gauss.amp;
// 	TabObj[Num_inTotal].ErrorRec = Eps_ErrorRec;
// 	if (CurrObj->SubObj == True) 
// 	{
//            info_obj2d *FatherObj = F_obj.get_obj_info(CurrObj->NumFatherObj);
// 	   TabObj[Num_inTotal].SubObj = True;
// 	   TabObj[Num_inTotal].FatherScale = FatherObj->ScaleMax+1;
// 	   TabObj[Num_inTotal].FatherNumberInScale  = FatherObj->ObjNumberInScale;
// 	   Nb_TotalSubObj ++;
// 	}
// 
// 
// // cout << "Xmax =  " <<  Gauss.mx+dx << " Ymax = " <<  Gauss.my+dy << endl;
// 
// 	// Flux Error calculation
// 	// search the 90% flux box (Sizebox), the signal standard deviation
//         // in this box and the noise standard deviation in this box
// 	// in the Poisson case, 
// 	// SNR = 90% Flux Source / sqrt(flux data in the box which contains  
// 	//                              the 90% of the flux) 
// 	// (and no SNR = (Nev_box - Nev_bgr) / sqrt(sigma^2_source +sigma^2_noise)
// 	// in the Gaussian case, SNR = sigma_signal / sigma_noise
// 	int Sizebox;
// 	float SigBox, SigNoise, FluxBox;
// 
//         // OUTPUT: Sizebox= size of the box containing 90% of the flux
//         //         SigBox=standard deviation in this box
//         //         SigNoise=standard deviation of the noise in this box
//         //                  calculated only if ModelData.UseRmsMap==True
//         // FluxBox = flux inside the box
//         box_obj_sig_and_noise(Im_rec, ModelData, Gauss,  Nlo, Nco,  Nl,  Nc,
//                         dx,  dy,  Flux, Sizebox,  SigBox, SigNoise, FluxBox);
// 
// //  cout << "Xmax =  " <<  Gauss.mx+dx << " Ymax = " <<  Gauss.my+dy << endl;
// 
//        switch (ModelData.which_noise())
//        {
//           case NOISE_GAUSSIAN: 
//           case NOISE_POISSON:
//           case NOISE_GAUSS_POISSON:
//           case NOISE_MULTI: 
//               // the error is proportional to the number of pixels
// 	      // of the box which contains 90% of the flux
// 	      // We consider that 90% of the flux calculation implies
// 	      // 90% of the error    
//               TabObj[Num_inTotal].ErrorFlux = 
//                     sqrt((double) Sizebox*Sizebox) * ModelData.SigmaNoise / 0.9;
//               
//               TabObj[Num_inTotal].ErrorFlux =
//                sqrt(TabObj[Num_inTotal].ErrorFlux*TabObj[Num_inTotal].ErrorFlux 
//                      + Eps_ErrorRec*Eps_ErrorRec);
//               TabObj[Num_inTotal].SNR_Obj = SigBox / ModelData.SigmaNoise;
//               break;
//           case NOISE_NON_UNI_ADD: 
// 	   if (ModelData.UseRmsMap == True) 
// 	   {
//  	      if (SigNoise > FLOAT_EPSILON)
//                   TabObj[Num_inTotal].SNR_Obj = SigBox / SigNoise;
// 	      else  TabObj[Num_inTotal].SNR_Obj = 10000.;
// 	      
// 	      // previous SigNoise was estimated only in the box 
// 	      // which contain 90% of the flux. We need here the 
// 	      // error in the box which contains the object
//  	      TabObj[Num_inTotal].ErrorFlux = Sizebox*Sizebox*SigNoise  / 0.9;
// 	   }
// 	   else 
// 	   {
// 	      TabObj[Num_inTotal].ErrorFlux = 0.;
// 	      TabObj[Num_inTotal].SNR_Obj = 0.;
// 	   }
// 	   break;
//         case NOISE_EVENT_POISSON:
//              { 
//               // the error is proportional to root square of the flux 
//               int ii,jj;
//               
// 	      SigNoise = 0.;
// 	      for (k = Sizebox/2-1; k <= Sizebox/2+1; k++)
// 	      for (l = Sizebox/2-1; l <= Sizebox/2+1; l++)
// 	      {   
// 	         ii = (int) (TabObj[Num_inTotal].PosY + k);
// 	         jj = (int) (TabObj[Num_inTotal].PosX + l);
// 	         if ( ind_ok(ii,Nl) && ind_ok(jj,Nc) )
// 	         {
// 	            SigNoise += Data(ii, jj);
// 	         }
// 	      }
// 	      SigNoise = sqrt(SigNoise);
// 	      TabObj[Num_inTotal].ErrorFlux = SigNoise / 0.9;
// 	      if (SigNoise > 0) 
//   	         TabObj[Num_inTotal].SNR_Obj = FluxBox / SigNoise;
//   	      else TabObj[Num_inTotal].SNR_Obj = 10000;
//  	    }
//  	   break;
//         case NOISE_UNDEFINED:
//         case NOISE_UNI_UNDEFINED:
//         case NOISE_NON_UNI_MULT:     
//           default:
// 	      TabObj[Num_inTotal].ErrorFlux = 0.;
// 	      TabObj[Num_inTotal].SNR_Obj = 0.;
//               break;
//        }
// 
// 
// 	// if a transform has been applied to an image, the inverse
// 	// transform must be performed on the reconstructed object,
// 	// This inverse transform can not be done in the case of
// 	// of Poisson noise with few events.
// 	// Several parameters must be reestimated:
// 	//   the flux, the error flux,  the moment, ...
// 	// SNR estimations don't have to be recalculated.
// 	if ((ModelData.TransImag == True) &&
//             (ModelData.which_noise() != NOISE_EVENT_POISSON))
// 	{
// 	   float ValIma;
// 	   float ValImaTrans;
// 
// 	   double FluxInvTrans=0.;
// 	   for (k=0;k<Nlo;k++)
// 	   for (l=0;l<Nco;l++)
// 	   {
// 	      ValIma = Data(k+dy,l+dx,I_CONT);
// 	      ValImaTrans = ModelData.val_transform(ValIma);
// 	      Im_rec(k,l) = ValIma - 
//                      ModelData.val_invtransform( ValImaTrans - Im_rec(k,l));
// 	      FluxInvTrans += Im_rec(k,l);
// 	   }
// 	   Estime_param_gauss_2D(Gauss, Im_rec);
// 	   TabObj[Num_inTotal].SigmaX = Gauss.sX;
// 	   TabObj[Num_inTotal].SigmaY = Gauss.sY;
// 	   TabObj[Num_inTotal].PosX = Gauss.mx + dx;
// 	   TabObj[Num_inTotal].PosY = Gauss.my + dy;
// 	   TabObj[Num_inTotal].Angle = Gauss.teta*RAD_DEG;
// 	   TabObj[Num_inTotal].ErrorFlux = FluxInvTrans * 
//                  TabObj[Num_inTotal].ErrorFlux / TabObj[Num_inTotal].Flux;
// 	   TabObj[Num_inTotal].Flux = FluxInvTrans;
// 	   TabObj[Num_inTotal].Magnitude = 
//                    (FluxInvTrans <= 0) ? 0. :-2.5 * (float) log10((double) FluxInvTrans);
// 	   TabObj[Num_inTotal].ValPixMax = Gauss.amp;
// 	}
// 
// if (Verbose == True)
// {
// cout<< "Object number: " << Num_inTotal+1 << " object:" << i+1 <<","<< j+1 << endl;
// cout << "   Last Scale = " <<  CurrObj->RootScale+1 << endl;
// cout<< "    Size = "<<"  "<< Nlo << "x" << Nco << endl;
// cout<< "    Depx = " << dx << " Depy = " << dy << endl;
// cout<< "    PosX = " << Gauss.mx + dx << " PosY = " << Gauss.my + dy << endl;
// cout<< "    SigmaX = " << Gauss.sX << " SigmaY = " << Gauss.sY << endl;
// cout<< "    Angle = " << Gauss.teta*RAD_DEG << endl;
//  if (AperPhot == True) 
// cout<< "    Aper. Phot. Flux = " << TabObj[Num_inTotal].Flux << endl;	     
// cout<< "    Reconstructed object Flux = " << TabObj[Num_inTotal].FluxObjRec << endl;
// cout<< "    ErrorFlux = " << TabObj[Num_inTotal].ErrorFlux << endl;
// cout<< "    Reconstruction Error " <<  TabObj[Num_inTotal].ErrorRec << endl;
// cout<< "    SNR Max Coef= " <<	TabObj[Num_inTotal].SNR_ValMaxCoef << endl; 
// cout<< "    Pos Max CoefX= " << TabObj[Num_inTotal].PosMaxCoef_X ;  
// cout<< "    Pos Max CoefY= " <<  TabObj[Num_inTotal].PosMaxCoef_Y << endl;
// cout<< "    SNR Obj= " <<	TabObj[Num_inTotal].SNR_Obj << endl; 
// cout<< "    SizeBox= " << Sizebox << endl;
// 
// if (TabObj[Num_inTotal].SubObj == True) 
//   {
//     cout<< "    Sub-object = True" << endl;
//     cout<< "    Father Number: " << TabObj[Num_inTotal].FatherScale << "-" <<  TabObj[Num_inTotal].FatherNumberInScale << endl;
//   }
// else cout<< "    Sub-object = False" << endl;
// } 
// 
// 	if (WriteAllObj == True)
// 	{
// 		sprintf(NameObj, "ima_obj_%d_%d.fits", i+1, j+1);
// 		io_write_ima_float(NameObj, Im_rec);
// 	}
// 
// 	// Additionne a l'image reconstruite finale si non sous-objet
// 	if (TabObj[Num_inTotal].SubObj == False) 
//                                add_image (Ima_SumObj, Im_rec, dx, dy);
// 
// 	if (KeepIma == True)
// 	{
// 	    TabObj[Num_inTotal].Image->alloc(Nlo, Nco, "TabObj");
// 	    *(TabObj[Num_inTotal].Image) = Im_rec;
// 	}
// 	Num_inTotal++;
//      }
//      else Nb_TotalSubObj ++;
//     }
//   }
//   Nb_TotalObj=Num_inTotal;
// }
// 
// /****************************************************************************/
// 
// void ListObj2D::create_list(Ifloat &Ima, Ifloat &Im_Segment, int Nb_Total, Bool WriteAllObj)
// {
//     int i,j,ind,Cpt=0;
// 
//    // initialization
//    Nl = Ima.nl();
//    Nc = Ima.nc();
//    Nb_TotalObj = Nb_Total;
// 
//    TabObj = new Object_2D[Nb_TotalObj];
// 
//    if (Verbose == True)
//    {
//       cout << "Total number of objects: " << Nb_TotalObj << endl;
//    }
// 
//    for (i = 0; i < Nl; i++)
//    for (j = 0; j < Nc; j++)
//       if (Im_Segment(i,j) != 0)
//       {
//          float Coef = Ima(i,j);
// 
//          Cpt ++;
//          ind = (int) (Im_Segment(i,j) + 0.5) - 1;
// 
//          // first passage
//          if (TabObj[ind].Surface == 0)
//          {
//             TabObj[ind].DepX = j;
//             TabObj[ind].DepY = i;
//             TabObj[ind].Nlo = 1;
//             TabObj[ind].Nco = 1;
//          }
//          else
//          {
//             if (j < TabObj[ind].DepX) TabObj[ind].DepX = j;
//             if (i < TabObj[ind].DepY) TabObj[ind].DepY = i;
//             if (j > TabObj[ind].Nco) TabObj[ind].Nco = j;
//             if (i > TabObj[ind].Nlo) TabObj[ind].Nlo = i;
//          }
//          TabObj[ind].Surface ++;
//          TabObj[ind].PosX  += Coef*j;
//          TabObj[ind].PosY  += Coef*i;
//          TabObj[ind].Mean += Coef;
//          TabObj[ind].Moment2XX += Coef*j*j;
//          TabObj[ind].Moment2YY  += Coef*i*i;
//          TabObj[ind].Moment2XY +=  Coef*i*j;
// 
//          if (ABS(Ima(i,j)) > TabObj[ind].ValPixMax)
//          {
//              TabObj[ind].ValPixMax = ABS(Ima(i,j));
//              TabObj[ind].Xmax = j;
//              TabObj[ind].Ymax = i;
//          }
//          if (i==0) TabObj[ind].Perimeter++;
//          if (i==Nl-1) TabObj[ind].Perimeter++;                 
//          if (j==0) TabObj[ind].Perimeter++;
//          if (j==Nc-1) TabObj[ind].Perimeter++;                 
// 
//          if ( (i>0) && (j>0) && (i<Nl-1) && (j<Nc-1)) 
//          {
//                 if (Im_Segment(i-1,j) < 0.5) TabObj[ind].Perimeter++;
//                 if (Im_Segment(i+1,j) < 0.5) TabObj[ind].Perimeter++;
//                 if (Im_Segment(i,j-1) < 0.5) TabObj[ind].Perimeter++;
//                 if (Im_Segment(i,j+1) < 0.5) TabObj[ind].Perimeter++;
//           }
//        }
// 
//    Signif =  (float) Cpt / (float)(Nl*Nc) * 100.;
// 
//    MeanMorpho = 0.;
//    for (i = 0; i < Nb_TotalObj; i++)
//    {
//       float Val, DiffVar,AddVar;
//       Val = (float) (TabObj[i].Perimeter*TabObj[i].Perimeter);
//       if (Val > FLOAT_EPSILON) 
//              TabObj[i].Morpho = 4. * PI * TabObj[i].Surface / Val;
//       else TabObj[i].Morpho = 0.;
//       MeanMorpho += TabObj[i].Morpho;
// 
//       TabObj[i].PosX /= TabObj[i].Mean;
//       TabObj[i].PosY /= TabObj[i].Mean;
//       TabObj[i].Moment2XX /= TabObj[i].Mean;
//       TabObj[i].Moment2YY /= TabObj[i].Mean;
//       TabObj[i].Moment2XY /= TabObj[i].Mean;
//       TabObj[i].Moment2XX -= TabObj[i].PosX*TabObj[i].PosX;
//       TabObj[i].Moment2YY -= TabObj[i].PosY*TabObj[i].PosY;
//       TabObj[i].Moment2XY -= TabObj[i].PosX*TabObj[i].PosY;
// 
//       DiffVar = TabObj[i].Moment2XX - TabObj[i].Moment2YY;
//       AddVar  = TabObj[i].Moment2XX + TabObj[i].Moment2YY;
//       if (ABS(DiffVar)<1E-7 || ABS(TabObj[i].Moment2XY)<1E-9)
//       {
//          TabObj[i].SigmaX = sqrt(ABS(TabObj[i].Moment2XX));
//          TabObj[i].SigmaY = sqrt(ABS(TabObj[i].Moment2YY));
//          TabObj[i].Angle = 0.;
//        }
//       else
//       {
//          TabObj[i].Angle = 0.5 * atan(2*TabObj[i].Moment2XY / DiffVar);
//          TabObj[i].SigmaX = sqrt( ABS(
//                    AddVar/2. + TabObj[i].Moment2XY/sin(2.*TabObj[i].Angle) ));
//          TabObj[i].SigmaY = sqrt( ABS(
//                     AddVar/2. - TabObj[i].Moment2XY/sin(2.*TabObj[i].Angle) ));
// 	TabObj[i].Angle *= RAD_DEG;
//        }
//       TabObj[i].Flux = TabObj[i].Mean;
//       if (TabObj[i].Surface) TabObj[i].Mean /= (float) TabObj[i].Surface;
//       else  TabObj[i].Mean = 0.;
//       TabObj[i].NumObj = i+1;
//       TabObj[i].Nli = Nl;
//       TabObj[i].Nci = Nc;
//       TabObj[i].Magnitude =  (TabObj[i].Flux <= 0) ? 0. :
//                                  -2.5 * (float) log10((double) TabObj[i].Flux);
//       TabObj[i].SubObj = False;
//       TabObj[i].Nlo = TabObj[i].Nlo - TabObj[i].DepY + 1;
//       TabObj[i].Nco = TabObj[i].Nco - TabObj[i].DepX + 1;
//    }
//    if (Nb_TotalObj > 0) MeanMorpho /= (float) Nb_TotalObj;
// 
//    // Write the results on the screen
//    if (Verbose == True)
//    {
//       for (i = 0; i < Nb_TotalObj; i++)
//       {
//       cout << "Object number: " << i << endl;
//       cout << "    Size = "<<"  "<< TabObj[i].Nlo << "x" << TabObj[i].Nco << endl;
//       cout << "    Depx = " << TabObj[i].DepX << " Depy = " << TabObj[i].DepY << endl;
//       cout << "    PosX = " << TabObj[i].PosX << " PosY = " << TabObj[i].PosY <<       endl;
//       cout<< "    SigmaX = " << TabObj[i].SigmaX << " SigmaY = " << TabObj[i].SigmaY << endl;
//       cout<< "    Angle = " << TabObj[i].Angle << endl;
//       cout<< "    Flux = " << TabObj[i].Flux << endl;
//       }
//    }
// 
//    // Creates the small images
//    if ((KeepIma == True) || (WriteAllObj == True))
//     for (i = 0; i < Nb_TotalObj; i++)
//     {
//        Ifloat *ImaObj = new Ifloat(TabObj[i].Nlo, TabObj[i].Nco, "TabObj");;
//        TabObj[i].get_box_obj(Ima, *ImaObj);
// 
//        if (KeepIma == True) TabObj[i].Image = ImaObj;
//        if (WriteAllObj == True)
//        {
//           char NameObj[80];
//           sprintf(NameObj, "ima_obj_%d.fits", i);
//           io_write_ima_float(NameObj, *ImaObj);
//        }
//        if (KeepIma == False) delete ImaObj;
//     }
// }

/****************************************************************************/

void Object_2D::get_box_obj(Ifloat &Ima_Big, Ifloat & ImaObj)
{
   int i,j;
   ImaObj.resize(Nlo,Nco);
   for(i=0;i<Nlo;i++)
   for(j=0;j<Nco;j++)
      ImaObj(i,j) = Ima_Big(DepY+i,DepX+j);
}

/****************************************************************************/

void Object_2D::draw (Ifloat &Ima_Big, Bool DotFlag, float Val)
{
    float Zoomx = (float) Ima_Big.nc() / (float) Nci;
    float Zoomy = (float) Ima_Big.nl() / (float) Nli;
    float Zoom = MIN(Zoomx, Zoomy);
    double a = SigmaX*3.*Zoom;
    double b = SigmaY*3.*Zoom;;
    double Theta0 = Angle;

    if (Val == 0.) Val = max(Ima_Big) * 1.5;
    im_draw_ellips(Ima_Big, PosX*Zoom, PosY*Zoom, a, b, Theta0, Val, DotFlag);
}

/****************************************************************************/

void Object_2D::simu (Ifloat &Ima_Big)
{
    float Zoomx = (float) Ima_Big.nc() / (float) Nci;
    float Zoomy = (float) Ima_Big.nl() / (float) Nli;
    float Zoom = MIN(Zoomx, Zoomy);
    double sx = SigmaX*Zoom;
    double sy = SigmaY*Zoom;
    double sX2 = sx*sx;
    double sY2 = sy*sy;
    double X,Y,x,y;
    double teta= Angle*PI/180.;
    int Nl=Ima_Big.nl();
    int Nc=Ima_Big.nc();
    int i,j, indi, indj;
    double dx = (PosX - DepX)*Zoom;
    double dy = (PosY - DepY)*Zoom;
    int Nl_obj = (int) (Nlo * Zoom);
    int Nc_obj = (int) (Nco * Zoom);
    int Z_DepX = (int) (DepX * Zoom);
    int Z_DepY = (int) (DepY * Zoom);

    for(i=0;i<Nl_obj;i++)
    for(j=0;j<Nc_obj;j++)
    {
        indi = i + Z_DepY;
        indj = j + Z_DepX;
        if ((indi >= 0) && (indj >= 0) && (indi < Nl) && (indj < Nc))
        {
           x = (float) j - dx;
           y = (float) i - dy;
           X =  x*cos(teta) + y*sin(teta);
           Y = -x*sin(teta) + y*cos(teta);
           if ((sX2 > FLOAT_EPSILON) && (sY2 > FLOAT_EPSILON)) 
                  Ima_Big(indi,indj) += ValPixMax*exp(-(X*X/sX2+Y*Y/sY2)/2.);
        }
    }
}

/****************************************************************************/


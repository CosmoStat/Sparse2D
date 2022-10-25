/******************************************************************************
**                   Copyright (C) 2008 by CEA  
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author:  Jean-Luc Starck  
**
**    Date:  25/10/08
**    
**    File:  mrs_inpainting.cc 
**
*******************************************************************************
**
**    DESCRIPTION:  Inpainting + Iterative wiener filtering
**    -------------------
**
******************************************************************************/

#include"HealpixClass.h"

char Name_Imag_In[512]; /* input file image */
char Name_Imag_Out[512]; /* output file name */
char Name_Mask[512]; 
char Name_RMSMap [512]; 
char SignalPowSpec_FileName [512];
char NoisePowSpec_FileName [512];
char Name_PowSpec_Out [512];

extern int  OptInd;
extern char *OptArg;
extern int  GetOpt(int argc, char **argv, char *opts);

Bool Verbose=False;
float SigmaNoise=0.;
#define  DEF_IWIENER_NBR_ITER  40
#define  DEF_ITER_GAUSSIANITY 40
 
 int Max_Iter=DEF_IWIENER_NBR_ITER;
Bool UseRMSMap = False;   // Apply a denoising using an iterative Wiener filtering using a RMS map.
Bool UseMask = True;          // Apply an inpainting
Bool ThresholdResi = True;   // Threshold the Alm of the residual instead of the Alm of the solution
Bool WienerOK = False;        // Apply a denoising using an iterative Wiener filtering
Bool HardThreshold = True;  // Apply an Alm iterative hard threshold for the inpainting
Bool SoftThreshold = False;
Bool Use_ZeroMeanConstraint = True; // Apply a zero mean constraint on the solution
Bool Use_WP_Constraint = True;   // Apply a variance constraint in the wavelet packet of the solution
Bool EstimPowSpec = False;
Bool All_WP_Band = False;

float  ZeroPadding=0.1;
int Lmax = 0.;
Bool OptN_PS = False;
Bool OptS_PS = False;

float Eps=1.;
Bool GaussianityConstraint = False;

inline double meyer_lowpass_fct(double x)
{
     double l,r;
     l = exp( 1-1 / (1-exp(1-1/(1-x))));
     r = exp(1-1/(1-exp(1-1/x)));
     return (l /= sqrt(l*l+r*r));
}

/***************************************************************************/
 
void make_meyer_filter_lowpass(fltarray &H,  int PixStart2Zero, int PixEnd2Zero)
{
	int i;
	if (PixEnd2Zero <= PixStart2Zero) PixEnd2Zero = PixStart2Zero + 1;
	int Npix2zero = PixEnd2Zero - PixStart2Zero + 1;
    H.init();
    for (i=0; i < PixStart2Zero; i++)  H(i) = 1.;    
 	for (i=0; i < MIN(H.nx(),PixEnd2Zero); i++)
	{
		 double x = i;
		 double r =  x  / ( (double) Npix2zero -1);
		 if ((PixStart2Zero+i) < H.nx())
		 {
		    if (r <=0) H(PixStart2Zero+i) = 1.;
 		    else if (r >= 1)  H(PixStart2Zero+i) = 0.;
		    else H(PixStart2Zero+i) = meyer_lowpass_fct(1.-r);
		 }
		 // cout <<  i+1 << " nx = " << H.nx() << " r= " << r  <<  "  H = " << H(i) << endl;
 	}
 	for (i=PixEnd2Zero; i < H.nx(); i++)  H(i) = 0;
}

/***************************************************************************/
#define PS_NB 8
int  PS_From[PS_NB] = {0, 2,  11,  31, 151, 421,  1201, 2501};
int  PS_To[PS_NB] = {1, 10, 30, 150, 420, 1200, 2500, 3000};
int  PS_Step[PS_NB] = {2, 1,   2,   5,  10,  20,  50,  100};

void get_planck_wp_meyer_filter(fltarray &TabH, fltarray &TabG, fltarray &Win, int Lmax, int LmaxT)
{
    int Nw=0;
    for (int s=0; s < PS_NB; s++)
    for (int f=PS_From[s] ; f < PS_To[s] ; f+=PS_Step[s])
    {
    	if (f <= Lmax) Nw++;
    }
  //  cout << "Nbr Win = " << Nw << endl;
    TabH.resize(LmaxT+1, Nw);
    Win.resize(Nw,2);
    fltarray H(TabH.nx());
    
    int w=0;
    for (int s=0; s < PS_NB; s++)
    for (int f=PS_From[s] ; f < PS_To[s] ; f+=PS_Step[s])
    {
    	int Beg = MAX(0,f - PS_Step[s]/2);
        int End = f + PS_Step[s]/2;
         if (f <= Lmax)
         {
          //  cout << w+1 << " " << Beg << "--> " << End << endl;
              make_meyer_filter_lowpass(H, f,  End);
              for (int l=0; l < H.nx(); l++)  TabH(l,w) = H(l);
              Win(w,0) = Beg;
              Win(w,1) = End;
              w++;
         }
    }
    
   // fits_write_fltarr("xxh.fits", TabH);
  // exit(0);
}

int TNbrWin[9] = {22,22,11,8,17,11,5,3,4};
double TabWin2048[22] ={64,128,192,256,320,384,512,768,1024,1280,1536,1664,1792,1920,2048,2176,2304,2432,2560,2688,2816,3000}; 
double TabWin1024[22] ={64,128,192,256,320,384,512,768,1024,1280,1536,1664,1792,1920,2048,2176,2304,2432,2560,2688,2816,3000}; 
double TabWin512[11] = {64,128,192,256,320,384,512,768,1024,1280,1536};
double TabWin256[8] = {64,128,192,256,320,384,512,768};
double TabWin128[17] = {64,80,96, 112,128,144,160, 176, 192,256, 272, 288,304, 320,336,352, 384};
double TabWin64[11] = {16,32,64,80,96, 112,128,144,160, 176, 192};
double TabWin32[5] = {8,16,32,64,92};
double TabWin16[3] = {8,16,48};
double TabWin8[4] = {4,8,16,24};

// Compute the filters for the Wavelet Packet decompositions.
void get_wp_meyer_filter(int nside, fltarray &TabH, fltarray &TabG, fltarray &Win, int  Lmax)
{
	int i,p,NbrBand=0;
	double *PtrTab=NULL;
	fltarray TabCl;
   // cout << "get_wp_meyer_filter: " <<nside << endl;
	
	if (nside == 2048)
	{
		NbrBand= TNbrWin[0];
		PtrTab = TabWin2048;
	}
	else if (nside == 1024)
	{
		NbrBand= TNbrWin[1];
		PtrTab = TabWin1024;
	}
	else if (nside == 512)
	{
		NbrBand= TNbrWin[2];
		PtrTab = TabWin512;
	}
	else if (nside == 256)
	{
		NbrBand= TNbrWin[3];
		PtrTab = TabWin256;
	}
	else if (nside == 128)
	{
		NbrBand= TNbrWin[4];
		PtrTab = TabWin128;
	}
	else if (nside == 64)
	{
		NbrBand= TNbrWin[5];
		PtrTab = TabWin64;
	}
	else if (nside == 32)
	{
		NbrBand= TNbrWin[6];
		PtrTab = TabWin32;
	}
	else if (nside == 16)
	{
		NbrBand= TNbrWin[7];
		PtrTab = TabWin16;
	}
	else if (nside == 8)
	{
		NbrBand= TNbrWin[8];
		PtrTab = TabWin8;
	}
	intarray TabB;
	
	int l=0;
	while ((l < NbrBand) && (PtrTab[l] <= Lmax)) l++;
	if (l < NbrBand) NbrBand = l;
	
	
 	if (PtrTab[NbrBand-1] == Lmax) TabB.alloc(NbrBand);
 	else TabB.alloc(NbrBand+1);
 	
	for (int i=0; i < NbrBand; i++)
	{
		TabB(i) = PtrTab[i];
        // cout << PtrTab[i] << endl;
	}
	if (PtrTab[NbrBand-1] != Lmax) 
	{
		TabB(NbrBand) = Lmax;
		NbrBand++;
	}
  	int LM = Lmax+1;
  	
 	fltarray H,Res;
 	Res.alloc(LM, NbrBand+1);
	Win.alloc(NbrBand,2);
    // cout << "get_wp_meyer_filter: " <<NbrBand << endl;
	
	H.resize(Lmax);
	int PixStart2Zero, PixEnd2Zero;
	int Np = TabB(0);
	float PercentPix2Zero=0.1;
	PixEnd2Zero = TabB(0);
	PixStart2Zero = TabB(0) - PercentPix2Zero * Np;
 	 make_meyer_filter_lowpass(H,  PixStart2Zero, PixEnd2Zero);
    // cout << "make_meyer_filter_lowpass: " <<  PixStart2Zero << " " << PixEnd2Zero<< endl;
 	// make_meyer_filters(TabB(0), Scale, H);
  	//  fits_write_fltarr((char *)"xx_0.fits", H);

	for (p=0; p < H.nx(); p++) Res(p,0) = H(p);
	Win(0,0) = 0;
	Win(0,1) = TabB(0);
	Win(1,0) = TabB(0);
	Win(1,1) = TabB(1);
	 //   fits_write_fltarr("xxres.fits", Res);
	 //   fits_write_fltarr("xxWin.fits", Win);

    // cout << " NbrBand = " << NbrBand << " " << H.nx() << endl;
	for (i=1; i < NbrBand; i++)
	{
	 // char FN[200];
	 //  cout << "   N = " << TabB(i-1) << " -> " << TabB(i) << endl;
 	  PixEnd2Zero = TabB(i);
      PixStart2Zero = TabB(i) - PercentPix2Zero * Np;
 	   make_meyer_filter_lowpass(H,  PixStart2Zero, PixEnd2Zero);

 	   // sprintf(FN, "xx_%d.fits", i);
 	   // fits_write_fltarr(FN, H);

 	   for (p=0; p < H.nx(); p++) Res(p,i) = H(p);
	}
 	for (p=0; p < Lmax; p++) Res(p, NbrBand) = 1.;
  	
 	fltarray Final;
   // cout << " Final = " << Lmax << " " << NbrBand << endl;
 	Final.alloc(LM, NbrBand+1);
 	for (p=0; p < LM; p++) Final(p, 0) = Res(p, 0);
 	for (i=0; i < NbrBand; i++)
	for (p=0; p < LM; p++) 
  	   Final(p, i+1) = Res(p, i+1) - Res(p, i);
 	
    // fltarray TabG;
    TabH.alloc(LM, NbrBand);
    TabG.alloc(LM, NbrBand);
    for (i=0; i < NbrBand; i++)
	for (p=0; p < LM; p++) 
	{
		TabH(p,i) = Res(p, i);
		if (i != NbrBand-1) TabG(p,i) = Res(p, i+1) - Res(p, i);
		else TabG(p,i) = 1. - Res(p, i);
	}
    // cout << " TabH = " << Lmax << " " << NbrBand << endl;
	
	for (i=1; i < NbrBand; i++)
	{
		Win(i,0) = TabB(i-1);
		Win(i,1) = TabB(i);	
	}
    // cout << " Win = " << Lmax << " " << NbrBand << endl;
     // fits_write_fltarr("xxh1.fits", TabH);
    // fits_write_fltarr("xxg1.fits", TabG);
}
 

/***************************************************************************/
 
static void usage(char *argv[])
{
    fprintf(OUTMAN, "Usage: %s options in_map in_mask out_map  [out_powspec]  \n\n", argv[0]);
    fprintf(OUTMAN, "   where options =  \n");

    fprintf(OUTMAN, "         [-i Nbr_of_Iter]\n");
    fprintf(OUTMAN, "             Default is %d.\n", DEF_IWIENER_NBR_ITER);
        
    fprintf(OUTMAN, "         [-l lmax]\n");
    fprintf(OUTMAN, "             Default is MIN(3000, 3*nside) \n");
    
    fprintf(OUTMAN, "         [-p Harmonics_ZeroPadding (in [0,1] )]\n");
    fprintf(OUTMAN, "             Default is 0.1.\n");

    fprintf(OUTMAN, "         [-e Relaxation_parameter (in ]0,1] )]\n");
    fprintf(OUTMAN, "             Default is 1.\n");
    
    fprintf(OUTMAN, "         [-G]\n");
    fprintf(OUTMAN, "             Add a Gaussianity constraint.\n");
    
    fprintf(OUTMAN, "         [-H]\n");
    fprintf(OUTMAN, "             Replace the hard-thresholding by a soft-thresholding.\n");
    exit(-1);
}
 
/*********************************************************************/

/* GET COMMAND LINE ARGUMENTS */
static void sinit(int argc, char *argv[])
{
    int c; 
    Bool OptG = False;

    /* get options */
    while ((c = GetOpt(argc,argv,(char *) "WHGae:l:p:g:i:vzZ")) != -1) 
    {
	switch (c) 
        { 
         	case 'W': Use_WP_Constraint = (Use_WP_Constraint == True) ? False: True;  
        	           break;
        	case 'H': SoftThreshold = (SoftThreshold == True) ? False: True;
        	          HardThreshold = (SoftThreshold == True) ? False: True;
        	          break;
         	case 'G': GaussianityConstraint = (GaussianityConstraint == True) ? False: True;  
        	           break;

        	case 'a': All_WP_Band = (All_WP_Band == True) ? False: True;  
        	           break;
           case 'e':
 	     	      if (sscanf(OptArg,"%f",&Eps) != 1) 
                  {
		            fprintf(OUTMAN, "Error: relaxation parameter  value: %s\n", OptArg);
		            exit(-1);
	            	}
 		          break;
           case 'p':
 	     	      if (sscanf(OptArg,"%f",&ZeroPadding) != 1) 
                  {
		            fprintf(OUTMAN, "Error: bad zero padding  value: %s\n", OptArg);
		            exit(-1);
	            	}
 		          break;
	        case 'l':
 	     	      if (sscanf(OptArg,"%d",&Lmax) != 1) 
                  {
		            fprintf(OUTMAN, "Error: bad lmax  value: %s\n", OptArg);
		            exit(-1);
	            	}
 		          break;
            case 'i':  
                   if (sscanf(OptArg,"%d",&Max_Iter) != 1) 
                   {
		            fprintf(OUTMAN, "Error: bad Max_Iter: %s\n", OptArg);
		    		exit(-1);
		           } 
                   if (Max_Iter <= 0)   Max_Iter =  DEF_IWIENER_NBR_ITER;
 	        	break;
	        case 'g':
 	     	      if (sscanf(OptArg,"%f",&SigmaNoise) != 1) 
                  {
		            fprintf(OUTMAN, "Error: bad sigma noise: %s\n", OptArg);
		            exit(-1);
	            	}
	            	OptG = True;
	            	WienerOK=True;
 		          break;
   	    case 'v': Verbose = True; break;
            case '?': usage(argv); break;
	    default: usage(argv); break;
 		}
	} 
      
        /* get optional input file names from trailing 
          parameters and open files */
    if (OptInd < argc) strcpy(Name_Imag_In, argv[OptInd++]);
         else usage(argv);

    if (OptInd < argc) strcpy(Name_Mask, argv[OptInd++]);
         else usage(argv);
           
	if (OptInd < argc) strcpy(Name_Imag_Out, argv[OptInd++]);
         else usage(argv);

	if (OptInd < argc) 
	{
		   strcpy(Name_PowSpec_Out, argv[OptInd++]);
		   EstimPowSpec = True;
	}

	/* make sure there are not too many parameters */
	if (OptInd < argc)
        {
		fprintf(OUTMAN, "Too many parameters: %s ...\n", argv[OptInd]);
		exit(-1);
	}
}
  
/*********************************************************************/

void alm_gaussianity_constraint(CAlmR &ALM, PowSpec & P, Bool UpdatePowSpec)
{
	// cout << "IN" << endl;
	int NiterG = 20;
    double Alpha = 0.05;
    double EpsM = Alpha;
    double EpsP = Alpha;
    double Ro = 1.;
    long NbrL = ALM.Lmax()+1;
    long NbrM = ALM.Lmax()+1;
    
    dblarray LM,LP;
    LM.alloc(NbrL,NbrM);
    LP.alloc(NbrL,NbrM);
    Bool Median=False;
    dblarray  StdPS;
    
    // xcomplex<REAL> A;
    fltarray A;
    A.alloc(NbrL, NbrM, 2);
    
/*    dblarray SpecM;
    SpecM.alloc(NbrL);
    for (int l=1; l <=  ALM.Lmax(); l++)
    {
    	int NbrM = 2*l+1;
    	dblarray SpecM;
    	SpecM.resize(NbrM);
    	SpecM(l) = norm ( ALM(l,0) );
    	for (int m=1; m <= l; m++)
    	{
    	   SpecM(l-m) = norm ( ALM(l,m) );
    	   SpecM(l+m) =	SpecM(l-m);
    	}
    	StdPS(l) = SpecM.sigma() ;
    }*/
    
   if (UpdatePowSpec == True)
   {
        if (Median == False) ALM.alm2powspec(P);
        else ALM.extract_median_powspec(P);
   }
    StdPS.alloc(NbrL);
    for (int l=1; l <=  ALM.Lmax(); l++) StdPS(l) =  P.tt(l);
    //StdPS.info("P");

    StdPS(0) = sqrt(P.tt(0));
 
    for (int l=1; l <=  ALM.Lmax(); l++) StdPS(l) = sqrt( P.tt(l)) * sqrt(1.+1./(2.*l+1.));
    // StdPS.info("StdPS");
    
    for (int IterG =0; IterG < NiterG; IterG ++)
    {
     	 for (int l=0; l <=  ALM.Lmax(); l++)
    	 {
    	 	double Mu = sqrt(P.tt(l));
            for (int m=0; m <= l; m++)
            {  
                // double m2 = sqrt( norm(ALM(l,m)));
                double T = (LP(l,m) - LM(l,m)) / 2.;
                xcomplex<REAL> SoftCoef = ALM.lm_soft_threshold(  l,   m, (float) T);
                A(l,m,0) = SoftCoef.real();
     		    A(l,m,1) = SoftCoef.imag();
     		    double m2n = sqrt(A(l,m,0)*A(l,m,0) + A(l,m,1)*A(l,m,1));
     		     LP(l,m) = MAX( LP(l,m) + Ro*(m2n - Mu - EpsP * StdPS(l) ), 0);
     		     LM(l,m) = MAX( LM(l,m) - Ro*(m2n - Mu - EpsM * StdPS(l)), 0);  
    	    }
    	    
    	 }
    	 // cout << IterG << ": MinA = " << A.min() << " MaxA = " << A.max() << " sigA = " << A.sigma() << endl;
    }
     
   for (int l=1; l <=  ALM.Lmax(); l++)
   for (int m=0; m <= l; m++) 
   {
   	   ALM(l,m) = xcomplex<REAL>(A(l,m,0), A(l,m,1));
   }
}

/*********************************************************************/

int main(int argc, char *argv[])
{
    int l,k;
    double Min, Max;
    fitsstruct Header, HD1;
    char Cmd[512];
    Cmd[0] = '\0';
    for (k =0; k < argc; k++) sprintf(Cmd, "%s %s", Cmd, argv[k]);
     
    /* Get command line arguments, open input file(s) if necessary */
    sinit(argc, argv);
    
   if (GaussianityConstraint == True)
   {
     	ZeroPadding = 0.;
        if (Max_Iter==DEF_IWIENER_NBR_ITER) Max_Iter = DEF_ITER_GAUSSIANITY;
   }


    if (Verbose == True)
    { 
        cout << "# PARAMETERS: " << endl;
        cout << "# File Name in = " << Name_Imag_In << endl;
        cout << "# InpMap  File Name Out = " << Name_Imag_Out << endl;   
        if (EstimPowSpec ==True)  cout << "# File Name Out = " << Name_Imag_Out << endl;   
        if (UseRMSMap ==  True) cout << "# File Name Noise RMS = " <<   Name_RMSMap << endl;  
        else if (SigmaNoise > 0)  cout << "Sigma Noise = " <<   SigmaNoise << endl; 
        if ( UseMask ==  True) cout << "# File Name Mask = " <<   Name_Mask << endl;  
        cout << "# Number of iterations = " <<  Max_Iter << endl;   
        if (Lmax > 0)  cout << "# Lmax = " <<  Lmax << endl;
        if (ZeroPadding > 0)  cout << "# ZeroPadding = " <<  ZeroPadding << endl;
        if (OptN_PS > 0)  cout << "# Noise power spectrum file name = " <<  NoisePowSpec_FileName << endl;
        if (OptS_PS > 0)  cout << "# Signal power spectrum file name = " <<  NoisePowSpec_FileName << endl;
        if (GaussianityConstraint == True) cout << "# Use Gaussianity constraint " << endl;
        if (HardThreshold == True) cout << "# Use HardThresholding " << endl;
        if (SoftThreshold == True) cout << "# Use SoftThresholding " << endl;
        if (Use_WP_Constraint == True) cout << "# Use WP constraint " << endl;
        if (ThresholdResi == True) cout << "# Use residal thresholding " << endl;
#ifdef MRSDEBUG                      
        cout << "DEBUG MODE " << endl;
#endif
    }
    
    
   Hdmap Map;
   Hdmap MaskMap;

   Map.read(Name_Imag_In);
   if (Verbose == True)
   {
   	   Map.minmax(Min,Max);
   	   if (Verbose == True) Map.info((char*) "Input MAP");
   }
   
    MaskMap.read(Name_Mask);
   
   if (Verbose == True) MaskMap.info((char*) "MaskMap");
   int NbrTotMask = 0; 
   double SigData = 0.;
   double MeanData = 0.;
   double FSky = 1.;
   if (UseMask  ==  True)  
   {
   	   for (int p=0; p<  Map.Npix(); p++)  
   	   {
    	   if (MaskMap[p] > 0)  
   	   	   {
   	   	   	   NbrTotMask ++;
   	   	   	   MeanData += Map[p];
   	   	   }
   	   }
   	   MeanData /= (double) NbrTotMask;
   	   
     // for (int p=0; p<  Map.Npix(); p++) 
    //       if (MaskMap[p] > 0) Map[p] -= MeanData;
           
   	   for (int p=0; p<  Map.Npix(); p++)  
   	   {
           if (MaskMap[p] > 0)   SigData += (Map[p]-MeanData)*(Map[p]-MeanData);
       }
       SigData = sqrt(SigData / (double) NbrTotMask);
       FSky = (double) Map.Npix() / NbrTotMask;
       if (Verbose == True) cout << "Percentage of good pixels: " << 1. / FSky * 100. << endl;
   }
   
   // We start with a zero mean image
   for (int p=0; p<  Map.Npix(); p++) if (MaskMap[p] > 0) Map[p] -= MeanData;
             
   int Nside = Map.Nside();
   int Lmax_Tot = mrs_get_lmax (Lmax,  Nside,   ZeroPadding);
   if (Verbose == True) cout << "==> Use Lmax = " << Lmax << ",  Lmax TOT (with zero padding)  = " <<  Lmax_Tot << endl;
 
   CAlmR  ALM;  
   ALM.alloc(Nside, Lmax_Tot);
   ALM.Norm = True;
   ALM.set_beam_eff(Lmax, Lmax_Tot);
   ALM.UseBeamEff = False;

   fltarray WP_H, WP_G, WP_W;
   int NbrWP_Band=0;
   CAlmR  ALM_Band;
   Hdmap Band;
   double Sig_InMap, Sig_Mask;
  
    if (Use_WP_Constraint == True)
    {
        ALM_Band.alloc(Nside, Lmax_Tot);
        ALM_Band.Norm = ALM.Norm;
        ALM_Band.set_beam_eff(Lmax, Lmax_Tot);
        ALM_Band.UseBeamEff =  ALM.UseBeamEff ;
       
        if (All_WP_Band == False) get_wp_meyer_filter(Nside, WP_H, WP_G, WP_W, Lmax_Tot);
        else get_planck_wp_meyer_filter(WP_H, WP_G, WP_W, Lmax, Lmax_Tot);
        NbrWP_Band = WP_H.ny();
        Band.SetNside ((int) Nside,  (Healpix_Ordering_Scheme) DEF_MRS_ORDERING);
        
// #ifdef MRSDEBUG                      
      fits_write_fltarr((char *)"xxh.fits", WP_H);
//      fits_write_fltarr((char *)"xxg.fits", WP_G);
//      fits_write_fltarr((char *)"xxw.fits", WP_W);
// #endif
     }

         
     //Signal noise power spectrum estimation.
    if (Verbose == True)  Map.info((char *) "INIT");
    ALM.alm_trans(Map);
    // ALM.write((char *) "xx_init_alm.fits", True);
    
   double MaxAlm =ALM.max_absalm();
   if (Verbose == True)  cout << "==> Max Alm = " << MaxAlm << " N = " << ALM.NormVal << endl;
   
         
    Hdmap Result;
    Hdmap Resi;
   	Result.SetNside ((int) Nside,  (Healpix_Ordering_Scheme) DEF_MRS_ORDERING);
   	Result.fill(0.);
   	Resi.SetNside ((int) Nside,  (Healpix_Ordering_Scheme) DEF_MRS_ORDERING);
     
  SigmaNoise = 1.;
  double FirstSoftThreshold = MaxAlm / SigmaNoise * 0.95;
  if (SoftThreshold == True) FirstSoftThreshold *= sqrt(2.);
  double LastSoftThreshold = 0.;
  double  DeltaThreshold = FirstSoftThreshold - LastSoftThreshold;
  // double StepL = DeltaThreshold / (float) (Max_Iter -1);
  double Lambda = FirstSoftThreshold;
   
   double Err=Map.rms();
   int MaxNonZeroL=0; 
   for (int i=0; i < Max_Iter; i++)
   { 
   	      // cout << FirstSoftThreshold << " " << LastSoftThreshold  << endl;
       	  Lambda =  LastSoftThreshold  + DeltaThreshold  *(1.-erf(2.8 * i/ Max_Iter));  // 4 versus 2.8  in IDL code
   	      if ((Lambda < LastSoftThreshold) || (i == Max_Iter-1)) Lambda = LastSoftThreshold;
   	      
   	      if (Verbose == True)  cout << "Iter  " << i+1 << ", Lambda = " << Lambda  <<  " Err = " << Err << endl;
   	      
   	      for (int p=0; p < Result.Npix(); p++) 
   	      {
   	      	   // Compute the residual
   	           Resi[p] =  Map[p] - Result[p];
     	       if (UseMask  ==  True)  Resi[p] *= MaskMap[p];
    	       Err += Resi[p]*Resi[p];
   	           if (ThresholdResi == False) Resi[p] += Result[p];
   	       }
             
            int CptAlm=0;
		      
		   // Resi.info((char*) "     --> BefW MAP");
           ALM.alm_trans(Resi);
           // ALM.info((char *) "ALM bef");
           float ThresholdLevel = SigmaNoise*Lambda;
           if ((HardThreshold == True) &&  (Lambda > 0)) CptAlm = ALM.hard_threshold(ThresholdLevel, MaxNonZeroL);
           else if ((SoftThreshold == True) &&  (Lambda > 0)) CptAlm = ALM.soft_threshold(ThresholdLevel*sqrt(2.), MaxNonZeroL);
           if (Verbose == True) cout << "  Number of unsignificant Alm = " << CptAlm <<  " MaxL = " << MaxNonZeroL << endl; 
           // ALM.info((char *) "ALM After");
           
           ALM.alm_rec(Resi);
 
   	       char NN[256];
   	    //   sprintf(NN, "xx_alm_%d.fits", i+1);
   	    //   ALM.write((char *) NN, True);
   	     
   	       // ALM.alm_rec(Resi);
           // sprintf(NN, "xx_it_%d.fits", i+1);
  		   // Resi.write(NN);
 
            // Resi.info((char*) "     -> AftW MAP");
		      // ALM.WienerFilter.info("WF: ");
		      // fits_write_fltarr(fitsname("xx_wien"), ALM.WienerFilter);
      
            if (ThresholdResi == False) 
            {
            	   for (int p=0; p<  Result.Npix(); p++)   Result[p] = (1-Eps)*Result[p] + Eps*Resi[p];
            }
            else for (int p=0; p<  Result.Npix(); p++)   Result[p]  +=  Eps*Resi[p];
            

   	       if (GaussianityConstraint ==  True) 
   	       {
   	       	   Bool UpdatePowSpec=False;
   	       	   PowSpec  P;
   	       	   ALM.alm_trans(Result);
   	       	   ALM.alm2powspec(P);
   	       	   alm_gaussianity_constraint(ALM, P, UpdatePowSpec);
   	       	                            // If Gaussianity constraint, we want to keep the ALM coeff
   	                                    // since we have fixed other spatial constraints to zero
   	                                    // and force the WP constraint, which requires also Alm coefficient
   	       }              
   	       
            //sprintf(NN, "xx_res_it_%d.fits", i+1);
           //Result.write(NN);
  		  // Result.info("  RES");
  		   
            if ((Use_WP_Constraint == True) &&  (i >= Max_Iter/2))
            {
            	// We compute the WP transform:
            	
            	//   ITER: Result contains the high resolution   (c_j)
            	//   We compute Band, which contains the low resolution (c_{j+1}
            	//   Coef =  Result - Band
            	//   Result = Band
            	
                  if (GaussianityConstraint ==  False) ALM.alm_trans(Result);
                  Resi.init();
                  // WP_H.info("WPH");
                  for (int b=0; b < NbrWP_Band; b++)
                  {  
                  	  int LMin = WP_W(NbrWP_Band-1-b,0) ;
   	                  int LMax = WP_W(NbrWP_Band-1-b,1) ;
   	                  if (Verbose == True) cout << "        WP:Band " << b+1 << ", Lmin = " << LMin << ", Lmax = " << LMax << endl; 
   	                  if  (LMin < MaxNonZeroL)
   	                  {
   	                  	
                      ALM_Band.Set(LMax,  LMax);
   	                  ALM_Band.SetToZero();
   	                  
   	                  // Compute the solution at a given resolution (convolution with H)
     	              for (int l=0; l <= MIN(ALM_Band.Lmax(),WP_H.nx()); l++)
                      for (int m=0; m <= l; m++)     ALM_Band(l,m) =  ALM(l,m) * xcomplex<REAL>(WP_H(l, NbrWP_Band-1-b),0);
    	                	                  
                        // cout << "        WP Band " << b+1 << ",  MAX ABS(ALM) = " << ALM.max_absalm() << endl;
                        // cout << "        WP Band " << b+1 << ",  MAX_ABS(ALM) = " << ALM_Band.max_absalm() << endl;
                        // ALM_Band.info("bandeALM BEF");
                        ALM_Band.alm_rec(Band);
                        // char FN[200];
                        // sprintf(FN, "        BandeBF_%d.fits",b+1);
                        // sprintf(NN, "xx_band_it_%d_%d.fits", i+1, b+1);
                        // Band.write(NN);
                        
                        // Compute the coefficients and standard deviations in and out the mask
                        for (int p=0; p < Map.Npix(); p++) 
                        {
                              double Val = Band[p];
                  	       	  Band[p] = Result[p] - Val;
                   	          Result[p] = Val;
                        }
                        double Mb=Band.average();
                        Sig_InMap = Sig_Mask = 0.;
                        for (int p=0; p < Map.Npix(); p++) 
                        {
                        	  Band[p] -= Mb;
           	                  if (MaskMap[p] == 0) Sig_Mask += Band[p]*Band[p];
                    	      else Sig_InMap  += Band[p]*Band[p];
                        }
                        
                        
                        // sprintf(FN, "        BandeAF_%d.fits",b+1);
                        // Band.info(FN);
                        // sprintf(NN, "xx_bandA_it_%d_%d.fits", i+1, b+1);
                        // Band.write(NN);
                        Sig_Mask = sqrt( Sig_Mask /  (double) (Map.Npix() - NbrTotMask));
                        Sig_InMap  = sqrt(Sig_InMap  /  (double) (NbrTotMask));
                        
                        // Compute the correction coefficient
                        float Coeff =  1.;
		                if ((Sig_Mask > 0) && (Sig_InMap > 0)) Coeff = Sig_InMap / Sig_Mask;
                        if (Verbose == True) cout << "        WP:Band " << b+1 << ",   Coeff " <<  Coeff << " Sig_Mask = " << Sig_Mask << ",  Sig_InMap = " << Sig_InMap <<  "   Npix =  " <<  Map.Npix() << endl;

                        for (int p=0; p < Map.Npix(); p++)  if (MaskMap[p] == 0)  Band[p] *= Coeff;
                        
                        // coadd the band to new new solution
                        Resi += Band;
                        // add the last smooth array
                        if (b == NbrWP_Band-1) Resi += Result;
    	            } //end else b == NbrWP_Band-1
                } // endfor b
                Result = Resi;
            }   // endif Use_WP_Constraint
           else  if (GaussianityConstraint ==  True)  ALM.alm_rec(Result);
           
           // We want a zero mean solution
           if (Use_ZeroMeanConstraint == True)
           {
              double MeanRes = Result.average();
              for (int p=0; p < Map.Npix(); p++)  Result[p] -= MeanRes;
           }
           // sprintf(NN, "xx_res_it_%d.fits", i+1);
           // Result.write(NN);
  		   Result.info("RES");
     }

 
     // We reinsert the correct pixel values 
     for (int p=0; p<  Map.Npix(); p++) if (MaskMap[p] > 0) Result[p] = Map[p];
     
     // The solution as the same mean value as the input not masked area
     for (int p=0; p<  Map.Npix(); p++) Result[p] += MeanData;
     
     Result.write(Name_Imag_Out);
   
    exit(0);
}


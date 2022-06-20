/******************************************************************************
**                   Copyright (C) 1998 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Jean-Luc Starck
**
**    Date:  15/02/99
**    
**    File:  mr1d_pix.cc
**
*******************************************************************************
**
**    DESCRIPTION  analyse one pixel by multiresolution. Default is the 
**    -----------  last point.
**                 
**
******************************************************************************/

  
#include "IM_Obj.h"
#include "IM_IO.h"
#include "IM_Noise.h"
#include "MR1D_Obj.h"
#include "IM1D_IO.h"
#include "MR1D_NoiseModel.h"
#include "Usage.h"
#include "Licence.h"

extern void polint(float xa[], float ya[], int n, float x, float *y, float *dy);


char Name_Imag_In[256];
char Name_Imag_Out[256];
char Name_Coeff[256];

int NbrScale=0;           // number of scales
float NSigma=DEF_NSIG;   // number of sigma (for the noise) 
float SigmaNoise=0.;       // noise standard deviation
type_noise Stat_Noise = DEFAULT_STAT_NOISE;   // type of noise
type_trans_1d Transform=TO1_PAVE_HAAR;    // type of transform
float CCD_Gain=1.;         // CCD gain 
float CCD_ReadOutSigma=0.; // CCD read-out noise standard deviation 
float CCD_ReadOutMean=0.;  // CCD read-out noise mean  

Bool UseNSigma=False;
Bool SetPlan = False;

extern int  OptInd;
extern char *OptArg;

extern int  GetOpt(int argc, char **argv, char *opts);

Bool Verbose = False;
int Pos = -1;
 
/*********************************************************************/

static int max_scale_number (int Nc)
{
    int Nmin = Nc;
    int ScaleMax;

    ScaleMax=iround(log((float)Nmin/(float)MAX_SIZE_LAST_SCALE) / log(2.)+ 1.);
    return (ScaleMax);
}
 

/***************************************************************************/
      
 
static void usage(char *argv[])
{
    int i;

    fprintf(OUTMAN, "Usage: %s options input output\n\n", argv[0]);
    fprintf(OUTMAN, "   where options =  \n");


    fprintf(OUTMAN, "\n");
    fprintf(OUTMAN, "         [-m type_of_noise]\n");
    for (i = 0; i < NBR_NOISE-2; i++)
    fprintf(OUTMAN, "              %d: %s \n",i+1,
                                            StringNoise((type_noise )i));
    fprintf(OUTMAN, "             Default is Gaussian noise.\n\n");    
   
    gauss_usage();
    manline();  
         
    ccd_usage();
    manline();    

    nsigma_usage(NSigma);
    manline();
  
    fprintf(OUTMAN, "         [-n number_of_scales]\n");
    fprintf(OUTMAN, "             number of scales used in the multiresolution transform.\n");
    manline();

    fprintf(OUTMAN, "         [-x Position]\n");
    fprintf(OUTMAN, "             Position to analyse. Default is the last point.\n");
    manline();
        
    verbose_usage();
    manline();
    exit(-1);
}
 
/*********************************************************************/

/* GET COMMAND LINE ARGUMENTS */
static void filtinit(int argc, char *argv[])
{
    int c;
   
    /* get options */
    while ((c = GetOpt(argc,argv,"x:m:g:c:n:s:v")) != -1) 
    {
	switch (c) 
        {
	   case 'x':
	   	if (sscanf(OptArg,"%d",&Pos ) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad position: %s\n", OptArg);
	            exit(-1);
		}
		if (Pos < 2)
		{
		    fprintf(OUTMAN, "Error: bad  position: %s\n", OptArg);
		    fprintf(OUTMAN, "        1 < x <= N_pix \n");
	            exit(-1);
		}
		break;
           case 'm':
		if (sscanf(OptArg,"%d",&c ) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad type of noise: %s\n", OptArg);
	            exit(-1);
		}
                if ((c > 0) && (c <= NBR_NOISE-1)) 
                                        Stat_Noise = (type_noise) (c-1);
                else  
                {
		    fprintf(OUTMAN, "Error: bad type of noise: %s\n", OptArg);
	            exit(-1);
 		}
		break;
	    case 'g':
		/* -g <sigma_noise> */
		if (sscanf(OptArg,"%f",&SigmaNoise) != 1) 
                {
		    fprintf(OUTMAN, "bad sigma noise: %s\n", OptArg);
		    usage(argv);
		}
                Stat_Noise = NOISE_GAUSSIAN;
		break;
            case 'c':
		/* -c <gain sigma mean> */
                printf("OptArg = %s\n", OptArg);
		if (sscanf(OptArg,"%f,%f,%f", &CCD_Gain,
                                     &CCD_ReadOutSigma, &CCD_ReadOutMean) <= 0) 
                {
		    fprintf(OUTMAN, "Error: bad noise parameter: %s\n", OptArg);
		    exit(-1);
		}
                Stat_Noise = NOISE_POISSON;
                SigmaNoise  = 1.;
		break;
	   case 'n':
		/* -n <NbrScale> */
		if (sscanf(OptArg,"%d",&NbrScale) != 1) 
                {
		    fprintf(OUTMAN, "bad number of scales: %s\n", OptArg);
		    usage(argv);
		}
                if ((NbrScale <= 1) || (NbrScale > MAX_SCALE_1D)) 
                 {
		    fprintf(OUTMAN, "bad number of scales: %s\n", OptArg);
		    fprintf(OUTMAN, "1 < Nbr Scales <= %d\n", MAX_SCALE_1D);
 		    usage(argv);
		}
                SetPlan=True;
		break;
	   case 's':
		/* -s <nsigma> */
		if (sscanf(OptArg,"%f",&NSigma) != 1) 
                {
		    fprintf(OUTMAN, "bad NSigma: %s\n", OptArg);
		    usage(argv);
		}
                if (NSigma <= 0.) NSigma = DEF_NSIG;
		UseNSigma = True;
		break;
	  case 'v': Verbose = True; break;
  	    case '?': usage(argv); break;
		}
	}

        /* get optional input file names from trailing 
          parameters and open files */
	if (OptInd < argc) strcpy(Name_Imag_In, argv[OptInd++]);
	if (OptInd < argc) strcpy(Name_Imag_Out, argv[OptInd++]);

        else usage(argv);

	/* make sure there are not too many parameters */
	if (OptInd < argc)
        {
		fprintf(OUTMAN, "Too many parameters: %s ...\n", argv[OptInd]);
		usage(argv);
	}
}

/*********************************************************************/

void mksimu(fltarray & Data)
{
   int Np = 20;
   float Min = 0.;
   float Max = 100.;

   Data.alloc(Np,2);
   Data(0,0) = Min;
   Data(Np-1,0) = Max;
   Data(0,1) = Min;
   Data(Np-1,1) = Max;

   for (int i = 1; i < Np-1; i++)
   {
       float x = (int) (get_random()*(Max-Min) + Min);
       Data(i,0) = x;
       Data(i,1) = x;
   }
   cout << "Input = ";
   for (int i = 0; i < Np; i++) cout << Data(i,1) << " " << endl;
   cout << endl;
}

/*********************************************************************/

class UnbalHarrTrans: public SubBand1D 
{
  public:
  void transform_one_step(int N, float *High, float *Iw, float *Low, 
                    float *Det, float *IwLow, int Step);
  void recons_one_step(int N, float *Low, float *Det, float *High, float *Iw, float *IwLow, int Step);
  void transform(fltarray &Data, fltarray &Iw, fltarray &Mallat, fltarray &TabIw, int Nbr_Plan);
  void recons (fltarray &Mallat, fltarray &Data, fltarray &Iw, int Nbr_Plan);
  int init(fltarray &Data, fltarray &Resol0, fltarray &Iw, int Step=1);

  UnbalHarrTrans (){}
  ~UnbalHarrTrans(){}
};

/*************************************************************/

void UnbalHarrTrans::transform_one_step(int N, float *High, float *Iw, float *Low, float *Det, float *IwLow, int Step)
{
   int i;
   
//    for (i = 0; i < N; i ++)
//    {
//        int i1 = test_index(i-Step, N);
//        float w2,w1 = Iw[i];
//        if (i-Step < 0) w2 = 0;
//        else w2 =  Iw[i1];
//        IwLow[i] = w1+w2;
//        if (IwLow[i] != 0) Low[i] =  (w1*High[i] + w2*High[i1]) / IwLow[i];
//        else Low[i] = 0;
//        if (Iw[i] != 0) Det[i] =  High[i] - Low[i];
//        else Det[i] = 0;
//    }
   
   for (i = 0; i < N; i ++)
   {
       float Tc[5] ={ 1./16., 1./4., 3./8., 1./4., 1./16};
       // float Tc[5] ={ 0., 1./4., 1./2., 1./4., 0.};

       Low[i] = IwLow[i] = 0.;
       float Norm = 0;
       for (int k=0; k < 5; k++)
       {
           int Ind = test_index(i+(k-2)*Step, N);
           IwLow[i] += Iw[Ind];
           Norm += Iw[Ind]*Tc[k];
           Low[i] += Tc[k]*Iw[Ind]*High[Ind];
       }
       if (Norm  != 0) Low[i] /= Norm;
       if (Iw[i] != 0) Det[i] =  High[i] - Low[i];
       else Det[i] = 0;
   }

}

void UnbalHarrTrans::recons_one_step(int N, float *Low, float *Det, float *High, float *Iw, float *IwLow, int Step)
{
    int i;

    for (i = 0; i < N; i ++)
    {
       High[i] = Low[i] + Det[i];
    }
}


int UnbalHarrTrans::init(fltarray &Data, fltarray &Resol0, fltarray &Iw, int Step)
{
    int i;
    Bool Weight = False;
    if (Data.naxis() > 2) Weight = True;
    int Nx = Data.nx();  

cout << "Init = " << Nx << endl;

    float Min = Data(0,0);
    float Max = Data(0,0);
    for (i=0; i < Nx; i++)
    {
       if (Min > Data(i,0)) Min = Data(i,0);
       if (Max < Data(i,0)) Max = Data(i,0);
    }

    int Np = (int) ((Max - Min) / Step + 1.);
    cout << "Step = " << Step << " Min = " << Min << " Max = " << Max << " Np = " << Np <<endl;

    Resol0.alloc(Np);
    Iw.alloc(Np);
    fltarray TabNpix(Np);

    Resol0.init();
    Iw.init();
    TabNpix.init();

    cout << "Min = " << Min << " Max = " << Max << endl;

    for (i=0; i < Nx; i++)
    {
        int NexPosx = (int) ((Data(i,0)+0.5 - Min) / Step);
        if ((NexPosx < 0) || (NexPosx >= Np))
        {
           cout << "Error: Pb is NexPosx = " << NexPosx << " Np = " << Np << endl;
           exit(-1);
        }
        TabNpix(NexPosx) += 1;
        if (Weight == False) 
        {
           Iw(NexPosx) += 1.;
           Resol0(NexPosx) += Data(i,1);
	}
        else 
        {
           Iw(NexPosx) += Data(i,2);
	   Resol0(NexPosx) +=  Data(i,2)*Data(i,1);
        }
    }
   cout << "Data = " << endl;
    for (i=0; i < Np; i++)
    {
      if (TabNpix(i) > 1)
      {
          Resol0(i) /= Iw(i);
      }
      cout << Resol0(i) << " ";
    }

   cout << endl << "END simu: Np = " << Np << endl;

    return Np;
}

void UnbalHarrTrans::transform(fltarray &Data, fltarray &Iw, fltarray &Trans, fltarray &TabIw, int Nbr_Plan)
{
    int i,s,Cpt;
    int Np = Data.nx();
    TabIw.alloc(Np, Nbr_Plan+1);
    fltarray ImagLow, ImagHigh, DataResol(Np);
 
   cout << "UBH : " << Nbr_Plan << " np = " << Np << endl;
     
    Trans.alloc(Np, Nbr_Plan);
    ImagHigh.alloc(Np);
    ImagLow.alloc(Np);
    DataResol = Data;

    for (i = 0; i < Np; i++) TabIw(i,0) = Iw(i); 
    for (s = 0; s < Nbr_Plan-1; s++)
    {  
       float *PtrIw = TabIw.buffer() + s*Np;
       float *PtrIwLow = TabIw.buffer() + (s+1)*Np;
       int Step = POW2(s);

    cout << " Band " << s+1 << endl;
//     cout << " D = ";
//         for (i=0; i< Np; i++) cout << " " <<   DataResol (i);
        transform_one_step(Np, DataResol.buffer(), PtrIw,
                          ImagLow.buffer(),ImagHigh.buffer(), PtrIwLow, Step);
 
//      cout << endl << " G = ";
//         for (i=0; i< Np; i++) cout << " " <<  ImagHigh (i);
//      cout << endl << " H = ";
//         for (i=0; i< Np; i++) cout << " " <<  ImagLow (i);
//       cout << endl;
cout << endl << "  PtrIwLow = ";
         for (i=0; i< Np; i++) cout << " " <<    PtrIwLow[i];
       cout << endl;
        for (i=0; i < Np; i++)  Trans(i,s) = ImagHigh(i); 
        Cpt=0;
         for (i=0; i < Np; i++)  if (PtrIwLow[i] == 0) Cpt ++;

cout << "Cpt = " << Cpt << endl;

  	if (s != Nbr_Plan-1)
             for (i = 0; i < Np; i++)  
	               DataResol(i) = ImagLow(i);
     }
     for (i=0; i< Np; i++) Trans(i, Nbr_Plan-1) = ImagLow(i);
  cout << "end transform " << endl;
}

/****************************************************************************/

void UnbalHarrTrans::recons (fltarray &Trans, fltarray &Data, fltarray &Iw, int Nbr_Plan)
{
    register int i,s,IndE,IndO;
    int Np = Trans.nx();
    fltarray  image_h1,  image_g1, image;
    float Interpol, Err;
    cout << "Rec UnbalHarrTrans  : " << Nbr_Plan << " np = " << Np << endl;
    
    if (Data.n_elem() != Np) Data.alloc(Np);

    /* Allocation */
    image.alloc(Np);
    image_h1.alloc(Np);
    image_g1.alloc(Np);
    fltarray TabInterXE(Np);
    fltarray TabInterValE(Np);
    fltarray TabInterXO(Np);
    fltarray TabInterValO(Np);
           
     /* initial image construction : image_h1_h1. */
    for (i=0; i< Np; i++) image_h1(i) = Trans(i,Nbr_Plan-1);

    for (s = Nbr_Plan-2; s >= 0; s--)
    {
       int Step = POW2(s);
       float *PtrIw = Iw.buffer() + s*Np;
       float *PtrIwLow = Iw.buffer() + (s+1)*Np;

       for (i=0; i< Np; i++) image_g1(i) = Trans(i, s);

      cout << " Band " << s+1 << endl;
//       cout << " H = ";
// 
//         for (i=0; i< Np; i++) cout << " " <<  image_h1(i);
//       cout << endl << " G = ";
//         for (i=0; i< Np; i++) cout << " " <<  image_g1(i);

        recons_one_step (Np, image_h1.buffer(), image_g1.buffer(), 
                         image.buffer(), PtrIw, PtrIwLow, Step);
 
       for (i=0; i< Np; i++)
       {
          int d=Step;
          if (PtrIw[i] == 0)
          {
            int i1m = test_index(i-d, Np);
            int i1p = test_index(i+d, Np);
            int i2m = test_index(i-2*d, Np);
            int i2p = test_index(i+2*d, Np);
            int i3m = test_index(i-3*d, Np);
            int i3p = test_index(i+3*d, Np);
            // image(i) =  0.5*(image(i1m) + image(i1p));
            image(i) =  0.5*(-1./8*image(i3m) +1./8*image(i2m) +
              image(i1m) + image(i1p)  
              -1./8*image_h1(i3p) +1./8*image_h1(i2p));
           //image(i) =   -1./16*image(i2m) +9./16.*image(i1m) -
           //              1./16*image(i2p) +9./16.*image(i1p);
           }
        }
       

//       cout << endl << " R = ";
//         for (i=0; i< Np; i++) cout << " " << image(i);
//       cout << endl << " Iw = ";
//         for (i=0; i< Np; i++) cout << " " <<  PtrIw[i];

        for (i = 0; i < Np; i++) image_h1(i) = image(i);
//       cout << endl << " RI = ";
//         for (i=0; i< Np; i++) cout << " " << image(i);
      cout << endl;
    }
    cout << "out" << endl;
    
    for (i=0; i< Np; i++) Data(i) = image(i);
}

/*********************************************************************/


int main(int argc, char *argv[])
{
    fltarray Data;
    char Cmd[256];
    int i,b;
    // extern softinfo Soft;

    // lic_test_date();  
    // lm_check(LIC_M1D);
    // shareware_info();

    Cmd[0] = '\0';
    for (int k =0; k < argc; k++) sprintf(Cmd, "%s %s", Cmd, argv[k]);
    /* Get command line arguments, open input file(s) if necessary */
     
    filtinit(argc, argv);
    
    // fits_read_fltarr(Name_Imag_In, Data);
    // io_1d_read_data(Name_Imag_In, Data);
    mksimu(Data);
    
    if (Data.naxis() < 2)
    {
       cerr << " Error: bad file format ... naxis = " << Data.naxis() <<  endl;
       
       exit(-1);
    }

    // FitsHeader.origin = Cmd;

    if (Data.naxis()  < 2)
    {
        cout << "Error: data coordinate must be given in UnbalHarrTrans ... " << endl;
        exit(0);
    }
    Bool Weight = False;
    if (Data.naxis() > 2) Weight = True;
    int Nx = Data.nx();  

    if (!SetPlan) NbrScale = max_scale_number(Nx);
    float Step = 1.;
 
    UnbalHarrTrans UHT;
    fltarray Signal, Iw;

    fits_write_fltarr("data.fits", Data); 

    int Np = UHT.init(Data, Signal, Iw, Step);
    fits_write_fltarr("signal.fits", Signal); 
    fits_write_fltarr("iw.fits", Iw); 

    fltarray Result (Np);
    
    if (Verbose == True)
    {
       cout << endl << endl << "PARAMETERS: " << endl << endl;
       cout << "File Name in = " << Name_Imag_In << endl;
       cout << "Number of scales = " << NbrScale << endl;
       if (Stat_Noise == NOISE_GAUSSIAN)
            cout << "Type of Noise = GAUSSIAN" << endl;
       else cout << "Type of Noise = POISSON" << endl;
       cout << "Sigma Noise = " << SigmaNoise << endl;
       cout << "NSigma = " << NSigma << endl;
       cout << "naxis = " << Data.naxis() << endl;
       cout << " Output number of samples = " << Np << endl;
    }
 
    
    fltarray  Resi,Trans, TabIw;
    UHT.Border = I_CONT;
    UHT.transform(Signal, Iw, Trans, TabIw, NbrScale);
    fits_write_fltarr("tabiw.fits", TabIw); 
    fits_write_fltarr("trans.fits", Trans); 

    cout << "Reconstruction : " << endl;
    UHT.recons (Trans, Result, TabIw, NbrScale);
    io_1d_write_data(Name_Imag_Out, Result);
     
    Resi = Signal - Result;
    cout << "Min( Resi) = " <<  Resi.min() << endl;
    cout << "Max( Resi) = " <<  Resi.max() << endl;
    cout << "Sigma( Resi) = " <<  Resi.sigma() << endl;

    fits_write_fltarr("resi.fits", Trans); 

  // for (b=0; b < NbrScale-1; b++)
  
   exit(0);
} 


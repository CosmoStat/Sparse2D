/******************************************************************************
**                   Copyright (C) 1995 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Jean-Luc Starck
**
**    Date:  96/07/03
**    
**    File:  mr_visu.cc
**
*******************************************************************************
**
**    DESCRIPTION  visualize the wavelet coef. of an image
**    ----------- 
**                 
**    PARAMETRES    
**    ----------    
** 
* 
**
******************************************************************************/
 
#include "IM_Obj.h"
#include "IM_IO.h"
#include "MR_Obj.h"
#include "IM_Graphics.h"
#include "IM_Sigma.h"

#define GRAY 1
#define CONT 2 
#define PERS 3

char Name_Imag_In[256]; /* input file image */
char Name_Imag_Out[256]; /* output file name */
char Name_Write_PS[256]; /* PS file */

int Nbr_Plan=DEFAULT_NBR_SCALE;  /* number of scales */
type_transform Transform=DEFAULT_TRANSFORM; /* type of transform */
 
extern int  OptInd;
extern char *OptArg;

extern int  GetOpt(int argc, char *const*argv, char *opts);
   
Bool NormCoef = True;
int  TypeVisu = GRAY;
Bool VisuColor = True;
Bool WritePS = False;
int LineInc=1;
float NSigma = 3.;
    
Bool OptMR=False;
Bool Verbose = False;

sb_type_norm Norm = NORM_L1;
type_sb_filter SB_Filter = F_MALLAT_7_9;

/****************************************************************************/

static void usage(char *argv[])
{

    fprintf(OUTMAN, "Usage: %s options image_or_MRfile out_image\n\n", argv[0]);
    fprintf(OUTMAN, "   where options =  \n");

    transform_usage(Transform);
    manline();
    sb_usage(SB_Filter);
    manline();    
    nbr_scale_usage(Nbr_Plan);
    manline();
 
    mrvisu_option_usage(NSigma);
    manline();
    vm_usage();
    manline();
    verbose_usage(); 
     manline();
     manline();
     exit(-1);
}
 
/*********************************************************************/

/* GET COMMAND LINE ARGUMENTS */
static void filtinit(int argc, char *argv[])
{
    int c;
    Bool OptL = False, Optf = False;
    Bool Opti=False;
    Bool Opts=False;
#ifdef LARGE_BUFF
    int VMSSize=-1;
    Bool OptZ = False;
    char VMSName[1024] = "";
#endif    
    /* get options */
    while ((c = GetOpt(argc,argv,"t:T:Ln:bcV:w:i:s:vzZ:")) != -1) 
    {
	switch (c) 
        {
	   case 'v': Verbose = True; break;
	   case 't':
		/* -d <type> type of transform */
		if (sscanf(OptArg,"%d",&c ) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad type of multiresolution transform: %s\n", OptArg);
	            exit(-1);
                    
		}
                if ((c > 0) && (c <= NBR_TRANSFORM)) 
                                        Transform = (type_transform) (c-1);
                else  
                {
		    fprintf(OUTMAN, "Error: bad type of transform: %s\n", OptArg);
	            exit(-1);
 		}
 		OptMR=True;
 		break;
	   case 'T': 
		Optf = True;
		SB_Filter = get_filter_bank(OptArg);
		break;
	   case 'L': Norm = NORM_L2; OptL = True; break;
	   case 'n':
		/* -n <Nbr_Plan> */
		if (sscanf(OptArg,"%d",&Nbr_Plan) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad number of scales: %s\n", OptArg);
		    exit(-1);
		}
                if ((Nbr_Plan <= 1) || (Nbr_Plan > 9)) 
                 {
		    fprintf(OUTMAN, "Error: bad number of scales: %s\n", OptArg);
		    fprintf(OUTMAN, "       1 < Nbr Scales <= %d\n", 9);
 		    exit(-1);
		}
		OptMR=True;
 		break;
 	   case 'b':
 	        VisuColor = False;
 	        break;
 	   case 'c':
 	        NormCoef = False;
 	        break;
 	   case 'V':
 	        if (sscanf(OptArg,"%d",&c ) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad type of visualization: %s\n", OptArg);
	            exit(-1);
 		}
                if (c == 1) TypeVisu = GRAY;
                else if (c == 2) TypeVisu = CONT;
                else if (c == 3) TypeVisu = PERS;
                else 
                {
                    fprintf(OUTMAN, "Error: bad type of visualization parameter: %s\n", OptArg);
	            exit(-1);
 		}
 		break;		   
 	   case 's':
  		/* -s <nsigma> */
		if (sscanf(OptArg,"%f",&NSigma) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad NSigma: %s\n", OptArg);
		    exit(-1);
		}
                if (NSigma <= 0.) NSigma = DEFAULT_N_SIGMA;
                Opts = True;
		break;
            case 'i':
 	        if (sscanf(OptArg,"%d",&LineInc ) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad increment parameter: %s\n", OptArg);
	            exit(-1);
 		}
                if (LineInc < 1)  
                {
                    fprintf(OUTMAN, "Error: bad increment parameter: %s\n", OptArg);
	            exit(-1);
 		}
 		Opti = True;
 		break;	   
 	    case 'w':
		/* -w < support file name> */
		if (sscanf(OptArg,"%s",Name_Write_PS) != 1) 
                {
		   fprintf(OUTMAN, "Error: bad file name: %s\n", OptArg);
		   exit(-1);
		}
                WritePS = True;
                break;
#ifdef LARGE_BUFF
	    case 'z':
	        if (OptZ == True)
		{
                   fprintf(OUTMAN, "Error: Z option already set...\n");
                   exit(-1);
                }
	        OptZ = True;
	        break;
            case 'Z':
	        if (sscanf(OptArg,"%d:%s",&VMSSize, VMSName) < 1)
		{
		   fprintf(OUTMAN, "Error: syntaxe is Size:Directory ... \n");
		   exit(-1);
		}
	        if (OptZ == True)
		{
                   fprintf(OUTMAN, "Error: z option already set...\n");
                   exit(-1);
                }
		OptZ = True;
                break;
#endif
            case '?': usage(argv); break;
	    default: usage(argv); break;
		}
	}

       /* get optional input file names from trailing 
          parameters and open files */
	if (OptInd < argc) strcpy(Name_Imag_In, argv[OptInd++]);
         else usage(argv);

	if (OptInd < argc) strcpy(Name_Imag_Out, argv[OptInd++]);
         else usage(argv);

	/* make sure there are not too many parameters */
	if (OptInd < argc)
        {
		fprintf(OUTMAN, "Error: too many parameters: %s ...\n", argv[OptInd]);
		exit(-1);
	}

 
  	// if the transform is no isotrop,   
 	if (isotrop(Transform) == False) 
 	{
  	   if ((TypeVisu == CONT) || (TypeVisu == PERS))  
 	   {
 	      cerr << "Error: this visualization is not possible with this transform." << endl;
 	      exit(-1);
 	   }
 	}
 	
 	if ((Opti == True) && (TypeVisu != PERS))
 	{
 	   fprintf(OUTMAN, "\nError: -i option is not valid with this visualization ...\n");
           exit(-1);
	}
	if ((Opts == True) && (TypeVisu == GRAY ))
 	{
 	   fprintf(OUTMAN, "\nError: -s option is not valid with this visualization ...\n");
           exit(-1);
	}
	if ((VisuColor == False) && (TypeVisu == GRAY ))
 	{
 	   fprintf(OUTMAN, "\nError: -b option is not valid with this visualization ...\n");
           exit(-1);
	}
	if ((NormCoef == False) && (TypeVisu != GRAY ))
 	{
 	   fprintf(OUTMAN, "\nError: -c option is not valid with this visualization ...\n");
           exit(-1);
	}
	if ((Transform != TO_UNDECIMATED_MALLAT) && (Transform != TO_MALLAT) && ((OptL == True) || (Optf == True)))
	{
	   fprintf(OUTMAN, "Error: option -T and -L are only valid with Mallat transform ... \n");
           exit(0);
	}

#ifdef LARGE_BUFF
    if (OptZ == True) vms_init(VMSSize, VMSName, Verbose);
#endif

}


/****************************************************************************/

void mr_coeff_in_pers (MultiResol & MR_Data, Ifloat & Result, int Increment, Bool Color,
        float KSigma)
{  
   int MAXX = 512;
   int MAXY = 512;
   int s,i,j,Size,Nl,Nc,Nls,Ncs,Nli,Nci;
   int NbrScale = MR_Data.nbr_band();
   int Inc = (Increment > 0) ? Increment: 1;
   int Mod_Visu = (Color == True) ? 1: 0;;
   float Zoom, Sigma;
   int Size_Plan;
   Ifloat Imag, Caval;
   int Offsx,Offsy;
   
   Nl = MR_Data.size_ima_nl();
   Nc = MR_Data.size_ima_nc();       
   Zoom = (float) MAXX / (float) Nc;
   Size_Plan = (int) ((float) MAXY / (float) NbrScale);       
   Caval.alloc(Nl,Nc, "Perspec");
   Size = MAXX * MAXY;    
   Result.alloc(MAXY, MAXX, "Result");
   
   for (s = 0; s < NbrScale ; s++)
   {
        Imag = MR_Data.band(s);
        Nls = MR_Data.size_band_nl(s);
        Ncs = MR_Data.size_band_nc(s);
        Sigma = KSigma * sigma(Imag);
 
        Nli = Size_Plan;
        Nci = (int) ((((float) Ncs * Zoom) > MAXX) ? MAXX : ((float)Ncs * Zoom));
        Caval.resize(Nli,Nci);
        if (s !=  NbrScale-1)
        for (j = 0; j < Nl*Nc; j++)
        {
            if (Imag (j) > Sigma) Imag (j) = Sigma;
            if (Imag (j) < - Sigma) Imag (j) = - Sigma;
        }
        im_3dview_2 (Imag, Caval, Inc, Mod_Visu);  
        Offsx = (MAXX - Nci) / 2;
        Offsy = Nli * s;
        for (i = 0; i < Nli; i++)
        for (j = 0; j < Nci; j++)  Result(Offsy + i, Offsx + j) = Caval(i,j);
     } 
}


/****************************************************************************/

void mr_coeff_in_cont (MultiResol & MR_Data, Ifloat & Result, Bool Color,
        float KSigma)
{  
   int s,i,j,Nl,Nc,Nls,Ncs;
   int NbrScale = MR_Data.nbr_band();
   int Mod_Visu;
   float Sigma=0;
   // float Mean;
   Ifloat Imag, Cont;
    
   Nl = MR_Data.size_ima_nl();
   Nc = MR_Data.size_ima_nc();       
   Cont.alloc(Nl,Nc, "Cont");
   Result.alloc(Nl,Nc,"Result");
   Ifloat Interp(Nl,Nc,"Interp");
   
   for (s = 0; s < NbrScale -1; s++)
   {
        Imag = MR_Data.band(s);
        if (s == 0) Sigma = sigma(Imag)*KSigma;
        else Sigma /= 4.;
        
        Nls = MR_Data.size_band_nl(s);
        Ncs = MR_Data.size_band_nc(s);
        // sigma_clip(Imag, Mean, Sigma);
        // if (s == 0) Sigma *= (KSigma+1); 
        // else Sigma *= KSigma;     
        for (j = 0; j < Nl*Nc; j++)
        {
            if (Imag (j) > Sigma) Imag (j) = Sigma;
            if (Imag (j) < - Sigma) Imag (j) = - Sigma;
        }
        if ((Nls != Nl) || (Ncs != Nl))
        {
           im_bilinear_interp(Imag, Interp);
           im_isophot (Interp, Cont, (float) 0., (float) (Sigma + Sigma/10.), (float) Sigma, (int) 0);
           // im_quick_isophot (Interp, Cont,   Sigma,  1.5*Sigma,  Sigma, 0);
        }
        else im_isophot (Imag, Cont, 0., Sigma + Sigma/10., Sigma, 0);
        Mod_Visu = (Color == False) ? 1: NbrScale-s;
        
        for (i = 0; i < Nl; i++)
        for (j = 0; j < Nc; j++) 
              if (Cont(i,j) != 0) Result(i,j) = Cont(i,j)* Mod_Visu;
      }
}

/****************************************************************************/

void mr_coeff_in_ima(MultiResol & MR_Data, Ifloat & Result)
{
   int s,i,j,Size,Nl,Nc,Nls,Ncs,Nlo,Nco;
   int NbrScale = MR_Data.nbr_band();
   int Depi=0,Depj=0;
   float Sq;
   
   Nl = MR_Data.size_ima_nl();
   Nc = MR_Data.size_ima_nc();   
   
   switch(MR_Data.Set_Transform)
    {
        case TRANSF_UNDECIMATED_MALLAT:
        case TRANSF_PAVE:
	case TRANSF_DIADIC_MALLAT:
	       Sq = sqrt( (float) NbrScale);
	       if (Sq - ((int) Sq) == 0) Size = (int) Sq;
	       else Size = (int) (Sq+1);
               Nlo = Nl * Size;
               Nco = Nc * Size;
               Result.alloc(Nlo,Nco, "Result");  
	       Depi = Depj;
               for (s = 0; s < NbrScale; s++)
               {  
                 for (i=0;i<Nl;i++)
                 for (j=0;j<Nc;j++) Result(Depi+i,Depj+j) = MR_Data(s,i,j);
		 Depj += Nc;
		 if (Depj >= Nco)
		 {
		    Depj = 0;
		    Depi += Nl;
		 }
               }
               /*
               for (i = 0; i < Nl; i++)
               {
                  Imag(Nlo/2,i) = 1.;
                  Imag(i,Nco/2) = 1.;
               } */
               break;
	case TRANSF_SEMIPYR:
        case TRANSF_PYR:
	       if (MR_Data.Set_Transform == TRANSF_SEMIPYR) Size = 3;
	       else Size = 2;
               Nlo = Nl * Size;
               Nco = Nc * Size;
               Result.alloc(Nlo,Nco, "Result");
               Depi = Depj = 0;
               for (s = 0; s < NbrScale; s++)
               {
                  Nls = MR_Data.size_band_nl(s);
                  Ncs = MR_Data.size_band_nc(s);
                  for (i=0;i<Nls;i++)
                  for (j=0;j<Ncs;j++) Result(Depi+i,Depj+j) = MR_Data(s,i,j);
                  Depi += Nls;
                  Depj += Ncs;
               }
               break;
        case TRANSF_MALLAT:
        case TRANSF_FEAUVEAU:
	        ortho_trans_to_ima(MR_Data, Result);
                break;
        default: break;
       }
}

/****************************************************************************/

int main(int argc, char *argv[])
{
    int k;
    Ifloat Data, Result;
      
    /* Get command line arguments, open input file(s) if necessary */
    lm_check(LIC_MR1);
    filtinit(argc, argv);

    if (Verbose == True)
    {
       cout << endl << endl << "PARAMETERS: " << endl << endl;
       cout << "File Name in = " << Name_Imag_In << endl;
       cout << "File Name Out = " << Name_Imag_Out << endl;
       cout << "Transform = " << StringTransform(Transform) << endl;
       cout << "Number of scales = " << Nbr_Plan << endl;
    }
    MultiResol MR_Data;
    if ( (strstr(Name_Imag_In, ".mr") != NULL) || (strstr(Name_Imag_In, ".MR") != NULL))
    {	
        if (OptMR == True)
 	{
 	   fprintf(OUTMAN, "Error: t,n options are not valid with multiresolution input files ...");
           exit(-1);
	}
       MR_Data.read (Name_Imag_In);
       Transform= MR_Data.Type_Transform;
       Nbr_Plan = MR_Data.nbr_scale();
       if ((Nbr_Plan > 9) && (isotrop(Transform) == True))
       {
          cerr << "Error: too many scales in the transform for visualization" << endl;
          exit(-1);
       }
        
  	// if the transform is no isotrop,   
 	if (isotrop(Transform) == False) 
 	{
  	   if ((TypeVisu == CONT) || (TypeVisu == PERS))  
 	   {
 	      cerr << "ERROR: this visualization is not possible with this transform." << endl;
 	      exit(-1);
 	   }
 	}
    }
    else
    {
       io_read_ima_float(Name_Imag_In, Data);
       check_scale(Data.nl(), Data.nc(), Nbr_Plan);
       FilterAnaSynt FAS;
       FilterAnaSynt *PtrFAS = NULL;
       if ((Transform == TO_MALLAT) || (Transform == TO_UNDECIMATED_MALLAT))
       {
           FAS.Verbose = Verbose;
           FAS.alloc(SB_Filter);
           PtrFAS = &FAS;
       }
       MR_Data.alloc (Data.nl(), Data.nc(), Nbr_Plan, Transform, PtrFAS, Norm);
       MR_Data.transform(Data);
    }
    
    if ((NormCoef == True) && (TypeVisu == GRAY)) MR_Data.norm((float) 255., True);
    switch(TypeVisu)
    {
       case GRAY:  
                  mr_coeff_in_ima( MR_Data,  Result); 
// cout << "Nl = " << Result.nl() << " Nc = " << Result.nc() <<  endl;
                  io_write_ima_float(Name_Imag_Out, Result);
		  break;
       case CONT: mr_coeff_in_cont (MR_Data, Result, VisuColor, NSigma); 
                  io_write_ima_float(Name_Imag_Out, Result, NULL, True);
		  break;
       case PERS: 
           mr_coeff_in_pers (MR_Data, Result, LineInc, VisuColor, NSigma);
	   io_write_ima_float(Name_Imag_Out, Result, NULL, True);
           break;
       default: break;
    }
     
    if (WritePS == True) k=im_create_ps_file (Result, Name_Write_PS);
    exit(0);
} 


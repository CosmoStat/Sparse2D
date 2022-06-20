/******************************************************************************
**                   Copyright (C) 1996 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 3.2
**
**    Author: Jean-Luc Starck 
**
**    Date:  96/06/13
**    
**    File:  mr_psupport.cc
**
*******************************************************************************
**
**    DESCRIPTION  perform the multi-resolution support using 
**		   the abaque.
**    ----------- 
**                 
**    PARAMETRES    
**    ----------    
** 
**    USAGE: mr_filter option image output
**        where options = 
**		[-a ascii_file]
**
**		[-i image_file]
**
**		[-w]
**
**              [-f abaque_file_name]
**                  default is Abaque.fits
**
**		a & i options can't be used together and one of both
**		must be set
**
**
**   creates a file (filename in option) which contains the significant 
**   reduced wavelet coefficients.
**   if -w  then creates 
**     Wavelet.mr : contains the wavelet transform of the image
**                  (but not the reduced coefficient).
**     
****************************************************************************/


// static char sccsid[] = "@(#)mr_psupport.cc 3.2 96/06/13 CEA 1995 @(#)";
 

#include <fstream>

#include "IM_Obj.h"
#include "IM_IO.h"
#include "MR_Obj.h"
#include "MR_Abaque.h"
#include "MR_NoiseModel.h"
#include "MR_Psupport.h"

char Name_Imag_In[80]; /* input file image */
static Bool setopt = False; /* make sure of mutual exclusivity of options incommand line */
static Bool ascii_opt = False;
char Name_Support[80];
char AnaFileName[80];
char TabAnaFileName[80];

extern int  OptInd;
extern char *OptArg;

extern int  GetOpt(int argc, char *const*argv, char *opts);

Bool WriteAll = False;
Bool KeepPositivSup = False;

char Name_Abaque[80];
int Abaque_N_Scale=0;
int NbrScale = DEF_N_SCALE;
int MinEvent = DEF_MINEVENT_NUMBER;
int FirstScale = DEF_FIRST_DETECT_SCALE;
Bool ObjAna = False;
Bool TabAsciiRes = False;
Bool AsciiRes = False;
Bool Verbose = False;

/***********************************************************************/

static void usage(char *argv[])
{

    fprintf(OUTMAN,"\n\nUsage: %s options in_image out_mr_file\n\n", argv[0]);
    fprintf(OUTMAN, "   where options =  \n");

    input_poisson_usage();
    manline();
    
    first_detect_scale_usage();
    manline();
    
    min_event_usage(MinEvent);   
    manline();
 
    write_wave_mr_usage();
    manline();
    
  
    signif_ana_usage();
    manline();
    
    ascii_signif_ana_usage();
    manline();
    
    detect_pos_usage();
    manline();

    abaque_file_usage();
    manline();

    nbr_scale_usage(NbrScale);
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
static void psupportinit(int argc, char *argv[])
{
    int c;
#ifdef LARGE_BUFF
    int VMSSize=-1;
    Bool OptZ = False;
    char VMSName[1024] = "";
#endif    
   
    strcpy(Name_Abaque, Name_Abaque_Default);

    /* get options */
    while ((c = GetOpt(argc,argv,"F:a:I:q:wn:e:ps:t:vzZ:")) != -1) 
    {
	switch (c) 
      	{
	   case 'v': Verbose = True; break;
	   case 'a': /* ascii_filename */
		if ( (sscanf(OptArg,"%s",Name_Imag_In) != 1) && (setopt == True) )
                {
		   fprintf(OUTMAN, "Error: bad file name: %s\n", OptArg);
		   exit(-1);
		}
		setopt = True;
		ascii_opt = True;
		break;
	 case 'p':
		 KeepPositivSup=True;
		break; 
	   case 'F': /* begin the detection */
		if (sscanf(OptArg,"%d", &FirstScale) != 1) 
                {
		   fprintf(OUTMAN, "Error: bad first detection scale number: %s\n", OptArg);
		   exit(-1);
		}
		FirstScale --;
		break;
	 case 'e':
 		/* -n <Nbr_Plan> */
		if (sscanf(OptArg,"%d",&MinEvent) != 1) 
                {
		    fprintf(OUTMAN, "Error: Min number of events: %s\n", OptArg);
		    exit(-1);
		}
                if (MinEvent < 1)  
                {
		    fprintf(OUTMAN, "Error: bad number of events min: %s\n", OptArg);
		    exit(-1);
		}
		break;
 	   case 'I': /* image_filename */
		if ( (sscanf(OptArg,"%s",Name_Imag_In) != 1) && (setopt == True) )
                {
		   fprintf(OUTMAN, "Error: bad file name: %s\n", OptArg);
		   exit(-1);
		}
		setopt = True;
		break;
 	   case 's': /* image_filename */
		if (sscanf(OptArg,"%s",AnaFileName) != 1)
                {
		   fprintf(OUTMAN, "Error: bad file name: %s\n", OptArg);
		   exit(-1);
		}
		ObjAna = True;
		AsciiRes = True;
		break;
 	   case 't': /* image_filename */
		if (sscanf(OptArg,"%s",TabAnaFileName) != 1)
                {
		   fprintf(OUTMAN, "Error: bad file name: %s\n", OptArg);
		   exit(-1);
		}
		ObjAna = True;
		TabAsciiRes = True;
		break;
	 case 'w':
		 WriteAll=True;
		break; 
	 case 'q':
		if (sscanf(OptArg,"%s",Name_Abaque) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad abaque parameter: %s\n", OptArg);
		    exit(-1);
		}
		break;
	 case 'n':
 		/* -n <Nbr_Plan> */
		if (sscanf(OptArg,"%d",&NbrScale) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad number of scales: %s\n", OptArg);
		    exit(-1);
		}
                if ((NbrScale <= 1) || (NbrScale > MAX_SCALE)) 
                {
		    fprintf(OUTMAN, "Error: bad number of scales: %s\n", OptArg);
		    fprintf(OUTMAN, "       1 < Nbr Scales <= %d\n", MAX_SCALE);
		    exit(-1);
		}
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


   /* get output filename for mr_support */
   	if (OptInd < argc)  strcpy(Name_Support,argv[OptInd++]);
        else usage(argv);

   /* make sure there are not too many parameters */
	if (OptInd < argc) {
		fprintf(OUTMAN, "Error: too many parameters: %s ...\n", argv[OptInd]);
		exit(-1);
	}

  if (setopt == False)
  {
      fprintf(OUTMAN, "\n\nError: option -a or -I has to be set\n");
      exit(-1);
  }
#ifdef LARGE_BUFF
    if (OptZ == True) vms_init(VMSSize, VMSName, Verbose);
#endif

}



/************************************************************************/

int main(int argc, char *argv[])
{
 int Nl,Nc,Nls,Ncs,i,j,s;	
 type_transform Transform = TO_PAVE_BSPLINE;
 Ifloat Image;
 Iint EI;
 Ifloat Event_Image;
 Ifloat I_Abaque;
 type_border Border = I_MIRROR;

 /* get command line */
 lm_check(LIC_MR1);
 psupportinit(argc,argv);

 io_read_ima_float(Name_Abaque, I_Abaque);
 // StatEventPoisson EvPoisson;
 // EvPoisson.compute_distrib(False);
 // EvPoisson.find_threshold(0.0001);      
 // EvPoisson.Threshold = I_Abaque;

 Abaque_N_Scale = I_Abaque.nl();
 fprintf(stdout, "Abaque: number max of events = 2^%d\n", Abaque_N_Scale);
 io_set_format(F_UNKNOWN);

 /* build the initial image of events and the initial image filtered */
 /* by the scale function */
 if (ascii_opt==True)  
 {
     building_imag_ascii(Name_Imag_In, EI, Image);
     for (i=0;i < Image.nl(); i++)
     Event_Image.alloc(Image.nl(), Image.nc(), "build_img");
     for (j=0;j < Image.nc(); j++) Event_Image(i,j) = EI(i,j);
 }
 else   
 {
    io_read_ima_float(Name_Imag_In, Image);
    building_imag_imag(Image, Event_Image);
 }

 /* compute the input signal support */
 Nl = Image.nl();
 Nc = Image.nc();
 check_scale(Nl, Nc, NbrScale);

 MultiResol MR_Data(Nl,Nc,NbrScale, Transform,"mr_psupport");
 MRNoiseModel ModelNoise(NOISE_EVENT_POISSON, Nl,Nc,NbrScale, Transform);

 /* compute the wavelet transform */
 MR_Data.Border = Border;
 MR_Data.transform(Image);
 
 ModelNoise.MinEventNumber = MinEvent;
 ModelNoise.MinEvent = True;
 ModelNoise.OnlyPositivDetect = KeepPositivSup;
 ModelNoise.TransImag = False;
 ModelNoise.SigmaNoise = 1.;
 ModelNoise.FirstDectectScale = FirstScale;
 ModelNoise.Event_Image = Event_Image;

 mr_psupport(Event_Image, MR_Data, I_Abaque, ModelNoise, Border, WriteAll);
 
 ModelNoise.threshold(MR_Data);
 /* writing on disk the support */
 MR_Data.write(Name_Support);

  /* segmentation */
  if (ObjAna == True)
  {
      int Cpt,ind,nmax;
      float Signif;
      ofstream Fic;
      ofstream Fic1;
      int *TabPeri, *TabSurf, *TabXmax, *TabYmax;
      float MeanMorpho, Morpho;
      float *TabMorpho, *TabValMax;
      float *TabMx, *TabMy, *TabMean, *TabVx;
      float *TabVy, *TabVxy, *TabSx, *TabSy, *TabAngle;
 
      if (AsciiRes == True)
      {
         Fic.open(AnaFileName);
         if(!Fic)
         {
             cerr<<"Error: Unable to open "<< AnaFileName << endl;
             exit(1);
         }
      }
      if (TabAsciiRes == True)
      {
         Fic1.open(TabAnaFileName);
         if(!Fic)
         {
             cerr<<"Error: Unable to open "<< TabAnaFileName << endl;
             exit(1);
         }
      }

      for (s = 0; s < MR_Data.nbr_band()-1; s++)
      {
          Nls = MR_Data.size_band_nl(s);
          Ncs = MR_Data.size_band_nc(s);
          Image.resize(Nls, Ncs);
          Image = MR_Data.band(s);
          im_segment (Image, MR_Data.band(s), nmax, FLOAT_EPSILON);

          if (AsciiRes == True)
          {
             Fic << endl;
             Fic << "SCALE " << s+1 << ": Number of structures = " << nmax << endl;
          }

          if (nmax > 0)
          {
          TabPeri = new int [nmax+1];
          TabSurf = new int [nmax+1];
          TabMorpho = new float [nmax+1];
          TabXmax = new int [nmax+1];
          TabYmax = new int [nmax+1];
          TabValMax  = new float [nmax+1];
          TabMx  = new float [nmax+1];
          TabMy  = new float [nmax+1];
          TabMean  = new float [nmax+1];
          TabVx  = new float [nmax+1];
          TabVy  = new float [nmax+1];
          TabVxy  = new float [nmax+1];
          TabSx  = new float [nmax+1];
          TabSy  = new float [nmax+1];
          TabAngle  = new float [nmax+1];

          Cpt = 0;
          for (i = 0; i <= nmax; i++) 
          {
              TabPeri[i] = TabSurf[i] = 0;
              TabMorpho[i] = 0.;
              TabXmax[i] = -1;
              TabYmax[i] = -1;
              TabValMax [i] = 0.;
              TabMx[i] = 0.;
              TabMy[i] = 0.;
              TabMean[i] = 0.;
              TabVx[i] = 0.;
              TabVy[i] = 0.;
              TabVxy[i] = 0.;
              TabSx[i] = 0.;
              TabSy[i] = 0.;
              TabAngle[i] = 0.;
          }
          for (i = 0; i < Nls; i++)
          for (j = 0; j < Ncs; j++)
          if (MR_Data(s,i,j) != 0)
          {
              float Coef;

              Coef = Image(i,j);
              Cpt ++;
              ind = (int) (MR_Data(s,i,j) + 0.5);
              TabSurf[ind] ++;
              TabMx[ind] += Coef*j;
              TabMy[ind] += Coef*i;
              TabMean[ind] += Coef;
              TabVx[ind] += Coef*j*j;
              TabVy[ind] += Coef*i*i;
              TabVxy[ind] +=  Coef*i*j;

              if (ABS(Image(i,j)) > TabValMax[ind])
              {
                 TabValMax[ind] = ABS(Image(i,j));
                 TabXmax[ind] = j;
                 TabYmax[ind] = i;
              }
              if (i==0) TabPeri[ind]++;
              if (i==Nl-1) TabPeri[ind]++;                 
              if (j==0) TabPeri[ind]++;
              if (j==Nc-1) TabPeri[ind]++;                 

              if ( (i>0) && (j>0) && (i<Nl-1) && (j<Nc-1)) 
              {
                if (MR_Data(s,i-1,j) < 0.5) TabPeri[ind]++;
                else if (MR_Data(s,i+1,j) < 0.5) TabPeri[ind]++;
                else if (MR_Data(s,i,j-1) < 0.5) TabPeri[ind]++;
                else if (MR_Data(s,i,j+1) < 0.5) TabPeri[ind]++;
              }
          }

          Signif =  (float) Cpt / (float)(Nls*Ncs) * 100.;
          if (AsciiRes == True)
          {
             Fic << "   Pourcentage of Significant Coefficients = " << Signif << endl;
          }

          MeanMorpho = 0.;
          for (i = 1; i <= nmax; i++)
          {
             float Val;
             Val = (float) (TabPeri[i]*TabPeri[i]);
             if (Val > FLOAT_EPSILON) Morpho = 4. * PI * TabSurf[i] / Val;
             else Morpho = 0.;
             MeanMorpho += Morpho;
             TabMorpho[i] = Morpho;

             TabMx[i] /= TabMean[i];
             TabMy[i] /= TabMean[i];
             TabVx[i] /= TabMean[i];
             TabVy[i] /= TabMean[i];
             TabVxy[i] /= TabMean[i];
             TabVx[i] -= TabMx[i]*TabMx[i];
             TabVy[i] -= TabMy[i]*TabMy[i];
             TabVxy[i] -= TabMx[i]*TabMy[i];
             if (ABS(TabVx[i]-TabVy[i])<1E-7 || ABS(TabVxy[i])<1E-9)
             {
                TabSx[i] = sqrt(ABS(TabVx[i]));
                TabSy[i] = sqrt(ABS(TabVy[i]));
                TabAngle[i] = 0.;
             }
             else
             {
                TabAngle[i] = 0.5 * atan(2*TabVxy[i]/(TabVx[i]-TabVy[i]));
                TabSx[i] = sqrt( ABS(
                   (TabVx[i]+TabVy[i])/2. + TabVxy[i]/sin(2.*TabAngle[i]) ));
                TabSy[i] = sqrt( ABS(
                    (TabVx[i]+TabVy[i])/2. - TabVxy[i]/sin(2.*TabAngle[i]) ));
             }
          }
          if (nmax > 0) MeanMorpho /= (float) nmax;

          if (AsciiRes == True)
          {
             Fic << "   Mean deviation of shape from sphericity = " << MeanMorpho << endl << endl;
          }

          for (i = 1; i <= nmax; i++)
          {
             if (AsciiRes == True)
             {

                Fic << "   S" << s+1 << ":" << i << " Surf = " << TabSurf[i] << "  Peri = " << TabPeri[i] << "  Morpho = " << TabMorpho[i] << endl;
                Fic << "       Angle = " << TabAngle[i] << "  SigmaX = " << TabSx[i] << "  SigmaY = " << TabSy[i] << endl;
             }
             if (TabAsciiRes == True)
             {
                Fic1 << s+1 <<" " << i <<" " << TabXmax[i] <<" " << TabYmax[i] <<" " << TabSurf[i] <<" " << TabPeri[i] <<" " << TabMorpho[i] << " " << TabAngle[i] << " " << TabSx[i] << " " << TabSy[i] << endl;
             }
          }
          delete [] TabPeri;
          delete [] TabSurf;
          delete [] TabMorpho;
          delete [] TabXmax;
          delete [] TabYmax;
          delete [] TabValMax;
          delete [] TabMx;
          delete [] TabMy;
          delete [] TabMean;
          delete [] TabVx;
          delete [] TabVy;
          delete [] TabVxy;
          delete [] TabSx;
          delete [] TabSy;
          delete [] TabAngle;
          }
         /* write results */
       }
       if (AsciiRes == True) Fic.close();
       if (TabAsciiRes == True) Fic1.close();
       MR_Data.write("xx_Segment.mr");
  }
  exit(0);
}

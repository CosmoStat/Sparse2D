/*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Jean-Luc Starck
**
**    Date:  20/03/00 
**    
**    File: Filter.cc
**
**    Modification history:
**
*******************************************************************************
**
**    DESCRIPTION  Filter design
**    -----------  
**                 
******************************************************************************/
 
#include "MSVST_Filter.h"

// Antoni 7/9 filters
/*
static float Filter7_9_h0[9] =
{
	 0.02674875741,
   	-0.0168641184 ,
	-0.0782232665 ,
	 0.26686411844,
	 0.60294901823,
	 0.26686411844,
	-0.0782232665 ,
	-0.0168641184 ,
	 0.02674875741
};

static float Filter7_9_g0[7] =
{
	 0.04563588155 ,
	-0.02877176311 ,
    -0.295635881557,
	 0.557543526229,
    -0.295635881557,
   	-0.02877176311 ,
	 0.04563588155 
};

static float Filter7_9_h1[7] = 
{
	-0.04563588155,
	-0.02877176311,
    0.295635881557,
    0.557543526229, 
    0.295635881557,
   	-0.02877176311,
	-0.04563588155
};

static float Filter7_9_g1[9] =
{
	 0.02674875741,
   	 0.0168641184 ,
	-0.0782232665 ,
	-0.26686411844,
	 0.60294901823,
	-0.26686411844,
	-0.0782232665 ,
	 0.0168641184 ,
	 0.02674875741
};
*/

// \bar{h} and \tilde{h}
// \bar{g} and \tilde{g}

static float HaarAnalysis[2] = {
    1. / 2.,
    1. / 2.
}; // same with HaarSynthesis

static float Haar2Synthesis[6] = {
   -1. / 16,
    1. / 16,
    1. / 2,   
    1. / 2,
    1. / 16,
   -1. / 16
};

static float Haar4Synthesis[10] = {
   3. / 256.,
  -3. / 256.,
  -22. / 256.,
   22. / 256.,
   1. / 2.,
   1. / 2.,
   22. / 256.,
  -22. / 256.,
  -3. / 256.,
   3. / 256.
};

static float U_HAAR_B3SAnalysis[2] = {
  0.5,
  0.5
};

static float U_HAAR_B3SAnalysis2[2] = {
  -0.5,
  0.5
};

static float U_HAAR_B3SSynthesis[6] = {
  1./16.,
  1./4.,
  3./8.,
  1./4.,
  1./16.,
  0
};

static float U_HAAR_B3SSynthesis2[6] = {
  1./ 16.,
  6./ 16.,
  1.,
  -6./ 16.,
  -1./ 16.,
  0
};

static float U_HAAR_B3S2Analysis[2] = {
  0.5,
  0.5
};

static float U_HAAR_B3S2Analysis2[2] = {
  -0.5,
  0.5
};

static float U_HAAR_B3S2Synthesis[6] = {
  -1./16., 
  -1./8., 
  1., 
  1./8., 
  1./16.,
  0
};

static float U_HAAR_B3S2Synthesis2[6] = {
  -1./ 16.,
  -1./4.,
  5./8.,
  -1./ 4., 
  -1./ 16.,
  0
};

static float F5_3Analysis[5] = {
   -0.125,
    0.25,
    0.75,
    0.25,
   -0.125
};

static float F5_3Synthesis[3] = {
  0.25,
  0.5,
  0.25
};

static float Daub2Analysis [4] = {
                 (1 - sqrt(3.)) / (4 * sqrt(2.)),
                 (3 - sqrt(3.)) / (4 * sqrt(2.)),
                 (3 + sqrt(3.)) / (4 * sqrt(2.)),
                 (1 + sqrt(3.)) / (4 * sqrt(2.))
                 };

static float Daub2Synthesis [4] = {
                 (1 + sqrt(3.)) / (4 * sqrt(2.)),
                 (3 + sqrt(3.)) / (4 * sqrt(2.)),
                 (3 - sqrt(3.)) / (4 * sqrt(2.)), 
                 (1 - sqrt(3.)) / (4 * sqrt(2.))
                 };

static float Daub4Analysis [8] = { 
       -0.01059740178500,
       0.03288301166698,
       0.03084138183599,
       -0.18703481171888,
       -0.02798376941698,
       0.63088076792959,
       0.71484657055254,
       0.23037781330886
       };

static float Daub4Synthesis [8] = {
       0.23037781330886,
       0.71484657055254,
       0.63088076792959,
       -0.02798376941698,
       -0.18703481171888,
       0.03084138183599,
       0.03288301166698,
       -0.01059740178500
       };

static float   AntoniniAnalysis[10] =  {
                  0,  
	              3.782845550699535e-02,
			     -2.384946501937986e-02,
			     -1.106244044184226e-01,
			      3.774028556126536e-01,
			      8.526986790094022e-01,
			      3.774028556126537e-01,
			     -1.106244044184226e-01,
			     -2.384946501937986e-02,
			      3.782845550699535e-02
                  };

static float  AntoniniSynthesis [10] = { 
                   0,
	              -6.453888262893856e-02,
			      -4.068941760955867e-02,
			       4.180922732222124e-01,
			       7.884856164056651e-01,
			       4.180922732222124e-01,
			      -4.068941760955867e-02,
			      -6.453888262893856e-02,
			       0,
                   0
                  };
			      
// J. E. Odegard and C. S. Burrus, 
//  "Smooth biorthogonal wavelets for applications in image compression," 
// in Proceedings of DSP Workshop, Loen, Norway, September
// 1996 (http://www-dsp.rice.edu/publications). 
				
static float OdegardAnalysis[10] = {
   0,
   5.2865768532960523e-02,
  -3.3418473279346828e-02,
  -9.3069263703582719e-02,
   3.8697186387262039e-01,
   7.8751377152779212e-01,
   3.8697186387262039e-01,
  -9.3069263703582719e-02,
  -3.3418473279346828e-02,
   5.2865768532960523e-02
};

static float OdegardSynthesis[10] = {
   0,
  -8.6748316131711606e-02,
  -5.4836926902779436e-02,
   4.4030170672498536e-01,
   8.1678063499210640e-01,
   4.4030170672498536e-01,
  -5.4836926902779436e-02,
  -8.6748316131711606e-02,
   0,
   0
};

char *UserFilterFileName = NULL;

/***********************************************************************/

type_sb_filter get_filter_bank(char *UserArg)
{
   int c1;
   char *ch = new char[256];
   type_sb_filter FilRet = DEF_SB_FILTER;

   int N = sscanf(UserArg,"%d,%s",&c1,ch);
   // cout << "N = " << N << endl;
   if (N < 1)
   {
      fprintf(OUTMAN, "bad type of filter: %s\n", UserArg);
      exit(-1);
   }
   if (N > 0) 
   {
       FilRet = (type_sb_filter) c1;
       if ((c1 < 1) || (c1 > NBR_SB_FILTER))
       {
	   fprintf(OUTMAN, "bad type of filter: %s\n", UserArg);
	   exit(-1);
       }
       // cout << "Filter bank = " << StringSBFilter(FilRet) << endl;
       if (N > 1) 
       {
          UserFilterFileName = ch;
          // cout << "User file name = " << UserFilterFileName  << endl;
       }
   }
   return FilRet;
}

/***********************************************************************/

char* StringSBFilter (type_sb_filter type)
{
  switch (type)
    {
    case  F_HAAR: 
      return ("Haar filter"); 
    case F_BI2HAAR:
      return ("Biorthogonal 2/6 Haar filters"); 
    case F_BI4HAAR:
      return ("Biorthogonal 2/10 Haar filters");
    case  U_HAAR_B3S:
      return ("Haar B3-Spline filters I (for undecimated transform)"); 
    case  U_HAAR_B3S2:
      return ("Haar B3-Spline filters II (for undecimated transform)");       
    case  F_DAUBE_4: 
      return ("Daubechies filter 4 points (order 2)"); 
    case  F_DAUBE_8: 
      return ("Daubechies filter 8 points (order 4)"); 
    case  F_3_5:
      return ("3/5 filter"); 
    case  F_5_3: 
      return ("5/3 filter"); 
    case  F_MALLAT_7_9:
      return ("Biorthogonal 7/9 filters"); 
    case  F_MALLAT_9_7: 
      return ("Biorthogonal 9/7 filter");
    case  F_ODEGARD_9_7:
      return ("Odegard 9/7 filters"); 
    case  F_USER:
      return ("User's filters");
      
    default: 
      return ("Undefined sub-band filters");
    }
}
/***********************************************************************/

void sb_usage(type_sb_filter Filter)
{
    fprintf(OUTMAN, "         [-T type_of_filters]\n");
    for (int i = 1; i <= NBR_SB_FILTER; i++)
        fprintf(OUTMAN, "              %d: %s \n", i, StringSBFilter((type_sb_filter) i));
    fprintf(OUTMAN, "             default is %s\n\n", StringSBFilter(Filter));
     
    fprintf(OUTMAN, "         [-L]\n");
    fprintf(OUTMAN, "              Use a L2 normalization. Default is L1.\n");
}


/***********************************************************************/
 
// return the filename with '.wvf' as extension if the extension '.wvf' is not available
char *filtername(char *NameStep)
{
    char Sreturn[256];
    char Name[256];
    char *ValRet;

    strcpy(Name, NameStep);
    if (strstr(Name, ".wvf") != NULL)  strcpy(Sreturn, Name);
    else sprintf(Sreturn, "%s.%s", Name, "wvf");
    ValRet = strdup(Sreturn);
    return (ValRet);
}

/***********************************************************************/

FilterAnaSynt::FilterAnaSynt(type_sb_filter T_Filter) 
{
    reset_param();
    alloc(T_Filter);
}

/***********************************************************************/

FilterAnaSynt::FilterAnaSynt(char *FileName) 
{
    reset_param();
    FilterFileName=FileName; 
    alloc(F_USER);
}

/***********************************************************************/

FilterAnaSynt::~FilterAnaSynt() 
{
    reset_param();
}

/***********************************************************************/

void FilterAnaSynt::reset_param()
{
   Analysis = Synthesis = NULL; Analysis2 = Synthesis2 = NULL;
   Size_Ana  = Size_Synt = 0; Size_Ana2  = Size_Synt2 = 0;
   Start_Ana = Start_Synt =0; Start_Ana2 = Start_Synt2 =0;   
   TypeFilter = SB_UNKNOWN;
   Verbose = False;
   FilterFileName = NULL;
}

/***********************************************************************/

void  FilterAnaSynt::alloc(type_sb_filter T_Filter)
{
    TypeFilter = T_Filter;
    TypeNorm =  NORM_L2;
    switch(T_Filter)
    {
	 case F_HAAR:
	    Size_Ana = 2;
  	    Size_Synt = 2;
 	    Analysis = HaarAnalysis;
	    Synthesis = HaarAnalysis;
	    TypeNorm =  NORM_L1;
     break;	

  	 case F_BI2HAAR:
	    Size_Ana = 2;
  	    Size_Synt = 6;
 	    Analysis = HaarAnalysis;
	    Synthesis = Haar2Synthesis;
	    TypeNorm =  NORM_L1;
     break;	
   	 
     case F_BI4HAAR:
	    Size_Ana = 2;
  	    Size_Synt = 10;
 	    Analysis = HaarAnalysis;
	    Synthesis = Haar4Synthesis;
	    TypeNorm =  NORM_L1;
     break;

     case  U_HAAR_B3S :
        Size_Ana = Size_Ana2 = 2;
        Size_Synt = Size_Synt2 = 6;
        Analysis = U_HAAR_B3SAnalysis;
	    Synthesis = U_HAAR_B3SSynthesis;
 	    Analysis2 = U_HAAR_B3SAnalysis2;
	    Synthesis2 = U_HAAR_B3SSynthesis2;
	    TypeNorm =  NORM_L1;
     break;        
	  
     case  U_HAAR_B3S2 :
        Size_Ana = Size_Ana2 = 2;
        Size_Synt = Size_Synt2 = 6;
 	    Analysis = U_HAAR_B3S2Analysis;
	    Synthesis = U_HAAR_B3S2Synthesis;
 	    Analysis2 = U_HAAR_B3S2Analysis2;
	    Synthesis2 = U_HAAR_B3S2Synthesis2;
	    TypeNorm =  NORM_L1;
     break;        

	 case F_3_5:
	    Size_Ana = 3; 
	    Size_Synt = 5;
	    Analysis = F5_3Synthesis;
 	    Synthesis = F5_3Analysis;
	    TypeNorm =  NORM_L1;
     break;

	 case F_5_3:
	    Size_Ana = 5; 
  	    Size_Synt = 3;
 	    Analysis = F5_3Analysis;
	    Synthesis = F5_3Synthesis;
	    TypeNorm =  NORM_L1;
     break;       

	 case F_DAUBE_4:
 	    Size_Synt = 4;
	    Size_Ana = 4;
 	    Analysis = Daub2Analysis;
 	    Synthesis = Daub2Synthesis;
     break;
     
	 case F_DAUBE_8:
 	    Size_Synt = 8;
	    Size_Ana = 8;
 	    Analysis = Daub4Analysis;
 	    Synthesis = Daub4Synthesis;
     break;

	 case F_MALLAT_7_9:
       Size_Ana = 10;
       Size_Synt = 10;
	   Analysis = AntoniniSynthesis;  
	   Synthesis = AntoniniAnalysis;
 	   break;

     case F_MALLAT_9_7:
       Size_Ana = 10;
       Size_Synt = 10; 
	   Analysis = AntoniniAnalysis;
	   Synthesis = AntoniniSynthesis;
 	   break;

     case F_ODEGARD_9_7 :
        Size_Ana = 10;
        Size_Synt = 10; 
 	    Analysis = OdegardAnalysis;
	    Synthesis = OdegardSynthesis;
     break;        

	 case F_USER:
            if (FilterFileName == NULL)
            {
               if (UserFilterFileName != NULL) 
                     FilterFileName = UserFilterFileName;
               else
               {
                  FilterFileName = DEF_USER_FILTER_FILE_NAME;
                  FILE *FileDes = fopen(FilterFileName, "r");
                  if (FileDes != NULL) fclose(FileDes);
                  else
                  {
                      FilterFileName = (char *) getenv(USER_FILER_FILE_NAME);
                      if (FilterFileName == NULL)
                      { 
                         cout << "Error: the filter bank is not defined ... " << endl;
                         exit(-1);
                      }
                  }
               }
            }
            read_from_file(FilterFileName);
            break;
         default:
           cerr << "Error: unknown filter ... " << endl;
           exit(-1);
         break;
    }   
}

/***********************************************************************/

void  FilterAnaSynt::read_from_file(char *FileName)
{
   //  [Range low] [Range high]
   //  [Analysis LP filter coefficients]
   //  .
   //  .
   //  .
   //  [Range low] [Range high]
   //  [Synthesis LP filter coefficients]
   FILE *input=NULL;
   int i,ind,AnaLow, AnaHigh, SyntLow, SyntHigh;
   float Val;
   char *FName = filtername(FileName);
   
   input = fopen(FName,"r");
   if (input == NULL) 
   {
        cout << "Error: cannot open file " <<   FileName  << " ... or file doesn't exist" << endl;
        exit(-1);
   }
   if (fscanf(input,"%d %d",&AnaLow, &AnaHigh) != 2) 
   {
      cout << "Error: bad filter file format ... " << endl;
      exit(-1);
   }
   int nLine = AnaHigh - AnaLow + 1;

   if (Verbose == True)
   {
      cout << "Read filters from file " << FileName << endl;
      cout << "  Analysis: Range low = " << AnaLow << "  Range high = " << AnaHigh << " Size = " << nLine <<  endl;
   }
   int M = MAX(AnaHigh, ABS(AnaLow));
   Size_Ana = 2*M+1;
   Start_Ana = - Size_Ana/2;
   Analysis = new float [Size_Ana];
   double Sum=0.;
   double Sum2=0.;
   for (i=0; i < Size_Ana; i++) Analysis[i] = 0.;
   for (i=0; i < nLine; i++)
   {
      if (fscanf(input,"%f",&Val) != 1) 
      {
         cout << "Error: bad filter file format ... " << endl;
         exit(-1);
      }
      ind = Size_Ana /2+i+AnaLow;
      if ((ind < 0) || (ind >= Size_Ana))
      {
        cout << "Error: bad index ind = " << ind << endl;
        exit(-1);
      }
      Analysis[ind] = Val;
      Sum += Val;
      Sum2 += Val*Val;
      // cout << "Val = " << Val << endl;
   }
   if (Verbose == True) cout << "  Sum_i H[i]^2 = " << Sum2 << endl;
   
   if (ABS(Sum2-1.) > ABS(Sum-1.)) TypeNorm = NORM_L1;
   
   if (fscanf(input,"%d %d",& SyntLow, &SyntHigh) != 2) 
   {
      cout << "Error: bad filter file format ... " << endl;
      exit(-1);
   }
   nLine = SyntHigh - SyntLow + 1;

   if (Verbose == True)
   {
       cout << "  Synthesis: Range low = " <<  SyntLow << "  Range high = " <<  SyntHigh << "  Size = " <<   nLine << endl;
   }
   M = MAX(SyntHigh, ABS(SyntLow));
   Size_Synt = 2*M+1;
   Start_Synt = - Size_Synt/2;
   Synthesis = new float [Size_Synt];
   for (i=0; i <  Size_Synt; i++) Synthesis[i] = 0.;
   for (i=0; i < nLine; i++)
   {
      if (fscanf(input,"%f",&Val) != 1) 
      {
         cout << "Error: bad filter file format ... " << endl;
         exit(-1);
      }
      ind =  Size_Synt/2+i+SyntLow;
      if ((ind < 0) || (ind >=  Size_Synt))
      {
        cout << "Error: bad index ind = " << ind << endl;
        exit(-1);
      }
       Synthesis[ind] = Val;
   }
   if (Verbose == True) cout << "  Sum_i H1[i]^2 = " << Sum2 << endl;
   fclose(input);
}

/***********************************************************************/

/******************************************************************************
**                   Copyright (C) 1999 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Jean-Luc Starck
**
**    Date:  04/09/99
**    
**    File:  mr3d_trans.cc
**
*******************************************************************************
**
**    DESCRIPTION  multiresolution transform of cube
**    ----------- 
**                 
**    Usage: mr3d_trans options cube output
**
******************************************************************************/

#include "IM_Obj.h"
#include "IM_IO.h"
#include "SB_Filter.h"
#include "IM3D_IO.h"
#include "MR3D_Obj.h"
 
/****************************************************************************/

char Name_Cube_In[256];
char Name_Cube_Out[256];

int Nbr_Plan=4;
type_trans_3d Transform=TO3_MALLAT;

extern int  OptInd;
extern char *OptArg;

extern int  GetOpt(int argc, char *const*argv, char *opts);

sb_type_norm Norm = NORM_L1;
type_sb_filter SB_Filter = F_MALLAT_7_9;
type_lift LiftingTrans = DEF_LIFT;
Bool Verbose=False;
Bool InfoBand=False;
char FileSimu[180];
Bool WithFileSimu=False;

/*********************************************************************/

/*static int max_scale_number (int Nc)
{
    int Nmin = Nc;
    int ScaleMax;

    ScaleMax=iround(log((float)Nmin/(float)MAX_SIZE_LAST_SCALE) / log(2.)+ 1.);
    return (ScaleMax);
}*/

/*********************************************************************/
 
static void usage(char *argv[])
{
    int i;

    fprintf(OUTMAN, "Usage: %s options image output\n\n", argv[0]);
    fprintf(OUTMAN, "   where options =  \n");

    fprintf(OUTMAN, "\n");
    fprintf(OUTMAN, "         [-t type_of_multiresolution_transform]\n");
    for (i = 0; i < NBR_TRANS_3D; i++)
    fprintf(OUTMAN, "              %d: %s \n",i+1,
                                            StringTransf3D((type_trans_3d)i));
    fprintf(OUTMAN, "              default is %s\n", StringTransf3D((type_trans_3d) Transform));
    
    manline();
    fprintf(OUTMAN, "         [-n number_of_scales]\n");
    fprintf(OUTMAN, "              number of scales used in the multiresolution transform\n");
    manline();
    
    sb_usage(SB_Filter);
    manline();
    
    lifting_usage(LiftingTrans);
    manline();
    fprintf(OUTMAN, "         [-i]\n");
    fprintf(OUTMAN, "              print statistical information about each band.\n");
    manline();
    
    verbose_usage();    

    manline();
    manline();
    exit(-1);
}

/*********************************************************************/

/* GET COMMAND LINE ARGUMENTS */

static void transinit(int argc, char *argv[])
{
    int c;
    Bool OptL = False, Optf = False;
#ifdef LARGE_BUFF
    int VMSSize=-1;
    Bool OptZ = False;
    char VMSName[1024] = "";
#endif   

    /* get options */
    while ((c = GetOpt(argc,argv,"it:n:S:Ll:T:vzZ:")) != -1) 
    {
	switch (c) 
        {
	   case 'i': InfoBand = True; break;
	   case 'v': Verbose = True; break;
	   case 'L': Norm = NORM_L2; OptL = True; break;
	   case 'l':
 		if (sscanf(OptArg,"%d",&c ) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad lifting method: %s\n", OptArg);
	            exit(-1);
                    
		}
                if ((c > 0) && (c <= NBR_LIFT)) LiftingTrans = (type_lift) (c);
                else  
                {
		    fprintf(OUTMAN, "Error: bad lifting method: %s\n", OptArg);
	            exit(-1);
 		}
  		break;
	   case 'T':
 		if (sscanf(OptArg,"%d",&c ) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad filter: %s\n", OptArg);
	            exit(-1);
                    
		}
		Optf = True;
		SB_Filter = get_filter_bank(OptArg);
		// SB_Filter = (type_sb_filter) (c);
		break;
           case 't':
		/* -d <type> type of transform */
		if (sscanf(OptArg,"%d",&c ) != 1) 
                {
		    fprintf(OUTMAN, "bad type of multiresolution transform: %s\n", OptArg);
	            exit(-1);
                    
		}
                if ((c > 0) && (c <= NBR_TRANS_3D)) 
                                        Transform = (type_trans_3d) (c-1);
                else  
                {
		    fprintf(OUTMAN, "bad type of transform: %s\n", OptArg);
	            exit(-1);
 		}
		break;
	   case 'n':
		/* -n <Nbr_Plan> */
		if (sscanf(OptArg,"%d",&Nbr_Plan) != 1) 
                {
		    fprintf(OUTMAN, "bad number of scales: %s\n", OptArg);
		    exit(-1);
		}
                if ((Nbr_Plan <= 1) || (Nbr_Plan > MAX_SCALE_1D)) 
                 {
		    fprintf(OUTMAN, "bad number of scales: %s\n", OptArg);
		    fprintf(OUTMAN, "1 < Nbr Scales <= %d\n", MAX_SCALE_1D);
 		    exit(-1);
		}
		break;
	   case 'S':
		/* -n <Nbr_Plan> */
		if (sscanf(OptArg,"%s",FileSimu) != 1) 
                {
		    fprintf(OUTMAN, "bad file name: %s\n", OptArg);
		    exit(-1);
		}
		WithFileSimu = True;
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
 	   case '?':
			usage(argv);
		}
	}
      
 
        if ((Transform != TO3_MALLAT) && ((OptL == True) || (Optf == True)))
	{
	   fprintf(OUTMAN, "Error: option -f and -L are only valid with transforms 14 and 16  ... \n");
           exit(0);
	}
 	if ((Transform != TO3_LIFTING) && ( LiftingTrans != DEF_LIFT))
	{
	   fprintf(OUTMAN, "Error: option -f and -L are only valid with transforms 14 and 16  ... \n");
           exit(0);
	}
	
       /* get optional input file names from trailing 
          parameters and open files */
	if (OptInd < argc) strcpy(Name_Cube_In, argv[OptInd++]);
         else usage(argv);

	if (OptInd < argc) strcpy(Name_Cube_Out, argv[OptInd++]);
        else usage(argv);

	/* make sure there are not too many parameters */
	if (OptInd < argc)
        {
		fprintf(OUTMAN, "Error: too many parameters: %s ...\n", argv[OptInd]);
		exit(-1);
	}
#ifdef LARGE_BUFF
    if (OptZ == True) vms_init(VMSSize, VMSName, Verbose);
#endif
}

/*********************************************************************/
 
int main(int argc, char *argv[])
{
    fltarray Dat;
    int b;
    /* Get command line arguments, open input file(s) if necessary */
    lm_check(LIC_MR3);
    transinit(argc, argv);
 
    io_3d_read_data(Name_Cube_In, Dat);
 
    int Nx = Dat.nx();
    int Ny = Dat.ny();
    int Nz = Dat.nz();
    if (Verbose == True)
    {
       cout << "Filename in = " << Name_Cube_In << endl;
       cout << "Filename out = " << Name_Cube_Out  << endl;
       cout << "Nx = " << Dat.nx() << " Ny = " << Dat.ny() << " Nz = " << Dat.nz() << endl;
    }
    
    
    MR_3D MR_Data;
    MR_Data.Verbose=Verbose;
    FilterAnaSynt FAS;
    FilterAnaSynt *PtrFAS = NULL;
    if (Transform == TO3_MALLAT) 
    {
        FAS.Verbose = Verbose;
        FAS.alloc(SB_Filter);
        PtrFAS = &FAS;
    }
    MR_Data.alloc (Nx, Ny, Nz, Transform, Nbr_Plan, PtrFAS, Norm);
    
    if (Transform == TO3_LIFTING)  MR_Data.LiftingTrans = LiftingTrans;
        
    if (Verbose == True)
    {
       cout << "Transform = " << StringTransf3D((type_trans_3d) Transform) << endl;
      if (Transform == TO3_MALLAT)  
      {
         cout << "Filter = " <<  StringSBFilter(MR_Data.SBFilter) << endl;
         if (MR_Data.TypeNorm == NORM_L1) cout << "Norm L1 " << endl;
         else cout << "Norm L2 " << endl;
      }
    
      if (Transform == TO3_LIFTING)  
       cout << "Lifting method = " <<  StringLSTransform(MR_Data.LiftingTrans) << endl;
      cout << "Number of scales = " << MR_Data.nbr_scale() << endl;
    }
    if (Verbose == True) MR_Data.info_pos_band();
    
    if (Verbose == True) cout << "3D transform ... " << endl;
    MR_Data.transform (Dat);
    
    if (Verbose == True) cout << "Write ... " << endl;
    // io_3d_write_data(Name_Cube_Out, MR_Data.cube());
    MR_Data.write(Name_Cube_Out);
     
    if (InfoBand == True) 
       for (b=0; b < MR_Data.nbr_band(); b++) 
	   MR_Data.info_band(b);
	 
    exit(0);
} 


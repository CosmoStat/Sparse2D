/******************************************************************************
**                   Copyright (C) 2000 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Jean-Luc Starck
**
**    Date:  03/02/00
**    
**    File:  ridgelet.cc
**
*******************************************************************************
**
**    DESCRIPTION  ridgelet transform program
**    ----------- 
**                 
******************************************************************************/

#include "IM_Obj.h"
#include "IM_IO.h"
#include "Ridgelet.h"
 
char Name_Imag_In[256];  // input file image  
char Name_Imag_Out[256]; //  output file name  
 
extern int  OptInd;
extern char *OptArg;

extern int  GetOpt(int argc, char **argv, char *opts);

Bool Verbose=False; 
Bool Reverse=False; // Reverse transform
int Nbr_Plan = DEF_RID_NBR_SCALE;
int BlockSize=0;
Bool BlockOverlap = False;
type_ridgelet_WTtrans RidTrans = DEF_RID_TRANS;
Bool StatInfo = False;
Bool ExtractBand = False;
Bool ColTrans = False;
Bool DebugRid = False;
Bool MadAngleNorm= False;
int NbrIter = DEF_RADON_FSS_NBR_ITER;

/***************************************/

static void usage(char *argv[])
{
    // int i;
    fprintf(OUTMAN, "Usage: %s options in_image result\n\n", argv[0]);
    fprintf(OUTMAN, "   where options =  \n");

    fprintf(OUTMAN, "         [-t type_of_ridgelet]\n");
    for (int i = 0; i < NBR_RID_TRANS; i++)
    fprintf(OUTMAN, "              %d: %s \n",i+1, StringRidTransform((type_ridgelet_WTtrans) (i+1)));
    fprintf(OUTMAN, "              Default is %s.\n",  StringRidTransform(RidTrans));
    manline();

    fprintf(OUTMAN, "         [-n number_of_scales]\n");
    fprintf(OUTMAN, "             Number of scales used in the ridgelet transform.\n");
    fprintf(OUTMAN, "             Default is automatically calculated. \n");    
    manline();

    fprintf(OUTMAN, "         [-b BlockSize]\n");
    fprintf(OUTMAN, "             Block Size. Default is image size. \n");    
    manline();

    fprintf(OUTMAN, "         [-i]\n");
    fprintf(OUTMAN, "             Print statistical information about.\n");
    fprintf(OUTMAN, "             each band. Default is no. \n");    
    manline();

    fprintf(OUTMAN, "         [-O]\n");
    fprintf(OUTMAN, "             Block overlapping. Default is no. \n");    
    manline();

    fprintf(OUTMAN, "         [-r]\n");
    fprintf(OUTMAN, "             Inverse RIDGELET transform.\n"); 
    fprintf(OUTMAN, "             NbrIter is the number of iteration used in the reconstruction of\n");  
    fprintf(OUTMAN, "             the Slant Stack Radon transform. Default is %d\n", NbrIter); 
    manline();  

    fprintf(OUTMAN, "         [-R NbrIter]\n");
    fprintf(OUTMAN, "             NbrIter is the number of iteration used in the reconstruction of\n");  
    fprintf(OUTMAN, "             the Slant Stack Radon transform. Default is %d\n", NbrIter); 
    manline();  
    
    write_scales_x_band_usage();
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
#ifdef LARGE_BUFF
    int VMSSize=-1;
    Bool OptZ = False;
    char VMSName[1024] = "";
#endif     

    /* get options */
    while ((c = GetOpt(argc,argv,"R:MXCxt:Ob:n:irvzZ")) != -1) 
    {
	switch (c) 
        {
	   case 'M': MadAngleNorm= True; break;
	   case 'X': DebugRid = True; break;
           case 'C':  ColTrans =(ColTrans  == True) ? False: True; break; 
           case 'x': ExtractBand = True; break;
           case 't': 
               if (sscanf(OptArg,"%d",&c ) != 1) 
                {
                    fprintf(OUTMAN, "Error: bad type of ridgelet transform: %s\n", OptArg);
                    exit(-1);
                }
                if ((c > 0) && (c <= NBR_RID_TRANS)) RidTrans = (type_ridgelet_WTtrans) (c);
                else  
                {
                   fprintf(OUTMAN, "Error: bad type of ridgelet transform: %s\n", OptArg);
                   exit(-1);
                }
                break;
           case 'i': StatInfo= (StatInfo ==True) ? False: True;break;
           case 'O': BlockOverlap= (BlockOverlap==True) ? False: True;break;
           case 'r': Reverse=True;break;
           case 'v': Verbose = True;break;
	   case 'b':
                if (sscanf(OptArg,"%d",&BlockSize) != 1) 
                {
                    fprintf(OUTMAN, "bad block size: %s\n", OptArg);
                    exit(-1);
                }
                break;
	   case 'R':
                /* -n <Nbr_Plan> */
                if (sscanf(OptArg,"%d",&NbrIter) != 1) 
                {
                    fprintf(OUTMAN, "bad number of iterations: %s\n", OptArg);
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
      
       if (Reverse == True)
       {
           if (StatInfo == True) 
           {
              fprintf(OUTMAN, "Error: i and r option are not compatible...\n");
              exit(-1);
           }
           if (ExtractBand == True) 
           {
              fprintf(OUTMAN, "Error: x and r option are not compatible...\n");
              exit(-1);
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
		fprintf(OUTMAN, "Too many parameters: %s ...\n", argv[OptInd]);
		exit(-1);
	}

       	   
#ifdef LARGE_BUFF
    if (OptZ == True) vms_init(VMSSize, VMSName, Verbose);
#endif  
}

/***************************************************************************/
 
int main(int argc, char *argv[])
{
    int s,k;
    fitsstruct Header;
    char Cmd[512];
    Cmd[0] = '\0';
    for (k =0; k < argc; k++) sprintf(Cmd, "%s %s", Cmd, argv[k]);
     
    // Get command line arguments, open input file(s) if necessary
    lm_check(LIC_MR4);
    filtinit(argc, argv);

    if ((Verbose == True) && (Reverse == False))
    { 
        cout << endl << endl << "PARAMETERS: " << endl << endl;
        cout << "File Name in = " << Name_Imag_In << endl;
        cout << "File Name Out = " << Name_Imag_Out << endl;  
        cout << "Ridgelet transform = " <<  StringRidTransform(RidTrans) << endl;  
        if (BlockOverlap == False) cout << "No overlapping " << endl;
	else cout << "Overlapping " << endl;
        if (Reverse == True) cout << "Inverse Rigelet transform  " <<  endl;  
        if (BlockSize > 0)  cout << "BlockSize = " << BlockSize <<endl;
    }

    Ifloat Data;
    Ifloat Result;
    if (Reverse == False) 
         io_read_ima_float(Name_Imag_In, Data, &Header);
    Header.origin = Cmd;
    if (RidTrans ==  RID_FINITE) 
    {    
        CPRIME_NUMBER CPN;

        if ((BlockSize > 0) && ( CPN.is_prime_number((unsigned long) BlockSize) == False))
        {
           cout << "Error: Block size must be a prime number ... " << endl;
           cout << "       Previous prime number = " << CPN.previous_prime_number((unsigned long) BlockSize) << endl;
           cout << "       Next prime number = " << CPN.next_prime_number((unsigned long) BlockSize) << endl;
           exit(-1);
        }
    }

    if (BlockSize > Data.nc())
    {
        cout << "Warning: Block size must lower than the image size ... " << endl;
        cout << "         Block size is set to image size " << endl;
        BlockSize = 0;
    }

    FilterAnaSynt SelectFilter(F_MALLAT_7_9);
    SubBandFilter SB1D(SelectFilter, NORM_L2);
    Ridgelet RG(SB1D);
    if (RidTrans == RID_UWT)
    {
        // UndecSubBandFilter *Ptr_Undec_SB1D
        // UndecSubBandFilter USF(*Ptr_SB1D);
    }


    RG.RidTrans = RidTrans;
    if (Nbr_Plan != DEF_RID_NBR_SCALE) 
    {
       RG.NbrScale = Nbr_Plan;
       RG.GetAutoNbScale = False;
       if (Nbr_Plan < 2)
       {
          RG.BlockOverlap = False;
       }
       else RG.BlockOverlap = BlockOverlap;
    }         
    else RG.BlockOverlap = BlockOverlap;

    RG.NbrScale = Nbr_Plan;
    RG.Verbose = Verbose;
    RG.StatInfo = StatInfo;
    RG.ColTrans = ColTrans;
    if (Reverse == False) 
    {
        RG.transform(Data,Result,BlockSize);
        if (MadAngleNorm == True) RG.mad_normalize_per_angle(Result);

        if (ExtractBand == True) 
        {
           char Prefix[256];
           char Name_Imag[256];
           io_strcpy_prefix(Prefix,  Name_Imag_Out);
           Ifloat ImaScale(Result.nl(), Result.nc(), "scale");
           for (s = 0; s < RG.NbrScale; s++)
           {
              sprintf (Name_Imag, "%s_band_%d", Prefix, s+1);
              RG.get_scale(Result,ImaScale,s);
              io_write_ima_float(Name_Imag,ImaScale);
           } 
        }
        // if (RG.VarianceStab == False) cout << "No VarianceStab " << endl;
        // else cout << "VarianceStab " << endl;
	RG.write(Name_Imag_Out, Result);
	
	// cout << "SIMU" << endl;
	// RG.set_tab_norm_from_simu(Data, RG.NbrScale);
	
       // RG.normalize_coef(Result, False, (float) (sqrt((double) (Data.nc()))));
       // RG.init_scale(Result, RG.NbrScale-1);  
       
       if (DebugRid == True)
       { 
         cout << endl << "RECONSTRUCTION " << endl;
         Ifloat Rec(Data.nl(), Data.nc());
         RG.RADON.NbrIter = NbrIter; 
	     RG.RADON.SSR.Verbose= True;
         RG.recons(Result,Rec); 
           cout << "Ridgelet transform Rec = " <<  StringRidTransform(RidTrans) << endl;  
         
           cout << "Inverse Rigelet transform  " <<  endl;  
           cout << "BlockSize = " << RG.block_size()  <<endl;
           cout << "NbrScale = " << RG.NbrScale <<endl;
           if (RG.BlockOverlap == False) cout << "No overlapping " << endl;
           else cout << "Overlapping " << endl;
           if (RG.VarianceStab == False) cout << "No VarianceStab " << endl;
           else cout << "VarianceStab " << endl;
           
         io_write_ima_float("xx_rec.fits",Rec); 
         INFO_X(Rec, "REC");
         Data -= Rec;
         INFO_X(Data, "RESI"); 
         io_write_ima_float("xx_resi.fits",Data);
       }
    }
    else 
    {
       RG.read(Name_Imag_In, Data);
       if (Verbose == True)   
       { 
           cout << endl << endl << "PARAMETERS: " << endl << endl;
           cout << "File Name in = " << Name_Imag_In << endl;
           cout << "File Name Out = " << Name_Imag_Out << endl;  
           cout << "Ridgelet transform = " <<  StringRidTransform(RG.RidTrans) << endl;  
           if (RG.BlockOverlap == False) cout << "No overlapping " << endl;
	else cout << "Overlapping " << endl;
        if (RG.block_size() > 0)  cout << "BlockSize = " << RG.block_size() <<endl;
       }       
       RG.RADON.NbrIter = NbrIter;
       // RG.RADON.SSR.Verbose= True;
       RG.recons(Data,Result);
       Header.width = Result.nc();
       Header.height = Result.nl();
       Header.npix = Header.width*Header.height;
       Header.bitpix = -32;
       type_format FormatData = io_which_format(Name_Imag_Out);
        if ((FormatData == F_UNKNOWN) && (RG.FormatInputImag != F_UNKNOWN))  
                                           io_set_format(RG.FormatInputImag);
       io_write_ima_float(Name_Imag_Out, Result, &Header);
    }

    exit(0);
}

/******************************************************************************
**                   Copyright (C) 2000 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Philippe
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
#include "Linelet3D.h"
#include "macro.h"
 
char Name_Dat_In[256];  // input file image  
char Name_Dat_Out[256]; //  output file name  
 
extern int  OptInd;
extern char *OptArg;

extern int  GetOpt(int argc, char **argv, char *opts);

Bool Verbose=False;         // verbose mode
Bool TrDebug=False;
Bool Control=False;
Bool WriteFile=False;
Bool StatInfo = False;      // StatInfo mode
Bool TakeOnlyCompletePlane=False; 
                           // if true take only planes with no zero
			   // should normaly not be used
Bool DoNotCorrectSigma=False;
                           // Multiply each plane in order to have the
			   // the noise = sqrt(N)
Bool DoNotWriteTransf=False;// do not write transform file
Bool Reverse=False;         // Reverse transform
Bool BlockOverlap=False;    // Block overlap
Bool TestingPurpose=False;  // DO not use it.
Bool Normalize=False;       // normalize input data
int BlockSize=0;            // default size block
char TraceFileName[80];
Bool WithTrace=False;           // no file trace
char FileSimu[180];
Bool WithFileSimu=False;
Bool ExtractBand = False;

/***************************************/

static void usage(char *argv[]) {
    
   manline();   
    
   fprintf(OUTMAN, "         [-b]\n");
   fprintf(OUTMAN, "             Block Size.\n");
   fprintf(OUTMAN, "             Default Cube Size.\n");
   manline();   
   
   fprintf(OUTMAN, "         [-O]\n");
   fprintf(OUTMAN, "             Block overlapping.\n");
   fprintf(OUTMAN, "             Default is no.\n");
   manline();

    
   fprintf(OUTMAN, "         [-r]\n");
   fprintf(OUTMAN, "             Inverse Beamlet3d transform.\n");  
   manline();    
   
   write_scales_x_band_usage();
   manline();
        
   vm_usage();
   manline();  
      
   verbose_usage();
   manline();         
    
  exit(-1);
}
  
/*********************************************************************/

/* GET COMMAND LINE ARGUMENTS */
static void filtinit(int argc, char *argv[]) {
   
   int c;
   int NewBlockSize=False;
#ifdef LARGE_BUFF
   int VMSSize=-1;
   Bool OptZ = False;
   char VMSName[1024] = "";
#endif     

   /* get options */
   while ((c = GetOpt(argc,argv,"xb:T:lunatOricwvdzZ:")) != -1) {
      
      switch (c) {
      case 'x': ExtractBand = True; break;
      case 'i': StatInfo= (StatInfo ==True) ? False: True;break;
      case 'v': Verbose = True;break;
      case 'r': Reverse = True;break;
      case 'c': Control = True; break; 
      case 'w': WriteFile = True; break; 
      case 'O':	BlockOverlap = True; break; 
      case 'd':	TrDebug = True; break; 
      case 't': TestingPurpose = True; break; 
      case 'n': Normalize = True; break;
      case 'u': TakeOnlyCompletePlane = True; break;
      case 'l': DoNotCorrectSigma = True; break;
      case 'a': DoNotWriteTransf = True; break;
      
      case 'b':
         if (sscanf(OptArg,"%d",&BlockSize) != 1) {
            fprintf(OUTMAN, "bad block size: %s\n", OptArg);
            exit(-1);
         }
	 NewBlockSize=True;
         break;
	 
      case 'T': /* abaque file */
         if (sscanf(OptArg,"%s", TraceFileName) != 1) {
            fprintf(stderr, "Error: bad trace file name: %s\n", OptArg); 
	    exit(-1);
         }
	 WithTrace = True;
	 break; 	 
#ifdef LARGE_BUFF
      case 'z':
         if (OptZ == True) {
            fprintf(OUTMAN, "Error: Z option already set...\n");
            exit(-1);
         }
	 OptZ = True;
	 break;
      case 'Z':
	 if (sscanf(OptArg,"%d:%s",&VMSSize, VMSName) < 1) {
	    fprintf(OUTMAN, "Error: syntaxe is Size:Directory ... \n");
	    exit(-1);
         }
	 if (OptZ == True) {
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
      
   /* get otrace_paramptional input file names from trailing parameters and open files */
   if (OptInd < argc) strcpy(Name_Dat_In, argv[OptInd++]);
   else usage(argv);

   if (OptInd < argc) strcpy(Name_Dat_Out, argv[OptInd++]);
   else usage(argv);

   /* make sure there are not too many parameters */
   if (OptInd < argc) {
      fprintf(OUTMAN, "Too many parameters: %s ...\n", argv[OptInd]);
      exit(-1);
   }
   
//    if (BlockOverlap && NewBlockSize) {
//       fprintf(OUTMAN, 
//         "BlockOverlap option (-O) and new block size (-b xx) are not allowed \n");
//       exit(-1);   
//    }

       	   
#ifdef LARGE_BUFF
   if (OptZ == True) vms_init(VMSSize, VMSName, Verbose);
#endif  
}
 
/***************************************************************************/
 
int main(int argc, char *argv[]) 
{
   
   // init local var
   //---------------
   fitsstruct Header;
   char Cmd[512];
   Cmd[0] = '\0';
    extern softinfo Soft;
    Soft.mr4();
    for (int k =0; k < argc; k++) sprintf(Cmd, "%s %s", Cmd, argv[k]);
     
   // Get command line arguments, open input file(s) if necessary
   lm_check(LIC_MR4);
   filtinit(argc, argv);
   
   if (Verbose) TRACEUR.TraceVerbose();
   if (TrDebug) TRACEUR.TraceDebug();
   TRACEUR.SetScreenFlux(true);
   if (WithTrace) TRACEUR.SetFileFlux(true, TraceFileName);
   TRACEUR.SetScreenFlux(true);
      
   char Msg[40];
   if (Verbose) 
      TRACEUR << new Trace (Trace::METIER, "Main", "info param main",
                            " -- [-v]  : Verbose" , NULL, Trace::VERBOSE);
   if (TrDebug)  	 
      TRACEUR << new Trace (Trace::METIER, "Main", "info param main",
                            " -- [-d]  : Debug trace" , NULL, Trace::VERBOSE);
   if (BlockSize != 0) {
      sprintf (Msg, " -- [-b:] : Block size = %d", BlockSize);
      TRACEUR << new Trace (Trace::METIER, "Main", "info param main",
                            Msg, NULL, Trace::VERBOSE);   
   }
   if (StatInfo)					    
      TRACEUR << new Trace (Trace::METIER, "Main", "info param main",
                            " -- [-i]  : Stat info", NULL, Trace::VERBOSE); 
   if (Reverse) 
      TRACEUR << new Trace (Trace::METIER, "Main", "info param main",
                            " -- [-r]  : Reverse transform", NULL, Trace::VERBOSE);        			    
   if (Control) 
      TRACEUR << new Trace (Trace::METIER, "Main", "info param main",
                            " -- [-c]  : Control data", NULL, Trace::VERBOSE);        	          
   if (WriteFile)
      TRACEUR << new Trace (Trace::METIER, "Main", "info param main",
                            " -- [-w]  : Write all file", NULL, Trace::VERBOSE);          
   if (BlockOverlap) 
      TRACEUR << new Trace (Trace::METIER, "Main", "info param main",
                            " -- [-O]  : Block overlap", NULL, Trace::VERBOSE);      
   if (TestingPurpose) 
      TRACEUR << new Trace (Trace::METIER, "Main", "info param main",
                            " -- [-t]  : Call reconstruction", NULL, Trace::VERBOSE);    
   if (Normalize)  
      TRACEUR << new Trace (Trace::METIER, "Main", "info param main",
                            " -- [-n]  : Normalize data in", NULL, Trace::VERBOSE);        
   if (TakeOnlyCompletePlane)  
      TRACEUR << new Trace (Trace::METIER, "Main", "info param main",
                            " -- [-u]  : Take only complete plane", NULL, Trace::VERBOSE);    
   if (DoNotCorrectSigma) 
      TRACEUR << new Trace (Trace::METIER, "Main", "info param main",
                            " -- [-l]  : Do not correct sigma value", NULL, Trace::VERBOSE);
   sprintf (Msg, "File Name in = %s", Name_Dat_In);			    
   TRACEUR << new Trace (Trace::METIER, "Main", "info param main",
                         Msg, NULL, Trace::VERBOSE);
   sprintf (Msg, "File Name out = %s", Name_Dat_Out);				        
   TRACEUR << new Trace (Trace::METIER, "Main", "info param main",
                         Msg, NULL, Trace::VERBOSE);			    
   
   
   // param control
   //--------------
   //if ((BlockSize > 0) && (is_power_of_2(BlockSize) == False)) {
   //   cout << "Error: Block size must be a power of two ... " << endl;
   //   exit(-1);
   //}
       
   // init Linelet classes
   //---------------------
   Linelet3D LL;    
    
   // direct linelet transform
   //-------------------------
   if (!Reverse){
   
      fltarray Data;
      fltarray Result;
      io_3d_read_data(Name_Dat_In, Data);
      int Nx = Data.nx();
      int Ny = Data.ny(); 
      int Nz = Data.nz();
      
      double Mean = Data.mean();
      double Sigma = Data.sigma();

      if (Normalize)
         for (int i=0;i<Nx;i++)
	 for (int j=0;j<Ny;j++)
	 for (int k=0;k<Nz;k++)
            Data(i,j,k) = (Data(i,j,k)-Mean)/Sigma;
       
      if (BlockSize > Data.nx()+1) {
         cout << "Warning: Block size must lower than the cube size ... " << endl;
         cout << "         Block size is set to cube size " << endl;
         BlockSize = 0;
      }
    
      FilterAnaSynt SelectFilter(F_MALLAT_7_9);
      SubBandFilter SB1D(SelectFilter, NORM_L2);
      
      LL.set_Verbose(Verbose);
      LL.set_WriteFile(WriteFile);
      LL.set_StatInfo(StatInfo);
      LL.set_TakeOnlyCompletePlane(TakeOnlyCompletePlane);
      LL.set_DoNotCorrectSigma(DoNotCorrectSigma);
      LL.set_LinTransf(LIN3D_OWT, &SB1D);
      LL.set_BlockOverlap(BlockOverlap);
      if (Verbose) LL.trace_param();
      LL.transform (Data, Result, BlockSize);
      if (StatInfo) LL.info_stat(Result, FileSimu, WithFileSimu);

      if (!DoNotWriteTransf)
         LL.write(Name_Dat_Out, Result);

      if (ExtractBand == True)
      {
         char Prefix[256];
         char Name_Imag[256]; 
         fltarray t;   
	 io_strcpy_prefix(Prefix,  Name_Dat_Out);
         
	 for (int i=0;i< LL.nbr_band();i++) 
	 {
 	    LL.get_band (Result,t,i);
 	    sprintf (Name_Imag, "%s_band_%d", Prefix, i+1);
            io_3d_write_data(Name_Imag, t);
	 }
     }
      
      // !!!!!!!!!!! PROV testing purpose !!!!!!!!!!!!
      if (TestingPurpose) {
         cout << "!!!!!!!!! PROV recons for testing purpose !!!!!!!!!" << endl;
	          	 
	 fltarray Recons (Nx,Ny,Nz);
         LL.recons (Result, Recons);

         if (Normalize)
            for (int i=0;i<Nx;i++)
	    for (int j=0;j<Ny;j++)
	    for (int k=0;k<Nz;k++)
               Recons(i,j,k) = Recons(i,j,k)*Sigma+Mean;	 
	 
         io_3d_write_data(Name_Dat_Out, Recons); 
      }
          
   // reverse linelet transform
   //--------------------------
   } else {
   
      fltarray Result;
      LL.mr_io_read  (Name_Dat_In, Result);      
      
      FilterAnaSynt SelectFilter(F_MALLAT_7_9);
      SubBandFilter SB1D(SelectFilter, NORM_L2);
      
      LL.set_Verbose(Verbose);
      LL.set_WriteFile(WriteFile);
      LL.set_StatInfo(StatInfo);
      LL.set_TakeOnlyCompletePlane(TakeOnlyCompletePlane);
      LL.set_DoNotCorrectSigma(DoNotCorrectSigma);
      LL.set_LinTransf(LIN3D_OWT, &SB1D);
      
      if (Verbose) LL.trace_param();
            
      //int Nx = Result.nx();
      //int Ny = Result.ny(); 
      //int Nz = Result.nz();
      fltarray Recons (LL.get_SizeX(),LL.get_SizeY(),LL.get_SizeZ());
          
      LL.recons (Result, Recons);
      io_3d_write_data(Name_Dat_Out, Recons);
   }
   
   TRACEUR.FermerFichierDeTraces();
   exit(0);
}

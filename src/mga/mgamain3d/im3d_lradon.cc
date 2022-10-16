/******************************************************************************
**                   Copyright (C) 1999 by CEA
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
**    File:  im3d_lradon.cc
**
*******************************************************************************
**
**    DESCRIPTION  Radon program
**    ----------- 
**                 
******************************************************************************/

#include "IM_Obj.h"
#include "IM_IO.h"
#include "IM3D_PartialRadon.h"
#include "IM3D_IO.h"
#include "macro.h"

char Name_Cube_In[256];  // input file image  
char Name_Imag_Out[256]; // output file name  
 
extern int  OptInd;
extern char *OptArg;

extern int  GetOpt(int argc, char **argv, char *opts);

Bool Verbose=False;
Bool Reverse=False; 
Bool Control=False;
Bool WriteFile=False;
Bool WriteCoord=False;
Bool TrDebug=False;
Bool StatInfo=False;
char TraceFileName[80];
Bool WithTrace=False;          

/******************************************************************************/
static void usage(char *argv[]) {
    
    fprintf(OUTMAN, "Usage: %s options in_image out_image\n\n", argv[0]);
    fprintf(OUTMAN, "   where options =  \n");
    manline();

    fprintf(OUTMAN, "         [-r]\n");
    fprintf(OUTMAN, "             Inverse RADON transform.\n");  
    manline();  

    vm_usage();
    manline(); 
       
    verbose_usage();
    manline();         
    manline();
    exit(-1);
}
  
/******************************************************************************/

/* GET COMMAND LINE ARGUMENTS */
static void filtinit(int argc, char *argv[]) {
   
   int c;
#ifdef LARGE_BUFF
   int VMSSize=-1;
   Bool OptZ = False;
   char VMSName[1024] = "";
#endif     

   /* get options */
   while ((c = GetOpt(argc,argv,"T:sdcwirvzZ")) != -1) {
	
      switch (c) { 
      case 'r': Reverse=True; break;
      case 'v': Verbose = True; break;
      case 'c': Control = True; break; 
      case 'w': WriteFile = True; break; 
      case 'i': WriteCoord = True; break; 
      case 'd': TrDebug = True; break; 
      case 's': StatInfo = True; break;  
      
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
       
   // get optional input file names from trailing parameters and open files 
   if (OptInd < argc) strcpy(Name_Cube_In, argv[OptInd++]);
   else usage(argv);

   if (OptInd < argc) strcpy(Name_Imag_Out, argv[OptInd++]);
   else usage(argv);

   // make sure there are not too many parameters */
   if (OptInd < argc) {
      fprintf(OUTMAN, "Too many parameters: %s ...\n", argv[OptInd]);
      exit(-1);
   }
 	   
#ifdef LARGE_BUFF
   if (OptZ == True) vms_init(VMSSize, VMSName, Verbose);
#endif  
}


/******************************************************************************/

int main(int argc, char *argv[]) 
{

   // local var
   //----------
   int k;
   fitsstruct Header;
   char Cmd[512];
   Cmd[0] = '\0';
   for (k =0; k < argc; k++) sprintf(Cmd, "%s %s", Cmd, argv[k]);
     
   // Get command line arguments, open input file(s) if necessary
   lm_check(LIC_MR4);
   filtinit(argc, argv);
   
   if (Verbose) TRACEUR.TraceVerbose();
   if (TrDebug) TRACEUR.TraceDebug();
   TRACEUR.SetScreenFlux(true);
   if (WithTrace) TRACEUR.SetFileFlux(true, TraceFileName);   
   
   if (Verbose == True) { 
      cout << " im3d_lradon parameters: " << endl;
      if (Reverse) cout << " -- [-r]  : Reverse transform" << endl;
      if (Verbose) cout << " -- [-v]  : Verbose" << endl;
      cout << " -- File Name in = " << Name_Cube_In << endl;
      cout << " -- File Name Out = " << Name_Imag_Out << endl;  
   }
 
   PartialRadon3D PRD;
   PRD.setVerbose(Verbose);
   PRD.setControl(Control);
   PRD.setWriteFile(WriteFile);
   PRD.setWriteCoord(WriteCoord); 
   PRD.setStatInfo(StatInfo);   
     
   if (!Reverse) {
       
      fltarray Data;
      io_3d_read_data (Name_Cube_In, Data);
      int Nx = Data.nx();
      int Ny = Data.ny();
      int Nz = Data.nz();
     
      
      if ((Nx != Ny) || (Nx != Nz)) {
         cout << "Error: Cube size must be equal ... " << endl;
         exit(-1);
      }   
      
      int SizeCube = Nx;  
      
      if (Nx%2 != 1) {
         cout << "Error: Cube size must be 2n+1 ... " << endl;
         exit(-1);      
      }
	 
      PRD.setParam(SizeCube,0);
      PRD.alloc();
        
      fltarray Result (SizeCube, SizeCube, PRD.getNbPlane());     
      PRD.transform (Data, Result);
      
      io_3d_write_data (Name_Imag_Out, Result); 
  
   } else {
   
      fltarray LRadon;
      io_3d_read_data (Name_Cube_In, LRadon);
      int Nx = LRadon.nx();
      int Ny = LRadon.ny();
      int Nz = LRadon.nz();
            
      int SizeCube = LRadon.ny();
      PRD.setParam(SizeCube,0);
      
      if (Nx!=Ny || Nz !=PRD.getNbPlane()) {
         cout << "Bad partial radon transform ... " << endl;
         exit(-1);      
      }     
      
      fltarray Recons (SizeCube, SizeCube,SizeCube);
            
      PRD.recons(LRadon, Recons);
      io_3d_write_data(Name_Imag_Out, Recons);
   }
   
   TRACEUR.FermerFichierDeTraces();
   
   exit(0);
}

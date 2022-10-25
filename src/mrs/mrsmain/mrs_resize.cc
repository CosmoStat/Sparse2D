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
**    Date:  28/11/08
**    
**    File:  mrs_resize.cc
**
*******************************************************************************
**
**    DESCRIPTION:  resize an image 
**    -----------
**
******************************************************************************/

#include"HealpixClass.h"

char Name_Imag_In[512]; /* input file image */
char Name_Imag_Out[512]; /* output file name */
char Filter_FileName [512];
 
extern int  OptInd;
extern char *OptArg;
extern int  GetOpt(int argc, char **argv, char *opts);

Bool Verbose=False;
int Nside = 0.;


/***************************************************************************/
 
static void usage(char *argv[])
{
    fprintf(OUTMAN, "Usage: %s options in_map out_map \n\n", argv[0]);
    fprintf(OUTMAN, "   where options =  \n");    
    fprintf(OUTMAN, "         [-n Nside]\n");
    fprintf(OUTMAN, "             Nside out the resized image.\n");    
    fprintf(OUTMAN, "         [-v]\n");
    fprintf(OUTMAN, "             Verbose. Default is no.\n");
    exit(-1);
}

/*********************************************************************/

/* GET COMMAND LINE ARGUMENTS */
static void sinit(int argc, char *argv[])
{
    int c; 
    /* get options */
    while ((c = GetOpt(argc,argv,(char *) "n:vzZ")) != -1) 
    {
	switch (c) 
        { 
 	           case 'n':
 	     	      if (sscanf(OptArg,"%d",&Nside) != 1) 
                  {
		            fprintf(OUTMAN, "Error: bad nside  value: %s\n", OptArg);
		            exit(-1);
	            	}
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

	if (OptInd < argc) strcpy(Name_Imag_Out, argv[OptInd++]);
         else usage(argv);

	/* make sure there are not too many parameters */
	if (OptInd < argc)
        {
		fprintf(OUTMAN, "Too many parameters: %s ...\n", argv[OptInd]);
		exit(-1);
	}
}
  
/*********************************************************************/

int main(int argc, char *argv[])
{
    int k;
    double Min, Max;
    fitsstruct Header, HD1;
    char Cmd[512];
    Cmd[0] = '\0';
    for (k =0; k < argc; k++) sprintf(Cmd, "%s %s", Cmd, argv[k]);
     
    /* Get command line arguments, open input file(s) if necessary */
    sinit(argc, argv);
    if (Verbose == True)
    { 
        cout << "# PARAMETERS: " << endl;
        cout << "# File Name in = " << Name_Imag_In << endl;
        cout << "# File Name Out = " << Name_Imag_Out << endl;   
        if (Nside > 0)  cout << "# Nside = " <<  Nside << endl;
    }
   
   Hdmap Map;
   Map.read(Name_Imag_In);
   if (Verbose == True)
   {
   	   Map.minmax(Min,Max);
   	   if (Verbose == True) Map.info((char*) "Input MAP");
   }

   int NsideIn = Map.Nside();
   if (Nside <=0) Nside = NsideIn;
 
   Hdmap outmap;
   outmap.alloc (Nside);
   outmap.Import(Map);

  //    {
  //    Healpix_Base::SetNside(nside, scheme);
  //    map.alloc(npix_);
  //    }(Map);
  
   if (Verbose == True)
   {
   	   if (Verbose == True) outmap.info((char*) "Output MAP");
   }

   outmap.write(Name_Imag_Out);
   exit(0);
}


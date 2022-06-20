/******************************************************************************
 **                   Copyright (C) 2010 by CEA  
 *******************************************************************************
 **
 **    UNIT
 **
 **    Version: 1.0
 **
 **    Author:  Jean-Luc Starck  
 **
 **    Date:  25/01/2010
 **    
 **    File:  mrs_matmask.cc 
 **
 *******************************************************************************
 **
 **    DESCRIPTION: create the matrix relative to a mask
 **    -------------------
 **
 ******************************************************************************/


#include"HealpixClass.h"
#include"MatMask.h"

extern int  OptInd;
extern char *OptArg;
extern int  GetOpt(int argc, char **argv, char *opts);

char Name_Mask_file[512]; /* input mask file */
char Name_MatMask_file[512]; /* output MatMask file */
int lmax = 0;
Bool Verbose=False;
Bool Reverse=False;

//--------------------------------------------------------------------------
 
static void usage(char *argv[])
{
    fprintf(OUTMAN, "Usage: %s options in_map in_mask out_map  [out_powspec]  \n\n", argv[0]);
    fprintf(OUTMAN, "   where options =  \n");

    fprintf(OUTMAN, "         [-l lmax]\n");
    fprintf(OUTMAN, "             Default is MIN(%d, 3*nside). \n", ALM_MAX_L);
    fprintf(OUTMAN, "         [-R]\n");
    fprintf(OUTMAN, "             Reverse the mask (i.e. Mask = 1 - Mask). \n");

    exit(-1);
}
 
/***************************************************************************/

static void sinit(int argc, char *argv[])
{
   int c;
   
    /* get options */
    while ((c = GetOpt(argc,argv,(char *) "Rl:vzZ")) != -1) 
    {
	switch (c) 
        {   
	        case 'R': Reverse=True;	break;        
            case 'l':
 	     	      if (sscanf(OptArg,"%d",&lmax) != 1) 
                  {
		            fprintf(OUTMAN, "Error: bad lmax  value: %s\n", OptArg);
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
    if (OptInd < argc) strcpy(Name_Mask_file, argv[OptInd++]);
         else usage(argv);

    if (OptInd < argc) strcpy(Name_MatMask_file, argv[OptInd++]);
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
	//Call "mrs_matmask file_in file_out"
    sinit(argc, argv);
    if (Verbose == True)
    { 
        cout << "# PARAMETERS: " << endl;
        cout << "# File Name in = " << Name_Mask_file << endl;
        cout << "# InpMap  File Name Out = " << Name_MatMask_file << endl;   
        if (lmax > 0) cout << "# Lmax = " << lmax << endl;
    }
	//cout << "Name_Mask_file" << Name_Mask_file << '\n';
	//cout << "Name_MatMask_file" << Name_MatMask_file << '\n';
	
    Hdmap Mask;
    
    Mask.read( Name_Mask_file );
    if (Reverse == True)
    for (int i=0; i < Mask.Npix(); i++) Mask[i] = 1. - Mask[i];
    
    //cout << "Mask read" << '\n';
    int Nside = Mask.nside();
    lmax = mrs_get_lmax (lmax,  Nside,   0.0);

	dblarray MatMask = mrs_matmask( Mask, lmax );

    if (Verbose == True)
    { 
        cout << " Matrix: Nl = " << MatMask.ny() << ", Nc = " << MatMask.nx() << endl;
        MatMask.info(" ");
    }
  //  fltarray Mat;
  //  Mat.alloc ( MatMask.nx(), MatMask.ny());
  //  for (int i=0; i < Mat.nx(); i++)
  //  for (int j=0; j < Mat.ny(); j++) Mat(i,j) = (float) MatMask(i,j);
  //  fits_write_fltarr( Name_MatMask_file, Mat);
    
 	fits_write_dblarr( Name_MatMask_file, MatMask);
    
    /*
    dblarray MatMaskn;
    fits_read_dblarr( "m1.fits", MatMaskn );
    
    PowSpec PS, PSout;
    dblarray Tab;
    fits_read_dblarr("clt1537.fits",Tab);
    Tab.info();
    mrs_alloc_powspec(PS, Tab.nx()-1);
    for (int i=0; i < Tab.nx(); i++) PS.tt(i) = Tab(i);
    
    cout << PS.Lmax() << endl;
    mrs_alloc_powspec(PSout, PS.Lmax());
 
    matmask_mult_cl(PS, MatMaskn, PSout);
    mrs_write_powspec("xx1.fits", PSout);
    matmask_mult_cl(PS, MatMaskn, PSout, true);
    mrs_write_powspec("xx2.fits", PSout);
    */
    
    exit(0);
}


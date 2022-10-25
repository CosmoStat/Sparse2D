/******************************************************************************
**                   Copyright (C) 2013 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Francois Lanusse
**
**    Date:  13/12/01
**
**    File:  pkrec.cpp
**
*******************************************************************************
**
**    DESCRIPTION  Primordial power spectrum reconstruction using sparsity
**    -----------
**
**    PARAMETRES
**    ----------
**
**    RESULTS
**    -------
**
**
******************************************************************************/

#include "IM_Obj.h"
#include "IM_IO.h"
#include <iostream>
#include "shearinversions.h"

#define DEF_ZERO_PAD  512

char Name_flexion[100];
char Name_kappaE[100];
char Name_kappaB[100];

int Zero_pad = DEF_ZERO_PAD;
double pixel_size = 1.0;

extern int  OptInd;
extern char *OptArg;

extern int  GetOpt(int argc, char *const*argv, char *opts);

/*********************************************************************/

static void usage(char *argv[])
{
    fprintf(OUTMAN, "Usage: %s options pixel_size flexion_in kappaE_out \n\n", argv[0]);
    fprintf(OUTMAN, "   where options =  \n");
    
    fprintf(OUTMAN, "         [-Z ZeroPad]\n");
    fprintf(OUTMAN, "              Size of the zero padding \n");
    fprintf(OUTMAN, "              default is %d\n", DEF_ZERO_PAD);
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
static void infinit(int argc, char *argv[])
{
    int c;
#ifdef LARGE_BUFF
    int VMSSize=-1;
    Bool OptZ = False;
    char VMSName[1024] = "";
#endif
    
    /* get options */
    while ((c = GetOpt(argc,argv,"?:Z:")) != -1)
    {
        switch (c)
        {
	  case 'Z':
                /* -Z < zero padding size > */
                if (sscanf(OptArg,"%d",&Zero_pad) != 1) {
                    fprintf(OUTMAN, "bad Zero padding argument: %s\n", OptArg);
                    usage(argv);
                }
                if (Zero_pad <= 0)  Zero_pad = DEF_ZERO_PAD;
	    break;
	    case '?': usage(argv); break;
            default: usage(argv); break;
        }
    }
    if (OptInd < argc) pixel_size = atof(argv[OptInd++]);
    else usage(argv);
    
    if (OptInd < argc) strcpy(Name_flexion, argv[OptInd++]);
    else usage(argv);
    
    if (OptInd < argc) strcpy(Name_kappaE, argv[OptInd++]);
    else usage(argv);
        
    /* make sure there are not too many parameters */
    if (OptInd < argc)
    {
        fprintf(OUTMAN, "Too many parameters: %s ...\n", argv[OptInd]);
        usage(argv);
    }
#ifdef LARGE_BUFF
    if (OptZ == True) vms_init(VMSSize, VMSName, Verbose);
#endif
}

int main(int argc, char **argv) {
  
    fltarray flexion_In;
    char Cmd[512];
    Cmd[0] = '\0';
    for (int k =0; k < argc; k++) sprintf(Cmd, "%s %s", Cmd, argv[k]);   
    
    /* Get command line arguments, open input file(s) if necessary */
    infinit(argc, argv);
    
    pixel_size *= (M_PI / 180. /60.);
    fits_read_fltarr(Name_flexion,flexion_In);
    
    int Nx = flexion_In.nx();
    int Ny = flexion_In.ny();
    
    dblarray F1(Nx,Ny);
    dblarray F2(Nx,Ny);
    dblarray kappaDbl(Nx,Ny);
    fltarray kappa_out(Nx,Ny);
    
    
    for(int x=0; x < Nx; x++)
    for(int y=0; y < Ny; y++){
	F1(x,y) = flexion_In(x,y,0);
	F2(x,y) = flexion_In(x,y,1);
    }
    
    ShearInversions inv(max(Zero_pad,2*max(Nx,Ny)));
    inv.flexion2kappa(pixel_size, F1,F2,kappaDbl);
    
    for(int x=0; x < Nx; x++)
    for(int y=0; y < Ny; y++){
	kappa_out(x,y) = kappaDbl(x,y);
    }
   
    fits_write_fltarr(Name_kappaE,kappa_out);
    
    return 0;
}

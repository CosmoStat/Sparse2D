/********************************************************************************
**                   Copyright (C) 2009 by CEA  
*********************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author:  Olivier Fourt  
**
**    Date:  01/07/09
**    
**    File:  mrsp_alm2spec.cc
**
*********************************************************************************
**
**    DESCRIPTION:  Compute spectrums for polarized maps from ALM coefficients
**    -------------------
**
********************************************************************************/

#include "PolaHealpixClass.h"

char Name_file_spec[512];
char Name_file_ALM[512];

extern int OptInd;
extern char *OptArg;
extern int GetOpt(int argc, char **argv, char *opts);

bool NormALM = false;
bool PolaFastALM = true;
float ZeroPadding=0;
int Lmax = 0;
int Nside_in;

/***************************************************************************/
 
static void usage(char *argv[])
{
    fprintf(OUTMAN, "Usage: %s options ALMPola_file_input spec_file_output Nside_in \n\n", argv[0]);
    fprintf(OUTMAN, "   where options =  \n");
    
    fprintf(OUTMAN, "         [-l lmax]\n");
    fprintf(OUTMAN, "             Default is MIN(3000, 3*Nside_in) \n");
    
    fprintf(OUTMAN, "         [-p Harmonics_ZeroPadding (in [0,1] )]\n");
    fprintf(OUTMAN, "             Default is 0.\n");

    fprintf(OUTMAN, "         [-n]\n");
    fprintf(OUTMAN, "             Normalization. Default is no.\n");
    
    fprintf(OUTMAN, "         [-s]\n");
    fprintf(OUTMAN, "             No Fast ALM. Default is yes.\n");
	
    exit(-1);
}

/*********************************************************************/

/* GET COMMAND LINE ARGUMENTS */
static void sinit( int argc, char *argv[] )
{
    int c;
    // get options
    while( (c = GetOpt( argc,argv,(char *) "n:l:p:s:" )) != -1 ) 
    {
		switch (c) 
        { 
			case 'n':
				OptInd--;
           		NormALM = true;
 		        break;
 		        
 		    case 's':
				OptInd--;
           		PolaFastALM = false;
 		        break;
 		        
			case 'p':
 	     	    if( sscanf(OptArg,"%f",&ZeroPadding) != 1 ) 
                {
                	fprintf(OUTMAN, "Error: bad zero padding value: %s\n", OptArg);
		            exit(-1);
	            }
 		        break;
 		        
			case 'l':
				if( sscanf(OptArg,"%d",&Lmax) != 1 ) 
				{
		            fprintf(OUTMAN, "Error: bad lmax value: %s\n", OptArg);
		            exit(-1);
				}
				break;
				
 	     	default:
 	     		usage(argv);
 	     		break;
 		}
	} 

	//get optional input file names from trailing parameters and open files
    if( OptInd < argc )
    {
    	strcpy( Name_file_ALM, argv[OptInd] );
    	OptInd++;
    }
    else
    {
    	usage(argv);
    }
    
    if( OptInd < argc )
    {
    	strcpy( Name_file_spec, argv[OptInd] );
    	OptInd++;
    }
    else
    {
    	usage(argv);
    }
    
    if( OptInd < argc )
    {
    	Nside_in = atoi( argv[OptInd] );
    	OptInd++;
    }
    else
    {
    	usage(argv);
    }

	// make sure there are not too many parameters
	if( OptInd < argc )
    {
		fprintf( OUTMAN, "Too many parameters: %s ...\n", argv[OptInd] );
		exit(-1);
	}
}
  
/*********************************************************************/

int main( int argc, char *argv[] )
{
	//Call "mrsp_alm2spec file_in file_out"
	
	// Get command line arguments, open input file(s) if necessary
    sinit( argc, argv );
	
	if( Lmax <= 0 )
	{
		Lmax = 3 * Nside_in;
	}
	int Lmax_Tot = (int) (Lmax + Lmax * ZeroPadding);
	
	PolaAlmR ALM;
	
	ALM.read( Name_file_ALM, Nside_in, Lmax_Tot, PolaFastALM );
	
	ALM.set_flag_NormALM( NormALM );
		
	Healpix_PowSpec SigPowSpec;
    
    ALM.pola_alm2powspec_all( SigPowSpec );
    
    SigPowSpec.write( Name_file_spec );
	
	exit(0);
}
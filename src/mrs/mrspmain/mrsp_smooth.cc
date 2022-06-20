/******************************************************************************
**                   Copyright (C) 2008 by CEA  
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author:  Olivier Fourt  
**
**    Date:  15/07/09
**    
**    File:  mrsp_smooth.cc
**
*******************************************************************************
**
**    DESCRIPTION:  Smooth a polarized image 
**    -----------
**
******************************************************************************/

#include "PolaHealpixClass.h"

char flag_tqu_teb[4];
char Name_file_Map_in[512];
char Name_file_Map_out[512];
char Filter_FileName[512];
 
extern int  OptInd;
extern char *OptArg;
extern int  GetOpt(int argc, char **argv, char *opts);
 
bool NormALM = false;
bool PolaFastALM = true;
float ZeroPadding=0;
int Lmax = 0;

float Fwhm=60.;

bool OptFilter = false;
fltarray Filter;

/***************************************************************************/
 
static void usage(char *argv[])
{
    fprintf(OUTMAN, "Usage: %s options TQU/TEB Name_file_Map_in Name_file_Map_out \n\n", argv[0]);
    fprintf(OUTMAN, "   where options =  \n");
    
    fprintf(OUTMAN, "         [-F Filter_FileName]\n");
    fprintf(OUTMAN, "             Filter file name . \n");

    fprintf(OUTMAN, "         [-l lmax]\n");
    fprintf(OUTMAN, "             Default is MIN(3000, 3*nside) \n");
    
    fprintf(OUTMAN, "         [-p Harmonics_ZeroPadding (in [0,1] )]\n");
    fprintf(OUTMAN, "             Default is 0.\n");

    fprintf(OUTMAN, "         [-f]\n");
    fprintf(OUTMAN, "             Full width at half maximum (in arc minute). Gaussian beam. Default is 60 arc minutes.\n");
	
    fprintf(OUTMAN, "         [-s]\n");
    fprintf(OUTMAN, "             No Fast ALM. Default is yes.\n");

    exit(-1);
}

/*********************************************************************/

/* GET COMMAND LINE ARGUMENTS */
static void sinit(int argc, char *argv[])
{
    int c; 
    /* get options */
    while ((c = GetOpt(argc,argv,(char *) "F:f:l:p:s:")) != -1) 
    {
	switch (c) 
        {
        	case 's':
				OptInd--;
           		PolaFastALM = false;
 		        break;
 		        
			case 'f':
				if( sscanf(OptArg,"%f",&Fwhm) != 1 ) 
				{
					fprintf( OUTMAN, "Error: bad FWHM value: %s\n", OptArg );
		            exit(-1);
				}
				break; 	

			case 'p':
				if( sscanf(OptArg,"%f",&ZeroPadding) != 1 ) 
				{
		            fprintf( OUTMAN, "Error: bad zero padding  value: %s\n", OptArg );
		            exit(-1);
				}
				break;

			case 'l':
				if( sscanf(OptArg,"%d",&Lmax) != 1 ) 
				{
		            fprintf( OUTMAN, "Error: bad lmax  value: %s\n", OptArg );
		            exit(-1);
				}
				break;

			case 'F':
				if( sscanf(OptArg,"%s", Filter_FileName) != 1 ) 
				{
					fprintf( OUTMAN, "Error: bad filter file name: %s\n", OptArg );
					exit(-1);
				}
                OptFilter = true;
                fits_read_fltarr( fitsname(Filter_FileName), Filter );
                if( Filter.nx() < 1 )
                {
                	fprintf( OUTMAN, "Error: bad filter data: %s\n", OptArg );
					exit(-1);
                }
  	        	break;

			default:
				usage(argv);
				break;
 		}
	} 

	// get optional input file names from trailing parameters and open files
	if( OptInd < argc )
	{
		strcpy( flag_tqu_teb, argv[OptInd] );
		OptInd++;
	}
	else
	{
		usage(argv);
	}
	
    if( OptInd < argc )
    {
    	strcpy( Name_file_Map_in, argv[OptInd] );
    	OptInd++;
    }
    else
    {
    	usage(argv);
    }
    
    if( OptInd < argc )
    {
    	strcpy( Name_file_Map_out, argv[OptInd] );
    	OptInd++;
    }
    else
    {
    	usage(argv);
    }

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
    // Get command line arguments, open input file(s) if necessary
    sinit(argc, argv);
    
    if( (strcmp( flag_tqu_teb, "TQU" ) != 0) && (strcmp( flag_tqu_teb, "TEB" ) != 0) )
	{
		fprintf(OUTMAN, "Wrong flag parameter in %s call \n\n", argv[0]);
		
		fprintf(OUTMAN, "Usage: %s TQU Name_file_Map_in Name_file_Map_out \n\n OR \n\n %s TEB Name_file_Map_in Name_file_Map_out \n\n", argv[0], argv[0]);
		
		exit(-1);
	}
    
    PolaHdmap Map_TQU;
    
    if( strcmp( flag_tqu_teb, "TQU" ) == 0 )
	{
		Map_TQU.read( Name_file_Map_in );//Read a map in TQU
	}
	
	if( strcmp( flag_tqu_teb, "TEB" ) == 0 )
	{
		Map_TQU.read( Name_file_Map_in, true );//Read a map in TEB
	}
			
	int Nside_map = Map_TQU.get_nside();
	
	if( Lmax <= 0 )
	{
		Lmax = 3 * Nside_map;
	}
	int Lmax_Tot = (int) (Lmax + Lmax * ZeroPadding);
	
	PolaAlmR ALM;
	
    ALM.alloc( Nside_map, Lmax_Tot, PolaFastALM );
    
    ALM.pola_alm_trans( Map_TQU );

	if( OptFilter == false )
	{
		ALM.convol( Fwhm );
	}
	else
	{
		if( Filter.ny() == 3 )
		{
			fltarray Filter_T, Filter_E, Filter_B;
			Filter_T.alloc( Filter.nx() );
			Filter_E.alloc( Filter.nx() );
			Filter_B.alloc( Filter.nx() );
			
			for( int k=0; k < Filter.nx(); k++ )
			{
				Filter_T(k) = Filter(0,k);
				Filter_E(k) = Filter(1,k);
				Filter_B(k) = Filter(2,k);
			}
			
			ALM.convol( Filter_T, Filter_E, Filter_B );
		}
		else
		{
			ALM.convol( Filter );
		}
	}
	     
	ALM.pola_alm_rec( Map_TQU );
    
    Map_TQU.write( Name_file_Map_out );

	exit(0);
}


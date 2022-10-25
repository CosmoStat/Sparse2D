/********************************************************************************************
**                   Copyright (C) 2009 by CEA  
*********************************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author:  Olivier Fourt  
**
**    Date:  22/06/09
**    
**    File:  mrsp_spec.cc
**
*********************************************************************************************
**
**    DESCRIPTION:  Power spectrums and cross spectrums computation for polarized maps
**    -------------------
**
*********************************************************************************************/

#include "PolaHealpixClass.h"

char flag_tqu_teb[4];
char Name_file_Map[512];
char Name_file_spec[512];

extern int OptInd;
extern char *OptArg;
extern int GetOpt(int argc, char **argv, char *opts);

bool NormALM = false;
bool PolaFastALM = true;
float ZeroPadding=0;
int Lmax = 0;

/***************************************************************************/
 
static void usage(char *argv[])
{
    fprintf(OUTMAN, "Usage: %s options TQU/TEB PolaMap_File_input spec_file_output \n\n", argv[0]);
    fprintf(OUTMAN, "   where options =  \n");
    
    fprintf(OUTMAN, "         [-l lmax]\n");
    fprintf(OUTMAN, "             Default is MIN(3000, 3*nside) \n");
    
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
    /* get options */
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
                	fprintf(OUTMAN, "Error: bad zero padding  value: %s\n", OptArg);
		            exit(-1);
	            }
 		        break;
 		        
			case 'l':
				if( sscanf(OptArg,"%d",&Lmax) != 1 ) 
				{
		            fprintf(OUTMAN, "Error: bad lmax  value: %s\n", OptArg);
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
    	strcpy( Name_file_Map, argv[OptInd] );
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

	/* make sure there are not too many parameters */
	if( OptInd < argc )
    {
		fprintf( OUTMAN, "Too many parameters: %s ...\n", argv[OptInd] );
		exit(-1);
	}
}
  
/*********************************************************************/

int main( int argc, char *argv[] )
{
	sinit(argc, argv);
		
	PolaHdmap Map_TQU;
	
	if( (strcmp( flag_tqu_teb, "TQU" ) != 0) && (strcmp( flag_tqu_teb, "TEB" ) != 0) )
	{
		fprintf(OUTMAN, "Wrong flag parameter in %s call \n\n", argv[0]);
		
		fprintf(OUTMAN, "Usage: %s TQU PolaMap_File_input spec_file_output \n\n OR \n\n %s TEB PolaMap_File_input spec_file_output \n\n", argv[0], argv[0]);
		
		exit(-1);
	}
	
	if( strcmp( flag_tqu_teb, "TQU" ) == 0 )
	{
		Map_TQU.read( Name_file_Map );// Read a map in TQU
	}
	
	if( strcmp( flag_tqu_teb, "TEB" ) == 0 )
	{
		Map_TQU.read( Name_file_Map, true );// Read a map in TEB
	}
	
	int Nside_map = Map_TQU.get_nside();
	
	if( Lmax <= 0 )
	{
		Lmax = 3 * Nside_map;
	}
	int Lmax_Tot = (int) (Lmax + Lmax * ZeroPadding);
	
	PolaAlmR ALM;
    
    ALM.alloc( Nside_map, Lmax_Tot, PolaFastALM );
    
    ALM.set_flag_NormALM( NormALM );
    
    ALM.pola_alm_trans( Map_TQU );
    
    Healpix_PowSpec SigPowSpec;
    
    ALM.pola_alm2powspec_all( SigPowSpec );
    
    SigPowSpec.write( Name_file_spec );
	
	exit(0);
}
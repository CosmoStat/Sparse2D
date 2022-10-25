/******************************************************************************
**                   Copyright (C) 2009 by CEA  
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Olivier Fourt
**
**    Date:  22/06/09
**    
**    File:  mrsp_teb2tqu.cc
**
*******************************************************************************
**
**    DESCRIPTION  convert Healpix polarized map from TQU fields to TEB fields
**    ----------- 
**                 
******************************************************************************/

#include "PolaHealpixClass.h"

char Name_file_Map1[512];
char Name_file_Map2[512];

extern int OptInd;
extern char *OptArg;
extern int GetOpt(int argc, char **argv, char *opts);

bool PolaFastALM = true;

/***************************************************************************/
 
static void usage(char *argv[])
{
    fprintf(OUTMAN, "Usage: %s options PolaMapTEB_File_input PolaMapTQU_file_output \n\n", argv[0]);
    fprintf(OUTMAN, "   where options =  \n");
    
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
    while( (c = GetOpt( argc,argv,(char *) "s:" )) != -1 ) 
    {
		switch (c) 
        { 
 		        
 		    case 's':
				OptInd--;
           		PolaFastALM = false;
 		        break;
 		         	     	
 	     	default:
 	     		usage(argv);
 	     		break;
 		}
	} 

	//get optional input file names from trailing parameters and open files
    
    if( OptInd < argc )
    {
    	strcpy( Name_file_Map1, argv[OptInd] );
    	OptInd++;
    }
    else
    {
    	usage(argv);
    }
    
    if( OptInd < argc )
    {
    	strcpy( Name_file_Map2, argv[OptInd] );
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
	//Call "mrsp_teb2tqu file_in file_out"
	
	sinit(argc, argv);
	
	PolaHdmap Map_TQU;
	
	Map_TQU.read( Name_file_Map1, true );// Read a map in TEB
	
	Map_TQU.swap_tqu_teb( PolaFastALM );
	
	Map_TQU.write( Name_file_Map2 );// Write map in TQU
	
	exit(0);
}


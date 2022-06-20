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
**    Date:  02/07/09
**    
**    File:  mrsp_resize.cc
**
*******************************************************************************
**
**    DESCRIPTION:  resize a polarized Healpix map 
**    -----------
**
******************************************************************************/

#include "PolaHealpixClass.h"

char Name_Imag_In[512]; /* input file image */
char Name_Imag_Out[512]; /* output file name */

extern int  OptInd;
extern char *OptArg;
extern int  GetOpt(int argc, char **argv, char *opts);

bool ViaAlm = false;
bool PolaFastALM = true;
int Nside = 0;

/***************************************************************************/
 
static void usage(char *argv[])
{
    fprintf( OUTMAN, "Usage: %s options in_map out_map Nside_out \n\n", argv[0] );
    fprintf(OUTMAN, "   where options =  \n");
    fprintf(OUTMAN, "         [-a]\n");
    fprintf(OUTMAN, "             Resizing via alm trans.\n");
    fprintf(OUTMAN, "         [-s]\n");
    fprintf(OUTMAN, "             No Fast ALM. Default is yes, no use if alm are not used.\n");
    
    exit(-1);
}
  
/****************************************************************************/

/* GET COMMAND LINE ARGUMENTS */
static void sinit( int argc, char *argv[] )
{
    int c;
    // get options
    while( (c = GetOpt( argc,argv,(char *) "a:s:" )) != -1 ) 
    {
		switch (c) 
        { 
			case 'a':
				OptInd--;
           		ViaAlm = true;
 		        break;
 		        
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
    	strcpy( Name_Imag_In, argv[OptInd] );
    	OptInd++;
    }
    else
    {
    	usage(argv);
    }
    
    if( OptInd < argc )
    {
    	strcpy( Name_Imag_Out, argv[OptInd] );
    	OptInd++;
    }
    else
    {
    	usage(argv);
    }
    
    if( OptInd < argc )
    {
    	Nside = atoi( argv[OptInd] );
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
    /* Get command line arguments, open input file(s) if necessary */
    
    sinit( argc, argv );
            
    PolaHdmap Map_TQU;
    	
	Map_TQU.read( Name_Imag_In );
	
	int NsideIn = Map_TQU.get_nside();

	if( Nside <= 0 )
	{
		Nside = NsideIn;
	}
 
	PolaHdmap outmap;
	outmap.alloc( Nside, Map_TQU.flag_teb(), Map_TQU.flag_nested() );
	
	if( ViaAlm == false )
	{
		outmap.import( Map_TQU );
	}
	else
	{
		outmap.import_via_alm( Map_TQU, PolaFastALM );
	}

	outmap.write( Name_Imag_Out );

	exit(0);
}


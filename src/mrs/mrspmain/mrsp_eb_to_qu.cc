/******************************************************************************************
**                   Copyright (C) 2009 by CEA  
*******************************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Olivier Fourt
**
**    Date:  26/06/09
**    
**    File:  mrsp_eb_to_qu.cc
**
********************************************************************************************
**
**    DESCRIPTION  convert Healpix polarized map from E and B fields to Q and U fields
**    ----------- 
**                 
*******************************************************************************************/

#include "HealpixClass.h"
// #include "map_alm_qu_eb.h"

char Name_file_Map1[512];
char Name_file_Map2[512];
char Name_file_Map3[512];
char Name_file_Map4[512];

extern int OptInd;
extern char *OptArg;
extern int GetOpt(int argc, char **argv, char *opts);

bool PolaFastALM = true;

/***********************************************************************************************/
 
static void usage(char *argv[])
{
    fprintf(OUTMAN, "Usage: %s options Map_E_File_input Map_B_File_input Map_Q_File_output Map_U_File_output \n\n", argv[0]);
    fprintf(OUTMAN, "   where options =  \n");
    
    fprintf(OUTMAN, "         [-s]\n");
    fprintf(OUTMAN, "             No Fast ALM. Default is yes.\n");
	
    exit(-1);
}

/************************************************************************************************/

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
    
    if( OptInd < argc )
    {
    	strcpy( Name_file_Map3, argv[OptInd] );
    	OptInd++;
    }
    else
    {
    	usage(argv);
    }
    
    if( OptInd < argc )
    {
    	strcpy( Name_file_Map4, argv[OptInd] );
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
	sinit(argc, argv);
	
	Hdmap Map_Q, Map_U;
	
	Map_Q.read( Name_file_Map1 );// Read map E
	Map_U.read( Name_file_Map2 );// Read map B
	
	if( Map_Q.Scheme() != Map_U.Scheme() )
	{
		fprintf( OUTMAN, "Maps are not in the same scheme in: %s ...\n", argv[0] );
		exit(-1);
	}
	
	if( Map_Q.nside() != Map_U.nside() )
	{
		fprintf( OUTMAN, "Maps have not the same Nside in: %s ...\n", argv[0] );
		exit(-1);
	}
	
	bool nested;
	
	if( Map_Q.Scheme() != DEF_MRS_ORDERING )
     {
     	nested = true;
     }
	
	int Nside = Map_Q.nside();
	int Lmax = 3*Nside;
	int Mmax = Lmax;
	
	arr<double> weight_T;
	weight_T.alloc( 2*Nside );
	if( PolaFastALM == false )
	{
		read_weight_ring( "/Users/Shared/software/Healpix_2.10/data/", Nside, weight_T );
		
		for(int m=0; m < weight_T.size(); ++m)
		{
			weight_T[m]+=1;
		}
	}
	else
	{
		weight_T.fill(1);
	}
	
	if( nested == true )
    {
    	Map_Q.swap_scheme();
    	Map_U.swap_scheme();
    }//ALM Trans Must be done in RING scheme
	
	Alm<xcomplex<double> > ALM_E(Lmax,Mmax), ALM_B(Lmax,Mmax);
	
	map2alm_iter( Map_Q, ALM_E, 0, weight_T );
    map2alm_iter( Map_U, ALM_B, 0, weight_T );
    	
    // alm2map_pol_QU( ALM_E, ALM_B, Map_Q, Map_U );
    alm2map_spin(ALM_E,ALM_B,Map_Q,Map_U,2);

	if( nested == true )
    {
    	Map_Q.swap_scheme();
    	Map_U.swap_scheme();
    }//ALM Trans Must be done in RING scheme
	
	Map_Q.write( Name_file_Map3 );// Write map Q
	Map_U.write( Name_file_Map4 );// Write map B
	
	exit(0);
}


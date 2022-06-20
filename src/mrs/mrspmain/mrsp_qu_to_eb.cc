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
**    File:  mrsp_qu_to_eb.cc
**
********************************************************************************************
**
**    DESCRIPTION  convert Healpix polarized map from Q and U fields to E and B fields
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
Bool Verbose=False;

float  ZeroPadding=0.0;
int Lmax = 0.;

/***********************************************************************************************/
 
static void usage(char *argv[])
{
    fprintf(OUTMAN, "Usage: %s options Map_Q_File_input Map_U_File_input Map_E_File_output Map_B_File_output \n\n", argv[0]);
    fprintf(OUTMAN, "   where options =  \n");
    
//  fprintf(OUTMAN, "         [-f]\n");
//  fprintf(OUTMAN, "             No Fast ALM. Default is yes.\n");
    fprintf(OUTMAN, "         [-l lmax]\n");
    fprintf(OUTMAN, "             Default is MIN(%d, 2*nside) \n",ALM_MAX_L);
	fprintf(OUTMAN, "         [-v Verbose]\n");
    exit(-1);
}

/************************************************************************************************/

/* GET COMMAND LINE ARGUMENTS */
static void sinit( int argc, char *argv[] )
{
    int c;
    // get options
    while( (c = GetOpt( argc,argv,(char *) "l:fv" )) != -1 ) 
    {
		switch (c) 
        { 
 		        
 		    case 'f':
				OptInd--;
           		PolaFastALM = false;
 		        break;
	        case 'l':
                if (sscanf(OptArg,"%d",&Lmax) != 1) 
                {
		            fprintf(OUTMAN, "Error: bad lmax  value: %s\n", OptArg);
		            exit(-1);
                }
                break;
     	    case 'v': Verbose = True; break;
		         	     	
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
    int k;
    fitsstruct Header, HD1;
    char Cmd[512];
    Cmd[0] = '\0';
    for (k =0; k < argc; k++) sprintf(Cmd, "%s %s", Cmd, argv[k]);

	sinit(argc, argv);

    if (Verbose == True)
    { 
        cout << "# PARAMETERS: " << endl;
        if (PolaFastALM == false )  cout << "# Do not use fast alm calculation. " <<  Lmax << endl;  
        if (Lmax > 0)  cout << "# Lmax = " <<  Lmax << endl;    
    }
            
	Hdmap Map_Q, Map_U;
	
	Map_Q.read( Name_file_Map1 );// Read map Q
	Map_U.read( Name_file_Map2 );// Read map U
	
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
	
 	
	int Nside = Map_Q.nside();
    int Lmax_Tot = mrs_get_lmax (Lmax,  Nside,   ZeroPadding);
 	int Mmax = Lmax_Tot;
	
	arr<double> weight_T;
	weight_T.alloc( 2*Nside );
	if( PolaFastALM == false )
	{
        cout << "Read weights in $HEALPIX/data ...  " << endl;
		read_weight_ring( "$HEALPIX/data/", Nside, weight_T );
 		for(int m=0; m < (int) weight_T.size(); ++m)   weight_T[m]+=1;
 	}
	else weight_T.fill(1);
 	
	 	
	Alm<xcomplex<double> > ALM_E(Lmax_Tot,Mmax), ALM_B(Lmax_Tot,Mmax);
	
	// map2alm_pol_iter_QU( Map_Q, Map_U, ALM_E, ALM_B, 0, weight_T );
    map2alm_spin(Map_Q,Map_U,ALM_E,ALM_B,2,weight_T,false);

    alm2map( ALM_E, Map_Q );
    alm2map( ALM_B, Map_U );
		
	Map_Q.write( Name_file_Map3 );// Write map E
	Map_U.write( Name_file_Map4 );// Write map B
	
	exit(0);
}


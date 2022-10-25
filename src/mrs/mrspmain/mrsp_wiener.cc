/*********************************************************************************************************************
**                   Copyright (C) 2009 by CEA  
**********************************************************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author:  Olivier Fourt  
**
**    Date:  15/07/09
**    
**    File:  mrsp_wiener.cc
**
***********************************************************************************************************************
**
**    DESCRIPTION:  Wiener filtering for healpix polarized maps, true signal spectrums and noise spectrums in imput
**    -------------------
**
***********************************************************************************************************************/

#include "PolaHealpixClass.h"

char flag_tqu_teb[4];
char Name_file_Map_in[512];
char Name_file_Map_out[512];
char Name_file_Signal_spec[512];
char Name_file_Noise_spec[512];

extern int OptInd;
extern char *OptArg;
extern int GetOpt(int argc, char **argv, char *opts);

bool input_spectrum = false;
bool NormALM = false;
bool PolaFastALM = true;
float ZeroPadding=0;
int Lmax = 0;

/***************************************************************************/
 
static void usage(char *argv[])
{
    fprintf(OUTMAN, "Usage: %s options TQU/TEB Name_file_Noise_spec Name_file_Map_in Name_file_Map_out \n\n", argv[0]);
    fprintf(OUTMAN, "   where options =  \n");
    
    fprintf(OUTMAN, "         [-l lmax]\n");
    fprintf(OUTMAN, "             Default is MIN(3000, 3*nside) \n");
    
    fprintf(OUTMAN, "         [-p Harmonics_ZeroPadding (in [0,1] )]\n");
    fprintf(OUTMAN, "             Default is 0.\n");

    fprintf(OUTMAN, "         [-n]\n");
    fprintf(OUTMAN, "             Normalization. Default is no.\n");
    
    fprintf(OUTMAN, "         [-s]\n");
    fprintf(OUTMAN, "             No Fast ALM. Default is yes.\n");
    
    fprintf(OUTMAN, "         [-i]\n");
    fprintf(OUTMAN, "             Input signal spectrum. Default is no, spectrum estimated from image in spectrum and noise spectrum.\n");
	
    exit(-1);
}

/*********************************************************************/

/* GET COMMAND LINE ARGUMENTS */
static void sinit( int argc, char *argv[] )
{
    int c;
    // get options
    while( (c = GetOpt( argc,argv,(char *) "n:l:p:s:i:" )) != -1 ) 
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
				
			case 'i':
           		input_spectrum = true;
           		if( sscanf(OptArg,"%s", Name_file_Signal_spec) != 1 ) 
				{
		            fprintf(OUTMAN, "Error: bad file name: %s\n", OptArg);
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
		strcpy( flag_tqu_teb, argv[OptInd] );
		OptInd++;
	}
	else
	{
		usage(argv);
	}
	 
    if( OptInd < argc )
    {
    	strcpy( Name_file_Noise_spec, argv[OptInd] );
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
	// Get command line arguments, open input file(s) if necessary
    sinit(argc, argv);
    
	PolaHdmap Map_TQU;
	
	if( (strcmp( flag_tqu_teb, "TQU" ) != 0) && (strcmp( flag_tqu_teb, "TEB" ) != 0) )
	{
		fprintf(OUTMAN, "Wrong flag parameter in %s call \n\n", argv[0]);
		
		fprintf(OUTMAN,"Usage: %s TQU Name_file_Noise_spec Name_file_Map_in Name_file_Map_out \n\n OR \n\n %s TEB Name_file_Noise_spec Name_file_Map_in Name_file_Map_out \n\n",argv[0],argv[0]);
		
		exit(-1);
	}
	
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
	
	Healpix_PowSpec signal_spec;
    Healpix_PowSpec noise_spec;
    
    if( input_spectrum == true )
    {
    	signal_spec.read( Name_file_Signal_spec );
    }
        
    noise_spec.read( Name_file_Noise_spec );
    
    if( signal_spec.Lmax() != noise_spec.Lmax() )
    {
    	fprintf(OUTMAN, "Signal spectrums and noise spectrums have not the same lmax parameter in %s call \n\n", argv[0]);
    	
    	exit(-1);
    }
    
    if( signal_spec.Lmax() != Lmax_Tot )
    {
    	fprintf(OUTMAN, "Input spectrums and lmax parameter for pola_almtrans are not the same in %s call \n\n", argv[0]);
    	
    	exit(-1);
    }
    
	PolaAlmR ALM;
	
    ALM.alloc( Nside_map, Lmax_Tot, PolaFastALM );
        
    ALM.set_flag_NormALM( NormALM );
    
    ALM.pola_alm_trans( Map_TQU );
    
    if( input_spectrum == false )
    {
    	ALM.pola_alm2powspec_all( signal_spec );
    	
    	for( int ind=0; ind <= Lmax_Tot; ind ++)
    	{
    		signal_spec.tt( ind ) = signal_spec.tt( ind ) - noise_spec.tt( ind );
    		if( signal_spec.tt( ind ) < 0 )
    		{
    			signal_spec.tt( ind ) = 0.0;
    		}
    		
    		signal_spec.ee( ind ) = signal_spec.ee( ind ) - noise_spec.ee( ind );
    		if( signal_spec.ee( ind ) < 0 )
    		{
    			signal_spec.ee( ind ) = 0.0;
    		}
    		
    		signal_spec.bb( ind ) = signal_spec.bb( ind ) - noise_spec.bb( ind );
    		if( signal_spec.bb( ind ) < 0 )
    		{
    			signal_spec.bb( ind ) = 0.0;
    		}
    	}
    }
    
    ALM.wiener( noise_spec, signal_spec );
    
    ALM.pola_alm_rec( Map_TQU );
    
    Map_TQU.write( Name_file_Map_out );
    	
	exit(0);
}


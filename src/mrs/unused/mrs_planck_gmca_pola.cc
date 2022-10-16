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
**    Date:  17/06/10
**    
**    File:  mrs_planck_gmca_pola.cc
**
*******************************************************************************
**
**    DESCRIPTION  Generalized Morphological Component Analysis for polarized planck datas
**    ----------- 
**                 
**    Usage: mrs_planck_gmca_pola options cube output
**
******************************************************************************/

#include "HealpixClass.h"
#include "map_alm_qu_eb.h"
#include "MRS_Sparse.h"
#include "GMCA_POLA.h"

/****************************************************************************/
/****************************************************************************/


char Name_Cube_In_Q[256];
char Name_Cube_In_U[256];
char Name_Out[256];
char BeamChannels_List[256];
char STD_Cube_In_Q[256];
char STD_Cube_In_U[256];

int NbrScale2d = 5;

// type_transform  Transform=TO_MALLAT;

extern int  OptInd;
extern char *OptArg;

extern int  GetOpt(int argc, char *const*argv, char *opts);

Bool Verbose=False;

int NbrCannels=0;  // Number of cannels 
int NbrSources=0;  // Number of Sources
float KMin = 3; // Last K-Mad Thresholding
int Max_GMCA_Iter = DEF_GMCA_MAX_ITER; // Maximum number of iterations in the GMCA algorithm

type_gmca_threshold_decrease TypeThresholdDecrease = GMCA_THRES_MAD; // Select the decreasing threshold strategy 

bool Nested = true;

bool FixCMB = false;

char list_channels_in[13];

bool WMAP = false;
bool HFI = false;
bool LFI = false;

bool ASCII = false;

bool estim_bias = false;
bool estim_std_src = false;

bool conv_EB = false;

int EstNside = 0;
int Nside = 0;

int channel_sel = 0;

/*********************************************************************/

/*static int max_scale_number (int Nc)
{
    int Nmin = Nc;
    int ScaleMax;

    ScaleMax=iround(log((float)Nmin/(float)MAX_SIZE_LAST_SCALE) / log(2.)+ 1.);
    return (ScaleMax);
}*/

/*********************************************************************/
 
static void usage(char *argv[])
{
    fprintf(OUTMAN, "Usage: %s options input_data_Q input_data_U output_source [output_channel] \n\n", argv[0]);
    fprintf(OUTMAN, "   where options =  \n");
    
    fprintf(OUTMAN, "         [-i NbrIter]\n");
    fprintf(OUTMAN, "             Number of iterations. Default is no %d. \n", Max_GMCA_Iter);
   
    fprintf(OUTMAN, "         [-N Number_of_Sources]\n");
    fprintf(OUTMAN, "             Number of sources. \n");    
        
    fprintf(OUTMAN, "         -b Mixing matrix CMB column fixed \n");
    
    fprintf(OUTMAN, "         -c Conversion of the datas input maps from QU to EB mode before separation. \n");
    fprintf(OUTMAN, "         The estimated sources maps will also be in EB mode. \n");
    fprintf(OUTMAN, "         Does not apply to the standard deviations maps (option -V) \n");
    
    //fprintf(OUTMAN, "         -W WMAP channels selected \n");
    
    //fprintf(OUTMAN, "         -L LFI channels selected \n");
    
    //fprintf(OUTMAN, "         -H HFI channels selected \n");
    
    fprintf(OUTMAN, "         -T for the input file, use a ASCII text with name and path of Healpix maps to be read \n");
    
    fprintf(OUTMAN, "         [-e EstNside]\n");
    fprintf(OUTMAN, "             Nside parameter used for the GMCA. \n");
    
    fprintf(OUTMAN, "         [-n Nside]\n");
    fprintf(OUTMAN, "              Nside parameter for recontruction when input maps have different sizes\n");
    fprintf(OUTMAN, "              only used if option -T is also used. \n");
    
    fprintf(OUTMAN, "         [-B BeamFile]\n");
    fprintf(OUTMAN, "             fits file with the beams of the channels for the estimation of bias on sources after GMCA. \n");
    
    fprintf(OUTMAN, "         [-V std_InputFile]\n");
    fprintf(OUTMAN, "             fits file with the standard deviations of the channels for the estimation of standard deviations on sources after GMCA. \n");
    fprintf(OUTMAN, "             Work only with Planck channels OR WMAP channels, warning message otherwise. \n");
    fprintf(OUTMAN, "             If the otpion -T is used, ASCII text file with name and path of Healpix maps to be read. \n");
    
    fprintf(OUTMAN, "         [-C channel_sel]\n");
    fprintf(OUTMAN, "             Channel selection, depends of the value: 1 Planck HFI, 2 Planck LFI, 3 WMAP.\n");
    fprintf(OUTMAN, "             4 Planck LFI+HFI, 5 Planck HFI+WMAP, 6 Planck LFI+WMAP, 7 Planck LFI+HFI+WMAP.\n");
    fprintf(OUTMAN, "             8 list_channel_in string with 1 character for each channel selected:\n");
    fprintf(OUTMAN, "             '1', '2', '3', '4', '5', '6', '7' for Planck channels\n");
    fprintf(OUTMAN, "             'K', 'A', 'Q', 'V', 'W' for WMAP channels, 'A' is for Ka channel\n");
    
    exit(-1);
}

/*********************************************************************/

/* GET COMMAND LINE ARGUMENTS */

static void transinit(int argc, char *argv[])
{
    int c;
#ifdef LARGE_BUFF
    int VMSSize=-1;
    Bool OptZ = False;
    char VMSName[1024] = "";
#endif   

    /* get options */
	while ((c = GetOpt(argc,argv,"i:n:N:e:K:B:V:C:cWLHTvbzZ:")) != -1) 
    {
		switch (c) 
        {  
	  		case 'N':
				if (sscanf(OptArg,"%d",&NbrSources) != 1) 
                {
		    		fprintf(OUTMAN, "bad number of sources: %s\n", OptArg);
		   			exit(-1);
				}
 				break;
 			
 			case 'e':
				if (sscanf(OptArg,"%d",&EstNside) != 1) 
                {
		    		fprintf(OUTMAN, "Bad EstNside value: %s\n", OptArg);
		   			exit(-1);
				}
				if( is_power_of_2( EstNside ) == False )
				{
					fprintf(OUTMAN, "EstNside value read is not a valid Nside value: %s\n", OptArg);
		   			exit(-1);
				}
 				break;
 				
			case 'i':
				if (sscanf(OptArg,"%d",&Max_GMCA_Iter) != 1) 
                {
		    		fprintf(OUTMAN, "bad number of iterations: %s\n", OptArg);
		    		exit(-1);
				}
 				break;
 				
			case 'n':
				if( sscanf( OptArg, "%d", &Nside) != 1 ) 
                {
		    		fprintf(OUTMAN, "Bad Nside value: %s\n", OptArg);
		    		exit(-1);
				}
								
				if( is_power_of_2( Nside ) == False )
				{
					fprintf(OUTMAN, "Nside value read is not a valid Nside value: %s\n", OptArg);
		   			exit(-1);
				}
 				break;
				
			case 'K':
				if (sscanf(OptArg,"%f",&KMin) < 0) 
                {
		    		fprintf(OUTMAN, "bad value of last K \n");
		   			exit(-1);
				}
 				break;
 			
 			case 'C':
 				if( sscanf( OptArg, "%d", &channel_sel) != 1 ) 
                {
		    		fprintf(OUTMAN, "Bad channel list value selected: %s\n", OptArg);
		    		exit(-1);
				}
				
				switch (channel_sel)
				{
					case 1:
						HFI = true;
						break;
					
					case 2:
						LFI = true;
						break;
					
					case 3:
						WMAP = true;
						break;
					
					case 4:
						LFI = true;
						HFI = true;
						break;
						
					case 5:
						HFI = true;
						WMAP = true;
						break;
					
					case 6:
						LFI = true;
						WMAP = true;
						break;
					
					case 7:
						LFI = true;
						HFI = true;
						WMAP = true;
						break;
					
					case 8:
						strcpy( list_channels_in, argv[OptInd++] );
						break;
					
					default:
						fprintf(OUTMAN, "Unknown channel list value selected: %s\n", OptArg);
						usage(argv);
				}
				break;
			
			case 'c':
				conv_EB = true;
				break;
 			
 			case 'b':
 				FixCMB = true;
 				//OptInd = OptInd - 1;
 				break;
				 			
 			case 'W':
 				WMAP = true;
 				//OptInd = OptInd - 1;
 				break;
 			
 			case 'L':
 				LFI = true;
 				//OptInd = OptInd - 1;
 				break;
 			
 			case 'H':
 				HFI = true;
 				//OptInd = OptInd - 1;
 				break;
 			
 			case 'T':
 				ASCII = true;
 				//OptInd = OptInd - 1;
 				break;
 			
			case 'v': Verbose = True;
				break;
			
			case 'B':
				estim_bias = true;
				strcpy( BeamChannels_List, argv[OptInd-1] );
 				break;
 			
 			case 'V':
				estim_std_src = true;
				strcpy( STD_Cube_In_Q, argv[OptInd-1] );
				strcpy( STD_Cube_In_U, argv[OptInd] );
				OptInd = OptInd + 1;
 				break;
				 				
#ifdef LARGE_BUFF
			case 'z':
	        	if (OptZ == True)
				{
					fprintf(OUTMAN, "Error: Z option already set...\n");
					exit(-1);
                }
				OptZ = True;
				break;
				
            case 'Z':
	        	if (sscanf(OptArg,"%d:%s",&VMSSize, VMSName) < 1)
				{
		   			fprintf(OUTMAN, "Error: syntaxe is Size:Directory ... \n");
		   			exit(-1);
				}
	        	
	        	if (OptZ == True)
				{
					fprintf(OUTMAN, "Error: z option already set...\n");
					exit(-1);
                }
				OptZ = True;
                break;
                
#endif
			case '?':
				usage(argv);
		}
	}
      
	
       /* get optional input file names from trailing 
          parameters and open files */
/*
    if( OptInd < argc )
	{
		strcpy( list_channels_in, argv[OptInd++] );
	}
    else
    {
    	usage(argv);
    }
*/
    
	if( OptInd < argc )
	{
		strcpy( Name_Cube_In_Q, argv[OptInd++] );
	}
    else
    {
    	usage(argv);
    }
    
    if( OptInd < argc )
	{
		strcpy( Name_Cube_In_U, argv[OptInd++] );
	}
    else
    {
    	usage(argv);
    }

	if( OptInd < argc )
	{
		strcpy( Name_Out, argv[OptInd++] );
	}
    else
    {
    	usage(argv);
    }
        
	/* make sure there are not too many parameters */
	if (OptInd < argc)
        {
		fprintf(OUTMAN, "Error: too many parameters: %s ...\n", argv[OptInd]);
		exit(-1);
	}
#ifdef LARGE_BUFF
    if (OptZ == True) vms_init(VMSSize, VMSName, Verbose);
#endif
}

/*********************************************************************/
  
int main(int argc, char *argv[]) 
{
	fltarray Dat_file_in_Q, Dat_file_in_U;
    /* Get command line arguments, open input file(s) if necessary */
	fitsstruct Header;
	char Cmd[512];
	Cmd[0] = '\0';
	for (int k =0; k < argc; k++) sprintf(Cmd, "%s %s", Cmd, argv[k]);
	
	cout << "argc: " << argc << '\n';
   
    transinit(argc, argv);
	
	//Fix selected channels
	intarray set_planck_channel(12);
	set_planck_channel.init(0);
	
	for( int k=0; k < 13; k++ )
	{
		switch( list_channels_in[k] )
		{
			case '1':
				set_planck_channel(0) = 1;
				break;
			
			case '2':
				set_planck_channel(1) = 1;
				break;
			
			case '3':
				set_planck_channel(2) = 1;
				break;
			
			case '4':
				set_planck_channel(3) = 1;
				break;
			
			case '5':
				set_planck_channel(4) = 1;
				break;
			
			case '6':
				set_planck_channel(5) = 1;
				break;
			
			case '7':
				set_planck_channel(6) = 1;
				break;
						
			case 'K':
				set_planck_channel(7) = 1;// K WMAP Band
				break;
			
			case 'A':
				set_planck_channel(8) = 1;// Ka WMAP Band
				break;
			
			case 'Q':
				set_planck_channel(9) = 1;// Q WMAP Band
				break;
			
			case 'V':
				set_planck_channel(10) = 1;// V WMAP Band
				break;
			
			case 'W':
				set_planck_channel(11) = 1;// W WMAP Band
				break;
			
			default:
				break;
		}
	}
	
	if( LFI == true )
	{
		set_planck_channel(0) = 1;
		set_planck_channel(1) = 1;
		set_planck_channel(2) = 1;
	}
	if( HFI == true )
	{
		set_planck_channel(3) = 1;
		set_planck_channel(4) = 1;
		set_planck_channel(5) = 1;
		set_planck_channel(6) = 1;
	}
	if( WMAP == true )
	{
		set_planck_channel(7) = 1;
		set_planck_channel(8) = 1;
		set_planck_channel(9) = 1;
		set_planck_channel(10) = 1;
		set_planck_channel(11) = 1;
	}
	//set_planck_channel.display( 12 );
	
	intarray list_channels_set;
	list_channels_set.alloc( set_planck_channel.total() );
	int nb_channel = 0;//final: number of selected channels

	for( int n=0; n < 12; n++ )
	{
		if( set_planck_channel(n) == 1 )
		{
			list_channels_set( nb_channel ) = n;
			nb_channel = nb_channel + 1;
		}
	}
	
	if( nb_channel == 0 )
	{
		fprintf( OUTMAN, "No Channels were selected" );
		usage(argv);
		exit(-1);
	}
	
	if (Verbose == True)
    {
		cout << "Filename Q in = " << Name_Cube_In_Q << endl;
		cout << "Filename U in = " << Name_Cube_In_U << endl;
		cout << "Filename out = " << Name_Out  << endl;
		cout << "NbrGMCA_Iter = " <<  Max_GMCA_Iter << endl;
		cout << "Mixing matrix CMB column fixed: " << FixCMB << '\n';
		cout << "Conversion QU to EB used for input data: " << conv_EB << '\n';
		list_channels_set.display( 12 );
		if( estim_bias == true )
		{
			cout << "Channels beams File: " << BeamChannels_List << '\n';
		}
		if( estim_std_src == true )
		{
			cout << "Channels Q maps of standard deviations: " << STD_Cube_In_Q << '\n';
			cout << "Channels U maps of standard deviations: " << STD_Cube_In_U << '\n';
		}
	}
	
	if (Verbose == True) cout << "\n Reading the data"<< endl;
    
    fltarray Dat_work_Q, Dat_list_Q, Dat_work_U, Dat_list_U;
    float Mean_Q, Mean_U;
	Hdmap ReadMap_Q, TempMap_Q, TempMap_EstNside_Q, ReadMap_U, TempMap_U, TempMap_EstNside_U;
    		
    if( ASCII == false )
    {     
		fits_read_fltarr( Name_Cube_In_Q, Dat_file_in_Q );
		fits_read_fltarr( Name_Cube_In_U, Dat_file_in_U );
		
		if( (Dat_file_in_Q.nx() != Dat_file_in_U.nx()) || (Dat_file_in_Q.ny() != Dat_file_in_U.ny()) )
		{
			fprintf( OUTMAN, "Error table Q and table U don't have the same size:\n" );
			Dat_file_in_Q.info("Table Q");
			Dat_file_in_U.info("Table U");
			usage(argv);
			exit(-1);
		}
		
		Nside = sqrt( Dat_file_in_Q.nx() / 12 );
		if( EstNside == 0 )
		{
			EstNside = Nside;
		}
		int Nz = Dat_file_in_Q.ny();
		
		//centering data and resizing to new EstNside for GMCA only !!!
		TempMap_Q.alloc( Nside, Nested );
		TempMap_U.alloc( Nside, Nested );
		
		if( EstNside != Nside)
		{
			TempMap_EstNside_Q.alloc( EstNside, Nested );
			TempMap_EstNside_U.alloc( EstNside, Nested );
		
			Dat_work_Q.alloc( TempMap_EstNside_Q.Npix(), nb_channel );
			Dat_work_U.alloc( TempMap_EstNside_U.Npix(), nb_channel );
		}
		else
		{
			Dat_work_Q.alloc( TempMap_Q.Npix(), nb_channel );
			Dat_work_U.alloc( TempMap_U.Npix(), nb_channel );
		}
		
		if( Nz == 12 )
		{//12 maps read selection of maps for GMCA must be done
			Dat_list_Q.alloc( TempMap_Q.Npix(), nb_channel );
			Dat_list_U.alloc( TempMap_U.Npix(), nb_channel );
		
			for( int k=0; k < nb_channel; k++ )
			{
				for( int i=0; i < Dat_list_Q.nx(); i++ )
				{
					TempMap_Q[i] = Dat_file_in_Q( i, list_channels_set(k) );
					Dat_list_Q(i,k) = TempMap_Q[i];
					TempMap_U[i] = Dat_file_in_U( i, list_channels_set(k) );
					Dat_list_U(i,k) = TempMap_U[i];
				}
			
				if( EstNside != Nside)
				{
					TempMap_EstNside_Q.import_via_alm( TempMap_Q );
					TempMap_EstNside_U.import_via_alm( TempMap_U );
		
					Mean_Q = TempMap_EstNside_Q.average();
					TempMap_EstNside_Q.add( -Mean_Q );
					Mean_U = TempMap_EstNside_U.average();
					TempMap_EstNside_U.add( -Mean_U );
		
					for( int i=0; i < TempMap_EstNside_Q.Npix(); i++ )
					{
						Dat_work_Q(i,k) = TempMap_EstNside_Q[i];
						Dat_work_U(i,k) = TempMap_EstNside_U[i];
					}
				}
				else
				{
					Mean_Q = TempMap_Q.average();
					TempMap_Q.add( -Mean_Q );
					Mean_U = TempMap_U.average();
					TempMap_U.add( -Mean_U );
		
					for( int i=0; i < TempMap_Q.Npix(); i++ )
					{
						Dat_work_Q(i,k) = TempMap_Q[i];
						Dat_work_U(i,k) = TempMap_U[i];
					}
				}//end if( EstNside != Nside)
			}//end for( int k=0; k < nb_channel; k++ )
		}
		else
		{
			if( Nz != nb_channel )
			{//error number of selected channels different from number of read maps
				fprintf( OUTMAN, "Bad number of selected channels, number different from read maps \n" );
				fprintf( OUTMAN, "Number of selected channels: %d\n", nb_channel );
				fprintf( OUTMAN, "Number of map read: %d\n", Nz );
				usage(argv);
				exit(-1);
			}
			else
			{
				Dat_list_Q = Dat_file_in_Q;
				Dat_list_U = Dat_file_in_U;
			
				for( int k=0; k < nb_channel; k++ )
				{
					for( int i=0; i < Dat_list_Q.nx(); i++ )
					{
						TempMap_Q[i] = Dat_list_Q(i,k);
						TempMap_U[i] = Dat_list_U(i,k);
					}
		
					if( EstNside != Nside)
					{
						TempMap_EstNside_Q.import_via_alm( TempMap_Q );
						TempMap_EstNside_U.import_via_alm( TempMap_U );
		
						Mean_Q = TempMap_EstNside_Q.average();
						TempMap_EstNside_Q.add( -Mean_Q );
						Mean_U = TempMap_EstNside_U.average();
						TempMap_EstNside_U.add( -Mean_U );
		
						for( int i=0; i < TempMap_EstNside_Q.Npix(); i++ )
						{
							Dat_work_Q(i,k) = TempMap_EstNside_Q[i];
							Dat_work_U(i,k) = TempMap_EstNside_U[i];
						}
					}
					else
					{
						Mean_Q = TempMap_Q.average();
						TempMap_Q.add( -Mean_Q );
						Mean_U = TempMap_U.average();
						TempMap_U.add( -Mean_U );
		
						for( int i=0; i < TempMap_Q.Npix(); i++ )
						{
							Dat_work_Q(i,k) = TempMap_Q[i];
							Dat_work_U(i,k) = TempMap_U[i];
						}
					}//end if( EstNside != Nside)
				}//end for( int k=0; k < nb_channel; k++ )
			}//end if Nz != nb_channel
		}
	
		Dat_file_in_Q.free();
		Dat_file_in_U.free();
    }
    else
    {//Read a ASCII file for localization of Healpix files 
    	char List_file_in_Q[12][256];
    	char str_Q[256];
    	int nb_Q = 0;
    	
    	FILE * pFile_Q;
    	
    	char List_file_in_U[12][256];
    	char str_U[256];
    	int nb_U = 0;
    	
    	FILE * pFile_U;
	
		pFile_Q = fopen( Name_Cube_In_Q, "r" );
		
		if( pFile_Q != NULL )
		{
			int count;
			while( !feof( pFile_Q ) )
			{
				count = fscanf( pFile_Q, "%s", str_Q );
				if( count <= 0 )
				{
					fprintf( OUTMAN, "Problem occurs while reading file %s \n", Name_Cube_In_Q );
					exit(-1);
				}
				for( int j=0; j <= strlen(str_Q); j++ )
				{
					List_file_in_Q[nb_Q][j] = str_Q[j];
				}
				nb_Q = nb_Q + 1;
			}

    		fclose( pFile_Q );
		}
		else
		{
			fprintf( OUTMAN, "Could not open input file %s in ASCII mode\n", Name_Cube_In_Q );
			usage(argv);
			exit(-1);
		}
		
		pFile_U = fopen( Name_Cube_In_U, "r" );
		
		if( pFile_U != NULL )
		{
			int count;
			while( !feof( pFile_U ) )
			{
				count = fscanf( pFile_U, "%s", str_U );
				if( count <= 0 )
				{
					fprintf( OUTMAN, "Problem occurs while reading file %s \n", Name_Cube_In_U );
					exit(-1);
				}
				for( int j=0; j <= strlen(str_Q); j++ )
				{
					List_file_in_U[nb_U][j] = str_Q[j];
				}
				nb_U = nb_U + 1;
			}

    		fclose( pFile_U );
		}
		else
		{
			fprintf( OUTMAN, "Could not open input file %s in ASCII mode\n", Name_Cube_In_U );
			usage(argv);
			exit(-1);
		}
		
		if( nb_Q != nb_U )
		{
			fprintf( OUTMAN, "List of maps Q and U don't have the same number of elements\n" );
			fprintf( OUTMAN, "Number of maps read in Q file list: %d \n", nb_Q );
			fprintf( OUTMAN, "Number of maps read in U file list: %d \n", nb_U );
			usage(argv);
			exit(-1);
		}
		
		if( nb_Q == 12 )
		{//12 maps names read selection of maps for GMCA must be done
			//centering data and resizing to new EstNside for GMCA only !!!
			if( Nested == false )
			{
				ReadMap_Q.read( List_file_in_Q[ list_channels_set(0) ], true );// we work in RING
				ReadMap_U.read( List_file_in_U[ list_channels_set(0) ], true );
			}
			else
			{
				ReadMap_Q.read( List_file_in_Q[ list_channels_set(0) ], false );
				ReadMap_U.read( List_file_in_U[ list_channels_set(0) ], false );
				if( ReadMap_Q.Scheme() == DEF_MRS_ORDERING )
				{
					ReadMap_Q.swap_scheme();// we work in NESTED
				}
				if( ReadMap_U.Scheme() == DEF_MRS_ORDERING )
				{
					ReadMap_U.swap_scheme();// we work in NESTED
				}
			}
			if( Nside == 0 )
			{
				Nside = ReadMap_Q.Nside();
				TempMap_Q.alloc( Nside, Nested );
				TempMap_Q = ReadMap_Q;
				TempMap_U.alloc( Nside, Nested );
				TempMap_U = ReadMap_U;
			}
			else
			{
				TempMap_Q.alloc( Nside, Nested );
				TempMap_U.alloc( Nside, Nested );
			
				if( Nside != ReadMap_Q.Nside() )
				{
					cout << "Input Q map " << List_file_in_Q[ list_channels_set(0) ] << " is in nside = " << ReadMap_Q.Nside() << " : resizing needed\n";
					TempMap_Q.import_via_alm( ReadMap_Q );
				}
				else
				{
					TempMap_Q = ReadMap_Q;
				}
				if( Nside != ReadMap_U.Nside() )
				{
					cout << "Input U map " << List_file_in_U[ list_channels_set(0) ] << " is in nside = " << ReadMap_U.Nside() << " : resizing needed\n";
					TempMap_U.import_via_alm( ReadMap_U );
				}
				else
				{
					TempMap_U = ReadMap_U;
				}
			}
			if( EstNside == 0 )
			{
				EstNside = Nside;
			}
		
			if( Verbose == True )
			{
				cout << "Map Q read: " << List_file_in_Q[ list_channels_set(0) ] << '\n';
				cout << "Map U read: " << List_file_in_U[ list_channels_set(0) ] << '\n';
			}
		
			if( EstNside != Nside)
			{
				TempMap_EstNside_Q.alloc( EstNside, Nested );
				TempMap_EstNside_U.alloc( EstNside, Nested );
		
				Dat_work_Q.alloc( TempMap_EstNside_Q.Npix(), nb_channel );
				Dat_work_U.alloc( TempMap_EstNside_U.Npix(), nb_channel );
			}
			else
			{
				Dat_work_Q.alloc( TempMap_Q.Npix(), nb_channel );
				Dat_work_U.alloc( TempMap_U.Npix(), nb_channel );
			}
		
			Dat_list_Q.alloc( TempMap_Q.Npix(), nb_channel );
			Dat_list_U.alloc( TempMap_U.Npix(), nb_channel );
		
			for( int i=0; i < Dat_list_Q.nx(); i++ )
			{
				Dat_list_Q(i,0) = TempMap_Q[i];
				Dat_list_U(i,0) = TempMap_U[i];
			}
			if( EstNside != Nside)
			{
				if( EstNside != ReadMap_Q.Nside() )
				{
					TempMap_EstNside_Q.import_via_alm( ReadMap_Q );
				}
				else
				{
					TempMap_EstNside_Q = ReadMap_Q;
				}
				if( EstNside != ReadMap_U.Nside() )
				{
					TempMap_EstNside_U.import_via_alm( ReadMap_U );
				}
				else
				{
					TempMap_EstNside_U = ReadMap_U;
				}
		
				Mean_Q = TempMap_EstNside_Q.average();
				TempMap_EstNside_Q.add( -Mean_Q );
				Mean_U = TempMap_EstNside_U.average();
				TempMap_EstNside_U.add( -Mean_U );
		
				for( int i=0; i < TempMap_EstNside_Q.Npix(); i++ )
				{
					Dat_work_Q(i,0) = TempMap_EstNside_Q[i];
					Dat_work_U(i,0) = TempMap_EstNside_U[i];
				}
			}
			else
			{
				Mean_Q = TempMap_Q.average();
				TempMap_Q.add( -Mean_Q );
				Mean_U = TempMap_U.average();
				TempMap_U.add( -Mean_U );
		
				for( int i=0; i < TempMap_Q.Npix(); i++ )
				{
					Dat_work_Q(i,0) = TempMap_Q[i];
					Dat_work_U(i,0) = TempMap_U[i];
				}
			}//end if( EstNside != Nside)
								
			for( int k=1; k < nb_channel; k++ )
			{
				if( Nested == false )
				{
					ReadMap_Q.read( List_file_in_Q[ list_channels_set(k) ], true );// we work in RING
					ReadMap_U.read( List_file_in_U[ list_channels_set(k) ], true );
				}
				else
				{
					ReadMap_Q.read( List_file_in_Q[ list_channels_set(k) ], false );
					ReadMap_U.read( List_file_in_U[ list_channels_set(k) ], false );
					if( ReadMap_Q.Scheme() ==  DEF_MRS_ORDERING )
					{
						ReadMap_Q.swap_scheme();// we work in NESTED
					}
					if( ReadMap_U.Scheme() ==  DEF_MRS_ORDERING )
					{
						ReadMap_U.swap_scheme();// we work in NESTED
					}
				}
				if( Nside != ReadMap_Q.Nside() )
				{
					cout << "Input Q map " << List_file_in_Q[ list_channels_set(k) ] << " is in nside = " << ReadMap_Q.Nside() << " : resizing needed\n";
					TempMap_Q.import_via_alm( ReadMap_Q );
				}
				else
				{
					TempMap_Q = ReadMap_Q;
				}
				if( Nside != ReadMap_U.Nside() )
				{
					cout << "Input U map " << List_file_in_U[ list_channels_set(k) ] << " is in nside = " << ReadMap_U.Nside() << " : resizing needed\n";
					TempMap_U.import_via_alm( ReadMap_U );
				}
				else
				{
					TempMap_U = ReadMap_U;
				}
				
				if( Verbose == True )
				{
					cout << "Map Q read: " << List_file_in_Q[ list_channels_set(k) ] << '\n';
					cout << "Map U read: " << List_file_in_U[ list_channels_set(k) ] << '\n';
				}
				
				for( int i=0; i < Dat_list_Q.nx(); i++ )
				{
					Dat_list_Q(i,k) = TempMap_Q[i];
					Dat_list_U(i,k) = TempMap_U[i];
				}
			
				if( EstNside != Nside)
				{
					if( EstNside != ReadMap_Q.Nside() )
					{
						TempMap_EstNside_Q.import_via_alm( ReadMap_Q );
					}
					else
					{
						TempMap_EstNside_Q = ReadMap_Q;
					}
					if( EstNside != ReadMap_U.Nside() )
					{
						TempMap_EstNside_U.import_via_alm( ReadMap_U );
					}
					else
					{
						TempMap_EstNside_U = ReadMap_U;
					}
		
					Mean_Q = TempMap_EstNside_Q.average();
					TempMap_EstNside_Q.add( -Mean_Q );
					Mean_U = TempMap_EstNside_U.average();
					TempMap_EstNside_U.add( -Mean_U );
		
					for( int i=0; i < TempMap_EstNside_Q.Npix(); i++ )
					{
						Dat_work_Q(i,k) = TempMap_EstNside_Q[i];
						Dat_work_U(i,k) = TempMap_EstNside_U[i];
					}
				}
				else
				{
					Mean_Q = TempMap_Q.average();
					TempMap_Q.add( -Mean_Q );
					Mean_U = TempMap_U.average();
					TempMap_U.add( -Mean_U );
		
					for( int i=0; i < TempMap_Q.Npix(); i++ )
					{
						Dat_work_Q(i,k) = TempMap_Q[i];
						Dat_work_U(i,k) = TempMap_U[i];
					}
				}//end if( EstNside != Nside)
			}//end for( int k=0; k < nb_channel; k++ )	
		}
		else
		{
			if( nb_Q != nb_channel )
			{//error number of selected channels different from number of read maps
				fprintf( OUTMAN, "Bad number of selected channels, number different from read maps \n" );
				fprintf( OUTMAN, "Number of selected channels: %d\n", nb_channel );
				fprintf( OUTMAN, "Number of map read: %d\n", nb_Q );
				usage(argv);
				exit(-1);
			}
			else
			{
				//centering data and resizing to new EstNside for GMCA only !!!
				if( Nested == false )
				{
					ReadMap_Q.read( List_file_in_Q[ 0 ], true );// we work in RING
					ReadMap_U.read( List_file_in_U[ 0 ], true );
				}
				else
				{
					ReadMap_Q.read( List_file_in_Q[ 0 ], false );
					ReadMap_U.read( List_file_in_U[ 0 ], false );
					if( ReadMap_Q.Scheme() ==  DEF_MRS_ORDERING )
					{
						ReadMap_Q.swap_scheme();// we work in NESTED
					}
					if( ReadMap_U.Scheme() ==  DEF_MRS_ORDERING )
					{
						ReadMap_U.swap_scheme();// we work in NESTED
					}
				}
				if( Nside == 0 )
				{
					Nside = ReadMap_Q.Nside();
					TempMap_Q.alloc( Nside, Nested );
					TempMap_Q = ReadMap_Q;
					TempMap_U.alloc( Nside, Nested );
					TempMap_U = ReadMap_U;
				}
				else
				{
					TempMap_Q.alloc( Nside, Nested );
					TempMap_U.alloc( Nside, Nested );
			
					if( Nside != ReadMap_Q.Nside() )
					{
						cout << "Input Q map " << List_file_in_Q[ 0 ] << " is in nside = " << ReadMap_Q.Nside() << " : resizing needed\n";
						TempMap_Q.import_via_alm( ReadMap_Q );
					}
					else
					{
						TempMap_Q = ReadMap_Q;
					}
					if( Nside != ReadMap_U.Nside() )
					{
						cout << "Input U map " << List_file_in_U[ 0 ] << " is in nside = " << ReadMap_U.Nside() << " : resizing needed\n";
						TempMap_U.import_via_alm( ReadMap_U );
					}
					else
					{
						TempMap_U = ReadMap_U;
					}
				}
				if( EstNside == 0 )
				{
					EstNside = Nside;
				}
		
				if( Verbose == True )
				{
					cout << "Map Q read: " << List_file_in_Q[ 0 ] << '\n';
					cout << "Map U read: " << List_file_in_U[ 0 ] << '\n';
				}
		
				if( EstNside != Nside)
				{
					TempMap_EstNside_Q.alloc( EstNside, Nested );
					TempMap_EstNside_U.alloc( EstNside, Nested );
		
					Dat_work_Q.alloc( TempMap_EstNside_Q.Npix(), nb_channel );
					Dat_work_U.alloc( TempMap_EstNside_U.Npix(), nb_channel );
				}
				else
				{
					Dat_work_Q.alloc( TempMap_Q.Npix(), nb_channel );
					Dat_work_U.alloc( TempMap_U.Npix(), nb_channel );
				}
		
				Dat_list_Q.alloc( TempMap_Q.Npix(), nb_channel );
				Dat_list_U.alloc( TempMap_U.Npix(), nb_channel );
		
				for( int i=0; i < Dat_list_Q.nx(); i++ )
				{
					Dat_list_Q(i,0) = TempMap_Q[i];
					Dat_list_U(i,0) = TempMap_U[i];
				}
				if( EstNside != Nside)
				{
					if( EstNside != ReadMap_Q.Nside() )
					{
						TempMap_EstNside_Q.import_via_alm( ReadMap_Q );
					}
					else
					{
						TempMap_EstNside_Q = ReadMap_Q;
					}
					if( EstNside != ReadMap_U.Nside() )
					{
						TempMap_EstNside_U.import_via_alm( ReadMap_U );
					}
					else
					{
						TempMap_EstNside_U = ReadMap_U;
					}
		
					Mean_Q = TempMap_EstNside_Q.average();
					TempMap_EstNside_Q.add( -Mean_Q );
					Mean_U = TempMap_EstNside_U.average();
					TempMap_EstNside_U.add( -Mean_U );
		
					for( int i=0; i < TempMap_EstNside_Q.Npix(); i++ )
					{
						Dat_work_Q(i,0) = TempMap_EstNside_Q[i];
						Dat_work_U(i,0) = TempMap_EstNside_U[i];
					}
				}
				else
				{
					Mean_Q = TempMap_Q.average();
					TempMap_Q.add( -Mean_Q );
					Mean_U = TempMap_U.average();
					TempMap_U.add( -Mean_U );
		
					for( int i=0; i < TempMap_Q.Npix(); i++ )
					{
						Dat_work_Q(i,0) = TempMap_Q[i];
						Dat_work_U(i,0) = TempMap_U[i];
					}
				}//end if( EstNside != Nside)
			
				for( int k=1; k < nb_channel; k++ )
				{
					if( Nested == false )
					{
						ReadMap_Q.read( List_file_in_Q[ k ], true );// we work in RING
						ReadMap_U.read( List_file_in_U[ k ], true );
					}
					else
					{
						ReadMap_Q.read( List_file_in_Q[ k ], false );
						ReadMap_U.read( List_file_in_U[ k ], false );
						if( ReadMap_Q.Scheme() ==  DEF_MRS_ORDERING )
						{
							ReadMap_Q.swap_scheme();// we work in NESTED
						}
						if( ReadMap_U.Scheme() ==  DEF_MRS_ORDERING )
						{
							ReadMap_U.swap_scheme();// we work in NESTED
						}
					}
					if( Nside != ReadMap_Q.Nside() )
					{
						cout << "Input Q map " << List_file_in_Q[ k ] << " is in nside = " << ReadMap_Q.Nside() << " : resizing needed\n";
						TempMap_Q.import_via_alm( ReadMap_Q );
					}
					else
					{
						TempMap_Q = ReadMap_Q;
					}
					if( Nside != ReadMap_U.Nside() )
					{
						cout << "Input U map " << List_file_in_U[ k ] << " is in nside = " << ReadMap_U.Nside() << " : resizing needed\n";
						TempMap_U.import_via_alm( ReadMap_U );
					}
					else
					{
						TempMap_U = ReadMap_U;
					}
					
					if( Verbose == True )
					{
						cout << "Map Q read: " << List_file_in_Q[ k ] << '\n';
						cout << "Map U read: " << List_file_in_U[ k ] << '\n';
					}
					
					for( int i=0; i < Dat_list_Q.nx(); i++ )
					{
						Dat_list_Q(i,k) = TempMap_Q[i];
						Dat_list_U(i,k) = TempMap_U[i];
					}
		
					if( EstNside != Nside)
					{
						if( EstNside != ReadMap_Q.Nside() )
						{
							TempMap_EstNside_Q.import_via_alm( ReadMap_Q );
						}
						else
						{
							TempMap_EstNside_Q = ReadMap_Q;
						}
						if( EstNside != ReadMap_U.Nside() )
						{
							TempMap_EstNside_U.import_via_alm( ReadMap_U );
						}
						else
						{
							TempMap_EstNside_U = ReadMap_U;
						}
		
						Mean_Q = TempMap_EstNside_Q.average();
						TempMap_EstNside_Q.add( -Mean_Q );
						Mean_U = TempMap_EstNside_U.average();
						TempMap_EstNside_U.add( -Mean_U );
		
						for( int i=0; i < TempMap_EstNside_Q.Npix(); i++ )
						{
							Dat_work_Q(i,k) = TempMap_EstNside_Q[i];
							Dat_work_U(i,k) = TempMap_EstNside_U[i];
						}
					}
					else
					{
						Mean_Q = TempMap_Q.average();
						TempMap_Q.add( -Mean_Q );
						Mean_U = TempMap_U.average();
						TempMap_U.add( -Mean_U );
		
						for( int i=0; i < TempMap_Q.Npix(); i++ )
						{
							Dat_work_Q(i,k) = TempMap_Q[i];
							Dat_work_U(i,k) = TempMap_U[i];
						}
					}//end if( EstNside != Nside)
				}//end for( int k=1; k < nb_channel; k++ )
			}//end if nb != nb_channel
		}
    }// end if ASCII == false
	
	if (Verbose == True) cout << "Nside = " << Nside << " EstNside = " << EstNside << " NbrChannels = " << nb_channel << endl;
	
	NbrScale2d = log( (double) EstNside )/log( (double) 2.0 ) - 4;
	cout << "GMCA - NbrScale : " << NbrScale2d << '\n';
	
	if( conv_EB == true )
	{// Convert Q and U data table in E and B ones. ONLY for datas, NOT for std
		Hdmap Map_work_Q, Map_work_U, Map_list_Q, Map_list_U;
		Map_work_Q.alloc( EstNside, Nested );
		Map_work_U.alloc( EstNside, Nested );
		Map_list_Q.alloc( Nside, Nested );
		Map_list_U.alloc( Nside, Nested );
		
		int Lmax_EB_work = 3*EstNside;
		int Mmax_EB_work = Lmax_EB_work;
		int Lmax_EB_list = 3*Nside;
		int Mmax_EB_list = Lmax_EB_list;
		
		arr<double> weight_T_work, weight_T_list;
		weight_T_work.alloc( 2*EstNside );
		weight_T_list.alloc( 2*Nside );
		weight_T_work.fill(1);
		weight_T_list.fill(1);
		
		Alm<xcomplex<double> > ALM_E_work(Lmax_EB_work,Mmax_EB_work), ALM_B_work(Lmax_EB_work,Mmax_EB_work), ALM_E_list(Lmax_EB_list,Mmax_EB_list), ALM_B_list(Lmax_EB_list,Mmax_EB_list);
		
		for( int k = 0; k < nb_channel; k++ )
		{
			for( int i = 0; i < Map_work_Q.Npix(); i++ )
			{
				Map_work_Q[i] = Dat_work_Q(i,k);
				Map_work_U[i] = Dat_work_U(i,k);
			}
			for( int i = 0; i < Map_list_Q.Npix(); i++ )
			{
				Map_list_Q[i] = Dat_list_Q(i,k);
				Map_list_U[i] = Dat_list_U(i,k);
			}
			
			if( Nested == true )
		    {
		    	Map_work_Q.swap_scheme();
    			Map_work_U.swap_scheme();
    			Map_list_Q.swap_scheme();
    			Map_list_U.swap_scheme();
		    }//ALM Trans Must be done in RING scheme
		    
		    map2alm_pol_iter_QU( Map_work_Q, Map_work_U, ALM_E_work, ALM_B_work, 0, weight_T_work );
				
		    alm2map( ALM_E_work, Map_work_Q );
		    alm2map( ALM_B_work, Map_work_U );
		    
		    map2alm_pol_iter_QU( Map_list_Q, Map_list_U, ALM_E_list, ALM_B_list, 0, weight_T_list );
				
		    alm2map( ALM_E_list, Map_list_Q );
		    alm2map( ALM_B_list, Map_list_U );
		    
		    if( Nested == true )
		    {
		    	Map_work_Q.swap_scheme();
    			Map_work_U.swap_scheme();
    			Map_list_Q.swap_scheme();
    			Map_list_U.swap_scheme();
		    }//ALM Trans Must be done in RING scheme
		    
		    for( int i = 0; i < Map_work_Q.Npix(); i++ )
			{
				Dat_work_Q(i,k) = Map_work_Q[i];
				Dat_work_U(i,k) = Map_work_U[i];
			}
			for( int i = 0; i < Map_list_Q.Npix(); i++ )
			{
				Dat_list_Q(i,k) = Map_list_Q[i];
				Dat_list_U(i,k) = Map_list_U[i];
			}
		}
	}
	  
	MRS_GMCA_POLA WT;
       
	if (Verbose == True) cout << "Alloc ...  " << endl;
	WT.alloc( EstNside, NbrScale2d, Nested ); // La Transformation 1D se fait ensuite !!
 
	WT.NbrCannels = nb_channel;
	if (NbrSources == 0) NbrSources = nb_channel;
	WT.NbrSources = NbrSources; 
	WT.Max_GMCA_Iter = Max_GMCA_Iter;
	WT.TypeThresholdDecrease = TypeThresholdDecrease;
	WT.Inpainting = False;
	WT.PositiveMatrix = False;
	WT.PositiveSource = False;
	WT.L1_Matrix_Constraint = False;
	WT.KMin = KMin;
	WT.DisjSpec = False;
	WT.OrthoSpectra = False;
	WT.GlobThrd = False;
	WT.BandThrd = false;//true pour mrs_allbandcleaning_perband
		
	int status;
	status = WT.set_planck( FixCMB, NbrScale2d, NbrSources, list_channels_set );
	if( status == -1 )
	{
		fprintf(OUTMAN, "Error: selected number of sources lower than number of known sources setted...\n");
		exit(-1);
	}
       
	fltarray TransDat;
	
	WT.transform_sources( Dat_work_Q, Dat_work_U, TransDat, false );
          
	// run the GMCA method
	if (Verbose == True) cout << "Running GMCA ... "<< endl;
	   
	fltarray TransTabSource;
	int NbrCoef = TransDat.nx();
	TransTabSource.alloc( NbrCoef, NbrSources );
	
	WT.GMCA::Verbose = Verbose;
	WT.run_gmca( TransDat, TransTabSource );

	// Reconstruction :
	if (Verbose == True) cout << "Reconstruction ... "<< endl;
	fltarray EstSources_Q, EstSources_U;
	WT.recons_sources( Dat_list_Q, EstSources_Q );
	WT.recons_sources( Dat_list_U, EstSources_U );

	if (Verbose == True) cout << "Write results ... "<< endl;

    char FN[512];
    sprintf( FN, "%s_Q_sources.fits", Name_Out );
    fits_write_fltarr( FN, EstSources_Q );
    
    sprintf( FN, "%s_U_sources.fits", Name_Out );
    fits_write_fltarr( FN, EstSources_U );
    
    sprintf( FN, "%s_mat.fits", Name_Out );
    fits_write_fltarr( FN, WT.MixingMat );
    
    sprintf( FN, "%s_invmat.fits", Name_Out );
    fits_write_fltarr( FN, WT.InvMixingMat );
    
    if( estim_bias == true )
	{
		dblarray MapBeam, Map_Beam_read, Bias;
		fits_read_dblarr( BeamChannels_List, Map_Beam_read );
		
		if( Map_Beam_read.nx() == 12 )
		{
			MapBeam.alloc( nb_channel, Map_Beam_read.ny() );
			
			for( int i=0; i < nb_channel; i++ )
			{
				for( int j=0; j < Map_Beam_read.ny(); j++ )
				{
					MapBeam( i, j ) = Map_Beam_read( list_channels_set(i), j );
				}
			}
		}
		else
		{
			if( Map_Beam_read.nx() != nb_channel )
			{//error number of selected channels different from number of read maps
				fprintf( OUTMAN, "Bad number of selected channels, number different from channel beams read\n" );
				fprintf( OUTMAN, "Number of selected channels: %d\n", nb_channel );
				fprintf( OUTMAN, "Number of channel beams read: %d\n", Map_Beam_read.nx() );
				usage(argv);
				exit(-1);
			}
			else
			{
				MapBeam = Map_Beam_read;
			}
		}
		Map_Beam_read.free();
		
		status = WT.planck_bias_estim( MapBeam, Bias );
		
		/*
		if( list_channels_set.min() >= 7 )
		{// WMAP only
			status = WT.planck_bias_estim( MapBeam, Bias, true );
		}
		else
		{
			status = WT.planck_bias_estim( MapBeam, Bias, false );
		}
		*/
		
		if( status == -1 )
		{
			fprintf(OUTMAN, "Error: option planck has not been setted for GMCA\n");
			exit(-1);
		}
		else
		{
			sprintf( FN, "%s_bias.fits", Name_Out );
    		fits_write_dblarr( FN, Bias );
		}
	}
	
	if( estim_std_src == true )
	{
		//Planck channels selected?
		bool planck_channel_set = false;
		//WMAP channels selected?
		bool wmap_channel_set = false;
		for( int i=0; i < nb_channel; i++ )
		{
			if( list_channels_set(i) < 7 )
			{
				planck_channel_set = true;
			}
			if( list_channels_set(i) >= 7 )
			{
				wmap_channel_set = true;
			}
		}
		
		if( ( planck_channel_set == true )&&( wmap_channel_set == true)&&( Nside > 512 ) )
		{
			fprintf(OUTMAN, "Warning: Channels from both Planck and WMAP selected with Nside parameter greater than 512: Nside = %d \n.", Nside );
			fprintf(OUTMAN, "impossible to estimate the standard deviation on the sources.\n");
		}
		else
		{
			fltarray std_file_in_Q, std_work_Q, std_estim_Q, std_file_in_U, std_work_U, std_estim_U;
			std_work_Q.alloc( 12*Nside*Nside, nb_channel );
			std_work_U.alloc( 12*Nside*Nside, nb_channel );
			//std_estim.alloc( 12*Nside*Nside, NbrSources );
			
			float Rfact;
			
			if( ASCII == false )
			{
				fits_read_fltarr( STD_Cube_In_Q, std_file_in_Q );
				fits_read_fltarr( STD_Cube_In_U, std_file_in_U );
				
				if( (std_file_in_Q.nx() != std_file_in_U.nx()) || (std_file_in_Q.ny() != std_file_in_U.ny()) )
				{
					fprintf( OUTMAN, "Error STD table Q and STD table U don't have the same size:\n" );
					std_file_in_Q.info("Table Q");
					std_file_in_U.info("Table U");
					usage(argv);
					exit(-1);
				}
		
				int Nside_std = sqrt( std_file_in_Q.nx() / 12 );
		
				int Nz_std = std_file_in_Q.ny();
				
				if( Nside_std != Nside)
				{
					TempMap_EstNside_Q.alloc( Nside, Nested );
					TempMap_EstNside_U.alloc( Nside, Nested );
				}
				
				TempMap_Q.alloc( Nside_std, Nested );
				TempMap_U.alloc( Nside_std, Nested );	
		
				if( Nz_std == 12 )
				{//12 std maps read selection of maps must be done
					
					for( int k=0; k < nb_channel; k++ )
					{
						for( int i=0; i < std_file_in_Q.nx(); i++ )
						{
							TempMap_Q[i] = std_file_in_Q( i, list_channels_set(k) );
							TempMap_U[i] = std_file_in_U( i, list_channels_set(k) );
						}
						
						if( list_channels_set(k) <= 6 )
						{
							Rfact = 2048 / ( (float) Nside );//Rfact = 2048.0*2048.0 / ( (float) Nside*Nside );
						}
						else
						{
							Rfact = 512 / ( (float) Nside );//Rfact = 512.0*512.0 / ( (float) Nside*Nside );
						}
			
						if( Nside_std != Nside)
						{
							TempMap_EstNside_Q.import_via_alm( TempMap_Q );
							TempMap_EstNside_U.import_via_alm( TempMap_U );
			
							for( int i=0; i < TempMap_EstNside_Q.Npix(); i++ )
							{
								std_work_Q(i,k) = TempMap_EstNside_Q[i]/Rfact;
								std_work_U(i,k) = TempMap_EstNside_U[i]/Rfact;
							}
						}
						else
						{	
							for( int i=0; i < TempMap_Q.Npix(); i++ )
							{
								std_work_Q(i,k) = TempMap_Q[i]/Rfact;
								std_work_U(i,k) = TempMap_U[i]/Rfact;
							}
						}//end if( Nside_std != Nside)
					}//end for( int k=0; k < nb_channel; k++ )
				}
				else
				{
					if( Nz_std != nb_channel )
					{//error number of selected channels different from number of read maps
						fprintf( OUTMAN, "Bad number of selected channels, number different from read std maps \n" );
						fprintf( OUTMAN, "Number of selected channels: %d\n", nb_channel );
						fprintf( OUTMAN, "Number of map read: %d\n", Nz_std );
						usage(argv);
						exit(-1);
					}
					else
					{			
						for( int k=0; k < nb_channel; k++ )
						{
							for( int i=0; i < std_file_in_Q.nx(); i++ )
							{
								TempMap_Q[i] = std_file_in_Q(i,k);
								TempMap_U[i] = std_file_in_U(i,k);
							}
							
							if( list_channels_set(k) <= 6 )
							{
								Rfact = 2048 / ( (float) Nside );//Rfact = 2048.0*2048.0 / ( (float) Nside*Nside );
							}
							else
							{
								Rfact = 512 / ( (float) Nside );//Rfact = 512.0*512.0 / ( (float) Nside*Nside );
							}
		
							if( Nside_std != Nside)
							{
								TempMap_EstNside_Q.import_via_alm( TempMap_Q );
								TempMap_EstNside_U.import_via_alm( TempMap_U );
			
								for( int i=0; i < TempMap_EstNside_Q.Npix(); i++ )
								{
									std_work_Q(i,k) = TempMap_EstNside_Q[i]/Rfact;
									std_work_U(i,k) = TempMap_EstNside_U[i]/Rfact;
								}
							}
							else
							{	
								for( int i=0; i < TempMap_Q.Npix(); i++ )
								{
									std_work_Q(i,k) = TempMap_Q[i]/Rfact;
									std_work_U(i,k) = TempMap_U[i]/Rfact;
								}
							}//end if( Nside_std != Nside)
						}//end for( int k=0; k < nb_channel; k++ )
					}//end if Nz != nb_channel
				}
	
				std_file_in_Q.free();
				std_file_in_U.free();				
			}
			else
			{
				//Read a ASCII file for localization of Healpix standard deviation files 
    			char List_STDfile_in_Q[12][256];
    			char str2_Q[256];
    			int nb_2_Q = 0;
    	
    			FILE * pFile_std_Q;
    			
    			char List_STDfile_in_U[12][256];
    			char str2_U[256];
    			int nb_2_U = 0;
    	
    			FILE * pFile_std_U;
	
				pFile_std_Q = fopen( STD_Cube_In_Q, "r" );
		
				if( pFile_std_Q != NULL )
				{
					int count_std;
					while( !feof( pFile_std_Q ) )
					{
						count_std = fscanf( pFile_std_Q, "%s", str2_Q );
						if( count_std <= 0 )
						{
							fprintf( OUTMAN, "Problem occurs while reading file %s \n", STD_Cube_In_Q );
							exit(-1);
						}
						for( int j=0; j <= strlen(str2_Q); j++ )
						{
							List_STDfile_in_Q[nb_2_Q][j] = str2_Q[j];
						}
						nb_2_Q = nb_2_Q + 1;
					}

    				fclose( pFile_std_Q );
				}
				else
				{
					fprintf( OUTMAN, "Could not open input file %s in ASCII mode\n", STD_Cube_In_Q );
					usage(argv);
					exit(-1);
				}
				
				pFile_std_U = fopen( STD_Cube_In_U, "r" );
		
				if( pFile_std_U != NULL )
				{
					int count_std;
					while( !feof( pFile_std_U ) )
					{
						count_std = fscanf( pFile_std_U, "%s", str2_U );
						if( count_std <= 0 )
						{
							fprintf( OUTMAN, "Problem occurs while reading file %s \n", STD_Cube_In_U );
							exit(-1);
						}
						for( int j=0; j <= strlen(str2_U); j++ )
						{
							List_STDfile_in_U[nb_2_U][j] = str2_U[j];
						}
						nb_2_U = nb_2_U + 1;
					}

    				fclose( pFile_std_U );
				}
				else
				{
					fprintf( OUTMAN, "Could not open input file %s in ASCII mode\n", STD_Cube_In_U );
					usage(argv);
					exit(-1);
				}
				
				if( nb_2_Q != nb_2_U )
				{
					fprintf( OUTMAN, "List of STD maps Q and U don't have the same number of elements\n" );
					fprintf( OUTMAN, "Number of STD maps read in Q file list: %d \n", nb_2_Q );
					fprintf( OUTMAN, "Number of STD maps read in U file list: %d \n", nb_2_U );
					usage(argv);
					exit(-1);
				}
				
				if( nb_2_Q == 12 )
				{//12 std maps names read selection of maps must be done
					if( Nested == false )
					{
						ReadMap_Q.read( List_STDfile_in_Q[ list_channels_set(0) ], true );// we work in RING
						ReadMap_U.read( List_STDfile_in_U[ list_channels_set(0) ], true );
					}
					else
					{
						ReadMap_Q.read( List_STDfile_in_Q[ list_channels_set(0) ], false );
						ReadMap_U.read( List_STDfile_in_U[ list_channels_set(0) ], false );
						if( ReadMap_Q.Scheme() ==  DEF_MRS_ORDERING )
						{
							ReadMap_Q.swap_scheme();// we work in NESTED
						}
						if( ReadMap_U.Scheme() ==  DEF_MRS_ORDERING )
						{
							ReadMap_U.swap_scheme();// we work in NESTED
						}
					}
					TempMap_Q.alloc( Nside, Nested );
					TempMap_U.alloc( Nside, Nested );
			
					if( Nside != ReadMap_Q.Nside() )
					{
						cout << "Input STD Q map " << List_STDfile_in_Q[ list_channels_set(0) ] << " is in nside = " << ReadMap_Q.Nside() << " : resizing needed\n";
						TempMap_Q.import_via_alm( ReadMap_Q );
					}
					else
					{
						TempMap_Q = ReadMap_Q;
					}
					if( Nside != ReadMap_U.Nside() )
					{
						cout << "Input STD U map " << List_STDfile_in_U[ list_channels_set(0) ] << " is in nside = " << ReadMap_U.Nside() << " : resizing needed\n";
						TempMap_U.import_via_alm( ReadMap_U );
					}
					else
					{
						TempMap_U = ReadMap_U;
					}
				
					if( Verbose == True )
					{
						cout << "STD Q map read: " << List_STDfile_in_Q[ list_channels_set(0) ] << '\n';
						cout << "STD U map read: " << List_STDfile_in_U[ list_channels_set(0) ] << '\n';
					}
				
					if( list_channels_set(0) <= 6 )
					{
						Rfact = 2048 / ( (float) Nside );//Rfact = 2048.0*2048.0 / ( (float) Nside*Nside );
					}
					else
					{
						Rfact = 512 / ( (float) Nside );//Rfact = 512.0*512.0 / ( (float) Nside*Nside );
					}
				
					for( int i=0; i < std_work_Q.nx(); i++ )
					{
						std_work_Q(i,0) = TempMap_Q[i]/Rfact;
						std_work_U(i,0) = TempMap_U[i]/Rfact;
					}
							
					for( int k=1; k < nb_channel; k++ )
					{
						if( Nested == false )
						{
							ReadMap_Q.read( List_STDfile_in_Q[ list_channels_set(k) ], true );// we work in RING
							ReadMap_U.read( List_STDfile_in_U[ list_channels_set(k) ], true );
						}
						else
						{
							ReadMap_Q.read( List_STDfile_in_Q[ list_channels_set(k) ], false );
							ReadMap_U.read( List_STDfile_in_U[ list_channels_set(k) ], false );
							if( ReadMap_Q.Scheme() ==  DEF_MRS_ORDERING )
							{
								ReadMap_Q.swap_scheme();// we work in NESTED
							}
							if( ReadMap_U.Scheme() ==  DEF_MRS_ORDERING )
							{
								ReadMap_U.swap_scheme();// we work in NESTED
							}
						}
						
						if( Nside != ReadMap_Q.Nside() )
						{
							cout << "Input STD Q map " << List_STDfile_in_Q[ list_channels_set(k) ] << " is in nside = " << ReadMap_Q.Nside() << " : resizing needed\n";
							TempMap_Q.import_via_alm( ReadMap_Q );
						}
						else
						{
							TempMap_Q = ReadMap_Q;
						}
						if( Nside != ReadMap_U.Nside() )
						{
							cout << "Input STD U map " << List_STDfile_in_U[ list_channels_set(k) ] << " is in nside = " << ReadMap_U.Nside() << " : resizing needed\n";
							TempMap_U.import_via_alm( ReadMap_U );
						}
						else
						{
							TempMap_U = ReadMap_U;
						}
				
						if( Verbose == True )
						{
							cout << "STD Q map read: " << List_STDfile_in_Q[ list_channels_set(k) ] << '\n';
							cout << "STD U map read: " << List_STDfile_in_U[ list_channels_set(k) ] << '\n';
						}
						
						if( list_channels_set(k) <= 6 )
						{
							Rfact = 2048 / ( (float) Nside );//Rfact = 2048.0*2048.0 / ( (float) Nside*Nside );
						}
						else
						{
							Rfact = 512 / ( (float) Nside );//Rfact = 512.0*512.0 / ( (float) Nside*Nside );
						}
				
						for( int i=0; i < std_work_Q.nx(); i++ )
						{
							std_work_Q(i,k) = TempMap_Q[i]/Rfact;
							std_work_U(i,k) = TempMap_U[i]/Rfact;
						}
					}
				}
				else
				{
					if( nb_2_Q != nb_channel )
					{//error number of selected channels different from number of read maps
						fprintf( OUTMAN, "Bad number of selected channels, number different from read std maps \n" );
						fprintf( OUTMAN, "Number of selected channels: %d\n", nb_channel );
						fprintf( OUTMAN, "Number of map read: %d\n", nb_2_Q );
						usage(argv);
						exit(-1);
					}
					else
					{
						if( Nested == false )
						{
							ReadMap_Q.read( List_STDfile_in_Q[ 0 ], true );// we work in RING
							ReadMap_U.read( List_STDfile_in_U[ 0 ], true );
						}
						else
						{
							ReadMap_Q.read( List_STDfile_in_Q[ 0 ], false );
							ReadMap_U.read( List_STDfile_in_U[ 0 ], false );
							if( ReadMap_Q.Scheme() ==  DEF_MRS_ORDERING )
							{
								ReadMap_Q.swap_scheme();// we work in NESTED
							}
							if( ReadMap_U.Scheme() ==  DEF_MRS_ORDERING )
							{
								ReadMap_U.swap_scheme();// we work in NESTED
							}
						}
						TempMap_Q.alloc( Nside, Nested );
						TempMap_U.alloc( Nside, Nested );
			
						if( Nside != ReadMap_Q.Nside() )
						{
							cout << "Input STD Q map " << List_STDfile_in_Q[ 0 ] << " is in nside = " << ReadMap_Q.Nside() << " : resizing needed\n";
							TempMap_Q.import_via_alm( ReadMap_Q );
						}
						else
						{
							TempMap_Q = ReadMap_Q;
						}
						if( Nside != ReadMap_U.Nside() )
						{
							cout << "Input STD U map " << List_STDfile_in_U[ 0 ] << " is in nside = " << ReadMap_U.Nside() << " : resizing needed\n";
							TempMap_U.import_via_alm( ReadMap_U );
						}
						else
						{
							TempMap_U = ReadMap_U;
						}
				
						if( Verbose == True )
						{
							cout << "STD Q map read: " << List_STDfile_in_Q[ 0 ] << '\n';
							cout << "STD U map read: " << List_STDfile_in_U[ 0 ] << '\n';
						}
				
						if( list_channels_set(0) <= 6 )
						{
							Rfact = 2048 / ( (float) Nside );//Rfact = 2048.0*2048.0 / ( (float) Nside*Nside );
						}
						else
						{
							Rfact = 512 / ( (float) Nside );//Rfact = 512.0*512.0 / ( (float) Nside*Nside );
						}
				
						for( int i=0; i < std_work_Q.nx(); i++ )
						{
							std_work_Q(i,0) = TempMap_Q[i]/Rfact;
							std_work_U(i,0) = TempMap_U[i]/Rfact;
						}
						
						for( int k=1; k < nb_channel; k++ )
						{
							if( Nested == false )
							{
								ReadMap_Q.read( List_STDfile_in_Q[ k ], true );// we work in RING
								ReadMap_U.read( List_STDfile_in_U[ k ], true );
							}
							else
							{
								ReadMap_Q.read( List_STDfile_in_Q[ k ], false );
								ReadMap_U.read( List_STDfile_in_U[ k ], false );
								if( ReadMap_Q.Scheme() ==  DEF_MRS_ORDERING )
								{
									ReadMap_Q.swap_scheme();// we work in NESTED
								}
								if( ReadMap_U.Scheme() ==  DEF_MRS_ORDERING )
								{
									ReadMap_U.swap_scheme();// we work in NESTED
								}
							}
							if( Nside != ReadMap_Q.Nside() )
							{
								cout << "Input STD Q map " << List_STDfile_in_Q[ k ] << " is in nside = " << ReadMap_Q.Nside() << " : resizing needed\n";
								TempMap_Q.import_via_alm( ReadMap_Q );
							}
							else
							{
								TempMap_Q = ReadMap_Q;
							}
							if( Nside != ReadMap_U.Nside() )
							{
								cout << "Input STD U map " << List_STDfile_in_U[ k ] << " is in nside = " << ReadMap_U.Nside() << " : resizing needed\n";
								TempMap_U.import_via_alm( ReadMap_U );
							}
							else
							{
								TempMap_U = ReadMap_U;
							}
					
							if( Verbose == True )
							{
								cout << "STD Q map read: " << List_STDfile_in_Q[ k ] << '\n';
								cout << "STD U map read: " << List_STDfile_in_U[ k ] << '\n';
							}
							
							
							if( list_channels_set(k) <= 6 )
							{
								Rfact = 2048 / ( (float) Nside );//Rfact = 2048.0*2048.0 / ( (float) Nside*Nside );
							}
							else
							{
								Rfact = 512 / ( (float) Nside );//Rfact = 512.0*512.0 / ( (float) Nside*Nside );
							}
					
							for( int i=0; i < std_work_Q.nx(); i++ )
							{
								std_work_Q(i,k) = TempMap_Q[i]/Rfact;
								std_work_U(i,k) = TempMap_U[i]/Rfact;
							}
						}//end for( int k=1; k < nb_channel; k++ )
					}//end if nb_2 != nb_channel
				}// end if( nb_2 == 14 )
			}// end if ASCII
			
			fltarray var_in_Q, var_out_Q, var_in_U, var_out_U, pinv_2;
			MatOper MAT;
			
			var_in_Q = std_work_Q;
			var_in_U = std_work_U;
			for(long i=0; i < var_in_Q.n_elem(); i++ )
			{
				var_in_Q(i) = (float) pow( (double) var_in_Q(i), 2 );
				var_in_U(i) = (float) pow( (double) var_in_U(i), 2 );
			}
			//var_in^2;
			pinv_2 = WT.InvMixingMat;
			for(long i=0; i < pinv_2.n_elem(); i++ )
			{
				pinv_2(i) = (float) pow( (double) pinv_2(i), 2 );
			}
			//pinv_2 = WT.InvMixingMat^2;
			
			fltarray Npinv_2=pinv_2;
			MAT.transpose(Npinv_2,pinv_2);
			MAT.mat_mult( pinv_2, var_in_Q, var_out_Q );
			MAT.mat_mult( pinv_2, var_in_U, var_out_U );
			
			std_estim_Q = var_out_Q;
			std_estim_U = var_out_U;
			for(long i=0; i < std_estim_Q.n_elem(); i++ )
			{
				std_estim_Q(i) = sqrt( std_estim_Q(i) );
				std_estim_U(i) = sqrt( std_estim_U(i) );
			}
			//std_estim = var_out^0.5;
			
			sprintf( FN, "%s_estim_src_Q_std.fits", Name_Out );
    		fits_write_fltarr( FN, std_estim_Q );
    		
    		sprintf( FN, "%s_estim_src_U_std.fits", Name_Out );
    		fits_write_fltarr( FN, std_estim_U );
		}// end if( ( planck_channel_set == true )&&( wmap_channel_set == true) )
	}// end if( estim_std_src == true )

    exit(0);
}

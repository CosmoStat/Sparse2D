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
**    Date:  09/12/09
**    
**    File:  mrs_planck_gmca.cc
**
*******************************************************************************
**
**    DESCRIPTION  Generalized Morphological Component Analysis for planck datas
**    ----------- 
**                 
**    Usage: mrs_planck_gmca options cube output
**
******************************************************************************/

#include "HealpixClass.h"
#include "GMCA.h"
#include "MRS_Sparse.h"

/****************************************************************************/
/****************************************************************************/
#define DEF_HASLAM_NAME "408MHz_512.fits"

char Name_Cube_In[256];
char Name_Out[256];
char Haslam_FileName[256];
char BeamChannels_List[256];
char STD_Cube_In[256];

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
bool FixSZ = false;
bool FixFreeFree = false;
bool UseSync = false;

char list_channels_in[15];

bool WMAP = false;
bool HFI = false;
bool LFI = false;

bool ASCII = false;

bool estim_bias = false;
bool estim_std_src = false;

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
    fprintf(OUTMAN, "Usage: %s options input_data output_source [output_channel] \n\n", argv[0]);
    fprintf(OUTMAN, "   where options =  \n");
    
    fprintf(OUTMAN, "         [-i NbrIter]\n");
    fprintf(OUTMAN, "             Number of iterations. Default is no %d. \n", Max_GMCA_Iter);
   
    fprintf(OUTMAN, "         [-N Number_of_Sources]\n");
    fprintf(OUTMAN, "             Number of sources. \n");    
        
    fprintf(OUTMAN, "         -b Mixing matrix CMB column fixed \n");
    
    fprintf(OUTMAN, "         -s Mixing matrix SZ column fixed \n");
    
    fprintf(OUTMAN, "         -f Mixing matrix Free-Free column fixed \n");
    
    fprintf(OUTMAN, "         -S Synchrotron spectrum estimated separately \n");
    
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
    fprintf(OUTMAN, "             list_channel_in string with 1 character for each channel selected:\n");
    fprintf(OUTMAN, "             '1', '2', '3', '4', '5', '6', '7', '8', '9' for Planck channels\n");
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
	while ((c = GetOpt(argc,argv,"i:n:N:e:K:B:V:C:WLHTvbsfSzZ:")) != -1) 
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
 			
 			case 'b':
 				FixCMB = true;
 				//OptInd = OptInd - 1;
 				break;
				
			case 's':
 				FixSZ = true;
 				//OptInd = OptInd - 1;
 				break;
 			
 			case 'f':
 				FixFreeFree = true;
 				//OptInd = OptInd - 1;
 				break;
 			
 			case 'S':
 				UseSync = true;
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
				strcpy( STD_Cube_In, argv[OptInd-1] );
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
		strcpy( Name_Cube_In, argv[OptInd++] );
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
	fltarray Dat_file_in;
    /* Get command line arguments, open input file(s) if necessary */
	fitsstruct Header;
	char Cmd[512];
	Cmd[0] = '\0';
	for (int k =0; k < argc; k++) sprintf(Cmd, "%s %s", Cmd, argv[k]);
	
	cout << "argc: " << argc << '\n';
   
    transinit(argc, argv);
	
	//Fix selected channels
	intarray set_planck_channel(14);
	set_planck_channel.init(0);
	
	for( int k=0; k < 15; k++ )
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
			
			case '8':
				set_planck_channel(7) = 1;
				break;
			
			case '9':
				set_planck_channel(8) = 1;
				break;
			
			case 'K':
				set_planck_channel(9) = 1;// K WMAP Band
				break;
			
			case 'A':
				set_planck_channel(10) = 1;// Ka WMAP Band
				break;
			
			case 'Q':
				set_planck_channel(11) = 1;// Q WMAP Band
				break;
			
			case 'V':
				set_planck_channel(12) = 1;// V WMAP Band
				break;
			
			case 'W':
				set_planck_channel(13) = 1;// W WMAP Band
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
		set_planck_channel(7) = 1;
		set_planck_channel(8) = 1;
	}
	if( WMAP == true )
	{
		set_planck_channel(9) = 1;
		set_planck_channel(10) = 1;
		set_planck_channel(11) = 1;
		set_planck_channel(12) = 1;
		set_planck_channel(13) = 1;
	}
	//set_planck_channel.display( 14 );
	
	intarray list_channels_set;
	list_channels_set.alloc( set_planck_channel.total() );
	int nb_channel = 0;//final: number of selected channels
	int nb_sync_channel = 0;//final: number of channels where synchrotron could be estimated
	for( int n=0; n < 14; n++ )
	{
		if( set_planck_channel(n) == 1 )
		{
			list_channels_set( nb_channel ) = n;
			nb_channel = nb_channel + 1;
			
			if( n <= 2 )
			{
				nb_sync_channel = nb_sync_channel + 1;
			}
			if( (n >= 9) && (n <= 12) )
			{
				nb_sync_channel = nb_sync_channel + 1;
			}
		}
	}
	
	if( nb_channel == 0 )
	{
		fprintf( OUTMAN, "No Channels were selected" );
		usage(argv);
		exit(-1);
	}
	
	if( (nb_sync_channel <= 2) && (UseSync == true) )
	{
		cout << "Not enough channels for synchrotron estimation, disabling option -S\n";
		UseSync = false;
	}
	
	if (Verbose == True)
    {
		cout << "Filename in = " << Name_Cube_In << endl;
		cout << "Filename out = " << Name_Out  << endl;
		cout << "NbrGMCA_Iter = " <<  Max_GMCA_Iter << endl;
		cout << "Mixing matrix CMB column fixed: " << FixCMB << '\n';
		cout << "Mixing matrix SZ column fixed: " << FixSZ << '\n';
		cout << "Mixing matrix FreeFree column fixed: " << FixFreeFree << '\n';
		cout << "Synchrotron spectrum estimated separately: " << UseSync << '\n';
		list_channels_set.display( 14 );
		if( estim_bias == true )
		{
			cout << "Channels beams File: " << BeamChannels_List << '\n';
		}
		if( estim_std_src == true )
		{
			cout << "Channels maps of standard deviations: " << STD_Cube_In << '\n';
		}
	}
	
	if (Verbose == True) cout << "\n Reading the data"<< endl;
    
    fltarray Dat_work, Dat_list;
    float Mean;
	Hdmap ReadMap, TempMap, TempMap_EstNside;
    		
    if( ASCII == false )
    {     
		fits_read_fltarr( Name_Cube_In, Dat_file_in );
		Dat_file_in.info();
		
		Nside = sqrt( Dat_file_in.nx() / 12 );
		if( EstNside == 0 )
		{
			EstNside = Nside;
		}
		int Nz = Dat_file_in.ny();
		
		//centering data and resizing to new EstNside for GMCA only !!!
		TempMap.alloc( Nside, Nested );
		
		if( EstNside != Nside)
		{
			TempMap_EstNside.alloc( EstNside, Nested );
		
			Dat_work.alloc( TempMap_EstNside.Npix(), nb_channel );
		}
		else
		{
			Dat_work.alloc( TempMap.Npix(), nb_channel );
		}
		
		if( Nz == 14 )
		{//14 maps read selection of maps for GMCA must be done
			Dat_list.alloc( TempMap.Npix(), nb_channel );
		
			for( int k=0; k < nb_channel; k++ )
			{
				for( int i=0; i < Dat_list.nx(); i++ )
				{
					TempMap[i] = Dat_file_in( i, list_channels_set(k) );
					Dat_list(i,k) = TempMap[i];
				}
			
				if( EstNside != Nside)
				{
					TempMap_EstNside.import_via_alm( TempMap );
		
					Mean = TempMap_EstNside.average();
					TempMap_EstNside.add( -Mean );
		
					for( int i=0; i < TempMap_EstNside.Npix(); i++ )
					{
						Dat_work(i,k) = TempMap_EstNside[i];
					}
				}
				else
				{
					Mean = TempMap.average();
					TempMap.add( -Mean );
		
					for( int i=0; i < TempMap.Npix(); i++ )
					{
						Dat_work(i,k) = TempMap[i];
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
				Dat_list = Dat_file_in;
			
				for( int k=0; k < nb_channel; k++ )
				{
					for( int i=0; i < Dat_list.nx(); i++ )
					{
						TempMap[i] = Dat_list(i,k);
					}
		
					if( EstNside != Nside)
					{
						TempMap_EstNside.import_via_alm( TempMap );
		
						Mean = TempMap_EstNside.average();
						TempMap_EstNside.add( -Mean );
		
						for( int i=0; i < TempMap_EstNside.Npix(); i++ )
						{
							Dat_work(i,k) = TempMap_EstNside[i];
						}
					}
					else
					{
						Mean = TempMap.average();
						TempMap.add( -Mean );
		
						for( int i=0; i < TempMap.Npix(); i++ )
						{
							Dat_work(i,k) = TempMap[i];
						}
					}//end if( EstNside != Nside)
				}//end for( int k=0; k < nb_channel; k++ )
			}//end if Nz != nb_channel
		}
	
		Dat_file_in.free();
    }
    else
    {//Read a ASCII file for localization of Healpix files 
    	char List_file_in[14][256];
    	char str[256];
    	int nb = 0;
    	
    	FILE * pFile;
	
		pFile = fopen( Name_Cube_In, "r" );
		
		if( pFile != NULL )
		{
			int count;
			while( !feof( pFile ) )
			{
				count = fscanf( pFile, "%s", str );
				if( count <= 0 )
				{
					fprintf( OUTMAN, "Problem occurs while reading file %s \n", Name_Cube_In );
					exit(-1);
				}
				for( int j=0; j <= strlen(str); j++ )
				{
					List_file_in[nb][j] = str[j];
				}
				nb = nb + 1;
			}

    		fclose( pFile );
		}
		else
		{
			fprintf( OUTMAN, "Could not open input file %s in ASCII mode\n", Name_Cube_In );
			usage(argv);
			exit(-1);
		}
		
		if( nb == 14 )
		{//14 maps names read selection of maps for GMCA must be done
			//centering data and resizing to new EstNside for GMCA only !!!
			if( Nested == false )
			{
				ReadMap.read( List_file_in[ list_channels_set(0) ], true );// we work in RING
			}
			else
			{
				ReadMap.read( List_file_in[ list_channels_set(0) ], false );
				if( ReadMap.Scheme() ==  DEF_MRS_ORDERING )
				{
					ReadMap.swap_scheme();// we work in NESTED
				}
			}
			if( Nside == 0 )
			{
				Nside = ReadMap.Nside();
				TempMap.alloc( Nside, Nested );
				TempMap = ReadMap;
			}
			else
			{
				TempMap.alloc( Nside, Nested );
			
				if( Nside != ReadMap.Nside() )
				{
					cout << "Input map " << List_file_in[ list_channels_set(0) ] << " is in nside = " << ReadMap.Nside() << " : resizing needed\n";
					TempMap.import_via_alm( ReadMap );
				}
				else
				{
					TempMap = ReadMap;
				}
			}
			if( EstNside == 0 )
			{
				EstNside = Nside;
			}
		
			if( Verbose == True ) cout << "Map read: " << List_file_in[ list_channels_set(0) ] << '\n';
		
			if( EstNside != Nside)
			{
				TempMap_EstNside.alloc( EstNside, Nested );
		
				Dat_work.alloc( TempMap_EstNside.Npix(), nb_channel );
			}
			else
			{
				Dat_work.alloc( TempMap.Npix(), nb_channel );
			}
		
			Dat_list.alloc( TempMap.Npix(), nb_channel );
		
			for( int i=0; i < Dat_list.nx(); i++ )
			{
				Dat_list(i,0) = TempMap[i];
			}
			if( EstNside != Nside)
			{
				if( EstNside != ReadMap.Nside() )
				{
					TempMap_EstNside.import_via_alm( ReadMap );
				}
				else
				{
					TempMap_EstNside = ReadMap;
				}
		
				Mean = TempMap_EstNside.average();
				TempMap_EstNside.add( -Mean );
		
				for( int i=0; i < TempMap_EstNside.Npix(); i++ )
				{
					Dat_work(i,0) = TempMap_EstNside[i];
				}
			}
			else
			{
				Mean = TempMap.average();
				TempMap.add( -Mean );
		
				for( int i=0; i < TempMap.Npix(); i++ )
				{
					Dat_work(i,0) = TempMap[i];
				}
			}//end if( EstNside != Nside)
								
			for( int k=1; k < nb_channel; k++ )
			{
				if( Nested == false )
				{
					ReadMap.read( List_file_in[ list_channels_set(k) ], true );// we work in RING
				}
				else
				{
					ReadMap.read( List_file_in[ list_channels_set(k) ], false );
					if( ReadMap.Scheme() ==  DEF_MRS_ORDERING )
					{
						ReadMap.swap_scheme();// we work in NESTED
					}
				}
				if( Nside != ReadMap.Nside() )
				{
					cout << "Input map " << List_file_in[ list_channels_set(k) ] << " is in nside = " << ReadMap.Nside() << " : resizing needed\n";
					TempMap.import_via_alm( ReadMap );
				}
				else
				{
					TempMap = ReadMap;
				}
				
				if( Verbose == True ) cout << "Map read: " << List_file_in[ list_channels_set(k) ] << '\n';
				
				for( int i=0; i < Dat_list.nx(); i++ )
				{
					Dat_list(i,k) = TempMap[i];
				}
			
				if( EstNside != Nside)
				{
					if( EstNside != ReadMap.Nside() )
					{
						TempMap_EstNside.import_via_alm( ReadMap );
					}
					else
					{
						TempMap_EstNside = ReadMap;
					}
		
					Mean = TempMap_EstNside.average();
					TempMap_EstNside.add( -Mean );
		
					for( int i=0; i < TempMap_EstNside.Npix(); i++ )
					{
						Dat_work(i,k) = TempMap_EstNside[i];
					}
				}
				else
				{
					Mean = TempMap.average();
					TempMap.add( -Mean );
		
					for( int i=0; i < TempMap.Npix(); i++ )
					{
						Dat_work(i,k) = TempMap[i];
					}
				}//end if( EstNside != Nside)
			}//end for( int k=0; k < nb_channel; k++ )	
		}
		else
		{
			if( nb != nb_channel )
			{//error number of selected channels different from number of read maps
				fprintf( OUTMAN, "Bad number of selected channels, number different from read maps \n" );
				fprintf( OUTMAN, "Number of selected channels: %d\n", nb_channel );
				fprintf( OUTMAN, "Number of map read: %d\n", nb );
				usage(argv);
				exit(-1);
			}
			else
			{
				//centering data and resizing to new EstNside for GMCA only !!!
				if( Nested == false )
				{
					ReadMap.read( List_file_in[ 0 ], true );// we work in RING
				}
				else
				{
					ReadMap.read( List_file_in[ 0 ], false );
					if( ReadMap.Scheme() ==  DEF_MRS_ORDERING )
					{
						ReadMap.swap_scheme();// we work in NESTED
					}
				}
				if( Nside == 0 )
				{
					Nside = ReadMap.Nside();
					TempMap.alloc( Nside, Nested );
					TempMap = ReadMap;
				}
				else
				{
					TempMap.alloc( Nside, Nested );
			
					if( Nside != ReadMap.Nside() )
					{
						cout << "Input map " << List_file_in[ 0 ] << " is in nside = " << ReadMap.Nside() << " : resizing needed\n";
						TempMap.import_via_alm( ReadMap );
					}
					else
					{
						TempMap = ReadMap;
					}
				}
				if( EstNside == 0 )
				{
					EstNside = Nside;
				}
		
				if( Verbose == True ) cout << "Map read: " << List_file_in[ 0 ] << '\n';
		
				if( EstNside != Nside)
				{
					TempMap_EstNside.alloc( EstNside, Nested );
		
					Dat_work.alloc( TempMap_EstNside.Npix(), nb_channel );
				}
				else
				{
					Dat_work.alloc( TempMap.Npix(), nb_channel );
				}
		
				Dat_list.alloc( TempMap.Npix(), nb_channel );
		
				for( int i=0; i < Dat_list.nx(); i++ )
				{
					Dat_list(i,0) = TempMap[i];
				}
				if( EstNside != Nside)
				{
					if( EstNside != ReadMap.Nside() )
					{
						TempMap_EstNside.import_via_alm( ReadMap );
					}
					else
					{
						TempMap_EstNside = ReadMap;
					}
		
					Mean = TempMap_EstNside.average();
					TempMap_EstNside.add( -Mean );
		
					for( int i=0; i < TempMap_EstNside.Npix(); i++ )
					{
						Dat_work(i,0) = TempMap_EstNside[i];
					}
				}
				else
				{
					Mean = TempMap.average();
					TempMap.add( -Mean );
		
					for( int i=0; i < TempMap.Npix(); i++ )
					{
						Dat_work(i,0) = TempMap[i];
					}
				}//end if( EstNside != Nside)
			
				for( int k=1; k < nb_channel; k++ )
				{
					if( Nested == false )
					{
						ReadMap.read( List_file_in[ k ], true );// we work in RING
					}
					else
					{
						ReadMap.read( List_file_in[ k ], false );
						if( ReadMap.Scheme() ==  DEF_MRS_ORDERING )
						{
							ReadMap.swap_scheme();// we work in NESTED
						}
					}
					if( Nside != ReadMap.Nside() )
					{
						cout << "Input map " << List_file_in[ k ] << " is in nside = " << ReadMap.Nside() << " : resizing needed\n";
						TempMap.import_via_alm( ReadMap );
					}
					else
					{
						TempMap = ReadMap;
					}
					
					if( Verbose == True ) cout << "Map read: " << List_file_in[ k ] << '\n';
					
					for( int i=0; i < Dat_list.nx(); i++ )
					{
						Dat_list(i,k) = TempMap[i];
					}
		
					if( EstNside != Nside)
					{
						if( EstNside != ReadMap.Nside() )
						{
							TempMap_EstNside.import_via_alm( ReadMap );
						}
						else
						{
							TempMap_EstNside = ReadMap;
						}
		
						Mean = TempMap_EstNside.average();
						TempMap_EstNside.add( -Mean );
		
						for( int i=0; i < TempMap_EstNside.Npix(); i++ )
						{
							Dat_work(i,k) = TempMap_EstNside[i];
						}
					}
					else
					{
						Mean = TempMap.average();
						TempMap.add( -Mean );
		
						for( int i=0; i < TempMap.Npix(); i++ )
						{
							Dat_work(i,k) = TempMap[i];
						}
					}//end if( EstNside != Nside)
				}//end for( int k=1; k < nb_channel; k++ )
			}//end if nb != nb_channel
		}
    }// end if ASCII == false
	
	if (Verbose == True) cout << "Nside = " << Nside << " EstNside = " << EstNside << " NbrChannels = " << nb_channel << endl;
	
	NbrScale2d = log( (double) EstNside )/log( (double) 2.0 ) - 2;
	cout << "GMCA - NbrScale : " << NbrScale2d << '\n';
	  
	MRS_GMCA WT;
       
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
	WT.BandThrd = true;
	
	// Read Haslam map
	Hdmap Haslam_Map, Haslam_Map_EstNside;
	char *MRS_ENV = getenv("ISAP");
	sprintf( Haslam_FileName, "%s/data/%s", MRS_ENV, DEF_HASLAM_NAME );
	Haslam_Map.read( Haslam_FileName );// Map is in RING scheme !!! "408MHz_512.fits"	 "/dsm/cosmo01/planck/wg2/ancillary/galactic/408MHz_512.fits"
	Haslam_Map_EstNside.alloc( EstNside, Haslam_Map.Scheme() );
	Haslam_Map_EstNside.import_via_alm( Haslam_Map );
	if( Nested == false )
	{//we work in RING
		if( Haslam_Map_EstNside.Scheme() !=  DEF_MRS_ORDERING )
		{
			Haslam_Map_EstNside.swap_scheme();
		}
	}
	else
	{//we work in NESTED
		if( Haslam_Map_EstNside.Scheme() ==  DEF_MRS_ORDERING )
		{
			Haslam_Map_EstNside.swap_scheme();
		}
	}
	
	C_OWT Haslam_OWT;
	Haslam_OWT.alloc( EstNside, NbrScale2d, Nested );
	Haslam_OWT.transform( Haslam_Map_EstNside );
	
	fltarray WT_Haslam( Haslam_OWT.n_elements() );
	for( int m=0; m < Haslam_OWT.n_elements(); m++ )
	{
		WT_Haslam(m) = Haslam_OWT.WTTrans[m];
	}
	
	int status;
	status = WT.set_planck( FixCMB, FixSZ, FixFreeFree, UseSync, WT_Haslam, NbrScale2d, NbrSources, list_channels_set );
	if( status == -1 )
	{
		fprintf(OUTMAN, "Error: selected number of sources lower than number of known sources setted...\n");
		exit(-1);
	}
       
	fltarray TransDat;
	
	TransDat.alloc( Dat_work.nx(), Dat_work.ny() );
	
	WT.transform_sources( Dat_work, TransDat, false );
          
	// run the GMCA method
	if (Verbose == True) cout << "Running GMCA ... "<< endl;
	   
	fltarray TransTabSource;
	int NbrCoef = TransDat.nx();
	TransTabSource.alloc( NbrCoef, NbrSources );
	
	WT.GMCA::Verbose = Verbose;
	WT.run_gmca( TransDat, TransTabSource );

	// Reconstruction :
	if (Verbose == True) cout << "Reconstruction ... "<< endl;
	fltarray EstSources;
	WT.recons_sources( Dat_list, EstSources );

	if (Verbose == True) cout << "Write results ... "<< endl;

    char FN[512];
    sprintf( FN, "%s_sources.fits", Name_Out );
    fits_write_fltarr( FN, EstSources );
    
    sprintf( FN, "%s_mat.fits", Name_Out );
    fits_write_fltarr( FN, WT.MixingMat );
    
    sprintf( FN, "%s_invmat.fits", Name_Out );
    fits_write_fltarr( FN, WT.InvMixingMat );
    
    if( estim_bias == true )
	{
		dblarray MapBeam, Map_Beam_read, Bias;
		fits_read_dblarr( BeamChannels_List, Map_Beam_read );
		
		if( Map_Beam_read.nx() == 14 )
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
		if( list_channels_set.min() >= 9 )
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
			if( list_channels_set(i) < 9 )
			{
				planck_channel_set = true;
			}
			if( list_channels_set(i) >= 9 )
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
			fltarray std_file_in, std_work, std_estim;
			std_work.alloc( 12*Nside*Nside, nb_channel );
			//std_estim.alloc( 12*Nside*Nside, NbrSources );
			
			float Rfact;
			
			if( ASCII == false )
			{
				fits_read_fltarr( STD_Cube_In, std_file_in );
		
				int Nside_std = sqrt( std_file_in.nx() / 12 );
		
				int Nz_std = std_file_in.ny();
				
				if( Nside_std != Nside)
				{
					TempMap_EstNside.alloc( Nside, Nested );
				}
				
				TempMap.alloc( Nside_std, Nested );	
		
				if( Nz_std == 14 )
				{//14 std maps read selection of maps must be done
					
					for( int k=0; k < nb_channel; k++ )
					{
						for( int i=0; i < std_file_in.nx(); i++ )
						{
							TempMap[i] = std_file_in( i, list_channels_set(k) );
						}
						
						if( list_channels_set(k) <= 8 )
						{
							Rfact = 2048 / ( (float) Nside );//Rfact = 2048.0*2048.0 / ( (float) Nside*Nside );
						}
						else
						{
							Rfact = 512 / ( (float) Nside );//Rfact = 512.0*512.0 / ( (float) Nside*Nside );
						}
			
						if( Nside_std != Nside)
						{
							TempMap_EstNside.import_via_alm( TempMap );
			
							for( int i=0; i < TempMap_EstNside.Npix(); i++ )
							{
								std_work(i,k) = TempMap_EstNside[i]/Rfact;
							}
						}
						else
						{	
							for( int i=0; i < TempMap.Npix(); i++ )
							{
								std_work(i,k) = TempMap[i]/Rfact;
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
							for( int i=0; i < std_file_in.nx(); i++ )
							{
								TempMap[i] = std_file_in(i,k);
							}
							
							if( list_channels_set(k) <= 8 )
							{
								Rfact = 2048 / ( (float) Nside );//Rfact = 2048.0*2048.0 / ( (float) Nside*Nside );
							}
							else
							{
								Rfact = 512 / ( (float) Nside );//Rfact = 512.0*512.0 / ( (float) Nside*Nside );
							}
		
							if( Nside_std != Nside)
							{
								TempMap_EstNside.import_via_alm( TempMap );
			
								for( int i=0; i < TempMap_EstNside.Npix(); i++ )
								{
									std_work(i,k) = TempMap_EstNside[i]/Rfact;
								}
							}
							else
							{	
								for( int i=0; i < TempMap.Npix(); i++ )
								{
									std_work(i,k) = TempMap[i]/Rfact;
								}
							}//end if( Nside_std != Nside)
						}//end for( int k=0; k < nb_channel; k++ )
					}//end if Nz != nb_channel
				}
	
				std_file_in.free();				
			}
			else
			{
				//Read a ASCII file for localization of Healpix standard deviation files 
    			char List_STDfile_in[14][256];
    			char str2[256];
    			int nb_2 = 0;
    	
    			FILE * pFile_std;
	
				pFile_std = fopen( STD_Cube_In, "r" );
		
				if( pFile_std != NULL )
				{
					int count_std;
					while( !feof( pFile_std ) )
					{
						count_std = fscanf( pFile_std, "%s", str2 );
						if( count_std <= 0 )
						{
							fprintf( OUTMAN, "Problem occurs while reading file %s \n", STD_Cube_In );
							exit(-1);
						}
						for( int j=0; j <= strlen(str2); j++ )
						{
							List_STDfile_in[nb_2][j] = str2[j];
						}
						nb_2 = nb_2 + 1;
					}

    				fclose( pFile_std );
				}
				else
				{
					fprintf( OUTMAN, "Could not open input file %s in ASCII mode\n", STD_Cube_In );
					usage(argv);
					exit(-1);
				}
				
				if( nb_2 == 14 )
				{//14 std maps names read selection of maps must be done
					if( Nested == false )
					{
						ReadMap.read( List_STDfile_in[ list_channels_set(0) ], true );// we work in RING
					}
					else
					{
						ReadMap.read( List_STDfile_in[ list_channels_set(0) ], false );
						if( ReadMap.Scheme() ==  DEF_MRS_ORDERING )
						{
							ReadMap.swap_scheme();// we work in NESTED
						}
					}
					TempMap.alloc( Nside, Nested );
			
					if( Nside != ReadMap.Nside() )
					{
						cout << "Input map " << List_STDfile_in[ list_channels_set(0) ] << " is in nside = " << ReadMap.Nside() << " : resizing needed\n";
						TempMap.import_via_alm( ReadMap );
					}
					else
					{
						TempMap = ReadMap;
					}
				
					if( Verbose == True ) cout << "Map read: " << List_STDfile_in[ list_channels_set(0) ] << '\n';
				
					if( list_channels_set(0) <= 8 )
					{
						Rfact = 2048 / ( (float) Nside );//Rfact = 2048.0*2048.0 / ( (float) Nside*Nside );
					}
					else
					{
						Rfact = 512 / ( (float) Nside );//Rfact = 512.0*512.0 / ( (float) Nside*Nside );
					}
				
					for( int i=0; i < std_work.nx(); i++ )
					{
						std_work(i,0) = TempMap[i]/Rfact;
					}
							
					for( int k=1; k < nb_channel; k++ )
					{
						if( Nested == false )
						{
							ReadMap.read( List_STDfile_in[ list_channels_set(k) ], true );// we work in RING
						}
						else
						{
							ReadMap.read( List_STDfile_in[ list_channels_set(k) ], false );
							if( ReadMap.Scheme() ==  DEF_MRS_ORDERING )
							{
								ReadMap.swap_scheme();// we work in NESTED
							}
						}
						
						if( Nside != ReadMap.Nside() )
						{
							cout << "Input map " << List_STDfile_in[ list_channels_set(k) ] << " is in nside = " << ReadMap.Nside() << " : resizing needed\n";
							TempMap.import_via_alm( ReadMap );
						}
						else
						{
							TempMap = ReadMap;
						}
				
						if( Verbose == True ) cout << "Map read: " << List_STDfile_in[ list_channels_set(k) ] << '\n';
						
						if( list_channels_set(k) <= 8 )
						{
							Rfact = 2048 / ( (float) Nside );//Rfact = 2048.0*2048.0 / ( (float) Nside*Nside );
						}
						else
						{
							Rfact = 512 / ( (float) Nside );//Rfact = 512.0*512.0 / ( (float) Nside*Nside );
						}
				
						for( int i=0; i < std_work.nx(); i++ )
						{
							std_work(i,k) = TempMap[i]/Rfact;
						}
					}
				}
				else
				{
					if( nb_2 != nb_channel )
					{//error number of selected channels different from number of read maps
						fprintf( OUTMAN, "Bad number of selected channels, number different from read std maps \n" );
						fprintf( OUTMAN, "Number of selected channels: %d\n", nb_channel );
						fprintf( OUTMAN, "Number of map read: %d\n", nb_2 );
						usage(argv);
						exit(-1);
					}
					else
					{
						if( Nested == false )
						{
							ReadMap.read( List_STDfile_in[ 0 ], true );// we work in RING
						}
						else
						{
							ReadMap.read( List_STDfile_in[ 0 ], false );
							if( ReadMap.Scheme() ==  DEF_MRS_ORDERING )
							{
								ReadMap.swap_scheme();// we work in NESTED
							}
						}
						TempMap.alloc( Nside, Nested );
			
						if( Nside != ReadMap.Nside() )
						{
							cout << "Input map " << List_STDfile_in[ 0 ] << " is in nside = " << ReadMap.Nside() << " : resizing needed\n";
							TempMap.import_via_alm( ReadMap );
						}
						else
						{
							TempMap = ReadMap;
						}
				
						if( Verbose == True ) cout << "Map read: " << List_STDfile_in[ 0 ] << '\n';
				
						if( list_channels_set(0) <= 8 )
						{
							Rfact = 2048 / ( (float) Nside );//Rfact = 2048.0*2048.0 / ( (float) Nside*Nside );
						}
						else
						{
							Rfact = 512 / ( (float) Nside );//Rfact = 512.0*512.0 / ( (float) Nside*Nside );
						}
				
						for( int i=0; i < std_work.nx(); i++ )
						{
							std_work(i,0) = TempMap[i]/Rfact;
						}
						
						for( int k=1; k < nb_channel; k++ )
						{
							if( Nested == false )
							{
								ReadMap.read( List_STDfile_in[ k ], true );// we work in RING
							}
							else
							{
								ReadMap.read( List_STDfile_in[ k ], false );
								if( ReadMap.Scheme() ==  DEF_MRS_ORDERING )
								{
									ReadMap.swap_scheme();// we work in NESTED
								}
							}
							if( Nside != ReadMap.Nside() )
							{
								cout << "Input map " << List_STDfile_in[ k ] << " is in nside = " << ReadMap.Nside() << " : resizing needed\n";
								TempMap.import_via_alm( ReadMap );
							}
							else
							{
								TempMap = ReadMap;
							}
					
							if( Verbose == True ) cout << "Map read: " << List_STDfile_in[ k ] << '\n';
							
							if( list_channels_set(k) <= 8 )
							{
								Rfact = 2048 / ( (float) Nside );//Rfact = 2048.0*2048.0 / ( (float) Nside*Nside );
							}
							else
							{
								Rfact = 512 / ( (float) Nside );//Rfact = 512.0*512.0 / ( (float) Nside*Nside );
							}
					
							for( int i=0; i < std_work.nx(); i++ )
							{
								std_work(i,k) = TempMap[i]/Rfact;
							}
						}//end for( int k=1; k < nb_channel; k++ )
					}//end if nb_2 != nb_channel
				}// end if( nb_2 == 14 )
			}// end if ASCII
			
			fltarray var_in, var_out, pinv_2;
			MatOper MAT;
			
			var_in = std_work;
			for(long i=0; i < var_in.n_elem(); i++ )
			{
				var_in(i) = (float) pow( (double) var_in(i), 2 );
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
			MAT.mat_mult( pinv_2, var_in, var_out );
			
			std_estim = var_out;
			for(long i=0; i < std_estim.n_elem(); i++ )
			{
				std_estim(i) = sqrt( std_estim(i) );
			}
			//std_estim = var_out^0.5;
			
			sprintf( FN, "%s_estim_src_std.fits", Name_Out );
    		fits_write_fltarr( FN, std_estim );
		}// end if( ( planck_channel_set == true )&&( wmap_channel_set == true) )
	}// end if( estim_std_src == true )

    exit(0);
}

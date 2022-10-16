/******************************************************************************
**                   Copyright (C) 2009 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Florent Sureau
**
**    Date:  09/12/09
**    
**    File:  mrs_planck_lgmca_inv.cc
**
*******************************************************************************
**
**    DESCRIPTION CMB estimation on the sphere using lGMCA inversion weights 
**    ----------- 
**                 
**    Usage: mrs_lgmca_inv options NameCubeIn Name_Beams Name_Weights Name_Struc Name_Out
**
******************************************************************************/

#include "HealpixClass.h"
#include "MRS_Sparse.h"
#include "mrs_planck_lGMCA.h"
#include <sys/time.h>
#include <omp.h>

char Name_Cube_In[1024];
char Name_Beams[1024];
char Name_Weights[1024];
char Name_Struc[1024];
char Name_Out[1024];
int NbrScale2d = 5;


extern int  OptInd;
extern char *OptArg;

extern int  GetOpt(int argc, char *const*argv, char *opts);

bool Verbose=false;
bool Nested_in=true;
bool set_threads=false;
bool save_sameres=false;
bool debug_mode=false;

int NsideRec;

int outer_loop_threads;
int inner_loop_threads;

bool weight_flag = false;

/*********************************************************************/
 
static void openmp_status()
  {
   #ifdef _OPENMP    
	    cout << "Application was compiled with OpenMP support," << endl;
    	if (omp_get_max_threads() == 1)
      	cout << "but running with one process only." << endl;
    	else
    	  cout << "running with up to " << omp_get_max_threads()
           << " processes." << endl;
           if (omp_get_nested()){
           		cout << "Nested Parallelism allowed" <<endl;
           		cout << "running outer loops with " << outer_loop_threads
           << " processes." << endl;
           	cout << "running inner loops with " << inner_loop_threads
           << " processes." << endl;

           } else cout << "running with" << outer_loop_threads
           << " processes." << endl;
    #else
     cout << "Application was compiled without OpenMP support;" << endl
         << "running in scalar mode." << endl;
    #endif
  }
/*********************************************************************/


static void usage(char *argv[])
{
    fprintf(OUTMAN, "Usage: %s [options] NameCubeIn Name_Beams Name_Weights Name_Struc Name_Out\n\n", argv[0]);
    fprintf(OUTMAN, "   where  :  \n");

    
    fprintf(OUTMAN, "         [NameCubeIn]\n");
    fprintf(OUTMAN, "             Name of fits file containing the input data [Npixels,Nchannel]. \n");
    fprintf(OUTMAN, "         [Name_Beams]\n");
    fprintf(OUTMAN, "             Name of fits file containing the beams [lmax,Nchannel]. \n");
    fprintf(OUTMAN, "         [Name_Weights]\n");
    fprintf(OUTMAN, "             Name of fits file containing the inversion weights for the wavelet coefficient for the band 0 [Npixels,Nscales,Nchannels]. \n");     
    fprintf(OUTMAN, "             This should be of the form \"*_Band0_*\" and the other resolution bands should be named \"*_BandY_*\". \n");     
    fprintf(OUTMAN, "         [Name_Struc]\n");
    fprintf(OUTMAN, "             Name of fits file containing the structure describing each required parameter [Nparam,NResBands]. \n"); 
    fprintf(OUTMAN, "             Each line should be of the form [Nchannels_used,RefChannel,used_Channels,Lmax].\n");
    fprintf(OUTMAN, "             For instance: [3 1 1 2 5 0 0 512] will use channel 1 as reference, process 1,2 and 5 (index in array: 0,1,4) degraded to the resolution of the ref channel, and use a Lmax of 512.\n");     
    fprintf(OUTMAN, "         [Name_Out]\n");
    fprintf(OUTMAN, "             Name of fits file for output of band 0. \n");     
    fprintf(OUTMAN, "             This should be of the form \"*_Band0_*\", \"_sameres_inverted\" will be added to this name, and the other resolution bands will be named accordingly \"*_BandY_*\". \n");     
    fprintf(OUTMAN, "   and [options] are  :  \n");
    fprintf(OUTMAN, "         [-S NsideRec]\n");
    fprintf(OUTMAN, "             Flag to indicate that we only want to put channels to the same resolution, with Nside Nsiderec. The inversion structure should only contain one line. \n");     
    fprintf(OUTMAN, "         [-s]\n");
    fprintf(OUTMAN, "             Flag to indicate that we want to save intermediate maps put to same resolution. \n");     
    fprintf(OUTMAN, "         [-r]\n");
    fprintf(OUTMAN, "             Flag to indicate that the input maps are in ring format, and not nested. Ouptut will however be in nested format. \n");     
    fprintf(OUTMAN, "         [-t oloop:iloop]\n");
    fprintf(OUTMAN, "             Use oloop outer threads and iloop inner threads . \n");     
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

	NsideRec=0;
    /* get options */
	while ((c = GetOpt(argc,argv,"?zrZ:t:S:sdw")) != -1) 
    {
		switch (c) 
        {  
        	case 't':
	        	if (sscanf(OptArg,"%d:%d",&outer_loop_threads,&inner_loop_threads) < 1)
				{
		   			fprintf(OUTMAN, "Error: syntaxe is outer_loop_threads:inner_loop_threads ... \n");
		   			exit(-1);
				}
	        	set_threads=true;
	            break;
	  		case 'r':
	        	Nested_in = false;
				break;	  		 					 				
 			case 'S':
	        	if (sscanf(OptArg,"%d",&NsideRec) != 1)
				{
		   			fprintf(OUTMAN, "Error: syntaxe is Nside... \n");
		   			exit(-1);
				}
                break;
            case 's':
	        	save_sameres=true;
                break;
	    case 'd':
	        	debug_mode=true;
                break;
		case 'w':
	         	weight_flag = true;
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
    if (OptInd < argc) strcpy(Name_Cube_In, argv[OptInd++]);
    else usage(argv);
    if (OptInd < argc) strcpy(Name_Beams, argv[OptInd++]);
    else usage(argv);
    if (OptInd < argc) strcpy(Name_Weights, argv[OptInd++]);
    else usage(argv);
    if (OptInd < argc) strcpy(Name_Struc, argv[OptInd++]);
    else usage(argv);
    if (OptInd < argc) strcpy(Name_Out, argv[OptInd++]);
    else usage(argv);

	/* make sure there are not too many parameters */
	if (OptInd < argc) {
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
	fltarray Dat_file_out;
	int index=0,f,nproc,nthr,mlt4,mlt6,mlt12;
    /* Get command line arguments, open input file(s) if necessary */
	fitsstruct Header;
	char Cmd[512];
	struct timeval start,end,diff;

	Cmd[0] = '\0';
	for (int k =0; k < argc; k++) sprintf(Cmd, "%s %s", Cmd, argv[k]);
	   
    transinit(argc, argv);
	gettimeofday(&start, (struct timezone *) NULL);

	lGMCA_inv Strucinv;
		
	#ifdef _OPENMP
		int nbprocs=omp_get_num_procs();
		omp_set_nested(1);
		if(nbprocs<4) omp_set_num_threads(nbprocs);	
		else omp_set_num_threads((int) (nbprocs/2));
		nthr=omp_get_max_threads();
		if(set_threads==true) {
			if(outer_loop_threads * inner_loop_threads > nbprocs) {
				printf("Too many threads assigned (%d * %d vs %d procs)\n",outer_loop_threads,inner_loop_threads,nbprocs);
			}
		} else {
			if(nbprocs<4) {
				outer_loop_threads=nbprocs;
				inner_loop_threads=1;	
			} else {
				mlt4=floor(nthr/4)*4;
				mlt6=floor(nthr/6)*6;
				mlt12=floor(nthr/12)*12;
				if((mlt12>=mlt6)&&(mlt12>=mlt4)) {
					outer_loop_threads=12;
					inner_loop_threads=floor(nthr/12);
				} else if (mlt6>=mlt4) {
					outer_loop_threads=6;
					inner_loop_threads=floor(nthr/6);
				} else {
					outer_loop_threads=4;
					inner_loop_threads=floor(nthr/4);
				}
			}
			omp_set_num_threads(outer_loop_threads*inner_loop_threads);
		}
		printf("Use %d outer threads and %d inner threads for nb procs=%d\n",outer_loop_threads,inner_loop_threads,nbprocs);
    #endif
    openmp_status();
  
    Strucinv.Timer=true;
    Strucinv.Verbose=false;
    Strucinv.Debug=debug_mode;
    Strucinv.set_nb_thr(outer_loop_threads,inner_loop_threads);
    printf("Process Bands, NsideRec=%d\n",NsideRec);
    if(NsideRec >0) { //just put to same resolution first channel
		Strucinv.SetOnlyResolution(Name_Cube_In,Name_Beams,Name_Struc,Name_Out,NsideRec); 
	} else {
		if (weight_flag==true) {
			Strucinv.Invert_Weights_Only(Name_Cube_In,Name_Weights,Name_Struc,Name_Out);
			
		} else {
			Strucinv.save_sameres=save_sameres;
			Strucinv.Process_AllBands(Name_Cube_In,Name_Beams,Name_Weights,Name_Struc,Name_Out);
		}    
    }
    
    gettimeofday(&end, (struct timezone *) NULL);
    timersub(&end,&start,&diff);
    printf("Total Time=%f\n",(double) diff.tv_sec+(double)diff.tv_usec/1000000.);



    exit(0);
}

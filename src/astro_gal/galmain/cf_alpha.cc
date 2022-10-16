/******************************************************************************
**                   Copyright (C) 1999 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**	  Author: Antoine Labatie
**
**    Date:  19/12/2011
**
**    File:  cf_ana.cc
**
*******************************************************************************
**
**    DESCRIPTION  Correlation function program with alpha belonging
**    -----------
**
******************************************************************************/

#include "Array.h"
#include "IM_IO.h"
#include "DefPoint.h"
#include "cf_alpha.h"
#include <time.h>

char Name_Imag_Out[256];		/* output file name */
char Name_Imag_Out_Prefix[256];		/* output file prefix name */
char Name_Imag_Out_Suffix[256];		/* output file suffix name */

char Name_Imag_In[256];			/* input file image */
char NameRndFile[256];			/* random catalogue file name */
char NameDataWeightFile[256];   /* data catalogue weights file name */
char NameRndWeightFile[256];    /* random catalogue weights file name */
char NameDataAlphaFile[256];    /* data catalogue alpha belonging file name */
char NameRndAlphaFile[256];	    /* random catalogue alpha belonging file name */

extern int  OptInd;
extern char *OptArg;

extern int  GetOpt(int argc, char **argv, char *opts);

Bool Verbose=False;

unsigned int InitRnd = time(NULL);
float DistMin=-1;
float DistMax=-1;

float Step=1.;
int NpRnd=-1;

Bool UseDataWeight=False;
Bool UseRndWeight=False;
Bool UseDataAlpha=False;
Bool UseRndAlpha=False;
Bool ReadSimu = False;

int nalpha;
double AlphaMin;
double AlphaMax;



//maximum number of procs used for the loops
int Nproc_max=40;

/****************************************************************************/

static void usage(char *argv[])
{
    fprintf(OUTMAN, "Usage: %s options in_catalogue result_suffix\n\n", argv[0]);
    fprintf(OUTMAN, "   where options =  \n");

    fprintf(OUTMAN, "         [-I InitRandomVal]\n");
    fprintf(OUTMAN, "             Value used for random value generator initialization.\n");
    fprintf(OUTMAN, "             default is 100. \n");
    manline();

    fprintf(OUTMAN, "         [-n NbrRnd]\n");
    fprintf(OUTMAN, "             Number of random points.\n");
    fprintf(OUTMAN, "             default is the number of data points. \n"); 
    manline();
	
    fprintf(OUTMAN, "         [-s BinStep]\n");
    fprintf(OUTMAN, "             Pair separations bin size.\n");
    fprintf(OUTMAN, "             default is %f. \n", Step);
    manline();

    fprintf(OUTMAN, "         [-m SepMin]\n");
    fprintf(OUTMAN, "             Pair separations min.\n");
    fprintf(OUTMAN, "             Must be set. \n");
    manline();

    fprintf(OUTMAN, "         [-M SepMax]\n");
    fprintf(OUTMAN, "             Pair separations max.\n");
    fprintf(OUTMAN, "             Must be set. \n");
    manline();

    fprintf(OUTMAN, "         [-s BinStep]\n");
    fprintf(OUTMAN, "             Pair separations bin size.\n");
    fprintf(OUTMAN, "             default is %f. \n", Step);
    manline();
	
    fprintf(OUTMAN, "         [-w FileName]\n");
    fprintf(OUTMAN, "             Use data weights in Filename.\n");
    fprintf(OUTMAN, "             Default is no. \n");
    manline();
	
    fprintf(OUTMAN, "         [-W FileName]\n");
    fprintf(OUTMAN, "             Use random weights in Filename.\n");
    fprintf(OUTMAN, "             Default is no. \n");
    manline();
	
    fprintf(OUTMAN, "         [-a FileName]\n");
    fprintf(OUTMAN, "             Use alpha belonging for data in FileName.\n");
    fprintf(OUTMAN, "             Default is no. \n");
    manline();
	
    fprintf(OUTMAN, "         [-A FileName]\n");
    fprintf(OUTMAN, "             Use alpha belonging for random catalogue in FileName.\n");
    fprintf(OUTMAN, "             Default is no. \n");
    manline();

    fprintf(OUTMAN, "         [-r FileName]\n");
    fprintf(OUTMAN, "             Read simulated random catalogue in FileName.\n");
    fprintf(OUTMAN, "             Default is no. \n");
    manline();


    vm_usage();
    manline();
    verbose_usage();
    manline();
    manline();
    exit(-1);
}

/*********************************************************************/

/* GET COMMAND LINE ARGUMENTS */
static void filtinit(int argc, char *argv[])
{

    Bool DMin = False;
    Bool DMax = False;

	
	/* Require arguments (need at least one) !! */
	if(argc == 1){
		usage(argv);
	}
	/* Start at i = 1 to skip the command name. */
	int i=1;
	
    /* Check for a switch (leading "-"). */	
	while(argv[i][0] == '-') {
		
		/* Use the next character to decide what to do. */
		
		switch (argv[i][1]) {
				
			case 'n': NpRnd=atoi(argv[++i]);
				break;
				
			case 'm': DistMin=atof(argv[++i]);
				DMin=True;
				break;
	
			case 'M': DistMax=atof(argv[++i]);
				DMax=True;
				break;

			case 's': Step = atof(argv[++i]);
				break;
				
			case 'r': strcpy(NameRndFile,argv[++i]);
				ReadSimu = True;
				break;
				
			case 'w': strcpy(NameDataWeightFile,argv[++i]);
				UseDataWeight = True;  
				break;
				
			case 'W': strcpy(NameRndWeightFile,argv[++i]);
				UseRndWeight = True;  
				break;
				
			case 'a': strcpy(NameDataAlphaFile,argv[++i]);
				UseDataAlpha = True;  
				break;
				
			case 'A': strcpy(NameRndAlphaFile,argv[++i]);
				UseRndAlpha = True;  
				break;
				
			case 'I': InitRnd  = atol(argv[++i]);
				break;
				
			case 'v': Verbose = True;
				break;
				
			case '?': usage(argv);
				break;
				
			default:  usage(argv);
				break;
		}
		i++;
		if(i>=(argc-1)) //there remains less than 2 parameters
		{
			usage(argv);
			exit(-1);
		}
	}

	strcpy(Name_Imag_In, argv[i++]);
	strcpy(Name_Imag_Out_Suffix, argv[i++]);
	
	if(i < argc){
		fprintf(stderr, "Too many parameters: %s ...\n", argv[i]);
		usage(argv);
		exit(-1);
	}
	
	if ((DMin == False) || (DMax == False))
	{
		fprintf(OUTMAN, "Error: -m and -M option must be set ...\n");
		exit(-1);
	}

}


/*********************************************************************/

/* GET PARAMETERS */

void get_param()
{
	char Name_Param_File[256];
	sprintf(Name_Param_File, "../param/cf_alpha.param");	
	FILE *File=fopen(Name_Param_File,"r");
	
    if (File == NULL)
    {
		cerr << "Error: cannot open file "  <<  Name_Param_File << endl;
		exit(-1);
    }
	char Temp[256];
	double temp_double;
	int ret;
	
	ret=fscanf(File, "%s\t%lf\n", Temp, &AlphaMin);			//AlphaMin
	ret=fscanf(File, "%s\t%lf\n", Temp, &AlphaMax);			//AlphaMax
	ret=fscanf(File, "%s\t%lf\n\n", Temp, &temp_double);	//nalpha
	nalpha=round(temp_double);
		
	ret=fscanf(File, "%s\t%s\n", Temp, Name_Imag_Out_Prefix);	//Name_Out prefix
	fclose(File);
}

/*********************************************************************/


int main(int argc, char *argv[])
{
	long int times=time(NULL);
	
    int  i,k,d, nbins;
    int Naxis,Np;
	fltarray Result;
    fltarray CF_DataData, CF_DataRnd, CF_RndRnd;
    fitsstruct Header;
    char Cmd[512];
    Cmd[0] = '\0';
    for (k =0; k < argc; k++) sprintf(Cmd, "%s %s", Cmd, argv[k]);

	
	
     /* Get command line arguments, open input file(s) if necessary */
    filtinit(argc, argv);
	get_param();

	float AlphaStep=(AlphaMax-AlphaMin)/double(nalpha-1.0);
	
	//Alpha table
	fltarray AlphaTable(nalpha);
	for(int i=0;i<nalpha;i++)
		AlphaTable(i)=AlphaMin+i*AlphaStep;
	
	strcpy(Name_Imag_Out, Name_Imag_Out_Prefix); strcat(Name_Imag_Out, Name_Imag_Out_Suffix);
    if (Verbose == True)
    {
        cout << endl << endl << "PARAMETERS: " << endl << endl;
        cout << "File Name in = " << Name_Imag_In << endl;
        cout << "File Name Out = " << Name_Imag_Out << endl;
		if(UseDataWeight)     cout << "File Name Data Weight in = " << NameDataWeightFile << endl;
		if(UseRndWeight)     cout << "File Name Rnd Weight in = " << NameRndWeightFile << endl;
        cout << "SeparationMin = " << DistMin;
        cout << " SeparationMax = " << DistMax;
        cout << " SeparationStep = " << Step << endl;

        if (ReadSimu == True) cout << "Read Random catalogue in " << NameRndFile <<  endl ;
		cout << "AlphaMin = "<< AlphaMin;
		cout << " AlphaMax = "<< AlphaMax;
		cout << " AlphaStep = " << AlphaStep << endl ;
		if(UseDataAlpha)     cout << "File Name Data Alpha in = " << NameDataAlphaFile << endl;
		if(UseRndAlpha)     cout << "File Name Rnd Alpha in = " << NameRndAlphaFile << endl << endl;
    }

	//read TabData
    ArrayPoint TabData;
    TabData.read(Name_Imag_In, Verbose);
    Naxis = TabData.dim(); Np = TabData.np();
    if (TabData.TCoord == 2 && TabData.dim() == 3) TabData.toxyz();
    if (NpRnd < Np) NpRnd = Np;
	Point Pmin(Naxis),Pmax(Naxis);
    Pmin = TabData.Pmin; Pmax = TabData.Pmax;

	
	//TabData weights
	ArrayPoint TabDataWeight(1,Np);
	Point P(1); P.axis(0)=1.0;
	for(k=0;k<Np;k++) TabDataWeight(k)=P;
	if(UseDataWeight==True) TabDataWeight.read(NameDataWeightFile, False);

	//TabData alpha belonging
	ArrayPoint TabDataAlpha(2,Np);
	Point P2(1); P2.axis(0)=0.0; P2.axis(1)=nalpha;
	for(k=0;k<Np;k++) TabDataAlpha(k)=P2;
	if(UseDataAlpha==True) 
	{
		TabDataAlpha.read(NameDataAlphaFile, False);
		for(k=0;k<Np;k++) TabDataAlpha(k).axis(0)=round((TabDataAlpha(k).axis(0)-AlphaMin)/AlphaStep);
		for(k=0;k<Np;k++) TabDataAlpha(k).axis(1)=round((TabDataAlpha(k).axis(1)-AlphaMin)/AlphaStep);
		for(k=0;k<Np;k++) TabDataAlpha(k).axis(0)=max(TabDataAlpha(k).axis(0),float(0.0));
		for(k=0;k<Np;k++) TabDataAlpha(k).axis(1)=max(TabDataAlpha(k).axis(1),float(0.0));
	}
	
	
	//Check number of points data
	if(TabDataAlpha.np()!=Np) 
	{ 		cerr << "Incorrect # points in alpha belonging for data catalogue" << endl; exit(-1); 	}
	if(TabDataWeight.np()!=Np) 
	{ 		cerr << "Incorrect # weights for data catalogue" << endl; exit(-1); 	}

    // Allocation of the correlation function CLASS
    CorrFunAnaAlpha CFA(DistMin, DistMax, Step);
    if (Verbose == True)
    {
		cout << endl ;
		cout << "Data correlation function ... " << endl;
		cout << "Number of data points = " << Np << endl;
		cout << "Number of random data points = " << NpRnd << endl;
		cout << "Number of separation bins = " << CFA.np() << endl;
    }
    nbins = CFA.np();
	
	// Allocate the result array
	//binning + DD, RR, DR
	int NbLineResult=1+3;
	Result.alloc(nbins, nalpha, NbLineResult);
	for (d=0; d < nbins; d++)  
	{
		for (i=0; i < nalpha; i++)
			Result(d,i,0) = CFA.coord(d)/AlphaTable(i);
	}

	// Allocate CF_DataData, CF_DataRnd, CF_RndRnd
    CF_DataData.alloc(nbins,nalpha); CF_DataRnd.alloc(nbins,nalpha); CF_RndRnd.alloc(nbins,nalpha); 
	
    // find the pairs histogram and put it in CF_DataData
    CFA.cf_find_pairs(TabData,TabDataWeight,TabDataAlpha,CF_DataData);
						
    // Random number generator initialization
    init_random (InitRnd);

	//Start simulation
	if (Verbose == True) cout << "Simulation" << endl;
		
    ArrayPoint TabRnd(Naxis,NpRnd);
    TabRnd.Pmin = Pmin; TabRnd.Pmax = Pmax; TabRnd.TCoord = TabData.TCoord;
    for (d=0; d < 3; d++) TabRnd.BootCoord[d] = TabData.BootCoord[d];
	
	if (ReadSimu == True)
	{
		if (Verbose == True) cout << " Read " << NameRndFile << endl;
		TabRnd.read(NameRndFile, Verbose);
		TabRnd.Pmin = Pmin; TabRnd.Pmax = Pmax; NpRnd=TabRnd.np();
		if (TabRnd.TCoord == 2 && TabRnd.dim() == 3) TabRnd.toxyz();
	}
	else
	{
		TabRnd.random(TabRnd);
		// Convert into rectangular coordinates if Random catalogue generated
		// in spherical coordinates from the original form of the catalogue
		if (TabRnd.TCoord == 2 && TabRnd.dim() == 3) TabRnd.toxyz(); 
	}
	
	//TabRandom Weight
	ArrayPoint TabRndWeight(1,NpRnd);
	for(k=0;k<NpRnd;k++) TabRndWeight(k)=P;
	if(UseRndWeight==True) TabRndWeight.read(NameRndWeightFile, False);
	
	//TabRandom alpha belonging
	ArrayPoint TabRndAlpha(2,NpRnd);
	for(k=0;k<NpRnd;k++) TabRndAlpha(k)=P2;
	if(UseRndAlpha==True) 
	{
		TabRndAlpha.read(NameRndAlphaFile, False);

		for(k=0;k<NpRnd;k++) TabRndAlpha(k).axis(0)=round((TabRndAlpha(k).axis(0)-AlphaMin)/AlphaStep);
		for(k=0;k<NpRnd;k++) TabRndAlpha(k).axis(1)=round((TabRndAlpha(k).axis(1)-AlphaMin)/AlphaStep);
		for(k=0;k<NpRnd;k++) TabRndAlpha(k).axis(0)=max(TabRndAlpha(k).axis(0),float(0.0));
		for(k=0;k<NpRnd;k++) TabRndAlpha(k).axis(1)=max(TabRndAlpha(k).axis(1),float(0.0));
	}

	//Check number of points random
	if(TabRndAlpha.np()!=NpRnd) 
	{ 		cerr << "Incorrect # points in alpha belonging for random catalogue" << endl; exit(-1); 	}
	if(TabRndWeight.np()!=NpRnd) 
	{ 		cerr << "Incorrect # weights for random catalogue" << endl; exit(-1); 	}
	
	// random-random pairs histogram calculation and put the result in CF_RndRnd
	CFA.cf_find_pairs(TabRnd,TabRndWeight,TabRndAlpha,CF_RndRnd);
	// data-random pairs histogram calculation and put the result in CF_DataRnd
	CFA.cf_find_pairs(TabData,TabRnd,TabDataWeight,TabRndWeight,TabDataAlpha,TabRndAlpha,CF_DataRnd);

    
	if (Verbose == True)
	{
		long int timee=time(NULL);
		cout << "Time in sec : " << timee-times << endl; 
	}

	//make pair histo DD, DR, RR
	make_histo(CF_DataData, CF_RndRnd,  CF_DataRnd,  Result); 

	//normalize by \sum w_i * \sum_w_j
	normalize_histo(TabDataWeight,TabRndWeight,TabDataAlpha,TabRndAlpha,Result); 


    // Write the results
    Header.hd_fltarray(Result, Cmd);
    fits_write_fltarr(Name_Imag_Out, Result, &Header);
    exit(0);
}

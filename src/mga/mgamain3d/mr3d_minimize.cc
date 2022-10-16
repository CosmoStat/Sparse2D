/******************************************************************************
**                   Copyright (C) 2008 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Arnaud Woiselle
**
**    Date:  10/09/08
**    
**    File:  mr3d_fista.cc
**
*******************************************************************************
**
**    DESCRIPTION  
**    ----------- 
**                 
******************************************************************************/

#include <time.h>
#include <map>

// General includes
#include "Array.h"
#include "IM_IO.h"
#include "GetLongOptions.h"

// Method includes
#include "Fista.h"
#include "Nesterov.h"

// Transforms includes
#include "bfct.h"
#include "FCur_float.h"
#include "SB_Filter_float.h"

char* Name_Imag_In; /* input file image */
char* Name_Mask; /* input mask name */
char* Name_Imag_Out; /* output file name */
extern int  OptInd;
extern char *OptArg;
extern int  GetOpt(int argc, char **argv, char *opts);
extern char** short_argv;
extern int short_argc;

/********************************/
//		Input variables			//
/********************************/

Bool Verbose=False;
int Niter=1;
bool Simu=false;

// Transforms parameters
int NbrScale=3;
int NbrDir2d = 16;

// Domains
enum {DOM_DWT, DOM_UWT, DOM_IWT, DOM_FCT, DOM_TV} typedef type_domain;
type_domain Type_domain = DOM_IWT;

// Algorithms
enum {ALG_ISTA, ALG_FISTA, ALG_MCA, ALG_NESTEROV} typedef type_algo;
type_algo Type_algo = ALG_FISTA;

// Operators
enum {OP_NONE, OP_CONV, OP_MASK} typedef type_op;
type_op Type_op = OP_CONV;

// Fista/ista/mca parameters
bool decrease=false;

// Nesterov parameters
float cont=1;

// Common parameters
float SigmaNoise=1;
bool Positivity=false;
float mu = -1;
float TolVar=-1;
bool no_coarse = false;
bool use_init=false;

fltarray Mask;

/***************************************/

static void usage(char *argv[])
{
	fprintf(OUTMAN, "Usage: %s options data mask result\n\n", argv[0]);
	fprintf(OUTMAN, "   where options =  \n");

	fprintf(OUTMAN, "         [--verbose][-v] Verbose\n");
	manline();

	fprintf(OUTMAN, "         [-a algorithm] Choose the algorithm\n");
	fprintf(OUTMAN, "             1 : ISTA (Iterative Soft Thresholding)\n");
	fprintf(OUTMAN, "            [2]: FISTA (Fast ISTA)\n");
	fprintf(OUTMAN, "             3 : MCA (Morphological Component Analysis)\n");
	fprintf(OUTMAN, "             4 : NESTEROV\n");
	manline();

	fprintf(OUTMAN, "         [-t transform] Choose the representation\n");
	fprintf(OUTMAN, "             1 : Decimated Orthogonal Wavelets\n");
	fprintf(OUTMAN, "             2 : Undecimated Orthogonal Wavelets\n");
	fprintf(OUTMAN, "            [3]: Isotropic Undecimated Wavelets\n");
	fprintf(OUTMAN, "             4 : Fast Curvelets\n");
	fprintf(OUTMAN, "             5 : Total Variation\n");
	manline();

	fprintf(OUTMAN, "         [--op operator] Choose the operator of degradation\n");
	fprintf(OUTMAN, "             0 : None (identity)\n");
	fprintf(OUTMAN, "            [1]: Convolution\n");
	fprintf(OUTMAN, "             2 : Missing data\n");
	manline();

	fprintf(OUTMAN, "         [-n NbrScale]\n");
	fprintf(OUTMAN, "             Number of scales of the transform\n");
	manline();

	fprintf(OUTMAN, "         [-g sigma][--sigmanoise sigma]\n");
	fprintf(OUTMAN, "             sigma = noise standard deviation\n");
	manline();

	fprintf(OUTMAN, "         [-i number][--Niter number]\n");
	fprintf(OUTMAN, "             Number of successive thresholdings\n");
	manline();

	fprintf(OUTMAN, "         [--simu] Simulation : generate the data \n");
	fprintf(OUTMAN, "             by applying the degradation operator to the input\n");
	manline();

	fprintf(OUTMAN, "         [-P][--positive] Positivity constraint\n");
	manline();

	fprintf(OUTMAN, "         [--mu mu] \n");
	fprintf(OUTMAN, "             The value of the step. default=1\n");
	manline();

	fprintf(OUTMAN, "         [--tolvar tv] \n");
	fprintf(OUTMAN, "             The stop criterion.\n");
	fprintf(OUTMAN, "             Stop when |x(n+1)-x(n)|_infty < tv\n");
	fprintf(OUTMAN, "             Default = 1e-6\n");
	manline();

	fprintf(OUTMAN, "         [--decrease] Set to use a decreasing threshold (MCA)\n");
	manline();

	fprintf(OUTMAN, "         [--use-init] Set to use the input as an initial point.\n");
	manline();

	fprintf(OUTMAN, "         [--no-coarse] To substract the coarse scale \n");
	fprintf(OUTMAN, "             before applying the algorithm\n");
	manline();

	exit(-1);
}
 
/*********************************************************************/

/* GET COMMAND LINE ARGUMENTS */
static void filtinit(int argc, char *argv[])
{
	int c,plop;

	/* get options */
	while ((c = GetOpt(argc,argv,(char*)"vn:g:i:Pt:a:")) != -1) 
	{
		switch (c) 
		{
			case 'v': 
				Verbose = True;
				break;
			case 'n': 
				if (sscanf(OptArg,"%d",&NbrScale) != 1) 
				{fprintf(OUTMAN, "Error: bad number of scales parameter: %s\n", OptArg);exit(-1); }
				break;
			case 'g': 
				if (sscanf(OptArg,"%f",&SigmaNoise) != 1) 
				{fprintf(OUTMAN, "Error: bad noise parameter: %s\n", OptArg);exit(-1);}
				break;
			case 'i': 
				if (sscanf(OptArg,"%d",&Niter) != 1) 
				{fprintf(OUTMAN, "Error: bad number of N_iteration: %s\n", OptArg);exit(-1);}
				break;
			case 'P': 
				Positivity = true;
				break;
			case 'a': 
				if (sscanf(OptArg,"%d",&plop) != 1) 
				{fprintf(OUTMAN, "Error: bad number of algorithm: %s\n", OptArg);exit(-1);}
				Type_algo = (type_algo) (plop-1);
				break;
			case 't': 
				if (sscanf(OptArg,"%d",&plop) != 1) 
				{fprintf(OUTMAN, "Error: bad number of transform: %s\n", OptArg);exit(-1);}
				Type_domain = (type_domain) (plop-1);
				break;
			default : 
				usage(argv);
				break;
		}
	} 


	/* get optional input file names from trailing 
	parameters and open files */
	if (OptInd < argc)
		Name_Imag_In = strdup(argv[OptInd++]);
	else
	{
		cerr<<"Not enough parameters : missing Name_In Name_mask Name_Out"<<endl;
		usage(argv);
	}

	if (OptInd < argc)
		Name_Mask = strdup(argv[OptInd++]);
	else
	{
		cerr<<"Not enough parameters : missing Name_mask Name_Out"<<endl;
		usage(argv);
	}

	if (OptInd < argc)
		Name_Imag_Out = strdup(argv[OptInd++]);
	else
	{
		cerr<<"Not enough parameters : missing Name_Out"<<endl;
		usage(argv);
	}

	/* make sure there are not too many parameters */
	if (OptInd < argc)
	{
		fprintf(OUTMAN, "Too many parameters: %s ...\n", argv[OptInd]);
		exit(-1);
	}
}


/***************/

void degrade_mask(fltarray &x, bool degrade, bool reverse)
{
	if(degrade) x *= Mask;
}

int mirror(int i, int n)
{
	if(i<0) return -i;
	if(i>n-1) return 2*n-2-i;
	return i;
}
void degrade_conv2d(fltarray &x, bool degrade, bool reverse)
{
	int Nl=x.nl();
	int Nc=x.nc();
	int nl0=Mask.nl()/2;
	int nc0=Mask.nc()/2;
	int nln=(Mask.nl()-1)/2;
	int ncn=(Mask.nc()-1)/2;
	
//	cerr<<nl0<<":"<<nln<<":"<<nc0<<":"<<ncn<<endl;
	if(degrade)
	{
		fltarray y(Nl,Nc);
		y.init(0.);
		for(int i=nl0;i<Nl-nln;i++)
		for(int j=nc0;j<Nc-ncn;j++)
			for(int di=-nl0;di<=nln;di++)
			for(int dj=-nc0;dj<=ncn;dj++)
				y(i,j) += x(i+di,j+dj)*Mask(nl0+di,nc0+dj);
// Periodic conditions
/*		for(int i=0;i<nl0;i++)
		for(int j=0;j<Nc;j++)
			for(int di=-nl0;di<=nln;di++)
			for(int dj=-nc0;dj<=ncn;dj++)
				y(i,j) += x((i+di)%Nl,(j+dj)%Nc)*Mask(nl0+di,nc0+dj);
		for(int i=Nc-ncn;i<Nc;i++)
		for(int j=0;j<Nc;j++)
			for(int di=-nl0;di<=nln;di++)
			for(int dj=-nc0;dj<=ncn;dj++)
				y(i,j) += x((i+di)%Nl,(j+dj)%Nc)*Mask(nl0+di,nc0+dj);
		for(int i=nl0;i<Nl-nl0;i++)
		for(int j=0;j<nc0;j++)
			for(int di=-nl0;di<=nln;di++)
			for(int dj=-nc0;dj<=ncn;dj++)
				y(i,j) += x((i+di)%Nl,(j+dj)%Nc)*Mask(nl0+di,nc0+dj);
		for(int i=nl0;i<Nl-nc0;i++)
		for(int j=Nc-ncn;j<Nc;j++)
			for(int di=-nl0;di<=nln;di++)
			for(int dj=-nc0;dj<=ncn;dj++)
				y(i,j) += x((i+di)%Nl,(j+dj)%Nc)*Mask(nl0+di,nc0+dj);
*/
// Mirror conditions
		for(int i=0;i<nl0;i++)
		for(int j=0;j<Nc;j++)
			for(int di=-nl0;di<=nln;di++)
			for(int dj=-nc0;dj<=ncn;dj++)
				y(i,j) += x(mirror(i+di,Nl),mirror(j+dj,Nc))*Mask(nl0+di,nc0+dj);
		for(int i=Nc-ncn;i<Nc;i++)
		for(int j=0;j<Nc;j++)
			for(int di=-nl0;di<=nln;di++)
			for(int dj=-nc0;dj<=ncn;dj++)
				y(i,j) += x(mirror(i+di,Nl),mirror(j+dj,Nc))*Mask(nl0+di,nc0+dj);
		for(int i=nl0;i<Nl-nl0;i++)
		for(int j=0;j<nc0;j++)
			for(int di=-nl0;di<=nln;di++)
			for(int dj=-nc0;dj<=ncn;dj++)
				y(i,j) += x(mirror(i+di,Nl),mirror(j+dj,Nc))*Mask(nl0+di,nc0+dj);
		for(int i=nl0;i<Nl-nc0;i++)
		for(int j=Nc-ncn;j<Nc;j++)
			for(int di=-nl0;di<=nln;di++)
			for(int dj=-nc0;dj<=ncn;dj++)
				y(i,j) += x(mirror(i+di,Nl),mirror(j+dj,Nc))*Mask(nl0+di,nc0+dj);
		x = y;
	}
}

void degrade_other(fltarray &x, bool degrade, bool reverse)
{
	int Nl=x.nl();
	int Nc=x.nc();
	int shift=10;
	
//	cerr<<nl0<<":"<<nln<<":"<<nc0<<":"<<ncn<<endl;
	if(degrade)
	{
		fltarray y(Nl,Nc);
		y.init(0.);
		for(int i=0;i<Nl;i++)
		{
			for(int j=0;j<Nc;j++)
				y(i,j) += (sqrt(2))*x(i,j);
			for(int j=0;j<Nc-shift;j++)
				y(i,j) -= (sqrt(2)/2)*x(i,j+shift);
			for(int j=shift;j<Nc;j++)
				y(i,j) -= (sqrt(2)/2)*x(i,j-shift);
		}
	}
	else
	{
		fltarray y(Nl,Nc);
		y.init(0.);
		for(int i=0;i<Nl;i++)
		{
			for(int j=0;j<Nc;j++)
				y(i,j) += (sqrt(2))*x(i,j);
			for(int j=0;j<Nc-shift;j++)
				y(i,j) -= (sqrt(2)/2)*x(i,j+shift);
			for(int j=shift;j<Nc;j++)
				y(i,j) -= (sqrt(2)/2)*x(i,j-shift);
		}
	}
}

void degrade_none(fltarray &x, bool degrade, bool reverse)
{
}


int main(int argc, char *argv[])
{
// Get long options
	map<string, string> opts;
	GetLongOptions(argv, argc, opts);
	int plop;
	
	GetLongOption(opts,"verbose",Verbose);
	GetLongOption(opts,"op",plop) && (Type_op=(type_op)plop);
	GetLongOption(opts,"simu",Simu);
	GetLongOption(opts,"SigmaNoise",SigmaNoise);
	GetLongOption(opts,"Niter",Niter);
	GetLongOption(opts,"positive",Positivity);
	GetLongOption(opts,"decrease",decrease);
	GetLongOption(opts,"mu",mu);
	GetLongOption(opts,"tolvar",TolVar);
	GetLongOption(opts,"use-init",use_init);
	GetLongOption(opts,"no-coarse",no_coarse);

// Get short options
	filtinit(short_argc, short_argv);

// Print parameters
	{
		cout << endl << endl << "PARAMETERS: " << endl << endl;
		if(use_init) cout << "Using input as initial point" << endl;

		cout << "File Name in = " << Name_Imag_In << endl;
		cout << "File Name Out = " << Name_Imag_Out << endl;   
		cout << "NbrScale = " <<  NbrScale    << endl;  

		// Filtering
		if(SigmaNoise>0)
		{
								cout << " Noise level = " <<  SigmaNoise    << endl;
			if(Niter>1)			cout << " Number of iterations : " <<  Niter    << endl;
			if(Positivity)	cout << " The reconstruction is positive" << endl;
//			if(use_min)		cout << " The reconstruction is > " <<P.min_value<< endl;
//			if(use_max)		cout << " The reconstruction is < " <<P.max_value<< endl;
		}
		cout << endl;  
	}

// Input data
	fltarray Data;
	fltarray Recons;
	Name_Imag_In = add_fits(Name_Imag_In);
	fits_read_fltarr(Name_Imag_In, Data);
	Name_Mask = add_fits(Name_Mask);
	fits_read_fltarr(Name_Mask, Mask);

/*
// UWT initialisation
	UOWT *Domain;
	FilterAnaSynt SelectFilter(F_MALLAT_7_9);
	SubBandFilter *SB1D = new SubBandFilter(SelectFilter, NORM_L2);
	Domain = new UOWT(SB1D);
	Domain->alloc(Data,NbrScale);

// FCurvelet initialisation
	BFCurvelet3D *Domain = new BFCurvelet3D;
	BFCurvelet3D_params P;
	P.Positivity = Positivity;
	P.NbrScale = NbrScale;
	P.SigmaNoise = SigmaNoise;
	Domain->alloc_from_coarse(Data,P);
*/	


// Generic transform
	FloatTrans *Domain;

// Specific transforms
	Ortho_2D_WT_float *Domain_dwt;
	PAVE_2D_WT_float *Domain_uwt;
	ATROUS_2D_WT_float *Domain_iwt;
	FCUR_float *Domain_fct;
	PAVE_2D_WT_float *Domain_tv;

// Filters fir the wavelet transforms and TV
	FilterAnaSynt SelectFilter;
	if(Type_domain == DOM_TV)
		SelectFilter.alloc(F_HAAR);
	else
		SelectFilter.alloc(F_MALLAT_7_9);
	SubBandFilter *SB1D = new SubBandFilter(SelectFilter, NORM_L2);

	switch(Type_domain)
	{
		case DOM_DWT:
			Domain_dwt = new Ortho_2D_WT_float(*SB1D);
			Domain_dwt->alloc(Data, NbrScale);
			Domain = Domain_dwt;
			break;
		case DOM_UWT:
			Domain_uwt = new PAVE_2D_WT_float(*SB1D);
			Domain_uwt->alloc(Data, NbrScale);
			Domain = Domain_uwt;
			break;
		case DOM_IWT:
			Domain_iwt = new ATROUS_2D_WT_float();
			Domain_iwt->alloc(Data, NbrScale);
			Domain_iwt->AdjointRec = True;
			//Domain_iwt->set_positive_coef(true);
			Domain = Domain_iwt;
			break;
		case DOM_FCT:
			Domain_fct = new FCUR_float();
			Domain_fct->alloc(Data, NbrScale, NbrDir2d);
			Domain_fct->get_norm_coeff(3);
			Domain = Domain_fct;
			break;
		case DOM_TV:
			Domain_tv = new PAVE_2D_WT_float(*SB1D);
			Domain_tv->alloc(Data, 2);
			Domain = Domain_tv;
			break;
	}

	// Degradation operator
	void (*degrade)(fltarray&,bool,bool);
	switch(Type_op)
	{
		case OP_CONV:
			degrade = &degrade_conv2d;
			break;
		case OP_MASK:
			degrade = &degrade_mask;
			break;
		default:
			degrade = &degrade_other;
			break;
	}
	
	if(Type_algo == ALG_ISTA || Type_algo == ALG_FISTA || Type_algo == ALG_MCA)
	{
	// Fista init
		Fista F(Domain);
		F.Decreasing = (Type_algo == ALG_MCA);
		F.Fast = (Type_algo == ALG_FISTA);
		F.Threshold = SigmaNoise;
		F.MaxNiter = Niter;
		F.No_coarse = no_coarse;
		F.Positivity = Positivity;
		if(mu>0) F.Mu = mu;
		if(TolVar>0) F.TolVar = TolVar;
		strcpy(F.NameOut,Name_Imag_Out);

	// Fista algorithm
		if(Simu) degrade(Data, !use_init, false);// degrade only if not using init
		char filename[128]; sprintf(filename,"%s_data.fits",Name_Imag_Out);writefltarr(filename, Data);
		F.run(Data, Recons, degrade);
	}
	else if(Type_algo == ALG_NESTEROV)
	{
	// Nesterov init
		Nesterov N(Domain);
		N.MaxNiter = Niter;
		N.Continuation = cont;
		N.Threshold = SigmaNoise;
		N.Positivity = Positivity;
		if(mu>0) N.Mu = mu;
		if(TolVar>0) N.TolVar = TolVar;
		strcpy(N.NameOut,Name_Imag_Out);

	// Nesterov algorithm
		if(Simu) degrade(Data, !use_init, false);// degrade only if not using init
		N.run(Data, Recons, degrade);
	}
	else 
		cerr<<"Unknown method"<<endl;
	
// Free all arrays
	for(int i=0;i<short_argc;i++) delete [] short_argv[i];
	delete [] short_argv;
	delete Domain;
	Data.~fltarray();
	Recons.~fltarray();
	exit(0);
}






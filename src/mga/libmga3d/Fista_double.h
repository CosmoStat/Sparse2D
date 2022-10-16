
#ifndef _FISTA_D_H
#define _FISTA_D_H

#include <fstream>
#include "GlobalInc.h"
#include "SparseTrans.h"

extern Bool Verbose;

class Fista_params
{
public:
	void reset();			// Initialise the class with default parameters
public:
	char NameOut[256];
	int MaxNiter;			// Maximum number of iterations
	double Threshold;		// Thresholding level
	double TolVar;			// Stopping criterium : ||xt+1 - xt||_infty < TolVar
	double Mu;				// Gradient regularization parameter
	bool Fast;				// (F)ISTA
	bool Decreasing;		// if true : linearily decreasing threshold
	bool No_coarse;			// if true : kills the coarse scale before apaplying fista, then thresholds it
};
void Fista_params::reset()
{
	MaxNiter=10;
	Threshold = 1;
	TolVar = 1e-1;
	Mu=2./15.;
	Fast = false;
	Decreasing = false;
	No_coarse = false;
}


class Fista
{
public:
	Fista(SparseTrans *domain);
	~Fista(){};
	
	Fista_params P;
	FloatTrans* Domain;
	
	void run(fltarray &data, fltarray &z, void (*_degrade)(fltarray &,bool,bool));
	// data : observed signal 
	// z : recovered signal
	// _degrade : degradation operator. bool to transpose
	void run(cfarray &data, cfarray &z, void (*_degrade)(cfarray &,bool,bool));
	// data : observed signal in fourier space
	// z : recovered signal in fourier space
	// _degrade(z, bool degrade, bool reverse) : degradation operator. bool to transpose
	// _degrade(z,*,true) is assumed real (imaginary part ignored)
};


Fista::Fista(FloatTrans *domain)
{
	P.reset();
	Domain = domain;
}

void Fista::run(fltarray &b, fltarray &z, void (*_degrade)(fltarray &,bool,bool))
{
	char filename[64];
	double *x, *xold, *xtmp; // current estimate, previous estimate, and temp variable

// Initialization
	double mu = P.Mu;
	double lvl = P.Threshold;
	double TolVar = P.TolVar;

// Other variables
	double tk=1,told;
	int i; bool done;
	double speed, old_speed;
	fltarray y;// previous reconstruction
	i=0; done = false; speed=0; old_speed=0;
	bool allocD=true; // Domain allocation
	
// Initial point z=b
	float* CS;
	if(P.No_coarse) // Substract the coarse scale
	{
		_degrade(b, false, true);// degraded to direct space, no degradation made
		Domain->transform(b, xtmp, allocD); allocD=false;// allocate xtmp and localTB
		Domain->substract_coarse_scale(xtmp,CS,true);
		Domain->recons(xtmp, b);// allocate xtmp and localTB
		_degrade(b, false, false);// direct to degraded space, no degradation made
	}
	z = b; // degraded space
	_degrade(z, false, true); // degraded to direct space, no degradation made
	Domain->transform(z, xtmp, allocD);// allocate xtmp and localTB
	int n=Domain->size_transform();
	
	
	x = new float[n]; for(int k=0;k<n;k++) x[k] = xtmp[k];
	xold = new float[n];

	cerr<<"##########\nBegin FISTA with mu="<<mu<<", and "<<P.MaxNiter<<" iterations, with fast="<<P.Fast<<" decrease="<<P.Decreasing<<".\n##########"<<endl;
	
// Fista algorithm. 
	while( i<P.MaxNiter && (!done || P.Decreasing) )
	{
	// Threshold update
		if(P.Decreasing) 
		{
			//lvl = P.Threshold * (P.Niter-i-1.)/(P.Niter-1.);
			lvl = 1-i*1.0/float(P.MaxNiter-1);			// 1..0
			lvl = (pow((double)lvl,5)+lvl/100.)/1.01;
			lvl = P.Threshold*(1+ lvl*(1000-1));
		}
		
	// Save current solution
		if(i!=0) Domain->recons(xtmp,z);
		sprintf(filename,"%s_%05d.fits",P.NameOut,i);writefltarr(filename, z);
		
	// Evaluate the evolution
		y = y-z; // error between the last two estimates
		speed = abs(y.maxfabs());
		done = (speed < TolVar) && (i>0) && (speed<old_speed);
		old_speed=speed;
		cerr<<" Step "<<i<<", lvl="<<lvl<<", || z - zt ||_infty = "<<abs(speed)<<", TolVar="<<TolVar<<endl;
		y=z; // Save the new solution
		
	// Gradient step
		_degrade(z, true, false);
		z = b - z ;
		_degrade(z, true, true);
		Domain->transform(z, xold);
		for(int k=0;k<n;k++) xtmp[k] = xtmp[k] + mu*xold[k];
		
	// Save estimate
		if(P.Fast) for(int k=0;k<n;k++) xold[k] = x[k];
		
	// Proximal operator : soft thresholding
		for(int k=0;k<n;k++) x[k] = xtmp[k];
		Domain->soft_threshold(x,lvl,P.No_coarse);
		
	// New point
		if(P.Fast) 
		{
			told=tk;
			tk = (1.+sqrt(1.+4.*told*told))/2.;
			for(int k=0;k<n;k++) xtmp[k] = x[k] + (told-1)/tk * (x[k]-xold[k]);
		}
		else
			for(int k=0;k<n;k++) xtmp[k] = x[k];
		
		i++;
	}
	if(P.No_coarse)
		Domain->add_coarse_scale(x,CS);
	Domain->recons(x,z);
	sprintf(filename,"%s_recons.fits",P.NameOut);writefltarr(filename, z);
}



#endif




#ifndef _FISTA_H
#define _FISTA_H

#include <fstream>
#include "GlobalInc.h"
#include "FloatTrans.h"

extern Bool Verbose;

class Fista_params
{
public:
	void reset();			// Initialise the class with default parameters
public:
	char NameOut[256];
	int MaxNiter;			// Maximum number of iterations
	float Threshold;		// Thresholding level
	float TolVar;			// Stopping criterium : ||xt+1 - xt||_infty < TolVar
	float Mu;				// Gradient regularization parameter
	bool Fast;				// (F)ISTA
	bool Decreasing;		// if true : linearily decreasing threshold
};
void Fista_params::reset()
{
	MaxNiter=10;
	Threshold = 1;
	TolVar = 1e-1;
	Mu=2./15.;
	Fast = false;
	Decreasing = false;
}


class Fista
{
public:
	Fista(FloatTrans *domain);
	~Fista(){};
	
	Fista_params P;
	FloatTrans* Domain;
	
	void run(fltarray &data, fltarray &z, void (*_degrade)(fltarray &,bool));
	// data : observed signal 
	// z : recovered signal
	// _degrade : degradation operator. bool to transpose
	void run(cfarray &data, cfarray &z, void (*_degrade)(cfarray &,bool));
	// data : observed signal in fourier space
	// z : recovered signal in fourier space
	// _degrade : degradation operator. bool to transpose
	// _degrade(z,true) is assumed real (imaginary part ignored)
};


Fista::Fista(FloatTrans *domain)
{
	P.reset();
	Domain = domain;
}

void Fista::run(fltarray &b, fltarray &z, void (*_degrade)(fltarray &,bool))
{
	char filename[64];
	float *x, *xold, *xtmp; // current estimate, previous estimate, and temp variable

// Initialization
	float mu = P.Mu;
	float lvl = P.Threshold;
	float TolVar = P.TolVar;

// Initial point
	z = b;
	Domain->transform(z, xtmp, true);// allocate xtmp and localTB
	int n=Domain->size_transform();
	x = new float[n]; for(int k=0;k<n;k++) x[k] = x[k];
	xold = new float[n];

// Other variables
	float tk=1,told;
	int i; bool done;
	float speed, old_speed;
	fltarray y;// previous reconstruction
	i=0; done = false; speed=0; old_speed=0;
	
	cerr<<"##########\nBegin FISTA with mu="<<mu<<", and "<<P.MaxNiter<<" iterations, with fast="<<P.Fast<<" decrease="<<P.Decreasing<<".\n##########"<<endl;
	
// Fista algorithm. 
	while( i<P.MaxNiter && (!done || P.Decreasing) )
	{
	// Threshold update
		if(P.Decreasing) 
		{
			//lvl = P.Threshold * (P.Niter-i-1.)/(P.Niter-1.);
			
			lvl = 1-i*1.0/float(P.MaxNiter-1);			// 1..0
//			lvl = (pow((double)lvl,3)+lvl/25.)/1.04;
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
		_degrade(z, false);
		z = b - z ;
		_degrade(z,true);
		Domain->transform(z, xold);
		for(int k=0;k<n;k++) xtmp[k] = xtmp[k] + mu*xold[k];
		
	// Save estimate
		if(P.Fast) for(int k=0;k<n;k++) xold[k] = x[k];
		
	// Proximal operator : soft thresholding
		for(int k=0;k<n;k++) x[k] = xtmp[k];
		Domain->soft_threshold(x,lvl);
		
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
	Domain->recons(x,z);
}

// b : degraded data, in degraded domain
// z : output data, in degraded domain
// _degrade(z,bool reverse) : degradation operator
void Fista::run(cfarray &b, cfarray &z, void (*_degrade)(cfarray&,bool))
{
	char filename[64];
	float *x, *xold, *xtmp; // current estimate, previous estimate, and temp variable
	fltarray zreal;
	
// Initialization
	float mu = P.Mu;
	float lvl = P.Threshold;
	float TolVar = P.TolVar;

// Initial point
	z = b; // degraded space
	_degrade(z, true); // degraded to direct space
	zreal.alloc(z.nx(),z.ny(),z.nz());
	if(z.nz()) for(int k=0;k<z.nz();k++) for(int j=0;j<z.ny();j++) for(int i=0;i<z.nx();i++) zreal(i,j,k) = z(i,j,k).real();
	else for(int j=0;j<z.ny();j++) for(int i=0;i<z.nx();i++) zreal(j,i) = z(j,i).real();
	Domain->transform(zreal, xtmp, true);// allocate xtmp and localTB
	int n=Domain->size_transform();
	x = new float[n]; for(int k=0;k<n;k++) x[k] = x[k];
	xold = new float[n];

// Other variables
	float tk=1,told;
	int iter; bool done;
	float speed, old_speed;
	fltarray yreal;// previous reconstruction
	iter=0; done = false; speed=0; old_speed=0;
	
	cerr<<"##########\nBegin FISTA with mu="<<mu<<", and "<<P.MaxNiter<<" iterations, with fast="<<P.Fast<<" decrease="<<P.Decreasing<<".\n##########"<<endl;
// Fista algorithm. 
	while( iter<P.MaxNiter && (!done || P.Decreasing) )
	{
	// Threshold update
		if(P.Decreasing) 
		{
			//lvl = P.Threshold * (P.Niter-i-1.)/(P.Niter-1.);
			
			lvl = 1-iter*1.0/float(P.MaxNiter-1);			// 1..0
			lvl = (pow((double)lvl,3)+lvl/25.)/1.04;
			lvl = P.Threshold*(1+ lvl*(1000-1));
		}
		
	// Save current solution
		if(iter!=0) Domain->recons(xtmp,zreal);
		if(z.nz()) for(int k=0;k<z.nz();k++) for(int j=0;j<z.ny();j++) for(int i=0;i<z.nx();i++) z(i,j,k) = complex_f(zreal(i,j,k),0.);
		else for(int j=0;j<z.ny();j++) for(int i=0;i<z.nx();i++) z(j,i) = complex_f(zreal(j,i),0.);
		sprintf(filename,"%s_%05d.fits",P.NameOut,iter);writefltarr(filename, zreal);
		
	// Evaluate the evolution
		yreal = yreal - zreal; // error between the last two estimates
		speed = abs(yreal.maxfabs());
		done = (speed < TolVar) && (iter>0) && (speed<old_speed);
		old_speed=speed;
		cerr<<" Step "<<iter<<", lvl="<<lvl<<", || z - zt ||_infty = "<<abs(speed)<<", TolVar="<<TolVar<<endl;
		yreal = zreal; // Save the new solution
		
	// Gradient step
		_degrade(z, false);// direct to degraded space
		z = b - z ;// degraded space
		_degrade(z,true); // degraded to direct space
		if(z.nz()) for(int k=0;k<z.nz();k++) for(int j=0;j<z.ny();j++) for(int i=0;i<z.nx();i++) zreal(i,j,k) = z(i,j,k).real();
		else for(int j=0;j<z.ny();j++) for(int i=0;i<z.nx();i++) zreal(j,i) = z(j,i).real();
		Domain->transform(zreal, xold);
		for(int k=0;k<n;k++) xtmp[k] = xtmp[k] + mu*xold[k];
		
	// Save estimate
		if(P.Fast) for(int k=0;k<n;k++) xold[k] = x[k];
		
	// Proximal operator : soft thresholding
		for(int k=0;k<n;k++) x[k] = xtmp[k];
		Domain->soft_threshold(x,lvl);
		
	// New point
		if(P.Fast) 
		{
			told=tk;
			tk = (1.+sqrt(1.+4.*told*told))/2.;
			for(int k=0;k<n;k++) xtmp[k] = x[k] + (told-1)/tk * (x[k]-xold[k]);
		}
		else
			for(int k=0;k<n;k++) xtmp[k] = x[k];
		
		iter++;
	}
	Domain->recons(x,zreal);
	for(int i=0;i<z.nx();i++) for(int j=0;j<z.ny();j++) for(int k=0;k<z.nz();k++) z(i,j,k) = complex_f(zreal(i,j,k),0.);
}


#endif



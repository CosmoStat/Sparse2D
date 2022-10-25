
#ifndef _MCA_H
#define _MCA_H

#include <fstream>
#include "GlobalInc.h"
#include "FloatTrans.h"
#include "writefits3d.h"

extern Bool Verbose;

class MCA_params
{
public:
	void reset();			// Initialise the class with default parameters
	char NameOut[256];
	int MaxNiter;			// Maximum number of iterations
	float Threshold;		// Thresholding level
	float TolVar;			// Stopping criterium : ||xt+1 - xt||_infty < TolVar
	float Mu;				// Gradient regularization parameter
	bool Positivity;		// The reconstruction must be >0
	bool No_coarse;			// if true : kills the coarse scale before apaplying fista, then thresholds it
    float ksigma;           // threshold level = k * SigmaNoise.
    double SigmaNoise;      // Noise standard deviation:  lambda = ksigma * SigmaNoise
};


void MCA_params::reset()
{
	MaxNiter=10;
	Threshold = 1;
	TolVar = 1e-4;
	Mu=1;
	Positivity = false;
	No_coarse = false;
    ksigma=3.;
    SigmaNoise=0.;
}


class MCA : public MCA_params
{
public:
	MCA(FloatTrans *domain);
	~MCA(){};
	
	FloatTrans* Domain;
	
	void run(cfarray &b, cfarray &z, void (*_degrade)(cfarray &,bool,bool));
	// b : observed signal, in degraded domain
	// z : recovered signal, in degraded domain
	// _degrade : degradation operator. 1bool to apply the degradation, 2bool to transpose
	


};

MCA::MCA(FloatTrans *domain)
{
	reset();
	Domain = domain;
}

void MCA::run(cfarray &b, cfarray &z, void (*_degrade)(cfarray&,bool,bool))
{
    float lvl;
    bool done, allocD;
    float speed, old_speed, top_level, array_max;
    int i;
    float *wavelet_coeff;
    float *CS;
    char filename[64];
    fltarray yreal, zreal, zr;
	
    lvl = Threshold;  TolVar = TolVar; 
    i = 0; done = false; speed = 0; old_speed = 0; allocD = true;
  
	if(Verbose) cout << "MCA params " <<  MaxNiter << " " << Threshold << " " << TolVar << " " <<  Mu << " " <<  Positivity << " " << No_coarse << endl;
	
    // Initial point
  
    zreal.alloc(b.nx(),b.ny(),b.nz());
    zr.alloc(Domain->size_transform());
     
	if(No_coarse)   // turn the FFT back into an image, then wavelet space, then fiddle with wavelets, then back to Fourier space
	{
		_degrade(b, false, true);	
		for(int k=0;k<b.n_elem();k++) zreal(k) = b(k).real();
		Domain->transform(zreal, wavelet_coeff, allocD); allocD=false; // allocate wavelet_coeff and localTB
		Domain->substract_coarse_scale(wavelet_coeff,CS,true);
		Domain->recons(wavelet_coeff, zreal);// allocate wavelet_coeff and localTB
		for(int k=0;k<b.n_elem();k++) b(k) = complex_f(zreal(k), 0.);
		_degrade(b, false, false);
	}

	z = b;
	_degrade(z, false, true);
 	for(int k=0;k<z.n_elem();k++) zreal(k) = z(k).real();

    // Find top level
	Domain->transform(zreal,wavelet_coeff,allocD); 
    Domain->normalize(wavelet_coeff);
	array_max = -1e39;
	for (int i=0; i<Domain->size_transform()-z.nx()*z.ny(); ++i) if ( wavelet_coeff[i] > array_max ) array_max = wavelet_coeff[i];
        top_level = 0.99*array_max;   
        cout << "Calculated top level: " << top_level  << endl << endl;   	

	
	if(Verbose) cerr<<"##########\nBegin ";
	if(Verbose) cerr<<"MCA";
	if(Verbose) cerr<<" with mu="<<Mu<<", and "<<MaxNiter<<" iterations.\n##########"<<endl;

// MCA   algorithm. 
    zreal.init();
	while( i<MaxNiter && !done )
	{
    // zreal.info();
    
	// Threshold update
		//lvl = P.Threshold * (P.Niter-i-1.)/(P.Niter-1.);
		lvl = 1-i*1.0/float(MaxNiter-1);			// 1..0
		// lvl = (pow((double)lvl,3)+lvl/25.)/1.04;
		// lvl = Threshold*(1+ lvl*(top_level*100-1));
	       
		
	 // Apply positivity and save
		
     if(Positivity) for(int k=0;k<zreal.n_elem();k++) zreal(k) = zreal(k) *(zreal(k)>0);
     // sprintf(filename,"%s_%05d.fits",NameOut,i); 
     // writefltarr(filename, zreal);

		 
	// Evaluate the evolution. Compare z to y, which holds the previous image
		yreal = yreal - zreal; // error between the last two estimates
		speed = abs(yreal.maxfabs()/zreal.maxfabs());
		done = (speed < TolVar) && (i>0) && (speed<old_speed);
		old_speed=speed;
		if(Verbose) cerr<<" Step "<<i<<", lvl="<<lvl<<", || z - zt ||_infty = "<<abs(speed)<<", TolVar="<<TolVar<<endl; 
		yreal = zreal; // Save the current solution image zreal to yreal
	
	// Gradient step
        for(int k=0;k<z.n_elem();k++) z(k) = complex_f(zreal(k),0);
        _degrade(z, true, false);
		 z = b - z;    
                  
    // Get Fourier deltas and create image deltas and add to image with mu
        _degrade(z,true, true);  
          
    // Calculate the MAD estimator per band
        if (Domain-> UseMad == true) 
        {
            // wavelets transform of the residual
            for (int k=0;k<zreal.n_elem();k++) { zr(k) = Mu * z(k).real(); }
            Domain->transform(zr, wavelet_coeff); 
            
            // cout << "MAD calc " << Domain->size_transform() << endl;
            Domain->mad_calculation(wavelet_coeff);
            // cout << "MAD ok calc " <<   endl;
        }                 
        // zreal.info();                                                                   		 
		for(int k=0;k<zreal.n_elem();k++) { zreal(k) += Mu*z(k).real(); }
	
	// Proximal operator : hard thresholding using wavelets
		Domain->transform(zreal, wavelet_coeff); 
        
        
        // zreal.info("SOL");
        if (Domain-> UseMad == true)
        {
           double Threshold = ksigma + lvl;
           Domain->mad_normalize(wavelet_coeff);
           Domain->hard_threshold(wavelet_coeff,Threshold,No_coarse);
           Domain->mad_unnormalize(wavelet_coeff);
           // Domain->recons(wavelet_coeff, zreal);
           Domain->adjoint_recons(wavelet_coeff, zreal);
        }
        else
        {
           double Threshold = (ksigma+lvl) * SigmaNoise;
           Domain->normalize(wavelet_coeff);
           Domain->hard_threshold(wavelet_coeff,Threshold, No_coarse);
           Domain->unnormalize(wavelet_coeff);
           Domain->recons(wavelet_coeff, zreal);
        }
		
		i++;   
	}
	
	if (No_coarse) {
	  Domain->transform(zreal,wavelet_coeff);
      Domain->normalize(wavelet_coeff);
      Domain->unnormalize(wavelet_coeff);
      Domain->add_coarse_scale(wavelet_coeff,CS);
	  Domain->recons(wavelet_coeff,zreal);
	}
	
	
	if(Positivity) for(int k=0;k<zreal.n_elem();k++) zreal(k) = zreal(k) *(zreal(k)>0);
	// sprintf(filename,"%s.fits",NameOut);
    // writefltarr(filename, zreal);
    fits_write_fltarr(NameOut, zreal);
	cerr<<"##########\nEnd.\n##########"<<endl;
	
	
}

#endif



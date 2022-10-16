//******************************************************************************
//**
//**    DESCRIPTION  Few Event Poisson Thresholding 
//**    -----------  
//**                 
//******************************************************************************
 

#include "Mr3d_FewEvent.h"

//******************************************************************************
FewEventPoisson::FewEventPoisson ( Bool InitOk, int NAutoConv, float epsilon) {
//	initialisation of datas
   _Epsilon=epsilon;
   _Param.alloc(2);
   _InitOk = InitOk;
   if(_InitOk == True){   
   	fits_read_fltarr(Name_Param, _Param);
	_NbAutoConv = (int)(_Param(1)+0.5);
   }else{
   	_NbAutoConv= NAutoConv;
   }
   _HistoConv.alloc (_NbAutoConv+1, MAXHISTONBPOINT);
   _HistoConv.init(0.);
   _HistoBound.alloc (_NbAutoConv+1, 2);
   _HistoBin.alloc (_NbAutoConv+1, 2); 
   _HistoDistrib.alloc (_NbAutoConv+1, MAXHISTONBPOINT); 
   _Threshold.alloc(_NbAutoConv+1, 2, ""); 
}



//******************************************************************************
void FewEventPoisson::compute_distribution (Bool WriteAllInfo, Bool WriteHisto, 
		Bool Verbose, int param) {
	//  compute histogramms distributions and  the values of threshold

   // compute Bspline Histogramm
   if  (Verbose==True) 
      cout << "Compute the histogram of the Wavelet ... " << endl;
   bspline_histo_3D (WriteAllInfo, Verbose); 
  
   // computes the auto-convolued histogramms
   if  (Verbose==True) cout << "Compute the autoconvolutions of the histogram ... " << endl;
   histo_convolution (WriteAllInfo,Verbose,param);
   
   
   // computes the distribution function of histogramms at each leve
   histo_distribution (WriteAllInfo,Verbose);
  
   _Param(0)=param;
   _Param(1)=_NbAutoConv;
   
   // writing debug result
   if ((WriteHisto)||(WriteAllInfo==True)) {
      fits_write_fltarr (Name_HistoConv, _HistoConv);
      fits_write_fltarr (Name_HistoBound, _HistoBound);      
      fits_write_fltarr (Name_HistoBin, _HistoBin);
      fits_write_fltarr(Name_HistoDistrib, _HistoDistrib);
      fits_write_fltarr(Name_Param, _Param);
   } 
   
   _InitOk = True;
}
 


//******************************************************************************
void FewEventPoisson::bspline_histo_3D (Bool WriteAllInfo, Bool Verbose) {
	//	create the histogramm of the wavelet B3Spline
   
   // init Phi cubic bspline
   // performs B(x) for x in [-2;+2] :
   // B(x) =1/12 (|x-2|^3 - 4 |x-1|^3 + 6 |x|^3 - 4 |x+1|^3 + |x+2| ^ 3)
   fltarray phi_bspline(BIN3D_2);
   float x;
   for (int i=-BIN3D/2;i<=BIN3D/2;i++) {
   
      x = 2 * ((float)i)/((float)BIN3D/2);
      phi_bspline(i+BIN3D) = (      CUBE((ABS(x-2))) - 4 * CUBE((ABS(x-1)))
			      + 6 * CUBE((ABS(x)))   - 4 * CUBE((ABS(x+1))) 
			      +     CUBE((ABS(x+2))) ) / 12;
   }      
   
   if (WriteAllInfo==True) {
	fits_write_fltarr (Name_Bspline_, (phi_bspline));
   	if (Verbose==True) cout << "File  "<< Name_Bspline_ << "  is created  (Bspline)"<< endl;  
   }
   
   // performs Psi wavelet 
   // psi(x,y,z) = B(x,y,z) - 1/8 B(x/2,y/2,z/2)
   // psi(x,y) = B(x) B(y) B(z) - 1/8 B(x/2) B(y/2) B(z/2)
   
   
   float MaxPsi = 0;
   float MinPsi = 0;
   float val;
   for (int i=-BIN3D;i<=BIN3D;i++) {
      for (int j=-BIN3D;j<= BIN3D;j++) 
      	for (int k=-BIN3D;k<= BIN3D;k++) {
		val=phi_bspline(i+BIN3D) * phi_bspline(j+BIN3D) * phi_bspline(k+BIN3D) 
			- 0.125 * phi_bspline(i/2+BIN3D) * phi_bspline(j/2+BIN3D) * phi_bspline(k/2+BIN3D);
         	if (val >= MaxPsi) 
	    		MaxPsi = val;
         	if (val <= MinPsi) 
	    		MinPsi = val;
      }  
   }
 
   // Fullfills the histogramm (first length = 1025 points)
   // 
   int index;   
   int LocalSizeSignal = LENGTH_FIRST_HISTO;
   for (int i=-BIN3D;i<=BIN3D;i++) {
      for (int j=-BIN3D;j<=BIN3D;j++) {
   	 for (int k=-BIN3D;k<=BIN3D;k++) {
            // we only use (HSP-1)/2 point in _HistoConv tab
	    val=phi_bspline(i+BIN3D) * phi_bspline(j+BIN3D) * phi_bspline(k+BIN3D) 
			- 0.125 * phi_bspline(i/2+BIN3D) * phi_bspline(j/2+BIN3D) * phi_bspline(k/2+BIN3D);
            index = (int) (  (val-MinPsi)*(LocalSizeSignal-1)
                        / (MaxPsi-MinPsi));
            _HistoConv(0,index)++;  // first histogram, no convolution
	 }
      } 
   }
   
   // bin
   float LocalHistoBin   = (MaxPsi-MinPsi) / (LocalSizeSignal-1);         
   
   // set min, max, bin and nb points of first level (histogram only)
   //
   _HistoBound(0,0) = MinPsi;
   _HistoBound(0,1) = MaxPsi;
   _HistoBin (0,0) = LocalHistoBin;
   _HistoBin (0,1) = LocalSizeSignal;   
   
   float sum_check=0.0;
   for(int i=0;i<_HistoBin(0,1);i++)
      	sum_check += _HistoConv(0,i)*_HistoBin (0,0);
   // it must be a probability density -> normalisation by sum_check
   for(int i=0;i<_HistoBin(0,1);i++) 
	_HistoConv(0,i) /= sum_check;
   
   
   
   // normalisation of the first histogram
   
   //Reduction to have a sigma = 1 and mean = 0
   //Dilatation on x axes
   // the histogramme must have a mean=0, and the bin is divided by sigma (std)
   // calculation of the mean and the sigma (std)
   float mean=0.0;
   float std=0.0;
      	for (int j=0; j<_HistoBin(0,1);j++) {
         	float real_value = ((float) j)*_HistoBin (0,0) + _HistoBound(0,0);
         	mean += _HistoConv(0,j)*_HistoBin (0,0)*real_value;
      	}

      	for(int j=0; j<_HistoBin(0,1);j++){
         	float real_value = ((float) j)*_HistoBin (0,0) + _HistoBound(0,0);
         	std += _HistoConv(0,j)*_HistoBin (0,0)
	 		*(real_value-mean)*(real_value-mean);
      	}
      	
      	if (std < 0.0)  std = 0.0;
      	else std = sqrt(std);
      	
	//if (Verbose==True) cout << "Old sigma = " << std << " and old mean = "
	//	<< mean << endl;
	
	//Update of min, max and bin
	//offset of min and max and dilatation of 1/std
	_HistoBound(0,0)=(_HistoBound(0,0)-mean)/std;
	_HistoBound(0,1)=(_HistoBound(0,1)-mean)/std;
	//dilatation of 1/std
	_HistoBin (0,0)/=std;
	
	//we do not have a density of probability anymore becauseof the new bin	
	//Reduction to have a density of probability
	// sum must be equal to 1
   	sum_check=0.0;
   	for(int i=0;i<_HistoBin(0,1);i++)
      		sum_check += _HistoConv(0,i)*_HistoBin (0,0);
   	// it must be a probability density -> normalisation by sum_check
   	for(int i=0;i<_HistoBin(0,1);i++) 
		_HistoConv(0,i) /= sum_check;
		
		
#if DEBUG_FEW
   if (Verbose==True) show_param ("End Convol Compute", 0, 
                  _HistoBound(0,0), _HistoBound(0,1),
		  _HistoBin(0,0), _HistoBin(0,1)); 		
   if (WriteAllInfo==True) {
      char FileName[256];
      fltarray Courbe ((int)(_HistoBin(0,1)));
      for (int i=0;i<_HistoBin(0,1);i++) Courbe(i)=_HistoConv(0,i);
      sprintf (FileName, "Histo_%d.fits", 0);
      fits_write_fltarr (FileName, Courbe);
      if (Verbose==True) cout << "File  "<< FileName << "  is created  .First Convolution."<< endl; 
   }
#endif   
		  
}



//******************************************************************************
void FewEventPoisson::histo_convolution (Bool WriteAllInfo, Bool Verbose, int param) {
//	compute all convolutions 


   //  param gives the rythm of autoconvolution.
   //  during 2*(2^param) all autoconvolution are computed: 
   // if param=2,  2*(2^param)=8
   // 1	   2       3         4          ...   8
   // S=S1 S*S=S2  S*S*S=S3  S*S*S*S=S4 ...   S*S*S*S*S*S*S*S=S8
   // the next 2^param convolutions are autoconvolutions of the 2^param before:
   // S9=S5*S5  S10=S6*S6  S11=S7*S7  S12=S8*S8
   // same thing for the 4 next...
   // S13=S9*S9  S14=S10*S10  S15=S11*S11 S16=S12*S12  ...
   
  

   // length of first histogram h(0,*)
   int HistoNbPoint= LENGTH_FIRST_HISTO;
   int param_2= 1 << param;// param_2=2^param
   
   int HistoNbPointBegin = HistoNbPoint*(param_2*2-1);
   if(param<0) param=0;
      
   // Reduced Histo used for computation
   // init with h0 (_HistoConv(0,*))
   int MaxHistoNbPoint2=MAXHISTONBPOINT;
   MaxHistoNbPoint2*=2;
   fltarray ReducedHisto (MaxHistoNbPoint2);
   fltarray ReducedHistoBegin (param_2*2-1,HistoNbPointBegin );
   fltarray ReducedHistoTransi (HistoNbPointBegin);
   ReducedHisto.init(0.);
   ReducedHistoBegin.init(0.);
   ReducedHistoTransi.init(0.);
   
   for (int i=0; i<HistoNbPoint; i++)
      ReducedHistoBegin(0,i)= _HistoConv(0,i);
      
   float LocalBin=_HistoBin(0,0);
 
   //Construction of the first histograms
   for (int nbConv=1, IndexFirstHisto=0, IndexSecHisto=0; nbConv<param_2*2-1 ; 
   		nbConv++, IndexFirstHisto=(nbConv)/2, IndexSecHisto=(nbConv-1)/2) {
		
	HistoNbPoint=(int)(_HistoBin(IndexFirstHisto,1)+_HistoBin(IndexSecHisto,1)-1+0.5);	 
	for (long i=0; i < HistoNbPoint ; i++) 
		for (long j=0; j < _HistoBin(IndexSecHisto,1) ; j++)
			if ((i-j < _HistoBin(IndexFirstHisto,1)) && (i-j>=0))
				ReducedHistoBegin(nbConv,i) += ReducedHistoBegin(IndexFirstHisto,i-j)
					*ReducedHistoBegin(IndexSecHisto,j);


      	_HistoBin (nbConv,1) = HistoNbPoint;
	
	//  new_min = min_first_histo + min_sec_histo
	_HistoBound(nbConv,0)=_HistoBound(IndexFirstHisto,0)+_HistoBound(IndexSecHisto,0);
	//  new_max = max_first_histo + max_sec_histo
      	_HistoBound(nbConv,1)=_HistoBound(IndexFirstHisto,1)+_HistoBound(IndexSecHisto,1);
	

	//the bin do not change	 
      	_HistoBin (nbConv,0) = LocalBin;
	if(Verbose==True){
		cout << "nbConv= " << nbConv << endl;
		cout << "Min= " << _HistoBound(nbConv,0) <<"  Max= " << _HistoBound(nbConv,1) << endl;
		cout << "Bin= " << _HistoBin (nbConv,0) <<"  Nb Points= " 
			<< _HistoBin (nbConv,1) << endl<<endl;
	}
	float sum=0.;
      	for (long i=0; i<HistoNbPointBegin; i++) 
         	sum += ReducedHistoBegin(nbConv,i) * LocalBin;
      	for (long i=0; i<HistoNbPointBegin; i++)
         	ReducedHistoBegin(nbConv,i) /= sum;
    }
   
   
   #if DEBUG_FEW      
      if (WriteAllInfo==True) {
	 char FileName[256]="HistoBeginInit.fits";
         fits_write_fltarr (FileName, ReducedHistoBegin);
	 if (Verbose==True) 
	 	cout << "File  "<< FileName << "  is created." <<endl; 
       }		  
   #endif
   
   if(Verbose==True) cout <<endl;
  
  
  
        
   //Reduction of the first histograms
   // no decimation in the first histogrammes
   for(int nbConv=1; nbConv<param_2*2-1 ; nbConv++){
   	
	LocalBin=_HistoBin(nbConv,0);
	HistoNbPoint=(int)(_HistoBin(nbConv,1)+0.5);

	//reduction of bounds
	long Index = 0;
	float IntegMin=0.;
      	while (    ((IntegMin+ReducedHistoBegin(nbConv,Index)*LocalBin) 
		<_Epsilon/2.) 
              		&& (Index < HistoNbPoint)) {
		IntegMin+=ReducedHistoBegin(nbConv,Index)*LocalBin;
		Index++;
	}
      	ReducedHistoBegin(nbConv,Index)+=IntegMin/LocalBin;
	long IndexMin = Index;
      	Index = 0;
	float IntegMax=0.;
      	while (   
	((IntegMax+ReducedHistoBegin(nbConv,(HistoNbPoint-1)-Index)*LocalBin) 
		<_Epsilon/2.) 
              		&& (Index < HistoNbPoint)) {
		IntegMax+=ReducedHistoBegin(nbConv,(HistoNbPoint-1)-Index)*LocalBin;
		Index++;
	}
   	ReducedHistoBegin(nbConv,(HistoNbPoint-1)-Index)+=IntegMax/LocalBin;
	long IndexMax = HistoNbPoint-Index-1;
	ReducedHistoTransi.init(0.);
      	
	for (long i=IndexMin; i<IndexMax+1; i++) 
        	ReducedHistoTransi(i-IndexMin) = ReducedHistoBegin(nbConv,i);
      	
	
      	HistoNbPoint = IndexMax-IndexMin+1;
	// new min, max. the bin do not change
      	_HistoBound(nbConv,0) = _HistoBound(nbConv,0)+IndexMin*LocalBin;
      	_HistoBound(nbConv,1) = _HistoBound(nbConv,1)-Index*LocalBin;     
	_HistoBin(nbConv,1) = HistoNbPoint;

	
	if(Verbose==True){
		cout << "nbConv= " << nbConv << endl;
		cout << "IntegMin= "<<IntegMin <<"  IntegMax="<<IntegMax<<endl;
		cout << "IndexMin= "<< IndexMin<<"  IndexMax= "<< IndexMax << "  Index= "<< Index<<endl; 
		cout << "Min= " << _HistoBound(nbConv,0) <<"  Max= " << _HistoBound(nbConv,1) << endl;
		cout << "Bin= " << _HistoBin (nbConv,0) <<"  Nb Points= " << _HistoBin (nbConv,1) << endl;
	}
	
	// decimation?
	while (HistoNbPoint > MAXHISTONBPOINT) { 
		Index=HistoNbPoint%2-1;
		if (HistoNbPoint+Index > HistoNbPointBegin) 
			cout << "Error of bounds"<< endl;
		
      		// new min, max
      		HistoNbPoint+=Index;
		_HistoBound(nbConv,1) = _HistoBound(nbConv,1)+Index*LocalBin;     
		_HistoBin(nbConv,1) = HistoNbPoint;
		
		
      		// new histogram at next scale will have 2*HistoNbPoint points
      		// if 2*HistoNbPoint > MAXHISTONBPOINT) we have to decimate 
      		// we used a filter with kernel [0.25, 0.5, 0.25] (function shape_signal)
  
        	shape_signal (ReducedHistoTransi);
	 	HistoNbPoint = (HistoNbPoint-1)/2 + 1;
		LocalBin = (_HistoBound(nbConv,1)-_HistoBound(nbConv,0)) / (HistoNbPoint-1);
		_HistoBin(nbConv,1) = HistoNbPoint;	
		_HistoBin(nbConv,0) = LocalBin;
		
		if(Verbose==True){
			cout << "Histo reshaped:"<<endl;
			cout << "New Bin= " << _HistoBin (nbConv,0) 
				<<"New Nb Points= " << _HistoBin (nbConv,1) << endl;
		}
   
   	}
		
		
		
		
	
	//Reduction to have a sigma = 1 and mean = 0
	//Dilatation on x axes
	// the histogramme must have a mean=0, and the bin is divided by std
	// calculation of the mean and the sigma (std)
      	float mean=0.0;
      	float std=0.0;
      	//float LocalBin = _HistoBin(nbConv,0);
	// cout<< _HistoBin(nbConv,1)<<endl;
      	for (int j=0; j<_HistoBin(nbConv,1);j++) {
         	float real_value = ((float) j)*LocalBin + _HistoBound(nbConv,0);
         	mean += ReducedHistoTransi(j)*LocalBin*real_value;
      	}
	// cout << endl;
      	for(int j=0; j<_HistoBin(nbConv,1);j++){
         	float real_value = ((float) j)*LocalBin + _HistoBound(nbConv,0);
         	std += ReducedHistoTransi(j)*LocalBin
	 		*(real_value-mean)*(real_value-mean);
      	}
      	
      	if (std < 0.0)  std = 0.0;
      	else std = sqrt(std);
      	
	//if (Verbose==True) cout << "Old sigma = " << std << " and old mean = "
	//	<< mean << endl;
	//Update of min, max and bin
	//offset of min and max and dilatation of 1/std
	_HistoBound(nbConv,0)=(_HistoBound(nbConv,0)-mean)/std;
	_HistoBound(nbConv,1)=(_HistoBound(nbConv,1)-mean)/std;
	//dilatation of 1/std
	LocalBin/=std;
	_HistoBin(nbConv,0)=LocalBin;
	
	//we do not have a density of probability anymore becauseof the new bin	
	//Reduction to have a density of probability
	// sum must be equal to 1
	float sum=0.;
	for(int i=0;i<HistoNbPoint;i++)
		sum += ReducedHistoTransi(i)*LocalBin;
	for(int i=0;i<HistoNbPoint;i++)
		ReducedHistoBegin(nbConv,i)=_HistoConv(nbConv,i)
			=ReducedHistoTransi(i)/sum;


	
	for(int i=HistoNbPoint;i<MAXHISTONBPOINT;i++)
		ReducedHistoBegin(nbConv,i)=_HistoConv(nbConv,i)=0;
	// for(int i=MAXHISTONBPOINT;i<HistoNbPointBegin;i++)
	// 	ReducedHistoBegin(nbConv,i)=10000;	

	
	
	
	
	#if DEBUG_FEW      
      		if (Verbose==True)
         		show_param ("End Convol Compute", nbConv, 
                     	_HistoBound(nbConv,0), _HistoBound(nbConv,1),
		     	_HistoBin (nbConv,0), _HistoBin (nbConv,1));
      		if (WriteAllInfo==True) {
	 		char FileName[256];
         		fltarray Courbe (HistoNbPoint);
         		for (int i=0;i<HistoNbPoint;i++) 
				Courbe(i)=_HistoConv(nbConv,i);
        		sprintf (FileName, "Histo_%d.fits", nbConv);
        		fits_write_fltarr (FileName, Courbe);
			if (Verbose==True) 
				cout << "File  "<< FileName << "  is created. "
					<< endl << endl; 
		}		  
	#endif 
        
   
   
   
   }
   
   #if DEBUG_FEW      
      if (WriteAllInfo==True) {
	 char FileName[256]="HistoBeginReduced.fits";
         fits_write_fltarr (FileName, ReducedHistoBegin);
	 if (Verbose==True) cout << "File  "<< FileName << "  is created." <<
	 	endl<<endl; 
      }		  
   #endif
   if(Verbose==True) cout <<endl;
   
   
   
   
   
   //building now the other autoconvolution from previous convolutions
   for(int nbConv=param_2*2-1, nbOldConv=param_2-1; 
   		nbConv<_NbAutoConv+1 ; nbConv++, nbOldConv++){ 
   	
	HistoNbPoint=(int)(2*_HistoBin(nbOldConv,1)-1+0.5);	 
	ReducedHisto.init(0.);
	for (long i=0; i < HistoNbPoint ; i++) 
		for (long j=0; j < _HistoBin(nbOldConv,1) ; j++)
			if ((i-j < _HistoBin(nbOldConv,1)) && (i-j>=0))
				ReducedHisto(i) += _HistoConv(nbOldConv,i-j)
					*_HistoConv(nbOldConv,j);
	
      	_HistoBound(nbConv,0) = 2.*_HistoBound(nbOldConv,0);
      	_HistoBound(nbConv,1) = 2.*_HistoBound(nbOldConv,1);
	

	
      	_HistoBin (nbConv,1) = HistoNbPoint;
      	////////LocalBin = _HistoBin(nbOldConv,0); !!!!!!!!!!!1
	LocalBin = (_HistoBound(nbConv,1)-_HistoBound(nbConv,0)) / (float)(HistoNbPoint);		 
      	_HistoBin(nbConv,0) = LocalBin;

	float sum=0.;
      	for (long i=0; i<HistoNbPoint; i++) 
         	sum += ReducedHisto(i) * LocalBin;
      	for (long i=0; i<HistoNbPoint; i++)
         	ReducedHisto(i) /= sum;
 	//reduction of bounds
	
	
	long Index = 0;
	float IntegMin=0.;
      	while (    (IntegMin+ReducedHisto(Index)*LocalBin < _Epsilon/2.) 
              		&& (Index < HistoNbPoint)) {
		IntegMin+=ReducedHisto(Index)*LocalBin;
		Index++;
	}
	
	ReducedHisto(Index)+=IntegMin/LocalBin;
      	long IndexMin = Index;
      	Index = 0;
      	float IntegMax=0.;
	while (    (IntegMax+ReducedHisto((HistoNbPoint-1)-Index)*LocalBin 
		<_Epsilon/2.) 
              		&& (Index < HistoNbPoint)) {
		IntegMax+=ReducedHisto((HistoNbPoint-1)-Index)*LocalBin;
		Index++;
	}

	ReducedHisto((HistoNbPoint-1)-Index)+=IntegMax/LocalBin;
   	long IndexMax = HistoNbPoint-Index-1;
	
      	
	for (long i=IndexMin; i<IndexMax+1; i++) 
        	ReducedHisto(i-IndexMin) = ReducedHisto(i);
      	for (long i=IndexMax-IndexMin+1; i<MaxHistoNbPoint2; i++) 
		ReducedHisto(i)=0.;

      	HistoNbPoint = IndexMax-IndexMin+1;
	// new min, max
      	_HistoBound(nbConv,0) = _HistoBound(nbConv,0)+IndexMin*LocalBin;
      	_HistoBound(nbConv,1) = _HistoBound(nbConv,1)-Index*LocalBin;     
	_HistoBin(nbConv,1) = HistoNbPoint;
	//!!!!!!!!!!!!!!!!!!!!1
	LocalBin = (_HistoBound(nbConv,1)-_HistoBound(nbConv,0)) / (HistoNbPoint-1);
	_HistoBin(nbConv,0) = LocalBin;
	
	if(Verbose==True){
		cout << "nbConv= " << nbConv << endl;
		cout << "IntegMin= "<<IntegMin <<"  IntegMax="<<IntegMax<<endl;
		cout << "IndexMin= "<< IndexMin<<"  IndexMax= "<< IndexMax << "  Index= "<< Index<<endl; 
		cout << "Min= " << _HistoBound(nbConv,0) <<"  Max= " << _HistoBound(nbConv,1) << endl;
		cout << "Bin= " << _HistoBin (nbConv,0) <<"  Nb Points= " << _HistoBin (nbConv,1) << endl;
	}
	
	
	// decimation?
	while (HistoNbPoint > MAXHISTONBPOINT) { 
		
		Index=HistoNbPoint%2-1;
		if (HistoNbPoint+Index > MaxHistoNbPoint2) 
			cout << "Error of bounds"<< endl;
		
      		// new min, max, bin
      		HistoNbPoint+=Index;
		_HistoBound(nbConv,1) = _HistoBound(nbConv,1)+Index*LocalBin;     
		_HistoBin(nbConv,1) = HistoNbPoint;
		
		
		
      		// new histogram at next scale will have 2*HistoNbPoint points
      		// if 2*HistoNbPoint > MAXHISTONBPOINT) we have to decimate 
      		// we used a filter with kernel [0.25, 0.5, 0.25] (function shape_signal)

        	shape_signal (ReducedHisto);
	 	HistoNbPoint = (HistoNbPoint-1)/2 + 1;
         	LocalBin = (_HistoBound(nbConv,1)-_HistoBound(nbConv,0)) / (HistoNbPoint-1);	
		_HistoBin(nbConv,0) = LocalBin;
		_HistoBin(nbConv,1) = HistoNbPoint;
		
		if(Verbose==True){
			cout << "Histo reshaped"<<endl;
			cout << "Min= " << _HistoBound(nbConv,0) <<"  Max= " << _HistoBound(nbConv,1) << endl;
			cout << "Bin= " << _HistoBin (nbConv,0) <<"  Nb Points= " << _HistoBin (nbConv,1) << endl;
		}
   
   	}
	

	
	//Reduction to have a sigma = 1 and mean = 0
	//Dilatation on x axes
	// the histogramme must have a mean=0, and the bin is divided by std
	// calculation of the mean and the sigma (std)
      	float mean=0.0;
      	float std=0.0;
      	//float LocalBin = _HistoBin(nbConv,0);
      	for (int j=0; j<_HistoBin(nbConv,1);j++) {
         	float real_value = ((float) j)*LocalBin + _HistoBound(nbConv,0);
         	mean += ReducedHisto(j)*LocalBin*real_value;
      	}

      	for(int j=0; j<_HistoBin(nbConv,1);j++){
         	float real_value = ((float) j)*LocalBin + _HistoBound(nbConv,0);
         	std += ReducedHisto(j)*LocalBin
	 		*(real_value-mean)*(real_value-mean);
      	}
      	if (std < 0.0)  std = 0.0;
      	else std = sqrt(std);
 	//Update of min, max and bin
	//offset of min and max and dilatation of 1/std
	if (std > 0)
	{
	   _HistoBound(nbConv,0)=(_HistoBound(nbConv,0)-mean)/std;
	   _HistoBound(nbConv,1)=(_HistoBound(nbConv,1)-mean)/std;
	  //dilatation of 1/std
	   LocalBin/=std;
	   _HistoBin(nbConv,0)=LocalBin;
	}
	else 
	{
	   _HistoBound(nbConv,0)= _HistoBound(nbConv,1)=0;
	}
 
	
	//we do not have a density of probability anymore becauseof the new bin	
	//Reduction to have a density of probability
	// sum must be equal to 1
	sum=0.;
	for(int i=0;i<HistoNbPoint;i++)
		sum += ReducedHisto(i)*LocalBin;
	for(int i=0;i<HistoNbPoint;i++)
		_HistoConv(nbConv,i)=ReducedHisto(i)/sum;

	//just to see the bugs
   	// for(long i=HistoNbPoint;i<MAXHISTONBPOINT;i++)
	// 	_HistoConv(nbConv,i)=10000;

		
		
		
	#if DEBUG_FEW      
      		if (Verbose==True)
         		show_param ("End Convol Compute", nbConv, 
                     	_HistoBound(nbConv,0), _HistoBound(nbConv,1),
		     	_HistoBin (nbConv,0), _HistoBin (nbConv,1));
      		if (WriteAllInfo==True) {
	 		char FileName[256];
         		fltarray Courbe (HistoNbPoint);
         		for (int i=0;i<HistoNbPoint;i++) 
				Courbe(i)=_HistoConv(nbConv,i);
        		sprintf (FileName, "Histo_%d.fits", nbConv);
        		fits_write_fltarr (FileName, Courbe);
			if (Verbose==True) {
				cout << "File  "<< FileName << "  is created. ";
				cout << "Convolution nbr "<< nbConv<< " ." << endl; 
			}
		}		  
	#endif 
   }

}




//******************************************************************************
void FewEventPoisson::histo_distribution (Bool WriteAllInfo, Bool Verbose) {
//compute the distribution's laws from the histogramms 

   for (int nbConv=0; nbConv<_NbAutoConv+1; nbConv++) {
   
      float LocalBin     = _HistoBin(nbConv,0);
      float LocalNbPoint = _HistoBin(nbConv,1);
      
      //Probability Reduced Values
      //
      // init first value
      _HistoDistrib(nbConv,0)   = _HistoConv(nbConv,0)*LocalBin;
      
      // integ other val (Dist(i) = Dist(i-1) + val(i)*bin)
      for (int nbp=1; nbp<LocalNbPoint; nbp++) 
	 	_HistoDistrib(nbConv,nbp) = _HistoDistrib(nbConv,nbp-1) 
	                    + _HistoConv(nbConv,nbp)*LocalBin;
   }
#if DEBUG_FEW   
   if (Verbose==True) 
      for (int nbConv=0; nbConv<_NbAutoConv+1; nbConv++) 
      {
         int HB1 = (int) _HistoBin(nbConv,1);
         cout << "Scale:" << nbConv << endl;
	 cout << "Min:" << _HistoBound(nbConv,0)
	      << ", Max:" << _HistoBound(nbConv,1)
	      << ", bin:" << _HistoBin(nbConv,0) << endl;
	 cout << "d(0):" << _HistoDistrib(nbConv,0)
	      << ", d(" << (HB1-1)/2 << "):" 
	                << _HistoDistrib(nbConv, (HB1-1)/2)
	      << ", d(" << _HistoBin(nbConv,1)-1 << "):" 
	                << _HistoDistrib(nbConv, HB1-1) << endl;		
      }
   if (WriteAllInfo==True){
      fits_write_fltarr("NewDistrib.fits", _HistoDistrib);   
      if(Verbose==True) cout << "File  NewDistrib.fits  is created." << endl; 
   }
#endif
}



//******************************************************************************
void FewEventPoisson::shape_signal (fltarray &Signal1d) {
//	the signal is filtered by the kernel { 1/4 ; 1/2 ; 1/4 }
//	and it is down sampled by 2


   // init var
   long LocalSignalLength = Signal1d.nx();
   fltarray Filter_Signal1d(LocalSignalLength);
  
   // filter signal
   for (int i=1;i<=(LocalSignalLength-1)/2-1;i++)
      Filter_Signal1d(2*i) = 0.5  * Signal1d(2*i) + 
                             0.25 * (Signal1d(2*i-1) + Signal1d(2*i+1));
   
   Filter_Signal1d(0)	= 0.5*(Signal1d(0) + Signal1d(1));		    
   Filter_Signal1d(LocalSignalLength-1) = 0.5*(  Signal1d(LocalSignalLength-2) 
                                               + Signal1d(LocalSignalLength-1));
						
   // Sampling and put zeros
   Signal1d.init(0.);
   for (int i=0;i<=(LocalSignalLength-1)/2;i++) 
      Signal1d(i)= Filter_Signal1d(2*i);
}


//******************************************************************************
void FewEventPoisson::show_param (char* text, float nbconv, float min, 
                                    float max, float bin, float nbechant) {
			    
   cout << "Histogram number <" << nbconv << "> // " << text << endl;
   cout	<< "  -min:" << min << ", -max:" << max << "' bin:" << bin
        << ", (NBP:" << nbechant << ")" << endl;
}

//******************************************************************************
void FewEventPoisson::histo_threshold (Bool WriteAllInfo, Bool Verbose) {
//	compute the threshold max and the threshold min for every 
//	distribution laws. 


   for (int nbConv=0; nbConv<_NbAutoConv+1; nbConv++) {
   
      // serach indice of max threshold (RepFunc >= 1-eps) 
      // and min threshold (RepFunc <= eps)
      //
      long IndexMax=0;
      long IndexMin=(long)(_HistoBin(nbConv,1)+0.5)-1;
      
      
      // _HistoDistrib is a distribution of probability
      // _HistoDistrib(min)=~0     _HistoDistrib(max)=~1
      // RepFuc = 1-eps in [IndexMax-1:IndexMax]
      while (    (_HistoDistrib(nbConv,IndexMax) < 1-_Epsilon)
              && (IndexMax<_HistoBin(nbConv,1)-1)) {IndexMax++;}
      // RepFuc = eps in [IndexMin:IndexMin+1]      
      while (    (_HistoDistrib(nbConv,IndexMin) > _Epsilon)
              && (IndexMin>0)) {IndexMin--;}
      //cout << "Ind min:" << IndexMin <<", IndMax : " << IndexMax << endl;	      
      // test	      
      if ((IndexMax==0) || (IndexMin==_HistoBin(nbConv,1)-1)) {
	  cout << "ERROR : bad sampling of the distribution function " << endl;
	  exit (-1);
      }

      // inter between IndexMax and IndexMax-1
      // X1-eps = Xmax + (1-eps-Ymax)/a,  a=(Ymax-1 - Ymax)/(Xmax-1 - Xmax)
      // 
      float LocalBin=_HistoBin(nbConv,0);
      float Min=_HistoBound(nbConv,0);
      
      float a =   _HistoConv(nbConv,IndexMax);//*Sigma;
	      
      //float b =        _HistoDistrib(nbConv+2,IndexMax-1) 
      //           - a * _HistoDistrib(nbConv+1,IndexMax-1);
         
      //_Threshold (nbConv, 1) =   _HistoDistrib(nbConv+1,IndexMax) 
      //                          + (1-_Epsilon-_HistoDistrib(nbConv+2,IndexMax))/a;
      
      _Threshold (nbConv, 1) =   (IndexMax*LocalBin+Min)///Sigma 
                 +(1-_Epsilon-_HistoDistrib(nbConv,IndexMax))/a;
       
      
      //_Threshold (1, nbConv) = (1-_Epsilon-b)/a;
      
      
      // inter between IndexMin and IndexMin+1
      // Xeps = Xmin + (eps-Ymin)/a,  a=(Ymin+1 - Ymin)/(Xmin+1 - Xmin)
      a =  _HistoConv(nbConv,IndexMin+1);//*Sigma;
      
      _Threshold (nbConv, 0) =   (IndexMin*LocalBin+Min)///Sigma 
	         + (_Epsilon-_HistoDistrib(nbConv,IndexMin))/a;  
      

      
      if (Verbose==True) {
         cout << "Nbr of histogram autoconv.:" << nbConv << ", SeuilMin := " << _Threshold (nbConv, 0)
                                    << ", SeuilMax := " << _Threshold (nbConv, 1)
      			            << endl; 
      }     
   }
   if (WriteAllInfo==True) {
      io_write_ima_float (Name_Threshold, _Threshold);
      if(Verbose==True) cout << "File  "<< Name_Threshold << "  is created (Threshold) ." << endl; 
     
   }
}

//******************************************************************************

void FewEventPoisson::find_threshold (Bool WriteAllInfo, Bool WriteHisto, 
       Bool Verbose,int param) {
//	test if the thresholds are computed or read.
   if(Verbose==True) cout << "Initalisation of all histograms" << endl;
   if (_InitOk == False) {
   	if(Verbose==True) cout <<"Computing Wavelet Distribution"<<endl;
   	compute_distribution(WriteAllInfo,WriteHisto,Verbose, param);
   	
   }
   else{
   	if(Verbose==True) cout <<"Reading Wavelet Distribution"<<endl;
   	fits_read_fltarr(Name_HistoConv, _HistoConv);
      	fits_read_fltarr(Name_HistoBound, _HistoBound);      
      	fits_read_fltarr(Name_HistoBin, _HistoBin);
	fits_read_fltarr(Name_HistoDistrib, _HistoDistrib);
   } 
      
   histo_threshold ( WriteAllInfo);
   if(Verbose==True) cout << "Histogram of the wavelet function" << endl;
}


//******************************************************************************

void FewEventPoisson::get_sigma_band(fltarray & Mean_Nb_Ev, fltarray & tab_sigma,
   		float N_Sigma,int Nb_Scale)
{
//give a approximation of the sigma of noise in each band. 
//_tab_sigma(*,0) is the sigma of 1
//_tab_sigma(*,1) is the sigma of N_Sigma
	// float epsilon1=(1.-erf(1/sqrt(2.)));//gaussian equivalent of 1
	// float epsilon2=(1.-erf(N_Sigma/sqrt(2.)));//gaussian equivalent of N_Sigma
   float sigma_wavelet=SIGMA_BSPLINE_WAVELET;
   float Alpha;
   int param=(int)(_Param(0)+0.5);
   int param_2=1 << param;
   int s;
   for(s=0, Alpha=1.;s<Nb_Scale-1;s++, Alpha*=8.)
   {
      // cout << "MMMMean_Nb_Ev " << s+1 << " = " << Mean_Nb_Ev(s) << endl;
      int NConvReal=(int)(Mean_Nb_Ev(s)+0.5)-1;
      int Power=-1;
      int NConv=1;
      int NConvIndex;
      if (NConvReal<0)
      { 
          NConvIndex=0; NConv=0;NConvReal=0;
      }
      else
      { 	
         while (NConv<=NConvReal) 
	 {
            Power++; 
            NConv *= 2;
      	 }
         NConv/=2;//now NConv=2^Power
      		
         if (param>Power) NConvIndex=NConv;	
         else
	 {
             int B_Step=(Power-param+1)*param_2;	//big step 
             int L_Step=(int)(( NConvReal - (1<<Power) ) / (1 << (Power - param)));
             NConvIndex=B_Step+L_Step;
      	 }
      }
      if (NConvIndex >= _NbAutoConv)  NConvIndex = _NbAutoConv;
	   
	   
      for(int k=0;k<2; k++)
      {
 		float Epsilon=1.-erf((k==0 ? 1:N_Sigma) / sqrt(2.));
		long IndexMax=0;
		long IndexMin=(long)(_HistoBin(NConvIndex,1)+0.5)-1;
      		while ((_HistoDistrib(NConvIndex,IndexMax) < 1-Epsilon)
              		&& (IndexMax<_HistoBin(NConvIndex,1)-1))      {IndexMax++;}
    
      		while ((_HistoDistrib(NConvIndex,IndexMin) > Epsilon)
              		&& (IndexMin>0))  {IndexMin--;}

       		float LocalBin=_HistoBin(NConvIndex,0);
      		float Min=_HistoBound(NConvIndex,0);
      		// float a =_HistoConv(NConvIndex,IndexMax);
 		float Threshold_Min=0,Threshold_Max=0;
		if (_HistoConv(NConvIndex,IndexMax) > 0)
		     Threshold_Max = (IndexMax*LocalBin+Min) 
                 	+(1-_Epsilon-_HistoDistrib(NConvIndex,IndexMax))
		 	/_HistoConv(NConvIndex,IndexMax);
       		if (_HistoConv(NConvIndex,IndexMin+1) > 0)
		     Threshold_Min = (IndexMin*LocalBin+Min)  
      			+ (_Epsilon-_HistoDistrib(NConvIndex,IndexMin))
			/_HistoConv(NConvIndex,IndexMin+1);
 			
		// normalisation of coeff:
		//cout <<"avant   TH MIN="<< Threshold_Min <<"	TH MAX="<< Threshold_Max<< " Alpha = " << Alpha << endl;
		//cout <<"NConvReal = " << NConvReal << " sigma_wavelet = " << sigma_wavelet << endl;
		Threshold_Min *= sqrt((float)(NConvReal+1)) * sigma_wavelet / Alpha;	
	 	Threshold_Max *= sqrt((float)(NConvReal+1)) * sigma_wavelet / Alpha;
		//cout <<NConvReal+1<<"  "<<sqrt((float)(NConvReal+1))<<"   "<<Alpha<<endl;
		
      		//cout <<"apresTH MIN="<< Threshold_Min <<"	TH MAX="<< Threshold_Max<<endl;
		//cout <<endl;
		tab_sigma(s,k)=(Threshold_Max-Threshold_Min)/2.;
		// cout << "s= " << s+1 << "tab_sigma = " << tab_sigma(s,k) << endl;
	}
    }
}

//******************************************************************************

void FewEventPoisson::get_threshold(int NbrEvent, int CurrentScale, float & SeuilMin, float & SeuilMax)
{
    int param=(int)(_Param(0)+0.5);
    int param_2=1 << param;
    int ABAQUE_N_SCALE = _NbAutoConv;
    float sigma_wavelet=SIGMA_BSPLINE_WAVELET;
    int NConvReal = NbrEvent-1;
    int Power=-1;
    int NConv=1;//NConv= index in _HistoConv to find histogramme of NConvReal
    Bool equality=False; //	equlity=true => no interpolation needed
    int NConvIndex;
    float Alpha=1.;
    //if (CurrentScale != 0) Alpha = pow((double)8., (double)(CurrentScale-1.));
    if (CurrentScale != 0) Alpha = pow((double)8., (double)(CurrentScale));   
    
   
    if (NConvReal<0)
    { 
	NConvIndex=0;
      	NConv=0;
	equality=True;
    }
    else
    { 	
        while (NConv<=NConvReal) 
	{
      	   Power++; 
	   NConv *= 2;
      	}
	NConv/=2; //now NConv=2^Power
      	
	if (param>Power) NConvIndex=NConv;	
	else
	{
      	   int B_Step=(Power-param+1)*param_2;	//big step 
	   int L_Step=(int)(( NConvReal - (1<<Power) ) / (1 << (Power - param)));
	   if ((( NConvReal - (1<<Power) ) % (1 << (Power - param)))==0)  
			equality=True;
	   else  equality=False;
	   NConvIndex=B_Step+L_Step;
      	}
    }
      
    // interpolation of thresholds using abaque
    SeuilMin = SeuilMax = 0.;
    if (NConvIndex >= ABAQUE_N_SCALE) 
    {
	 // nb scale in abaque to small, take the max value in abaque     
         NConvIndex = ABAQUE_N_SCALE;
	 equality=True;
	 //cout << "Abaque too small" << endl;
	 SeuilMin = _Threshold(NConvIndex,0);
	 SeuilMax = _Threshold(NConvIndex,1);     
     
    } 
    else if (equality) 
    {
          // NEventReal is a power of two = NEvent, read value in _Threshold tab
         SeuilMin = _Threshold(NConvIndex,0);
         SeuilMax = _Threshold(NConvIndex,1);
    } 
    else 
    {
          // interpolate in Abaque tab between level power-1 and power
        if (NConv==0) 
        {
	    SeuilMin = _Threshold(0,0);
	    SeuilMax = _Threshold(0,1);
	 } 
	 else 
	 {
	    // linear interp 
	    
	    SeuilMin = -2*(_Threshold(NConvIndex+1,0) - _Threshold(NConvIndex,0))
                          *(NConvReal - NConv) / NConvReal 
			  + _Threshold(NConvIndex+1,0);
            SeuilMax = -2*(_Threshold(NConvIndex+1,1) - _Threshold(NConvIndex,1))
                          *(NConvReal - NConv) / NConvReal 
			  + _Threshold(NConvIndex+1,1);
         }
      }	
      
      // sigma of wavelet IN (function psy of bspline_histo_1D)
      //
      //float Sigma=0.121070;   
      int NEventReal=NConvReal+1;
      if (Power > 0) 
      {
         SeuilMin *= sqrt((float)(NEventReal)) * sigma_wavelet / Alpha;	
	 SeuilMax *= sqrt((float)(NEventReal)) * sigma_wavelet / Alpha;
      } 
      else 
      {
	 SeuilMin *= sigma_wavelet  / Alpha;
	 SeuilMax *= sigma_wavelet  / Alpha;
      }
}


//******************************************************************************

void FewEventPoisson::event_set_support(
//	create the support with the threshold values and the wavelet coeff


		       fltarray *&	 Data_In,
		       int               CurrentScale, 
		       int		 FirstDectectScale,
		       int		 NEGFirstDectectScale,
		       int 		 MinEventNumber,
		       Bool		 OnlyPositivDetect,
		       type_border       Border,
		       //number of events
		       intarray *& 	 Nb_Event,
		       //cube of support of the current scale
		       intarray *& 	 Support,
		       Bool              WriteAllInfo) 
{
   int FDS = FirstDectectScale;
   int NFDS = (NEGFirstDectectScale < 0) ? FirstDectectScale: NEGFirstDectectScale;   
   int Nx = Data_In[CurrentScale].nx();
   int Ny = Data_In[CurrentScale].ny();
   int Nz = Data_In[CurrentScale].nz();

   for (int i=0;i<Nx;i++)
   for (int j=0;j<Ny;j++)
   for (int k=0;k<Nz;k++)
   {
      // reading the number of events
      float SeuilMin, SeuilMax;
      // number of convolution by the first histogramme needed
      // compute "power" : 2^(power) <= NConvReal < 2^(power+1) 
      get_threshold(Nb_Event[CurrentScale](i,j,k), CurrentScale, SeuilMin, SeuilMax);
      
      // detection of signal : comparison with thresholds
      //

      Support[CurrentScale](i,j,k) = VAL_SupNull;
      if ((Data_In[CurrentScale](i,j,k)<=SeuilMin) || (Data_In[CurrentScale](i,j,k)>=SeuilMax)) 
      {
	     
         Support[CurrentScale](i,j,k) = VAL_SupOK;
	 	 
	 // 
	 if (Nb_Event[CurrentScale](i,j,k) < MinEventNumber) 
                     Support[CurrentScale](i,j,k) = VAL_SupMinEv;
	 
	 //
	 if ((OnlyPositivDetect) && (Data_In[CurrentScale](i,j,k)<0))
	    Support[CurrentScale](i,j,k) = VAL_SupNull;
	 
	 //    
	 if ((FDS > CurrentScale) && (Data_In[CurrentScale](i,j,k) >= 0))
	                    Support[CurrentScale](i,j,k) = VAL_SupFirstScale;

         if ((NFDS > CurrentScale) && (Data_In[CurrentScale](i,j,k) < 0))
	                    Support[CurrentScale](i,j,k) = VAL_SupFirstScale;
      }
   }
}

//****************************************************************************** 
int FewEventPoisson::get_pix(int i,int j,int k, intarray & DataIn, type_border Border){
//	give the value of a pixel depending to the border
	int Nx = DataIn.nx();
   	int Ny = DataIn.ny();
   	int Nz = DataIn.nz();
	int ind_i=0, ind_j=0, ind_k=0;
	switch(Border){	
		case I_CONT:
			ind_i=(i < 0 ? 0 : i);
			if(ind_i>=Nx) ind_i=Nx-1;
			ind_j=(j < 0 ? 0 : j);
			if(ind_j>=Ny) ind_j=Ny-1;
			ind_k=(k < 0 ? 0 : k);
			if(ind_k>=Nz) ind_k=Nz-1;
			break;
		case I_MIRROR:
			ind_i=(i < 0 ? -i : i);
			if(ind_i>=Nx) ind_i=2*(Nx-1)-ind_i;
			ind_j=(j < 0 ? -j : j);
			if(ind_j>=Ny) ind_j=2*(Ny-1)-ind_j;
			ind_k=(k < 0 ? -k : k);
			if(ind_k>=Nz) ind_k=2*(Nz-1)-ind_k;
			break;
		case I_ZERO:
			ind_i=(i < 0 ? 0 : i);
			if(ind_i>=Nx) ind_i=0;
			ind_j=(j < 0 ? 0 : j);
			if(ind_j>=Ny) ind_j=0;
			ind_k=(k < 0 ? 0 : k);
			if(ind_k>=Nz) ind_k=0;
			break;
		case I_PERIOD:
			ind_i=(i < 0 ? Nx+i : i);
			if(ind_i>=Nx) ind_i=i-Nx;
			ind_j=(j < 0 ? Ny+j : j);
			if(ind_j>=Ny) ind_j=j-Ny;
			ind_k=(k < 0 ? Nz+k : k);
			if(ind_k>=Nz) ind_k=k-Nz;
			break;
	}
	return( DataIn(ind_i, ind_j, ind_k));	
}
//******************************************************************************
int FewEventPoisson::get_ev(int i,int j,int k,intarray & nb_event, type_border  Border){
	//get the number of events is a box already computed depending to the border
	
	int Nx = nb_event.nx();
   	int Ny = nb_event.ny();
   	int Nz = nb_event.nz();
	int ind_i, ind_j, ind_k;
	int Return=0;
	switch(Border){	
		case I_CONT:
			cout << "I_CONT border is not implemented in case of no approximation "<<endl;
			exit(0);
			break;
		case I_MIRROR:
			ind_i=(i < 0 ? -i+1 : i);
			if(ind_i>=Nx) ind_i=2*Nx-1-ind_i;
			ind_j=(j < 0 ? -j+1 : j);
			if(ind_j>=Ny) ind_j=2*Ny-1-ind_j;
			ind_k=(k < 0 ? -k+1 : k);
			if(ind_k>=Nz) ind_k=2*Nz-1-ind_k;
			Return=nb_event(ind_i, ind_j, ind_k);
			break;
		case I_ZERO:
			if(i<0 || i>=Nx || j<0 || j>=Ny || k<0 || k>=Nz)
				Return=0;
			else 
				Return=nb_event(i, j, k);
			break;
		case I_PERIOD:
			ind_i=(i < 0 ? Nx+i : i);
			if(ind_i>=Nx) ind_i=i-Nx;
			ind_j=(j < 0 ? Ny+j : j);
			if(ind_j>=Ny) ind_j=j-Ny;
			ind_k=(k < 0 ? Nz+k : k);
			if(ind_k>=Nz) ind_k=k-Nz;
			Return=nb_event(ind_i, ind_j, ind_k);
			break;
	}
	return( Return);	
}





//******************************************************************************
void FewEventPoisson::event_scale_3D_Approx (	
					intarray & last_nb_event, 
					intarray & nb_event,
					int Current_Scale, 
					type_border Border,
					Bool WriteAllInfo) {
//compute the number of event in boxes as event_scale_3D but with the approxination
// of box = 2^n * 2^n * 2^n  and not box=(2^n+1) * (2^n+1) * (2^n+1) 
// the convention is nb_event[n](i,j,k)=sum(	cur_i from i-2^(n+2) to i+2^(n+2)-1, 
//						cur_j from j-2^(n+2) to j+2^(n+2)-1,
//						cur_k from k-2^(n+2) to k+2^(n+2)-1,
//							last_nb_event(cur_i, cur_j, cur_k)
	int Nx = last_nb_event.nx();
   	int Ny = last_nb_event.ny();
   	int Nz = last_nb_event.nz();
	if (Current_Scale==-2)
	{
 		int * tab;
		tab= new int[2];
		//*****************computation of the first point: (0, 0, 0)****************
		//(0,0,0)
		nb_event(0,0,0)=0;
		for (int cur_i=-1;cur_i<=0;cur_i++)
		for (int cur_j=-1;cur_j<=0;cur_j++)
		for (int cur_k=-1;cur_k<=0;cur_k++)
			nb_event(0,0,0)+=get_pix(cur_i, cur_j, cur_k, last_nb_event, Border);
			
		//*****************computation of the first ligne:  (1:Nx-1 ,0 ,0)**********
		//(1:2,0,0)
		for (int i=1;i<=2;i++)
		{
			tab[i-1]=0;
			nb_event(i,0,0)=nb_event(i-1,0,0);
			for (int cur_j=-1;cur_j<=0;cur_j++)
			for (int cur_k=-1;cur_k<=0;cur_k++){
				tab[i-1]+=get_pix(i, cur_j, cur_k, last_nb_event, Border);
				nb_event(i,0,0)-=get_pix(i-2, cur_j, cur_k, last_nb_event, Border);
		
			}
			nb_event(i,0,0)+=tab[i-1];
		}
		//(3:Nx-1,0,0)
		for (int i=3;i<Nx;i++)
		{
			nb_event(i,0,0)=nb_event(i-1,0,0) - tab[(i-1) % 2];
			tab[(i-1) % 2]=0;
			for (int cur_j=-1;cur_j<=0;cur_j++)
			for (int cur_k=-1;cur_k<=0;cur_k++)
				tab[(i-1) % 2]+=get_pix(i, cur_j, cur_k, last_nb_event, Border);
			nb_event(i,0,0)+=tab[(i-1) % 2];	
		}
	
		
	
		//********************computation of the first plan:  (0:Nx-1 ,1:Ny-1 ,0)*********
		for (int i=0;i<Nx;i++)
		{
			//(0:Nx-1, 1:2, 0)
			for (int j=1;j<=2;j++)
			{
				tab[j-1]=0;
				nb_event(i,j,0)=nb_event(i,j-1,0);
				for (int cur_i=-1;cur_i<=0;cur_i++)
				for (int cur_k=-1;cur_k<=0;cur_k++)
				{
					tab[j-1]+=get_pix(i+cur_i, j, cur_k, last_nb_event, Border);
					nb_event(i,j,0)-=get_pix(i+cur_i, j-2, cur_k, last_nb_event, 
						Border);
				}
				nb_event(i,j,0) += tab[j-1];
			}
			//(0:Nx-1 ,3:Ny, 0)
			for (int j=3;j<Ny;j++){
				nb_event(i,j,0)=nb_event(i,j-1,0) - tab[(j-1) % 2];
				tab[(j-1) % 2]=0;
				for (int cur_i=-1;cur_i<=0;cur_i++)
				for (int cur_k=-1;cur_k<=0;cur_k++)
					tab[(j-1) % 2] += get_pix(i+cur_i, j, cur_k, last_nb_event, Border);   
				nb_event(i,j,0) += tab[(j-1) % 2];
			}
		}	
		
		//********************computation of the cube:  (0:Nx-1 ,0:Ny-1 ,1:Nz-1)********
		for (int i=0;i<Nx;i++){
		for (int j=0;j<Ny;j++){
			//(0:Nx-1 ,0:Ny-1 ,1:2)
			for (int k=1;k<=2;k++){
				nb_event(i,j,k)=nb_event(i,j,k-1);
				tab[k-1]=0;
				for (int cur_i=-1;cur_i<=0;cur_i++)
				for (int cur_j=-1;cur_j<=0;cur_j++){
					tab[k-1] += get_pix(i+cur_i, j+cur_j, k, last_nb_event, Border);
					nb_event(i,j,k) -= 
						get_pix(i+cur_i, j+cur_j, k-2, last_nb_event, Border);
				}
				nb_event(i,j,k) += tab[k-1];
			}
			//(0:Nx-1 ,0:Ny-1 ,3:Nz)
			for (int k=3;k<Nz;k++){
				nb_event(i,j,k)=nb_event(i,j,k-1) - tab[(k-1) % 2];
				tab[(k-1) % 2]=0;
				for (int cur_i=-1;cur_i<=0;cur_i++)
				for (int cur_j=-1;cur_j<=0;cur_j++)
		   			tab[(k-1) % 2] += 
						get_pix(i+cur_i, j+cur_j, k, last_nb_event, Border);   
				nb_event(i,j,k) += tab[(k-1) % 2];
			}
		
		}
		if(WriteAllInfo==True) {
			cout << flush;
			float pc=i*100/(float)(Nx);
			pc=((int)(pc*10))/10.;
			cout << '\r' << "Scale number -2, " << pc << "% done.       ";
		}
		}
		if(WriteAllInfo==True) cout << '\r';
	delete[] tab;		
	} 
	else 
	{
		int ind=(int) (pow((double)2., (double)(Current_Scale+1)) + 0.5);
		for(int i=0;i<Nx;i++)
		{
		for(int j=0;j<Ny;j++)
		for(int k=0;k<Nz;k++)  nb_event(i,j,k)=
			 get_ev(i-ind, j-ind, k-ind, last_nb_event, Border)
			+get_ev(i-ind, j-ind, k+ind, last_nb_event, Border)
			+get_ev(i-ind, j+ind, k-ind, last_nb_event, Border)
			+get_ev(i-ind, j+ind, k+ind, last_nb_event, Border)
			+get_ev(i+ind, j-ind, k-ind, last_nb_event, Border)
			+get_ev(i+ind, j-ind, k+ind, last_nb_event, Border)
			+get_ev(i+ind, j+ind, k-ind, last_nb_event, Border)
			+get_ev(i+ind, j+ind, k+ind, last_nb_event, Border);
			
		  if(WriteAllInfo==True) 
		  {		
			cout << flush;
			float pc=i*100/(float)(Nx);
			pc=((int)(pc*10))/10.;
			cout << '\r' << "Scale number "<<Current_Scale << ", " << pc << "% done.       ";
		  }
		}
		if(WriteAllInfo==True) cout << '\r';
	}		
	// int total_ind=(int) (pow((double)2., (double)(Current_Scale+3)) + 0.5);
	// if(WriteAllInfo==True) 
	// cout << "Table of number of events in cube "
	//    << total_ind << '*' << total_ind << '*' << total_ind << " is created" <<endl;
	// cout << "Event_scale_3D : " << nb_event.min() << endl;

}



//******************************************************************************

void FewEventPoisson::event_scale_3D (	intarray & data_event, 
					intarray *& nb_event,
					int Current_Scale, 
					type_border Border,
					Bool WriteAllInfo) {
	//number of events in each cube of total_ind*total_ind*total_ind points
	// in case of no approximation 
	
	int ind= (int) (pow((double)2., (double)(Current_Scale+2)) + 0.5);
	int total_ind=ind*2+1;
	
	int*tab;
	tab= new int[total_ind];
	int Nx = data_event.nx();
   	int Ny = data_event.ny();
   	int Nz = data_event.nz();
	

	
	//********************computation of the first point: (0, 0, 0)***********
	nb_event[Current_Scale](0,0,0)=0;
	for (int cur_i=-ind;cur_i<=ind;cur_i++)
	for (int cur_j=-ind;cur_j<=ind;cur_j++)
	for (int cur_k=-ind;cur_k<=ind;cur_k++)
		nb_event[Current_Scale](0,0,0)+=get_pix(cur_i, cur_j, cur_k, 
			data_event, Border);
	// cout << "Event_scale_3D 1: " << nb_event[Current_Scale].min() << endl;
	
	//********************computation of the first ligne:  (0:Nx-1 ,0 ,0)*****
	
	//(1:total_ind,0,0)
	for (int i=1;i<=total_ind;i++){
		tab[i-1]=0;
		nb_event[Current_Scale](i,0,0)=nb_event[Current_Scale](i-1,0,0);
		for (int cur_j=-ind;cur_j<=ind;cur_j++)
		for (int cur_k=-ind;cur_k<=ind;cur_k++){
			tab[i-1]+=get_pix(i+ind, cur_j, cur_k, data_event, Border);
			nb_event[Current_Scale](i,0,0)-=get_pix(i-ind-1, cur_j, 
				cur_k, data_event, Border);
		}
		nb_event[Current_Scale](i,0,0)+=tab[i-1];
	}
        // cout << "Event_scale_3D 2: " << nb_event[Current_Scale].min() << endl;

	//(total_ind+1:Nx-1,0,0)
	for (int i=total_ind+1;i<Nx;i++){
		nb_event[Current_Scale](i,0,0)=nb_event[Current_Scale](i-1,0,0) 
			- tab[(i-1) % total_ind];
		tab[(i-1) % total_ind]=0;
		for (int cur_j=-ind;cur_j<=ind;cur_j++)
		for (int cur_k=-ind;cur_k<=ind;cur_k++)
			tab[(i-1) % total_ind]+=get_pix(i+ind, cur_j, cur_k, 
				data_event, Border);
		nb_event[Current_Scale](i,0,0)+=tab[(i-1) % total_ind];	
	}
	
       // cout << "Event_scale_3D 3: " << nb_event[Current_Scale].min() << endl;

	
	
	//********************computation of the first plan:  (0:Nx-1 ,0:Ny-1 ,0)****
	
	
	for (int i=0;i<Nx;i++){
		//(0:Nx-1, 1:total_ind, 0)
		for (int j=1;j<=total_ind;j++){
			tab[j-1]=0;
			nb_event[Current_Scale](i,j,0)=
				nb_event[Current_Scale](i,j-1,0);
			for (int cur_i=-ind;cur_i<=ind;cur_i++)
			for (int cur_k=-ind;cur_k<=ind;cur_k++){
				tab[j-1]+=get_pix(i+cur_i, j+ind, cur_k, data_event, 
					Border);
				nb_event[Current_Scale](i,j,0)-=get_pix(i+cur_i, 
					j-ind-1, cur_k, data_event, Border);
			}
			nb_event[Current_Scale](i,j,0) += tab[j-1];
		}
		//(0:Nx-1 ,total_ind+1:Ny, 0)
		for (int j=total_ind+1;j<Ny;j++){
			nb_event[Current_Scale](i,j,0)=nb_event[Current_Scale](i,j-1,0)
				- tab[(j-1) % total_ind];
			tab[(j-1) % total_ind]=0;
			for (int cur_i=-ind;cur_i<=ind;cur_i++)
			for (int cur_k=-ind;cur_k<=ind;cur_k++)
				tab[(j-1) % total_ind] += get_pix(i+cur_i, j+ind, cur_k,
					data_event, Border);   
			nb_event[Current_Scale](i,j,0) += tab[(j-1) % total_ind];
		}
	}	
	// cout << "Event_scale_3D 4: " << nb_event[Current_Scale].min() << endl;	
	
	//********************computation of the cube:  (0:Nx-1 ,0:Ny-1 ,0:Nz-1)********
	
	
	for (int i=0;i<Nx;i++){
	for (int j=0;j<Ny;j++){
		//(0:Nx-1 ,0:Ny-1 ,1:total_ind)
		for (int k=1;k<=total_ind;k++){
			nb_event[Current_Scale](i,j,k)=nb_event[Current_Scale](i,j,k-1);
			tab[k-1]=0;
			for (int cur_i=-ind;cur_i<=ind;cur_i++)
			for (int cur_j=-ind;cur_j<=ind;cur_j++){
				tab[k-1] += get_pix(i+cur_i, j+cur_j, k+ind, data_event, 
					Border);
				nb_event[Current_Scale](i,j,k) -= get_pix(i+cur_i, 
					j+cur_j, k-ind-1, data_event, Border);
			}
			nb_event[Current_Scale](i,j,k) += tab[k-1];
		}
		//(0:Nx-1 ,0:Ny-1 ,total_ind+1:Nz)
		for (int k=total_ind+1;k<Nz;k++){
			nb_event[Current_Scale](i,j,k)=nb_event[Current_Scale](i,j,k-1) 
				- tab[(k-1) % total_ind];
			tab[(k-1) % total_ind]=0;
			for (int cur_i=-ind;cur_i<=ind;cur_i++)
			for (int cur_j=-ind;cur_j<=ind;cur_j++)
		   		tab[(k-1) % total_ind] += get_pix(i+cur_i, j+cur_j, 
					k+ind, data_event, Border);   
			nb_event[Current_Scale](i,j,k) += tab[(k-1) % total_ind];
		}
	}
	// cout << "Event_scale_3D 5: " << nb_event[Current_Scale].min() << endl;	
	if(WriteAllInfo==True) 
	{
	   cout << flush;
	   float pc=i*100/(float)(Nx);
	   pc=((int)(pc*10))/10.;
	   cout << '\r' << "Scale number " << Current_Scale << ", " << pc 
		<< "% done.       ";
	}
	/*switch(i%4){
		case 0: cout <<'\r'<< "-"; break;
		case 1: cout <<'\r'<< "\\"; break;
		case 2: cout <<'\r'<< "|"; break;
		case 3: cout <<'\r'<< "/"; break;
	}*/
	}
	if(WriteAllInfo==True) cout << '\r';
	if(WriteAllInfo==True) cout << "Table of number of events in cube "
	    << total_ind << '*' << total_ind << '*' << total_ind << " is created" <<endl;
	delete[] tab;
}

// ******************************************************************************

/*
long FewEventPoisson::event_one_scale_3D (int i, int j, int k, int s, 
					intarray & data_event, 
					intarray & nb_event, 
					type_border Border) {

   if(s==0) 
   	return(nb_event(i,j,k));
   else{
   	long total;
	int ind= (int) (pow((double)2., (double)(s+1)) + 0.5);
	int ind_2=ind*2;
	
	total =	 event_one_scale_3D(i-ind,j-ind,k-ind,s-1,data_event,nb_event,Border)
		+event_one_scale_3D(i-ind,j-ind,k+ind,s-1,data_event,nb_event,Border)
		+event_one_scale_3D(i-ind,j+ind,k-ind,s-1,data_event,nb_event,Border)
		+event_one_scale_3D(i-ind,j+ind,k+ind,s-1,data_event,nb_event,Border)
		+event_one_scale_3D(i+ind,j-ind,k-ind,s-1,data_event,nb_event,Border)
		+event_one_scale_3D(i+ind,j-ind,k+ind,s-1,data_event,nb_event,Border)
		+event_one_scale_3D(i+ind,j+ind,k-ind,s-1,data_event,nb_event,Border)
		+event_one_scale_3D(i+ind,j+ind,k+ind,s-1,data_event,nb_event,Border);
	
	for(int cur_i=i-ind_2; cur_i <= i+ind_2; cur_i++)
	for(int cur_j=j-ind_2; cur_j <= j+ind_2; cur_j++)
		total -= get_pix(cur_i,cur_j,k,data_event, Border);

	for(int cur_i=i-ind_2; cur_i <= i+ind_2; cur_i++)
	for(int cur_k=k-ind_2; cur_k <= k+ind_2; cur_k++)
		total -= get_pix(cur_i,j,cur_k,data_event, Border);

	for(int cur_j=j-ind_2; cur_j <= j+ind_2; cur_j++)
	for(int cur_k=k-ind_2; cur_k <= k+ind_2; cur_k++)
		total -= get_pix(i,cur_j,cur_k,data_event, Border);
	
	for(int cur_i=i-ind_2; cur_i <= i+ind_2; cur_i++)
		total -= get_pix(cur_i,j,k,data_event, Border);
	
	for(int cur_j=j-ind_2; cur_j <= j+ind_2; cur_j++)
		total -= get_pix(i,cur_j,k,data_event, Border);
	
	for(int cur_k=k-ind_2; cur_k <= k+ind_2; cur_k++)
		total -= get_pix(i,j,cur_k,data_event, Border);
		
	total -= get_pix(i,j,k,data_event, Border);
		
   	return(total);
   }
}*/

// ******************************************************************************
/*
void FewEventPoisson::init_event_all_scales_3D (intarray & data_event, intarray *& nb_event,
	int nb_scale, type_border Border,Bool WriteAllInfo){
	
 	
   init_event_first_scale_3D(data_event, nb_event[0], Border, WriteAllInfo);	
   int Nx = data_event.nx();
   int Ny = data_event.ny();
   int Nz = data_event.nz();
   
   for(int s=1;s<nb_scale;s++){
   	int ind= (int) (pow((double)2., (double)(s+1)) + 0.5);
   	int ind_2=ind*2;
	if(WriteAllInfo==True) cout << "Scale " << s <<" computing :) ..." << endl;
   	for(int i=0;i<Nx;i++){ cout << i << endl;
   	for(int j=0;j<Ny;j++)
   	for(int k=0;k<Nz;k++){
   	    long total=0;
   	    if( i-ind>=0 && j-ind>=0 && k-ind>=0)
	    	total += nb_event[s-1](i-ind,j-ind,k-ind)
		   +nb_event[s-1](i-ind,j-ind,k+ind)
		   +nb_event[s-1](i-ind,j+ind,k-ind)
		   +nb_event[s-1](i-ind,j+ind,k+ind)
		   +nb_event[s-1](i+ind,j-ind,k-ind)
		   +nb_event[s-1](i+ind,j-ind,k+ind)
		   +nb_event[s-1](i+ind,j+ind,k-ind)
		   +nb_event[s-1](i+ind,j+ind,k+ind);
		    
	    for(int cur_i=i-ind_2; cur_i <= i+ind_2; cur_i++)
	    for(int cur_j=j-ind_2; cur_j <= j+ind_2; cur_j++)
		total -= get_pix(cur_i,cur_j,k,data_event, Border);

	    for(int cur_i=i-ind_2; cur_i <= i+ind_2; cur_i++)
	    for(int cur_k=k-ind_2; cur_k <= k+ind_2; cur_k++)
		total -= get_pix(cur_i,j,cur_k,data_event, Border);

	    for(int cur_j=j-ind_2; cur_j <= j+ind_2; cur_j++)
	    for(int cur_k=k-ind_2; cur_k <= k+ind_2; cur_k++)
		total -= get_pix(i,cur_j,cur_k,data_event, Border);
	
	    for(int cur_i=i-ind_2; cur_i <= i+ind_2; cur_i++)
		total -= get_pix(cur_i,j,k,data_event, Border);
	
	    for(int cur_j=j-ind_2; cur_j <= j+ind_2; cur_j++)
		total -= get_pix(i,cur_j,k,data_event, Border);
	
	    for(int cur_k=k-ind_2; cur_k <= k+ind_2; cur_k++)
		total -= get_pix(i,j,cur_k,data_event, Border);
		
	    total -= get_pix(i,j,k,data_event, Border);
	    nb_event[s](i,j,k) = total;
   	}
   	}	
   }
   
}*/

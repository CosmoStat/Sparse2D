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
**    File:  mrs_planck_mgmca.cc
**
*******************************************************************************
**
**    DESCRIPTION  Multiscale Generalized Morphological Component Analysis for planck data
**    ----------- 
**                 
**    Usage: mrs_planck_mgmca options cube output
**
******************************************************************************/

#include "HealpixClass.h"
#include "MRS_Sparse.h"

#include "DefMath.h"
#include "MatrixOper.h"
#include "Array.h"
#include "NR.h"
#include "IM_IO.h"
#include <cmath>
#include <sys/time.h>

 
#define ALM_ITER 0
#define ALM_NRM false
#define ALM_BEAM false

extern int NitALM;
extern bool Nested_in;
extern void  error(int num, char *msg1, char *msg2);

class lGMCA_inv {
	int Nband; //Number of resolution bands
	int NchanMax; //Max Number of observations
	int Nbeams;//Number of beams > NchanMax
	int Nchan; //Number of observations for the current resolution
	int Nside; //Nside for the maps
	unsigned long int Npix; //current number of pixels -> Nside*Nside*12
	int NsideRec; //Nside for the current resolution
	int NsideFace; //Nside used to calculate face indices
	unsigned long int NpixRec; //current number of pixels -> Nside*Nside*12
	double Prec;
	
	//For Threading
	int inner_loop_threads;
	int outer_loop_threads;
	int NumThread;
	
	 /*For WT Transform*/
	 int Nscale; //Nscale for the current resolution
	
	 /*For Alm Transform*/
	int lLmax;
	int weight_Nside;
    int Alm_Niter;
	bool Alm_Norm;
	bool Alm_Beam;
	bool Alm_Fast;
	bool AllocALMMem;
	intarray AlmXComp;
	intarray LmaxChan;
	arr<double> weight_TRing;	
	
	void Set_almtrans(int iter=ALM_ITER, bool Norm=ALM_NRM, bool beam=ALM_BEAM,bool fast=DEF_ALM_FAST) {
    	Alm_Niter=iter;
    	Alm_Norm = Norm;
   		Alm_Beam = beam;
   		Alm_Fast=fast;
   		AlmData.Niter=iter;
   		AlmData.Norm=Norm;
   		AlmData.UseBeamEff=beam; 
   		AllocALMMem=false;  		
    }
	void Set_Weight(int Nside) { 
    	 string DirWeight = string(getenv("HEALPIX")) + string("/data");
    	 weight_TRing.alloc(2*Nside);
    	 if (Alm_Fast == false){
             read_weight_ring ( DirWeight, Nside, weight_TRing);
             for (int m=0; m< (int) weight_TRing.size(); ++m) weight_TRing[m]+=1;
          } else weight_TRing.fill(1);
          weight_Nside=Nside;
	}
	
public:
	fltarray InputMaps; //[Npix,NchanMax] -> Input Maps
	fltarray Beams; //[Npix,Nchan] -> Beams 
	fltarray SameResMaps; //[Npix,Nchan] -> Maps put to the current resolution
	fltarray InversionVectors; //[Nside*Nside*Nface,Nchan,Nscale] -> Weighting for the inversion in wavelet domain
	fltarray ResolStruc; //[NBand,Nstruct]-> Nstruct=[Nchan,ListChan,lmax] eg. [8 2 3 4 5 6 7 8 9 0 1024] -> 8 channels, from indices 1->8, with lmax 1024
	fltarray OutputMap; //[Npix,Nchan] -> Final maps after convolution + weighting in wavelet domain

    CAlmR AlmData; //will contain the SPH. HARM. TRANS routines
	Alm<xcomplex<double> > *AlmX;//Will contain the ALM transform of the input maps 

	C_UWT2D_ATROUS WT_Trans; //will contain the wavelet transform routines for the Nchannels 
	
 	to_array<unsigned long,true> FaceIndices;//Index of each pixel in a face in Healpix Ring format
 	to_array<unsigned long,true> Nest2Ring_indices;//Index of each pixel in a face in Healpix Ring format
 	to_array<unsigned long,true> Ring2Nest_indices;//Index of each pixel in a face in Healpix Ring format
 	
	bool Verbose;
	bool Debug;
	bool Timer;
	bool save_sameres;
	lGMCA_inv() {Nchan=Nside=Npix=NchanMax=0;Nscale=4;Verbose=false;AllocALMMem=false;WT_Trans.UWT2D = new ATROUS_2D_WT(); 
			WT_Trans.AllocMem=1;Set_almtrans();Timer=false;Prec=1e-8;inner_loop_threads=outer_loop_threads=NumThread=1;NsideFace=0;Debug=false;Nbeams=0;save_sameres=false;};
	
	//I/O Routines
	void Read_input_maps(char *Name);	
	void Read_Beams(char *Name);		
	void Read_inversion_vectors(char *Name);		
	void Read_inversion_structure(char *Name);		
	unsigned long int get_Npix() const {return Npix;};
	void set_nb_thr(int oloop,int iloop) {
		inner_loop_threads=iloop;
		outer_loop_threads=oloop;
		NumThread=iloop*oloop;
	};
	
	 //ALM Space routines
    double Get_ALMNormval() const { return AlmData.NormVal;	}
    void Set_lLmax(int lL){ lLmax=lL;}
    int Get_Lmax() const { return lLmax;}
    void Set_lNrm(int Nside){ 
    	double Nl=sqrt( (double) (Nside)* (double) (Nside)*12.);
    	AlmData.NormVal = sqrt( (double)(Nl*(Nl+1.)/(4.* PI)));
    }
	void Set_Same_Resolution(int index_Array,int CChan, int RefBeam, int ChanLmax);
	void create_nest2ringindices(int lNside);
	void create_ring2nestindices(int lNside);
	
	//Only set to same resolution
	void SetOnlyResolution(char *Name_in, char *Name_beams, char *Name_inv_struc, char *Name_Out, int lNside);

	//only inversion
	void Invert_Weights_Only(char *Name_in,char *Name_inv_vectors,char *Name_inv_struc, char *Name_Out);

	//WT routines
	int get_Nscale() const {return Nscale;};
	void WT_invert_single_face(int nf) ;
	void WT_transform_single_map(int Chan);
	void create_faceindices(bool cNest=Nested_in);
	void WT_transform_single_face(fltarray &OutMap,int nf) ;
	
	//Process bands	
	void Process_AllBands(char *Name_in, char *Name_beams,char *Name_inv_vectors,char *Name_inv_struc, char *Name_Out);
	void Process_Band(int band, intarray &listchannel);
	
	~lGMCA_inv() { 	if(AllocALMMem==true) delete [] AlmX; 
 };
};





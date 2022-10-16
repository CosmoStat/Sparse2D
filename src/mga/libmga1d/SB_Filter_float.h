//
// Floating version of the FCur, adapted to FloatTrans class, for
// algorithms like Fista/Nesterov
// This transform is fully normalized
//

#ifndef _SB_FILTER_FLOAT_H
#define _SB_FILTER_FLOAT_H

#include <vector>
#include "SB_Filter.h"
#include "FloatTrans.h"
#include "writefits3d.h"

// ***********************************************************************
// *   Orthogonal 2D wavelet transform
// ***********************************************************************
 
class Ortho_2D_WT_float: public Ortho_2D_WT, public FloatTrans
{
	int DataNx;
	int DataNy;
	int NbrScale;
public:
	Ortho_2D_WT_float(SubBand1D &SB1D):Ortho_2D_WT(SB1D),DataNx(0),DataNy(0),NbrScale(0){}
	~Ortho_2D_WT_float(){}
	
	void alloc(fltarray& in, int _NbrScale);
	void alloc(cfarray& in, int _NbrScale);
	
// FloatTrans methods
	int size_transform();
	void transform(fltarray& in, float* &out, bool alloc=false);
	void recons(float* in, fltarray& out);
	void soft_threshold(float* in, float lvl, bool threshold_coarse=false);
        void hard_threshold(float* in, float lvl, bool threshold_coarse=false);
	void substract_coarse_scale(float* in, float* &coarse, bool alloc_coarse=false);
	void add_coarse_scale(float* in, float* coarse);
};
// Global functions for MEX
void dwt2d_clear(vector< vector< Ifloat* > > &vTabBand);
int dwt2d_transform(Ifloat &Data, Ifloat &vTabBand, int _Nbr_Scale, int wavelet_type);
int dwt2d_transform(Ifloat &Data, vector< vector< Ifloat* > > &vTabBand, int _Nbr_Scale, int wavelet_type);
int dwt2d_recons(Ifloat &vTabBand, Ifloat &Data, int _Nbr_Scale, int wavelet_type);
int dwt2d_recons(vector< vector< Ifloat* > > &vTabBand, Ifloat &Data, int wavelet_type);
int dwt2d_threshold(Ifloat &vTabBand, int _NbrScale, float lvl, int filtering_type);
int dwt2d_threshold(vector< vector< Ifloat* > > &vTabBand, float lvl, int filtering_type);
int dwt2d_filter(Ifloat &Data, Ifloat &Recons, int _Nbr_Scale, float lvl, int filtering_type, int wavelet_type);

// ***********************************************************************
// *   Orthogonal Undecimated 2D wavelet transform
// ***********************************************************************

class PAVE_2D_WT_float: public PAVE_2D_WT, public FloatTrans
{
	int DataNx;
	int DataNy;
	int NbrScale;
    dblarray TabMADNorm;
public:
	PAVE_2D_WT_float(SubBand1D &SB1D):PAVE_2D_WT(SB1D),DataNx(0),DataNy(0),NbrScale(0){}
	~PAVE_2D_WT_float(){}
	
	void alloc(fltarray& in, int _NbrScale);
	void alloc(cfarray& in, int _NbrScale);
	
// FloatTrans methods
	int size_transform();
	void transform(fltarray& in, float* &out, bool alloc=false);
	void recons(float* in, fltarray& out);
    void adjoint_recons(float* in, fltarray& out) {recons(in, out);};
	void soft_threshold(float* in, float lvl, bool threshold_coarse=false);
    void hard_threshold(float* in, float lvl, bool threshold_coarse=false);
	void substract_coarse_scale(float* in, float* &coarse, bool alloc_coarse=false);
	void add_coarse_scale(float* in, float* coarse);
    
    void normalize(float* in){};
    void unnormalize(float* in){};
    void mad_calculation(float* in);
    void mad_normalize(float* in);
    void mad_unnormalize(float* in);
};

// Global functions for MEX
void uwt2d_clear(vector< vector< Ifloat* > > &C);
int uwt2d_transform(Ifloat &Data, vector< vector< Ifloat* > > &vTabBand, int _Nbr_Scale, int wavelet_type);
int uwt2d_recons(vector< vector< Ifloat* > > &vTabBand, Ifloat &Data, int wavelet_type);
int uwt2d_threshold(vector< vector< Ifloat* > > &vTabBand, float lvl, int filtering_type);
int uwt2d_filter(Ifloat &Data, Ifloat &Recons, int _Nbr_Scale, float lvl, int filtering_type, int wavelet_type);


/************************************************************************/
//    Isotropic a trous wavelet transform
/************************************************************************/

class ATROUS_2D_WT_float: public ATROUS_2D_WT, public FloatTrans
{
	int DataNx;
	int DataNy;
	int NbrScale;
	bool positive_coef;
    dblarray TabMADNorm;
public:
	fltarray locations;
	ATROUS_2D_WT_float():ATROUS_2D_WT(),DataNx(0),DataNy(0),NbrScale(0),positive_coef(false){ModifiedAWT=True;}
	~ATROUS_2D_WT_float(){}
	
	void alloc(fltarray& in, int _NbrScale);
	void alloc(cfarray& in, int _NbrScale);
	void set_positive_coef(bool pc){positive_coef = pc;}
	void set_use_adjoint(bool ad=true) {ModifiedAWT=(Bool)ad;}
// FloatTrans methods
	int size_transform();
	void transform(fltarray& in, float* &out, bool alloc=false);
	void recons(float* in, fltarray& out);
    void adjoint_recons(float* in, fltarray& out);
	void soft_threshold(float* in, float lvl, bool threshold_coarse=false);
        void hard_threshold(float* in, float lvl, bool threshold_coarse=false);
	void substract_coarse_scale(float* in, float* &coarse, bool alloc_coarse=false);
	void add_coarse_scale(float* in, float* coarse);
    void normalize(float* in);
    void unnormalize(float* in);
    void mad_calculation(float* in);
    void mad_normalize(float* in);
    void mad_unnormalize(float* in);
};

// Global functions for MEX
void iwt2d_clear(vector< Ifloat* > &C);
int iwt2d_transform(Ifloat &Data, vector< Ifloat* > &vTabBand, int _Nbr_Scale);
int iwt2d_recons(vector< Ifloat* > &vTabBand, Ifloat &Data);
int iwt2d_threshold(vector< Ifloat* > &vTabBand, float lvl, int filtering_type);
int iwt2d_filter(Ifloat &Data, Ifloat &Recons, int _Nbr_Scale, float lvl, int filtering_type);

#endif

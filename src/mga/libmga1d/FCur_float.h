//
// Floating version of the FCur, adapted to FloatTrans class, for
// algorithms like Fista/Nesterov
//

#ifndef _FASTCUR_FLOAT_H
#define _FASTCUR_FLOAT_H

#include <vector>
#include "FCur.h"
#include "FloatTrans.h"

class FCUR_float_params
{
public:
	FCUR_float_params() {reset_params();}
	~FCUR_float_params() {}//{if(TabNorm!=NULL) delete TabNorm;}
	void reset_params();
	
	int Nl;
	int Nc;
	int NbrScale;
	int NbrDir;
	bool Undecimated;
//	fltarray *TabNorm;
}; 

class FCUR_float: public FCUR, public FloatTrans
{
    dblarray TabMADNorm;

public:
	FCUR_float():FCUR(){}
	~FCUR_float(){};//{clear();} // clear is set in FCUR

	void export_trans(float* out);
	void import_trans(float* in);
	void alloc(fltarray& in, int _NbrScale, int _NbrDir);
//	void clear() {for(int i=0;i<nbr_scale();i++) delete [] TabCF_Band[i];}
	
// FloatTrans methods
	int size_transform();
	void transform(fltarray& in, float* &out, bool alloc=false);
	void recons(float* in, fltarray& out);
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
//int fct_transform(Ifloat &Data, vector< vector< Ifloat* > > &vTabBand, int _Nbr_Scale, int _NbrDir);
//int fct_recons(vector< vector< Ifloat* > > &vTabBand, Ifloat &Data, int Nl, int Nc, int _Nbr_Scale, int _NbrDir);
//int fct_threshold(vector< vector< Ifloat* > > &vTabBand, float lvl, int filtering_type);
//int fct_filter(Ifloat &Data, Ifloat &Recons, int _Nbr_Scale, int _NbrDir, float lvl, int filtering_type);
void fct_clear(vector< vector< Ifloat* > > &C);
int fct_get_norm(Ifloat &Data, FCUR_float_params &P);
int fct_transform(Ifloat &Data, vector< vector< Ifloat* > > &vTabBand, FCUR_float_params &P);
int fct_recons(vector< vector< Ifloat* > > &vTabBand, Ifloat &Data, FCUR_float_params &P);
int fct_threshold(vector< vector< Ifloat* > > &vTabBand, float lvl, int filtering_type);
int fct_filter(Ifloat &Data, Ifloat &Recons, FCUR_float_params &P, float lvl, int filtering_type);


/***********************************************************************/

#endif

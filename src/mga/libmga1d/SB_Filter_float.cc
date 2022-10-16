

#include "SB_Filter_float.h"

// ***********************************************************************
// *   Orthogonal 2D wavelet transform
// ***********************************************************************
 
void Ortho_2D_WT_float::alloc(fltarray& in, int _NbrScale)
{
	DataNx = in.nx();
	DataNy = in.ny();
	NbrScale = _NbrScale;
}

void Ortho_2D_WT_float::alloc(cfarray& in, int _NbrScale)
{
        DataNx = in.nx();
        DataNy = in.ny();
        NbrScale = _NbrScale;
}

int Ortho_2D_WT_float::size_transform()
{	
	return DataNx*DataNy;
}

void Ortho_2D_WT_float::transform(fltarray& in, float* &out, bool alloc)
{
	if(alloc) out = new float [size_transform()];
	Ifloat f_in; f_in.alloc(in.buffer(),DataNy,DataNx);
	Ifloat f_out; f_out.alloc(out,DataNy,DataNx);
	Ortho_2D_WT::transform (f_in, f_out, NbrScale);
}
void Ortho_2D_WT_float::recons(float* in, fltarray& out)
{
	out.resize(DataNx,DataNy,0);
	Ifloat f_in; f_in.alloc(in,DataNy,DataNx);
	Ifloat f_out; f_out.alloc(out.buffer(),DataNy,DataNx);
	Ortho_2D_WT::recons(f_in, f_out, NbrScale);
}
void Ortho_2D_WT_float::soft_threshold(float* in, float lvl, bool threshold_coarse)
{// threshold_coarse=false not implemented yet
	Ifloat f_in; f_in.alloc(in,DataNy,DataNx);
	int Nl = DataNy;
	int Nc = DataNx;
	for(int s=0;s<NbrScale-1;s++)
	{
		Nl = (Nl+1)/2;
		Nc = (Nc+1)/2;
	}
	for(int i=0;i<DataNy;i++)
	for(int j=0;j<DataNx;j++)
		if(i>=Nl*int(!threshold_coarse) || j>=Nc*int(!threshold_coarse))
			f_in(i,j) = ::soft_threshold(f_in(i,j), lvl);
//	for(int i=0;i<size_transform();i++)
//		in[i] = ::soft_threshold(in[i], lvl);
}
void Ortho_2D_WT_float::hard_threshold(float* in, float lvl, bool threshold_coarse)
{// threshold_coarse=false not implemented yet
	Ifloat f_in; f_in.alloc(in,DataNy,DataNx);
	int Nl = DataNy;
	int Nc = DataNx;
	for(int s=0;s<NbrScale-1;s++)
	{
		Nl = (Nl+1)/2;
		Nc = (Nc+1)/2;
	}
	for(int i=0;i<DataNy;i++)
	for(int j=0;j<DataNx;j++)
		if(i>=Nl*int(!threshold_coarse) || j>=Nc*int(!threshold_coarse))
			f_in(i,j) = ::hard_threshold(f_in(i,j), lvl);
//	for(int i=0;i<size_transform();i++)
//		in[i] = ::soft_threshold(in[i], lvl);
}
void Ortho_2D_WT_float::substract_coarse_scale(float* in, float* &coarse, bool alloc_coarse)
{
	
}
void Ortho_2D_WT_float::add_coarse_scale(float* in, float* coarse)
{
	
}

int dwt2d_transform(Ifloat &Data, Ifloat &vTabBand, int _Nbr_Scale, int wavelet_type)
{
// Wavelet initialisation
	FilterAnaSynt SelectFilter;
	SelectFilter.alloc((type_sb_filter)wavelet_type);
	SubBandFilter *SB1D = new SubBandFilter(SelectFilter, NORM_L2);
	Ortho_2D_WT_float *DataW = new Ortho_2D_WT_float(*SB1D);

// Forward Transform
	DataW->Ortho_2D_WT::transform(Data,vTabBand,_Nbr_Scale);

	delete DataW;
	return 1;
}

void dwt2d_clear(vector< vector< Ifloat* > > &C)
{
	for(int i=0;i<C.size();i++)
	{
		for(int j=0;j<C[i].size();j++)
			delete C[i][j];
		C[i].clear();
	}
	C.clear();
}

int dwt2d_transform(Ifloat &Data, vector< vector< Ifloat* > > &vTabBand, int _Nbr_Scale, int wavelet_type)
{
	Ifloat TabBand;
	dwt2d_transform(Data, TabBand, _Nbr_Scale, wavelet_type);
	vTabBand.resize(_Nbr_Scale);
	int Nl = Data.nl(); int nl;
	int Nc = Data.nc(); int nc;
	for(int s=0;s<_Nbr_Scale;s++)
	{
		// 3 directions per scale
		int nb = s == _Nbr_Scale-1 ? 1 : 3;
		vTabBand[s].resize(nb);
		// next scale's size
		nl = (Nl+1)/2;
		nc = (Nc+1)/2;
		if(s<_Nbr_Scale-1)
		{
			vTabBand[s][0] = new Ifloat(Nl-nl,nc);
			for (int i=0; i < Nl-nl; i++)
			for (int j=0; j < nc; j++)
				(*vTabBand[s][0])(i,j) = TabBand(i,j+nl);
			vTabBand[s][1] = new Ifloat(nl,Nc-nc);
			for (int i=0; i < nl; i++)
			for (int j=0; j < Nc-nc; j++)
				(*vTabBand[s][1])(i,j) = TabBand(i+nc,j);
			vTabBand[s][2] = new Ifloat(Nl-nl,Nc-nc);
			for (int i=0; i < Nl-nl; i++)
			for (int j=0; j < Nc-nc; j++)
				(*vTabBand[s][2])(i,j) = TabBand(i+nc,j+nl);
		}
		else
		{
			vTabBand[s][0] = new Ifloat(Nl,Nc);
			for (int i=0; i < Nl; i++)
			for (int j=0; j < Nc; j++)
				(*vTabBand[s][0])(i,j) = TabBand(i,j);
		}
		Nl = nl; Nc = nc;
	}
	return 1;
}

int dwt2d_recons(Ifloat &vTabBand, Ifloat &Recons, int _Nbr_Scale, int wavelet_type)
{
// Wavelet initialisation
	FilterAnaSynt SelectFilter;
	SelectFilter.alloc((type_sb_filter)wavelet_type);
	SubBandFilter *SB1D = new SubBandFilter(SelectFilter, NORM_L2);
	Ortho_2D_WT_float *DataW = new Ortho_2D_WT_float(*SB1D);

// Backward Transform
	DataW->Ortho_2D_WT::recons(vTabBand, Recons, _Nbr_Scale);

	delete DataW;
	return 1;
}

int dwt2d_recons(vector< vector< Ifloat* > > &vTabBand, Ifloat &Recons, int wavelet_type)
{
	int _Nbr_Scale = vTabBand.size();
// Data size calculation
	int Nl=0, nl, Nc=0, nc;
	for(int s=0;s<_Nbr_Scale-1;s++)
	{
		Nl += vTabBand[s][2]->nl();
		Nc += vTabBand[s][2]->nc();
	}
	Nl += vTabBand[_Nbr_Scale-1][0]->nl();
	Nc += vTabBand[_Nbr_Scale-1][0]->nc();
		
// Data recopy
	Ifloat TabBand(Nl,Nc);
	for(int s=0;s<_Nbr_Scale;s++)
	{
		// next scale's size
		nl = (Nl+1)/2;
		nc = (Nc+1)/2;
		if(s<_Nbr_Scale-1)
		{
			for (int i=0; i < Nl-nl; i++)
			for (int j=0; j < nc; j++)
				TabBand(i,j+nl) = (*vTabBand[s][0])(i,j);
			for (int i=0; i < nl; i++)
			for (int j=0; j < Nc-nc; j++)
				TabBand(i+nc,j) = (*vTabBand[s][1])(i,j);
			for (int i=0; i < Nl-nl; i++)
			for (int j=0; j < Nc-nc; j++)
				TabBand(i+nc,j+nl) = (*vTabBand[s][2])(i,j);
		}
		else
		{
			for (int i=0; i < Nl; i++)
			for (int j=0; j < Nc; j++)
				TabBand(i,j) = (*vTabBand[s][0])(i,j);
		}
		Nl = nl; Nc = nc;
	}
	dwt2d_recons(TabBand, Recons, _Nbr_Scale, wavelet_type);

	return 1;
}

int dwt2d_threshold(Ifloat &vTabBand, int _Nbr_Scale, float lvl, int filtering_type)
{
// Thresholding
	int nl = vTabBand.nl();
	int nc = vTabBand.nc();
	int i0 = nl;
	int j0 = nc;
	
	if(filtering_type==1)// Hard Thresholding
		for(int s=0;s<_Nbr_Scale-1;s++)
		{
			int in = i0;
			int jn = j0;
			i0 = (in+1)/2;
			j0 = (jn+1)/2;
			
			for(int b=0;b<3;b++)
				for (int i =  (b != 0)? i0 : 0 ; i < ((b != 0) ? in : i0) ; i++)
				for (int j = (b % 2 == 0) ? j0 : 0 ; j < ((b % 2 == 0) ? jn : j0) ; j++)
				{
					if(abs(vTabBand(i,j)) < lvl)
						vTabBand(i,j) = 0.F;
				}
		}
	else if(filtering_type==2) // soft thresholding
		for(int s=0;s<_Nbr_Scale-1;s++)
		{
			int in = i0;
			int jn = j0;
			i0 = (in+1)/2;
			j0 = (jn+1)/2;
			
			for(int b=0;b<3;b++)
				for (int i =  (b != 0)? i0 : 0 ; i < ((b != 0) ? in : i0) ; i++)
				for (int j = (b % 2 == 0) ? j0 : 0 ; j < ((b % 2 == 0) ? jn : j0) ; j++)
					vTabBand(i,j) = ::soft_threshold(vTabBand(i,j),lvl);
		}
	else cerr<<"Filtering method not implemented yet"<<endl;

	return 1;
}

int dwt2d_threshold(vector< vector< Ifloat* > > &vTabBand, float lvl, int filtering_type)
{
	int _Nbr_Scale = vTabBand.size();
// Data size calculation
	int Nl=0, nl, Nc=0, nc;
	for(int s=0;s<_Nbr_Scale-1;s++)
	{
		Nl += vTabBand[s][2]->nl();
		Nc += vTabBand[s][2]->nc();
	}
	Nl += vTabBand[_Nbr_Scale-1][0]->nl();
	Nc += vTabBand[_Nbr_Scale-1][0]->nc();
	
	if(filtering_type==1)// Hard Thresholding
		for(int s=0;s<_Nbr_Scale-1;s++)
		{
			// next scale's size
			nl = (Nl+1)/2;
			nc = (Nc+1)/2;
			for(int b=0;b<3;b++)
			for (int i=0; i < ((b%2) ? nl : Nl-nl); i++)
			for (int j=0; j < (b ? Nc-nc : nc); j++)
				if(abs((*vTabBand[s][b])(i,j)) < lvl)
					(*vTabBand[s][b])(i,j) = 0.F;
			
			Nl = nl; Nc = nc;
		}
	else if(filtering_type==2) // soft thresholding
		for(int s=0;s<_Nbr_Scale-1;s++)
		{
			// next scale's size
			nl = (Nl+1)/2;
			nc = (Nc+1)/2;
			for(int b=0;b<3;b++)
			for (int i=0; i < ((b%2) ? nl : Nl-nl); i++)
			for (int j=0; j < (b ? Nc-nc : nc); j++)
					(*vTabBand[s][b])(i,j) = ::soft_threshold((*vTabBand[s][b])(i,j),lvl);
			Nl = nl; Nc = nc;
		}
	else cerr<<"Filtering method not implemented yet"<<endl;	
	return 1;
}

int dwt2d_filter(Ifloat &Data, Ifloat &Recons, int _Nbr_Scale, float lvl, int filtering_type, int wavelet_type)
{
// Wavelet initialisation
	FilterAnaSynt SelectFilter;
	SelectFilter.alloc((type_sb_filter)wavelet_type);
	SubBandFilter *SB1D = new SubBandFilter(SelectFilter, NORM_L2);
	Ortho_2D_WT_float *DataW = new Ortho_2D_WT_float(*SB1D);

	Ifloat vTabBand;
	
// Filtering
	DataW->Ortho_2D_WT::transform(Data,vTabBand,_Nbr_Scale);
	dwt2d_threshold(vTabBand, _Nbr_Scale, lvl, filtering_type);
	DataW->Ortho_2D_WT::recons(vTabBand, Recons, _Nbr_Scale);
	
	delete DataW;
	return 1;
}


// ***********************************************************************
// *   Orthogonal Undecimated 2D wavelet transform
// ***********************************************************************

void PAVE_2D_WT_float::alloc(fltarray& in, int _NbrScale)
{
	DataNx = in.nx();
	DataNy = in.ny();
	NbrScale = _NbrScale;
}

void PAVE_2D_WT_float::alloc(cfarray& in, int _NbrScale)
{
        DataNx = in.nx();
        DataNy = in.ny();
        NbrScale = _NbrScale;
}

int PAVE_2D_WT_float::size_transform()
{
	return DataNx*DataNy*(3*NbrScale-2);
}

void PAVE_2D_WT_float::transform(fltarray& in, float* &out, bool alloc)
{
	if(alloc) out = new float [size_transform()];
	int nplan = 3*NbrScale-2;
	Ifloat f_in; f_in.alloc(in.buffer(),DataNy,DataNx);
	Ifloat * tab_f_out = new Ifloat[nplan];
	for(int i=0;i<nplan;i++)
		tab_f_out[i].alloc(out+i*DataNy*DataNx,DataNy,DataNx);
	
	PAVE_2D_WT::transform(f_in, tab_f_out, NbrScale);
	delete [] tab_f_out;
}

void PAVE_2D_WT_float::recons(float* in, fltarray& out)
{
	int nplan = 3*NbrScale-2;
	out.resize(DataNx,DataNy,0);
	Ifloat * tab_f_in = new Ifloat[nplan];
	for(int i=0;i<nplan;i++)
		tab_f_in[i].alloc(in+i*DataNy*DataNx,DataNy,DataNx);
	Ifloat f_out; f_out.alloc(out.buffer(),DataNy,DataNx);
	
	PAVE_2D_WT::recons(tab_f_in, f_out, NbrScale);
	delete [] tab_f_in;
}

void PAVE_2D_WT_float::mad_calculation(float* in)
{
     int nplan = 3*NbrScale-2;
     if (TabMADNorm.nx() != nplan) TabMADNorm.alloc(nplan);
     Ifloat * tab_f_in = new Ifloat[nplan];
     for(int i=0;i<nplan;i++)
         tab_f_in[i].alloc(in+i*DataNy*DataNx,DataNy,DataNx);
          
     for(int i=0;i<nplan;i++)
     {
         TabMADNorm(i) = get_sigma_mad(tab_f_in[i].buffer(), DataNy*DataNx);
         cout << "Mad scale " << i+1 << " = " << TabMADNorm(i) << ", sig = " << tab_f_in[i].sigma() << endl;
     }
     delete [] tab_f_in;
 }
 
 
void PAVE_2D_WT_float::mad_normalize(float* in)
{
    int nplan = 3*NbrScale-2;
    Ifloat * tab_f_in = new Ifloat[nplan];
    for(int i=0;i<nplan;i++)
         tab_f_in[i].alloc(in+i*DataNy*DataNx,DataNy,DataNx);
     
    for(int i=0;i<nplan;i++)
    {
        for (int i1=0; i1 < tab_f_in[i].nl(); i1++)
        for (int i2=0; i2 < tab_f_in[i].nc(); i2++)
            tab_f_in[i](i1,i2) /= TabMADNorm(i);
    }
    delete [] tab_f_in;
}
 
void PAVE_2D_WT_float::mad_unnormalize(float* in)
{
    int nplan = 3*NbrScale-2;
    Ifloat * tab_f_in = new Ifloat[nplan];
    for(int i=0;i<nplan;i++)
        tab_f_in[i].alloc(in+i*DataNy*DataNx,DataNy,DataNx);
    
    for(int i=0;i<nplan;i++)
    {
        for (int i1=0; i1 < tab_f_in[i].nl(); i1++)
            for (int i2=0; i2 < tab_f_in[i].nc(); i2++)
                tab_f_in[i](i1,i2) *= TabMADNorm(i);
    }
    delete [] tab_f_in;
}


void PAVE_2D_WT_float::soft_threshold(float* in, float lvl, bool threshold_coarse)
{
	int nplan = 3*NbrScale-2;
	for(int i=0;i<(nplan-int(!threshold_coarse))*DataNx*DataNy;i++)
		in[i] = ::soft_threshold(in[i], lvl);
}


void PAVE_2D_WT_float::hard_threshold(float* in, float lvl, bool threshold_coarse)
{
	int nplan = 3*NbrScale-2;
	for(int i=0;i<(nplan-int(!threshold_coarse))*DataNx*DataNy;i++)
		in[i] = ::hard_threshold(in[i], lvl);
}



void PAVE_2D_WT_float::substract_coarse_scale(float* in, float* &coarse, bool alloc_coarse)
{
	int start = 3*(NbrScale-1)*DataNx*DataNy;
	if(alloc_coarse) coarse = new float[DataNx*DataNy];
	for(int i=start;i<size_transform();i++)
	{
		coarse[i-start] = in[i];
		in[i] = 0.;
	}
}
void PAVE_2D_WT_float::add_coarse_scale(float* in, float* coarse)
{
	int start = 3*(NbrScale-1)*DataNx*DataNy;
	for(int i=start;i<size_transform();i++)
		in[i] += coarse[i-start];
}

// Global functions for MEX
int uwt2d_transform(Ifloat &Data, vector< vector< Ifloat* > > &vTabBand, int _Nbr_Scale, int wavelet_type)
{
	int nl = Data.nl();
	int nc = Data.nc();

// Wavelet initialisation
	FilterAnaSynt SelectFilter;
	SelectFilter.alloc((type_sb_filter)wavelet_type);
	SubBandFilter *SB1D = new SubBandFilter(SelectFilter, NORM_L2);
	PAVE_2D_WT_float *DataW = new PAVE_2D_WT_float(*SB1D);
	Ifloat * tab_f_out;
	DataW->PAVE_2D_WT::alloc(tab_f_out, nl, nc, _Nbr_Scale);

// Forward Transform
	DataW->PAVE_2D_WT::transform(Data, tab_f_out, _Nbr_Scale);

// vTabBand allocation
	vTabBand.resize(_Nbr_Scale);
	for(int s=0;s<_Nbr_Scale;s++)
	{
		// 3 directions per scale
		int nb = s == _Nbr_Scale-1 ? 1 : 3;
		vTabBand[s].resize(nb);
		for(int b=0;b<nb;b++)
			vTabBand[s][b] = &tab_f_out[3*s+b];
	}
	delete DataW;
	return 1;
}

void uwt2d_clear(vector< vector< Ifloat* > > &C)
{
	// coresponds to PAVE_WT_2D::alloc
	delete [] C[0][0];
}	

int uwt2d_recons(vector< vector< Ifloat* > > &vTabBand, Ifloat &Recons, int wavelet_type)
{
	int nl = vTabBand[0][0]->nl();
	int nc = vTabBand[0][0]->nc();
	int _Nbr_Scale = vTabBand.size();
	
// Wavelet initialisation
	FilterAnaSynt SelectFilter;
	SelectFilter.alloc((type_sb_filter)wavelet_type);
	SubBandFilter *SB1D = new SubBandFilter(SelectFilter, NORM_L2);
	PAVE_2D_WT_float *DataW = new PAVE_2D_WT_float(*SB1D);
	int nplan = 3*_Nbr_Scale-2;

	Ifloat * tab_f_out = new Ifloat[nplan];
	for(int s=0;s<_Nbr_Scale;s++)
	{
		// 3 directions per scale
		int nb = s == _Nbr_Scale-1 ? 1 : 3;
		for(int b=0;b<nb;b++)
			tab_f_out[3*s+b] = *vTabBand[s][b]; // recopy as vTabBand are not contiguous
	}

// Backward Transform
	Recons.alloc(nl,nc);
	DataW->PAVE_2D_WT::recons(tab_f_out, Recons, _Nbr_Scale);
	
	delete [] tab_f_out;
	delete DataW;
	return 1;
}
		
int uwt2d_threshold(vector< vector< Ifloat* > > &vTabBand, float lvl, int filtering_type)
{
// Thresholding
	int ns = vTabBand.size();
	
	if(filtering_type==1)// Hard Thresholding
		for(int s=0;s<ns-1;s++)
			for(int b=0;b<3;b++)
			{
				Ifloat *A = vTabBand[s][b];
				for (int i=0; i < A->n_elem(); i++)
					if(abs((*A)(i)) < lvl)
						(*A)(i) = 0.F;
			}
	else if(filtering_type==2) // soft thresholding
		for(int s=0;s<ns-1;s++)
			for(int b=0;b<3;b++)
			{
				Ifloat *A = vTabBand[s][b];
				for (int i=0; i < A->n_elem(); i++)
					(*A)(i) = ::soft_threshold((*A)(i),lvl);
			}
	else cerr<<"Filtering method not implemented yet"<<endl;
	
	return 1;
}
		
int uwt2d_filter(Ifloat &Data, Ifloat &Recons, int _Nbr_Scale, float lvl, int filtering_type, int wavelet_type)
{
	int nl = Data.nl();
	int nc = Data.nc();

// Wavelet initialisation
	FilterAnaSynt SelectFilter;
	SelectFilter.alloc((type_sb_filter)wavelet_type);
	SubBandFilter *SB1D = new SubBandFilter(SelectFilter, NORM_L2);
	PAVE_2D_WT_float *DataW = new PAVE_2D_WT_float(*SB1D);
	Ifloat * tab_f_out;
	DataW->PAVE_2D_WT::alloc(tab_f_out, nl, nc, _Nbr_Scale);

// Forward Transform
	DataW->PAVE_2D_WT::transform(Data, tab_f_out, _Nbr_Scale);

// thresholding
	if(filtering_type==1)// Hard Thresholding
		for(int s=0;s<_Nbr_Scale-1;s++)
			for(int b=0;b<3;b++)
			{
				Ifloat *A = tab_f_out + 3*s+b;
				for (int i=0; i < A->n_elem(); i++)
					if(abs((*A)(i)) < lvl)
						(*A)(i) = 0.F;
			}
	else if(filtering_type==2) // soft thresholding
		for(int s=0;s<_Nbr_Scale-1;s++)
			for(int b=0;b<3;b++)
			{
				Ifloat *A = tab_f_out + 3*s+b;
				for (int i=0; i < A->n_elem(); i++)
					(*A)(i) = ::soft_threshold((*A)(i),lvl);
			}
	else cerr<<"Filtering method not implemented yet"<<endl;

// Backward Transform
	Recons.alloc(nl,nc);
	DataW->PAVE_2D_WT::recons(tab_f_out, Recons, _Nbr_Scale);

	delete DataW;
	delete [] tab_f_out;
	return 1;
}


// ***********************************************************************
// *   Isotropic a trous wavelet transform
// ***********************************************************************

void ATROUS_2D_WT_float::alloc(fltarray& in, int _NbrScale)
{
	DataNx = in.nx();
	DataNy = in.ny();
	NbrScale = _NbrScale;
}

void ATROUS_2D_WT_float::alloc(cfarray& in, int _NbrScale)
{
        DataNx = in.nx();
        DataNy = in.ny();
        NbrScale = _NbrScale;
}


int ATROUS_2D_WT_float::size_transform()
{
	return DataNx*DataNy*NbrScale;
}

void ATROUS_2D_WT_float::transform(fltarray& in, float* &out, bool alloc)
{
	if(alloc) out = new float [size_transform()];
	Ifloat f_in; f_in.alloc(in.buffer(),DataNy,DataNx);
	Ifloat * tab_f_out = new Ifloat[NbrScale];
	for(int i=0;i<NbrScale;i++)
		tab_f_out[i].alloc(out+i*DataNy*DataNx,DataNy,DataNx);
	
	ATROUS_2D_WT::transform(f_in, tab_f_out, NbrScale);

// Normalization	
//	for(int i=0;i<size_transform();i++)
//	{
//		float norm = float(norm_band(int(floor(i/(DataNx*DataNy)))));
//		out[i] /= norm;
//	}
	
	delete [] tab_f_out;
}
void ATROUS_2D_WT_float::recons(float* in, fltarray& out)
{
	out.resize(DataNx,DataNy,0);
	Ifloat * tab_f_in = new Ifloat[NbrScale];
	for(int i=0;i<NbrScale;i++)
		tab_f_in[i].alloc(in+i*DataNy*DataNx,DataNy,DataNx);
	Ifloat f_out; 
    f_out.alloc(out.buffer(),DataNy,DataNx);
	
// Normalization	
//	for(int i=0;i<size_transform();i++)
//	{
//		float norm = float(norm_band(int(floor(i/(DataNx*DataNy)))));
//		in[i] *= norm;
//	}
	
	ATROUS_2D_WT::recons(tab_f_in, f_out, NbrScale);

// re-normalization	(because in may be used again ouside this function)
//	for(int i=0;i<size_transform();i++)
//	{
//		float norm = float(norm_band(int(floor(i/(DataNx*DataNy)))));
//		in[i] /= norm;
//	}
	delete [] tab_f_in;
}

void ATROUS_2D_WT_float::adjoint_recons(float* in, fltarray& out)
{
	out.resize(DataNx,DataNy,0);
	Ifloat * tab_f_in = new Ifloat[NbrScale];
	for(int i=0;i<NbrScale;i++)
		tab_f_in[i].alloc(in+i*DataNy*DataNx,DataNy,DataNx);
	Ifloat f_out; 
    f_out.alloc(out.buffer(),DataNy,DataNx);
	
    // Normalization	
    //	for(int i=0;i<size_transform();i++)
    //	{
    //		float norm = float(norm_band(int(floor(i/(DataNx*DataNy)))));
    //		in[i] *= norm;
    //	}
	
	ATROUS_2D_WT::adjoint_recons(tab_f_in, f_out, NbrScale);
    
    // re-normalization	(because in may be used again ouside this function)
    //	for(int i=0;i<size_transform();i++)
    //	{
    //		float norm = float(norm_band(int(floor(i/(DataNx*DataNy)))));
    //		in[i] /= norm;
    //	}
	delete [] tab_f_in;
}


void ATROUS_2D_WT_float::unnormalize(float* in)
{
 	
    // Normalization	
    for(int i=0;i<size_transform();i++)
    {
        float norm = float(norm_band(int(floor(i/(DataNx*DataNy)))));
    	in[i] *= norm;
    }
}

void ATROUS_2D_WT_float::normalize(float* in)
{
    // re-normalization	(because in may be used again ouside this function)
	for(int i=0;i<size_transform();i++)
	{
		float norm = float(norm_band(int(floor(i/(DataNx*DataNy)))));
		in[i] /= norm;
	}
}


void ATROUS_2D_WT_float::soft_threshold(float* in, float lvl, bool threshold_coarse)
{
	// There seems to be problems when thresholding the coarse scale
	threshold_coarse=false;
	for(int i=0;i<(NbrScale-int(!threshold_coarse))*DataNx*DataNy;i++)
		in[i] = ::soft_threshold(in[i], lvl);
	if(positive_coef)
	for(int i=0;i<size_transform();i++)
		if(in[i]<0) in[i] = 0;
}
void ATROUS_2D_WT_float::hard_threshold(float* in, float lvl, bool threshold_coarse)
{
        // There seems to be problems when thresholding the coarse scale - AW.
        threshold_coarse=false;
        for(int i=0;i<(NbrScale-int(!threshold_coarse))*DataNx*DataNy;i++) 
                in[i] = ::hard_threshold(in[i], lvl);
        if(positive_coef)
        for(int i=0;i<size_transform();i++)
                if(in[i]<0) in[i] = 0;
}
void ATROUS_2D_WT_float::substract_coarse_scale(float* in, float* &coarse, bool alloc_coarse)
{
//	int start = (NbrScale-1)*DataNx*DataNy;
//	if(alloc_coarse) coarse = new float[DataNx*DataNy];
//	for(int i=start;i<size_transform();i++)
//	{
//		coarse[i-start] = in[i];
//		in[i] = 0.;
//	}
}
void ATROUS_2D_WT_float::add_coarse_scale(float* in, float* coarse)
{
//	int start = (NbrScale-1)*DataNx*DataNy;
//	for(int i=start;i<size_transform();i++)
//		in[i] += coarse[i-start];
}

void ATROUS_2D_WT_float::mad_calculation(float* in)
{
    if (TabMADNorm.nx() != NbrScale) TabMADNorm.alloc(NbrScale);
    Ifloat * tab_f_in = new Ifloat[NbrScale];
   //  cout << "Nscale = " << NbrScale << endl;
	for(int i=0;i<NbrScale;i++)
		tab_f_in[i].alloc(in+i*DataNy*DataNx,DataNy,DataNx);
   //  cout << "Nscale1 = " << NbrScale << endl;
  
	for(int i=0;i<NbrScale;i++)
    {
        TabMADNorm(i) = get_sigma_mad(tab_f_in[i].buffer(), DataNy*DataNx);
        // cout << "Mad scale " << i+1 << " = " << TabMADNorm(i) << ", sig = " << tab_f_in[i].sigma() << endl;
    }
	delete [] tab_f_in;
}

void ATROUS_2D_WT_float::mad_normalize(float* in)
{
    Ifloat * tab_f_in = new Ifloat[NbrScale];
	for(int i=0;i<NbrScale;i++)
		tab_f_in[i].alloc(in+i*DataNy*DataNx,DataNy,DataNx);
    
	for(int i=0;i<NbrScale;i++)
    {
        for (int i1=0; i1 < tab_f_in[i].nl(); i1++)
        for (int i2=0; i2 < tab_f_in[i].nc(); i2++)
           tab_f_in[i](i1,i2) /= TabMADNorm(i);
    }
 	delete [] tab_f_in;
}

void ATROUS_2D_WT_float::mad_unnormalize(float* in)
{
    Ifloat * tab_f_in = new Ifloat[NbrScale];
	for(int i=0;i<NbrScale;i++)
		tab_f_in[i].alloc(in+i*DataNy*DataNx,DataNy,DataNx);
    
	for(int i=0;i<NbrScale;i++)
    {
        for (int i1=0; i1 < tab_f_in[i].nl(); i1++)
            for (int i2=0; i2 < tab_f_in[i].nc(); i2++)
                tab_f_in[i](i1,i2) *= TabMADNorm(i);
    }
 	delete [] tab_f_in;
}




// Global functions for MEX
void iwt2d_clear(vector< Ifloat* > &C)
{
	// coresponds to PAVE_WT_2D::alloc
	delete [] C[0];
}	

int iwt2d_transform(Ifloat &Data, vector< Ifloat* > &vTabBand, int _Nbr_Scale)
{
	int nl = Data.nl();
	int nc = Data.nc();

// Wavelet initialisation
	ATROUS_2D_WT_float *DataW = new ATROUS_2D_WT_float();
	DataW->AdjointRec = True;
	Ifloat * tab_f_out = new Ifloat[_Nbr_Scale];
	for(int i=0;i<_Nbr_Scale;i++)
		tab_f_out[i].alloc(nl,nc);

// Forward Transform
	DataW->ATROUS_2D_WT::transform(Data, tab_f_out, _Nbr_Scale);

// Normalization	
	for(int s=0;s<_Nbr_Scale-1;s++)
	{
		float norm = float(DataW->norm_band(s));
		float * a = tab_f_out[s].buffer();
		for(int i=0;i<nl*nc;i++)
			 *(a+i) /= norm;
	}

// vTabBand allocation
	vTabBand.resize(_Nbr_Scale);
	for(int s=0;s<_Nbr_Scale;s++)
		vTabBand[s] = &tab_f_out[s];

	delete DataW;
	return 1;
}
		
int iwt2d_recons(vector< Ifloat* > &vTabBand, Ifloat &Recons)
{
	int nl = vTabBand[0]->nl();
	int nc = vTabBand[0]->nc();
	int _Nbr_Scale = vTabBand.size();
	
// Wavelet initialisation
	ATROUS_2D_WT_float *DataW = new ATROUS_2D_WT_float();
	DataW->AdjointRec = True;
	Ifloat * tab_f_out = new Ifloat[_Nbr_Scale];

	for(int s=0;s<_Nbr_Scale;s++)
		tab_f_out[s] = *vTabBand[s];// recopy as vTabBand is not necessaryly contiguous

// Normalization	
	for(int s=0;s<_Nbr_Scale-1;s++)
	{
		float norm = float(DataW->norm_band(s));
		float * a = tab_f_out[s].buffer();
		for(int i=0;i<nl*nc;i++)
			 *(a+i) *= norm;
	}

// Backward Transform
	Recons.alloc(nl,nc);
	DataW->ATROUS_2D_WT::recons(tab_f_out, Recons, _Nbr_Scale);
	
	delete [] tab_f_out;
	delete DataW;
	return 1;
}
		
int iwt2d_threshold(vector< Ifloat* > &vTabBand, float lvl, int filtering_type)
{
	int ns = vTabBand.size();
	
// Thresholding
	if(filtering_type==1)// Hard Thresholding
		for(int s=0;s<ns-1;s++)
		{
			Ifloat *A = vTabBand[s];
			for (int i=0; i < A->n_elem(); i++)
				if(abs((*A)(i)) < lvl)
					(*A)(i) = 0.F;
		}
	else if(filtering_type==2) // soft thresholding
		for(int s=0;s<ns-1;s++)
		{
			Ifloat *A = vTabBand[s];
			for (int i=0; i < A->n_elem(); i++)
				(*A)(i) = ::soft_threshold((*A)(i),lvl);
		}
	else cerr<<"Filtering method not implemented yet"<<endl;
	
	return 1;
}
		
int iwt2d_filter(Ifloat &Data, Ifloat &Recons, int _Nbr_Scale, float lvl, int filtering_type)
{
	int nl = Data.nl();
	int nc = Data.nc();

// Wavelet initialisation
	ATROUS_2D_WT_float *DataW = new ATROUS_2D_WT_float();
	DataW->AdjointRec = True;
	Ifloat * tab_f_out = new Ifloat[_Nbr_Scale];
	for(int i=0;i<_Nbr_Scale;i++)
		tab_f_out[i].alloc(nl,nc);

// Forward Transform
	DataW->ATROUS_2D_WT::transform(Data, tab_f_out, _Nbr_Scale);

// Filtering
	if(filtering_type==1)// Hard Thresholding
		for(int s=0;s<_Nbr_Scale-1;s++)
		{
			float nlvl = lvl * float(DataW->norm_band(s));
			Ifloat *A = tab_f_out+s;
			for (int i=0; i < A->n_elem(); i++)
				if(abs((*A)(i)) < nlvl)
					(*A)(i) = 0.F;
		}
	else if(filtering_type==2) // soft thresholding
		for(int s=0;s<_Nbr_Scale-1;s++)
		{
			float nlvl = lvl * float(DataW->norm_band(s));
			Ifloat *A = tab_f_out+s;
			for (int i=0; i < A->n_elem(); i++)
				(*A)(i) = ::soft_threshold((*A)(i),nlvl);
		}
	else cerr<<"Filtering method not implemented yet"<<endl;

// Backward Transform
	Recons.alloc(nl,nc);
	DataW->ATROUS_2D_WT::recons(tab_f_out, Recons, _Nbr_Scale);

	delete [] tab_f_out;
	delete DataW;
	return 1;
}




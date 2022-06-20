
#include "FCur_float.h"

void FCUR_float_params::reset_params()
{
	Nl=0;
	Nc=0;
	NbrScale=3;
	NbrDir=16;
//	TabNorm=NULL;
}

int FCUR_float::size_transform()
{
	int cnt=0;
	for (int s=0; s < nbr_scale(); s++)
	for (int b=0; b < nbr_band(s); b++) 
		cnt += size_band_nl(s,b)*size_band_nc(s,b);
	
	return cnt;
}

void FCUR_float::alloc(fltarray& in, int _NbrScale, int _NbrDir)
{
	alloc_from_coarse(_NbrScale, in.ny(), in.nx(), _NbrDir, False, False, True);
}

void FCUR_float::transform(fltarray& in, float* &out, bool alloc)
{
	if(alloc) out = new float[size_transform()];
	Ifloat inFrame;
	inFrame.alloc(in.buffer(),in.ny(),in.nx());
	cur_trans(inFrame);
	export_trans(out);
}

void FCUR_float::recons(float* in, fltarray& out)
{
	import_trans(in);
	Ifloat outFrame;
	outFrame.alloc(out.buffer(),out.ny(),out.nx());
	cur_recons(outFrame);
}


void FCUR_float::mad_calculation(float* in)
{
    float *PtrBand = in;
    if (TabMADNorm.nx() != nbr_tot_band()) TabMADNorm.alloc(nbr_tot_band());
    int ind=0;
    
    for (int s=0; s < nbr_scale(); s++)
    for (int b=0; b < nbr_band(s); b++) 
    {
        TabMADNorm(ind++) = get_sigma_mad(PtrBand, size_band_nl(s,b)*size_band_nc(s,b));
        PtrBand += size_band_nl(s,b)*size_band_nc(s,b);
    }
}


void FCUR_float::mad_normalize(float* in)
{
    float *PtrBand = in;
    if (TabMADNorm.nx() != nbr_tot_band()) TabMADNorm.alloc(nbr_tot_band());
    int ind=0;
    
    for (int s=0; s < nbr_scale(); s++)
    for (int b=0; b < nbr_band(s); b++) 
    {
        float Mad = TabMADNorm(ind++);
        for (int c=0; c < size_band_nl(s,b)*size_band_nc(s,b); c++) PtrBand[c] /= Mad;
        PtrBand += size_band_nl(s,b)*size_band_nc(s,b);
    }
}

void FCUR_float::mad_unnormalize(float* in)
{
    float *PtrBand = in;
    if (TabMADNorm.nx() != nbr_tot_band()) TabMADNorm.alloc(nbr_tot_band());
    int ind=0;
    
    for (int s=0; s < nbr_scale(); s++)
        for (int b=0; b < nbr_band(s); b++) 
        {
            float Mad = TabMADNorm(ind++);
            for (int c=0; c < size_band_nl(s,b)*size_band_nc(s,b); c++) PtrBand[c] *= Mad;
            PtrBand += size_band_nl(s,b)*size_band_nc(s,b);
        }
}

void FCUR_float::soft_threshold(float* in, float lvl, bool threshold_coarse)
{
	import_trans(in);
	int i,j;
	int th = int(!threshold_coarse);
	for (int s=0; s < nbr_scale()-th; s++)
	for (int b=0; b < nbr_band(s); b++) 
	{
		for (i=0; i < size_band_nl(s,b); i++)
		for (j=0; j < size_band_nc(s,b); j++)
			(*this)(s,b,i,j) = ::soft_threshold((*this)(s,b,i,j),lvl);
	}
	export_trans(in);
}
void FCUR_float::hard_threshold(float* in, float lvl, bool threshold_coarse)
{
	import_trans(in);
	int i,j;
	int th = int(!threshold_coarse);
	for (int s=0; s < nbr_scale()-th; s++)
	for (int b=0; b < nbr_band(s); b++) 
	{
		for (i=0; i < size_band_nl(s,b); i++)
		for (j=0; j < size_band_nc(s,b); j++)
			(*this)(s,b,i,j) = ::hard_threshold((*this)(s,b,i,j),lvl);
	}
	export_trans(in);
}

void FCUR_float::substract_coarse_scale(float* in, float* &coarse, bool alloc_coarse)
{
	import_trans(in);
	int s = nbr_scale()-1;
	int b = 0;
	if(alloc_coarse) coarse = new float[size_band_nl(s,b)*size_band_nc(s,b)];
	
	int cnt=0;
	for (int i=0; i < size_band_nl(s,b); i++)
	for (int j=0; j < size_band_nc(s,b); j++)
	{
		coarse[cnt++] = (*this)(s,b,i,j);
		(*this)(s,b,i,j) = 0.;
	}
	export_trans(in);
}
void FCUR_float::add_coarse_scale(float* in, float* coarse)
{
	import_trans(in);
	int s = nbr_scale()-1;
	int b = 0;
	int cnt=0;
	for (int i=0; i < size_band_nl(s,b); i++)
	for (int j=0; j < size_band_nc(s,b); j++)
		(*this)(s,b,i,j) += coarse[cnt++];
	export_trans(in);
}

void FCUR_float::export_trans(float* out)
{
	int i,j,cnt=0;
	for (int s=0; s < nbr_scale(); s++)
	for (int b=0; b < nbr_band(s); b++) 
	{
		for (i=0; i < size_band_nl(s,b); i++)
		for (j=0; j < size_band_nc(s,b); j++)
			out[cnt++] = (*this)(s,b,i,j);
	}
}
void FCUR_float::import_trans(float* in)
{
	int i,j,cnt=0;
	for (int s=0; s < nbr_scale(); s++)
	for (int b=0; b < nbr_band(s); b++) 
	{
		for (i=0; i < size_band_nl(s,b); i++)
		for (j=0; j < size_band_nc(s,b); j++)
			(*this)(s,b,i,j) = in[cnt++];
	}
} 

// ########################################################
// #########             MEX functions            #########
// ########################################################

void fct_clear(vector< vector< Ifloat* > > &C)
{
	for(int i=0;i<(int)C.size();i++)
	{
		for(int j=0;j<(int)C[i].size();j++)
			delete C[i][j];
		C[i].clear();
	}
	C.clear();
}	

/*
int fct_get_norm(Ifloat &Data, FCUR_float_params &P)
{
// FCUR initialisation
	FCUR_float *DataC = new FCUR_float;
	P.Nl = Data.nl();
	P.Nc = Data.nc();
	DataC->Undecimated=P.Undecimated;
	DataC->alloc_from_coarse(P.NbrScale, Data.nl(), Data.nc(), P.NbrDir, False, False, True);//Bool ExtendWT, Bool IsotropWT, Bool Real
	
// Normalizing coefficients estimation and export
	DataC->get_norm_coeff(3); 
//	DataC->export_norm_coef(P.TabNorm);
	
	delete DataC;
	return 1;
}
*/
		
int fct_transform(Ifloat &Data, vector< vector< Ifloat* > > &vTabBand, FCUR_float_params &P)
{
//fstream cx;
//cx.open ("log", fstream::out);
// FCUR initialisation
	FCUR_float *DataC = new FCUR_float;
	DataC->Undecimated = P.Undecimated;
	DataC->FCUR::alloc_from_coarse(P.NbrScale, Data.nl(), Data.nc(), P.NbrDir, False, False, True);
//	if(P.TabNorm!=NULL) DataC->import_norm_coef(P.TabNorm);
	
// Forward Transform
	DataC->cur_trans(Data);
	
// vTabBand allocation
	vTabBand.resize(DataC->nbr_scale());
	for(int s=0;s<DataC->nbr_scale();s++)
	{
		vTabBand[s].resize(DataC->nbr_band(s));
		for(int b=0;b<DataC->nbr_band(s);b++)
		{
			//cerr<<"b nonorm"<<DataC->nbr_scale()<<" "<<DataC->nbr_band(s)<<endl;
			vTabBand[s][b] = new Ifloat;
			int quad = b/max(DataC->nbr_band(s)/4,1);
			int Nl = DataC->size_band_nl(s,b);
			int Nc = DataC->size_band_nc(s,b);
			// rotation of the bands in different quadrants
			if(quad%2)
			{
				vTabBand[s][b]->alloc(Nc,Nl);
				for (int i=0; i < Nl; i++)
				for (int j=0; j < Nc; j++)
					(*vTabBand[s][b])(Nc-1-j,i) = (*DataC)(s,b,i,j);
			}
			else 
			{
				vTabBand[s][b]->alloc(Nl,Nc);
				for (int i=0; i < Nl; i++)
				for (int j=0; j < Nc; j++)
					(*vTabBand[s][b])(i,j) = (*DataC)(s,b,i,j);
			}
		}
	}
	delete DataC;
	return 1;
}

int fct_recons(vector< vector< Ifloat* > > &vTabBand, Ifloat &Data, FCUR_float_params &P)
{
// FCUR initialisation
	FCUR_float *DataC = new FCUR_float;
	DataC->Undecimated=P.Undecimated;
	DataC->alloc_from_coarse(P.NbrScale, P.Nl, P.Nc, P.NbrDir, False, False, True);//Bool ExtendWT, Bool IsotropWT, Bool Real
//	if(P.TabNorm!=NULL) DataC->import_norm_coef(P.TabNorm);

// TabCF_Band allocation
	for(int s=0;s<P.NbrScale;s++)
	{
		DataC->TabCF_Band[s] = new Icomplex_f[DataC->nbr_band(s)];
		for(int b=0;b<DataC->nbr_band(s);b++)
		{
			Ifloat * A = vTabBand[s][b];
			int quad = b/max(DataC->nbr_band(s)/4,1);
			int Nl = A->nl();
			int Nc = A->nc();
		// export to DataC->TabCF_Band
		// rotation of the bands in different quadrants
			if(quad%2)
			{
				DataC->TabCF_Band[s][b].alloc(Nc, Nl);
				for (int i=0; i < A->nl(); i++)
				for (int j=0; j < A->nc(); j++)
					(*DataC)(s,b,j,Nl-1-i) = (*vTabBand[s][b])(i,j);
			}
			else
			{
				DataC->TabCF_Band[s][b].alloc(Nl, Nc);
				for (int i=0; i < A->nl(); i++)
				for (int j=0; j < A->nc(); j++)
					(*DataC)(s,b,i,j) = (*vTabBand[s][b])(i,j);
			}
		}
	}

// Backward Transform
	DataC->cur_recons(Data);

	delete DataC;
	return 1;
}


int fct_threshold(vector< vector< Ifloat* > > &vTabBand, float lvl, int filtering_type)
{
// Thresholding
	int ns = vTabBand.size();

	if(filtering_type==1)// Hard Thresholding
		for(int s=0;s<ns-1;s++)
			for(int b=0;b<(int)vTabBand[s].size();b++)
			{
				Ifloat *A = vTabBand[s][b];
				for (int i=0; i < A->n_elem(); i++)
					if(abs((*A)(i)) < lvl)
						(*A)(i) = 0.F;
			}
	else if(filtering_type==2) // soft thresholding
		for(int s=0;s<ns-1;s++)
			for(int b=0;b<(int)vTabBand[s].size();b++)
			{
				Ifloat *A = vTabBand[s][b];
				for (int i=0; i < A->n_elem(); i++)
					(*A)(i) = ::soft_threshold((*A)(i),lvl);
			}
	else cerr<<"Filtering method not implemented yet"<<endl;
		
	return 1;
}

int fct_filter(Ifloat &Data, Ifloat &Recons, FCUR_float_params &P, float lvl, int filtering_type)
{
// FCUR initialisation
	FCUR_float *DataC = new FCUR_float;
	DataC->Undecimated=P.Undecimated;
	DataC->alloc_from_coarse(P.NbrScale, Data.nl(), Data.nc(), P.NbrDir, False, False, True);//Bool ExtendWT, Bool IsotropWT, Bool Real
//	if(P.TabNorm!=NULL) DataC->import_norm_coef(P.TabNorm);

// Filtering
	DataC->cur_trans(Data);
	
	int ns = DataC->nbr_scale();
	if(filtering_type==1)// Hard Thresholding
	{
		for(int s=0;s<ns-1;s++)
		{
			int nb = DataC->nbr_band(s);
			for(int b=0;b<nb;b++)
			{
				for (int i=0; i < DataC->size_band_nl(s,b); i++)
				for (int j=0; j < DataC->size_band_nc(s,b); j++)
					if(abs((*DataC)(s,b,i,j)) < lvl ) 
						(*DataC)(s,b,i,j) = 0.F;
			}
		}
	}
	else if(filtering_type==2) // soft thresholding
	{
		for(int s=0;s<ns-1;s++)
		{
			int nb = DataC->nbr_band(s);
			for(int b=0;b<nb;b++)
			{
				for (int i=0; i < DataC->size_band_nl(s,b); i++)
				for (int j=0; j < DataC->size_band_nc(s,b); j++)
					(*DataC)(s,b,i,j) = ::soft_threshold((*DataC)(s,b,i,j),lvl);
			}
		}
	}
	else cerr<<"Filtering method not implemented yet"<<endl;

	DataC->cur_recons(Recons);
	
	delete DataC;
	return 1;
}


/*
int fct_transform(Ifloat &Data, vector< vector< Ifloat* > > &vTabBand, int _Nbr_Scale, int _NbrDir)
{
//fstream cx;
//cx.open ("log", fstream::out);
// FCUR initialisation
	FCUR_float *DataC = new FCUR_float;
	DataC->alloc_from_coarse(_Nbr_Scale, Data.nl(), Data.nc(), _NbrDir, False, False, True);//Bool ExtendWT, Bool IsotropWT, Bool Real
	if(!DataC->isset_norm) DataC->get_norm_coeff(3);

// Forward Transform
	DataC->cur_trans(Data);
	
// vTabBand allocation
	vTabBand.resize(DataC->nbr_scale());
	for(int s=0;s<DataC->nbr_scale();s++)
	{
		vTabBand[s].resize(DataC->nbr_band(s));
		for(int b=0;b<DataC->nbr_band(s);b++)
		{
//cx<<"b nonorm"<<DataC->nbr_scale()<<" "<<DataC->nbr_band(s)<<endl;
			vTabBand[s][b] = new Ifloat;
			vTabBand[s][b]->alloc(DataC->size_band_nl(s,b),DataC->size_band_nc(s,b));

			for (int i=0; i < DataC->size_band_nl(s,b); i++)
			for (int j=0; j < DataC->size_band_nc(s,b); j++)
				(*vTabBand[s][b])(i,j) = (*DataC)(s,b,i,j);
		}
	}
	delete DataC;
	return 1;
}

*/


/*
int fct_recons(vector< vector< Ifloat* > > &vTabBand, Ifloat &Data, int Nl, int Nc, int _Nbr_Scale, int _NbrDir)
{
// FCUR initialisation
	FCUR_float *DataC = new FCUR_float;
	DataC->alloc_from_coarse(_Nbr_Scale, Nl, Nc, _NbrDir, False, False, True);//Bool ExtendWT, Bool IsotropWT, Bool Real
	if(!DataC->isset_norm()) DataC->get_norm_coeff(3);

// TabCF_Band allocation
	for(int s=0;s<_Nbr_Scale;s++)
	{
		DataC->TabCF_Band[s] = new Icomplex_f[DataC->nbr_band(s)];
		for(int b=0;b<DataC->nbr_band(s);b++)
		{
			Ifloat * A = vTabBand[s][b];
			DataC->TabCF_Band[s][b].alloc(A->nl(), A->nc());
			for (int i=0; i < A->nl(); i++)
			for (int j=0; j < A->nc(); j++)
				(*DataC)(s,b,i,j) = (*vTabBand[s][b])(i,j);
		}
	}

// Backward Transform
	DataC->cur_recons(Data);

	delete DataC;
	return 1;
}
*/

/*
int fct_filter(Ifloat &Data, Ifloat &Recons, int _Nbr_Scale, int _NbrDir, float lvl, int filtering_type)// bool Undecimated
{
// FCUR initialisation
	FCUR_float *DataC = new FCUR_float;
	DataC->alloc_from_coarse(_Nbr_Scale, Data.nl(), Data.nc(), _NbrDir, False, False, True);//Bool ExtendWT, Bool IsotropWT, Bool Real
	if(!DataC->isset_norm()) DataC->get_norm_coeff(3);

// Filtering
	DataC->cur_trans(Data);
	
	int ns = DataC->nbr_scale();
	if(filtering_type==1)// Hard Thresholding
	{
		for(int s=0;s<ns-1;s++)
		{
			int nb = DataC->nbr_band(s);
			for(int b=0;b<nb;b++)
			{
				for (int i=0; i < DataC->size_band_nl(s,b); i++)
				for (int j=0; j < DataC->size_band_nc(s,b); j++)
					if(abs((*DataC)(s,b,i,j)) < lvl ) 
						(*DataC)(s,b,i,j) = 0.F;
			}
		}
	}
	else if(filtering_type==2) // soft thresholding
	{
		for(int s=0;s<ns-1;s++)
		{
			int nb = DataC->nbr_band(s);
			for(int b=0;b<nb;b++)
			{
				for (int i=0; i < DataC->size_band_nl(s,b); i++)
				for (int j=0; j < DataC->size_band_nc(s,b); j++)
					(*DataC)(s,b,i,j) = ::soft_threshold((*DataC)(s,b,i,j),lvl);
			}
		}
	}
	else cerr<<"Filtering method not implemented yet"<<endl;

	DataC->cur_recons(Recons);
	
	delete DataC;
	return 1;
}*/


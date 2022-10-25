#include "C_OMP.h"
#include <vector>
#include <MatrixOper.h>
using namespace std;

#include "IM_IO.h"
#include "IM1D_IO.h"





C_OMP::C_OMP(dblarray &dictionary)
{

	TabDico = dictionary;
	Na = TabDico.ny();
	Npix = TabDico.nx();
	eps = 2.220446049250313e-16;
	rescale_factor.alloc(Na);
	rescale_factor.init(1.0);
}

void C_OMP::update_scaling(dblarray &rescale){
	rescale_factor = rescale;
}

C_OMP::~C_OMP()
{
}

void C_OMP::update_dictionary(dblarray &new_dictionary)
{
	TabDico = new_dictionary;
	Na = TabDico.ny();
	Npix = TabDico.nx();
}
dblarray C_OMP::omp(dblarray &input_vector,int SparsityTarget, double ErrorTarget, bool Verb)
{
	dblarray residual(Npix);
	dblarray Data; // Data exctracted from input_vector, usually equals to it unless input_vector had a single line
	vector<int> indx;
	double residual_norm2;
	dblarray proj(Na);
	int Nv = input_vector.ny(); // number of vectors
	if (Nv == 0) // creating a second line for Data filled with zeros
	{
		Nv = 1;
		Data.alloc(Npix,2);
		Data.init(0);
		for (int i=0;i<Npix;i++)
		Data(i,0) = input_vector(i);
	}
	else
	Data = input_vector;
	int new_pos;
	dblarray subDico;
	dblarray pinvsubDico;
	dblarray approx_coeff;
	dblarray approx_vector;
	dblarray coefficient_vector(Na,Nv);
	coefficient_vector.init(0);
	dblarray tsubDico;
	dblarray tpinvsubDico;
	double average_sparsity = 0;
	MatOper robocop;
	double minT;
	
	if (Npix != Data.nx())
	{
		cout << "Error:Input vector and dictionary dimensions do not match. " << endl;
		cout << "   Dico.nx = " << Npix << ", InputVector.nx = " <<Data.nx() << endl;
		exit(1);
	}
	if (Verb == true)
	cout << "Running OMP with SparsityTarget = "<< SparsityTarget << " and ErrorTarget = " << ErrorTarget << endl;
	// Iterating for every vector
	for (int v=0;v<Nv;v++)
	{
		if (Verb == true)
		{
			if (Nv>100)
			{
				if (((v+1))%(Nv/100)  == 1)
				{
					cout << "Processing sample " << 1+v << " / " << Nv << ", average sparsity " << average_sparsity << "\r" ;
					cout.flush();
				}
			}
			else
			{
				cout << "Processing sample " << 1+v << " / " << Nv << ", average sparsity " << average_sparsity << "\r" ;
				cout.flush();
			}
		}
		for (int k=0;k<Npix;k++)
		residual(k) = Data(k,v);
		residual_norm2=residual.energy();
		indx.clear();
		for (int j=0;j < SparsityTarget && residual_norm2 > (1+1e-5)*ErrorTarget;j++)
		{
			// cout << "Residual " << residual_norm2 << " j " << j << endl;
			// Compute residual projection in dictionary
			proj.init(0);
			for (int k=0;k < Na;k++)
			{
				for (int l=0;l < Npix;l++)
				proj(k) += TabDico(l,k)*residual(l);
				proj(k) = proj(k) / rescale_factor(k);
			}
			proj.maxfabs(new_pos);
			indx.push_back(new_pos); // storing new index in indx
			// Extracting sub-dictionary
			subDico.alloc(Npix,indx.size());
			for (int m=0;m<indx.size();m++)
			for (int n=0;n<Npix;n++)
			subDico(n,m) = TabDico(n,indx[m]);
			//robocop.lin_eq_svd(subDico,input_vector,approx);
			if (subDico.nx()>subDico.ny())
			minT = subDico.nx()*eps;
			else
			minT = subDico.ny()*eps;
			if (j>0)
			{

				// Transposing subDico to invert it properly
				tsubDico.alloc(subDico.ny(),subDico.nx());
				for (int k=0;k<subDico.nx();k++)
				for (int l=0;l<subDico.ny();l++)
				tsubDico(l,k) = subDico(k,l);
				// Computing subDico pseudo inverse
				robocop.inv_mat_svd(tsubDico, tpinvsubDico,minT);
				tsubDico.free();
				// Transposing back the pseudo inverse of tsubDico
				pinvsubDico.alloc(tpinvsubDico.ny(),tpinvsubDico.nx());
				for (int k=0;k<tpinvsubDico.nx();k++)
				for (int l=0;l<tpinvsubDico.ny();l++)
				pinvsubDico(l,k) = tpinvsubDico(k,l);
			}
			else robocop.inv_mat_svd(subDico, pinvsubDico,minT);
			approx_coeff.alloc(indx.size());
			// Computing coefficient approximation
			for (int k=0;k<pinvsubDico.nx();k++)
			{
				approx_coeff(k) = 0;
				for (int p=0;p<pinvsubDico.ny();p++)
				{
					approx_coeff(k) += pinvsubDico(k,p)*Data(p,v);
				}
			}
			pinvsubDico.free();
			// Computing vector approximation
			approx_vector.alloc(Data.nx());
			for (int k=0;k<subDico.nx();k++)
			{
				approx_vector(k) = 0;
				for (int p=0;p<subDico.ny();p++)
				{
					approx_vector(k) += subDico(k,p)*approx_coeff(p);
				}
			}
			subDico.free();
			// Updating residual
			for (int k=0;k<Npix;k++)
			residual(k) = Data(k,v) - approx_vector(k);
			// Updating residual energy
			residual_norm2=residual.energy();
			/*residual_norm2=0;
			for (int k=0;k<residual.nx();k++)
				residual_norm2 += residual(k) * residual(k);*/
		}
		// Filling indx indices with coefficients
		for (int k=0;k < indx.size();k++)
		coefficient_vector(indx[k],v) = approx_coeff(k);
		approx_coeff.free();
		approx_vector.free();
		// Updating average sparsity
		average_sparsity = 0;
		for (int i=0;i<coefficient_vector.nx();i++)
		for (int j=0;j<coefficient_vector.ny();j++)
		if (coefficient_vector(i,j) !=0)
		average_sparsity++;
		average_sparsity = average_sparsity / (v + 1);
	}
	if (Verb == true)
	{
		cout << "OMP processing complete" << endl;
	}
	return coefficient_vector;


}

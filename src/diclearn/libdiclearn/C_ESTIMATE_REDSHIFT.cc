#include "C_DL1D.h"
#include <vector>
#include <MatrixOper.h>
#include "IM_IO.h"
#include "C_OMP.h"
#include "C_ESTIMATE_REDSHIFT.h"
using namespace std;
#include <math.h>


C_ESTIMATE_REDSHIFT::C_ESTIMATE_REDSHIFT(dblarray &Dictionary)
{
	Na = Dictionary.ny();
	Npix = Dictionary.nx();
	Dico = Dictionary;
}
C_ESTIMATE_REDSHIFT::~C_ESTIMATE_REDSHIFT()
{
}
dblarray C_ESTIMATE_REDSHIFT::estimate_redshift(dblarray &Data,intarray &Shift, int SparsityTarget, double ErrorTarget, dblarray &rms,bool Verb)
{
	int Nshift = Shift.ny(); // number of shift to be tested
	int Nspec = Data.ny(); // number of vector to test
	if (Nspec == 0)
		Nspec = 1;
	int shift_value; // current shift value
	dblarray blue_spectrum(Npix); // spectrum truncated at the current shift value
	double blue_spectrum_mean; // truncated spectrum mean
	dblarray rescale_factor(Na);
	dblarray X; // sparse coding coefficients at the current shift value
	dblarray approx_blue_spectrum(Npix); // sparse approximation of blue_spectrum
	dblarray reconstruction_error(Nshift,Nspec); // error between vectors and their sparse approximation
	double residual; // residual on the current vector at the current shift
	dblarray score(Nshift,Nspec); // score of each vector truncated at each shift value
	dblarray optimal_shift(Nspec); // shift value estimated as optimal, for each vector processed
	int optimal_index; // position of optimal shift of current vector
	dblarray vect_score(Nshift); // list of the score value of a given vector for each shift value
	dblarray rms_trunc(Npix); // rms curve truncated at the current shift

	if (Verb == true)
	{
		cout << "Running SHIELD with SparsityTarget = "<< SparsityTarget << " and ErrorTarget = " << ErrorTarget << endl;
		cout << "Processing " << Nspec << " samples and testing " << Nshift << " shift values per spectrum" << endl;
	}


	for (int v=0;v<Nspec;v++) // looping on the spectrum index
	{
		if (Nspec>10)
		{
			{
				cout << "Processing spectrum " << 1+v << " / " << Nspec << ", search centered on " << Shift(v,(Nshift - 1)/2) << endl;
			}
		}
		else
		{
			cout << "Processing spectrum " << 1+v << " / " << Nspec << endl ;
		}

		for (int s=0;s<Nshift;s++) // Computing sparse aproximation error for each vectors at a given shift
		{
			shift_value = Shift(v,s);
			// Truncating data and RMS curve at current shift value
			for (int p=0;p<Npix;p++)
			{
				blue_spectrum(p) = Data(p+shift_value-1,v);
				rms_trunc(p) = rms(p+shift_value-1);
			}
			// Computing and removing mean of each truncated spectrum
			blue_spectrum_mean = blue_spectrum.mean();
			for (int p=0;p<Npix;p++)
				blue_spectrum(p) = blue_spectrum(p) - blue_spectrum_mean;
			// Computing rescaling factor given rms curve and dictionary
			for (int a=0;a<Na;a++)
			{
				rescale_factor(a) = 0;
				for (int p=0;p<Npix;p++)
					rescale_factor(a) += rms_trunc(p)*Dico(p,a)*Dico(p,a);
			}
			// Building sparse encoder given rescaling factors
			C_OMP coder(Dico);
			coder.update_scaling(rescale_factor);
			// Sparse coding blue shifted spectrum in the dictionary
			X = coder.omp(blue_spectrum,SparsityTarget,ErrorTarget,False);
			// Computing sparse approximation spectrum given sparse coefficient
			for (int k=0;k<Npix;k++)
			{
				approx_blue_spectrum(k) = 0;
				for (int p=0;p<Na;p++)
				{
					approx_blue_spectrum(k) += Dico(k,p)*X(p);
				}
			}
			// Computing score (residual energy normilized by rms) of current shifted spectrum
			residual = 0;
			for (int k=0;k<Npix;k++)
				//  error_temp = sqrt(sum((res.^2)./rms)) / norm(D*X);
				residual += ((blue_spectrum(k)- approx_blue_spectrum(k))*(blue_spectrum(k)- approx_blue_spectrum(k))) / rms_trunc(k);
			residual = sqrt(residual) / sqrt(approx_blue_spectrum.energy());
			score(s,v) = residual;
		}
	}
	for (int i=0;i<Nspec;i++) // finding for each vector which shift values gives the lowest score
	{
		for (int k=0;k<Nshift;k++)
			vect_score(k) = score(k,i);
		optimal_index = 0;
		for (int k=0;k<Nshift;k++)
			if (vect_score(k) < vect_score(optimal_index))
				optimal_index = k;
		optimal_shift(i) = Shift(i,optimal_index);
	}
	return optimal_shift;
}













/***********************************************************
**	Copyright (C) 1999 CEA
************************************************************
**
**    UNIT
**
**    Version:  1.0
**
**    Author: 	J.L. Starck
**
**    Date: 	27/08/99
**    
**    File:  	cf_alpha_obj.cc
**
************************************************************
**
**  1D,2D,3D correlation function analysis 
**  
************************************************************/

#include "cf_tools.h"
#include "DefPoint.h"
#include "cf.h"
#include <omp.h>

extern int Nproc_max;

/****************************************************************************/

void CorrFunAna::cf_find_pairs(ArrayPoint & Data, ArrayPoint & DataWeight, fltarray &CF_DataData)
{
	int N = Data.np();
	int i,j;
	double weight;
	
	init();
	
	int nbins=np();
	
	if (CF_DataData.nx() != nbins)
		CF_DataData.alloc(nbins);

	dblarray TempHisto(nbins);
	
	#ifdef _OPENMP
		int Nproc=omp_get_num_procs();
		if(Nproc>Nproc_max) Nproc=Nproc_max;
		printf("\n %u processors used for the loop \n\n",Nproc);
	#endif
   
   #pragma omp parallel default(shared)  shared(N) private(i,j) firstprivate(TempHisto) num_threads(Nproc)
	{
		if (Data.TCoord == 2 && Data.dim() == 2) 
		{
			#pragma omp for schedule(dynamic)
			for (i=0; i < N-1; i++)
			{
				for (j=i+1; j < N; j++)
				{ 
					float r = squaresphdist(Data(i), Data(j));
					int Ind = index_dist(r);
					
					if(Ind >=0)
					{
						weight=DataWeight(i).x()*DataWeight(j).x();
						TempHisto(Ind) += weight;
					}
				}       
			}
		}
		else
		{
			#pragma omp for schedule(dynamic)
			for (i=0; i < N-1; i++)
			{
				for (j=i+1; j < N; j++)
				{ 
					float r = squaredist(Data(i), Data(j),SquareDistMax);
					int Ind = index_dist(r);
					if(Ind >=0)
					{
						weight=DataWeight(i).x()*DataWeight(j).x();
						TempHisto(Ind) += weight;
					}
				}
			}
		}		
		#pragma omp critical
		{
			for(i=0; i < nbins; i++) 
			{
				CF_DataData(i)+=TempHisto(i); //Join TempHisto
			}
		}
	}

}


/****************************************************************************/

void CorrFunAna::cf_find_pairs(ArrayPoint & Data1, ArrayPoint & Data2, ArrayPoint & DataWeight1, ArrayPoint & DataWeight2, fltarray &CF_Data1Data2)
{
	int N1 = Data1.np();
	int N2 = Data2.np();
	int i,j;
	double weight;
	
	init();
	
   	int nbins=np();
	if (CF_Data1Data2.nx() != nbins)
		CF_Data1Data2.alloc(nbins);
	dblarray TempHisto(nbins);
  

	#ifdef _OPENMP
		int Nproc=omp_get_num_procs();
		if(Nproc>Nproc_max) Nproc=Nproc_max;
		printf("\n %u processors used for the loop \n\n",Nproc);
	#endif
	
	#pragma omp parallel default(shared)  shared(N1,N2) private(i,j) firstprivate(TempHisto) num_threads(Nproc)
	{
		if (Data1.TCoord == 2 && Data1.dim() == 2) 
		{
			#pragma omp for schedule(dynamic)
			for (i=0; i < N1; i++)
			{
				for (j=0; j < N2; j++)
				{ 
					float r = squaresphdist(Data1(i), Data2(j));
					int Ind = index_dist(r);
					if(Ind >=0)
					{
						weight=DataWeight1(i).x()*DataWeight2(j).x();
						TempHisto(Ind) += weight;
					}
				}       
			}
		}
		else
		{
			#pragma omp for schedule(dynamic)
			for (i=0; i < N1; i++)
			{
				for (j=0; j < N2; j++)
				{ 
					float r = squaredist(Data1(i), Data2(j),SquareDistMax);
					int Ind = index_dist(r);
					if(Ind >= 0)
					{
						weight=DataWeight1(i).x()*DataWeight2(j).x();	
						TempHisto(Ind) += weight;
					}
				}
			}
		}
		#pragma omp critical
		{
			for(i=0; i < nbins; i++) 
			{
				CF_Data1Data2(i)+=TempHisto(i); //Join TempHisto
			}
		}
	}

}
	       


/****************************************************************************/

int CorrFunAna::index_dist(float r)
{
   int Ind=-1;	
   if ((r > SquareDistMin) && (r < SquareDistMax))
   {
	  r=sqrt(r);
      Ind = (int) ((r - DistMin)/Step);
	  if ((Ind < 0) || (Ind >= Nc))
	  {
	    //cerr << "Error: bad distance index = " << Ind << " Dist = " << r << endl;
		Ind=-1;
      }
   } 
   return Ind;  
}

/****************************************************************************/

CorrFunAna::CorrFunAna(float Dmin, float Dmax, int Nbin)
{
   int i;
   DistMin = Dmin;
   DistMax = Dmax;
	
   Step = (DistMax - DistMin) / (float) Nbin;
   Nc = Nbin;
   PairHisto.alloc(2,Nc);
      
   for (i=0; i < Nc; i++)
     PairHisto(0,i) = DistMin+ i*Step + Step/2.;

	SquareDistMin = DistMin*DistMin;
	SquareDistMax = DistMax*DistMax;
}

/****************************************************************************/

CorrFunAna::CorrFunAna(float Dmin, float Dmax, float StepVal)
{
   int i;
   DistMin = Dmin;
   DistMax = Dmax;

   Step = StepVal;

   if((DistMax - DistMin)/Step == int((DistMax - DistMin) / Step) ) 
	   Nc = int((DistMax - DistMin) / Step);
   else
	   Nc = int((DistMax - DistMin) / Step) + 1;
   DistMax=DistMin+Nc*Step;
   PairHisto.alloc(2,Nc);
      
   DistMax = DistMin+ Nc*Step;
   // Must be the effective Distmax which is not neccessarily the case if we fix
   // Dmin, Dmax and StepVal. Because when we add pair we reject the pairs which
   // are not between Distmin and Distmax in index_dist

   for (i=0; i < Nc; i++)
     PairHisto(0,i) = DistMin+ i*Step + Step/2.;

	SquareDistMin = DistMin*DistMin;
	SquareDistMax = DistMax*DistMax;
}

/****************************************************************************/
/*

void make_histo(fltarray &dd, fltarray &rr, fltarray &dr, fltarray &Result)
{
	int nbins=dd.nx();
	int nalpha=dd.ny();
	
	
	for (int i=0; i< nbins;i++) 
	{
		Result(i,1) = dd(i);
		Result(i,2) = rr(i);
		Result(i,3) = dr(i);
	}
		
}*/

/****************************************************************************/


void normalize_histo(ArrayPoint & DataWeight, ArrayPoint & RndWeight, fltarray &Result)
{

	
	int Np = DataWeight.np();
	int NpRnd = RndWeight.np();
	int i;
	
	int nbins=Result.nx();
	
	//\sum_wi for data, rnd
	double TotalWeightData1=0.0;
	double TotalWeightRnd1=0.0;
	
	//\sum_wi^2 for data, rnd
	double TotalWeightData2=0.0;
	double TotalWeightRnd2=0.0;
	
	double TotalWeightDD;
	double TotalWeightRR;
	double TotalWeightDR;
	
		

	for(i=0;i<Np;i++)
	{
			TotalWeightData1+=double(DataWeight(i).x());
			TotalWeightData2+=double(DataWeight(i).x())*double(DataWeight(i).x());
			
	}
	for(i=0;i<NpRnd;i++)
	{
			TotalWeightRnd1+=double(RndWeight(i).x());
			TotalWeightRnd2+=double(RndWeight(i).x())*double(RndWeight(i).x());
	}
		
	TotalWeightDD=0.5*(TotalWeightData1*TotalWeightData1-TotalWeightData2);
	TotalWeightRR=0.5*(TotalWeightRnd1*TotalWeightRnd1-TotalWeightRnd2);
	TotalWeightDR=TotalWeightData1*TotalWeightRnd1;		

	
	for(i=0;i<nbins;i++)
	{
		Result(i,1) /= TotalWeightDD ;
		Result(i,2) /= TotalWeightRR ;
		Result(i,3) /= TotalWeightDR ;
	}
	
}





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
#include "cf_alpha.h"
#include <omp.h>

extern int Nproc_max;

/****************************************************************************/

void CorrFunAnaAlpha::cf_find_pairs(ArrayPoint & Data, ArrayPoint & DataWeight, ArrayPoint & DataAlpha, fltarray &CF_DataData)
{
	int N = Data.np();
	int i,j;
	double weight;
	
	init();
	
	int nbins=np();
	
	if (CF_DataData.nx() != nbins)
		CF_DataData.alloc(nbins);
	int nalpha=CF_DataData.ny();
	int minAlpha[N]; int maxAlpha[N];

	for(i=0;i<N;i++)
	{
		minAlpha[i]=round(DataAlpha(i).x());
		maxAlpha[i]=round(DataAlpha(i).y());
	}

	dblarray TempHisto(nbins,nalpha);
	
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
						int IndAlphaMin=max(minAlpha[i],minAlpha[j]);
						int IndAlphaMax=min(maxAlpha[i],maxAlpha[j]);
					
						if (IndAlphaMin<IndAlphaMax && IndAlphaMin<nalpha) 
						{
							weight=DataWeight(i).x()*DataWeight(j).x();
							TempHisto(Ind,IndAlphaMin) += weight;
							if(IndAlphaMax<nalpha) 
								TempHisto(Ind,IndAlphaMax) -= weight;
						}
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
					//float r = squaredist(Data(i), Data(j));
					int Ind = index_dist(r);
					if(Ind >=0)
					{
						int IndAlphaMin=max(minAlpha[i],minAlpha[j]);
						int IndAlphaMax=min(maxAlpha[i],maxAlpha[j]);
					
						if (IndAlphaMin<IndAlphaMax && IndAlphaMin<nalpha) 
						{
							weight=DataWeight(i).x()*DataWeight(j).x();
							TempHisto(Ind,IndAlphaMin) += weight;
							if(IndAlphaMax<nalpha) 
								TempHisto(Ind,IndAlphaMax) -= weight;
						}
					}
				}
			}
		}		
		#pragma omp critical
		{
			for(i=0; i < nbins; i++) 
			{
				for(j=0; j < nalpha; j++) 
					CF_DataData(i,j)+=TempHisto(i,j); //Join TempHisto
			}
		}
	}

}


/****************************************************************************/

void CorrFunAnaAlpha::cf_find_pairs(ArrayPoint & Data1, ArrayPoint & Data2, ArrayPoint & DataWeight1, ArrayPoint & DataWeight2,
							   ArrayPoint & DataAlpha1, ArrayPoint & DataAlpha2, fltarray &CF_Data1Data2)
{
	int N1 = Data1.np();
	int N2 = Data2.np();
	int i,j;
	double weight;
	
	init();
	
   	int nbins=np();
	int nalpha=CF_Data1Data2.ny();
	if (CF_Data1Data2.nx() != nbins)
		CF_Data1Data2.alloc(nbins,nalpha);
	dblarray TempHisto(nbins,nalpha);
  
	int minAlpha1[N1]; int maxAlpha1[N1];
	int minAlpha2[N2]; int maxAlpha2[N2];
	for(i=0;i<N1;i++) {minAlpha1[i]=round(DataAlpha1(i).x()); maxAlpha1[i]=round(DataAlpha1(i).y());}
	for(i=0;i<N2;i++) {minAlpha2[i]=round(DataAlpha2(i).x()); maxAlpha2[i]=round(DataAlpha2(i).y());}
	

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
						int IndAlphaMin=max(minAlpha1[i],minAlpha2[j]);
						int IndAlphaMax=min(maxAlpha1[i],maxAlpha2[j]);
					
						if (IndAlphaMin<IndAlphaMax && IndAlphaMin<nalpha) 
						{
							weight=DataWeight1(i).x()*DataWeight2(j).x();
							TempHisto(Ind,IndAlphaMin) += weight;
							if(IndAlphaMax<nalpha) 
								TempHisto(Ind,IndAlphaMax) -= weight;
						}
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
					//float r = squaredist(Data1(i), Data2(j));
					int Ind = index_dist(r);
					if(Ind >= 0)
					{
						int IndAlphaMin=max(minAlpha1[i],minAlpha2[j]);
						int IndAlphaMax=min(maxAlpha1[i],maxAlpha2[j]);
					
						if (IndAlphaMin<IndAlphaMax && IndAlphaMin<nalpha) 
						{
							weight=DataWeight1(i).x()*DataWeight2(j).x();	
							TempHisto(Ind,IndAlphaMin) += weight;
							if(IndAlphaMax<nalpha) 
								TempHisto(Ind,IndAlphaMax) -= weight;
						}
					}
				}
			}
		}
		#pragma omp critical
		{
			for(i=0; i < nbins; i++) 
			{
				for(j=0; j < nalpha; j++) 
					CF_Data1Data2(i,j)+=TempHisto(i,j); //Join TempHisto
			}
		}
	}

}
	       


/****************************************************************************/

int CorrFunAnaAlpha::index_dist(float r)
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

CorrFunAnaAlpha::CorrFunAnaAlpha(float Dmin, float Dmax, int Nbin)
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

CorrFunAnaAlpha::CorrFunAnaAlpha(float Dmin, float Dmax, float StepVal)
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


void make_histo(fltarray &dd, fltarray &rr, fltarray &dr, fltarray &Result)
{
	int nbins=dd.nx();
	int nalpha=dd.ny();
	
	
	for (int i=0; i< nbins;i++) 
	{
		Result(i,0,1) = dd(i,0);
		Result(i,0,2) = rr(i,0);
		Result(i,0,3) = dr(i,0);
		for (int j=1; j< nalpha;j++) 
		{
			Result(i,j,1) = Result(i,j-1,1)+dd(i,j);
			Result(i,j,2) = Result(i,j-1,2)+rr(i,j);
			Result(i,j,3) = Result(i,j-1,3)+dr(i,j);
		}
	}
		
}

/****************************************************************************/


void normalize_histo(ArrayPoint & DataWeight, ArrayPoint & RndWeight,
					 ArrayPoint & DataAlpha,  ArrayPoint & RndAlpha,  fltarray &Result)
{
	int Np = DataWeight.np();
	int NpRnd = RndWeight.np();
	int i,j;
	
	int nbins=Result.nx();
	int nalpha=Result.ny();
	
	//\sum_wi for data, rnd
	dblarray TotalWeightData1(nalpha);
	dblarray TotalWeightRnd1(nalpha);
	
	//\sum_wi^2 for data, rnd
	dblarray TotalWeightData2(nalpha);
	dblarray TotalWeightRnd2(nalpha);
	
	dblarray TotalWeightDD(nalpha);
	dblarray TotalWeightRR(nalpha);
	dblarray TotalWeightDR(nalpha);
	
		
	for(j=0;j<nalpha;j++)
	{
		for(i=0;i<Np;i++)
		{
			if( (round(DataAlpha(i).x()) <=j) && (round(DataAlpha(i).y())>j) )
			{
				TotalWeightData1(j)+=double(DataWeight(i).x());
				TotalWeightData2(j)+=double(DataWeight(i).x())*double(DataWeight(i).x());
			}
		}
		for(i=0;i<NpRnd;i++)
		{
			if( (round(RndAlpha(i).x()) <=j) && (round(RndAlpha(i).y())>j) )
			{
				TotalWeightRnd1(j)+=double(RndWeight(i).x());
				TotalWeightRnd2(j)+=double(RndWeight(i).x())*double(RndWeight(i).x());
			}
		}
		
		TotalWeightDD(j)=0.5*(TotalWeightData1(j)*TotalWeightData1(j)-TotalWeightData2(j));
		TotalWeightRR(j)=0.5*(TotalWeightRnd1(j)*TotalWeightRnd1(j)-TotalWeightRnd2(j));
		TotalWeightDR(j)=TotalWeightData1(j)*TotalWeightRnd1(j);		
	}
		
	for(j=0;j<nalpha;j++)
	{
		for(i=0;i<nbins;i++)
		{
			Result(i,j,1) /= TotalWeightDD(j) ;
			Result(i,j,2) /= TotalWeightRR(j) ;
			Result(i,j,3) /= TotalWeightDR(j) ;
		}
	}
	
}





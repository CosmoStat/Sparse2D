/*******************************************************************************
**
**    DESCRIPTION  tools for filtering using the A trous algorithm
**    -----------  
**                 
******************************************************************************/
 
#include "Atrou3DFil.h"

 
//****************************************************************************** 

void ATROUS_3D_FIL::Filtering_Structure(intarray *& Support, int LastScale,Bool InfoBand,
	Bool RemoveBorder, int FistScaleDetect)
{
   //remove the structure on the border of the cube 
   if((InfoBand == True)|| (RemoveBorder == True))
   {
       int Nx = (Support[0]).nx();
       int Ny = (Support[0]).ny();
       int Nz = (Support[0]).nz();
       intarray seg(Nx, Ny,Nz);
       for(int s=FistScaleDetect; s<LastScale-1;s++)
       {
          int nb_seg;
          cube_segment(Support[s], seg ,nb_seg, 0., RemoveBorder, 1);
          if (RemoveBorder ==True)
	  {
             for (int i=0; i<Nx;i++)
             for (int j=0; j<Ny;j++)
             for (int k=0; k<Nz;k++)
		 if(seg(i,j,k)==0) (Support[s])(i,j,k)=0;
	  }
	  if (InfoBand == True) cout << "Scale nbr "<< s+1 <<":    "<< nb_seg <<
			" different structures are detected"<<endl;
      }
   }
}

 
/****************************************************************************/

void ATROUS_3D_FIL::threshold(intarray *& Support, int NbrScale)
{
  //when the support is computed, it thresholds the differents scales.
  //it doesn't threshold the last scale (the smooth scale)
	
   for (int s = 0; s <  NbrScale-1; s++)
   for (int i=0;i<Support[s].nx(); i++)
   for (int j=0;j<Support[s].ny(); j++)  
   for (int k=0;k<Support[s].nz(); k++)
   	if((Support[s])(i,j,k)!=1)  (Wavelet_Coef[s])(i,j,k) =0.;
   
   
}


/****************************************************************************/

void ATROUS_3D_FIL::No_Single_Point(int Nx,int Ny, int Nz,
	intarray *& Support,FewEventPoisson FEP, int Nbr_Scale, Bool Verbose)
	//remove the single points of the cubes of the support
{
 
    Bool zero;   
    int Cpt = 0; 
    for(int s=0;s<Nbr_Scale;s++)
    {
    for(int i=0;i<Nx;i++)
    for(int j=0;j<Ny;j++)
    for(int k=0;k<Nz;k++)
    {
    	zero=True;
	int cur_i=-1;
	int cur_j=-1;
	int cur_k=-1;
    	if ((Support[s])(i,j,k) !=0)
	{
	   while ((zero==True) && (cur_i<2))
	   {
	   	while ((zero==True) && (cur_j<2))
		{
	   	     while ((zero==True) && (cur_k<2) && 
		            !((cur_i==0) && (cur_j==0) && (cur_k==0)))
		     {
	   		  if(FEP.get_pix(i+cur_i,j+cur_j,k+cur_k,Support[s],I_ZERO)!=0)
		    		zero=False;
			  cur_k++;
	   	     }
	   	     cur_k=-1;
	             cur_j++;
	   	}
	        cur_j=-1;
	        cur_i++;
	   }
	}
	if ((zero==True) && ((Support[s])(i,j,k) !=0))
	{ 
	    Cpt ++;
	    (Support[s])(i,j,k)=0;
	}
     }}
     if (Verbose == True)
           cout << "Nbr isolated pixels = " << Cpt << endl;
}
 
/****************************************************************************/
 
//	A trous transformation 
void  ATROUS_3D_FIL::wttransform(fltarray &Cube, int NbrScale,
	Bool WriteRes,Bool InfoBand,char Name_Cube_Out[])
{
   char Name_Cube_Out2 [50];
   Wavelet_Coef[0] = Cube;
   
   transform(Cube, Wavelet_Coef, NbrScale);
   for (int s=0;s< NbrScale;s++)
   {
    	if (WriteRes == True)
	{
		sprintf(Name_Cube_Out2,"%s%d.fits",Name_Cube_Out,s);
		fits_write_fltarr(Name_Cube_Out2, Wavelet_Coef[s]);
	}
	sprintf(Name_Cube_Out2,"Scale %d :",s);
	//if(InfoBand == True)  (Wavelet_Coef[s]).info(Name_Cube_Out2);
	if(InfoBand == True)  (Wavelet_Coef[s]).info();
   }
}

/****************************************************************************/

void ATROUS_3D_FIL::wtrecons(fltarray & Data_Out, int NbrScale, Bool InfoBand, 
		Bool AddLastScale, Bool adjoint)
{
	Adjoint=adjoint;
   recons(Wavelet_Coef,  Data_Out, NbrScale, AddLastScale);
   if(InfoBand== True) cout <<"Sigma:	"<< Data_Out.sigma()<<endl;
}

/****************************************************************************/

void free(fltarray * TabBand, int Nbr_Plan)
{
    if (Nbr_Plan != 0) delete [] TabBand;
}

/****************************************************************************/

void free(intarray * TabBand, int Nbr_Plan)
{
    if (Nbr_Plan != 0) delete [] TabBand;
}

/****************************************************************************/

void alloc_array (fltarray * & TabBand, int Nx, int Ny, int Nz, int NbrBand)
{
    TabBand = new fltarray [NbrBand];
    for (int s = 0; s < NbrBand; s++)
        TabBand[s].alloc (Nx, Ny, Nz);   
}

/****************************************************************************/

void alloc_array (intarray * & Nb_Event, int Nx, int Ny, int Nz, int NbrBand)
{
    Nb_Event = new intarray [NbrBand];
    for (int s = 0; s < NbrBand; s++)
        Nb_Event[s].alloc (Nx, Ny, Nz);   
}

/****************************************************************************/

void init (fltarray * & TabBand, int NbrBand, float val)
{
    for (int s = 0; s < NbrBand; s++)
        TabBand[s].init(val);   
}

/****************************************************************************/

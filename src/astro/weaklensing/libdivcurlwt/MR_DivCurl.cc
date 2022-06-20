/******************************************************************************
**                   Copyright (C) 2009 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Erwan Deriaz & Jean-Luc Starck
**
**    Date:  08/09/08
**    
**    File:  MR_DivCurl.cc
**
*******************************************************************************
**
**    DESCRIPTION  isotropic and anisotropic div and curl decomposition class implementation  
**    ----------- 
**                 
******************************************************************************/

#include "MR_DivCurl.h"

#include "helmhb_neu.c"

/*********************************************************************/





/*********************************************************************/
// N must be a power of 2
void HelmCatalog::allocn(fltarray & TabXYG, int N)
{
 	CatNp = TabXYG.nx();
    TabX.alloc(CatNp);
    TabY.alloc(CatNp);
    TabSigmaX.alloc(CatNp);
    TabSigmaY.alloc(CatNp);
    DataG.alloc(CatNp, 2);
	int i=0;
	MinX = TabXYG(i,0);
	MinY = TabXYG(i,1);
	MaxX = TabXYG(i,0);
	MaxY = TabXYG(i,1);
	for (i=0; i < CatNp; i++)
	{
        TabX(i) = TabXYG(i,0);
        TabY(i) = TabXYG(i,1);
        DataG(i,0) = TabXYG(i,2);
        DataG(i,1) = TabXYG(i,3);
        TabSigmaX(i) = TabXYG(i,4);
        TabSigmaY(i) = TabXYG(i,5);
        
		if (MinX > TabXYG(i,0)) MinX = TabXYG(i,0);
		if (MinY > TabXYG(i,1)) MinY = TabXYG(i,1);
		if (MaxX < TabXYG(i,0)) MaxX = TabXYG(i,0);
		if (MaxY < TabXYG(i,1)) MaxY = TabXYG(i,1);
	}
    
    float xDyn = MaxX - MinX;
    float yDyn = MaxY - MinY;
    CatNx = N;
    CatNy = N;
    BinCatX = (float) xDyn / (float) (N);
    BinCatY = (float) yDyn / (float) (N);
    BinCat = BinCatX;
    if (BinCatY > BinCat) BinCat = BinCatY;
    
    // if (Verbose == True)
    {
        cout << "CAT: Minx = " << MinX <<  ", Maxx = " << MaxX << ", Miny = " << MinY <<  ", MaxY = " << MaxY << "=> Nx = " << CatNx << ", Ny = " << CatNy << endl;
    }
    TabNp.alloc(CatNx,CatNy);
    for(int i=0; i<CatNp; i++) TabNp (posx(i),posy(i)) ++; 
    DataG.info("G");   
}

/*********************************************************************/

void HelmCatalog::alloc(fltarray & TabXYG, float BCat)
{
    BinCat = BinCatX = BinCatY = BCat;
	CatNp = TabXYG.nx();
    TabX.alloc(CatNp);
    TabY.alloc(CatNp);
    TabSigmaX.alloc(CatNp);
    TabSigmaY.alloc(CatNp);
    DataG.alloc(CatNp, 2);
	int i=0;
	MinX = TabXYG(i,0);
	MinY = TabXYG(i,1);
	MaxX = TabXYG(i,0);
	MaxY = TabXYG(i,1);
	for (i=0; i < CatNp; i++)
	{
        TabX(i) = TabXYG(i,0);
        TabY(i) = TabXYG(i,1);
        DataG(i,0) = TabXYG(i,2);
        DataG(i,1) = TabXYG(i,3);
        TabSigmaX(i) = TabXYG(i,4);
        TabSigmaY(i) = TabXYG(i,5);
        
		if (MinX > TabXYG(i,0)) MinX = TabXYG(i,0);
		if (MinY > TabXYG(i,1)) MinY = TabXYG(i,1);
		if (MaxX < TabXYG(i,0)) MaxX = TabXYG(i,0);
		if (MaxY < TabXYG(i,1)) MaxY = TabXYG(i,1);
	}

    float xDyn = MaxX - MinX;
    float yDyn = MaxY - MinY;
    CatNx = int( xDyn / BinCat ) + 1;
    CatNy = int( yDyn / BinCat ) + 1;
    
    
   // if (Verbose == True)
    {
        cout << "CAT: Minx = " << MinX <<  ", Maxx = " << MaxX << ", Miny = " << MinY <<  ", MaxY = " << MaxY << "=> Nx = " << CatNx << ", Ny = " << CatNy << endl;
    }
    TabNp.alloc(CatNx,CatNy);
    for(int i=0; i<CatNp; i++) TabNp (posx(i),posy(i)) ++;
}

/*********************************************************************/

void HelmCatalog::cat2ima(fltarray & Ima)
{
    cat2ima(DataG, Ima);
}

/*********************************************************************/

void HelmCatalog::cat2ima(fltarray & Shear_or_Pol1D, fltarray & Ima)
{
   Ima.resize( nx(), ny(), 2 );
   Ima.init();
   for(int i=0; i< np() ; i++)
   { 
       Ima (posx(i),posy(i), 0) += Shear_or_Pol1D(i,0);
       Ima (posx(i),posy(i), 1) += Shear_or_Pol1D(i,1);
   }

   for(int j=0; j < nx() ; j++)
   for(int i=0; i< ny() ; i++)
   { 
       if (np_in_xy(j,i) > 0)
       {
          Ima (j,i, 0) /= (float) np_in_xy(j,i);
          Ima (j,i, 1) /= (float) np_in_xy(j,i);
       }
   }
}

/*********************************************************************/

void HelmCatalog::cat2imasigma(fltarray & ImaSigmaXY)
{
    ImaSigmaXY.resize( nx(), ny(), 2 );
    ImaSigmaXY.init();
    for(int i=0; i< np() ; i++)
    { 
        ImaSigmaXY (posx(i),posy(i), 0) += sigmax(i)*sigmax(i);
        ImaSigmaXY (posx(i),posy(i), 1) += sigmay(i)*sigmay(i);
    }
    
    for(int j=0; j < nx() ; j++)
    for(int i=0; i< ny() ; i++)
    { 
        if (np_in_xy(j,i) > 0)
        {
            ImaSigmaXY (j,i, 0) = sqrt( ImaSigmaXY (j,i, 0) / (float) np_in_xy(j,i));
            ImaSigmaXY (j,i, 1) = sqrt( ImaSigmaXY (j,i, 1) / (float) np_in_xy(j,i));
        }
    }
}

/*********************************************************************/

void div_quasi_interpol(fltarray & Tabp, fltarray & Tabc, Bool BorderWavelet)
{
  int Nx = Tabp.nx();
  int Ny = Tabp.ny();
  int N=Nx-1;
  if (BorderWavelet == False) 
  {
     Tabc.resize(Nx,Ny,2);
      int i1,i2;
      float f[4]={-1./8.,5./8.,5./8.,-1./8.};
      
      for (i2=0;i2<Ny;i2++) 
          for (i1=0;i1<Nx;i1++) 
              Tabc(i1,i2,0) = f[0] * Tabp(i1-1,i2,0,DEF_DIVROT_BORD) + f[1] * Tabp(i1,i2,0) + f[2] * Tabp(i1+1,i2,0,DEF_DIVROT_BORD) + f[3] * Tabp(i1+2,i2,0,DEF_DIVROT_BORD);
      
      for (i2=0;i2<Ny;i2++) 
          for (i1=0;i1<Nx;i1++)
              Tabc(i1,i2,1) = f[0] * Tabp(i1,i2-1,1,DEF_DIVROT_BORD) + f[1] * Tabp(i1,i2,1) + f[2] * Tabp(i1,i2+1,1,DEF_DIVROT_BORD) + f[3] * Tabp(i1,i2+2,1,DEF_DIVROT_BORD);
  }
  else  
  {
     //  cout << " div_quasi_interpol Border " << endl;
     Tabc.resize(Nx+1,Ny+1,2);
     float *u1 = Tabp.buffer();
     float *u2 = Tabp.buffer() + Nx*Ny;
     float *v1 = Tabc.buffer();
     float *v2 = Tabc.buffer() + (Nx+1)*(Ny+1);
     qi10(u1,v1, N);
     qi01(u2,v2, N); 
  }
}

/*********************************************************************/

void curl_quasi_interpol(fltarray & Tabp, fltarray & Tabc, Bool BorderWavelet)
{
  int Nx = Tabp.nx();
  int Ny = Tabp.ny();
  int N=Nx-1;
  
  if (BorderWavelet == False) 
  {
     Tabc.resize(Nx,Ny,2);
    
     int i1,i2;
     float f[4]={-1./8.,5./8.,5./8.,-1./8.};
   
     for (i2=0;i2<Ny;i2++) 
     for (i1=0;i1<Nx;i1++) 
        Tabc(i1,i2,0)  = f[0] * Tabp(i1,i2-1,0, DEF_DIVROT_BORD) + f[1] * Tabp(i1,i2,0) 
                   + f[2] *  Tabp(i1, i2+1, 0, DEF_DIVROT_BORD) + f[3] * Tabp(i1, i2+2,0, DEF_DIVROT_BORD);
   
     for (i2=0;i2<Ny;i2++) 
     for (i1=0;i1<Nx;i1++)
        Tabc(i1,i2,1) = f[0] * Tabp(i1-1,i2,1,DEF_DIVROT_BORD) + f[1] * Tabp(i1,i2,1) + f[2] * Tabp(i1+1,i2,1,DEF_DIVROT_BORD) + f[3] * Tabp(i1+2,i2, 1,DEF_DIVROT_BORD);
  }
  else  
  {
      // cout << " curl_quasi_interpol Border " << endl;
      Tabc.resize(Nx+1,Ny+1,2);
      float *u1 = Tabp.buffer();
      float *u2 = Tabp.buffer() + Nx*Ny;
      float *v1 = Tabc.buffer();
      float *v2 = Tabc.buffer() + (Nx+1)*(Ny+1);
      qi01(u1,v1, N);
      qi10(u2,v2, N); 
  }    
}
 
/*********************************************************************/

void div_point_value(fltarray & Tabc, fltarray & Tabp, Bool BorderWavelet)
{
  int Nx = Tabc.nx();
  int Ny = Tabc.ny();
   
    if (BorderWavelet == False) 
    {
      Tabp.resize(Nx,Ny,2);
      int i1,i2;
      float f[2]={1./2.,1./2.};
  
      for (i2=0;i2<Ny;i2++) 
      for (i1=0;i1<Nx;i1++) 
      Tabp(i1,i2,0) = f[0] *  Tabc(i1-1,i2,0,DEF_DIVROT_BORD) + f[1] * Tabc(i1,i2,0);
  
      for (i2=0;i2<Ny;i2++) 
      for (i1=0;i1<Nx;i1++)
      Tabp(i1,i2,1) = f[0]  * Tabc(i1,i2-1, 1,DEF_DIVROT_BORD) + f[1] * Tabc(i1,i2,1);
    }
    else  
    {
       //  cout << " div_point_value Border " << endl;
        int N=Nx-2;
        Tabp.resize(Nx-1,Ny-1,2);
        float *w1 = Tabp.buffer();
        float *w2 = Tabp.buffer() + (Nx-1)*(Ny-1);
        float *v1 = Tabc.buffer();
        float *v2 = Tabc.buffer() + Nx*Ny;
        rec10(w1,v1,N);
        rec01(w2,v2,N);
    }        
}

/*********************************************************************/

void curl_point_value(fltarray & Tabc, fltarray & Tabp, Bool BorderWavelet)
{
  int Nx = Tabc.nx();
  int Ny = Tabc.ny();
  
    if (BorderWavelet == False) 
    {
      Tabp.resize(Nx,Ny,2);
    
       int i1,i2;
      float f[2]={1./2.,1./2.};
  
      for (i2=0;i2<Ny;i2++) 
      for (i1=0;i1<Nx;i1++) 
         Tabp(i1,i2,0) = f[0]*  Tabc(i1,i2-1,0,DEF_DIVROT_BORD) + f[1] *Tabc(i1,i2,0);
  
      for (i2=0;i2<Ny;i2++) 
      for (i1=0;i1<Nx;i1++)
         Tabp(i1,i2,1) = f[0] * Tabc(i1-1,i2,1,DEF_DIVROT_BORD) + f[1] * Tabc(i1,i2,1);
    }
    else  
    {
        // cout << " curl_point_value Border " << endl;
        int N=Nx-2;
        Tabp.resize(Nx-1,Ny-1,2);
        float *w1 = Tabp.buffer();
        float *w2 = Tabp.buffer() + (Nx-1)*(Ny-1);
        float *v1 = Tabc.buffer();
        float *v2 = Tabc.buffer() + Nx*Ny;
        rec01 (w1,v1,N);
        rec10 (w2,v2,N);
    }          
}

/*********************************************************************/

static void set_tab_ani(int NScale, int Np, intarray &TabPos, intarray &TabN)
{
   int j,Nx=Np;
   TabPos.alloc(NScale);
   TabN.alloc(NScale);
   
   for (j=0; j < NScale-1; j++)
   {
      int Nx2 = (Nx + 1) / 2;
	  TabN(j) = Nx / 2;
      TabPos(j) = Nx2;
	  Nx = Nx2;
   }
    TabN(NScale-1) = Nx;
	TabPos(NScale-1) = 0;
}

/***************************************/

static void indima(int Scale, int i,int j, int Numband, int Nl, int Nc, int &ind_i, int &ind_j)
// return the pixel position  (ind_i, ind_j) in an image, for a wavelet coefficient at scale Scale and position (i,j)
// This is only valid for fully decimated wavelet transforms
{
	ind_i = i;
	ind_j = j;
	int s,Nls,Ncs;
	
	Nls = Nl;
	Ncs = Nc;
	for (s=0; s <= Scale; s++)
	{
		Nls = (Nls+1) / 2;
		Ncs = (Ncs+1) / 2;
	}
    
    switch (Numband)
    {
        case 1:
			ind_j += Ncs;
			break;
		case 2:  
			ind_i += Nls;
			break;
		case 3:  
			ind_i += Nls;
			ind_j += Ncs;
			break;
		default: 
			break;
	}
}


/***************************************/

void MR_ISO_DIVCURL::read_from_tab(fltarray &Trans, Bool IsDivTrans)
{
	DivTrans = IsDivTrans;
 	fltarr2band(Trans);
}

/***************************************/

void MR_ISO_DIVCURL::band2fltarr(fltarray &Data, Bool DivRotTrans)
{
	// cout << "band2fltarr " << endl;
	// Copy the two components in a fltarr table
	// if DivRotTrans = true, we copy the div-free or the rot-free transform
	// while if DivRotTrans = false, we copy the two wavelet transform components
	
	// if the number of undecimated scales equals to 0 then 
	// the two transformation has the same size as the input data
	// and we create two images Dat(*,*,0) and Dat(*,*,1) which
	// will contain the two wavelet components.
	//  if the number of undecimated scales > 0 then we create
	// a cube for each component, and a frame of the cube will
	// contains a band
	
	if (NbrUndecimatedScale == 0)
	{
		Data.resize(Nc,Nl,2);
		for (int b=0; b < nbr_band(); b++)
		{
			int ind_i, ind_j, s = b/3;
			int TypeBand = b % 3 + 1;
			if (b == nbr_band() -1) TypeBand = 0;
			for (int i=0; i < size_band_nl(b); i++)
				for (int j=0; j < size_band_nc(b); j++)
				{
					// int jj = size_band_nc(b) - j - 1;
					// int ii = size_band_nl(b) - i - 1;
					indima(s,i,j, TypeBand,Nl,Nc,ind_i,ind_j);
					if (DivRotTrans == True)
					{
						Data(ind_j,ind_i,0) = (TabTransDivRot_0[b])(i,j);
						Data(ind_j,ind_i,1) = (TabTransDivRot_1[b])(i,j);
					}
					else
					{
						Data(ind_j,ind_i,0) = (TabWT_Trans_0[b])(i,j);
						Data(ind_j,ind_i,1) = (TabWT_Trans_1[b])(i,j);
					}
				}
		}
	}
	else
	{
		int CubeDimNz = nbr_band() * 2;
		Data.resize(Nc,Nl,CubeDimNz);
		for (int b=0; b < nbr_band(); b++)
			for (int i=0; i < size_band_nl(b); i++)
				for (int j=0; j < size_band_nc(b); j++)
				{
					if (DivRotTrans == True)
					{
						Data(j,i,b) = (TabTransDivRot_0[b])(i,j);
						Data(j,i,b+nbr_band()) = (TabTransDivRot_1[b])(i,j);
					}
					else
					{
						Data(j,i,b) = (TabWT_Trans_0[b])(i,j);
						Data(j,i,b+nbr_band()) = (TabWT_Trans_1[b])(i,j);
					}
				}
	}
	//  cout << "endband2fltarr " << endl;
}
/***************************************/

void MR_ISO_DIVCURL::fltarr2band(fltarray &Data, Bool DivRotTrans)
{
	if (NbrUndecimatedScale == 0)
	{
		for (int b=0; b < nbr_band(); b++)
		{
			int ind_i, ind_j, s = b/3;
			int TypeBand = b % 3 + 1;
			if (b == nbr_band() -1) TypeBand = 0;
			for (int i=0; i < size_band_nl(b); i++)
				for (int j=0; j < size_band_nc(b); j++)
				{
					// int jj = size_band_nc(b) - j - 1;
					// int ii = size_band_nl(b) - i - 1;
					indima(s,i,j, TypeBand,Nl,Nc,ind_i,ind_j);
					if (DivRotTrans == True)
					{
						(TabTransDivRot_0[b])(i,j) = Data(ind_j,ind_i,0);
						(TabTransDivRot_1[b])(i,j) = Data(ind_j,ind_i,1);
					}
					else
					{
						(TabWT_Trans_0[b])(i,j) = Data(ind_j,ind_i,0);
						(TabWT_Trans_1[b])(i,j) = Data(ind_j,ind_i,1);
					}
				}
		}
	}
	else
	{
		for (int b=0; b < nbr_band(); b++)
		{
			if (DivRotTrans == True)
			{
				(TabTransDivRot_0[b]).init();
				(TabTransDivRot_1[b]).init();
			}
			else
			{
				( TabWT_Trans_0[b]).init();
				( TabWT_Trans_1[b]).init();
			}
			
			for (int i=0; i < size_band_nl(b); i++)
				for (int j=0; j < size_band_nc(b); j++)
				{
					if (DivRotTrans == True)
					{
						(TabTransDivRot_0[b])(i,j) = Data(j,i,b);
						(TabTransDivRot_1[b])(i,j) = Data(j,i,b+nbr_band());
					}
					else
					{
						(TabWT_Trans_0[b])(i,j) = Data(j,i,b);
						(TabWT_Trans_1[b])(i,j) = Data(j,i,b+nbr_band());
					}
				}
		}  
	}
}

/***************************************/

inline float & MR_ISO_DIVCURL::operator() (int Num, int b, int i, int j) const
{
	if (Num == 0) return TabTransDivRot_0[b] (i,j);
	else return TabTransDivRot_1[b] (i,j);
}


/***************************************/

void MR_ISO_DIVCURL::alloc(int Nli, int Nci, int NbrUndec, int NbrPlan)
{
    init_filter_bank();
	Nl=Nli;
	Nc=Nci;
 	DivTrans=True;
	NbrUndecimatedScale=NbrUndec;
	if (NbrPlan >= 2) NbrScale=NbrPlan;
	else 
	{
		float N1=MIN(Nl,Nc);
		NbrScale = (int) (log(N1) / log(2.)+1);
	}
	if (Verbose == True) 
		cout << "NbrScale = " << NbrScale << " " << " NbrUndecimatedScale = " << NbrUndecimatedScale << endl;
        
	NbrBand = (WT_1->alloc)(TabWT_Trans_0, Nl, Nc, NbrScale, NbrUndecimatedScale);
	NbrBand = (WT_2->alloc)(TabWT_Trans_1, Nl, Nc, NbrScale, NbrUndecimatedScale);
	NbrBand = (WT_1->alloc)(TabTransDivRot_0, Nl, Nc, NbrScale, NbrUndecimatedScale);
	NbrBand = (WT_2->alloc)(TabTransDivRot_1, Nl, Nc, NbrScale, NbrUndecimatedScale);
} 

/***************************************/

void MR_ISO_DIVCURL::init_filter_bank()
{
    FilterAnaSynt *PtrFAS_1 = NULL;
    FAS_1.Verbose = Verbose;
    type_sb_filter SB_1 = F_5_3;
    FAS_1.alloc(SB_1);
    PtrFAS_1 = &FAS_1;
    SBF_1 = new SubBandFilter(FAS_1, Norm);
    // SBF_1 ->Border = DEF_DIVROT_BORD;
    SBF_1 ->setBorder(DEF_DIVROT_BORD);

	
    FilterAnaSynt *PtrFAS_2 = NULL;
    FAS_2.Verbose = Verbose;
    type_sb_filter SB_2 = F_4_4;
    FAS_2.alloc(SB_2);
    PtrFAS_2 = &FAS_2;
    SBF_2 = new SubBandFilter(FAS_2, Norm);
    // SBF_2 ->Border = DEF_DIVROT_BORD;
    SBF_2 ->setBorder(DEF_DIVROT_BORD);

    WT_1 = new HALF_DECIMATED_2D_WT(*SBF_2, *SBF_1);
    WT_2 = new HALF_DECIMATED_2D_WT(*SBF_1, *SBF_2);
	// HALF_DECIMATED_2D_WT WT_Curl(SubBand1D &SB1D_LINE, SubBand1D &SB1D_COL);
	
	NbrScale=0;
} 

/***************************************/

void MR_ISO_DIVCURL::transform(fltarray &Data, Bool Div)
{
    fltarray Tabc;
	
	DivTrans = Div;
    if (Verbose == True)
    {
        if (DivTrans == True) cout << " Anisotropic DIV transform ... " << endl;
 	    else cout << " Anisotropic CURL transform ... " << endl;
    }

   if (DivTrans == True) 
   {
      if (V0_Proj == True) 
	  {
	      if (Verbose == True) cout << " V0 projection ... " << endl;
	      div_quasi_interpol(Data,Tabc);
		  div_isotrop_transform (Tabc, RotDiv_Trans);
	  }
	  else  div_isotrop_transform (Data, RotDiv_Trans);
   }
   else 
   {
      if (V0_Proj == True) 
	  {
	      if (Verbose == True) cout << " V0 projection ... " << endl;
	      curl_quasi_interpol(Data,Tabc);
          rot_isotrop_transform (Tabc, RotDiv_Trans);
	   }
	   else rot_isotrop_transform (Data, RotDiv_Trans);
   }
}

/***************************************/

// void MR_ISO_DIVCURL::component_recons(fltarray &Data)
// {
// }

/***************************************/

void MR_ISO_DIVCURL::recons(fltarray &Data, Bool Div)
{
    fltarray Tabc;
	DivTrans = Div;
    Data.resize(Nc, Nl, 2);
    if (Verbose == True)
    {
        if (DivTrans == True)  cout << " Anisotropic DIV reconstruction ... " << endl;
 	    else cout << " Anisotropic CURL reconstruction ... " << endl;
    }
   
   if (DivTrans == True)  
   {
      div_isotrop_recons (RotDiv_Trans, Data);
	  if (PointValueRec == True)  
	  {
	       Tabc = Data;
	       div_point_value(Tabc, Data);
	  }
   }
   else 
   {
      rot_isotrop_recons (RotDiv_Trans, Data);
      if (PointValueRec == True) 
	  {
	       Tabc = Data;
	       curl_point_value(Tabc, Data);
	  }
   }
}

/***************************************/

void MR_ISO_DIVCURL::info()
{
  // For each band, write statistical information
  for (int b=0; b < nbr_band(); b++)
  {
      cout << " STAT BAND " << b+1 << " Nl = " <<  size_band_nl(b) << " Nc = " <<  size_band_nc(b) <<endl;
//     (TabTransDivRot_0[b]).info(" Component 1 : ");
//     (TabTransDivRot_1[b]).info(" Component 2 : ");
  }
  

}

/***************************************/

void MR_ISO_DIVCURL::div_isotrop_transform (fltarray &Data, fltarray &Trans)
{
	int b;
    Ifloat Imag;
    float *BuffData = Data.buffer();
    Imag.alloc (BuffData,Nl,Nc);
	(WT_1->transform)(Imag, TabWT_Trans_0, NbrScale, NbrUndecimatedScale);
	
	BuffData += Nl*Nc;
	Imag.alloc (BuffData,Nl,Nc);
	(WT_2->transform)(Imag, TabWT_Trans_1, NbrScale, NbrUndecimatedScale);
	// cout << " nbr = " <<  NbrScale << " " <<  nbr_band() << endl;
	
	//   Ifloat Buf;
	//    Buf = Imag;
	//    Buf.init();
	//    (WT_2->recons)(TabWT_Trans_1, Buf, NbrScale, NbrUndecimatedScale);
	//    Buf = Imag - Buf;
	//    Buf.info();
	//    exit(0);
	
	// fltarray Buf;
	// band2fltarr(Buf,False);
	// fits_write_fltarr("xx_wt1", Data);
	
	
	
	// Parameter renormalization
	int Step=1;
	if (NbrUndecimatedScale > 0) Step = 2;
	for (b=0; b < nbr_band()-1; b++)
	{
		int Nlb=size_band_nl(b);
		int Ncb=size_band_nc(b);
		// cout << b << ": Band size " << size_band_nl(b) << " " << size_band_nc(b) << endl;
		
		int Ind, s = b/3;
		int j1 = Nc / POW(2,s+1);
		int j2 = Nl / POW(2,s+1);
		int TypeBand = b % 3; 
		float Den = 1. / (float)(j1*j1+j2*j2);
		for (int i=0; i < Nlb; i++)
			for (int j=0; j < Ncb; j++)
			{
				// if (( i == 1) && ( j == 1)) cout << " j = " << s << " j1 = " << j1 << " j2 = " << j2 << " " << float(j1) / (j1*j1+j2*j2) << " " << float(j2) / (j1*j1+j2*j2) << endl;
				
				switch (TypeBand) // 0 for Horizontal, 1 for vertical and 2 for diagonal
				{
					case 0: (TabTransDivRot_0[b])(i,j) =  - (TabWT_Trans_1[b])(i,j) /4.;
						Ind = (i-Step+Nlb) % Nlb;
						(TabTransDivRot_1[b])(i,j) =  (TabWT_Trans_0[b])(i,j) +  ((TabWT_Trans_1[b])(i,j) - (TabWT_Trans_1[b])(Ind,j))/4.;	   
						break;
					case 1: (TabTransDivRot_0[b])(i,j) =  - (TabWT_Trans_0[b])(i,j) /4.;
						Ind = (j-Step+Ncb) % Ncb;
						(TabTransDivRot_1[b])(i,j) =  (TabWT_Trans_1[b])(i,j) +  ((TabWT_Trans_0[b])(i,j) - (TabWT_Trans_0[b])(i,Ind))/4.;
						break;
					case 2: (TabTransDivRot_0[b])(i,j) = (j2*(TabWT_Trans_0[b])(i,j) - j1*(TabWT_Trans_1[b])(i,j)) * Den;
						(TabTransDivRot_1[b])(i,j) = (j1*(TabWT_Trans_0[b])(i,j) + j2*(TabWT_Trans_1[b])(i,j)) * Den;
						break;
					default: cout << "Error: bad band type ... " << endl;
						exit(-1);
						break;
				}
			}
		if ((s > 0) && (b % 3 == 0) && (NbrUndecimatedScale > s)) Step = Step * 2;
	}
	b = nbr_band()-1;
	// cout << b << ": Last Band size " << size_band_nl(b) << " " << size_band_nc(b) << endl;
	int Nlb=size_band_nl(b);
	int Ncb=size_band_nc(b);
	for (int i=0; i < Nlb; i++)
		for (int j=0; j < Ncb; j++)
		{
			// cout << "BB" << endl;
			(TabTransDivRot_0[b])(i,j) = (TabWT_Trans_0[b])(i,j);
		//	(TabTransDivRot_0[b-3])(i,j) = (TabWT_Trans_1[b-3])(i,j);
		//	(TabTransDivRot_0[b-2])(i,j) = (TabWT_Trans_0[b-2])(i,j);
		//	(TabTransDivRot_0[b-1])(i,j) = (TabWT_Trans_0[b-1])(i,j);
			
			(TabTransDivRot_1[b])(i,j) = (TabWT_Trans_1[b])(i,j);
		//	(TabTransDivRot_1[b-3])(i,j) = (TabWT_Trans_0[b-3])(i,j);
		//	(TabTransDivRot_1[b-2])(i,j) = (TabWT_Trans_1[b-2])(i,j); 
		//	(TabTransDivRot_1[b-1])(i,j) = (TabWT_Trans_1[b-1])(i,j); 
		}

	band2fltarr(Trans);
}


/***************************************/
		   
void MR_ISO_DIVCURL::div_isotrop_recons (fltarray &Trans, fltarray &Data)
{ 
	int b;
	Ifloat Imag;
	Data.resize(Nc,Nl,2);
	
	fltarr2band(Trans);
	
	b = nbr_band()-1;
	int Nlb=size_band_nl(b);
	int Ncb=size_band_nc(b);
	
	// cout << "Last Band size " << size_band_nl(b) << " " << size_band_nc(b) << endl;
	for (int i=0; i < Nlb; i++)
		for (int j=0; j < Ncb; j++)
		{
			(TabWT_Trans_0[b])(i,j) = (TabTransDivRot_0[b])(i,j);
		//	(TabWT_Trans_1[b-3])(i,j) = (TabTransDivRot_0[b-3])(i,j);
		//	(TabWT_Trans_0[b-2])(i,j) = (TabTransDivRot_0[b-2])(i,j);
		//	(TabWT_Trans_0[b-1])(i,j) = (TabTransDivRot_0[b-1])(i,j);
			
			(TabWT_Trans_1[b])(i,j) = (TabTransDivRot_1[b])(i,j);
		//	(TabWT_Trans_0[b-3])(i,j) = (TabTransDivRot_1[b-3])(i,j);
		//	(TabWT_Trans_1[b-2])(i,j) = (TabTransDivRot_1[b-2])(i,j); 
		//	(TabWT_Trans_1[b-1])(i,j) = (TabTransDivRot_1[b-1])(i,j); 
		}
	
	int Step=1;
	if (NbrUndecimatedScale > 0) Step = 2;  
	for (b=0; b < nbr_band()-1; b++)
	{
		int Nlb=size_band_nl(b);
		int Ncb=size_band_nc(b);
		// cout << b << ": Band size " << size_band_nl(b) << " " << size_band_nc(b) << endl;
		
		int Ind, s = b/3;
		int j1 = Nc / POW(2,s+1);
		int j2 = Nl / POW(2,s+1);
		int TypeBand = b % 3; 
		//      float Den = 1. / (float)(j1*j1+j2*j2);
		for (int i=0; i < Nlb; i++)
			for (int j=0; j < Ncb; j++)
			{
				// if (( i == 1) && ( j == 1)) cout << " j = " << s << " j1 = " << j1 << " j2 = " << j2 << " " << float(j1) / (j1*j1+j2*j2) << " " << float(j2) / (j1*j1+j2*j2) << endl;
				
				switch (TypeBand) // 0 for Horizontal, 1 for vertical and 2 for diagonal
				{
					case 0: (TabWT_Trans_1[b])(i,j) = -4. * (TabTransDivRot_0[b])(i,j) ;
						Ind = (i-Step+Nlb) % Nlb;
						(TabWT_Trans_0[b])(i,j) =  (TabTransDivRot_1[b])(i,j) + (TabTransDivRot_0[b])(i,j) - (TabTransDivRot_0[b])(Ind,j);
						// ddn[1][i1][i2]+ddn[0][i1][i2]-ddn[0][i1][(i2-1+j)%(j)];
						// (TabWT_Trans_0[b])(i,j) +  ((TabWT_Trans_1[b])(i,j) - (TabWT_Trans_1[b])(Ind,j))/4.;	   
						break;
					case 1: (TabWT_Trans_0[b])(i,j) = -4. * (TabTransDivRot_0[b])(i,j);
						Ind = (j-Step+Ncb) % Ncb;
						(TabWT_Trans_1[b])(i,j) = (TabTransDivRot_1[b])(i,j) + (TabTransDivRot_0[b])(i,j) - (TabTransDivRot_0[b])(i,Ind);
						// c[1][i2][i1]=ddn[1][i2][i1]+ddn[0][i2][i1]-ddn[0][(i2-1+j)%(j)][i1];
						// (TabWT_Trans_1[b])(i,j) =  (TabTransDivRot_1[b])(i,j) +  ((TabWT_Trans_0[b])(i,j) - (TabWT_Trans_0[b])(i,Ind))/4.;
						break;
					case 2: (TabWT_Trans_0[b])(i,j) =  j2*(TabTransDivRot_0[b])(i,j) + j1*(TabTransDivRot_1[b])(i,j);
						(TabWT_Trans_1[b])(i,j) = -j1*(TabTransDivRot_0[b])(i,j) + j2*(TabTransDivRot_1[b])(i,j);
						break;
						
					default: cout << "Error: bad band type ... " << endl;
						exit(-1);
						break;
				}
			}
		if ((s > 0) && (b % 3 == 0) && (NbrUndecimatedScale > s)) Step = Step * 2;
	}
	
	
	// fltarray Buf;
	// band2fltarr(Buf,False);
	// fits_write_fltarr((char *) "xx_wt2", Data);
	
	float *BuffData = Data.buffer();
	Imag.alloc (BuffData,Nl,Nc);
	(WT_1->recons)(TabWT_Trans_0, Imag, NbrScale, NbrUndecimatedScale);
	BuffData += Nl*Nc;
	Imag.alloc (BuffData,Nl,Nc);
	(WT_2->recons)(TabWT_Trans_1, Imag, NbrScale, NbrUndecimatedScale);
}

/***************************************/

void MR_ISO_DIVCURL::rot_isotrop_transform (fltarray &Data, fltarray &Trans)
{
	int b;
    Ifloat Imag;
    float *BuffData = Data.buffer();
    Imag.alloc (BuffData,Nl,Nc);
	(WT_2->transform)(Imag, TabWT_Trans_0, NbrScale, NbrUndecimatedScale);
	BuffData += Nl*Nc;
	Imag.alloc (BuffData,Nl,Nc);
	(WT_1->transform)(Imag, TabWT_Trans_1, NbrScale, NbrUndecimatedScale);
	
	// Parameter renormalization
	int Step=1;
	if (NbrUndecimatedScale > 0) Step = 2;
	
	for (b=0; b < nbr_band()-1; b++)
	{
		int Nlb=size_band_nl(b);
		int Ncb=size_band_nc(b);
		int Ind, s = b/3;
		int j1 = Nc / POW(2,s+1);
		int j2 = Nl / POW(2,s+1);
		int TypeBand = b % 3; 
		float Den = 1. / (float)(j1*j1+j2*j2);
		for (int i=0; i < Nlb; i++)
			for (int j=0; j < Ncb; j++)
			{
				// if (( i == 1) && ( j == 1)) cout << " j = " << s << " j1 = " << j1 << " j2 = " << j2 << " " << float(j1) / (j1*j1+j2*j2) << " " << float(j2) / (j1*j1+j2*j2) << endl;
				// ddn[0][i1][i2]=d[1][i1][i2]-d[0][i1][i2]/4.+d[0][i1][(i2-1+j)%(j)]/4.;
				// ddn[1][i1][i2]=d[0][i1][i2]/4.;
				// ddn[0][i2][i1]=d[0][i2][i1]-d[1][i2][i1]/4.+d[1][(i2-1+j)%(j)][i1]/4.;
				// ddn[1][i2][i1]=d[1][i2][i1]/4.;
				
				// ddn[0][j1+i1][j2+i2]=(j2*d[0][j1+i1][j2+i2]-j1*d[1][j1+i1][j2+i2])/(j1*j1+j2*j2);
				//  ddn[1][j1+i1][j2+i2]=(j1*d[0][j1+i1][j2+i2]+j2*d[1][j1+i1][j2+i2])/(j1*j1+j2*j2);
				
				switch (TypeBand)   // 0 for Horizontal, 1 for vertical and 2 for diagonal
				{
					case 0:  Ind = (i-Step+Nlb) % Nlb;
						(TabTransDivRot_0[b])(i,j) = (TabWT_Trans_1[b])(i,j) - (TabWT_Trans_0[b])(i,j)/ 4. + (TabWT_Trans_0[b])(Ind, j)/4.;
						(TabTransDivRot_1[b])(i,j) = (TabWT_Trans_0[b])(i,j) / 4;  
						break;
					case 1: Ind = (j-Step+Ncb) % Ncb;
						(TabTransDivRot_0[b])(i,j) = (TabWT_Trans_0[b])(i,j) - (TabWT_Trans_1[b])(i,j)/ 4. + (TabWT_Trans_1[b])(i,Ind)/4.;
						(TabTransDivRot_1[b])(i,j) = (TabWT_Trans_1[b])(i,j)/ 4;
						break;
					case 2: (TabTransDivRot_0[b])(i,j) = (j2*(TabWT_Trans_0[b])(i,j) - j1*(TabWT_Trans_1[b])(i,j)) * Den;
						(TabTransDivRot_1[b])(i,j) = (j1*(TabWT_Trans_0[b])(i,j) + j2*(TabWT_Trans_1[b])(i,j)) * Den;
						break;
						
					default: cout << "Error: bad band type ... " << endl;
						exit(-1);
						break;
				}
			}
		if ((s > 0) && (b % 3 == 0) && (NbrUndecimatedScale > s)) Step = Step * 2;
	}
	b = nbr_band()-1;
	// cout << "Last Band size " << size_band_nl(b) << " " << size_band_nc(b) << endl;
	int Nlb=size_band_nl(b);
	int Ncb=size_band_nc(b);
	for (int i=0; i < Nlb; i++)
		for (int j=0; j < Ncb; j++)  
		{
			// cout << "BB" << endl;
			(TabTransDivRot_0[b])(i,j) = (TabWT_Trans_0[b])(i,j);
	//		(TabTransDivRot_0[b-3])(i,j) = (TabWT_Trans_1[b-3])(i,j);
	//		(TabTransDivRot_0[b-2])(i,j) = (TabWT_Trans_0[b-2])(i,j);
	//		(TabTransDivRot_0[b-1])(i,j) = (TabWT_Trans_0[b-1])(i,j);
			
			(TabTransDivRot_1[b])(i,j) = (TabWT_Trans_1[b])(i,j);
	//		(TabTransDivRot_1[b-3])(i,j) = (TabWT_Trans_0[b-3])(i,j);
	//		(TabTransDivRot_1[b-2])(i,j) = (TabWT_Trans_1[b-2])(i,j); 
	//		(TabTransDivRot_1[b-1])(i,j) = (TabWT_Trans_1[b-1])(i,j); 
		}
	band2fltarr(Trans);
	
}


/***************************************/

void MR_ISO_DIVCURL::rot_isotrop_recons (fltarray &Trans, fltarray &Data)
{ 
	int b;
	Ifloat Imag;
	Data.resize(Nc,Nl,2);
	
 	fltarr2band(Trans);
	b = nbr_band()-1;
	int Nlb=size_band_nl(b);
	int Ncb=size_band_nc(b);
	// cout << "Last Band size " << size_band_nl(b) << " " << size_band_nc(b) << endl;
	for (int i=0; i < Nlb; i++)
		for (int j=0; j < Ncb; j++)  
		{
			(TabWT_Trans_0[b])(i,j) = (TabTransDivRot_0[b])(i,j);
		//	(TabWT_Trans_1[b-3])(i,j) = (TabTransDivRot_0[b-3])(i,j);
		//	(TabWT_Trans_0[b-2])(i,j) = (TabTransDivRot_0[b-2])(i,j);
		//	(TabWT_Trans_0[b-1])(i,j) = (TabTransDivRot_0[b-1])(i,j);
			
			(TabWT_Trans_1[b])(i,j) = (TabTransDivRot_1[b])(i,j);
		//	(TabWT_Trans_0[b-3])(i,j) = (TabTransDivRot_1[b-3])(i,j);
		//	(TabWT_Trans_1[b-2])(i,j) = (TabTransDivRot_1[b-2])(i,j); 
		//	(TabWT_Trans_1[b-1])(i,j) = (TabTransDivRot_1[b-1])(i,j); 
		}
	
	int Step=1;
	if (NbrUndecimatedScale > 0) Step = 2;    
	for (b=0; b < nbr_band()-1; b++)
	{
		int Nlb=size_band_nl(b);
		int Ncb=size_band_nc(b);
		int Ind, s = b/3;
		int j1 = Nc / POW(2,s+1);
		int j2 = Nl / POW(2,s+1);
		int TypeBand = b % 3; 
		//      float Den = 1. / (float)(j1*j1+j2*j2);
		for (int i=0; i < Nlb; i++)
			for (int j=0; j < Ncb; j++)
			{
				// if (( i == 1) && ( j == 1)) cout << " j = " << s << " j1 = " << j1 << " j2 = " << j2 << " " << float(j1) / (j1*j1+j2*j2) << " " << float(j2) / (j1*j1+j2*j2) << endl;
				
				switch (TypeBand)
				{
					case 0: (TabWT_Trans_0[b])(i,j) = 4. * (TabTransDivRot_1[b])(i,j) ;
						Ind = (i-Step+Nlb) % Nlb;
						(TabWT_Trans_1[b])(i,j) =  (TabTransDivRot_0[b])(i,j) + (TabTransDivRot_1[b])(i,j) - (TabTransDivRot_1[b])(Ind,j);
						// (TabTransDivRot_1[b])(i,j) + (TabTransDivRot_0[b])(i,j) - (TabTransDivRot_0[b])(Ind,j);
						break;
					case 1: (TabWT_Trans_1[b])(i,j) = 4. * (TabTransDivRot_1[b])(i,j);
						Ind = (j-Step+Ncb) % Ncb;
						(TabWT_Trans_0[b])(i,j) = (TabTransDivRot_0[b])(i,j) + (TabTransDivRot_1[b])(i,j) - (TabTransDivRot_1[b])(i,Ind);
						// c[1][i2][i1]=ddn[1][i2][i1]+ddn[0][i2][i1]-ddn[0][(i2-1+j)%(j)][i1];
						// (TabWT_Trans_1[b])(i,j) =  (TabTransDivRot_1[b])(i,j) +  ((TabWT_Trans_0[b])(i,j) - (TabWT_Trans_0[b])(i,Ind))/4.;
						break;
					case 2: (TabWT_Trans_0[b])(i,j) =  j2*(TabTransDivRot_0[b])(i,j) + j1*(TabTransDivRot_1[b])(i,j);
						(TabWT_Trans_1[b])(i,j) = -j1*(TabTransDivRot_0[b])(i,j) + j2*(TabTransDivRot_1[b])(i,j);
						break;
						
					default: cout << "Error: bad band type ... " << endl;
						exit(-1);
						break;
				}
			}
		if ((s > 0) && (b % 3 == 0) && (NbrUndecimatedScale > s)) Step = Step * 2;
	}
	
	float *BuffData = Data.buffer();
	Imag.alloc (BuffData,Nl,Nc);
	(WT_2->recons)(TabWT_Trans_0, Imag, NbrScale, NbrUndecimatedScale);
	BuffData += Nl*Nc;
	Imag.alloc (BuffData,Nl,Nc);
	(WT_1->recons)(TabWT_Trans_1, Imag, NbrScale, NbrUndecimatedScale);
}

/***************************************/

void MR_ISO_DIVCURL::component_recons(fltarray &Data)
{ 
	Ifloat Imag;
	Data.resize(Nc,Nl,2);
	
	float *BuffData = Data.buffer();
	Imag.alloc (BuffData,Nl,Nc);
	(WT_2->recons)(TabTransDivRot_1, Imag, NbrScale, NbrUndecimatedScale);
	// INFO(Imag, (char *) "Comp 1");
	BuffData += Nl*Nc;
	Imag.alloc (BuffData,Nl,Nc);
	(WT_1->recons)(TabTransDivRot_0, Imag, NbrScale, NbrUndecimatedScale);
	// INFO(Imag, (char *) "Comp 2");
}

/***************************************/
/***************************************/
/***************************************/

void MR_ANI_DIVCURL::alloc(int Nli, int Nci, int NbrPlan)
{
	
	init_filter_bank();
	Nl=Nli;
	Nc=Nci;
	DivTrans=True;
	if (NbrPlan >= 2) NbrScaleX=NbrScaleY=NbrPlan;
	else 
	{
 		NbrScaleX = (int) (log(Nc) / log(2.)+1);
		NbrScaleY = (int) (log(Nl) / log(2.)+1);
	}
	if (Verbose == True) 
		cout << "NbrScaleX = " << NbrScaleX << " " <<  ", NbrScaleY = " << NbrScaleY <<endl;
	
  	NbrBand = NbrScaleX*NbrScaleY;
	set_tab_ani(NbrScaleX, Nc, TabPosX, TabNx);
	set_tab_ani(NbrScaleY, Nl, TabPosY, TabNy);
	
} 

/***************************************/

void MR_ANI_DIVCURL::init_filter_bank()
{
	FilterAnaSynt *PtrFAS_1 = NULL;
    FAS_1.Verbose = Verbose;
    type_sb_filter SB_1 = F_5_3;
    FAS_1.alloc(SB_1);
    PtrFAS_1 = &FAS_1;
    SBF_1 = new SubBandFilter(FAS_1, Norm);
    SBF_1 ->setBorder(DEF_DIVROT_BORD);
	
    FilterAnaSynt *PtrFAS_2 = NULL;
    FAS_2.Verbose = Verbose;
    type_sb_filter SB_2 = F_4_4;
    FAS_2.alloc(SB_2);
    PtrFAS_2 = &FAS_2;
    SBF_2 = new SubBandFilter(FAS_2, Norm);
    // SBF_2 ->Border = DEF_DIVROT_BORD;
    SBF_2 ->setBorder(DEF_DIVROT_BORD);
	
    Bool UseMirror=False;
	WT_1 = new BordLineCol;
	WT_2 = new BordLineCol;
	(WT_1->alloc)(*SBF_2, *SBF_1, UseMirror);
    (WT_2->alloc)(*SBF_1, *SBF_2, UseMirror);     
    NbrScale=0;
}

/***************************************/

void MR_ANI_DIVCURL::transform(fltarray &Data, Bool Div)
{
    fltarray Tabc;
	
	DivTrans = Div;
    if (Verbose == True)
    {
        if (DivTrans == True) cout << " Anisotropic DIV transform ... " << endl;
 	    else cout << " Anisotropic CURL transform ... " << endl;
    }
	
	if (DivTrans == True) 
	{
		if (V0_Proj == True) 
		{
	      if (Verbose == True) cout << " V0 projection ... " << endl;
			div_quasi_interpol(Data,Tabc,BorderWavelet);
            fits_write_fltarr("xx_v0div.fits", Tabc);
			if (BorderWavelet == False) div_transform (Tabc, RotDiv_Trans);
            else div_bord_transform(Tabc, RotDiv_Trans);
		}
		else  
        {
           if (BorderWavelet == False) div_transform (Data, RotDiv_Trans);
           else div_bord_transform(Data, RotDiv_Trans);
        }
	}
	else 
	{
		if (V0_Proj == True) 
	 	{
	       if (Verbose == True) cout << " V0 projection ... " << endl;
			curl_quasi_interpol(Data,Tabc,BorderWavelet);
           if (BorderWavelet == False) curl_transform (Tabc, RotDiv_Trans);
           else curl_bord_transform(Tabc, RotDiv_Trans);
		}
		else 
        {
           if (BorderWavelet == False) curl_transform (Data, RotDiv_Trans);
           else curl_bord_transform(Data, RotDiv_Trans);
        }
	}
}

/***************************************/

// void MR_ANI_DIVCURL::component_recons(fltarray &Data)
// {
// }

/***************************************/

void MR_ANI_DIVCURL::recons(fltarray &Data, Bool Div)
{
    fltarray Tabc;
	DivTrans = Div;
    Data.resize(Nc, Nl, 2);
    if (Verbose == True)
    {
        if (DivTrans == True)  cout << " Anisotropic DIV reconstruction ... " << endl;
 	    else cout << " Anisotropic CURL reconstruction ... " << endl;
    }
	
	if (DivTrans == True)  
	{
		if (BorderWavelet == False) div_recons (RotDiv_Trans, Data);
        else div_bord_recons(RotDiv_Trans, Data);
		if (PointValueRec == True)  
		{
			Tabc = Data;
            if (Verbose == True)  cout << " Point Value ... " << endl;
			div_point_value(Tabc, Data, BorderWavelet);
		}
	}
	else 
	{
		if (BorderWavelet == False) curl_recons (RotDiv_Trans, Data);
        else curl_bord_recons(RotDiv_Trans, Data);
		if (PointValueRec == True) 
		{
			Tabc = Data;
            if (Verbose == True)  cout << " Point Value ... " << endl;
			curl_point_value(Tabc, Data, BorderWavelet);
		}
	}
}

/***************************************/

void MR_ANI_DIVCURL::info()
{
	// For each band, write statistical information
	for (int b=0; b < nbr_band(); b++)
	{
		cout << " STAT BAND " << b+1 << " Nl = " <<  size_band_nl(b) << " Nc = " <<  size_band_nc(b) <<endl;
		//     (TabTransDivRot_0[b]).info(" Component 1 : ");
		//     (TabTransDivRot_1[b]).info(" Component 2 : ");
	}
	
	for (int bx=0; bx < nbr_band_nx(); bx++)
		for (int by=0; by < nbr_band_ny(); by++)
		{
			cout << " STAT BAND " << bx+1 << "," << by+1 << ":  Nx = " <<  size_band_nx(bx) << " Ny = " <<  size_band_nx(by) << ",  Pos = " <<  TabPosX(bx) << "," << TabPosX(by) <<  endl;
			//     (TabTransDivRot_0[b]).info(" Component 1 : ");
			//     (TabTransDivRot_1[b]).info(" Component 2 : ");
		}
}

  
/***************************************/

void MR_ANI_DIVCURL::div_bord_transform (fltarray &Data1, fltarray &Data2, fltarray &Trans1, fltarray &Trans2)
{
cout << "DDDD " << endl;
   int N = Data1.nx() - 2;
   int Nx1= Data1.nx();
   int Ny1 = Data1.ny();
   int Nx2= Data2.nx();
   int Ny2 = Data2.ny();
    Trans1.alloc(Nx1,Ny2);
    Trans2.alloc(Nx2,Ny1);
    // tfodab2d(Nx1-2, Data1.buffer(),Data2.buffer(),Trans1.buffer(),Trans2.buffer());
   tfodab2d(Data1.buffer(),Data2.buffer(),Trans1.buffer(),Trans2.buffer(), N);
}

void MR_ANI_DIVCURL::div_bord_transform (fltarray &Data, fltarray &Trans)
{
    cout << "2DDDD " << endl;
    int N = Data.nx() - 2;
    Trans.alloc(N+2,N+2,2);
    float *Data1 = Data.buffer();
    float *Data2 = Data.buffer() + Data.nx() * Data.ny();
    float *Trans1 = Trans.buffer();
    float *Trans2 = Trans.buffer() + Trans.nx() * Trans.ny();
    tfodab2d(Data1,Data2,Trans1,Trans2, N);
}

/***************************************/


void MR_ANI_DIVCURL::div_bord_recons ( fltarray &Data1, fltarray &Data2, fltarray &Trans1, fltarray &Trans2)
{
    int Nx1= Trans1.nx();
    int Ny1 = Trans1.ny();
    int Nx2= Trans2.nx();
    int Ny2 = Trans2.ny();
    Data1.alloc(Nx1,Ny2);
    Data2.alloc(Nx2,Ny1);
   //  tfodab2dinv(Nx1-2,Data1.buffer(),Data2.buffer(), Trans1.buffer(),Trans2.buffer());
    tfodab2dinv(Data1.buffer(),Data2.buffer(),Trans1.buffer(),Trans2.buffer(), Nx1-2);
}

void MR_ANI_DIVCURL::div_bord_recons (fltarray &Trans, fltarray &Data)
{
    int N = Trans.nx() - 2;
    Data.alloc(N+2,N+2,2);
    float *Data1 = Data.buffer();
    float *Data2 = Data.buffer() + Data.nx() * Data.ny();
    float *Trans1 = Trans.buffer();
    float *Trans2 = Trans.buffer() + Trans.nx() * Trans.ny();
    tfodab2dinv(Data1, Data2, Trans1, Trans2, N);
}

/***************************************/

void MR_ANI_DIVCURL::curl_bord_transform (fltarray &Data1, fltarray &Data2, fltarray &Trans1, fltarray &Trans2)
{
    int Nx1= Data1.nx();
    int Ny1 = Data1.ny();
    int Nx2= Data2.nx();
    int Ny2 = Data2.ny();
    Trans1.alloc(Nx1,Ny2);
    Trans2.alloc(Nx2,Ny1);
    cout << "T1: " << Trans1.nx() << " " << Trans1.ny() << endl;
    cout << "T2: " << Trans2.nx() << " " << Trans2.ny() << endl;
    
    cout << "Nx1 = " << Nx1 << endl;
    tfogab2d(Data1.buffer(),Data2.buffer(),Trans1.buffer(),Trans2.buffer(), Nx1-2);
}

void MR_ANI_DIVCURL::curl_bord_transform (fltarray &Data, fltarray &Trans)
{
    int N = Data.nx() - 2;
    Trans.alloc(N+2,N+2,2);
    float *Data1 = Data.buffer();
    float *Data2 = Data.buffer() + Data.nx() * Data.ny();
    float *Trans1 = Trans.buffer();
    float *Trans2 = Trans.buffer() + Trans.nx() * Trans.ny();
    tfogab2d(Data1,Data2,Trans1,Trans2, N);
}

/***************************************/

void MR_ANI_DIVCURL::curl_bord_recons ( fltarray &Data1, fltarray &Data2, fltarray &Trans1, fltarray &Trans2)
{
    int Nx1= Trans1.nx();
    int Ny1 = Trans1.ny();
    int Nx2= Trans2.nx();
    int Ny2 = Trans2.ny();
    Data1.alloc(Nx1,Ny2);
    Data2.alloc(Nx2,Ny1);
    cout << "R1: " << Data1.nx() << " " << Data1.ny() << endl;
    cout << "RR2: " << Data2.nx() << " " << Data2.ny() << endl;
    
    //  tfodab2dinv(Nx1-2,Data1.buffer(),Data2.buffer(), Trans1.buffer(),Trans2.buffer());
    tfogab2dinv(Data1.buffer(),Data2.buffer(),Trans1.buffer(),Trans2.buffer(), Nx1-2);
    cout << "EE div_bord_transform " << endl;
    
}

void MR_ANI_DIVCURL::curl_bord_recons (fltarray &Trans, fltarray &Data)
{
    int N = Trans.nx() - 2;
    Data.alloc(N+2,N+2,2);
    float *Data1 = Data.buffer();
    float *Data2 = Data.buffer() + Data.nx() * Data.ny();
    float *Trans1 = Trans.buffer();
    float *Trans2 = Trans.buffer() + Trans.nx() * Trans.ny();
    tfogab2dinv(Data1, Data2, Trans1, Trans2, N);
}


/***************************************/

/***************************************/


void MR_ANI_DIVCURL::div_transform (fltarray &Data, fltarray &Trans)
{
    fltarray WT;
 	WT.alloc(Data.nx(), Data.ny(), 2);
 	Trans.alloc(Nc,Nl,2);
	// cout << Nc << " " << Nc << endl;
	for (int x=0; x < Data.nx(); x++)
		for (int y=0; y < Data.ny(); y++) 
		{
			WT(x,y,0) = Data(x,y,0);
			WT(x,y,1) = Data(x,y,1);
		}
    Ifloat Imag;
    float *BuffData = WT.buffer();
    Imag.alloc (BuffData,Nl,Nc);
	(WT_1->transform)(Imag, NbrScaleX, NbrScaleY);
	
    BuffData += Nl*Nc;
    Imag.alloc (BuffData,Nl,Nc);
    (WT_2->transform)(Imag, NbrScaleX, NbrScaleY);
 	
	//   fits_write_fltarr((char *) "a1div.fits", WT);
	
 
    
    if (Debug == 1)
    {
    // cout << "old" << endl;
	for (int j1=Nc/2; j1>0; j1=j1/2) 
		for (int j2=Nl/2; j2>0; j2=j2/2) 
		{
			double Omega = 1. / (j1*j1+j2*j2);
            // if (Verbose == True) cout << " j1 = " << j1 << ", j2 = " << j2 << ", nx = " << j1 << ", ny = " << j2 << endl;
			for (int i1=0;i1<j1;i1++) 
				for (int i2=0;i2<j2;i2++) 
				{
					Trans(j1+i1,j2+i2,0) = (j2 * WT(j1+i1,j2+i2,0) - j1 * WT(j1+i1,j2+i2,1)) / (j1*j1+j2*j2);
					Trans(j1+i1,j2+i2,1) = (j1 * WT(j1+i1,j2+i2,0) + j2 * WT(j1+i1,j2+i2,1)) / (j1*j1+j2*j2);
				}
		}
    }
    else
    {
      //  cout << "new" << endl;
        int j1;
        int j2;
        for (int bx=0, j1=Nc/2; bx < nbr_band_nx(); bx++,j1=j1/2)
        for (int by=0, j2=Nl/2; by < nbr_band_ny(); by++,j2=j2/2)
        {
			double Omega = 1. / (j1*j1+j2*j2);
            //if (Verbose == True) cout << " j1 = " << j1 << ", j2 = " << j2 << ", nx = " << size_band_nx(bx) << ", ny = " << size_band_ny(by) << endl;
            // if (Verbose == True) cout << " bx = " << TabPosX(bx) << ", by = " << TabPosY(by) << endl;
            if ((bx !=  nbr_band_nx()-1) || (by != nbr_band_ny()-1))
            {
                for (int i1=0;i1<size_band_nx(bx);i1++) 
                for (int i2=0;i2<size_band_ny(by);i2++) 
                {
                    int x = TabPosX(bx) + i1;
                    int y = TabPosY(by) + i2;
                    Trans(x,y,0) = (j2 * WT(x,y,0) - j1 * WT(x,y,1)) / (j1*j1+j2*j2);
                    Trans(x,y,1) = (j1 * WT(x,y,0) + j2 * WT(x,y,1)) / (j1*j1+j2*j2);
                }
            }
        }
       
    }
        
	// fits_write_fltarr((char *) "a1divbe.fits", Trans);
	// cout << "Div trans " << WT.naxis() << " " << WT.nx() <<  " " << WT.ny() << " " << WT.nz() << endl;
	// cout << "Div trans " << Trans.naxis() << " " << Trans.nx() <<  " " << Trans.ny() << " " << Trans.nz() << " " << Nc << endl;
	
	
	if (Debug == 1)
    {
     for (int i1=1; i1 < Nc; i1++) 
	 {
		Trans(i1,0,1) = WT(i1,0, 0);
		Trans(i1,0,0) = WT(i1,0,1);
		Trans(0,i1,1) = WT(0,i1,1);
		Trans(0,i1,0) = WT(0,i1,0);
	 }
	}
	 

	// cout << "Div trans1 " << endl;
	
	if (Debug == 1)
    {
		Trans(0,0,0) = WT(0,0,0);
	    Trans(0,0,1) = WT(0,0,1);
	}
	else {
		int bx = nbr_band_nx() -1;
		int by = nbr_band_ny() -1;
		for (int i1=0;i1<size_band_nx(bx);i1++) 
		for (int i2=0;i2<size_band_ny(by);i2++) 
		{
			int x = TabPosX(bx) + i1;
			int y = TabPosY(by) + i2;
			Trans(x,y,0) =  WT(x,y,0);
			Trans(x,y,1) =  WT(x,y,1);
		}
	}

	// fits_write_fltarr((char *) "a1dive.fits", Trans);
	// cout << "Div trans2 " << endl;
}

/***************************************/

void MR_ANI_DIVCURL::div_recons (fltarray &Trans, fltarray &Data)
{
	
	/*    c[0][0][0]=ddn[0][0][0];
	 c[1][0][0]=ddn[1][0][0];
	 for (i1=1;i1<N;i1++) {
	 c[0][i1][0]=ddn[0][i1][0];
	 c[1][i1][0]=ddn[1][i1][0];
	 c[1][0][i1]=ddn[0][0][i1];
	 c[0][0][i1]=ddn[1][0][i1];
	 }
	 for (j1=N/2;j1>0;j1=j1/2) {
	 for (j2=N/2;j2>0;j2=j2/2) {
	 for (i1=0;i1<j1;i1++) {
	 for (i2=0;i2<j2;i2++) {
     c[0][j1+i1][j2+i2]=j2*ddn[0][j1+i1][j2+i2]+j1*ddn[1][j1+i1][j2+i2];
     c[1][j1+i1][j2+i2]=-j1*ddn[0][j1+i1][j2+i2]+j2*ddn[1][j1+i1][j2+i2];
	 }
	 }
	 }
	 } */
    fltarray WT;
 	WT.alloc(Data.nx(), Data.ny(), 2);
 	Data.alloc(Nc,Nl,2);
	
	if (Debug == 1)
    {
	WT(0,0,0) = Trans(0,0,0);
	WT(0,0,1) = Trans(0,0,1);
	}
	else {
		int bx = nbr_band_nx() -1;
		int by = nbr_band_ny() -1;
		for (int i1=0;i1<size_band_nx(bx);i1++) 
			for (int i2=0;i2<size_band_ny(by);i2++) 
			{
				int x = TabPosX(bx) + i1;
				int y = TabPosY(by) + i2;
				WT(x,y,0) = Trans(x,y,0);
				WT(x,y,1) = Trans(x,y,1);
			}
	}

	
	
	if (Debug == 1)
    {	for (int i1=1; i1 < Nc; i1++) 
	{
		WT(i1,0,1) = Trans(i1,0,0);
		WT(i1,0,0) = Trans(i1,0,1);
		WT(0,i1,0) = Trans(0,i1,0);
		WT(0,i1,1) = Trans(0,i1,1);
	}
	}

	if (Debug == 1)
    {
		// cout << "old" << endl;
		for (int j1=Nc/2;j1>0;j1=j1/2) 
		for (int j2=Nc/2;j2>0;j2=j2/2) 
		{
           // if (Verbose == True) cout << " j1 = " << j1 << ", j2 = " << j2 << ", nx = " << j1 << ", ny = " << j2 << endl;
			for (int i1=0;i1<j1;i1++) 
				for (int i2=0;i2<j2;i2++) 
				{
					WT(j1+i1,j2+i2,0) =  j2*Trans(j1+i1,j2+i2,0)  +  j1*Trans(j1+i1,j2+i2,1);
					WT(j1+i1,j2+i2,1) =  -j1*Trans(j1+i1,j2+i2,0)  +  j2*Trans(j1+i1,j2+i2,1);
				}
		}
	}
	else
    {
      //  cout << "new" << endl;
        int j1;
        int j2;
        for (int bx=0, j1=Nc/2; bx < nbr_band_nx(); bx++,j1=j1/2)
		for (int by=0, j2=Nl/2; by < nbr_band_ny(); by++,j2=j2/2)
		{
			// if (Verbose == True) cout << " j1 = " << j1 << ", j2 = " << j2 << ", nx = " << size_band_nx(bx) << ", ny = " << size_band_ny(by) << endl;
			// if (Verbose == True) cout << " bx = " << TabPosX(bx) << ", by = " << TabPosY(by) << endl;
            if ((bx !=  nbr_band_nx()-1) || (by != nbr_band_ny()-1))
		    {
               for (int i1=0;i1<size_band_nx(bx);i1++) 
               for (int i2=0;i2<size_band_ny(by);i2++) 
			   {
					int x = TabPosX(bx) + i1;
					int y = TabPosY(by) + i2;
					if ((x < 0) || (x >= WT.nx()) || (y < 0) || (y >= WT.ny()))
					{
					    cout << "Error " << x << " " << y << " " << WT.nx() << " " << WT.ny() << endl;
					    cout << TabPosX(bx) << " " <<  TabPosY(by) << " " << i1 << " " << i2 << endl;
					    exit(-1);
					}
					WT(x,y,0) =  j2*Trans(x,y,0)  +  j1*Trans(x,y,1);
					WT(x,y,1) = -j1*Trans(x,y,0)  +  j2*Trans(x,y,1);
			   }
			}
	    }
	}
	
	float *BuffData = WT.buffer();
	Ifloat Imag;
	Imag.alloc (BuffData,Nl,Nc);
	(WT_1->recons)(Imag, NbrScaleX, NbrScaleY);
	BuffData += Nl*Nc;
	Imag.alloc (BuffData,Nl,Nc);
	(WT_2->recons)(Imag, NbrScaleX, NbrScaleY);
	Data = WT;
}

/***************************************/

void MR_ANI_DIVCURL::curl_transform (fltarray &Data, fltarray &Trans)
{
    fltarray WT;
 	WT.alloc(Data.nx(), Data.ny(), 2);
 	Trans.alloc(Nc,Nl,2);
	
	for (int x=0; x < Data.nx(); x++)
		for (int y=0; y < Data.ny(); y++) 
		{
			WT(x,y,0) = Data(x,y,0);
			WT(x,y,1) = Data(x,y,1);
		}
    Ifloat Imag;
    float *BuffData = WT.buffer();
    Imag.alloc (BuffData,Nl,Nc);
	(WT_2->transform)(Imag, NbrScaleX, NbrScaleY);
	
    BuffData += Nl*Nc;
    Imag.alloc (BuffData,Nl,Nc);
    (WT_1->transform)(Imag, NbrScaleX, NbrScaleY);
	
	// fits_write_fltarr((char *) "a1cur.fits", WT);
	
	/* for (j1=N/2;j1>0;j1=j1/2) {
	 for (j2=N/2;j2>0;j2=j2/2) {
	 for (i1=0;i1<j1;i1++) {
	 for (i2=0;i2<j2;i2++) {
     dg[1][j1+i1][j2+i2]=(j1*d[0][j1+i1][j2+i2]+j2*d[1][j1+i1][j2+i2])/(j1*j1+j2*j2);
     dg[0][j1+i1][j2+i2]=(j2*d[0][j1+i1][j2+i2]-j1*d[1][j1+i1][j2+i2])/(j1*j1+j2*j2);
	 }
	 }
	 }
	 }
	 for (i1=1;i1<N;i1++) {
	 dg[0][i1][0]=d[0][i1][0];
	 dg[1][i1][0]=d[1][i1][0];
	 dg[0][0][i1]=d[1][0][i1];
	 dg[1][0][i1]=d[0][0][i1];
	 }
	 dg[0][0][0]=d[0][0][0];
	 dg[1][0][0]=d[1][0][0];
	 */
	if (Debug == 1)
    {
		// cout << "old" << endl;
		for (int j1=Nc/2;j1>0;j1=j1/2) 
		for (int j2=Nl/2;j2>0;j2=j2/2) 
		{
			for (int i1=0;i1<j1;i1++) 
				for (int i2=0;i2<j2;i2++) 
				{
					Trans(j1+i1,j2+i2,1)=(j1*WT(j1+i1,j2+i2,0)+j2*WT(j1+i1,j2+i2,1))/(j1*j1+j2*j2);
					Trans(j1+i1,j2+i2,0)=(j2*WT(j1+i1,j2+i2,0)-j1*WT(j1+i1,j2+i2,1))/(j1*j1+j2*j2);
				}
		}
	}
	else 
	{
       //  cout << "new" << endl;
        int j1;
        int j2;
        for (int bx=0, j1=Nc/2; bx < nbr_band_nx(); bx++,j1=j1/2)
		for (int by=0, j2=Nl/2; by < nbr_band_ny(); by++,j2=j2/2)
		{
			double Omega = 1. / (j1*j1+j2*j2);
			// if (Verbose == True) cout << " j1 = " << j1 << ", j2 = " << j2 << ", nx = " << size_band_nx(bx) << ", ny = " << size_band_ny(by) << endl;
			// if (Verbose == True) cout << " bx = " << TabPosX(bx) << ", by = " << TabPosY(by) << endl;
			if ((bx !=  nbr_band_nx()-1) || (by != nbr_band_ny()-1))
			{
				for (int i1=0;i1<size_band_nx(bx);i1++) 
				for (int i2=0;i2<size_band_ny(by);i2++) 
				{
					int x = TabPosX(bx) + i1;
					int y = TabPosY(by) + i2;
					Trans(x,y,1) = (j1 * WT(x,y,0) + j2 * WT(x,y,1)) / (j1*j1+j2*j2);
					Trans(x,y,0) = (j2 * WT(x,y,0) - j1 * WT(x,y,1)) / (j1*j1+j2*j2);
				}
			}
		}
		
	}

	if (Debug == 1)
    {	for (int i1=1; i1 < Nc; i1++) 
	{
		Trans(i1,0,1) =  WT(i1,0, 0);
		Trans(i1,0,0) = WT(i1,0,1);
		Trans(0,i1,1) = WT(0,i1,1);
		Trans(0,i1,0) = WT(0,i1,0);
	}
	}
	// cout << "Cur trans1 " << endl;
	
	
	if (Debug == 1)
    {
		Trans(0,0,0) = WT(0,0,0);
	    Trans(0,0,1) = WT(0,0,1);
	}
	else {
		int bx = nbr_band_nx() -1;
		int by = nbr_band_ny() -1;
		for (int i1=0;i1<size_band_nx(bx);i1++) 
			for (int i2=0;i2<size_band_ny(by);i2++) 
			{
				int x = TabPosX(bx) + i1;
				int y = TabPosY(by) + i2;
				Trans(x,y,0) =  WT(x,y,0);
				Trans(x,y,1) =  WT(x,y,1);
			}
	}
	
	// fits_write_fltarr((char *) "a1cure.fits", Trans);
	// cout << "Cur trans2 " << endl;
	
}


/***************************************/

void MR_ANI_DIVCURL::curl_recons (fltarray &Trans, fltarray &Data)
{
	// cout << "cur recons " << endl;
    fltarray WT;
 	WT.alloc(Data.nx(), Data.ny(), 2);
 	Data.alloc(Nc,Nl,2);
	
	// cout << "cur rec 1   " << endl;
 	if (Debug == 1)
    {
		WT(0,0,0) = Trans(0,0,0);
		WT(0,0,1) = Trans(0,0,1);
	}
	else {
		int bx = nbr_band_nx() -1;
		int by = nbr_band_ny() -1;
		for (int i1=0;i1<size_band_nx(bx);i1++) 
			for (int i2=0;i2<size_band_ny(by);i2++) 
			{
				int x = TabPosX(bx) + i1;
				int y = TabPosY(by) + i2;
				WT(x,y,0) = Trans(x,y,0);
				WT(x,y,1) = Trans(x,y,1);
			}
	}
	
	
 	
	if (Debug == 1)
    {
		// cout << "old" << endl;
		for (int i1=1; i1 < Nc; i1++) 
	{
		WT(i1,0,1) = Trans(i1,0,0);
		WT(i1,0,0) = Trans(i1,0,1);
		WT(0,i1,0) = Trans(0,i1,0);
		WT(0,i1,1) = Trans(0,i1,1);
	}
	for (int j1=Nc/2;j1>0;j1=j1/2) 
		for (int j2=Nc/2;j2>0;j2=j2/2) 
		{
			for (int i1=0;i1<j1;i1++) 
				for (int i2=0;i2<j2;i2++) 
				{
					WT(j1+i1,j2+i2,0) =  j1*Trans(j1+i1,j2+i2,1)  +  j2*Trans(j1+i1,j2+i2,0);
					WT(j1+i1,j2+i2,1) =  j2*Trans(j1+i1,j2+i2,1)  -  j1*Trans(j1+i1,j2+i2,0);
				}
		}
	}
	else
    {
      //   cout << "new" << endl;
        int j1;
        int j2;
        for (int bx=0, j1=Nc/2; bx < nbr_band_nx(); bx++,j1=j1/2)
			for (int by=0, j2=Nl/2; by < nbr_band_ny(); by++,j2=j2/2)
			{
				// if (Verbose == True) cout << " j1 = " << j1 << ", j2 = " << j2 << ", nx = " << size_band_nx(bx) << ", ny = " << size_band_ny(by) << endl;
				// if (Verbose == True) cout << " bx = " << TabPosX(bx) << ", by = " << TabPosY(by) << endl;
				if ((bx !=  nbr_band_nx()-1) || (by != nbr_band_ny()-1))
				{
					for (int i1=0;i1<size_band_nx(bx);i1++) 
						for (int i2=0;i2<size_band_ny(by);i2++) 
						{
							int x = TabPosX(bx) + i1;
							int y = TabPosY(by) + i2;
							WT(x,y,0) =  j1*Trans(x,y,1)  +  j2*Trans(x,y,0);
							WT(x,y,1) =  j2*Trans(x,y,1)  -  j1*Trans(x,y,0);
						}
				}
			}
	}
	// cout << "cur wt rec " << endl;
	
	float *BuffData = WT.buffer();
	Ifloat Imag;
	Imag.alloc (BuffData,Nl,Nc);
	(WT_2->recons)(Imag, NbrScaleX, NbrScaleY);
	BuffData += Nl*Nc;
	Imag.alloc (BuffData,Nl,Nc);
	(WT_1->recons)(Imag, NbrScaleX, NbrScaleY);
	Data = WT;
}

/*%%%%%%%%%%%%%%%%%%%%%%% Erwan %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

/*********************** en periodique ************************************

void MR_ISO_DIVCURL::div_quasi_interpol(fltarray & Shear, fltarray & Ima2D)
{
 fltarray ponder;
 ponder.resize(nx(),ny(),2);
 Ima2D.resize(nx(),ny(),2);
 for (int e=0;e<2;e++) for (int i2=0;i2<ny();i2++) for (int i1=0;i1<nx();i1++) { Ima2D(i1,i2,e)=0.; ponder(i1,i2,e)=0.;}
 for (int i=0;i<np();i++) {
  i1 = posx(i); i2 = posy(i);
  float x1 = nx()*(TabX(i) -MinX) / (MaxX-MinX); float x2 = ny()*(TabY(i) -MinY) / (MaxY-MinY);
  
  /* dans V1 x V0 *
  Ima2D((i1+1)%nx(),i2,0)+=0.5*(x1-i1)*(x1-i1)*(1-fabs(x2-i2))*Shear(i,0);
  Ima2D((i1+1)%nx(),(i2+1)%ny(),0)+=0.5*(x1-i1)*(x1-i1)*(1-fabs(x2-i2-1))*Shear(i,0);
  Ima2D(i1,i2,0)+=(0.75-(x1-i1-0.5)*(x1-i1-0.5))*(1-fabs(x2-i2))*Shear(i,0);
  Ima2D(i1,(i2+1)%ny(),0)+=(0.75-(x1-i1-0.5)*(x1-i1-0.5))*(1-fabs(x2-i2-1))*Shear(i,0);
  Ima2D((i1-1+nx())%nx(),i2,0)+=0.5*(x1-i1-1)*(x1-i1-1)*(1-fabs(x2-i2))*Shear(i,0);
  Ima2D((i1-1+nx())%nx(),(i2+1)%ny(),0)+=0.5*(x1-i1-1)*(x1-i1-1)*(1-fabs(x2-i2-1))*Shear(i,0);
  
  ponder((i1+1)%nx(),i2,0)+=0.5*(x1-i1)*(x1-i1)*(1-fabs(x2-i2));
  ponder((i1+1)%nx(),(i2+1)%ny(),0)+=0.5*(x1-i1)*(x1-i1)*(1-fabs(x2-i2-1));
  ponder(i1,i2,0)+=(0.75-(x1-i1-0.5)*(x1-i1-0.5))*(1-fabs(x2-i2));
  ponder(i1,(i2+1)%ny(),0)+=(0.75-(x1-i1-0.5)*(x1-i1-0.5))*(1-fabs(x2-i2-1));
  ponder((i1-1+nx())%nx(),i2,0)+=0.5*(x1-i1-1)*(x1-i1-1)*(1-fabs(x2-i2));
  ponder((i1-1+nx())%nx(),(i2+1)%ny(),0)+=0.5*(x1-i1-1)*(x1-i1-1)*(1-fabs(x2-i2-1));
  
  /* dans V0 x V1 *
  Ima2D(i1,(i2+1)%ny(),1)+=0.5*(x2-i2)*(x2-i2)*(1-fabs(x1-i1))*Shear(i,1);
  Ima2D((i1+1)%nx(),(i2+1)%ny(),1)+=0.5*(x2-i2)*(x2-i2)*(1-fabs(x1-i1-1))*Shear(i,1);
  Ima2D(i1,i2,1)+=(0.75-(x2-i2-0.5)*(x2-i2-0.5))*(1-fabs(x1-i1))*Shear(i,1);
  Ima2D((i1+1)%nx(),i2,1)+=(0.75-(x2-i2-0.5)*(x2-i2-0.5))*(1-fabs(x1-i1-1))*Shear(i,1);
  Ima2D(i1,(i2-1+ny())%ny(),1)+=0.5*(x2-i2-1)*(x2-i2-1)*(1-fabs(x1-i1))*Shear(i,1);
  Ima2D((i1+1)%nx(),(i2-1+ny())%ny(),1)+=0.5*(x2-i2-1)*(x2-i2-1)*(1-fabs(x1-i1-1))*Shear(i,1);
  
  ponder(i1,(i2+1)%ny(),1)+=0.5*(x2-i2)*(x2-i2)*(1-fabs(x1-i1));
  ponder((i1+1)%nx(),(i2+1)%ny(),1)+=0.5*(x2-i2)*(x2-i2)*(1-fabs(x1-i1-1));
  ponder(i1,i2,1)+=(0.75-(x2-i2-0.5)*(x2-i2-0.5))*(1-fabs(x1-i1));
  ponder((i1+1)%nx(),i2,1)+=(0.75-(x2-i2-0.5)*(x2-i2-0.5))*(1-fabs(x1-i1-1));
  ponder(i1,(i2-1+ny())%ny(),1)+=0.5*(x2-i2-1)*(x2-i2-1)*(1-fabs(x1-i1));
  ponder((i1+1)%nx(),(i2-1+ny())%ny(),1)+=0.5*(x2-i2-1)*(x2-i2-1)*(1-fabs(x1-i1-1));
 }
 for (int e=0;e<2;e++) for (int i2=0;i2<ny();i2++) for (int i1=0;i1<nx();i1++) if (ponder(i1,i2,e)==0.) ponder(i1,i2,e)=1.;
 for (int e=0;e<2;e++) for (int i2=0;i2<ny();i2++) for (int i1=0;i1<nx();i1++) Ima2D(i1,i2,e)=Ima2D(i1,i2,e)/ponder(i1,i2,e);
}

void MR_ISO_DIVCURL::curl_quasi_interpol(fltarray & Shear, fltarray & Ima2D)
{
 fltarray ponder;
 ponder.resize(nx(),ny(),2);
 Ima2D.resize(nx(),ny(),2);
 for (int e=0;e<2;e++) for (int i2=0;i2<ny();i2++) for (int i1=0;i1<nx();i1++) { Ima2D(i1,i2,e)=0.; ponder(i1,i2,e)=0.;}
 for (int i=0;i<np();i++) {
  i1 = posx(i); i2 = posy(i);
  float x1 = nx()*(TabX(i) -MinX) / (MaxX-MinX); float x2 = ny()*(TabY(i) -MinY) / (MaxY-MinY);
  
  /* dans V0 x V1 *
  Ima2D(i1,(i2+1)%ny(),0)+=0.5*(x2-i2)*(x2-i2)*(1-fabs(x1-i1))*Shear(i,0);
  Ima2D((i1+1)%nx(),(i2+1)%ny(),0)+=0.5*(x2-i2)*(x2-i2)*(1-fabs(x1-i1-1))*Shear(i,0);
  Ima2D(i1,i2,0)+=(0.75-(x2-i2-0.5)*(x2-i2-0.5))*(1-fabs(x1-i1))*Shear(i,0);
  Ima2D((i1+1)%nx(),i2,0)+=(0.75-(x2-i2-0.5)*(x2-i2-0.5))*(1-fabs(x1-i1-1))*Shear(i,0);
  Ima2D(i1,(i2-1+ny())%ny(),0)+=0.5*(x2-i2-1)*(x2-i2-1)*(1-fabs(x1-i1))*Shear(i,0);
  Ima2D((i1+1)%nx(),(i2-1+ny())%ny(),0)+=0.5*(x2-i2-1)*(x2-i2-1)*(1-fabs(x1-i1-1))*Shear(i,0);
  
  ponder(i1,(i2+1)%ny(),0)+=0.5*(x2-i2)*(x2-i2)*(1-fabs(x1-i1));
  ponder((i1+1)%nx(),(i2+1)%ny(),0)+=0.5*(x2-i2)*(x2-i2)*(1-fabs(x1-i1-1));
  ponder(i1,i2,0)+=(0.75-(x2-i2-0.5)*(x2-i2-0.5))*(1-fabs(x1-i1));
  ponder((i1+1)%nx(),i2,0)+=(0.75-(x2-i2-0.5)*(x2-i2-0.5))*(1-fabs(x1-i1-1));
  ponder(i1,(i2-1+ny())%ny(),0)+=0.5*(x2-i2-1)*(x2-i2-1)*(1-fabs(x1-i1));
  ponder((i1+1)%nx(),(i2-1+ny())%ny(),0)+=0.5*(x2-i2-1)*(x2-i2-1)*(1-fabs(x1-i1-1));
  
  /* dans V1 x V0 *
  Ima2D((i1+1)%nx(),i2,1)+=0.5*(x1-i1)*(x1-i1)*(1-fabs(x2-i2))*Shear(i,1);
  Ima2D((i1+1)%nx(),(i2+1)%ny(),1)+=0.5*(x1-i1)*(x1-i1)*(1-fabs(x2-i2-1))*Shear(i,1);
  Ima2D(i1,i2,1)+=(0.75-(x1-i1-0.5)*(x1-i1-0.5))*(1-fabs(x2-i2))*Shear(i,1);
  Ima2D(i1,(i2+1)%ny(),1)+=(0.75-(x1-i1-0.5)*(x1-i1-0.5))*(1-fabs(x2-i2-1))*Shear(i,1);
  Ima2D((i1-1+nx())%nx(),i2,1)+=0.5*(x1-i1-1)*(x1-i1-1)*(1-fabs(x2-i2))*Shear(i,1);
  Ima2D((i1-1+nx())%nx(),(i2+1)%ny(),1)+=0.5*(x1-i1-1)*(x1-i1-1)*(1-fabs(x2-i2-1))*Shear(i,1);
  
  ponder((i1+1)%nx(),i2,1)+=0.5*(x1-i1)*(x1-i1)*(1-fabs(x2-i2));
  ponder((i1+1)%nx(),(i2+1)%ny(),1)+=0.5*(x1-i1)*(x1-i1)*(1-fabs(x2-i2-1));
  ponder(i1,i2,1)+=(0.75-(x1-i1-0.5)*(x1-i1-0.5))*(1-fabs(x2-i2));
  ponder(i1,(i2+1)%ny(),1)+=(0.75-(x1-i1-0.5)*(x1-i1-0.5))*(1-fabs(x2-i2-1));
  ponder((i1-1+nx())%nx(),i2,1)+=0.5*(x1-i1-1)*(x1-i1-1)*(1-fabs(x2-i2));
  ponder((i1-1+nx())%nx(),(i2+1)%ny(),1)+=0.5*(x1-i1-1)*(x1-i1-1)*(1-fabs(x2-i2-1));
 }
 for (int e=0;e<2;e++) for (int i2=0;i2<ny();i2++) for (int i1=0;i1<nx();i1++) if (ponder(i1,i2,e)==0.) ponder(i1,i2,e)=1.;
 for (int e=0;e<2;e++) for (int i2=0;i2<ny();i2++) for (int i1=0;i1<nx();i1++) Ima2D(i1,i2,e)=Ima2D(i1,i2,e)/ponder(i1,i2,e);
}

void MR_ISO_DIVCURL::div_point_value(fltarray & Ima2D, fltarray & Shear)
{
 Shear.resize(nx(),ny(),2);
 for (int e=0;e<2;e++) for (int i2=0;i2<ny();i2++) for (int i1=0;i1<nx();i1++) { Shear(i1,i2,e)=0.;}
 for (int i=0;i<np();i++) {
  i1 = posx(i); i2 = posy(i);
  float x1 = nx()*(TabX(i) -MinX) / (MaxX-MinX); float x2 = ny()*(TabY(i) -MinY) / (MaxY-MinY);
  
  /* dans V1 x V0 *
  Shear(i,0)+=Ima2D((i1+1)%nx(),i2,0)*0.5*(x1-i1)*(x1-i1)*(1-fabs(x2-i2));
  Shear(i,0)+=Ima2D((i1+1)%nx(),(i2+1)%ny(),0)*0.5*(x1-i1)*(x1-i1)*(1-fabs(x2-i2-1));
  Shear(i,0)+=Ima2D(i1,i2,0)*(0.75-(x1-i1-0.5)*(x1-i1-0.5))*(1-fabs(x2-i2));
  Shear(i,0)+=Ima2D(i1,(i2+1)%ny(),0)*(0.75-(x1-i1-0.5)*(x1-i1-0.5))*(1-fabs(x2-i2-1));
  Shear(i,0)+=Ima2D((i1-1+nx())%nx(),i2,0)*0.5*(x1-i1-1)*(x1-i1-1)*(1-fabs(x2-i2));
  Shear(i,0)+=Ima2D((i1-1+nx())%nx(),(i2+1)%ny(),0)*0.5*(x1-i1-1)*(x1-i1-1)*(1-fabs(x2-i2-1));
  
  /* dans V0 x V1 *
  Shear(i,1)+=Ima2D(i1,(i2+1)%ny(),1)*0.5*(x2-i2)*(x2-i2)*(1-fabs(x1-i1));
  Shear(i,1)+=Ima2D((i1+1)%nx(),(i2+1)%ny(),1)*0.5*(x2-i2)*(x2-i2)*(1-fabs(x1-i1-1));
  Shear(i,1)+=Ima2D(i1,i2,1)*(0.75-(x2-i2-0.5)*(x2-i2-0.5))*(1-fabs(x1-i1));
  Shear(i,1)+=Ima2D((i1+1)%nx(),i2,1)*(0.75-(x2-i2-0.5)*(x2-i2-0.5))*(1-fabs(x1-i1-1));
  Shear(i,1)+=Ima2D(i1,(i2-1+ny())%ny(),1)*0.5*(x2-i2-1)*(x2-i2-1)*(1-fabs(x1-i1));
  Shear(i,1)+=Ima2D((i1+1)%nx(),(i2-1+ny())%ny(),1)*0.5*(x2-i2-1)*(x2-i2-1)*(1-fabs(x1-i1-1));
 }
}

void MR_ISO_DIVCURL::curl_point_value(fltarray & Ima2D, fltarray & Shear)
{
 Shear.resize(nx(),ny(),2);
 for (int e=0;e<2;e++) for (int i2=0;i2<ny();i2++) for (int i1=0;i1<nx();i1++) { Shear(i1,i2,e)=0.;}
 for (int i=0;i<np();i++) {
  i1 = posx(i); i2 = posy(i);
  float x1 = nx()*(TabX(i) -MinX) / (MaxX-MinX); float x2 = ny()*(TabY(i) -MinY) / (MaxY-MinY);
  
  /* dans V0 x V1 *
  Shear(i,0)+=Ima2D(i1,(i2+1)%ny(),0)*0.5*(x2-i2)*(x2-i2)*(1-fabs(x1-i1));
  Shear(i,0)+=Ima2D((i1+1)%nx(),(i2+1)%ny(),0)*0.5*(x2-i2)*(x2-i2)*(1-fabs(x1-i1-1));
  Shear(i,0)+=Ima2D(i1,i2,0)*(0.75-(x2-i2-0.5)*(x2-i2-0.5))*(1-fabs(x1-i1));
  Shear(i,0)+=Ima2D((i1+1)%nx(),i2,0)*(0.75-(x2-i2-0.5)*(x2-i2-0.5))*(1-fabs(x1-i1-1));
  Shear(i,0)+=Ima2D(i1,(i2-1+ny())%ny(),0)*0.5*(x2-i2-1)*(x2-i2-1)*(1-fabs(x1-i1));
  Shear(i,0)+=Ima2D((i1+1)%nx(),(i2-1+ny())%ny(),0)*0.5*(x2-i2-1)*(x2-i2-1)*(1-fabs(x1-i1-1));
  
  /* dans V1 x V0 *
  Shear(i,1)+=Ima2D((i1+1)%nx(),i2,1)*0.5*(x1-i1)*(x1-i1)*(1-fabs(x2-i2));
  Shear(i,1)+=Ima2D((i1+1)%nx(),(i2+1)%ny(),1)*0.5*(x1-i1)*(x1-i1)*(1-fabs(x2-i2-1));
  Shear(i,1)+=Ima2D(i1,i2,1)*(0.75-(x1-i1-0.5)*(x1-i1-0.5))*(1-fabs(x2-i2));
  Shear(i,1)+=Ima2D(i1,(i2+1)%ny(),1)*(0.75-(x1-i1-0.5)*(x1-i1-0.5))*(1-fabs(x2-i2-1));
  Shear(i,1)+=Ima2D((i1-1+nx())%nx(),i2,1)*0.5*(x1-i1-1)*(x1-i1-1)*(1-fabs(x2-i2));
  Shear(i,1)+=Ima2D((i1-1+nx())%nx(),(i2+1)%ny(),1)*0.5*(x1-i1-1)*(x1-i1-1)*(1-fabs(x2-i2-1));
 }
}


/**************************************************************************************/
/*********************** avec conditions aux bords ************************************/
/**************************************************************************************/


void HelmCatalog::div_quasi_interpol(fltarray & Shear, fltarray & Ima2D)
{
 fltarray ponder;
 ponder.resize(nx(),ny(),2);
 Ima2D.resize(nx(),ny(),2);
 
 for (int e=0;e<2;e++) for (int i2=0;i2<ny();i2++) for (int i1=0;i1<nx();i1++) { Ima2D(i1,i2,e)=0.; ponder(i1,i2,e)=0.;}
 for (int i=0;i<np();i++) {
  int i1 = posx(i); int i2 = posy(i);
  float x1 =(TabX(i)-MinX)/BinCat;
  float x2 =(TabY(i)-MinY)/BinCat;

  /* dans V1 x V0 */
  Ima2D.setxyz(i1+1,i2,0,DEF_DIVROT_BORD)   +=0.5*(x1-i1)*(x1-i1)*(1-fabs(x2-i2))*Shear(i,0);
  Ima2D.setxyz(i1+1,i2+1,0,DEF_DIVROT_BORD) +=0.5*(x1-i1)*(x1-i1)*(1-fabs(x2-i2-1))*Shear(i,0);
  Ima2D.setxyz(i1,i2,0,DEF_DIVROT_BORD)     +=(0.75-(x1-i1-0.5)*(x1-i1-0.5))*(1-fabs(x2-i2))*Shear(i,0);
  Ima2D.setxyz(i1,i2+1,0,DEF_DIVROT_BORD)   +=(0.75-(x1-i1-0.5)*(x1-i1-0.5))*(1-fabs(x2-i2-1))*Shear(i,0);
  Ima2D.setxyz(i1-1,i2,0,DEF_DIVROT_BORD)   +=0.5*(x1-i1-1)*(x1-i1-1)*(1-fabs(x2-i2))*Shear(i,0);
  Ima2D.setxyz(i1-1,i2+1,0,DEF_DIVROT_BORD) +=0.5*(x1-i1-1)*(x1-i1-1)*(1-fabs(x2-i2-1))*Shear(i,0);
  
  ponder.setxyz(i1+1,i2,0,DEF_DIVROT_BORD)   +=0.5*(x1-i1)*(x1-i1)*(1-fabs(x2-i2));
  ponder.setxyz(i1+1,i2+1,0,DEF_DIVROT_BORD) +=0.5*(x1-i1)*(x1-i1)*(1-fabs(x2-i2-1));
  ponder.setxyz(i1,i2,0,DEF_DIVROT_BORD)     +=(0.75-(x1-i1-0.5)*(x1-i1-0.5))*(1-fabs(x2-i2));
  ponder.setxyz(i1,i2+1,0,DEF_DIVROT_BORD)   +=(0.75-(x1-i1-0.5)*(x1-i1-0.5))*(1-fabs(x2-i2-1));
  ponder.setxyz(i1-1,i2,0,DEF_DIVROT_BORD)   +=0.5*(x1-i1-1)*(x1-i1-1)*(1-fabs(x2-i2));
  ponder.setxyz(i1-1,i2+1,0,DEF_DIVROT_BORD) +=0.5*(x1-i1-1)*(x1-i1-1)*(1-fabs(x2-i2-1));
  
  /* dans V0 x V1 */
  Ima2D.setxyz(i1,i2+1,1,DEF_DIVROT_BORD)    +=0.5*(x2-i2)*(x2-i2)*(1-fabs(x1-i1))*Shear(i,1);
  Ima2D.setxyz(i1+1,i2+1,1,DEF_DIVROT_BORD)  +=0.5*(x2-i2)*(x2-i2)*(1-fabs(x1-i1-1))*Shear(i,1);
  Ima2D.setxyz(i1,i2,1,DEF_DIVROT_BORD)      +=(0.75-(x2-i2-0.5)*(x2-i2-0.5))*(1-fabs(x1-i1))*Shear(i,1);
  Ima2D.setxyz(i1+1,i2,1,DEF_DIVROT_BORD)    +=(0.75-(x2-i2-0.5)*(x2-i2-0.5))*(1-fabs(x1-i1-1))*Shear(i,1);
  Ima2D.setxyz(i1,i2-1,1,DEF_DIVROT_BORD)    +=0.5*(x2-i2-1)*(x2-i2-1)*(1-fabs(x1-i1))*Shear(i,1);
  Ima2D.setxyz(i1+1,i2-1,1,DEF_DIVROT_BORD)  +=0.5*(x2-i2-1)*(x2-i2-1)*(1-fabs(x1-i1-1))*Shear(i,1);
  
  ponder.setxyz(i1,i2+1,1,DEF_DIVROT_BORD)+=0.5*(x2-i2)*(x2-i2)*(1-fabs(x1-i1));
  ponder.setxyz(i1+1,i2+1,1,DEF_DIVROT_BORD)+=0.5*(x2-i2)*(x2-i2)*(1-fabs(x1-i1-1));
  ponder.setxyz(i1,i2,1,DEF_DIVROT_BORD)+=(0.75-(x2-i2-0.5)*(x2-i2-0.5))*(1-fabs(x1-i1));
  ponder.setxyz(i1+1,i2,1,DEF_DIVROT_BORD)+=(0.75-(x2-i2-0.5)*(x2-i2-0.5))*(1-fabs(x1-i1-1));
  ponder.setxyz(i1,i2-1,1,DEF_DIVROT_BORD)+=0.5*(x2-i2-1)*(x2-i2-1)*(1-fabs(x1-i1));
  ponder.setxyz(i1+1,i2-1,1,DEF_DIVROT_BORD)+=0.5*(x2-i2-1)*(x2-i2-1)*(1-fabs(x1-i1-1));
 }
 for (int e=0;e<2;e++) for (int i2=0;i2<ny();i2++) for (int i1=0;i1<nx();i1++) if (ponder(i1,i2,e)==0.) ponder(i1,i2,e)=1.;
 for (int e=0;e<2;e++) for (int i2=0;i2<ny();i2++) for (int i1=0;i1<nx();i1++) Ima2D(i1,i2,e)=Ima2D(i1,i2,e)/ponder(i1,i2,e);
}


void HelmCatalog::div_bord_quasi_interpol(fltarray & Shear, fltarray & Ima2D,int N)
{
    float phi1,phi2;
    fltarray ponder;
    ponder.resize(N+2,N+2,2);
    Ima2D.resize(N+2,N+2,2);
     
   // float BinCatX=(maxX-minX)/N;  float BinCatY=(maxY-minY)/N;
    float BinCatX=BinCat;
    float BinCatY=BinCat;
    
    for (int e=0;e<2;e++) 
    for (int i2=0;i2<N+2;i2++) 
    for (int i1=0;i1<N+2;i1++) 
    { 
        Ima2D(i1,i2,e)=0.; 
        ponder(i1,i2,e)=0.;
    }
   
     for (int i=0;i<np();i++) 
     {
        float x1 =(TabX(i)-MinX)/BinCat;
        float x2 =(TabY(i)-MinY)/BinCat;
        int i1=(int) x1; 
        int i2=(int) x2; 
        i1-=(i1==N); 
        i2-=(i2==N);
        
        /* dans V1 x V0 */
        for (int e1=0;e1<3;e1++) 
        {
            if (i1==0) phi1=(e1==0)*(1.-x1)*(1.-x1)+(e1==1)*x1*(2.-3./2.*x1)+(e1==2)*x1*x1/2.;
            else if (i1>=(N-1)) 
            { 
               float X=N-x1; 
               phi1=(e1==0)*X*X/2.+(e1==1)*X*(2.-3./2.*X)+(e1==2)*(1.-X)*(1.-X);
            }
            else phi1=3./4.*(e1==1)+(1./2.-3./2.*(e1==1))*(x1-i1-(2.-e1)/2.)*(x1-i1-(2.-e1)/2.);
            for (int e2=0;e2<2;e2++) 
            {
                phi2=1-fabs(x2-i2-e2);
                Ima2D(i1+e1,i2+e2,0)   += phi1*phi2*Shear(i,0);
                ponder(i1+e1,i2+e2,0)  += phi1*phi2;
            }
        }
        
        /* dans V0 x V1 */
        for (int e1=0;e1<2;e1++) 
        {
            phi1=1-fabs(x1-i1-e1);
            for (int e2=0;e2<3;e2++) 
            {
                if (i2==0) phi2=(e2==0)*(1.-x2)*(1.-x2)+(e2==1)*x2*(2.-3./2.*x2)+(e2==2)*x2*x2/2.;
                else if (i2>=(N-1)) 
                { 
                   float X=N-x2; phi2=(e2==0)*X*X/2.+(e2==1)*X*(2.-3./2.*X)+(e2==2)*(1.-X)*(1.-X);
                }
                else phi2=3./4.*(e2==1)+(1./2.-3./2.*(e2==1))*(x2-i2-(2.-e2)/2.)*(x2-i2-(2.-e2)/2.);
                
                Ima2D(i1+e1,i2+e2,1)   += phi1*phi2*Shear(i,1);
                ponder(i1+e1,i2+e2,1)  += phi1*phi2;
            }
        }
    }
    for (int e=0;e<2;e++) for (int i2=0;i2<N+2;i2++) for (int i1=0;i1<N+2;i1++) Ima2D(i1,i2,e)=Ima2D(i1,i2,e)/(ponder(i1,i2,e)+(ponder(i1,i2,e)==0.));
}






void HelmCatalog::curl_quasi_interpol(fltarray & Shear, fltarray & Ima2D)
{
 fltarray ponder;
 ponder.resize(nx(),ny(),2);
 Ima2D.resize(nx(),ny(),2);
 for (int e=0;e<2;e++) for (int i2=0;i2<ny();i2++) for (int i1=0;i1<nx();i1++) { Ima2D(i1,i2,e)=0.; ponder(i1,i2,e)=0.;}
 for (int i=0;i<np();i++) {
  int i1 = posx(i); int i2 = posy(i);
  float x1 =(TabX(i)-MinX)/BinCat;
  float x2 =(TabY(i)-MinY)/BinCat;
 
  /* dans V0 x V1 */
  Ima2D.setxyz(i1,i2+1,0,DEF_DIVROT_BORD)+=0.5*(x2-i2)*(x2-i2)*(1-fabs(x1-i1))*Shear(i,0);
  Ima2D.setxyz(i1+1,i2+1,0,DEF_DIVROT_BORD)+=0.5*(x2-i2)*(x2-i2)*(1-fabs(x1-i1-1))*Shear(i,0);
  Ima2D.setxyz(i1,i2,0,DEF_DIVROT_BORD)+=(0.75-(x2-i2-0.5)*(x2-i2-0.5))*(1-fabs(x1-i1))*Shear(i,0);
  Ima2D.setxyz(i1+1,i2,0,DEF_DIVROT_BORD)+=(0.75-(x2-i2-0.5)*(x2-i2-0.5))*(1-fabs(x1-i1-1))*Shear(i,0);
  Ima2D.setxyz(i1,i2-1,0,DEF_DIVROT_BORD)+=0.5*(x2-i2-1)*(x2-i2-1)*(1-fabs(x1-i1))*Shear(i,0);
  Ima2D.setxyz(i1+1,i2-1,0,DEF_DIVROT_BORD)+=0.5*(x2-i2-1)*(x2-i2-1)*(1-fabs(x1-i1-1))*Shear(i,0);
  
  ponder.setxyz(i1,i2+1,0,DEF_DIVROT_BORD)+=0.5*(x2-i2)*(x2-i2)*(1-fabs(x1-i1));
  ponder.setxyz(i1+1,i2+1,0,DEF_DIVROT_BORD)+=0.5*(x2-i2)*(x2-i2)*(1-fabs(x1-i1-1));
  ponder.setxyz(i1,i2,0,DEF_DIVROT_BORD)+=(0.75-(x2-i2-0.5)*(x2-i2-0.5))*(1-fabs(x1-i1));
  ponder.setxyz(i1+1,i2,0,DEF_DIVROT_BORD)+=(0.75-(x2-i2-0.5)*(x2-i2-0.5))*(1-fabs(x1-i1-1));
  ponder.setxyz(i1,i2-1,0,DEF_DIVROT_BORD)+=0.5*(x2-i2-1)*(x2-i2-1)*(1-fabs(x1-i1));
  ponder.setxyz(i1+1,i2-1,0,DEF_DIVROT_BORD)+=0.5*(x2-i2-1)*(x2-i2-1)*(1-fabs(x1-i1-1));
  
  /* dans V1 x V0 */
  Ima2D.setxyz(i1+1,i2,1,DEF_DIVROT_BORD)+=0.5*(x1-i1)*(x1-i1)*(1-fabs(x2-i2))*Shear(i,1);
  Ima2D.setxyz(i1+1,i2+1,1,DEF_DIVROT_BORD)+=0.5*(x1-i1)*(x1-i1)*(1-fabs(x2-i2-1))*Shear(i,1);
  Ima2D.setxyz(i1,i2,1,DEF_DIVROT_BORD)+=(0.75-(x1-i1-0.5)*(x1-i1-0.5))*(1-fabs(x2-i2))*Shear(i,1);
  Ima2D.setxyz(i1,i2+1,1,DEF_DIVROT_BORD)+=(0.75-(x1-i1-0.5)*(x1-i1-0.5))*(1-fabs(x2-i2-1))*Shear(i,1);
  Ima2D.setxyz(i1-1,i2,1,DEF_DIVROT_BORD)+=0.5*(x1-i1-1)*(x1-i1-1)*(1-fabs(x2-i2))*Shear(i,1);
  Ima2D.setxyz(i1-1,i2+1,1,DEF_DIVROT_BORD)+=0.5*(x1-i1-1)*(x1-i1-1)*(1-fabs(x2-i2-1))*Shear(i,1);
  
  ponder.setxyz(i1+1,i2,1,DEF_DIVROT_BORD)+=0.5*(x1-i1)*(x1-i1)*(1-fabs(x2-i2));
  ponder.setxyz(i1+1,i2+1,1,DEF_DIVROT_BORD)+=0.5*(x1-i1)*(x1-i1)*(1-fabs(x2-i2-1));
  ponder.setxyz(i1,i2,1,DEF_DIVROT_BORD)+=(0.75-(x1-i1-0.5)*(x1-i1-0.5))*(1-fabs(x2-i2));
  ponder.setxyz(i1,i2+1,1,DEF_DIVROT_BORD)+=(0.75-(x1-i1-0.5)*(x1-i1-0.5))*(1-fabs(x2-i2-1));
  ponder.setxyz(i1-1,i2,1,DEF_DIVROT_BORD)+=0.5*(x1-i1-1)*(x1-i1-1)*(1-fabs(x2-i2));
  ponder.setxyz(i1-1,i2+1,1,DEF_DIVROT_BORD)+=0.5*(x1-i1-1)*(x1-i1-1)*(1-fabs(x2-i2-1));
 }
 for (int e=0;e<2;e++) for (int i2=0;i2<ny();i2++) for (int i1=0;i1<nx();i1++) if (ponder(i1,i2,e)==0.) ponder(i1,i2,e)=1.;
 for (int e=0;e<2;e++) for (int i2=0;i2<ny();i2++) for (int i1=0;i1<nx();i1++) Ima2D(i1,i2,e)=Ima2D(i1,i2,e)/ponder(i1,i2,e);
}


void HelmCatalog::curl_bord_quasi_interpol(fltarray & Shear, fltarray & Ima2D,int N)
{
    float phi1,phi2;
    fltarray ponder;
    ponder.resize(N+2,N+2,2);
    Ima2D.resize(N+2,N+2,2);
    for (int i1=0;i1<N+2;i1++) 
    for (int i2=0;i2<N+2;i2++) Ima2D(i1,i2,0) = Ima2D(i1,i2,1) = 0;
    
    // float BinCatX=(maxX-minX)/N;  float BinCatY=(maxY-minY)/N;
    float BinCatX=BinCat;
    float BinCatY=BinCat;
    
    for (int e=0;e<2;e++) 
    for (int i2=0;i2<N+2;i2++) 
    for (int i1=0;i1<N+2;i1++) 
    { 
        Ima2D(i1,i2,e)=0.; 
        ponder(i1,i2,e)=0.;
    }
    
    for (int i=0;i<np();i++) 
    {
        float x1 =(TabX(i)-MinX)/BinCat;
        float x2 =(TabY(i)-MinY)/BinCat;
        int i1=(int) x1; int i2=(int) x2; i1-=(i1==N); i2-=(i2==N);
        
        /* dans V0 x V1 */
        for (int e1=0;e1<2;e1++) 
        {
            phi1=1-fabs(x1-i1-e1);
            for (int e2=0;e2<3;e2++) 
            {
                if (i2==0) phi2=(e2==0)*(1.-x2)*(1.-x2)+(e2==1)*x2*(2.-3./2.*x2)+(e2==2)*x2*x2/2.;
                else if (i2>=(N-1)) 
                { 
                   float X=N-x2; 
                   phi2=(e2==0)*X*X/2.+(e2==1)*X*(2.-3./2.*X)+(e2==2)*(1.-X)*(1.-X);
                }
                else phi2=3./4.*(e2==1)+(1./2.-3./2.*(e2==1))*(x2-i2-(2.-e2)/2.)*(x2-i2-(2.-e2)/2.);
                
                Ima2D(i1+e1,i2+e2,0)   += phi1*phi2*Shear(i,0);
                ponder(i1+e1,i2+e2,0)  += phi1*phi2;
            }
        }
        
        /* dans V1 x V0 */
        for (int e1=0;e1<3;e1++) 
        {
            if (i1==0) phi1=(e1==0)*(1.-x1)*(1.-x1)+(e1==1)*x1*(2.-3./2.*x1)+(e1==2)*x1*x1/2.;
            else if (i1>=(N-1)) 
            { 
               float X=N-x1; 
               phi1=(e1==0)*X*X/2.+(e1==1)*X*(2.-3./2.*X)+(e1==2)*(1.-X)*(1.-X);
            }
            else phi1=3./4.*(e1==1)+(1./2.-3./2.*(e1==1))*(x1-i1-(2.-e1)/2.)*(x1-i1-(2.-e1)/2.);
            for (int e2=0;e2<2;e2++) 
            {
                phi2=1-fabs(x2-i2-e2);
                
                Ima2D(i1+e1,i2+e2,1)   += phi1*phi2*Shear(i,1);
                ponder(i1+e1,i2+e2,1)  += phi1*phi2;
            }
        }
    }
    for (int e=0;e<2;e++) for (int i2=0;i2<N+2;i2++) 
    for (int i1=0;i1<N+2;i1++) Ima2D(i1,i2,e)=Ima2D(i1,i2,e)/(ponder(i1,i2,e)+(ponder(i1,i2,e)==0.));
}






void HelmCatalog::div_point_value(fltarray & Ima2D, fltarray & Shear)
{
 Shear.resize(np(), 2);
 for (int i=0;i<np();i++) Shear(i,0) = Shear(i,1) = 0;
 
 for (int i=0;i<np();i++) 
 {
  int i1 = posx(i); int i2 = posy(i);
  float x1 =(TabX(i)-MinX)/BinCat;
  float x2 =(TabY(i)-MinY)/BinCat;
  
  /* dans V1 x V0 */
  Shear(i,0)+=Ima2D(i1+1,i2,0,DEF_DIVROT_BORD)*0.5*(x1-i1)*(x1-i1)*(1-fabs(x2-i2));
  Shear(i,0)+=Ima2D(i1+1,i2+1,0,DEF_DIVROT_BORD)*0.5*(x1-i1)*(x1-i1)*(1-fabs(x2-i2-1));
  Shear(i,0)+=Ima2D(i1,i2,0,DEF_DIVROT_BORD)*(0.75-(x1-i1-0.5)*(x1-i1-0.5))*(1-fabs(x2-i2));
  Shear(i,0)+=Ima2D(i1,i2+1,0,DEF_DIVROT_BORD)*(0.75-(x1-i1-0.5)*(x1-i1-0.5))*(1-fabs(x2-i2-1));
  Shear(i,0)+=Ima2D(i1-1,i2,0,DEF_DIVROT_BORD)*0.5*(x1-i1-1)*(x1-i1-1)*(1-fabs(x2-i2));
  Shear(i,0)+=Ima2D(i1-1,i2+1,0,DEF_DIVROT_BORD)*0.5*(x1-i1-1)*(x1-i1-1)*(1-fabs(x2-i2-1));
  
  /* dans V0 x V1 */
  Shear(i,1)+=Ima2D(i1,i2+1,1,DEF_DIVROT_BORD)*0.5*(x2-i2)*(x2-i2)*(1-fabs(x1-i1));
  Shear(i,1)+=Ima2D(i1+1,i2+1,1,DEF_DIVROT_BORD)*0.5*(x2-i2)*(x2-i2)*(1-fabs(x1-i1-1));
  Shear(i,1)+=Ima2D(i1,i2,1,DEF_DIVROT_BORD)*(0.75-(x2-i2-0.5)*(x2-i2-0.5))*(1-fabs(x1-i1));
  Shear(i,1)+=Ima2D(i1+1,i2,1,DEF_DIVROT_BORD)*(0.75-(x2-i2-0.5)*(x2-i2-0.5))*(1-fabs(x1-i1-1));
  Shear(i,1)+=Ima2D(i1,i2-1,1,DEF_DIVROT_BORD)*0.5*(x2-i2-1)*(x2-i2-1)*(1-fabs(x1-i1));
  Shear(i,1)+=Ima2D(i1+1,i2-1,1,DEF_DIVROT_BORD)*0.5*(x2-i2-1)*(x2-i2-1)*(1-fabs(x1-i1-1));
 }
}


void HelmCatalog::div_bord_point_value(fltarray & Ima2D, fltarray & Shear,int N)
{
    float phi1,phi2;
    Shear.resize(np(), 2);
    for (int i=0;i<np();i++) Shear(i,0) = Shear(i,1) = 0;
    
    // float BinCatX=(MaxX-MinX)/N;  float BinCatY=(MaxY-MinY)/N;
    float BinCatX=BinCat;
    float BinCatY=BinCat;    
    for (int i=0;i<np();i++) 
    {
        float x1=(TabX(i)-MinX)/BinCatX; 
        float x2=(TabY(i)-MinY)/BinCatY;
        int i1=(int) x1; 
        int i2=(int) x2; 
        i1-=(i1==N); 
        i2-=(i2==N);
        
        /* dans V1 x V0 */
        for (int e1=0;e1<3;e1++) 
        {
            if (i1==0) phi1=(e1==0)*(1.-x1)*(1.-x1)+(e1==1)*x1*(2.-3./2.*x1)+(e1==2)*x1*x1/2.;
            else if (i1==(N-1)) 
            { 
               float X=N-x1; phi1=(e1==0)*X*X/2.+(e1==1)*X*(2.-3./2.*X)+(e1==2)*(1.-X)*(1.-X);
            }
            else phi1=3./4.*(e1==1)+(1./2.-3./2.*(e1==1))*(x1-i1-(2.-e1)/2.)*(x1-i1-(2.-e1)/2.);
            for (int e2=0;e2<2;e2++) 
            {
                phi2=1-fabs(x2-i2-e2);
                Shear(i,0)+=phi1*phi2*Ima2D(i1+e1,i2+e2,0);
            }
        }
        
        /* dans V0 x V1 */
        for (int e1=0;e1<2;e1++) 
        {
            phi1=1-fabs(x1-i1-e1);
            for (int e2=0;e2<3;e2++) 
            {
                if (i2==0) phi2=(e2==0)*(1.-x2)*(1.-x2)+(e2==1)*x2*(2.-3./2.*x2)+(e2==2)*x2*x2/2.;
                else if (i2>=(N-1)) 
                { 
                   float X=N-x2; 
                   phi2=(e2==0)*X*X/2.+(e2==1)*X*(2.-3./2.*X)+(e2==2)*(1.-X)*(1.-X);
                }
                else phi2=3./4.*(e2==1)+(1./2.-3./2.*(e2==1))*(x2-i2-(2.-e2)/2.)*(x2-i2-(2.-e2)/2.);
                
                Shear(i,1)+=phi1*phi2*Ima2D(i1+e1,i2+e2,1);
            }
        }
    }
}
//



void HelmCatalog::curl_point_value(fltarray & Ima2D, fltarray & Shear)
{
 Shear.resize(np(), 2);
 for (int i=0;i<np();i++) Shear(i,0) = Shear(i,1) = 0;

 for (int i=0;i<np();i++) {
  int i1 = posx(i); int i2 = posy(i);
  float x1 =(TabX(i)-MinX)/BinCat;
  float x2 =(TabY(i)-MinY)/BinCat;
  
  /* dans V0 x V1 */
  Shear(i,0)+=Ima2D(i1,i2+1,0,DEF_DIVROT_BORD)*0.5*(x2-i2)*(x2-i2)*(1-fabs(x1-i1));
  Shear(i,0)+=Ima2D(i1+1,i2+1,0,DEF_DIVROT_BORD)*0.5*(x2-i2)*(x2-i2)*(1-fabs(x1-i1-1));
  Shear(i,0)+=Ima2D(i1,i2,0,DEF_DIVROT_BORD)*(0.75-(x2-i2-0.5)*(x2-i2-0.5))*(1-fabs(x1-i1));
  Shear(i,0)+=Ima2D(i1+1,i2,0,DEF_DIVROT_BORD)*(0.75-(x2-i2-0.5)*(x2-i2-0.5))*(1-fabs(x1-i1-1));
  Shear(i,0)+=Ima2D(i1,i2-1,0,DEF_DIVROT_BORD)*0.5*(x2-i2-1)*(x2-i2-1)*(1-fabs(x1-i1));
  Shear(i,0)+=Ima2D(i1+1,i2-1,0,DEF_DIVROT_BORD)*0.5*(x2-i2-1)*(x2-i2-1)*(1-fabs(x1-i1-1));
  
  /* dans V1 x V0 */
  Shear(i,1)+=Ima2D(i1+1,i2,1,DEF_DIVROT_BORD)*0.5*(x1-i1)*(x1-i1)*(1-fabs(x2-i2));
  Shear(i,1)+=Ima2D(i1+1,i2+1,1,DEF_DIVROT_BORD)*0.5*(x1-i1)*(x1-i1)*(1-fabs(x2-i2-1));
  Shear(i,1)+=Ima2D(i1,i2,1,DEF_DIVROT_BORD)*(0.75-(x1-i1-0.5)*(x1-i1-0.5))*(1-fabs(x2-i2));
  Shear(i,1)+=Ima2D(i1,i2+1,1,DEF_DIVROT_BORD)*(0.75-(x1-i1-0.5)*(x1-i1-0.5))*(1-fabs(x2-i2-1));
  Shear(i,1)+=Ima2D(i1-1,i2,1,DEF_DIVROT_BORD)*0.5*(x1-i1-1)*(x1-i1-1)*(1-fabs(x2-i2));
  Shear(i,1)+=Ima2D(i1-1,i2+1,1,DEF_DIVROT_BORD)*0.5*(x1-i1-1)*(x1-i1-1)*(1-fabs(x2-i2-1));
 }
}

void HelmCatalog::curl_bord_point_value(fltarray & Ima2D, fltarray & Shear,int N)
{
    // int N = CatNx;
    float phi1,phi2;
    Shear.resize(np(), 2);
    for (int i=0;i<np();i++) Shear(i,0) = Shear(i,1) = 0;
    
    // float BinCatX=(MaxX-MinX)/N;  float BinCatY=(MaxY-MinY)/N;
    float BinCatX=BinCat;
    float BinCatY=BinCat;    
    
    for (int i=0;i<np();i++) 
    {
        float x1=(TabX(i)-MinX)/BinCatX; 
        float x2=(TabY(i)-MinY)/BinCatY;
        int i1=(int) x1; 
        int i2=(int) x2; i1-=(i1==N); i2-=(i2==N);
        
        /* dans V0 x V1 */
        for (int e1=0;e1<2;e1++) 
        {
            phi1=1-fabs(x1-i1-e1);
            for (int e2=0;e2<3;e2++) 
            {
                if (i2==0) phi2=(e2==0)*(1.-x2)*(1.-x2)+(e2==1)*x2*(2.-3./2.*x2)+(e2==2)*x2*x2/2.;
                else if (i2>=(N-1)) 
                { 
                   float X=N-x2; 
                   phi2=(e2==0)*X*X/2.+(e2==1)*X*(2.-3./2.*X)+(e2==2)*(1.-X)*(1.-X);
                }
                else phi2=3./4.*(e2==1)+(1./2.-3./2.*(e2==1))*(x2-i2-(2.-e2)/2.)*(x2-i2-(2.-e2)/2.);
                
                Shear(i,0)+=phi1*phi2*Ima2D(i1+e1,i2+e2,0);
            }
        }
        
        /* dans V1 x V0 */
        for (int e1=0;e1<3;e1++) 
        {
            if (i1==0) phi1=(e1==0)*(1.-x1)*(1.-x1)+(e1==1)*x1*(2.-3./2.*x1)+(e1==2)*x1*x1/2.;
            else if (i1==(N-1)) 
            { 
               float X=N-x1; 
               phi1=(e1==0)*X*X/2.+(e1==1)*X*(2.-3./2.*X)+(e1==2)*(1.-X)*(1.-X);
            }
            else phi1=3./4.*(e1==1)+(1./2.-3./2.*(e1==1))*(x1-i1-(2.-e1)/2.)*(x1-i1-(2.-e1)/2.);
            for (int e2=0;e2<2;e2++) 
            {
                phi2=1-fabs(x2-i2-e2);
                
                Shear(i,1)+=phi1*phi2*Ima2D(i1+e1,i2+e2,1);
            }
        }
    }
}

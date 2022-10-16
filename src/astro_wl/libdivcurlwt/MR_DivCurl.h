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
**    File:  MR_DivCurl.ch
**
*******************************************************************************
**
**    DESCRIPTION  div and curl anisotropic decomposition class definition
**    ----------- 
**                 
******************************************************************************/

#ifndef _DIVCURL_H_
#define _DIVCURL_H_

#include "GlobalInc.h"
#include "IM_Obj.h"
#include "IM_IO.h"
#include "Filter.h"
#include "SB_Filter1D.h"
#include "SB_Filter.h"
#include "WT_Bord.h"

#define DEF_DIVROT_BORD I_PERIOD   // default border manadgement. Other are I_PERIOD, I_MIRROR, I_ZERO, I_CONT


// Catalog for a vectorial data set 
// TabXYG(i,0) = x
// TabXYG(i,1) = y
// TabXYG(i,2) = X-component
// TabXYG(i,3) = Y-component
// TabXYG(i,4) = Sigma-X
// TabXYG(i,5) = Sigma-Y

class HelmCatalog { 
    fltarray TabX;
    fltarray TabY;
    fltarray TabSigmaX;
    fltarray TabSigmaY;
    
    int CatNp;
    int CatNx;
    int CatNy;
    float BinCat,BinCatX,BinCatY;
    float MinX,MaxX,MinY,MaxY;
    
   public: 
    intarray TabNp;  // TabNp(0:CatNx-1, 0:CatNx-1) = image of nbr of points per pixel position
    fltarray DataG;  // DataG(i, 0) = shear-gamma1 or q-pol for the ith point
                     // DataG(i, 1) = shear-gamma2 or u-pol for the ith point
    
   // return the number of points in the catalog
    int np() const {return CatNp;}  
    // return the x image size
    int nx() const {return CatNx;}   
    // return the y image size
    int ny() const {return CatNy;}   
    // return the number of points in pixel (j,i)
    int np_in_xy(int j, int i)  {return TabNp(j,i);}  
   
    // return the x position of the ith point of the catalog
    int posx(int i) {
        int pos_x=(int)((TabX(i) -MinX) / BinCat); 
        if (pos_x==nx()) pos_x=(int)(nx()-1.e-06);
        return pos_x; }
    // return the y position of the ith point of the catalog
    int posy(int i) {
        int pos_y=(int)((TabY(i) -MinY) / BinCat); 
        if (pos_y==ny()) pos_y=(int)(ny()-1.e-06);
        return pos_y; } 
    // return the noise standard deviation of the ith point
    float sigmax(int i) {return TabSigmaX(i);} 
    float sigmay(int i) {return TabSigmaY(i);} 
     
    Bool Verbose;

    void alloc(fltarray & TabXYG, float BCat=1.);
    void allocn(fltarray & TabXYG, int N);
    
    // Create an image from the catalog by averaging all points falling in a given pixel
    // Shear_or_Pol1D = input fltarr[0:np()-1,0:1] (same size as DataG)
    // Ima = outout 3D fltarr [0:nx()-1, 0:ny()-1, 0:1] = image of average shear of pola.
    void cat2ima(fltarray & Shear_or_Pol1D, fltarray & Ima);
    
    // same as before but input is DataG
    void cat2ima(fltarray & Ima);
    // Create the standard deviation map related to the average image = outout 3D fltarr [0:nx()-1, 0:ny()-1, 0:1]
    void cat2imasigma(fltarray & ImaSigmaXY);
    
    // V0 projection and point value interpolation for the divergence WT
    void div_quasi_interpol(fltarray & Shear_or_Pol1D, fltarray & Ima2D);
    void div_bord_quasi_interpol(fltarray & Shear, fltarray & Ima2D,int N);
    
    void div_point_value(fltarray & Ima2D, fltarray & Shear_or_Pol1D);
    void div_bord_point_value(fltarray & Ima2D, fltarray & Shear,int N);

    // V0 projection and point value interpolation for the curl WT
    void curl_point_value(fltarray & Ima2D, fltarray & Shear_or_Pol1D);
    void curl_bord_point_value(fltarray & Ima2D, fltarray & Shear,int N);

    void curl_quasi_interpol(fltarray & Shear_or_Pol1D, fltarray & Ima2D);
    void curl_bord_quasi_interpol(fltarray & Shear, fltarray & Ima2D,int N);

    // We may want a final solution of the grid
    void div_grid_point_value(fltarray & Ima2D, fltarray & ImaGrid2D);
    void curl_grid_point_value(fltarray & Ima2D, fltarray & ImaGrid2D);

     HelmCatalog() {BinCat=1; Verbose=False;CatNp=0;CatNx=0;CatNy=0;}
     ~HelmCatalog() {};
};



class DivCurl {
         public:
           Bool BorderWavelet;
 	       Bool PointValueRec;
	       Bool V0_Proj;
		   fltarray RotDiv_Trans;
           HelmCatalog *Cat;
           
 	       virtual void transform(fltarray &Data, Bool Div=True){};
           virtual void recons(fltarray &Data, Bool Div=True){};

		   // virtual void div_transform (fltarray &Data, fltarray &Trans){};
           // virtual void curl_transform (fltarray &Data, fltarray &Trans){};
		   
           // virtual void div_recons (fltarray &Trans, fltarray &Data){};
		   // virtual void curl_recons (fltarray &Trans, fltarray &Data){};
           virtual ~DivCurl(){};
         
		
};


void div_point_value(fltarray & Tabc, fltarray & Tabp, Bool BorderWavelet=False);
// compute the point value on a regular grid
// Tabc = fltarray(0..N-1,0..N-1,2) = input scaling function coefficients at the finest scale
//                                    V1xV0, V0xV1
// Tabp= fltarray(0..N-1,0..N-1,2) =  output point values on a regular grid

void curl_point_value(fltarray & Tabc, fltarray & Tabp, Bool BorderWavelet=False);
// compute the point value on a regular grid
// Tabc = fltarray(0..N-1,0..N-1,2) = input scaling function coefficients at the finest scale
//                                    V0xV1, V1xV0
// Tabp= fltarray(0..N-1,0..N-1,2) =  output point values on a regular grid


void div_quasi_interpol(fltarray & Tabp, fltarray & Tabc, Bool BorderWavelet=False);
void curl_quasi_interpol(fltarray & Tabp, fltarray & Tabc, Bool BorderWavelet=False);


/***************************************/
//
//  CLASS definition for DIV-free and ROT-free transform and reconstruction

class MR_ANI_DIVCURL: public DivCurl {
    FilterAnaSynt FAS_1;    // First filter  
    FilterAnaSynt FAS_2;    // Second filter  
    SubBandFilter *SBF_1;   // First filter bank class
    SubBandFilter *SBF_2;   // Second filter bank class
    BordLineCol *WT_1; // First 2D Wavelet Transform Class
    BordLineCol *WT_2; // Second 2D Wavelet Transform Class
	
	void init_filter_bank(); // initalize the two filter banks
    int NbrScale;            // Number of scales used in the decomposition 
    int NbrBand;             // Number of bands used in the decomposition = 3 * (NbrScale-1)  + 1
	
 	int NbrScaleX;           
	int NbrScaleY;
	intarray TabPosX;
	intarray TabPosY;
	intarray TabNy;
	intarray TabNx;
 	
	int Nl;       // input image number of lines
    int Nc;       // input image number of columns
    Bool DivTrans;     // True if we use a dif-free transform, and false for a rot-free transform

	int indexi(int by, int i) const {return TabPosY(by)+i;}   
	int indexj(int bx, int j) const {return TabPosX(bx)+j;}
    int indexij(int by, int bx, int i, int j) const {return indexi(by, i)*Nc + indexj(bx,j);}

    void div_transform (fltarray &Data, fltarray &Trans);
    void div_bord_transform (fltarray &Data, fltarray &Trans);
    void curl_transform (fltarray &Data, fltarray &Trans);
    void curl_bord_transform (fltarray &Data, fltarray &Trans);
	void div_recons (fltarray &Trans, fltarray &Data);
    void div_bord_recons (fltarray &Trans, fltarray &Data);
    void curl_recons (fltarray &Trans, fltarray &Data);
    void curl_bord_recons (fltarray &Trans, fltarray &Data);
	
   public:
    sb_type_norm Norm;
    
    int Debug;
 
    int nbr_band_nx() const { return NbrScaleX;}   // return the number of bands on the x-axis
    int nbr_band_ny() const { return NbrScaleY;}   // return the number of bands on the x-axis
    int size_band_ny(int by) const { return TabNy(by);} // return the number of lines of a given band bt in [0 ... NbrScaleY-1]
	int size_band_nx(int bx) const { return TabNx(bx);} // return the number of column of a given band
    float & operator() (int Num, int bx, int by, int j, int i) const { return RotDiv_Trans(TabPosX(bx)+j, TabPosY(by)+i, Num);}  
       // return the value for component Num (Num = 0 or 1), the band (bx,by), and pixel position (j,i)
	
	
	int nbr_band() const { return nbr_band_nx()*nbr_band_ny();}   // return the number of bands
    int size_band_nl(int b) const { return size_band_ny(b / nbr_band_ny());} // return the number of lines of a given band
    int size_band_nc(int b) const { return  size_band_nx(b % nbr_band_ny()); } // return the number of column of a given band
 
    Bool Verbose;
 	
    MR_ANI_DIVCURL() {Norm = NORM_L1;PointValueRec=True;V0_Proj=True;Cat=NULL;Debug=0;BorderWavelet=False;}
    void alloc(int Nl, int Nc, int NbrPlan=0); // Allocation of the class

    void div_bord_transform (fltarray &Data1, fltarray &Data2, fltarray &Trans1, fltarray &Trans2);
    void div_bord_recons ( fltarray &Trans1, fltarray &Trans2, fltarray &Data1, fltarray &Data2);
    void curl_bord_transform (fltarray &Data1, fltarray &Data2, fltarray &Trans1, fltarray &Trans2);
    void curl_bord_recons ( fltarray &Trans1, fltarray &Trans2, fltarray &Data1, fltarray &Data2);
    
    void transform(fltarray &Data, Bool Div=True);
    // Main routine to perform a decomposition
    // Data = input data (i.e. fltarray(*,*,2))
    // By default, an isotropic  div-free transform is done
    //    if Div == False then a rot-free transform is done
    //    if Anisotrop == True then an anisotropic transform is done

    void info();  // print statistical information relative to each band
    
    void recons(fltarray &Data, Bool Div=True);
    // Reconstruction of the vector field from its decomposition

    // void component_recons(fltarray &Data);
    // Apply an inverse wavelet coefficients on each of the two wavelet components
    // and store the result in Data

    void set2zero();
    // Set to zero the decomposition
    
   ~MR_ANI_DIVCURL() {}; 
};


class MR_ISO_DIVCURL : public DivCurl {
    FilterAnaSynt FAS_1;    // First filter  
    FilterAnaSynt FAS_2;    // Second filter  
    SubBandFilter *SBF_1;   // First filter bank class
    SubBandFilter *SBF_2;   // Second filter bank class
    HALF_DECIMATED_2D_WT *WT_1; // First 2D Wavelet Transform Class
    HALF_DECIMATED_2D_WT *WT_2; // Second 2D Wavelet Transform Class
    Ifloat * TabWT_Trans_0;  // Array of bands containing the first wavelet transform
    Ifloat * TabWT_Trans_1;  // Array of bands containing the second wavelet transform
    Ifloat * TabTransDivRot_0;  // Array of bands containing the first component of the div-free transform or the rot-free transform
    Ifloat * TabTransDivRot_1;  // Array of bands containing the second component of the div-free transform or the rot-free transform
    void init_filter_bank(); // initalize the two filter banks
    int NbrScale;            // Number of scales used in the decomposition 
    int NbrBand;             // Number of bands used in the decomposition = 3 * (NbrScale-1)  + 1
	
 
	
    int NbrUndecimatedScale; // Number of undecimated scales (default is 0)
    int Nl;       // input image number of lines
    int Nc;       // input image number of columns
	Bool DivTrans;     // True if we use a dif-free transform, and false for a rot-free transform

	void div_isotrop_transform(fltarray &Data, fltarray &Trans);  // isotropic div-free transform
    void div_isotrop_recons(fltarray &Trans, fltarray &Data);     // isotropic div-free  reconstruction
    void rot_isotrop_transform(fltarray &Data, fltarray &Trans);  // isotropic rot-free transform
    void rot_isotrop_recons(fltarray &Trans, fltarray &Data);     // isotropic rot-free  reconstruction
	
public:
	sb_type_norm Norm;
    
    int nbr_band() const { return NbrBand;}   // return the number of bands
    int nbr_scale() const { return NbrScale;} // return the number of scales
    int size_band_nl(int b) const { return TabWT_Trans_0[b].nl();} // return the number of lines of a given band
    int size_band_nc(int b) const { return TabWT_Trans_0[b].nc();} // return the number of column of a given band
    float & operator() (int Num, int b, int i, int j) const;  // return the value for component Num (Num = 0 or 1), the band b, and pixel position (i,j)
	
	Bool Verbose;
	MR_ISO_DIVCURL() {Norm = NORM_L1;PointValueRec=True;V0_Proj=True;NbrUndecimatedScale=0;Cat=NULL;}
	
    void alloc(int Nl, int Nc, int NbrUndec=0, int NbrPlan=0); // Allocation of the class
    
     
	//    void quasi_interpol_div(fltarray &Data);  // quasi-interpolation of data in spline space (quadratic splines in the x-direction)
	//    void quasi_interpol_curl(fltarray &Data);  // quasi-interpolation of data in spline space (quadratic splines in the y-direction)
	//    void point_value_div(fltarray &Data);  // compute the point values from spline space (quadratic splines in the x-direction)
	//    void point_value_curl(fltarray &Data);  // compute the point values from spline space (quadratic splines in the y-direction)
    
    void band2fltarr(fltarray &Data, Bool DivRotTrans=True); 
    // store the decomposition in a fltarray class
    // For a DECIMATED transform, Data is a fltarray(*,*,2) object,
    //                            where Data(*,*,0) is the first wavelet component
    //                              and Data(*,*,1) is the first wavelet component
    // For an unDECIMATED transform, Data is a fltarray(*,*, nbr_of_band * 2) object
    //                           where Data(*,*,0:nbr_band-1) is the first wavelet component
    //                             and Data(*,*, nbr_band:2*nbr_band-1) is the second wavelet component
    // if DivRotTrans == True then the div-free (resp rot-free) is stored
    // if DivRotTrans == False then the wavelet transforms of the x and y field are stored
    
    void fltarr2band(fltarray &Data, Bool DivRotTrans=True);
    //  Set the decomposition from data given in Data
    //  if DivRotTrans == True then the div-free (resp rot-free) is set
    //  if DivRotTrans == False then the wavelet transforms of the x and y field are set
	
    void transform(fltarray &Data, Bool Div=True);
    // Main routine to perform a decomposition
    // Data = input data (i.e. fltarray(*,*,2))
    // By default, an isotropic  div-free transform is done
    //    if Div == False then a rot-free transform is done
 	
    void info();  // print statistical information relative to each banf
    
    void recons(fltarray &Data, Bool Div=True);
    // Reconstruction of the vector field from its decomposition
	
    void component_recons(fltarray &Data);
    // Apply an inverse wavelet coefficients on each of the two wavelet components
    // and store the result in Data
    
    void read_from_tab(fltarray &Trans, Bool IsDivTrans=True);
    // Set the decomposition to data contained in Trans, and initialize the paramter
    // DivTrans  of the class.
    
    void set2zero();
    // Set to zero the decomposition
    
    ~MR_ISO_DIVCURL() {}; 
};

enum type_helmholtz {H_ANI_WT, H_ISO_WT, H_FFT};
#define DEF_HELMHOLTZ H_ANI_WT
#define DEF_ITER_HELMHOLTZ 10

const char * StringHelmholtzTransform (type_helmholtz type);

enum type_constraints_helmholtz {C_NO, C_HARD, C_SOFT, C_ZERO_BORDER, C_ZERO};
const char * StringHelmholtzConstraint (type_constraints_helmholtz type);

class Helmholtz {
	MR_ANI_DIVCURL *MAD;
	MR_ISO_DIVCURL *MID;
  	type_helmholtz THelm;
	DivCurl *DC;
	HelmCatalog *HCat;
	Bool UseCat;
	void fft_trans(fltarray &Data);
	void iter_trans(fltarray &Data);
    void cat_iter_trans(fltarray &Data);
public:
    Bool BorderWavelet;
	Bool Verbose;
	fltarray U;
	fltarray V;
	fltarray EB;
	Bool PointValueRec;
	Bool V0_Proj;
	int NbrUndec;
	int Niter;
    type_constraints_helmholtz THCons;
    Bool UseConstraint_on_B;  // if true,  All common modes in E mode 
    Bool UseConstraint_on_E;  // if true,  All common modes in E mode 
	Helmholtz() {Niter=DEF_ITER_HELMHOLTZ;PointValueRec=True;V0_Proj=True;DC=NULL;NbrUndec=0; HCat=NULL;UseCat=False;BorderWavelet=False;UseConstraint_on_B=True;
                 UseConstraint_on_E=False;THCons=C_NO;}
	void decomposition(fltarray &Data);
	void alloc(type_helmholtz TypeHel, int Nx, int Ny);
	void alloc(type_helmholtz TypeHel, HelmCatalog & Cat);
	~Helmholtz(){};
};


#endif


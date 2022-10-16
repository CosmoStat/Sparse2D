/******************************************************************************
**                   Copyright (C) 2000 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Philippe
**
**    Date:  5/03/2000 
**    
**    File:  linelet3d.cc
**
**    Modification history:
**
*******************************************************************************
**
**    DESCRIPTION linelet3d transform and reconstruction
**    -----------  
**                 
******************************************************************************/
#include "Linelet3D.h"
#include "SB_Filter.h"
#include "macro.h"
#include "IM3D_IO.h"
#include "MR_Obj.h"

extern Bool Verbose;

Linelet3D::Linelet3D()
{
	reset();
}  

Linelet3D::~Linelet3D()
{
	if(Verbose) cerr<<"Linelet3D::~Linelet3D()"<<endl;
	dealloc();
}

void Linelet3D::dealloc()
{
	if(TabBandNy!=NULL) delete [] TabBandNy;
	if(TabBandNx!=NULL) delete [] TabBandNx;
	if(TabBandNz!=NULL) delete [] TabBandNz;
	if(TabDepx!=NULL) delete [] TabDepx;
	if(TabDepy!=NULL) delete [] TabDepy;
	if(TabDepz!=NULL) delete [] TabDepz;
	AllocClass=False;
}



/******************************************************************************/
// reset function
/******************************************************************************/
 
void Linelet3D::reset()
{
    _Control=False;
    _WriteFile=False;
    _StatInfo=False;
    _TakeOnlyCompletePlane=False;
    _LinTransf=LIN3D_UNKNOWN;
    _Verbose=False;
    _Nx=_Ny=_Nz=0;
    _BlockOverlap=False;
	_OverlapFactor=0;
    _Nxb=_Nyb=_Nzb=0;
    _BlockSize=0;
    _NbrBlock=0;
    _NbrScale=0;
    _RadNbPlane=0;
    _LinNbPlane=0;
    _LinPlaneSize=0;
    _NbPlaneUsedInStat=0;
    _NbPlaneUsedInStatByBlock=0;
    NbrBand = 0;
    TabBandNx = NULL;
    TabBandNy = NULL;
    TabBandNz = NULL;
    TabDepx = NULL;
    TabDepy = NULL;
    TabDepz = NULL;
	AllocClass=False;
	no_recons=false;
	keep_energy=false;
}

/******************************************************************************/
// linelet3d transform info
/******************************************************************************/

const char * StringLin3DTransform (type_linelet3d_WTtrans  type) {
   switch (type) {
   case  LIN3D_OWT:     return ("Standard bi-orthogonal WT"); 
   case  LIN3D_PYR_FFT: return ("FFT based Pyramidal WT"); 
   default:             return ("Undefined transform");
   }
}

type_linelet3d_WTclass lin3d_class(type_linelet3d_WTtrans Type) {
   switch (Type) {
   case  LIN3D_OWT:     return  LIN3D_CL_ORTHO; 
   case  LIN3D_PYR_FFT: return  LIN3D_CL_PYR; 
   default:             return  LIN3D_CL_UNKNOW;
  }
}
  

/******************************************************************************/
// max_scale_number
/******************************************************************************/

static int max_scale_number (int N)
{
    int ScaleMax;
    ScaleMax  = iround((float)log((float) (N / 4. * 3.) / log(2.)));
    return ScaleMax;
} 

/******************************************************************************/
// info, verbose, stat on linelet3d classes
/******************************************************************************/
   
void Linelet3D::info_stat(fltarray &Transf, char* FileSimu,              
                          Bool WithFileSimu) 
{
   fltarray StatLinelet;
   fltarray Simu;
   if (WithFileSimu) {
      StatLinelet.alloc(_NbrScale,7);
      fits_read_fltarr(FileSimu, Simu);
   } else StatLinelet.alloc(_NbrScale,6);
   fltarray MaxVal(_NbrScale, 100);
   
   char Name[50];
   for (int s=0;s<_NbrScale;s++) {
  
      float SigmaSimu=0;
      if (WithFileSimu) {
         SigmaSimu = Simu(s);
         cout << "SigmaSimu:" << SigmaSimu << endl;
      }
      sprintf (Name, "Scale%d stat:", s);
      //(_pScale[s]).info(Name);
      
      StatLinelet(s,0)=(_pScale[s]).mean();
      StatLinelet(s,1)=(_pScale[s]).sigma();
      StatLinelet(s,2)= skewness((_pScale[s]).buffer(), 
                                 (_pScale[s]).n_elem());
      StatLinelet(s,3)= curtosis((_pScale[s]).buffer(), 
                                 (_pScale[s]).n_elem());
      StatLinelet(s,4)=0;
	for (int i=0;i<(_pScale[s]).n_elem();i++)
	   if (fabs((_pScale[s])(i)-StatLinelet(s,0)) > 3*StatLinelet(s,1))
	      StatLinelet(s,4) += ((_pScale[s])(i)*(_pScale[s])(i));				 
      StatLinelet(s,5)=(_pScale[s]).energy();
	          
      if (WithFileSimu) {
         StatLinelet(s,6)=0;
	    for (int i=0;i<(_pScale[s]).n_elem();i++) 
	       if (fabs((_pScale[s])(i)- StatLinelet(s,0))> 3*SigmaSimu)
	          StatLinelet(s,6) += ((_pScale[s])(i)*(_pScale[s])(i));				 
      }	 
      
      //!!!!!!!!!!!!!!! modify the pScale[s] value !!!!!!!!!!!!!!!!!!!!!!!
     fltarray Prov=_pScale[s];
     for (int i=0;i<100;i++) {
        int IndMax=0;
	MaxVal(s,i) = Prov.maxfabs(IndMax);
	Prov(IndMax)=0.0;
     }  			 
       			    		    
   }   
   
   fits_write_fltarr ((char*)"StatLinTransf.fits", StatLinelet);   
   fits_write_fltarr ((char*)"MaxLinelet.fits", MaxVal);
}


void Linelet3D::trace_param() {

}                  
 
 
void Linelet3D::stat_alloc () {
   
   intarray* MemNotUsed = PRad3D.getMemNotUsed();
   // if (_Verbose) MemNotUsed->info("Plane information"); 
   
   if (_TakeOnlyCompletePlane) {
      for (int k=0;k<_RadNbPlane;k++) 
         if ((*MemNotUsed)(k) == 0) _NbPlaneUsedInStatByBlock++;
   } else {
      _NbPlaneUsedInStatByBlock = _RadNbPlane;
   }
   _NbPlaneUsedInStat = _NbrBlock*_NbPlaneUsedInStatByBlock;
      
   _pScale =  new fltarray [_NbrScale];
   
   int NbPoint = _BlockSize;
   for (int s=0;s<_NbrScale-1;s++) {
   
      int ScaleNbPt = NbPoint/2;
      int BegInd    = NbPoint-ScaleNbPt;
      
      _pScale[s].alloc (ScaleNbPt*ScaleNbPt + 2*(NbPoint-ScaleNbPt)*ScaleNbPt, 
                       _NbPlaneUsedInStat);
      
      NbPoint = BegInd;
   }
   
   _pScale[_NbrScale-1].alloc (NbPoint*NbPoint, _NbPlaneUsedInStat);
      
}

void Linelet3D::stat_add (int Bi, int Bj, int Bk, fltarray& CubeBlock) {

   int NumBlock = num_block(Bi,Bj,Bk);

   fltarray* CoefSigma = PRad3D.getCoefSigma();
   if (_DoNotCorrectSigma) CoefSigma->init(1.);
   // if (_Verbose) CoefSigma->info("Corrected SigmaCoef");
   
   intarray* MemNotUsed = PRad3D.getMemNotUsed();
   // if (_Verbose) MemNotUsed->info("Plane information"); 
/*      
fltarray prov(_BlockSize,_BlockSize,_RadNbPlane);
for (int k=0;k<_RadNbPlane;k++) 
for (int i=0;i<_BlockSize;i++)
for (int j=0;j<_BlockSize;j++)  
   prov(i,j,k) = (*CoefSigma)(k)*CubeBlock (i,j,k);
char Name[50];
sprintf(Name,"CorBlckLin%d.fits",NumBlock);
fits_write_fltarr(Name, prov);
CubeBlock.info("block transf:");    
prov.info("Corrected block transf:");   
*/   
   
   fltarray SigmaInf(_NbrScale,2);
   
   int NbPoint = _BlockSize;
   for (int s=0;s<_NbrScale-1;s++) {
      
      int ScaleNbPt = NbPoint/2;           // nb point in current scale
      int BegInd    = NbPoint-ScaleNbPt;   // begin indice in CubeBlock
                                           // for current scale (for 
					   // diagonal block)   
      int currentPoint=0;                
      int currentPlane=0;
      
      fltarray Val ( ScaleNbPt*ScaleNbPt +              // diag block
                     2*(NbPoint-ScaleNbPt)*ScaleNbPt);  // hor and vert
			                                // block
      fltarray SigmaofBlock (_NbPlaneUsedInStat/_NbrBlock); // sigma on each plane   
      
      for (int k=0;k<_RadNbPlane;k++) {
      
         if (    !_TakeOnlyCompletePlane
	     || (_TakeOnlyCompletePlane && (*MemNotUsed)(k) == 0)) {
        
            currentPoint=0;
    
            for (int i=0;i<ScaleNbPt;i++)
            for (int j=0;j<ScaleNbPt;j++) 
	       (_pScale[s]) (currentPoint, 
	                     currentPlane + NumBlock*_NbPlaneUsedInStatByBlock) 
	                   = Val(currentPoint++)
		           = (*CoefSigma)(k)*CubeBlock(BegInd+i,BegInd+j,k);
	 
	    for (int i=0;i<ScaleNbPt;i++)
            for (int j=0;j<NbPoint-ScaleNbPt;j++) 
	       (_pScale[s]) (currentPoint, 
	                     currentPlane + NumBlock*_NbPlaneUsedInStatByBlock) 
	                   = Val(currentPoint++)
		           = (*CoefSigma)(k)*CubeBlock(BegInd+i,j,k);
	 
	    for (int i=0;i<NbPoint-ScaleNbPt;i++)
            for (int j=0;j<ScaleNbPt;j++) 	 
	       (_pScale[s]) (currentPoint, 
	                     currentPlane + NumBlock*_NbPlaneUsedInStatByBlock) 
	                   = Val(currentPoint++)
		           = (*CoefSigma)(k)*CubeBlock(i,BegInd+j,k);
	    
	    SigmaofBlock(currentPlane) = Val.sigma(); 
	    currentPlane++;
	 }
      }
      	
      NbPoint = BegInd;                   //  point number for next scales
      
      SigmaInf(s,0) = SigmaofBlock.mean();
      SigmaInf(s,1) = SigmaofBlock.sigma();
   }
   
   fltarray LastScale (NbPoint, NbPoint, _NbPlaneUsedInStatByBlock);
   int currentPoint = 0;
   int currentPlane = 0;
   for (int k=0;k<_RadNbPlane;k++) { 
      
      if (    !_TakeOnlyCompletePlane
	   || (_TakeOnlyCompletePlane && (*MemNotUsed)(k) == 0)) {
	      
         currentPoint=0;
         for (int i=0;i<NbPoint;i++)
         for (int j=0;j<NbPoint;j++) {
      
            LastScale(i,j,currentPlane) 
	          = (_pScale[_NbrScale-1]) (currentPoint, 
	                       currentPlane + NumBlock*_NbPlaneUsedInStatByBlock)
                  = (*CoefSigma)(k)*CubeBlock(i,j,k);
	    currentPoint++;
         }
	 currentPlane++;
      }
   }
   
   //LastScale.info("Last Scale"); 
   
   //fits_write_fltarr ("StatLinTransf.fits", SigmaInf);   
}



/******************************************************************************/
// set_block_param - init block param
/******************************************************************************/

void  Linelet3D::set_tab_band()
{
   int Ny_s=_BlockSize;
   int Nx_s=_BlockSize;
   int s,ind_i,ind_j;
   
   NbrBand = 3*(_NbrScale-1)+1;
   TabBandNy = new int [NbrBand];
   TabBandNx = new int [NbrBand];
   TabBandNz = new int [NbrBand];
   TabDepx = new int [NbrBand];
   TabDepy = new int [NbrBand];
   TabDepz = new int [NbrBand];
   for (s = 0; s < NbrBand-1; s+=3)
   {
      TabBandNy[s] = (Ny_s+1)/2;
      TabBandNx[s] = Nx_s/2;
      TabBandNz[s] = _LinNbPlane;
      ind_orthog_transf(s/3,0,0,D_HORIZONTAL, _BlockSize, _BlockSize, ind_i,  ind_j);
      TabDepx[s] = ind_j;
      TabDepy[s] = ind_i;
      TabDepz[s] = 0;
      
      TabBandNy[s+1] = Ny_s/2;
      TabBandNx[s+1] = (Nx_s+1)/2;
      TabBandNz[s+1] = _LinNbPlane;
      ind_orthog_transf(s/3,0,0,D_VERTICAL, _BlockSize, _BlockSize, ind_i,  ind_j);
      TabDepx[s+1] = ind_j;
      TabDepy[s+1] = ind_i;
      TabDepz[s+1] = 0;
      
      TabBandNy[s+2] = Ny_s/2;
      TabBandNx[s+2] = Nx_s/2;
      TabBandNz[s+2] = _LinNbPlane;
      ind_orthog_transf(s/3,0,0,D_DIAGONAL, _BlockSize, _BlockSize, ind_i,  ind_j);
      TabDepx[s+2] = ind_j;
      TabDepy[s+2] = ind_i;
      TabDepz[s+2] = 0;
      Ny_s = (Ny_s+1)/2;
      Nx_s = (Nx_s+1)/2;
   }
   s = NbrBand-1;
   TabBandNy[s] = Ny_s;
   TabBandNx[s] = Nx_s;
   TabBandNz[s] = Nx_s;
   ind_orthog_transf(s,0,0,I_SMOOTH, _BlockSize, _BlockSize, ind_i,  ind_j);
   TabDepx[s] = ind_j;
   TabDepy[s] = ind_i;
   TabDepz[s] = 0;
}

/******************************************************************************/

void Linelet3D::set_block_param() 
{
   if (_BlockSize <= 0) _BlockOverlap = False;
   B3D.BlockOverlap = _BlockOverlap;
   B3D.Overlap = _OverlapFactor;
   B3D.Verbose = _Verbose;
   B3D.alloc(_Nx,_Ny,_Nz,_BlockSize);
   _BlockSize = B3D.block_size();
   _BlockOverlap = B3D.BlockOverlap;
   _Nxb = B3D.nbr_block_nx();
   _Nyb = B3D.nbr_block_ny();
   _Nzb = B3D.nbr_block_nz();
   
   // PRad3D.setVerbose(_Verbose);
   PRad3D.setControl(_Control);
   PRad3D.setWriteFile(_WriteFile);
   PRad3D.setStatInfo(_StatInfo);
   PRad3D.setTakeOnlyCompletePlane(_TakeOnlyCompletePlane);
   PRad3D.setDoNotCorrectSigma(_DoNotCorrectSigma);   
   _NbrScale = max_scale_number(_BlockSize);
	       
   PRad3D.setParam(_BlockSize,_NbrScale);         //!!!!!add type of WT...!!! 
   PRad3D.alloc();
   
   _RadNbPlane = PRad3D.getNbPlane();
   _NbrBlock = _Nxb*_Nyb*_Nzb;
   _LinNbPlane = _RadNbPlane * _NbrBlock;
   _LinPlaneSize = _BlockSize*_BlockSize;
   set_tab_band();
   
   if (_StatInfo) stat_alloc();
   		     
   _PRadTransf = fltarray (_BlockSize, _BlockSize, _RadNbPlane);
    
   switch (lin3d_class(_LinTransf)) {
   case LIN3D_CL_PYR: 
      cout << "Not implemented yet!" << endl; exit(-1);
   case LIN3D_CL_ORTHO:break;
   default:break;
   }
   
   if (_Verbose) 
   {
      printf (" \n\n INITIALIZATION:\n");
      printf (" Nbr block x direction  : %d\n", _Nxb);
      printf (" Nbr block y direction  : %d\n", _Nyb);
      printf (" Nbr block z direction  : %d\n", _Nzb);
      printf (" Block size             : %d\n", _BlockSize);
      printf (" Nbr blocks             : %d\n", _NbrBlock);
      printf (" Nbr of scales          : %d\n", _NbrScale);
      printf (" Nbr of radon planes    : %d\n", _RadNbPlane);
      printf (" Nbr of beamlet planes  : %d\n", _LinNbPlane);
      printf (" Size of each plane     : %d\n", _LinPlaneSize);
   }  
}

void Linelet3D::set_recons_block_param() {

 
   // PRad3D.setVerbose(_Verbose);
   PRad3D.setControl(_Control);
   PRad3D.setWriteFile(_WriteFile);
   PRad3D.setStatInfo(_StatInfo);
   PRad3D.setTakeOnlyCompletePlane(_TakeOnlyCompletePlane);
   PRad3D.setDoNotCorrectSigma(_DoNotCorrectSigma);   
   _NbrScale = max_scale_number(_BlockSize);
   NbrBand = 3*(_NbrScale-1)+1;
   PRad3D.setParam(_BlockSize,_NbrScale);         //!!!!!add type of WT...!!! 
   PRad3D.alloc();
     
   _RadNbPlane = PRad3D.getNbPlane();
   _NbrBlock = _Nxb*_Nyb*_Nzb;
   _LinNbPlane = _RadNbPlane * _NbrBlock;
   _LinPlaneSize = _BlockSize*_BlockSize;
   set_tab_band();
   
   if (_StatInfo) stat_alloc();
   _PRadTransf = fltarray (_BlockSize, _BlockSize, _RadNbPlane);
    
   switch (lin3d_class(_LinTransf)) {
   case LIN3D_CL_PYR: 
      cout << "LIN3D_CL_PYR Not implemented yet!" << endl; exit(-1);
   case LIN3D_CL_ORTHO:
	   break;
   default:
      cout << "Wrong Transform type!" << endl; exit(-1);
   }
   
   
   

}
/******************************************************************************/
// set_attribut
/******************************************************************************/
void Linelet3D::set_attribut(fltarray &Cube, int & BlockSize) {

   _Nx = Cube.nx();
   _Ny = Cube.ny();  
   _Nz = Cube.nz();  
   _BlockSize = BlockSize;
      
}
void Linelet3D::set_attribut(int nx,int ny,int nz, int & BlockSize) {

   _Nx = nx;
   _Ny = ny;
   _Nz = nz;
   _BlockSize = BlockSize;
      
}

/******************************************************************************/
// get band
/******************************************************************************/

void Linelet3D::get_band(fltarray& Transf, fltarray& Band, int NumBand) 
{
    int Sx = size_band_nx(NumBand);
    int Sy = size_band_ny(NumBand);
    int Sz = size_band_nz(NumBand);
    int Depx = pos_band_nx(NumBand);
    int Depy = pos_band_ny(NumBand);
    int Depz = pos_band_nz(NumBand);
    int i,j,k;
    Band.alloc(Sx,Sy,Sz);

    // printf("BAND %2d: Size = (%2d,%2d,%2d), Dep = (%2d,%2d,%2d) \n",
    //    NumBand+1,Sx,Sy,Sz,Depx,Depy,Depz);
    
    for (i=0; i < Sx; i++)
    for (j=0; j < Sy; j++)
    for (k=0; k < Sz; k++) Band(i,j,k) = Transf(Depx+i,Depy+j,Depz+k);
}

/******************************************************************************/
// linelet3d transform and recons
/******************************************************************************/

void Linelet3D::alloc(int NxCube, int NyCube, int NzCube, int BlockSize)
{
	if(AllocClass) dealloc();
   _Nx = NxCube;
   _Ny = NyCube;  
   _Nz = NzCube;  
   _BlockSize = BlockSize;
   set_block_param();
   AllocClass=True;
}

void Linelet3D::transform(fltarray &Cube, fltarray &Transf)
{
   if ((Cube.nx() != _Nx) || (Cube.ny() != _Ny) || (Cube.nz() != _Nz))
   {
       cout << "Error: the Beamlet class is not allocated for this cube size ... " << endl;
       cout << "       Class cube size = " << _Nx << " " << _Ny << " " << _Nz << endl;
       cout << "       Cube size = " << Cube.nx() << " " << Cube.ny() << " " <<  Cube.nz() << endl;
       exit(-1);
   }
   block_transform(Cube, Transf);
}

void Linelet3D::transform(fltarray &Cube, fltarray &Transf, int BlockSize) 
{
   // linelet3d_check_size_cube (Cube, BlockSize);
   set_attribut (Cube, BlockSize);
   set_block_param();
   block_transform(Cube, Transf);
   //if (_StatInfo == True) info_stat(Transf, SigmaSimu);
}

void Linelet3D::filter(fltarray &Cube, fltarray &Recons, int BlockSize, 
		int s3, float SigmaNoise, float NSigma, filter_type FilterType, float Alpha, int WienerBS, bool force4, bool UseCubeSigma, 
		fltarray *TabSigma, fltarray *CubeSigma, dblarray * TabStat, int NbrScale3D)
{
	if(Verbose) cerr<<"LIN::BlockFilter()"<<endl;
	set_attribut (Cube, BlockSize);
	set_block_param();
	Recons.alloc(Cube.nx(),Cube.ny(),Cube.nz());
	
	if(TabStat!=NULL) 
	{
		TabStat[s3].alloc(_NbrBlock,NbrBand,6);
		if(s3==0) TabStat[NbrScale3D].alloc(NbrScale3D,NbrBand,6);
	}
	
	// Loops on the blocks
	#if USE_OMP_L
	omp_set_nested(1);
	#pragma omp parallel for 
	#endif
	for (int xblock=0; xblock<_Nxb; xblock++)
	#if USE_OMP_L
	#pragma omp parallel for 
	#endif
	for (int yblock=0; yblock<_Nyb; yblock++)
	#if USE_OMP_L
	#pragma omp parallel for 
	#endif
	for (int zblock=0; zblock<_Nzb; zblock++)
	{
//cerr<<"block "<<xblock<<","<<yblock<<","<<zblock<<endl;
		fltarray CurrentBlock (_BlockSize,_BlockSize,_BlockSize);
		fltarray BlockTransf (_BlockSize,_BlockSize,PRad3D.getNbPlane());
		B3D.get_block_cube (xblock, yblock, zblock, Cube, CurrentBlock);

		transform_one_block (CurrentBlock, BlockTransf);

		if(TabStat!=NULL)
		{
			float val;
			int n=num_block(xblock, yblock, zblock);
			int depth = 3*_BlockSize*_BlockSize;
			for(int s=0;s<NbrBand;s++)
			{
				fltarray temp(s2d_xn(s),s2d_yn(s),depth);
				// Normalization and recopy
				for(int i = 0; i < s2d_xn(s); i++)
					for(int j = 0; j < s2d_yn(s); j++)
						for(int k = 0; k < depth; k++)
						{
							if(UseCubeSigma) val = BlockTransf(i+s2d_x0(s),j+s2d_y0(s),k)/(CubeSigma[s3])(i+s2d_x0(s),j+s2d_y0(s),k);
							else val = BlockTransf(i+s2d_x0(s),j+s2d_y0(s),k)/(*TabSigma)(s3,s);
							temp(i,j,k)=val;
						}
				double Mean, Sig, Skew, Curt;
				float Mini, Maxi;
				moment4(temp.buffer(), temp.n_elem(), Mean, Sig, Skew, Curt, Mini, Maxi, true); // not centered
				TabStat[s3](n,s,0)=Mean;
				TabStat[s3](n,s,4)=Mini;
				TabStat[s3](n,s,5)=Maxi;
				int Ntot = depth * s2d_xn(s)*s2d_yn(s) ;
				moment4_center(Ntot, Mean, Sig, Skew, Curt, TabStat[s3](n,s,1), TabStat[s3](n,s,2), TabStat[s3](n,s,3));

				if(Verbose)
				{
					cerr<<" ("<<s2d_x0(s)<<":"<<s2d_x0(s)+s2d_xn(s)-1<<")x("<<s2d_y0(s)<<":"<<s2d_y0(s)+s2d_yn(s)-1<<")";
					cerr<< "\tStat echelle ("<<s3<<","<<s<<") : (mean,std) = ("<<Mean<<", \t"<<Sig<<")"<<endl;
					cerr << " MaxCoef Value = "<<Mini<<","<<Maxi << endl;
				}

				// Global stats
				TabStat[NbrScale3D](s3,s,0)+=Mean;
				TabStat[NbrScale3D](s3,s,1)+=Sig;
				TabStat[NbrScale3D](s3,s,2)+=Skew;
				TabStat[NbrScale3D](s3,s,3)+=Curt;
				if(n==0)
				{
					TabStat[NbrScale3D](s3,s,4) = Mini;
					TabStat[NbrScale3D](s3,s,5) = Maxi;
				}
				else
				{
					TabStat[NbrScale3D](s3,s,4) = TabStat[NbrScale3D](s3,s,4) < Mini ? TabStat[NbrScale3D](s3,s,4) : Mini;
					TabStat[NbrScale3D](s3,s,5) = TabStat[NbrScale3D](s3,s,5) > Maxi ? TabStat[NbrScale3D](s3,s,5) : Maxi;
				}
			}
		}

		if(SigmaNoise>0)
			if(FilterType==FT_HARD || FilterType==FT_SOFT) threshold_single(BlockTransf, BlockSize, s3, SigmaNoise, NSigma, FilterType, 
				force4, UseCubeSigma, TabSigma, CubeSigma);
			else if(FilterType==FT_WIENER) wiener_single(BlockTransf, BlockSize, s3, SigmaNoise*NSigma/3., WienerBS, 
				UseCubeSigma, TabSigma, CubeSigma);
			else cerr<<"Filtering method '"<<string_filter_type(FilterType)<<"'not implemented yet";

		if(!no_recons)
		{
			if(keep_energy)
			{
				if(FilterType==FT_WIENER) // if wiener, hard at 1 sigma
					threshold_single(BlockTransf, BlockSize, s3, SigmaNoise, 1., FT_HARD, force4, UseCubeSigma, TabSigma, CubeSigma);
				fltarray BlockTransfNew (_BlockSize,_BlockSize,PRad3D.getNbPlane());
				recons_one_block (BlockTransf, CurrentBlock);
				transform_one_block (CurrentBlock, BlockTransfNew);

				float lvl;
				for(int s=0;s<nbr_band();s++)
					for(int i = s2d_x0(s); i < s2d_x0(s)+s2d_xn(s); i++)
						for(int j = s2d_y0(s); j < s2d_y0(s)+s2d_yn(s); j++)
							for(int k = 0; k < BlockTransf.nz(); k++)
							{
								lvl = SigmaNoise * (*TabSigma)(s3,s) ; // sigma_noise
								
								float max_amplif = 7.;
								if ( abs(BlockTransf(i,j,k)) > lvl/2. )
								{
									if( abs(BlockTransfNew(i,j,k)) > abs(BlockTransf(i,j,k)/max_amplif) )
										BlockTransf(i,j,k) = BlockTransf(i,j,k)*BlockTransf(i,j,k)/BlockTransfNew(i,j,k);
									//else BlockTransf(i,j,k) *= 10.;
								}
							}
			}
		
			recons_one_block (BlockTransf, CurrentBlock);

			if (_BlockOverlap) 
				B3D.add_block_cube(xblock, yblock, zblock, Recons, CurrentBlock,True);
			else              
				B3D.put_block_cube(xblock, yblock, zblock, Recons, CurrentBlock);
		}
	}
	if(Verbose) cerr<<"end LIN::BlockFilter()"<<endl;
}



void Linelet3D::recons(fltarray& Transf, fltarray& Cube) {
   
   set_recons_block_param();
   block_recons (Transf, Cube);
}


void Linelet3D::block_transform (fltarray& Cube, fltarray& Linelet)
{
	if (_NbrBlock==1)
	{
		fltarray CurrentBlock (_BlockSize,_BlockSize,_BlockSize);   
		Linelet.alloc (_BlockSize,_BlockSize,_LinNbPlane);

		B3D.get_block_cube (0, 0, 0, Cube, CurrentBlock);
		transform_one_block (CurrentBlock, Linelet);
		//sim_transform_one_block (CurrentBlock, BlockTransf);

		if (_StatInfo)
		stat_add(0, 0, 0, Linelet);	        

	}
	else
	{
		Linelet.alloc (_BlockSize,_BlockSize,_LinNbPlane);
		int CpmtBlck=0;
		
		// Loops on the block
		#if USE_OMP_L
		#pragma omp parallel for 
		#endif
		for (int xblock=0; xblock<_Nxb; xblock++)
			#if USE_OMP_L
			#pragma omp parallel for 
			#endif
			for (int yblock=0; yblock<_Nyb; yblock++)
				#if USE_OMP_L
				#pragma omp parallel for 
				#endif
				for (int zblock=0; zblock<_Nzb; zblock++)
				{
					fltarray LOCALCurrentBlock (_BlockSize,_BlockSize,_BlockSize);
					fltarray LOCALBlockTransf (_BlockSize,_BlockSize,PRad3D.getNbPlane());
					CpmtBlck++;

					B3D.get_block_cube (xblock, yblock, zblock, Cube, LOCALCurrentBlock);

					transform_one_block (LOCALCurrentBlock, LOCALBlockTransf);

					put_block_trans (xblock, yblock, zblock, Linelet, LOCALBlockTransf);
					
					if (_StatInfo)
						stat_add(xblock, yblock, zblock, LOCALBlockTransf);
				}
	}
}



//		#if USE_OMP_L
//		#pragma omp parallel for 
//		#endif
void Linelet3D::block_filter (fltarray& Cube, fltarray& Recons)
{
} 

int Linelet3D::s2d_y0(int s)
{
	int NbPoint=_BlockSize;
	int ScaleNbPt,BegInd;
	
	// previous scales
	for(int ss=0;ss<s/3;ss++)
	{
		//pos_band_nx size_band_nx
		ScaleNbPt = NbPoint/2;
		BegInd    = NbPoint-ScaleNbPt;
		NbPoint = BegInd;
	}

	// if not the last scale
	if( s != 3*_NbrScale-3 )
	{
		ScaleNbPt = NbPoint/2;
		BegInd    = NbPoint-ScaleNbPt;
		if(s%3==0) return BegInd;
		else if(s%3==1) return BegInd;
		else return 0;//(s%3==2)
	}
	else // last scale
		return 0;
}

int Linelet3D::s2d_yn(int s)
{
	int NbPoint=_BlockSize;
	int ScaleNbPt,BegInd;

	// previous scales
	for(int ss=0;ss<s/3;ss++)
	{
		//pos_band_nx size_band_nx
		ScaleNbPt = NbPoint/2;
		BegInd    = NbPoint-ScaleNbPt;
		NbPoint = BegInd;
	}

	// if not the last scale
	if( s != 3*_NbrScale-3 )
	{
		ScaleNbPt = NbPoint/2;
		BegInd    = NbPoint-ScaleNbPt;
		if(s%3==0) return ScaleNbPt;
		else if(s%3==1) return ScaleNbPt;
		else return BegInd;//(s%3==2)
	}
	else // last scale
		return BegInd;
}
int Linelet3D::s2d_x0(int s)
{
	int NbPoint=_BlockSize;
	int ScaleNbPt,BegInd;

	// previous scales
	for(int ss=0;ss<s/3;ss++)
	{
		//pos_band_nx size_band_nx
		ScaleNbPt = NbPoint/2;
		BegInd    = NbPoint-ScaleNbPt;
		NbPoint = BegInd;
	}

	// if not the last scale
	if( s != 3*_NbrScale-3 )
	{
		ScaleNbPt = NbPoint/2;
		BegInd    = NbPoint-ScaleNbPt;
		if(s%3==0) return BegInd;
		else if(s%3==1) return 0;
		else return BegInd;//(s%3==2)
	}
	else // last scale
		return 0;
}

int Linelet3D::s2d_xn(int s)
{
	int NbPoint=_BlockSize;
	int ScaleNbPt,BegInd;

	// previous scales
	for(int ss=0;ss<s/3;ss++)
	{
		//pos_band_nx size_band_nx
		ScaleNbPt = NbPoint/2;
		BegInd    = NbPoint-ScaleNbPt;
		NbPoint = BegInd;
	}

	// if not the last scale
	if( s != 3*_NbrScale-3 )
	{
		ScaleNbPt = NbPoint/2;
		BegInd    = NbPoint-ScaleNbPt;
		if(s%3==0) return ScaleNbPt;
		else if(s%3==1) return BegInd;
		else return ScaleNbPt;//(s%3==2)
	}
	else // last scale
		return BegInd;
}

void Linelet3D::threshold_single(fltarray &Band, int BlockSize, int s3, float SigmaNoise, float NSigma, filter_type FilterType, 
		bool force4, bool UseCubeSigma, fltarray *TabSigma, fltarray *CubeSigma)
{
	if(Verbose) cerr<<"Linelet3D::threshold_single(.,"<<SigmaNoise<<","<<NSigma<<","<<FilterType<<","<<force4<<","<<UseCubeSigma<<")"<<endl;

	float lvl;
	for(int s=0;s<nbr_band();s++)
	{
		float cnt=0,cnt2=0;
		for(int i = s2d_x0(s); i < s2d_x0(s)+s2d_xn(s); i++)
			for(int j = s2d_y0(s); j < s2d_y0(s)+s2d_yn(s); j++)
				for(int k = 0; k < Band.nz(); k++)
				{
					float Nsig=NSigma;
					if(s3==0) Nsig=(force4 ? (NSigma+1) : NSigma);
					if(UseCubeSigma) lvl = SigmaNoise * Nsig * (*CubeSigma)(i,j,k);
					else lvl = SigmaNoise * Nsig * (*TabSigma)(s3,s);

					cnt+=1;
					if( abs(Band(i,j,k)) < lvl )
					{// Hard and soft put theese to 0
						cnt2+=1;
						Band(i,j,k)=0;
					}// soft lowers the remaining by lvl
					else if(FilterType==FT_SOFT) Band(i,j,k) -= (2*int( Band(i,j,k)>0 )-1)*lvl;
				}
		if(Verbose) cerr<<"Band "<<s3<<","<<s<<") n#proportion non seuillee ("<<lvl<<")="<<cnt-cnt2<<"#"<<(cnt-cnt2)/cnt<<endl;
	}
	
	if(Verbose) cerr<<"...End BCurvelet3D::threshold_single"<<endl;
}

int Linelet3D::get_BC_pos(int kx,int ky,int bx,int by,int LocalBS,int N)
{
	int Px = ( bx + kx*LocalBS ) % (2*N) ;
	int Py = ( by + ky*LocalBS ) ;
	int s,f; // slow and fast components in terms of exterior and interior loop
	
	if(Py >= N)
	{
		if(Px < N)
		{
			f = Px ;
			s = (2*N - 1 - Py)*N ;
		}
		else // Px >= N
		{
			f = 2*N - 1 - Py ;
			s = (2*N - 1 - Px)*N ;
		}
	//cerr<<" kv kh by bx -> h v y = "<<kv<<","<<kh<<","<<by<<","<<bx<<" -> "<<h<<" "<<v<<" "<<y<<"   "<<endl;
	}
	else 
	{
		if(Px < N)
		{
			f = Py ;
			s = (Px + 2*N)*N ;
		}
		else // Px >= N
		{
			f = Py ;
			s = (3*N - Px)*N - 1 ;
		}
	}
	return s+f;
}
void Linelet3D::wiener_single(fltarray &Band, int BlockSize, int s3, float noise_lvl, int LocalBS, 
		bool UseCubeSigma, fltarray *TabSigma, fltarray *CubeSigma)
{
	if(Verbose) cerr<<"Lin3D::wiener("<<BlockSize<<","<<s3<<","<<noise_lvl<<","<<LocalBS<<","<<UseCubeSigma<<")..."<<endl;

	float val;
	float noise2 = noise_lvl*noise_lvl;
	

	int N = BlockSize; // for height : number of angles 

	// Linelet scales
	for(int s=0;s<nbr_band();s++)
	{
		int Nax = ceil(float(2*N)/float(LocalBS)); // horizontal angle
		int Nay = ceil(float(N)/float(LocalBS)); // vertical angle
		int Npx = ceil(float(s2d_xn(s))/float(LocalBS)); // horizontal space position
		int Npy = ceil(float(s2d_yn(s))/float(LocalBS)); // vertical space position
		//cerr<<"s N nax Nay xn yn Npx Npy = ("<<s<<","<<N<<","<<Nax<<","<<Nay<<","<<s2d_xn(s)<<","<<s2d_yn(s)<<","<<Npx<<","<<Npy<<")"<<endl;

		// Spatial blocks B
			fltarray coef_wiener(s2d_xn(s),s2d_yn(s),3*N*N); // 3*N² angles and s2d_xn(s)*s2d_yn(s) positions

			coef_wiener.init(-2);

			// Wiener blocks (angle^2 and space^2) for mainly horizontal Beamlets
			for(int ky=0 ; ky < Nay ; ky++)
				for(int kx=0 ; kx < Nax ; kx++)
					for(int kpy = 0; kpy < Npy; kpy++)
						for(int kpx = 0; kpx < Npx; kpx++)
						{
							// cerr<<endl<<"  B ky kx kpy kpx = "<<B<<" ("<<ky<<","<<kx<<","<<kpy<<","<<kpx<<") z ";
							double sigma2 = 0.0;
							float cnt = 0;

						// Sigma calculation
							// Pixels in a wiener block (angle)
							for(int by = 0; by < LocalBS; by++)
								for(int bx = 0; bx < LocalBS; bx++)
								{
									// angular position
									int z = get_BC_pos(kx,ky,bx,by,LocalBS,N);

									// pixels in a wiener block (space)
									for(int bpy = 0; bpy < LocalBS; bpy++)
										for(int bpx = 0; bpx < LocalBS; bpx++)
										{
											// if(B==0) cerr<<bx<<" "<<by<<" "<<bpx<<" "<<bpy<<endl;
											cnt+=1;
											// spatial position = position in the block + block_position
											int x = ((bpx + kpx*LocalBS) % s2d_xn(s)) ;
											int y = ((bpy + kpy*LocalBS) % s2d_yn(s)) ;
											if(UseCubeSigma)	val = Band(x+s2d_x0(s), y+s2d_y0(s), z) 
																  / (*CubeSigma)(x+s2d_x0(s), y+s2d_y0(s), z);
											else 				val = Band(x+s2d_x0(s), y+s2d_y0(s), z) / (*TabSigma)(s3,s);

											sigma2+=pow(val,2);
										}
								}
							float sig2 = max( 0.0, sigma2/cnt - noise2 );
							float norm = sig2 / (sig2+noise2);

						// Store the coef in the table
							for(int by = 0; by < LocalBS; by++)
								for(int bx = 0; bx < LocalBS; bx++)
								{
									int z = get_BC_pos(kx,ky,bx,by,LocalBS,N);

									// pixels in a wiener block (space)
									for(int bpy = 0; bpy < LocalBS; bpy++)
										for(int bpx = 0; bpx < LocalBS; bpx++)
										{
											// spatial position = position in the block + block_position
											int x = ((bpx + kpx*LocalBS) % s2d_xn(s)) ;
											int y = ((bpy + kpy*LocalBS) % s2d_yn(s)) ;
											if( coef_wiener(x,y,z) < -1 )
												coef_wiener(x,y,z) = norm;
										}
								}
						}// end wiener blocks

			// Wiener blocks (angle^2 and space^2) for mainly vertical Beamlets
			for(int ky=0 ; ky < Nay ; ky++)
				for(int kx=0 ; kx < Nay ; kx++) // Nay becuse it's a square
					for(int kpy = 0; kpy < Npy; kpy++)
						for(int kpx = 0; kpx < Npx; kpx++)
						{
							double sigma2 = 0.0;
							float cnt = 0;

						// Sigma calculation
							// Pixels in a wiener block (angle)
							for(int by = 0; by < LocalBS; by++)
								for(int bx = 0; bx < LocalBS; bx++)
								{
									// angular position
									int z = ((by + ky*LocalBS) % (Nay*LocalBS))*N + ((bx + kx*LocalBS) % (Nay*LocalBS));
									//z=0;
									// pixels in a wiener block (space)
									for(int bpy = 0; bpy < LocalBS; bpy++)
										for(int bpx = 0; bpx < LocalBS; bpx++)
										{
											cnt+=1;
											// spatial position = position in the block + block_position
											int x = ((bpx + kpx*LocalBS) % s2d_xn(s)) ;
											int y = ((bpy + kpy*LocalBS) % s2d_yn(s)) ;
											if(UseCubeSigma)	val = Band(x+s2d_x0(s), y+s2d_y0(s), z) 
																  / (*CubeSigma)(x+s2d_x0(s), y+s2d_y0(s), z);
											else				val = Band(x+s2d_x0(s), y+s2d_y0(s), z) / (*TabSigma)(s3,s);

											sigma2+=pow(val,2);
										}
								}
							float sig2 = max( 0.0, sigma2/cnt - noise2 );
							float norm = sig2 / (sig2+noise2);

						// Store the coef in the table
							for(int by = 0; by < LocalBS; by++)
								for(int bx = 0; bx < LocalBS; bx++)
								{
									int z = ((by + ky*LocalBS) % (Nay*LocalBS))*N + ((bx + kx*LocalBS) % (Nay*LocalBS));

									// pixels in a wiener block (space)
									for(int bpy = 0; bpy < LocalBS; bpy++)
										for(int bpx = 0; bpx < LocalBS; bpx++)
										{
											// spatial position = position in the block + block_position
											int x = ((bpx + kpx*LocalBS) % s2d_xn(s)) ;
											int y = ((bpy + kpy*LocalBS) % s2d_yn(s)) ;
											if( coef_wiener(x,y,z) < -1 )
												coef_wiener(x,y,z) = norm;
										}
								}
						}// end wiener blocks

//		char filename[64];
//		sprintf(filename,"coef_w_s%d_b%d.fits",s,B);
//		writefltarr(filename, coef_wiener);

		// Apply the wiener coefficient
			// Wiener blocks (int angle2 and space) for mainly horizontal beamlets
			for(int ky=0 ; ky < Nay ; ky++)
				for(int kx=0 ; kx < Nax ; kx++)
					for(int kpy = 0; kpy < Npy; kpy++)
						for(int kpx = 0; kpx < Npx; kpx++)
							for(int by = 0; by < LocalBS; by++)
								for(int bx = 0; bx < LocalBS; bx++)
								{
									// angular position
									int z = get_BC_pos(kx,ky,bx,by,LocalBS,N);

									// pixels in a wiener block (space)
									for(int bpy = 0; bpy < LocalBS; bpy++)
										for(int bpx = 0; bpx < LocalBS; bpx++)
										{
											// spatial position = position in the block + block_position
											int x = ((bpx + kpx*LocalBS) % s2d_xn(s)) ;
											int y = ((bpy + kpy*LocalBS) % s2d_yn(s)) ;
											if( coef_wiener(x,y,z) > -1 )
											{
												Band(x+s2d_x0(s), y+s2d_y0(s), z) *= coef_wiener(x,y,z);
												coef_wiener(x,y,z) = -2;
											}
										}
								}

			// Wiener blocks (int angle2 and space) for mainly vertical bebamlets
			for(int ky=0 ; ky < Nay ; ky++)
				for(int kx=0 ; kx < Nay ; kx++)
					for(int kpy = 0; kpy < Npy; kpy++)
						for(int kpx = 0; kpx < Npx; kpx++)
							for(int by = 0; by < LocalBS; by++)
								for(int bx = 0; bx < LocalBS; bx++)
								{
									// angular position
									int z = ((by + ky*LocalBS) % (Nay*LocalBS))*N + ((bx + kx*LocalBS) % (Nay*LocalBS));

									// pixels in a wiener block (space)
									for(int bpy = 0; bpy < LocalBS; bpy++)
										for(int bpx = 0; bpx < LocalBS; bpx++)
										{
											// spatial position = position in the block + block_position
											int x = ((bpx + kpx*LocalBS) % s2d_xn(s)) ;
											int y = ((bpy + kpy*LocalBS) % s2d_yn(s)) ;
											if( coef_wiener(x,y,z) > -1 )
											{
												Band(x+s2d_x0(s), y+s2d_y0(s), z) *= coef_wiener(x,y,z);
												coef_wiener(x,y,z) = -2;
											}
										}
								}
	}// end 2D scale
			
	if(Verbose) cerr<<"...End Lin3D::wiener"<<endl;	
}

void Linelet3D::block_recons (fltarray& Transf, fltarray& Cube)
{
// Loop on the blocks
	#if USE_OMP_L
	#pragma omp parallel for 
	#endif
	for (int xblock=0; xblock<_Nxb; xblock++)
		#if USE_OMP_L
		#pragma omp parallel for 
		#endif
		for (int yblock=0; yblock<_Nyb; yblock++)
			#if USE_OMP_L
			#pragma omp parallel for 
			#endif
			for (int zblock=0; zblock<_Nzb; zblock++)
			{
				fltarray LOCALCurrentBlock (_BlockSize,_BlockSize,_BlockSize);
				fltarray LOCALBlockTransf (_BlockSize,_BlockSize,PRad3D.getNbPlane());

				get_block_trans(xblock, yblock, zblock, Transf,  LOCALBlockTransf);
				recons_one_block (LOCALBlockTransf, LOCALCurrentBlock);

				if (_BlockOverlap) 
					B3D.add_block_cube(xblock, yblock, zblock, Cube, LOCALCurrentBlock,True);
				else              
					B3D.put_block_cube(xblock, yblock, zblock, Cube, LOCALCurrentBlock);
			}
} 

void Linelet3D::sim_transform_one_block (fltarray& Cube, fltarray& Transf) {

   for (int curPlane=0;curPlane<_RadNbPlane;curPlane=curPlane+3*_BlockSize)
   for (int i=0;i<_BlockSize;i++)
   for (int j=0;j<_BlockSize;j++) 
   Transf(i,j,curPlane) = Cube(i,j,curPlane/3/_BlockSize);
}

void Linelet3D::transform_one_block (fltarray& Cube, fltarray& Transf) 
{
	switch (lin3d_class(_LinTransf))
	{
		case LIN3D_CL_ORTHO:
		{
			fltarray _PRadTransf_local(_BlockSize, _BlockSize, _RadNbPlane);

			PartialRadon3D PRad3D_local;
			PRad3D_local.setControl(_Control);
			PRad3D_local.setWriteFile(_WriteFile);
			PRad3D_local.setStatInfo(_StatInfo);
			PRad3D_local.setTakeOnlyCompletePlane(_TakeOnlyCompletePlane);
			PRad3D_local.setDoNotCorrectSigma(_DoNotCorrectSigma);   
			PRad3D_local.setParam(_BlockSize,_NbrScale);         //!!!!!add type of WT...!!! 
			PRad3D_local.alloc();
			PRad3D_local.transform (Cube, _PRadTransf_local);

			if (_Control) fits_write_fltarr((char*)"xx_rad.fits", _PRadTransf_local);
			if (_ptrSB1D == (SubBand1D*)NULL) {
				cout << "Error: filter bank class not initialized ... " << endl; exit(-1);}

			Ifloat Plane (_BlockSize,_BlockSize);
			Ifloat OWTPlane (_BlockSize,_BlockSize);
			Ortho_2D_WT WT2D (*_ptrSB1D);
			for (int curPlane=0;curPlane<_RadNbPlane;curPlane++)
			{

				for (int i=0;i<_BlockSize;i++)
					for (int j=0;j<_BlockSize;j++) 
						Plane(i,j)=_PRadTransf_local(i,j,curPlane);

				WT2D.transform (Plane, OWTPlane, _NbrScale);

				for (int i=0;i<_BlockSize;i++)
					for (int j=0;j<_BlockSize;j++) 
						Transf(i,j,curPlane) = OWTPlane(i,j);
			}
			break;
		}

		case LIN3D_CL_PYR:
			cout << "Not yet implemented ..." << endl;
			break;
			
		default:
			break;
	}
}
      
void Linelet3D::sim_recons_one_block (fltarray& Transf, fltarray& Cube)
{
	for (int curPlane=0;curPlane<_RadNbPlane;curPlane=curPlane+3*_BlockSize)
	for (int i=0;i<_BlockSize;i++)
	for (int j=0;j<_BlockSize;j++) 
		Cube(i,j,curPlane/3/_BlockSize) = Transf(i,j,curPlane);
}   
 
void Linelet3D::recons_one_block (fltarray& Transf, fltarray& Cube)
{	
	switch (lin3d_class(_LinTransf))
	{
		case LIN3D_CL_ORTHO:
		{
			fltarray _PRadTransf_local(_BlockSize, _BlockSize, _RadNbPlane);

			PartialRadon3D PRad3D_local;
			PRad3D_local.setControl(_Control);
			PRad3D_local.setWriteFile(_WriteFile);
			PRad3D_local.setStatInfo(_StatInfo);
			PRad3D_local.setTakeOnlyCompletePlane(_TakeOnlyCompletePlane);
			PRad3D_local.setDoNotCorrectSigma(_DoNotCorrectSigma);   
			PRad3D_local.setParam(_BlockSize,_NbrScale);         //!!!!!add type of WT...!!! 
			PRad3D_local.alloc();
			PRad3D_local.transform (Cube, _PRadTransf_local);

			Ifloat Plane (_BlockSize,_BlockSize);
			Ifloat OWTPlane (_BlockSize,_BlockSize);
			Ortho_2D_WT WT2D (*_ptrSB1D);
			for (int curPlane=0;curPlane<_RadNbPlane;curPlane++)
			{

				for (int i=0;i<_BlockSize;i++)
				for (int j=0;j<_BlockSize;j++) 
				Plane(i,j)=Transf(i,j,curPlane);

				WT2D.recons (Plane, OWTPlane, _NbrScale);

				for (int i=0;i<_BlockSize;i++)
				for (int j=0;j<_BlockSize;j++) 
				_PRadTransf_local(i,j,curPlane) = OWTPlane(i,j);
			}

			if (_Control)
				fits_write_fltarr((char*)"xx_wtinv.fits", _PRadTransf_local);
			
			if (_ptrSB1D == (SubBand1D*)NULL)
			{
				cout << "Error: filter bank class not initialized ... " << endl;
				exit(-1);	 
			}	

			PRad3D_local.recons (_PRadTransf_local, Cube);
			break;   
		}

		case LIN3D_CL_PYR:
			cout << "Not yet implemented ..." << endl;
			break;
			
		default:
				break;
	}
}
                    	
/******************************************************************************/
// read write functions
/******************************************************************************/
void Linelet3D::write (char* Name, fltarray& Cube) {
   
   //fits_write_fltarr(Name, Cube);
  
   mr_io_write (Name, Cube); 
   
   if (_StatInfo) {
      char FileName[80];
      for (int s=0;s<_NbrScale;s++) {
         sprintf (FileName,"%s_scale_%d", Name, s);
	 fits_write_fltarr(FileName, _pScale[s]);
      }
   }
}


/******************************************************************************/
// get_block_cube, get_block_trans
/******************************************************************************/

void Linelet3D::get_block_trans (int Bi, int Bj, int Bk, fltarray& Transf, 
                                 fltarray& CurrentBlock) {
				 
   int NumBlock = num_block(Bi,Bj,Bk);
   
   switch (lin3d_class(_LinTransf)) {
   
      case LIN3D_CL_ORTHO: {   				

         for (int k=0; k<_BlockSize; k++)
         for (int l=0; l<_BlockSize; l++)
         for (int m=0; m<_RadNbPlane; m++) {
        
            CurrentBlock(k,l,m) = Transf(k,l,m + NumBlock*_RadNbPlane);  
         }
	 break;
      }
	
      case LIN3D_CL_PYR: cout << "Not yet implemented ..." << endl; break;
      default: break;	 
   }   
   
   if (_WriteFile) {
      char FileName[50];
      sprintf (FileName, "getBlkTrans_%d.fits", num_block(Bi,Bj,Bk));
      fits_write_fltarr (FileName, CurrentBlock);
   }				
}  
  
/******************************************************************************/
// put_block_trans, 
/******************************************************************************/
void Linelet3D::put_block_trans(int Bi, int Bj, int Bk,  fltarray& Transf,  
                                fltarray& CurrentBlock) {

   int NumBlock = num_block(Bi,Bj,Bk);
   
   switch (lin3d_class(_LinTransf)) {
   
      case LIN3D_CL_ORTHO: {   				

         for (int k=0; k<_BlockSize; k++)
         for (int l=0; l<_BlockSize; l++)
         for (int m=0; m<_RadNbPlane; m++) {
        
            Transf(k,l,m + NumBlock*_RadNbPlane) = CurrentBlock(k,l,m);  
         }
	 break;
      }
	
      case LIN3D_CL_PYR: cout << "Not yet implemented ..." << endl; break;
      default: break;	 
   }   
   
   if (_WriteFile) {
      char FileName[50];
      sprintf (FileName, "putBlkTrans_%d.fits", NumBlock);
      fits_write_fltarr (FileName, CurrentBlock);
   }											
}
  
  
static void PrintError( int status)
{
    /*****************************************************/
    /* Print out cfitsio error messages and exit program */
    /*****************************************************/

    char status_str[FLEN_STATUS], errmsg[FLEN_ERRMSG];
  
    if (status)
      fprintf(stderr, "\n*** Error occurred during program execution ***\n");

    ffgerr(status, status_str);        /* get the error status description */
    fprintf(stderr, "\nstatus = %d: %s\n", status, status_str);

    if ( ffgmsg(errmsg) )  /* get first message; null if stack is empty */
    {
         fprintf(stderr, "\nError message stack:\n");
         fprintf(stderr, " %s\n", errmsg);

         while ( ffgmsg(errmsg) )  /* get remaining messages */
             fprintf(stderr, " %s\n", errmsg);
    }

    exit( status );       /* terminate the program, returning error status */
}  
static void mr_io_name (char *File_Name_In, char *File_Name_Out)
{
     int L;

    strcpy (File_Name_Out, File_Name_In);

    L = strlen (File_Name_In);
    if ((L < 4) || (File_Name_In[L-1] != 't')
                || (File_Name_In[L-2] != 'e')
                || (File_Name_In[L-3] != 'b')
                || (File_Name_In[L-4] != '.'))
    {
        strcat (File_Name_Out, ".bet");
    }
} 

void Linelet3D::mr_io_write(char* Name, fltarray& Cube) {
    
  
  //---------------------------------------------------------
  // take in Ridgelet3d (with minor modifications) !!!!!!!!!!
  //--------------------------------------------------------- 
 char filename[256];
 fitsfile *fptr;    
 int status;
 //int i,j,s;
 //float *Ptr;
 int simple;
 int bitpix;
 long naxis=0;
 long naxes[3];
 long nelements;
 long group = 1;  /* group to write in the fits file, 1= first group */
 long firstpixel = 1;    /* first pixel to write (begin with 1) */
 Ifloat Aux;
 //long fpixels[3];
 //long lpixels[3];

/* we keep mr as extension even if its fits ! */
 mr_io_name (Name, filename);

#if DEBUG_IO  
    cout << "Write on " << filename << endl;
#endif

 FILE *FEXIST = fopen(filename, "rb");
 if (FEXIST)
 {
    fclose(FEXIST);
    remove(filename);               /* Delete old file if it already exists */
 }

 status = 0;         /* initialize status before calling fitsio routines */

    /* open the file */
 if ( ffinit(&fptr, filename, &status) )     /* create the new FITS file */
     PrintError( status );           /* call PrintError if error occurs */
                                                                              
/* write  the header */
 simple   = True;
 bitpix   =  -32;   /* 32-bit real pixel values      */
 long pcount   =   0;  /* no group parameters */
 long gcount   =   1;  /* only a single image/group */
 int  extend   =   False;
 
 naxis = 3;
 naxes[0] = _BlockSize; 
 naxes[1] = _BlockSize;
 naxes[2] = _LinNbPlane;
 
 // write first header part (parameters)
 if ( ffphpr(fptr,simple,bitpix,naxis,naxes,pcount,gcount,extend,&status) )
     PrintError( status );          /* call PrintError if error occurs */

  // write the header of the multiresolution file
  mr_io_fill_header(fptr);

  nelements = naxes[0] * naxes[1] * naxes[2];
  if ( ffppre(fptr, group, firstpixel, nelements, Cube.buffer(), &status) )
              PrintError( status );  
 
 /* close the FITS file */
 if ( ffclos(fptr, &status) )  PrintError( status );  
// cout << " end of write fits " << endl;
  
} 
   
void Linelet3D::mr_io_fill_header(fitsfile *fptr) {

  int status = 0; // this means OK for cfitsio !!!
    /*****************************************************/
     /* write optional keyword to the header */
    /*****************************************************/
  if ( ffpkyj(fptr, (char*)"Nx", (long)_Nx,(char*)"Nx",&status))
     PrintError( status );  
  if ( ffpkyj(fptr,(char*)"Ny",(long)_Ny,(char*)"Ny",&status))
     PrintError( status );  
  if ( ffpkyj(fptr,(char*)"Nz",(long)_Nz,(char*)"Nz",&status))
     PrintError( status ); 
  if ( ffpkyj(fptr,(char*)"BlockSiz",(long)_BlockSize,(char*)"BlockSize",&status))
     PrintError( status ); 
  if ( ffpkyj(fptr, (char*)"Nxb", (long)_Nxb,(char*)"Number of blocks Nxb",&status))
     PrintError( status );  
  if ( ffpkyj(fptr,(char*)"Nyb",(long)_Nyb,(char*)"Number of blocks Nyb",&status))
     PrintError( status );  
  if ( ffpkyj(fptr,(char*)"Nzb",(long)_Nzb,(char*)"Number of blocks Nzb",&status))
     PrintError( status ); 
  int Val = (_BlockOverlap == True) ? 1: 0;
  if ( ffpkyj(fptr,(char*)"Overlap",(long) Val,(char*)"Block overlap",&status))
     PrintError( status );
     
     
     
  if ( ffpkyj(fptr,(char*)"NbrScale",(long)_NbrScale,(char*)"Number of scales",&status))
     PrintError( status );      
           
  if ( ffpkyj(fptr, (char*)"Type_Tra", (long)  _LinTransf, 
                          (char*)StringLin3DTransform(_LinTransf), &status))
     PrintError( status );
}


void Linelet3D::mr_io_read (char *Name, fltarray& Cube) {
    // for fits
    char filename[256];
    fitsfile *fptr;           /* pointer to the FITS file */
    int status=0, hdutype ;
    long hdunum;
    char comment[FLEN_COMMENT];
    int naxis;
    long naxes[3];
    long mon_long;
    int anynul = 0;
    long nulval = 0;
    long inc[3];
    void PrintError( int status);
    long nelements = 0 ; // naxes[0] * naxes[1] in the image
    //long fpixels[3];
    //long int lpixels[3];
    
     // for multiresol
    float *Ptr;
    //int my_logical; // sais pas...

     mr_io_name (Name, filename);

    inc[0]=1;  inc[1]=1; inc[2]=1;
   
    /* open the file */
    status = 0;         /* initialize status before calling fitsio routines */
    if ( ffopen(&fptr, filename, (int) READONLY, &status) ) 
         PrintError( status );
                                    
    hdunum = 1;  /*read  table */
    if ( ffmahd(fptr, hdunum, &hdutype, &status) ) /* move to the HDU */
           PrintError( status );

    int simple, bitpix, extend;
    long pcount, gcount;
    if ( ffghpr(fptr, 3, &simple, &bitpix, &naxis, naxes, &pcount,
            &gcount, &extend, &status)) /* move to the HDU */
           PrintError( status );

    nelements = naxes[0] * naxes[1] * naxes[2];
    // cout << " begin to read " << endl;
    Cube.alloc(naxes[0], naxes[1], naxes[2],(char*)"read io");
     
    if (ffgkyj(fptr,(char*)"Nx", &mon_long, comment, &status)) PrintError( status );
    _Nx = (int) mon_long;  
    if (ffgkyj(fptr,(char*)"Ny", &mon_long, comment, &status)) PrintError( status );
    _Ny = (int) mon_long; 
    if (ffgkyj(fptr,(char*)"Nz", &mon_long, comment, &status)) PrintError( status );
    _Nz = (int) mon_long;  
    if (ffgkyj(fptr,(char*)"BlockSiz", &mon_long, comment, &status)) PrintError( status );
    _BlockSize = (int) mon_long;
    if (ffgkyj(fptr,(char*)"NbrScale", &mon_long, comment, &status)) PrintError( status );
    _NbrScale = (int) mon_long;
    if (ffgkyj(fptr,(char*)"Nxb", &mon_long, comment, &status)) PrintError( status );
    _Nxb = (int) mon_long;  
    if (ffgkyj(fptr,(char*)"Nyb", &mon_long, comment, &status)) PrintError( status );
    _Nyb = (int) mon_long; 
    if (ffgkyj(fptr,(char*)"Nzb", &mon_long, comment, &status)) PrintError( status );
    _Nzb = (int) mon_long;      
    
    if (ffgkyj(fptr,(char*)"Overlap", &mon_long, comment, &status)) PrintError( status );
    _BlockOverlap = (mon_long == 1) ? True: False;
    if (ffgkyj(fptr,(char*)"Type_Tra", &mon_long, comment, &status)) PrintError( status );
    _LinTransf = (type_linelet3d_WTtrans) mon_long;

    Ptr = Cube.buffer();
    if ( ffgpve(fptr, 1, 1, nelements, nulval, Ptr, &anynul, &status))
         PrintError( status );
  
    if ( ffclos(fptr, &status) ) PrintError( status );
    // cout << " end of read fits file " << endl;

    if (_BlockSize <= 0) _BlockOverlap = False;
    B3D.BlockOverlap = _BlockOverlap;
    B3D.Verbose = _Verbose;
    B3D.alloc(_Nx,_Ny,_Nz,_BlockSize);
    _BlockSize = B3D.block_size();
    _BlockOverlap = B3D.BlockOverlap;
    _Nxb = B3D.nbr_block_nx();
    _Nyb = B3D.nbr_block_ny();
    _Nzb = B3D.nbr_block_nz();
}

 

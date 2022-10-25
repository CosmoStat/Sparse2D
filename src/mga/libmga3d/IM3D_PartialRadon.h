/******************************************************************************
**                   Copyright (C) 2001 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Jean-Luc Starck
**
**    Date:  10/10/2001 
**    
**    File:  IM3D_Radon.h
**
**    Modification history :
**
*******************************************************************************
**
**    DESCRIPTION  
**    ----------- 
**                 
**    PARAMETRES    
**    ----------    
** 
**    RESULTS      
**    -------  
**
**
******************************************************************************/

#ifndef _PARTRADON3D_H_
#define _PARTRADON3D_H_

#include "Array.h"
#include "IM_Obj.h"
#include "IM3D_IO.h"
#include "IM_IO.h"
#include "FFTN_1D.h"
#include "FFTN_2D.h"
#include "FFTN_3D.h"
#include "macro.h"


/*****************************************************************************/
class PartialRadon3D; 
typedef int  (PartialRadon3D::*tp_mfunc)(int c1, int c2, 
                                        double a, double b, double c);
                                        
typedef void (PartialRadon3D::*tp_setfunc)(cfarray& InCmpArray,
                                           cfarray& OutCmpArrayint,
                                           int n, int l1, int l2, 
                                           int x, int y, int z, Bool TakeSym);
                                           // Bool TakeSym=False);                                        

class PartialRadon3D {

/*****************************************************************************/
// private attributs
/*****************************************************************************/
private: 
public: 
   int ProjectNbr;            // Number of line in 3d radon transform
   int SizeCube;              // size of data cube  
   Bool Verbose;              // verbose mode
   Bool Control;              // control mode
   Bool WriteFile;            // write internal file
   Bool WriteCoord;           // write ccord <=> Plane
   Bool Alloc;                // alloc control
   Bool Plane;                // flag make plane   
   Bool RemoveEqualPlane;  // do not remove equal plane in PRad transf
   Bool StatInfo;             // Statistic information
   Bool TakeOnlyCompletePlane;// Take only ......
   Bool DoNotCorrectSigma;    // Do not ...
   
   FFTN_3D CFFT3D;            // FFT Class
    
   cfarray F3DCube;           // inter array used in (inv-)partial_slice 
   cfarray FPartRad3DImage;   // inter array used in (inv-)partial_slice 
   
   intarray xPos;             // x position of pixel for angle calculation
   intarray yPos;             // y position of pixel for angle calculation
   intarray zPos;             // z position of pixel for angle calculation
   //intarray DataCount1;       // pixel count for plane
   intarray DataCount2;       // pixel count for cube out
   intarray DataCount3;       // pixel count for cube out   
   fltarray MemPlane;         // memorize trnasf position (l1,l2,plane number)
                              // in function of cube position f(x,y,z)
                              
                              // used only if StatInfo is true
   intarray MemNotUsed;       // number of point not used in each Plane
   fltarray CoefSigmaPlane;   // coeff sigma (sigma real/sigma plane)   
   
/*****************************************************************************/
// private members
/*****************************************************************************/
private: 
public:
   //int radon_nl()  {return ProjectNbr;}
   
   void reset ();           // reset internal param
   
   void make_table();       
   void make_plane();                            
        
   void control_cube_size (fltarray& Cube); 
                            // control sise cube in (allocated, same size)
   void control_pradon_size (fltarray& Cube); 
                            // control sise cube in (power of two, same size)                            
                            
   void control_even_sym (cfarray& CmpArray);
                            // control plane symetrie 
   void control_odd_sym (cfarray& CmpArray);
                            // control plane symetrie                             
                              
   void partial_slice (fltarray &Cube, fltarray &PartRad3DImage);
                            // compute partial slice
                            // 1) 3d FFT transform of data in
                            // 2) compute all plane by center of cube
                            // 3) FFT inv on all plane
   
   void inv_partial_slice (fltarray& PartRad3DImage, fltarray& Cube);
                            // inverse partial slice
                            // 1) 2d FFT on all plane
                            // 2) compute cube from plane
                            // 3) inv 3d FFT 
                                 
   void comp_Plane (cfarray& F3DCube, cfarray& FPartRad3DImage,
                    Bool Inverse=False);
                            // compute all plane =>
                            // for each line by center of cube (theta, phi)
                            // we compute the orthogonal plane
                                 
   virtual void comp_odd_Plane (cfarray& F3DCube, cfarray& FPartRad3DImage,
                        Bool Inverse=False);
                                                                                 
   void comp_InvFFT2D (cfarray& FPartRad3DImage , fltarray& PRadImag);
                            // inverse 2d fft of all plane                                                                

   virtual void comp_FFT2D (fltarray& PartRad3DImage, cfarray& FPartRad3DImage);
                            // 2d fft of all plane
                            
   void comp_FFT3D (fltarray& Cube, cfarray& F3DCube);
                            // 3d fft of cube in
                            
   void comp_InvFFT3D (cfarray& F3DCube, fltarray& Cube);
                            // 3d inv fft                                                                          
                                                                                            
   int double2int (double point);
                            // iround function              
                            

   void even_make_plane (cfarray& OutPlane, cfarray& InPlane,
                    tp_mfunc f1, tp_mfunc f2, tp_mfunc f3, 
                    double a, double b, double c, int NumLine,
                    tp_setfunc set_f, Bool Inverse=False);
   void odd_make_plane (cfarray& OutPlane, cfarray& InPlane,
                    tp_mfunc f1, tp_mfunc f2, tp_mfunc f3, 
                    double a, double b, double c, int NumLine,
                    tp_setfunc set_f, Bool Inverse=False);                    
                    
   Bool even_border_manager (int x, int y, int z, int& xPix, int& yPix, int& zPix);
   
   void comp_xyz (int l1, int l2, int n, int& xPix, int& yPix, int& zPix);
                                                 
   void WriteCmpArr (char* FileName, cfarray&  CmpArray); 
   void InfoCmp (cfarray& F3DCube, string Name=""); 
   void EnergyCtrl  (cfarray& CmpArray);
   
   void direct_make_plane  (cfarray& F3DCube, cfarray& FPartRad3DImage,
                            int n, int l1, int l2, int x, int y, int z,
                            Bool TakeSym=False);
   void inverse_make_plane (cfarray& FPartRad3DImage, cfarray& F3DCube,
                            int n, int l1, int l2, int x, int y, int z,
                            Bool TakeSym=False);
   void memorize_make_plane (cfarray& F3DCube, cfarray& FPartRad3DImage,
                             int n, int l1, int l2, int x, int y, int z,
                             Bool TakeSym=False);
                             
   void remove_equal_Plane ();
                            // remove equal plane, onky for odd plane
                            // called in make_plane for odd size.
                                 
   void stat_alloc();       // alloc all array used in stat
   void stat_init();        // init all array used in stat 
                            // (+MemPlane if odd) 
   void stat_comp_coef (cfarray& FPartRad3DImage);  
                            // compute coeff = sima real/sigma mesured
                            // after FFT3D (before inv FFT2D) 
                            
   void stat_info (fltarray& PartRad3DImage);
                            // compute stat on Radon transform data
                            // mean, sigma.... (mean on al plane)
                            
   float stat_comp_true_sigma (int NbPtsUsed,int NbZero, 
                               float Sigma, float Mean);
                            // compute sigma for point used in a plane
                            // (all point wiht MemPlane(i,j) != -1)                             
                                  
   inline int fx_uXv8w8 (int c1,int c2, double a, double b, double c);                                  
   inline int fy_uXv8w8 (int c1,int c2, double a, double b, double c);
   inline int fz_uXv8w8 (int c1,int c2, double a, double b, double c);

   inline int fx1_uXvXw8 (int c1,int c2, double a, double b, double c);
   inline int fy1_uXvXw8 (int c1,int c2, double a, double b, double c);
   inline int fz1_uXvXw8 (int c1,int c2, double a, double b, double c);        
   
   inline int fx2_uXvXw8 (int c1,int c2, double a, double b, double c);
   inline int fy2_uXvXw8 (int c1,int c2, double a, double b, double c);
   inline int fz2_uXvXw8 (int c1,int c2, double a, double b, double c);          
   
   inline int fx1_uXvXwX (int c1,int c2, double a, double b, double c);
   inline int fy1_uXvXwX (int c1,int c2, double a, double b, double c);
   inline int fz1_uXvXwX (int c1,int c2, double a, double b, double c);     
   
   inline int fx2_uXvXwX (int c1,int c2, double a, double b, double c);
   inline int fy2_uXvXwX (int c1,int c2, double a, double b, double c);
   inline int fz2_uXvXwX (int c1,int c2, double a, double b, double c); 
   
   inline int fx3_uXvXwX (int c1,int c2, double a, double b, double c);
   inline int fy3_uXvXwX (int c1,int c2, double a, double b, double c);
   inline int fz3_uXvXwX (int c1,int c2, double a, double b, double c);    
                       
                                   
/*****************************************************************************/
// public attributs
/*****************************************************************************/
public: 


/***********MemNotUsed******************************************************************/
// public members
/*****************************************************************************/
public: 
   
   int sizecube()    {return SizeCube;}
   int getNbPlane() {return ProjectNbr;}
   void setVerbose (Bool Flag) {Verbose=Flag;}
   void setControl (Bool Flag) {Control=Flag;}
   void setWriteFile (Bool Flag) {WriteFile=Flag;}
   void setWriteCoord (Bool Flag) {WriteCoord=Flag;}   
   void setParam (int BlockSize, int NbrScale=0);
   void setRemovePlane (Bool Flag) {RemoveEqualPlane=Flag;}
   void setStatInfo (Bool Flag) {StatInfo=Flag;}
   void setTakeOnlyCompletePlane(Bool Flag) {TakeOnlyCompletePlane=Flag;}
   void setDoNotCorrectSigma(Bool Flag) {DoNotCorrectSigma=Flag;}
   fltarray* getCoefSigma () {return &CoefSigmaPlane;}
   intarray* getMemNotUsed () {return &MemNotUsed;}
   
   PartialRadon3D () {reset();}
   ~PartialRadon3D() {reset();}
    
   void alloc ();           // allocate the image dimension  
   void cleanDirect ();     // clean fltarr used in recons
                            // when call from a loop as in Linelet3d....
   void cleanInverse ();                            
                               
   void transform(fltarray &Cube, fltarray &PartRad3DImage);
   // Partial Radon transform of a cube
   //              RadonImage  should be first allocated before calling
   //              if it is not alloc is called with default parameter
   //              i.e. number of plan = 6*cube.nx*cube.nx
   //                   size of plan   = cube.nx, cube.nx
   
   void recons(fltarray &PartRad3DImage, fltarray &CubeRec);
   // Radon inverse transform of an image
   // if the allocation has not been done, it calls the alloc routine
   // with the following parameter:
   //      Nx_Cube = RadonImage.nx
   //      Ny_Cube = RadonImage.nx
   //      Nz_Cube = RadonImage.nx
                	
};

#endif

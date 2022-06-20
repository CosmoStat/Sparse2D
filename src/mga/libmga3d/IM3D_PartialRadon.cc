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
**    Date:  15/10/2001 
**    
**    File:  IM_Radon3D.cc
**
**    Modification history:
**
*******************************************************************************
**
**    DESCRIPTION  Radon3D transform and reconstruction
**    -----------  
**                 
******************************************************************************/

#include "IM3D_PartialRadon.h"


#define OUT_CUBE 1000000
#define EPS_IMAGPART 5.e-05
  
static Bool PTakeSym=False;


//------------------------------------------------------------------------------
// 
//

/******************************************************************************/
// control size param
/******************************************************************************/
void PartialRadon3D::control_cube_size (fltarray& Cube) {
} 
        
void PartialRadon3D::control_pradon_size (fltarray& PRadImag) {
}

/******************************************************************************/
// write complex in two array
/******************************************************************************/ 
void PartialRadon3D::WriteCmpArr (char* FileName,  cfarray &CmpArray) {

   char Name[100];
   int Nx = CmpArray.nx(), Ny = CmpArray.ny(), Nz = CmpArray.nz();
   fltarray PRadReal(Nx, Ny, Nz);
   fltarray PRadImag(Nx, Ny, Nz);
   for (int i=0;i<Nx;i++)
   for (int j=0;j<Ny;j++) 
   for (int n=0;n<Nz;n++)
   {
      PRadReal(i,j,n) = CmpArray(i,j,n).real();
      PRadImag(i,j,n) = CmpArray(i,j,n).imag();  
   }
   sprintf (Name, "%sReal.fits", FileName);
   fits_write_fltarr(Name, PRadReal);
   sprintf (Name, "%sImag.fits", FileName);       
   fits_write_fltarr(Name, PRadImag);  
}  

/******************************************************************************/
// write complex in two array
/******************************************************************************/ 
void PartialRadon3D::InfoCmp (cfarray& CmpArray, string Name) {
  
   int Nx = CmpArray.nx(), Ny = CmpArray.ny(), Nz = CmpArray.nz();
   fltarray PRadReal(Nx, Ny, Nz);
   fltarray PRadImag(Nx, Ny, Nz);
   for (int i=0;i<Nx;i++)
   for (int j=0;j<Ny;j++) 
   for (int n=0;n<Nz;n++)
   {
      PRadReal(i,j,n) = CmpArray(i,j,n).real();
      PRadImag(i,j,n) = CmpArray(i,j,n).imag();  
   }
   
   if (Name=="") cout << "  Name:" << Name;
   else cout << "  " << Name << ", Name:" << Name;
   cout << ", mean = (" << PRadReal.mean() << "," 
                        << PRadImag.mean() << ") ,";
   cout << ", sigma = (" << PRadReal.sigma() << "," 
                         << PRadImag.sigma() << ") ,";
   cout << ", min = (" << PRadReal.min() << ","
                       << PRadImag.min() << ") ,";
   
   cout << ", max = (" << PRadReal.max() << ","
                       << PRadImag.max() << ") " << endl;
}

/******************************************************************************/
// control energy, symetry for 
/******************************************************************************/ 
void PartialRadon3D::EnergyCtrl (cfarray &CmpArray) {
   
   int Nx = CmpArray.nx(), Ny = CmpArray.ny(), Nz = CmpArray.nz();   
   fltarray ImagPart (Ny,Nz);
   for (int n=0;n<Nz;n++) {
      for (int i=0;i<Ny;i++) 
      for (int j=0;j<Nx;j++)
         ImagPart(i,j) =  (CmpArray(i,j,n)).imag();
      float Energy = ImagPart.energy();
      if (Energy > 1e-04)
         cout << "Cmp Energy plane(" << n << ") = " << Energy << endl;
   }
}      

void PartialRadon3D::control_even_sym (cfarray& CmpArray) {

   for (int n=0;n<ProjectNbr;n++) {
   
      //dont test row 0, line SizeCube-1
      for (int i=1; i<SizeCube; i++)
      for (int j=i; j<SizeCube; j++) {
      
         int iSym = (SizeCube-1)-(i-1);
         int jSym = (SizeCube-1)-(j-1);
         if (   fabs((CmpArray(i,j,n)).real()-(CmpArray(iSym,jSym,n)).real()) > 1e-03
             || fabs((CmpArray(i,j,n)).imag()+(CmpArray(iSym,jSym,n)).imag()) > 1e-03) {
         
            cout << "Plane number " << n << ", pos(" 
                 << i << "," << j << ")/(" << iSym 
                 << "," << jSym << ") is not sym...."  
                 << (CmpArray(i,j,n)).real() << "/" <<
                 (CmpArray(iSym,jSym,n)).real() << " , "
                 << (CmpArray(i,j,n)).imag() << "/" << 
                 -(CmpArray(iSym,jSym,n)).imag() << endl;
         
         }
      }
   }
}

 
void PartialRadon3D::control_odd_sym (cfarray& CmpArray) {
   
   for (int n=0;n<ProjectNbr;n++) {
   
      for (int i=0; i<SizeCube; i++)
      for (int j=i; j<SizeCube; j++) {
      
         int iSym = (SizeCube-1)-(i);
         int jSym = (SizeCube-1)-(j);
         if (   fabs((CmpArray(i,j,n)).real()-(CmpArray(iSym,jSym,n)).real()) > 1e-03
             || fabs((CmpArray(i,j,n)).imag()+(CmpArray(iSym,jSym,n)).imag()) > 1e-03) {
         
            cout << fabs((CmpArray(i,j,n)).real()-(CmpArray(iSym,jSym,n)).real())
                 << "/"
                 << fabs((CmpArray(i,j,n)).imag()+(CmpArray(iSym,jSym,n)).imag())
                 << endl;
            cout << "Plane number " << n << ", pos(" 
                 << i << "," << j << ")/(" << iSym 
                 << "," << jSym << ") is not sym...."  
                 << (CmpArray(i,j,n)).real() << "/" <<
                 (CmpArray(iSym,jSym,n)).real() << " , "
                 << (CmpArray(i,j,n)).imag() << "/" << 
                 (CmpArray(iSym,jSym,n)).imag() << endl;
         
         }
      }
   }
}



/******************************************************************************/
// reset function
/******************************************************************************/
void PartialRadon3D::reset () {
     
   Verbose=False;
   Control=False;
   WriteCoord=False;
   WriteFile=False;   
   StatInfo=False;
   xPos.free();
   yPos.free();
   zPos.free();
   Alloc = False;
   Plane = False;
   RemoveEqualPlane = False;
   TakeOnlyCompletePlane = False;
   DoNotCorrectSigma = False;
}

/******************************************************************************/
// clean function
/******************************************************************************/
void PartialRadon3D::cleanDirect () {
     
   //DataCount1.init (0.0); 
   if (SizeCube%2==0) {
      DataCount2.init (0);  
      MemPlane.init (-1.0);
   }   
}
void PartialRadon3D::cleanInverse () {

   F3DCube.init (complex_f(0));
   DataCount3.init (0);
} 


/******************************************************************************/
// set param
/******************************************************************************/
void PartialRadon3D::setParam (int BlockSize, int NbrScale) {

   SizeCube = BlockSize;
   if (SizeCube%2==1) ProjectNbr = SizeCube*SizeCube*3;
   else               ProjectNbr = (SizeCube+1)*(SizeCube+1)*3;
   
   make_table();
   if (SizeCube%2==1) make_plane();
   if (StatInfo) stat_alloc();
}

/******************************************************************************/
// alloc. make_table function
/******************************************************************************/
void PartialRadon3D::alloc () {
   Alloc = True;
   
   CFFT3D.CenterZeroFreq = True;

   F3DCube.alloc (SizeCube, SizeCube, SizeCube);
   FPartRad3DImage.alloc (SizeCube, SizeCube, ProjectNbr);
                                                                           
   //DataCount1.alloc (SizeCube,SizeCube,ProjectNbr); 

   if (SizeCube%2==0) {
      DataCount2.alloc (SizeCube,SizeCube,SizeCube);
      MemPlane.alloc (SizeCube,SizeCube,ProjectNbr);
      MemPlane.init (-1.0);
   }   
                                
   DataCount3.alloc (SizeCube, SizeCube, SizeCube);   
}

void PartialRadon3D::make_table() {
   xPos.alloc (ProjectNbr, (char*)"x");
   yPos.alloc (ProjectNbr, (char*)"y");
   zPos.alloc (ProjectNbr, (char*)"z");
      
   int MaxSize;
   int Decal;
   int Cmpt=0;
   if (SizeCube%2==1) {
      MaxSize = SizeCube-1;
      Decal   = SizeCube*SizeCube;
   } else {
      MaxSize = SizeCube;
      Decal   = (SizeCube+1)*(SizeCube+1);   
   }
      
   for (int i=0;i<=MaxSize;i++)
   for (int j=0;j<=MaxSize;j++) {
      
      // first face
      xPos(Cmpt)=i;
      yPos(Cmpt)=j;
      zPos(Cmpt)=MaxSize;
         
      // second face
      xPos(Cmpt+Decal)=i;
      yPos(Cmpt+Decal)=MaxSize;
      zPos(Cmpt+Decal)=j;   
         
      // third face
      xPos(Cmpt+2*Decal)=MaxSize;
      yPos(Cmpt+2*Decal)=i;
      zPos(Cmpt+2*Decal)=j;  
         
      Cmpt++;
   }   
}

void PartialRadon3D::make_plane() {
   Plane = True;
   MemPlane.alloc (SizeCube,SizeCube,ProjectNbr);
   MemPlane.init(-1);
   
   cfarray LocalF3DCube, LocalFPartRad3DImage; // not used....                         
   comp_Plane (LocalF3DCube, LocalFPartRad3DImage, False);
   if (WriteFile) fits_write_fltarr((char*)"MemPlane.fits", MemPlane);
   
   if (RemoveEqualPlane) remove_equal_Plane(); 
   //if (StatInfo) memorize_stat();       
}
   
void PartialRadon3D::remove_equal_Plane () {
   int Compt=0;
   Bool FlagNotEqual=False;
   for (int n=0;n<ProjectNbr;n++) {
   
      for (int m=n+1;m<ProjectNbr;m++) {
      
         FlagNotEqual=False;
         for (int i=0;i<SizeCube;i++) {
            for (int j=0;j<SizeCube;j++) {
      
               if (MemPlane(i,j,n) != MemPlane(i,j,m)) {
                 FlagNotEqual=True; 
                 break;
               }
            }
            if (FlagNotEqual) break;
         }
         if (!FlagNotEqual) {
            //cout << n << "," << m << "," << (int)Flag << endl;
            break;
         }
      }
      if (FlagNotEqual) {
         xPos(Compt) = xPos(n);
         yPos(Compt) = yPos(n);
         zPos(Compt) = zPos(n);
         Compt++; //one more plane
      }
   }
   ProjectNbr = Compt;                            
}


//------------------------------------------------------------------------------
// stat function
//

/******************************************************************************/
// stat_alloc
/******************************************************************************/
void PartialRadon3D::stat_alloc() {
   CoefSigmaPlane.alloc (ProjectNbr);
   MemNotUsed.alloc (ProjectNbr);
   
   // if SizeCube even, MemPlane is initialized in direct_make_plane                         
   if (SizeCube%2==1) stat_init();
}


/******************************************************************************/
// stat_init
/******************************************************************************/
void PartialRadon3D::stat_init() {
   CoefSigmaPlane.init(1.);   
   
   for (int n=0;n<ProjectNbr;n++) {
  
      int Compt=0;
      for (int i=0;i<SizeCube;i++) 
      for (int j=0;j<SizeCube;j++) {
         if (MemPlane(i,j,n) == -1) Compt++; 
      }
      MemNotUsed(n) = Compt;
   } 
}


/******************************************************************************/
// stat_comp_true_sigma
/******************************************************************************/
float PartialRadon3D::stat_comp_true_sigma (int NbPtsUsed, int NbZero, 
                       float Sigma, float Mean) {

   float Variance =  (((float)NbZero+NbPtsUsed)/NbPtsUsed * Sigma*Sigma)
         - ((float)NbZero*(NbZero+NbPtsUsed)/NbPtsUsed/NbPtsUsed * Mean*Mean);
   return (sqrt(Variance));
}


/******************************************************************************/
// stat_comp_coef
/******************************************************************************/
void PartialRadon3D::stat_comp_coef (cfarray& FPartRad3DImage) {
/*
   sprintf (PRadMsg, "Alloc cmp[%d], int[%d]...", ProjectNbr, ProjectNbr);
   TRACEUR << new Trace (Trace::MEMORY, "PartialRadon3D", "stat_comp_coef",
                         PRadMsg, NULL, Trace::DEBUG); 
                         
   cfarray MeanPlane (ProjectNbr);
   fltarray SigmaPlane (ProjectNbr);
   
   TRACEUR << new Trace (Trace::MEMORY, "PartialRadon3D", "stat_alloc",
                         "End Alloc", NULL, Trace::DEBUG); 
                            
   for (int n=0;n<ProjectNbr;n++) {  
   
      for (int i=0;i<SizeCube;i++) 
      for (int j=0;j<SizeCube;j++) { 
         if (MemPlane(i,j,n) >= 0)    
            MeanPlane(n)  += FPartRad3DImage(i,j,n);
      }
   }
*/  
   
   for (int n=0;n<ProjectNbr;n++) {
      
      int NbPointsUsed = SizeCube*SizeCube-MemNotUsed(n);
      
/*
      MeanPlane(n)/=(float)(NbPointsUsed);
      float Som=0.;
      
      for (int i=0;i<SizeCube;i++) 
      for (int j=0;j<SizeCube;j++) { 
         if (MemPlane(i,j,n) >= 0)    
            Som += norm (FPartRad3DImage(i,j,n)-MeanPlane(n));
      }
      SigmaPlane(n) =   sqrt((float)(Som / (float)(NbPointsUsed)));
*/
                    
      complex_f m=complex_f(0.,0.);           
      for (int i=0;i<SizeCube;i++) 
      for (int j=0;j<SizeCube;j++)
         m += FPartRad3DImage(i,j,n);
      m /= (SizeCube*SizeCube);
      float s=0;   
      for (int i=0;i<SizeCube;i++) 
      for (int j=0;j<SizeCube;j++)   
         s += norm (FPartRad3DImage(i,j,n)-m);
      s = sqrt(s/(SizeCube*SizeCube));              
                      
      float TrueSigma = stat_comp_true_sigma (NbPointsUsed, MemNotUsed(n),
                                              s, m.real());
      CoefSigmaPlane(n) = TrueSigma/s;

/*
      float Norm = SizeCube*sqrt((float)SizeCube);
      cout << "[pn:" << n << "," << NbPointsUsed << "pts]" 
           << ", m1:" << MeanPlane(n)/Norm
           << ", m2:" << m/Norm << endl;
      cout << "[pn:" << n 
           << ", s1:" << SigmaPlane(n)/Norm  
           << ", s2:" << s/Norm  
           << ", s3:" << TrueSigma/Norm 
           << ", coef:" << CoefSigmaPlane(n) << endl;
*/
   }   
}

/******************************************************************************/
// stat_info
/******************************************************************************/
void PartialRadon3D::stat_info (fltarray& PartRad3DImage) {
   int NbPlaneUsed = 0;
   if (TakeOnlyCompletePlane) {
      for (int k=0;k<ProjectNbr;k++) 
         if (MemNotUsed(k) == 0) NbPlaneUsed++;
      cout << "NbPlaneComplete:" << NbPlaneUsed << endl;
   } else { 
      NbPlaneUsed = ProjectNbr;
   }
                                      
   fltarray Plane(SizeCube,SizeCube);
   fltarray Mean(NbPlaneUsed);
   fltarray Sigma(NbPlaneUsed);   
   fltarray EstimSigma(NbPlaneUsed);
   
   int ComptPlaneUsed=0;
   for (int n=0;n<ProjectNbr;n++) {
           
      for (int i=0;i<SizeCube;i++) 
      for (int j=0;j<SizeCube;j++) Plane(i,j) =  PartRad3DImage(i,j,n);
      
      int NbPointsUsed = SizeCube*SizeCube-MemNotUsed(n);
      
      if (TakeOnlyCompletePlane) {
         if (NbPointsUsed == SizeCube*SizeCube) {
            Mean(ComptPlaneUsed) =  Plane.mean();
            Sigma(ComptPlaneUsed) = Plane.sigma();         
            EstimSigma(ComptPlaneUsed) = Plane.sigma();
            ComptPlaneUsed++;
         }
      } else {
         Mean(n) =  Plane.mean();
         Sigma(n) = Plane.sigma();
         if (NbPointsUsed != SizeCube*SizeCube && !DoNotCorrectSigma) {
            EstimSigma(n) = CoefSigmaPlane(n)*Sigma(n);
         } else {
            EstimSigma(n) = Sigma(n);
         }
      }
        
      //sprintf (PRadMsg, "[plane:%d,pts:%d]  Sig:%f  EstSig:%f", 
      //         n, NbPointsUsed, Sigma(ComptPlaneUsed-1), 
        //       EstimSigma(ComptPlaneUsed-1));
      //TRACEUR << new Trace (Trace::METIER, "PartialRadon3D", "stat_info",
      //                      PRadMsg, NULL, Trace::DEBUG); 
   }
   fltarray StatInf(2);
   StatInf(0) = EstimSigma.mean(); StatInf(1) = EstimSigma.sigma();
   //fits_write_fltarr ("StatPRadTransf.fits", StatInf);
}

/******************************************************************************/
// iround mrethod
/******************************************************************************/
int PartialRadon3D::double2int (double point) {
   
   //int result;
   //if (point >= (double) 0.0) result = (int) (point);
   //else result = (int) (point);
   //return result;
   int result;
   if (point >= (double) 0.0) result = (int) ( point + (double) 0.5001);
   else result = (int) (point - 0.5001);
   // printf ("iiround: %F ==> %F, %d \n", point, point+0.5, result);
   return result;
}



//------------------------------------------------------------------------------
// 
//

/******************************************************************************/
// transform function
/******************************************************************************/
void PartialRadon3D::transform (fltarray& Cube, fltarray& PartRad3DImage) {
   if (!Alloc) {
      cout << "Classe must be initialized!" << endl; exit(-1);
   }
   
   if (Control) {
      control_cube_size  (Cube);
      control_pradon_size (PartRad3DImage);  
   }
   partial_slice (Cube, PartRad3DImage);
   
   if (StatInfo) stat_info (PartRad3DImage);
}

/******************************************************************************/
// recons function
/******************************************************************************/
void PartialRadon3D::recons(fltarray& PartRad3DImage, fltarray& CubeRec) {
   if (Control) {
      control_cube_size  (CubeRec);
     control_pradon_size (PartRad3DImage); 
   }
   

   int RecNbPlane=PartRad3DImage.nz();
   int RecSizeCube=PartRad3DImage.nx();
   
   // recons for even size if alloc true only
   
   if (!Alloc) {
cerr<<"realloc"<<endl;
      setParam (RecSizeCube);
      alloc();
   }
   
   inv_partial_slice (PartRad3DImage, CubeRec);
}



//-----------f (Control) CmpControl(numPlane,i,j)=CmpCube(-------------------------------------------------------------------
//
//

/******************************************************************************/
// partial_slice function
/******************************************************************************/
void PartialRadon3D::partial_slice (fltarray& Cube, 
                                    fltarray& PartRad3DImage) {
   cleanDirect();
   
   comp_FFT3D (Cube, F3DCube);
   if (SizeCube%2==0 || !Plane) comp_Plane (F3DCube, FPartRad3DImage);
   else comp_odd_Plane (F3DCube, FPartRad3DImage);
   comp_InvFFT2D (FPartRad3DImage, PartRad3DImage);
}

/******************************************************************************/
// inv_partial_slice function
/******************************************************************************/
void PartialRadon3D::inv_partial_slice (fltarray& PartRad3DImage, 
                                        fltarray& Cube) {
   cleanInverse();
                            
   comp_FFT2D (PartRad3DImage, FPartRad3DImage);
   if (SizeCube%2==0 || !Plane) comp_Plane(F3DCube, FPartRad3DImage, True);
   else comp_odd_Plane (F3DCube, FPartRad3DImage, True);
   comp_InvFFT3D(F3DCube, Cube);
}


//------------------------------------------------------------------------------
//
//

/******************************************************************************/
// comp_FFT3D function
/******************************************************************************/
void PartialRadon3D::comp_FFT3D (fltarray& Cube, cfarray& F3DCube) {

   if (WriteFile)
      fits_write_fltarr((char*)"InFFT3D", Cube);
//FFTN_3D CFFT3D_local;
//CFFT3D_local.CenterZeroFreq = True;
//CFFT3D_local.fftn3d(Cube, F3DCube);
   CFFT3D.fftn3d(Cube, F3DCube); 
   if (WriteFile) 
      WriteCmpArr ((char*)"OutFFT3D", F3DCube);
}


/******************************************************************************/
// comp_InvFFT3D function
/******************************************************************************/
void PartialRadon3D::comp_InvFFT3D (cfarray& F3DCube, fltarray& Cube) {

   if (WriteFile)
      WriteCmpArr ((char*)"InInvFFT3D", F3DCube);
//FFTN_3D CFFT3D_local;
//CFFT3D_local.CenterZeroFreq = True;
//CFFT3D_local.fftn3d(F3DCube, True);
   CFFT3D.fftn3d(F3DCube, True);
   if (WriteFile) 
      WriteCmpArr ((char*)"OutInvFFT3D", F3DCube);
   for (int i=0;i<SizeCube;i++)
   for (int j=0;j<SizeCube;j++)
   for (int k=0;k<SizeCube;k++)
      Cube(i,j,k) = (float)(F3DCube(i,j,k)).real(); 
}



//------------------------------------------------------------------------------
//
//

/******************************************************************************/
// comp_Plane function
/******************************************************************************/
void PartialRadon3D::comp_odd_Plane (cfarray& F3DCube, cfarray& FPartRad3DImage,
                                     Bool Inverse) { 
//   cerr<<"PRAD::Comp odd plane"<<endl;
   if (!Inverse) {   
      if (WriteFile) WriteCmpArr ((char*)"InPlane", F3DCube);
   } else {
      if (WriteFile) WriteCmpArr ((char*)"InInvPlane", FPartRad3DImage);
   }                                      

   for (int numLine=0; numLine<ProjectNbr; numLine++) 
   for (int l1=0; l1<SizeCube; l1++)
   for (int l2=0; l2<SizeCube; l2++) {   
      
      int xPix,yPix,zPix;
      if (!Inverse && MemPlane(l1,l2,numLine)>=0) {
         
         comp_xyz (l1,l2,numLine,xPix,yPix,zPix);
         FPartRad3DImage(l1,l2,numLine) = F3DCube(xPix,yPix,zPix);
         
      } 
      else 
      if (Inverse && MemPlane(l1,l2,numLine)>=0) {
      
         comp_xyz (l1,l2,numLine,xPix,yPix,zPix);
         F3DCube(xPix,yPix,zPix) += FPartRad3DImage(l1,l2,numLine); 
         DataCount3(xPix,yPix,zPix)++;  
      }
   }       
   
   if (!Inverse) {
      if (Control) control_odd_sym (FPartRad3DImage);
      if (WriteFile) WriteCmpArr ((char*)"OutPlane", FPartRad3DImage);
      if (StatInfo) stat_comp_coef (FPartRad3DImage);
   }
   
   if (Inverse) {
 
      for (int i=0;i<SizeCube;i++)
      for (int j=0;j<SizeCube;j++) 
      for (int k=0;k<SizeCube;k++) {
  
         if (DataCount3(i,j,k) != 0) F3DCube(i,j,k)/=(double)DataCount3(i,j,k);
         else F3DCube(i,j,k)=complex<float>(0.,0.);
         
      }        
      if (WriteFile) {
         WriteCmpArr ((char*)"OutInvPlane", F3DCube);
         fits_write_intarr((char*)"DCount3.fits", DataCount3);   
      }
   } 
}                

void PartialRadon3D::comp_Plane (cfarray& F3DCube, cfarray& FPartRad3DImage,
                                 Bool Inverse) {  
   tp_setfunc set_f; 
   if (SizeCube%2==0) {
      if (!Inverse) {   
         if (WriteFile) WriteCmpArr ((char*)"InPlane", F3DCube);
         //FPartRad3DImage.init(complex_f(0));
         set_f = &PartialRadon3D::direct_make_plane;
      } else {
         if (WriteFile) WriteCmpArr ((char*)"InInvPlane", FPartRad3DImage);
         //F3DCube.init(complex_f(0));
         set_f = &PartialRadon3D::inverse_make_plane;
      } 
   } else {
      set_f = &PartialRadon3D::memorize_make_plane;
   }
     
   int MidPixVal;
   if (SizeCube%2==0) MidPixVal = SizeCube/2;        // even signal size
   else               MidPixVal = (SizeCube-1)/2;    // odd signal size
   
   int MaxBorder;
   if (SizeCube%2==0) MaxBorder = SizeCube;         // even signal size
   else               MaxBorder = SizeCube-1;       // odd signal size
    
   
   for (int numLine=0; numLine<ProjectNbr; numLine++) {
        
      double a = xPos(numLine)-SizeCube/2;
      double b = yPos(numLine)-SizeCube/2;
      double c = zPos(numLine)-SizeCube/2;
      
      if (WriteCoord) 
         cout << "["<<xPos(numLine)<<","<<yPos(numLine)<<","<<zPos(numLine)
              << ", numPlane:" << numLine 
              << ",   {" << a << "," << b << "," << c << "}" << endl;
 
      // c=0 => w=8
      if (zPos(numLine) == MidPixVal) { // (c==0) vertival plane
      
         //b=0 => v=8 (plane c=0,b=0,a)
         if (yPos(numLine) == MidPixVal) {

            tp_mfunc f_x = &PartialRadon3D::fx_uXv8w8;
            tp_mfunc f_y = &PartialRadon3D::fy_uXv8w8;
            tp_mfunc f_z = &PartialRadon3D::fz_uXv8w8;
            
            if (SizeCube%2==0)
               even_make_plane (F3DCube, FPartRad3DImage, f_x, f_y, f_z,
                                a, b, c, numLine, set_f, Inverse);
            else
               odd_make_plane (F3DCube, FPartRad3DImage, f_x, f_y, f_z,
                                a, b, c, numLine, set_f, Inverse);
            
         } else {
            // cout << numLine << endl;
            
            //theta in [PI/2,3PI/2] 
            if (yPos(numLine) == 0 || yPos(numLine) == MaxBorder) {
            
               tp_mfunc f_x = &PartialRadon3D::fx1_uXvXw8;
               tp_mfunc f_y = &PartialRadon3D::fy1_uXvXw8;
               tp_mfunc f_z = &PartialRadon3D::fz1_uXvXw8;
         
               if (SizeCube%2==0)
                  even_make_plane (F3DCube, FPartRad3DImage, f_x, f_y, f_z,
                                   a, b, c, numLine, set_f, Inverse);
               else
                  odd_make_plane (F3DCube, FPartRad3DImage, f_x, f_y, f_z,
                                   a, b, c, numLine, set_f, Inverse);             
            
            //theta in [0,PI/2[ U ]3PI/2,PI]
            } else {
                           
                tp_mfunc f_x = &PartialRadon3D::fx2_uXvXw8;
               tp_mfunc f_y = &PartialRadon3D::fy2_uXvXw8;
               tp_mfunc f_z = &PartialRadon3D::fz2_uXvXw8;
         
               if (SizeCube%2==0)
                  even_make_plane (F3DCube, FPartRad3DImage, f_x, f_y, f_z,
                                   a, b, c, numLine, set_f, Inverse);
               else
                  odd_make_plane (F3DCube, FPartRad3DImage, f_x, f_y, f_z,
                                   a, b, c, numLine, set_f, Inverse);                        
            }
         }
         
      } else {

         // a != 0
         if (xPos(numLine) == 0 || xPos(numLine) == MaxBorder) {
                  
            tp_mfunc f_x = &PartialRadon3D::fx1_uXvXwX;
            tp_mfunc f_y = &PartialRadon3D::fy1_uXvXwX;
            tp_mfunc f_z = &PartialRadon3D::fz1_uXvXwX;
            
            if (SizeCube%2==0)
               even_make_plane (F3DCube, FPartRad3DImage, f_x, f_y, f_z,
                                a, b, c, numLine, set_f, Inverse);
            else
               odd_make_plane (F3DCube, FPartRad3DImage, f_x, f_y, f_z,
                                a, b, c, numLine, set_f, Inverse); 
                        
         } else 
         // b != 0
         if (yPos(numLine) == 0 || yPos(numLine) == MaxBorder) {
         
            tp_mfunc f_x = &PartialRadon3D::fx2_uXvXwX;
            tp_mfunc f_y = &PartialRadon3D::fy2_uXvXwX;
            tp_mfunc f_z = &PartialRadon3D::fz2_uXvXwX;
            
            if (SizeCube%2==0)
               even_make_plane (F3DCube, FPartRad3DImage, f_x, f_y, f_z,
                                a, b, c, numLine, set_f, Inverse);
            else
               odd_make_plane (F3DCube, FPartRad3DImage, f_x, f_y, f_z,
                                a, b, c, numLine, set_f, Inverse); 
                  
         } else
         // others
         if (zPos(numLine) == 0 || zPos(numLine) == MaxBorder) {
         
            tp_mfunc f_x = &PartialRadon3D::fx3_uXvXwX;
            tp_mfunc f_y = &PartialRadon3D::fy3_uXvXwX;
            tp_mfunc f_z = &PartialRadon3D::fz3_uXvXwX;
            
            if (SizeCube%2==0)
               even_make_plane (F3DCube, FPartRad3DImage, f_x, f_y, f_z,
                                a, b, c, numLine, set_f, Inverse);
            else
               odd_make_plane (F3DCube, FPartRad3DImage, f_x, f_y, f_z,
                                a, b, c, numLine, set_f, Inverse);
         }         
      }
   }    
   
   
   if (!Inverse) { 
   
      if (Control && SizeCube%2==0)
         control_even_sym (FPartRad3DImage);
      //if (Control && SizeCube%2==1)
      //   control_odd_sym (FPartRad3DImage);
         
      if (WriteFile && SizeCube%2==0) {
         WriteCmpArr ((char*)"OutPlane", FPartRad3DImage);   
         //fits_write_intarr ("DCount1.fits", DataCount1); 
         fits_write_intarr ((char*)"DCount2.fits", DataCount2);
      }
      
      if (StatInfo && SizeCube%2==0) {
         stat_init ();
         stat_comp_coef (FPartRad3DImage);
      }
   } 

   if (Inverse) {
     
      for (int i=0;i<SizeCube;i++)
      for (int j=0;j<SizeCube;j++) 
      for (int k=0;k<SizeCube;k++) {
  
         if (DataCount3(i,j,k) != 0) F3DCube(i,j,k)/=(double)DataCount3(i,j,k);
         else F3DCube(i,j,k)=complex<float>(0.,0.);
      }
      
      if (WriteFile && SizeCube%2==0) {
         WriteCmpArr ((char*)"OutInvPlane", F3DCube);
         fits_write_intarr((char*)"DCount3.fits", DataCount3);      
      }
   }  
}


//------------------------------------------------------------------------------
//
//

/******************************************************************************/
// comp_InvFFT2D function
/******************************************************************************/
void PartialRadon3D::comp_InvFFT2D (cfarray& FPartRad3DImage , 
                                    fltarray& PRadImag) {
   cfarray CmpControl;
   if (WriteFile || Control) {
      WriteCmpArr ((char*)"InInvFFT2D", FPartRad3DImage);
      CmpControl.alloc (SizeCube,SizeCube,ProjectNbr);
   }
                                               
   FFTN_2D FFT2D;
   FFT2D.CenterZeroFreq = True;
   Icomplex_f CmpCube(SizeCube,SizeCube);  

   for (int numPlane=0; numPlane<ProjectNbr ; numPlane++) {
   
      for (int i=0;i<SizeCube;i++)
      for (int j=0;j<SizeCube;j++) {
         CmpCube(j,i) = FPartRad3DImage(i,j,numPlane);  
      }
      FFT2D.fftn2d (CmpCube, True);
      for (int i=0;i<SizeCube;i++)
      for (int j=0;j<SizeCube;j++) { 
         if (WriteFile || Control) CmpControl(i,j,numPlane)=CmpCube(j,i); 
         PRadImag(i,j,numPlane) = (float) (CmpCube(j,i)).real();
      }    
   }
   
   if (WriteFile) WriteCmpArr ((char*)"OutInvFFT2D", CmpControl);
   if (WriteFile) fits_write_fltarr ((char*)"OutInvFFT2D", PRadImag);
   if (Control) EnergyCtrl (CmpControl);
}

/******************************************************************************/
// comp_FFT2D function
/******************************************************************************/
void PartialRadon3D::comp_FFT2D (fltarray& PartRad3DImage, 
                                 cfarray& FPartRad3DImage) {
   if (Control) fits_write_fltarr ((char*)"InFFT2D", PartRad3DImage);
                                                  
   FFTN_2D FFT2D;
   FFT2D.CenterZeroFreq = True;
   Icomplex_f CmpCube(SizeCube,SizeCube);  
   
   for (int numPlane=0; numPlane<ProjectNbr ; numPlane++) {
      
      for (int i=0;i<SizeCube;i++)
      for (int j=0;j<SizeCube;j++) {
         CmpCube(j,i) = PartRad3DImage(i,j,numPlane);
      }
      FFT2D.fftn2d (CmpCube, False);
      for (int i=0;i<SizeCube;i++)
      for (int j=0;j<SizeCube;j++) {  
         FPartRad3DImage(i,j,numPlane) = CmpCube(j,i);
      }    
   }
   if (Control) WriteCmpArr ((char*)"OutFFT2D", FPartRad3DImage);
}   

 






//------------------------------------------------------------------------------
//
//

/******************************************************************************/
// internal function
/******************************************************************************/


void PartialRadon3D::direct_make_plane  (cfarray& F1, 
                                         cfarray& F2,
                                         int n, int l1, int l2, 
                                         int x, int y, int z,
                                         Bool TakeSym) {                                 

    F2(l1,l2,n) = F1(x,y,z);
    //DataCount1(l1,l2,n)++; 
    if (SizeCube%2==0) {
       DataCount2(x,y,z)++; 
       MemPlane (l1,l2,n)=x+SizeCube*y+z*SizeCube*SizeCube;
    }
    if (TakeSym) {
       F2(l1,l2,n)=complex_f (F1(x,y,z).real(),-F1(x,y,z).imag());
   }
}
                            
void PartialRadon3D::inverse_make_plane (cfarray& F1, 
                                         cfarray& F2,
                                         int n, int l1, int l2, 
                                         int x, int y, int z,
                                         Bool TakeSym) {

   //cout <<"["<<n<<","<<l1<<","<<l2<<"]--{"<<x<<","<<y<<","<<z<<"}"<<endl;
   DataCount3(x,y,z)++;         
   if (!TakeSym) {
      F1(x,y,z) += F2(l1,l2,n);
   } else {
      F1(x,y,z) += complex_f(F2(l1,l2,n).real(),-F2(l1,l2,n).imag());
   }                                                           
}

void PartialRadon3D::memorize_make_plane  (cfarray& F1, 
                                           cfarray& F2,
                                           int n, int l1, int l2, 
                                           int x, int y, int z,
                                           Bool TakeSym) {
                                           
    //cout <<"["<<n<<","<<l1<<","<<l2<<"]--{"<<x<<","<<y<<","<<z<<"}"<<endl;                                                                          
    MemPlane(l1,l2,n)=x+SizeCube*y+z*SizeCube*SizeCube;
}

Bool PartialRadon3D::even_border_manager (int x, int y, int z, 
                                     int& xPix, int& yPix, int& zPix) {

   if (x==SizeCube/2) {
      xPix=0; yPix=-y+SizeCube/2; zPix=-z+SizeCube/2;
   } else if (y==SizeCube/2) {    
      xPix=-x+SizeCube/2; yPix=0; zPix=-z+SizeCube/2;   
   } else if (z==SizeCube/2) { 
      xPix=-x+SizeCube/2; yPix=-y+SizeCube/2; zPix=0;   
   } 
   
   if (xPix != SizeCube && yPix != SizeCube && zPix != SizeCube)
      return True;
   else
      return False;
}


void PartialRadon3D::comp_xyz (int l1, int l2, int n,
                               int& xPix, int& yPix, int& zPix) {

  // called only in even case...                               
  int Inter = int(MemPlane(l1,l2,n))%(SizeCube*SizeCube);
  xPix = Inter%SizeCube;
  yPix = (Inter-xPix)/SizeCube;
  zPix = (int(MemPlane(l1,l2,n))-Inter)/SizeCube/SizeCube;                                                                                
}


void PartialRadon3D::even_make_plane (cfarray& InPlane, cfarray& OutPlane,
                                 tp_mfunc f1, tp_mfunc f2, tp_mfunc f3, 
                                 double a, double b, double c,
                                 int NumLine, tp_setfunc set_f, 
                                 Bool Inverse)  {

   // local var
   intarray PlaneMask (SizeCube, SizeCube);
   intarray CubeMask  (SizeCube, SizeCube, SizeCube);
           
   // symetrical part       
   for (int l1=-SizeCube/2+1;l1<SizeCube/2;l1++)
   for (int l2=l1+1;l2<SizeCube/2;l2++) {
   
      int x = (this->*f1)(l1,l2,a,b,c);
      int y = (this->*f2)(l1,l2,a,b,c);
      int z = (this->*f3)(l1,l2,a,b,c);               
      int xPix = x + SizeCube/2;
      int yPix = y + SizeCube/2;
      int zPix = z + SizeCube/2;

      if (x == OUT_CUBE || y == OUT_CUBE || z == OUT_CUBE) {
         PlaneMask(l1+SizeCube/2,l2+SizeCube/2)=0;
      } else 
      if (x == SizeCube/2 || y == SizeCube/2 || z == SizeCube/2) {
         if (even_border_manager (x,y,z,xPix,yPix,zPix)) {
            PlaneMask(l1+SizeCube/2,l2+SizeCube/2)=1;
            CubeMask(xPix,yPix,zPix)++;
            (this->*set_f)(InPlane, OutPlane, NumLine, l1+SizeCube/2,
                           l2+SizeCube/2,xPix,yPix,zPix,True);
         } else {
            PlaneMask(l1+SizeCube/2,l2+SizeCube/2)=0;
         }
      } else {
         PlaneMask(l1+SizeCube/2,l2+SizeCube/2)=1;
         CubeMask(xPix,yPix,zPix)++;
         (this->*set_f)(InPlane, OutPlane, NumLine, l1+SizeCube/2,
                        l2+SizeCube/2,xPix,yPix,zPix, PTakeSym);
         //OutPlane(NumLine,l1+SizeCube/2,l2+SizeCube/2) = InPlane(xPix,yPix,zPix);
      }

         xPix = -x + SizeCube/2;
         yPix = -y + SizeCube/2;
         zPix = -z + SizeCube/2;    

         if (x == OUT_CUBE || y == OUT_CUBE || z == OUT_CUBE) {
            PlaneMask(-l1+SizeCube/2,-l2+SizeCube/2)=0;
         } else if (-x == SizeCube/2 || -y == SizeCube/2 || -z == SizeCube/2) {
            if (even_border_manager (-x,-y,-z,xPix,yPix,zPix)) {
               PlaneMask(-l1+SizeCube/2,-l2+SizeCube/2)=1;
               CubeMask(xPix,yPix,zPix)++;
               (this->*set_f)(InPlane, OutPlane, NumLine, -l1+SizeCube/2,
                              -l2+SizeCube/2,xPix,yPix,zPix,True);                   
            } else {
              PlaneMask(-l1+SizeCube/2,-l2+SizeCube/2)=0;
            }
         } else {
            PlaneMask(-l1+SizeCube/2,-l2+SizeCube/2)=1;
            CubeMask(xPix,yPix,zPix)++;
            (this->*set_f)(InPlane, OutPlane, NumLine, -l1+SizeCube/2,
                           -l2+SizeCube/2,xPix,yPix,zPix, PTakeSym);
            //OutPlane(NumLine,-l1+SizeCube/2,-l2+SizeCube/2) = InPlane(xPix,yPix,zPix);
         }
   }
   
   // symetrical diagonal
   for (int l1=-SizeCube/2+1;l1<1;l1++) {
      int l2=l1;
      
      int x = (this->*f1)(l1,l2,a,b,c);
      int y = (this->*f2)(l1,l2,a,b,c);
      int z = (this->*f3)(l1,l2,a,b,c);               
      int xPix = x + SizeCube/2;
      int yPix = y + SizeCube/2;
      int zPix = z + SizeCube/2;
      if (x == OUT_CUBE || y == OUT_CUBE || z == OUT_CUBE) {
         PlaneMask(l1+SizeCube/2,l2+SizeCube/2)=0;
      } else 
      if (x == SizeCube/2 || y == SizeCube/2 || z == SizeCube/2) {
         if (even_border_manager (x,y,z,xPix,yPix,zPix)) {
            PlaneMask(l1+SizeCube/2,l2+SizeCube/2)=1;
            CubeMask(xPix,yPix,zPix)++;
            (this->*set_f)(InPlane, OutPlane, NumLine, l1+SizeCube/2,
                           l2+SizeCube/2,xPix,yPix,zPix,True);
         } else {
            PlaneMask(l1+SizeCube/2,l2+SizeCube/2)=0;
         }
      } else {
         PlaneMask(l1+SizeCube/2,l2+SizeCube/2)=1;
         CubeMask(xPix,yPix,zPix)++;
         (this->*set_f)(InPlane, OutPlane, NumLine, l1+SizeCube/2,
                        l2+SizeCube/2,xPix,yPix,zPix, PTakeSym);
         //OutPlane(NumLine,l1+SizeCube/2,l2+SizeCube/2) = InPlane(xPix,yPix,zPix);
      }

      if ( l1!=0) { // central point do not have a sym !!
         xPix = -x + SizeCube/2;
         yPix = -y + SizeCube/2;
         zPix = -z + SizeCube/2;      
         if (x == OUT_CUBE || y == OUT_CUBE || z == OUT_CUBE) {
            PlaneMask(-l1+SizeCube/2,-l2+SizeCube/2)=0;
         } else if (-x == SizeCube/2 || -y == SizeCube/2 || -z == SizeCube/2) {
            if (even_border_manager (-x,-y,-z,xPix,yPix,zPix)) {
               PlaneMask(-l1+SizeCube/2,-l2+SizeCube/2)=1;
               CubeMask(xPix,yPix,zPix)++;
               (this->*set_f)(InPlane, OutPlane, NumLine, -l1+SizeCube/2,
                              -l2+SizeCube/2,xPix,yPix,zPix,True);                   
            } else {
              PlaneMask(-l1+SizeCube/2,-l2+SizeCube/2)=0;
            }
         } else {
            PlaneMask(-l1+SizeCube/2,-l2+SizeCube/2)=1;
            CubeMask(xPix,yPix,zPix)++;
            (this->*set_f)(InPlane, OutPlane, NumLine, -l1+SizeCube/2,
                           -l2+SizeCube/2,xPix,yPix,zPix,PTakeSym);
            //OutPlane(NumLine,-l1+SizeCube/2,-l2+SizeCube/2) = InPlane(xPix,yPix,zPix);
         }
      }
   }
   
   int Inf,Sup;
   //if (NumLine%2 ==0) {
      Inf=-SizeCube/2+1; Sup=0;
   //} else {
   //   Inf=1;Sup=SizeCube/2;
   //}
   // first row
   //for (int l1=-SizeCube/2;l1<1;l1++) {
   //for (int l1=-SizeCube/2;l1<SizeCube/2;l1++) {
   for (int l1=Inf;l1<Sup;l1++) {
   
      int l2 = -SizeCube/2;
      int x = (this->*f1)(l1,l2,a,b,c);
      int y = (this->*f2)(l1,l2,a,b,c);
      int z = (this->*f3)(l1,l2,a,b,c);
      int xPix = x + SizeCube/2;
      int yPix = y + SizeCube/2;
      int zPix = z + SizeCube/2;  

      int xSym = (this->*f1)(-l1,l2,a,b,c);
      int ySym = (this->*f2)(-l1,l2,a,b,c);
      int zSym = (this->*f3)(-l1,l2,a,b,c);
      int xPixSym = xSym + SizeCube/2;
      int yPixSym = ySym + SizeCube/2;
      int zPixSym = zSym + SizeCube/2;  
  
      if (x == OUT_CUBE || y == OUT_CUBE || z == OUT_CUBE ||
          x==SizeCube/2 || y==SizeCube/2 || z==SizeCube/2 ||
          xSym == OUT_CUBE || ySym == OUT_CUBE || zSym == OUT_CUBE ||
          xSym==SizeCube/2 || ySym==SizeCube/2 || zSym==SizeCube/2) {
            
            PlaneMask(l1+SizeCube/2,l2+SizeCube/2)=0;
            PlaneMask(-l1+SizeCube/2,l2+SizeCube/2)=0;
         
      } else {
                        
            if (!Inverse) {
            
               if (DataCount2(xPix,yPix,zPix) > DataCount2(xPixSym,yPixSym,zPixSym)) {
                  xPix = xPixSym;
                  yPix = yPixSym;
                  zPix = zPixSym;
               }              
            
               PlaneMask(l1+SizeCube/2,l2+SizeCube/2)=1;
               CubeMask(xPix,yPix,zPix)++;
               (this->*set_f)(InPlane, OutPlane, NumLine, l1+SizeCube/2,
                              l2+SizeCube/2,xPix,yPix,zPix,PTakeSym);
                              
               PlaneMask(-l1+SizeCube/2,l2+SizeCube/2)=1;
               CubeMask(xPix,yPix,zPix)++;
               (this->*set_f)(InPlane, OutPlane, NumLine, -l1+SizeCube/2,
                              l2+SizeCube/2,xPix,yPix,zPix,True);
            } else {
            
               if (MemPlane(l1+SizeCube/2,l2+SizeCube/2,NumLine)>=0) {
                  comp_xyz (l1+SizeCube/2,l2+SizeCube/2,NumLine,xPix,yPix,zPix);
                  PlaneMask(l1+SizeCube/2,l2+SizeCube/2)=1;
                  CubeMask(xPix,yPix,zPix)++;
                  (this->*set_f)(InPlane, OutPlane, NumLine, l1+SizeCube/2,
                                 l2+SizeCube/2,xPix,yPix,zPix,PTakeSym);
               }
               
               if (MemPlane(-l1+SizeCube/2,l2+SizeCube/2,NumLine)>=0) {
                  comp_xyz (-l1+SizeCube/2,l2+SizeCube/2,NumLine,xPix,yPix,zPix);              
                  PlaneMask(-l1+SizeCube/2,l2+SizeCube/2)=1;
                  CubeMask(xPix,yPix,zPix)++;
                  (this->*set_f)(InPlane, OutPlane, NumLine, -l1+SizeCube/2,
                                 l2+SizeCube/2,xPix,yPix,zPix,True);
               }            
            }                
         }      
   }
   
                  
   // first line
   //   if (NumLine%2 ==1) {
  // Inf=-SizeCube/2+1; Sup=0;
   //} else {
      Inf=1;Sup=SizeCube/2;
   //}
   
   //for (int l2=-SizeCube/2+1;l2<1;l2++) {
   //for (int l2=-SizeCube/2+1;l2<SizeCube/2;l2++) {
   for (int l2=Inf;l2<Sup;l2++) {
      int l1 = -SizeCube/2;
      int x = (this->*f1)(l1,l2,a,b,c);
      int y = (this->*f2)(l1,l2,a,b,c);
      int z = (this->*f3)(l1,l2,a,b,c);      
      int xPix = x + SizeCube/2;
      int yPix = y + SizeCube/2;
      int zPix = z + SizeCube/2;        

      int xSym = (this->*f1)(l1,-l2,a,b,c);
      int ySym = (this->*f2)(l1,-l2,a,b,c);
      int zSym = (this->*f3)(l1,-l2,a,b,c);
      int xPixSym = xSym + SizeCube/2;
      int yPixSym = ySym + SizeCube/2;
      int zPixSym = zSym + SizeCube/2; 
      
              
      
      if (x == OUT_CUBE || y == OUT_CUBE || z == OUT_CUBE ||
          x==SizeCube/2 || y==SizeCube/2 || z==SizeCube/2 ||     
          xSym == OUT_CUBE || ySym == OUT_CUBE || zSym == OUT_CUBE ||
          xSym==SizeCube/2 || ySym==SizeCube/2 || zSym==SizeCube/2) {
            
            PlaneMask(l1+SizeCube/2,l2+SizeCube/2)=0;
            PlaneMask(l1+SizeCube/2,-l2+SizeCube/2)=0;
            
      } else {
  
              if (!Inverse) {
            
               if (DataCount2(xPix,yPix,zPix) > DataCount2(xPixSym,yPixSym,zPixSym)) {
                  xPix = xPixSym;
                  yPix = yPixSym;
                  zPix = zPixSym;
               }            
             
               PlaneMask(l1+SizeCube/2,l2+SizeCube/2)=1;
               CubeMask(xPix,yPix,zPix)++;
               (this->*set_f)(InPlane, OutPlane, NumLine, l1+SizeCube/2,
                              l2+SizeCube/2,xPix,yPix,zPix,PTakeSym);
                              
               PlaneMask(l1+SizeCube/2,-l2+SizeCube/2)=1;
               CubeMask(xPix,yPix,zPix)++;
               (this->*set_f)(InPlane, OutPlane, NumLine, l1+SizeCube/2,
                      -l2+SizeCube/2,xPix,yPix,zPix,True);
            } else {
               if (MemPlane(l1+SizeCube/2,l2+SizeCube/2,NumLine)>=0) {
                  comp_xyz (l1+SizeCube/2,l2+SizeCube/2,NumLine,xPix,yPix,zPix);
                  PlaneMask(l1+SizeCube/2,l2+SizeCube/2)=1;
                  CubeMask(xPix,yPix,zPix)++;
                  (this->*set_f)(InPlane, OutPlane, NumLine, l1+SizeCube/2,
                                 l2+SizeCube/2,xPix,yPix,zPix,PTakeSym);
               }
               
               if (MemPlane(l1+SizeCube/2,-l2+SizeCube/2,NumLine)>=0) {              
                  comp_xyz (l1+SizeCube/2,-l2+SizeCube/2,NumLine,xPix,yPix,zPix);              
                  PlaneMask(l1+SizeCube/2,-l2+SizeCube/2)=1;
                  CubeMask(xPix,yPix,zPix)++;
                  (this->*set_f)(InPlane, OutPlane, NumLine, l1+SizeCube/2,
                         -l2+SizeCube/2,xPix,yPix,zPix,True);        
               }     
            }                
         }     
   }
   
   
   int l1=-SizeCube/2;
   int l2=-SizeCube/2;
   {
      int x = (this->*f1)(l1,l2,a,b,c);
      int y = (this->*f2)(l1,l2,a,b,c);
      int z = (this->*f3)(l1,l2,a,b,c);      
      int xPix = x + SizeCube/2;
      int yPix = y + SizeCube/2;
      int zPix = z + SizeCube/2;           

      if (x == OUT_CUBE || y == OUT_CUBE || z == OUT_CUBE) {
         PlaneMask(l1+SizeCube/2,l2+SizeCube/2)=0;
      } else 
      if (x == SizeCube/2 || y == SizeCube/2 || z == SizeCube/2) {
         PlaneMask(l1+SizeCube/2,l2+SizeCube/2)=0;
      } else {
         if (!Inverse) {
            if (fabs (InPlane(xPix,yPix,zPix).imag()) <= EPS_IMAGPART) {
               PlaneMask(l1+SizeCube/2,l2+SizeCube/2)=1;
               CubeMask(xPix,yPix,zPix)++;
               (this->*set_f)(InPlane, OutPlane, NumLine, l1+SizeCube/2,
                              l2+SizeCube/2,xPix,yPix,zPix,PTakeSym);
            } else
            PlaneMask(l1+SizeCube/2,l2+SizeCube/2)=0;
         } else {
            if (MemPlane(l1+SizeCube/2,l2+SizeCube/2,NumLine)>=0) {
               comp_xyz (l1+SizeCube/2,l2+SizeCube/2,NumLine,xPix,yPix,zPix);
               (this->*set_f)(InPlane, OutPlane, NumLine, l1+SizeCube/2,
                              l2+SizeCube/2,xPix,yPix,zPix,PTakeSym);
            }
         }
      }
   }   
      
         
   l1=0;
   l2=-SizeCube/2;
   {
      int x = (this->*f1)(l1,l2,a,b,c);
      int y = (this->*f2)(l1,l2,a,b,c);
      int z = (this->*f3)(l1,l2,a,b,c);      
      int xPix = x + SizeCube/2;
      int yPix = y + SizeCube/2;
      int zPix = z + SizeCube/2;
                 
      if (x == OUT_CUBE || y == OUT_CUBE || z == OUT_CUBE) {
         PlaneMask(l1+SizeCube/2,l2+SizeCube/2)=0;
      } else 
      if (x == SizeCube/2 || y == SizeCube/2 || z == SizeCube/2) {
         PlaneMask(l1+SizeCube/2,l2+SizeCube/2)=0;
      } else {
         if (!Inverse) {
            if (fabs (InPlane(xPix,yPix,zPix).imag()) <= EPS_IMAGPART ) {
               PlaneMask(l1+SizeCube/2,l2+SizeCube/2)=1;
               CubeMask(xPix,yPix,zPix)++;
               (this->*set_f)(InPlane, OutPlane, NumLine, l1+SizeCube/2,
                              l2+SizeCube/2,xPix,yPix,zPix,PTakeSym);
            } else
            PlaneMask(l1+SizeCube/2,l2+SizeCube/2)=0; 
         } else {
            if (MemPlane(l1+SizeCube/2,l2+SizeCube/2,NumLine)>=0) {               
               comp_xyz (l1+SizeCube/2,l2+SizeCube/2,NumLine,xPix,yPix,zPix);
               (this->*set_f)(InPlane, OutPlane, NumLine, l1+SizeCube/2,
                              l2+SizeCube/2,xPix,yPix,zPix,PTakeSym);
            }
         }         
      } 
   }  
      
   l1=-SizeCube/2;
   l2=0;
   {
      int x = (this->*f1)(l1,l2,a,b,c);
      int y = (this->*f2)(l1,l2,a,b,c);
      int z = (this->*f3)(l1,l2,a,b,c);      
      int xPix = x + SizeCube/2;
      int yPix = y + SizeCube/2;
      int zPix = z + SizeCube/2;           
  
      if (x == OUT_CUBE || y == OUT_CUBE || z == OUT_CUBE) {
         PlaneMask(l1+SizeCube/2,l2+SizeCube/2)=0;
      } else 
      if (x == SizeCube/2 || y == SizeCube/2 || z == SizeCube/2) {
         PlaneMask(l1+SizeCube/2,l2+SizeCube/2)=0;
      } else {
         if (!Inverse) {
            if (fabs (InPlane(xPix,yPix,zPix).imag()) <= EPS_IMAGPART ) {
               PlaneMask(l1+SizeCube/2,l2+SizeCube/2)=1;
               CubeMask(xPix,yPix,zPix)++;
               (this->*set_f)(InPlane, OutPlane, NumLine, l1+SizeCube/2,
                              l2+SizeCube/2,xPix,yPix,zPix,PTakeSym);
            } else
            PlaneMask(l1+SizeCube/2,l2+SizeCube/2)=0;
         } else {
            if (MemPlane(l1+SizeCube/2,l2+SizeCube/2,NumLine)>=0) {               
               comp_xyz (l1+SizeCube/2,l2+SizeCube/2,NumLine,xPix,yPix,zPix);
               (this->*set_f)(InPlane, OutPlane, NumLine, l1+SizeCube/2,
                              l2+SizeCube/2,xPix,yPix,zPix,PTakeSym);
            }
         }          
      }   
   }
   
   if (WriteFile) {
      char Name[20];
      sprintf (Name, "Plane%d", NumLine);
      //fits_write_intarr(Name, PlaneMask);   
      sprintf (Name, "Cube%d", NumLine);
      //fits_write_intarr(Name, CubeMask);   
   }  
}



void PartialRadon3D::odd_make_plane (cfarray& InPlane, cfarray& OutPlane,
                                     tp_mfunc f1, tp_mfunc f2, tp_mfunc f3, 
                                     double a, double b, double c,
                                     int NumLine, tp_setfunc set_f, 
                                     Bool Inverse)  {

   // local var
   intarray PlaneMask (SizeCube, SizeCube);
   intarray CubeMask  (SizeCube, SizeCube, SizeCube);
   
   int MidPixVal = (SizeCube-1)/2;
           
   // symetrical part       
   for (int l1=-MidPixVal;l1<=MidPixVal;l1++)
   for (int l2=l1+1;l2<=MidPixVal;l2++) {
   
      int x = (this->*f1)(l1,l2,a,b,c);
      int y = (this->*f2)(l1,l2,a,b,c);
      int z = (this->*f3)(l1,l2,a,b,c);               
      int xPix = x + MidPixVal;
      int yPix = y + MidPixVal;
      int zPix = z + MidPixVal;

      if (x == OUT_CUBE || y == OUT_CUBE || z == OUT_CUBE) {
         PlaneMask(l1+MidPixVal,l2+MidPixVal)=0;
         PlaneMask(-l1+MidPixVal,-l2+MidPixVal)=0;
      } else {
         PlaneMask(l1+MidPixVal,l2+MidPixVal)=1;
         PlaneMask(-l1+MidPixVal,-l2+MidPixVal)=1;
         
         CubeMask(xPix,yPix,zPix)++;
         (this->*set_f)(InPlane, OutPlane, NumLine, l1+MidPixVal,
                        l2+MidPixVal,xPix,yPix,zPix,PTakeSym);
                              
         xPix = -x + MidPixVal;
         yPix = -y + MidPixVal;
         zPix = -z + MidPixVal; 
                   
         CubeMask(xPix,yPix,zPix)++;
         (this->*set_f)(InPlane, OutPlane, NumLine, -l1+MidPixVal,
                           -l2+MidPixVal,xPix,yPix,zPix,PTakeSym);                        
      }
   }
   
   // symetrical diagonal
   for (int l1=-MidPixVal;l1<0;l1++) {
      int l2=l1;
      
      int x = (this->*f1)(l1,l2,a,b,c);
      int y = (this->*f2)(l1,l2,a,b,c);
      int z = (this->*f3)(l1,l2,a,b,c);               
      int xPix = x + MidPixVal;
      int yPix = y + MidPixVal;
      int zPix = z + MidPixVal;
      
      // central point used two times....
      if (x == OUT_CUBE || y == OUT_CUBE || z == OUT_CUBE) {
         PlaneMask(l1+MidPixVal,l2+MidPixVal)=0;
         PlaneMask(-l1+MidPixVal,-l2+MidPixVal)=0;
      } else {
         PlaneMask(l1+MidPixVal,l2+MidPixVal)=1;
         PlaneMask(-l1+MidPixVal,-l2+MidPixVal)=1;
         
         CubeMask(xPix,yPix,zPix)++;
         (this->*set_f)(InPlane, OutPlane, NumLine, l1+MidPixVal,
                        l2+MidPixVal,xPix,yPix,zPix,PTakeSym);
                        
         xPix = -x + MidPixVal;
         yPix = -y + MidPixVal;
         zPix = -z + MidPixVal;                    
         
         CubeMask(xPix,yPix,zPix)++;
         (this->*set_f)(InPlane, OutPlane, NumLine, -l1+MidPixVal,
                        -l2+MidPixVal,xPix,yPix,zPix,PTakeSym);                
      }
   }
   
   //central point x,y,z must be zero!!!
   int l1=0,l2=0;
                  
   int xPix = 0 + MidPixVal;
   int yPix = 0 + MidPixVal;
   int zPix = 0 + MidPixVal;   
  
   PlaneMask(l1+MidPixVal,l2+MidPixVal)=1;
         
   CubeMask(xPix,yPix,zPix)++;
   (this->*set_f)(InPlane, OutPlane, NumLine, l1+MidPixVal,
                  l2+MidPixVal,xPix,yPix,zPix,PTakeSym); 
     
     
   if (WriteFile) {
      char Name[20];
      sprintf (Name, "Plane%d", NumLine);
      //fits_write_intarr(Name, PlaneMask);   
      sprintf (Name, "Cube%d", NumLine);
      //fits_write_intarr(Name, CubeMask);   
   }  
}



/******************************************************************************/
// plane c=0, b=0
/******************************************************************************/
inline int PartialRadon3D::fx_uXv8w8 (int l1,int l2, double a, double b, double c) {
   return (0);                                  
}
inline int PartialRadon3D::fy_uXv8w8 (int l1,int l2, double a, double b, double c) {
   return (l1);                                  
}
inline int PartialRadon3D::fz_uXv8w8 (int l1,int l2, double a, double b, double c) {
   return (l2);                                  
}



/******************************************************************************/
// plane c=0
/******************************************************************************/
inline int PartialRadon3D::fx1_uXvXw8 (int l1,int l2, double a, double b, double c) {
   return (l1);                                  
}
inline int PartialRadon3D::fy1_uXvXw8 (int l1,int l2, double a, double b, double c) {
   int Ret = double2int((-(a/b)*l1));        
   if      (Ret>SizeCube/2) Ret=OUT_CUBE;
   else if (Ret<-SizeCube/2) Ret=OUT_CUBE;
   return Ret;                                  
}
inline int PartialRadon3D::fz1_uXvXw8 (int l1,int l2, double a, double b, double c) {
   return (l2);                  
}

inline int PartialRadon3D::fx2_uXvXw8 (int l1,int l2, double a, double b, double c) {
   int Ret = double2int((-b/a)*l1);        
   if      (Ret>SizeCube/2) Ret=OUT_CUBE;
   else if (Ret<-SizeCube/2) Ret=OUT_CUBE;
   return Ret;                                    
}
inline int PartialRadon3D::fy2_uXvXw8 (int l1,int l2, double a, double b, double c) {
   return (l1);                                  
}
inline int PartialRadon3D::fz2_uXvXw8 (int l1,int l2, double a, double b, double c) {
   return (l2);                          
}



/******************************************************************************/
// plane 
/******************************************************************************/
inline int PartialRadon3D::fx1_uXvXwX (int l1,int l2, double a, double b, double c) {                            
   int Ret = double2int((-b*l1-c*l2)/a);        
   if      (Ret>SizeCube/2) Ret=OUT_CUBE;
   else if (Ret<-SizeCube/2) Ret=OUT_CUBE;
   return Ret;
}
inline int PartialRadon3D::fy1_uXvXwX (int l1,int l2, double a, double b, double c) {
   return (l1);                                  
}
inline int PartialRadon3D::fz1_uXvXwX (int l1,int l2, double a, double b, double c) {
   return (l2);                          
}

inline int PartialRadon3D::fx2_uXvXwX (int l1,int l2, double a, double b, double c) {                              
   return (l1);
}
inline int PartialRadon3D::fy2_uXvXwX (int l1,int l2, double a, double b, double c) {
   int Ret = double2int((-a*l1-c*l2)/b);
   if      (Ret>SizeCube/2) Ret=OUT_CUBE;
   else if (Ret<-SizeCube/2) Ret=OUT_CUBE;
   return Ret;                                    
}
inline int PartialRadon3D::fz2_uXvXwX (int l1,int l2, double a, double b, double c) {
   return (l2);                          
}

inline int PartialRadon3D::fx3_uXvXwX (int l1,int l2, double a, double b, double c) {                              
   return (l1);
}
inline int PartialRadon3D::fy3_uXvXwX (int l1,int l2, double a, double b, double c) {
   return (l2);                                   
}
inline int PartialRadon3D::fz3_uXvXwX (int l1,int l2, double a, double b, double c) {
   int Ret = double2int((-a*l1-b*l2)/c);        
   if      (Ret>SizeCube/2) Ret=OUT_CUBE;
   else if (Ret<-SizeCube/2) Ret=OUT_CUBE;
   return Ret;                              
}




/******************************************************************************
**                   Copyright (C) 2001 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Arnaud Woiselle
**
**    Date:  27/05/2009
**    
**    File:  Project_3d2d.h
**
**    Modification history :
**
*******************************************************************************
**
**    DESCRIPTION  see .h
**    -----------  
**
******************************************************************************/

#include "Project_3d2d.h"


void Project_3d2d::comp_odd_Plane (cfarray& F3DCube, cfarray& FPartRad3DImage,
                                     Bool Inverse) { 
//   cerr<<"Project_3d2d::Comp odd plane"<<endl;

// Estimate the numLine corresponding to a,b,c_ref
	int numLine;
	for (int i=0; i<ProjectNbr; i++)
		if((xPos(i)-SizeCube/2 == a_ref)
		&& (yPos(i)-SizeCube/2 == b_ref)
		&& (zPos(i)-SizeCube/2 == c_ref))
			numLine = i;

//cerr<<"Numline "<<numLine<<"/"<<ProjectNbr<<"xyz="<<xPos(numLine)<<","<<yPos(numLine)<<","<<zPos(numLine)<<endl;
/*	if(Inverse)
	{
		for(int i=0;i<SizeCube;i++)
		for(int j=0;j<SizeCube;j++)
			FPartRad3DImage(i,j,0) = FPartRad3DImage(i,j,0);
	}
*/



//   for (int numLine=0; numLine<ProjectNbr; numLine++) 
   for (int l1=0; l1<SizeCube; l1++)
   for (int l2=0; l2<SizeCube; l2++) {   
      
      int xPix,yPix,zPix;
      if (!Inverse && MemPlane(l1,l2,numLine)>=0) {
         
         comp_xyz (l1,l2,numLine,xPix,yPix,zPix);
         FPartRad3DImage(l1,l2,0) = F3DCube(xPix,yPix,zPix);
         
      } 
      else 
      if (Inverse && MemPlane(l1,l2,numLine)>=0) {
      
         comp_xyz (l1,l2,numLine,xPix,yPix,zPix);
         F3DCube(xPix,yPix,zPix) += FPartRad3DImage(l1,l2,0); 
         DataCount3(xPix,yPix,zPix)++;  
      }
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
   
   
/*   
	if(!Inverse)
	{
		for(int i=0;i<SizeCube;i++)
		for(int j=0;j<SizeCube;j++)
			FPartRad3DImage(i,j,0) = FPartRad3DImage(i,j,numLine);
	}
*/
//   cerr<<" end Project_3d2d::Comp odd plane"<<endl;
}


void Project_3d2d::alloc () {
   Alloc = True;
   
   CFFT3D.CenterZeroFreq = True;

   F3DCube.alloc (SizeCube, SizeCube, SizeCube);
   FPartRad3DImage.alloc (SizeCube, SizeCube, 1);
                                                                           
   //DataCount1.alloc (SizeCube,SizeCube,ProjectNbr); 

   if (SizeCube%2==0) {
      DataCount2.alloc (SizeCube,SizeCube,SizeCube);
      MemPlane.alloc (SizeCube,SizeCube,ProjectNbr);
      MemPlane.init (-1.0);
   }   
                                
   DataCount3.alloc (SizeCube, SizeCube, SizeCube);   
}

void Project_3d2d::comp_FFT2D (fltarray& PartRad3DImage, 
                                 cfarray& FPartRad3DImage) {
	if (Control) fits_write_fltarr ((char*)"InFFT2D", PartRad3DImage);

	FFTN_2D FFT2D;
	FFT2D.CenterZeroFreq = True;
	Icomplex_f CmpCube(SizeCube,SizeCube);  

	int numPlane=0; 
	for (int i=0;i<SizeCube;i++)
	for (int j=0;j<SizeCube;j++)
		CmpCube(j,i) = PartRad3DImage(i,j,numPlane);
	FFT2D.fftn2d (CmpCube, False);
	for (int i=0;i<SizeCube;i++)
	for (int j=0;j<SizeCube;j++)
		FPartRad3DImage(i,j,numPlane) = CmpCube(j,i);

	if (Control) WriteCmpArr ((char*)"OutFFT2D", FPartRad3DImage);
}   


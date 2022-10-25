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
**    Date:  03/02/00
**    
**    File:  ridgelet.cc
**
*******************************************************************************
**
**    DESCRIPTION  ridgelet transform program
**    ----------- 
**                 
******************************************************************************/

#include "DefPoint.h"
#include "IM3D_IO.h"
#include<iostream>
#include "macro.h"

#define  CUBE_FABS(x)    (fabs(x) * fabs(x) * fabs(x)) 

char Name_Dat_In[256];  // input file image  
char Name_Dat_Out[256]; //  output file name  
 
extern int  OptInd;
extern char *OptArg;

extern int  GetOpt(int argc, char **argv, char *opts);

float BinCat = 1.;
Bool MakeProjV0 = True;
Bool Verbose=False;         // verbose mode
Bool TrVerbose=False;
Bool TrDebug=False;
char TraceFileName[80];
Bool WithTrace=False;       // no file trace
int NbrPoint=0;

int NoX=-1;
int NoY=-1;
int NoZ=-1;
Bool OptDim=False;
/***************************************/

static void usage(char *argv[]) {
    
   fprintf(OUTMAN, "Usage: %s options catalog.cat catalog.fits\n\n", argv[0]);
   fprintf(OUTMAN, "   where options =  \n");    
    
   fprintf(OUTMAN, "         [-b bin] (only for catalogue)\n");
   fprintf(OUTMAN, "              bin gives resolution ( default=1 ).\n");    
   
   fprintf(OUTMAN, "         [-p ]\n");
   fprintf(OUTMAN, "              Do not project in V0.\n");    

   fprintf(OUTMAN, "         [-n NbrPoints]\n");
   fprintf(OUTMAN, "              Keep randomly only NbrPoints of the input catalog.\n"); 
   fprintf(OUTMAN, "              By default, all points in the catalog are used.\n");  
   
   fprintf(OUTMAN, "         [-D Nx,Ny,Nz]\n");
   fprintf(OUTMAN, "              Force the output cube to have this given dimension (Nx,Ny,Nz) .\n"); 
   fprintf(OUTMAN, "              Points outside this volume are not taken into account.\n"); 
   fprintf(OUTMAN, "              By default, the cube dimensions are automatically calculated from the input catalog.\n");
   // manline();   
   vm_usage();
  //  manline();  
      
   verbose_usage();
   // manline();         
    
  exit(-1);
}
  
/*********************************************************************/

/* GET COMMAND LINE ARGUMENTS */
static void filtinit(int argc, char *argv[]) {
   
   int c;
#ifdef LARGE_BUFF
   int VMSSize=-1;
   Bool OptZ = False;
   char VMSName[1024] = "";
#endif     

   /* get options */
   while ((c = GetOpt(argc,argv,"D:n:b:dpvzZ:")) != -1) {
      
      switch (c) {
      case 'D':
            if (sscanf(OptArg,"%d,%d,%d",&NoX, &NoY, &NoZ) != 3){
            fprintf(OUTMAN, "bad cube dimension : %s\n", OptArg);
            exit(-1);
         } 
	 OptDim=True;
         break;	
      case 'n':
            if (sscanf(OptArg,"%d",&NbrPoint) != 1){
            fprintf(OUTMAN, "bad number of points : %s\n", OptArg);
            exit(-1);
         }
         break;	
      case 'v': Verbose = True;break;
      case 'd':	TrDebug = True; TrVerbose = True; break; 
      case 'p':	MakeProjV0 = False; break; 
      case 'b': /* -B <bin> */
         if ((sscanf(OptArg,"%f",&BinCat) != 1)|| (BinCat < 0)){
            fprintf(OUTMAN, "bad value of bin: %s\n", OptArg);
            exit(-1);
         }
         break;	
      case 'T': /* abaque file */
         if (sscanf(OptArg,"%s", TraceFileName) != 1) {
            fprintf(stderr, "Error: bad trace file name: %s\n", OptArg); 
	    exit(-1);
         }
	 WithTrace = True;
      break;  
         	 
#ifdef LARGE_BUFF
      case 'z':
         if (OptZ == True) {
            fprintf(OUTMAN, "Error: Z option already set...\n");
            exit(-1);
         }
	 OptZ = True;
	 break;
      case 'Z':
	 if (sscanf(OptArg,"%d:%s",&VMSSize, VMSName) < 1) {
	    fprintf(OUTMAN, "Error: syntaxe is Size:Directory ... \n");
	    exit(-1);
         }
	 if (OptZ == True) {
            fprintf(OUTMAN, "Error: z option already set...\n");
            exit(-1);
         }
	 OptZ = True;
         break;
#endif
      case '?': usage(argv); break;
      default: usage(argv); break;
      }
   } 
      
   /* get otrace_paramptional input file names from trailing parameters and open files */
   if (OptInd < argc) strcpy(Name_Dat_In, argv[OptInd++]);
   else usage(argv);

   if (OptInd < argc) strcpy(Name_Dat_Out, argv[OptInd++]);
   else usage(argv);

   /* make sure there are not too many parameters */
   if (OptInd < argc) {
      fprintf(OUTMAN, "Too many parameters: %s ...\n", argv[OptInd]);
      exit(-1);
   }

       	   
#ifdef LARGE_BUFF
   if (OptZ == True) vms_init(VMSSize, VMSName, Verbose);
#endif  
}
 
/***************************************************************************/
 
int main(int argc, char *argv[]) {
   
   // init local var
   //---------------
   fitsstruct Header;
   char Cmd[512];
   Cmd[0] = '\0';
   for (int k =0; k < argc; k++) sprintf(Cmd, "%s %s", Cmd, argv[k]);
   
   // Get command line arguments, open input file(s) if necessary
   lm_check(LIC_MR4);
   filtinit(argc, argv);
   
   if (TrVerbose) TRACEUR.TraceVerbose();
   if (TrDebug) TRACEUR.TraceDebug();
   if (TrVerbose) TRACEUR.SetScreenFlux(true);
   if (TrVerbose)  if (WithTrace) TRACEUR.SetFileFlux(true, TraceFileName);
   if (TrVerbose) TRACEUR.SetScreenFlux(true);
      
   char Msg[140];
   if (TrVerbose)
   { 
      TRACEUR << new Trace (Trace::METIER, "Main", "info param main",
                            " -- [-v]  : Verbose" , NULL, Trace::VERBOSE);
       sprintf (Msg, " -- [-B :] : %f ", BinCat);			    
      TRACEUR << new Trace (Trace::METIER, "Main", "info param main",
                            Msg, NULL, Trace::VERBOSE);
     sprintf (Msg, "File Name in = %s", Name_Dat_In);			    
     TRACEUR << new Trace (Trace::METIER, "Main", "info param main",
                         Msg, NULL, Trace::VERBOSE);
     sprintf (Msg, "File Name out = %s", Name_Dat_Out);				        
    TRACEUR << new Trace (Trace::METIER, "Main", "info param main",
                         Msg, NULL, Trace::VERBOSE);			    
   }
   
   // read file in, must be *.cat
   //----------------------------
   if (     strstr(Name_Dat_In,".fits") 
         || strstr(Name_Dat_In,".Fits")
         || strstr(Name_Dat_In,".FITS")) {
      cout << "file must be catalog (*.cat)" << endl;
      exit(-1);
   }
   
   // init local var
   //---------------
   ArrayPoint Catalogue;
   Catalogue.read(Name_Dat_In, TrVerbose);    
   if(Catalogue.dim()!=3){
      cout << "catalog must have 3 dimentions" << endl;
      exit(-1);
   }		
   
   if (TrVerbose)  Catalogue.trace(); 
   
   //compute size of fits files
   //--------------------------
   float xDyn = ( Catalogue.Pmax.x()-Catalogue.Pmin.x() );
   int xSize = int( xDyn / BinCat ) + 1;
   float yDyn = ( Catalogue.Pmax.y()-Catalogue.Pmin.y() );
   int ySize = int( yDyn / BinCat ) + 1;
   float zDyn = ( Catalogue.Pmax.z()-Catalogue.Pmin.z() );
   int zSize = int( zDyn / BinCat ) + 1;
   
   float xRealMin = Catalogue.Pmin.x();	  
   float yRealMin = Catalogue.Pmin.y();	  
   float zRealMin = Catalogue.Pmin.z();	  
   
   fltarray xCenter( xSize );
   for( int i=0; i<xSize; ++i ) xCenter(i) = ( xRealMin + BinCat/2. ) + i*BinCat;
   fltarray yCenter( ySize );
   for( int i=0; i<ySize; ++i ) yCenter(i) = ( yRealMin + BinCat/2. ) + i*BinCat;
   fltarray zCenter( zSize );
   for( int i=0; i<zSize; ++i ) zCenter(i) = ( zRealMin + BinCat/2. ) + i*BinCat;
   
   if(TrVerbose) {   
      sprintf (Msg, "dynamique [%f,%f,%f]", xDyn, yDyn, zDyn);			    
      TRACEUR << new Trace (Trace::METIER, "main", "info",
                            Msg , NULL, Trace::VERBOSE);
      sprintf (Msg, "(xMin,xMin+bin*Nx,(xMax)) [%f,%f,(%f)]", 
               Catalogue.Pmin.x(), Catalogue.Pmin.x()+BinCat*xDyn, Catalogue.Pmax.x());		    
      TRACEUR << new Trace (Trace::METIER, "main", "info",
                            Msg , NULL, Trace::VERBOSE);
      sprintf (Msg, "(yMin,yMin+bin*Ny,(yMax)) [%f,%f,(%f)]", 
               Catalogue.Pmin.y(), Catalogue.Pmin.y()+BinCat*yDyn, Catalogue.Pmax.y());			    
      TRACEUR << new Trace (Trace::METIER, "main", "info",
                            Msg , NULL, Trace::VERBOSE);
      sprintf (Msg, "(zMin,zMin+bin*Nz,(zMax)) [%f,%f,(%f)]", 
               Catalogue.Pmin.z(), Catalogue.Pmin.z()+BinCat*zDyn, Catalogue.Pmax.z());		    
      TRACEUR << new Trace (Trace::METIER, "main", "info",
                            Msg , NULL, Trace::VERBOSE);
      sprintf (Msg, "New size cube out [%d,%d,%d]", xSize, ySize, zSize );			    
      TRACEUR << new Trace (Trace::METIER, "main", "info",
                            Msg , NULL, Trace::VERBOSE);
   } 
   
   
      int Npx=(int)((Catalogue.Pmax.x()-Catalogue.Pmin.x())/BinCat+0.5);
      int Npy=(int)((Catalogue.Pmax.y()-Catalogue.Pmin.y())/BinCat+0.5);
      int Npz=(int)((Catalogue.Pmax.z()-Catalogue.Pmin.z())/BinCat+0.5);
      
      if(TrVerbose) {
         sprintf (Msg, "Old size cube out [%d,%d,%d]", Npx, Npy, Npz);			    
         TRACEUR << new Trace (Trace::METIER, "main", "info",
                               Msg , NULL, Trace::VERBOSE);
      } 
      
    
      
   fltarray EventCatalogue;
   fltarray ProjV0;   
   if (OptDim == False)	 
   {
      EventCatalogue.resize( xSize, ySize, zSize );
      ProjV0.resize( xSize, ySize, zSize );
   }
   else
    {
      EventCatalogue.resize( NoX, NoY, NoZ );
      ProjV0.resize( NoX, NoY, NoZ );
   }
   if (Verbose == True)
   {
       cout << "Filename in = " << Name_Dat_In << endl;
       cout << "Filename out = " << Name_Dat_Out  << endl;
       if (MakeProjV0 == False) cout << "Do not apply a projection in V0" << endl;
       cout << "BinCat = " << BinCat << endl;
       if (NbrPoint > 0) cout << "Keep only " << NbrPoint << " in the catalog " <<  endl;
       cout << "Catalog cube size = " << xSize << " , " << ySize << " , " << zSize  << endl;
       if (OptDim==True)
          cout << "Output cube size = " << NoX << " , " << NoY << " , " << NoZ  << endl;
   }
 
   ProjV0.init();
   if (NbrPoint < 0) NbrPoint = 0;
   if (NbrPoint > Catalogue.np()) NbrPoint = Catalogue.np();
   float PerCentKill = (NbrPoint == 0) ? 0.: 1. - (float) NbrPoint / (float) Catalogue.np();
   if (PerCentKill > 1) PerCentKill = 1.;

   init_random ();
     
   if (MakeProjV0) {
      
      // projection on V0...
      //--------------------
      
      // compute bspline on SLINESIZE pts
      //---------------------------------
      int SLINESIZE = 401;
      dblarray bspline( SLINESIZE );
      for( int i=-SLINESIZE/2; i<=SLINESIZE/2; i++) {
         float x = 2. * float(i) / ( float(SLINESIZE/2) ) ;
         bspline( i+SLINESIZE/2 ) = (           CUBE_FABS( x-2. ) 
	                                  - 4.* CUBE_FABS( x-1. )
				          + 6.* CUBE_FABS( x )
				          - 4.* CUBE_FABS( x+1.)
				          +     CUBE_FABS( x+2. ) ) / 12.; 
      } 
      
      
      // read events on catalog and project them
      for(int i=0;i<Catalogue.np();i++) {
      
         float RND = get_random();
	 
	 if(TrVerbose)
	  {
	    sprintf( Msg, "RND:%f, percentkill:%f", RND, PerCentKill );
	     TRACEUR << new Trace (Trace::METIER, "main", "info",
                               Msg , NULL, Trace::VERBOSE);
	 }
	 
	 // current point is kept ?
	 //------------------------
	 if (RND > PerCentKill) {
	 
	    float xRealVal  = Catalogue(i).x();
	    int   xPixelInd = int( ( xRealVal - xRealMin ) / BinCat );
	    float yRealVal  = Catalogue(i).y();
	    int   yPixelInd = int( ( yRealVal - yRealMin ) / BinCat );
	    float zRealVal  = Catalogue(i).z();
	    int   zPixelInd = int( ( zRealVal - zRealMin ) / BinCat );
	    
	    if(TrVerbose) 
	    {
	    sprintf( Msg, "x val:%f, pixel indice:%d, pixel lim:[%f,(c:%f),%f]",
	             xRealVal, xPixelInd, xRealMin+xPixelInd*BinCat, 
		     xCenter(xPixelInd), xRealMin+(xPixelInd+1)*BinCat );
	    TRACEUR << new Trace (Trace::METIER, "main", "info",
                                  Msg , NULL, Trace::VERBOSE);
	    sprintf( Msg, "y val:%f, pixel indice:%d, pixel lim:[%f,(c:%f),%f]",
	             yRealVal, yPixelInd, yRealMin+yPixelInd*BinCat, 
		     yCenter(yPixelInd), yRealMin+(yPixelInd+1)*BinCat );
	    TRACEUR << new Trace (Trace::METIER, "main", "info",
                                  Msg , NULL, Trace::VERBOSE);
	    sprintf( Msg, "z val:%f, pixel indice:%d, pixel lim:[%f,(c:%f),%f]",
	             zRealVal, zPixelInd, zRealMin+zPixelInd*BinCat, 
		     zCenter(zPixelInd), zRealMin+(zPixelInd+1)*BinCat );
	    TRACEUR << new Trace (Trace::METIER, "main", "info",
                                  Msg , NULL, Trace::VERBOSE);
	   }
	   
           // add event in event catalog
            if (OptDim == False) EventCatalogue( xPixelInd, yPixelInd, zPixelInd )++;   
	    else
	    {
	       int pos_x = xPixelInd;
	       int pos_y = yPixelInd;
	       int pos_z = zPixelInd;
	       if ((pos_x >= 0) && (pos_x < NoX) && (pos_y >= 0) && (pos_y < NoY )  && (pos_z>= 0) && (pos_z < NoZ) ) 
	       {
	          EventCatalogue(pos_x, pos_y, pos_z) ++; 
 	       }
	    }
	      
	         
           
           // compute offset from center
	   //---------------------------
	   int xBeg=-1,xEnd=-1;
	   float xOffset = xRealVal - xCenter(xPixelInd);
	   if( fabs( xOffset ) <= 1e-06 ) {xBeg = -1; xEnd = 1;}
	   else if( xOffset > 0 )         {xBeg = -2; xEnd = 1;}
	   else if( xOffset < 0 )         {xBeg = -1; xEnd = 2;}
	   
	    if(TrVerbose) 
	    {
	   sprintf( Msg, "x Offset:%f, [xBeg,xEnd] : [%d;%d]", 
	            xOffset, xBeg, xEnd );
	   TRACEUR << new Trace (Trace::METIER, "main", "info",
                                 Msg , NULL, Trace::VERBOSE);
	   }
	      
	   int yBeg=-1,yEnd=-1;
	   float yOffset = yRealVal - yCenter(yPixelInd);
	   if( fabs( yOffset ) <= 1e-06 ) {yBeg = -1; yEnd = 1;}
	   else if( yOffset > 0 )         {yBeg = -2; yEnd = 1;}
	   else if( yOffset < 0 )         {yBeg = -1; yEnd = 2;}
	   
	   if(TrVerbose) 
	    {
	   sprintf( Msg, "y Offset:%f, [yBeg,yEnd] : [%d;%d]", 
	            yOffset, yBeg, yEnd );
	   TRACEUR << new Trace (Trace::METIER, "main", "info",
                                 Msg , NULL, Trace::VERBOSE);
 	   }
	   
	   int zBeg=-1,zEnd=-1;
	   float zOffset = zRealVal - zCenter(zPixelInd);
	   if( fabs( zOffset ) <= 1e-06 ) {zBeg = -1; zEnd = 1;}
	   else if( zOffset > 0 )         {zBeg = -2; zEnd = 1;}
	   else if( zOffset < 0 )         {zBeg = -1; zEnd = 2;}
	   
	    if(TrVerbose) 
	    {
	   sprintf( Msg, "z Offset:%f, [zBeg,zEnd] : [%d;%d]", 
	            zOffset, zBeg, zEnd );
	   TRACEUR << new Trace (Trace::METIER, "main", "info",
                                 Msg , NULL, Trace::VERBOSE);
	   }
	   //bounds tests
	   while( xPixelInd+xBeg < 0 ) xBeg++;
	   while( xPixelInd+xEnd > xSize-1 ) xEnd--;
	   while( yPixelInd+yBeg < 0 ) yBeg++;
	   while( yPixelInd+yEnd > ySize-1 ) yEnd--;
	   while( zPixelInd+zBeg < 0 ) zBeg++;
	   while( zPixelInd+zEnd > zSize-1 ) zEnd--;
	   
	   if(TrVerbose) 
	    {
	   sprintf( Msg, "corrected bound x[%d,%d], y[%d,%d], z[%d,%d]", 
	            xBeg, xEnd, yBeg, yEnd, zBeg, zEnd );
	   TRACEUR << new Trace (Trace::METIER, "main", "info",
                                 Msg , NULL, Trace::VERBOSE);
	   }
	   
	   // proj in V0
           for( int xInd=xBeg; xInd<=xEnd; ++xInd )
           for( int yInd=yBeg; yInd<=yEnd; ++yInd )
           for( int zInd=zBeg; zInd<=zEnd; ++zInd ) {
	   
	     if(TrVerbose) 
	     {
	      sprintf( Msg, "== point %d, pixel pos [%d+%d,%d+%d,%d+%d]",
	               i, xPixelInd, xInd, yPixelInd, yInd, zPixelInd, zInd );
	      TRACEUR << new Trace (Trace::METIER, "main", "info",
                                    Msg , NULL, Trace::VERBOSE);
              }
	      int xCur = int( SLINESIZE/2+1 + xInd*( SLINESIZE/4 ) 
	                      + xOffset*( SLINESIZE/4 ) );
	      int yCur = int( SLINESIZE/2+1 + yInd*( SLINESIZE/4 ) 
	                      + yOffset*( SLINESIZE/4 ) );
	      int zCur = int( SLINESIZE/2+1 + zInd*( SLINESIZE/4 ) 
	                      + zOffset*( SLINESIZE/4 ) );
	      
	      if(TrVerbose) 
	     {
	      sprintf( Msg, "xInd:%d, in [%f,%f], val=%d",
	               xInd, SLINESIZE/2+1 + xInd*( SLINESIZE/4 ) - SLINESIZE/8.,
		       SLINESIZE/2+1 + xInd*( SLINESIZE/4 ) + SLINESIZE/8., xCur );   
	      TRACEUR << new Trace (Trace::METIER, "main", "info",
                                    Msg , NULL, Trace::VERBOSE);
	      sprintf( Msg, "yInd:%d, in [%f,%f], val=%d",
	               yInd, SLINESIZE/2+1 + yInd*( SLINESIZE/4 ) - SLINESIZE/8.,
		       SLINESIZE/2+1 + yInd*( SLINESIZE/4 ) + SLINESIZE/8., yCur );   
	      TRACEUR << new Trace (Trace::METIER, "main", "info",
                                    Msg , NULL, Trace::VERBOSE);
	      sprintf( Msg, "zInd:%d, in [%f,%f], val=%d",
	               zInd, SLINESIZE/2+1 + zInd*( SLINESIZE/4 ) - SLINESIZE/8.,
		       SLINESIZE/2+1 + zInd*( SLINESIZE/4 ) + SLINESIZE/8., zCur );   
	      TRACEUR << new Trace (Trace::METIER, "main", "info",
                                    Msg , NULL, Trace::VERBOSE);
             }   
	      double xspline = bspline(xCur);
	      double yspline = bspline(yCur);
	      double zspline = bspline(zCur);
	      
	      if(TrVerbose) 
	     {
	      sprintf( Msg, "xSpline:%lf, ySpline:%lf, zSpline:%lf", 
	               xspline, yspline, zspline );
	      TRACEUR << new Trace (Trace::METIER, "main", "info",
                                    Msg , NULL, Trace::VERBOSE);
	      
	      }
	      if (OptDim == False)  ProjV0( xPixelInd+xInd, yPixelInd+yInd, zPixelInd+zInd )  += float( xspline*yspline*zspline ); 
	      else
	      {
	          int pos_x = xPixelInd+xInd;
		  int pos_y = yPixelInd+yInd;
		  int pos_z = zPixelInd+zInd;
	          if ((pos_x >= 0) && (pos_x < NoX) && (pos_y >= 0) && (pos_y < NoY )  && (pos_z>= 0) && (pos_z < NoZ) ) 
	          {
	             ProjV0(pos_x, pos_y, pos_z) += float( xspline*yspline*zspline); 
 	          }
	      }
	      
	       if(TrVerbose) 
	     {
	      sprintf( Msg, "ProjV0[%d,%d,%d]=%f", xPixelInd+xInd, yPixelInd+yInd, 
	               zPixelInd+zInd, 
		       ProjV0( xPixelInd+xInd, yPixelInd+yInd, zPixelInd+zInd )
	              );
	      TRACEUR << new Trace (Trace::METIER, "main", "info",
                                    Msg , NULL, Trace::VERBOSE);
	     }
	    } 
         } 
      }
      
   } else {
      // no projection in V0
      //--------------------
      
      for(int i=0;i<Catalogue.np();i++)
      {
         float RND = get_random();
	 if (RND > PerCentKill)
	 {
            int pos_x=(int)((Catalogue(i).x()-Catalogue.Pmin.x())/BinCat);
            int pos_y=(int)((Catalogue(i).y()-Catalogue.Pmin.y())/BinCat);
            int pos_z=(int)((Catalogue(i).z()-Catalogue.Pmin.z())/BinCat);
	 
	    // a little bad trick.....
	    if (pos_x==Npx) pos_x=(int)(Npx-1.e-06);
	    if (pos_y==Npy) pos_y=(int)(Npy-1.e-06);
	    if (pos_z==Npz) pos_z=(int)(Npz-1.e-06);
	 
	      	
            if (    pos_x<0 || pos_x>=Npx 
	      || pos_y<0 || pos_y>=Npy 
	      || pos_z<0 || pos_z>=Npz ) 
	   {
	    cout << "event number : " << i 
	         << "[" << pos_x << "," << pos_y << "," << pos_z << "]" << endl;
	    cout << "Error in coordinate... " <<endl; exit(-1);
	  }
	  
	  if (OptDim == False)
	  {
             ProjV0 (pos_x,pos_y,pos_z)++;
             EventCatalogue (pos_x,pos_y,pos_z)++;
	  }
	  else
	  {
	       if (    (pos_x >= 0) && (pos_x < NoX)
	      && (pos_y >= 0) && (pos_y < NoY )
	      && (pos_z>= 0) && (pos_z < NoZ) ) 
	   {
	      ProjV0 (pos_x,pos_y,pos_z)++;
              EventCatalogue (pos_x,pos_y,pos_z)++;
	   }
	  }
	    
        }
     }    
   }
   // ProjV0.info("proj");
   Header.hd_fltarray(ProjV0);
   Header.origin = Cmd;
   // io_3d_write_data("EventCat.fits", EventCatalogue);   
   io_3d_write_data(Name_Dat_Out, ProjV0, &Header); 
   
   if(TrVerbose) TRACEUR.FermerFichierDeTraces();
   exit(0);
}

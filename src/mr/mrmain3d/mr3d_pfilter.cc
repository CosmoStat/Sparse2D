/******************************************************************************
**                   Copyright (C) 2001 by CEA
*******************************************************************************
*******************************************************************************
**
**    DESCRIPTION  multiresolution iterative transform  for cubes
**    ----------- 
**                 
******************************************************************************/
#include "mr3d_pfilter.h"
#include "Memory.h"

Bool Verbose=False;

/*****************************************************************************/
void Projection_V0(fltarray & Event_of_V0, intarray & Catalogue_Event, 
		FewEventPoisson & FEP,User_Param & Param){
   //		Projection in V0 if input data is a fits

   
   int Nx = Event_of_V0.nx();
   int Ny = Event_of_V0.ny();
   int Nz = Event_of_V0.nz();
   
      
   
   char Name_Cube_Out2[256];
   if (Param.ProjV0){
	//FEP.building_V0_proj_3D_fltarr(Event_of_V0,Catalogue_Event);
	//Catalogue_Event=Event_of_V0
	for (int i=0;i<Nx;i++) 
   	for (int j=0;j<Ny;j++)
   	for (int k=0;k<Nz;k++)
   		Catalogue_Event(i,j,k) = (int) (Event_of_V0(i,j,k)+.5);

	//Event_of_V0=0.
   	Event_of_V0.init();
   	float tab [3]={0.1666666666,0.6666666666666,0.16666666666};
	// for the number of events in Event_of_V0 at pos i,j,k
      	// projection in V0 (scalar product with function phi)
   	for (int i=0;i<Nx;i++) 
   	for (int j=0;j<Ny;j++)
   	for (int k=0;k<Nz;k++) 	
	  for (int x_cur=-1;x_cur<=1;x_cur++)       
	  for (int y_cur=-1;y_cur<=1;y_cur++) 
	  for (int z_cur=-1;z_cur<=1;z_cur++) 
		if ((i+x_cur >= 0)    &&  (i+x_cur < Nx) 
			&& (j+y_cur >= 0) && (j+y_cur < Ny) 
			&& (k+z_cur >= 0) && (k+z_cur < Nz)) 
		   Event_of_V0(i+x_cur,j+y_cur,k+z_cur) += tab[x_cur+1]*tab[y_cur+1]
				*tab[z_cur+1]*Catalogue_Event(i,j,k);
         		
	if (Param.WriteRes){
        	sprintf(Name_Cube_Out2,"%s_V0.fits",Param.Name_Cube_Out);
		fits_write_fltarr(Name_Cube_Out2, (Event_of_V0));
	}
    }else
    	for(int i=0;i<Event_of_V0.nx();i++)
	for(int j=0;j<Event_of_V0.ny();j++)
	for(int k=0;k<Event_of_V0.nz();k++)
    		Catalogue_Event(i,j,k)=(int)(Event_of_V0(i,j,k) + 0.5);
}

/*********************************************************************/
//		Projection in V0 if input data is a catalogue
void Projection_V0(ArrayPoint & Catalogue, fltarray & ProjV0, 
		intarray & EventCatalogue, User_Param & Param){
	//initialisation of Event_of_V0 and Catalogue_Event
 
   float BinCat = Param.bin_cat;

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
   
 
   
    int Npx=(int)((Catalogue.Pmax.x()-Catalogue.Pmin.x())/BinCat+0.5);
    int Npy=(int)((Catalogue.Pmax.y()-Catalogue.Pmin.y())/BinCat+0.5);
    int Npz=(int)((Catalogue.Pmax.z()-Catalogue.Pmin.z())/BinCat+0.5);
 
   EventCatalogue.alloc( xSize, ySize, zSize );
   ProjV0.alloc( xSize, ySize, zSize );
   
 	if(Param.ProjV0)
	{
		//projection in V0
		// B(x) =1/12 (|x-2|^3 - 4 |x-1|^3 + 6 |x|^3 - 4 |x+1|^3 + |x+2| ^ 3)
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
      
	    float xRealVal  = Catalogue(i).x();
	    int   xPixelInd = int( ( xRealVal - xRealMin ) / BinCat );
	    float yRealVal  = Catalogue(i).y();
	    int   yPixelInd = int( ( yRealVal - yRealMin ) / BinCat );
	    float zRealVal  = Catalogue(i).z();
	    int   zPixelInd = int( ( zRealVal - zRealMin ) / BinCat );
	    
           // add event in event catalog
           EventCatalogue( xPixelInd, yPixelInd, zPixelInd )++;   
   
           
           // compute offset from center
	   //---------------------------
	   int xBeg=-1,xEnd=-1;
	   float xOffset = xRealVal - xCenter(xPixelInd);
	   if( fabs( xOffset ) <= 1e-06 ) {xBeg = -1; xEnd = 1;}
	   else if( xOffset > 0 )         {xBeg = -2; xEnd = 1;}
	   else if( xOffset < 0 )         {xBeg = -1; xEnd = 2;}
	      
	   int yBeg=-1,yEnd=-1;
	   float yOffset = yRealVal - yCenter(yPixelInd);
	   if( fabs( yOffset ) <= 1e-06 ) {yBeg = -1; yEnd = 1;}
	   else if( yOffset > 0 )         {yBeg = -2; yEnd = 1;}
	   else if( yOffset < 0 )         {yBeg = -1; yEnd = 2;}
 	   
	   int zBeg=-1,zEnd=-1;
	   float zOffset = zRealVal - zCenter(zPixelInd);
	   if( fabs( zOffset ) <= 1e-06 ) {zBeg = -1; zEnd = 1;}
	   else if( zOffset > 0 )         {zBeg = -2; zEnd = 1;}
	   else if( zOffset < 0 )         {zBeg = -1; zEnd = 2;}
	  
	   //bounds tests
	   while( xPixelInd+xBeg < 0 ) xBeg++;
	   while( xPixelInd+xEnd > xSize-1 ) xEnd--;
	   while( yPixelInd+yBeg < 0 ) yBeg++;
	   while( yPixelInd+yEnd > ySize-1 ) yEnd--;
	   while( zPixelInd+zBeg < 0 ) zBeg++;
	   while( zPixelInd+zEnd > zSize-1 ) zEnd--;
	   
	   // proj in V0
           for( int xInd=xBeg; xInd<=xEnd; ++xInd )
           for( int yInd=yBeg; yInd<=yEnd; ++yInd )
           for( int zInd=zBeg; zInd<=zEnd; ++zInd ) {
	      int xCur = int( SLINESIZE/2+1 + xInd*( SLINESIZE/4 ) 
	                      + xOffset*( SLINESIZE/4 ) );
	      int yCur = int( SLINESIZE/2+1 + yInd*( SLINESIZE/4 ) 
	                      + yOffset*( SLINESIZE/4 ) );
	      int zCur = int( SLINESIZE/2+1 + zInd*( SLINESIZE/4 ) 
	                      + zOffset*( SLINESIZE/4 ) );
	      
	      double xspline = bspline(xCur);
	      double yspline = bspline(yCur);
	      double zspline = bspline(zCur);
	      
              ProjV0( xPixelInd+xInd, yPixelInd+yInd, zPixelInd+zInd ) 
	         += float( xspline*yspline*zspline ); 
	    } 
         }   
	
   
	}
	else{		
		for(int i=0;i<Catalogue.np();i++)
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
                ProjV0 (pos_x,pos_y,pos_z)++;
                EventCatalogue (pos_x,pos_y,pos_z)++;  
		}
	
	}
}

/*********************************************************************/
//	read input data and test if it is a catalogue or a fits
void read_data(fltarray & Event_of_V0, intarray & Catalogue_Event ,
		FewEventPoisson & FEP,int & Nx,int & Ny,int & Nz,
		User_Param & Param, fitsstruct *HD){
		
	// reading input data 
	//test if input data is a fltarray or a catalogue
	if(strstr(Param.Name_Cube_In,".fits") || strstr(Param.Name_Cube_In,".Fits")
			|| strstr(Param.Name_Cube_In,".FITS")){ 
		Param.FitsMod=True;
		io_3d_read_data(Param.Name_Cube_In, Event_of_V0, HD);
                Catalogue_Event.alloc(Event_of_V0.nx(),Event_of_V0.ny(),Event_of_V0.nz());
		Projection_V0(Event_of_V0,Catalogue_Event,FEP,Param);//Projection in V0
 	}else if (strstr(Param.Name_Cube_In,".cat") || strstr(Param.Name_Cube_In,".Cat")
			|| strstr(Param.Name_Cube_In,".CAT")){
		Param.FitsMod=False;
		ArrayPoint Catalogue;//initialisation of catalogue
    		Catalogue.read(Param.Name_Cube_In);
		if(Catalogue.dim()!=3){
			cout <<"The catalogue do not have 3 dimentions"<<endl;
			exit(0);
		}
		Projection_V0(Catalogue,Event_of_V0,Catalogue_Event,Param);
		HD->hd_fltarray(Event_of_V0);
	}
	
	//	give a good Nbr_Plan if necessary
	Nx = Event_of_V0.nx();
    	Ny = Event_of_V0.ny();
    	Nz = Event_of_V0.nz();
	int min=MIN3(Nx, Ny, Nz);
	if ( min<9){
		cout << "All dimentions of the input cube must be >9"<<endl;
		exit(0);
	}
        if ( (min < 32) && (Param.Nbr_Plan>2 || !Param.Nbr_Plan)) 
		Param.Nbr_Plan=2;
      	else if ( min < 64 && (Param.Nbr_Plan>3 || !Param.Nbr_Plan)) 
		Param.Nbr_Plan=3;
      	else if ( min < 128 && (Param.Nbr_Plan>4 || !Param.Nbr_Plan)) 
		Param.Nbr_Plan=4;     
      	else if ( min < 256 && (Param.Nbr_Plan>5 || !Param.Nbr_Plan)) 
		Param.Nbr_Plan=5; 
      	else if ( min < 512 && (Param.Nbr_Plan>6 || !Param.Nbr_Plan)) 
		Param.Nbr_Plan=6;
      	else if ( min < 1024 && (Param.Nbr_Plan>7 || !Param.Nbr_Plan)) 
		Param.Nbr_Plan=7;
      	else if ( min < 2048 && (Param.Nbr_Plan>8 || !Param.Nbr_Plan)) 
		Param.Nbr_Plan=8;
      	else if ( min < 2*2048 && (Param.Nbr_Plan>9 || !Param.Nbr_Plan)) 
		Param.Nbr_Plan=9;
      	else if ( min < 2*4096 && (Param.Nbr_Plan>10 || !Param.Nbr_Plan)) 
		Param.Nbr_Plan=10;
		
}


/*********************************************************************/

void init_Nb_Event(intarray & Catalogue_Event,intarray *& Nb_Event, 
		fltarray & Mean_Nb_Ev, FewEventPoisson & FEP,
		ATROUS_3D_FIL & WT_tools, User_Param & Param)
{
	//	compute the cubes Nb_Event
	int Nx = Catalogue_Event.nx();
    	int Ny = Catalogue_Event.ny();
    	int Nz = Catalogue_Event.nz();
	char Name_Cube_Out2[50];
	if(Param.Verbose) cout << "Initalisation of Nb_Event"<<endl;
    	for(int s=-2;s<Param.Nbr_Plan-1;s++)
	{
		if (Param.Approx)
		{//we do not compute nb of events in boxes of 
			// (2^n+1)^3 but boxes of (2^n)^3
			if (s==-2)  FEP.event_scale_3D_Approx( Catalogue_Event, 
				Nb_Event[0],s, Param.Border, Param.Verbose );
			else if(s==-1)  FEP.event_scale_3D_Approx( Nb_Event[0], 
				Nb_Event[1],s, Param.Border, Param.Verbose );
			else if(s==0)  FEP.event_scale_3D_Approx( Nb_Event[1], 
				Nb_Event[0],s, Param.Border, Param.Verbose );	
    			else FEP.event_scale_3D_Approx( Nb_Event[s-1], Nb_Event[s], 
				s, Param.Border, Param.Verbose );  
		}
		else //we compute nb of events in boxes of (2^n+1)^3  
			if(s>=0)  FEP.event_scale_3D( Catalogue_Event, Nb_Event, 
					s, Param.Border, Param.Verbose); 

    	        // if (s >=0 ) Nb_Event[s].info("FIRST INFO EVENT");
	}
	if(Param.Approx)
	{
	   intarray Buff(Nx,Ny,Nz);
	   for(int s=0;s<Param.Nbr_Plan-1;s++)
	   {
	   	//for a better appoximation of counting boxes 
		//it computes the equivalant of a box "9*9*9"
		//so it multiply by 9*9*9/(8*8*8)
		//and it makes an average with the 7 neighbors of right-hand side
	   	float mult=pow((double)(1+ pow((double)2.,-s-3.)),3.)/8.;
		Buff=Nb_Event[s];
		for(int i=0;i<Nx-1;i++)
		for(int j=0;j<Ny-1;j++)
		for(int k=0;k<Nz-1;k++)
			Nb_Event[s](i,j,k)=(int)((Buff(i,j,k)+Buff(i,j,k+1)+
			   Buff(i,j+1,k)+Buff(i,j+1,k+1)+Buff(i+1,j,k)+
			   Buff(i+1,j,k+1)+Buff(i+1,j+1,k)+Buff(i+1,j+1,k+1))
			   *mult +0.5);
		
	        mult=pow((double)(1+ pow((double)2.,-s-3.)),3.)/4.;
		
		
		//for the plan z=z_max-1
	   	for(int i=0;i<Nx-1;i++)
	   	for(int j=0;j<Ny-1;j++)
	   		Nb_Event[s](i,j,Nz-1)=(int)((Buff(i,j,Nz-1)
				+Buff(i,j+1,Nz-1)+Buff(i+1,j,Nz-1)
				+Buff(i+1,j+1,Nz-1))*mult +0.5);
				
		//for the plan y=y_max-1		
		for(int i=0;i<Nx-1;i++)
		for(int k=0;k<Nz-1;k++)
			Nb_Event[s](i,Ny-1,k)=(int)((Buff(i,Ny-1,k)+Buff(i,Ny-1,k+1)
			   +Buff(i+1,Ny-1,k)+Buff(i+1,Ny-1,k+1))*mult +0.5);
			   
		//for the plan x=x_max-1	   
		for(int j=0;j<Ny-1;j++)
		for(int k=0;k<Nz-1;k++)
			Nb_Event[s](Nx-1,j,k)=(int)((Buff(Nx-1,j,k)+Buff(Nx-1,j,k+1)
			   +Buff(Nx-1,j+1,k)+Buff(Nx-1,j+1,k+1))*mult +0.5);
	   }
 	} 
	// for(int s=0;s<Param.Nbr_Plan-1;s++) Nb_Event[s].info("INFO EVENT");
	for(int s=0;s<Param.Nbr_Plan-1;s++)
	{
		if(Param.WriteRes)
		{
			sprintf(Name_Cube_Out2,"%s_nb_ev_%d.fits",
				Param.Name_Cube_Out,s+1);
			fits_write_intarr(Name_Cube_Out2, 
				Nb_Event[s]);
		}
		Mean_Nb_Ev(s)=Nb_Event[s].mean();
	}
	if(Param.WriteRes||Param.WriteSupport){
		sprintf(Name_Cube_Out2,"%s_nb_ev_mean",
				Param.Name_Cube_Out);
		fits_write_fltarr(Name_Cube_Out2, Mean_Nb_Ev);
	}
}

/*********************************************************************/

void init_support(intarray *& Nb_Event,intarray *& Support, 
		fltarray & Mean_Nb_Ev,  ATROUS_3D_FIL & WT_tools, 
		FewEventPoisson & FEP,int Nx, int Ny, int Nz, 
		User_Param & Param){
		//	read or compute the support 
    char Name_Cube_Out2[50];
    if(Param.ReadSupport)
    {
    	for(int s=0;s<Param.Nbr_Plan-1 ;s++)
	{
		sprintf(Name_Cube_Out2,"%s_support_%d",Param.Name_Cube_Out,s+1);
		fits_read_intarr(Name_Cube_Out2, Support[s]);
		if((Support[s].nx()!=Nx) || (Support[s].ny()!=Ny)
				|| (Support[s].nz()!=Nz)){
			cout <<"Error, the support read do not corespond to "
				<<"the input cube"<<endl;
			exit(0);
		}else 
			if(Param.Verbose) cout << "The support of scale "
				<<s<<" is read."<<endl; 
		    
	}
	sprintf(Name_Cube_Out2,"%s_nb_ev_mean",Param.Name_Cube_Out);
	fits_read_fltarr(Name_Cube_Out2, Mean_Nb_Ev);
    }
    else
    {
    	for(int s=0;s<Param.Nbr_Plan-1 ;s++)
	{
           FEP.event_set_support(WT_tools.Wavelet_Coef, s, 
			Param.FirstScaleDetect, Param.NEGFirstScaleDetect, Param.MinEvent, 
			Param.OnlyPosVal, I_MIRROR, Nb_Event, Support, 
			Param.WriteRes);
           if(Param.Verbose) cout << "Support of scale number " 
			<< s+1 << " completed"<< endl; 
		if (Param.WriteRes||Param.WriteSupport){
    			sprintf(Name_Cube_Out2,"%s_support_%d",
				Param.Name_Cube_Out,s+1);
    			fits_write_intarr(Name_Cube_Out2, Support[s]);
    		}
    	}
    }
    //filtring of single points in the support 
    if(!Param.SinglePoint) {
    	WT_tools.No_Single_Point(Nx, Ny,Nz, Support, FEP, Param.Nbr_Plan-1);
    	if(Param.Verbose) cout << "Single points have been cleared"<<endl;
    }

  
}


/*********************************************************************/ 
void write_result_inter(fltarray & Output_Cube,int s, User_Param & Param,char text[])
{
    char Name_Cube_Out2[50];
    if(Param.Verbose) cout << "Cube is reconstructed"<< endl;
    sprintf(Name_Cube_Out2,"%s_%s_%d",Param.Name_Cube_Out,text,s);
    fits_write_fltarr(Name_Cube_Out2, Output_Cube);
}
 
 
/*********************************************************************/ 
void write_result(fltarray & Output_Cube, User_Param & Param, fitsstruct *HD){
    char Name_Cube_Out2[50];
    if(Param.Verbose) cout << "Cube is reconstructed"<< endl;
    sprintf(Name_Cube_Out2,"%s",Param.Name_Cube_Out);
    fits_write_fltarr(Name_Cube_Out2, Output_Cube, HD);
}   


 
/*********************************************************************/
//			Class User_Param
/*********************************************************************/
void User_Param::put_default_val(){
	Nbr_Plan=0;  
	FirstScaleDetect=0;
	NEGFirstScaleDetect=-1;
	MinEvent=5;
	Rythm=0;
	NAutoConv=30;
	NbIter=10;
	Border=I_ZERO;
	epsilon=1e-3;
	bin_cat=1.;
	N_Sigma=0.5;
	Verbose = False;
	ProjV0 = True;
	WriteRes = False;
	InfoBand = False; 
	Read_Histo = False;
	WriteHisto = False;
	SinglePoint = True;
	OnlyPosVal = False;
	PositiveImage = True;
	WriteOpt=False;
	AddLastScale=True;
	ReadSupport=False;
	WriteSupport=False;
	Approx=True;
	Adjoint=True;
	RemoveBorder=False;
	CFAIter=True;
	WriteDecomp=False;
}

void User_Param::usage(char *argv[]){
    fprintf(OUTMAN, "Usage: %s options in_cube out_cube\n\n", argv[0]);
    fprintf(OUTMAN, "   where options =  \n");

    fprintf(OUTMAN, "         [-n number_of_scales]\n");
    fprintf(OUTMAN, "             Number of scales used in the multiresolution transform.\n");
    fprintf(OUTMAN, "             Default number is automatically calculated.\n");

    fprintf(OUTMAN, "         [-F first_detection_scale]\n");
    fprintf(OUTMAN, "              First scale used for the detection.\n");
    fprintf(OUTMAN, "              Default is 1.\n");
    
    fprintf(OUTMAN, "         [-N first_detection_scale_neg_coeff]\n");
    fprintf(OUTMAN, "              First scale used for the detection of negative coefficients.\n");
    fprintf(OUTMAN, "              Default is first_detection_scale.\n");
    
    fprintf(OUTMAN, "         [-E epsilon]\n");
    fprintf(OUTMAN, "              Epsilon = precision for computing thresholds\n");
    fprintf(OUTMAN, "              Default is 1.00e-03.\n");
        
    fprintf(OUTMAN, "         [-m min_event]\n");
    fprintf(OUTMAN, "              Minimum number of event to detect ( default=5 ).\n");
      
//     fprintf(OUTMAN, "         [-b type_border = 1 to 4]\n");
//     fprintf(OUTMAN, "              give the type of border {1:I_CONT|2:I_MIRROR|3:I_ZERO|4:I_PERIOD}\n");
//     fprintf(OUTMAN, "              ( default=2:I_MIRROR ).\n");
  
    fprintf(OUTMAN, "         [-i number_iteration]\n");
    fprintf(OUTMAN, "              Number of Iteration for the reconstruction (default=10 ).\n");
    
    fprintf(OUTMAN, "         [-B bin] (only for catalogue)\n");
    fprintf(OUTMAN, "              If the input data is a catalogue, bin gives resolution ( default=1 ).\n");
    
    fprintf(OUTMAN, "         [-R]\n");
    fprintf(OUTMAN, "              Remove the structures near the borders.\n");
    
//     fprintf(OUTMAN, "         [-A]\n");
//     fprintf(OUTMAN, "              do not use adjoint operator for the reconstruction.\n");
//     

    fprintf(OUTMAN, "         [-C]\n");
    fprintf(OUTMAN, "              Remove the Smoothness constraint. \n");
    
    fprintf(OUTMAN, "         [-p]\n");
    fprintf(OUTMAN, "              Detect only positive structure. Default is no.\n");
  
    fprintf(OUTMAN, "         [-a]\n");
    fprintf(OUTMAN, "              No approximation for calculating the number of events in boxes.\n");
    fprintf(OUTMAN, "              By default, an approximation is made which speed up the algorithm.\n");
// 	      
    fprintf(OUTMAN, "         [-w]\n");
    fprintf(OUTMAN, "              Write the multiresolution support to the disk.\n");
        
//    fprintf(OUTMAN, "         [-S]\n");
//    fprintf(OUTMAN, "              try to read the support already computed.\n");
    
    fprintf(OUTMAN, "         [-k]\n");
    fprintf(OUTMAN, "              Suppress isolated pixels in the support. Default is no.\n");

    fprintf(OUTMAN, "         [-K]\n");
    fprintf(OUTMAN, "              Suppress the last scale. Default is no.\n");
    
    fprintf(OUTMAN, "         [-H]\n");
    fprintf(OUTMAN, "              Read pre-computed histogram table (produced by mr3d_histo). \n");
    fprintf(OUTMAN, "               Default is no.\n");
    verbose_usage(); 
        
    manline();
    vm_usage();
    manline();    
    exit(-1);
}

/*********************************************************************/

/* GET COMMAND LINE ARGUMENTS */
User_Param::User_Param(int argc, char *argv[]){
    int c;

#ifdef LARGE_BUFF
    int VMSSize=-1;
    Bool OptZ = False;
    char VMSName[1024] = "";
#endif  
    put_default_val();//put default values for each parametres 
    /* get options */
    while ((c = GetOpt(argc,argv,"HN:s:pwCWASRakKi:B:m:n:F:E:vzZ:")) != -1) 
    {
	switch (c) 
        {
	  case 'H': Read_Histo = True; break;
          case 's':
	       {
	        float N_Sigma;
		/* -s <nsigma> */
		if (sscanf(OptArg,"%f",&N_Sigma) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad N_Sigma: %s\n", OptArg);
		    exit(-1);
		}
                if (N_Sigma <= 0.)   
		{
		    fprintf(OUTMAN, "Error: bad N_Sigma: %s\n", OptArg);
		    exit(-1);
		}
		epsilon = (1. - erf((double) N_Sigma / sqrt((double) 2.)));
		}
		break;
	   case 'v': Verbose = True; WriteOpt=True; break;
 	   case 'k': SinglePoint=False; break;
 	   case 'p': OnlyPosVal=True;break;
	   case 'S': ReadSupport=True;break;
 	   case 'K': AddLastScale=False; break;
	   case 'W': WriteRes=True; break;
	   case 'a': Approx=False; break;
	   case 'A': Adjoint=False; break;
 	   case 'C': CFAIter=False;break;
	   case 'w': WriteSupport=True; break;
	   case 'R': RemoveBorder=True; break;
	   case 'n':
		/* -n <Nbr_Plan> */
		if (sscanf(OptArg,"%d",&Nbr_Plan) != 1) 
                {
		    fprintf(OUTMAN, "bad number of scales: %s\n", OptArg);
		    exit(-1);
		}
                if ((Nbr_Plan <= 1)) //|| (Nbr_Plan > MAX_SCALE_1D)) 
                {
		    fprintf(OUTMAN, "bad number of scales: %s\n", OptArg);
 		    fprintf(OUTMAN, "nbr Scales > 1");
		    exit(-1);
		}
		break;
	   case 'F':
		/* -s <First_Plan_Detect> */
		if ((sscanf(OptArg,"%d",&FirstScaleDetect) != 1) || (FirstScaleDetect <= 0)){
		    fprintf(OUTMAN, "bad first scale number: %s\n", OptArg);
		    exit(-1);
		}
		FirstScaleDetect -=1;
		break;
	  case 'N':
		/* -s <First_Plan_Detect> */
		if ((sscanf(OptArg,"%d",&NEGFirstScaleDetect) != 1) || (NEGFirstScaleDetect <= 0)){
		    fprintf(OUTMAN, "bad first scale number: %s\n", OptArg);
		    exit(-1);
		}
		NEGFirstScaleDetect  -=1;
		break;
	  case 'm':
		/* -m <Min_Event> */
		if ((sscanf(OptArg,"%d",&MinEvent) != 1) || (MinEvent < 0)){
		    fprintf(OUTMAN, "bad first scale detected: %s\n", OptArg);
		    exit(-1);
		}
		break;	
	   case 'b':
		/* -b <border> */
		int i;
		if ((sscanf(OptArg,"%d",&i) != 1)||(i<=0)||(i>4)){
		    fprintf(OUTMAN, "bad number for border: %s\n", OptArg);
		    exit(-1);
		}
		switch(i){
			case 1:Border=I_CONT;break;
			case 2:Border=I_MIRROR;break;
			case 3:Border=I_ZERO;break;
			case 4:Border=I_PERIOD;break;
		}
		break;	
	   case 'E':
		/* -e <epsilon> */
		if ((sscanf(OptArg,"%f",&epsilon) != 1)|| (epsilon < 0)){
		    fprintf(OUTMAN, "bad value of epsilon: %s\n", OptArg);
		    exit(-1);
		}
		break;
 	   case 'B':
		/* -B <bin> */
		if ((sscanf(OptArg,"%f",&bin_cat) != 1)|| (bin_cat < 0)){
		    fprintf(OUTMAN, "bad value of bin: %s\n", OptArg);
		    exit(-1);
		}
		break;	
	   case 'i':
		/* -i <NbIter> */
		if ((sscanf(OptArg,"%d",&NbIter) != 1)|| (NbIter < 0)){
		    fprintf(OUTMAN, "bad number of iteration: %s\n", OptArg);
		    exit(-1);
		}
		break;	
#ifdef LARGE_BUFF
	    case 'z':
	        if (OptZ == True)
		{
                   fprintf(OUTMAN, "Error: Z option already set...\n");
                   exit(-1);
                }
	        OptZ = True;
	        break;
            case 'Z':
	        if (sscanf(OptArg,"%d:%s",&VMSSize, VMSName) < 1)
		{
		   fprintf(OUTMAN, "Error: syntaxe is Size:Directory ... \n");
		   exit(-1);
		}
	        if (OptZ == True)
		{
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

       /* get optional input file names from trailing 
          parameters and open files */
	if (OptInd < argc) strcpy(Name_Cube_In, argv[OptInd++]);
        else usage(argv);

	if (OptInd < argc) strcpy(Name_Cube_Out, argv[OptInd++]);
        else usage(argv);

	/* make sure there are not too many parameters */
	if (OptInd < argc)
        {
		fprintf(OUTMAN, "Too many parameters: %s ...\n", argv[OptInd]);
		usage(argv);
	}

	
#ifdef LARGE_BUFF
    if (OptZ == True) vms_init(VMSSize, VMSName, Verbose);
#endif

}

/*********************************************************************/
void User_Param::WriteOptProc(){
	cout << "####################################################"
		<< endl;
	cout << "STATE OF OPTIONS: "<< endl;
	cout << "Number of Scales:		"<<Nbr_Plan<<endl;
	cout << "First Detection Scale:	"<<FirstScaleDetect+1<<endl;
	if (NEGFirstScaleDetect > 0) 
	cout << "First Detection Scale for Negative Coefficients:	"<<NEGFirstScaleDetect+1<<endl;
	cout << "Epsilon:			"<<epsilon<<endl;
	cout << "Min_event:			"<<MinEvent<<endl;
	if(!Read_Histo){ 
		cout << "rythm:				"<<Rythm<<endl;
		cout << "number_of_convolution:		"<<NAutoConv<<endl;
	}
	cout << "Number of Iteration:		"<<NbIter<<endl;
	if(!FitsMod)
		cout << "bin:				"<<bin_cat<<endl;
	if(CFAIter) cout << "N_Sigma:			"<<N_Sigma<<endl;
	// if(Adjoint) cout << "adjoint operator is used for reconstruction"<<endl;
	// else cout << "adjoint operator is not used for reconstruction"<<endl;
	if(CFAIter) cout << "Use the CFA algorithm for the iterations"<<endl;
	if(OnlyPosVal) cout << "Detect only in positive wavelet coef" << endl;
	else cout << "detect in positive and negative wavelet coef" << endl;
	if(RemoveBorder) cout << "Remove structures near the borders."<<endl;
	else cout << "Do not remove structures near the borders."<<endl;
	if(InfoBand){ 
		cout << "Print statistical information about each band";
		cout <<endl;}
	if(Verbose) cout << "Verbose mode"<<endl;
	if(ProjV0) cout <<"Projection in V0 space before wavelet transform";
	else cout <<"No projection in V0 space";
	cout <<endl;
	if(WriteRes) cout << "Write all results in files"<<endl;
 	if(WriteDecomp) cout <<"Write the wavelet's coefficients of the result"<< endl;
	cout << "The border type is ";
	switch(Border){
		case I_CONT: cout <<"continuous";break;
		case I_MIRROR: cout <<"mirror";break;
		case I_ZERO: cout <<"zero";break;
		case I_PERIOD: cout <<"periodic";break;
	}
	cout <<endl;
	if(WriteHisto) cout<<"Write histograms"<<endl;
	if(Read_Histo) cout << "Try to read the histograms already computed"; 
	else cout<< "Compute histograms";
	cout << endl;
	if(SinglePoint) cout<< "filter the single points"<<endl; 
	if (!AddLastScale) cout <<"Remove the last scale"<<endl;
 	if(PositiveImage) cout << "Positivity Constraint"<<endl;
 	if(!Approx) cout << "Number of event in boxes is computed exactely"<<endl;
	else cout << "Number of event in boxes is approximated"<<endl;
	if(ReadSupport) cout << "Try to read the support already computed"<<endl;
	if(WriteSupport) cout << "write the support"<<endl;
	cout << "INPUT Filename = " << Name_Cube_In  << endl;
        cout << "OUTPUT Filename = " << Name_Cube_Out  << endl;
	cout << "####################################################"
		<< endl;
	cout << endl;
}


/*********************************************************************/ 
/*			Class Iteration	 CFA			     */
/*********************************************************************/ 

//	initialisation of the class Iteration with CFA algorithm
Iteration_CFA::Iteration_CFA(fltarray & Event_of_V0,fltarray & Output_Cube, 
			intarray *& Support,intarray *& Nb_Event_In_Boxes, 
			fltarray & Mean_Nb_Ev, ATROUS_3D_FIL  & WT_tools,FewEventPoisson & FEP, 
			User_Param & Param)
{
    _Nx=Event_of_V0.nx();
    _Ny=Event_of_V0.ny();
    _Nz=Event_of_V0.nz();		
    _NbIter=Param.NbIter;
    _Param= &Param;
    _Output_Cube= &Output_Cube;
    _Mean_Nb_Ev=&Mean_Nb_Ev;
    _Event_of_V0= &Event_of_V0;
    _Nb_Event_In_Boxes=Nb_Event_In_Boxes;
    _WT_tools= &WT_tools;
    _Support=Support;
    //Allocation of first Wavelet Coefficient  
    WT_tools.alloc(_Wavelet_Coef_Init,_Nx,_Ny,_Nz,_Param->Nbr_Plan);
    //Allocation of cube for the results 
    _FEP=&FEP;
    _tab_sigma=fltarray(_Param->Nbr_Plan-1,2);
    
}

void Iteration_CFA::init(){ //	iter()
//initialisation of the iterator CFA
	

    _CurrentIter=1;
    _Lmax=1;
    _delta=(float)_Lmax/_NbIter;


    // wavelet 3D transformation
    // wavelet 3D transformation
    if(_NbIter==1)
    	_WT_tools->wttransform(*_Event_of_V0, _Param->Nbr_Plan,_Param->WriteDecomp,
    		_Param->InfoBand, _Param->Name_Cube_Out);
    else 
    	_WT_tools->wttransform(*_Event_of_V0, _Param->Nbr_Plan,False,
    		_Param->InfoBand, _Param->Name_Cube_Out);
    
    _lambda=_Lmax;

    //creation of the support
    init_support(_Nb_Event_In_Boxes,_Support, *_Mean_Nb_Ev,*_WT_tools,*_FEP,_Nx,_Ny,_Nz,*_Param);
    
    _FEP->get_sigma_band(*_Mean_Nb_Ev,_tab_sigma,_Param->N_Sigma,_Param->Nbr_Plan);
    // _tab_sigma.display();
    // if(_Param->ReadSupport == False) free( (intarray *) _Nb_Event_In_Boxes, (int) _Param->Nbr_Plan);
    
    _WT_tools->Filtering_Structure(_Support, _Param->Nbr_Plan,_Param->InfoBand,_Param->RemoveBorder,
    		                   _Param->FirstScaleDetect);
    //_WT_tools->threshold(_Support,_Param->Nbr_Plan);
    (*_Output_Cube).init();
    
    _WT_tools->ModifiedAWT = True;
    if (_WT_tools->ModifiedAWT == True)
       _WT_tools->wttransform(*_Event_of_V0, _Param->Nbr_Plan,False,
    		              _Param->InfoBand, _Param->Name_Cube_Out);
    
    for(int s=0; s<_Param->Nbr_Plan;s++)
    {
    	_Wavelet_Coef_Init[s]=(_WT_tools->Wavelet_Coef)[s];
        // cout << "INITBAND: " << s << " "  << _Wavelet_Coef_Init[s].min() << " " << _Wavelet_Coef_Init[s].max()<< endl;
    }
}


void Iteration_CFA::operator++(){  //	++iter
	(*this)++;
}

void Iteration_CFA::operator++(int x)
{
   //	iter++
   // iteration
	
  if(_Param->Verbose) cout << endl << "Computing "<< _CurrentIter  <<" iteration"<<endl;
  // if (_WT_tools->ModifiedAWT == True) cout << "Mod. AWT " << endl;
  // if (_Param->Adjoint== True) cout << "Use Adjoint " << endl;
  // else cout << "no Adjoint " << endl;
  
  _WT_tools->wttransform(*_Output_Cube, _Param->Nbr_Plan, False, False, 
		        _Param->Name_Cube_Out);

  //_WT_tools->threshold(_Support,_Param->Nbr_Plan);
  for(int s=0;s<_Param->Nbr_Plan-1;s++)
  {
     // _local_threshold=_tab_sigma(s,1);
     float SoftThresholdLevel; //  = _lambda*_tab_sigma(s,0);
     //cout << "Scale " << s+1 << endl;
     //cout << " _Nb_Event_In_Boxes " << _Nb_Event_In_Boxes[s].nx() <<  " " << _Nb_Event_In_Boxes[s].ny() <<  " " << _Nb_Event_In_Boxes[s].nz() << endl;
     
     for(int i=0;i<_Nx;i++)
     for(int j=0;j<_Ny;j++)
     for(int k=0;k<_Nz;k++)
     {
         float SeuilMin, SeuilMax;
         float coef_init=_Wavelet_Coef_Init[s](i,j,k);
	 int Nev = _Nb_Event_In_Boxes[s](i,j,k);
	 _FEP->get_threshold(Nev, s, SeuilMin, SeuilMax);
         float DetectLevel = (coef_init > 0) ? SeuilMax: ABS(SeuilMin);
	 float ErrTolerance = DetectLevel  / 5;
         float coef=_WT_tools->Wavelet_Coef[s](i,j,k);
         // if ((_Support[s](i,j,k)==1) && (fabs(coef_init-coef)> ErrTolerance)) coef=coef_init;
	 if (_Support[s](i,j,k)==1)  coef=coef_init;
	 else
	 {
	     if ((coef_init > 0) && (coef_init < SeuilMax) && (coef > SeuilMax)) coef= SeuilMax;
	     else if ((coef_init < 0) && (coef_init > SeuilMin) && (coef < SeuilMin)) coef= SeuilMin;
	 }
	 SoftThresholdLevel = _lambda * DetectLevel / 2;
         coef=soft_threshold(coef,SoftThresholdLevel);
         _WT_tools->Wavelet_Coef[s](i,j,k)=coef;
     }
  }
	
  _WT_tools->Wavelet_Coef[_Param->Nbr_Plan-1]=
			_Wavelet_Coef_Init[_Param->Nbr_Plan-1];
	

  if(_CurrentIter >= _NbIter-1)
  {
     _WT_tools->wtrecons(*_Output_Cube, _Param->Nbr_Plan, False, 
			_Param->AddLastScale, _Param->Adjoint);
     if(_Param->WriteDecomp)
        for(int s=_Param->FirstScaleDetect;s<_Param->Nbr_Plan;s++)
				write_result_inter(_WT_tools->Wavelet_Coef[s]
					,s,(*_Param),"Band");
   }
   else  _WT_tools->wtrecons(*_Output_Cube, _Param->Nbr_Plan, False, True, _Param->Adjoint);
			
  if(_Param->PositiveImage) (*_Output_Cube).inf_threshold(0.);

// 		for(int i=0;i<_Nx;i++)
// 		for(int j=0;j<_Ny;j++)
// 		for(int k=0;k<_Nz;k++)
// 			if((*_Output_Cube)(i,j,k)<0) (*_Output_Cube)(i,j,k)=0;
// 	
  //if(_Param->InfoBand == True)  (*_Output_Cube).info("Result :");
  if(_Param->InfoBand == True)  (*_Output_Cube).info();

  if(_Param->WriteRes == True) write_result_inter((*_Output_Cube),_CurrentIter,(*_Param),"iter");
   _CurrentIter++;
   _lambda-=_delta;
}

Iteration_CFA::~Iteration_CFA(){
	free(_Wavelet_Coef_Init,_Param->Nbr_Plan);
}


/*********************************************************************/ 
/*			Class Iteration				     */
/*********************************************************************/ 
//	initialisation of the class Iteration
Iteration::Iteration(fltarray & Event_of_V0,fltarray & Output_Cube, 
			intarray *& Support,intarray *& Nb_Event_In_Boxes, 
			fltarray & Mean_Nb_Ev,ATROUS_3D_FIL & WT_tools,FewEventPoisson & FEP, 
			User_Param & Param){
    _Nx=Event_of_V0.nx();
    _Ny=Event_of_V0.ny();
    _Nz=Event_of_V0.nz();		
    _NbIter=Param.NbIter;
    _CurrentIter=1;
    _Param= &Param;
    _Output_Cube= &Output_Cube;
    _Event_of_V0= &Event_of_V0;
    _Nb_Event_In_Boxes=Nb_Event_In_Boxes;
    _WT_tools= &WT_tools;
    _Support=Support;
    //Allocation of cube for the residues 
    _Residue.alloc(_Nx,_Ny,_Nz);
    //Allocation of cube for the results 
    _Result.alloc(_Nx,_Ny,_Nz);
    _FEP=&FEP;
    _Mean_Nb_Ev=&Mean_Nb_Ev;;
}

void Iteration::init(){ //	iter()
//initialisation of the iterator
	

    // wavelet 3D transformation
    _WT_tools->wttransform(*_Event_of_V0, _Param->Nbr_Plan,_Param->WriteRes,
    	_Param->InfoBand, _Param->Name_Cube_Out);		
	
    //creation of the support
    init_support(_Nb_Event_In_Boxes,_Support,*_Mean_Nb_Ev,*_WT_tools,*_FEP,_Nx,_Ny,_Nz,*_Param);
    
    if(_Param->ReadSupport == False) free(_Nb_Event_In_Boxes,_Param->Nbr_Plan);
    
    _WT_tools->Filtering_Structure(_Support, _Param->Nbr_Plan,_Param->InfoBand,_Param->RemoveBorder,
    		_Param->FirstScaleDetect);
		
		
    _WT_tools->ModifiedAWT = True;
    if (_WT_tools->ModifiedAWT == True)
       _WT_tools->wttransform(*_Event_of_V0, _Param->Nbr_Plan,False,
    		              _Param->InfoBand, _Param->Name_Cube_Out);
    
    _WT_tools->threshold(_Support,_Param->Nbr_Plan);
	
    //reconstruction 
    if(_CurrentIter >= _NbIter-1)
    {
    	_WT_tools->wtrecons(*_Output_Cube ,_Param->Nbr_Plan,False, _Param->AddLastScale,_Param->Adjoint);
	if(_Param->WriteDecomp)
		for(int s=_Param->FirstScaleDetect;s<_Param->Nbr_Plan;s++)
		 write_result_inter(_WT_tools->Wavelet_Coef[s],s,(*_Param),"Band"); 
    }
    else 
	_WT_tools->wtrecons(*_Output_Cube ,_Param->Nbr_Plan,False, True,_Param->Adjoint); 
	
    if(_Param->PositiveImage) (*_Output_Cube).inf_threshold(0.);

// 		for(int i=0;i<_Nx;i++)
// 		for(int j=0;j<_Ny;j++)
// 		for(int k=0;k<_Nz;k++)
// 			if((*_Output_Cube)(i,j,k)<0) (*_Output_Cube)(i,j,k)=0;
    if(_Param->NbIter>0 && _Param->WriteRes) write_result_inter(*_Output_Cube,0,*_Param,"iter");


    for(int i=0;i<_Nx;i++)
    for(int j=0;j<_Ny;j++)
    for(int k=0;k<_Nz;k++)
		_Residue(i,j,k)=(*_Event_of_V0)(i,j,k) - (*_Output_Cube)(i,j,k);
    if(_Param->Verbose == True) cout<<endl<<"Computing 1 iteration"<<endl;
    if(_Param->InfoBand == True)
    {	
		cout <<"sigma residue    "<< _Residue.sigma()<<endl;
		cout <<"mean  residue    "<< _Residue.mean()<<endl;
		cout <<"sigma new output "
			<< (*_Output_Cube).sigma()<<endl;
		cout <<"mean new output  "
			<< (*_Output_Cube).mean()<<endl;
    }

}



void Iteration::operator++(){  //	++iter
	(*this)++;
}

void Iteration::operator++(int x){	//	iter++
	// iteration
	
    	_WT_tools->wttransform(_Residue, _Param->Nbr_Plan, False, False, 
		_Param->Name_Cube_Out);
	_WT_tools->threshold(_Support,_Param->Nbr_Plan);
	if(_CurrentIter >= _NbIter-1)
		_WT_tools->wtrecons(_Result, _Param->Nbr_Plan, False, _Param->AddLastScale, _Param->Adjoint);
	else 
		_WT_tools->wtrecons(_Result, _Param->Nbr_Plan, False, True, _Param->Adjoint);
	for(int i=0;i<_Nx;i++)
	for(int j=0;j<_Ny;j++)
	for(int k=0;k<_Nz;k++)
		(*_Output_Cube)(i,j,k)+=_Result(i,j,k);
		
	if(_Param->PositiveImage == True) (*_Output_Cube).inf_threshold(0.);
//      for(int i=0;i<_Nx;i++)
// 		for(int j=0;j<_Ny;j++)
// 		for(int k=0;k<_Nz;k++)
// 			if((*_Output_Cube)(i,j,k)<0) (*_Output_Cube)(i,j,k)=0;
		 	
	for(int i=0;i<_Nx;i++)
	for(int j=0;j<_Ny;j++)
	for(int k=0;k<_Nz;k++)
		_Residue(i,j,k)=(*_Event_of_V0)(i,j,k) - (*_Output_Cube)(i,j,k);
	if(_Param->Verbose) cout<<endl<<"Computing "<<_CurrentIter+1 
		<<" iteration"<<endl;

        // if (_WT_tools->ModifiedAWT == True) cout << "Mod. AWT " << endl;
	// if (_Param->Adjoint== True) cout << "Use Adjoint " << endl;
        // else cout << "no Adjoint " << endl;
	if(_Param->InfoBand){	
		cout <<"sigma residue    "<< _Residue.sigma()<<endl;
		cout <<"mean  residue    "<< _Residue.mean()<<endl;
		cout <<"sigma new output "
			<< (*_Output_Cube).sigma()<<endl;
		cout <<"mean new output  "
			<< (*_Output_Cube).mean()<<endl;
	}
        if(_Param->WriteRes) write_result_inter((*_Output_Cube),_CurrentIter,(*_Param),"iter");
	_CurrentIter++;
}

Iteration::~Iteration(){
	_Result.free();
	_Residue.free();
}
 /********************************************************************/   
 /***************************P P**************************************/    
 /********************************************************************/    

int main(int argc, char *argv[])
{
    lm_check(LIC_MR3);
    User_Param Param(argc, argv);// Get command line arguments
    int k,Nx,Ny,Nz;        	// dimention of input signal 
    intarray Catalogue_Event;	//input catalogue
    fltarray Event_of_V0;	//same as Catalogue_Event but projected 
    				//in V0 (if necessary)
    fitsstruct HD;
    char Cmd[512];
    extern softinfo Soft;
   	    
    Soft.mr3();
    Cmd[0] = '\0';
    for (k =0; k < argc; k++) sprintf(Cmd, "%s %s", Cmd, argv[k]);


    //creating Atrous_3D class 
    ATROUS_3D_FIL WT_tools; 
    
    //creating FewEventPoisson class
    FewEventPoisson FEP(Param.Read_Histo,Param.NAutoConv,Param.epsilon);
      
    // reading input data 
    read_data(Event_of_V0,Catalogue_Event,FEP,Nx,Ny,Nz,Param, &HD);
    if(Param.WriteOpt) Param.WriteOptProc();
    if (Param.Verbose == True)
    {
       cout << endl << endl << "PARAMETERS: " << endl << endl;
       cout << "File Name in = " << Param.Name_Cube_In << endl;
       cout << "File Name out = " <<  Param.Name_Cube_Out << endl; 
       cout <<"Dimentions of the cube are: "<<Nx<<"  "<<Ny<<"  "<<Nz<<"  "<<endl;
       cout <<"The Number of Scale is: "<<Param.Nbr_Plan<<endl;
       if (Param.SinglePoint == False) cout << "Kill isolated pixels" << endl;
       else cout << "Keep isolated pixels" << endl;
       cout <<"Epsilon = " << Param.epsilon << endl;
       if (Param.OnlyPosVal == False) cout << "Detect both positive and negative coeff." << endl;
       else cout << "Detect only positive coeff." << endl;
       if (Param.RemoveBorder == True) cout << "Remove structures on the border " << endl;
       else cout << "Keep structures on the border " << endl;
       cout << endl << endl;       
    }
    //if(Param.InfoBand) Event_of_V0.info("Input cube after the projection");
    if(Param.InfoBand) Event_of_V0.info();
  
    fltarray Mean_Nb_Ev;
    intarray *Nb_Event_In_Boxes; //cubes of number of event per box
    if (Param.ReadSupport == False)
    {
    	alloc_array(Nb_Event_In_Boxes,Nx,Ny,Nz,Param.Nbr_Plan);
	Mean_Nb_Ev.alloc(Param.Nbr_Plan);
    	//initalisation of Nb_Event
    	init_Nb_Event(Catalogue_Event,Nb_Event_In_Boxes,Mean_Nb_Ev,FEP,WT_tools,Param);
    }
    
    intarray *Support;//Allocation of the support for each plan
    alloc_array(Support,Nx,Ny,Nz,Param.Nbr_Plan);
 
    //FEP.find_threshold(Param.WriteRes, Param.WriteHisto, 
    //	              Param.Verbose, Param.Rythm);
    FEP.find_threshold(True, Param.WriteHisto,Param.Verbose, Param.Rythm);
		      
		      
    //Allocation of cube for the first multiresolution transformation 
    WT_tools.alloc(WT_tools.Wavelet_Coef, Nx, Ny, Nz, Param.Nbr_Plan);
    
    fltarray Output_Cube(Nx,Ny,Nz);//Allocation of output cube 
         
    
    //creation of the iterator
    if(Param.CFAIter == True)
    { 
    	Iteration_CFA Iter(Event_of_V0, Output_Cube, Support,Nb_Event_In_Boxes, 
    		           Mean_Nb_Ev, WT_tools, FEP, Param);
	for(Iter.init();Iter.test_end();++Iter){}
    }
    else
    { 
    	Iteration Iter(Event_of_V0, Output_Cube, Support,Nb_Event_In_Boxes, 
    		Mean_Nb_Ev,WT_tools, FEP, Param); 
	//iterator
    	for(Iter.init();Iter.test_end();Iter++){}
    }
    HD.origin = Cmd;
    if ((Param.NbIter==0 || Param.WriteRes == False)) 
                                       write_result(Output_Cube,Param, &HD);

    WT_tools.free(WT_tools.Wavelet_Coef, Param.Nbr_Plan);
    free(Support,Param.Nbr_Plan);
    if(Param.Verbose) cout << "Done."<<endl;
    exit(0);
}

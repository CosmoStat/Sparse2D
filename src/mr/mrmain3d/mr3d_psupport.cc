/******************************************************************************
**                   Copyright (C) 1996 by CEA
*******************************************************************************
*******************************************************************************
**
**    DESCRIPTION  perform the multi-resolution support using 
**		   the abaque.
****************************************************************************/

#include "mr3d_psupport.h"

/*****************************************************************************/
void Projection_V0(fltarray & Event_of_V0, intarray & Catalogue_Event, 
		FewEventPoisson & FEP){
   //		Projection in V0 if input data is a fits

   
   int Nx = Event_of_V0.nx();
   int Ny = Event_of_V0.ny();
   int Nz = Event_of_V0.nz();
   
      
   
   // char Name_Cube_Out2[256];
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
         		
	
}

/*********************************************************************/
//		Projection in V0 if input data is a catalogue
void Projection_V0(ArrayPoint & Catalogue, fltarray & Event_of_V0, 
		intarray & Catalogue_Event){
	//initialisation of Event_of_V0 and Catalogue_Event
	Event_of_V0.init(0.);
	Catalogue_Event.init(0);
	int Nx=	Event_of_V0.nx();
	int Ny=	Event_of_V0.ny();
	int Nz=	Event_of_V0.nz();
	
	
		//projection in V0
		// B(x) =1/12 (|x-2|^3 - 4 |x-1|^3 + 6 |x|^3 - 4 |x+1|^3 + |x+2| ^ 3)
		fltarray bspline(401);
		for(int i=-200;i<=200;i++)//bspline define beetwin -2 and 2
			bspline(i+200)=(CUBE_FABS(i/100.-2.)-4.*CUBE_FABS(i/100.-1.)
				+6.*CUBE_FABS(i/100.)-4.*CUBE_FABS(i/100.+1.)
				+CUBE_FABS(i/100.+2.))/12.;
		for(int i=0;i<Catalogue.np();i++){
			int pos_x=(int)((Catalogue(i).x()-Catalogue.Pmin.x())
				/bin_cat);
			int pos_y=(int)((Catalogue(i).y()-Catalogue.Pmin.y())
				/bin_cat);
			int pos_z=(int)((Catalogue(i).z()-Catalogue.Pmin.z())
				/bin_cat);
			int offset_x=(int)(((Catalogue(i).x()-Catalogue.Pmin.x())
				/bin_cat-pos_x)*100.);
			int offset_y=(int)(((Catalogue(i).y()-Catalogue.Pmin.y())
				/bin_cat-pos_y)*100.);
			int offset_z=(int)(((Catalogue(i).z()-Catalogue.Pmin.z())
				/bin_cat-pos_z)*100.);
			Catalogue_Event(pos_x,pos_y,pos_z)++;
			for(int cur_x=-1;cur_x<3;cur_x++)
			for(int cur_y=-1;cur_y<3;cur_y++)
			for(int cur_z=-1;cur_z<3;cur_z++)
				if (!(pos_x+cur_x<0 || pos_x+cur_x>=Nx) && !(pos_y+cur_y<0 || pos_y+cur_y>=Ny)
					&& !(pos_z+cur_z<0 || pos_z+cur_z>=Nz)) 
					Event_of_V0(pos_x+cur_x,pos_y+cur_y,pos_z+cur_z)+=
						bspline(200+cur_x*100-offset_x)*
						bspline(200+cur_y*100-offset_y)*
						bspline(200+cur_z*100-offset_z);	
		}
	
}

/*********************************************************************/
//	read input data and test if it is a catalogue or a fits
void read_data(fltarray & Event_of_V0, intarray & Catalogue_Event,
		fltarray & Filtred_Input, Bool Filter_Def,
		FewEventPoisson & FEP,int & Nx,int & Ny,int & Nz)
{
		
	// reading input data 
	//test if input data is a fltarray or a catalogue
	if(strstr(Name_Cube_In,".fits") || strstr(Name_Cube_In,".Fits")
			|| strstr(Name_Cube_In,".FITS")){ 
		FitsMod=True;
		io_3d_read_data(Name_Cube_In, Event_of_V0);
		Catalogue_Event.alloc(Event_of_V0.nx(),Event_of_V0.ny(),Event_of_V0.nz());
		Projection_V0(Event_of_V0,Catalogue_Event,FEP);//Projection in V0
		
	}else if (strstr(Name_Cube_In,".cat") || strstr(Name_Cube_In,".Cat")
			|| strstr(Name_Cube_In,".CAT")){
		FitsMod=False;
		ArrayPoint Catalogue;//initialisation of catalogue
    		Catalogue.read(Name_Cube_In);
		if(Catalogue.dim()!=3){
			cout <<"The catalogue do not have 3 dimentions"<<endl;
			exit(0);
		}
		int npx=(int)((Catalogue.Pmax.x()-Catalogue.Pmin.x())/bin_cat+0.5);
		int npy=(int)((Catalogue.Pmax.y()-Catalogue.Pmin.y())/bin_cat+0.5);
		int npz=(int)((Catalogue.Pmax.z()-Catalogue.Pmin.z())/bin_cat+0.5);
		Event_of_V0.alloc(npx,npy,npz);
		Catalogue_Event.alloc(npx,npy,npz);
		Projection_V0(Catalogue,Event_of_V0,Catalogue_Event);
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
	if ( (min < 32) && (Nbr_Plan>2 || !Nbr_Plan)) 
		Nbr_Plan=2;
      	else if ( min < 64 && (Nbr_Plan>3 || !Nbr_Plan)) 
		Nbr_Plan=3;
      	else if ( min < 128 && (Nbr_Plan>4 || !Nbr_Plan)) 
		Nbr_Plan=4;     
      	else if ( min < 256 && (Nbr_Plan>5 || !Nbr_Plan)) 
		Nbr_Plan=5; 
      	else if ( min < 512 && (Nbr_Plan>6 || !Nbr_Plan)) 
		Nbr_Plan=6;
      	else if ( min < 1024 && (Nbr_Plan>7 || !Nbr_Plan)) 
		Nbr_Plan=7;
      	else if ( min < 2048 && (Nbr_Plan>8 || !Nbr_Plan)) 
		Nbr_Plan=8;
      	else if ( min < 2*2048 && (Nbr_Plan>9 || !Nbr_Plan)) 
		Nbr_Plan=9;
      	else if ( min < 2*4096 && (Nbr_Plan>10 || !Nbr_Plan)) 
		Nbr_Plan=10;
	
	if(Filter_Def){ 
		io_3d_read_data(Filtred_FileName, Filtred_Input);
		if(Event_of_V0.nx() !=Filtred_Input.nx() || Event_of_V0.ny() !=Filtred_Input.ny() ||
		 		Event_of_V0.nz() !=Filtred_Input.nz()){
			cout << "Filtered Input do not correspond with Input Data"<<endl;
			exit(-1);
		}
	}	
	
	
}


/*********************************************************************/
void init_Nb_Event(intarray & Catalogue_Event,intarray *& Nb_Event,
		FewEventPoisson & FEP,
		ATROUS_3D_FIL & WT_tools){
	//	compute the cubes Nb_Event
	int Nx = Catalogue_Event.nx();
    	int Ny = Catalogue_Event.ny();
    	int Nz = Catalogue_Event.nz();
	// char Name_Cube_Out2[50];
	if(Verbose) cout << "Initalisation of Nb_Event"<<endl;
    	for(int s=-2;s<Nbr_Plan-1;s++){
		if (Approx){//we do not compute nb of events in boxes of 
			// (2^n+1)^3 but boxes of (2^n)^3
			if(s==-2)  FEP.event_scale_3D_Approx( Catalogue_Event, 
				Nb_Event[0],s, Border, Verbose );
			else if(s==-1)  FEP.event_scale_3D_Approx( Nb_Event[0], 
				Nb_Event[1],s, Border, Verbose );
			else if(s==0)  FEP.event_scale_3D_Approx( Nb_Event[1], 
				Nb_Event[0],s, Border, Verbose );	
    			else FEP.event_scale_3D_Approx( Nb_Event[s-1], Nb_Event[s], 
				s, Border, Verbose );  
		}else //we compute nb of events in boxes of (2^n+1)^3  
			if(s>=0)  FEP.event_scale_3D( Catalogue_Event, Nb_Event, 
					s, Border, Verbose); 

    	}
	if(Approx){
	   intarray Buff(Nx,Ny,Nz);
	   for(int s=0;s<Nbr_Plan-1;s++){
	   	//for a better appoximation of counting boxes 
		//it computes the equivalant of a box "9*9*9"
		//so it multiply by 9*9*9/(8*8*8)
		//and it makes an average with the 7 neighbors of right-hand side
	   	float mult=pow((double)(1+ pow(2.,(double)-s-3.)),3.)/8.;
		Buff=Nb_Event[s];
		for(int i=0;i<Nx-1;i++)
		for(int j=0;j<Ny-1;j++)
		for(int k=0;k<Nz-1;k++)
			Nb_Event[s](i,j,k)=(int)((Buff(i,j,k)+Buff(i,j,k+1)+
			   Buff(i,j+1,k)+Buff(i,j+1,k+1)+Buff(i+1,j,k)+
			   Buff(i+1,j,k+1)+Buff(i+1,j+1,k)+Buff(i+1,j+1,k+1))
			   *mult +0.5);
		
	        mult=pow((double)(1+ pow(2.,(double)-s-3.)),3.)/4.;
		
		
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
	
}

/*********************************************************************/
void init_support(intarray *& Nb_Event,intarray *& Support, 
		ATROUS_3D_FIL & WT_tools, FewEventPoisson & FEP){
		//	read or compute the support 
    char Name_Cube_Out2[50];
    if (ReadSupport)
    {
    	for(int s=0;s<Nbr_Plan-1 ;s++)
	{
	    sprintf(Name_Cube_Out2,"%s_support_%d",Name_Support,s+1);
	    fits_read_intarr(Name_Cube_Out2, Support[s]);
	    if(Verbose) cout << "The support of scale "
				<<s<<" is read."<<endl;   
	}
    }
    else
    {
    	for(int s=0;s<Nbr_Plan-1 ;s++)
	{
    	    FEP.event_set_support(WT_tools.Wavelet_Coef, s,FirstScale, FirstScale, MinEvent, 
			          KeepPositivSup, I_MIRROR,Nb_Event, Support, 
			          False);
            if (Verbose) cout << "Support of scale number " << s+1 << " completed"<< endl; 
    	}
    }

    //filtring of single points in the support 
    if(!SinglePoint) 
    {
    	WT_tools.No_Single_Point(Support[0].nx(),Support[0].ny(),
		Support[0].nz(), Support, FEP, Nbr_Plan-1,Verbose);
    	if(Verbose) cout << "Single points have been cleared"<<endl;
    }
        
}


/***********************************************************************/
 

static void usage(char *argv[])
{

    fprintf(OUTMAN,"\n\nUsage: %s options in_cube out_mr_file\n\n", argv[0]);
    fprintf(OUTMAN, "   where options =  \n");

    
    first_detect_scale_usage();
    
    min_event_usage(MinEvent);   
   
//     fprintf(OUTMAN, "         [-s SignifStructureAnalysis_FileName]\n");
//     fprintf(OUTMAN, "             Analyse the detected wavelet coefficients,\n");
//     fprintf(OUTMAN, "             and write in the file:\n");
//     fprintf(OUTMAN, "               Number of detected structures per scale\n");
//     fprintf(OUTMAN, "               Percentage of significant wavelet coefficents\n");
//     fprintf(OUTMAN, "               Mean deviation of shape from sphericity\n");
//     fprintf(OUTMAN, "               For each detected structure, its volume, its surface, and\n");
//     fprintf(OUTMAN, "               its deviation of shape from sphericity, \n");
//     fprintf(OUTMAN, "               its angles, its elongation in the three axis directions.\n");
// 
// 
//     fprintf(OUTMAN, "         [-t SignifStructureAnalysis_FileName]\n");
//     fprintf(OUTMAN, "             Same as -s option, but results are stored\n");
//     fprintf(OUTMAN, "             in an ascii table format.\n");
//     fprintf(OUTMAN, "             The table contains: scale number, structure number, \n");
//     fprintf(OUTMAN, "             Max_x, Max_y, Max_z, Volume, Surface, Morpho, \n");
//     fprintf(OUTMAN, "             Angle_Theta, Angle_Phi, Sigma_1, Sigma_2, Sigma_3. \n"); 
//     
    detect_pos_usage();
    
//     fprintf(OUTMAN, "         [-f filename_filtered_data]\n");
//     fprintf(OUTMAN, "             File with filtered input data.\n");

    fprintf(OUTMAN, "         [-n number_of_scales]\n");
    fprintf(OUTMAN, "             Number of scales used in the multiresolution transform.\n");
        
    verbose_usage();    
    
    fprintf(OUTMAN, "         [-E epsilon]\n");
    fprintf(OUTMAN, "             Value of epsilon used for the threshold ( default=1e-3 ).\n");
    
    fprintf(OUTMAN, "         [-k]\n");
    fprintf(OUTMAN, "             Suppress isolated pixels in the support. Default is no.\n");
    
    fprintf(OUTMAN, "         [-R]\n");
    fprintf(OUTMAN, "             Remove the structures near the borders.\n");
    
//     fprintf(OUTMAN, "         [-b type_border = 1 to 4]\n");
//     fprintf(OUTMAN, "              give the type of border {1:I_CONT|2:I_MIRROR|3:I_ZERO|4:I_SurfVolOD}\n");
//     fprintf(OUTMAN, "              ( default=3:I_ZERRO ).\n");
    
    fprintf(OUTMAN, "         [-B bin] (only for catalogue)\n");
    fprintf(OUTMAN, "             If the input data is a catalogue give the bin of it ( default=1 ).\n");
   
//     fprintf(OUTMAN, "         [-a]\n");
//     fprintf(OUTMAN, "              No approximation in computing the number of events in boxes.\n");
// 	
//    fprintf(OUTMAN, "         [-w]\n");
//    fprintf(OUTMAN, "             Write the multiscale support to the disk.\n");

//     fprintf(OUTMAN, "         [-S]\n");
//     fprintf(OUTMAN, "             Read the support already computed.\n");
    vm_usage();
    manline();    
    
    exit(-1);
}

 
/*********************************************************************/

/* GET COMMAND LINE ARGUMENTS */
static void psupportinit(int argc, char *argv[])
{
    int c;
#ifdef LARGE_BUFF
    int VMSSize=-1;
    Bool OptZ = False;
    char VMSName[1024] = "";
#endif    
   
    strcpy(Name_Abaque, Name_Abaque_Default);

    /* get options */
    while ((c = GetOpt(argc,argv,"F:RSwkan:e:f:B:E:ps:t:vzZ:b:")) != -1) 
    {
	switch (c) 
      	{
 	 case 'k': SinglePoint=False; break;
	 case 'v': Verbose = True; break;
	 case 'R': RemoveBorder=True; break;
	 case 'a': Approx = (Approx == True) ? False: True; break;
	 case 'w': WriteSupport=True; break;
	 case 'S': ReadSupport=True; break;
	 case 'E':
		/* -E <epsilon> */
		if ((sscanf(OptArg,"%f",&epsilon) != 1)|| (epsilon < 0)){
		    fprintf(OUTMAN, "bad value of epsilon: %s\n", OptArg);
		    exit(-1);
		}
		break;
	 case 'p':
		 KeepPositivSup=True;
		break; 
	   case 'F': /* begin the detection */
		if (sscanf(OptArg,"%d", &FirstScale) != 1) 
                {
		   fprintf(OUTMAN, "Error: bad first detection scale number: %s\n", OptArg);
		   exit(-1);
		}
		FirstScale --;
		break;
	 case 'e':
 		/* -n <Nbr_Plan> */
		if (sscanf(OptArg,"%d",&MinEvent) != 1) 
                {
		    fprintf(OUTMAN, "Error: Min number of events: %s\n", OptArg);
		    exit(-1);
		}
                if (MinEvent < 1)  
                {
		    fprintf(OUTMAN, "Error: bad number of events min: %s\n", OptArg);
		    exit(-1);
		}
		break;
 	   
 	   case 's': /* info_filename */
		if (sscanf(OptArg,"%s",AnaFileName) != 1)
                {
		   fprintf(OUTMAN, "Error: bad file name: %s\n", OptArg);
		   exit(-1);
		}
		ObjAna = True;
		AsciiRes = True;
		break;
 	   case 't': /* info_ascii_filename */
		if (sscanf(OptArg,"%s",TabAnaFileName) != 1)
                {
		   fprintf(OUTMAN, "Error: bad file name: %s\n", OptArg);
		   exit(-1);
		}
		ObjAna = True;
		TabAsciiRes = True;
		break;
	   case 'f': /* Filtered_FileName */
		if (sscanf(OptArg,"%s",Filtred_FileName) != 1)
                {
		   fprintf(OUTMAN, "Error: bad file name: %s\n", OptArg);
		   exit(-1);
		}
		Filter_Def = True;
		break;
 	 case 'n':
 		/* -n <Nbr_Plan> */
		if (sscanf(OptArg,"%d",&Nbr_Plan) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad number of scales: %s\n", OptArg);
		    exit(-1);
		}
                if ((Nbr_Plan <= 1) || (Nbr_Plan > MAX_SCALE)) 
                {
		    fprintf(OUTMAN, "Error: bad number of scales: %s\n", OptArg);
		    fprintf(OUTMAN, "       1 < Nbr Scales <= %d\n", MAX_SCALE);
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
	case 'B':
		/* -B <bin> */
		if ((sscanf(OptArg,"%f",&bin_cat) != 1)|| (bin_cat < 0)){
		    fprintf(OUTMAN, "bad value of bin: %s\n", OptArg);
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
   
   /* get input filename for input cube */
   	if (OptInd < argc)  strcpy(Name_Cube_In,argv[OptInd++]);
        else usage(argv);

   /* get output filename for mr_support */
   	if (OptInd < argc)  strcpy(Name_Support,argv[OptInd++]);
        else usage(argv);

   /* make sure there are not too many parameters */
	if (OptInd < argc) {
		fprintf(OUTMAN, "Error: too many parameters: %s ...\n", argv[OptInd]);
		exit(-1);
	}

#ifdef LARGE_BUFF
    if (OptZ == True) vms_init(VMSSize, VMSName, Verbose);
#endif

}



/************************************************************************/

int main(int argc, char *argv[])
{
 int s,Nx,Ny,Nz;	

 // get command line 
    lm_check(LIC_MR3);
    psupportinit(argc,argv);
    

    intarray Catalogue_Event;	//input catalogue
    fltarray Event_of_V0;	//same as Catalogue_Event but projected 
    				//in V0 (if necessary)
	
    fltarray Filtred_Input;     //Filtred_Input MUST be the filtered of 
    				//Event_of_V0 of Catalogue_Event
    // char Name_Cube_Out2[256];
    //creating Atrous_3D class 
    ATROUS_3D_FIL WT_tools; 
    //creating FewEventPoisson class
    FewEventPoisson FEP(True,0,epsilon);
      
    
    read_data(Event_of_V0,Catalogue_Event,Filtred_Input,Filter_Def,
    	FEP,Nx,Ny,Nz);// reading input data 
    if(Verbose)
    {       
       cout << endl << endl << "PARAMETERS: " << endl << endl;
       cout << "File Name in = " << Name_Cube_In << endl;
       cout << "File Name out = " <<  Name_Support << endl; 
       cout <<"Dimentions of the cube are: "<<Nx<<"  "<<Ny<<"  "<<Nz<<"  "<<endl;
       cout <<"The Number of Scale is: "<<Nbr_Plan<<endl;
       if (SinglePoint == False) cout << "Kill isolated pixels" << endl;
       else cout << "Keep isolated pixels" << endl;
       cout <<"Epsilon = " << epsilon << endl;
       if (KeepPositivSup == False) cout << "Detect both positive and negative coeff." << endl;
       else cout << "Detect only positive coeff." << endl;
       if (ObjAna == True)
           cout << "Wavelet Significant Structures Analysis" << endl;
       else cout << "No Wavelet Significant Structures Analysis" << endl;
       if (RemoveBorder == True) cout << "Remove structures on the border " << endl;
       else cout << "Keep structures on the border " << endl;
       cout << endl << endl;       
    }
    intarray *Nb_Event_In_Boxes;//cubes of number of event per box
    if(!ReadSupport){
    	alloc_array(Nb_Event_In_Boxes,Nx,Ny,Nz,Nbr_Plan);
    	//initalisation of Nb_Event
    	init_Nb_Event(Catalogue_Event,Nb_Event_In_Boxes,FEP,WT_tools);
    }
    
    intarray *Support;//Allocation of the support for each plan
    alloc_array(Support,Nx,Ny,Nz,Nbr_Plan);
 
    FEP.find_threshold( False,False,Verbose, 0);
  
    //Allocation of cube for the first multiresolution transformation 
    //WT_tools.Init_Wavelet_Coef(Nx, Ny, Nz, Nbr_Plan);
    WT_tools.alloc(WT_tools.Wavelet_Coef, Nx, Ny, Nz, Nbr_Plan);

    if(Filter_Def) Event_of_V0=Filtred_Input;
    WT_tools.wttransform(Event_of_V0, Nbr_Plan, False, False, "titi");
    
    init_support(Nb_Event_In_Boxes,Support,WT_tools,FEP);

    fltarray Cube(Nx, Ny, Nz);
 
    // if (RemoveBorder == True)
    for(s=0;s<Nbr_Plan-1 ;s++)
    {
//        char Name_Cube_Out2[50];
//        char Name_Cube_Out3[50];
//        sprintf(Name_Cube_Out2,"%s_seg_%d",Name_Support,s+1);
//        sprintf(Name_Cube_Out3,"%s_xx_%d",Name_Support,s+1);
//        fits_write_intarr(Name_Cube_Out3, Support[s]);
       int i,j,k,nmax;
       for (i=0;i<Nx; i++)
       for (j=0;j<Ny; j++)  
       for (k=0;k<Nz; k++) 
         if (Support[s](i,j,k) == 1) Cube(i,j,k) = 1;
	 else Cube(i,j,k) = 0;
       if(Verbose) cout<< "SCALE " << s+1 << endl;
       cube_segment(Cube, Support[s], nmax, 1e-10, RemoveBorder,1);
//        fits_write_intarr(Name_Cube_Out2, Support[s]);
       if (Verbose) cout<< "SCALE " << s+1 << ": Number of structures = " 
	     		                                   << nmax << endl;
    }
    

 
    for(s=0;s<Nbr_Plan-1 ;s++)
    {
       int i,j,k;
       for (i=0;i<Nx; i++)
       for (j=0;j<Ny; j++)  
       for (k=0;k<Nz; k++)
         if (Support[s](i,j,k) != 0) Support[s](i,j,k) = 1;
       char Name_Cube_Out2[50];
       sprintf(Name_Cube_Out2,"%s_support_%d",Name_Support,s+1);
       fits_write_intarr(Name_Cube_Out2, Support[s]);
    }
    
    WT_tools.threshold(Support,Nbr_Plan);
    io_set_format(F_UNKNOWN);//?????????????
 

  // segmentation 
 if (ObjAna == True)
  {
      int Cpt,ind,nmax;
      float Signif;
      ofstream Fic;
      ofstream Fic1;
      int  *TabXmax, *TabYmax, *TabZmax;
      float *TabSurfVol, *TabVolume;
      float MeanMorpho, Morpho;
      float *TabMorpho, *TabValMax;
      float *TabMx, *TabMy, *TabMz, *TabMean, *TabVx;
      float *TabVy, *TabVz, *TabVxy, *TabVxz, *TabVyz, *TabSx;
      float *TabSy, *TabSz, *TabTeta, *TabPhi;
      fltarray dia(3,3);
      fltarray valp(3);
      fltarray vecp(3,3);
      fltarray mat(3,3);
      if (AsciiRes == True)
      {
         Fic.open(AnaFileName);
         if(!Fic)
         {
             cerr<<"Error: Unable to open "<< AnaFileName << endl;
             exit(1);
         }
      }
      if (TabAsciiRes == True)
      {
         Fic1.open(TabAnaFileName);
         if(!Fic)
         {
             cerr<<"Error: Unable to open "<< TabAnaFileName << endl;
             exit(1);
         } 
      }

      for (int s = FirstScale; s < Nbr_Plan-1; s++)
      {
          Cube = WT_tools.Wavelet_Coef[s];
          cube_segment (Cube, Support[s], nmax, 1e-10, RemoveBorder,1);
	  if(Verbose) cout<< "SCALE " << s+1 << ": Number of structures = " 
	     		<< nmax << endl;
          if (AsciiRes == True)
          {
             Fic << "SCALE " << s+1 << ": Number of structures = " 
	     		<< nmax << endl;
          }

          if (nmax > 0)
          {
          TabSurfVol = new float [nmax+1];
          TabVolume = new float [nmax+1];
          TabMorpho = new float [nmax+1];
          TabXmax = new int [nmax+1];
          TabYmax = new int [nmax+1];
	  TabZmax = new int [nmax+1];
          TabValMax  = new float [nmax+1];
          TabMx  = new float [nmax+1];
          TabMy  = new float [nmax+1];
	  TabMz  = new float [nmax+1];
          TabMean  = new float [nmax+1];
          TabVx  = new float [nmax+1];
          TabVy  = new float [nmax+1];
	  TabVz  = new float [nmax+1];
          TabVxy  = new float [nmax+1];
	  TabVxz  = new float [nmax+1];
	  TabVyz  = new float [nmax+1];
          TabSx  = new float [nmax+1];
          TabSy  = new float [nmax+1];
	  TabSz  = new float [nmax+1];
          TabTeta  = new float [nmax+1];
	  TabPhi   = new float [nmax+1];
          Cpt = 0;
          for (int i = 0; i <= nmax; i++) 
          {
              TabSurfVol[i] = TabVolume[i] = 0;
              TabMorpho[i] = 0.;
              TabXmax[i] = -1;
              TabYmax[i] = -1;
	      TabZmax[i] = -1;
              TabValMax [i] = 0.;
              TabMx[i] = 0.;
              TabMy[i] = 0.;
	      TabMz[i] = 0.;
              TabMean[i] = 0.;
              TabVx[i] = 0.;
              TabVy[i] = 0.;
	      TabVz[i] = 0.;
              TabVxy[i] = 0.;
	      TabVxz[i] = 0.;
	      TabVyz[i] = 0.;
              TabSx[i] = 0.;
              TabSy[i] = 0.;
	      TabSz[i] = 0.;
              TabTeta[i] = 0.;
	      TabPhi[i] = 0.;
          }
          for (int i = 0; i < Nx; i++)
          for (int j = 0; j < Ny; j++)
	  for (int k = 0; k < Nz; k++)
          if (Support[s](i,j,k) != 0)
          {

              float Coef = Cube(i,j,k);
              Cpt ++;
              ind = (int) (Support[s](i,j,k) + 0.5);
              TabVolume[ind] ++;
              TabMx[ind] += Coef*i;
              TabMy[ind] += Coef*j;
	      TabMz[ind] += Coef*k;
              TabMean[ind] += Coef;
              TabVx[ind] += Coef*i*i;
              TabVy[ind] += Coef*j*j;
	      TabVz[ind] += Coef*k*k;
              TabVxy[ind] +=  Coef*i*j;
	      TabVyz[ind] +=  Coef*j*k;
	      TabVxz[ind] +=  Coef*i*k;
              
	      if (ABS(Cube(i,j,k)) > TabValMax[ind])
              {
                 TabValMax[ind] = ABS(Cube(i,j,k));
                 TabXmax[ind] = i;
                 TabYmax[ind] = j;
		 TabZmax[ind] = k;
              }
              if ((i==0)||(i==Nx-1)) TabSurfVol[ind]++;
	      if ((j==0)||(j==Ny-1)) TabSurfVol[ind]++;
	      if ((k==0)||(k==Nz-1)) TabSurfVol[ind]++;                
		
              if ( (i>0) && (j>0) && (k>0) && (i<Nx-1) && (j<Ny-1) && (k<Nz-1)){
	      	float sq2=sqrt(2.)-1;
	      	float sq3=sqrt(3.)-1;
		
                if (Support[s](i-1,j,k) < 0.5) TabSurfVol[ind]++;
		else if (Support[s](i,j-1,k) < 0.5) TabSurfVol[ind]++;
		else if (Support[s](i,j,k-1) < 0.5) TabSurfVol[ind]++;
		else if (Support[s](i+1,j,k) < 0.5) TabSurfVol[ind]++;
		else if (Support[s](i,j+1,k) < 0.5) TabSurfVol[ind]++;
		else if (Support[s](i,j,k+1) < 0.5) TabSurfVol[ind]++;
		
		else if (Support[s](i-1,j-1,k) < 0.5) TabSurfVol[ind]+=sq2;
		else if (Support[s](i-1,j,k-1) < 0.5) TabSurfVol[ind]+=sq2;
		else if (Support[s](i-1,j+1,k) < 0.5) TabSurfVol[ind]+=sq2;
		else if (Support[s](i-1,j,k+1) < 0.5) TabSurfVol[ind]+=sq2;
		
		else if (Support[s](i+1,j-1,k) < 0.5) TabSurfVol[ind]+=sq2;
		else if (Support[s](i+1,j,k-1) < 0.5) TabSurfVol[ind]+=sq2;
		else if (Support[s](i+1,j+1,k) < 0.5) TabSurfVol[ind]+=sq2;
		else if (Support[s](i+1,j,k+1) < 0.5) TabSurfVol[ind]+=sq2;
		
		else if (Support[s](i-1,j-1,k) < 0.5) TabSurfVol[ind]+=sq2;
		else if (Support[s](i,j-1,k-1) < 0.5) TabSurfVol[ind]+=sq2;
		else if (Support[s](i+1,j-1,k) < 0.5) TabSurfVol[ind]+=sq2;
		else if (Support[s](i,j-1,k+1) < 0.5) TabSurfVol[ind]+=sq2;
		
                else if (Support[s](i-1,j+1,k) < 0.5) TabSurfVol[ind]+=sq2;
		else if (Support[s](i,j+1,k-1) < 0.5) TabSurfVol[ind]+=sq2;
		else if (Support[s](i+1,j+1,k) < 0.5) TabSurfVol[ind]+=sq2;
		else if (Support[s](i,j+1,k+1) < 0.5) TabSurfVol[ind]+=sq2;
		
		else if (Support[s](i-1,j,k-1) < 0.5) TabSurfVol[ind]+=sq2;
		else if (Support[s](i,j-1,k-1) < 0.5) TabSurfVol[ind]+=sq2;
		else if (Support[s](i+1,j,k-1) < 0.5) TabSurfVol[ind]+=sq2;
		else if (Support[s](i,j+1,k-1) < 0.5) TabSurfVol[ind]+=sq2;
		
		else if (Support[s](i-1,j,k+1) < 0.5) TabSurfVol[ind]+=sq2;
		else if (Support[s](i,j-1,k+1) < 0.5) TabSurfVol[ind]+=sq2;
		else if (Support[s](i+1,j,k+1) < 0.5) TabSurfVol[ind]+=sq2;
		else if (Support[s](i,j+1,k+1) < 0.5) TabSurfVol[ind]+=sq2;
		
		else if (Support[s](i-1,j-1,k-1) < 0.5) TabSurfVol[ind]+=sq3;
		else if (Support[s](i-1,j-1,k+1) < 0.5) TabSurfVol[ind]+=sq3;
		else if (Support[s](i-1,j+1,k-1) < 0.5) TabSurfVol[ind]+=sq3;
		else if (Support[s](i-1,j+1,k+1) < 0.5) TabSurfVol[ind]+=sq3;
		else if (Support[s](i+1,j-1,k-1) < 0.5) TabSurfVol[ind]+=sq3;
		else if (Support[s](i+1,j-1,k+1) < 0.5) TabSurfVol[ind]+=sq3;
		else if (Support[s](i+1,j+1,k-1) < 0.5) TabSurfVol[ind]+=sq3;
		else if (Support[s](i+1,j+1,k+1) < 0.5) TabSurfVol[ind]+=sq3;
		
              }
          }

          Signif =  (float) Cpt / (float)(Nx*Ny*Nz) * 100.;
          if (AsciiRes == True)
          {
             Fic << "   Pourcentage of Significant Coefficients = " << Signif << endl;
          }

          MeanMorpho = 0.;
          for (int i = 1; i <= nmax; i++)
          {
             
	     Morpho=float(TabVolume[i]) * 6 / float(TabSurfVol[i]) * sqrt(PI/float(TabSurfVol[i]));
	     if(Morpho<1e-10) Morpho=0.;
	     
	     MeanMorpho += Morpho;
             TabMorpho[i] = Morpho;

             TabMx[i]  /= TabMean[i];
             TabMy[i]  /= TabMean[i];
	     TabMz[i]  /= TabMean[i];
	     
             TabVx[i]  /= TabMean[i];
             TabVy[i]  /= TabMean[i];
	     TabVz[i]  /= TabMean[i];
             TabVxy[i] /= TabMean[i];
	     TabVxz[i] /= TabMean[i];
	     TabVyz[i] /= TabMean[i];
	     
             TabVx[i]  -= TabMx[i]*TabMx[i];
             TabVy[i]  -= TabMy[i]*TabMy[i];
	     TabVz[i]  -= TabMz[i]*TabMz[i];
             TabVxy[i] -= TabMx[i]*TabMy[i];
	     TabVyz[i] -= TabMy[i]*TabMz[i];
	     TabVxz[i] -= TabMx[i]*TabMz[i];
	     
	     
	   
	     
	     mat(0,0)=TabVx[i];
	     mat(1,1)=TabVy[i];
	     mat(2,2)=TabVz[i];
	     mat(0,1)=mat(1,0)=TabVxy[i];
	     mat(0,2)=mat(2,0)=TabVxz[i];
	     mat(1,2)=mat(2,1)=TabVyz[i];
	     
	     int  nrot;
	     
	     diag_matrix(mat, dia, valp, vecp, nrot, True);
	     TabTeta[i]=-atan(vecp(0,2)/vecp(0,1));
	     TabPhi[i]=(fabs(vecp(0,0))<1e-7 ? 
	     	acos(-vecp(1,2)/vecp(1,0)):asin(-vecp(0,1)/vecp(0,0)));
	     TabSx[i]=sqrt(fabs(valp(0)));
	     TabSy[i]=sqrt(fabs(valp(1)));
	     TabSz[i]=sqrt(fabs(valp(2)));
	    
		
		
		
          }
          if (nmax > 0) MeanMorpho /= (float) nmax;

          if (AsciiRes == True)
          {
             Fic << "   Mean deviation of shape from sphericity = " << MeanMorpho << endl << endl;
          }

          for (int i = 1; i <= nmax; i++)
          {
             if (AsciiRes == True)
             {
                Fic << "   S" << s+1 << ":" << i << " Volume = " << TabVolume[i] << "  SurfVol = " 
			<< TabSurfVol[i] << "  Morpho = " << TabMorpho[i] << endl;
                Fic << "       Angle Teta = " << TabTeta[i] << "       Angle Phi = " << TabPhi[i]<<endl; 
		Fic << " Sigma1 = " <<TabSx[i] << " Sigma2 = " <<TabSy[i]<< " Sigma3 = " <<TabSz[i]<<endl; 
             }
             if (TabAsciiRes == True)
             {  
                Fic1 << s+1 <<" " << i <<" " << TabXmax[i] <<" " << TabYmax[i] <<" "<< TabZmax[i] <<" "
			<< TabVolume[i] <<" " << TabSurfVol[i] <<" " << TabMorpho[i]
			<< " " << TabTeta[i]<<" " << TabPhi[i] << " " << TabSx[i] << " " << TabSy[i] <<" " << TabSz[i] << endl;
             }
          }
          delete [] TabSurfVol;
          delete [] TabVolume;
          delete [] TabMorpho;
          delete [] TabXmax;
          delete [] TabYmax;
	  delete [] TabZmax;
          delete [] TabValMax;
          delete [] TabMx;
          delete [] TabMy;
	  delete [] TabMz;
          delete [] TabMean;
          delete [] TabVx;
          delete [] TabVy;
	  delete [] TabVz;
          delete [] TabVxy;
	  delete [] TabVxz;
	  delete [] TabVyz;
          delete [] TabSx;
          delete [] TabSy;
	  delete [] TabSz;
          delete [] TabTeta;
	  delete [] TabPhi;
	  
          }
         // write results 
       }
      if (AsciiRes == True) Fic.close();
      if (TabAsciiRes == True) Fic1.close();
       //MR_Data.write("xx_Segment.mr");
  }
  

  exit(0);
}

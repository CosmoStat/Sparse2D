/******************************************************************************
**                   Copyright (C) 2009 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Florent Sureau
**
**    Date:  09/12/09
**    
**    File:  mrs_planck_mgmca.cc
**
*******************************************************************************
**
**    DESCRIPTION  Multiscale Generalized Morphological Component Analysis for planck data
**    ----------- 
**                 
**    Usage: mrs_planck_mgmca options cube output
**
******************************************************************************/

#include "mrs_planck_lGMCA.h"
#include <omp.h>
extern int inner_loop_threads;
extern int outer_loop_threads;


/*********************************************************************/
void lGMCA_inv::Read_input_maps(char *Name) {
	fits_read_fltarr(Name,InputMaps);
	//Check 2D array
	if(InputMaps.naxis() > 2) {
		InputMaps.free();
		error(EXIT_FAILURE,"lGMCA->Size of input maps is not 2D",":Read_input_maps()");
	}
	//Check NchanMax consistent
	NchanMax=InputMaps.ny();
	if(Nbeams > 0) {
		if(NchanMax > Nbeams) {
			InputMaps.free();
			Beams.free();
			error(EXIT_FAILURE,"lGMCA->Not enough beams",":Read_input_maps()");
		}
	} 
	Npix=InputMaps.nx();
	Nside=floor(sqrt(Npix / 12));
	if (Debug) printf("Allocate memory for ALM transform (Nchanmax=%d)\n",NchanMax);
	AlmX=new Alm<xcomplex<double> >[NchanMax];//Will contain the ALMs of the input maps
	AlmXComp.alloc(NchanMax);//Alm X already performed if =1
	AllocALMMem=true; 
}

/*********************************************************************/
void lGMCA_inv::Read_Beams(char *Name) {
	
	fits_read_fltarr(Name,Beams);
	Nbeams=Beams.ny();

	//Check NchanMax consistent
	if(NchanMax > 0) {
		if(Nbeams < NchanMax) {
			InputMaps.free();
			Beams.free();
			error(EXIT_FAILURE,"lGMCA->NchanMax does not agree",":Read_input_maps()");
		}
	} 
}


/*********************************************************************/
void lGMCA_inv::Read_inversion_vectors(char *Name) {
	
	if(InversionVectors.n_elem() >0 ) {
		if(Verbose) printf("Deallocate inversion vector\n");
		InversionVectors.free();
	}
	
	fits_read_fltarr(Name,InversionVectors);
	
	//Check 2D or 3D array
	if((InversionVectors.naxis() != 3)&&(InversionVectors.naxis() != 2)) {
		printf("Reading %s: InversionVectors.naxis()=%d\n",Name,InversionVectors.naxis());
		InversionVectors.free();
		error(EXIT_FAILURE,"lGMCA->Size of inversion vectors is neither 2D nor 3D",":Read_inversion_vectors()");
	}
	
	//Check Nchan consistent
	if(InversionVectors.naxis()==2) { //Case a single channel with different weightings for the planes
		if(Nchan > 0) {
			if(Nchan != 1) {
				InversionVectors.free();
				error(EXIT_FAILURE,"lGMCA->Nchan does not agree",":Read_inversion_vectors()");
			} 
		} else Nchan=1;
		Nscale=InversionVectors.ny();
	} else {
		if(Nchan > 0) {
			if(Nchan != InversionVectors.nz()) {
				printf("Nchan=%d,InversionVectors.ny()=%d\n",Nchan, InversionVectors.ny());
				InversionVectors.free();
				error(EXIT_FAILURE,"lGMCA->Nchan does not agree",":Read_inversion_vectors()");
			}
		} else Nchan=InversionVectors.nz();
		Nscale=InversionVectors.ny();
	}
	
	//Check nside consistent
	NsideRec=floor(sqrt(InversionVectors.nx() / 12));
	NpixRec=(unsigned long) InversionVectors.nx();
	if(Debug) printf("Use for this band NSide =%d\n",NsideRec);
}

/*********************************************************************/
void lGMCA_inv::Read_inversion_structure(char *Name) {
	
	fits_read_fltarr(Name,ResolStruc);
	
	//Check 1D or 2D array
	if((ResolStruc.naxis() != 1)&&(ResolStruc.naxis() != 2)) {
		printf("%s: ResolStruc.naxis()=%d\n",Name,ResolStruc.naxis());
		ResolStruc.free();
		error(EXIT_FAILURE,"lGMCA->Size of resolution structure is neither 1D nor 2D",":Read_inversion_structure()");
	}
	
	if(ResolStruc.naxis()==1) Nband=1;
	else Nband=ResolStruc.ny();	

}

/*********************************************************************/
void lGMCA_inv::SetOnlyResolution(char *Name_in, char *Name_beams, char *Name_inv_struc, char *Name_Out, int lNside) {
	int kc;
	intarray listchannel;
	
	struct timeval start,end,diff;
	NsideRec=lNside;
	NpixRec=NsideRec*NsideRec*12l;
	if(Timer==True) gettimeofday(&start, (struct timezone *) NULL);

	Read_input_maps(Name_in);//Read input maps
	Read_Beams(Name_beams);//Read beams
	Read_inversion_structure(Name_inv_struc); //read structure containing the list of channel per resolution band
	
	//Get Lmax per channel
	LmaxChan.alloc(NchanMax);
	LmaxChan.init(0.);
	if(Nband==1) {
		if(ResolStruc(0) > NchanMax) error(EXIT_FAILURE,"lGMCA->Nchan does not agree",":Process_AllBands()");
		if(ResolStruc(0) > ResolStruc.nx()-3) error(EXIT_FAILURE,"lGMCA->Not enough channels in structure",":Process_AllBands()");
		for(kc=0;kc<NchanMax;kc++) {
			LmaxChan(kc)=ResolStruc(ResolStruc.nx()-1);
			AlmX[kc].Set(LmaxChan(kc),LmaxChan(kc));
 			AlmX[kc].SetToZero ();
		}
	} else error(EXIT_FAILURE,"lGMCA-> Function to only set to same resolution only works for a single band",":SetOnlyResolution()");

	if(Debug==true) for(kc=0;kc<NchanMax;kc++) printf("LmaxChan[%d]=%d\n",kc,LmaxChan(kc));
	
	create_nest2ringindices(Nside);
	
	Nchan=ResolStruc(0);
	listchannel.alloc(Nchan+1);
	for (kc=0;kc<=Nchan;kc++) listchannel(kc)=ResolStruc(kc+1)-1;
	lLmax=ResolStruc(ResolStruc.nx()-1);
	if(Verbose==true) printf("Set to same resolution for band %d\n",0);
	for (kc=0;kc<Nchan;kc++) {
		if(Verbose==true) printf("Process channel %d (ref %d)\n",	listchannel(kc+1),listchannel(0));
		Set_Same_Resolution(kc,listchannel(kc+1),listchannel(0),LmaxChan(listchannel(kc+1)));	//Set to same resolution all channels in a Band	
	}	
	fits_write_fltarr(Name_Out,SameResMaps);

	if(Timer==True) {
		gettimeofday(&end, (struct timezone *) NULL);
		timersub(&end,&start,&diff);
    	printf("Time for Process_AllBands=%f\n",(double) diff.tv_sec+(double)diff.tv_usec/1000000.);
	}	
}

/*********************************************************************/
void lGMCA_inv::Invert_Weights_Only(char *Name_in,char *Name_inv_vectors,char *Name_inv_struc, char *Name_Out) {

	int kc;
	intarray listchannel;
	int f;
	
	struct timeval start,end,diff;

	if(Timer==True) gettimeofday(&start, (struct timezone *) NULL);

	//Process output names
	
	Read_input_maps(Name_in);//Read input maps
	Read_inversion_structure(Name_inv_struc); //read structure containing the list of channel per resolution band
	
	//Get Lmax per channel
	if(ResolStruc(0) > NchanMax) error(EXIT_FAILURE,"lGMCA->Nchan does not agree",":Process_AllBands()");
	if(ResolStruc(0) > ResolStruc.nx()-3) error(EXIT_FAILURE,"lGMCA->Not enough channels in structure",":Process_AllBands()");
	
	Nchan=ResolStruc(0);
	listchannel.alloc(Nchan+1);
	for (kc=0;kc<=Nchan;kc++) listchannel(kc)=ResolStruc(kc+1)-1;
	Read_inversion_vectors(Name_inv_vectors);
	OutputMap.alloc(NpixRec); //Allocate memory for output
	printf("Copy %d Maps\n",Nchan);
	SameResMaps=InputMaps;
	//Perform wavelet transform and inversion	
	//if(Verbose==true) 

	create_faceindices(Nested_in); 

	printf("Invert each face, NsideRec=%d, Nside=%d\n",NsideRec,Nside);
	#ifdef _OPENMP
    		#pragma omp parallel for private(f) num_threads(outer_loop_threads)
    	#endif
    	for(f=0;f<12;f++) {
    		 WT_invert_single_face(f);
    	}


	fits_write_fltarr(Name_Out,OutputMap);
	
	if(Timer==True) {
		gettimeofday(&end, (struct timezone *) NULL);
		timersub(&end,&start,&diff);
    		printf("Time for Process_AllBands=%f\n",(double) diff.tv_sec+(double)diff.tv_usec/1000000.);
	}	
}


/*********************************************************************/
void lGMCA_inv::Process_AllBands(char *Name_in, char *Name_beams,char *Name_inv_vectors,char *Name_inv_struc, char *Name_Out) {
	int kc,kb;
	intarray listchannel;
	char fname_inv[1024];
	char fname_out[1024];
	char prefix[1024],suffix[1024],prefix_out[1024],suffix_out[1024], * pos_substr;	
	int len_substr,len_str;
	
	struct timeval start,end,diff;

	if(Timer==True) gettimeofday(&start, (struct timezone *) NULL);

	//Process input file names for inversion vectors
	len_str=strlen(Name_inv_vectors);
	pos_substr=strstr(Name_inv_vectors,"_Band0");
	if(pos_substr==NULL) error(EXIT_FAILURE,"lGMCA->Cannot find \"_Band0\" in the name of the inversion vector file provided ",":Process_AllBands()");
	len_substr=pos_substr-Name_inv_vectors;
	strncpy(prefix,Name_inv_vectors,len_substr);
	prefix[len_substr]='\0';
	strncpy(suffix,&Name_inv_vectors[len_substr+6],len_str-len_substr-6);
	suffix[len_str-len_substr-6]='\0';
	if(Debug == true) printf("prefix= %s, suffix = %s\n",prefix,suffix);

	//Process output names
	len_str=strlen(Name_Out);
	pos_substr=strstr(Name_Out,"_Band0");
	if(pos_substr==NULL) error(EXIT_FAILURE,"lGMCA->Cannot find \"_Band0\" in the name of the output file provided ",":Process_AllBands()");
	len_substr=pos_substr-Name_Out;
	strncpy(prefix_out,Name_Out,len_substr);
	prefix_out[len_substr]='\0';
	strncpy(suffix_out,&Name_Out[len_substr+6],len_str-len_substr-6);
	suffix_out[len_str-len_substr-6]='\0';
	if(Debug == true) printf("prefix_out= %s, suffix_out = %s\n",prefix_out,suffix_out);

	Read_input_maps(Name_in);//Read input maps
	Read_Beams(Name_beams);//Read beams
	Read_inversion_structure(Name_inv_struc); //read structure containing the list of channel per resolution band
	
	//Get Lmax per channel
	LmaxChan.alloc(NchanMax);
	LmaxChan.init(0.);
	if(Nband==1) {
		if(ResolStruc(0) > NchanMax) error(EXIT_FAILURE,"lGMCA->Nchan does not agree",":Process_AllBands()");
		if(ResolStruc(0) > ResolStruc.nx()-3) error(EXIT_FAILURE,"lGMCA->Not enough channels in structure",":Process_AllBands()");
		for(kc=0;kc<NchanMax;kc++) {
			LmaxChan(kc)=ResolStruc(ResolStruc.nx()-1);
			AlmX[kc].Set(LmaxChan(kc),LmaxChan(kc));
 			AlmX[kc].SetToZero ();
		}
	} else {
		for(kb=0;kb<Nband;kb++) {
			if(ResolStruc(0,kb) > NchanMax) error(EXIT_FAILURE,"lGMCA->Nchan does not agree",":Process_AllBands()");
			if(ResolStruc(0,kb) > ResolStruc.nx()-3) error(EXIT_FAILURE,"lGMCA->Not enough channels in structure",":Process_AllBands()");
			for(kc=0;kc<ResolStruc(0,kb);kc++) LmaxChan(ResolStruc(kc+2,kb)-1)=max((int) ResolStruc(ResolStruc.nx()-1,kb),LmaxChan(ResolStruc(kc+2,kb)-1)) ;
		}
		for(kc=0;kc<NchanMax;kc++) {
			AlmX[kc].Set(LmaxChan(kc),LmaxChan(kc));
 			AlmX[kc].SetToZero ();
		}
	}	
	if(Debug==true) for(kc=0;kc<NchanMax;kc++) printf("LmaxChan[%d]=%d\n",kc,LmaxChan(kc));
	
	create_nest2ringindices(Nside);
	
	if(Nband==1) {//Single Band
		Nchan=ResolStruc(0);
		listchannel.alloc(Nchan+1);
		for (kc=0;kc<=Nchan;kc++) listchannel(kc)=ResolStruc(kc+1)-1;
		lLmax=ResolStruc(ResolStruc.nx()-1);
		sprintf(fname_inv,"%s_Band%d%s",prefix,0,suffix);
		Read_inversion_vectors(fname_inv);
		OutputMap.alloc(NpixRec); //Allocate memory for output
		Process_Band(0,listchannel);
		if(save_sameres==true) {
			sprintf(fname_out,"%s_Band%d_sameres%s",prefix_out,0,suffix_out);
			fits_write_fltarr(fname_out,SameResMaps);
		}
		sprintf(fname_out,"%s_Band%d_sameres_inverted%s",prefix_out,0,suffix_out);
		fits_write_fltarr(fname_out,OutputMap);
	} else {//More than one band
		for(kb=0;kb<Nband;kb++) {
			Nchan=ResolStruc(0,kb);
			listchannel.alloc(Nchan+1);
			for (kc=0;kc<=Nchan;kc++) {
				listchannel(kc)=ResolStruc(kc+1,kb)-1;
				if(Verbose==true) printf("listchannel[%d]=%d\n",kc,listchannel(kc));
			}
			lLmax=ResolStruc(ResolStruc.nx()-1,kb);
			sprintf(fname_inv,"%s_Band%d%s",prefix,kb,suffix);
			Read_inversion_vectors(fname_inv);
			if(kb==0) OutputMap.alloc(NpixRec);
			if((kb >0)&& (NpixRec!=(unsigned long)OutputMap.nx())) {
				 if(Debug==true) printf("Deallocate memory (%d ->%lu)\n",OutputMap.nx(),NpixRec); //Allocate memory for output
				 OutputMap.free();
				 if(Debug==true) printf("Reallocate memory (%lu)\n",NpixRec); //Allocate memory for output
				 OutputMap.alloc(NpixRec);				 
			}
			Process_Band(kb,listchannel);
			if(save_sameres==true) {
				sprintf(fname_out,"%s_Band%d_sameres%s",prefix_out,kb,suffix_out);
				fits_write_fltarr(fname_out,SameResMaps);
			}
			sprintf(fname_out,"%s_Band%d_sameres_inverted%s",prefix_out,kb,suffix_out);
			fits_write_fltarr(fname_out,OutputMap);
		}
	}
	if(Timer==True) {
		gettimeofday(&end, (struct timezone *) NULL);
		timersub(&end,&start,&diff);
    	printf("Time for Process_AllBands=%f\n",(double) diff.tv_sec+(double)diff.tv_usec/1000000.);
	}	
}

/*********************************************************************/
void lGMCA_inv::Process_Band(int band, intarray &listchannel) {

	int kc,f;
	struct timeval start,end,diff;

	if(Timer==True) gettimeofday(&start, (struct timezone *) NULL);
	
	//Set to same resolution
	if(Verbose==true) printf("Set to same resolution for band %d\n",band);
	for (kc=0;kc<Nchan;kc++) {
		if(Verbose==true) printf("Process channel %d (ref %d)\n",	listchannel(kc+1),listchannel(0));
		Set_Same_Resolution(kc,listchannel(kc+1),listchannel(0),LmaxChan(listchannel(kc+1)));	//Set to same resolution all channels in a Band	
	}
	//if(Debug==true) fits_write_fltarr("SameRes.fits",SameResMaps);
	//Compute face indices if necessary	
	if(Verbose==true) printf("Create face indices if necessary, %d\n",band);
	create_faceindices(true); 

	if(Debug) printf("Max of face indices=%lu\n",FaceIndices.max());

	//Perform wavelet transform and inversion	
	if(Verbose==true) printf("Invert each face, %d\n",band);
	#ifdef _OPENMP
    	#pragma omp parallel for private(f) num_threads(outer_loop_threads)
    #endif
    for(f=0;f<12;f++) {
    	 WT_invert_single_face(f);
    }
    
    if(Timer==True) {
		gettimeofday(&end, (struct timezone *) NULL);
		timersub(&end,&start,&diff);
    	printf("Time for Process_Band Band %d=%f\n",band,(double) diff.tv_sec+(double)diff.tv_usec/1000000.);
	}
    
}

/*********************************************************************/
void lGMCA_inv::Set_Same_Resolution(int index_Array,int CChan, int RefBeam, int ChanLmax) {
	int l,m,iter;
	unsigned long kp,ind;
	dblarray Beameq(ChanLmax+1);
	Hmap<double> Map;
	double avg=0.;
	struct timeval start,end,diff;
	Alm<xcomplex<double> > AlmSameRes;
	
	if(Timer==True) gettimeofday(&start, (struct timezone *) NULL);
	
	if(Beams.nx() < ChanLmax+1) error(EXIT_FAILURE,"lGMCA->Not enough l's in beam file",":Set_Same_Resolution()");
	if(AlmXComp.buffer() == NULL) error(EXIT_FAILURE,"lGMCA->Input Maps not read",":Set_Same_Resolution()");

    if (Debug == True) cout << "==> Use Lmax = " <<   ChanLmax << " NumThread = " << NumThread << endl ;
   		
	//Allocate memory if necessary
	if(NchanMax >1) {
		if((SameResMaps.buffer()==NULL)||(NpixRec!=(unsigned long)SameResMaps.nx())||(SameResMaps.ny()!=Nchan)) SameResMaps.alloc(NpixRec,Nchan);
	} else {
		if((SameResMaps.buffer()==NULL)||(NpixRec!=(unsigned long) SameResMaps.nx())) SameResMaps.alloc(NpixRec);
	}
    if (Debug == True) cout << "==> Use SameResMaps size = " <<    SameResMaps.nx() << " , " << SameResMaps.ny() << endl ;
	
	if(((CChan == RefBeam)||(NchanMax==1))&&(Nside==NsideRec)) {//Nothing to be done, just copy the map
		if(NchanMax >1) {
			#ifdef _OPENMP
				#pragma omp parallel for private(kp) shared(CChan) schedule(static) num_threads(NumThread)  
			#endif
			for(kp=0;kp<Npix;kp++) SameResMaps(kp,index_Array)=InputMaps(kp,CChan);
		} else {
			#ifdef _OPENMP
				#pragma omp parallel for private(kp) schedule(static) num_threads(NumThread)  
			#endif
			for(kp=0;kp<Npix;kp++) SameResMaps(kp)=InputMaps(kp);
		}
		if(Debug) printf("No need to transform the map %d\n",CChan);
		
	} else { // Need to put to same resolution

		//Compute equivalent beam
		#ifdef _OPENMP
			#pragma omp parallel for private(l) schedule(static) num_threads(NumThread)  
		#endif
		for(l=0;l<=ChanLmax;l++) {
			if(Beams(l,CChan) > Prec ) Beameq(l)=Beams(l,RefBeam)/Beams(l,CChan);
			else Beameq(l)=0.;
		}
	
		if(AlmXComp(CChan) !=1 ){//Compute alm Xform	
		
			AlmX[CChan].Set(ChanLmax,ChanLmax);//Here we just go to lLmax
 			AlmX[CChan].SetToZero ();
			
			if (Debug==true) printf("Need to transform chan %d\n",CChan);
			//Set the properties of the transform  
			Set_almtrans(Alm_Niter,Alm_Norm,Alm_Beam,Alm_Fast); 	
			if(weight_Nside!=Nside) Set_Weight(Nside);
		 	Map.alloc(Nside,RING);//Must be in ring format 
			if(Debug==true) printf("Nside=%d\n",Nside);
	
			if(Timer==True) {
				gettimeofday(&end, (struct timezone *) NULL);
				timersub(&end,&start,&diff);
    			printf("Time for Allocate Map Chan %d=%f\n",CChan,(double) diff.tv_sec+(double)diff.tv_usec/1000000.);
			}
	
			//Compute monopole and assign to Healpix Map 
			if(NchanMax >1) {
				#ifdef _OPENMP
					#pragma omp parallel for private(kp) shared(CChan) schedule(static) reduction(+: avg) num_threads(NumThread)  
				#endif
   				for(kp=0;kp<Npix;kp++) avg+=InputMaps(kp,CChan);//Zero freq not estimated thru ALMs
   				avg/=(double) Npix;
				if(Nested_in == false ) { //Direct copy
					#ifdef _OPENMP
						#pragma omp parallel for private(kp) shared(avg,CChan) schedule(static) num_threads(NumThread)  
					#endif 
   					for(kp=0;kp<Npix;kp++) {
   						Map[kp]=InputMaps(kp,CChan)-avg;//Copy input data into a map
   					}
				} else { //Use previously compute indices
					if(Npix == (unsigned long) Nest2Ring_indices.nx()) {
						#ifdef _OPENMP
							#pragma omp parallel for private(kp,ind) shared(avg,CChan) schedule(static) num_threads(NumThread)  
						#endif 
   						for(kp=0;kp<Npix;kp++) {
   							ind=Nest2Ring_indices(kp);
   							Map[ind]=InputMaps(kp,CChan)-avg;//Copy input data into a map
   						} 
					} else error(EXIT_FAILURE,"lGMCA->Not enough pixels in Nest2Ring_indices",":Set_Same_Resolution()");
				}
			} else {
				#ifdef _OPENMP
					#pragma omp parallel for private(kp) reduction(+: avg) schedule(static) num_threads(NumThread)  
				#endif
   				for(kp=0;kp<Npix;kp++) avg+=InputMaps(kp);//Zero freq not estimated thru ALMs
   				avg/=(double) Npix;
   				if(Nested_in == false ) { //Direct copy
					#ifdef _OPENMP
						#pragma omp parallel for private(kp) shared(avg) schedule(static) num_threads(NumThread)  
					#endif 
   					for(kp=0;kp<Npix;kp++) Map[kp]=InputMaps(kp)-avg;//Copy input data into a map
   				} else {
   					if(Npix == (unsigned long) Nest2Ring_indices.nx()) {
						#ifdef _OPENMP
							#pragma omp parallel for private(kp,ind) shared(avg,CChan) schedule(static) num_threads(NumThread)  
						#endif 
   						for(kp=0;kp<Npix;kp++) {
   							ind=Nest2Ring_indices(kp);
   							Map[ind]=InputMaps(kp)-avg;//Copy input data into a map
   						} 
					} else error(EXIT_FAILURE,"lGMCA->Not enough pixels in Nest2Ring_indices",":Set_Same_Resolution()");	
   				}
			}
			if(Timer==True) {
				gettimeofday(&end, (struct timezone *) NULL);
				timersub(&end,&start,&diff);
    			printf("Time for Average subtraction Chan %d=%f\n",CChan,(double) diff.tv_sec+(double)diff.tv_usec/1000000.);
			}
//			if(Nested_in == true) Map.swap_scheme();//Make Sure we are in RING order -> required for alm xform map2alm
			if(Timer==True) {
				gettimeofday(&end, (struct timezone *) NULL);
				timersub(&end,&start,&diff);
    			printf("Time for swapping scheme Chan %d=%f\n",CChan,(double) diff.tv_sec+(double)diff.tv_usec/1000000.);
			}


			if(Debug==true) printf("Num_Alms=%ld, InputMaps.nx()=%d,Npix=%lu\n",AlmX[CChan].Num_Alms(ChanLmax,ChanLmax),InputMaps.nx(),Npix); 
			//Get estimate of ALMS
			map2alm(Map,AlmX[CChan],weight_TRing);
			if(Timer==True) {
				gettimeofday(&end, (struct timezone *) NULL);
				timersub(&end,&start,&diff);
    			printf("Time for alm Xform strictly Chan %d=%f\n",CChan,(double) diff.tv_sec+(double)diff.tv_usec/1000000.);
			}
			
		
			if(Alm_Niter>0) { //Case iterative recons
				if(Verbose==true) printf("Iterative ALM Xform (iter=%d)\n",Alm_Niter); 
				Hmap<double> Map2;
				Map2.alloc(Nside,RING);
				for (iter=1; iter<=Alm_Niter; ++iter)
    			{
    				alm2map(AlmX[CChan],Map2);
					#ifdef _OPENMP
						#pragma omp parallel for default(none) private(kp) shared(Map,Map2) num_threads(NumThread) 
    				#endif
    				for (kp=0; kp<Npix; ++kp) Map2[kp] = Map[kp]-Map2[kp];
    				map2alm(Map2,AlmX[CChan],weight_TRing,true);
    			}
			}
		   	AlmX[CChan](0,0) = xcomplex<REAL> (avg*sqrt(fourpi), 0.);//Zero freq not estimated thru ALMs
		   	AlmXComp(CChan)=1;
		} else if (Debug==true) printf("No need to transform again chan %d\n",CChan);
		
		if(Timer==True) {
			gettimeofday(&end, (struct timezone *) NULL);
			timersub(&end,&start,&diff);
    		printf("Time for ALM Xform Chan %d=%f\n",CChan,(double) diff.tv_sec+(double)diff.tv_usec/1000000.);
		}
		
		if((AlmX[CChan].Lmax() < lLmax)||(AlmX[CChan].Mmax() < lLmax)) {
		 	  printf("Lmax=%d, Mmax=%d, lLmax=%d\n",AlmX[CChan].Lmax(),AlmX[CChan].Mmax(),lLmax);
			  error(EXIT_FAILURE,"lGMCA->Not enough memory allocated for Alm transform",":Set_Same_Resolution()");
		}
  		//Multiply in alm space with equivalent beam
		Set_almtrans(Alm_Niter,Alm_Norm,Alm_Beam,Alm_Fast); 
		AlmSameRes.Set(lLmax,lLmax);//Here we just go to lLmax
 		AlmSameRes.SetToZero ();	
 		
 		if(Debug) printf("Lmax=%d, Mmax=%d, lLmax=%d, Beameq.nx()=%d\n",AlmX[CChan].Lmax(),AlmX[CChan].Mmax(),lLmax,Beameq.nx());

		#ifdef _OPENMP
			#pragma omp parallel for default(none) private(m,l) shared(CChan,Beameq,AlmSameRes) num_threads(NumThread) 
		#endif
		for(m=0;m<=lLmax;m++) {
   	   		for(l=m;l<=lLmax;l++) {
       			// AlmSameRes(l,m).real()= AlmX[CChan](l,m).real()*Beameq(l);
       			// AlmSameRes(l,m).imag()= AlmX[CChan](l,m).imag()*Beameq(l);
                AlmSameRes(l,m) *= (REAL) Beameq(l);
   	   		}
    	}
    	avg=real( AlmSameRes(0,0)) /sqrt(fourpi);
		AlmSameRes(0,0) = xcomplex<REAL> (0.,0.);
		
		if(Debug) printf("Set everything to NsideRec=%d\n",NsideRec);
		
		//Reconstruct
    	Map.SetNside(NsideRec, RING); //alm2map outputs a ring map
		if(weight_Nside!=NsideRec) Set_Weight(NsideRec);
		if(Debug) printf("size of Map= %d / size of SameresMaps=%d\n",Map.Npix(),SameResMaps.nx());
		alm2map(AlmSameRes,Map);//output ring map
    	
    	if(Timer==True) {
			gettimeofday(&end, (struct timezone *) NULL);
			timersub(&end,&start,&diff);
    		printf("Time for ALM recons Chan %d=%f\n",CChan,(double) diff.tv_sec+(double)diff.tv_usec/1000000.);
		}
		
		create_ring2nestindices(NsideRec); //Make sure we go back to RING before computing WT
    	
//		if(Nested_in == true)  Map.swap_scheme();//Make Sure we are back to right order
		
		if(Timer==True) {
			gettimeofday(&end, (struct timezone *) NULL);
			timersub(&end,&start,&diff);
    		printf("Time for swapping back Chan %d=%f\n",CChan,(double) diff.tv_sec+(double)diff.tv_usec/1000000.);
		}
		
		//Assign result to output array
		#ifdef _OPENMP
			#pragma omp parallel for private(kp,ind) shared(avg,CChan) num_threads(NumThread) 
		#endif
		for(kp=0;kp<NpixRec;kp++) {
			 ind=Ring2Nest_indices(kp);
			 SameResMaps(ind,index_Array)=(float) Map[kp]+avg;//ORIGINAL FORMAT
		}
		if(Timer==True) {
			gettimeofday(&end, (struct timezone *) NULL);
			timersub(&end,&start,&diff);
    		printf("Time for Set_Same_Resolution Chan %d=%f\n",CChan,(double) diff.tv_sec+(double)diff.tv_usec/1000000.);
		}
		
	}
		
}

/*********************************************************************/
void lGMCA_inv::create_faceindices(bool cNest) { 
	int kx,ky,kf;

	struct timeval start,end,diff;
	
	if(Timer==True) gettimeofday(&start, (struct timezone *) NULL);
	if (NsideRec > 0) {
		if(NsideFace != NsideRec) {
			Hmap<REAL> Map;
			Map.alloc(NsideRec,cNest);
			FaceIndices.alloc(NsideRec,NsideRec,12);
			if(cNest==true) { //NESTED INPUT FORMAT
				unsigned long offset[12];
				unsigned long NsideRec2=(unsigned long) NsideRec * (unsigned long) NsideRec;
				for(kf=0;kf<12;kf++) offset[kf]= NsideRec2* (unsigned long) kf;
				#ifdef _OPENMP
					#pragma omp parallel for private(kx,ky,kf) shared(offset) num_threads(NumThread)
				#endif
				for (kx=0;kx<NsideRec;kx++){
					for (ky=0;ky<NsideRec;ky++) {
						FaceIndices(kx,ky,0)= (unsigned long) Map.xyf2ind (kx,ky,0); 
						for (kf=1;kf<12;kf++) FaceIndices(kx,ky,kf)=FaceIndices(kx,ky,0)+offset[kf];
					}
				}
			} else {
				for (kf=0;kf<12;kf++){
					#ifdef _OPENMP
						#pragma omp parallel for private(kx,ky) shared(kf) num_threads(NumThread)
					#endif
					for (kx=0;kx<Nside;kx++) {
						for (ky=0;ky<Nside;ky++) FaceIndices(kx,ky,kf)= (unsigned long) Map.xyf2ind (kx,ky,kf);	
					}
				}
			}
			NsideFace=NsideRec;
		} 
	} else {
		error(EXIT_FAILURE,"lGMCA->Nside not strictly greater than 0 ",":create_faceindices()");
	}
	if(Timer==True) {
		gettimeofday(&end, (struct timezone *) NULL);
		timersub(&end,&start,&diff);
    	printf("Time for create_faceindices=%f\n",(double) diff.tv_sec+(double)diff.tv_usec/1000000.);
	}	
}

/*********************************************************************/
void lGMCA_inv::create_nest2ringindices(int lNside) { 
	unsigned long kp,lNpix;
	Hmap<REAL> Map;
	Map.alloc(lNside,true);
	struct timeval start,end,diff;
	
	if(Timer==True) gettimeofday(&start, (struct timezone *) NULL);

	
	lNpix=(unsigned long)lNside*(unsigned long)lNside*12ul;
	
	if(lNpix!= (unsigned long) Nest2Ring_indices.nx()) {
		Nest2Ring_indices.alloc(lNpix);
		#ifdef _OPENMP
			#pragma omp parallel for private(kp) num_threads(NumThread)
		#endif
		for(kp=0;kp<lNpix;kp++) {
			Nest2Ring_indices(kp)=Map.nest2ring(kp);
		}

		if(Timer==True) {
			gettimeofday(&end, (struct timezone *) NULL);
			timersub(&end,&start,&diff);
    		printf("Time for create_nest2ringindices=%f\n",(double) diff.tv_sec+(double)diff.tv_usec/1000000.);
		}
	}	
}

/*********************************************************************/
void lGMCA_inv::create_ring2nestindices(int lNside) { 
	unsigned long kp,lNpix;
	Hmap<REAL> Map;
	Map.alloc(lNside,true);
	struct timeval start,end,diff;
	
	if(Timer==True) gettimeofday(&start, (struct timezone *) NULL);

	
	lNpix=(unsigned long)lNside*(unsigned long)lNside*12ul;
	
	if(lNpix!= (unsigned long) Ring2Nest_indices.nx()) {
		Ring2Nest_indices.alloc(lNpix);
		#ifdef _OPENMP
			#pragma omp parallel for private(kp) num_threads(NumThread)
		#endif
		for(kp=0;kp<lNpix;kp++) {
			Ring2Nest_indices(kp)=Map.ring2nest(kp);
		}

		if(Timer==True) {
			gettimeofday(&end, (struct timezone *) NULL);
			timersub(&end,&start,&diff);
    		printf("Time for create_ring2nestindices=%f\n",(double) diff.tv_sec+(double)diff.tv_usec/1000000.);
		}
	}	
}

/*********************************************************************/
void lGMCA_inv::WT_invert_single_face(int nf) {
	int kc,kx,ky,kp,ind;
	Ifloat Face,*AllChanFaces;
	struct timeval start,end,diff;
	if(Timer==True) gettimeofday(&start, (struct timezone *) NULL);

	Face.alloc(NsideRec,NsideRec);
	AllChanFaces=new Ifloat[Nscale*Nchan];
	for (kp=0;kp<Nscale*Nchan;kp++) AllChanFaces[kp].alloc(Face.ny(), Face.nx());
	
	if(Nchan == 1) { //Case of a single channel: just transform + weighting // in that case, InversionVectors are 2D
		#ifdef _OPENMP
			#pragma omp parallel for default(none) private(kx,ky) shared(Face,nf) num_threads(inner_loop_threads)
		#endif
		for(kx=0;kx<NsideRec;kx++) {
			for(ky=0;ky<NsideRec;ky++) Face(ky,kx)=SameResMaps( FaceIndices(kx,ky,nf) ) ;
		}
		(WT_Trans.UWT2D)->transform(Face,AllChanFaces,Nscale);
		for(kp=0;kp<Nscale;kp++) {
			#ifdef _OPENMP
				#pragma omp parallel for private(kx,ky,ind) shared(Face,kp,AllChanFaces,nf) num_threads(inner_loop_threads)
			#endif
			for(kx=0;kx<NsideRec;kx++) {
				for(ky=0;ky<NsideRec;ky++) {
					ind= FaceIndices(kx,ky,nf);
					(AllChanFaces[kp])(ky,kx)*=InversionVectors(ind,kp);
				}
			}
		}
		(WT_Trans.UWT2D)->recons(AllChanFaces,Face,Nscale,True); //perform reconstruction of face	
	} else {
		for(kc=0;kc<Nchan;kc++) {
			#ifdef _OPENMP
				#pragma omp parallel for default(none) private(kx,ky) shared(Face,kc,nf) num_threads(inner_loop_threads)
			#endif
			for(kx=0;kx<NsideRec;kx++) {
				for(ky=0;ky<NsideRec;ky++) Face(ky,kx)=SameResMaps( FaceIndices(kx,ky,nf),kc);
			}
			(WT_Trans.UWT2D)->transform(Face,&(AllChanFaces[kc*Nscale]),Nscale);
			for(kp=0;kp<Nscale;kp++) {
				#ifdef _OPENMP
					#pragma omp parallel for private(kx,ky,ind) shared(Face,kp,kc,AllChanFaces,nf) num_threads(inner_loop_threads)
				#endif
				for(kx=0;kx<NsideRec;kx++) {
					for(ky=0;ky<NsideRec;ky++) {
						ind= FaceIndices(kx,ky,nf);
						if(kc ==0 ) AllChanFaces[kp](ky,kx)*=InversionVectors(ind,kp,0);
						else AllChanFaces[kp](ky,kx)+=(AllChanFaces[kc*Nscale+kp])(ky,kx)*InversionVectors(ind,kp,kc);
					}
				}
			}
		}
		(WT_Trans.UWT2D)->recons(AllChanFaces,Face,Nscale,True); //perform reconstruction of face
		delete [] AllChanFaces;
	}
	#ifdef _OPENMP
		#pragma omp parallel for default(none) private(kx,ky) shared(Face,nf) num_threads(inner_loop_threads)
	#endif
	for(kx=0;kx<NsideRec;kx++) {
		for(ky=0;ky<NsideRec;ky++) {
			OutputMap(FaceIndices(kx,ky,nf))=Face(ky,kx);
		}
	}						
	Face.free();
	if(Timer==True) {
		gettimeofday(&end, (struct timezone *) NULL);
		timersub(&end,&start,&diff);
    	printf("Time for WT_invert_single_face face %d=%f\n",nf,(double) diff.tv_sec+(double)diff.tv_usec/1000000.);
	}
};

/*********************************************************************/
//Not Used Wavelet transform routines
void lGMCA_inv::WT_transform_single_face(fltarray &OutMap,int nf) {
	int kc,kx,ky,kp;
	Ifloat Face,*AllChanFaces;

	Face.alloc(Nside,Nside);
	AllChanFaces=new Ifloat[Nscale*Nchan];
	for (kp=0;kp<Nscale*Nchan;kp++) AllChanFaces[kp].alloc(Face.ny(), Face.nx());
	
	if(Nchan == 1) {
		#ifdef _OPENMP
			#pragma omp parallel for default(none) private(kx,ky) shared(Face,nf) num_threads(inner_loop_threads)
		#endif
		for(kx=0;kx<Nside;kx++) {
			for(ky=0;ky<Nside;ky++) Face(ky,kx)=InputMaps( FaceIndices(kx,ky,nf) ) ;
		}
		(WT_Trans.UWT2D)->transform(Face,AllChanFaces,Nscale);
	} else {
		for(kc=0;kc<Nchan;kc++) {
			#ifdef _OPENMP
				#pragma omp parallel for default(none) private(kx,ky) shared(Face,kc,nf) num_threads(inner_loop_threads)
			#endif
			for(kx=0;kx<Nside;kx++) {
				for(ky=0;ky<Nside;ky++) Face(ky,kx)=InputMaps( FaceIndices(kx,ky,nf),kc);
			}
			(WT_Trans.UWT2D)->transform(Face,&(AllChanFaces[kc*Nscale]),Nscale);	
		}
	}
	for(kp=0;kp<Nscale;kp++) {
		#ifdef _OPENMP
			#pragma omp parallel for private(kx,ky) shared(Face,kp,AllChanFaces,nf) num_threads(inner_loop_threads)
		#endif
		for(kx=0;kx<Nside;kx++) {
			for(ky=0;ky<Nside;ky++) 
				OutMap(FaceIndices(kx,ky,nf),kp)=(AllChanFaces[kp])(ky,kx);
		}
	}				
	delete [] AllChanFaces;
	Face.free();
};


void lGMCA_inv::WT_transform_single_map(int Chan) {
	unsigned long kul;
	Hmap<REAL> Map;
	
	if (Chan < Nchan) {
		Map.alloc(Nside,Nested_in);
		WT_Trans.alloc(Nside, Nscale,Nested_in,0); //default values
		if(InputMaps.naxis()>1) for (kul=0;kul< Npix;kul++) Map[kul]=InputMaps(kul,Chan);//copy value
		else for (kul=0;kul< Npix;kul++) Map[kul]=InputMaps(kul);//copy value
		WT_Trans.transform(Map);		
	} else {
		error(EXIT_FAILURE,"lGMCA->Chan greater than Nchan",":WT_transform_single_map()");
	}	
};





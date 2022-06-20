#include "mrsp_lib_matmask.cc"

//****************************************************************************************************************************//
//****************************************************************************************************************************//
//****************************************************************************************************************************//
//IDL/C++ Wrappers
//Main Philosophy: Avoid copying pointer elements by elements as much as possible
//Use directly the Healpix arr class to perform so.
//NOTE: ALL INPUT/OUTPUT MEMORY SHOULD BE PREVIOUSLY ALLOCATED IN IDL
//****************************************************************************************************************************//
//****************************************************************************************************************************//
//****************************************************************************************************************************//




//****************************************************************************************************************************//
//****************************************************************************************************************************//
//SET INPUT OR INTERMEDIARY DATA ROUTINES
//****************************************************************************************************************************//
//****************************************************************************************************************************//
template <class T>
void MasterPola <T>::idl_set_map_TQU(T *cMap, long Npix, bool cflag_tqu,bool onlyT, bool cflag_nested, bool PolaFastALM ) {
	Healpix_Ordering_Scheme NFlag;
	if(cflag_nested) NFlag=NEST; 
	else NFlag=RING;
	long NsideIn=sqrt(Npix/12l);
	if(Npix!= NsideIn*NsideIn*12l) printf("Input Map size not compatible with Healpix Map (%ld vs %ld for Nside = %ld)\n. No Map allocated\n", Npix,NsideIn*NsideIn*12l,NsideIn);
	else {
		arr<T> MapT(cMap, Npix);
		if(onlyT) {
			this->Map_TQU.map_T.Set(MapT, NFlag);	
			this->Nside_map = Map_TQU.map_T.Nside();
			if(Verbose) printf("Nside Map onlyT:%d\n", this->Nside_map);			
			this->Nmaps=1;
		} else {
			long offset=Npix;
			arr<T> MapQ(&cMap[offset], Npix);
			offset+=Npix;
			arr<T> MapU(&cMap[offset], Npix);
			this->Map_TQU.alloc(MapT,MapQ,MapU, !cflag_tqu, cflag_nested);
			this->Nside_map = this->Map_TQU.map_T.Nside();
			this->Nmaps=3;
			if((!Flag_Work_EB) && (cflag_tqu == false))  this->Map_TQU.swap_tqu_teb(PolaFastALM);
			if(( Flag_Work_EB) && (cflag_tqu == true))  this->Map_TQU.swap_tqu_teb(PolaFastALM);
		}
	}
}


//****************************************************************************************************************************//
template <class T>
void MasterPola <T>::idl_set_Noisemap_TQU(T *cNoiseMap, long Npix, bool cflag_tqu,bool onlyT, bool cflag_nested, bool PolaFastALM ) {
	Healpix_Ordering_Scheme NFlag;
	if(cflag_nested) NFlag=NEST; 
	else NFlag=RING;
	long NsideIn=sqrt(Npix/12l);
	if(NsideIn != this->Nside_map) printf("BEWARE MISMATCH IN NSIDE FROM MAPS (%ld) AND NOISE (%ld)\n",(long) Nside_map, NsideIn);
	if(Npix!= NsideIn*NsideIn*12l) printf("Input Noise Map size not compatible with Healpix Map (%ld vs %ld for Nside = %ld)\n. No Map allocated\n", Npix,NsideIn*NsideIn*12l,NsideIn);
	else {
		arr<T> NoiseMapT(cNoiseMap, Npix);
		if(onlyT) this->NoiseMap_TQU.map_T.Set(NoiseMapT, NFlag);	
		else {
			long offset=Npix;
			arr<T> NoiseMapQ(&cNoiseMap[offset], Npix);
			offset+=Npix;
			arr<T> NoiseMapU(&cNoiseMap[offset], Npix);
			this->NoiseMap_TQU.alloc(NoiseMapT, NoiseMapQ, NoiseMapU, !cflag_tqu, cflag_nested);
			if((!Flag_Work_EB) && (cflag_tqu == false))  this->NoiseMap_TQU.swap_tqu_teb(PolaFastALM);
			if(( Flag_Work_EB) && (cflag_tqu == true)) this->NoiseMap_TQU.swap_tqu_teb(PolaFastALM);
		}
	}
}

//****************************************************************************************************************************//
template <class T>
void MasterPola <T>:: idl_set_MaskT(T *MaskT, long Npix, bool cflag_nested, bool set_ring) {
	arr<T> MMT(MaskT,Npix);
	Healpix_Ordering_Scheme NFlag;
	if(cflag_nested) NFlag=NEST; 
	else NFlag=RING;
	this->MaskT.Set(MMT, NFlag); 
	if((set_ring == true) && (NFlag ==NEST)) this->MaskT.swap_scheme();
	this->MaskTFlag=true;
}

//****************************************************************************************************************************//
template <class T>
void MasterPola <T>::idl_set_MaskP(T *MaskP, long Npix, bool cflag_nested, bool set_ring) {
	arr<T> MMP(MaskP,Npix);
	Healpix_Ordering_Scheme NFlag;
	if(cflag_nested) NFlag=NEST; 
	else NFlag=RING;
	this->MaskP.Set(MMP, NFlag); 
	if((set_ring == true) && (NFlag ==NEST)) this->MaskP.swap_scheme();
	this->MaskPFlag=true;
}

//****************************************************************************************************************************//
template <class T>
void MasterPola <T>::idl_set_spec_TEB(double* PSPEC_MSKMAP, bool onlyT) { 
	long Nls= (long) this->lLmax +1;
	arr<double> PspecTT(PSPEC_MSKMAP, Nls);
	if(onlyT) {
		this->Nmaps=1;
		this->MskMapPowSpec.Set(PspecTT);
		printf("Nmaps=%d, %d\n",Nmaps,this->MskMapPowSpec.num_specs);
	} else {
		this->Nmaps=3;
		long offset= Nls;
		arr<double> PspecEE(&PSPEC_MSKMAP[offset], Nls);
		offset+= Nls;
		arr<double> PspecBB(&PSPEC_MSKMAP[offset], Nls);
		offset+= Nls;
		arr<double> PspecTE(&PSPEC_MSKMAP[offset], Nls);
		offset+= Nls;
		arr<double> PspecTB(&PSPEC_MSKMAP[offset], Nls);
		offset+= Nls;
		arr<double> PspecEB(&PSPEC_MSKMAP[offset], Nls);
		this->MskMapPowSpec.Set(PspecTT, PspecEE, PspecBB, PspecTE, PspecTB, PspecEB);
		printf("Nmaps=%d, %d\n",Nmaps,this->MskMapPowSpec.num_specs);
	}
	this->PMapFlag=true;						
}


//****************************************************************************************************************************//
template <class T>
void MasterPola <T>::idl_set_spec_NoiseTEB(double* PSPEC_MSKNOISE, bool onlyT) { 
	long Nls= (long) this->lLmax +1;
	arr<double> PspecTT(PSPEC_MSKNOISE, Nls);
	if(onlyT) this->MskNoisePowSpec.Set(PspecTT);
	else {
		long offset= Nls;
		arr<double> PspecEE(& PSPEC_MSKNOISE[offset], Nls);
		offset+= Nls;
		arr<double> PspecBB(& PSPEC_MSKNOISE[offset], Nls);
		offset+= Nls;
		arr<double> PspecTE(& PSPEC_MSKNOISE[offset], Nls);
		offset+= Nls;
		arr<double> PspecTB(& PSPEC_MSKNOISE[offset], Nls);
		offset+= Nls;
		arr<double> PspecEB(& PSPEC_MSKNOISE[offset], Nls);
		this->MskNoisePowSpec.Set(PspecTT, PspecEE, PspecBB, PspecTE, PspecTB, PspecEB);
	}
	this->PNoiseMapFlag =true;						
}

//****************************************************************************************************************************//
template <class T>
void MasterPola <T>::idl_set_spec_mask(double* PSPEC_MASKS, bool onlyT) { 
	long Nls= (long) this->lLmax +1;
	arr<double> PspecTT(PSPEC_MASKS, Nls);
	this->MaskTT_Powspec.transfer(PspecTT);
	if(!onlyT) {
		long offset=Nls;
		arr<double> PspecPP(&PSPEC_MASKS[offset], Nls);
		this->MaskPP_Powspec.transfer(PspecPP);
		offset+=Nls;
		arr<double> PspecTP(&PSPEC_MASKS[offset], Nls);		
		this->MaskTP_Powspec.transfer(PspecTP);
	} 
	PMaskFlag =true;						
}

//****************************************************************************************************************************//
template <class T>
void MasterPola <T>::idl_set_coupling_mat(double* KMat, bool onlyT) { 
	long Nls= (long) this->lLmax +1, Nls2=Nls*Nls;
	long offset,offset2=0;

	this->MAT_TT_TT.alloc(Nls,Nls);
	for(long kl2=0;kl2<	Nls;kl2++) {
		offset= Nls*kl2;
		for(long kl1=0;kl1<Nls;kl1++) this->MAT_TT_TT(kl1,kl2)= KMat[kl1+ offset];
	}
	if(!onlyT) {
		offset2+= Nls2;
		this->MAT_EE_EE.alloc(Nls,Nls);
		for(long kl2=0;kl2<	Nls;kl2++) {
			offset= Nls*kl2+offset2;
			for(long kl1=0;kl1<Nls;kl1++) this->MAT_EE_EE(kl1,kl2)= KMat[kl1+ offset];
		}
		offset2+= Nls2;
		this->MAT_EE_BB.alloc(Nls,Nls);
		for(long kl2=0;kl2<	Nls;kl2++) {
			offset= Nls*kl2+offset2;
			for(long kl1=0;kl1<Nls;kl1++) this->MAT_EE_BB(kl1,kl2)= KMat[kl1+ offset];
		}
		offset2+= Nls2;
		this->MAT_TE_TE.alloc(Nls,Nls);
		for(long kl2=0;kl2<	Nls;kl2++) {
			offset= Nls*kl2+offset2;
			for(long kl1=0;kl1<Nls;kl1++) this->MAT_TE_TE(kl1,kl2)= KMat[kl1+ offset];
		}
		offset2+= Nls2;
		this->MAT_EB_EB.alloc(Nls,Nls);
		for(long kl2=0;kl2<	Nls;kl2++) {
			offset= Nls*kl2+offset2;
			for(long kl1=0;kl1<Nls;kl1++) this->MAT_EB_EB(kl1,kl2)= KMat[kl1+ offset];
		}
	}
}

//****************************************************************************************************************************//
template <class T>
void MasterPola <T>::idl_set_spec_radii(double* SpecRad) { 
	this->gamma_TT_TT=SpecRad[0];
	if(Nmaps==3) {
		this->gamma_EE_EE=SpecRad[1];
		this->gamma_EE_BB =SpecRad[2];
		this->gamma_TE_TE =SpecRad[3];
		this->gamma_EB_EB =SpecRad[4];
		this->gamma_EB_4BLOCKS =SpecRad[5];
	}	
}

//****************************************************************************************************************************//
template <class T>
void MasterPola <T>::idl_set_inv_coupling_mat(double* InvMat) { 
	long Nls= (long) this->lLmax+1, Nls2=Nls*Nls;
	long offset,offset2=0,offset3;

	this->inv_MAT_TT_TT.alloc(Nls,Nls);
	for(long kl2=0;kl2<	Nls;kl2++) {
		offset= Nls*kl2;
		for(long kl1=0;kl1<Nls;kl1++) this->inv_MAT_TT_TT(kl1,kl2)= InvMat[kl1+ offset];
	}
	offset2+= Nls2;
	this->inv_MAT_TE_TE.alloc(Nls,Nls);
	for(long kl2=0;kl2<	Nls;kl2++) {
		offset= Nls*kl2+offset2;
		for(long kl1=0;kl1<Nls;kl1++) this->inv_MAT_TE_TE(kl1,kl2)= InvMat[kl1+ offset];
	}
	offset2+= Nls2;
	this->inv_MAT_EB_EB.alloc(Nls,Nls);
	for(long kl2=0;kl2<	Nls;kl2++) {
		offset= Nls*kl2+offset2;
		for(long kl1=0;kl1<Nls;kl1++) this->inv_MAT_EB_EB(kl1,kl2)= InvMat[kl1+ offset];
	}
	offset2+= Nls2;
	this->inv_MAT_EB_4Blocks.alloc(2*Nls,2*Nls);
	for(long kl2=0;kl2<	Nls;kl2++) {
		offset= Nls*kl2+offset2;
		offset3=offset+Nls2;
		for(long kl1=0;kl1<Nls;kl1++) { //Assign all 4 blocks
			this->inv_MAT_EB_4Blocks(kl1,kl2)= InvMat[kl1+ offset];
			this->inv_MAT_EB_4Blocks(kl1+Nls,kl2+Nls)= InvMat[kl1+ offset];
			this->inv_MAT_EB_4Blocks(kl1+Nls,kl2)= InvMat[kl1+ offset3];
			this->inv_MAT_EB_4Blocks(kl1,kl2+Nls)= InvMat[kl1+ offset3];
		}
	}
}



//****************************************************************************************************************************//
//****************************************************************************************************************************//
//GET INTERMEDIARY AND FINAL RESULTS
//****************************************************************************************************************************//
//****************************************************************************************************************************//
template <class T>
void MasterPola <T>::idl_get_TEBmap(T *EBMap, bool PolaFastALM) { 
	if(this->Nmaps==3) {
		PolaHmap<T> Map_Save;
		Map_Save.alloc(this->Map_TQU.get_nside(), this->Map_TQU.flag_teb(), this->Map_TQU.flag_nested());
		Map_Save.set_map_T(this->Map_TQU.map_T);
		Map_Save.set_map_Q(this->Map_TQU.map_Q);
		Map_Save.set_map_U(this->Map_TQU.map_U);
		if(! Flag_Work_EB) Map_Save.swap_tqu_teb(PolaFastALM);

		for(long kl=0;kl<Map_Save.map_T.Npix();kl++) EBMap[kl]=Map_Save.map_T[kl];
		long offset=Map_Save.map_T.Npix();
		for(long kl=0;kl<Map_Save.map_Q.Npix();kl++) EBMap[kl+ offset]=Map_Save.map_Q[kl];
		offset+=Map_Save.map_Q.Npix();
		for(long kl=0;kl<Map_Save.map_U.Npix();kl++) EBMap[kl+ offset]=Map_Save.map_U[kl];
	} else printf("Only Temperature : no need to save E/B maps\n");			
}

//****************************************************************************************************************************//
template <class T>
void MasterPola <T>::idl_get_spec_masks(T *PSPEC_MASKS, bool onlyT) { 
	long Nls= this->lLmax + 1l ;
	long offset= 0;
	if(this->PMaskFlag) {
		for(long kl=0;kl<Nls;kl++) PSPEC_MASKS[kl]= this->MaskTT_Powspec[kl];
		if(!onlyT) {
		printf("POLA MASK ALSO\n");
			offset+=Nls;
			for(long kl=0;kl<Nls;kl++) PSPEC_MASKS[kl+offset]= this->MaskPP_Powspec[kl];
			offset+=Nls;
			for(long kl=0;kl<Nls;kl++) PSPEC_MASKS[kl+offset]= this->MaskTP_Powspec[kl];
		}
	}
}

//****************************************************************************************************************************//
template <class T>
void MasterPola <T>::idl_get_spec_mskmap(T *PSPEC_MSKMAP) { 
	long Nls= this->lLmax + 1l ;
	if(this->PMapFlag) {
		long offset= 0;
		if( this->MskMapPowSpec.num_specs >0 ) {
			for(long kl=0;kl<Nls;kl++) PSPEC_MSKMAP[kl]=this-> MskMapPowSpec.tt_[kl];
			if( this->MskMapPowSpec.num_specs == 6 ) {
			printf("SAVE ALSO POLARIZED MASKED MAPS\n");
				offset+=Nls;
				for(long kl=0;kl<Nls;kl++) PSPEC_MSKMAP[kl+offset]= this->MskMapPowSpec.ee_[kl];
				offset+=Nls;
				for(long kl=0;kl<Nls;kl++) PSPEC_MSKMAP[kl+offset]= this->MskMapPowSpec.bb_[kl];
				offset+=Nls;
				for(long kl=0;kl<Nls;kl++) PSPEC_MSKMAP[kl+offset]= this->MskMapPowSpec.te_[kl];
				offset+=Nls;
				for(long kl=0;kl<Nls;kl++) PSPEC_MSKMAP[kl+offset]= this->MskMapPowSpec.tb_[kl];
				offset+=Nls;
				for(long kl=0;kl<Nls;kl++) PSPEC_MSKMAP[kl+offset]= this->MskMapPowSpec.eb_[kl];
			}
		}
	}
}

//****************************************************************************************************************************//
template <class T>
void MasterPola <T>::idl_get_spec_msknoisemap(T *PSPEC_MSKNOISE) { 
	long Nls= this->lLmax + 1l ;
	if(this-> PNoiseMapFlag) {
		long offset= 0;
		if( this->MskNoisePowSpec.num_specs >0 ) {
			for(long kl=0;kl<Nls;kl++) PSPEC_MSKNOISE[kl]= this->MskNoisePowSpec.tt_[kl];
			if( MskNoisePowSpec.num_specs == 6 ) {
				offset+=Nls;
				for(long kl=0;kl<Nls;kl++) PSPEC_MSKNOISE[kl+offset]= this->MskNoisePowSpec.ee_[kl];
				offset+=Nls;
				for(long kl=0;kl<Nls;kl++) PSPEC_MSKNOISE[kl+offset]= this->MskNoisePowSpec.bb_[kl];
				offset+=Nls;
				for(long kl=0;kl<Nls;kl++) PSPEC_MSKNOISE[kl+offset]= this->MskNoisePowSpec.te_[kl];
				offset+=Nls;
				for(long kl=0;kl<Nls;kl++) PSPEC_MSKNOISE[kl+offset]= this->MskNoisePowSpec.tb_[kl];
				offset+=Nls;
				for(long kl=0;kl<Nls;kl++) PSPEC_MSKNOISE[kl+offset]= this->MskNoisePowSpec.eb_[kl];
			}
		}
	}
}

//****************************************************************************************************************************//
template <class T>
void MasterPola <T>::idl_get_coupling_mat(double* KMat) { 
	long naxes[2],offset,offset2=0;
  	naxes[0] = (long) this->MAT_TT_TT.axis(1);   
   	naxes[1] = (long) this->MAT_TT_TT.axis(2);
	for(long kl2=0;kl2<	naxes[1];kl2++) {
		offset=naxes[0]*kl2;
		for(long kl1=0;kl1<naxes[0];kl1++) KMat[kl1+ offset]=this->MAT_TT_TT(kl1,kl2);
	}
	if(Nmaps==3) {
		offset2+=naxes[0]*naxes[1];
	   	naxes[0] = (long) this->MAT_EE_EE.axis(1);  
	   	naxes[1] = (long) this->MAT_EE_EE.axis(2);
		for(long kl2=0;kl2<	naxes[1];kl2++) {
			offset=naxes[0]*kl2+offset2;
			for(long kl1=0;kl1<naxes[0];kl1++) KMat[kl1+ offset]=this->MAT_EE_EE(kl1,kl2);
		}
		offset2+=naxes[0]*naxes[1];
	   	naxes[0] = (long) this->MAT_EE_BB.axis(1); 
	   	naxes[1] = (long) this->MAT_EE_BB.axis(2);
		for(long kl2=0;kl2<	naxes[1];kl2++) {
			offset=naxes[0]*kl2+offset2;
			for(long kl1=0;kl1<naxes[0];kl1++) KMat[kl1+ offset]=this->MAT_EE_BB(kl1,kl2);
		}
		offset2+=naxes[0]*naxes[1];
	   	naxes[0] = (long) this->MAT_TE_TE.axis(1);   
	   	naxes[1] = (long) this->MAT_TE_TE.axis(2);
		for(long kl2=0;kl2<	naxes[1];kl2++) {
			offset=naxes[0]*kl2+offset2;
			for(long kl1=0;kl1<naxes[0];kl1++) KMat[kl1+ offset]=this->MAT_TE_TE(kl1,kl2);
		}
		offset2+=naxes[0]*naxes[1];
   		naxes[0] = (long) this->MAT_EB_EB.axis(1); 
   		naxes[1] = (long) this->MAT_EB_EB.axis(2);
		for(long kl2=0;kl2<	naxes[1];kl2++) {
			offset=naxes[0]*kl2+offset2;
			for(long kl1=0;kl1<naxes[0];kl1++) KMat[kl1+ offset]=this->MAT_EB_EB(kl1,kl2);
		}
	}	
}

//****************************************************************************************************************************//
template <class T>
void MasterPola <T>::idl_get_spec_radii(double* SpecRad) { 
	SpecRad[0]=this->gamma_TT_TT;
	if(Nmaps==3) {
		SpecRad[1]=this->gamma_EE_EE;
		SpecRad[2]=this->gamma_EE_BB;
		SpecRad[3]=this->gamma_TE_TE ;
		SpecRad[4]=this->gamma_EB_EB;
		SpecRad[5]=this->gamma_EB_4BLOCKS;
	}	
}

//****************************************************************************************************************************//
template <class T>
void MasterPola <T>::idl_get_inv_coupling_mat(double* invMat) { 
	long naxes[2],offset,offset2=0;
  	naxes[0] = (long) this->inv_MAT_TT_TT.axis(1);   
   	naxes[1] = (long) this->inv_MAT_TT_TT.axis(2);
	for(long kl2=0;kl2<	naxes[1];kl2++) {
		offset=naxes[0]*kl2;
		for(long kl1=0;kl1<naxes[0];kl1++) invMat[kl1+ offset]=this->inv_MAT_TT_TT(kl1,kl2);
	}
	offset2+=naxes[0]*naxes[1];
   	naxes[0] = (long) this->inv_MAT_TE_TE.axis(1);   
   	naxes[1] = (long) this->inv_MAT_TE_TE.axis(2);
	for(long kl2=0;kl2<	naxes[1];kl2++) {
		offset=naxes[0]*kl2+offset2;
		for(long kl1=0;kl1<naxes[0];kl1++) invMat[kl1+ offset]=this->inv_MAT_TE_TE(kl1,kl2);
	}
	offset2+=naxes[0]*naxes[1];
   	naxes[0] = (long) this->inv_MAT_EB_EB.axis(1); 
   	naxes[1] = (long) this->inv_MAT_EB_EB.axis(2);
	for(long kl2=0;kl2<	naxes[1];kl2++) {
		offset=naxes[0]*kl2+offset2;
		for(long kl1=0;kl1<naxes[0];kl1++) invMat[kl1+ offset]=this->inv_MAT_EB_EB(kl1,kl2);
	}	
	offset2+=naxes[0]*naxes[1];
   	naxes[0] = (long) this->inv_MAT_EB_4Blocks.axis(1)/2; 
   	naxes[1] = (long) this->inv_MAT_EB_4Blocks.axis(2)/2;
	//Diagonal Blocks
	for(long kl2=0;kl2<	naxes[1];kl2++) {
		offset=naxes[0]*kl2+offset2;
		for(long kl1=0;kl1<naxes[0];kl1++) invMat[kl1+ offset]=this->inv_MAT_EB_4Blocks(kl1,kl2);
	}	
	offset2+=naxes[0]*naxes[1];
	//OffDiagonal Blocks
	for(long kl2=0;kl2<	naxes[1];kl2++) {
		offset=naxes[0]*kl2+offset2;
		for(long kl1=0;kl1<naxes[0];kl1++) invMat[kl1+ offset]=this->inv_MAT_EB_4Blocks(kl1+naxes[0],kl2);
	}	
}


//****************************************************************************************************************************//
template <class T>
void MasterPola <T>::idl_get_master_pspec(double* MasterPspec) { 
	long offset= 0;
	long Nls=this->Master_TT_Powspec.size();
	if(Nls >0) {
		for(long kl=0;kl<Nls;kl++) MasterPspec[kl]= this->Master_TT_Powspec[kl];
		if(this->Master_EE_Powspec.size()== (unsigned int) Nls) {
			offset+=Nls;
			for(long kl=0;kl<Nls;kl++) MasterPspec[kl+offset]= this->Master_EE_Powspec[kl];
			if(this->Master_BB_Powspec.size()== (unsigned int) Nls) {
				offset+=Nls;
				for(long kl=0;kl<Nls;kl++) MasterPspec[kl+offset]= this->Master_BB_Powspec[kl];
				if(this->Master_TE_Powspec.size()== (unsigned int)Nls) {
					offset+=Nls;
					for(long kl=0;kl<Nls;kl++) MasterPspec[kl+offset]= this->Master_TE_Powspec[kl];
					if(this->Master_TB_Powspec.size()== (unsigned int) Nls) {
						offset+=Nls;
						for(long kl=0;kl<Nls;kl++) MasterPspec[kl+offset]= this->Master_TB_Powspec[kl];
						if(this->Master_EB_Powspec.size()== (unsigned int) Nls) {
							offset+=Nls;
							for(long kl=0;kl<Nls;kl++) MasterPspec[kl+offset]= this->Master_EB_Powspec[kl];
						}	
					}
				}
			}
		}
	}
}



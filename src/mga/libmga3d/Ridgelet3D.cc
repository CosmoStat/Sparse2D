/******************************************************************************
**                   Copyright (C) 2000 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Jean-Luc Starck
**
**    Date:  5/03/2000 
**    
**    File:  Ridgelet3D.cc
**
**    Modification history:
**
*******************************************************************************
**
**    DESCRIPTION Ridgelet3D transform and reconstruction
**    -----------  
**                 
******************************************************************************/

#include "Ridgelet3D.h"
// #include "LineCol.cc"


#define DEGUG_RID 0


/*************************************************************/
// Ridgelet3D transform
/*************************************************************/

const char * StringRid3DTransform (type_ridgelet3d_WTtrans  type)
{
    switch (type)
    {
        case  RID3D_OWT:
	      return ("Standard bi-orthogonal WT"); 
        case  RID3D_PYR_FFT: 
              return ("FFT based Pyramidal WT"); 
        default: 
	      return ("Undefined transform");
     }
}

type_ridgelet3d_WTclass rid3d_class(type_ridgelet3d_WTtrans Type)
{
    switch (Type)
    {
         case  RID3D_OWT:
	      return  RID3D_CL_ORTHO; 
        case  RID3D_PYR_FFT: 
              return  RID3D_CL_PYR; 
         default: 
	      return RID3D_CL_UNKNOW;
     }
}

/*************************************************************************/

static int max_scale_number (int N)
{
    int ScaleMax;
    ScaleMax  = iround((float)log((float) (N / 4. * 3.) / log(2.)));
    return ScaleMax;
}

/****************************************************************************/

void Ridgelet3D::reset(SubBand1D *SB1D)
{
    Ptr_SB1D = SB1D; 
    VarianceStab=False;
    MadThreshold=False;
    StatInfo=False;
    NbrScale=DEF_RID3D_NBR_SCALE;
    Verbose=False;
    BlockCubeOverLap=0;
    BlockOverlap=False;
    GetAutoNbScale=True;
    RidTrans=DEF_RID3D_TRANS; 
    RidClass = rid3d_class(RidTrans);
    WTinFFT=True;
    FirstDetectScale=0;
    KillLastScale=False;
    SetTable = True;
    // NbrBlock=Nxb=Nyb=Nzb=0;
    UserSizeBlockParam=0;
}

/*************************************************************************/

void Ridgelet3D::set_block_param_from_transform(int Nlt, int Nct, int BlockSize)
// initialize the variables for block reconstruction from
// the input ridgelet transform size
{
   int N3D;

   if (BlockSize == 0) BlockOverlap = False;

   // We need to calculate the size of the original cube
   // for setting the different parameters

   // The radon transform multiply be 2 the number of lines
   N3D = Nct;
 
   // For pyramidal transform, the redundancy of the WT is 2         
   if ((GetAutoNbScale == True) || (NbrScale > 1))
        if (rid3d_class(RidTrans) == RID3D_CL_PYR) N3D /= 2;

   // When overlapping, there is another redundancy of 2 in all directions
   if (BlockOverlap == True) N3D /= 2;
 
   // initialize the block parameters
//   cout << "N3d = " << N3D << " BS = " << BlockSize << endl;
   set_block_param(N3D, N3D, N3D, BlockSize);
}

/*************************************************************************/

void Ridgelet3D::set_block_param(int NxCube, int NyCube, int NzCube, int BlockSize)
{
    UserSizeBlockParam=BlockSize;
    RidClass = rid3d_class(RidTrans);

    if (BlockSize <= 0) BlockOverlap = False;
	B3D.BlockOverlap = (BlockOverlap>0)?True:False;
	B3D.Overlap = OverlapFactor;
    B3D.alloc(NxCube,NyCube,NzCube,BlockSize);
    int BS = block_size();
	
   // Determine the number of scales to use
   if (GetAutoNbScale == True)  NbrScale = max_scale_number(BS);
 
   if (Verbose == True) cout << "Number of scales (auto="<<GetAutoNbScale<<") = " << NbrScale << endl;
   if (NbrScale > 1) RidPixNbr = set_tab_pos(NbrScale, BS); // set the tables

    // Ridgelet3D image size
    // The FFT based Radon transform multiply by 2 the number of lines
    // If a pyramidal WT is used instead of an orthogonal one, 
    // the number of columns are multiplied by 2.
    // Overlapping increase the ridgelet image size by the overlapping ratio
    NlBlockTransSize = 3*BS*BS; 
    NcBlockTransSize = BS;

    NlRid = NlBlockTransSize;
    NcRid = NcBlockTransSize*nbr_block();

    if (RidClass == RID3D_CL_PYR) 
    {
        if ((GetAutoNbScale == True) || (NbrScale > 1))
        {
           NcRid *=2;
           NcBlockTransSize *=2;
        }
    }
 
    if (Verbose == True)
    {
       cout << " Ridgelet image size: Nl  = " <<  NlRid << " Nc = " <<  NcRid << endl;
       cout << " Number of block: NbrBx = " << nbr_block_nx() << " NbrBy = " << nbr_block_ny() <<  " NbrBz = " <<  nbr_block_nz() <<endl;
       cout << " BlockCubeSize = " << BS << endl;
       if (OverlapFactor>0)
          cout << " Overlapping size = " << BlockSize-ceil(BlockSize*(1-OverlapFactor)) << endl;
       cout << " Ridgelet Block image size: Nbrb = " << NlBlockTransSize << " Ncb = " << NcBlockTransSize << endl;
    }
}

/*************************************************************************/

void Ridgelet3D::block_transform(fltarray &Cube, Ifloat &Transf)
// Ridgelet3D transform of image block per block
{
	Bool UserVerbose = Verbose;
	
	int BS = block_size();
	//int Nc = Ima.nc();
	//int Nl = Ima.nl();
	//set_block_param(Nl, Nc, BlockSize);

	if ((Transf.nl() != NlRid) || (Transf.nc() != NcRid))
		Transf.resize(NlRid, NcRid);
	
	// Loops on the block
	// cout << "Block BS = " <<  BS << endl;
	// cout << "Block ImaTransBlock  = " <<  NlBlockTransSize << " " <<  NcBlockTransSize  << endl;
	//  cout << "Min InDat = " << Cube.min() <<  " Max InDat = " << Cube.max() << endl;
#if USE_OMP_R
	omp_set_nested(1);
	#pragma omp parallel for 
#endif
	for (int i = 0; i < nbr_block_nx(); i++)
	{
#if USE_OMP_R
	#pragma omp parallel for 
#endif
		for (int j = 0; j < nbr_block_ny(); j++)
		{
#if USE_OMP_R
	#pragma omp parallel for 
#endif
			for (int k = 0; k < nbr_block_nz(); k++)
			{
		//		cout << "Block " << i << " " << j << " " << k << endl;

				fltarray CubeBlock(BS, BS,BS);
				Ifloat ImaTransBlock(NlBlockTransSize,NcBlockTransSize);
				get_block_cube(i,j,k, Cube, CubeBlock);
				// CubeBlock.info("Block");
				// cout << " Transform_one_block " << CubeBlock.nx() << " " << CubeBlock.ny() << " " << CubeBlock.nz() << endl;

		//char filename[64];
		//sprintf(filename,"CubeBlock_%d.fits",i);
		//cerr<<"write file "<<filename<<endl;
		//fits_write_fltarr(filename,CubeBlock);		



				transform_one_block(CubeBlock, ImaTransBlock);
				// cout << " Put_block_trans  " << endl;
				//  cout << "Min trans = " <<  min(ImaTransBlock) << endl;
				// cout << "Max trans = " <<  max(ImaTransBlock) << endl;

		//sprintf(filename,"ImaTransBlock_%d.fits",i);
		//io_write_ima_float(filename,ImaTransBlock);		


				put_block_trans(i,j,k, Transf, ImaTransBlock);
		//		if ((i == 0) && (j == 0) && (k == 0)) Verbose = False;
			}
		}
	}
	// cout << "Min trans = " <<  min(Transf) << endl;
	// cout << "Max trans = " <<  max(Transf) << endl;
	Verbose = UserVerbose;
}

/***************************************/

void Ridgelet3D::block_recons(Ifloat &Transf, fltarray &Cube)
// Cube reconstruction from its ridgelet transform
// Block per block
{
    int BS = block_size();

    if ((Cube.nx() != nx()) || (Cube.ny() !=  ny())|| (Cube.nz() !=  nz())) 
      Cube.reform(nx(),ny(),nz());
    Cube.init();
    Bool UserVerbose = Verbose;
#if USE_OMP_R
		#pragma omp parallel for 
#endif
	for (int i = 0; i < nbr_block_nx(); i++)
	{
#if USE_OMP_R
		#pragma omp parallel for 
#endif
		for (int j = 0; j < nbr_block_ny(); j++)
		{
//#pragma omp parallel for 
			for (int k = 0; k < nbr_block_nz(); k++)
			{
		//		cout << "Block " << i << " " << j << " " << k << endl;

				fltarray CubeBlock(BS, BS,BS);
				Ifloat ImaTransBlock(NlBlockTransSize,NcBlockTransSize);
				
    		   // cout << "Block " << i << " " << j << endl;
    		   get_block_trans(i, j, k, Transf,  ImaTransBlock);
    		   // cout << " get_block_trans " << endl;

        		recons_one_block(ImaTransBlock, CubeBlock);
    		   // cout << "  recons_one_block " << endl;

        		if (BlockOverlap == True) add_block_cube(i, j, k, Cube, CubeBlock, True);
        		else put_block_cube(i, j, k, Cube, CubeBlock);
    		   // cout << "  add_block_ima   " << endl;
    		   Verbose = False;
		   }
	   }
    }
    Verbose = UserVerbose;
}

/***************************************/

void Ridgelet3D::variance_stab(Ifloat &Transf)
{
    int i,j;
    double C1 = 2.;
    double C2 = 3. / 8.;
    for (i = 0; i <  Transf.nl(); i++)
    for (j = 0; j <  Transf.nc(); j++)
    { 
       if (Transf(i,j) > 0.)
               Transf(i,j) = (float)(C1*sqrt((double) Transf(i,j) + C2));
       else Transf(i,j) = 0.;
    }
}

/***************************************/

void Ridgelet3D::inv_variance_stab(Ifloat &Transf)
{
   int i,j;
   double C1 = 4.;
   double C2 = 3. / 8.;
   for (i = 0; i <  Transf.nl(); i++)
   for (j = 0; j <  Transf.nc(); j++)
            Transf(i,j) = (float)((double) Transf(i,j)*Transf(i,j)/C1 - C2);
}


/****************************************************************************/
  
int Ridgelet3D::set_tab_pos(int NbrScale, int N)
// Set the table position, for band position
{
    int s,Ns=N;
    int TotSize=0;
    if (TabSize.nx() != NbrScale) TabSize.reform(NbrScale);
    if (TabPos.nx() != NbrScale)  TabPos.reform(NbrScale);
 
    if (RidClass != RID3D_CL_ORTHO)
    {
      for (s = 0; s < NbrScale; s++)
      {
         TotSize += Ns;
         TabSize(s) = Ns;
         Ns = (Ns+1)/2;
      }
    }
    else
    {
      for (s = 0; s < NbrScale; s++)
      {
         if (s != NbrScale-1) 
         {
            TabSize(s) = Ns/2;
            Ns = (Ns+1)/2;
         }
         else TabSize(s) = Ns;
         TotSize += TabSize(s);
      }
    }
    TabPos(NbrScale-1) = 0;
    for (s =  NbrScale-2; s >= 0; s--) 
        TabPos(s) = TabPos(s+1) + TabSize(s+1);

    // cout << "N = " << N << endl;
    // for (s = 0; s < NbrScale; s++) 
    //   cout << "Scale " << s+1 << " Pos = " << TabPos(s) << "Size = " << TabSize(s) << " End = " << TabPos(s)+TabSize(s)-1 << endl;

    return TotSize;
}

/****************************************************************************/
//int numero_du_block=0;
void Ridgelet3D::transform_one_block(fltarray &Cube, Ifloat &Transf)
{ 
	bool LocVerbose = false && Verbose;
	int i,j,Nlr,Ncr;
	
	Ifloat TransRadon;
	cfarray TransRadonCF;

	if (LocVerbose == true)  
		cout << "Ridgelet3D transform of Cube " << Cube.nx() << " " << Cube.ny() <<  " " << Cube.nz() <<endl;

//	if ((GetAutoNbScale == False) && (NbrScale < 2))
//	{
//		cout << "Error: number of scales must be larger than 2 ... " << endl;
//		exit(-1);
//	}

	switch(RidTrans)
	{
		case RID3D_OWT:
			WTinFFT=False;
			if (Ptr_SB1D == NULL)
			{
				cout << "Error: filter bank class not initialized ... " << endl;
				exit(-1);
			}
			break;
		case RID3D_PYR_FFT: 
			WTinFFT=True; 
			break;
		default: 
			cout << "Error: unknown ridgelet transform ... " << endl; exit(-1);
			break;
	}
	RidClass = rid3d_class(RidTrans);

	//char filename[64];
	//sprintf(filename,"Cube_%d.fits",numero_du_block);
	//fits_write_fltarr(filename, Cube);
	// identiques
	Radon3D RADON_local;

	// Radon transform of the input image
	if (WTinFFT == False) 
	{
		RADON_local.transform(Cube, TransRadon);
		Nlr = TransRadon.nl();
		Ncr = TransRadon.nc();
	}	
	else
	{
		RADON_local.fft_trans_cf(Cube, TransRadonCF);
		Nlr = TransRadonCF.ny();
		Ncr = TransRadonCF.nx();

		//fltarray prov(Ncr,Nlr);
		//for (int i=0;i<Ncr;i++)
		//for (int j=0;j<Nlr;j++)
		//   prov(i,j)=TransRadonCF(i,j).real();
		//sprintf(filename,"RadT_%d.fits",numero_du_block);
		//fits_write_fltarr(filename, prov);
		// pas identiques
	}
	//numero_du_block++;

#if DEGUG_RID
	if (WTinFFT == False) io_write_ima_float("xx_rad.fits", TransRadon);
#endif

	// Variance stabilization of the Radon transform
	if (VarianceStab == True) 
	{
		if (WTinFFT == False) variance_stab(TransRadon);
		else
		{
			if ((TransRadon.nl() != Nlr) || (TransRadon.nl() != Nlr))
				TransRadon.resize(Nlr,Ncr);
			RADON_local.invfft_line(TransRadonCF, TransRadon);
			variance_stab(TransRadon);
			RADON_local.fft_line(TransRadon, TransRadonCF);
		}
	}


	// Apply the 1D wavelet transform on each line of the Radon transform 
	if (NbrScale > 1)
	{
		int BNlRid = Nlr;
		// For non orthogonal WT transform, the redundancy is 2
		int BNcRid;
		if  (RidClass == RID3D_CL_ORTHO)
			BNcRid = Ncr;
		else BNcRid = 2*Ncr;

		if (LocVerbose == true) 
		{
			cout << "Radon image size: Nl = " <<  Nlr<< " Nc = " <<  Ncr << endl;
			cout << "Ridgelet3D image size: Nl = " << BNlRid  << " Nc = " <<  BNcRid  << endl;
		}


		if ((Transf.nl() != BNlRid) || (Transf.nc() != BNlRid))
			Transf.resize(BNlRid,BNcRid);

		// Orthogonal wavelet transform
		if (RidTrans == RID3D_OWT)  
		{
			fltarray Line(Ncr);      // Radon line
			fltarray WTLine(BNcRid); // Ridgelet3D line
			Ortho_1D_WT WT1D(*Ptr_SB1D);
			for (i= 0; i < Nlr; i++)
			{
				for (j=0; j < Ncr; j++)
					Line(j) = TransRadon(i,j);
				WT1D.transform(Line, WTLine, NbrScale);
				for (j=0; j < BNcRid; j++)
					Transf(i,j) = WTLine(j);
			} 
		}
		else // FFT-based pyramidal wavelet transform
		{
			WT1D_FFT WT1D(WT_PYR_FFT_DIFF_RESOL, NbrScale, Ncr);
			int NPixTrans = WT1D.pyr_np();
//cerr << "NPixTrans = " << NPixTrans << " Ncr = " << Ncr << endl;
			fltarray WTLine(NPixTrans);
			cfarray LineCF(Ncr,(char*)"line CF");         
			// Radon line of size Ncr

			for (i= 0; i < Nlr; i++)
			{
				// cout << "i = " << i << endl;
				// for (j=0; j < Nc; j++) Line(j) = TransRadon(i,j);
				// WT1D.transform(Line, WTLine, NbrScale);
				for (j=0; j < Ncr; j++)
					LineCF(j) = TransRadonCF(j,i);
				WT1D.transform_cf(LineCF, WTLine, NbrScale);
				for (j=0; j < NPixTrans; j++)
					Transf(i,j) = WTLine(j);
				for (j=NPixTrans; j < BNcRid; j++)
					Transf(i,j) = 0.;
			}
		}
	}
	else 
	{
		if (RidTrans != RID3D_OWT)  
		{
			TransRadon.resize(Nlr,Ncr);
			RADON_local.invfft_line(TransRadonCF, TransRadon);
			INFO_X(TransRadon, (char*)"TR");
		}
		Transf = TransRadon;
	}
}


/****************************************************************************/
 
void Ridgelet3D::transform(fltarray &Cube, Ifloat &Transf, int BlockSize)
{
   set_block_param(Cube.nx(), Cube.ny(), Cube.nz(), BlockSize);
   block_transform(Cube, Transf);
   if (StatInfo == True) info_stat(Transf);
}

/****************************************************************************/

void Ridgelet3D::recons_one_block(Ifloat &Transf, fltarray &Cube)
{
    int i,j;
    int Nl = Transf.nl();
    int Nc = Transf.nc();

Radon3D RADON_local;

   // cout << "NbrScale = " << NbrScale << endl;
    if (NbrScale < 2)  RADON_local.recons(Transf, Cube);
    else
    {
    switch(RidTrans)
    {
        case RID3D_OWT:     
              WTinFFT=False; 
              if (Ptr_SB1D == NULL)
              {
                     cout << "Error: filter bank class not initialized ... " << endl;
                     exit(-1);
              }
              break;
        case RID3D_PYR_FFT: WTinFFT=True; break;
        default: 
         cout << "Error: unknown ridgelet transform ... " << endl; exit(-1);
    }
    RidClass = rid3d_class(RidTrans);

    if (RidTrans == RID3D_OWT) 
    {
       fltarray WTLine(Nc);
       fltarray RecLine(Nc);
    
       if (Verbose == True) 
       {
          cout << "Ridgelet3D reconstruction ..." << endl;
          cout << "Input image size: Nl = " << Nl << " Nc = " << Nc << endl;
       }
       // if (GetAutoNbScale == True) NbrScale = max_scale_number(Nc);

       if (NbrScale < 2)
       {
           cout << "Error: number of scales must be larger than 2 ... " << endl;
           exit(-1);
       }
 
       if (Verbose == True) cout << "Number of scales = " << NbrScale << endl;
#if DEGUG_RID
       io_write_ima_float("xx_trans1.fits", Transf);
#endif
      if (NbrScale > 1)
      {
         Ortho_1D_WT WT1D(*Ptr_SB1D);
         for (i= 0; i < Nl; i++)
         {
            for (j=0; j < Nc; j++) WTLine(j) = Transf(i,j);
            WT1D.recons(WTLine, RecLine, NbrScale);
            for (j=0; j < Nc; j++)  Transf(i,j) = RecLine(j);
         }
      }
       if (VarianceStab == True) inv_variance_stab(Transf);
#if DEGUG_RID
       RADON_local.Verbose = Verbose;
       io_write_ima_float("xx_rad_rec.fits", Transf);
#endif
       RADON_local.recons(Transf, Cube);
    }
    else
    {
       Ifloat TransRadon;
       // fltarray RecLine(Nc/2); 

       cfarray TransRadonCF(Nc/2,Nl);
       cfarray RecLineCF(Nc/2) ;
                               // RecLine is one line of the Radon transform
                               // FFT WT1D has a factor 2 of redundancy
                               // The radon transform has half the number
			       // of columns than the Ridglet transform
       // if (GetAutoNbScale == True) NbrScale = max_scale_number(Nc/2);
       if (Verbose == True) 
       {
          cout << "Ridgelet3D reconstruction ..." << endl;
          cout << "Input image size: Nl = " << Transf.nl() << " Nc = " << Transf.nc() << endl;
          cout << "Number of scales = " << NbrScale << endl;
       }
       WT1D_FFT WT1D(WT_PYR_FFT_DIFF_RESOL, NbrScale, Nc/2);
       int NPixTrans = WT1D.pyr_np();
       fltarray WTLine(NPixTrans);
       
//        cout << "1D WT" << endl;
//        i = 10;
//        for (j=0; j < NPixTrans; j++) WTLine(j) = Transf(i,j);
//        cout << "Max = " << WTLine.max() << endl;
//        WT1D.recons_cf(WTLine, RecLineCF, NbrScale);
//        WTLine.init();
//        for (j=0; j < Nc/2; j++)  WTLine(j) = RecLineCF(j).real();
//        cout << "Max rec = " << WTLine.max() << endl;

       for (i= 0; i < Nl; i++)
       {
          for (j=0; j < NPixTrans; j++) WTLine(j) = Transf(i,j);
          // WT1D.recons(WTLine, RecLine, NbrScale);
          // for (j=0; j < Nc/2; j++) TransRadon(i,j) = RecLine(j);
	  // RecLineCF.resize (1, N);

	  WT1D.recons_cf(WTLine, RecLineCF, NbrScale);
          for (j=0; j < Nc/2; j++) TransRadonCF(j,i) = RecLineCF(j);
       }
       if (VarianceStab == True) 
       {
           TransRadon.alloc(Nl,Nc/2,(char*)"trans radon");;
 	   RADON_local.invfft_line(TransRadonCF, TransRadon);
           inv_variance_stab(TransRadon);
	   RADON_local.fft_line(TransRadon, TransRadonCF);
       }
#if DEGUG_RID
       RADON_local.Verbose = Verbose;
       TransRadon.resize(Nl,Nc/2);
       for (i= 0; i < Nl; i++)
       for (j= 0; j < Nc/2; j++) TransRadon(i,j) = TransRadonCF(j,i).real();
       io_write_ima_float("xx_rad_rec.fits", TransRadon);
#endif
       // cout << "Rec radon " << TransRadonCF.nl() << " " << TransRadonCF.nc() << endl;
        // RADON.recons(TransRadon, Image);
       RADON_local.fft_rec_cf(TransRadonCF, Cube);
    }
   }
}

/****************************************************************************/

void Ridgelet3D::recons(Ifloat &Transf, fltarray &Cube, int BlockSize)
{
   // set_block_param_from_transform(Transf.nl(), Transf.nc(), BlockSize);
   block_recons(Transf, Cube);
}

/****************************************************************************/

void Ridgelet3D::recons(Ifloat &Transf, fltarray &Cube, int BlockSize, 
                        int InputNx, int InputNy, int InputNz)
{
   set_block_param(InputNx, InputNy, InputNz, BlockSize);
   block_recons(Transf, Cube);
}

/****************************************************************************/

void Ridgelet3D::info_stat(Ifloat &Transf)
{
    int s;
    Ifloat Buff(rid_nl(), rid_nc(), (char*)"Scale");

    if (NbrScale > 1)
    {
        for (s=0; s < NbrScale; s++)
        {
           get_scale(Transf, Buff,s);
           // cout << "WT1D Band " << s+1 << " Sigma " << Buff.sigma() << endl;
           cout << "Band " << s+1;
           cout << "   Mean = " << average(Buff)  << " Sigma " << sigma(Buff) << endl;
           cout << "   Min  = " << min(Buff) << " Max = " <<  max(Buff) << endl;
           cout << "   Skewness  = " << skewness(Buff.buffer(),Buff.n_elem()) 
	        << " Kurtosis = " <<  curtosis(Buff.buffer(),Buff.n_elem()) << endl << endl;
        
	
	}
    }
}

 
/****************************************************************************/
//  RIDGELET DATA Structure management
/****************************************************************************/

// float Ridgelet3D::get_coef(Ifloat &Trans, int i, int j, int s, int Bi, int Bj, int Bk)
// {
//    int Indi = 0;
//    int Indj = 0;
//    int Indk = 0;     
//    if (NbrScale <= 1)
//    {
//        cout << "Error in get_coef: NbrScale = " << NbrScale << endl;
//        exit(-1);
//    }
//    else
//    {
//        Indi = ipos(s, Bi) + i;
//        Indj = jpos(s, Bj+Bk) + j;
//    }
//    return  Trans(Indi,Indj);
// }
// 
// /****************************************************************************/
// 
// void Ridgelet3D::put_coef(Ifloat &Trans, int i, int j, float Val, int s, int Bi, int Bj, int Bk)
// {
//    int Indi = 0;
//    int Indj = 0;
//      
//    if (NbrScale <= 1)
//    {
//       cout << "Error in get_coef: NbrScale = " << NbrScale << endl;
//       exit(-1);
//    }
//    else
//    {
//       Indi = ipos(s,Bi) + i;  
//       Indj = jpos(s,Bj+Bk) + j;
//    }
//    Trans(Indi,Indj) = Val;
// }  
// 
// /****************************************************************************/

void Ridgelet3D::get_scale(Ifloat &ImaTrans, Ifloat &ImaScale, int s)
{
   int i,j;
   int Nls = size_scale_nl(s);
   int Ncs = size_scale_nc(s);
   int Posi = ipos(s);
   int Posj = jpos(s);

   if ((ImaScale.nl() != Nls) || (ImaScale.nc() != Ncs)) 
                                             ImaScale.resize(Nls,Ncs);

   for (i=0; i < Nls; i++)
   for (j=0; j < Ncs; j++) ImaScale(i,j) = ImaTrans(Posi+i,Posj+j);
}

/****************************************************************************/

void Ridgelet3D::put_scale(Ifloat &ImaTrans, Ifloat &ImaScale, int s)
{
   int i,j;
   int Nls = size_scale_nl(s);
   int Ncs = size_scale_nc(s);
   int Posi = ipos(s);
   int Posj = jpos(s);

   if ((ImaScale.nl() != Nls) || (ImaScale.nc() != Ncs)) 
   {
       cout << "Error: the scale image has not the correct size ... " << endl;
       exit(-1);
   }
   for (i=0; i < Nls; i++)
   for (j=0; j < Ncs; j++) ImaScale(i,j) = ImaTrans(i+Posi,Posj+j);
}

/****************************************************************************/
   

/****************************************************************************/

//               Ridgelet3D Block manadgement

/*************************************************************************/

// float lap_weight(int BlockSize, int Overlap, int k, int B, int Nb)
// // BlockSize = block image size
// // Overlap = overlap size
// // int k = pixel position in the block B
// // Nb = Number of blocks
// // return a weight value, following pixel position in the block
// {
//      float ValReturn=1.;
//      if (Overlap != 0)
//      {
//         int Center = Overlap + BlockSize / 2;
//         // int Size = BlockSize+2*Overlap;
//         if ((B != 0) && (k < Center)) 
//           ValReturn = (float) pow(sin((double) k / (double) Center*PI/2.), 2.);
//           // ValReturn = (float) k / (float) Center;  
//         else if ((B != Nb-1) && (k >= Center)) 
//           ValReturn = (float) pow(cos((double)(k-Center)/ (double) Center *PI/2.), 2.);
//           // ValReturn =  (float) (Size - k) / (float) Center;  
//      }
//      return ValReturn;
// }
// 
// /*************************************************************************/
// 
// float Ridgelet3D::get_weight(int Bk, int Bl, int Bm, 
//                              int k, int l, int m)
// // return the weight value for pixel position k,l in the block Bi,Bj
// {
//     float ValReturn=1;
//     
//     ValReturn = lap_weight(BlockCubeSize, BlockCubeOverLap, k, Bk, Nxb)* 
//                 lap_weight(BlockCubeSize, BlockCubeOverLap, l, Bl, Nyb) * 
//                 lap_weight(BlockCubeSize, BlockCubeOverLap, m, Bm, Nzb);
//  
//     return ValReturn;
// }
// 
// /*************************************************************************/
// 
// void Ridgelet3D::get_block_cube(int Bi, int Bj, int Bk, fltarray &Cube, fltarray &CubeBlock)
// // Extract the block (Bi,Bj,Bk) from cube and put it in CubeBlock
// {
//    int k,l,m;
//    int Depi = BlockCubeSize*Bi - BlockCubeOverLap;
//    int Depj = BlockCubeSize*Bj - BlockCubeOverLap;
//    int Depk = BlockCubeSize*Bk - BlockCubeOverLap;
// 
//    for (k = 0; k < BS; k++)
//    for (l = 0; l < BS; l++)
//    for (m = 0; m < BS; m++)
//    {
//       int Indk = test_index_mirror(Depi+k,Nx);
//       int Indl = test_index_mirror(Depj+l,Ny);
//       int Indm = test_index_mirror(Depk+m,Nz);
//       CubeBlock(k,l,m) = Cube(Indk,Indl,Indm);  
//    }
// }
// 
// /*************************************************************************/
// 
// void Ridgelet3D::put_block_cube(int Bi, int Bj, int Bk, fltarray &Cube, fltarray &CubeBlock)
// // Put the block (Bi,Bj,Bk) CubeBlock in Cube  
// {
//    int k,l,m;
//    int Depi = BlockCubeSize*Bi - BlockCubeOverLap;
//    int Depj = BlockCubeSize*Bj - BlockCubeOverLap;
//    int Depk = BlockCubeSize*Bk - BlockCubeOverLap;
// 
//    for (k = 0; k < BS; k++)
//    for (l = 0; l < BS; l++)
//    for (m = 0; m < BS; m++)
//          if ((Depi+k >= 0) && (Depi+k < Nx) 
//              && (Depj+l >= 0) && (Depj+l < Ny)
// 	     && (Depk+m >= 0) && (Depk+m < Nz))
//                                Cube(Depi+k,Depj+l,Depk+m) = CubeBlock(k,l,m);  
// }
// 
// /*************************************************************************/
// 
// void Ridgelet3D::add_block_cube(int Bi, int Bj, int Bk, fltarray &Cube, fltarray &CubeBlock)
// // Add the block (Bi,Bj,Bk) CubeBlock in Cube  with weighted values
// 
// {
//    int k,l,m;
//    int Depi = BlockCubeSize*Bi - BlockCubeOverLap;
//    int Depj = BlockCubeSize*Bj - BlockCubeOverLap;
//    int Depk = BlockCubeSize*Bk - BlockCubeOverLap;
// 
//    for (k = 0; k < BS; k++)
//    for (l = 0; l < BS; l++)
//    for (m = 0; m < BS; m++)
//          if ((Depi+k >= 0) && (Depi+k < Nx) 
//              && (Depj+l >= 0) && (Depj+l < Ny)
// 	     && (Depk+m >= 0) && (Depk+m < Nz))
//             Cube(Depi+k,Depj+l,Depk+m) += CubeBlock(k,l,m)*get_weight(Bi,Bj,Bk,k,l,m);  
// }

/*************************************************************************/

void Ridgelet3D::get_block_trans(int Bi, int Bj, int Bk, Ifloat & Transf, Ifloat &  ImaTransBlock)
// Get a block from the ridgelet transformed image 
{
   int s,i,j,k,l;
   int NumBlock = num_block(Bi,Bj,Bk);
//cerr<<" start getblock "<<Bi<<","<<Bj<<","<<Bk<<endl;
   if (NbrScale <= 1)
   {
      int Depi = 0;
      int Depj =  NumBlock*NcBlockTransSize;
      for (k = 0; k < NlBlockTransSize; k++)
      for (l = 0; l < NcBlockTransSize; l++)  
             ImaTransBlock(k,l) = Transf(Depi+k, Depj+l);
   }
   else
   {
      for (s=0; s < NbrScale; s++)
      {
          // cout << "Scale = " << s+1 << endl;
          int NFirst= rid_pos(s);
          int Nsize = rid_size(s);
          int Posi = ipos(s,NumBlock);
          int NFirstTrans= jpos(s,NumBlock);

          // cout << " ImaTransBlock.nl() = " << ImaTransBlock.nl() << endl;
          // cout << "  NFirst = " << NFirst << endl;
          // cout << "  NFirstTrans  = " <<  NFirstTrans << endl;

          for (i=0; i < ImaTransBlock.nl(); i++)
          for (j=0; j < Nsize; j++)
                ImaTransBlock(i,j+NFirst) = Transf(i+Posi,NFirstTrans+j);
      }
   }
//cerr<<" end getblock "<<Bi<<","<<Bj<<","<<Bk<<endl;
}

/*************************************************************************/

void Ridgelet3D::put_block_trans(int Bi, int Bj, int Bk, Ifloat & Transf, Ifloat & ImaTransBlock)
// Put a block from the ridgelet transformed image
{
   int s,i,j,k,l;
   int NumBlock = num_block(Bi,Bj,Bk);
//cerr<<" start putblock "<<Bi<<","<<Bj<<","<<Bk<<endl;
   if ((Transf.nl() < ImaTransBlock.nl()) || (Transf.nc() < ImaTransBlock.nc()))
   {
       cout << "Block size < image size " <<  endl;
       cout << "Nbr line: " << ImaTransBlock.nl() << " " << Transf.nl() << endl;
       cout << "Nbr Col: " << ImaTransBlock.nc() << " " << Transf.nc() << endl;
       exit(-1);
   }
  // cout << "Nbr line block and ima: " << ImaTransBlock.nl() << " " << Transf.nl() << endl;
  // cout << "Nbr Col  block and ima: " << ImaTransBlock.nc() << " " << Transf.nc() << endl;

   if (NbrScale <= 1)
   {
      int Depi = 0;
      int Depj =  NumBlock*NcBlockTransSize;
      for (k = 0; k < NlBlockTransSize; k++)
      	for (l = 0; l < NcBlockTransSize; l++)  
            Transf(Depi+k, Depj+l) = ImaTransBlock(k,l);
   }
   else
   {
      for (s=0; s < NbrScale; s++)
      {
          int NFirst= rid_pos(s);
          int Nsize = rid_size(s);
          int Posi = ipos(s,NumBlock);
          int NFirstTrans= jpos(s,NumBlock);
          // cout << "Scale " << s+1 << " NFirst = " << NFirst << endl;
          // cout << "NFirst end  = " << NFirst  + Nsize -1  << endl;
          // cout << "NFirstTrans = " <<  NFirstTrans << endl;
          // cout << " END NFirstTrans " << NFirstTrans + Nsize -1 << endl;

          for (i=0; i < ImaTransBlock.nl(); i++)
          	for (j=0; j < Nsize; j++)
                 Transf(i+Posi,NFirstTrans+j) = ImaTransBlock(i,j+NFirst);
      }
   } 
//cerr<<" end putblock "<<Bi<<","<<Bj<<","<<Bk<<endl;
}

/*************************************************************************/

static void mr_io_name (char *File_Name_In, char *File_Name_Out)
{
    int L;

    strcpy (File_Name_Out, File_Name_In);

    L = strlen (File_Name_In);
    if ((L < 4) || (File_Name_In[L-1] != 'd')
                || (File_Name_In[L-2] != 'i')
                || (File_Name_In[L-3] != 'r')
                || (File_Name_In[L-4] != '.'))
    {
        strcat (File_Name_Out, ".rid");
    }
}

/****************************************************************************/
/****************************************************************************/
/****************************************************************************/

/*--------------------------------------------------------------------------*/
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

/*--------------------------------------------------------------------------*/

void Ridgelet3D::mr_io_fill_header(fitsfile *fptr)
{
  int status = 0; // this means OK for cfitsio !!!
    /*****************************************************/
     /* write optional keyword to the header */
    /*****************************************************/
  if ( ffpkyj(fptr, (char*)"Nx", (long) nx(),(char*)"Nx",&status))
     PrintError( status );  
  if ( ffpkyj(fptr,(char*)"Ny",(long) ny(),(char*)"Ny",&status))
     PrintError( status );  
  if ( ffpkyj(fptr,(char*)"Nz",(long) nz(),(char*)"Nz",&status))
     PrintError( status ); 
  if ( ffpkyj(fptr,(char*)"BlockSiz",(long)UserSizeBlockParam,(char*)"BlockSize",&status))
     PrintError( status ); 
  if ( ffpkyj(fptr, (char*)"Nxb", (long) nbr_block_nx(),(char*)"Number of blocks Nxb",&status))
     PrintError( status );  
  if ( ffpkyj(fptr,(char*)"Nyb",(long)nbr_block_ny(),(char*)"Number of blocks Nyb",&status))
     PrintError( status );  
  if ( ffpkyj(fptr,(char*)"Nzb",(long) nbr_block_nz(),(char*)"Number of blocks Nzb",&status))
     PrintError( status ); 
  int Val = (BlockOverlap == True) ? 1: 0;
  if ( ffpkyj(fptr,(char*)"Overlap",(long) Val,(char*)"Block overlap",&status))
     PrintError( status );
  
  if ( ffpkyj(fptr,(char*)"NbrScale",(long)NbrScale,(char*)"Number of scales",&status))
     PrintError( status );      
           
  if ( ffpkyj(fptr, (char*)"Type_Tra", (long) RidTrans , 
                          (char*)StringRid3DTransform(RidTrans), &status))
     PrintError( status );  
     
} 


/****************************************************************************/

void Ridgelet3D::write (char *Name, Ifloat &Ima)
/* new version with fits */
{
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
 naxis = 2;
 naxes[0] = rid_nc(); 
 naxes[1] = rid_nl();

 // write first header part (parameters)
 if ( ffphpr(fptr,simple,bitpix,naxis,naxes,pcount,gcount,extend,&status) )
     PrintError( status );          /* call PrintError if error occurs */

  // write the header of the multiresolution file
  mr_io_fill_header(fptr);

  nelements = naxes[0] * naxes[1];
  if ( ffppre(fptr, group, firstpixel, nelements, Ima.buffer(), &status) )
              PrintError( status );  
 
 /* close the FITS file */
 if ( ffclos(fptr, &status) )  PrintError( status );  
// cout << " end of write fits " << endl;

}
/****************************************************************************/

void Ridgelet3D::read (char *Name, Ifloat & Ima)
{
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
 
#if DEBUG_IO
    cout << "Read in " << filename << endl;
#endif
   
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

     nelements = naxes[0] * naxes[1];
     // cout << " begin to read " << endl;
     Ima.alloc(naxes[1], naxes[0], (char*)"read io");
     
    if (ffgkyj(fptr,(char*)"Nx", &mon_long, comment, &status)) PrintError( status );
    int Nxi = (int) mon_long;  
    if (ffgkyj(fptr,(char*)"Ny", &mon_long, comment, &status)) PrintError( status );
    int Nyi = (int) mon_long; 
    if (ffgkyj(fptr,(char*)"Nz", &mon_long, comment, &status)) PrintError( status );
    int Nzi = (int) mon_long;  
    if (ffgkyj(fptr,(char*)"BlockSiz", &mon_long, comment, &status)) PrintError( status );
    int Bi = (int) mon_long;
    if (ffgkyj(fptr,(char*)"NbrScale", &mon_long, comment, &status)) PrintError( status );
    NbrScale = (int) mon_long;
    if (ffgkyj(fptr,(char*)"Type_Tra", &mon_long, comment, &status)) PrintError( status );
    RidTrans = (type_ridgelet3d_WTtrans) mon_long;
    if (ffgkyj(fptr,(char*)"Overlap", &mon_long, comment, &status)) PrintError( status );
    BlockOverlap = (mon_long == 1) ? True: False;
    set_block_param(Nxi,Nyi,Nzi,Bi);

#if DEBUG_IO 
    cout << "Read in " << filename << endl;
    cout << "Nx = " << size_ima_nl () << endl;
    cout << "Ny = " << size_ima_nc () << endl;
    cout << "Nbr_Plan = " << nbr_scale () << endl;
    cout << "Type_Transform = " << RidTrans << " " <<
             StringRid3DTransform(RidTrans) << endl;
 #endif

    Ptr = Ima.buffer();
    if ( ffgpve(fptr, 1, 1, nelements, nulval, Ptr, &anynul, &status))
             PrintError( status );
  
  if ( ffclos(fptr, &status) ) PrintError( status );
// cout << " end of read fits file " << endl;

#if DEBUG_IO
    cout << "Read out " << filename << endl;
#endif
}

/****************************************************************************/

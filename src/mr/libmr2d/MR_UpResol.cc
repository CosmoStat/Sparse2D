
// pyramidal median transform with integer
// #include "NR.h"
#include "IM_Obj.h"
#include "IM_IO.h"
#include "IM_BIO.h"
#include "MR_Obj.h"
#include "IM_CompTool.h"
#include "IM_Comp.h"
#include "IM_Noise.h"
#include "MR_Noise.h"
#include "MR_Comp.h"

/***************************************************/

void  BlockMRInfo::write(char *FileName)
{
   FILE *FileDes;
   FileDes = fopen(FileName, "w");
   fprintf(FileDes, "%d %d %d %d\n", L_Nl, L_Nc, Nl, Nc);
   fprintf(FileDes, "%d %d %d\n", Resol, BlockSize, KeepResi);
   fprintf(FileDes, "%d %d %d %d\n", FirstBlockNl, LastBlockNl, FirstBlockNc, LastBlockNc);
   fclose(FileDes);
}

/***************************************************/

void  BlockMRInfo::read(char *FileName)
{
   FILE *FileDes;
   FileDes = fopen(FileName, "r");
   if (FileDes == NULL)
   {
       cerr << "Error: cannot open file " << FileName << endl;
       exit(-1);
   }
   fscanf(FileDes, "%d %d %d %d\n", &L_Nl, &L_Nc, &Nl, &Nc);
   fscanf(FileDes, "%d %d %d\n", &Resol, &BlockSize, &KeepResi);
   fscanf(FileDes, "%d %d %d %d\n", &FirstBlockNl, &LastBlockNl, &FirstBlockNc, &LastBlockNc);
   fclose(FileDes);
}

/***************************************************/

BlockMRInfo::BlockMRInfo (char *FileName)
{
   int L_Nli, L_Nci, Nli, Nci, Scale, Bs;
   int Fbi,Lbi,Fbj,Lbj,Kr;
   FILE *FileDes;
   FileDes = fopen(FileName, "r");
   if (FileDes == NULL)
   {
       cerr << "Error: cannot open file " << FileName << endl;
       exit(-1);
   }
   if ((fscanf(FileDes, "%d %d %d %d\n", &L_Nli, &L_Nci, &Nli, &Nci) != 4)
     || (fscanf(FileDes, "%d %d %d\n", &Scale, &Bs, &Kr) != 3)
     || (fscanf(FileDes, "%d %d %d %d\n", &Fbi, &Lbi, &Fbj, &Lbj) != 4))
   {
      cerr << "Error: cannot read in file " << FileName << endl;
       exit(-1);
   }
   fclose(FileDes);
   
   init_param(L_Nli, L_Nci, Bs);
   Resol = Scale;
   FirstBlockNl = Fbi;
   FirstBlockNc = Fbj;
   LastBlockNl = Lbi;
   LastBlockNc = Lbj;
   KeepResi = Kr;
   set_pos();
   if ((Nl != Nli) || (Nc != Nci))
   {
      cerr << "Error: values are not correct in file " <<  FileName << endl;
      exit(-1);
   }
}

/***************************************************/

BlockMRInfo::BlockMRInfo(int Nli, int Nci, int Bs)
{
  init_param(Nli, Nci, Bs);
}

/***************************************************/

void BlockMRInfo::init_param(int Nli, int Nci, int Bs)
{
  Nl = L_Nl = Nli;
  Nc = L_Nc = Nci;
  BlockSize = Bs;
  Verbose = False;
  KeepResi = 0;
  
  if (BlockSize < 0) BlockSize = MAX(Nl,Nc);
  
  NbrBlockNl = L_Nl/BlockSize;
  // if (NbrBlockNl * BlockSize != L_Nl)  NbrBlockNl++;
  NbrBlockNc  = L_Nc/BlockSize;
  // if (NbrBlockNc * BlockSize != L_Nc)  NbrBlockNc++;
  NbrBlock = NbrBlockNl*NbrBlockNc;

  TabNl.alloc(NbrBlockNl,NbrBlockNc);
  TabNc.alloc(NbrBlockNl,NbrBlockNc);
  TabDepNl.alloc(NbrBlockNl,NbrBlockNc);
  TabDepNc.alloc(NbrBlockNl,NbrBlockNc);  
  for (int i=0; i < NbrBlock; i++)
  {
      int Posi = i / NbrBlockNc;
      int Posj = i % NbrBlockNc;

      if (Posi <  NbrBlockNl-1)  TabNl(Posi,Posj) = BlockSize;
      else  TabNl(Posi,Posj) = L_Nl - BlockSize * Posi;
      
      if (Posj <  NbrBlockNc-1) TabNc(Posi,Posj) = BlockSize;
      else TabNc(Posi,Posj) =  L_Nc - BlockSize * Posj;
      
      TabDepNl(Posi,Posj) = BlockSize*Posi;  
      TabDepNc(Posi,Posj) = BlockSize*Posj;
      //cout << "Block " << i+1 << ": " << TabNl(Posi,Posj) << "x" << TabNc(Posi,Posj) <<
      //        " DepNl = " << TabDepNl(Posi,Posj) << " DepNc = " << TabDepNc(Posi,Posj) << endl;
  }
  FirstBlockNl= FirstBlockNc=0;
  LastBlockNl= NbrBlockNl-1;
  LastBlockNc=NbrBlockNc-1;
  Resol = 0;
}

/***************************************************/

void BlockMRInfo::set_pos()
{
   int i,j;
   
   for (i=FirstBlockNl; i <= LastBlockNl; i++)
   for (j=FirstBlockNc; j <= LastBlockNc; j++)
   {
      int Re = MAX(Resol, 0);
      TabNl(i,j) = size_resol(Re, TabNl(i,j));
      TabNc(i,j) = size_resol(Re, TabNc(i,j));
   }
   
   for (i=FirstBlockNl; i <= LastBlockNl; i++)
   for (j=FirstBlockNc; j <= LastBlockNc; j++)
   {
     if (i == FirstBlockNl) TabDepNl(i,j) = 0;
     else TabDepNl(i,j) = TabDepNl(i-1,j) + TabNl(i-1,j);
    
     if (j == FirstBlockNc)  TabDepNc(i,j) = 0;
     else TabDepNc(i,j) = TabDepNc(i,j-1) + TabNc(i,j-1);
   }
   Nl = 0;
   Nc = 0;

   for (i=FirstBlockNl; i <= LastBlockNl; i++) Nl += TabNl(i, FirstBlockNc);
   for (j=FirstBlockNc; j <= LastBlockNc; j++) Nc += TabNc(FirstBlockNl,j);
   if (Verbose == True) cout << "Image: Nl = " << Nl << " Nc = " << Nc << endl;
}

/***************************************************/

BlockMRInfo::BlockMRInfo (BlockMRInfo &BI, int Scale)
{
   init_param(BI.ima_nl(), BI.ima_nc(), BI.block_size());
   Resol = Scale;
   KeepResi = BI.KeepResi;
   FirstBlockNl = BI.fb_nl();
   FirstBlockNc = BI.fb_nc();
   LastBlockNl = BI.lb_nl();
   LastBlockNc = BI.lb_nc();
   set_pos();
}
/***************************************************/

BlockMRInfo::BlockMRInfo(int Zoomi, int Zoomj, BlockMRInfo &BV,  
                          int NlVisu, int NcVisu)
{
  int HalfZoomSizei,HalfZoomSizej,BegZoomi,EndZoomi,BegZoomj,EndZoomj;
  int K, BB, BE, BB_NEW, BE_NEW;

  init_param(BV.ima_nl(), BV.ima_nc(), BV.block_size());
  Resol = BV.resol()-1;
  KeepResi = BV.KeepResi;
  
  if ((Zoomi < 0) || (Zoomj < 0) || (Zoomi >= BV.nl()) || (Zoomj >= BV.nc()))
  {
     cerr << "Error: zooming parameters are incorrect ... " << endl;
     cerr << "       " << Zoomi << "  " << Zoomj << endl;
     exit(-1);
  }
  // cout << " NlVisu  = " << NlVisu  << " NcVisu = " << NcVisu  << endl;
  
   // if the double resolution can be contained in the window, then ok
  if (((Nl <= NlVisu) && (Nc <= NcVisu)) || (Resol < 0))
  {
     if (Verbose == True)
         cout << "out: Nl = " << Nl << "Nc = " << Nc << endl;
      // if (Resol < 0) Resol = 0;
     FirstBlockNl = BV.fb_nl();
     FirstBlockNc = BV.fb_nc();
     LastBlockNl = BV.lb_nl();
     LastBlockNc = BV.lb_nc();
     set_pos();
  }
  else
  {
    //  Find the aera to zoom
    HalfZoomSizei = NlVisu / 4;
    HalfZoomSizej = NcVisu / 4;

    BegZoomi = Zoomi - HalfZoomSizei;
    EndZoomi = Zoomi + HalfZoomSizei - 1;
    BegZoomj = Zoomj -  HalfZoomSizej;
    EndZoomj = Zoomj +  HalfZoomSizej - 1;

    if (BegZoomi < 0)
    {
       BegZoomi = 0;
       EndZoomi = NlVisu/2 - 1;
    }
 
    if (BegZoomj < 0)  
    {
      BegZoomj = 0;
      EndZoomj = NcVisu/2 - 1;
    }

    if (EndZoomi >= BV.nl())
    {
       EndZoomi = BV.nl() - 1;
       BegZoomi = BV.nl() - 1 - HalfZoomSizei*2;
    }
    
    if (EndZoomj >= BV.nc())
    {
       EndZoomj = BV.nc() - 1;
       BegZoomj = BV.nc() - 1 - HalfZoomSizej*2;
    }
 
    if (BegZoomi < 0) BegZoomi = 0;
    if (BegZoomj < 0) BegZoomj = 0;
  //cout << "Zoomi  = " <<   Zoomi << " Zoomj = " << Zoomj << endl;
  //cout << "ZNL: " << BegZoomi << " --> " << EndZoomi << endl;
  //cout << "ZNC: " << BegZoomj << " --> " << EndZoomj << endl;
 
    if (Nl > NlVisu)
    {
        K = BV.fb_nc();
        BB = BV.fb_nl();
        BE = BV.lb_nl();
        BB_NEW = BB;
        while (BegZoomi > BV.block_nl(BB_NEW,K)+BV.pos_nl(BB_NEW,K)) BB_NEW++;
   
        BE_NEW = BB_NEW;
        while ((BE_NEW < BV.lb_nl()) &&
           ((BV.pos_nl(BE_NEW,K) + BV.block_nl(BE_NEW,K)  < EndZoomi )  
	  || (BV.pos_nl(BE_NEW,K) + BV.block_nl(BE_NEW,K) < NlVisu/2)))
        {
           //cout << EndZoomi << BV.pos_nl(BE_NEW,K) + BV.block_nl(BE_NEW,K) <<
	   //         BV.pos_nl(BE_NEW,K) + BV.block_nl(BE_NEW,K) << endl;
           BE_NEW ++;
	}
    
        FirstBlockNl = BB_NEW;
        LastBlockNl  = BE_NEW;
        if (Verbose == True)
	 cout << "NL: New first block   = " <<   BB_NEW << " New last block = " << BE_NEW << endl;
    }

    if (Nc > NcVisu)
    {
        K = BV.fb_nl();
        BB = BV.fb_nc();
        BE = BV.lb_nc();
        BB_NEW = BB;
        while (BegZoomj >  BV.block_nc(K,BB_NEW)+BV.pos_nc(K,BB_NEW)) BB_NEW++;
        BE_NEW = BB_NEW;
        while ((BE_NEW < BV.lb_nc()) &&
           ((BV.pos_nc(K,BE_NEW) + BV.block_nc(K,BE_NEW)  < EndZoomj )  
	  || (BV.pos_nc(K,BE_NEW) + BV.block_nl(BE_NEW,K) < NcVisu/2)))
        {
           //cout << EndZoomi << BV.pos_nc(K,BE_NEW) + BV.block_nc(K,BE_NEW) <<
	   //         BV.pos_nc(K,BE_NEW) + BV.block_nc(K,BE_NEW) << endl;
           BE_NEW ++;
	}        
        FirstBlockNc = BB_NEW;
        LastBlockNc  = BE_NEW ;  
        if (Verbose == True)
	 cout << "NC: New first block   = " <<   BB_NEW << " New last block = " << BE_NEW << endl;
    }
    set_pos();
   }

}

/***************************************************/

int BlockMRInfo::check_size(int NbrScale, int MaxSizeLastScale)
{
   int Nm = MIN(TabNl(0,0), TabNc(0,0));
   int IndMini = 0;
   int IndMinj = 0;
   for (int i=0; i < NbrBlockNl; i++)
   for (int j=0; j < NbrBlockNc; j++)
   {
      int Nmin = MIN(TabNl(i,j), TabNc(i,j));
      if (Nmin < Nm) 
      {
         IndMini = i;
	 IndMinj = j;
	 Nm = Nmin;
      }
   }   
   int ScaleMax=iround(log((float)Nm/(float) MaxSizeLastScale) / log(2.)+ 1.);
   if (NbrScale > ScaleMax)
   {
      cerr << endl << endl;
      cerr << "Error: the number of scales is too high " << endl;
      cerr << "       Block number " << IndMini << IndMinj  << " size = " 
             << TabNl(IndMini,IndMinj) << "x" << TabNc (IndMini,IndMinj) << endl;
      cerr << "       the maximum allowed number of scales is " << ScaleMax << endl;
      if (ScaleMax < 2)
             cerr << "       the block size must be modified "  << endl;
      else
             cerr << "       the block size or the number of scales must be modified "  << endl;
      exit(-10);
   }
   return 0;
}


/***************************************************/

void MR_CompData::upblock(Ifloat &Low)
{
   int Nls = Low.nl();
   int Ncs = Low.nc();
   int i,j,Nls1,Ncs1;
   long  BufSize=0;
   int NbrSkipResol, NbrBand, NbrCoef;
   int *data;
   float Level;    
   Bool VerboseUp = False;
   
    if (NbrBlock > 1) BufSize = readint(InFile); 

    // read the original block image size
    Nl =readint( InFile);  
    Nc =readint( InFile);  
    Nbr_Plan=readint( InFile);    
   
            
    // calculate the image size from Nbr_Plan,Nl,Nc,Transform => Nls1,Ncs1
    int s= 0;
    while ((s < Nbr_Plan) && (Nls != size_ima_resol(s, Nl))
                  && (Ncs != size_ima_resol(s, Nc))) s++;
    if (s == Nbr_Plan)
    {
       cerr << "Error: the input image does not correspond to a given scale ... " << endl;
       exit(-1);
    }         
    s--; 
    Resol = s;
    if (VerboseUp == True)
    {
          cerr << "Resolution " << s+1 << " ==> " << s << endl;
	  if (NbrBlock > 1)  cout << "BufSize = " << BufSize << endl;
	  cout << "Nl = " << Nl << " Nc = " << Nc << " Np = " <<  Nbr_Plan << endl;
    }
    	  
    // extract multiresolution coefficients if the given scale s is > 0 
    if (s >= 0)
    {
       Nls1 = size_ima_resol(s, Nl);
       Ncs1 = size_ima_resol(s, Nc);
    
       // create the multiresolution object with 2 resolutions
        MultiResol MR_Data (Nls1, Ncs1, 2, Transform, "MR_Imag");
        NbrBand = MR_Data.nbr_band_per_resol(); 
       MR_Data.LiftingTrans=LiftingTrans;
       
       // transform the input low resolution image if needed
       if ((Stat_Noise == NOISE_POISSON) && 
             (Comp_Method != COMP_O_MIN_MAX) &&
             (Comp_Method != COMP_INT_LIFTING))
                        noise_poisson_transform (Low,Low); 
                         
       // Haar renormalization
       if (Comp_Method == COMP_HAAR)
       {
          float Coef = 1.;
          for (i=0; i <= s; i++) Coef *= 2.;
          for (i=0; i < Nls*Ncs; i++)   Low(i) *= Coef;
       }
       
       // insert the low resolution image in the multiresolution object       
       MR_Data.band(MR_Data.nbr_band()-1) = Low;
     
       // allocate output image
       Dat.alloc(Nls1, Ncs1,"High resol");
    
       // calculate the number of resolution to skip
       NbrSkipResol = Nbr_Plan  - s - 1;
    
       // go to the good resolution
       seek_resol(NbrSkipResol);
       
       // read the bands corresponding to the good resolution
       float *LevelQ = new float[NbrBand];
   
       for (int b = NbrBand-1; b  >= 0; b--) 
       {
          NbrCoef = MR_Data.nbr_coeff_in_band(b);
          
          // decode the multiresolution coefficients
          decode (InFile, &data, &Nls, &Ncs, &Level);
          LevelQ[b] = Level;
          if (VerboseUp == True)
          {
             extern long size_enc;
             cerr << "Set " << b+1 << " " << Nls << "x" << Ncs << 
                 " Quantif = " << Level  <<  " Size = " <<  size_enc << endl;
          }
          
          
          if (Level > FLOAT_EPSILON)
          {
            for (i = 0; i <  Nls; i++)
            for (j = 0; j <  Ncs; j++)
                              MR_Data(b,i,j) = ((float) data[i*Ncs+j]) * Level;
          }
          else 
          {
             for (i = 0; i <  Nls; i++)
             for (j = 0; j <  Ncs; j++)
                      MR_Data(b,i,j) =  data[0];
          }
          i_free(data);
       }
       
       // reconstruct the high resolution image
       if ((IterRec > 1) && (MR_Data.Set_Transform == TRANSF_MALLAT))
               mr_ortho_regul_rec (MR_Data, Dat, LevelQ, IterRec);
       else    MR_Data.recons (Dat);
       delete [] LevelQ;
       
       // Haar renormalization
       if (Comp_Method == COMP_HAAR)
       {
          float Coef = 1.;
          for (i=0; i < s; i++) Coef *= 2.;
          for (i=0; i <  Dat.nl()* Dat.nc(); i++)    Dat(i) /= Coef;
       }            
          
       // inverse transform    
       if ((Stat_Noise == NOISE_POISSON) && 
             (Comp_Method != COMP_O_MIN_MAX) &&
	     (Comp_Method !=COMP_INT_LIFTING))
              noise_inverse_poisson_transform (Dat, Dat);
	             
       //skip the noise
          if (Resol > 0)
	  {  
	     int NB_seek = number_band_per_resol(Transform)* Resol;
	     // cout << " NB_seek = " << NB_seek << endl;
	     seek_band(NB_seek);
	  }

        // cout << "decode noise" << endl;          
          if (NoiseThresholding == True)
	  {
	      // read info parameters
	      read_noise_info();
 	      if (Resol < 0)  
              {
                  decode_residu();
                  if ((IntTrans == True) && (IResidual != NULL)) i_free(IResidual);
              }
              else if (KeepResi)
	      {
                 // cout << "fseek noise" << endl;
		 seek_band(1);
	      }
	  }
    }
    else   // the scale < 0 ==> we want to add the residu 
    {
       // output image has the same size as the input image
       Nls1 = Nl;
       Ncs1 = Nc;
       // test if the input image size is compatible with the size of
       // compressed image
       if ((Nls1 != Nls) && (Ncs1 != Ncs))
       {
          cerr << "Error: the input image does not correspond to compressed image size " << endl;
          exit(-1); 
       }
       
       // allocate output image
       Dat.alloc(Nls1, Ncs1, "ImagHighResol");
    
       // calculate the number of resolution to skip
       NbrSkipResol = Nbr_Plan;
    
       // go to the good resolution
       seek_resol(NbrSkipResol);
       
       if (NoiseThresholding == True)
       {
          // get the residu
          read_noise_info();
          if (VerboseUp == True)   cerr << "Noise ima = " << Noise_Ima << endl;

          // test if the residu is stored and decode it
          if  (KeepResi!=0)
          {        
              decode(InFile, &data, &Nls, &Ncs, &Level);
              if (Verbose == True)   cerr << "Sigma residu = " << Level << endl;
              for (j = 0; j < Nls*Ncs; j++) 
                         Dat (j) = Low(j) + data[j] * Level;
          }
          else 
          {
          cerr << "Warning: residual part cannot be read ... " << endl;
          cerr << "         output image is set to the input image" << endl;
          Dat = Low;
          }
       }
       else 
       {
          cerr << "Warning: residual part cannot be read ... " << endl;
          cerr << "         output image is set to the input image" << endl;
          Dat = Low;
       }
    }
}
   
   
/***************************************************/

void MR_CompData::up_mrresol(Ifloat &ImagLowResol, BlockMRInfo &BV, BlockMRInfo &BVup)
{
    int i, j;
    int Nls,Ncs;
    Nls = ImagLowResol.nl();
    Ncs = ImagLowResol.nc();
    Ifloat HighResol(BVup.nl(),BVup.nc(), "High");
     
    // if NoBscale is True, then we need to rescale to input data
    if (NoBscale == True)
    {
        for (i=0; i <  ImagLowResol.nl()* ImagLowResol.nc(); i++)
           ImagLowResol(i) =  (ImagLowResol(i) - Header.bzero) / Header.bscale;
    }                    
                          
    int NumBlock=BVup.fb_nl()*BVup.nbr_block_nc()+BVup.fb_nc();
    seek_block(NumBlock);
    NumBlock =  NbrBlocNc - (BVup.lb_nc() - BVup.fb_nc() + 1);
    Ifloat Low;
        
    // the decompression is performed using float
    IntTrans = False;
    
    // test compression method: must be based on multiresolution.
    if ((Comp_Method == COMP_MORPHO) || (Comp_Method == COMP_O_MORPHO))
    {
       cerr << "Error: this compression method is not multiresolution based. " << endl;
       exit(-1);
    }
    
    if ((BVup.nbr_block_nl() != NbrBlocNl) 
                   || ( BVup.nbr_block_nc() != NbrBlocNc))
    {
	cerr << "Error: Block nl and nc in info file not compatible with MRC File ... " << endl;
	cerr << "  NbrBlocNl = " << NbrBlocNl << " NbrBlocNc = " << NbrBlocNc << endl;
	cerr << "  nbr_block_nl= " << BVup.nbr_block_nl() << "  nbr_block_nc = " << BVup.nbr_block_nc() << endl;
	exit(-1);
    }
    if (BVup.block_size() != BlockSize)
    {
	cerr << "Error: Block size in info file not compatible with MRC File ... " << endl;
	exit(-1);
    }
	   
    // seek to the good block number

    for (int Bi=BVup.fb_nl(); Bi <= BVup.lb_nl(); Bi++)
    {
      for (int Bj=BVup.fb_nc(); Bj <= BVup.lb_nc(); Bj++)
      {
        if (Verbose == True)
	{
	   cout << "Block " << Bi << " " << Bj << endl;
	}
        int Posi, Posj;
        // extract a block from the low resolution image
        Nls = BV.block_nl(Bi,Bj);
        Ncs = BV.block_nc(Bi,Bj);
        Low.alloc(Nls, Ncs,"block");
        Posi = BV.pos_nl(Bi,Bj);
        Posj = BV.pos_nc(Bi,Bj);
        for (i = 0; i < Nls; i++)
        for (j = 0; j < Ncs; j++) Low(i,j) = ImagLowResol(Posi+i,Posj+j);
	//if ((Bi == 1) && (Bj==1)) io_write_ima_float("bl.fits", Low);
        upblock(Low);
        //if ((Bi == 1) && (Bj==1)) io_write_ima_float("bh.fits", Dat);

        
        // put the decompressed block in the high resolution image
        Posi = BVup.pos_nl(Bi,Bj);
        Posj = BVup.pos_nc(Bi,Bj);
        Nls = BVup.block_nl(Bi,Bj);
        Ncs = BVup.block_nc(Bi,Bj);
	//cout << "Nl = " << Dat.nl() << " Nc = " << Dat.nc() << endl;
 	//cout << "Pos i = " << Posi << " Pos j = " << Posj << endl;
	//cout << "Nls = " << Nls << " Ncs = " << Ncs << endl;
	if ((Nls != Dat.nl()) || (Ncs != Dat.nc()))
	{
	   cerr << "Error: Block size in information file not compatible with MRC File ... " << endl;
	   exit(-1);
        }
        for (i = 0; i < Nls; i++)
        for (j = 0; j < Ncs; j++) HighResol(Posi+i,Posj+j) = Dat(i,j);
       }
       if (Bi != BVup.lb_nl())
       {
           // skip the first blocks
           i =0;
           while ((i < NumBlock) && (i < NbrBlock)) 
           {
              long BufSize = readint(InFile);
	      if (Verbose == True) cerr << " BufSize block = " << BufSize << endl;
              if ( fseek(InFile, BufSize, SEEK_CUR) < 0)
              {
	         cerr << "Error in fseek: request shift is " <<  BufSize << endl;
                 exit(-1);
              }       
	      i++;
           }      
       }
     }
    // io_write_ima_float("bb.fits",  HighResol);
     Dat = HighResol;
     
     if (ReadHd > 0) 
     {
         Header.width = Dat.nc();
         Header.height = Dat.nl();
         Header.npix = Dat.n_elem();
	 BlockMRInfo BI( L_Nl, L_Nc, BlockSize);
         int Nx = BI.pos_nc(BVup.fb_nl(), BVup.fb_nc());
	 int Ny = BI.pos_nl(BVup.fb_nl(), BVup.fb_nc());
         if (Resol > 0)  
         {
          for (i = 0; i < Resol; i++)
          {
              Nx = (Nx+1) / 2;
              Ny = (Ny+1) / 2;
          } 
         }
         Header.crpixx -= Nx;
         Header.crpixy -= Ny;  	 
     }
    if (InFile != NULL) fclose (InFile);
    dec_data();
}


/***************************************************/
// increase the resolution of a single block
void MR_CompData::up_mrresol(Ifloat &ImagLowResol, int Nb)
{
    int i, NumBlock=Nb-1;
    long  BufSize;
    int Nls,Ncs;
    Nls = ImagLowResol.nl();
    Ncs = ImagLowResol.nc();

    
    // if NoBscale is True, then we need to rescale to input data
    if (NoBscale == True)
    {
        for (i=0; i <  ImagLowResol.nl()* ImagLowResol.nc(); i++)
           ImagLowResol(i) =  (ImagLowResol(i) - Header.bzero) / Header.bscale;
    }                    
                          
    // open the file and read the header
    // init_decompress();
    
    // the decompression is performed using float
    IntTrans = False;
    
    // test compression method: must be based on multiresolution.
    if ((Comp_Method == COMP_MORPHO) || (Comp_Method == COMP_O_MORPHO))
    {
       cerr << "Error: this compression method is not multiresolution based. " << endl;
       exit(-1);
    }
    
    // seek to the good block number
    seek_block(NumBlock);
    
    if (UseBlock == True)   // (NbrBlock > 1)
               BufSize = readint(InFile); 

    upblock(ImagLowResol);

    // correct the fits header
    if (Resol > 0)  header_resol();
 
    // Modify the fits header 
    if (NbrBlock > 1)
    {
       int Nx = TabBlockDepNc[NumBlock];
       int Ny = TabBlockDepNl[NumBlock];
       if ((Resol > 0) && (ReadHd > 0))
       {
          for (i = 0; i < Resol; i++)
          {
              Nx = (Nx+1) / 2;
              Ny = (Ny+1) / 2;
          } 
       }
       Header.crpixx -= Nx;
       Header.crpixy -= Ny;       
    }    
    dec_data();
     
    if (InFile != NULL) fclose (InFile);
}


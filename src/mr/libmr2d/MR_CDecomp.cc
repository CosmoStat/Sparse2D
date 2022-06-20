
// pyramidal median transform with integer
// #include "NR.h"
#include "IM_Obj.h"
#include "IM_Math.h"
#include "IM_IO.h"
#include "IM_BIO.h"
#include "MR_Obj.h"
#include "IM_CompTool.h"
#include "IM_Noise.h"
#include "MR_Noise.h"
#include "MR_Comp.h"
#include "IM_Lut.h"
#include "IM_Comp.h"

extern float PasCodeur;  /* CCD gain */
extern float SigmaGauss; /* CCD read-out noise standard deviation */
extern float MeanGauss;  /* CCD read-out noise mean */
extern float BadPixalVal; 
extern Bool BadPixel;

Bool InputFromStdin = False;
static Iint ImagInt;
static IOInfoData InfoData;

/***************************************************/
// Decompression part
/***************************************************/

void MR_CompData::read_tinfo()
{
     extern softinfo Soft;
   
    short ValIdent[3];
    fread(ValIdent, sizeof(short), 3, InFile);
//    cout << "V1 = " << ValIdent[0] << endl;
//    cout << "V2 = " << ValIdent[1] << endl;
//    cout << "V3 = " << ValIdent[2] << endl;
    
    float SoftVersion;
    SoftVersion = readfloat (InFile);
    if (SoftVersion > Soft.release())
    {
       cerr << "Warning: the image has been compressed with a more recent version " << endl;
       cerr << "         than the decompression program version. " << endl;
       cerr << "         You should update the decompression program." << endl;
    }
    Comp_Method  = (type_comp) readint (InFile);
    switch (Comp_Method)
    {       
       case COMP_HAAR: Transform = TO_HAAR; break;
       case COMP_MALLAT: 
                Transform = TO_MALLAT;
		if (SoftVersion > 2)
		{
                   type_sb_filter TF;
	           sb_type_norm TN;
	           TF = (type_sb_filter) readint(InFile);
		   TN = (sb_type_norm) readint(InFile);
		   if (readint(InFile) == 1) UseLiftingInsteadOfFilter = True;
		   else UseLiftingInsteadOfFilter = False;
		   if (readint(InFile) == 1) UnifQuant = True;
		   else UnifQuant = False;
 		}
		break;
       case COMP_FEAUVEAU: Transform = TO_FEAUVEAU;break;
       case COMP_MIN_MAX: Transform = TM_MIN_MAX; break;
       case COMP_MRMEDIAN: Transform = TM_PYR_MEDIAN; break;
       case COMP_WT_PMT: Transform = TM_TO_PYR; break;
       case COMP_O_MRMEDIAN: Transform = TM_PYR_MEDIAN; break;
       case COMP_O_MIN_MAX: Transform = TM_MIN_MAX;break;
       case COMP_O_MORPHO:
       case COMP_MORPHO: break;
       case COMP_INT_LIFTING: 
             LiftingTrans = (type_lift) readint (InFile);
	     Transform = TO_LIFTING;
             break;
       default: break;
    }  
    Nbr_Plan = readint (InFile);
    if (SoftVersion > 1.)
    {
       L_Nl = readint (InFile);
       L_Nc = readint (InFile);
       KeepResi = readint (InFile);
    }
    if (Resol >= Nbr_Plan) Resol = Nbr_Plan-1;
    
    IntTrans = work_with_int(Comp_Method);     
    NoiseThresholding = noise_threshold(Comp_Method);
        
    FormatInput = (type_format) readint (InFile); 
    // test for compatibility software version     
    if (SoftVersion > 1.)
    {
       int NCol = readint (InFile);
       if (NCol > 0)
       {
      cbyte *r = new cbyte[NCol];
	  cbyte *g = new cbyte[NCol];
      cbyte *b = new cbyte[NCol];
	  fread(r, 1, NCol, InFile);
	  fread(g, 1, NCol, InFile);
	  fread(b, 1, NCol, InFile);
	  IO_RGB.alloc(NCol,r,g,b);
       }
    }
    NbrImag = readint (InFile);
    if (NbrImag > 1)
    {
        cout << "Error: multichannel data cannot be decompressed by this program ... " << endl;
        exit(-1);
    }
    TypeQuant = (type_quant) readint (InFile);
    TypeCoder  = (type_coding) readint (InFile);
    
    Stat_Noise = (type_noise) readint (InFile);
    if (Stat_Noise == NOISE_POISSON)
    {
        PasCodeur = readfloat(InFile);
        SigmaGauss = readfloat(InFile);
        MeanGauss = readfloat(InFile);
    }
    
    if (readint (InFile) == 1)   
    {
       UseBlock  = True;
       // NbrImag = readint (InFile);
       NbrBlock = readint (InFile);
       BlockSize = readint (InFile);
       BlockOverlap= readint (InFile);
    }
     
    if (readint (InFile) == 1)
    { 
       BadPixel = True;
       BadPixalVal = readfloat (InFile);
       MinDat = readfloat (InFile);
       MaxDat = readfloat (InFile);
    }
    
    if (Verbose == True) 
    {
        cerr << endl;
        cerr << Soft.banner() << endl;
        cerr <<  StringComp(Comp_Method) << endl << endl;
        cerr << "Input File = " <<  File_Name_Transform << endl;
        cerr << "Output File = " << File_Name_Imag << endl;
        if (NbrBlock > 1) cerr <<  "Number of blocks = " << NbrBlock << endl;
        if (NbrBlock > 1) cerr <<  "block size = " << BlockSize << endl;
        if (NbrBlock > 1) cerr <<  "block overlapping = " << BlockOverlap<< endl;
        if (Resol >= 0) cerr << "Resolution = " << Resol << endl;
        if (BadPixel == True)
        {
            cerr << "MinDat = " << MinDat << endl;
            cerr << "MaxDat = " << MaxDat << endl;
        }
    }        
}
/*****************************************************************/

void MR_CompData::header_resol()
{
  if (Resol > 0)
  {
       for (int k = 0; k < Resol; k++)
       {
            Header.cdeltx *= 2.;
            Header.cdelty *= 2.;
            Header.crpixx /= 2.;
            Header.crpixy /= 2.;
            Header.width = (Header.width+1)/2;
            Header.height = (Header.height+1)/2;
        } 
        Header.npix = Header.width*Header.height;   
    }
}

/*****************************************************************/

void MR_CompData::read_header()
{
  initfield(&Header);

  ReadHd = readint (InFile);
    
    // cout << "ReadHd = " << ReadHd << endl;
    if (ReadHd > 0)
    {
        int Sname;
        char *Buff;
        char *Buff_Head;

        if (ReadHd > 1)
        {
            Buff_Head = new char [ReadHd];
            fread((char *) Buff_Head, sizeof(char), ReadHd, InFile);
            Header.fitshead = Buff_Head;
            Header.fitsheadsize = ReadHd;
        }
        else Header.fitshead = NULL;

        Sname = readint (InFile);
        Buff = new char [Sname+1];
        fread(Buff, sizeof(char), Sname, InFile);
        Buff[Sname] = '\0';
        Header.origin = Buff;
        strcpy(Header.rident, Cmd);

        Header.bitpix = readint (InFile);
        Header.bytepix = readint (InFile);
        Header.naxis = 2;
        Header.width = readint (InFile);
        Header.height = readint (InFile);
        Header.crpixx = (double) readfloat (InFile);
        Header.crpixy = (double) readfloat (InFile);
        Header.crvalx = (double) readfloat (InFile);
        Header.crvaly = (double) readfloat (InFile);
        Header.cdeltx = (double) readfloat (InFile);
        Header.cdelty = (double) readfloat (InFile);
        Header.crotax = (double) readfloat (InFile);
        Header.crotay = (double) readfloat (InFile);
        Header.bscale = (double) readfloat (InFile);
        Header.bzero = (double) readfloat (InFile);
        Header.ngamma = (double) readfloat (InFile);
        Header.epoch = (double) readfloat (InFile);
        Header.pixscale = (double) readfloat (InFile);
        strcpy(Header.ident, Header.origin);
        Sname = readint (InFile);
        Buff = new char [Sname+1];
        fread(Buff, sizeof(char), Sname, InFile);
        Buff[Sname] = '\0';
        strcpy(Header.ctypex, Buff);
        Sname = readint (InFile);
        Buff = new char [Sname+1];
        fread(Buff, sizeof(char), Sname, InFile);
        Buff[Sname] = '\0';
        strcpy(Header.ctypey,Buff);
        header_resol();
        Header.npix = Header.width*Header.height;   
        if (Verbose == True) 
        {
      cerr << "cmd = " <<  Header.origin << endl;
      cerr << "bitpix = " <<    Header.bitpix;
      cerr << " bytepix = " <<    Header.bytepix << endl;
      cerr << "width = " <<    Header.width;
      cerr << " height = " <<    Header.height << endl;
      if ((Header.crpixx > 0) && (Header.crpixy > 0))
      {
         cerr << "crpixx = " <<    Header.crpixx;
         cerr << " crpixy = " <<    Header.crpixy << endl;
      }
      if ((Header.crvalx > 0) && (Header.crvaly> 0))
      {
         cerr << "crvalx = " <<    Header.crvalx;
         cerr << " crvaly = " <<    Header.crvaly << endl;
      }
      if ((Header.cdeltx > 0) && (Header.cdelty > 0))
      {
         cerr << "cdeltx = " <<    Header.cdeltx;
         cerr << " cdelty = " <<   Header.cdelty << endl;
      }
      if ((Header.crotax> 0) || (Header.crotay > 0))
      {
         cerr << "crotax = " <<    Header.crotax;
         cerr << " crotay = " <<    Header.crotay << endl;
      }
      
      cerr << "bscale = " <<    Header.bscale;
      cerr << " bzero = " <<    Header.bzero << endl;
      if (Header.ngamma> 0)
         cerr << "ngamma = " <<    Header.ngamma << endl;
      if (Header.epoch> 0)
         cerr << "epoch = " <<    Header.epoch << endl;
      if (Header.pixscale> 0)
         cerr << "pixscale = " <<    Header.pixscale << endl;
      if (strlen(Header.ctypex) > 0)
          cerr << "ctypex = " <<    Header.ctypex;
      if (strlen(Header.ctypey) > 0)   
          cerr << " ctypey = " <<    Header.ctypex << endl << endl;
        }
    }

    if (IntTrans == True)
    {
       if (readint (InFile) == 0) NoBscale= False;
       else NoBscale = True;
    }
    if ((Comp_Method == COMP_O_MORPHO) || (Comp_Method == COMP_MORPHO))
    {
       Elstr_Size = readint (InFile);
       Elstr_Shape = readint (InFile);
    }
  
    if (Verbose == True) 
    {
      if (NoBscale == True) cerr << "Bscaling not applied " << endl;
      if ((Comp_Method == COMP_O_MORPHO) || (Comp_Method == COMP_MORPHO))
      {
         cerr << "Struct. element size = " << Elstr_Size << endl;
         if (Elstr_Shape == 0) cerr << "Struct. element = square " <<  endl;
         else cerr << "Struct. element = circle " <<  endl;
      }
    }


   // modify the fits structure following user request
   type_format FormatData = io_which_format(File_Name_Imag);
   if (FormatData == F_UNKNOWN)
   {  
      if (FormatInput != F_UNKNOWN) io_set_format(FormatInput);
      else FormatInput = FormatData;
   }

   // in case of fits format, set the output type of data
   if ((OutputType == 'i') && (ReadHd > 0))
   {
       Header.bitpix = 32;
       Header.bytepix = 4;
   }
   else if ((OutputType == 's') && (ReadHd > 0))
   {
      Header.bitpix = 16;
      Header.bytepix = 2;
   }
   else if ((OutputType == 'f') && (ReadHd > 0))
   {
      Header.bitpix = -32;
      Header.bytepix = 4;
   }    
}

/****************************************************************/

void MR_CompData::decode_data()
{
   int i, NScale,Nli,Nci;
   // cout << "DEC" << Comp_Method << endl;
    switch (Comp_Method)
    {
       case COMP_HAAR:  
       case COMP_MALLAT:
       case COMP_FEAUVEAU:
       case COMP_MIN_MAX:
       case COMP_WT_PMT:
          // if (Comp_Method != COMP_HAAR) IterRec = 0;
	  if ((Comp_Method == COMP_MALLAT)
	         && (UseLiftingInsteadOfFilter == True))
		 lifting_decompress(*BDat, Resol,InFile,Nli,Nci,Verbose);
          else
	  {
          if (Resol == -1)   
	       mr_buddecompress(*BDat, Transform, 0,InFile,Nli,Nci,IterRec,Verbose);
          else mr_buddecompress(*BDat, Transform, Resol,InFile,Nli,Nci,IterRec,Verbose);
          Nl =  BDat->nl();
          Nc =  BDat->nc(); 
             // Haar scaling function is not normalized by 1.
             // It means that pixels values of the image at a lower 
             // resolution are higher than the original image
             // ===> we do a renormalization
             if ((Resol > 0) && (Comp_Method == COMP_HAAR))
             {
                float Coef = 1.;
                for (i=0; i < Resol; i++) Coef *= 2.;
                for (i=0; i < Dat.nl()*Dat.nc(); i++) (*BDat)(i) /= Coef;
             }
	  }
          break;
       case COMP_MRMEDIAN:
             if (Resol == -1) mr_decompress (*BDat,0,InFile,Nli,Nci,Verbose);
             else mr_decompress (*BDat,Resol,InFile,Nli,Nci,Verbose);
             Nl =  BDat->nl();
             Nc =  BDat->nc();
	     break;
       case COMP_O_MRMEDIAN:
             if (Resol == -1) 
                      mri_decompress (&BIDat,&Nl,&Nc,& NScale,0,InFile,Nli,Nci,Verbose);
             else mri_decompress (&BIDat,&Nl,&Nc,&NScale, Resol,InFile,Nli,Nci,Verbose);
	     ImagInt.alloc(BIDat,Nl,Nc);
             break;
       case COMP_MORPHO:
             morpho_decompress (*BDat, N_Iter,InFile,
                                         Elstr_Size,Elstr_Shape,Verbose);
             Nl = Dat.nl();
             Nc = Dat.nc();
	     Nli = Nl;
	     Nci = Nc;
            break;
       case COMP_O_MIN_MAX:
            {
 	     Iint GDat;
             if (Resol == -1)  
	          mr_i_ortho_decompress (GDat,0, NScale, InFile,Nli,Nci,Verbose);
             else mr_i_ortho_decompress (GDat, Resol, NScale, InFile,Nli,Nci,Verbose);
             Nl = GDat.nl();
             Nc = GDat.nc();
  	     /* image reconstruction */
	     ImagInt.alloc(Nl,Nc,"rec");
             inverse_transform_g (GDat, ImagInt, NScale-1);
	     BIDat = ImagInt.buffer();
             }
            break;  
       case COMP_INT_LIFTING:
            {
	     // cout << "lift" << endl;
             if (Resol == -1)  
	          mr_i_ortho_decompress (ImagInt,0,  NScale, InFile, Nli,Nci,Verbose);
             else mr_i_ortho_decompress (ImagInt, Resol, NScale, InFile, Nli,Nci,Verbose);
             Nl = ImagInt.nl();
             Nc = ImagInt.nc();
	     Lifting Clift1D(LiftingTrans);
             Ortho_2D_WT WT2D(Clift1D);
             WT2D.recons(ImagInt, NScale);
	     BIDat = ImagInt.buffer();
	    }
             break;        
       case COMP_O_MORPHO:
             Nl =readint(InFile); /* x size of image */
             Nc =readint(InFile); /* y size of image */
             ImagInt.alloc(Nl,Nc,"ImagInt");
             BIDat= ImagInt.buffer();
             morphoi_decompress (BIDat, N_Iter, Nl, Nc,
                                          InFile,Elstr_Size,Elstr_Shape,Verbose);
             Nli = Nl;
	     Nci = Nc;
             break;
       default:
             cerr << "Error: unknown compression method ... " << endl;
             exit(-1);
             break;
    } 

    // inverse transform    
    if ((Stat_Noise == NOISE_POISSON) && 
            (Comp_Method != COMP_O_MIN_MAX) &&
            (Comp_Method != COMP_INT_LIFTING))
    {
        if (IntTrans == False)
             noise_inverse_poisson_transform (*BDat, *BDat);
        else noise_inverse_poisson_transform (BIDat, BIDat, Nl, Nc);
    }
    
    // if it is not a block compression, the original image size is
    // equal to Nli,Nci
    if (NbrBlock <= 1)
    {
       L_Nl = Nli;
       L_Nc = Nci;
    }
}


/***************************************************************/

void MR_CompData::read_noise_info()
{
 
    KeepResi = readint(InFile);
    Noise_Ima = readfloat(InFile);
}

/*****************************************************************/

void MR_CompData::decode_residu()
{
   float Sigmaf;
   int i,j,k=0;
        
   if (IntTrans == False)
   {
    if ((Resol==-1) && (KeepResi!=0))
    {        
        decode(InFile, &IResidual, &Nl, &Nc, &Sigmaf);

        // if (Verbose == True)   cerr << "Sigma residu = " << Sigmaf << endl;
        for (i = 0; i < Nl; i++)
        for (j = 0; j < Nc; j++) (*BDat)(i,j) += IResidual[k++] * Sigmaf;
    }
    else if (SimuNoise == True) 
         {
            if (Stat_Noise == NOISE_POISSON)
                           noise_poisson_transform (*BDat, *BDat);
            *BDat +=im_noise(Noise_Ima, BDat->nl(), BDat->nc());
            if (Stat_Noise == NOISE_POISSON)
               noise_inverse_poisson_transform (*BDat, *BDat);
        else noise_inverse_poisson_transform (BIDat, BIDat, Nl, Nc);
        }
   }
   else
   {
         if ((Resol==-1) && (KeepResi!=0))
         {        
            decode(InFile, &IResidual, &Nl, &Nc, &Sigmaf);

            // if (Verbose == True)  cerr << "Sigma residu =  " << Sigmaf << endl;
            for (i = 0; i < Nl*Nc; i++) 
                 BIDat[i] += (int) (IResidual[i] * Sigmaf);
         }
         else if (SimuNoise == True) 
         {
             Ifloat DatN(Nl,Nc,"im_noise");
             for (i = 0; i < Nl*Nc; i++)  DatN(i) = BIDat[i];
             if (Stat_Noise == NOISE_POISSON)
                           noise_poisson_transform (DatN, DatN);
             DatN += im_noise(Noise_Ima, Nl, Nc);             
             if (Stat_Noise == NOISE_POISSON)
                     noise_inverse_poisson_transform (DatN, DatN);
             for (i = 0; i < Nl*Nc; i++) BIDat[i] = quant(DatN(i), 1.);
         }
   }
}

/*****************************************************************/

void MR_CompData::put_block(int NumBlock)
{
    int Depi, Depj;
    Depj = TabBlockDepNc[NumBlock];
    Depi = TabBlockDepNl[NumBlock];
    
   if (IntTrans == False)
   { 
    io_write_block_ima(File_Name_Imag, *BDat, Depi, Depj, InfoData, NoBscale);          
   }
   else
   {
    io_write_block_ima(File_Name_Imag, ImagInt, Depi, Depj, InfoData, NoBscale);          
   }
       
}

/*****************************************************************/

void MR_CompData::dec_data()
{
   int i;
   
    // BSCALE operation and send to the output the result              
   if (IntTrans == False)
   {
      if (NoBscale == True)
           for (i=0; i < Dat.nl()*Dat.nc(); i++)
                 Dat(i) = Dat(i) * Header.bscale + Header.bzero;
   }
    else
    {
         Dat.resize(Nl,Nc);
         if (NoBscale == True)
         {
            for (i=0; i < Nl*Nc; i++)
               Dat(i) = IDat[i] * Header.bscale + Header.bzero;
         }
         else for (i=0; i < Nl*Nc; i++)  Dat(i) = IDat[i];
    }
    
    // if bad pixels must be be replaced, do it.
    if ((BadPixel == True) && (Resol < 0))
    {
        for (i=0; i < Dat.nl()*Dat.nc(); i++)
        {
            if   (Dat(i)-Noise_Ima*5 < MinDat)  Dat(i) = BadPixalVal;
            else if (Dat(i)+Noise_Ima*5 > MaxDat) Dat(i) = BadPixalVal;
        }
    }    
}
   
 

   
/*****************************************************************/

void MR_CompData::write()
{
   // write to ouput file the decompressed data
   // cout << "ReadHd = " << ReadHd << endl;
   
   if (ReadHd > 0)  io_write_ima_float(File_Name_Imag, Dat, &Header);
   else io_write_ima_float(File_Name_Imag, Dat);
 }

/*****************************************************************/

void MR_CompData::set_block_info_data()
{
   int i;
   int Nli=L_Nl,Nci=L_Nc;
  // cerr << "set_block_info_data:  L_Nl " << L_Nl  << " L_Nc =  " <<  L_Nc << endl;
  // cerr << "StringTransform = " << StringTransform(Transform) << endl;
  
     if (BlockSize < 0) BlockSize = MAX(Nli,Nci);
   set_tabblock(); 
  
   // find final image size
   if (Resol > 0)    
   {
       // BlockSize  = BlockSize >> Resol;
       int Depi=0;
       int Depj=0; 
       Nli = 0;
       Nci = 0;
       for (i = 0; i < NbrBlock; i++)
       {  
	     TabBlockNl[i] = size_ima_resol(Resol, TabBlockNl[i]);
	     TabBlockNc[i] = size_ima_resol(Resol, TabBlockNc[i]);
             if (i <  NbrBlocNc) Nci += TabBlockNc[i];
             if (i %  NbrBlocNc == 0)
             {
                 Depj = 0;
                 Nli += TabBlockNl[i];
                 if (i!= 0) Depi += TabBlockNl[i-1];
             }
             TabBlockDepNl[i] = Depi;
             TabBlockDepNc[i] = Depj;
             Depj += TabBlockNc[i];
//      cerr << "Block " << i+1 << ": " << TabBlockNl[i] << "x" << TabBlockNc[i] <<
//              " DepNl = " << TabBlockDepNl[i] << " DepNc = " << TabBlockDepNc[i] << endl;
   
          }              
          Header.width = Nci;
          Header.height = Nli;
          Header.npix = Nci*Nli;
   }
   // cout << "FINAL Image Size= " << Nli  << "x" << Nci  <<  endl;
       
   InfoData.Nl = Nli;
   InfoData.Nc = Nci;
   InfoData.Format = FormatInput;
   InfoData.Nima = 1;
   
   if (ReadHd > 0) InfoData.PtrFits = &Header;
   else
   {
       InfoData.PtrFits = new fitsstruct;
       init_fits_struct(InfoData.PtrFits, Nli, Nci);
   }
   if (InfoData.Format == F_GIF) InfoData.PtrGif = new PICINFO;
   InfoData.put_info(File_Name_Imag);
   if (IntTrans == False) BDat = new Ifloat;
}

/*****************************************************************/

void MR_CompData::decompress()
{ 
    int i;
    long BufSize;
    
    // intialize the decompression
    init_decompress();
        
    if ((Resol > 0) && (    (Comp_Method == COMP_MORPHO) 
                         || (Comp_Method == COMP_O_MORPHO)))
    {
       cerr << "Error: morphology based method is not compatible with -r option ... " << endl;
       exit(-1);
       
    }
    
    // decompress all blocks
    if (UseBlock  == True) // (NbrBlock > 1) 
    {
       L_Nl = readint(InFile);
       L_Nc = readint(InFile);  
       if (Verbose == True)
           cerr << "Compressed image size: Nl = " << L_Nl << " Nc = " << L_Nc << endl;
       set_block_info_data();

       // decode each block
       for (i =0; i < NbrBlock; i++)
       {
          BufSize = readint(InFile);
          if (Verbose == True)
          {  
          cerr << "Block " << i+1 << ": " << TabBlockNl[i] << "x" << TabBlockNc[i] <<
                  " Size = " << BufSize << endl;
          }
          decode_data(); 
          
        // cout << "skip to next block" << endl;
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
          
       //  cout << "put_block" << endl; 
          put_block(i);
          if (IntTrans == True) ImagInt.free();
	  else (*BDat).free();
       }
       io_close_block_ima(File_Name_Imag, InfoData);
    }   
    else
    {
        BDat = &Dat;
        decode_data(); 
        IDat = BIDat;
        // decode noise if all scale are decompressed (Resol == -1)
       if ((NoiseThresholding == True) && (Resol < 0))
       {
          read_noise_info();
          decode_residu();
       }
       dec_data();
    }

  if ((NoiseThresholding == True)
     && (Verbose == True) && (Stat_Noise == NOISE_GAUSSIAN) && (Resol < 0) )
    {
        fprintf(stderr, "Sigma estimated noise = %f\n", Noise_Ima);
    }
    
   if ((InFile != NULL) && (InFile != stdin)) fclose (InFile);
}


/*****************************************************************/

void MR_CompData::init_decompress()
{
   if (InputFromStdin == False)
    {
        InFile = fopen (File_Name_Transform, "rb");
        if (InFile == NULL)
        {
         cerr << "Error: cannot open input file: "<<File_Name_Transform<< endl;
         exit(0);
        }
    }
    else InFile = stdin;

    read_tinfo();
    read_header();
}

/*****************************************************************/

void MR_CompData::decompress(int NB)
{ 
    int i, NumBlock=NB-1;
       
    init_decompress();
    long BufSize;

    if ((NumBlock<0) || (NumBlock >= NbrBlock)) NumBlock = 0;
    
    // skip the first blocks
    seek_block(NumBlock);
    
    // read the total size 
    if (UseBlock  == True)   // (NbrBlock > 1) 
      BufSize = readint(InFile);
    
    BDat = &Dat;
    decode_data();
    IDat = BIDat;
    //L_Nl = Nl;
    //L_Nc = Nc;
    
    // Modify the fits header 
    if (NbrBlock > 1)
    {
       Header.width = Nc;
       Header.height = Nl;
       Header.npix = Nc*Nl;
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
          
    if (Verbose == True)
         cerr << "Block " <<  NumBlock+1 << ": " <<  Nl << "x" <<  Nc << endl;
         
    // decode noise
    if ((NoiseThresholding == True) && (Resol < 0))
    {
        read_noise_info();
        decode_residu();
    }
 
   if ((NoiseThresholding == True)
     && (Verbose == True) && (Stat_Noise == NOISE_GAUSSIAN) && (Resol < 0) )
         fprintf(stderr, "Sigma estimated noise = %f\n", Noise_Ima);
   
    dec_data();
    if (InFile != NULL) fclose (InFile);
}

/***************************************************/

void MR_CompData::get_resol(int RN, int &SizeBufResol, char * & BuffResol, int NB)
{
    int i,NumBlock=NB-1;
    int NbrSkipResol;
    long BufSize;
    int ResolNumber=RN;
    init_decompress();
    if (ResolNumber >= Nbr_Plan) ResolNumber = Nbr_Plan-1;
    
    // test compression method: must be based on multiresolution.
    if ((Comp_Method == COMP_MORPHO) || (Comp_Method == COMP_O_MORPHO))
    {
       cerr << "Error: this compression method is not multiresolution based. " << endl;
       exit(-1);
    }    
    seek_block(NumBlock);
    if (UseBlock  == True) // (NbrBlock > 1) 
       BufSize = readint(InFile); 

    // read the original block image size
    Nl =readint( InFile);  
    Nc =readint( InFile);  

    // read the number of scales 
    Nbr_Plan=readint( InFile);    
    
    // calculate the number of resolution to skip
    NbrSkipResol = Nbr_Plan  - ResolNumber - 1;
     
    if (NbrSkipResol >= Nbr_Plan) 
    {
       seek_resol(Nbr_Plan);
       read_noise_info();
       if (Verbose == True)   cerr << "Noise ima = " << Noise_Ima << endl;
       // test if the residu is stored 
       if ((ResolNumber==-1) && (KeepResi!=0)) 
       {
           SizeBufResol  = readint(InFile);
           if (Verbose == True)  cerr << "Size residual = " << SizeBufResol << endl;
           BuffResol = new char [SizeBufResol];
           fread((void *) BuffResol, sizeof(char), SizeBufResol, InFile);
       }
       else 
       {
           cerr << "Error: no residual part in the compressed file ... " << endl;
           exit(-1);
       } 
    }
    else
    {   
        seek_resol(NbrSkipResol);
        
        // read the number of bands (contained in a single resolution) to read 
        int NbrB;
        if (NbrSkipResol == 0) NbrB = 1;
        else NbrB = number_band_per_resol(Transform);
        
        // read the bands
        SizeBufResol = 0;
        long SizeOneBand;
        for (i=0; i < NbrB; i++)
        {
           char *PtrBuff;
           
           SizeOneBand = readint(InFile);
           if (Verbose == True)   
             cerr << "Resol " << ResolNumber <<  " Band " <<  i+1 << " Size = " << SizeOneBand << endl;
           if (SizeBufResol == 0) BuffResol = new char [SizeOneBand];
           else BuffResol = (char *) realloc(BuffResol, SizeBufResol+SizeOneBand);
           PtrBuff = BuffResol;
           PtrBuff += SizeBufResol;
           SizeBufResol += SizeOneBand;
           fread((void *) PtrBuff, sizeof(char), SizeOneBand, InFile);
        }
     }
    
    if (InFile != NULL) fclose (InFile);
}

/***************************************************/
// seek the first N blocks
void MR_CompData::seek_block(int NumBlock)
{
    long BufSize;
    int i;
    
    if (UseBlock  == True) // (NbrBlock > 1) 
    {
       L_Nl = readint(InFile);
       L_Nc = readint(InFile);  
       set_tabblock();
    }
    
    // skip the first blocks
    i =0;
    while ((i < NumBlock) && (i < NbrBlock)) 
    {
       BufSize = readint(InFile);
       if ( fseek(InFile, BufSize, SEEK_CUR) < 0)
       {
	  cerr << "Error in fseek: request shift is " <<  BufSize << endl;
          exit(-1);
       }       i++;
    }      
}

/***************************************************/

void MR_CompData::seek_resol(int NbrSkipResol)   
{
    int NumberBand_per_Resol = number_band_per_resol(Transform);
    int NbrSkipBand = NumberBand_per_Resol*(NbrSkipResol-1)+1;
//     cout << "NumberBand_per_Resol = " << NumberBand_per_Resol << endl;
//     cout << "NbrSkipBand = " <<NbrSkipBand<< endl;
//     cout << "NbrSkipResol= " <<NbrSkipResol<< endl;
    seek_band(NbrSkipBand);
}

/***************************************************/

void MR_CompData::seek_band(int NbrSkipBand)   
{
    long BufSize;
      
     // skip the resolutions
    if (NbrSkipBand >  3*Nbr_Plan+1)  
    {
       cerr << "Error: the number of band to skip is too high ... " << endl;
       cerr << "       NbrSkipBand = " << NbrSkipBand << " Nbr_Plan = " << Nbr_Plan  << endl;
       exit(-1);
    }
      
    for (int r=0; r < NbrSkipBand; r++)
    {
         BufSize = readint(InFile);
         // cout << "size skip " << BufSize << endl;
         if (fseek(InFile, BufSize, SEEK_CUR))
 	 {
	      cerr << "Error in fseek: request shift is " <<  BufSize  << endl;
              exit(-1);
	 }
    }       
}

/***************************************************/


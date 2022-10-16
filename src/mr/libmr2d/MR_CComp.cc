

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
#include "IM_Lut.h"

extern float PasCodeur;  /* CCD gain */
extern float SigmaGauss; /* CCD read-out noise standard deviation */
extern float MeanGauss;  /* CCD read-out noise mean */
extern float BadPixalVal; 
extern Bool BadPixel;

static Iint ImagInt;
static Iint *BImagInt;
static int InputNbByte;

static IOInfoData InfoData;

/***************************************************/


// static float ecart_type_i (int *Imag, int Nl, int Nc)
// {
//     int j;
//     double Moy = 0., Ecart = 0.;
// 
//     for (j = 0; j < Nl * Nc; j++)
//     {
//         Moy += Imag [j];
//         Ecart += Imag [j] * Imag [j];
//     }
//     Moy /= (float) (Nl * Nc);
//     Ecart = Ecart / (float) (Nl * Nc) - Moy*Moy;
//     if (Ecart >= 0.) Ecart = sqrt(Ecart);
//     else Ecart = 0.;
//     return ((float) Ecart);
// }

/***************************************************/

void MR_CompData::encode_residu()
{
    float Sigmaf=0;
    int *Resi;
    int i;

    Resi = i_alloc(Nl*Nc);
    if (IntTrans == False)
    {
       Sigmaf =  Noise_Ima * NoiseQuantif;
       if (KeepResi == 2) Sigmaf = 0.5;

       if (Stat_Noise == NOISE_POISSON)
       {
           noise_inverse_poisson_transform (Residual,Residual);
           noise_inverse_poisson_transform (*BDat, *BDat);
       }
   
       for (i = 0; i < Nl*Nc; i++) 
             Resi[i] = quant ((*BDat)(i) - Residual(i), Sigmaf);
      }
     else
     {  
          Sigmaf = Noise_Ima * NoiseQuantif;
          if (KeepResi == 2) Sigmaf = 0.5;

          if (Stat_Noise == NOISE_POISSON)  
          {
              noise_inverse_poisson_transform (IResidual,IResidual,Nl,Nc);
              noise_inverse_poisson_transform ( BIDat, BIDat,Nl,Nc);
          }
          for (i = 0; i < Nl*Nc; i++) 
               Resi[i] = quant ((float)(BIDat[i] - IResidual[i]), Sigmaf);
    }
    if (Verbose == True) cerr <<  "Residual coding"<< endl;

    encode (OutFile, Resi, Nl, Nc, Sigmaf);
    i_free(Resi);
}


/***************************************************/

void MR_CompData::store_noise_info()
{
 

    writeint (OutFile, KeepResi);
    writefloat(OutFile, Noise_Ima);
}

/***************************************************/


void MR_CompData::encode_data()
{
    Bool Rec = (KeepResi > 0) ? True : False;
     
    InputNbByte = Nl*Nc;
    switch ( which_data_type() )
    {
       case T_SHORT:    InputNbByte *=2;break;
       case T_INT:      InputNbByte *=4;break;
       case T_FLOAT:    InputNbByte *=4;break;
       case T_DOUBLE:   InputNbByte *=8;break;
       default: break;
    }
    if (UseBudget == True)
    {       
       TargetNbrByte =  (int) ( (float)InputNbByte / CompressionRatio + 0.5);
       if (IntTrans == True)
       {
       cerr << "Error: fixed budget compression method must work with float ..." << endl;
       exit(-1);
       }
       if (TargetNbrByte  < 1)
       {
       cerr << "Error: budget amount is not initialized ..." << endl;
       exit(-1);
       }
    }
    else TargetNbrByte = 0;
    
    
     switch (Comp_Method)
     {
        case COMP_MRMEDIAN:
            mr_compress (*BDat, Nbr_Plan, Noise_Ima, N_Sigma, 
                    SignalQuantif, (int)TargetNbrByte, OutFile, Residual,  
                    SupIsol, NoiseInData, MedianWindowSize, SupNeg, Verbose);
          break;
       case COMP_O_MRMEDIAN:
             mri_compress (BIDat, Nl, Nc, Nbr_Plan, &Noise_Ima,  
                           SignalQuantif, N_Sigma,   
                           OutFile, IResidual, SupIsol, NoiseInData,
                           Rec, MedianWindowSize, SupNeg, Verbose);
             break;       
       case COMP_MORPHO:            
              morpho_compress (*BDat, OutFile, Noise_Ima, 
                              N_Sigma, Elstr_Size,Elstr_Shape, 
                              UseMR_for_BGR, NpixBgr, Verbose);
             break;
       case COMP_O_MORPHO:
             morphoi_compress (BIDat, Nl, Nc, OutFile, Noise_Ima, 
                               N_Sigma, Elstr_Size,
                               Elstr_Shape, UseMR_for_BGR, NpixBgr,Verbose);
             break; 
      case COMP_INT_LIFTING: 
            {
             Lifting Clift1D(LiftingTrans);
             Ortho_2D_WT WT2D(Clift1D);
             WT2D.transform(*BImagInt, Nbr_Plan);
	     mr_i_ortho_compress(*BImagInt, Nbr_Plan, OutFile);
            } 
 	     break;
      case COMP_O_MIN_MAX:
            {
               Iint GImag(Nl,Nc,"GDat");
               transform_g (*BImagInt,GImag,Nbr_Plan-1);
               mr_i_ortho_compress(GImag, Nbr_Plan,OutFile);
	    }
            break;
       case COMP_HAAR:
       case COMP_MALLAT:
       case COMP_MIN_MAX:
       case COMP_FEAUVEAU:
       case COMP_WT_PMT:
             if ((Comp_Method == COMP_MALLAT)
	         && (UseLiftingInsteadOfFilter == True))
		 lifting_compress (*BDat, Nbr_Plan, OutFile, N_Sigma,  
		                   SignalQuantif, UnifQuant, Verbose, Rec);
             else
             mr_budcompress (*BDat, Transform, Nbr_Plan, (int) TargetNbrByte,
                            OutFile, Residual, NoiseInData, Noise_Ima, 
                            N_Sigma,  SignalQuantif, SupIsol, True, 
                            SupNeg, MedianWindowSize, MaxIter, Verbose, Rec);
             break;
       default:
             cerr << "Error: unknown compression method ... " << endl;
             exit(-1);
             break;
    }
  
}

/***************************************************/

MR_CompData::~MR_CompData()
{
     if (IntTrans == True)
     {
          if (IResidual != NULL) i_free(IResidual);
          // if (IDat != NULL) delete IDat;
     }
     if (TabBlockNl != NULL) delete  [] TabBlockNl;
     if (TabBlockNc != NULL) delete [] TabBlockNc;
     if (TabBlockDepNl  != NULL) delete [] TabBlockDepNl;
     if (TabBlockDepNc != NULL) delete [] TabBlockDepNc;
}

/***************************************************/

void MR_CompData::variance_stab()
{
    /* Anscombe's transform */
    if (IntTrans == False)
    {
        noise_poisson_transform (*BDat, *BDat);
    }
    else
    {
        noise_poisson_transform (BIDat, BIDat, Nl, Nc);
    } 
    Noise_Ima = 1.;
}

/***************************************************/

void MR_CompData::read_data()
{
    initfield(&Header);
    if (UseBlock == False) 
    {
      if (IntTrans == False)
      {
          io_read_ima_float((char *) File_Name_Imag, Dat, &Header);
          Nl = Dat.nl();
          Nc = Dat.nc();
      }
      else
      {
          io_read_ima_int((char *) File_Name_Imag, ImagInt, &Header, NoBscale);
           Nl = ImagInt.nl();
           Nc = ImagInt.nc();
           IResidual = i_alloc (Nl*Nc);
           IDat = ImagInt.buffer();
      }
      L_Nl = Nl;
      L_Nc = Nc;
      FormatInput = which_format();
    }
    else
    {
       InfoData.get_info((char *) File_Name_Imag);
       Nl = L_Nl = InfoData.Nl;
       Nc = L_Nc = InfoData.Nc;
       FormatInput = InfoData.Format;
       if (FormatInput == F_FITS) Header = *(InfoData.PtrFits);          
    }
}

/***************************************************/

void MR_CompData::replace_badpix()
{
   int i,j;
   float Val;
   
    if (BadPixel == True) 
    {
       if (IntTrans == False)
       {
              for (i=0; i < Nl; i++)
              for (j=0; j < Nc; j++)
              {
                 Val = Dat(i,j);
                 if (ABS(Val-BadPixalVal) > FLOAT_EPSILON)
                 {
                    if (MinDat > Val) MinDat = Val;
                    else if (MaxDat < Val) MaxDat = Val;
                 }
              }
       }
       else 
       {
              for (i=0; i < Nl*Nc; i++)
              {
                 Val = IDat[i];
                 if (ABS(Val-BadPixalVal) > FLOAT_EPSILON)
                 {
                    if (MinDat > Val) MinDat = Val;
                    else if (MaxDat < Val) MaxDat = Val;
                 }
              }
         }
    }
}

/***************************************************/

void MR_CompData::print_info()
{
   extern softinfo Soft;
   cerr << endl;
   cerr <<  Soft.banner() << endl;
   cerr <<  StringComp(Comp_Method) << endl << endl;

   cerr <<  "Input File = "<< File_Name_Imag << endl;
   cerr <<   "Image size (" << Nl << "," << Nc << ") N_Sigma " << N_Sigma << endl;
   cerr <<  "Output File = "<< File_Name_Transform << endl;
#if WRITE_PARAM
    cerr << endl << endl << "PARAMETERS: " << endl << endl;
    cerr << "File Name in = " << File_Name_Imag << endl;
    cerr << "File Name Out = " << File_Name_Transform << endl;
    cerr << "Number of scales = " << Nbr_Plan << endl;
    if (Stat_Noise == NOISE_GAUSSIAN)
                             cerr << "Type of Noise = GAUSSIAN" << endl;
    else
    {
           cerr << "Type of Noise = POISSON" << endl;
           cerr << "  Gain = " << PasCodeur << endl;
           cerr << "  Read-out Noise Sigma  = " << SigmaGauss << endl;
           cerr << "  Read-out Mean = " << MeanGauss << endl;
    }
    cerr << "Sigma Noise = " << Noise_Ima << endl;
    cerr << "N_Sigma = " << N_Sigma << endl;
    // cerr << "Max_Iter = " << MaxIter << endl;
    cerr << "SignalQuantif = " << SignalQuantif << endl;
    if (KeepResi)
      cerr << "Keep residual : NoiseQuantif = " << NoiseQuantif << endl;
    else 
      cerr << "residual not kept : " << endl;
#endif
}


/***************************************************/

void MR_CompData::store_header()
{
    extern softinfo Soft;
    fwrite(ValIdent, sizeof(short), 3, OutFile);

    writefloat (OutFile, (float)  Soft.release());
    writeint (OutFile, (int) Comp_Method);
    if (Comp_Method == COMP_INT_LIFTING) 
         writeint (OutFile, (int) LiftingTrans);

    if (Comp_Method == COMP_MALLAT)
    {
       type_sb_filter TF = F_MALLAT_7_9;
       sb_type_norm TN = NORM_L1;
       if (UseLiftingInsteadOfFilter == True) TN = NORM_L2;
       writeint (OutFile, (int) TF);
       writeint (OutFile, (int) TN);
       if (UseLiftingInsteadOfFilter == True) writeint (OutFile,1);
       else writeint (OutFile,0);
       if (UnifQuant == True) writeint (OutFile,1);
       else writeint (OutFile,0);
    }
    writeint (OutFile, (int) Nbr_Plan);
    writeint (OutFile, (int) L_Nl);
    writeint (OutFile, (int) L_Nc);
    writeint (OutFile, (int) KeepResi);
    writeint (OutFile, (int) FormatInput);
    if (IO_RGB.lut_size() != 0)
    {
       writeint (OutFile, (int) IO_RGB.lut_size());
       fwrite(IO_RGB.red(), 1, IO_RGB.lut_size(), OutFile);
       fwrite(IO_RGB.green(), 1, IO_RGB.lut_size(), OutFile);
       fwrite(IO_RGB.blue(), 1, IO_RGB.lut_size(), OutFile);
    }
    else writeint (OutFile, (int) 0);
    
    writeint (OutFile, (int) NbrImag);
    writeint (OutFile, (int) TypeQuant);
    writeint (OutFile, (int) TypeCoder);
    
    if (Stat_Noise == NOISE_GAUSSIAN) writeint (OutFile, (int) Stat_Noise);
    else 
    {
        writeint (OutFile, (int) Stat_Noise);
        writefloat(OutFile, PasCodeur);
        writefloat(OutFile, SigmaGauss);
        writefloat(OutFile, MeanGauss);
    }
    
    if (UseBlock  == True)
    {
       writeint (OutFile, (int) 1);
       writeint (OutFile, NbrBlock);
       writeint (OutFile, BlockSize);
       writeint (OutFile, BlockOverlap);
    }
    else writeint (OutFile, (int) 0);
    
    if (BadPixel == True)
    {
       writeint (OutFile, (int) 1);
       writefloat (OutFile, BadPixalVal);
       writefloat (OutFile, MinDat);
       writefloat (OutFile, MaxDat);
    }
    else writeint (OutFile, (int) 0);
 
 
    /* Save Header Information */
    if (Header.npix != 0) 
    {
        int Sname = strlen(Cmd);

        if (KeepFitsHeader == 1)
        {
            writeint (OutFile, Header.fitsheadsize);
            fwrite((char *) Header.fitshead, sizeof(char), 
                   Header.fitsheadsize, OutFile);
        }
        else writeint (OutFile, 1);

        writeint (OutFile, Sname);
        fwrite(Cmd, sizeof(char), Sname, OutFile);
        writeint (OutFile, Header.bitpix);
        writeint (OutFile, Header.bytepix);
        writeint (OutFile, Header.width);
        writeint (OutFile, Header.height);
        writefloat (OutFile, (float) Header.crpixx);
        writefloat (OutFile, (float) Header.crpixy);
        writefloat (OutFile, (float) Header.crvalx);
        writefloat (OutFile, (float) Header.crvaly);
        writefloat (OutFile, (float) Header.cdeltx);
        writefloat (OutFile, (float) Header.cdelty);
        writefloat (OutFile, (float) Header.crotax);
        writefloat (OutFile, (float) Header.crotay);
        writefloat (OutFile, (float) Header.bscale);
        writefloat (OutFile, (float) Header.bzero);
        Header.ngamma=0;
        writefloat (OutFile, (float) Header.ngamma);
        writefloat (OutFile, (float) Header.epoch);
        Header.pixscale=0;
        writefloat (OutFile, (float) Header.pixscale);
        Sname = strlen(Header.ctypex);
        writeint (OutFile, Sname);
        fwrite(Header.ctypex, sizeof(char), Sname, OutFile);
        Sname = strlen(Header.ctypey);
        writeint (OutFile, Sname);
        fwrite(Header.ctypey, sizeof(char), Sname, OutFile);
    }
    else writeint (OutFile, 0);

    if (IntTrans == True)
    {
       if (NoBscale == True) writeint (OutFile, 1);
       else writeint (OutFile, 0);
    }
    if ((Comp_Method == COMP_O_MORPHO) || (Comp_Method == COMP_MORPHO))
    {
        writeint (OutFile, Elstr_Size);
        writeint (OutFile, Elstr_Shape);
    }
}


/***************************************************/
 
void MR_CompData::set_tabblock()
{
   NbrBlocNl = L_Nl/BlockSize;
   // if (NbrBlocNl * BlockSize != L_Nl)  NbrBlocNl++;
   NbrBlocNc  = L_Nc/BlockSize;
   // if (NbrBlocNc * BlockSize != L_Nc)  NbrBlocNc++;
   NbrBlock = NbrBlocNl*NbrBlocNc;
   
   TabBlockNl = new int [NbrBlock];
   TabBlockNc = new int [NbrBlock];
   TabBlockDepNl = new int [NbrBlock];
   TabBlockDepNc = new int [NbrBlock];
   
   if (Verbose == True)
   {
   cerr << "BLOCK PARAM" << endl;
   cerr << "L_Nl = " << L_Nl << " L_Nc = " << L_Nc << endl;
   cerr << "NbrBlocNl = " << NbrBlocNl  << " NbrBlocNc = " <<  NbrBlocNc << endl;
   }
   
   for (int i=0; i < NbrBlock; i++)
   {
      int PosBNl = i / NbrBlocNc;
      int PosBNc = i % NbrBlocNc;

      if (PosBNl <  NbrBlocNl-1)  TabBlockNl[i] = BlockSize;
      else  TabBlockNl[i] = L_Nl - BlockSize * PosBNl;
      
      if (PosBNc <  NbrBlocNc-1) TabBlockNc[i] = BlockSize;
      else TabBlockNc[i] =  L_Nc - BlockSize * PosBNc;
      
      TabBlockDepNl[i] = BlockSize*PosBNl;  
      TabBlockDepNc[i] = BlockSize*PosBNc;
      // cerr << "Block " << i+1 << ": " << TabBlockNl[i] << "x" << TabBlockNc[i] <<
      //        " DepNl = " << TabBlockDepNl[i] << " DepNc = " << TabBlockDepNc[i] << endl;
   }

   if ((Comp_Method != COMP_MORPHO) && (Comp_Method != COMP_O_MORPHO))
   {
   int Nm = MIN(TabBlockNl[0], TabBlockNc[0]);
   int IndMin = 0;
   for (int i=1; i < NbrBlock; i++)
   {
      int Nmin = MIN(TabBlockNl[i], TabBlockNc[i]);
      if (Nmin < Nm) 
      {
         IndMin = i;
	 Nm = Nmin;
      }
   }
   int ScaleMax=(int) (log((float)Nm/(float)MAX_SIZE_LAST_SCALE) / log(2.)+ 1.);
   if (Nbr_Plan > ScaleMax)
   {
      cerr << endl << endl;
      if (ScaleMax < 2)
      {
         cerr << "Error: the number of scales is too high " << endl;
         cerr << "       Block number " << IndMin  << " size = " << TabBlockNl[IndMin] << "x" << TabBlockNc[IndMin] << endl;
         cerr << "       the maximum allowed number of scales is " << ScaleMax << endl;
         cerr << "       the block size must be modified "  << endl;
	 exit(-1);
      }
      else 
      {
         cerr << "Warning: the number of scales is too high " << endl;
         cerr << "       Block number " << IndMin  << " size = " << TabBlockNl[IndMin] << "x" << TabBlockNc[IndMin] << endl;
         cerr << "       the maximum allowed number of scales is " << ScaleMax << endl;
         cerr << "       The number of scales is modified: new values is "  << ScaleMax << endl;
         Nbr_Plan = ScaleMax; 
      }
   }
   }
}

/***************************************************/

void MR_CompData::get_block(int NumBlock)
{
   // int i,j;
   int Depi, Depj;
    if (UseBlock == False)
    {
        if (IntTrans == False)  BDat = &Dat;
        else
        {
            BImagInt = &ImagInt;
            BIDat = BImagInt->buffer();
        }
    }
    else
    {
       Nl = TabBlockNl[NumBlock];
       Nc = TabBlockNc[NumBlock];
       Depi = TabBlockDepNl[NumBlock];
       Depj = TabBlockDepNc[NumBlock];
       if (IntTrans == False)
       {
         BDat->resize(Nl,Nc);
	 io_read_block_ima((char *) File_Name_Imag, *BDat, Depi, Depj, InfoData);
         // for (i=0; i < Nl; i++)
         // for (j=0; j < Nc; j++) (*BDat)(i,j) = Dat(i+Depi, j+Depj);
       }
       else
       {
         BImagInt->resize(Nl,Nc);
         BIDat = BImagInt->buffer();
	 io_read_block_ima((char *) File_Name_Imag, *BImagInt, Depi, Depj, InfoData, NoBscale);
         //for (i=0; i < Nl; i++)
         //for (j=0; j < Nc; j++)  (*BImagInt)(i,j) =  ImagInt(i+Depi, j+Depj);
       }
    }
       
}

/***************************************************/

void MR_CompData::compress()
{
    int NumBlock;
    long Bef0=0;
    float UserNoiseIma = Noise_Ima;
    
    // controle if the method works with int or float
    IntTrans = work_with_int(Comp_Method);
     
    // controle for methods eliminating the noise
    NoiseThresholding = noise_threshold(Comp_Method);
     
    noise_compute (MAX_SCALE,Transform, DEFAULT_NL, 
                   DEFAULT_NC, MedianWindowSize);

    read_data();
    if (UseBlock == True)
    {
       if (BlockSize >= L_Nl)
       {
          cerr << "Warning: block option has no effect ... " << endl;
	  cerr << "         BlockSize is larger than image size," << endl;
       }
    }
    if (UseBlock == True) set_tabblock();
    // if ((UseBlock == True) && (NbrBlock == 1)) UseBlock = False;
    
    int TotalInputNbByte = L_Nl*L_Nc;
    switch (which_data_type())
    {
       case T_SHORT:    TotalInputNbByte *=2;break;
       case T_INT:      TotalInputNbByte *=4;break;
       case T_FLOAT:    TotalInputNbByte *=4;break;
       case T_DOUBLE:   TotalInputNbByte *=8;break;
       default: break;
    }
    
    // if  UseBlock = true, data are not in memory   
    if (UseBlock == False) replace_badpix();
    

    if (Verbose == True) print_info();

    OutFile = fopen (File_Name_Transform, "wb");
    store_header();
    if (UseBlock == True)   
    {
        if (IntTrans == False) 
                  BDat = new Ifloat(BlockSize, BlockSize, "block");
        else BImagInt = new Iint(BlockSize, BlockSize, "block");
        writeint(OutFile, L_Nl);
        writeint(OutFile, L_Nc);
     }
       
    for (NumBlock = 0; NumBlock < NbrBlock; NumBlock++)
    {
       get_block(NumBlock);
       Noise_Ima = UserNoiseIma;
       if ((Stat_Noise == NOISE_POISSON) && 
                 (Comp_Method != COMP_O_MIN_MAX) &&
		 (Comp_Method != COMP_INT_LIFTING)) variance_stab();

       if ((Verbose == True) && (UseBlock == True))
          cerr <<  "    Encoding block number " 
                 << NumBlock+1 << " Size " << Nl << "x" << Nc << endl;
       
                  
       if (UseBlock == True)  writeint(OutFile, 0); // write size
       long Bef=ftell(OutFile);
       if (NumBlock == 0) Bef0 = Bef;
       
       if (IntTrans == False) Residual.resize(Nl, Nc);
       else
       {
          if (IResidual != NULL) i_free(IResidual);
          IResidual = i_alloc (Nl*Nc);
       }
       encode_data();
     
     
       if (NoiseThresholding == True)
       {
          store_noise_info();
          // compress the residual
          // if KeepResi=1, the noise is compressed with a lost of information
          // if KeepResi=2, the noise is keeped without loosing anything
          if (KeepResi) encode_residu();
	  if (IntTrans == False) Residual.free();
        } 
        EncodedNbrByte = ftell(OutFile) - Bef;
        if (UseBlock == True)
        {	   
	   long Val = (-1) * (long) EncodedNbrByte - 4;
           int nf = fseek( OutFile,  Val, SEEK_CUR);
	   if (nf < 0)
	   {
	      cerr << "Error in fseek: request shift is " << Val  << endl;
	      exit(-1);
	   }           
	   writeint(OutFile,  EncodedNbrByte);
           nf = fseek(OutFile, EncodedNbrByte,SEEK_CUR);
	   if (nf < 0)
	   {
	      cerr << "Error in fseek: request shift is " << EncodedNbrByte  << endl;
              exit(-1);
	   }           
	   if (Verbose == True) cerr <<  "    Size block = " << EncodedNbrByte << endl;
        }
     }   
     
     
     if (Verbose == True)
     {
          if (Stat_Noise == NOISE_GAUSSIAN)
                            cerr << "Estimated Noise = " << Noise_Ima << endl;
          if (KeepResi) 
          cerr <<  "    RESIDUAL coding" << endl  << "    NoiseQuantif = " << NoiseQuantif<< endl;
              
          EncodedNbrByte = ftell(OutFile) - Bef0;
          float RealCompRatio = (float) TotalInputNbByte / (float) EncodedNbrByte;
          if (UseBudget == True)
          {
            cerr << "Target Compression Ratio = " <<  CompressionRatio <<
                "  Compression Ratio = " << RealCompRatio  << endl;
          }
          else cerr << "Compression Ratio (without header) = " << RealCompRatio  << endl;
     }      
     fclose (OutFile);
     if (IResidual != NULL) i_free(IResidual);
     if (IntTrans == False)  Residual.free();
}


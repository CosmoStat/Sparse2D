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
**    Date:  14/03/00 
**    
**    File:  MC_Comp.h
**
*******************************************************************************
**
**    DESCRIPTION   multichannel image compression
**    ----------- 
**                 
**    PARAMETRES    
**    ----------    
** 
**    RESULTS      
**    -------  
**
**
******************************************************************************/


#include "MC_Comp.h"
#include "IM3D_IO.h"

extern float PasCodeur;  /* CCD gain */
extern float SigmaGauss; /* CCD read-out noise standard deviation */
extern float MeanGauss;  /* CCD read-out noise mean */
 
static int InputNbByte;

static IO3DInfoData InfoData; 

#define NO_COLOR 0
#define COLOR_Y  1
#define COLOR_Cr 2
#define COLOR_Cb 3

/***************************************************/

/***************************************************/

void wt_compress(Ifloat &Imag, int Nbr_Plan, FILE *FileDes,  
                      float &Noise, float NSigma, float SignalQuantif,
                      Bool Verbose, Bool RecImag, int TColor)
{
    int Nl = Imag.nl();
    int Nc = Imag.nc();
    int i,j,s;
    int Nl_s=Nl,Nc_s=Nc;
    int NbrBand = (Nbr_Plan-1)*3 + 1;
    float *TabSigma = new float [NbrBand];
    float *TabParamQuant = new float [NbrBand];
    float *TabRate = new float [NbrBand];
    float *TabQuantif = new float [NbrBand];
    float *TabMean = new float [NbrBand];
    int *TabBandNl = new int [NbrBand];
    int *TabBandNc = new int [NbrBand];
    int *TabResolNl = new int [NbrBand];
    int *TabResolNc = new int [NbrBand];
    // FilterAnaSynt SelectFilter(F_MALLAT_7_9);
    // SubBandFilter SB1D(SelectFilter, NORM_L2);
    // SB1D.Border = I_MIRROR;
    // Ortho_2D_WT WT2D(SB1D);
    type_lift LiftingTrans=TL_F79;
    Lifting Clift1D(LiftingTrans);
    Clift1D.Border = I_MIRROR;
    Ortho_2D_WT WT2D(Clift1D);

    writeint(FileDes, Nl);
    writeint(FileDes, Nc);
    writeint(FileDes, Nbr_Plan); 
    // if (Noise < FLOAT_EPSILON) Noise = 0.5;
    // cout << "WT2D " << endl;
    WT2D.transform(Imag, Nbr_Plan);
    // cout << "WT2D OK" << endl;
    // Band image size computation
    if (Noise < FLOAT_EPSILON)
    {
        int Nlb = (Nl_s+1)/2;
        int Ncb = Nc_s/2;
        int Depj = (Nc_s+1)/2;
        int It, Nit=3;
        double S0,S1,S2,Sm=0,x;
        double Sigma=0.;

        for (It = 0; It < Nit; It++)
        {
            S0 = S1 = S2 = 0.;
            for (i = 0; i < Nlb; i++)
            for (j = 0; j < Ncb; j++)
            {
                x = Imag(i,j+Depj);
                if (ABS(x) > FLOAT_EPSILON)
                {
	           if ((It == 0) || (ABS(x) < Sm))
	           { 
	              S0 ++;
	              S1 += x;
	              S2 += x*x;
	           }
                }
            }
            if (S0 == 0) S0=1;
            Sigma = sqrt(S2/S0);
            Sm = 3. * Sigma;       
        }
        Noise = (float) Sigma;
        if (Verbose == True) cout << "Estimated noise = " << Noise << endl;
    }
    // else cout << "  noise = " << Noise << endl;

    for (s=0; s < NbrBand; s++)
    {
        TabParamQuant[s] = SignalQuantif;
        TabSigma[s]= Noise; 
        TabQuantif[s] = TabParamQuant[s]*TabSigma[s];
	if (TColor == 2) TabQuantif[s] *= 2;
	else if (TColor == 3) TabQuantif[s] *= 4;
        TabMean[s] = 0.;
        // cout << "Band s+1 : Noise = " << TabSigma[s] << endl;
    }
  
    for (s=0; s < NbrBand-1; s+=3)
    {
        TabBandNl[s] = (Nl_s+1)/2;
        TabBandNc[s] = Nc_s/2;
        TabBandNl[s+1] = Nl_s/2;
        TabBandNc[s+1] = (Nc_s+1)/2;
        TabBandNl[s+2] = Nl_s/2;
        TabBandNc[s+2] = Nc_s/2;
        Nl_s = (Nl_s+1)/2;
        Nc_s = (Nc_s+1)/2;
        TabResolNl[s] = Nl_s;
        TabResolNc[s] = Nc_s;
        TabResolNl[s+1] = Nl_s;
        TabResolNc[s+1] = Nc_s;
        TabResolNl[s+2] = Nl_s;
        TabResolNc[s+2] = Nc_s;
   }     
    s = NbrBand-1;
    TabBandNl[s] = Nl_s;
    TabBandNc[s] = Nc_s;
    TabResolNl[s] = Nl_s;
    TabResolNc[s] = Nc_s;

    // extern long size_enc;
    int NbrCoef;
    float Quantif;
    int *qdata;
    for (s = NbrBand-1; s >= 0; s--)
    {      
       int Ind=0, Depi=0, Depj=0;
       details Detail;
       int Scale;
       Nl_s = TabBandNl[s];
       Nc_s =  TabBandNc[s];
       NbrCoef =  Nl_s*Nc_s;
       qdata = i_alloc(NbrCoef);
       band2scale(s,TO_MALLAT, NbrBand, Scale, Detail);
       switch (Detail)
       {
        case D_HORIZONTAL:
               Depj += TabResolNc[s];
               break;
        case D_DIAGONAL:  
               Depi += TabResolNl[s];
               Depj += TabResolNc[s];
	       if (TColor > 0)
	       {
	          if (Scale == 0) TabQuantif[s] *= 2;
	          else TabQuantif[s] *= 1.5;
	       }
             break;
        case D_VERTICAL:  
               Depi += TabResolNl[s];
               break;
        case I_SMOOTH: break;
        case D_NULL:
        default: cerr << "Error: unknown detail" << endl;
                 exit(0);
                 break;
       }         
       if (Scale == 0) TabQuantif[s] *= 1.5;
       Quantif =  TabQuantif[s];

       if (Quantif > 0) 
       {
         float LevelDetect = TabSigma[s]*NSigma;

     // cerr << "Set " << s+1 << "  LevelDetect  = " <<  LevelDetect 
     //     << " Depi = " << Depi << "   Depj = " << Depj << endl;


         for (i = 0; i < Nl_s; i++)
         for (j = 0; j < Nc_s; j++)
         {
             float Val = Imag(i+Depi,j+Depj);
             if (ABS(Val) < LevelDetect)  qdata[Ind++] = 0;
             else qdata[Ind++] = quant(Val, Quantif);
         }
         encode(FileDes,  qdata,  Nl_s,  Nc_s, Quantif);

         if (RecImag == True)
         {
            Ind=0; 
            for (i = 0; i < Nl_s; i++)
            for (j = 0; j < Nc_s; j++) 
                  Imag(i+Depi,j+Depj) = (float) qdata[Ind++] * Quantif;
         }
       }
       else
       { 
          if (s !=  NbrBand-1) qdata[0] = 0;
          else qdata[0] =  quant(TabMean[s], 1.);
          encode( FileDes,  qdata, 1,1, (float) 0.);
          if (RecImag == True)
          {
             for (i = 0; i < Nl_s; i++)
           for (j = 0; j < Nc_s; j++) 
                  Imag(i+Depi,j+Depj) = (float) qdata[0];
          }
       }
       if (Verbose == True)
       {
           //extern long size_enc;
           // cerr << "Set " << s+1 << "  NbrCoef = " <<  NbrCoef 
           //     << " Quant = " << Quantif <<  " Size = " << size_enc << endl;
         }
        i_free(qdata);
     }
   if (RecImag == True) WT2D.recons(Imag, Nbr_Plan);
    
   delete [] TabSigma;
   delete [] TabParamQuant;
   delete [] TabRate;
   delete [] TabMean;
   delete [] TabQuantif;
   delete [] TabBandNl;
   delete [] TabBandNc;
   delete [] TabResolNl;
   delete [] TabResolNc;
}

/*************************************/

void MR_Comp3DData::encode_residu()
{
    float Sigmaf=0;
    int *Resi;
    int i,j;
    
    if (Verbose == True) cerr <<  "Residual coding"<< endl;

    if (RGBImage == True) yuv_to_rgb(*BDat);
    Resi = i_alloc(Nl*Nc);
    for (int f=0; f < Nf; f++) 
    {
        Sigmaf = TabSigmaNoise[f] * NoiseQuantif;
        if (KeepResi == 2) Sigmaf = 1.;
        for (i = 0; i < Nl; i++) 
	for (j = 0; j < Nc; j++)
	{
            Resi[i*Nc+j] = quant ((*BDat)(j,i,f) - Residual(j,i,f), Sigmaf); 
	    Residual(j,i,f) = Resi[i*Nc+j] * Sigmaf;
	}    
        encode (OutFile, Resi, Nl, Nc, Sigmaf);
    }
    io_3d_write_data("resi.fits", Residual);
    i_free(Resi);
}
 

/***************************************************/

void MR_Comp3DData::store_noise_info()
{
    writeint (OutFile, KeepResi);
    for (int f=0; f < Nf; f++)
    {
       writefloat(OutFile, TabSigmaNoise[f]);
       // cout << "Band " << f+1 << " Noise = " << TabSigmaNoise[f] << endl;
    }
}

/***************************************************/


void MR_Comp3DData::encode_data()
{
    int i,j;
    Bool SupIsol=False;
    int MaxIter=1;
    Bool Rec = (KeepResi > 0) ? True : False;
    int TColor=0;
    
    InputNbByte = Nl*Nc*Nf;
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
    
    // Buffer alloction in which we insert each frame
    // if (KeepResi > 0) we need the original data to calculate the 
    // residual, otherwise we can use the same buffer for the WT calculation
    Ifloat Buffer;
    if (KeepResi > 0) Buffer.alloc(Nl, Nc, "buffer");

    // if (TargetNbrByte > 0) TargetNbrByte /= Nf;
    for (int f=0; f < Nf; f++)
    {
       int BudgetFrame=0;
       if (TargetNbrByte > 0) 
       {
          if ((RGBImage == False) || (RGBChromDownSample == False)) 
	                               BudgetFrame = TargetNbrByte / Nf;
          else if (f == 0) BudgetFrame = TargetNbrByte / 2;
	  else BudgetFrame = TargetNbrByte / 4;
       }
       
       TabSigmaNoise[f] = Noise_Ima;
       if (RGBImage == False) TColor=0;
       else TColor=f+1;
       
       
       if ((f == 0) || (RGBImage == False) || (RGBChromDownSample == False))
       {
          if (KeepResi > 0) 
          {
             for (i=0; i < Nl; i++)
             for (j=0; j < Nc; j++) Buffer(i,j) = (*BDat)(j,i,f);
          }
          else
          {
             float *PtrData = (*BDat).buffer() + f*Nl*Nc;
             Buffer.alloc(PtrData, Nl, Nc, "buffer");
          }
       }
       else // Color images using downsampled chrominance maps
       {
          if (KeepResi > 0) Buffer.resize((Nl+1)/2, (Nc+1)/2);
          else  if (f == 1) Buffer.alloc((Nl+1)/2, (Nc+1)/2, "buffer");
          for (i=0; i < Buffer.nl(); i++)
          for (j=0; j < Buffer.nc(); j++) Buffer(i,j) = (*BDat)(2*j,2*i,f);
          TabSigmaNoise[f] = TabSigmaNoise[f-1];
       }
       // INFO(Buffer, "buffer");
       switch (Comp_Method)
       {
        case COMP_MRMEDIAN:
             mr_compress (Buffer, Nbr_Plan, TabSigmaNoise[f], N_Sigma, 
                     SignalQuantif, (int)BudgetFrame, OutFile, Buffer,  
                     SupIsol, NoiseInData, MedianWindowSize, SupNeg, Verbose);
          break;
  
         case COMP_INT_LIFTING: 
//             {
//              Lifting Clift1D(LiftingTrans);
//              Ortho_2D_WT WT2D(Clift1D);
//              WT2D.transform(*BImagInt, Nbr_Plan);
// 	     mr_i_ortho_compress(*BImagInt, Nbr_Plan, OutFile);
//             } 
 	     break;
         case COMP_MALLAT:
          // cout << " COMP_MALLAT " <<  TabSigmaNoise[f] << endl;

            if (TargetNbrByte > 0)
              mr_budcompress (Buffer, Transform, Nbr_Plan, (int) BudgetFrame,
                            OutFile, Buffer, NoiseInData, TabSigmaNoise[f], 
                            N_Sigma,  SignalQuantif, SupIsol, True, 
                            SupNeg, MedianWindowSize, MaxIter, Verbose, Rec); 
            else wt_compress(Buffer, Nbr_Plan, OutFile, TabSigmaNoise[f], N_Sigma, 
                            SignalQuantif, Verbose, Rec, TColor);
  	     break;
         case COMP_HAAR:
         case COMP_MIN_MAX:
         case COMP_FEAUVEAU:
         case COMP_WT_PMT:
          cout << "mr_budcompress " << endl;
             mr_budcompress (Buffer, Transform, Nbr_Plan, (int) TargetNbrByte,
                            OutFile, Buffer, NoiseInData, TabSigmaNoise[f], 
                            N_Sigma,  SignalQuantif, SupIsol, True, 
                            SupNeg, MedianWindowSize, MaxIter, Verbose, Rec);
             break;
         default:
             cerr << "Error: unknown compression method ... " << endl;
             exit(-1);
             break;
      } 
       if (KeepResi > 0) 
       {       
          if ((f == 0) || (RGBImage == False) || (RGBChromDownSample == False))
          {
             for (i=0; i < Nl; i++)
             for (j=0; j < Nc; j++) Residual(j,i,f) = Buffer(i,j);
          }
          else
          {
             for (i=0; i < Nl; i++)
             for (j=0; j < Nc; j++) Residual(j,i,f) = Buffer(i/2,j/2);
          }
       }    
   }
}

/***************************************************/

MR_Comp3DData::~MR_Comp3DData()
{
     if (TabBlockNl != NULL) delete  [] TabBlockNl;
     if (TabBlockNc != NULL) delete [] TabBlockNc;
     if (TabBlockDepNl  != NULL) delete [] TabBlockDepNl;
     if (TabBlockDepNc != NULL) delete [] TabBlockDepNc;
}

/***************************************************/
 
void MR_Comp3DData::read_data()
{
    initfield(&Header);
    if (UseBlock == False) 
    {
       extern type_3d_format IO_3D_Format;
       io_3d_read_data(File_Name_Imag, Dat);
       if (RGBImage == True) 
       {
          if (Dat.nz() != 3)
          {
             cout << "Error: input image is not 3-channel RBG image ... " << endl;
             exit(-1);
          }
          rgb_to_yuv(Dat);
       }
       Nl = Dat.ny();
       Nc = Dat.nx();
       Nf = NbrImag = Dat.nz();
       L_Nl = Nl;
       L_Nc = Nc;
       FormatInput = IO_3D_Format;
    }
    else
    {
       InfoData.init_reading(File_Name_Imag, ReadInPseudo);
       Nl = L_Nl = InfoData.Nl;
       Nc = L_Nc = InfoData.Nc;
       Nf = NbrImag = 3;
       FormatInput = InfoData.Format;
    }
    TabSigmaNoise = new float [Nf];
}

/***************************************************/

void MR_Comp3DData::print_info()
{
   extern softinfo Soft;
   cerr << endl;
   cerr <<  Soft.banner() << endl;
   cerr <<  StringComp(Comp_Method) << endl << endl;

   cerr <<  "Input File = "<< File_Name_Imag << endl;
   cerr <<   "Image size (" << Nc << "," << Nl << "," << Nf << ") N_Sigma " << N_Sigma << endl;
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

void MR_Comp3DData::store_header()
{
    extern softinfo Soft;
    fwrite(ValIdent, sizeof(short), 3, OutFile);

    writefloat (OutFile, (float)  Soft.release());
    writeint (OutFile, (int) Comp_Method);
    if (Comp_Method == COMP_INT_LIFTING) 
         writeint (OutFile, (int) LiftingTrans);
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
    if (NbrImag < 2)
    {
        cout << "Error: input data are not a multichannel data set ... " << endl;
        exit(-1);
    }
    if (RGBImage == True) 
    {
         writeint (OutFile, (int) 1);
         if (NbrImag != 3)
         {
            cout << "Error: input data are not a RGB image ... " << endl;
            exit(-1);
         }
         if (RGBChromDownSample == True) writeint (OutFile, (int) 1);
         else  writeint (OutFile, (int) 0);
    }
    else 
    {
       writeint (OutFile, (int) 0);
       cout << "Error: not implemented ... " << endl;
       exit(-1);
    }
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
    
    // BadPixel = false 
    writeint (OutFile, (int) 0);
 
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
}


/***************************************************/

void MR_Comp3DData::get_block(int NumBlock)
{
   // int i,j;
   int Depi, Depj;
    if (UseBlock == False)
    {
        BDat = &Dat;
    }
    else
    {
       Nl = TabBlockNl[NumBlock];
       Nc = TabBlockNc[NumBlock];
       Depi = TabBlockDepNl[NumBlock];
       Depj = TabBlockDepNc[NumBlock];
       BDat->reform(Nc,Nl,Nf);
       io_3d_read_block_ima(File_Name_Imag, *BDat, Depi, Depj, InfoData); 
       if (RGBImage == True) rgb_to_yuv(*BDat);
    }
}


/***************************************************/


void MR_Comp3DData::compress()
{
    int NumBlock;
    long Bef0=0;
    float UserNoiseIma = Noise_Ima;
              
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
    int TotalInputNbByte = L_Nl*L_Nc*Nf;
    switch (which_data_type())
    {
       case T_SHORT:    TotalInputNbByte *=2;break;
       case T_INT:      TotalInputNbByte *=4;break;
       case T_FLOAT:    TotalInputNbByte *=4;break;
       case T_DOUBLE:   TotalInputNbByte *=8;break;
       default: break;
    }
    
    if (Verbose == True) print_info();

    OutFile = fopen (File_Name_Transform, "wb");
    store_header();
    if (UseBlock == True)   
    {
        BDat = new fltarray(BlockSize, BlockSize, Nf);
        writeint(OutFile, L_Nl);
        writeint(OutFile, L_Nc);
	writeint(OutFile, Nf);
     }
       
    for (NumBlock = 0; NumBlock < NbrBlock; NumBlock++)
    {
       get_block(NumBlock);
       Noise_Ima = UserNoiseIma;
 
       if ((Verbose == True) && (UseBlock == True))
          cerr <<  "    Encoding block number " 
                 << NumBlock+1 << " Size " << Nl << "x" << Nc << endl;
       
       if (NbrBlock > 1)  writeint(OutFile, 0); // write size

       long Bef=ftell(OutFile);
       if (NumBlock == 0) Bef0 = Bef;
       
       if (KeepResi) Residual.reform(Nc, Nl, Nf);
       encode_data();
     
       if ((KeepResi) && (NoiseThresholding == True))
       {
          store_noise_info();
          // compress the residual
          // if KeepResi=1, the noise is compressed with a lost of information
          // if KeepResi=2, the noise is keeped without loosing anything
          if (KeepResi) encode_residu();
	  if (IntTrans == False) Residual.free();
        } 
        EncodedNbrByte = ftell(OutFile) - Bef;
        if (NbrBlock > 1)
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
     Residual.free();
}

/***************************************************/
 
void MR_Comp3DData::set_tabblock()
{
   NbrBlocNl = L_Nl/BlockSize;
   if (NbrBlocNl * BlockSize != L_Nl)  NbrBlocNl++;
   NbrBlocNc  = L_Nc/BlockSize;
   if (NbrBlocNc * BlockSize != L_Nc)  NbrBlocNc++;
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
      //cout << "Block " << i+1 << ": " << TabBlockNl[i] << "x" << TabBlockNc[i] <<
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
   int ScaleMax=iround(log((float)Nm/(float)MAX_SIZE_LAST_SCALE) / log(2.)+ 1.);
   if (Nbr_Plan > ScaleMax)
   {
      cerr << endl << endl;
      cerr << "Error: the number of scales is too high " << endl;
      cerr << "       Block number " << IndMin  << " size = " << TabBlockNl[IndMin] << "x" << TabBlockNc[IndMin] << endl;
      cerr << "       the maximum allowed number of scales is " << ScaleMax << endl;
      if (ScaleMax < 2)
             cerr << "       the block size must be modified "  << endl;
      else
             cerr << "       the block size or the number of scales must be modified "  << endl;
      exit(-10);
   }}
}

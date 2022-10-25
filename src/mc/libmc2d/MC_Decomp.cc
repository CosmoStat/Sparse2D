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
**    DESCRIPTION  multichannel image decompression
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
// extern float BadPixalVal; 
// extern Bool BadPixel;

extern Bool InputFromStdin;// = False;

static IO3DInfoData InfoData;

/***************************************************/
// Decompression part
/***************************************************/

/*********************************************************************/

// void param_visual_quant(int Orient, int NumCol, int Level, double Dist)
void param_visual_quant()
{
   int Orient;
   int NumCol;
   int Level;
   double Dist=32.;
   double Ampli[4][6] =
            { {0.62171, 0.34537, 0.18004, 0.091401, 0.045943, 0.023013},
              {0.67234, 0.41317, 0.22727, 0.11792, 0.059758, 0.030018},
              {0.72709, 0.49428, 0.28688, 0.15214, 0.077727, 0.039156},
              {0.67234, 0.41317, 0.22727, 0.11792, 0.059758, 0.030018}};

   double  A[3] = {0.495, 0.9, 2.};
   double  K[3] = {0.466, 0.53, 0.4};
   double F[3] = {0.401, 0.38, 0.269};
   double G[4][3] = {{1.5, 1.9, 1.6},
                     {1., 1., 1.},
                     {0.535, 0.59, 0.59},
                     {1., 1., 1.}};

   for ( NumCol = 0;  NumCol < 3;  NumCol++)
   {
   for (Orient = 0; Orient < 4; Orient++)
   {
   for (Level = 1; Level <= 4; Level++)
   {
   double Val;
   double Num;
   double LogVal= log(pow(2., (double) Level) * F[NumCol]*G[Orient][NumCol]/ Dist);
   Num = K[NumCol]*LogVal*LogVal;
   Val = 2. / Ampli[Orient][NumCol] * A[NumCol]  * exp(Num);
   printf("%5.2f ", (float ) sqrt(Val));
   }
   printf("\n");
   }
   printf("\n");
   }



}

/*********************************************************************/

void wt_decompress (Ifloat &Imag, type_transform Transform, 
                    int Resol, FILE *FileDes, int &Nli, int &Nci, Bool Verbose)
{
    int i,j,s;
    int Nl,Nc,Nl_s,Nc_s;
    int Nbr_Plan;
 
    /* read the original image size */
    Nl =readint(FileDes); /* x size of image */
    Nc =readint(FileDes); /* y size of image */
    Nli = Nl;
    Nci = Nc;
    Imag.resize (Nl, Nc);

    /* read the number of scales */
    Nbr_Plan=readint(FileDes); 
     
    if (Verbose == True)
    {
       fprintf (stderr,"DECODING wavelet transform: nbr_scale = %d\n", Nbr_Plan);
       fprintf (stderr,"Original image: Nl = %d, Nc = %d\n", Nl, Nc);
     }

    /* Resol must be < Nbr_Plan and > 0 */
    if ((Resol < 0) || (Resol >= Nbr_Plan))
    {
      fprintf (stderr, "Error: bad parameter Resol: Resol = %d\n", Resol);
      fprintf (stderr, "       Resol must verify 0 <= Resol < %d\n", Nbr_Plan);
      exit (0);
    }

    /* if we don't want the full resolution */
    if (Resol > 0)
    {
       for (s=0; s < Resol; s++)
       {
          Nl = (Nl+1)/2;
          Nc = (Nc+1)/2;
       }
       Nbr_Plan -= Resol;
    }     
    if (Verbose == True) printf ("\n Create an image of size (%d, %d)\n", Nl, Nc);
    Imag.resize (Nl, Nc);
    int NbrBand = (Nbr_Plan-1)*3 + 1;

//     FilterAnaSynt SelectFilter(F_MALLAT_7_9);
//     SubBandFilter SB1D(SelectFilter, NORM_L2);
//     SB1D.Border = I_MIRROR;
//     Ortho_2D_WT WT2D(SB1D);
    type_lift LiftingTrans=TL_F79;
    Lifting Clift1D(LiftingTrans);
    Clift1D.Border = I_MIRROR;
    Ortho_2D_WT WT2D(Clift1D);
    int *TabBandNl = new int [NbrBand];
    int *TabBandNc = new int [NbrBand];
    int *TabResolNl = new int [NbrBand];
    int *TabResolNc = new int [NbrBand];
    Nl_s = Nl;
    Nc_s = Nc;
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

    //extern long size_enc;
    float *LevelQ = new float[NbrBand];
    for (s =  NbrBand -1 ; s  >= 0; s--) 
    {
       int NbrCoef;
       int Ind=0, Depi=0, Depj=0;
       Nl_s = TabBandNl[s];
       Nc_s =  TabBandNc[s];
       NbrCoef = Nl_s*Nc_s;
       int *data;
       float Level;
       int Nls,Ncs;
       decode (FileDes, &data, &Nls, &Ncs, &Level);
       LevelQ[s] = Level;
       //int SizeScale = size_enc;
       // if (Verbose == True)
       //      cerr << "Set " << s+1 << " " << Nls << "x" << Ncs << 
       //          " Quantif = " << Level  <<  " Size = " << SizeScale << endl;
       if (s != NbrBand -1)
       {
          details Detail;
          int Scale;
          band2scale(s,TO_MALLAT, NbrBand, Scale, Detail);
          switch (Detail)
          {
             case D_HORIZONTAL:
               Depj += TabResolNc[s];
               break;
             case D_DIAGONAL:  
               Depi += TabResolNl[s];
               Depj += TabResolNc[s];
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
       }
       if (Level > FLOAT_EPSILON)
       {
           if (NbrBand > 1)
           {
              for (i = 0; i < Nl_s; i++)
              for (j = 0; j < Nc_s; j++)
                  Imag(i+Depi,j+Depj) = ((float) data[Ind++]) * Level;
           }
           else
           {
              for (i = 0; i < Nl_s; i++)
              for (j = 0; j < Nc_s; j++) Imag(i,j) = ((float) data[Ind++]) * Level;
           }
       }
       else
       {
           if (NbrBand > 1)
           {
              for (i = 0; i < Nl_s; i++)
              for (j = 0; j < Nc_s; j++) Imag(i+Depi,j+Depj) =  data[0];
           }
           else  
           {
              for (i = 0; i < Nl_s; i++)
              for (j = 0; j < Nc_s; j++)  Imag(i,j) = data[0];
           }
       }
       i_free(data);
    }        

   if (NbrBand > 1) WT2D.recons(Imag, Nbr_Plan);
 
   // In case of sub resolution reconstruction, the reconstructued
   // image must be renormalized because we used a L2 normalization
   float NormCoef=1.;
   if (Resol > 0)
   {
      for (s=0; s < Resol; s++) NormCoef *=2;
      for (i = 0; i < Nl; i++)
      for (j = 0; j < Nc; j++) Imag(i,j) /= NormCoef;
  
   }     

   delete [] LevelQ;
   delete [] TabBandNl;
   delete [] TabBandNc;
   delete [] TabResolNl;
   delete [] TabResolNc;
}
    

/****************************************************************************/

 
void MR_Comp3DData::read_tinfo()
{
     extern softinfo Soft;
   
    short ValIdent[3];
    fread(ValIdent, sizeof(short), 3, InFile);
//    cout << "V1 = " << ValIdent[0] << endl;
//    cout << "V2 = " << ValIdent[1] << endl;
//    cout << "V3 = " << ValIdent[2] << endl;
    
    float SoftVersion;
    SoftVersion = readfloat (InFile);
    if ((SoftVersion > 1) && (Soft.release() == 1))
    {
       cerr << "Error: MR/1 version 1.0 cannot decompress files " << endl;
       cerr << "       which have been compressed with more recent release of MR/1 ... " << endl;
    }
    Comp_Method  = (type_comp) readint (InFile);
    switch (Comp_Method)
    {       
       case COMP_HAAR: Transform = TO_HAAR; break;
       case COMP_MALLAT: Transform = TO_MALLAT;break;
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
    // if (SoftVersion > 1.)
    {
       L_Nl = readint (InFile);
       L_Nc = readint (InFile);
       KeepResi = readint (InFile);
    }
    if (Resol >= Nbr_Plan) Resol = Nbr_Plan-1;
    
    IntTrans = work_with_int(Comp_Method);     
    NoiseThresholding = noise_threshold(Comp_Method);
        
    FormatInput = (type_3d_format) readint (InFile); 
    // test for compatibility software version     
    // if (SoftVersion > 1.)
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
    Nf = NbrImag = readint (InFile);
    if (NbrImag < 2)
    {
        cout << "Error: input data are not a multichannel data set ... " << endl;
        exit(-1);
    }
    RGBImage = (readint (InFile) == 1) ?  True: False;
    if (RGBImage == True)
    {
       if (Verbose == True) cout << "RGB image " << endl;
       RGBChromDownSample = (readint (InFile) == 1) ?  True: False;
       if ((Verbose == True) && (RGBChromDownSample == True))
                            cout << "Chrominance map downsampled. " << endl;
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
       cout << "Problem: BadPixel must be false " << endl;
       exit(0);
       // BadPixel = True;
       // BadPixalVal = readfloat (InFile);
       // MinDat = readfloat (InFile);
       // MaxDat = readfloat (InFile);
    }
    
    if (Verbose == True) 
    {
        cerr << endl;
        cerr << Soft.banner() << endl;
        cerr <<  StringComp(Comp_Method) << endl << endl;
        cerr << "Input File = " <<  File_Name_Transform << endl;
        cerr << "Nx = " << L_Nc << " Ny = " << L_Nl << " Nz = " << Nf << endl;
        cerr << "Output File = " << File_Name_Imag << endl;
        if (NbrBlock > 1) cerr <<  "Number of blocks = " << NbrBlock << endl;
        if (NbrBlock > 1) cerr <<  "block size = " << BlockSize << endl;
        if (NbrBlock > 1) cerr <<  "block overlapping = " << BlockOverlap<< endl;
        if (Resol >= 0) cerr << "Resolution = " << Resol << endl;
//         if (BadPixel == True)
//         {
//             cerr << "MinDat = " << MinDat << endl;
//             cerr << "MaxDat = " << MaxDat << endl;
//         }
    }
}

/*****************************************************************/
 
void MR_Comp3DData::header_resol()
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

void MR_Comp3DData::read_header()
{
  initfield(&Header);

  extern type_3d_format IO_3D_Format;
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
   type_3d_format FormatData = io_which_3d_format(File_Name_Imag);
   if (FormatData == F3D_UNKNOWN)
   {  
      if (FormatInput != F3D_UNKNOWN) IO_3D_Format = FormatInput;
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

void MR_Comp3DData::decode_data()
{
   int i,j,Nli,Nci;
   Ifloat Buffer;
   // int NScale;
   // cout << "DEC" << Comp_Method << endl;
   // param_visual_quant();

    for (int f=0; f < Nf; f++)
    {
    switch (Comp_Method)
    {
       case COMP_MALLAT:
         if (Resol == -1) 
              wt_decompress(Buffer,Transform,0,InFile,Nli,Nci,Verbose);
         else wt_decompress(Buffer,Transform,Resol,InFile,Nli,Nci,Verbose);
         Nl = Buffer.nl();
         Nc = Buffer.nc(); 
         break;
       case COMP_HAAR:  
       case COMP_FEAUVEAU:
       case COMP_MIN_MAX:
       case COMP_WT_PMT:
          // if (Comp_Method != COMP_HAAR) IterRec = 0;
          if (Resol == -1)   
	       mr_buddecompress(Buffer, Transform, 0,InFile,Nli,Nci,IterRec,Verbose);
          else mr_buddecompress(Buffer, Transform, Resol,InFile,Nli,Nci,IterRec,Verbose);
              Nl = Buffer.nl();
              Nc = Buffer.nc(); 

             // Haar scaling function is not normalized by 1.
             // It means that pixels values of the image at a lower 
             // resolution are higher than the original image
             // ===> we do a renormalization
             if ((Resol > 0) && (Comp_Method == COMP_HAAR))
             {
                float Coef = 1.;
                for (i=0; i < Resol; i++) Coef *= 2.;
                for (i=0; i < Buffer.nl(); i++)
                for (j=0; j < Buffer.nc(); j++) Buffer(i,j) /= Coef;
             }
          break;
       case COMP_MRMEDIAN:
             if (Resol == -1) mr_decompress (Buffer,0,InFile,Nli,Nci,Verbose);
             else mr_decompress (Buffer,Resol,InFile,Nli,Nci,Verbose);
             Nl =  Buffer.nl();
             Nc =  Buffer.nc();
	     break;
      case COMP_MORPHO:
             morpho_decompress (Buffer, N_Iter,InFile,
                                         Elstr_Size,Elstr_Shape,Verbose);
             Nl = Dat.ny();
             Nc = Dat.nx();
	     Nli = Nl;
	     Nci = Nc;
            break;
       case COMP_INT_LIFTING:
//            {
//	     // cout << "lift" << endl;
//             if (Resol == -1)  
// 	          mr_i_ortho_decompress (ImagInt,0,  NScale, InFile, Nli,Nci,Verbose);
//              else mr_i_ortho_decompress (ImagInt, Resol, NScale, InFile, Nli,Nci,Verbose);
//              Nl = ImagInt.nl();
//              Nc = ImagInt.nc();
// 	     Lifting Clift1D(LiftingTrans);
//              Ortho_2D_WT WT2D(Clift1D);
//              WT2D.recons(ImagInt, NScale);
// 	     BIDat = ImagInt.buffer();
//	    }
//             break;        
       case COMP_O_MORPHO:
       case COMP_O_MRMEDIAN:
             cout << "Error: not implemented method " << endl;
             exit(-1);
             break;
       default:
             cerr << "Error: unknown compression method ... " << endl;
             exit(-1);
             break;
    } 
     if (f == 0) (*BDat).reform(Nc,Nl,Nf);
     if ((f == 0) || (RGBImage == False) || (RGBChromDownSample == False))
     {
        for (i=0; i < Buffer.nl(); i++)
        for (j=0; j < Buffer.nc(); j++) (*BDat)(j,i,f) = Buffer(i,j);
     }
     else
     {
        Nl = Nli = (*BDat).ny();
        Nc = Nci = (*BDat).nx();
        for (i=0; i < Nl; i++)
        for (j=0; j < Nc; j++)  (*BDat)(j,i,f) = Buffer(i/2,j/2);
     }
     // cout << "skip to next block" << endl;
     if (Resol > 0)
     {  
	 int NB_seek = number_band_per_resol(Transform)*Resol;
	     // cout << " NB_seek = " << NB_seek << endl;
	 seek_band(NB_seek);
     }
   }

    // inverse transform    
//     if ((Stat_Noise == NOISE_POISSON) && 
//             (Comp_Method != COMP_O_MIN_MAX) &&
//             (Comp_Method != COMP_INT_LIFTING))
//     {
//         if (IntTrans == False)
//              noise_inverse_poisson_transform (*BDat, *BDat);
//         else noise_inverse_poisson_transform (BIDat, BIDat, Nl, Nc);
//     }
    
    // if it is not a block compression, the original image size is
    // equal to Nli,Nci
    if (NbrBlock <= 1)
    {
       L_Nl = Nli;
       L_Nc = Nci;
    }
}


/***************************************************************/

void MR_Comp3DData::read_noise_info()
{
    KeepResi = readint(InFile);
    TabSigmaNoise = new float [Nf];
    for (int f=0; f < Nf; f++)
    {
        TabSigmaNoise[f] = readfloat(InFile);
        // cout << "Band " << f+1 << " Noise = " << TabSigmaNoise[f] << endl;
    }
}

/*****************************************************************/

void MR_Comp3DData::decode_residu()
{
   float Sigmaf;
   int i,j,k;
   int *IResidual=NULL;  
   for (int f=0; f < Nf; f++) 
   {
     k=0;
     if ((Resol==-1) && (KeepResi!=0))
     {        
        decode(InFile, &IResidual, &Nl, &Nc, &Sigmaf);
        // if (Verbose == True)   cerr << "Sigma residu = " << Sigmaf << endl;
        // if ((f == 0) || (RGBImage == False) || (RGBChromDownSample == False))
	{
           for (i = 0; i < Nl; i++)
           for (j = 0; j < Nc; j++) (*BDat)(j,i,f) += IResidual[k++] * Sigmaf;
	}
// 	else
// 	{
//            for (i=0; i < (*BDat).ny(); i++)
//            for (j=0; j < (*BDat).nx(); j++)  
// 	     (*BDat)(j,i,f) += IResidual[i/2*Nc+j/2] * Sigmaf;
// 	}
     }
   }
   if (IResidual != NULL) delete IResidual;
}

/*****************************************************************/

void MR_Comp3DData::put_block(int NumBlock)
{
    int Depi, Depj;
    Depj = TabBlockDepNc[NumBlock];
    Depi = TabBlockDepNl[NumBlock];
    
    io_3d_write_block_ima(File_Name_Imag, *BDat, Depi, Depj, InfoData, NoBscale);          
}

/*****************************************************************/

void MR_Comp3DData::write()
{
   // write to ouput file the decompressed data
   // cout << "ReadHd = " << ReadHd << endl;
   
   if (ReadHd > 0)  io_3d_write_data(File_Name_Imag, Dat, &Header);
   else  io_3d_write_data(File_Name_Imag, Dat);
}

/*****************************************************************/

void MR_Comp3DData::set_block_info_data()
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
    
   InfoData.Format = FormatInput;
   if (ReadHd > 0) InfoData.PtrFits = &Header;
   else
   {
       InfoData.PtrFits = new fitsstruct;
       (InfoData.PtrFits)->hd_init(Nci, Nli, Nf);
   }
   InfoData.init_writing(File_Name_Imag, Nli, Nci, Nf);
   BDat = new fltarray;
}

/*****************************************************************/

void MR_Comp3DData::decompress()
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
    if (NbrBlock > 1) 
    {
       L_Nl = readint(InFile);
       L_Nc = readint(InFile); 
       Nf = readint(InFile);

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
          


        // cout << "decode noise" << endl;          
          if (RGBImage == True) yuv_to_rgb(*BDat);
          if ((KeepResi) && (NoiseThresholding == True))
	  {
	      // read info parameters
	      read_noise_info();
 	      if (Resol < 0)  
              {
                  decode_residu();
               }
              else if (KeepResi)
	      {
                 // cout << "fseek noise" << endl;
		 seek_band(1);
	      }
	  }
          // if (RGBImage == True) yuv_to_rgb(*BDat);
          
       //  cout << "put_block" << endl; 
          put_block(i);
          (*BDat).free();
       }
       InfoData.end_writing();
    }   
    else
    {
        BDat = &Dat;
        decode_data();
         // decode noise if all scale are decompressed (Resol == -1)
       if ((KeepResi) && (NoiseThresholding == True) && (Resol < 0))
       {
          read_noise_info();
          decode_residu();
       }
       if (RGBImage == True) yuv_to_rgb(*BDat);
       // dec_data();
    }

  if ((NoiseThresholding == True)
     && (Verbose == True) && (Stat_Noise == NOISE_GAUSSIAN) && (Resol < 0) )
    {
        fprintf(stderr, "Sigma estimated noise = %f\n", Noise_Ima);
    }
    
   if ((InFile != NULL) && (InFile != stdin)) fclose (InFile);
}


/*****************************************************************/

void MR_Comp3DData::init_decompress()
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

void MR_Comp3DData::decompress(int NB)
{ 
    int NumBlock=NB-1;
       
    init_decompress();
    long BufSize;

    if ((NumBlock<0) || (NumBlock >= NbrBlock)) NumBlock = 0;
    
    // skip the first blocks
    seek_block(NumBlock);
    
    // read the total size 
    if (NbrBlock > 1) BufSize = readint(InFile);
    
    BDat = &Dat;
    decode_data();
    //L_Nl = Nl;
    //L_Nc = Nc;
          
    if (Verbose == True)
         cerr << "Block " <<  NumBlock+1 << ": " <<  Nl << "x" <<  Nc << endl;
         
    // decode noise
    if ((KeepResi) && (NoiseThresholding == True) && (Resol < 0))
    {
        read_noise_info();
        decode_residu();
    }
    if (RGBImage == True) yuv_to_rgb(*BDat);
 
   if ((NoiseThresholding == True)
     && (Verbose == True) && (Stat_Noise == NOISE_GAUSSIAN) && (Resol < 0) )
         fprintf(stderr, "Sigma estimated noise = %f\n", Noise_Ima);
   
    // dec_data();
    if (InFile != NULL) fclose (InFile);
}

/***************************************************/

void MR_Comp3DData::get_resol(int RN, int &SizeBufResol, char * & BuffResol, int NB)
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
    if (NbrBlock > 1) BufSize = readint(InFile); 

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
void MR_Comp3DData::seek_block(int NumBlock)
{
    long BufSize;
    int i;
    
    if (NbrBlock > 1) 
    {
       L_Nl = readint(InFile);
       L_Nc = readint(InFile);  
       Nf = readint(InFile);
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

void MR_Comp3DData::seek_resol(int NbrSkipResol)   
{
    int NumberBand_per_Resol = number_band_per_resol(Transform);
    int NbrSkipBand = NumberBand_per_Resol*(NbrSkipResol-1)+1;
//     cout << "NumberBand_per_Resol = " << NumberBand_per_Resol << endl;
//     cout << "NbrSkipBand = " <<NbrSkipBand<< endl;
//     cout << "NbrSkipResol= " <<NbrSkipResol<< endl;
    seek_band(NbrSkipBand);
}

/***************************************************/

void MR_Comp3DData::seek_band(int NbrSkipBand)   
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


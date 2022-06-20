/******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Jean-Luc Starck
**
**    Date:  98/06/17 
**    
**    File:  IM_IO.cc
**
*******************************************************************************
**
**    DESCRIPTION  Block input output routines
**    ----------- 
**                 
**    PARAMETRES    
**    ----------    
** 
**    RESULTS      
**    -------  
**
**
**   
******************************************************************************/

#include"GlobalInc.h"
#include"IM_Obj.h"
#include "IM_IO.h"
#include "IM_BIO.h"

#if IO_DISP_OK
extern io_disp IO_DISP;
#endif

#if IO_MIDAS_OK
extern io_midas IO_MIDAS;
#endif

extern type_format Format_Imag;
extern type_data TypeInputData;

#if IO_JPEG_OK
#include "IM_JPEG.h"
#endif

#ifdef USELM
extern DemoLic MRDEMO;
#endif

 /**************************************************************************/

void IOInfoData::get_info(char *File_Name)
{  
    if (Format_Imag == F_UNKNOWN)
                     Format_Imag = io_detect_format (File_Name);
    Format = Format_Imag;
    
    switch (Format)
    {
        case F_DISP: 
#if IO_DISP_OK
              IO_DISP.read_info_ima(File_Name, Nl, Nc, Type);
  #else
              fprintf (stderr, "Error: DISP is not active\n");
              exit (-1); 
#endif
              break;
        case F_MIDAS:
#if IO_MIDAS_OK
              IO_MIDAS.read_info(File_Name, (*this));
#else
              fprintf (stderr, "Error: MIDAS is not active\n");
              exit (-1); 
#endif
              break;
        case F_FITS:
#if IO_FITS_OK
              PtrFits = new fitsstruct;
	      initfield(PtrFits);
              fits_read_header(File_Name, PtrFits);
	      if (PtrFits->naxis > 0) Nc = (PtrFits->TabAxis)[0];
	      if (PtrFits->naxis > 1) Nl = (PtrFits->TabAxis)[1];
 	      Nima = 1;
	      if (PtrFits->naxis == 3) Nima = (PtrFits->TabAxis)[2];
	      switch (PtrFits->bitpix)
	      {
	         case BP_BYTE: Type = T_BYTE; break;
		 case BP_SHORT: Type = T_SHORT; break;        
                 case BP_INT: Type = T_INT; break;           
                 case BP_FLOAT: Type = T_FLOAT; break;
                 case BP_DOUBLE: Type = T_DOUBLE; break;
		 default:
		    cerr << "Error: unknown format ... " << endl;
		    exit(-1);
		    break;
	      }
#else
              fprintf (stderr, "Error: FITS is not active\n");
              exit (-1); 
#endif
              break;
         case F_GIF: 
#if IO_GIF_OK
               readhdgif(File_Name, Nl, Nc, PtrGif);
               Type = T_BYTE;
#else
              fprintf (stderr, "Error: GIF is not active\n");
              exit (-1); 
#endif	       
               break;
         case F_PGM: 
#if IO_PGM_OK
               readhdpgm  (File_Name, Nl, Nc);
               Type = T_BYTE;
#else
              fprintf (stderr, "Error: PGM is not active\n");
              exit (-1); 
#endif
               break;        
	 case F_JPEG:
#if IO_JPEG_OK
               readhdjpeg(File_Name, Nl, Nc);
               Type = T_BYTE;
#else
              fprintf (stderr, "Error: JPEG is not active\n");
              exit (-1); 
#endif	       
               break;	 
         default:
              fprintf (stderr, "Error: bad image format. cannot read ...\n");
              exit (-1);
              break;
    }
    TypeInputData = Type;
#ifdef USELM
    MRDEMO.test(Nc, Nl);
#endif
}

/**************************************************************************/

void IOInfoData::put_info(char *File_Name)
{
    if (Format_Imag == F_UNKNOWN)
                     Format_Imag = io_detect_format (File_Name);
    Format = Format_Imag;
    
    switch (Format)
    {
        case F_DISP: 
#if IO_DISP_OK
              IO_DISP.Nl_io = Nl;
	      IO_DISP.Nc_io = Nc;
	      IO_DISP.Type = Type;
              IO_DISP.create_file(File_Name);
#else
              fprintf (stderr, "Error: DISP is not active\n");
              exit (-1); 
#endif
              break;
        case F_MIDAS:
#if IO_MIDAS_OK
              // IO_MIDAS.read_info(File_Name, (*this));
#else
              fprintf (stderr, "Error: MIDAS is not active\n");
              exit (-1); 
#endif
              break;
        case F_FITS:
#if IO_FITS_OK
             fits_write_header(File_Name, PtrFits);
#else
              fprintf (stderr, "Error: FITS is not active\n");
              exit (-1); 
#endif
              break;
         case F_GIF: 
#if IO_GIF_OK
                writehdgif (File_Name, Nl, Nc, PtrGif);
#else
              fprintf (stderr, "Error: GIF is not active\n");
              exit (-1); 
#endif
                break;
         case F_PGM: 
#if IO_PGM_OK
              writehdpgm(File_Name, Nl, Nc);
#else
              fprintf (stderr, "Error: PGM is not active\n");
              exit (-1); 
#endif
               break;        
	 case F_JPEG:
#if IO_JPEG_OK
              writehdjpeg (File_Name, Nl, Nc);
#else
              fprintf (stderr, "Error: JPEG is not active\n");
              exit (-1); 
#endif
              break;	 
         default:
              fprintf (stderr, "Error: bad image format. cannot read ...\n");
              exit (-1);
              break;
    }
}

/**************************************************************************/

void io_read_block_ima(char *File_Name, Ifloat &Image, 
                       int Indi, int Indj, IOInfoData &InfoDat, Bool NoBscale)
{
    if ((Image.nl() + Indi > InfoDat.Nl) ||
        (Image.nc() + Indj > InfoDat.Nc))
    {
       cerr << "Error: this block cannot be extracted from file: "  <<  File_Name << endl;
       cerr << "       Xs  = " <<        Indj << " Ys  = " <<  Indi       << endl;
       cerr << "       Ncb = " <<  Image.nc() << " Nlb = " <<  Image.nl() << endl;
       cerr << "       Nc  = " <<  InfoDat.Nc << " Nl  = " <<  InfoDat.Nl << endl;
    }
    
    switch (InfoDat.Format)
    {
        case F_DISP: 
#if IO_DISP_OK
              IO_DISP.read_block_ima(File_Name, Image, Indi, Indj);              
 #else
              fprintf (stderr, "Error: DISP is not active\n");
              exit (-1); 
#endif
              break;
        case F_MIDAS:
              cerr << "Error:  MIDAS block read is not implemented ... " << endl;
	      exit(-1);
              break;
        case F_FITS:
#if IO_FITS_OK
              fits_read_block(File_Name, Image, Indi, Indj, NoBscale);
#else
              fprintf (stderr, "Error: FITS is not active\n");
              exit (-1); 
#endif
              break;
         case F_GIF: 
#if IO_GIF_OK	 
               read_block_gif (Image, Indi, Indj);
#else
              fprintf (stderr, "Error: GIF is not active\n");
              exit (-1); 
#endif
               break;
         case F_PGM:
#if IO_PGM_OK	  
               read_block_pgm (File_Name, Image, Indi, Indj);
#else
              fprintf (stderr, "Error: PGM is not active\n");
              exit (-1); 
#endif
               break; 
	case F_JPEG:   
#if IO_JPEG_OK	
	       read_block_jpeg (Image, Indi, Indj);
#else
              fprintf (stderr, "Error: JPEG is not active\n");
              exit (-1); 
#endif	       
               break;
        default:
              fprintf (stderr, "Error: bad image format. cannot read ...\n");
              exit (-1);
              break;
    }
}

/**************************************************************************/

void io_read_block_ima(char *File_Name, Iint &Image, 
                       int Indi, int Indj, IOInfoData &InfoDat, Bool NoBscale)
{
    if ((Image.nl() + Indi > InfoDat.Nl) ||
        (Image.nc() + Indj > InfoDat.Nc))
    {
       cerr << "Error: this block cannot be extracted from file: "  <<  File_Name << endl;
       cerr << "       Xs  = " <<        Indj << " Ys  = " <<  Indi       << endl;
       cerr << "       Ncb = " <<  Image.nc() << " Nlb = " <<  Image.nl() << endl;
       cerr << "       Nc  = " <<  InfoDat.Nc << " Nl  = " <<  InfoDat.Nl << endl;
    }
    
    switch (InfoDat.Format)
    {
        case F_DISP: 
#if IO_DISP_OK
              IO_DISP.read_block_ima(File_Name, Image, Indi, Indj);              
 #else
              fprintf (stderr, "Error: DISP is not active\n");
              exit (-1); 
#endif
              break;
        case F_MIDAS:
              cerr << "Error:  MIDAS block read is not implemented ... " << endl;
	      exit(-1);
              break;
        case F_FITS:
#if IO_FITS_OK
              fits_read_block(File_Name, Image, Indi, Indj, NoBscale);
#else
              fprintf (stderr, "Error: FITS is not active\n");
              exit (-1); 
#endif
              break;
         case F_GIF: 
#if IO_GIF_OK	 
               read_block_gif (Image, Indi, Indj);
#else
              fprintf (stderr, "Error: GIF is not active\n");
              exit (-1); 
#endif	       
               break;
         case F_PGM: 
#if IO_PGM_OK	 
               read_block_pgm (File_Name, Image, Indi, Indj);
#else
              fprintf (stderr, "Error: PGM is not active\n");
              exit (-1); 
#endif	       
               break;      
        case F_JPEG: 
#if IO_JPEG_OK	
	       read_block_jpeg (Image, Indi, Indj);
#else
              fprintf (stderr, "Error: JPEG is not active\n");
              exit (-1); 
#endif
               break;
	default:
              fprintf (stderr, "Error: bad image format. cannot read ...\n");
              exit (-1);
              break;
    }
}

/**************************************************************************/

void io_write_block_ima(char *File_Name, Ifloat &Image, 
                       int Indi, int Indj, IOInfoData &InfoDat, Bool NoBscale)
{
    if ((Image.nl() + Indi > InfoDat.Nl) ||
        (Image.nc() + Indj > InfoDat.Nc))
    {
       cerr << "Error: this block cannot be inserted in file: "  <<  File_Name << endl;
       cerr << "       Xs  = " <<        Indj << " Ys  = " <<  Indi       << endl;
       cerr << "       Ncb = " <<  Image.nc() << " Nlb = " <<  Image.nl() << endl;
       cerr << "       Nc  = " <<  InfoDat.Nc << " Nl  = " <<  InfoDat.Nl << endl;
    }
    
    switch (InfoDat.Format)
    {
        case F_DISP: 
#if IO_DISP_OK
              IO_DISP.write_block_ima(File_Name, Image, Indi, Indj);              
 #else
              fprintf (stderr, "Error: DISP is not active\n");
              exit (-1); 
#endif
              break;
        case F_MIDAS:
              cerr << "Error:  MIDAS block read is not implemented ... " << endl;
	      exit(-1);
              break;
        case F_FITS:
#if IO_FITS_OK
              fits_write_block(File_Name, Image, Indi, Indj, NoBscale);
#else
              fprintf (stderr, "Error: FITS is not active\n");
              exit (-1); 
#endif
              break;
         case F_GIF: 
#if IO_GIF_OK
               write_block_gif (Image, Indi, Indj);
#else
              fprintf (stderr, "Error: GIF is not active\n");
              exit (-1); 
#endif
               break;
         case F_PGM: 
#if IO_PGM_OK
               write_block_pgm (File_Name, InfoDat.Nl, InfoDat.Nc, 
	                        Image, Indi, Indj);
#else
              fprintf (stderr, "Error: PGM is not active\n");
              exit (-1); 
#endif
               break; 
	case F_JPEG:  
#if IO_JPEG_OK
	       write_block_jpeg (Image, Indi, Indj);
#else
              fprintf (stderr, "Error: JPEG is not active\n");
              exit (-1); 
#endif
	       break;
        default:
              fprintf (stderr, "Error: bad image format. cannot read ...\n");
              exit (-1);
              break;
    }
}

/**************************************************************************/

void io_write_block_ima(char *File_Name, Iint &Image, 
                       int Indi, int Indj, IOInfoData &InfoDat, Bool NoBscale)
{
    if ((Image.nl() + Indi > InfoDat.Nl) ||
        (Image.nc() + Indj > InfoDat.Nc))
    {
       cerr << "Error: this block cannot be extracted from file: "  <<  File_Name << endl;
       cerr << "       Xs  = " <<        Indj << " Ys  = " <<  Indi       << endl;
       cerr << "       Ncb = " <<  Image.nc() << " Nlb = " <<  Image.nl() << endl;
       cerr << "       Nc  = " <<  InfoDat.Nc << " Nl  = " <<  InfoDat.Nl << endl;
    }
    
    switch (InfoDat.Format)
    {
        case F_DISP: 
#if IO_DISP_OK
              IO_DISP.write_block_ima(File_Name, Image, Indi, Indj);              
 #else
              fprintf (stderr, "Error: DISP is not active\n");
              exit (-1); 
#endif
              break;
        case F_MIDAS:
              cerr << "Error:  MIDAS block read is not implemented ... " << endl;
	      exit(-1);
              break;
        case F_FITS:
#if IO_FITS_OK
              fits_write_block(File_Name, Image, Indi, Indj, NoBscale);
#else
              fprintf (stderr, "Error: FITS is not active\n");
              exit (-1); 
#endif
              break;
         case F_GIF: 
#if IO_GIF_OK	 
                write_block_gif (Image, Indi, Indj);
#else
              fprintf (stderr, "Error: GIF is not active\n");
              exit (-1); 
#endif		
               break;
         case F_PGM: 
#if IO_PGM_OK	 
               write_block_pgm (File_Name, InfoDat.Nl, InfoDat.Nc, Image, Indi, Indj);
#else
              fprintf (stderr, "Error: PGM is not active\n");
              exit (-1); 
#endif
               break;      
        case F_JPEG: 
#if IO_JPEG_OK	
	       write_block_jpeg (Image, Indi, Indj);
#else
              fprintf (stderr, "Error: JPEG is not active\n");
              exit (-1); 
#endif	       
               break;
	default:
              fprintf (stderr, "Error: bad image format. cannot read ...\n");
              exit (-1);
              break;
    }
}


/**************************************************************************/

void io_close_block_ima(char *File_Name, IOInfoData &InfoDat)
{
    switch (InfoDat.Format)
    {
        case F_DISP: 
        case F_MIDAS:
        case F_FITS:
        case F_PGM:
	     break;
        case F_JPEG:
#if IO_JPEG_OK		
	     closejpeg(File_Name, InfoDat.PtrGif);
#else
              fprintf (stderr, "Error: JPEG is not active\n");
              exit (-1); 
#endif
	   break;
         case F_GIF: 
#if IO_GIF_OK	 
             closegif (File_Name, InfoDat.PtrGif);
#else
              fprintf (stderr, "Error: GIF is not active\n");
              exit (-1); 
#endif		
               break;
 	default:
              fprintf (stderr, "Error: bad image format. cannot read ...\n");
              exit (-1);
              break;
    }
}

/**************************************************************************/

void LIfloat::close()
{
  io_close_block_ima(FileName, Info);
}

void LIfloat::init(char *FileImageName, int BSize)
{
   Info.get_info(FileImageName);
   BlockSize = BSize;
   L_Nl = Info.Nl;
   L_Nc = Info.Nc;
   set_tabblock();
   FileName = strdup(FileImageName);
}

void LIfloat::create(char *FileImageName, int NbrL, int NbrC, int BSize)
{
   L_Nl = NbrL;
   L_Nc = NbrC;
   BlockSize = BSize;
   set_tabblock();
   Info.Nl = L_Nl;
   Info.Nc = L_Nc;
   if (Info.Format == F_UNKNOWN)
              Info.Format = io_detect_format (FileImageName);
   if (Info.Type == UNKNOWN) Info.Type = T_FLOAT;
   if ((Info.Format == F_FITS) && (Info.PtrFits == NULL))
   {
       Info.PtrFits = new fitsstruct;
       init_fits_struct(Info.PtrFits, L_Nl, L_Nc);
   }
   else if ((Info.Format == F_GIF) && (Info.PtrGif == NULL))
   {
       Info.PtrGif = new PICINFO;
   }
   Info.put_info(FileImageName);
   FileName = strdup(FileImageName);
}

void LIfloat::create(char *FileImageName, IOInfoData &IData, int BSize)
{
   L_Nl =  IData.Nl;
   L_Nc =  IData.Nc;
   BlockSize = BSize;
   set_tabblock();
   Info.Format = IData.Format;
   Info.Type = IData.Type;
   Info.Nl = IData.Nl;
   Info.Nc = IData.Nc;
   Info.Nima = IData.Nima;
   if ((Info.Format == F_FITS) && (Info.PtrFits == NULL))
   {
       Info.PtrFits = new fitsstruct;
       init_fits_struct(Info.PtrFits, L_Nl, L_Nc);
   }
   else if ((Info.Format == F_GIF) && (Info.PtrGif == NULL))
   {
       Info.PtrGif = new PICINFO;
   }   
   Info.put_info(FileImageName);
   FileName = strdup(FileImageName);
}

LIfloat::~LIfloat()
{
   if (TabBlockNl != NULL) delete  [] TabBlockNl;
   if (TabBlockNc != NULL) delete [] TabBlockNc;
   if (TabBlockDepNl  != NULL) delete [] TabBlockDepNl;
   if (TabBlockDepNc != NULL) delete [] TabBlockDepNc;
   if (FileName != NULL) free ((char *) FileName);
   reset_param();
}

void LIfloat::set_block(int Nb)
{
   if (NbrBlock < 1)
   {
      cerr << "Error: image not initialized ... " << endl;
      exit(-1);
   }
   if ((Nb < 0) || (Nb >=  NbrBlock))
   {
      cerr << "Error: bad block number ... " << endl;
      exit(-1);
   }
   
   if (Nb != NumBlock)
   {
       NumBlock = Nb;
       Nl = TabBlockNl[NumBlock];
       Nc = TabBlockNc[NumBlock];
       Depi = TabBlockDepNl[NumBlock];
       Depj = TabBlockDepNc[NumBlock];
       Block.resize(Nl,Nc);
   }
}

void LIfloat::get_block(int Nb)
{
   if (NbrBlock < 1)
   {
      cerr << "Error: image not initialized ... " << endl;
      exit(-1);
   }
   if ((Nb < 0) || (Nb >=  NbrBlock))
   {
      cerr << "Error: bad block number ... " << endl;
      exit(-1);
   }
   
   if (Nb != NumBlock)
   {
       NumBlock = Nb;
       Nl = TabBlockNl[NumBlock];
       Nc = TabBlockNc[NumBlock];
       Depi = TabBlockDepNl[NumBlock];
       Depj = TabBlockDepNc[NumBlock];
       Block.resize(Nl,Nc);
       io_read_block_ima(FileName, Block, Depi, Depj, Info, NoBscale);
   }
}
       

void LIfloat::put_block()
{
   if (NbrBlock < 1)
   {
      cerr << "Error: image not initialized ... " << endl;
      exit(-1);
   }
   io_write_block_ima(FileName, Block, Depi, Depj, Info, NoBscale);
}


void LIfloat::set_tabblock()
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
      //      " DepNl = " << TabBlockDepNl[i] << " DepNc = " << TabBlockDepNc[i] << endl;
              
   }
}

/**************************************************************************/


void LIint::close()
{
  io_close_block_ima(FileName, Info);
}

void LIint::init(char *FileImageName, int BSize)
{
   Info.get_info(FileImageName);
   BlockSize = BSize;
   L_Nl = Info.Nl;
   L_Nc = Info.Nc;
   set_tabblock();
   FileName = strdup(FileImageName);
}

void LIint::create(char *FileImageName, int NbrL, int NbrC, int BSize)
{
   L_Nl = NbrL;
   L_Nc = NbrC;
   BlockSize = BSize;
   set_tabblock();
   Info.Nl = L_Nl;
   Info.Nc = L_Nc;
   if (Info.Format == F_UNKNOWN)
              Info.Format = io_detect_format (FileImageName);
   if (Info.Type == UNKNOWN) Info.Type = T_FLOAT;
   if ((Info.Format == F_FITS) && (Info.PtrFits == NULL))
   {
       Info.PtrFits = new fitsstruct;
       init_fits_struct(Info.PtrFits, L_Nl, L_Nc);
   }
   else if ((Info.Format == F_GIF) && (Info.PtrGif == NULL))
   {
       Info.PtrGif = new PICINFO;
   }
   Info.put_info(FileImageName);
   FileName = strdup(FileImageName);
}

void LIint::create(char *FileImageName, IOInfoData &IData, int BSize)
{
   L_Nl =  IData.Nl;
   L_Nc =  IData.Nc;
   BlockSize = BSize;
   set_tabblock();
   Info.Format = IData.Format;
   Info.Type = IData.Type;
   Info.Nl = IData.Nl;
   Info.Nc = IData.Nc;
   Info.Nima = IData.Nima;
   if ((Info.Format == F_FITS) && (Info.PtrFits == NULL))
   {
       Info.PtrFits = new fitsstruct;
       init_fits_struct(Info.PtrFits, L_Nl, L_Nc);
   }
   else if ((Info.Format == F_GIF) && (Info.PtrGif == NULL))
   {
       Info.PtrGif = new PICINFO;
   }   
   Info.put_info(FileImageName);
   FileName = strdup(FileImageName);
}

LIint::~LIint()
{
   if (TabBlockNl != NULL) delete  [] TabBlockNl;
   if (TabBlockNc != NULL) delete [] TabBlockNc;
   if (TabBlockDepNl  != NULL) delete [] TabBlockDepNl;
   if (TabBlockDepNc != NULL) delete [] TabBlockDepNc;
   if (FileName != NULL) free ((char *) FileName);
   reset_param();
}

void LIint::set_block(int Nb)
{
   if (NbrBlock < 1)
   {
      cerr << "Error: image not initialized ... " << endl;
      exit(-1);
   }
   if ((Nb < 0) || (Nb >=  NbrBlock))
   {
      cerr << "Error: bad block number ... " << endl;
      exit(-1);
   }
   
   if (Nb != NumBlock)
   {
       NumBlock = Nb;
       Nl = TabBlockNl[NumBlock];
       Nc = TabBlockNc[NumBlock];
       Depi = TabBlockDepNl[NumBlock];
       Depj = TabBlockDepNc[NumBlock];
       Block.resize(Nl,Nc);
   }
}

void LIint::get_block(int Nb)
{
   if (NbrBlock < 1)
   {
      cerr << "Error: image not initialized ... " << endl;
      exit(-1);
   }
   if ((Nb < 0) || (Nb >=  NbrBlock))
   {
      cerr << "Error: bad block number ... " << endl;
      exit(-1);
   }
   
   if (Nb != NumBlock)
   {
       NumBlock = Nb;
       Nl = TabBlockNl[NumBlock];
       Nc = TabBlockNc[NumBlock];
       Depi = TabBlockDepNl[NumBlock];
       Depj = TabBlockDepNc[NumBlock];
       Block.resize(Nl,Nc);
       io_read_block_ima(FileName, Block, Depi, Depj, Info, NoBscale);
   }
}
       

void LIint::put_block()
{
   if (NbrBlock < 1)
   {
      cerr << "Error: image not initialized ... " << endl;
      exit(-1);
   }
   io_write_block_ima(FileName, Block, Depi, Depj, Info, NoBscale);
}


void LIint::set_tabblock()
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
      //      " DepNl = " << TabBlockDepNl[i] << " DepNc = " << TabBlockDepNc[i] << endl;
              
   }
}

/*********************************************************************/

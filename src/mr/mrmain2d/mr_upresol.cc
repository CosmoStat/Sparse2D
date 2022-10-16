/******************************************************************************
**                   Copyright (C) 1995 CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 3.3
**
**    Author: Jean-Luc Starck
**
**    Date:  96/05/07 
**    
**    File:  mr_decomp.cc
**
*******************************************************************************
**
**    DECRIPTION  decompression program
**    ---------- 
**
*******************************************************************************
**
**    PARAMETRES
**    ----------
**
**        Input_File: image file name to compress
**        [Output_File]: compressed image file
**                       by default, "Input_File.P"
**           [-v] : verbose
**
**		"Usage: %s [-v] [-r resolution] input output\n",
**        [-v] : verbose
**        [-r resolution] : we are always interested to decompress
**                          the image at full resolution. This 
**                          parameter allows to extract from the
**                          compressed file, the image at worse 
**                          resolution, and produces by this way 
**                          a smaller image than the original
**                          Even is the compressed file contained
**                          the noise of the input image, the noise
**                          will not be decompressed.
**                          resolution must be >= 0 and < number of scales
**                          of the transform
**                          be default, the image is reconstructed at
**                          full resolution with its noise if it exists
**                          in the compressed file.
**         input : compressed file name
**        output : decompressed file name
**      
**        [-s] : the decompressed image is send to the 
**                standard output (screen). 
**
**        [-i] : if this option is not set, input 
**                data are read from the standard input. 
**
**        [-o] : if this option is set, the decompressed
**                image is written to a file. If not, the
**                -s option has to be set\n");
**
**         [-t] output type\n");
**                if the input image was a fits image, 
**                the image output type can be fixed to 
**                by the user to 'i' for integer, or 's' 
**                for short. By default, the output type 
**                is the same as the type of the original image
**
**         [-g] \n");
**                add a simulated noise to the decompressed
**                image with the same properties as in the
**                original image. So they look very similar.
**
******************************************************************************/

// static char sccsid[] = "@(#)mr_decomp.cc 3.3 96/05/07 CEA 1995 @(#)";

#include "IM_Obj.h"
#include "IM_IO.h"
#include "MR_Obj.h"
#include "IM_CompTool.h"
#include "IM_Noise.h"
#include "IM_Comp.h"
#include "MR_Comp.h"

int Resol=-1;
char File_Name_Imag[80], File_Name_Transform[80];
char File_Name_Resol[80];
char InfoFile[80];
char InfoFileOut[80];

extern int  OptInd;
extern char *OptArg;

extern int  GetOpt(int argc, char *const*argv, char *opts);

extern Bool InputFromStdin;
int IterRec=0;
Bool LowResol = False;
char OutputType='u';
int Verbose=0;
int Zy=0;
int Zx=0;
int NlVisu=256;
int NcVisu=256;
Bool RealUnit = False;
int RY=0;
int RX=0;

/*********************************************************************/

static void usage(char *argv[])
{
    fprintf(OUTMAN, "\n\nUsage: %s options CompressedFileName OutputFileName\n\n", argv[0]);
    fprintf(OUTMAN, "   where options are = \n");

    manline();
    fprintf(OUTMAN, "         [-t] output type\n");
    fprintf(OUTMAN, "              if the input image was a fits image, \n");
    fprintf(OUTMAN, "              the image output type can be fixed by the user \n");
    fprintf(OUTMAN, "              to 'f' for float, to 'i' for integer, or 's' \n");
    fprintf(OUTMAN, "              for short. By default, the output type \n");
    fprintf(OUTMAN, "              is the same as the type of the original image. \n");
    manline();  
    
    fprintf(OUTMAN, "         [-I IterRecNbr]\n");
    fprintf(OUTMAN, "              Use an iterative reconstruction. \n");
    fprintf(OUTMAN, "              Only used with orthogonal transform.\n");
    fprintf(OUTMAN, "              Default is no iteration.\n");
    
    manline();
    fprintf(OUTMAN, "         [-l LowResol_Image_FileName]\n");
    fprintf(OUTMAN, "              Image at a lower resolution. \n");
    manline();
    
    fprintf(OUTMAN, "         [-x PosX]\n");
    fprintf(OUTMAN, "              Center pixel position (x) of the region to zoom \n");
    fprintf(OUTMAN, "              in the input low resolution image.\n");
    fprintf(OUTMAN, "              Default value is 0.\n");
    manline();
    
    fprintf(OUTMAN, "         [-y PosY]\n");
    fprintf(OUTMAN, "              Center pixel position (y) of the region to zoom \n");
    fprintf(OUTMAN, "              in the input low resolution image.\n");
    fprintf(OUTMAN, "              Default value is 0.\n");
    manline();

    fprintf(OUTMAN, "         [-X PosX]\n");
    fprintf(OUTMAN, "              Center pixel position (x) of the region to zoom \n");
    fprintf(OUTMAN, "              in the full resolution image (real pixel coordinate).\n");
    fprintf(OUTMAN, "              Default value is 0.\n");
    manline();
    
    fprintf(OUTMAN, "         [-Y PosY]\n");
    fprintf(OUTMAN, "              Center pixel position (y) of the region to zoom \n");
    fprintf(OUTMAN, "              in the full resolution image (real pixel coordinate).\n");
    fprintf(OUTMAN, "              Default value is 0.\n");
    manline();
            
    fprintf(OUTMAN, "         [-W WindowSize]\n");
    fprintf(OUTMAN, "              Window size. Default is %d.\n", NlVisu);
    
    manline();
    verbose_usage();
    manline();
    manline();
    manline(); 
    exit(-1);
}

/*********************************************************************/

/* GET COMMAND LINE ARGUMENTS */
static void hcinit(int argc, char *argv[])
{
    int c, L;
    char *Ptr;
    Bool Optxy = False;
     
    Verbose = 0;
    strcpy(File_Name_Transform, "stdin");
    InputFromStdin = True;

    /* get options */
    while ((c = GetOpt(argc,argv,"vl:I:t:W:x:y:X:Y:")) != -1) 
    {
	switch (c) 
        {
           case 'I': if (sscanf(OptArg,"%d",&IterRec) != 1) 
                     {
		    fprintf(OUTMAN, "Error: bad iteration number parameter: %s\n", OptArg);
		    usage(argv);
		     } 
		     break;
           case 'W': if (sscanf(OptArg,"%d",&NlVisu) != 1) 
                     {
		        fprintf(OUTMAN, "Error: bad window size: %s\n", OptArg);
		        exit(-1);
		     }
		     if (NlVisu < 16)
		     {
 		        fprintf(OUTMAN, "Error: bad window size: %d\n", NlVisu);
			fprintf(OUTMAN, "       minimum value is 16.\n");
		        exit(-1);
		     }
		     NcVisu = NlVisu;
		     break;
	   case 'x': if (sscanf(OptArg,"%d",&Zx) != 1) 
                     {
		        fprintf(OUTMAN, "Error: bad x position: %s\n", OptArg);
		        exit(-1);
		     }
		     if (Zx < 0)
		     {
		       fprintf(OUTMAN, "Error: bad x position: %d\n", Zx);
		       exit(-1);
		     }
 		     Optxy = True;
		     break;
	   case 'y': if (sscanf(OptArg,"%d",&Zy) != 1) 
                     {
		        fprintf(OUTMAN, "Error: bad y position: %s\n", OptArg);
		         usage(argv);
		     }
		     if (Zx < 0)
		     {
		       fprintf(OUTMAN, "Error: bad y position: %d\n", Zy);
		       exit(-1);
		     }
		     Optxy = True;
 		     break;
	  case 'X': if (sscanf(OptArg,"%d",&RX) != 1) 
                     {
		        fprintf(OUTMAN, "Error: bad X position: %s\n", OptArg);
		         usage(argv);
		     }
		     if (RX < 0)
		     {
		       fprintf(OUTMAN, "Error: bad X position: %d\n", RX);
		       exit(-1);
		     }
		     RealUnit = True;
 		     break;
	  case 'Y': if (sscanf(OptArg,"%d",&RY) != 1) 
                     {
		        fprintf(OUTMAN, "Error: bad Y position: %s\n", OptArg);
		         usage(argv);
		     }
		     if (RY < 0)
		     {
		       fprintf(OUTMAN, "Error: bad Y position: %d\n", RY);
		       exit(-1);
		     }
		     RealUnit = True;
 		     break;	     
 	  case 'v':
		/* verbose flag -v */
		Verbose = 1;
		break;
           case 'l':  
		if (sscanf(OptArg,"%s", File_Name_Resol) != 1) 
                {
		    fprintf(OUTMAN, "Error: bad input image: %s\n", OptArg);
		    exit(-1);
		} 
		LowResol = True;
		break;
	   case 't':
               /* output can be integer or short */
		if (sscanf(OptArg,"%c",&OutputType) != 1) 
                {
		   fprintf(OUTMAN, "Error: bad output type: %s\n", OptArg);
		   exit(-1);
		}
		if ((OutputType != 'i') && (OutputType != 's') && (OutputType != 'f'))
		{
		    fprintf(OUTMAN, "Error: bad output type: %s\n", OptArg);
		    fprintf(OUTMAN, "       output type has to equal to 'i','s' or 'f'\n");
		    exit(-1);
		}
		break;
	   case '?':
		usage(argv);
		break;
	}
    }
    if ((Optxy == True) && (RealUnit == True))
    {
      fprintf(OUTMAN, "Error: x,y options are not compatible with X and Y options...\n");
      exit(-1);
    }
    if ((LowResol == False) && (RealUnit == True))
    {
        fprintf(OUTMAN, "Warning: X,Y parameters will not be used ...\n");
    }
    if ((LowResol == False) && (Optxy == True))
    {
        fprintf(OUTMAN, "Warning: x,y parameters will not be used ...\n");
    }
    
   /*  MRC file name */
    if (OptInd < argc) strcpy(File_Name_Transform, argv[OptInd++]);
    else usage(argv);
    L = strlen (File_Name_Transform);
    if ((L == 1) && (File_Name_Transform[0] == '-')) InputFromStdin = True;
    else 
    {
       InputFromStdin = False;
       Ptr = File_Name_Transform;
       if (L < 5) strcat (File_Name_Transform, ".MRC");
       else if ((Ptr[L-1] != 'C') || (Ptr[L-2] != 'R') || (Ptr[L-3] != 'M')) 
                         strcat (File_Name_Transform, ".MRC");
    }		
    if (InputFromStdin == True)
    {
       cout << " Error: input file cannot be read from stdin ..." << endl;
       exit(-1);
    }
    /* output file name */		
    if (OptInd < argc) strcpy(File_Name_Imag, argv[OptInd++]);
    else usage(argv);   
    
    /* make sure there are not too many parameters */
    if (OptInd < argc)
    {
	fprintf(OUTMAN, "Too many parameters: %s ...\n", argv[OptInd]);
	usage(argv);
    }
}

/***************************************************/

int main (int argc, char *argv[])
{
    int k;
    char Cmd[512];
    fitsstruct Header;
	
    Cmd[0] = '\0';
    for (k =0; k < argc; k++) sprintf(Cmd, "%s %s", Cmd, argv[k]);
    /* Get command line arguments*/

    hcinit(argc, argv);

 MR_CompData DecIma;
    DecIma.IterRec = IterRec;
    DecIma.Cmd = Cmd;
    DecIma.File_Name_Imag=File_Name_Imag;
    DecIma.File_Name_Transform=File_Name_Transform;
    DecIma.OutputType=OutputType;
    if (Verbose) DecIma.Verbose = True;

    char *ptr;
    strcpy(InfoFileOut, File_Name_Imag);
    ptr = strrchr(InfoFileOut, '.');
    if (ptr != NULL) *ptr = '\0';
    strcat(InfoFileOut,".inf");
      
     Ifloat ImagLowResol;
     if (LowResol == True)
     {
        strcpy(InfoFile, File_Name_Resol);
	ptr = strrchr(InfoFile, '.');
        if (ptr != NULL) *ptr = '\0';
        strcat(InfoFile,".inf");
        io_read_ima_float(File_Name_Resol, ImagLowResol, &Header);
        BlockMRInfo BV(InfoFile);
 	if ((BV.nl() != ImagLowResol.nl()) || (BV.nc() != ImagLowResol.nc()))
	{
	   cerr << "Error: low resolution image size is not compatible with the description file ... " << endl;
           cerr << "IF.Nl = " << BV.nl() << " IF.Nc = " << BV.nc() << endl;
	   cerr << "LR.Nl = " <<  ImagLowResol.nl() << " LR.Nc = " <<  ImagLowResol.nc() << endl;
 	   exit(-1);
	}
	
	// Test if the full resolution has already been achived
	if (BV.resol() < 0)
	{
	   cerr << "Warning: the low resolution image is already at " << endl;
           cerr << "         the full resolution ..." << endl;
	   cerr << "         output image is set to the input image." << endl;
           BV.write(InfoFileOut);
           io_write_ima_float(File_Name_Imag, ImagLowResol, &Header);	   
 	}
	else
	{ 
	   if (RealUnit == True)
	   {
	      if ((RX < 0) || (RX >= BV.ima_nc()) ||
	          (RY < 0) || (RY >= BV.ima_nl()))
	      {
	         cerr << "Error: bad coordinates ... " << endl;
		 cerr << "       image size is : " << endl;
		 cerr << "              Nl = " << BV.ima_nl() << endl;
		 cerr << "              Nc = " << BV.ima_nc() << endl;
		 exit(0);
	      }
              BlockMRInfo BI(BV.ima_nl(), BV.ima_nc(), BV.block_size());
  	      int Bi=0;
	      int Bj=0;
 	      while ((RY > BI.block_nl(Bi,Bj)) && (Bi < BI.nbr_block_nl()))
	      {
	         RY -= BI.block_nl(Bi,Bj);
	         Bi++;
	      } 
 	      while ((RX > BI.block_nc(Bi,Bj)) && (Bj < BI.nbr_block_nc()))  
	      {
	         RX -= BI.block_nc(Bi,Bj);
	         Bj++;
	      } 
	      if (BV.resol() > 0)
	      {
	         RY = BV.pos_nl(Bi,Bj) + (int) (RY / pow(2., (double) BV.resol()) + 0.5);
		 RX = BV.pos_nc(Bi,Bj) + (int) (RX / pow(2., (double) BV.resol()) + 0.5);
	      }
	      Zx = RX;
	      Zy = RY;
	      if (Verbose == True)
	      {
	         cout << "Pixel position in the low resolution image: " << endl;
		 cout << "              X = " << Zx << endl;
		 cout << "              Y = " << Zy << endl;
              }
  	   }
	   
           BlockMRInfo BVup(Zy,Zx, BV, NlVisu, NcVisu);  
	   // open the file and read the header
	   DecIma.Resol = BVup.resol();
           DecIma.init_decompress();
	   
	   // Test if the residual exists
           if ((BV.resol() == 0) &&
 	       ((noise_threshold(DecIma.Comp_Method) == False) || (DecIma.KeepResi==0)))
	   {
	      cerr << "Warning: residual part cannot be read ... " << endl;
              cerr << "         output image is set to the input image." << endl;
	      BV.write(InfoFileOut);
              io_write_ima_float(File_Name_Imag, ImagLowResol, &Header);
	   }
 	   else
	   {
   	      if (BV.nbr_block() > 1) DecIma.up_mrresol(ImagLowResol, BV, BVup);
	      else DecIma.up_mrresol(ImagLowResol);
              BVup.write(InfoFileOut);
	      DecIma.write();
 	   }    
	   // if (DecIma.InFile != NULL) fclose (DecIma.InFile);
  	}
     }
     else 
     {
        // decompress the last scale
        DecIma.Resol = 20;
        DecIma.decompress();
	
        if (DecIma.NbrBlock <= 1) DecIma.write();
	//cout << "Block info:  " << DecIma.L_Nl << endl;

        BlockMRInfo BI(DecIma.L_Nl,  DecIma.L_Nc, DecIma.BlockSize);
// cout << "Block info " << DecIma.KeepResi <<  endl;
        BI.KeepResi = DecIma.KeepResi;
        BlockMRInfo BV(BI,DecIma.Nbr_Plan-1);
// cout << "Block info " <<  BV.KeepResi <<  endl;

        BV.write(InfoFileOut);
     }
     
     exit(0);
}



 



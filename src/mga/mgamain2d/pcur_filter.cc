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
**    Date:  03/02/00
**    
**    File:  cur_filter.cc
**
*******************************************************************************
**
**    DESCRIPTION  curvelet filtering program
**    ----------- 
**                 
******************************************************************************/

#include "IM_Obj.h"
#include "IM_IO.h"
#include "Usage.h"
#include "PCur.h"
#include "FCur.h"

char Name_Imag_In[256];  // input file image  
char Name_Imag_Out[256]; //  output file name  
 
extern int  OptInd;
extern char *OptArg;

extern int  GetOpt(int argc, char **argv, char *opts);

Bool Verbose=False; 
int NbrScale2D = DEF_CUR_NBR_SCALE_2D;
float Noise_Ima=0.;
float N_Sigma=DEFAULT_N_SIGMA;
Bool Simu=False;
int BlockSize=DEF_CUR_BLOCK_SIZE;
int NbrScale1D = -1;
Bool BlockOverlap=True;
type_ridgelet_WTtrans RidTrans = DEF_RID_TRANS;
Bool PositivSol=True;
Bool ColTrans = False;
type_curvelet_block TypeBlock = DEF_TYPE_CUR_BLOCK;

char NoiseFile[256];  /*  list of input image names  */
Bool NoiseIma = False;
#ifdef DO_SIMU
#include"Simu.cc"
#endif

type_transform Transform = TO_PYR_FFT_DIFF_RESOL; // TM_PYR_MEDIAN;
Bool MIX = False;
Bool FastCur = False;

Bool ExtendWT=False;
Bool IsotropWT=False;
Bool RecWT=False;

/***************************************/

static void usage(char *argv[])
{
    fprintf(OUTMAN, "Usage: %s options in_image result\n\n", argv[0]);
    fprintf(OUTMAN, "   where options =  \n");

    fprintf(OUTMAN, "         [-t type_of_ridgelet]\n");
    for (int i = 0; i < NBR_RID_DECIMATED; i++)
    fprintf(OUTMAN, "              %d: %s \n",i+1, StringRidTransform((type_ridgelet_WTtrans) (i+1)));
    fprintf(OUTMAN, "              Default is %s.\n",  StringRidTransform(RidTrans));
    manline();

    fprintf(OUTMAN, "         [-n number_of_scales]\n");
    fprintf(OUTMAN, "             Number of scales used in the wavelet transform.\n");
    fprintf(OUTMAN, "             default is %d. \n", NbrScale2D);    
    manline();

    fprintf(OUTMAN, "         [-N number_of_scales]\n");
    fprintf(OUTMAN, "             Number of scales used in the ridgelet transform.\n");
    fprintf(OUTMAN, "             default is automatically calculated. \n");    
    manline();
 
    fprintf(OUTMAN, "         [-b BlockSize]\n");
    fprintf(OUTMAN, "             Block Size.\n");
    fprintf(OUTMAN, "             default is %d. \n", BlockSize);    
    manline();

//     fprintf(OUTMAN, "         [-B type_of_curvelet_block]\n");
//     for (int i = 0; i < NBR_CUR_TYPE_BLOCK; i++)
//     fprintf(OUTMAN, "              %d: %s \n",i+1,  StringBlockType((type_curvelet_block) (i)));
//     fprintf(OUTMAN, "              Default is %s.\n",   StringBlockType(TypeBlock));
//     manline();
     
    gauss_usage();
    manline();  
     
    nsigma_usage(N_Sigma);
    manline();

    fprintf(OUTMAN, "         [-O]\n");
    fprintf(OUTMAN, "             Do not apply block overlapping. \n"); 
    fprintf(OUTMAN, "             By default, block overlapping is used. \n");    
    manline();
    
    fprintf(OUTMAN, "         [-P]\n");
    fprintf(OUTMAN, "             Supress the positivity constraint. Default is no. \n");    
    manline();
    
    fprintf(OUTMAN, "         [-I NoiseFileName]\n");
    manline();
    
    vm_usage();
    manline();    
    verbose_usage();
    manline();         
    manline();
    exit(-1);
}
  
/*********************************************************************/

/* GET COMMAND LINE ARGUMENTS */
static void filtinit(int argc, char *argv[])
{
    int c;
#ifdef LARGE_BUFF
    int VMSSize=-1;
    Bool OptZ = False;
    char VMSName[1024] = "";
#endif     

    /* get options */
    while ((c = GetOpt(argc,argv,"reifxT:I:B:CPt:Ob:Ss:g:n:N:vzZ")) != -1) 
    {
	switch (c) 
        { 
	   case 'r':  RecWT = (RecWT== True) ? False: True; break;
	   case 'e':  ExtendWT = (ExtendWT== True) ? False: True; break;
	   case 'i':  IsotropWT = (IsotropWT== True) ? False: True; break;
	   case 'f': FastCur= True; break;
	   case 'x': MIX = True; break;
	   case 'T':
		/* -d <type> type of transform */
		if (sscanf(OptArg,"%d",&c ) != 1) 
                {
		    fprintf(OUTMAN, "bad type of multiresolution transform: %s\n", OptArg);
	            exit(-1);
                    
		}
                if ((c > 0) && (c <= NBR_TRANSFORM+1)) 
                                        Transform = (type_transform) (c-1);
                else  
                {
		    fprintf(OUTMAN, "bad type of transform: %s\n", OptArg);
	            exit(-1);
 		}                
 		break;
          case 'I':
	        if (sscanf(OptArg,"%s",NoiseFile) != 1) 
                {
                    fprintf(OUTMAN, "bad file name: %s\n", OptArg);
                    exit(-1);
                }
		NoiseIma = True;
	        break;
	   case 'C':  ColTrans =(ColTrans  == True) ? False: True; break; 
           case 'P': PositivSol = (PositivSol == True) ? False: True; break;
            case 't': 
              if (sscanf(OptArg,"%d",&c ) != 1) 
                {
                    fprintf(OUTMAN, "Error: bad type of ridgelet transform: %s\n", OptArg);
                    exit(-1);
                }
                if ((c > 0) && (c <= NBR_RID_TRANS)) RidTrans = (type_ridgelet_WTtrans) (c);
                else  
                {
                   fprintf(OUTMAN, "Error: bad type of ridgelet transform: %s\n", OptArg);
                   exit(-1);
                }
                break;
	  case 'B': 
              if (sscanf(OptArg,"%d",&c ) != 1) 
                {
                    fprintf(OUTMAN, "Error: bad type of block: %s\n", OptArg);
                    exit(-1);
                }
                if ((c > 0) && (c <= NBR_CUR_TYPE_BLOCK)) TypeBlock = (type_curvelet_block) (c-1);
                else  
                {
                   fprintf(OUTMAN, "Error: bad type of block: %s\n", OptArg);
                   exit(-1);
                }
                break;
           case 'O': BlockOverlap =(BlockOverlap == True) ? False: True; break; 
           case 'S': Simu=True;break;
	   case 'b':
                if (sscanf(OptArg,"%d",&BlockSize) != 1) 
                {
                    fprintf(OUTMAN, "bad block size: %s\n", OptArg);
                    exit(-1);
                }
                break;
           case 'g':
                /* -g <sigma_noise> */
                if (sscanf(OptArg,"%f",&Noise_Ima) != 1) 
                {
                    fprintf(OUTMAN, "bad sigma noise: %s\n", OptArg);
                    usage(argv);
                }
                break;
          case 's':
                /* -s <nsigma> */
                if (sscanf(OptArg,"%f",&N_Sigma) != 1) 
                {
                    fprintf(OUTMAN, "bad N_Sigma: %s\n", OptArg);
                    usage(argv);
                }
                if (N_Sigma <= 0.) N_Sigma = DEFAULT_N_SIGMA;
                break;
           case 'v': Verbose = True;break;
	   case 'n':
                /* -n <NbrScale2D> */
                if (sscanf(OptArg,"%d",&NbrScale2D) != 1) 
                {
                    fprintf(OUTMAN, "bad number of scales: %s\n", OptArg);
                    exit(-1);
                }
                if (NbrScale2D > MAX_SCALE)
                 {
                    fprintf(OUTMAN, "bad number of scales: %s\n", OptArg);
                    fprintf(OUTMAN, " Nbr Scales <= %d\n", MAX_SCALE);
                    exit(-1);
                }
                break;
	   case 'N':
                 if (sscanf(OptArg,"%d",&NbrScale1D) != 1) 
                {
                    fprintf(OUTMAN, "bad number of scales: %s\n", OptArg);
                    exit(-1);
                }
                break;
#ifdef LARGE_BUFF
	    case 'z':
	        if (OptZ == True)
		{
                   fprintf(OUTMAN, "Error: Z option already set...\n");
                   exit(-1);
                }
	        OptZ = True;
	        break;
            case 'Z':
	        if (sscanf(OptArg,"%d:%s",&VMSSize, VMSName) < 1)
		{
		   fprintf(OUTMAN, "Error: syntaxe is Size:Directory ... \n");
		   exit(-1);
		}
	        if (OptZ == True)
		{
                   fprintf(OUTMAN, "Error: z option already set...\n");
                   exit(-1);
                }
		OptZ = True;
                break;
#endif
            case '?': usage(argv); break;
	    default: usage(argv); break;
 		}
	} 
      

       /* get optional input file names from trailing 
          parameters and open files */
       if (OptInd < argc) strcpy(Name_Imag_In, argv[OptInd++]);
         else usage(argv);

	if (OptInd < argc) strcpy(Name_Imag_Out, argv[OptInd++]);
         else usage(argv);

	/* make sure there are not too many parameters */
	if (OptInd < argc)
        {
		fprintf(OUTMAN, "Too many parameters: %s ...\n", argv[OptInd]);
		exit(-1);
	}

       	   
#ifdef LARGE_BUFF
    if (OptZ == True) vms_init(VMSSize, VMSName, Verbose);
#endif  
}


/***************************************************************************/

void PyrCurvelet::mix_threshold_one_block(int s2d, Ifloat &ImaBlock, 
                                   Ifloat &RidIma, float Nsig, float SigmaNoise)
{
    int BS = TabRidgelet[s2d].block_size();
    TabRidgelet[s2d].transform_one_block(ImaBlock,RidIma);	   
    for (int s1d=0; s1d < TabRidgelet[s2d].NbrScale; s1d++)
    {
       int NFirst,NLast; 
       NFirst= TabRidgelet[s2d].rid_pos(s1d);
       NLast = TabRidgelet[s2d].rid_pos(s1d) + TabRidgelet[s2d].rid_size(s1d);
       for (int k=0; k < RidIma.nl(); k++)
       for (int l=NFirst; l < NLast; l++)
       {
          float Level = SigmaNoise * Nsig * (TabCurSigma[s2d])(s1d,k%BS);
 	  if (ABS(RidIma(k,l)) < Level) RidIma(k,l) = 0.;
       }
    }
    //TabRidgelet[s2d].put_block_trans(i, j, TabRidTrans[s2d], RidIma);
    TabRidgelet[s2d].recons_one_block(RidIma,ImaBlock);
}

/***************************************************************************/
/*
Use wavelet at the first scale and curvelet for the other scales.
Results are not better than curvelet only
void PyrCurvelet::mix_filtering(Ifloat &Image)
{
    int s2d,i,j,k,l;
    int Nl = Image.nl();
    int Nc = Image.nc();
    if (AllocClass == False)
    {
       cout << "Error: the curvelet class is not initialized ... " << endl;
       exit(-1);
    }
    if ((Nl != NlIma) || (Nc != NcIma))
    {
       cout << "Error: the curvelet class is not initialized for this given image size ... " << endl;
       cout << "       Image size = " << Nl << " " << Nc << endl;
       cout << "       Expected size = " << NlIma << " " << NcIma << endl;
       exit(-1);
    }
    if (SizeTabRid != NbrBand2D-1) 
    {
       cout << "Error: the curvelet class is not initialized for this number of scales ... " << endl;
       cout << "       Scale number  = " << NbrScale2D <<  endl;
       cout << "       Expected Scale number = " << SizeTabRid+1  << endl;
       exit(-1);
    }
 
    MR_Data.transform(Image);
    HALF_DECIMATED_2D_WT WT(*Ptr_SB1D);
    Ifloat * WTTabBand; 
    int NbrUndec = NbrBand2D;
    int NbrWTBand = WT.alloc(WTTabBand, Nl, Nc, NbrBand2D, NbrUndec);

    for (s2d = 0; s2d < NbrBand2D-1; s2d++)
    {
 	int Nlb = MR_Data.size_band_nl(s2d);
	int Ncb = MR_Data.size_band_nc(s2d);
        float Nsig = (s2d == 0) ? N_Sigma+1: N_Sigma;

        if (Verbose == True) 
        {
           // cout << endl;
           cout << " BAND " << s2d+1 << " Block size = " << TabBlockSize(s2d) << " NbrScale WT1D = " << TabRidgelet[s2d].NbrScale << endl;
        }
	if (s2d == 0)
	{ 
	   WT.transform(MR_Data.band(0), WTTabBand,NbrBand2D, NbrUndec);
	   for (int b= 0; b < NbrWTBand-1; b++)
	   {
	      for (k = 0; k < WTTabBand[b].nl(); k++)
              for (l = 0; l < WTTabBand[b].nc(); l++)
	         if (ABS(WTTabBand[b](k,l)) < Nsig*Noise_Ima) WTTabBand[b](k,l) = 0.;
	   }
	   WT.recons(WTTabBand, MR_Data.band(0), NbrBand2D, NbrUndec);
	}
	else
	{   
	 TabRidgelet[s2d].alloc(Nlb, Ncb, TabBlockSize(s2d));
         int BS = TabRidgelet[s2d].block_size();
         Ifloat ImaBlock(BS, BS,"block");
         Ifloat Filter(Nlb,Ncb, "Filter");
 	 cout << " Nlb X Ncb = " << TabRidgelet[s2d].rid_block_nl() << " " << TabRidgelet[s2d].rid_block_nc() << endl;
         for (i = 0; i < TabRidgelet[s2d].rid_block_nl(); i++)
         for (j = 0; j < TabRidgelet[s2d].rid_block_nc(); j++)
         {
            Ifloat RidIma;
           //if (Verbose == True) 
           //    cout << TabBlockSize(s2d) << ": Filter block " << i + 1 << " " << j+1 << endl;
           TabRidgelet[s2d].get_block_ima(i, j, MR_Data.band(s2d), ImaBlock);
           mix_threshold_one_block(s2d, ImaBlock, RidIma,Nsig,Noise_Ima);
  	   TabRidgelet[s2d].add_block_ima(i,j,Filter,ImaBlock);
           if ((i == 0) && (j == 0)) Verbose = False;
         }
	 MR_Data.band(s2d)=Filter;
	}
   }
 	
        // INFO(MR_Data.band(b), "BAND");
   MR_Data.recons(Image);
}
*/

/*************************************************************************/

static int max_scale_number (int N)
{
    int ScaleMax;
    ScaleMax  = iround((float)log((float) (N / 4. * 3.) / log(2.)));
    return ScaleMax;
}

/*************************************************************************/

void PyrCurvelet::mix_filtering(Ifloat &Image)
{
    int s2d,i,j,k,l;
    int Nl = Image.nl();
    int Nc = Image.nc();
    if (AllocClass == False)
    {
       cout << "Error: the curvelet class is not initialized ... " << endl;
       exit(-1);
    }
    if ((Nl != NlIma) || (Nc != NcIma))
    {
       cout << "Error: the curvelet class is not initialized for this given image size ... " << endl;
       cout << "       Image size = " << Nl << " " << Nc << endl;
       cout << "       Expected size = " << NlIma << " " << NcIma << endl;
       exit(-1);
    }
    if (SizeTabRid != NbrBand2D-1) 
    {
       cout << "Error: the curvelet class is not initialized for this number of scales ... " << endl;
       cout << "       Scale number  = " << NbrScale2D <<  endl;
       cout << "       Expected Scale number = " << SizeTabRid+1  << endl;
       exit(-1);
    }
    Ifloat ImSimu(Nl,Nc,"ImSimu");
    unsigned int InitRnd = 10;
    im_noise_gaussian (ImSimu, 1., InitRnd);
    fltarray TabStat[NbrBand2D];
    MR_Data.transform(ImSimu);
    for (s2d = 0; s2d < NbrBand2D-1; s2d++)
    {
 	int Nlb = MR_Data.size_band_nl(s2d);
	int Ncb = MR_Data.size_band_nc(s2d);
	HALF_DECIMATED_2D_WT WT(*Ptr_SB1D);
        Ifloat * WTTabBand; 
 	int NbrBand2DScale = max_scale_number(MIN(Nlb,Ncb));
        int NbrUndec = NbrBand2DScale;
        int NbrWTBand = WT.alloc(WTTabBand, Nlb, Ncb, NbrBand2DScale, NbrUndec);
        TabStat[s2d].alloc(NbrWTBand);
        WT.transform(MR_Data.band(s2d), WTTabBand, NbrBand2DScale, NbrUndec);
	for (int b= 0; b < NbrWTBand; b++)
	        TabStat[s2d](b) =  WTTabBand[b].sigma();
    }		    

    MR_Data.transform(Image);


    for (s2d = 0; s2d < NbrBand2D-1; s2d++)
    {
 	int Nlb = MR_Data.size_band_nl(s2d);
	int Ncb = MR_Data.size_band_nc(s2d);
        float Nsig = (s2d == 0) ? N_Sigma+1: N_Sigma;

    
        if (Verbose == True) 
        {
           // cout << endl;
           cout << " Scale " << s2d+1 << " Block size = " << TabBlockSize(s2d) << " NbrScale WT1D = " << TabRidgelet[s2d].NbrScale << endl;
        }
 
	 TabRidgelet[s2d].alloc(Nlb, Ncb, TabBlockSize(s2d));
         int BS = TabRidgelet[s2d].block_size();
         Ifloat ImaBlock(BS, BS,"block");
         HALF_DECIMATED_2D_WT WT(*Ptr_SB1D);
         Ifloat * WTTabBand; 
	 
 	 int NbrBand2DScale = max_scale_number(MIN(Nlb,Ncb));
         int NbrUndec = NbrBand2DScale;
         int NbrWTBand = WT.alloc(WTTabBand, Nlb, Ncb, NbrBand2DScale, NbrUndec);

         Ifloat WTFilter(Nlb,Ncb, "Filter");
         WT.transform(MR_Data.band(s2d), WTTabBand, NbrBand2DScale, NbrUndec);
	 for (int b= 0; b < NbrWTBand; b++)
	 {
	    // cout << "     Band " << b+1 << " Norm = " << TabStat[s2d](b) << 
	    float Level = Noise_Ima*Nsig*TabStat[s2d](b);
            for (k = 0; k < WTTabBand[b].nl(); k++)
            for (l = 0; l < WTTabBand[b].nc(); l++)
	            if (ABS(WTTabBand[b](k,l)) < Level) WTTabBand[b](k,l) = 0.;
	 }
         WT.recons(WTTabBand, WTFilter, NbrBand2DScale, NbrUndec);
	 // MR_Data.band(s2d)=WTFilter;
	       	
         Ifloat Filter(Nlb,Ncb, "Filter");
 	 cout << " NbrBand2DScale = " << NbrBand2DScale << " Nlb X Ncb = " << TabRidgelet[s2d].rid_block_nl() << " " << TabRidgelet[s2d].rid_block_nc() << endl;
         for (i = 0; i < TabRidgelet[s2d].rid_block_nl(); i++)
         for (j = 0; j < TabRidgelet[s2d].rid_block_nc(); j++)
         {
            Ifloat RidIma;
           //if (Verbose == True) 
           //    cout << TabBlockSize(s2d) << ": Filter block " << i + 1 << " " << j+1 << endl;
           TabRidgelet[s2d].get_block_ima(i, j, MR_Data.band(s2d), ImaBlock);
	   if ( ((i % 2 == 0) && (j% 2 == 0)) ||
	        ((i % 2 == 1) && (j% 2 == 1))) 
	   {
	      //cout << "RIDRUN" << i+1 << " " << j+1 << endl;
              mix_threshold_one_block(s2d, ImaBlock, RidIma,Nsig,Noise_Ima);
	      //cout << "RIDEND" << endl;
	   }
	   else
	   {
	      //cout << "WTRUN" << i+1 << " " << j+1  << endl;
	      TabRidgelet[s2d].get_block_ima(i, j, WTFilter, ImaBlock);
// 	      WT.transform(ImaBlock, WTTabBand, NbrBand2DScale, NbrUndec);
// 	      for (int b= 0; b < NbrWTBand; b++)
// 	      {
// 	          for (k = 0; k < WTTabBand[b].nl(); k++)
//                   for (l = 0; l < WTTabBand[b].nc(); l++)
// 	            if (ABS(WTTabBand[b](k,l)) < Nsig*Noise_Ima) WTTabBand[b](k,l) = 0.;
// 	      }
// 	      WT.recons(WTTabBand, ImaBlock, NbrBand2DScale, NbrUndec);
	      // cout << "WTEND" << endl;
	   } 
  	   TabRidgelet[s2d].add_block_ima(i,j,Filter,ImaBlock);
           // if ((i == 0) && (j == 0)) Verbose = False;
	   
         }
	 MR_Data.band(s2d)=Filter;
   }
 	
        // INFO(MR_Data.band(b), "BAND");
   MR_Data.recons(Image);
}
/****************************************************************************/

int main(int argc, char *argv[])
{
    int k;
    fitsstruct Header;
    char Cmd[512];
    Cmd[0] = '\0';
    for (k =0; k < argc; k++) sprintf(Cmd, "%s %s", Cmd, argv[k]);
     
    // Get command line arguments, open input file(s) if necessary
    lm_check(LIC_MR4);
    filtinit(argc, argv);

    if (Verbose == True)
    { 
        cout << endl << endl << "PARAMETERS: " << endl << endl;
        cout << "File Name in = " << Name_Imag_In << endl;
        cout << "File Name Out = " << Name_Imag_Out << endl;  
        cout << "Ridgelet transform = " <<  StringRidTransform(RidTrans) << endl;  
        cout << "Transform WT2D = " << StringTransform(Transform) << endl;
        if (BlockSize > 0)  cout << "BlockSize = " << BlockSize <<endl;
        else  cout << "No partitionning  " << endl;
        cout << "Nbr Scale = " <<   NbrScale2D << endl;  
        cout << "NSigma = " << N_Sigma  << endl;  
	cout << "Type block = " << StringBlockType(TypeBlock) << endl;
	if (BlockOverlap == True) cout << "Block overlapping " << endl;
	else cout << "No block overlapping " << endl;
    }

    Ifloat Data;
    io_read_ima_float(Name_Imag_In, Data, &Header);
    Header.origin = Cmd;
 
    if (BlockSize > Data.nc())
    {
        cout << "Warning: Block size must lower than the image size ... " << endl;
        cout << "         Block size is set to image size " << endl;
        BlockSize = 0;
    }
#ifdef DO_SIMU
    if (Simu == True) make_ima(Data);
#endif

    if (Noise_Ima < FLOAT_EPSILON)
    {
       Noise_Ima = detect_noise_from_med (Data);
       if (Verbose == True) cout << "Sigma Noise = " << Noise_Ima << endl;
    }
    Ifloat NoiseData;
    Ifloat Filter(Data.nl(), Data.nc(), "filter");

    if (NoiseIma == True)
    {
        io_read_ima_float(NoiseFile, NoiseData);
        if ((NoiseData.nl() != Data.nl()) || (NoiseData.nc() != Data.nc()))
       {
         cout << "Error: the noise map must have the same size" << endl;
         cout << "       as the input data ... " << endl;
	 exit(-1);
       }             
    }
        
    if (FastCur == False)
    {
    //cout << "MK SelectFilter " << endl;

    FilterAnaSynt SelectFilter(F_MALLAT_7_9);
    SubBandFilter SB1D(SelectFilter, NORM_L2);
    SB1D.Border = I_MIRROR;
    //cout << "MK Cur(SB1D) " << endl;
    PyrCurvelet Cur(SB1D);
    //cout << "DONE Cur(SB1D) " << endl;
    Cur.RidTrans = RidTrans;
    Cur.NbrScale2D = NbrScale2D;
    Cur.ColTrans = ColTrans;
    Cur.TypeBlock = TypeBlock;
    if ((RidTrans != RID_PYR_FFT) && (RidTrans !=  RID_PYR_DIRECT) && (RidTrans != RID_OWT)) Cur.OddBlockSize=False;
    Cur.WPTrans = False;
    Cur.TransformWT2D = Transform;
    if (NbrScale1D >= 0)
    {  
        Cur.NbrScale1D = NbrScale1D;
        Cur.GetAutoNbScale = False;
        if (NbrScale1D < 2) Cur.BlockOverlap = False;
        else Cur.BlockOverlap = BlockOverlap;
    }
    else Cur.BlockOverlap = BlockOverlap;
    Cur.Verbose = Verbose;     
    Cur.alloc(Data.nl(), Data.nc(), BlockSize);

    if (NoiseIma == True) Cur.get_norm_coeff(NoiseData, N_Sigma);
    else Cur.get_norm_coeff(N_Sigma);
    
    // fltarray TabStat;
    // Cur.get_stat(TabStat);
    
    if (Verbose == True) cout << "Filtering ... " << endl;
    if (MIX == True)
    {
       // Cur.mix_filtering(Data);
       Cur.transform(Data);
       Cur.wiener_filter(Noise_Ima);
       Cur.recons(Filter);
       if (PositivSol == True) threshold(Data);
       io_write_ima_float(Name_Imag_Out, Data, &Header);
    }
    else
    {
    Cur.transform(Data);
    Cur.threshold(Noise_Ima, N_Sigma);
    Cur.recons(Filter);
     }
    }
    
    else // FastCur
    {
      FCUR  Cur;
      Cur.Verbose = Verbose;   
      Cur.alloc_from_fine(NbrScale2D, Data.nl(), Data.nc(), BlockSize, ExtendWT,IsotropWT, True);
      if (NoiseIma == True) Cur.get_norm_coeff(NoiseData, N_Sigma);
      else Cur.get_norm_coeff(N_Sigma);
     
      if (Verbose == True) cout << "Filtering ... " << Noise_Ima <<  N_Sigma << endl;
      Cur.cur_trans(Data);
      for (int s=0; s < Cur.nbr_scale()-1; s++)
      for (int b=0; b < Cur.nbr_band(s); b++) 
      {
         Ifloat BandR;
 	 Cur.get_band(s, b, BandR);
 	 cout << "Scale " << s+1 << " Band " << b+1 << "Sig = " << BandR.sigma() << " sqrt(Ener) = " <<  sqrt(BandR.energy()) << endl;
      }
      if (RecWT == True)
      {
         Cur.cur_recons(Filter);
         io_write_ima_float("xx_rec.fits", Filter);
         Cur.cur_trans(Data);
      }
      if (MIX == False) Cur.threshold(Noise_Ima, N_Sigma);
      else Cur.wiener_filter(Noise_Ima);
      Cur.cur_recons(Filter);
      }
     if (PositivSol == True) threshold(Filter);
      io_write_ima_float(Name_Imag_Out, Filter, &Header);
    exit(0);
}


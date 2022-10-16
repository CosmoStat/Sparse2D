/******************************************************************************
**                   Copyright (C) 1994 by CEA
*******************************************************************************
*******************************************************************************
**
**    DESCRIPTION  compute the histograms for the few event noise model 3D
**    ----------- 
**                 
******************************************************************************/

#include "mr3d_phisto.h"


/*********************************************************************/
//			Class User_Param
/*********************************************************************/
void User_Param::put_default_val(){
	Rythm=0;
	NAutoConv=30;
	Verbose = False;
	WriteOpt=False;	
	WriteRes=False;
	DefaultOpt=False;
}

void User_Param::usage(char *argv[]){
    fprintf(OUTMAN, "Usage: %s options\n\n", argv[0]);
    fprintf(OUTMAN, "   where options =  \n");
    
    fprintf(OUTMAN, "         [-n number_of_convolution]\n");
    fprintf(OUTMAN, "              Total number of calculated histograms. Default is 30.\n");
    
    fprintf(OUTMAN, "         [-r rythm]\n");
    fprintf(OUTMAN, "              parameter relative to the rythm of autoconvolution ( default=0 )\n");
    fprintf(OUTMAN, "              the program computes the 2^n autoconvolutions of the histogram of\n");
    fprintf(OUTMAN, "              the wavelet. It computes also 2^r convolutions between\n");
    fprintf(OUTMAN, "              2^n and 2^(n+1).\n");

    fprintf(OUTMAN, "         [-d]\n");
    fprintf(OUTMAN, "              Use all default options.\n");
        
//    fprintf(OUTMAN, "         [-w ]\n");
//    fprintf(OUTMAN, "              write each histogram an independent file\n");
    
    verbose_usage();
    
    manline();

    exit(-1);
}

/*********************************************************************/

/* GET COMMAND LINE ARGUMENTS */
User_Param::User_Param(int argc, char *argv[]){
    int c;
    put_default_val();//put default values for each parametres 
    /* get options */
    while ((c = GetOpt(argc,argv,"dvwr:n:")) != -1) 
    {
	switch (c) 
        {
	   case 'd': DefaultOpt = True; break;
	   case 'v': Verbose = True; WriteOpt=True; break;
	   case 'w': WriteRes=True; break;
	   case 'r':
		/* -r <rythm> */
		if ((sscanf(OptArg,"%d",&Rythm) != 1) || (Rythm < 0)){
		    fprintf(OUTMAN, "bad parameter: %s\n", OptArg);
		    exit(-1);
		}
		break;	
	   case 'n':
		/* -c <NAutoConv> */
		if ((sscanf(OptArg,"%d",&NAutoConv) != 1) || (NAutoConv < 0)){
		    fprintf(OUTMAN, "bad number of convolution: %s\n", OptArg);
		    exit(-1);
		}
		break;	
	   
            case 'h': usage(argv); break;
	    default: usage(argv); break;
  	}
    }
    if (DefaultOpt == False) usage(argv);
    
    /* make sure there are not too many parameters */
    if (OptInd < argc)
    {
	fprintf(OUTMAN, "Too many parameters: %s ...\n", argv[OptInd]);
	usage(argv);
    }

}

/*********************************************************************/
void User_Param::WriteOptProc(){
	cout << "####################################################"
		<< endl;
	cout << "STATE OF OPTIONS FOR mr3D_histo: "<< endl;
	cout << "rythm:				"<<Rythm<<endl;
	cout << "number_of_convolution:		"<<NAutoConv<<endl;
	if(WriteRes)cout << "Writing histograms in different files"<<endl;
	if(Verbose) cout << "verbose mode"<<endl;
	cout << "####################################################"
		<< endl;
	cout << endl;
}



 /********************************************************************/   
 /***************************P P**************************************/    
 /********************************************************************/    
int main(int argc, char *argv[])
{
    double EpsVal = 1e-10; // The threshold computation will be 
                           // correct for espsilon values lower than 
			   // 1e-10
    lm_check(LIC_MR3);
    User_Param Param(argc, argv);// Get command line arguments

    //creating FewEventPoisson class
    FewEventPoisson FEP(False,Param.NAutoConv,EpsVal);
      
    if(Param.WriteOpt) Param.WriteOptProc();

    if(Param.Verbose){
    	cout <<"Initalisation of all histogrammes" << endl;
    	cout <<"Computing Wavelet Distribution"<<endl;
    }
    FEP.compute_distribution(Param.WriteRes, True, Param.Verbose, Param.Rythm);

}

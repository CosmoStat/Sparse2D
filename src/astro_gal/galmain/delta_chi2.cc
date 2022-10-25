/******************************************************************************
**                   Copyright (C) 2011 by CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 1.0
**
**    Author: Antoine Labatie
**
**    Date:  25/07/2011
**    
**    File:  delta_chi2.cc
**
*******************************************************************************
**
**    DESCRIPTION  Delta chi2 method 
**    ----------- 
**                 
******************************************************************************/
 

#include "Array.h"
#include "IM_IO.h"
#include <omp.h>


int h=0;
Bool Verbose=False;
Bool VarCov=False;
Bool Data=False;
char Name_Xi_In[256]; //Name of input correlation file for analyzing data

// (TO BE SET IN DELTA_CHI2.PARAM FILE)
double o_min1,o_max1,delta_o,o_min2,o_max2;				//Omega_m h^2
double a_min1,a_max1,delta_a,a_min2,a_max2;				//alpha
double B_min1,B_max1,delta_B,B_min2,B_max2,B_model;		//B

char Name_o_table[256]; char Name_a_table[256]; char Name_B_table[256]; 

char Name_Inverse_Cov[256];     //name for inverse convariance matrix
char Name_Sqrt_Cov[256];        //name for square root convariance matrix
char Name_Sqrt_Cov_all[256];        //name for square root convariance matrix
char Name_Model_BAO_all[256];   //name for BAO model correlations
char Name_Model_noBAO_all[256]; //name for noBAO model correlations
double n_simu;

char Name_Histo_Out_Prefix[256];		/* output file prefix name */


//maximum number of procs used for the loops
int Nproc_max=40;

/****************************************************************************/

static void usage(char *argv[])
{
    fprintf(OUTMAN, "Usage: %s options \n\n", argv[0]);
    fprintf(OUTMAN, "   where options =  \n");
  
    fprintf(OUTMAN, "         [-h h_value]\n");
    fprintf(OUTMAN, "             Hypothesis sampled (h=0 for noBAO, h=1 for BAO).\n"); 
    fprintf(OUTMAN, "             Default is 0."); 
    manline();

    fprintf(OUTMAN, "         [-c ]\n");
    fprintf(OUTMAN, "             Use non constant covariance matrix.\n"); 
    manline();
	
	fprintf(OUTMAN, "         [-d Name_Xi_In]\n");
    fprintf(OUTMAN, "             Apply method to data in Name_Xi_In.\n");
	manline();

	
    vm_usage();
    manline();    
    verbose_usage();
    manline();         
    manline();
    exit(-1);
}
 
/*********************************************************************/

void gauss(double sigma, double &x, double &y)    /* Recipes */
{
	//Create a 2D gaussian independent in x,y  => f(x,y)=1/(2*pi*sigma) e^-(x^2+y^2)/(2*sigma^2)
    double v1,v2,r,fac;
	
    do {
        v1=2.0*drand48()-1; 
        v2=2.0*drand48()-1; 
        r=v1*v1+v2*v2;
	}
    while(r>=1.0);
    fac=sqrt(-2*sigma*log((double) r)/r);
    x=v1*fac;
    y=v2*fac;
}


/*********************************************************************/

/* GET PARAMETERS */

void get_param()
{
	char Name_Param_File[256];
	sprintf(Name_Param_File, "../param/delta_chi2.param");	
	FILE *File=fopen(Name_Param_File,"r");
	
    if (File == NULL)
    {
		cerr << "Error: cannot open file "  <<  Name_Param_File << endl;
		exit(-1);
    }
	char Temp[256];
	int ret;
	ret=fscanf(File, "%s\t%lf\t%lf\t%lf\n", Temp, &o_min1, &o_max1, &delta_o);	//min, max, step in Omega_m h^2 for hypotheses
	ret=fscanf(File, "%s\t%lf\t%lf\t%lf\n", Temp, &a_min1, &a_max1, &delta_a);	//min, max, step in alpha for hypotheses
	ret=fscanf(File, "%s\t%lf\t%lf\t%lf\n", Temp, &B_min1, &B_max1, &delta_B);	//min, max, step in B for hypotheses
	ret=fscanf(File, "%s\t%lf\t%lf\n", Temp, &o_min2, &o_max2);					//min, max, in Omega_m h^2 for fitting
	ret=fscanf(File, "%s\t%lf\t%lf\n", Temp, &a_min2, &a_max2);					//min, max, in alpha for fitting
	ret=fscanf(File, "%s\t%lf\t%lf\n", Temp, &B_min2, &B_max2);					//min, max, in B for fitting
	ret=fscanf(File, "%s\t%lf\n\n", Temp, &B_model);								//B model
	
	ret=fscanf(File, "%s\t%s\n", Temp, Name_a_table);	
	ret=fscanf(File, "%s\t%s\n", Temp, Name_o_table);	
	ret=fscanf(File, "%s\t%s\n\n", Temp, Name_B_table);
	
	ret=fscanf(File, "%s\t%s\n", Temp, Name_Model_BAO_all);	
	ret=fscanf(File, "%s\t%s\n", Temp, Name_Model_noBAO_all);	
	ret=fscanf(File, "%s\t%s\n", Temp, Name_Sqrt_Cov_all);
	ret=fscanf(File, "%s\t%s\n", Temp, Name_Inverse_Cov);	
	ret=fscanf(File, "%s\t%s\n", Temp, Name_Sqrt_Cov);
	ret=fscanf(File, "%s\t%lf\n\n", Temp, &n_simu);

	ret=fscanf(File, "%s\t%s\n", Temp, Name_Histo_Out_Prefix);
	
	fclose(File);
}

/*********************************************************************/

/* GET COMMAND LINE ARGUMENTS */
static void filtinit(int argc, char *argv[])
{
	int i=1;
	while(i<argc) 
	{
		if(argv[i][0] != '-' ) break;
		
		switch (argv[i][1]) 
        {     
			case 'h': 
				if(i+1==argc) 
				{
					fprintf(OUTMAN, "Error: no argument for -h.\n"); exit(-1);
				}
				h=atoi(argv[++i]);
				if(strcmp(argv[i],"0")!=0 && strcmp(argv[i],"1")!=0)
					{
						fprintf(OUTMAN, "Error: bad number value for h: %s\n",argv[i] );
						exit(-1);
					}
					break;
			case 'd': 
				if(i+1==argc) 
				{
					fprintf(OUTMAN, "Error: no argument for -d.\n"); exit(-1);
				}
				strcpy(Name_Xi_In, argv[++i]);
				Data=True;
				h=-1;
				break;
			case 'c': VarCov = True;break;
			case 'v': Verbose = True;break;
			case '?': usage(argv); break;
			default: usage(argv); break;
			}
		i++;
	}

       
	/* make sure there are not too many parameters */
	if (i < argc)
        {
		fprintf(OUTMAN, "Too many parameters: %s ...\n", argv[i]);
		exit(-1);
	}


}


/***************/


int main(int argc, char *argv[])
{
     /* Get command line arguments, open input file(s) if necessary */
    filtinit(argc, argv);

	/* Get parameters */	
	get_param();
	
	if(o_min1 > o_max1 || a_min1>a_max1 || B_min1 >B_max1 || o_min2 > o_max2 || a_min2>a_max2 || B_min2 >B_max2 ) 
	{
		fprintf(OUTMAN,"minimum values of parameters must be less than maximum values.\n");
		exit(-1);
	}
	

	//Names
	char Name_Histo_Out[256];       //name for writing out histrogram
	strcpy(Name_Histo_Out, Name_Histo_Out_Prefix); 
	
	if(h==0 && VarCov==False)	strcat(Name_Histo_Out, "Dchi2_h0.fits");
	if(h==0 && VarCov==True)	strcat(Name_Histo_Out, "Dchi2_varcov_h0.fits");
	if(h==1 && VarCov==False)	strcat(Name_Histo_Out, "Dchi2_h1.fits");
	if(h==1 && VarCov==True)	strcat(Name_Histo_Out, "Dchi2_varcov_h1.fits");
	
	
	//Read arrays
	fltarray iC,sC,iC_all,sC_all,model_BAO_all, model_noBAO_all;
	fits_read_fltarr(Name_Sqrt_Cov_all, sC_all);
	fits_read_fltarr(Name_Inverse_Cov, iC);
	fits_read_fltarr(Name_Sqrt_Cov, sC);
	fits_read_fltarr(Name_Model_BAO_all, model_BAO_all);	
	fits_read_fltarr(Name_Model_noBAO_all, model_noBAO_all);

	
	int nr=model_BAO_all.nz(); 
	if(nr != model_noBAO_all.nz()) 
	{
		fprintf(OUTMAN,"Incompatible dimensions between models.\n");
		exit(-1);
	}
	if(nr != iC.nx() || nr != iC.ny()) 
	{
		fprintf(OUTMAN,"Incompatible dimensions between models and covariance matrix.\n");
		exit(-1);
	}

	//read parameter table
	fltarray o_table,a_table,B_table;
	fits_read_fltarr(Name_o_table, o_table); fits_read_fltarr(Name_a_table, a_table); fits_read_fltarr(Name_B_table, B_table);
	
	int no_table=o_table.nx(); 	int na_table=a_table.nx(); 	int nB_table=B_table.nx();
	if(no_table != model_BAO_all.nx() || no_table != model_noBAO_all.nx() ||  \
	   na_table != model_BAO_all.ny() || na_table != model_noBAO_all.ny())
	{
		cout <<"Incompatible dimensions between models and parameter tables" << endl;
		exit(-1);
	}
	double delta_o_table=o_table(1)-o_table(0); double delta_a_table=a_table(1)-a_table(0); double delta_B_table=B_table(1)-B_table(0);
	
	
	//Find indices for hypotheses model correlations according to requested limits
	int oind_min1,oind_max1,aind_min1,aind_max1,Bind_min1,Bind_max1;
	int delta_oind,delta_aind,delta_Bind;
	delta_oind=round(delta_o/delta_o_table); delta_aind=round(delta_a/delta_a_table); delta_Bind=round(delta_B/delta_B_table);

 	//readjust delta parameters
	delta_o=delta_oind*delta_o_table; delta_a=delta_aind*delta_a_table; delta_B=delta_Bind*delta_B_table;
	
	oind_min1=round((o_min1-o_table.min())/(delta_o_table)); oind_max1=round((o_max1-o_table.min())/delta_o_table);
	aind_min1=round((a_min1-a_table.min())/(delta_a_table)); aind_max1=round((a_max1-a_table.min())/delta_a_table);
	Bind_min1=round((B_min1-B_table.min())/(delta_B_table)); Bind_max1=round((B_max1-B_table.min())/delta_B_table);
	if(oind_min1<0 || oind_min1 >=no_table || oind_max1<0 || oind_max1 >=no_table || \
	   aind_min1<0 || aind_min1 >=na_table || aind_max1<0 || aind_max1 >=na_table || \
	   Bind_min1<0 || Bind_min1 >=nB_table || Bind_max1<0 || Bind_max1 >=nB_table)
	{
		cout << "Hypotheses range out of limits" << endl; exit(-1);
	}
	if(delta_oind<1 || delta_aind<1 || delta_Bind<1)  
	{
		cout << "Too small binning of the parameters (smaller than the model grid)" << endl; exit(-1);
	}
	oind_min1=round(double(oind_min1)/double(delta_oind)); oind_max1=round(double(oind_max1)/double(delta_oind));
	aind_min1=round(double(aind_min1)/double(delta_aind)); aind_max1=round(double(aind_max1)/double(delta_aind));
	Bind_min1=round(double(Bind_min1)/double(delta_Bind)); Bind_max1=round(double(Bind_max1)/double(delta_Bind));
	int nind_o1=oind_max1-oind_min1+1; 	int nind_a1=aind_max1-aind_min1+1; 	int nind_B1=Bind_max1-Bind_min1+1;
	
	//Readjust parameter limits
	o_min1=o_table(oind_min1*delta_oind); o_max1=o_table(oind_max1*delta_oind);
	a_min1=a_table(aind_min1*delta_aind); a_max1=a_table(aind_max1*delta_aind);
	B_min1=B_table(Bind_min1*delta_Bind); B_max1=B_table(Bind_max1*delta_Bind);

	
	
	//Find indices for model correlation functions to fit according to requested limits
	int oind_min2,oind_max2,aind_min2,aind_max2;

	oind_min2=round((o_min2-o_table.min())/(delta_o_table)); oind_max2=round((o_max2-o_table.min())/delta_o_table);
	aind_min2=round((a_min2-a_table.min())/(delta_a_table)); aind_max2=round((a_max2-a_table.min())/delta_a_table);
	if(oind_min2<0 || oind_min2 >=no_table || oind_max2<0 || oind_max2 >=no_table || \
	   aind_min2<0 || aind_min2 >=na_table || aind_max2<0 || aind_max2 >=na_table )
	{
		cout << "Fitting range out of limits" << endl; exit(-1);
	}
	oind_min2=round(double(oind_min2)/double(delta_oind)); oind_max2=round(double(oind_max2)/double(delta_oind));
	aind_min2=round(double(aind_min2)/double(delta_aind)); aind_max2=round(double(aind_max2)/double(delta_aind));
	
	//Readjust parameter limitis
	o_min2=o_table(oind_min2*delta_oind); o_max2=o_table(oind_max2*delta_oind);
	a_min2=a_table(aind_min2*delta_aind); a_max2=a_table(aind_max2*delta_aind);

	
	
	//Print parameters
    if (Verbose == True)
    {
		cout << endl << endl << "DELTA CHI2 (difference of best-fits)" << endl << endl;
		if(h==0) cout << "-Sample H0 (noBAO) hypothesis " << endl;
		if(h==1) cout << "-Sample H1 (BAO) hypothesis " << endl;
		if(Data==True) cout << "-Delta chi2 for single xi data in " << Name_Xi_In << endl;
		if(VarCov==True) cout << "-non constant covariance matrix ";		
		
		cout << endl << endl << "RANGE FOR HYPOTHESES: " << endl;
		cout << "Omega_m h^2 grid used (min,max,delta) = ("<< o_min1 << " , " << o_max1 << " , " << delta_o << ")"<< endl ;
		cout << "alpha grid used (min,max,delta) = ("<< a_min1 << " , " << a_max1 << " , " << delta_a << ")"<< endl;
		cout << "B (min,max,delta) used (min,max,delta) = ("<< B_min1 << " , " << B_max1 << " , " << delta_B << ")"<< endl << endl;
		
		cout << endl << endl << "FITTING RANGE: " << endl;
		cout << "Omega_m h^2 grid used (min,max,delta) = ("<< o_min2 << " , " << o_max2 << " , " << delta_o << ")"<< endl ;
		cout << "alpha grid used (min,max,delta) = ("<< a_min2 << " , " << a_max2 << " , " << delta_a << ")"<< endl;
		cout << "B grid used (min,max) = ("<< B_min2 << " , " << B_max2 << ")"<< endl << endl;
		
		cout << "Number of simulations for each point = " << n_simu << endl << endl;
		
		cout << "Write histo in file " << Name_Histo_Out << endl << endl << endl;
    }
	
	//Delta chi2 for single xi data
	if(Data==True)
	{
		fltarray model_BAO(nr);
		fltarray model_noBAO(nr);
		fltarray xi(nr);
		fits_read_fltarr(Name_Xi_In, xi);
		
		dblarray chi2_noBAO(oind_max2-oind_min2+1,aind_max2-aind_min2+1) ;
		dblarray chi2_BAO(oind_max2-oind_min2+1,aind_max2-aind_min2+1) ;	
		double B;
		
		for(int oind2=oind_min2; oind2<=oind_max2; oind2++)
		{
			for(int aind2=aind_min2; aind2<=aind_max2; aind2++)
			{
				for(int i=0;i<nr;i++) model_BAO(i)=model_BAO_all(oind2*delta_oind,aind2*delta_aind,i);
				for(int i=0;i<nr;i++) model_noBAO(i)=model_noBAO_all(oind2*delta_oind,aind2*delta_aind,i);
				
				B=(xi*mult(iC,model_BAO)).total()/(model_BAO*mult(iC,model_BAO)).total() ; //bias giving best-fit \chi^2
				if(B>B_max2/B_model) B=B_max2/B_model;
				if(B<B_min2/B_model) B=B_min2/B_model;
				chi2_BAO(oind2-oind_min2,aind2-aind_min2)=((xi-B*model_BAO)*mult(iC,xi-B*model_BAO)).total() ;

				B=(xi*mult(iC,model_noBAO)).total()/(model_noBAO*mult(iC,model_noBAO)).total()  ; //bias giving best-fit \chi^2
				if(B>B_max2/B_model) B=B_max2/B_model;
				if(B<B_min2/B_model) B=B_min2/B_model;
				chi2_noBAO(oind2-oind_min2,aind2-aind_min2)=((xi-B*model_noBAO)*mult(iC,xi-B*model_noBAO)).total();

			}
		}

		fltarray Dchi2_data(1);
		Dchi2_data(0)=chi2_noBAO.min()-chi2_BAO.min() ;
		char Name_Data_Out[256];
		sprintf(Name_Data_Out, "../output_files/delta_chi2/Dchi2_data.fits");
		fits_write_fltarr(Name_Data_Out,Dchi2_data);
		cout << "Delta chi2: " << Dchi2_data(0) << endl << endl;
		exit(0);
	}
	
	//histogram under H0 or H1
	fltarray Dchi2_histo(nind_o1*nind_a1*nind_B1*n_simu);
	long int count=0;

	#ifdef _OPENMP
		int Nproc=omp_get_num_procs();
		if(Nproc>Nproc_max) Nproc=Nproc_max;
		printf("\n %u processors used for the loop \n\n",Nproc);
	#endif	
	
	
	//Main loop for simu histogram
    #pragma omp parallel default(shared) num_threads(Nproc)
    {
		int index_z;
		
		
		//Create variables for the loop
		long int histo_ind=0;
		fltarray model0(nr,1);
		fltarray model_BAO(nr,1);
		fltarray model_noBAO(nr,1);

		fltarray g(nr); //gaussian standard variable
		fltarray xi(nr,1);

		dblarray chi2_noBAO(oind_max2-oind_min2+1,aind_max2-aind_min2+1) ;
		dblarray chi2_BAO(oind_max2-oind_min2+1,aind_max2-aind_min2+1) ;

		double B,B1;
	
		//Start loop
		#pragma omp for schedule(dynamic)
		for(int j=0;j<nind_o1*nind_a1*nind_B1;j++)
		{
			int oind1=j/(nind_a1*nind_B1)+oind_min1;
			int temp=j-(oind1-oind_min1)*(nind_a1*nind_B1);
			int aind1=temp/nind_B1+aind_min1;
			int Bind1=temp-(aind1-aind_min1)*nind_B1+Bind_min1;
			
			B1=B_table(Bind1*delta_Bind)/B_model;
			index_z=na_table*oind1*delta_oind+aind1*delta_aind;
			
			#pragma omp critical
			{
				count++;
				if(Verbose==True) cout << count << '/' << nind_o1*nind_a1*nind_B1 << endl ;
			}
			
			//model for hypothesis H0 or H1
			if(h==0)
				for(int i=0;i<nr;i++) model0(i,0)=model_noBAO_all(oind1*delta_oind,aind1*delta_aind,i);
			if(h==1)
				for(int i=0;i<nr;i++) model0(i,0)=model_BAO_all(oind1*delta_oind,aind1*delta_aind,i); 
			for(int sind=0; sind<n_simu; sind++)
			{
				g.init();
				int k=0;
				double temp_x,temp_y;
				while(k<=nr-1)
				{
					gauss(1.0,temp_x,temp_y);
					g(k)=temp_x; k++;
					if(k<=nr-1) { g(k)=temp_y; k++;}
				}
				
				if(VarCov==False) xi=mult(sC,g);   //Gaussian of mean 0 and covariance C
				if(VarCov==True)  xi=mult(sC_all,g,index_z);
					
				if(VarCov==False) xi+=B1*model0;      //Constant covariance matrix
				if(VarCov==True)  xi=B1*(xi+model0);  //Varying covariance matrix					
				
				
				for(int oind2=oind_min2; oind2<=oind_max2; oind2++)
				{
					for(int aind2=aind_min2; aind2<=aind_max2; aind2++)
					{
						for(int i=0;i<nr;i++) model_BAO(i,0)=model_BAO_all(oind2*delta_oind,aind2*delta_aind,i);
						for(int i=0;i<nr;i++) model_noBAO(i,0)=model_noBAO_all(oind2*delta_oind,aind2*delta_aind,i);

						B=(xi*mult(iC,model_BAO)).total()/(model_BAO*mult(iC,model_BAO)).total(); //bias giving best-fit \chi^2 
						if(B>B_max2/B_model) B=B_max2/B_model;
						if(B<B_min2/B_model) B=B_min2/B_model;
						chi2_BAO(oind2-oind_min2,aind2-aind_min2)=((xi-B*model_BAO)*mult(iC,xi-B*model_BAO)).total() ;

						B=(xi*mult(iC,model_noBAO)).total()/(model_noBAO*mult(iC,model_noBAO)).total(); //bias giving best-fit \chi^2
						if(B>B_max2/B_model) B=B_max2/B_model;
						if(B<B_min2/B_model) B=B_min2/B_model;
						chi2_noBAO(oind2-oind_min2,aind2-aind_min2)=((xi-B*model_noBAO)*mult(iC,xi-B*model_noBAO)).total();
					}
				}
				histo_ind=sind+n_simu*(Bind1-Bind_min1)+n_simu*nind_B1*(aind1-aind_min1)+n_simu*nind_B1*nind_a1*(oind1-oind_min1);
				#pragma omp critical
				{
					Dchi2_histo(histo_ind)=chi2_noBAO.min()-chi2_BAO.min(); 
					//cout << chi2_noBAO.min() << " " << chi2_BAO.min() << "  ";
					//if(chi2_noBAO.min()>chi2_BAO.min()) cout << sqrt(chi2_noBAO.min()-chi2_BAO.min())<< endl;
					//else cout << endl;
				}
			}
		}
	}
						

    // Write the histogram
    fits_write_fltarr(Name_Histo_Out, Dchi2_histo);
    exit(0);
}

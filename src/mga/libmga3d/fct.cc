
#include "fct.h"

extern bool Verbose;

static const double PI4 = M_PI/4.0;
static const double PI2 = M_PI/2.0;

static inline int nearest_odd(int i)
{
	return 2*int(floor(i/2))+1;
}

static inline bool isnan(std::complex<float> a)
{
	return (a!=a);
}
static inline double prob_noise(float Val)
{
	return (abs(Val) < FLOAT_EPSILON) ? 1.:erfc (abs(Val));
}

// ********************************************************************

FCurvelet3D::FCurvelet3D() : MEYER_WT3D()
{
	RealBand=False;
	tabsigma=false;
	extract_type=DEF_TYPE_EXTRACTION;
	
// Filtering properties
	no_norm=false;
	no_fine=false;
	no_coarse=false;
	lowmem=false;
	eval_stat=false;
	use_min=false; use_max=false;
	min_value=0; max_value=0;
	TabBand_norm = NULL;
	th_coarse=false;
	
	pointer=false;
	address=NULL;
}

FCurvelet3D::~FCurvelet3D()
{
	delete [] TabSizeNx;
	delete [] TabSizeNy;
	delete [] TabSizeNz;
}


// ********************************************************************

double FCurvelet3D::fw(double x)
{
	double r,l;
	if(x<=0)
		r = 1;
	else if(x>=1)
		r = 0;
	else
	{
		l = exp(1-1/(1-exp(1-1/(1-x))));
		l *= l;
		r = exp(1-1/(1-exp(1-1/x)));
		r *= r;
		double norm = l+r;
		r /= norm;
//		r=(1.0-x);
//		r=1.0*(x<0.45);//+0.5*(x==0.5);
	}
	return r;
}

// ********************************************************************

void FCurvelet3D::scale_into_single_wedge(cfarray& Scale, cfarray & FBand, int b, int nd, double DNX, double DNY, double DNZ, bool backward)
{
	bool LocVerbose = false && Verbose;
	if(Verbose)
	{
		if(backward)cerr<<"Wedge_into_single_scale(.,.,"<<nd<<","<<DNX<<","<<DNY<<","<<DNZ<<")"<<endl;
		else cerr<<"Scale_into_single_wedge(.,.,"<<nd<<","<<DNX<<","<<DNY<<","<<DNZ<<")"<<endl;
	}
	
	int normalize_wedge = (backward && extract_type==EXTRACT_B) + (!backward && extract_type==EXTRACT_F) + 2*(extract_type==EXTRACT_FB);

// band numbers
	int g2 = b%nd;
	int g3 = ((b%(nd*nd))-(b%nd))/nd;
	int face = (b-(b%(nd*nd)))/(nd*nd);
	
// Size of the scale
	int S1;
	int S2;
	int S3;

	switch(face)
	{
		case 0:
			S1 = nearest_odd(DNX); //cerr<<"x=s1";
			S2 = nearest_odd(DNY); //cerr<<"y=s2";
			S3 = nearest_odd(DNZ); //cerr<<"z=s3";
			break;
		case 1:
			S1 = nearest_odd(DNY); //cerr<<"y=s1";
			S2 = nearest_odd(DNZ); //cerr<<"z=s2";
			S3 = nearest_odd(DNX); //cerr<<"x=s3";
			break;
		case 2:
			S1 = nearest_odd(DNZ); //cerr<<"z=s1";
			S2 = nearest_odd(DNX); //cerr<<"x=s2";
			S3 = nearest_odd(DNY); //cerr<<"y=s3";
			break;
		case 3:
			S1 = nearest_odd(DNX);
			S2 = nearest_odd(DNY);
			S3 = nearest_odd(DNZ);
			break;
		case 4:
			S1 = nearest_odd(DNY);
			S2 = nearest_odd(DNZ);
			S3 = nearest_odd(DNX);
			break;
		case 5:
			S1 = nearest_odd(DNZ);
			S2 = nearest_odd(DNX);
			S3 = nearest_odd(DNY);
			break;
		default:
			break;
	}	

// Offset coord (Center coord)
	int O1 = floor(S1/2);
	int O2 = floor(S2/2);
	int O3 = floor(S3/2);
	
// Wedge maximal half size
	float W1 = 0.5*S1/nd;
	float W2 = 0.5*S2/nd;
	float W3 = 0.5*S3/nd;

// corner semiblock angle
	float da = PI4-atan2(1.0-1.0/float(nd),1.0);
	
// Main coord : Starting and ending (centered)
	float s1 = 1./4.*O1; // theory : 3/8, but not for MeyerW
	float e1 = (float) O1;
	int n1 = int(ceil(e1-s1)+1);
	int x1s = int(ceil(s1));
	int x1e = int(floor(e1));

// Other coord : Starting and ending (not centered)
	float s3,e3;
	int n3;

	s3 = (2*g3-1)*W3;
	e3 = s3 + 4*W3;
	n3 = int(ceil(4*W3));

	float s2,e2;
	int n2;

	s2 = (2*g2-1)*W2;
	e2 = s2 + 4*W2;
	n2 = int(ceil(4*W2));

	int nx = n1*(face%3==0) + n3*(face%3==1) + n2*(face%3==2);
	int ny = n2*(face%3==0) + n1*(face%3==1) + n3*(face%3==2);
	int nz = n3*(face%3==0) + n2*(face%3==1) + n1*(face%3==2);

	if(!backward) FBand.alloc(nx,ny,nz);

// limit angles (start, middle, end) 
	float a2s,a2m,a2e,a3s,a3m,a3e;
	if(g2==0)
	{
		a2s = atan2( -1.0 , 1.-1./nd );
		a2m = atan2( (2.*g2+1.)/nd -1., 1.0 );
		a2e = atan2( (2.*g2+3.)/nd -1., 1.0 );
	}
	else if(g2==nd-1)
	{
		a2s = atan2( (2.*g2-1.)/nd -1., 1.0 );
		a2m = atan2( (2.*g2+1.)/nd -1., 1.0 );
		a2e = atan2( 1.0 , 1.-1./nd );
	}
	else
	{
		a2s = atan2( (2.*g2-1.)/nd -1., 1.0 );
		a2m = atan2( (2.*g2+1.)/nd -1., 1.0 );
		a2e = atan2( (2.*g2+3.)/nd -1., 1.0 );
	}
	if(g3==0)
	{
		a3s = atan2( -1.0 , 1.-1./nd );
		a3m = atan2( (2.*g3+1.)/nd -1., 1.0 );
		a3e = atan2( (2.*g3+3.)/nd -1., 1.0 );
	}
	else if(g3==nd-1)
	{
		a3s = atan2( (2.*g3-1.)/nd -1., 1.0 );
		a3m = atan2( (2.*g3+1.)/nd -1., 1.0 );
		a3e = atan2( 1.0 , 1.-1./nd );
	}
	else
	{
		a3s = atan2( (2.*g3-1.)/nd -1., 1.0 );
		a3m = atan2( (2.*g3+1.)/nd -1., 1.0 );
		a3e = atan2( (2.*g3+3.)/nd -1., 1.0 );
	}

	float aa2m = abs(a2m);
	float aa2s = min(abs(a2s),abs(a2e));
	float aa2e = max(abs(a2e),abs(a2s));
	float aa3m = abs(a3m);
	float aa3s = min(abs(a3s),abs(a3e));
	float aa3e = max(abs(a3e),abs(a3s));

	if(LocVerbose) cerr<<" g2 g3="<<g2<<" "<<g3<<" => a2="<<a2s<<":"<<a2e<<" a3="<<a3s<<":"<<a3e<<"; x1t="<<x1s<<":"<<x1e<<endl;
	
	for(int x1t=x1s;x1t<=x1e;x1t++) // x1temp
	{
		// ranges of tangential directions
		int x2s = max( 0 ,	(int) ceil ( O2 + W2/W1*x1t*tan(a2s) ) );
		int x2e = min( S2-1,(int) floor( O2 + W2/W1*x1t*tan(a2e) ) );
		int x3s = max( 0 ,	(int) ceil ( O3 + W3/W1*x1t*tan(a3s) ) );
		int x3e = min( S3-1,(int) floor( O3 + W3/W1*x1t*tan(a3e) ) );

		for(int x2t=x2s;x2t<=x2e;x2t++)
			for(int x3t=x3s;x3t<=x3e;x3t++)
			{
				int x1,x2,x3;
				double coefw=1;
				float w2=1,w3=1;

				if( normalize_wedge )
				{
					float a2 = atan2( (x2t-O2)/W2, x1t/W1 );// angular position, squeezed to a cube of size x1^3
					float a3 = atan2( (x3t-O3)/W3, x1t/W1 );
					float aa2 = abs(a2);
					float aa3 = abs(a3);
					float aa1 = atan2( abs(x2t-O2)/W2, abs(x3t-O3)/W3 );

					if(a2>a2m)
						w2 = fw((a2-a2m)/(a2e-a2m));
					else // a2<a2m
						w2 = fw((a2m-a2)/(a2m-a2s));

					if(a3>a3m)
						w3 = fw((a3-a3m)/(a3e-a3m));
					else
						w3 = fw((a3m-a3)/(a3m-a3s));

					bool corner=false, inside=false, Face=false, Up=false, Side=false;
					corner = abs(aa2-PI4)<da && abs(aa3-PI4)<da;
					inside = corner && abs(aa1-PI4)<da;
					Face = x1t/W1 >= abs(x2t-O2)/W2 && x1t/W1 >= abs(x3t-O3)/W3 ;
					Up = abs(x3t-O3)/W3 > x1t/W1 && abs(x3t-O3)/W3 > abs(x2t-O2)/W2 ;
					Side = abs(x2t-O2)/W2 > x1t/W1 && abs(x2t-O2)/W2 > abs(x3t-O3)/W3 ;

					if(corner)
					{
						if(inside || Face)	// Being inside Face is not possible but must be included 
											// for rare cases on the surface between inside and face
							coefw = fw((aa2-aa2m)/(aa2e-aa2m)) * fw((aa3-aa3m)/(aa3e-aa3m)) +
									fw((PI2-aa1-aa2m)/(aa2e-aa2m)) * fw((PI2-aa2-aa3m)/(aa3e-aa3m)) +
									fw((PI2-aa3-aa2m)/(aa2e-aa2m)) * fw((aa1-aa3m)/(aa3e-aa3m));
						else if(Up)
						{
							coefw = fw((aa2-aa2m)/(aa2e-aa2m)) * fw((aa3-aa3m)/(aa3e-aa3m)) +
									fw((PI2-aa3-aa2m)/(aa2e-aa2m)) * fw((aa3m-aa1)/(aa3m-aa3s)) +
									fw((PI2-aa3-aa2m)/(aa2e-aa2m)) * fw((aa1-aa3s)/(aa3m-aa3s));
						}
						else if(Side)
							coefw = fw((aa2-aa2m)/(aa2e-aa2m)) * fw((aa3-aa3m)/(aa3e-aa3m)) +
									fw((aa2m-(PI2-aa1))/(aa2m-aa2s)) * fw((PI2-aa2-aa3m)/(aa3e-aa3m)) +
									fw((PI2-aa1-aa2s)/(aa2m-aa2s)) * fw((PI2-aa2-aa3m)/(aa3e-aa3m));
					}
					else // not corner
					{
						// corner in the other faces' referential
						if( abs(aa1-PI4)<da )
						{// the triangles on the current face
							if(abs(aa2-PI4)<da)
							{
								// from the edge point
								if( ((g2==1 || g2==nd-2) && (g3==0 || g3==nd-1)) || ((g3==1 || g3==nd-2) && (g2==0 || g2==nd-1)) )
									coefw = fw((PI2-aa1-aa2m)/(aa2e-aa2m)) * fw((PI2-aa2-aa3e)/(PI2-2*aa3e)) +
											fw((aa2-aa2m)/(aa2e-aa2m)) * fw((aa3e-aa3)/(aa3e-aa3m)) +
											fw((aa2-aa2m)/(aa2e-aa2m)) * fw((aa3-aa3m)/(aa3e-aa3m)) ;
								// from the corner point
								else if( ((g2==0 || g2==nd-1) && (g3==0 || g3==nd-1)) )
									coefw = fw((PI2-aa1-aa2m)/(aa2e-aa2m)) * fw((PI2-aa2-aa3m)/(aa3e-aa3m)) +
											fw((aa2-aa2m)/(aa2e-aa2m)) * fw((aa3m-aa3)/(aa3m-aa3s)) +
											fw((aa2-aa2m)/(aa2e-aa2m)) * fw((aa3-aa3s)/(aa3m-aa3s));
							}
							else if( abs(aa3-PI4)<da )
							{
								// from the edge point
								if( ((g2==1 || g2==nd-2) && (g3==0 || g3==nd-1)) || ((g3==1 || g3==nd-2) && (g2==0 || g2==nd-1)) )
									coefw = fw((PI2-aa3-aa2e)/(PI2-2*aa2e)) * fw((aa1-aa3m)/(aa3e-aa3m)) +
											fw((aa2e-aa2)/(aa2e-aa2m)) * fw((aa3-aa3m)/(aa3e-aa3m)) +
											fw((aa2-aa2m)/(aa2e-aa2m)) * fw((aa3-aa3m)/(aa3e-aa3m)) ;
								// from the corner point
								else if( ((g2==0 || g2==nd-1) && (g3==0 || g3==nd-1)) )
									coefw = fw((PI2-aa3-aa2m)/(aa2e-aa2m)) * fw((aa1-aa3m)/(aa3e-aa3m)) +
											fw((aa2m-aa2)/(aa2m-aa2s)) * fw((aa3-aa3m)/(aa3e-aa3m)) +
											fw((aa2-aa2s)/(aa2m-aa2s)) * fw((aa3-aa3m)/(aa3e-aa3m)) ;
							}
						}// end corner in other referential
					}// end corner or not
				}// end normalize wedge

				// Rotation and symetry for the 6 faces
				x1=x1t+O1;
				x2=x2t;
				x3=x3t;
				int x,y,z,sx,sy,sz, Sx,Sy,Sz;
				switch(face)
				{
					case 0:
						sx=x1; sy=x2; sz=x3;	Sx=S1; Sy=S2; Sz=S3; 
						x=(sx-x1s-O1+n1/2)%nx; y=sy%ny; z=sz%nz; 
						break;
					case 1:
						sx=x3; sy=x1; sz=x2;	Sx=S3; Sy=S1; Sz=S2; 
						x=sx%nx; y=(sy-x1s-O1+n1/2)%ny; z=sz%nz; 
						break;
					case 2:
						sx=x2; sy=x3; sz=x1;	Sx=S2; Sy=S3; Sz=S1; 
						x=sx%nx; y=sy%ny; z=(sz-x1s-O1+n1/2)%nz; 
						break;
					case 3:
						sx=S1-1-x1; sy=S2-1-x2; sz=S3-1-x3; 
						x=(sx-n1/2)%nx; y=sy%ny; z=sz%nz; 
						break;
					case 4:
						sx=S3-1-x3; sy=S1-1-x1; sz=S2-1-x2; 
						x=sx%nx; y=(sy-n1/2)%ny; z=sz%nz; 
						break;
					case 5:
						sx=S2-1-x2; sy=S3-1-x3; sz=S1-1-x1; 
						x=sx%nx; y=sy%ny; z=(sz-n1/2)%nz; 
						break;
					default:
						break;
				}

				double coef;
				if(normalize_wedge==2) coef = coefw <= 1e-8 ? 0 : sqrt(w2*w3/coefw);
				else coef = coefw <= 1e-8 ? 0 : w2*w3/coefw;

				if(!backward) FBand(x,y,z) = coef*complex<double>(Scale(sx,sy,sz));
				else
				{
					Scale(sx,sy,sz) += coef*complex<double>(FBand(x,y,z));
					if(RealBand) Scale(Sx-1-sx,Sy-1-sy,Sz-1-sz) += coef*conj(complex<double>(FBand(x,y,z)));
				}

			}// end for x2t x3t
	}// end of wedge browsing

	if(Verbose)
	{
		if(backward)cerr<<"...End Wedge_into_single_scale"<<endl;
		else cerr<<"...End Scale_into_single_wedge"<<endl;
	}
}


// ********************************************************************

void FCurvelet3D::get_size()
{
	bool LocVerbose = false && Verbose;
	if(Verbose) cerr<<"get size()..."<<endl;
	

	double DNX = nxd();
	double DNY = nyd();
	double DNZ = nzd();
	if(LocVerbose) cerr<<"DNXYZ "<<DNX<<","<<DNY<<","<<DNZ<<endl;
	
	for  (int s=0; s < nbr_scale()-1; s++) 
	{
		int nd = sqrt(TabNbrAnglePerScale(s)/6);
		
	// Size of the scale
		int S1;
		int S2;
		int S3;

		
		int cnt=0;
		for(int face = 0;face<6;face++)
		{
			switch(face)
			{
				case 0:
					S1 = nearest_odd(DNX); //cerr<<"x=s1";
					S2 = nearest_odd(DNY); //cerr<<"y=s2";
					S3 = nearest_odd(DNZ); //cerr<<"z=s3";
					break;
				case 1:
					S1 = nearest_odd(DNY); //cerr<<"y=s1";
					S2 = nearest_odd(DNZ); //cerr<<"z=s2";
					S3 = nearest_odd(DNX); //cerr<<"x=s3";
					break;
				case 2:
					S1 = nearest_odd(DNZ); //cerr<<"z=s1";
					S2 = nearest_odd(DNX); //cerr<<"x=s2";
					S3 = nearest_odd(DNY); //cerr<<"y=s3";
					break;
				case 3:
					S1 = nearest_odd(DNX);
					S2 = nearest_odd(DNY);
					S3 = nearest_odd(DNZ);
					break;
				case 4:
					S1 = nearest_odd(DNY);
					S2 = nearest_odd(DNZ);
					S3 = nearest_odd(DNX);
					break;
				case 5:
					S1 = nearest_odd(DNZ);
					S2 = nearest_odd(DNX);
					S3 = nearest_odd(DNY);
					break;
				default:
					break;
			}	

			// Offset coord (Center coord)
			int O1 = floor(S1/2);
			// Wedge maximal half size
			float W2 = 0.5*S2/nd;
			float W3 = 0.5*S3/nd;

			// Main coord : Starting and ending (centered)
			float s1 = 1./4.*O1; // theory : 3/8, but not for MeyerW
			float e1 = (float) O1;
			int n1 = int(ceil(e1-s1)+1);

			int n2 = int(ceil(4*W2));
			int n3 = int(ceil(4*W3));
			
			for(int h=0; h<nd; h++)
				for(int g=0; g<nd; g++)
				{
					// switch from (x,y,z) to (1,2,3) with a rotation for 6-face symetry
					int nx = n1*(face%3==0) + n3*(face%3==1) + n2*(face%3==2);
					int ny = n2*(face%3==0) + n1*(face%3==1) + n3*(face%3==2);
					int nz = n3*(face%3==0) + n2*(face%3==1) + n1*(face%3==2);
					TabSizeNx[s](cnt) = int(nx);
					TabSizeNy[s](cnt) = int(ny);
					TabSizeNz[s](cnt) = int(nz);
					if(!real()) // complex cube : 2 bands per angle
					{
						TabSizeNx[s](cnt+nbr_angle(s)) = int(nx);
						TabSizeNy[s](cnt+nbr_angle(s)) = int(ny);
						TabSizeNz[s](cnt+nbr_angle(s)) = int(nz);
					}
					cnt++;
				}
		}
		
		DNX =  DNX/2.;
		DNY =  DNY/2.;
		DNZ =  DNZ/2.;
	}//end for scale
	
	int s = nbr_scale()-1;
	TabSizeNx[s](0) = nearest_odd(DNX);
	TabSizeNy[s](0) = nearest_odd(DNY);
	TabSizeNz[s](0) = nearest_odd(DNZ);

	if(Verbose) cerr<<"...End get size"<<endl;
}

// ********************************************************************

void FCurvelet3D::set_tabsigma(fltarray &TS)
{
	if(Verbose) cerr<<"Heritate TabSigma"<<endl;
	TabSigma=TS;
	tabsigma=true;
}

// ********************************************************************

void FCurvelet3D::alloc_with_tab(int Nbr_Scale, int Nx, int Ny, int Nz, intarray & TabDir, Bool ExtendWT, Bool IsotropWT, Bool RealCur)
{
	bool LocVerbose = false & Verbose;
	if(Verbose) cerr<<"alloc_with_tab("<<Nbr_Scale<<","<<Nx<<","<<Ny<<","<<Nz<<",.,"<<ExtendWT<<","<<IsotropWT<<","<<RealCur<<")..."<<endl;
	DataNx = Nx;
	DataNy = Ny;
	DataNz = Nz;
	
	if ((Nx %2 == 0) || (Ny %2 == 0) || (Nz %2 == 0)) ModifSize = True;
	else ModifSize = False;
	if(LocVerbose) cerr<<"ModifSize(even2odd) "<<ModifSize<<endl;
	
	if (Nx %2 == 0) NewNx = Nx + 1;
	else NewNx = Nx;
	if (Ny %2 == 0) NewNy = Ny + 1;
	else NewNy = Ny;
	if (Nz %2 == 0) NewNz = Nz + 1;
	else NewNz = Nz;
	if(LocVerbose) cerr<<"NewNx "<<NewNx<<","<<NewNy<<","<<NewNz<<endl;
	
/*	if(NewNx!=NewNy || NewNx!=NewNz )
	{
		cerr<<"The fast curvelet transform works only on square data."<<endl;
		exit(0);
	}
*/
	// Wavelet init
	init(Nbr_Scale,NewNx,NewNy,NewNz,ExtendWT,IsotropWT,True);
	RealBand = RealCur;

	TabNbrBandPerScale.alloc(Nbr_Scale);
	TabNbrAnglePerScale.alloc(Nbr_Scale);
	if (LocVerbose) cout << " Number of Scales = " << Nbr_Scale << " " << nbr_scale() << endl;

	for (int s=0; s < Nbr_Scale-1; s++)
	{
		TabNbrAnglePerScale(s) = TabDir(s);
		if (real() == False) TabNbrBandPerScale(s) = TabNbrAnglePerScale(s)*2;
		else TabNbrBandPerScale(s) = TabNbrAnglePerScale(s);
		if (LocVerbose) cout << "  Scale " << s << " NbrAnglePerScale = " << TabNbrAnglePerScale(s) << " NbrBandPerScale = " << TabNbrBandPerScale(s) << endl;
	}
	TabNbrAnglePerScale(Nbr_Scale-1) = 1;
	TabNbrBandPerScale(Nbr_Scale-1) = 1;

	TabSizeNx = new intarray [nbr_scale()];
	TabSizeNy = new intarray [nbr_scale()];
	TabSizeNz = new intarray [nbr_scale()];
	
	for(int s=0; s < nbr_scale(); s++) 
	{
		TabSizeNx[s].alloc(TabNbrBandPerScale(s));
		TabSizeNy[s].alloc(TabNbrBandPerScale(s));
		TabSizeNz[s].alloc(TabNbrBandPerScale(s));
	}
	
// Scale/Band size calculation
	get_size();
	
// TabSigma initialisation
	if(LocVerbose) cerr<<" TabSigma initialisation"<<endl;
	TabSigma.resize(NbrScale,TabDir(0));
	char filename[128];
	char extract_string[3];
	if(extract_type==EXTRACT_B) strcpy(extract_string,"BB");
	else if(extract_type==EXTRACT_F) strcpy(extract_string,"FF");
	else strcpy(extract_string,"FB"); // EXTRACT_FB
	
	if(NewNx==NewNy && NewNx==NewNz)
		sprintf(filename,"%s/FastCur3D_norm/TabSigma_N%04d_n%02d_d%02d.fits", getenv("SAP_COM"), 
						NewNx, NbrScale, NbrDir2d);
	else 
		sprintf(filename,"%s/FastCur3D_norm/TabSigma_N%04d_%04d_%04d_n%02d_d%02d.fits", getenv("SAP_COM"), 
						NewNx, NewNy, NewNz, NbrScale, NbrDir2d);
	
//cerr<<filename<<endl;
	if(!tabsigma)
	{
		if(access(filename, F_OK)==0)
		{
			if(Verbose) cerr<<" Using previously calculated TabSigma"<<endl;
			fits_read_fltarr(filename, TabSigma);
//cerr<<TabSigma.nx()<<","<<TabSigma.ny()<<","<<TabSigma.nz()<<endl;
			tabsigma=true;
		}
		else if(Verbose) cerr<<" File "<<filename<<" not found"<<endl;
	}

// TabSigmaNoise init
	TabSigmaNoise.alloc(NbrScale);
	TabSigmaNoise.init(-1.);
	
	if(Verbose) cerr<<"...End alloc_with_tab"<<endl;
}

// ********************************************************************

// Nbr_Scale : nombre de direction total en 2D => NDir3DTotal = 6*(NbrDir/4)^2
void FCurvelet3D::alloc_from_coarse(int Nbr_Scale, int Nx, int Ny, int Nz, int _NbrDir2d, Bool ExtendWT, Bool IsotropWT, Bool Real, type_extraction _extract_type)
{
	bool LocVerbose = true && Verbose ;
	if(Verbose) cerr<<"alloc_from_coarse("<<Nbr_Scale<<","<<Nx<<","<<Ny<<","<<Nz<<","<<_NbrDir2d<<","<<ExtendWT<<","<<IsotropWT<<","<<Real<<")..."<<endl;
	
	if(LocVerbose) cerr<<"ExtendWT "<<ExtendWT<<endl;
	if(LocVerbose) cerr<<"Real "<<Real<<endl;
	
	extract_type = _extract_type;
	NbrScale=Nbr_Scale;
	NbrDir2d = max(_NbrDir2d,12);
	// NbrDir2d cannot be 8 (ie. nd=2)
	// because then the formulas in wedge_scale are not correct,
	// as aX becomes negative in the corner then we cannot use aaX (see corner extraction in scale_to_wedge)
	// anyway 12 is useful as it contains horizontals and verticals, while 8 doesn't.
	
	if(LocVerbose) cerr<<" NbrDir2d="<<NbrDir2d<<endl;
	intarray TabDir;
	TabDir.alloc(Nbr_Scale);
	
	int nd=NbrDir2d/4;
	TabDir(Nbr_Scale-2) = 6*nd*nd;
	if(LocVerbose) cerr<<" tabdir("<<Nbr_Scale-2<<")="<<TabDir(Nbr_Scale-2)<<endl;
	
	int Pas=1;
	for (int s=Nbr_Scale-3; s >=0; s--) 
	{
		nd = nd*(1+Pas);
		Pas = 1-Pas;
		TabDir(s) = 6*nd*nd;
		if(LocVerbose) cerr<<" tabdir("<<s<<")="<<TabDir(s)<<endl;
	}
	alloc_with_tab(Nbr_Scale, Nx, Ny, Nz, TabDir, ExtendWT, IsotropWT, Real);
	
	TabStat.alloc(NbrScale,nbr_band(0),6);
	
	if(Verbose) cerr<<"...End alloc from coarse"<<endl;
}

// ********************************************************************

void FCurvelet3D::alloc_from_fine(int Nbr_Scale, int Nx, int Ny, int Nz, int _NbrDir2d, Bool ExtendWT, Bool IsotropWT, Bool Real)
{
	bool LocVerbose = false && Verbose ;
	if(Verbose) cerr<<"alloc_from_fine("<<Nbr_Scale<<","<<Nx<<","<<Ny<<","<<Nz<<","<<_NbrDir2d<<","<<ExtendWT<<","<<IsotropWT<<","<<Real<<")..."<<endl;
	
	NbrScale=Nbr_Scale;
	NbrDir2d = max(_NbrDir2d,12);
	
	intarray TabDir;
	TabDir.alloc(Nbr_Scale);
	
	int nd=NbrDir2d/4;
	TabDir(Nbr_Scale-2) = 6*nd*nd;
	if(LocVerbose) cerr<<" tabdir("<<Nbr_Scale-2<<")="<<TabDir(Nbr_Scale-2)<<endl;
	
	int Pas=0;
	for (int s=Nbr_Scale-3; s >=0; s--) 
	{
		nd=nd*(1+Pas);
		Pas = 1-Pas;
		TabDir(s) = 6*nd*nd;
		if(LocVerbose) cerr<<" tabdir("<<s<<")="<<TabDir(s)<<endl;
	}
	alloc_with_tab(Nbr_Scale, Nx, Ny, Nz, TabDir, ExtendWT, IsotropWT, Real);
	
	if(Verbose) cerr<<"...End alloc from fine"<<endl;
}

// ********************************************************************

void FCurvelet3D::get_wedges(cfarray * &TabFWT, fltarray** TabBand)
{
	if (Verbose) cerr<<"Get Wedges(.,.)..."<<endl;
	
	bool LVerbose = Verbose;
	
//int x=0;
//cerr<<"get(0,"<<x<<")="<<get_pointer_of_band(0,x)<<", get(0,"<<x+1<<")="<<get_pointer_of_band(0,x+1)<<" diff="<<get_pointer_of_band(0,x+1)-get_pointer_of_band(0,x)<<", size(0,"<<x<<")="<<size_band_nx(0,x)*size_band_ny(0,x)*size_band_nz(0,x)<<"="<<float(size_band_nx(0,x))*float(size_band_ny(0,x))*float(size_band_nz(0,x))<<"("<<size_band_nx(0,x)<<","<<size_band_ny(0,x)<<","<<size_band_nz(0,x)<<")"<<endl;
//x=40;
//cerr<<"get(0,"<<x<<")="<<get_pointer_of_band(0,x)<<", get(0,"<<x+1<<")="<<get_pointer_of_band(0,x+1)<<" diff="<<get_pointer_of_band(0,x+1)-get_pointer_of_band(0,x)<<", size(0,"<<x<<")="<<size_band_nx(0,x)*size_band_ny(0,x)*size_band_nz(0,x)<<"="<<float(size_band_nx(0,x))*float(size_band_ny(0,x))*float(size_band_nz(0,x))<<"("<<size_band_nx(0,x)<<","<<size_band_nx(0,x)<<","<<size_band_nx(0,x)<<")"<<endl;
	for  (int s=0; s < NbrScale; s++)
	{
		#if USE_OMP_FC3D
			#pragma omp parallel for 
		#endif
		for  (int b=0; b < use_NAngle(s); b++)
		{
			Verbose = LVerbose & (b==0);
			fltarray **Band;
			Band = new fltarray*[2];
			Band[0] = &(TabBand[s][b]);
			Band[1] = &(TabBand[s][b+use_NAngle(s)]);
			get_single_wedge(TabFWT, Band, s, b);
//cerr<<"->band"<<Band[0][0](0,0,0)<<"Y";
			delete [] Band;
		}
		Verbose = LVerbose;
	}
	if(Verbose) cerr<<"...End Get Wedge"<<endl;
}

// ********************************************************************

int FCurvelet3D::get_size_of_transform()
{
	int a=0;
	for  (int s=0; s <NbrScale; s++)
		for  (int b=0; b < nbr_band(s); b++)
			a += size_band_nx(s,b)*size_band_ny(s,b)*size_band_nz(s,b);
	return a;
}
float* FCurvelet3D::get_pointer_of_band(int scale, int band)
{
	int a=0;
	for  (int s=0; s <= scale; s++)
	{
		int nband = nbr_band(s);
		if(s==scale) nband = band;
		for  (int b=0; b < nband; b++)
			a += size_band_nx(s,b)*size_band_ny(s,b)*size_band_nz(s,b);
	}
	return &(address[a]);
}
void FCurvelet3D::get_single_wedge(cfarray * &TabFWT, fltarray ** Band, int s, int b)
{
	bool LocVerbose = true && Verbose;
	if (Verbose) cerr<<"get_single_wedge(.,.,"<<s<<","<<b<<")..."<<endl;
	
	cfarray FBand;
	
	int nd;
	double DNX = nxd();
	double DNY = nyd();
	double DNZ = nzd();
	if(LocVerbose) cerr<<" DNXYZ "<<DNX<<","<<DNY<<","<<DNZ<<endl;

	for  (int is=0; is < s; is++) 
	{
		DNX =  DNX/2.;
		DNY =  DNY/2.;
		DNZ =  DNZ/2.;
	}
	
	if(LocVerbose) cerr<<"Scale "<<s<<", UseNAngle="<<use_NAngle(s)<<", real="<<real()<<endl;

// Coarse : copy
// Fine : scale (TabFWT) into wedges (FBand)
	if (nbr_angle(s) > 1)
	{
		if (LocVerbose) cout << " Get Wedges Scale " << s << ", NbrAngle = " << nbr_angle(s) << ", UseNAngle = " << use_NAngle(s) <<  endl;
		nd=(int)sqrt(nbr_angle(s)/6);
		scale_into_single_wedge(TabFWT[s], FBand, b, nd, DNX, DNY, DNZ, false);
	}
	else FBand = TabFWT[s];

// IFFT of the wedges + separate real and imaginary parts if RealTrans
	int nx = FBand.nx();
	int ny = FBand.ny();
	int nz = FBand.nz();

	float Norm = sqrt(float(nx * ny * nz));
	for (int i=0; i < nx; i++)
	for (int j=0; j < ny; j++)
	for (int k=0; k < nz; k++) FBand(i,j,k) *= Norm;

	FFTN_3D FFT3D_local;
	FFT3D_local.CenterZeroFreq = True;
	FFT3D_local.ifftn3d(FBand);



//if real : Separate Real and imaginary parts of the reconstructed bands
// then recopy to TabBand
	if(real())
	{
		if (s != nbr_scale()-1)
		{
			if(pointer)
			{
				Band[0]->alloc(get_pointer_of_band(s,b),nx,ny,nz);
				Band[1]->alloc(get_pointer_of_band(s,b+use_NAngle(s)),nx,ny,nz);
//				float*a=Band[0]->buffer();
//				int n = size_band_nx(s,b)*size_band_ny(s,b)*size_band_nz(s,b);
//				double *c = (double*)a;
//cerr<<"Band ("<<s<<","<<b<<"), pointer="<<a<<":"<<&(a[n-1])<<", nelem="<<n<<"("<<size_band_nx(s,b)<<","<<size_band_nx(s,b)<<","<<size_band_nx(s,b)<<")"<<endl;
//cerr<<"Band ("<<s<<","<<b+use_NAngle(s)<<"), pointer="<<Band[1]->buffer()<<":"<<&(Band[1]->buffer()[n-1])<<", nelem="<<n<<"("<<size_band_nx(s,b+use_NAngle(s))<<","<<size_band_nx(s,b+use_NAngle(s))<<","<<size_band_nx(s,b+use_NAngle(s))<<")"<<endl;
//cerr<<":adiff="<<(int)(&(a[n-1])-a)<<",nelem="<<n;
//cerr<<" ; a "<<a<<" ; a+1 "<<a+1<<" &(a[0])="<<&(a[0])<<",&(a[1])="<<&(a[1])<<", (int)(&(a[1])-&(a[0]))="<<(int)(&(a[1])-&(a[0]))<<", (&(a[1])-&(a[0]))="<<(&(a[1])-&(a[0]))<<", (int(&(a[1]))-int(&(a[0])))="<<(int(&(a[1]))-int(&(a[0])))<<endl;
//cerr<<" ; c "<<c<<" ; c+1 "<<c+1<<" &(c[0])="<<&(c[0])<<",&(c[1])="<<&(c[1])<<", (int)(&(c[1])-&(c[0]))="<<(int)(&(c[1])-&(c[0]))<<", (&(c[1])-&(c[0]))="<<(&(c[1])-&(c[0]))<<", (int(&(c[1]))-int(&(c[0])))="<<(int(&(c[1]))-int(&(c[0])))<<endl;
			}
			else
			{
				Band[0]->alloc(nx,ny,nz);
				Band[1]->alloc(nx,ny,nz);
			}

			for (int i=0; i < nx; i++)
				for (int j=0; j < ny; j++)
					for (int k=0; k < nz; k++)
					{
						float Re = FBand(i,j,k).real()*sqrt(2.);
						float Im = FBand(i,j,k).imag()*sqrt(2.);
						(*Band[0])(i,j,k) = Re ;
						(*Band[1])(i,j,k) = Im ;
					}
		}
		else if (s == nbr_scale()-1)
		{
			if(pointer)
			{
				Band[0]->alloc(get_pointer_of_band(s,b),nx,ny,nz);
//cerr<<"Band ("<<s<<"), pointer:"<<Band[0]->buffer();
			}
			else
				Band[0]->alloc(nx,ny,nz);
			
			for (int i=0; i < nx; i++)
				for (int j=0; j < ny; j++)
					for (int k=0; k < nz; k++)
					{
						float Re = FBand(i,j,k).real();
						(*Band[0])(i,j,k) = Re ;
					}
		}
	}
	else // Complex : recopy to TabBand
	{
		// not done yet, as TabBand is only fltarray for now
	}
	
// Normalization of the wedge
	normalize_band(*Band[0],s,b);
	if (s != nbr_scale()-1)
		normalize_band(*Band[1],s,b+use_NAngle(s));
	
	if(Verbose) cerr<<"...End get_single_wedge"<<endl;
}

// ********************************************************************

void FCurvelet3D::meyer_transform(fltarray & Data)
{
	bool LocVerbose = true & Verbose ;
	if (Verbose) cout << "FCurvelet3D::meyer_transform(.)... " << endl;
	
	if (ModifSize == False) transform(Data);
	else
	{
		fltarray ModCube(NewNx, NewNy, NewNz);

		// Coeur du cube
		for (int i=0; i < DataNx; i++)
			for (int j=0; j < Data.ny(); j++)
				for (int k=0; k < Data.nz(); k++)
					ModCube(i,j,k) = Data(i,j,k);
		
		// Ajout des plans X,Y Z d'allongement
		if (NewNx != DataNx)
			for (int j=0; j < DataNy; j++)
				for (int k=0; k < DataNz; k++)
					ModCube(DataNx,j,k) = Data(DataNx-1,j,k);
		
		if (NewNy != DataNy)
			for (int i=0; i < DataNx; i++)
				for (int k=0; k < DataNz; k++)
					ModCube(i,DataNy,k) = Data(i,DataNy-1,k);
		
		if (NewNz != DataNz)
			for (int i=0; i < DataNx; i++)
				for (int j=0; j < DataNy; j++)
					ModCube(i,j,DataNz) = Data(i,j,DataNz-1);
		
		// Ajout des lignes de raccord entre les plans d'allongement
		if ((NewNx != DataNx) && (NewNy != DataNy))
			for (int k=0; k < DataNz; k++)
				ModCube(DataNx,DataNy,k) = Data(DataNx-1,DataNy-1,k);
		
		if ((NewNx != DataNx) && (NewNz != DataNz))
			for (int j=0; j < DataNy; j++)
				ModCube(DataNx,j,DataNz) = Data(DataNx-1,j,DataNz-1);
		
		if ((NewNy != DataNy) && (NewNz != DataNz))
			for (int i=0; i < DataNx; i++)
				ModCube(i,DataNy,DataNz) = Data(i,DataNy-1,DataNz-1);
		
		// Ajout du point de raccord si le cube était pair trois fois
		if ((NewNx != DataNx) && (NewNy != DataNy) && (NewNz != DataNz))
			ModCube(DataNx,DataNy,DataNz) = Data(DataNx-1,DataNy-1,DataNz-1);
		
		// Wavelet Transform in Fourier domain into Tabcf_WT_Band (MeyerWT::)
		transform(ModCube);
		free_fft_data();
	}
	
	if(LocVerbose)
	{
		// save real part of wavelet transform (in euclidian space and in fourier space)
		char filename[64];
		for  (int s=0; s < NbrScale; s++)
		{
			fltarray TabWavelet;
			cfarray TabWaveletF;
			TabWaveletF=Tabcf_WT_Band[s];

			float Norm =  sqrt((float)(TabNx(s) * TabNy(s) * TabNz(s)));
			TabWavelet.alloc(TabNx(s), TabNy(s), TabNz(s));

			// save module in fourier space
			for (int i=0; i < TabNx(s); i++)
			for (int j=0; j < TabNy(s); j++)
			for (int k=0; k < TabNz(s); k++) TabWavelet(i,j,k) = (TabWaveletF)(i,j,k).real();
			sprintf(filename,"%s_FW%d.fits","out",s);
			writefltarr(filename, TabWavelet);

			// save in euclidian space
			for (int i=0; i < TabNx(s); i++)
			for (int j=0; j < TabNy(s); j++)
			for (int k=0; k < TabNz(s); k++) TabWavelet(i,j,k) = (TabWaveletF)(i,j,k).real()*Norm;
			FFT3D.ifftn3d(TabWaveletF);
			sprintf(filename,"%s_W%d.fits","out",s);
			writefltarr(filename, TabWavelet);
		}
	}
	if (Verbose) cout << " End FCurvelet3D::meyer_transform(.) " << endl;
}

// ********************************************************************

void FCurvelet3D::cur_trans(fltarray &Data, fltarray ** &TabBand, bool allocTB)
{
	bool LocVerbose = true & Verbose;
	if (Verbose) cout << "Cur_trans... " << endl;
	
	if(LocVerbose) cerr<<"Transform, modifsize ="<<(int)ModifSize<<endl;
	if(LocVerbose) cerr<<" NxNyNz, New xyz = "<<DataNx<<" "<<DataNy<<" "<<DataNz<<" "<<NewNx<<" "<<NewNy<<" "<<NewNz<<endl;
	
// Transforms data (may extend it to NewN.) into Tabcf_WT_Band
	meyer_transform(Data);
	
// TabBand alloc, only done in real case, else TabBand not used yet
	if(allocTB && RealBand)
	{
		TabBand = new fltarray*[NbrScale];
		for(int s=0;s<nbr_scale()-1; s++)
			TabBand[s] = new fltarray[nbr_angle(s)];
		TabBand[nbr_scale()-1] = new fltarray[1];
	}

// get wedge doesn't modify Tabcf_WT_Banf, it reads it and creates TabBand 
//cerr<<Tabcf_WT_Band[0](0,0,0).imag()<<"X";	
	get_wedges(Tabcf_WT_Band, TabBand);
//cerr<<TabBand[0][0](0,0,0)<<endl;	

	if (Verbose) cout << "...End Cur_trans" << endl;
}

// ********************************************************************

void FCurvelet3D::put_wedges(cfarray * &TabWT, fltarray ** TabBand)
{
	if (Verbose) cout << "Put wedges()" << endl;
	
	bool LVerbose = Verbose;
	
	for  (int s=0; s < NbrScale; s++)
	{
		TabWT[s].init();
		
		#if USE_OMP_FC3D
			#pragma omp parallel for 
		#endif
		for  (int b=0; b < use_NAngle(s); b++)
		{
			Verbose = LVerbose & (b==0);
			fltarray **Band;
			Band = new fltarray*[2];
			Band[0] = &(TabBand[s][b]);
			Band[1] = &(TabBand[s][b+use_NAngle(s)]);
			put_single_wedge(TabWT, Band, s, b);
			delete [] Band;
		}
		Verbose = LVerbose;
	}

	if (Verbose) cout << "...End Put wedges" << endl;
}

// ********************************************************************

void FCurvelet3D::put_single_wedge(cfarray * &TabWT, fltarray ** Band, int s, int b)
{
	bool LocVerbose = false && Verbose;
	if (Verbose) cout << "Put_single_wedges()" << endl;
	
	int nx = Band[0]->nx();
	int ny = Band[0]->ny();
	int nz = Band[0]->nz();
	
// Normalization of the wedge; inverse
	normalize_band(*Band[0],s,b,true);
	if (s != nbr_scale()-1)
		normalize_band(*Band[1],s,b+use_NAngle(s),true);
	
	cfarray FBand(nx,ny,nz);
	
	double DNX =  nxd();
	double DNY =  nyd();
	double DNZ =  nzd();
	if(LocVerbose) cerr<<" DNXYZ "<<DNX<<","<<DNY<<","<<DNZ<<endl;
	for(int is=0;is<s;is++)
	{
		DNX = DNX/2.;
		DNY = DNY/2.;
		DNZ = DNZ/2.;
	}

	if(LocVerbose) cerr<<"  Echelle s="<<s<<endl;

	// Unif the real and imaginary parts of the bands in euclidian space
	if ( s != (nbr_scale()-1) )
	{
		for (int i=0; i < nx; i ++)
			for (int j=0; j < ny; j ++)
				for (int k=0; k < nz; k ++)
				{
					float Re = (*Band[0])(i,j,k)*sqrt(0.5);
					float Im = (*Band[1])(i,j,k)*sqrt(0.5);
					FBand(i,j,k) =  complex_f(Re,Im);
				}
	}

	// keep only the real part of the coarse scale
	if (s == (nbr_scale()-1))
	{
		for (int i=0; i < nx; i++)
			for (int j=0; j < ny; j++)
				for (int k=0; k < nz; k++)
				{
					float Re = (*Band[0])(i,j,k) ;
					FBand(i,j,k) =  complex_f(Re,0.);
				}
	}

// Fourier transform the wedge
	float Norm = sqrt(float(nx * ny * nz));
	FFTN_3D FFT3D_local;
	FFT3D_local.CenterZeroFreq = True;
	FFT3D_local.fftn3d(FBand);
	
	for (int i=0; i < nx; i ++)
		for (int j=0; j < ny; j ++)
			for (int k=0; k < nz; k ++)
				FBand(i,j,k) /= Norm;

	if (LocVerbose) cout << " Put Wedges Scale " << s << ", NbrDir = " << TabNbrAnglePerScale(s) << endl;
	if( s != (nbr_scale()-1) )
	{
		int nd=sqrt(nbr_angle(s)/6);
		scale_into_single_wedge(TabWT[s], FBand, b, nd, DNX, DNY, DNZ, true);
	}
	else
		TabWT[s]=FBand;

// Normalization of the wedge; forward again (in case of reuse of the coefficients)
	normalize_band(*Band[0],s,b,false);
	if (s != nbr_scale()-1)
		normalize_band(*Band[1],s,b+use_NAngle(s),false);
	

	if (Verbose) cout << "...End Put_single_wedges" << endl;
}

// ********************************************************************

// Reconstruit dans Data, à partir de TabBand (curvelet), en passant par Tabcf_WT_Band, puis IFFT
void FCurvelet3D::cur_recons(fltarray ** TabBand, fltarray &Data)
{
	bool LocVerbose = false && Verbose;
	if (Verbose) cout << "cur_recons()" << endl;

	if (LocVerbose) cout << " DataNx " <<DataNx<<" "<<DataNy<<" "<<DataNz<<" ModifSize "<<(int)ModifSize<< endl;
	Data.resize(DataNx, DataNy, DataNz);

	// From wedges to wavelet bands in fourier
	put_wedges(Tabcf_WT_Band, TabBand); // réinitialise tabcf, lit tabband, ecrit dans tabcf

	if(LocVerbose)
	{
		// Detection des NAN
		for(int s=0;s<nbr_scale();s++)
		{
			//cerr<<" Size TabWTBand "<<(Tabcf_WT_Band[s]).nx()<<" "<<(Tabcf_WT_Band[s]).ny()<<" "<<(Tabcf_WT_Band[s]).nz()<<endl;;
			for (int i=0; i < (Tabcf_WT_Band[s]).nx(); i++)
			for (int j=0; j < (Tabcf_WT_Band[s]).ny(); j++)
			for (int k=0; k < (Tabcf_WT_Band[s]).nz(); k++)
				if(isnan((Tabcf_WT_Band[s])(i,j,k))) cerr<<" CT nan:"<<s<<" "<<i<<" "<<j<<" "<<k<<endl;
		}

		// save real part of wavelet transform (in euclidian and fourier space)
		char filename[64];
		for  (int s=0; s < NbrScale; s++)
		{
			fltarray TabWavelet;
			cfarray TabWaveletF;
			TabWaveletF=Tabcf_WT_Band[s];

			float Norm = sqrt((float)(TabNx(s) * TabNy(s) * TabNz(s)));
			TabWavelet.alloc(TabNx(s), TabNy(s), TabNz(s));

			// save module in fourier space
			for (int i=0; i < TabNx(s); i++)
			for (int j=0; j < TabNy(s); j++)
			for (int k=0; k < TabNz(s); k++) TabWavelet(i,j,k) = (TabWaveletF)(i,j,k).real();
			sprintf(filename,"%s_rFW%d.fits","out",s);
			writefltarr(filename, TabWavelet);

			// save in euclidian space
			FFT3D.ifftn3d(TabWaveletF);
			for (int i=0; i < TabNx(s); i++)
			for (int j=0; j < TabNy(s); j++)
			for (int k=0; k < TabNz(s); k++) TabWavelet(i,j,k) = (TabWaveletF)(i,j,k).real()*Norm;
			sprintf(filename,"%s_rW%d.fits","out",s);
			writefltarr(filename, TabWavelet);
		}
	}

//cerr<<Tabcf_WT_Band[0].nx()<<","<<Tabcf_WT_Band[0].ny()<<","<<Tabcf_WT_Band[0].nz()<<endl;
//cerr<<Tabcf_WT_Band[1].nx()<<","<<Tabcf_WT_Band[1].ny()<<","<<Tabcf_WT_Band[1].nz()<<endl;
//cerr<<Tabcf_WT_Band[2].nx()<<","<<Tabcf_WT_Band[2].ny()<<","<<Tabcf_WT_Band[2].nz()<<endl;
	// Wavelet reconstruction
	if (ModifSize == False)
		recons(Data);
	else // extract the smaller original cube
	{
		fltarray ModCube(NewNx, NewNy, NewNz, (char*)"New");
		recons(ModCube);
		for (int i=0; i < Data.nx(); i++)
			for (int j=0; j < Data.ny(); j++)
				for (int k=0; k < Data.nz(); k++)
					Data(i,j,k) =  ModCube(i,j,k);
	}
	
	if(use_min)
		for (int i=0; i < Data.nx(); i++)
			for (int j=0; j < Data.ny(); j++)
				for (int k=0; k < Data.nz(); k++)
					if( Data(i,j,k) < min_value ) Data(i,j,k) = min_value;

	if(use_max)
		for (int i=0; i < Data.nx(); i++)
			for (int j=0; j < Data.ny(); j++)
				for (int k=0; k < Data.nz(); k++)
					if( Data(i,j,k) > max_value ) Data(i,j,k) = max_value;

	if (Verbose) cout << "...End cur_recons" << endl;
}

// ********************************************************************

void FCurvelet3D::get_band(int s, int b, fltarray **TabBand, fltarray &band)
{
//	if(Verbose) cerr << "FCurvelet3D::get_band("<<s<<","<<b<<".,.)..." <<  endl;
	band=TabBand[s][b];
//	if(Verbose) cerr << "End FCurvelet3D::get_band" <<  endl;
}

// ********************************************************************

void FCurvelet3D::estim_normalization(char* Name_Imag_Out)
{
	if(Verbose) cerr << "BFCurvelet3D::estim_normalization("<<Name_Imag_Out<<")..."<< endl;
	
// Transform iditentically 0
	fltarray CubeNoise(DataNx,DataNy,DataNz);
	fltarray **TabBand;
	bool allocTB = true;
	cur_trans(CubeNoise, TabBand, allocTB);

	no_norm=true;
	
// we must estimate the normalization coefficients at all scales, for about 1/8 of one face
// the rest of the bands will be copied
	fltarray Recons(DataNx, DataNy, DataNz);
	for(int S=0;S<NbrScale;S++)
	{
		bool ok=false;
		//int nd = NbrDir2d/4;
		int nd = sqrt(TabNbrAnglePerScale(S)/6);
		int nd2 = (nd+1)/2;
		int Bx = 0;
		int By = 0;
		
	// Norm of the current scale
		float norm = get_norm(S)*sqrt( (float(nxs(S))*float(nys(S))*float(nzs(S))) / (float(nbr_band(S))*float(size_band_nx(S,0))*float(size_band_ny(S,0))*float(size_band_nz(S,0))) ); 
		
		while(!ok)
		{
			int B = Bx+nd*By;
			Recons.init();
			
		// Normalized Dirac
			TabBand[S][B](0,0,0) = sqrt(float(DataNx)*float(DataNy)*float(DataNz));
			
		// From wedges to wavelet bands in fourier
			for  (int s=0; s < NbrScale; s++)
				Tabcf_WT_Band[s].init();

			fltarray **Band;
			Band = new fltarray*[2];
			Band[0] = &(TabBand[S][B]);
			Band[1] = &(TabBand[S][B+use_NAngle(S)]);
			put_single_wedge(Tabcf_WT_Band, Band, S, B);
			delete [] Band;

		// Wavelet reconstruction
			if (ModifSize == False)
				recons(Recons);
			else // extract the smaller original cube
			{
				fltarray ModCube(NewNx, NewNy, NewNz, (char*)"New");
				recons(ModCube);
				for (int i=0; i < Recons.nx(); i++)
					for (int j=0; j < Recons.ny(); j++)
						for (int k=0; k < Recons.nz(); k++)
							Recons(i,j,k) =  ModCube(i,j,k);
			}

			double sig=0.;
			for (int i=0; i < Recons.nx(); i++)
			for (int j=0; j < Recons.ny(); j++)
			for (int k=0; k < Recons.nz(); k++)
				sig+=Recons(i,j,k)*Recons(i,j,k);
			TabSigma(S,B) = sqrt(sig/(Recons.nx()*Recons.ny()*Recons.nz()));
			
			if(Verbose) cerr<<"Norm("<<S<<","<<B<<"("<<Bx<<","<<By<<")) = "<<TabSigma(S,B)<<" / "<<norm<<" => "<<TabSigma(S,B)/norm<<endl;
			TabSigma(S,B) /= norm;
			
		// apply the symmetries : 6 faces, 8 symmetries inside each face
			for(int f=0; f<6;f++) 
			{
				TabSigma(S,Bx		+ By		*nd+f*(nbr_band(S)/6)) = TabSigma(S,B);
				TabSigma(S,By		+ Bx		*nd+f*(nbr_band(S)/6)) = TabSigma(S,B);
				TabSigma(S,(nd-1-Bx)+ By		*nd+f*(nbr_band(S)/6)) = TabSigma(S,B);
				TabSigma(S,(nd-1-By)+ Bx		*nd+f*(nbr_band(S)/6)) = TabSigma(S,B);
				TabSigma(S,Bx		+ (nd-1-By)	*nd+f*(nbr_band(S)/6)) = TabSigma(S,B);
				TabSigma(S,By		+ (nd-1-Bx)	*nd+f*(nbr_band(S)/6)) = TabSigma(S,B);
				TabSigma(S,(nd-1-Bx)+ (nd-1-By) *nd+f*(nbr_band(S)/6)) = TabSigma(S,B);
				TabSigma(S,(nd-1-By)+ (nd-1-Bx) *nd+f*(nbr_band(S)/6)) = TabSigma(S,B);
			}
			
		// delete the normalized dirac
			TabBand[S][B](0,0,0) = 0.;
			
		// Next band if applicable
			if(Bx<nd2-1)
				Bx++;
			else if(By<nd2-1)
			{
				By++;
				Bx=By;
			}
			else
				ok=true;
		}
	}
	
	if(Name_Imag_Out!=NULL)
	{
		char filename[64];
		sprintf(filename,"%s_norm.fits",Name_Imag_Out);
		Ifloat toto;
		toto.alloc(TabSigma.buffer(), TabSigma.ny(),TabSigma.nx());
		io_write_ima_float(filename, toto);
//		writefltarr(filename, TabSigma);
//		fits_write_image(filename, TabSigma);
	}
	
	for(int s=0;s<NbrScale-1;s++)
		delete [] TabBand[s];
	
	if(Verbose) cerr << "End BFCurvelet3D::estim_normalization"<<endl;
}

// ********************************************************************

void FCurvelet3D::get_norm_coeff(float N_Sigma, char* Name_Imag_Out)
{
	if(Verbose) cerr << "FCurvelet3D::get_norm_coeff("<<N_Sigma<<","<<Name_Imag_Out<<")..."<< endl;
	float Sig=1.;
	unsigned int InitRnd = 10;

	bool oldverbose=Verbose;
	Verbose=false;
	
	// generate a random gaussian cube
	fltarray CubeNoise(DataNx,DataNy,DataNz);
	randomn(CubeNoise, DataNx*DataNy*DataNz, Sig, InitRnd);
	tabsigma = false; // to prevent normalization if existing coefficients
	fltarray **TabBand;
	bool allocTB = true;
	cur_trans(CubeNoise, TabBand, allocTB);
	calib_noise_nsig(TabBand, N_Sigma, Name_Imag_Out);
	tabsigma = true;
	for(int s=0;s<NbrScale;s++)
		delete [] TabBand[s];
	delete [] TabBand;
	
	Verbose=oldverbose;
	
	if(Verbose) cerr << "...End FCurvelet3D::get_norm_coeff" << endl;
}

// ********************************************************************

void FCurvelet3D::calib_noise_nsig(fltarray ** TabBand, float N_Sigma, char* Name_Imag_Out)
{
	if(Verbose) cerr << "FCurvelet3D::calib_noise_nsig..." <<  endl;
	bool LocVerbose = false & Verbose;
	double LMax;
	fltarray Band;
	
	TabSigma.init(0); // useless when there is no bug
	
	for(int s=0;s<NbrScale-1;s++)
		for(int b=0;b<nbr_band(s);b++)
		{
			float Nsig = N_Sigma;
			get_band(s,b,TabBand, Band);
			
			CArrProb *TabCP;
			TabCP = new CArrProb;
			
			TabCP->set(Band,true);
			TabCP->find_gthreshold(Nsig, LMax);
			TabSigma(s,b) = LMax / N_Sigma;
			delete TabCP;
			if(LocVerbose) cerr << "Band (" << s<<","<<b << ") Lmax = " << TabSigma(s,b)<< endl;
		}
		
	if(Name_Imag_Out!=NULL)
	{
		char filename[64];
		sprintf(filename,"%s_nsig.fits",Name_Imag_Out);
		writefltarr(filename, TabSigma);
	}
	
	if(Verbose) cerr << "End FCurvelet3D::calib_noise_nsig" <<  endl;
}

// ********************************************************************
/*
// The same as extract_stat except that there is no normalisation 
void FCurvelet3D::noise_calibration(fltarray **TabBand, char* Outname)
{
	if(Verbose) cerr<<"FCurvelet3D::Noise_calibration(.,"<<Outname<<")..."<<endl;
	bool LocVerbose = false & Verbose;

// Output stat file
	char Statname[250];
	strcpy(Statname, Outname);
	strcat(Statname, "_noise_calib.dat");
//	fstream cstat;
//	cstat.open (Statname, fstream::out);

	double var;
	double mean;

// fine scales : 0..N-2
	for(int s=0;s<NbrScale-1;s++)
	{
		for(int b=0;b<nbr_band(s);b++)
		{
			mean=0;	// mean
			var=0;	// variance
			
			int nx = size_band_nx(s,b);
			int ny = size_band_ny(s,b);
			int nz = size_band_nz(s,b);
			for(int i=0;i<nx;i++)
				for(int j=0;j<ny;j++)
					for(int k=0;k<nz;k++)
					{
						mean += (TabBand[s][b])(i,j,k);
						var += pow( (TabBand[s][b])(i,j,k) , 2 );
					}
			mean /= nx*ny*nz;
			var  /= nx*ny*nz;
			TabSigma(s,b) = sqrt(var);

		// Write output
			if(LocVerbose)	cerr<< "\tStat band ("<<s<<","<<b<<") : (mean,std) = ("<<mean<<",\t"<<TabSigma(s,b)<<")"<<endl;
//			cstat <<mean<<"\t"<<TabSigma(s,b)<<endl;
		}
	}

// coarse scale : N-1
	// No calibration for the coarsest scale : no "/TabSigma(N-1,0)" defined, only stats

	char filename[64];
	sprintf(filename,"%s_nsig.fits",Outname);
	writefltarr(filename, TabSigma);
	
//	cstat.close();
	if(Verbose) cerr<<"...End Noise_calibration"<<endl;
}
*/
// ********************************************************************

void FCurvelet3D::extract_stat(fltarray ** TabBand, char* Outname, bool not_centered)
{
//	bool LocVerbose = false & Verbose;
	if(Verbose) cerr<<"FCurvelet3D::Extract_stat..."<<endl;

// Output stat file
	char Statname[250];
	fstream cstat;
	if(Outname!=NULL)
	{
		strcpy(Statname, Outname);
		strcat(Statname, "_stat.dat");
		cstat.open (Statname, fstream::out);
	}
	
	for(int s=0;s<NbrScale;s++)
		for(int b=0;b<nbr_band(s);b++)
		{
			double Mean, Sig, Skew, Curt;
			float Mini, Maxi;
			moment4(TabBand[s][b].buffer(), TabBand[s][b].n_elem(), Mean, Sig, Skew, Curt, Mini, Maxi, not_centered);

		// Output
			TabStat(s,b,0)=Mean; TabStat(s,b,1)=Sig; TabStat(s,b,2)=Skew; TabStat(s,b,3)=Curt;
			TabStat(s,b,4)=Mini; TabStat(s,b,5)=Maxi; 
			if(Verbose)
			{
				cerr<<" Stat scale("<<s<<") : " <<Mean<<","<<Sig<<","<<Skew<<","<<Curt<<","<<endl;
				cerr << "  MinMaxCoef Values = "<< Mini << ":"<< Maxi<<endl;
			}
			if(Outname!=NULL)
				cstat << s << "\t" << b <<"\t"<< Mean <<"\t"<< Sig <<"\t"<< Skew<<"\t"<<Curt<<"\t"<<endl;
		}

	if(Outname!=NULL)
		cstat.close();

	if(Verbose) cerr<<"...End FCurvelet3D::Extract_stat"<<endl;
}

// ********************************************************************

void FCurvelet3D::normalize_band(fltarray &TabBand, int s, int b, bool inverse)
{
	if(!no_norm)
	{
		float norm,coef=1.;
		if(tabsigma) // use exact, measured coefficients
			coef = TabSigma(s,b);

		if(s==NbrScale-1)//coef*
			norm = coef*get_norm(s);
		else // sqrt( N_band/N_scale * 1/Ndir )
			norm = coef*get_norm(s)*sqrt( (float(nxs(s))*float(nys(s))*float(nzs(s))) / (float(nbr_band(s))*float(size_band_nx(s,b))*float(size_band_ny(s,b))*float(size_band_nz(s,b))) );
//cerr<<s<<","<<b<<":"<<get_norm(s)<<"*"<<nxs(s)<<"*"<<nys(s)<<"*"<<nzs(s)<<"/"<<
//		nbr_band(s)<<"*"<<size_band_nx(s,b)<<"*"<<size_band_ny(s,b)<<"*"<<size_band_nz(s,b)<<"="<<norm<<" : "<<coef<<endl;
//cerr<<norm<<endl;
		for(int i=0;i<size_band_nx(s,b);i++)
			for(int j=0;j<size_band_ny(s,b);j++)
				for(int k=0;k<size_band_nz(s,b);k++)
					if(inverse)
						TabBand(i,j,k) *= norm ; 
					else
						TabBand(i,j,k) /= norm ; 
	}
}

// ********************************************************************
float xtnorm(float a)
{
	return (abs(a)*30. + 1.);
}

void FCurvelet3D::extern_band_norm(fltarray & Band, fltarray & Band_norm, bool backward)
{
	if(!backward)
		for(int i=0;i<Band.nx();i++)
			for(int j=0;j<Band.ny();j++)
				for(int k=0;k<Band.nz();k++)
					Band(i,j,k) /= xtnorm(Band_norm(i,j,k));
	else
		for(int i=0;i<Band.nx();i++)
			for(int j=0;j<Band.ny();j++)
				for(int k=0;k<Band.nz();k++)
					Band(i,j,k) *= xtnorm(Band_norm(i,j,k));
}

// ********************************************************************

void FCurvelet3D::filter(fltarray &Data, fltarray &Recons, float SigmaNoise, float NSigma, filter_type FilterType, 
		float Alpha, // FDR
		int LocalBS, // WIENER
		float Rho, float Kmin, float Kmean, float Lmax, // Contrast enhancement
		bool force4sigma)
{
	bool LocVerbose = true & Verbose;
	if (Verbose) cout << "FCurvelet3D::XXfilter(.,.,"<<SigmaNoise<<","<<NSigma<<","<<FilterType<<","<<Alpha<<","<<LocalBS<<","<<force4sigma<<")"<< endl;
	if(LocVerbose) cerr<<"Transform, modifsize ="<<(int)ModifSize<<endl;
	if(LocVerbose) cerr<<" NxNyNz, New xyz = "<<DataNx<<" "<<DataNy<<" "<<DataNz<<" "<<NewNx<<" "<<NewNy<<" "<<NewNz<<endl;

	meyer_transform(Data);
	
	cfarray* TabFWBand = new cfarray[NbrScale];
	bool LVerbose = Verbose;
	
	// Band extraction, filtering and reconstruction
	int NS = NbrScale-1 + th_coarse;
	for  (int s=0; s < NS; s++)
	{
		TabFWBand[s].alloc(Tabcf_WT_Band[s].nx(),Tabcf_WT_Band[s].ny(),Tabcf_WT_Band[s].nz());
		#if USE_OMP_FC3D
			#pragma omp parallel for 
		#endif
		for  (int b=0; b < max(nbr_angle(s)/2,1); b++)
		{
		//	Verbose = LVerbose & (b==0);
			fltarray *Band;
			Band = new fltarray[2];
			fltarray **Band_address;
			Band_address = new fltarray*[2];
			Band_address[0]=&(Band[0]); Band_address[1]=&(Band[1]);	
			get_single_wedge(Tabcf_WT_Band, Band_address, s, b);
//		cerr<<"TabBand_norm="<<b<<endl;
			if(TabBand_norm!=NULL) extern_band_norm(Band[0], TabBand_norm[s][b]);
			if(FilterType==FT_HARD || FilterType==FT_SOFT ) threshold_single(Band, s, b, SigmaNoise, NSigma, FilterType, force4sigma);
			else if(FilterType==FT_FDR) fdr_single(Band, s, b, Alpha, SigmaNoise);
			else if(FilterType==FT_WIENER) wiener_single(Band, s, b, SigmaNoise*NSigma/3., LocalBS);
			else if(FilterType==FT_SBT) stein_block_threshold_single(Band, s, b, SigmaNoise);
			else if(FilterType==FT_CONTRAST && s!=NbrScale-1) enhance_single(Band, s, b, Rho, SigmaNoise*Kmin, SigmaNoise*Kmean, Lmax);
			if(TabBand_norm!=NULL) extern_band_norm(Band[0], TabBand_norm[s][b], true);
			put_single_wedge(TabFWBand, Band_address, s, b);
			delete [] Band;
		}
	}
	for  (int s=0; s < NS; s++)
		Tabcf_WT_Band[s] = TabFWBand[s];
	// if NS<<NbrScale, Tabcf_WT_Band[NbrScale-1] is unchanged, no copy needed
	delete [] TabFWBand;
	
	Verbose = LVerbose;
	
	Recons.alloc(Data.nx(),Data.ny(),Data.nz());

	// Wavelet reconstruction
	alloc_fft_data();
	if (ModifSize == False)
		recons(Recons);
	else // extract the smaller original cube
	{
		fltarray ModCube(NewNx, NewNy, NewNz);
		recons(ModCube);
		for (int i=0; i < Data.nx(); i++)
			for (int j=0; j < Data.ny(); j++)
				for (int k=0; k < Data.nz(); k++)
					Recons(i,j,k) =  ModCube(i,j,k);
	}
	
	if(use_min)
		for (int i=0; i < Data.nx(); i++)
			for (int j=0; j < Data.ny(); j++)
				for (int k=0; k < Data.nz(); k++)
					if( Data(i,j,k) < min_value ) Data(i,j,k) = min_value;

	if(use_max)
		for (int i=0; i < Data.nx(); i++)
			for (int j=0; j < Data.ny(); j++)
				for (int k=0; k < Data.nz(); k++)
					if( Data(i,j,k) > max_value ) Data(i,j,k) = max_value;

	if (Verbose) cout << "...End filter" << endl;
}

// ********************************************************************

float FCurvelet3D::get_noise_lvl(float sigma, int s)
{
	if(TabSigmaNoise(0)<-0.5)
		return sigma;
	else
		return TabSigmaNoise(s);
}

// ********************************************************************
// Enhancement between Tmin and Tmax
float coef_enhance(float val, float Rho, float Tmin, float Tmean, float Tmax)
{
	Tmax = max(2*Tmin,Tmax);
	if(val < Tmin) return 1.F;
	if(val > Tmax) return 1.F;
	if(val < 2*Tmin) return (val-Tmin)/Tmin * pow(Tmax/(2.F*Tmin),Rho) + (2*Tmin-val)/Tmin;
	return pow(Tmax/val,Rho);
}

// ********************************************************************
// Rho : non linearity parameter
// enhance from Tmin to max(coef)*Lmax
void FCurvelet3D::enhance(fltarray **TabBand, float Rho, float Tmin, float Tmean, float Lmax)
{
	if(Verbose) cerr<<"FCurvelet3D::enhance(.,"<<Rho<<","<<Tmin<<","<<Tmean<<","<<Lmax<<")"<<endl;

// Fine scales
	for(int s=0;s<NbrScale-1;s++)
	{
		for(int b=0;b<max(nbr_angle(s)/2,1);b++)
		{
			int nx = TabBand[s][b].nx();
			int ny = TabBand[s][b].ny();
			int nz = TabBand[s][b].nz();
			fltarray *Band = new fltarray[2];
			Band[0].alloc(TabBand[s][b].buffer(),nx,ny,nz);
			Band[1].alloc(TabBand[s][b+use_NAngle(s)].buffer(),nx,ny,nz);
			enhance_single(Band, s, b, Rho, Tmin, Tmean, Lmax);
		}
	}
	
	if(no_coarse)
		TabBand[NbrScale-1][0].init(0.0);
	
	if(no_fine)
		for(int b=0;b<nbr_band(0);b++)
			TabBand[0][b].init(0.0);
	
	if(Verbose) cerr<<"...End FCurvelet3D::enhance"<<endl;
}

// *********************************************************************

void FCurvelet3D::enhance_single(fltarray *Band, int s, int b, float Rho, float Tmin, float Tmean, float Lmax)
{// Lmax >0 : max = Lmax*max(band) ; Lmax<0 : max = -Lmax (->Lmax must be -K*sigma)
	if(Verbose) cerr<<"FCurvelet3D::enhance_single(.,"<<Rho<<","<<Tmin<<","<<Tmean<<","<<Lmax<<")"<<endl;
	int nx = Band[0].nx();
	int ny = Band[0].ny();
	int nz = Band[0].nz();
	
	float m0,m1;
	if(Lmax>0)
	{
		m0 = Band[0].max();
		m1 = Band[1].max();
	}
	else// Lmax is -Kmax*SigmaNoise
	{
		m0 = -1.;
		m1 = -1.;
	}

	for(int k = 0; k < nz; k++)
		for(int j = 0; j < ny; j++)
			for(int i = 0; i < nx; i++)
				(Band[0])(i,j,k) *= coef_enhance(abs((Band[0])(i,j,k)), Rho, Tmin, Tmean, m0*Lmax);
	
	if(s != NbrScale-1)
	for(int k = 0; k < nz; k++)
		for(int j = 0; j < ny; j++)
			for(int i = 0; i < nx; i++)
				(Band[1])(i,j,k) *= coef_enhance(abs((Band[1])(i,j,k)), Rho, Tmin, Tmean, m1*Lmax);
	
	if(Verbose) cerr<<"...End FCurvelet3D::enhance_single"<<endl;
}

// *********************************************************************

void FCurvelet3D::threshold(fltarray **TabBand, float SigmaNoise, float NSigma, filter_type FilterType, bool force4)
{ //!!!!! only applicable for real transform
	if(Verbose) cerr<<"FCurvelet3D::threshold(.,"<<SigmaNoise<<","<<NSigma<<","<<FilterType<<","<<force4<<")"<<endl;

	int NS = NbrScale-1 + th_coarse;

// Fine scales
	for(int s=0;s<NS;s++)
	{
		for(int b=0;b<max(nbr_angle(s)/2,1);b++)
		{
			int nx = TabBand[s][b].nx();
			int ny = TabBand[s][b].ny();
			int nz = TabBand[s][b].nz();
			fltarray *Band = new fltarray[2];
			Band[0].alloc(TabBand[s][b].buffer(),nx,ny,nz);
			Band[1].alloc(TabBand[s][b+use_NAngle(s)].buffer(),nx,ny,nz);
			threshold_single(Band, s, b, SigmaNoise, NSigma, FilterType, force4);
		}
	}
	
	if(no_coarse)
		TabBand[NbrScale-1][0].init(0.0);
	
	if(no_fine)
		for(int b=0;b<nbr_band(0);b++)
			TabBand[0][b].init(0.0);
	
	if(Verbose) cerr<<"...End FCurvelet3D::threshold"<<endl;
}

// *********************************************************************

void FCurvelet3D::threshold_single(fltarray * Band, int s, int b, float SigmaNoise, float NSigma, filter_type FilterType, bool force4)
{ //!!!!! only applicable for real transform
	bool LocVerbose = true & Verbose;
	if(LocVerbose) cerr<<"FCurvelet3D::threshold_single(.,"<<s<<","<<b<<","<<SigmaNoise<<","<<NSigma<<","<<FilterType<<","<<force4<<")"<<endl;
	float lvl;
	float cnt=0;
	int nx = Band[0].nx();
	int ny = Band[0].ny();
	int nz = Band[0].nz();
	int tot = 2*nx*ny*nz;
	float Nsig=NSigma;
	if(s==0) Nsig = (force4 ? (NSigma+1) : NSigma);
	lvl = get_noise_lvl(SigmaNoise,s) * Nsig;
//cerr<<"threshold band "<<s<<","<<b<<" : "<<lvl<<endl;
//cerr<<"BAND.nx()"<<Band[0].nx()<<endl;
//cerr<<"BAND.nx()"<<Band[1].nx()<<endl;
	for(int i = 0; i < nx; i++)
		for(int j = 0; j < ny; j++)
			for(int k = 0; k < nz; k++)
			{
				if( abs((Band[0])(i,j,k)) < lvl )
				{
					cnt++;
					(Band[0])(i,j,k)=0; // hard
				}
				else if(FilterType==FT_SOFT) (Band[0])(i,j,k) -= (2*int( (Band[0])(i,j,k)>0 )-1)*lvl;
			}
	if(s != NbrScale-1)
	for(int i = 0; i < nx; i++)
		for(int j = 0; j < ny; j++)
			for(int k = 0; k < nz; k++)
			{
				if( abs((Band[1])(i,j,k)) < lvl )
				{
					cnt++;
					(Band[1])(i,j,k)=0; // hard
				}
				else if(FilterType==FT_SOFT) (Band[1])(i,j,k) -= (2*int( (Band[1])(i,j,k)>0 )-1)*lvl;
			}
			
//cerr<<"threshold band "<<s<<","<<b<<" : "<<get_noise_lvl(SigmaNoise,s) << ","<< Nsig << ","<<lvl<<endl;
	if(LocVerbose) cerr<<" Scale "<<s<<", n#proportion non seuillee ("<<lvl<<")="<<tot-cnt<<"#"<<(tot-cnt)/tot<<endl;
	
	if(LocVerbose) cerr<<"...End FCurvelet3D::threshold_single"<<endl;
}

// *********************************************************************

void FCurvelet3D::soft_threshold(fltarray **TabBand, float lvl, bool threshold_coarse)
{ 
	bool LocVerbose = false & Verbose;
	
	int NS = NbrScale-1 + threshold_coarse;

// All scales
	for(int s=0;s<NS;s++)
	{
		for(int b=0;b<nbr_band(s);b++)
		{
			float cnt=0;
			int nx = TabBand[s][b].nx();
			int ny = TabBand[s][b].ny();
			int nz = TabBand[s][b].nz();
			int tot = nx*ny*nz;

			for(int i = 0; i < nx; i++)
				for(int j = 0; j < ny; j++)
					for(int k = 0; k < nz; k++)
					{
						//(TabBand[s][b])(i,j,k) = soft_threshold((TabBand[s][b])(i,j,k),lvl);
						if( abs((TabBand[s][b])(i,j,k)) < lvl )
						{
							cnt++;
							(TabBand[s][b])(i,j,k)=0; // hard
						}
						else (TabBand[s][b])(i,j,k) -= (2*( (TabBand[s][b])(i,j,k)>0 )-1)*lvl;
						
					}
			if(LocVerbose) cerr<<" Scale "<<s<<", n#proportion non seuillee ("<<lvl<<")="<<tot-cnt<<"#"<<(tot-cnt)/tot<<endl;
		}
	}
	
	if(no_coarse)
		TabBand[NbrScale-1][0].init(0.0);
	
	if(no_fine)
		for(int b=0;b<nbr_band(0);b++)
			TabBand[0][b].init(0.0);
	
	if(Verbose) cerr<<"...End FCurvelet3D::threshold"<<endl;
}

 
void FCurvelet3D::hard_threshold(fltarray **TabBand, float lvl, bool threshold_coarse)
{ 
	bool LocVerbose = false & Verbose;
	
	int NS = NbrScale-1 + threshold_coarse;
    
    // All scales
	for(int s=0;s<NS;s++)
	{
		for(int b=0;b<nbr_band(s);b++)
		{
			float cnt=0;
			int nx = TabBand[s][b].nx();
			int ny = TabBand[s][b].ny();
			int nz = TabBand[s][b].nz();
			int tot = nx*ny*nz;
            
			for(int i = 0; i < nx; i++)
				for(int j = 0; j < ny; j++)
					for(int k = 0; k < nz; k++)
					{
						//(TabBand[s][b])(i,j,k) = soft_threshold((TabBand[s][b])(i,j,k),lvl);
						if( abs((TabBand[s][b])(i,j,k)) < lvl )
						{
							cnt++;
							(TabBand[s][b])(i,j,k)=0; // hard
						}
 					}
			if(LocVerbose) cerr<<" Scale "<<s<<", n#proportion non seuillee ("<<lvl<<")="<<tot-cnt<<"#"<<(tot-cnt)/tot<<endl;
		}
	}
	
	if(no_coarse)
		TabBand[NbrScale-1][0].init(0.0);
	
	if(no_fine)
		for(int b=0;b<nbr_band(0);b++)
			TabBand[0][b].init(0.0);
	
	if(Verbose) cerr<<"...End FCurvelet3D::threshold"<<endl;
}

// *********************************************************************

void FCurvelet3D::wiener(fltarray ** &TabBand, float noise_lvl, int LocalBS)
{
	if(Verbose) cerr<<"FCurvelet3D::wiener("<<noise_lvl<<","<<LocalBS<<")..."<<endl;
	
	// Fine wavelet scales
	for(int s=0;s<NbrScale-1;s++)
	{
		for(int b=0;b<max(nbr_angle(s)/2,1);b++)
		{
			int nx = size_band_nx(s,b);
			int ny = size_band_ny(s,b);
			int nz = size_band_nz(s,b);
			fltarray *Band = new fltarray[2];
			Band[0].alloc(TabBand[s][b].buffer(),nx,ny,nz);
			Band[1].alloc(TabBand[s][b+use_NAngle(s)].buffer(),nx,ny,nz);
			wiener_single(Band, s, b, noise_lvl, LocalBS);
		}
	}// end scale

	// Coarse scale
	// Do not denoise the coarsest scale

	if(Verbose) cerr<<"...End FCurvelet3D::wiener"<<endl;
}

// *********************************************************************

void FCurvelet3D::wiener_single(fltarray * Band, int s, int b, float noise_lvl, int LocalBS)
{ //!!!!! only applicable for real transform
	if(Verbose) cerr<<"FCurvelet3D::wiener_single("<<noise_lvl<<","<<LocalBS<<")..."<<endl;
	float val;
	float noise2 = get_noise_lvl(noise_lvl,s)*get_noise_lvl(noise_lvl,s);
	
	int nx = size_band_nx(s,b);
	int ny = size_band_ny(s,b);
	int nz = size_band_nz(s,b);
	int Nx = ceil(float(nx)/float(LocalBS));
	int Ny = ceil(float(ny)/float(LocalBS));
	int Nz = ceil(float(nz)/float(LocalBS));

	fltarray coef_wiener(nx,ny,nz);
	for(int RI=0;RI<2;RI++)
	{
		coef_wiener.init(-2);
		// Wiener blocks
		for(int kx=0 ; kx < Nx ; kx++)
			for(int ky=0 ; ky < Ny ; ky++)
				for(int kz = 0; kz < Nz; kz++)
				{
					double sigma2 = 0.0;
					float cnt=pow((float)LocalBS,3);

				// Sigma calculation
					// Pixels in a wiener block (angle)
					for(int bx = 0; bx < LocalBS; bx++)
						for(int by = 0; by < LocalBS; by++)
							for(int bz = 0; bz < LocalBS; bz++)
							{
								//cnt+=1;
								int x = (bx + kx*LocalBS) % nx;
								int y = (by + ky*LocalBS) % ny;
								int z = (bz + kz*LocalBS) % nz;
								val = (Band[RI])(x,y,z) ;
								sigma2+=pow(val,2);
							}

					float sig2 = max( 0.0, sigma2/cnt - noise2 );
					float norm = sig2 / (sig2+noise2);


				// Store the coef in the table
					for(int bx = 0; bx < LocalBS; bx++)
						for(int by = 0; by < LocalBS; by++)
							for(int bz = 0; bz < LocalBS; bz++)
							{
								int x = (bx + kx*LocalBS) % nx;
								int y = (by + ky*LocalBS) % ny;
								int z = (bz + kz*LocalBS) % nz;
								if( coef_wiener(x,y,z) < -1 )
									coef_wiener(x,y,z) = norm;
							}

	//cerr<<kx<<","<<ky<<","<<kz<<", norm = "<<norm;
				// Apply the coefficients
					for(int bx = 0; bx < LocalBS; bx++)
						for(int by = 0; by < LocalBS; by++)
							for(int bz = 0; bz < LocalBS; bz++)
							{
								int x = (bx + kx*LocalBS) % nx;
								int y = (by + ky*LocalBS) % ny;
								int z = (bz + kz*LocalBS) % nz;
								(Band[RI])(x,y,z) *= coef_wiener(x,y,z);
							}
				}
	}

	if(Verbose) cerr<<"...End FCurvelet3D::wiener_single"<<endl;
}

// *********************************************************************

void FCurvelet3D::fdr(fltarray ** &TabBand, float Alpha, float SigmaNoise)
{
	if(Verbose) cerr<<"RCurvelet3D::fdr(.,"<<Alpha<<","<<SigmaNoise<<")"<<endl;
	bool LocVerbose = false & Verbose;
	
	for(int s3=0;s3<NbrScale-1;s3++)
		for(int b=0;b<nbr_band(s3);b++)
		{
		//getband
			fltarray Band;
			get_band(s3, b, TabBand, Band);

			int nx = Band.nx();
			int ny = Band.ny();
			int nz = Band.nz();
			
			float lvl = get_noise_lvl(SigmaNoise,s3);
			
		// P-Values
			fltarray PVal(Band.nx(),Band.ny(),Band.nz());
			for (int i=0; i < nx; i++)
				for (int j=0; j < ny; j++)
					for (int k=0; k < nz; k++)
						PVal(i,j,k) = prob_noise(Band(i,j,k)/lvl);
	
			double PDet = fdr_pvalue(PVal.buffer(), PVal.n_elem(), Alpha);
			
			if(LocVerbose) cerr<<" band ("<<s3<<","<<b<<") : PDet="<<PDet<<endl;

		// threshold
			for (int i=0; i < nx; i++)
				for (int j=0; j < ny; j++)
					for (int k=0; k < nz; k++)
						if( PVal(i,j,k) > PDet) (TabBand[s3][b])(i,j,k) = 0;
		}
	
	if(Verbose) cerr<<"...End RCurvelet3D::fdr"<<endl;
}

// **************************************************************************

void FCurvelet3D::stein_block_threshold(fltarray ** &TabBand, float SigmaNoise)
{
	if(Verbose) cerr<<"FCurvelet3D::stein_block_threshold(.)"<<endl;
	//bool LocVerbose = true & Verbose;
	float val;
	
	bool use_theorical_scales = true;
	
	int d = 3; 						// number of dimensions of the space
	float delta = 0.;				// parameter (0 for denoising)
	float lambda = 4.50524;			// parameter (root of x-log(x)=3)
	int n = max(NewNx,max(NewNy,NewNz));	// size of the data
	float r = d;						// ?? in [1,d]
	int L = (int) floor( pow(r*log(n),1./float(d)) );
	float cnt=pow((float)L,3);
	float mu1=1, mu2=0.5, mu3=0.5;	// anisotropy, in R+
	float dstar = mu1+mu2+mu3;		//
	float nu = 1;					// ?? in [0,1] (dstar+nu=d in general)
	int J = (int) floor(log2(n)); 		// maximum number of scales
	int j0 = (int) floor( 1/min(mu1,min(mu2,mu3)) * log2(L) );	// coarse scale limit
	int Jstar = (int) floor( r/(dstar+delta+nu)*log2(n) );		// scale limit
	
	// change of reference : now the coarsest is NbrScale-1, and 0 si the finest
	// 0 -> J ; j0 -> J-j0 ; Jstar -> J-Jstar
	// from j0 to Jstar -> from J-Jstar to J-j0
	if(!use_theorical_scales)
	{
		Jstar=J; 
		j0=J+2-NbrScale;
	}
L=3;
	for(int s=max(0,J-Jstar);s<=min(NbrScale-1,J-j0);s++)
	{
//		float factor = lambda*pow((float)2,(float)delta)*pow(get_noise_lvl(SigmaNoise,s),(float)2)/pow((float)n,r) ;
//		float factor = lambda*pow((float)2,(float)delta)*pow(get_noise_lvl(SigmaNoise,s),(float)2)/pow((float)L,d) ;
		float factor = lambda*get_noise_lvl(SigmaNoise,s)*get_noise_lvl(SigmaNoise,s) ;

cerr<<"scales "<<max(0,J-Jstar)<<" to "<<min(NbrScale-1,J-j0)<<", BlockSize = "<<L<<", factor="<<factor<<endl;
		for(int b=0;b<nbr_band(s);b++)
		{
			int nx = TabBand[s][b].nx();
			int ny = TabBand[s][b].ny();
			int nz = TabBand[s][b].nz();
			int Nx = ceil(float(nx)/float(L));
			int Ny = ceil(float(ny)/float(L));
			int Nz = ceil(float(nz)/float(L));
			
			fltarray coef_SB(nx,ny,nz);
			coef_SB.init(-2);
			
			// Local blocks 
			for(int kx=0 ; kx < Nx ; kx++)
				for(int ky=0 ; ky < Ny ; ky++)
					for(int kz = 0; kz < Nz; kz++)
					{
						double sigma2 = 0.0;

					// Threshold calculation
						// Pixels in a wiener block (angle)
						for(int bx = 0; bx < L; bx++)
							for(int by = 0; by < L; by++)
								for(int bz = 0; bz < L; bz++)
								{
									int x = (bx + kx*L) % nx;
									int y = (by + ky*L) % ny;
									int z = (bz + kz*L) % nz;
									val = (TabBand[s][b])(x,y,z);
									sigma2+=val*val;
								}
						sigma2 /= cnt;
						//float norm = max( 0.,1-factor/sigma2);
						float norm = max( 0.,1-factor/sigma2);

					// Store the coef in the table
						for(int bx = 0; bx < L; bx++)
							for(int by = 0; by < L; by++)
								for(int bz = 0; bz < L; bz++)
								{
									int x = (bx + kx*L) % nx;
									int y = (by + ky*L) % ny;
									int z = (bz + kz*L) % nz;
									if( coef_SB(x,y,z) < -1 )
										coef_SB(x,y,z) = norm;
								}

					// Apply the coefficients
						for(int bx = 0; bx < L; bx++)
							for(int by = 0; by < L; by++)
								for(int bz = 0; bz < L; bz++)
								{
									int x = (bx + kx*L) % nx;
									int y = (by + ky*L) % ny;
									int z = (bz + kz*L) % nz;
									(TabBand[s][b])(x,y,z) *= coef_SB(x,y,z);
								}
					}
		}
	}
	
// coarse scale : N-1
	// Do not denoise the coarse scale
	
	if(Verbose) cerr<<"...End FCurvelet3D::stein_block_threshold"<<endl;
}

// ********************************************************************

void FCurvelet3D::fdr_single(fltarray * Band, int s, int b, float Alpha, float SigmaNoise)
{ //!!!!! only applicable for real transform
	if(Verbose) cerr<<"RCurvelet3D::fdr(.,"<<Alpha<<","<<SigmaNoise<<")"<<endl;
	bool LocVerbose = false & Verbose;
	
	int nx = Band[0].nx();
	int ny = Band[0].ny();
	int nz = Band[0].nz();

	float lvl = get_noise_lvl(SigmaNoise,s);

// REAL and IMAG parts
	for(int RI=0;RI<2;RI++)
	{
	// P-Values
		fltarray PVal(Band[0].nx(),Band[0].ny(),Band[0].nz());
		for (int i=0; i < nx; i++)
			for (int j=0; j < ny; j++)
				for (int k=0; k < nz; k++)
					PVal(i,j,k) = prob_noise((Band[RI])(i,j,k)/lvl);

		double PDet = fdr_pvalue(PVal.buffer(), PVal.n_elem(), Alpha);

		if(LocVerbose) cerr<<" band ("<<s<<","<<b<<"R) : PDet="<<PDet<<endl;

	// threshold
		for (int i=0; i < nx; i++)
			for (int j=0; j < ny; j++)
				for (int k=0; k < nz; k++)
					if( PVal(i,j,k) > PDet) (Band[RI])(i,j,k) = 0;
	}
	
	if(Verbose) cerr<<"...End RCurvelet3D::fdr"<<endl;
}

// ***************************************************************************

void FCurvelet3D::stein_block_threshold_single(fltarray * Band, int s, int b, float SigmaNoise)
{ //!!!!! only applicable for real transform
	if(Verbose) cerr<<"FCurvelet3D::stein_block_threshold(.)"<<endl;
	bool LocVerbose = false & Verbose;
	float val;
	
	bool use_theorical_scales = true;
	
	int d = 3; 								// number of dimensions of the space
	float delta = 0.;						// parameter (0 for denoising)
	float lambda = 4.50524;					// parameter (root of x-log(x)=3)
	int n = max(NewNx,max(NewNy,NewNz));	// size of the data
	float r = d;							// ?? in [1,d]
	int L = (int) floor( pow(r*log(n),1./float(d)) );
	float cnt=pow((float)L,3);
	float mu1=1, mu2=0.5, mu3=0.5;			// anisotropy, in R+
	float dstar = mu1+mu2+mu3;				//
	float nu = 1;							// ?? in [0,1] (dstar+nu=d in general)
	int J = (int) floor(log2(n)); 			// maximum number of scales
	int j0 = (int) floor( 1/min(mu1,min(mu2,mu3)) * log2(L) );	// coarse scale limit
	int Jstar = (int) floor( r/(dstar+delta+nu)*log2(n) );		// scale limit
	
	float factor = lambda*get_noise_lvl(SigmaNoise,s)*get_noise_lvl(SigmaNoise,s);

	// change of reference : now the coarsest is NbrScale-1, and 0 si the finest
	// 0 -> J ; j0 -> J-j0 ; Jstar -> J-Jstar
	// from j0 to Jstar -> from J-Jstar to J-j0
	if(!use_theorical_scales)
	{
		Jstar=J; 
		j0=J+2-NbrScale;
	}
	if(LocVerbose) cerr<<"scales "<<max(0,J-Jstar)<<" to "<<min(NbrScale-1,J-j0)<<", BlockSize = "<<L<<", factor="<<factor<<endl;
	
	if(s>max(0,J-Jstar) && s<=min(NbrScale-1,J-j0))
	{
		for(int RI=0;RI<2;RI++)
		{
			int nx = Band[RI].nx();
			int ny = Band[RI].ny();
			int nz = Band[RI].nz();
			int Nx = ceil(float(nx)/float(L));
			int Ny = ceil(float(ny)/float(L));
			int Nz = ceil(float(nz)/float(L));

			fltarray coef_SB(nx,ny,nz);
			coef_SB.init(-2);

			// Local blocks 
			for(int kx=0 ; kx < Nx ; kx++)
				for(int ky=0 ; ky < Ny ; ky++)
					for(int kz = 0; kz < Nz; kz++)
					{
						double sigma2 = 0.0;

					// Threshold calculation
						// Pixels in a wiener block (angle)
						for(int bx = 0; bx < L; bx++)
							for(int by = 0; by < L; by++)
								for(int bz = 0; bz < L; bz++)
								{
									int x = (bx + kx*L) % nx;
									int y = (by + ky*L) % ny;
									int z = (bz + kz*L) % nz;
									val = (Band[RI])(x,y,z) ; 
									sigma2+=val*val;
								}
						sigma2 /= cnt;
						float norm = max( 0.,1-factor/sigma2);

					// Store the coef in the table
						for(int bx = 0; bx < L; bx++)
							for(int by = 0; by < L; by++)
								for(int bz = 0; bz < L; bz++)
								{
									int x = (bx + kx*L) % nx;
									int y = (by + ky*L) % ny;
									int z = (bz + kz*L) % nz;
									if( coef_SB(x,y,z) < -1 )
										coef_SB(x,y,z) = norm;
								}

					// Apply the coefficients
						for(int bx = 0; bx < L; bx++)
							for(int by = 0; by < L; by++)
								for(int bz = 0; bz < L; bz++)
								{
									int x = (bx + kx*L) % nx;
									int y = (by + ky*L) % ny;
									int z = (bz + kz*L) % nz;
									(Band[RI])(x,y,z) *= coef_SB(x,y,z);
								}
					}
		}
	}
	
// coarse scale : N-1
	// Do not denoise the coarse scale
	
	if(Verbose) cerr<<"...End FCurvelet3D::stein_block_threshold"<<endl;
}

// *********************************************************************

static inline void PrintError( int status)
{
    // ***************************************************** 
    // * Print out cfitsio error messages and exit program * 
    // ***************************************************** 

    char status_str[FLEN_STATUS], errmsg[FLEN_ERRMSG];

    if (status)
      fprintf(stderr, "\n*** Error occurred during program execution ***\n");

    ffgerr(status, status_str);        // get the error status description 
    fprintf(stderr, "\nstatus = %d: %s\n", status, status_str);

    if ( ffgmsg(errmsg) )  // get first message; null if stack is empty 
    {
         fprintf(stderr, "\nError message stack:\n");
         fprintf(stderr, " %s\n", errmsg);

         while ( ffgmsg(errmsg) )  // get remaining messages 
             fprintf(stderr, " %s\n", errmsg);
    }

    exit( status );
}

// *********************************************************************

static inline void mr_io_name (char *File_Name_In, char *File_Name_Out)
{
    int L;

    strcpy (File_Name_Out, File_Name_In);

    L = strlen (File_Name_In);
    if ((L < 3) || (File_Name_In[L-1] != 'r')
                || (File_Name_In[L-2] != 'm')
                || (File_Name_In[L-3] != '.'))
    {
        strcat (File_Name_Out, ".mr");
    }
}

// *********************************************************************

void FCurvelet3D::write(char *Name, fltarray ** TabBand)
{
	if(Verbose) cerr<<"FCurvelet3D::write("<<Name<<",.,"<<")"<<endl;
	
	char filename[256];
	fitsfile *fptr;    
	int status;
	int simple;
	int bitpix;
	long naxis=0;
	long naxes[3];
	long group = 1; 

	// .mr extention
	mr_io_name (Name, filename);

	FILE *FEXIST = fopen(filename, "rb");
	if (FEXIST)
	{
		fclose(FEXIST);
		remove(filename);               // Delete old file if it already exists 
	}
	status = 0;         // initialize status before calling fitsio routines 

// open the file
	if ( ffinit(&fptr, filename, &status) )	// create the new FITS file 
		PrintError( status );					// call PrintError if error occurs 
	
// write  the header 
	simple   = True;
	bitpix   =  -32;   // 32-bit real pixel values      
	long pcount   =   0;  // no group parameters 
	long gcount   =   1;  // only a single image/group 
	int  extend   =   False;

// write first header part (parameters)
	naxis=0;
	if (ffphps(fptr, bitpix, naxis, naxes, &status))
		PrintError( status );  
	// write optional keyword to the header 
		if ( ffpkyj(fptr, (char*)"Type_Tra", (long) 0, (char*)"3D FastCurvelet", &status))
			PrintError( status );  
		if ( ffpkyj(fptr, (char*)"NbrScale", (long) NbrScale, (char*)"Number of 3D scales", &status))
			PrintError( status );  
		if ( ffpkyj(fptr, (char*)"NbrDir2d", (long) NbrDir2d, (char*)"Number of 2D directions (on a circle)", &status))
			PrintError( status );  
		if ( ffpkyj(fptr, (char*)"NbrBands", (long) nbr_tot_band(), (char*)"Total number of bands", &status))
			PrintError( status );  
		if ( ffpkyj(fptr, (char*)"ModifSiz", (long) ModifSize, (char*)"1 if size is modified (ie. was not odd)", &status))
			PrintError( status );  
		if ( ffpkyj(fptr, (char*)"ExtendWT", (long) Extend, (char*)"Number of bands 3D", &status))
			PrintError( status );  
		if ( ffpkyj(fptr, (char*)"IsotropW", (long) Isotrop, (char*)"Number of bands 3D", &status))
			PrintError( status );  
		if ( ffpkyj(fptr, (char*)"RealCur", (long) RealBand, (char*)"Number of bands 3D", &status))
			PrintError( status );  
		if ( ffpkyj(fptr, (char*)"Nx_Cube", (long) DataNx, (char*)"x size of the original cube", &status))
			PrintError( status );  
		if ( ffpkyj(fptr, (char*)"Ny_Cube", (long) DataNy, (char*)"y size of the original cube", &status))
			PrintError( status );  
		if ( ffpkyj(fptr, (char*)"Nz_Cube", (long) DataNz, (char*)"z size of the original cube", &status))
			PrintError( status );  
	
// write other headers and associated data	
// Fine scales
	for (int s=0; s < NbrScale; s++)
	{
		for (int b=0; b < nbr_angle(s); b++)
		{
			naxis=3;
			naxes[0] = TabBand[s][b].nx();
			naxes[1] = TabBand[s][b].ny();
			naxes[2] = TabBand[s][b].nz();
//	cerr<<"save : "<<s<<" "<<naxes[0]<<" "<<naxes[1]<<" "<<naxes[2]<<endl;
			if(ffcrhd(fptr,&status))
				PrintError( status );
			if ( ffphpr(fptr,simple,bitpix,naxis,naxes,pcount,gcount,extend,&status) )
				PrintError( status );

		// save the data
			if ( ffppre(fptr, group, 1, naxes[0]*naxes[1]*naxes[2], (TabBand[s][b]).buffer(), &status) )
				PrintError( status );
		}
	}
	
// close the FITS file 
	if ( ffclos(fptr, &status) )  PrintError( status );  
	
	if(Verbose) cerr<<"...end FCurvelet3D::write"<<endl;
}

// ********************************************************************

void FCurvelet3D::read(char *Name, fltarray ** &TabBand)
{
	if(Verbose) cerr<<"BCurvelet3D::read("<<Name<<",.,.,.)"<<endl;
	bool LocVerbose = false & Verbose;
	
	char filename[256];
	fitsfile *fptr;           // pointer to the FITS file 
	int status=0, hdutype ;
	char comment[FLEN_COMMENT];
	long mon_long;
	int anynul = 0;
	long nulval = 0;
	void PrintError( int status);

	mr_io_name (Name, filename);

// open the file 
	status = 0;         // initialize status before calling fitsio routines 
	if ( ffopen(&fptr, filename, (int) READONLY, &status) ) 
		PrintError( status );

// get number of Headers
	int nhead;
	fits_get_num_hdus(fptr, &nhead, &status);
	
// read primary header
	if ( ffmahd(fptr, 1, &hdutype, &status) ) PrintError( status );
	
	// Read params
	int _NbrScale, _Nx, _Ny, _Nz, _NbrDir2d;
	Bool _ExtendWT, _IsotropWT, _RealCur;
	if (ffgkyj(fptr,(char*)"Type_Tra", &mon_long, comment, &status)) PrintError( status );
	if (ffgkyj(fptr,(char*)"NbrScale", &mon_long, comment, &status)) PrintError( status );
	_NbrScale = (int)mon_long;
	if (ffgkyj(fptr,(char*)"NbrDir2d", &mon_long, comment, &status)) PrintError( status );
	_NbrDir2d = (int)mon_long;
	if (ffgkyj(fptr,(char*)"ExtendWT", &mon_long, comment, &status)) PrintError( status );
	_ExtendWT = (Bool)mon_long;
	if (ffgkyj(fptr,(char*)"IsotropW", &mon_long, comment, &status)) PrintError( status );
	_IsotropWT = (Bool)mon_long;
	if (ffgkyj(fptr,(char*)"RealCur", &mon_long, comment, &status)) PrintError( status );
	_RealCur = (Bool)mon_long;
	
// Check the number of bands
//	if(nhead!=_NbrScale+1) { cerr<<"Wrong header number in fits file : Hdr="<<nhead<<", NbrScale="<<_NbrScale<<endl; exit(0); }
	
	if (ffgkyj(fptr,(char*)"Nx_Cube", &mon_long, comment, &status)) PrintError( status );
	_Nx = (int)mon_long;
	if (ffgkyj(fptr,(char*)"Ny_Cube", &mon_long, comment, &status)) PrintError( status );
	_Ny = (int)mon_long;
	if (ffgkyj(fptr,(char*)"Nz_Cube", &mon_long, comment, &status)) PrintError( status );
	_Nz = (int)mon_long;
	
//	init the structure
	alloc_from_coarse(_NbrScale,_Nx, _Ny, _Nz, _NbrDir2d, _ExtendWT, _IsotropWT, _RealCur);
	
// TabBand allocation
	TabBand = new fltarray*[NbrScale];
	for(int s=0;s<nbr_scale()-1; s++)
	{
		TabBand[s] = new fltarray[nbr_angle(s)];
	}
	TabBand[nbr_scale()-1] = new fltarray[1];
	
// Read data
	int cnt=1;
	for(int s=0;s<_NbrScale;s++)
	{
		for (int b=0; b < nbr_angle(s); b++)
		{
			if(LocVerbose) cerr<<" read : "<<s<<" "<<b <<endl;
			int NX,NY,NZ;
			if (fits_movabs_hdu(fptr, ++cnt, NULL, &status)) PrintError( status );

			if (ffgkyj(fptr,(char*)"NAXIS1", &mon_long, comment, &status)) PrintError( status );
			NX = (int)mon_long;
			if (ffgkyj(fptr,(char*)"NAXIS2", &mon_long, comment, &status)) PrintError( status );
			NY = (int)mon_long;
			if (ffgkyj(fptr,(char*)"NAXIS3", &mon_long, comment, &status)) PrintError( status );
			NZ = (int)mon_long;

			if(LocVerbose) cerr<<" read TB("<<s<<","<<b<<") : "<<NX<<" "<<NY<<" "<<NZ<<endl;
			
			TabBand[s][b].alloc(NX,NY,NZ);
			if (ffgpve(fptr, 1, 1, NX*NY*NZ, nulval, (TabBand[s][b]).buffer(), &anynul, &status)) PrintError( status );
		}
	}
	
	
// close the FITS file 
	if ( ffclos(fptr, &status) ) PrintError( status );
	
	if(Verbose) cerr<<"...end BCurvelet3D::read"<<endl;
}

// ********************************************************************

void FCurvelet3D::write_nohead(char *Name, fltarray ** TabBand)
{
	if(Verbose) cerr<<"FCurvelet3D::write("<<Name<<",.,"<<")"<<endl;
	
	char filename[256];
	fitsfile *fptr;    
	int status;
	int simple;
	int bitpix;
	long naxis=0;
	long naxes[3];
	long group = 1; 

	// .mr extention
	mr_io_name (Name, filename);

	FILE *FEXIST = fopen(filename, "rb");
	if (FEXIST)
	{
		fclose(FEXIST);
		remove(filename);               // Delete old file if it already exists 
	}
	status = 0;         // initialize status before calling fitsio routines 

// open the file
	if ( ffinit(&fptr, filename, &status) )	// create the new FITS file 
		PrintError( status );					// call PrintError if error occurs 
	
// write  the header 
	simple   = True;
	bitpix   =  -32;   // 32-bit real pixel values      
//	long pcount   =   0;  // no group parameters 
//	long gcount   =   1;  // only a single image/group 
//	int  extend   =   False;

// write first header part (parameters)
	naxis=0;
int totalsize=0;
for (int s=0; s < NbrScale; s++)
for (int b=0; b < nbr_angle(s); b++)
totalsize += TabBand[s][b].nx()*TabBand[s][b].ny()*TabBand[s][b].nz() ;
naxis=1;
naxes[0]=totalsize;
naxes[1]=0;
naxes[2]=0;
	if (ffphps(fptr, bitpix, naxis, naxes, &status))
		PrintError( status );  
	// write optional keyword to the header 
		if ( ffpkyj(fptr, (char*)"Type_Tra", (long) 0, (char*)"3D FastCurvelet", &status))
			PrintError( status );  
		if ( ffpkyj(fptr, (char*)"NbrScale", (long) NbrScale, (char*)"Number of 3D scales", &status))
			PrintError( status );  
		if ( ffpkyj(fptr, (char*)"NbrDir2d", (long) NbrDir2d, (char*)"Number of 2D directions (on a circle)", &status))
			PrintError( status );  
		if ( ffpkyj(fptr, (char*)"NbrBands", (long) nbr_tot_band(), (char*)"Total number of bands", &status))
			PrintError( status );  
		if ( ffpkyj(fptr, (char*)"ModifSiz", (long) ModifSize, (char*)"1 if size is modified (ie. was not odd)", &status))
			PrintError( status );  
		if ( ffpkyj(fptr, (char*)"ExtendWT", (long) Extend, (char*)"Number of bands 3D", &status))
			PrintError( status );  
		if ( ffpkyj(fptr, (char*)"IsotropW", (long) Isotrop, (char*)"Number of bands 3D", &status))
			PrintError( status );  
		if ( ffpkyj(fptr, (char*)"RealCur", (long) RealBand, (char*)"Number of bands 3D", &status))
			PrintError( status );  
		if ( ffpkyj(fptr, (char*)"Nx_Cube", (long) DataNx, (char*)"x size of the original cube", &status))
			PrintError( status );  
		if ( ffpkyj(fptr, (char*)"Ny_Cube", (long) DataNy, (char*)"y size of the original cube", &status))
			PrintError( status );  
		if ( ffpkyj(fptr, (char*)"Nz_Cube", (long) DataNz, (char*)"z size of the original cube", &status))
			PrintError( status );  
	
// write other headers and associated data	
// Fine scales
	int cnt=1;
	for (int s=0; s < NbrScale; s++)
	{
		for (int b=0; b < nbr_angle(s); b++)
		{
			naxis=3;
			naxes[0] = TabBand[s][b].nx();
			naxes[1] = TabBand[s][b].ny();
			naxes[2] = TabBand[s][b].nz();
//	cerr<<"save ("<<s<<","<<b<<") "<<naxes[0]<<" "<<naxes[1]<<" "<<naxes[2]<<endl;

		// save the data
			if ( ffppre(fptr, group, cnt, naxes[0]*naxes[1]*naxes[2], (TabBand[s][b]).buffer(), &status) )
				PrintError( status );
			cnt += naxes[0]*naxes[1]*naxes[2];
//			if ( fits_write_subset(fptr, group, ) )
//				PrintError( status );
			
//	cerr<<"saved"<<endl;
		}
	}
	
// close the FITS file 
	if ( ffclos(fptr, &status) )  PrintError( status );  
	
	if(Verbose) cerr<<"...end FCurvelet3D::write"<<endl;
}

// ********************************************************************

void FCurvelet3D::read_nohead(char *Name, fltarray ** &TabBand)
{
	if(Verbose) cerr<<"BCurvelet3D::read("<<Name<<",.,.,.)"<<endl;
	bool LocVerbose = false & Verbose;
	
	char filename[256];
	fitsfile *fptr;           // pointer to the FITS file 
	int status=0, hdutype ;
	char comment[FLEN_COMMENT];
	long mon_long;
	int anynul = 0;
	long nulval = 0;
	void PrintError( int status);

	mr_io_name (Name, filename);

// open the file 
	status = 0;         // initialize status before calling fitsio routines 
	if ( ffopen(&fptr, filename, (int) READONLY, &status) ) 
		PrintError( status );

// get number of Headers
	int nhead;
	fits_get_num_hdus(fptr, &nhead, &status);
	
// read primary header
	if ( ffmahd(fptr, 1, &hdutype, &status) ) PrintError( status );
	
	// Read params
	int _NbrScale, _Nx, _Ny, _Nz, _NbrDir2d;
	Bool _ExtendWT, _IsotropWT, _RealCur;
	if (ffgkyj(fptr,(char*)"Type_Tra", &mon_long, comment, &status)) PrintError( status );
	if (ffgkyj(fptr,(char*)"NbrScale", &mon_long, comment, &status)) PrintError( status );
	_NbrScale = (int)mon_long;
	if (ffgkyj(fptr,(char*)"NbrDir2d", &mon_long, comment, &status)) PrintError( status );
	_NbrDir2d = (int)mon_long;
	if (ffgkyj(fptr,(char*)"ExtendWT", &mon_long, comment, &status)) PrintError( status );
	_ExtendWT = (Bool)mon_long;
	if (ffgkyj(fptr,(char*)"IsotropW", &mon_long, comment, &status)) PrintError( status );
	_IsotropWT = (Bool)mon_long;
	if (ffgkyj(fptr,(char*)"RealCur", &mon_long, comment, &status)) PrintError( status );
	_RealCur = (Bool)mon_long;
	
// Check the number of bands
//	if(nhead!=_NbrScale+1) { cerr<<"Wrong header number in fits file : Hdr="<<nhead<<", NbrScale="<<_NbrScale<<endl; exit(0); }
	
	if (ffgkyj(fptr,(char*)"Nx_Cube", &mon_long, comment, &status)) PrintError( status );
	_Nx = (int)mon_long;
	if (ffgkyj(fptr,(char*)"Ny_Cube", &mon_long, comment, &status)) PrintError( status );
	_Ny = (int)mon_long;
	if (ffgkyj(fptr,(char*)"Nz_Cube", &mon_long, comment, &status)) PrintError( status );
	_Nz = (int)mon_long;
	
//	init the structure
	alloc_from_coarse(_NbrScale,_Nx, _Ny, _Nz, _NbrDir2d, _ExtendWT, _IsotropWT, _RealCur);
	
// TabBand allocation
	TabBand = new fltarray*[NbrScale];
	for(int s=0;s<nbr_scale()-1; s++)
	{
		TabBand[s] = new fltarray[nbr_angle(s)];
	}
	TabBand[nbr_scale()-1] = new fltarray[1];
	
// Read data
	int cnt=1;
	for(int s=0;s<_NbrScale;s++)
	{
		for (int b=0; b < nbr_angle(s); b++)
		{
			int NX,NY,NZ;
			NX = TabSizeNx[s](b);
			NY = TabSizeNy[s](b);
			NZ = TabSizeNz[s](b);
			if(LocVerbose) cerr<<" read TB("<<s<<","<<b<<") : "<<NX<<" "<<NY<<" "<<NZ<<endl;
			
			TabBand[s][b].alloc(NX,NY,NZ);
			if (ffgpve(fptr, 1, cnt, NX*NY*NZ, nulval, (TabBand[s][b]).buffer(), &anynul, &status)) PrintError( status );
			cnt += NX*NY*NZ;
		}
	}
	
	
// close the FITS file 
	if ( ffclos(fptr, &status) ) PrintError( status );
	
	if(Verbose) cerr<<"...end BCurvelet3D::read"<<endl;}

// ********************************************************************

void FCurvelet3D::redundancy()
{
	double tot=DataNx*DataNy*DataNz;
	double ftot=0;
	
	for(int s=0;s<NbrScale;s++)
		for (int b=0; b < nbr_angle(s); b++)
			ftot += size_band_nx(s,b)*size_band_ny(s,b)*size_band_nz(s,b);
	
	if(Verbose) cerr<<"Redundancy factor = "<<ftot/tot<<endl;
}

// ********************************************************************

void FCurvelet3D::temp(fltarray **TabBand)
{

// Set the transformto 0
	for(int s=0;s<NbrScale;s++)
		for(int b=0;b<nbr_band(s);b++)
			for(int i=0;i<size_band_nx(s,b);i++)
			for(int j=0;j<size_band_ny(s,b);j++)
			for(int k=0;k<size_band_nz(s,b);k++)
				(TabBand[s][b])(i,j,k)=0;
	
	
	
// Put a coefficient
	(TabBand[1][0])(0,0,0)=1;
	(TabBand[1][1])(35,15,15)=1;
//	(TabBand[0][0])(25,16,16)=1;

}



#include "dirdct3d_single.h"

extern Bool Verbose;

DirDCT3D::DirDCT3D()
{
	
}

DirDCT3D::~DirDCT3D()
{
	
}




void DirDCT3D::init(int BS)
{
	BlockSize = BS;
}

void DirDCT3D::alloc()
{
	TabSize.alloc(BlockSize,BlockSize,BlockSize);
}

extern void cosft2(float y[], int n, int isign);


void DirDCT3D::transform(fltarray &Cube, fltarray &Trans)
{
	if(Verbose) cerr<<"DirDCT3D::transform(.,.)..."<<endl;
	
	PRad3D.setTakeOnlyCompletePlane(False);
	PRad3D.setDoNotCorrectSigma(True);   
	PRad3D.setParam(BlockSize,0);  	   // unused second parameter
	PRad3D.alloc();
	
	int Nplane = PRad3D.getNbPlane();
	
	
//cerr<<"TRANS "<<Nplane<<endl;

	Trans.alloc(BlockSize,BlockSize,Nplane);
//cerr<<Cube.nx()<<endl;
//cerr<<Trans.nx()<<endl;
	PRad3D.transform (Cube, Trans);
//cerr<<"TRANS"<<endl;
	


// no odd-size dct for now
BlockSize-=1;



//char filename[64];
//sprintf(filename,"%s_trans.fits","out");
//writefltarr(filename, Trans);

/*
	Cube.init(1.);
	// Slice of the transformed cube
	float *Ptr = Trans.buffer();
	Ifloat F,Ft;
	//F.alloc(Ptr,Trans.ny(),Trans.nx());
	F.alloc(Trans.ny(),Trans.nx());
	for(int i=0;i<Trans.nx();i++)
		for(int j=0;j<Trans.ny();j++)
			F(j,i)=Trans(i,j,0);
	Ft.alloc(Trans.ny(),Trans.nx());
//	Ft=F;
	
	im_dct(F,Ft);
	
	for(int i=0;i<Trans.nx();i++)
		for(int j=0;j<Trans.ny();j++)
			Trans(i,j,0)=Ft(j,i);
	
sprintf(filename,"%s_transT.fits","out");
writefltarr(filename, Trans);
	
//	Ft.init(2.);
	F.init(2.);
	im_dct(Ft,F,True);
	
	for(int i=0;i<Trans.nx();i++)
		for(int j=0;j<Trans.ny();j++)
			Trans(i,j,0)=F(j,i);

sprintf(filename,"%s_transR.fits","out");
writefltarr(filename, Trans);
*/
	float *d = new float[BlockSize];
	
	// Transform on the three dimensions
	for(int j=0;j<Nplane;j++)
	{
		for(int k=0;k<BlockSize;k++)
		{
			for(int i=0;i<BlockSize;i++)
				d[i]=Trans(i,k,j);
			cosft2(d-1, BlockSize, 1);
			for(int i=0;i<BlockSize;i++)
				Trans(i,k,j)=d[i];
		}
		for(int k=0;k<BlockSize;k++)
		{
			for(int i=0;i<BlockSize;i++)
				d[i]=Trans(k,i,j);
			cosft2(d-1, BlockSize, 1);
			for(int i=0;i<BlockSize;i++)
				Trans(k,i,j)=d[i];
		}
	}

//sprintf(filename,"%s_transT.fits","out");
//writefltarr(filename, Trans);
	
// no odd-size dct for now
BlockSize++;

	if(Verbose) cerr<<"end DirDCT3D::transform(.,.)"<<endl;
}

void DirDCT3D::recons(fltarray &Trans, fltarray &Recons)
{
	if(Verbose) cerr<<"DirDCT3D::recons(.,.)..."<<endl;
	
	float *d = new float[BlockSize];
	Recons.alloc(BlockSize,BlockSize,BlockSize);
//	cerr<<BlockSize<<","<<Recons.nx()<<","<<Recons.ny()<<","<<Recons.nz()<<endl;


//char filename[64];

// no odd-size dct for now
BlockSize--;

	int Nplane = PRad3D.getNbPlane();
	for(int j=0;j<Nplane;j++)
	{
		for(int k=0;k<BlockSize;k++)
		{
			for(int i=0;i<BlockSize;i++)
				d[i]=Trans(i,k,j);
			cosft2(d-1, BlockSize, -1);
			for(int i=0;i<BlockSize;i++)
				Trans(i,k,j)=d[i];
		}
		for(int k=0;k<BlockSize;k++)
		{
			for(int i=0;i<BlockSize;i++)
				d[i]=Trans(k,i,j);
			cosft2(d-1, BlockSize, -1);
			for(int i=0;i<BlockSize;i++)
				Trans(k,i,j)=d[i];
		}
	}
// Normalisation
	for(int k=0;k<Nplane;k++)
		for(int j=0;j<BlockSize;j++)
			for(int i=0;i<BlockSize;i++)
				Trans(i,j,k) /= (BlockSize*BlockSize/4);

// no odd-size dct for now
BlockSize++;

//sprintf(filename,"%s_transR.fits","out");
//writefltarr(filename, Trans);
	
	

	PRad3D.recons (Trans, Recons);
	

//sprintf(filename,"%s_recons.fits","out");
//writefltarr(filename, Cube);

	if(Verbose) cerr<<"end DirDCT3D::recons(.,.)"<<endl;
}

void DirDCT3D::threshold(fltarray &Trans, float SigmaNoise, float NSigma, filter_type FilterType)
{
//	34.6 & 4.1
	float lvl;
	float cnt=0,cnt2=0;
	int Nplane = PRad3D.getNbPlane();
	for(int k=0;k<Nplane;k++)
		for(int j=0;j<BlockSize;j++)
			for(int i=0;i<BlockSize;i++)
			{
						lvl = SigmaNoise * NSigma * 34.6;

						cnt+=1;
						if( abs((Trans)(i,j,k)) < lvl )
						{// Hard and soft put these to 0
							cnt2+=1;
							(Trans)(i,j,k)=0;
						}// soft lowers the remaining by lvl
						else if(FilterType==FT_SOFT) (Trans)(i,j,k) -= (2*int( (Trans)(i,j,k)>0 )-1)*lvl;
			}
	if(Verbose) cerr<<"n#proportion non seuillee ("<<lvl<<")="<<cnt-cnt2<<"#"<<(cnt-cnt2)/cnt<<endl;
}

void DirDCT3D::write(char *Name, fltarray &Trans, bool Normalize)
{
	writefltarr(Name,Trans);
}

































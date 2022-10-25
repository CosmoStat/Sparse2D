/*******************************************************************************
**                   Copyright (C) 2003 by CEA
*******************************************************************************
**    UNIT
**
**    Version: 1.0
**
**    Author: Jean-Luc Starck
**
**    Date:  03/04/03 
**    
**    File:  LineCol.cc
**
**    Modification history:
**
*******************************************************************************
**
**    DESCRIPTION  Line Column Multiscale  decomposition
**    -----------  
**                 
******************************************************************************/
 

#include "IM_Obj.h"
#include "WT_Bord.h"
 #include "IM_IO.h"


#define maxf 7

static float g1d[2*maxf+1]={0.,0.,0.,0.,0.,0.,0.,-1./4.,1./2.,-1./4.,0.,0.,0.,0.,0.};
static float h1i[2*maxf+1]={0.,0.,0.,0.,0.,0.,1./2.,1.,1./2.,0.,0.,0.,0.,0.,0.};
static float h1d[2*maxf+1]={0.,0.,0.,0.,0.,-1./8.,1./4.,3./4.,1./4.,-1./8.,0.,0.,0.,0.,0.};
static float g1i[2*maxf+1]={0.,0.,0.,0.,0.,0.,-1./4.,-1./2.,3./2.,-1./2.,-1./4.,0.,0.,0.,0.};

static float g2d[2*maxf+1]={0.,0.,0.,0.,0.,0.,1./8.,-3./8.,3./8.,-1./8.,0.,0.,0.,0.,0.};
static float h2i[2*maxf+1]={0.,0.,0.,0.,0.,0.,1./4.,3./4.,3./4.,1./4.,0.,0.,0.,0.,0.};
static float h2d[2*maxf+1]={0.,0.,0.,0.,0.,0.,-1./4.,3./4.,3./4.,-1./4.,0.,0.,0.,0.,0.};
static float g2i[2*maxf+1]={0.,0.,0.,0.,0.,0.,-1./2.,-3./2.,3./2.,1./2.,0.,0.,0.,0.,0.};

/* Coiflets 12 */
static float gCoid[2*maxf+1]={0.,0.000720549445000,-0.001823208870999,-0.005611434818997,0.023680171946988,0.059434418645969,-0.076488599077960,-0.417005184423784,0.812723635449579,-0.386110066822800,-0.067372554721965,0.041464936781979,0.016387336462992,0.,0.};
static float hCoid[2*maxf+1]={0.,0.,0.,0.016387336462992,-0.041464936781979,-0.067372554721965,0.386110066822800,0.812723635449579,0.417005184423784,-0.076488599077960,-0.059434418645969,0.023680171946988,0.005611434818997,-0.001823208870999,-0.000720549445000};
static float hCoii[2*maxf+1]={0.,0.,0.,0.016387336462992,-0.041464936781979,-0.067372554721965,0.386110066822800,0.812723635449579,0.417005184423784,-0.076488599077960,-0.059434418645969,0.023680171946988,0.005611434818997,-0.001823208870999,-0.000720549445000};
static float gCoii[2*maxf+1]={0.,0.000720549445000,-0.001823208870999,-0.005611434818997,0.023680171946988,0.059434418645969,-0.076488599077960,-0.417005184423784,0.812723635449579,-0.386110066822800,-0.067372554721965,0.041464936781979,0.016387336462992,0.,0.};



void tfd1D_bl_V0(float *Input,float *Output, int Np) /* transformee avec bords libres avec splines lineaires */
{
   float *c = Input;
   float *d = Output;
   float g1[2*maxf+1]={0.,0.,0.,0.,0.,0.,0.,-1./4.,1./2.,-1./4.,0.,0.,0.,0.,0.};
   float h1[2*maxf+1]={0.,0.,0.,0.,0.,-1./8.,1./4.,3./4.,1./4.,-1./8.,0.,0.,0.,0.,0.};
   int N=Np-1;
   
    int j,i,l;
    float *k;
    k=(float*) malloc(sizeof(float)*(N+1));
    
    for (i=0;i<=N;i++) {
        k[i]=c[i];
    }
    for (j=N/2;j>1;j=j/2) {
        for (i=2;i<(j-1);i++) {
            d[i]=0.;
            for (l=-maxf;l<=maxf;l++) d[i]+=h1[l+maxf]*k[(2*i+l+maxf*(2*j+1))%(2*j+1)];
        }
        for (i=1;i<(j-1);i++) {
            d[j+1+i]=0.;
            for (l=-maxf;l<=maxf;l++) d[j+1+i]+=g1[l+maxf]*k[(2*i+l+maxf*(2*j+1))%(2*j+1)];
        }
        // exception pour les bords ds V0 !
        d[j+1]=-0.25*k[0]+0.5*k[1]-0.25*k[2];
        d[0]=k[0]+0.75*d[j+1];
        d[1]=0.25*k[1]+0.75*k[2]+0.25*k[3]-1./8.*k[4]-1./8.*d[0]-1./32.*d[j+1];
        
        d[2*j]=-0.25*k[2*j]+0.5*k[2*j-1]-0.25*k[2*j-2];
        d[j]=k[2*j]+0.75*d[2*j];
        d[j-1]=0.25*k[2*j-1]+0.75*k[2*j-2]+0.25*k[2*j-3]-1./8.*k[2*j-4]-1./8.*d[j]-1./32.*d[2*j];
        if (j==2) d[1]=0.25*k[2*j-1]+0.75*k[2*j-2]+0.25*k[2*j-3]-1./8.*d[0]-1./32.*d[j+1]-1./8.*d[j]-1./32.*d[2*j];
        
        for (i=0;i<=j;i++) k[i]=d[i];
    }
    d[0]=1./8.*k[0]+0.75*k[1]+1./8.*k[2];
    d[1]=7./12.*k[0]-0.5*k[1]-1./12.*k[2];
    d[2]=-1./12.*k[0]-0.5*k[1]+7./12.*k[2];
    free(k);
}

/****************************************************************************/

void tfd1Dinv_bl_V0(float *Input,float *Output, int Np) /* transformee avec bords libres avec splines quadratique */
{
    float *c = Output;
    float *d = Input;
    int N=Np-1;
    float h1[2*maxf+1]={0.,0.,0.,0.,0.,0.,1./2.,1.,1./2.,0.,0.,0.,0.,0.,0.};
    float g1[2*maxf+1]={0.,0.,0.,0.,0.,0.,-1./4.,-1./2.,3./2.,-1./2.,-1./4.,0.,0.,0.,0.};
    
    int j,i,l;
    float *cw;
    cw=(float*) malloc(sizeof(float)*(N+1));
    
    for (i=0;i<=N;i++) {
        cw[i]=d[i];
        c[i]=d[i];
    }
    c[0]=d[0]+1.5*d[1];
    c[1]=d[0]-0.25*d[1]-0.25*d[2];
    c[2]=d[0]+1.5*d[2];
    for (j=2;j<N;j=2*j) {
        for (i=0;i<=j;i++) cw[i]=c[i];
        for (i=0;i<=2*j;i++) c[i]=0.;
        for (i=1;i<j;i++) {
            for (l=-maxf;l<=maxf;l++) {
                c[(2*i+l+maxf*(2*j+1))%(2*j+1)]+=h1[l+maxf]*cw[i];
            }
        }
        for (i=1;i<(j-1);i++) {
            for (l=-maxf;l<=maxf;l++) {
                c[(2*i+l+maxf*(2*j+1))%(2*j+1)]+=g1[l+maxf]*d[j+1+i];
            }
        }
        c[0]+=-0.75*d[j+1]+cw[0];
        c[1]+=23./16.*d[j+1]+0.5*cw[0];
        c[2]+=-3./8.*d[j+1];
        c[3]+=-3./16.*d[j+1];
        c[2*j]+=-0.75*d[2*j]+cw[j];
        c[2*j-1]+=23./16.*d[2*j]+0.5*cw[j];
        c[2*j-2]+=-3./8.*d[2*j];
        c[2*j-3]+=-3./16.*d[2*j];
    }
    free(cw);
}


void tfd1D_bl_V1(float *Input,float *Output, int Np) /* transformee avec bords libres avec splines quadratique */
{
    float *c = Input;
    float *d = Output;
    int N=Np-2;
    float h1[2*maxf+1]={0.,0.,0.,0.,0.,0.,-1./4.,3./4.,3./4.,-1./4.,0.,0.,0.,0.,0.};
    float g1[2*maxf+1]={0.,0.,0.,0.,0.,0.,1./8.,-3./8.,3./8.,-1./8.,0.,0.,0.,0.,0.};
    int j,i,l;
    float *k;
    k=(float*) malloc(sizeof(float)*(N+2));
    
    for (i=0;i<=(N+1);i++) {
        k[i]=c[i];
    }
    for (j=N/2;j>1;j=j/2) {
        for (i=2;i<j;i++) {
            d[i]=0.;
            for (l=-maxf;l<=maxf;l++) d[i]+=h1[l+maxf]*k[(2*i-1+l+maxf*(2*j+2))%(2*j+2)];
        }
        for (i=1;i<(j-1);i++) {
            d[j+2+i]=0.;
            for (l=-maxf;l<=maxf;l++) d[j+2+i]+=g1[l+maxf]*k[(2*i+1+l+maxf*(2*j+2))%(2*j+2)];
        }
        // exception pour les bords ds V1 !
        d[j+2]=-0.25*k[0]+0.5*k[1]-3./8.*k[2]+1./8.*k[3];
        d[0]=k[0]-d[j+2];
        d[1]=1.5*k[2]-0.5*k[3]+1.5*d[j+2];
        
        d[2*j+1]=-0.25*k[2*j+1]+0.5*k[2*j]-3./8.*k[2*j-1]+1./8.*k[2*j-2];
        d[j+1]=k[2*j+1]-d[2*j+1];
        d[j]=1.5*k[2*j-1]-0.5*k[2*j-2]+1.5*d[2*j+1];
        
        for (i=0;i<=(j+1);i++) k[i]=d[i];
    }
    d[0]=1./6.*k[0]+1./3.*k[1]+1./3.*k[2]+1./6.*k[3];
    d[1]=-1./8.*k[0]-0.25*k[1]+0.25*k[2]+1./8.*k[3];
    d[2]=7./12.*k[0]-5./6.*k[1]+1./6.*k[2]+1./12.*k[3];
    d[3]=1./12.*k[0]+1./6.*k[1]-5./6.*k[2]+7./12.*k[3];
    free(k);
}


/****************************************************************************/


void tfd1Dinv_bl_V1(float *Input,float *Output, int Np) /* transformee avec bords libres avec splines quadratique */
{
    float *c = Output;
    float *d = Input;
    int N=Np-2;
    float h1[2*maxf+1]={0.,0.,0.,0.,0.,0.,1./4.,3./4.,3./4.,1./4.,0.,0.,0.,0.,0.};
    float g1[2*maxf+1]={0.,0.,0.,0.,0.,0.,-1./2.,-3./2.,3./2.,1./2.,0.,0.,0.,0.,0.};
    int j,i,l;
    float *cw;
    cw = new float[N+2];
    
    for (i=0;i<=N;i++) {
        cw[i]=d[i];
        c[i]=d[i];
    }
    c[0]=d[0]-2.*d[1]+d[2];
    c[1]=d[0]-d[1]-0.5*d[2];
    c[2]=d[0]+d[1]-0.5*d[3];
    c[3]=d[0]+2.*d[1]+d[3];
    for (j=2;j<N;j=2*j) {
        for (i=0;i<=(j+1);i++) cw[i]=c[i];
        for (i=0;i<=(2*j+1);i++) c[i]=0.;
        for (i=2;i<j;i++) {
            for (l=-maxf;l<=maxf;l++) {
                c[(2*i-1+l+maxf*(2*j+2))%(2*j+2)]+=h1[l+maxf]*cw[i];
            }
        }
        for (i=1;i<(j-1);i++) {
            for (l=-maxf;l<=maxf;l++) {
                c[(2*i+1+l+maxf*(2*j+2))%(2*j+2)]+=g1[l+maxf]*d[j+2+i];
            }
        }
        c[0]+=d[j+2]+cw[0];
        c[1]+=7./4.*d[j+2]+0.5*cw[1]+0.5*cw[0];
        c[2]+=-9./8.*d[j+2]+0.75*cw[1];
        c[3]+=-3./8.*d[j+2]+0.25*cw[1];
        c[2*j+1]+=d[2*j+1]+cw[j+1];
        c[2*j]+=7./4.*d[2*j+1]+0.5*cw[j]+0.5*cw[j+1];
        c[2*j-1]+=-9./8.*d[2*j+1]+0.75*cw[j];
        c[2*j-2]+=-3./8.*d[2*j+1]+0.25*cw[j];
    }
    delete cw;
}



/****************************************************************************/


/***********************   Ondelettes a div nulle   **********************************************/


void tfodab2d(int N, float *c1,float *c2,float *ddn1,float *ddn2) /* transformee directe */
{
    int e,e1,e2,j,j1,j2,i,i1,i2,k,k1,k2;
    float *c,*d,*d1,*d2;
    c=(float*) malloc (sizeof(float)*(N+2));
    d=(float*) malloc (sizeof(float)*(N+2));
    d1=(float*) malloc (sizeof(float)*(N+2)*(N+2));
    d2=(float*) malloc (sizeof(float)*(N+2)*(N+2));
    
    float g1d[2*maxf+1]={0.,0.,0.,0.,0.,0.,0.,-1./4.,1./2.,-1./4.,0.,0.,0.,0.,0.};
    float h1d[2*maxf+1]={0.,0.,0.,0.,0.,-1./8.,1./4.,3./4.,1./4.,-1./8.,0.,0.,0.,0.,0.};
    
    float g2d[2*maxf+1]={0.,0.,0.,0.,0.,0.,1./8.,-3./8.,3./8.,-1./8.,0.,0.,0.,0.,0.};
    float h2d[2*maxf+1]={0.,0.,0.,0.,0.,0.,-1./4.,3./4.,3./4.,-1./4.,0.,0.,0.,0.,0.};
    
    /* transformee en ondelettes classiques de bord de c1 et c2 */
    for (i2=0;i2<(N+1);i2++) {
        for (i1=0;i1<(N+2);i1++) c[i1]=c1[(N+2)*i2+i1];
        tfd1D_bl_V1(c,d,N+2);
        for (i1=0;i1<(N+2);i1++) d1[(N+2)*i2+i1]=d[i1]; 
    }
    for (i1=0;i1<(N+2);i1++) {
        for (i2=0;i2<(N+1);i2++) c[i2]=d1[(N+2)*i2+i1];
        tfd1D_bl_V0(c,d,N+1);
        for (i2=0;i2<(N+1);i2++) d1[(N+2)*i2+i1]=d[i2];
    }
    
    for (i2=0;i2<(N+2);i2++) {
        for (i1=0;i1<(N+1);i1++) c[i1]=c2[(N+1)*i2+i1];
        tfd1D_bl_V0(c,d,N+1);
        for (i1=0;i1<(N+1);i1++) d2[(N+1)*i2+i1]=d[i1];
    }
    for (i1=0;i1<(N+1);i1++) {
        for (i2=0;i2<(N+2);i2++) c[i2]=d2[(N+1)*i2+i1];
        tfd1D_bl_V1(c,d,N+2);
        for (i2=0;i2<(N+2);i2++) d2[(N+1)*i2+i1]=d[i2];
    }
   // cout << "OK 2 " << endl;
    
    /* transformee en ondelettes a divergence nulle */
    for (j1=N;j1>0;j1=j1) { j1=j1/2;
        for (j2=N;j2>0;j2=j2) { j2=j2/2;
            for (i1=0;i1<(j1+(j1<=1));i1++) {
                for (i2=0;i2<(j2+(j2<=1));i2++) 
                {
                    e1=(i1==0)-(j1==0); e2=(i2==0)-(j2==0); e=1-2*e1-2*e2+4*e1*e2; k1=j1+(j1==0); k2=j2+(j2==0);
                    long x = (1+(j2>1)+j2+i2)*(N+2)+1+(j1>1)+j1+i1;
                    long y = ((j2>1)+j2+i2)*(N+1)+(j1>1)+j1+i1;
                   //  if ((x > (N+
                    ddn1[(1+(j2>1)+j2+i2)*(N+2)+1+(j1>1)+j1+i1]=(k2*d1[((j2>1)+j2+i2)*(N+2)+1+(j1>1)+j1+i1]-e*k1*d2[(1+(j2>1)+j2+i2)*(N+1)+(j1>1)+j1+i1])/(k1*k1+k2*k2);
                    ddn2[((j2>1)+j2+i2)*(N+1)+(j1>1)+j1+i1]=(k1*d1[((j2>1)+j2+i2)*(N+2)+1+(j1>1)+j1+i1]+e*k2*d2[(1+(j2>1)+j2+i2)*(N+1)+(j1>1)+j1+i1])/(k1*k1+k2*k2);
                }
            }
        }
    }
   // cout << "OK 3 " << endl;

    
    for (i=0;i<(N+1);i++) {
        ddn1[(i+1)*(N+2)+0]=d1[i*(N+2)+0];
        ddn1[0*(N+2)+i+1]=d2[0*(N+1)+i];
    }
  //  cout << "OK 4 " << endl;

    
}

void tfodab2dinv(int N, float *c1,float *c2,float *ddn1,float *ddn2) /* transformee inverse */
{
cout << "tfodab2dinv - N = " << N << endl;
    int e,e1,e2,j,j1,j2,i,i1,i2,k,k1,k2;
    float *c,*d,*d1,*d2;
    
    float h1i[2*maxf+1]={0.,0.,0.,0.,0.,0.,1./2.,1.,1./2.,0.,0.,0.,0.,0.,0.};
    float g1i[2*maxf+1]={0.,0.,0.,0.,0.,0.,-1./4.,-1./2.,3./2.,-1./2.,-1./4.,0.,0.,0.,0.};
    
    float h2i[2*maxf+1]={0.,0.,0.,0.,0.,0.,1./4.,3./4.,3./4.,1./4.,0.,0.,0.,0.,0.};
    float g2i[2*maxf+1]={0.,0.,0.,0.,0.,0.,-1./2.,-3./2.,3./2.,1./2.,0.,0.,0.,0.,0.};
    
    c=(float*) malloc (sizeof(float)*(N+2));
    d=(float*) malloc (sizeof(float)*(N+2));
    d1=(float*) malloc (sizeof(float)*(N+2)*(N+2));
    d2=(float*) malloc (sizeof(float)*(N+2)*(N+2));
    
    /* On recalcule les coefficients d'ondelettes classiques a partir des coeff d'ond a div nulle + complementaires */
    for (i=0;i<(N+1);i++) {
        d1[i*(N+2)+0]=ddn1[(i+1)*(N+2)+0];
        d2[0*(N+1)+i]=ddn1[0*(N+2)+i+1];
    }
    for (j1=N;j1>0;j1=j1) { j1=j1/2;
        for (j2=N;j2>0;j2=j2) { j2=j2/2;
            for (i1=0;i1<(j1+(j1<=1));i1++) {
                for (i2=0;i2<(j2+(j2<=1));i2++) {
                    e1=(i1==0)-(j1==0); e2=(i2==0)-(j2==0); e=1-2*e1-2*e2+4*e1*e2; k1=j1+(j1==0); k2=j2+(j2==0);
                    d1[((j2>1)+j2+i2)*(N+2)+1+(j1>1)+j1+i1]=k2*ddn1[(1+(j2>1)+j2+i2)*(N+2)+1+(j1>1)+j1+i1]+k1*ddn2[((j2>1)+j2+i2)*(N+1)+(j1>1)+j1+i1];
                    d2[(1+(j2>1)+j2+i2)*(N+1)+(j1>1)+j1+i1]=-e*k1*ddn1[(1+(j2>1)+j2+i2)*(N+2)+1+(j1>1)+j1+i1]+e*k2*ddn2[((j2>1)+j2+i2)*(N+1)+(j1>1)+j1+i1];
                }
            }
        }
    }
    
    /* Transformee en ondelettes inverse classique (dans ddn1) */
    for (i2=0;i2<(N+1);i2++) {
        for (i1=0;i1<(N+2);i1++) d[i1]=d1[(N+2)*i2+i1];
        tfd1Dinv_bl_V1(c,d,N+2);
        for (i1=0;i1<(N+2);i1++) c1[(N+2)*i2+i1]=c[i1];
    }
    for (i1=0;i1<(N+2);i1++) {
        for (i2=0;i2<(N+1);i2++) d[i2]=c1[(N+2)*i2+i1];
        tfd1Dinv_bl_V0(c,d,N+1);
        for (i2=0;i2<(N+1);i2++) c1[(N+2)*i2+i1]=c[i2];
    }
    
    for (i2=0;i2<(N+2);i2++) {
        for (i1=0;i1<(N+1);i1++) d[i1]=d2[(N+1)*i2+i1];
        tfd1Dinv_bl_V0(c,d,N+1);
        for (i1=0;i1<(N+1);i1++) c2[(N+1)*i2+i1]=c[i1];
    }
    for (i1=0;i1<(N+1);i1++) {
        for (i2=0;i2<(N+2);i2++) d[i2]=c2[(N+1)*i2+i1];
        tfd1Dinv_bl_V1(c,d,N+2);
        for (i2=0;i2<(N+2);i2++) c2[(N+1)*i2+i1]=c[i2];
    }
}


/***********************   Ondelettes gradient   **********************************************/

void tfogab2d(int N, float *c1,float *c2,float *dg1,float *dg2) /* transformee directe */
{
    int e,e1,e2,j,j1,j2,i,i1,i2,k,k1,k2;
    float *c,*d,*d1,*d2;
    c=(float*) malloc (sizeof(float)*(N+2));
    d=(float*) malloc (sizeof(float)*(N+2));
    d1=(float*) malloc (sizeof(float)*(N+1)*(N+2));
    d2=(float*) malloc (sizeof(float)*(N+2)*(N+1));
    
    float g1d[2*maxf+1]={0.,0.,0.,0.,0.,0.,0.,-1./4.,1./2.,-1./4.,0.,0.,0.,0.,0.};
    float h1d[2*maxf+1]={0.,0.,0.,0.,0.,-1./8.,1./4.,3./4.,1./4.,-1./8.,0.,0.,0.,0.,0.};
    
    float g2d[2*maxf+1]={0.,0.,0.,0.,0.,0.,1./8.,-3./8.,3./8.,-1./8.,0.,0.,0.,0.,0.};
    float h2d[2*maxf+1]={0.,0.,0.,0.,0.,0.,-1./4.,3./4.,3./4.,-1./4.,0.,0.,0.,0.,0.};
    
    /* transformee en ondelettes classiques de bord de c1 et c2 */
    for (i2=0;i2<(N+2);i2++) {
        for (i1=0;i1<(N+1);i1++) c[i1]=c1[(N+1)*i2+i1];
        tfd1D_bl_V0(c,d,N+1);
        for (i1=0;i1<(N+1);i1++) d1[(N+1)*i2+i1]=d[i1];
    }
    for (i1=0;i1<(N+1);i1++) {
        for (i2=0;i2<(N+2);i2++) c[i2]=d1[(N+1)*i2+i1];
        tfd1D_bl_V1(c,d,N+2);
        for (i2=0;i2<(N+2);i2++) d1[(N+1)*i2+i1]=d[i2];
    }
    
    for (i2=0;i2<(N+1);i2++) {
        for (i1=0;i1<(N+2);i1++) c[i1]=c2[(N+2)*i2+i1];
        tfd1D_bl_V1(c,d,N+2);
        for (i1=0;i1<(N+2);i1++) d2[(N+2)*i2+i1]=d[i1];
    }
    for (i1=0;i1<(N+2);i1++) {
        for (i2=0;i2<(N+1);i2++) c[i2]=d2[(N+2)*i2+i1];
        tfd1D_bl_V0(c,d,N+1);
        for (i2=0;i2<(N+1);i2++) d2[(N+2)*i2+i1]=d[i2];
    }
    
    /* transformee en ondelettes gradient (dans dg2) */
    for (j1=N;j1>0;j1=j1) { j1=j1/2;
        for (j2=N;j2>0;j2=j2) { j2=j2/2;
            for (i1=0;i1<(j1+(j1<=1));i1++) {
                for (i2=0;i2<(j2+(j2<=1));i2++) {
                    e1=(i1==0)-(j1==0); e2=(i2==0)-(j2==0); e=1-2*e1-2*e2+4*e1*e2; k1=j1+(j1==0); k2=j2+(j2==0);
                    dg1[((j2>1)+j2+i2)*(N+1)+(j1>1)+j1+i1]=(k2*d1[(1+(j2>1)+j2+i2)*(N+1)+(j1>1)+j1+i1]-e*k1*d2[((j2>1)+j2+i2)*(N+2)+1+(j1>1)+j1+i1])/(k1*k1+k2*k2);
                    dg2[(1+(j2>1)+j2+i2)*(N+2)+1+(j1>1)+j1+i1]=(k1*d1[(1+(j2>1)+j2+i2)*(N+1)+(j1>1)+j1+i1]+e*k2*d2[((j2>1)+j2+i2)*(N+2)+1+(j1>1)+j1+i1])/(k1*k1+k2*k2);
                }
            }
        }
    }
    for (i=0;i<(N+1);i++) {
        dg2[(i+1)*(N+2)+0]=d2[i*(N+2)+0];
        dg2[0*(N+2)+i+1]=d1[0*(N+1)+i];
    }
}

void tfogab2dinv(int N, float *c1,float *c2,float *dg1,float *dg2) /* transformee inverse */
{
    int e,e1,e2,j,j1,j2,i,i1,i2,k,k1,k2;
    float *c,*d,*d1,*d2;
    
    float h1i[2*maxf+1]={0.,0.,0.,0.,0.,0.,1./2.,1.,1./2.,0.,0.,0.,0.,0.,0.};
    float g1i[2*maxf+1]={0.,0.,0.,0.,0.,0.,-1./4.,-1./2.,3./2.,-1./2.,-1./4.,0.,0.,0.,0.};
    
    float h2i[2*maxf+1]={0.,0.,0.,0.,0.,0.,1./4.,3./4.,3./4.,1./4.,0.,0.,0.,0.,0.};
    float g2i[2*maxf+1]={0.,0.,0.,0.,0.,0.,-1./2.,-3./2.,3./2.,1./2.,0.,0.,0.,0.,0.};
    
    c=(float*) malloc (sizeof(float)*(N+2));
    d=(float*) malloc (sizeof(float)*(N+2));
    d1=(float*) malloc (sizeof(float)*(N+1)*(N+2));
    d2=(float*) malloc (sizeof(float)*(N+2)*(N+1));
    
    /* On recalcule les coefficients d'ondelettes classiques a partir des coeff d'ond grad + complementaires */
    for (i=0;i<(N+1);i++) {
        d2[i*(N+2)+0]=dg2[(i+1)*(N+2)+0];
        d1[0*(N+1)+i]=dg2[0*(N+2)+i+1];
    }
    for (j1=N;j1>0;j1=j1) { j1=j1/2;
        for (j2=N;j2>0;j2=j2) { j2=j2/2;
            for (i1=0;i1<(j1+(j1<=1));i1++) {
                for (i2=0;i2<(j2+(j2<=1));i2++) {
                    e1=(i1==0)-(j1==0); e2=(i2==0)-(j2==0); e=1-2*e1-2*e2+4*e1*e2; k1=j1+(j1==0); k2=j2+(j2==0);
                    d1[(1+(j2>1)+j2+i2)*(N+1)+(j1>1)+j1+i1]=k2*dg1[((j2>1)+j2+i2)*(N+1)+(j1>1)+j1+i1]+k1*dg2[(1+(j2>1)+j2+i2)*(N+2)+1+(j1>1)+j1+i1];
                    d2[((j2>1)+j2+i2)*(N+2)+1+(j1>1)+j1+i1]=-e*k1*dg1[((j2>1)+j2+i2)*(N+1)+(j1>1)+j1+i1]+e*k2*dg2[(1+(j2>1)+j2+i2)*(N+2)+1+(j1>1)+j1+i1];
                }
            }
        }
    }
    
    /* Transformee en ondelettes inverse classique */
    for (i2=0;i2<(N+2);i2++) {
        for (i1=0;i1<(N+1);i1++) d[i1]=d1[(N+1)*i2+i1];
        tfd1Dinv_bl_V0(c,d,N+1);
        for (i1=0;i1<(N+1);i1++) c1[(N+1)*i2+i1]=c[i1];
    }
    for (i1=0;i1<(N+1);i1++) {
        for (i2=0;i2<(N+2);i2++) d[i2]=c1[(N+1)*i2+i1];
        tfd1Dinv_bl_V1(c,d,N+2);
        for (i2=0;i2<(N+2);i2++) c1[(N+1)*i2+i1]=c[i2];
    }
    
    for (i2=0;i2<(N+1);i2++) {
        for (i1=0;i1<(N+2);i1++) d[i1]=d2[(N+2)*i2+i1];
        tfd1Dinv_bl_V1(c,d,N+2);
        for (i1=0;i1<(N+2);i1++) c2[(N+2)*i2+i1]=c[i1];
    }
    for (i1=0;i1<(N+2);i1++) {
        for (i2=0;i2<(N+1);i2++) d[i2]=c2[(N+2)*i2+i1];
        tfd1Dinv_bl_V0(c,d,N+1);
        for (i2=0;i2<(N+1);i2++) c2[(N+2)*i2+i1]=c[i2];
    }
}

/*

void qi10(float *u,float *v, int N)
{
    int i1,i2;
    for (i2=0;i2<(N+1);i2++) {
        v[i2*(N+2)+0]=u[i2*(N+1)+0];
        v[i2*(N+2)+1]=u[i2*(N+1)+1]+0.5*(u[i2*(N+1)+1]-u[i2*(N+1)+2]);
        for (i1=2;i1<N;i1++) {
            v[i2*(N+2)+i1]=0.5*(u[i2*(N+1)+i1-1]+u[i2*(N+1)+i1]);
        }
        v[i2*(N+2)+N]=u[i2*(N+1)+N-1]+0.5*(u[i2*(N+1)+N-1]-u[i2*(N+1)+N-2]);
        v[i2*(N+2)+N+1]=u[i2*(N+1)+N];
    }
}

void qi01(float *u,float *v, int N)
{
    int i1,i2;
    for (i1=0;i1<(N+1);i1++) {
        v[0*(N+1)+i1]=u[0*(N+1)+i1];
        v[1*(N+1)+i1]=u[1*(N+1)+i1]+0.5*(u[1*(N+1)+i1]-u[2*(N+1)+i1]);
        for (i2=2;i2<N;i2++) {
            v[i2*(N+1)+i1]=0.5*(u[(i2-1)*(N+1)+i1]+u[i2*(N+1)+i1]);
        }
        v[N*(N+1)+i1]=u[(N-1)*(N+1)+i1]+0.5*(u[(N-1)*(N+1)+i1]-u[(N-2)*(N+1)+i1]);
        v[(N+1)*(N+1)+i1]=u[N*(N+1)+i1];
    }
}

 void rec10(float *u,float *v, int N)
{
    int i1,i2;
    for (i2=0;i2<(N+1);i2++) {
        u[i2*(N+1)+0]=v[i2*(N+2)+0];
        for (i1=1;i1<N;i1++) {
            u[i2*(N+1)+i1]=0.5*(v[i2*(N+2)+i1]+v[i2*(N+2)+i1+1]);
        }
        u[i2*(N+1)+N]=v[i2*(N+2)+N+1];
    }
}

void rec01(float *u,float *v, int N)
{
    int i1,i2;
    for (i1=0;i1<(N+1);i1++) {
        u[0*(N+1)+i1]=v[0*(N+1)+i1];
        for (i2=1;i2<N;i2++) {
            u[i2*(N+1)+i1]=0.5*(v[i2*(N+1)+i1]+v[(i2+1)*(N+1)+i1]);
        }
        u[N*(N+1)+i1]=v[(N+1)*(N+1)+i1];
    }
}
*/


 
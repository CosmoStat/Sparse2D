#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#define N 256
#define pi 3.1415926535898
#define maxf 7  /* le filtre est contenu dans [-maxf,+maxf]  */

#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr

void tfd1D(float *c,float *d,float *h1,float *g1) /* transformee en ondelettes 1D classique */
{
 int j,i,l;
 float *k;
 k=(float*) malloc(sizeof(float)*N);
 
 for (i=0;i<N;i++) {
  k[i]=c[i];
 }
 for (j=N/2;j>0;j=j/2) {
  for (i=0;i<j;i++) {
   d[i]=0.; d[j+i]=0.;
   for (l=-maxf;l<=maxf;l++) {
    d[i]+=h1[l+maxf]*k[(2*i+l+2*maxf*j)%(2*j)];
    d[j+i]+=g1[l+maxf]*k[(2*i+l+2*maxf*j)%(2*j)];
   }
  }
  for (i=0;i<j;i++) {
   k[i]=d[i];
  }
 }
 free(k);
}

void tfd1Dinv(float *c,float *d,float *h1,float *g1) /* transformee 1D inverse classique */
{
 int j,i,l;
 float *cw;
 cw=(float*) malloc(sizeof(float)*N);
 
 for (i=0;i<N;i++) {
  cw[i]=d[i];
  c[i]=d[i];
 }
 for (j=1;j<N;j=2*j) {
  for (i=0;i<2*j;i++) {
   cw[i]=c[i];
   c[i]=0.;
  }
  for (i=0;i<j;i++) {
   for (l=-maxf;l<=maxf;l++) {
    c[(2*i+l+2*maxf*j)%(2*j)]+=h1[l+maxf]*cw[i]+g1[l+maxf]*cw[j+i];
   }
  }
 }
 free(cw);
}

void tfd1D_bl_V0(float *c,float *d,float *h1,float *g1) /* transformee avec bords libres avec splines lineaires */
{
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

void tfd1Dinv_bl_V0(float *c,float *d,float *h1,float *g1) /* transformee inverse bords libres */
{
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

/***********************   Ondelettes dans V1   **********************************************/

void tfd1D_bl_V1(float *c,float *d,float *h1,float *g1) /* transformee avec bords libres avec splines quadratiques */
{
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

void tfd1Dinv_bl_V1(float *c,float *d,float *h1,float *g1) /* transformee inverse bords libres */
{
 int j,i,l;
 float *cw;
 cw=(float*) malloc(sizeof(float)*(N+2));
 
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
 free(cw);
}

void tfodab2d(float *c1,float *c2,float *ddn1,float *ddn2) /* transformee en ondelettes a div nulle avec bords libres */
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
 for (i2=0;i2<(N+1);i2++) {
  for (i1=0;i1<(N+2);i1++) c[i1]=c1[(N+2)*i2+i1];
  tfd1D_bl_V1(c,d,h2d,g2d);
  for (i1=0;i1<(N+2);i1++) d1[(N+2)*i2+i1]=d[i1]; 
 }
 for (i1=0;i1<(N+2);i1++) {
  for (i2=0;i2<(N+1);i2++) c[i2]=d1[(N+2)*i2+i1];
  tfd1D_bl_V0(c,d,h1d,g1d);
  for (i2=0;i2<(N+1);i2++) d1[(N+2)*i2+i1]=d[i2];
 }
 
 for (i2=0;i2<(N+2);i2++) {
  for (i1=0;i1<(N+1);i1++) c[i1]=c2[(N+1)*i2+i1];
  tfd1D_bl_V0(c,d,h1d,g1d);
  for (i1=0;i1<(N+1);i1++) d2[(N+1)*i2+i1]=d[i1];
 }
 for (i1=0;i1<(N+1);i1++) {
  for (i2=0;i2<(N+2);i2++) c[i2]=d2[(N+1)*i2+i1];
  tfd1D_bl_V1(c,d,h2d,g2d);
  for (i2=0;i2<(N+2);i2++) d2[(N+1)*i2+i1]=d[i2];
 }
 
 /* transformee en ondelettes a divergence nulle */
 for (j1=N;j1>0;j1=j1) { j1=j1/2;
  for (j2=N;j2>0;j2=j2) { j2=j2/2;
   for (i1=0;i1<(j1+(j1<=1));i1++) {
    for (i2=0;i2<(j2+(j2<=1));i2++) {
     e1=(i1==0)-(j1==0); e2=(i2==0)-(j2==0); e=1-2*e1-2*e2+4*e1*e2; k1=j1+(j1==0); k2=j2+(j2==0);
     ddn1[(1+(j2>1)+j2+i2)*(N+2)+1+(j1>1)+j1+i1]=(k2*d1[((j2>1)+j2+i2)*(N+2)+1+(j1>1)+j1+i1]-e*k1*d2[(1+(j2>1)+j2+i2)*(N+1)+(j1>1)+j1+i1])/(k1*k1+k2*k2);
     ddn2[((j2>1)+j2+i2)*(N+1)+(j1>1)+j1+i1]=(k1*d1[((j2>1)+j2+i2)*(N+2)+1+(j1>1)+j1+i1]+e*k2*d2[(1+(j2>1)+j2+i2)*(N+1)+(j1>1)+j1+i1])/(k1*k1+k2*k2);
    }
   }
  }
 }
 for (i=0;i<(N+1);i++) {
  ddn1[(i+1)*(N+2)+0]=d1[i*(N+2)+0];
  ddn1[0*(N+2)+i+1]=d2[0*(N+1)+i];
 }
}

void tfodab2dinv(float *c1,float *c2,float *ddn1,float *ddn2) /* transformee inverse ond a div nulle avec bords libres */
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
  tfd1Dinv_bl_V1(c,d,h2i,g2i);
  for (i1=0;i1<(N+2);i1++) c1[(N+2)*i2+i1]=c[i1];
 }
 for (i1=0;i1<(N+2);i1++) {
  for (i2=0;i2<(N+1);i2++) d[i2]=c1[(N+2)*i2+i1];
  tfd1Dinv_bl_V0(c,d,h1i,g1i);
  for (i2=0;i2<(N+1);i2++) c1[(N+2)*i2+i1]=c[i2];
 }
 
 for (i2=0;i2<(N+2);i2++) {
  for (i1=0;i1<(N+1);i1++) d[i1]=d2[(N+1)*i2+i1];
  tfd1Dinv_bl_V0(c,d,h1i,g1i);
  for (i1=0;i1<(N+1);i1++) c2[(N+1)*i2+i1]=c[i1];
 }
 for (i1=0;i1<(N+1);i1++) {
  for (i2=0;i2<(N+2);i2++) d[i2]=c2[(N+1)*i2+i1];
  tfd1Dinv_bl_V1(c,d,h2i,g2i);
  for (i2=0;i2<(N+2);i2++) c2[(N+1)*i2+i1]=c[i2];
 }
}


void tfogab2d(float *c1,float *c2,float *dg1,float *dg2) /* transformee en ondelettes gradient avec bords libres */
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
  tfd1D_bl_V0(c,d,h1d,g1d);
  for (i1=0;i1<(N+1);i1++) d1[(N+1)*i2+i1]=d[i1];
 }
 for (i1=0;i1<(N+1);i1++) {
  for (i2=0;i2<(N+2);i2++) c[i2]=d1[(N+1)*i2+i1];
  tfd1D_bl_V1(c,d,h2d,g2d);
  for (i2=0;i2<(N+2);i2++) d1[(N+1)*i2+i1]=d[i2];
 }
 
 for (i2=0;i2<(N+1);i2++) {
  for (i1=0;i1<(N+2);i1++) c[i1]=c2[(N+2)*i2+i1];
  tfd1D_bl_V1(c,d,h2d,g2d);
  for (i1=0;i1<(N+2);i1++) d2[(N+2)*i2+i1]=d[i1];
 }
 for (i1=0;i1<(N+2);i1++) {
  for (i2=0;i2<(N+1);i2++) c[i2]=d2[(N+2)*i2+i1];
  tfd1D_bl_V0(c,d,h1d,g1d);
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

void tfogab2dinv(float *c1,float *c2,float *dg1,float *dg2) /* transformee inverse ond gradient avec bords libres */
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
  tfd1Dinv_bl_V0(c,d,h1i,g1i);
  for (i1=0;i1<(N+1);i1++) c1[(N+1)*i2+i1]=c[i1];
 }
 for (i1=0;i1<(N+1);i1++) {
  for (i2=0;i2<(N+2);i2++) d[i2]=c1[(N+1)*i2+i1];
  tfd1Dinv_bl_V1(c,d,h2i,g2i);
  for (i2=0;i2<(N+2);i2++) c1[(N+1)*i2+i1]=c[i2];
 }
 
 for (i2=0;i2<(N+1);i2++) {
  for (i1=0;i1<(N+2);i1++) d[i1]=d2[(N+2)*i2+i1];
  tfd1Dinv_bl_V1(c,d,h2i,g2i);
  for (i1=0;i1<(N+2);i1++) c2[(N+2)*i2+i1]=c[i1];
 }
 for (i1=0;i1<(N+2);i1++) {
  for (i2=0;i2<(N+1);i2++) d[i2]=c2[(N+2)*i2+i1];
  tfd1Dinv_bl_V0(c,d,h1i,g1i);
  for (i2=0;i2<(N+1);i2++) c2[(N+2)*i2+i1]=c[i2];
 }
}

 /* quasi-interpolation */
void qi10(float *u,float *v)
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
 
void qi01(float *u,float *v)
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

 /* reconstruction */
void rec10(float *u,float *v)
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
 
void rec01(float *u,float *v)
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
 
float g1d[2*maxf+1]={0.,0.,0.,0.,0.,0.,0.,-1./4.,1./2.,-1./4.,0.,0.,0.,0.,0.};
float h1i[2*maxf+1]={0.,0.,0.,0.,0.,0.,1./2.,1.,1./2.,0.,0.,0.,0.,0.,0.};
float h1d[2*maxf+1]={0.,0.,0.,0.,0.,-1./8.,1./4.,3./4.,1./4.,-1./8.,0.,0.,0.,0.,0.};
float g1i[2*maxf+1]={0.,0.,0.,0.,0.,0.,-1./4.,-1./2.,3./2.,-1./2.,-1./4.,0.,0.,0.,0.};

float g2d[2*maxf+1]={0.,0.,0.,0.,0.,0.,1./8.,-3./8.,3./8.,-1./8.,0.,0.,0.,0.,0.};
float h2i[2*maxf+1]={0.,0.,0.,0.,0.,0.,1./4.,3./4.,3./4.,1./4.,0.,0.,0.,0.,0.};
float h2d[2*maxf+1]={0.,0.,0.,0.,0.,0.,-1./4.,3./4.,3./4.,-1./4.,0.,0.,0.,0.,0.};
float g2i[2*maxf+1]={0.,0.,0.,0.,0.,0.,-1./2.,-3./2.,3./2.,1./2.,0.,0.,0.,0.,0.};

/* Coiflets 12 */
float gCoid[2*maxf+1]={0.,0.000720549445000,-0.001823208870999,-0.005611434818997,0.023680171946988,0.059434418645969,-0.076488599077960,-0.417005184423784,0.812723635449579,-0.386110066822800,-0.067372554721965,0.041464936781979,0.016387336462992,0.,0.};
float hCoid[2*maxf+1]={0.,0.,0.,0.016387336462992,-0.041464936781979,-0.067372554721965,0.386110066822800,0.812723635449579,0.417005184423784,-0.076488599077960,-0.059434418645969,0.023680171946988,0.005611434818997,-0.001823208870999,-0.000720549445000};
float hCoii[2*maxf+1]={0.,0.,0.,0.016387336462992,-0.041464936781979,-0.067372554721965,0.386110066822800,0.812723635449579,0.417005184423784,-0.076488599077960,-0.059434418645969,0.023680171946988,0.005611434818997,-0.001823208870999,-0.000720549445000};
float gCoii[2*maxf+1]={0.,0.000720549445000,-0.001823208870999,-0.005611434818997,0.023680171946988,0.059434418645969,-0.076488599077960,-0.417005184423784,0.812723635449579,-0.386110066822800,-0.067372554721965,0.041464936781979,0.016387336462992,0.,0.};

int main () {
 float *w,*c,*d,*u1,*u2,*v1,*v2,*w1,*w2,*c1,*c2,*ddn,*dg,*dn;
 int i,j,k,k1,k2,i1,i2,it,It=40;
 double x,y,x1,x2,mx,f1=2324.53,f2=426.63,h=1./N,p=1.2;
 FILE *data1=fopen("data1.sci","w");
 FILE *data2=fopen("data2.sci","w");
 FILE *res1=fopen("res1.sci","w");
 FILE *res2=fopen("res2.sci","w");
 
 // allocate memory 
 w=(float*) malloc (sizeof(float)*(N+2));
 c=(float*) malloc (sizeof(float)*(N+2));
 d=(float*) malloc (sizeof(float)*(N+2));
 u1=(float*) malloc (sizeof(float)*(N+1)*(N+1));
 u2=(float*) malloc (sizeof(float)*(N+1)*(N+1));
 v1=(float*) malloc (sizeof(float)*(N+1)*(N+2));
 v2=(float*) malloc (sizeof(float)*(N+2)*(N+1));
 w1=(float*) malloc (sizeof(float)*(N+1)*(N+1));
 w2=(float*) malloc (sizeof(float)*(N+1)*(N+1));
 ddn=(float*) malloc (sizeof(float)*(N+2)*(N+2));
 dg=(float*) malloc (sizeof(float)*(N+2)*(N+2));
 dn=(float*) malloc (sizeof(float)*(N+1)*(N+1));
 fprintf(data1,"v1=[\n",N,N);
 fprintf(data2,"v2=[\n",N,N);
 fprintf(res1,"w1=[\n",N,N);
 fprintf(res2,"w2=[\n",N,N);
 
 for (i2=0;i2<(N+1);i2++) {
  for (i1=0;i1<(N+1);i1++) {
   x2=((float) i2)/N;
   x1=((float) i1)/N;
   u1[(N+1)*i2+i1]=f2*cos(f1*pi*x1)*sin(f2*pi*x2)+exp(x1)*sin(x2);
   u2[(N+1)*i2+i1]=f1*sin(f1*pi*x1)*cos(f2*pi*x2)+exp(x1)*cos(x2);
  }
 }
 
  mx=0.;
  for (i=0;i<((N+1)*(N+1));i++) mx+=((u2[i])*(u2[i])+(u1[i])*(u1[i]))/(N+1)/(N+2);
  printf("erreur=%e\n",sqrt(mx));
 for (it=0;it<It;it++) {
  qi10(u1,v1);
  qi01(u2,v2);
  tfodab2d(v1,v2,ddn,dn);
  for (i=0;i<((N+1)*(N+1));i++)  dn[i]=0.;
  tfodab2dinv(v1,v2,ddn,dn);
  rec10(w1,v1);
  rec01(w2,v2);
  for (i=0;i<((N+1)*(N+1));i++) { u1[i]-=p*w1[i]; u2[i]-=p*w2[i];}
 
  qi01(u1,v1);
  qi10(u2,v2);
  tfogab2d(v1,v2,dn,dg);
  for (i=0;i<((N+1)*(N+1));i++)  dn[i]=0.;
  tfogab2dinv(v1,v2,dn,dg);
  rec01(w1,v1);
  rec10(w2,v2);
  for (i=0;i<((N+1)*(N+1));i++) { u1[i]-=p*w1[i]; u2[i]-=p*w2[i];}
 
  mx=0.;
  for (i=0;i<((N+1)*(N+1));i++) mx+=(u2[i]*u2[i]+u1[i]*u1[i])/(N+1)/(N+1);
  printf("erreur=%e\n",sqrt(mx));
 
 }
 
 /*
 tfogab2d(v1,v2,dn,dg);
 
 for (i=0;i<((N+1)*(N+1));i++)  dn[i]=0.;
// for (i=0;i<((N+2)*(N+2));i++)  dg[i]=0.;
 
 tfogab2dinv(w1,w2,dn,dg);
 */
 
 mx=0.;
 for (i=0;i<((N+1)*(N+1));i++) mx+=((u2[i])*(u2[i])+(u1[i])*(u1[i]))/(N+1)/(N+2);
// printf("norme2u=%e\n",sqrt(mx));
 mx=0.;
 for (i=0;i<((N+1)*(N+2));i++) mx+=((v2[i])*(v2[i])+(v1[i])*(v1[i]))/(N+1)/(N+2);
// printf("norme2v=%e\n",sqrt(mx));
 mx=0.;
 for (i=0;i<((N+1)*(N+1));i++) mx+=((w2[i])*(w2[i])+(w1[i])*(w1[i]))/(N+1)/(N+1);
// printf("norme2w=%e\n",sqrt(mx));
 mx=0.;
 for (i=0;i<((N+1)*(N+1));i++) mx+=((w2[i]-u2[i])*(w2[i]-u2[i])+(w1[i]-u1[i])*(w1[i]-u1[i]))/(N+1)/(N+1);
// printf("erreur=%e\n",sqrt(mx));
}



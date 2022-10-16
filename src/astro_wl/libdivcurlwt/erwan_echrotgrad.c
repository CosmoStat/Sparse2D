#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>

#define N 512
#define pi 3.1415926535898

#define TWOPI (6.2831853071795864769252867665590057683943387987502) /* 2 * pi */
#define RAND (rand()+1)/((double) RAND_MAX +2)
#define RANDN (sqrt(-2.0*log(RAND))*cos(TWOPI*RAND))

#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr

/* Transformee en ondelettes a divergence nulle (projection oblique sur Hdiv) */

void tfodanp2d(float c[2][N][N],float ddn[2][N][N])
{
 int j,j1,j2,i1,i2,e1,e2,l1,l2;
 float k1[N],k2[N];
 float d[2][N][N];
 float h[2][5]={{-1./8.,1./4.,3./4.,1./4.,-1./8.},{0.,-1./4.,3./4.,3./4.,-1./4.}};
 float g[2][5]={{0.,0.,-1./4.,1./2.,-1./4},{0.,1./8.,-3./8.,3./8.,-1./8.}};

 for (i2=0;i2<N;i2++) {
  for (i1=0;i1<N;i1++) {
   d[0][i1][i2]=c[0][i1][i2];
   d[1][i1][i2]=c[1][i1][i2];
  }
 }
  
 for (j=N/2;j>0;j=j/2) {
  for (i2=0;((i2<4*j)&&(i2<N));i2++) {
   for (i1=0;i1<2*j;i1++) {
    k1[i1]=d[0][i1][i2];
    k2[i1]=d[1][i1][i2];
   }
   for (i1=0;i1<j;i1++) {
    d[0][i1][i2]=h[1][1]*k1[(2*i1-1+2*j)%(2*j)]+h[1][2]*k1[2*i1]+h[1][3]*k1[2*i1+1]+h[1][4]*k1[(2*i1+2)%(2*j)];
    d[0][j+i1][i2]=g[1][1]*k1[(2*i1-1+2*j)%(2*j)]+g[1][2]*k1[2*i1]+g[1][3]*k1[2*i1+1]+g[1][4]*k1[(2*i1+2)%(2*j)];
    d[1][i1][i2]=h[0][0]*k2[(2*i1-2+2*j)%(2*j)]+h[0][1]*k2[(2*i1-1+2*j)%(2*j)]+h[0][2]*k2[2*i1]+h[0][3]*k2[2*i1+1]+h[0][4]*k2[(2*i1+2)%(2*j)];
    d[1][j+i1][i2]=g[0][2]*k2[2*i1]+g[0][3]*k2[2*i1+1]+g[0][4]*k2[(2*i1+2)%(2*j)];
   }
  }
  for (i1=0;((i1<4*j)&&(i1<N));i1++) {
   for (i2=0;i2<2*j;i2++) {
    k1[i2]=d[0][i1][i2];
    k2[i2]=d[1][i1][i2];
   }
   for (i2=0;i2<j;i2++) {
    d[0][i1][i2]=h[0][0]*k1[(2*i2-2+2*j)%(2*j)]+h[0][1]*k1[(2*i2-1+2*j)%(2*j)]+h[0][2]*k1[2*i2]+h[0][3]*k1[2*i2+1]+h[0][4]*k1[(2*i2+2)%(2*j)];
    d[0][i1][j+i2]=g[0][2]*k1[2*i2]+g[0][3]*k1[2*i2+1]+g[0][4]*k1[(2*i2+2)%(2*j)];
    d[1][i1][i2]=h[1][1]*k2[(2*i2-1+2*j)%(2*j)]+h[1][2]*k2[2*i2]+h[1][3]*k2[2*i2+1]+h[1][4]*k2[(2*i2+2)%(2*j)];
    d[1][i1][j+i2]=g[1][1]*k2[(2*i2-1+2*j)%(2*j)]+g[1][2]*k2[2*i2]+g[1][3]*k2[2*i2+1]+g[1][4]*k2[(2*i2+2)%(2*j)];
   }
  }
 }
 
 for (j=N/2;j>1;j=j/2) {
  for (j1=j;j1>j/4;j1=j1/2) {
   for (j2=j;j2>j/4;j2=j2/2) {
    for (i1=0;i1<j1;i1++) {
     for (i2=0;i2<j2;i2++) {
      ddn[0][j1+i1][j2+i2]=(j2*d[0][j1+i1][j2+i2]-j1*d[1][j1+i1][j2+i2])/(j1*j1+j2*j2);
      ddn[1][j1+i1][j2+i2]=(j1*d[0][j1+i1][j2+i2]+j2*d[1][j1+i1][j2+i2])/(j1*j1+j2*j2);
     }
    }
   }
  }
  for (i1=j;i1<2*j;i1++) {
   for (i2=0;i2<j/2;i2++) {
    ddn[0][i1][i2]=-d[1][i1][i2]/8.;
    ddn[1][i1][i2]=d[0][i1][i2]+d[1][i1][i2]/8.-d[1][i1][(i2-1+j/2)%(j/2)]/8.;
    ddn[0][i2][i1]=-d[0][i2][i1]/8.;
    ddn[1][i2][i1]=d[1][i2][i1]+d[0][i2][i1]/8.-d[0][(i2-1+j/2)%(j/2)][i1]/8.;
   }
  }
 }
 
 ddn[0][1][0]=d[1][1][0];
 ddn[1][1][0]=d[0][1][0];
 ddn[0][0][1]=d[0][0][1];
 ddn[1][0][1]=d[1][0][1];
 ddn[0][0][0]=d[0][0][0];
 ddn[1][0][0]=d[1][0][0];
}

/* Transformee inverse partant des ondelettes a divergence nulle */

void tfodanpinv2d(float c[2][N][N],float ddn[2][N][N])
{
 int j,j1,j2,i1,i2,e1,e2,l;
 float k1[N],k2[N];
 float h[2][5]={{1./2.,1.,1./2.,0.,0.},{1./4.,3./4.,3./4.,1./4.,0.}};
 float g[2][5]={{-1./4.,-1./2.,3./2.,-1./2.,-1./4.},{-1./2.,-3./2.,3./2.,1./2.,0.}};
 
 c[0][0][0]=ddn[0][0][0];
 c[1][0][0]=ddn[1][0][0];
 c[1][1][0]=ddn[0][1][0];
 c[0][1][0]=ddn[1][1][0];
 c[0][0][1]=ddn[0][0][1];
 c[1][0][1]=ddn[1][0][1];
 for (j=N/2;j>1;j=j/2) {
  for (j1=j;j1>j/4;j1=j1/2) {
   for (j2=j;j2>j/4;j2=j2/2) {
    for (i1=0;i1<j1;i1++) {
     for (i2=0;i2<j2;i2++) {
      c[0][j1+i1][j2+i2]=j2*ddn[0][j1+i1][j2+i2]+j1*ddn[1][j1+i1][j2+i2];
      c[1][j1+i1][j2+i2]=-j1*ddn[0][j1+i1][j2+i2]+j2*ddn[1][j1+i1][j2+i2];
     }
    }
   }
  }
  for (i1=j;i1<2*j;i1++) {
   for (i2=0;i2<j/2;i2++) {
    c[1][i1][i2]=-8.*ddn[0][i1][i2];
    c[0][i1][i2]=ddn[1][i1][i2]+ddn[0][i1][i2]-ddn[0][i1][(i2-1+j/2)%(j/2)];
    c[0][i2][i1]=-8.*ddn[0][i2][i1];
    c[1][i2][i1]=ddn[1][i2][i1]+ddn[0][i2][i1]-ddn[0][(i2-1+j/2)%(j/2)][i1];
   }
  }
 }
 
 for (i1=0;i1<N;i1++) {
  k1[i1]=0.;
  k2[i1]=0.;
 }
 
 for (j=1;j<N;j=2*j) {
  for (i2=0;((i2<4*j)&&(i2<N));i2++) {
   for (i1=0;i1<j;i1++) {
    for (l=-1;l<4;l++) {
     k1[(2*i1+l+2*j)%(2*j)]=h[1][l+1]*c[0][i1][i2]+g[1][l+1]*c[0][j+i1][i2]+k1[(2*i1+l+2*j)%(2*j)];
     k2[(2*i1+l+2*j)%(2*j)]=h[0][l+1]*c[1][i1][i2]+g[0][l+1]*c[1][j+i1][i2]+k2[(2*i1+l+2*j)%(2*j)];
    }
   }
   for (i1=0;i1<2*j;i1++) {
    c[0][i1][i2]=k1[i1];
    k1[i1]=0.;
    c[1][i1][i2]=k2[i1];
    k2[i1]=0.;
   }
  }
  
  for (i1=0;((i1<4*j)&&(i1<N));i1++) {
   for (i2=0;i2<j;i2++) {
    for (l=-1;l<4;l++) {
     k1[(2*i2+l+2*j)%(2*j)]=h[0][l+1]*c[0][i1][i2]+g[0][l+1]*c[0][i1][j+i2]+k1[(2*i2+l+2*j)%(2*j)];
     k2[(2*i2+l+2*j)%(2*j)]=h[1][l+1]*c[1][i1][i2]+g[1][l+1]*c[1][i1][j+i2]+k2[(2*i2+l+2*j)%(2*j)];
    }
   }
   for (i2=0;i2<2*j;i2++) {
    c[0][i1][i2]=k1[i2];
    k1[i2]=0.;
    c[1][i1][i2]=k2[i2];
    k2[i2]=0.;
   }
  }
 }

}

/* Transformee en ondelettes gradients (projection oblique sur Hrot) */

void tfoganp2d(float c[2][N][N],float ddn[2][N][N])
{
 int j,j1,j2,i1,i2,e1,e2,l1,l2;
 float k1[N],k2[N];
 float d[2][N][N];
 float h[2][5]={{-1./8.,1./4.,3./4.,1./4.,-1./8.},{0.,-1./4.,3./4.,3./4.,-1./4.}};
 float g[2][5]={{0.,0.,-1./4.,1./2.,-1./4},{0.,1./8.,-3./8.,3./8.,-1./8.}};

 for (i2=0;i2<N;i2++) {
  for (i1=0;i1<N;i1++) {
   d[0][i1][i2]=c[0][i1][i2];
   d[1][i1][i2]=c[1][i1][i2];
  }
 }
  
 for (j=N/2;j>0;j=j/2) {
  for (i2=0;((i2<4*j)&&(i2<N));i2++) {
   for (i1=0;i1<2*j;i1++) {
    k1[i1]=d[0][i1][i2];
    k2[i1]=d[1][i1][i2];
   }
   for (i1=0;i1<j;i1++) {
    d[0][i1][i2]=h[0][0]*k1[(2*i1-2+2*j)%(2*j)]+h[0][1]*k1[(2*i1-1+2*j)%(2*j)]+h[0][2]*k1[2*i1]+h[0][3]*k1[2*i1+1]+h[0][4]*k1[(2*i1+2)%(2*j)];
    d[0][j+i1][i2]=g[0][1]*k1[(2*i1-1+2*j)%(2*j)]+g[0][2]*k1[2*i1]+g[0][3]*k1[2*i1+1]+g[0][4]*k1[(2*i1+2)%(2*j)];
    d[1][i1][i2]=h[1][0]*k2[(2*i1-2+2*j)%(2*j)]+h[1][1]*k2[(2*i1-1+2*j)%(2*j)]+h[1][2]*k2[2*i1]+h[1][3]*k2[2*i1+1]+h[1][4]*k2[(2*i1+2)%(2*j)];
    d[1][j+i1][i2]=g[1][1]*k2[(2*i1-1+2*j)%(2*j)]+g[1][2]*k2[2*i1]+g[1][3]*k2[2*i1+1]+g[1][4]*k2[(2*i1+2)%(2*j)];
   }
  }
  for (i1=0;((i1<4*j)&&(i1<N));i1++) {
   for (i2=0;i2<2*j;i2++) {
    k1[i2]=d[0][i1][i2];
    k2[i2]=d[1][i1][i2];
   }
   for (i2=0;i2<j;i2++) {
    d[0][i1][i2]=h[1][0]*k1[(2*i2-2+2*j)%(2*j)]+h[1][1]*k1[(2*i2-1+2*j)%(2*j)]+h[1][2]*k1[2*i2]+h[1][3]*k1[2*i2+1]+h[1][4]*k1[(2*i2+2)%(2*j)];
    d[0][i1][j+i2]=g[1][1]*k1[(2*i2-1+2*j)%(2*j)]+g[1][2]*k1[2*i2]+g[1][3]*k1[2*i2+1]+g[1][4]*k1[(2*i2+2)%(2*j)];
    d[1][i1][i2]=h[0][0]*k2[(2*i2-2+2*j)%(2*j)]+h[0][1]*k2[(2*i2-1+2*j)%(2*j)]+h[0][2]*k2[2*i2]+h[0][3]*k2[2*i2+1]+h[0][4]*k2[(2*i2+2)%(2*j)];
    d[1][i1][j+i2]=g[0][1]*k2[(2*i2-1+2*j)%(2*j)]+g[0][2]*k2[2*i2]+g[0][3]*k2[2*i2+1]+g[0][4]*k2[(2*i2+2)%(2*j)];
   }
  }
 }
 
 for (j=N/2;j>1;j=j/2) {
  for (j1=j;j1>j/4;j1=j1/2) {
   for (j2=j;j2>j/4;j2=j2/2) {
    for (i1=0;i1<j1;i1++) {
     for (i2=0;i2<j2;i2++) {
      ddn[0][j1+i1][j2+i2]=(j2*d[0][j1+i1][j2+i2]-j1*d[1][j1+i1][j2+i2])/(j1*j1+j2*j2);
      ddn[1][j1+i1][j2+i2]=(j1*d[0][j1+i1][j2+i2]+j2*d[1][j1+i1][j2+i2])/(j1*j1+j2*j2);
     }
    }
   }
  }
  for (i1=j;i1<2*j;i1++) {
   for (i2=0;i2<j/2;i2++) {
    ddn[0][i1][i2]=d[1][i1][i2]-d[0][i1][i2]/8.+d[0][i1][(i2-1+j/2)%(j/2)]/8.;
    ddn[1][i1][i2]=d[0][i1][i2]/8.;
    ddn[0][i2][i1]=d[0][i2][i1]-d[1][i2][i1]/8.+d[1][(i2-1+j/2)%(j/2)][i1]/8.;
    ddn[1][i2][i1]=d[1][i2][i1]/8.;
   }
  }
 }
 
 ddn[0][1][0]=d[1][1][0];
 ddn[1][1][0]=d[0][1][0];
 ddn[0][0][1]=d[0][0][1];
 ddn[1][0][1]=d[1][0][1];
 ddn[0][0][0]=d[0][0][0];
 ddn[1][0][0]=d[1][0][0];
}

/* Transformee inverse partant des ondelettes gradients */

void tfoganpinv2d(float c[2][N][N],float ddn[2][N][N])
{
 int j,j1,j2,i1,i2,e1,e2,l;
 float k1[N],k2[N];
 float h[2][5]={{1./2.,1.,1./2.,0.,0.},{1./4.,3./4.,3./4.,1./4.,0.}};
 float g[2][5]={{-1./4.,-1./2.,3./2.,-1./2.,-1./4.},{-1./2.,-3./2.,3./2.,1./2.,0.}};
 
 c[0][0][0]=ddn[0][0][0];
 c[1][0][0]=ddn[1][0][0];
 c[1][1][0]=ddn[0][1][0];
 c[0][1][0]=ddn[1][1][0];
 c[0][0][1]=ddn[0][0][1];
 c[1][0][1]=ddn[1][0][1];
 for (j=N/2;j>1;j=j/2) {
  for (j1=j;j1>j/4;j1=j1/2) {
   for (j2=j;j2>j/4;j2=j2/2) {
    for (i1=0;i1<j1;i1++) {
     for (i2=0;i2<j2;i2++) {
      c[0][j1+i1][j2+i2]=j2*ddn[0][j1+i1][j2+i2]+j1*ddn[1][j1+i1][j2+i2];
      c[1][j1+i1][j2+i2]=-j1*ddn[0][j1+i1][j2+i2]+j2*ddn[1][j1+i1][j2+i2];
     }
    }
   }
  }
  for (i1=j;i1<2*j;i1++) {
   for (i2=0;i2<j/2;i2++) {
    c[0][i1][i2]=8.*ddn[1][i1][i2];
    c[1][i1][i2]=ddn[0][i1][i2]+ddn[1][i1][i2]-ddn[1][i1][(i2-1+j/2)%(j/2)];
    c[0][i2][i1]=ddn[0][i2][i1]+ddn[1][i2][i1]-ddn[1][(i2-1+j/2)%(j/2)][i1];
    c[1][i2][i1]=8.*ddn[1][i2][i1];
   }
  }
 }
 
 for (i1=0;i1<N;i1++) {
  k1[i1]=0.;
  k2[i1]=0.;
 }
 
 for (j=1;j<N;j=2*j) {
  for (i2=0;((i2<4*j)&&(i2<N));i2++) {
   for (i1=0;i1<j;i1++) {
    for (l=-1;l<4;l++) {
     k1[(2*i1+l+2*j)%(2*j)]=h[0][l+1]*c[0][i1][i2]+g[0][l+1]*c[0][j+i1][i2]+k1[(2*i1+l+2*j)%(2*j)];
     k2[(2*i1+l+2*j)%(2*j)]=h[1][l+1]*c[1][i1][i2]+g[1][l+1]*c[1][j+i1][i2]+k2[(2*i1+l+2*j)%(2*j)];
    }
   }
   for (i1=0;i1<2*j;i1++) {
    c[0][i1][i2]=k1[i1];
    k1[i1]=0.;
    c[1][i1][i2]=k2[i1];
    k2[i1]=0.;
   }
  }
  
  for (i1=0;((i1<4*j)&&(i1<N));i1++) {
   for (i2=0;i2<j;i2++) {
    for (l=-1;l<4;l++) {
     k1[(2*i2+l+2*j)%(2*j)]=h[1][l+1]*c[0][i1][i2]+g[1][l+1]*c[0][i1][j+i2]+k1[(2*i2+l+2*j)%(2*j)];
     k2[(2*i2+l+2*j)%(2*j)]=h[0][l+1]*c[1][i1][i2]+g[0][l+1]*c[1][i1][j+i2]+k2[(2*i2+l+2*j)%(2*j)];
    }
   }
   for (i2=0;i2<2*j;i2++) {
    c[0][i1][i2]=k1[i2];
    k1[i2]=0.;
    c[1][i1][i2]=k2[i2];
    k2[i2]=0.;
   }
  }
 }

}

/* Decomposition de Helmholtz par ondelettes splines a divergence nulle et gradient */

void hodge(float vn[2][N][N]/*donnee*/,float d[2][N][N]/*resultat*/,float vn0[2][N][N]/*first guess*/,float d0[2][N][N]/*first guess*/,int It/*nb d'it*/)
{
 int j,i1,i2;
 float c[2][N][N],cn[2][N][N],dn[2][N][N];
 for (i2=0;i2<N;i2++) {
  for (i1=0;i1<N;i1++) {
   d[0][i1][i2]=d0[0][i1][i2];
   d[1][i1][i2]=d0[1][i1][i2];
   vn[0][i1][i2]=vn[0][i1][i2]-vn0[0][i1][i2];
   vn[1][i1][i2]=vn[1][i1][i2]-vn0[1][i1][i2];
  }
 }
 for (j=0;j<It;j++) {
  for (i2=0;i2<N;i2++) {
   for (i1=0;i1<N;i1++) {
   c[0][i1][i2]=5./8.*(vn[0][i1][i2]+vn[0][(i1+1)%N][i2])-1./8.*(vn[0][(i1-1+N)%N][i2]+vn[0][(i1+2)%N][i2]);
   c[1][i1][i2]=5./8.*(vn[1][i1][i2]+vn[1][i1][(i2+1)%N])-1./8.*(vn[1][i1][(i2-1+N)%N]+vn[1][i1][(i2+2)%N]);
   }
  }
  tfodanp2d(c,dn);
  for (i2=0;i2<N;i2++) {
   for (i1=0;i1<N;i1++) {
    d[0][i1][i2]=d[0][i1][i2]+dn[0][i1][i2];
    d[1][i1][i2]=d[1][i1][i2]+dn[1][i1][i2];
    if ((i1+i2)!=0) {
     d[1][i1][i2]=d[1][i1][i2]-dn[1][i1][i2];
     dn[1][i1][i2]=0.;
    }
   }
  }
  tfodanpinv2d(c,dn);
  for (i2=0;i2<N;i2++) {
   for (i1=0;i1<N;i1++) {
    vn[0][i1][i2]=vn[0][i1][i2]-1./2.*(c[0][(i1-1+N)%N][i2]+c[0][i1][i2]);
    vn[1][i1][i2]=vn[1][i1][i2]-1./2.*(c[1][i1][(i2-1+N)%N]+c[1][i1][i2]);
   }
  }
  for (i2=0;i2<N;i2++) {
   for (i1=0;i1<N;i1++) {
   c[0][i1][i2]=5./8.*(vn[0][i1][i2]+vn[0][i1][(i2+1)%N])-1./8.*(vn[0][i1][(i2-1+N)%N]+vn[0][i1][(i2+2)%N]);
   c[1][i1][i2]=5./8.*(vn[1][i1][i2]+vn[1][(i1+1)%N][i2])-1./8.*(vn[1][(i1-1+N)%N][i2]+vn[1][(i1+2)%N][i2]);
   }
  }
  tfoganp2d(c,dn);
  for (i2=0;i2<N;i2++) {
   for (i1=0;i1<N;i1++) {
    if ((i1+i2)!=0) {
     d[1][i1][i2]=d[1][i1][i2]+dn[1][i1][i2];
    }
    dn[0][i1][i2]=0.;
    if ((i1+i2)==0) dn[1][i1][i2]=0.;
   }
  }
  tfoganpinv2d(c,dn);
  for (i2=0;i2<N;i2++) {
   for (i1=0;i1<N;i1++) {
    vn[0][i1][i2]=vn[0][i1][i2]-1./2.*(c[0][i1][(i2-1+N)%N]+c[0][i1][i2]);
    vn[1][i1][i2]=vn[1][i1][i2]-1./2.*(c[1][(i1-1+N)%N][i2]+c[1][i1][i2]);
   }
  }
 }
 /*actualisation de la donnée initiale*/
 for (i2=0;i2<N;i2++) {
  for (i1=0;i1<N;i1++) {
   d0[0][i1][i2]=d[0][i1][i2];
   d0[1][i1][i2]=0.;
   dn[1][i1][i2]=d[1][i1][i2];
   dn[0][i1][i2]=0.;
  }
 }
 d0[1][0][0]=d[1][0][0];
 dn[1][0][0]=0.;
 tfodanpinv2d(c,d0);
 tfoganpinv2d(cn,dn);
 for (i2=0;i2<N;i2++) {
  for (i1=0;i1<N;i1++) {
   vn0[0][i1][i2]=0.5*(c[0][(i1-1+N)%N][i2]+c[0][i1][i2])+0.5*(cn[0][i1][(i2-1+N)%N]+cn[0][i1][i2]);
   vn0[1][i1][i2]=0.5*(c[1][i1][(i2-1+N)%N]+c[1][i1][i2])+0.5*(cn[1][(i1-1+N)%N][i2]+cn[1][i1][i2]);
   d0[0][i1][i2]=d[0][i1][i2];
   d0[1][i1][i2]=d[1][i1][i2];
  }
 }
}

float v[2][N][N],vn[2][N][N],u[2][N][N],c[2][N][N],cn[2][N][N],d[2][N][N],dn[2][N][N],mx,x,y,dx=1./N,mxo,mxd;
float omega[N][N],diver[N][N];
int i,i1,i2,j,l1,l2,k1,k2,seed=101,R=10,t=12,R1,R2,im,cr,cb,cg;
int dim[3]={0,N,N};
/* regularite de la partie rot :*/ float ar=3.1;
/* regularite de la partie div :*/ float ad=3.1;

int main()
{
 FILE *s=fopen("echant.h","w");
 FILE *sr=fopen("echantr.h","w");
 FILE *sd=fopen("echantd.h","w");
 FILE *dr=fopen("omega.ppm","w");
 FILE *dd=fopen("diver.ppm","w");
 
 fprintf(dr,"P6\n%d %d\n255\n",N,N);
 fprintf(dd,"P6\n%d %d\n255\n",N,N);
 
 srandom(seed);
 for (i2=0;i2<N;i2++) {
  for (i1=0;i1<N;i1++) {
   d[0][i1][i2]=RANDN;
   d[1][i1][i2]=0.;
   dn[1][i1][i2]=RANDN;
   dn[0][i1][i2]=0.;
  }
 }
 d[0][0][0]=0.;
 dn[1][i1][i2]=0.;
 
 /* reglage de la regularite */
 for (j=N/2;j>1;j=j/2) {
  for (l1=j;l1>j/4;l1=l1/2) {
   for (l2=j;l2>j/4;l2=l2/2) {
    for (i1=0;i1<l1;i1++) {
     for (i2=0;i2<l2;i2++) {
      d[0][l1+i1][l2+i2]=((float) pow((double) j,(double) -ar))*(l1+l2)*d[0][l1+i1][l2+i2];
      dn[1][l1+i1][l2+i2]=((float) pow((double) j,(double) -ad))*(l1+l2)*dn[1][l1+i1][l2+i2];
     }
    }
   }
  }
  for (i1=j;i1<2*j;i1++) {
   for (i2=0;i2<j/2;i2++) {
    d[0][i1][i2]=((float) pow((double) j,(double) -ar))*d[0][i1][i2]/9;
    dn[1][i1][i2]=((float) pow((double) j,(double) -ad))*dn[1][i1][i2]/9;
    d[0][i2][i1]=((float) pow((double) j,(double) -ar))*d[0][i2][i1]/9;
    dn[1][i2][i1]=((float) pow((double) j,(double) -ad))*dn[1][i2][i1]/9;
   }
  }
 }
 
 tfodanpinv2d(c,d);
 tfoganpinv2d(cn,dn);

 for (i2=0;i2<N;i2++) {
  for (i1=0;i1<N;i1++) {
   v[0][i1][i2]=0.5*(c[0][(i1-1+N)%N][i2]+c[0][i1][i2]);
   v[1][i1][i2]=0.5*(c[1][i1][(i2-1+N)%N]+c[1][i1][i2]);
   vn[0][i1][i2]=0.5*(cn[0][i1][(i2-1+N)%N]+cn[0][i1][i2]);
   vn[1][i1][i2]=0.5*(cn[1][(i1-1+N)%N][i2]+cn[1][i1][i2]);
  }
 }

 /* Representation */
 mx=0.;
 for (i2=0;i2<N;i2++) {
  for (i1=0;i1<N;i1++) {
   omega[i1][i2]=(v[1][(i1+1)%N][i2]-v[1][i1][i2])/dx-(v[0][i1][(i2+1)%N]-v[0][i1][i2])/dx;
   mx=(fabs(mx-fabs(omega[i1][i2]))+mx+fabs(omega[i1][i2]))/2;
  }
 }
 mxo=mx;
 printf("mx(omega)=%f\n",mx);
 mx=0.;
 for (i2=0;i2<N;i2++) {
  for (i1=0;i1<N;i1++) {
   diver[i1][i2]=(vn[0][(i1+1)%N][i2]-vn[0][i1][i2])/dx+(vn[1][i1][(i2+1)%N]-vn[1][i1][i2])/dx;
   mx=(fabs(mx-fabs(diver[i1][i2]))+mx+fabs(diver[i1][i2]))/2;
  }
 }
 mxd=mx;
 printf("mx(diver)=%f\n",mx);
 
  for (j=0;j<N;j++) {
   for (i=0;i<N;i++) {
    im=(int) ((-omega[i][N-1-j])/mxo*255);
    cr=(int) (255+(im-abs(im))/2);
    cb=(int) (255-(im+abs(im))/2);
    cg=(cr+cb-abs(cr-cb))/2;
    cg=(255+cg-abs(255-cg))/2;
    fprintf(dr,"%c%c%c", (unsigned char) cr,(unsigned char) cg,(unsigned char) cb);
   }
  }
  fclose(dr);

  for (j=0;j<N;j++) {
   for (i=0;i<N;i++) {
    im=(int) ((diver[i][N-1-j])/mxo*255);
    cr=(int) (255+(im-abs(im))/2);
    cb=(int) (255-(im+abs(im))/2);
    cg=(cr+cb-abs(cr-cb))/2;
    cg=(255+cg-abs(255-cg))/2;
    fprintf(dd,"%c%c%c", (unsigned char) cr,(unsigned char) cg,(unsigned char) cb);
   }
  }
  fclose(dd);

 for (i2=0;i2<N;i2++)  /* ou (i2=0;i2<N/2;i2++) pour briser la periodicite */{
  for (i1=0;i1<N;i1++) /* ou (i2=0;i2<N/2;i2++) pour briser la periodicite */{
   u[0][i1][i2]=v[0][i1][i2]+vn[0][i1][i2];
   u[1][i1][i2]=v[1][i1][i2]+vn[1][i1][i2];
   fprintf(s,"%f %f\n",u[0][i1][i2],u[1][i1][i2]);
   fprintf(sr,"%f %f\n",v[0][i1][i2],v[1][i1][i2]);
   fprintf(sd,"%f %f\n",vn[0][i1][i2],vn[1][i1][i2]);
  }
 }

 
}






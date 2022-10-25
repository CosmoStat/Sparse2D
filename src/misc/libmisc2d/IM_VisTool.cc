

#include "IM_Obj.h"
#include "IM_VisTool.h"

//----------------------------------------------------------
//	decime_2 : decime un pixel sur 2 de Im
//		et place le resultat dans Im2 allouee
//----------------------------------------------------------
Ifloat *		decime_2 (Ifloat & Im)
{
	int	i, j;
	float	*l, *lr;
	Ifloat	*rep = new Ifloat (Im.nl() >> 1, Im.nc() >> 1, "decime 2");

	lr = rep->buffer();
	for (j=0 ; j<rep->nl() ; j++)
	{
                l = Im.buffer() + Im.nc()*j*2;
		for (i=0 ; i<rep->nc() ; i++, l+=2, lr++)
			*lr = *l;
	}
	return rep;
}


//----------------------------------------------------------
//	zoom_2
//		grossi l'image Im par 2
//----------------------------------------------------------
Ifloat *		zoom_2 (Ifloat & Im, int dx, int dy)
{
	Ifloat *rep;
	int		i, j;
	float		*l1, *l2, *l;

	rep = new Ifloat ((Im.nl()<<1) + dy, (Im.nc()<<1) + dx, "zoom 2");
	l1 = rep->buffer();
	l2 = rep->buffer() + rep->nc();
	l = Im.buffer();
	for (j=0 ; j<Im.nl() ; j++)
	{
		for (i=0 ; i<Im.nc() ; i++, l++, l1++, l2++)
		{
			*l1 = *l2 = *l;
			*++l1 = *++l2 = *l;
		}
		l1 += dx;
		l2 += dx;
		l1 += rep->nc();
		l2 += rep->nc();
	}

	complete_bord (*rep, dx, dy);
	return rep;
}


//----------------------------------------------------------
//	complete_bord
//----------------------------------------------------------

void complete_bord (Ifloat & Ima, int dx, int dy)
{
	float		*f, *b;
	int		i;

	if (dx)
	{
		f = Ima.buffer() + Ima.nc() - 2;
		for (i=0 ; i<Ima.nl() ; i++, f+=Ima.nc())
			f[1] = f[0];
	}
	if (dy)
	{
		f = Ima.buffer() + Ima.nc()*(Ima.nl()-1);
		b = Ima.buffer() + Ima.nc()*(Ima.nl()-2);
		for (i=0 ; i<Ima.nc() ; i++, f++, b++)
			*f = *b;
	}
}

//----------------------------------------------------------
//	masque
//----------------------------------------------------------
void masque (Ifloat & Ima, Ifloat& Im)
{
   float *data = Ima.buffer();
   float *ptr = Im.buffer();

	assert (Ima.nl()==Im.nl() && Ima.nc()==Im.nc());

	for (int i=0 ; i<Ima.nl()*Ima.nc() ; i++)

		if (ptr[i] == 0.)
			data[i] = 0.;
}

//----------------------------------------------------------
//	masque
//----------------------------------------------------------
void masque (Ifloat &Ima, Ifloat& Im, int dx, int dy)
{
        int Nl = Ima.nl();
        int Nc = Ima.nc();
	Ifloat *temp = new Ifloat (Nl, Nc, "masque 2");
	temp->init (0.0);
	int	i, j, x1, y1, x2, y2, c, l;
	x1 = dx;
	y1 = dy;
	x2 = dx + Im.nc();
	y2 = dy + Im.nl();

	x1 = MIN (Nc, x1);
	x1 = MAX (0, x1);
	y1 = MIN (Nl, y1);
	y1 = MAX (0, y1);

	x2 = MIN (Nc, x2);
	x2 = MAX (0, x2);
	y2 = MIN (Nl, y2);
	y2 = MAX (0, y2);

	for (j=y1, l=0 ; j<y2 ; j++, l++)
		for (i=x1, c=0 ; i<x2 ; i++, c++)
			(*temp)(j,i) = Im(l,c);
	masque (Ima, *temp);
	delete temp;
}

//----------------------------------------------------------
//	first_pixel:
//----------------------------------------------------------
int first_pixel (Ifloat &Ima, int & x, int & y)
{
        float *data = Ima.buffer();
	int Nc = Ima.nc();

	x = y = -1;
	for (int i=0 ; i<Ima.nl()*Ima.nc() ; i++)
		if (data[i])
		{
			x = i % Nc;
			y = i / Nc;
			return i;
		}
	return -1;
}

//----------------------------------------------------------
//	add_image
//----------------------------------------------------------
void add_image (Ifloat & Ima1, Ifloat& Im, int dx, int dy)
{
	int	i, j, x1, y1, x2, y2, c, l;
        int Nl = Ima1.nl();
        int Nc = Ima1.nc();
	x1 = dx;
	y1 = dy;
	x2 = dx + Im.nc();
	y2 = dy + Im.nl();
        
	x1 = MIN (Nc, x1);
	x1 = MAX (0, x1);
	y1 = MIN (Nl, y1);
	y1 = MAX (0, y1);

	x2 = MIN (Nc, x2);
	x2 = MAX (0, x2);
	y2 = MIN (Nl, y2);
	y2 = MAX (0, y2);

	for (j=y1, l=y1-dy ; j<y2 ; j++, l++)
		for (i=x1, c=x1-dx ; i<x2 ; i++, c++)
			Ima1(j,i) += Im(l,c);
}

//----------------------------------------------------------
//	moments: calcul les moments d'ordre 1 et 2  
//		de l'image exprimes en pixels
//----------------------------------------------------------
t_Moment_image moments(Ifloat &Ima)
{
	int	x,y;
        double data;
	t_Moment_image mom;
	double	M=0.;
	double	mx=0.;
	double	my=0.;
	double	vx=0.;
	double	vy=0.;
	double	vxy=0.;
        int Nc = Ima.nc();
        int Nl = Ima.nl();

	for(x=0;x<Nc;x++)
	{
		for(y=0;y<Nl;y++)
		{
			data=Ima(y,x);
			M += data;
			
			mx+=x*data;
			my+=y*data;

			vx+=x*x*data;
			vy+=y*y*data;

			vxy+=x*y*data;
		}
	}

	mom.mx=(mx/=M);
	mom.my=(my/=M);
	vx/=M;
	vy/=M;
	vxy/=M;

	mom.vx=vx-=mx*mx;
	mom.vy=vy-=my*my;
	mom.vxy=vxy-=mx*my;

	return mom;
}

//----------------------------------------------------------
//	Convol_5_trou
//		convolution d'une image par le filtre m ^ m
//		par l'algo a trou
//		plan est le niveau du dilate
//		m est la matrice du filtre en ligne
//		exemple:
//			m[0] = m[4] = 1/16
//			m[1] = m[3] = 4/16;
//			m[2] = 6/16
//----------------------------------------------------------

void convol_5_trou (Ifloat & Im, Ifloat & Im_out, int plan)
{
	int	pas, pas2, x, y, i, j;
	float	*l, *lr, f;
        float m[5];
	m[0] = m[4] = 1./16;
        m[1] = m[3] = 4./16;
	m[2] = 6./16;
	pas = 1 << plan;
	pas2 = pas << 1;
        int Nl = Im.nl();
        int Nc = Im.nc();
	Ifloat rep1( Nl, Nc, "convol 5 trou");

	//-----------------------------------------
	//	Traitement en ligne de l'image
	//-----------------------------------------
	lr = rep1.buffer();
	for (j=0 ; j< Nl ; j++)
	{
		l = Im.buffer() + j*Nc; 
		for (i=0 ; i<pas ; i++, lr++)
		{
			x = i -pas2;
			f = m[0] * l[-x];
			x += pas;
			f += m[1] * l[-x];
			x += pas;
			f += m[2] * l[x];
			x += pas;
			f += m[3] * l[x];
			x += pas;
			f += m[4] * l[x];
			*lr = f;
		}
		for ( ; i<pas2 ; i++, lr++)
		{
			x = i -pas2;
			f = m[0] * l[-x];
			x += pas;
			f += m[1] * l[x];
			x += pas;
			f += m[2] * l[x];
			x += pas;
			f += m[3] * l[x];
			x += pas;
			f += m[4] * l[x];
			*lr = f;
		}
		for ( ; i<Nc-pas2 ; i++, lr++)
		{
			x = i -pas2;
			f = m[0] * l[x];
			x += pas;
			f += m[1] * l[x];
			x += pas;
			f += m[2] * l[x];
			x += pas;
			f += m[3] * l[x];
			x += pas;
			f += m[4] * l[x];
			*lr = f;
		}
		for ( ; i<Nc-pas ; i++, lr++)
		{
			x = i - pas2;
			f = m[0] * l[x];
			x += pas;
			f += m[1] * l[x];
			x += pas;
			f += m[2] * l[x];
			x += pas;
			f += m[3] *l[x];
			x += pas;
			f += m[4] * l[x-2*(x%(Nc-1))-1];
			*lr = f;
		}
		for ( ; i<Nc ; i++, lr++)
		{
			x = i - pas2;
			f = m[0] * l[x];
			x += pas;
			f += m[1] * l[x];
			x += pas;
			f += m[2] * l[x];
			x += pas;
			f += m[3] * l[x-2*(x%(Nc-1))-1];
			x += pas;
			f += m[4] * l[x-2*(x%( Nc-1))-1];
			*lr = f;
		}
	}


	//-----------------------------------------
	//	Traitement en colomne de l'image
	//-----------------------------------------
	lr = Im_out.buffer();
	for (j=0 ; j<pas ; j++)
		for (i=0 ; i< Nc ; i++, lr++)
		{
			y = j - pas2;
			f = m[0] * rep1(-y,i);
			y += pas;
			f += m[1] * rep1(-y,i);
			y += pas;
			f += m[2] * rep1(y,i);
			y += pas;
			f += m[3] * rep1(y,i);
			y += pas;
			f += m[4] * rep1(y,i);
			*lr = f;
		}
	for ( ; j<pas2 ; j++)
		for (i=0 ; i< Nc ; i++, lr++)
		{
			y = j - pas2;
			f = m[0] * rep1(-y,i);
			y += pas;
			f += m[1] * rep1(y,i);
			y += pas;
			f += m[2] * rep1(y,i);
			y += pas;
			f += m[3] * rep1(y,i);
			y += pas;
			f += m[4] * rep1(y,i);
			*lr = f;
		}
	for ( ; j< Nl-pas2 ; j++)
		for (i=0 ; i< Nc ; i++, lr++)
		{
			y = j - pas2;
			f = m[0] * rep1(y,i);
			y += pas;
			f += m[1] * rep1(y,i);
			y += pas;
			f += m[2] * rep1(y,i);
			y += pas;
			f += m[3] * rep1(y,i);
			y += pas;
			f += m[4] * rep1(y,i);
			*lr = f;
		}
	for ( ; j< Nl-pas ; j++)
		for (i=0 ; i< Nc ; i++, lr++)
		{
			y = j - pas2;
			f = m[0] * rep1(y,i);
			y += pas;
			f += m[1] * rep1(y,i);
			y += pas;
			f += m[2] * rep1(y,i);
			y += pas;
			f += m[3] * rep1(y,i);
			y += pas;
			f += m[4] * rep1(y-2*(y%(Nl-1))-1,i);
			*lr = f;
		}
	for ( ; j< Nl ; j++)
		for (i=0 ; i< Nc ; i++, lr++)
		{
			y = j - pas2;
			f = m[0] * rep1(y,i);
			y += pas;
			f += m[1] * rep1(y,i);
			y += pas;
			f += m[2] * rep1(y,i);
			y += pas;
			f += m[3] * rep1(y-2*(y%(Nl-1))-1,i);
			y += pas;
			f += m[4] * rep1(y-2*(y%(Nl-1))-1,i);
			*lr = f;
		}
}




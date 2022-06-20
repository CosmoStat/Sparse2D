
/******************************************************************************
**                   Copyright (C) 1995 CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 3.1
**
**    Author: J.L. Starck
**
**    Date:  96/05/02 
**    
**    File:  IM_Code.cc
**
*******************************************************************************
**
**    DECRIPTION    ROUTINES for image compression
**    ---------- 
**
*****************************************************************************/

// static char sccsid[] = "@(#)IM_Code.cc 3.1 96/05/02 CEA 1995 @(#)";

#include"IM_Obj.h"
#include"IM_CompTool.h"

#define input_nybble(infile)	input_nbits(infile,4)
#define output_nybble(outfile,c) output_nbits(outfile,c,4)

static char code_magic[2] = { (char)0xDD, (char)0x99 };
extern int bitcount;
static int DebugCode=0;
long size_enc;
/********************************************************************/


void encode(FILE *outfile, int *a, int nx, int ny, float scale)
{
/* outfile=	Output file	
   a[]=  input  array (nx,ny)	
   nx,ny= size of H-transform array		
   scale= scale factor for digitization	*/

	int nel, i, j, k, vmax, nsign, bits_to_go;
	unsigned char nbitplanes;
	unsigned char *signbits=NULL;

	nel = nx*ny;
	writeint(outfile, 0);
        int Bef = ftell(outfile);        
	writeint(outfile, nx);	/* size of image */
	writeint(outfile, ny);
	writefloat(outfile, scale); /* scale factor for digitization	*/

	/* allocate array for sign bits and save values, 8 per byte */
	signbits = new unsigned char [(nel+7)/8+1];
	if (signbits == (unsigned char *) NULL) 
	{
		fprintf(stderr, "encode: insufficient memory\n");
		exit(-1);
	}
	nsign = 0;
	bits_to_go = 8;
	signbits[0] = 0;
	for (i=0; i<nel; i++) 
        {
		if (a[i] > 0)
                {
			/*
			 * positive element, put zero at end of buffer
			 */
			signbits[nsign] <<= 1;
			bits_to_go -= 1;
		} 
		else if (a[i] < 0) 
                {
			/*
			 * negative element, shift in a one
			 */
			signbits[nsign] <<= 1;
			signbits[nsign] |= 1;
			bits_to_go -= 1;
			/*
			 * replace a by absolute value
			 */
			a[i] = -a[i];
		}
		if (bits_to_go == 0) 
		{
			/*
			 * filled up this byte, go to the next one
			 */
			bits_to_go = 8;
			nsign += 1;
			signbits[nsign] = 0;
		}
	}
	if (bits_to_go != 8) 
        {
		/*
		 * some bits in last element
		 * move bits in last byte to bottom and increment nsign
		 */
		signbits[nsign] <<= bits_to_go;
		nsign += 1;
	}
	/* calculate number of bit planes image	 */
	vmax = 0;

	/* get maximum absolute value in image */
	j=0;	/* column counter	*/
	k=0;	/* row counter		*/
	for (i=0; i<nel; i++) 
        {
		if (vmax < a[i]) vmax = a[i];
		if (++j >= ny) 
                {
			j = 0;
			k += 1;
		}
	}
	/* now calculate number of bits for image */
	nbitplanes = (unsigned char) (log((float) (vmax+1))/log(2.0)+0.5);
	if ( (vmax+1) > (1<<nbitplanes) ) 
        {
	    nbitplanes += 1;
	}

	/* write nbitplanes */
	qwrite(outfile, (char *) &nbitplanes, sizeof(nbitplanes));

	/* Initialize bit output */
	start_outputing_bits();
	qtree_encode(outfile, &a[0], ny, nx, ny, nbitplanes);
	/* Add zero as an EOF symbol */
	output_nybble(outfile, 0);
	done_outputing_bits(outfile);

	/* write sign bits */
	if (nsign > 0) qwrite(outfile, (char *) signbits, nsign);
	if (DebugCode) 
        {
                // fprintf(stderr, "nbitplanes = %d\n", nbitplanes);
		/* total number of bits written to file */
		i = bitcount  + 8 * (nsign + sizeof(code_magic) 
                                   +sizeof(nbitplanes) + 4*sizeof(int));
               //fprintf(stderr, "%6.3f bits/pixel, compression factor %5.1f\n",
		//	((float) i)/nel,
		//	16.0*nel/((float) i));
	}
	if (signbits != NULL) delete [] signbits;        
	int Aft = ftell(outfile);
	size_enc =  Aft - Bef;
	long Val = -size_enc-4;
	int nf = fseek(outfile, Val, SEEK_CUR);
        writeint(outfile, size_enc);
        nf = fseek(outfile, size_enc,SEEK_CUR);
}

/********************************************************************/

void dodecode(FILE *infile, int *a, int nx, int ny, 
                     unsigned char nbitplanes)
/*FILE *infile;
int a[];  Array of values to decode	
int nx,ny;  Array dimensions are [nx][ny] 
unsigned char nbitplanes;  Number of bit planes  */
{
    int i, nel;

    nel = nx*ny;

    /* initialize a to zero */
    for (i=0; i<nel; i++) a[i] = 0;

    /* Initialize bit input */
    start_inputing_bits();

    /* read bit planes for each quadrant */
    qtree_decode(infile, &a[0], ny, nx, ny,  nbitplanes);

    /* make sure there is an EOF symbol (nybble=0) at end */
    if (input_nybble(infile) != 0) 
    {
        fprintf(stderr, "dodecode: bad bit plane values\n");
        exit(-1);
    }

    /* now get the sign bits Re-initialize bit input */
    start_inputing_bits();
    for (i=0; i<nel; i++) 
    {
        if (a[i] != 0) 
        {
            if (input_bit(infile) != 0) a[i] = -a[i];
	}
    }
}

/********************************************************************/

void decode(FILE *infile, int  **a, int  *nx, int  *ny, float  *scale)
/*FILE *infile;  input file 
int  **a;  address of output array [nx][ny]	
int  *nx,*ny;	 size of output array	
int  *scale;  scale factor for digitization */
{
    int nel;
    unsigned char nbitplanes;
    size_enc = readint(infile);
    *nx =readint(infile); /* x size of image */
    *ny =readint(infile); /* y size of image */
    *scale=readfloat(infile); /* scale factor for digitization */

    /* allocate memory for array */
    nel = (*nx) * (*ny);
    *a = i_alloc(nel);
    if (*a == (int *) NULL)
    {
        fprintf(stderr, "decode: insufficient memory\n");
        exit(-1);
    }

    /* # bits in quadrants */
    qread(infile, (char *) &nbitplanes, sizeof(nbitplanes));
    dodecode(infile, *a, *nx, *ny, nbitplanes);
}

/********************************************************************/

/******************************************************************************
**                   Copyright (C) 1994 CEA
*******************************************************************************
**
**    UNIT
**
**    Version: 3.1
**
**    Author: R. White & J.L. Starck
**
**    Date:  96/05/02 
**    
**    File:  IM_QuadTree.cc
**
*******************************************************************************
**
**    DECRIPTION  input-output routines for reading and witing coded values
**    ---------- 
**
*******************************************************************************
**
** void qtree_encode(FILE *outfile, int a[], int n, int nqx, 
**                   iny nqy,int nbitplanes)
**
** outfile: IN = file where we put the coded data
** a: IN = image to code
** n: IN = physical dimension of row in a
** nqx: IN =  length of row 
** nqy: IN = length of column (<=n)
** nbitplanes: IN =  number of bit planes to output
**
** encodes an image by doing first a quad tree, and then applying
** a huffmann coding
**
*********************************************************************************
** void qtree_decode(FILE *infile, int a[], int n, int nqx, 
**                   int nqy, int nbitplanes)
** 
** infile: IN = file where we get the coded data
** a: OUT = decoded image  
** n: IN =  length of full row in a
** nqx: IN = partial length of row to decode   
** nqy: IN = partial length of column (<=n) 
** nbitplanes: IN =  number of bit planes to decode
** 
** decodes an image
**
*****************************************************************************/

// static char sccsid[] = "@(#)IM_QuadTree.cc 3.1 96/05/02 CEA 1994 @(#)";

#include "IM_Obj.h"
#include "IM_CompTool.h"

/* Huffman code values and number of bits in each code */
static int code[16] =
	{
	0x3e, 0x00, 0x01, 0x08, 0x02, 0x09, 0x1a, 0x1b,
	0x03, 0x1c, 0x0a, 0x1d, 0x0b, 0x1e, 0x3f, 0x0c
	};
static int ncode[16] =
	{
	6,    3,    3,    4,    3,    4,    5,    5,
	3,    5,    4,    5,    4,    5,    6,    4
	};

/* variables for bit output to buffer when Huffman coding */
static int bitbuffer, bits_to_go;

/* macros to write out 4-bit nybble, Huffman code for this value */
#define output_nybble(outfile,c)	output_nbits(outfile,c,4)
#define output_huffman(outfile,c)	output_nbits(outfile,code[c],ncode[c])
#define input_nybble(infile)    input_nbits(infile,4)


static void qtree_reduce(unsigned char a[], int n, int nx,
                         int ny, unsigned char b[]);

static void write_bdirect(FILE *outfile, int a[], int n, int nqx,
                          int nqy, unsigned char scratch[], int bit);

static int bufcopy(unsigned char a[], int n, unsigned char buffer[],
                   int *b, int bmax);
static void qtree_onebit(int a[], int n, int nx, int ny, 
                         unsigned char b[], int bit);

static void qtree_expand(FILE *infile, unsigned char a[],
                         int nx, int ny, unsigned char b[]);

static void qtree_bitins(unsigned char a[], int nx, int ny, 
                         int b[], int n, int bit);

static void qtree_copy(unsigned char a[], int nx, int ny, unsigned char b[],
                       int n);

static void read_bdirect(FILE *infile, int a[], int n, int nqx,
                         int nqy, unsigned char scratch[], int bit);

static int input_huffman(FILE *infile);

/********************************************************************/

void qtree_encode(FILE *outfile, int a[], int n, 
                         int nqx, int nqy, int nbitplanes)
/* FILE *outfile;
int a[];
int n;		physical dimension of row in a		
int nqx;	 length of row				
int nqy;         length of column (<=n)			
int nbitplanes;	 number of bit planes to output		*/
{
int log2n, i, k, bit, b, bmax, nqmax, nqx2, nqy2, nx, ny;
unsigned char *scratch, *buffer;

	/* log2n is log2 of max(nqx,nqy) rounded up to next power of 2 */
	nqmax = (nqx>nqy) ? nqx : nqy;
	log2n = (int) (log((float) nqmax)/log(2.0)+0.5);
	if (nqmax > (1<<log2n))  log2n += 1;

	/* initialize buffer point, max buffer size */
	nqx2 = (nqx+1)/2;
	nqy2 = (nqy+1)/2;
	bmax = (nqx2*nqy2+1)/2;

	/*
	 * We're indexing A as a 2-D array with dimensions (nqx,nqy).
	 * Scratch is 2-D with dimensions (nqx/2,nqy/2) rounded up.
	 * Buffer is used to store string of codes for output.
	 */
	scratch = (unsigned char *) malloc(2*bmax);
	buffer = (unsigned char *) malloc(bmax);
	if ((scratch == (unsigned char *) NULL) ||
	    (buffer  == (unsigned char *) NULL)) 
        {
	    fprintf(stderr, "qtree_encode: insufficient memory\n");
            exit(-1);
	}

	/* now encode each bit plane, starting with the top */
	for (bit=nbitplanes-1; bit >= 0; bit--) 
        {
		/* initial bit buffer */
		b = 0;
		bitbuffer = 0;
		bits_to_go = 0;
		/* on first pass copy A to scratch array */
		qtree_onebit(a,n,nqx,nqy,scratch,bit);
		nx = (nqx+1)>>1;
		ny = (nqy+1)>>1;
		/*
		 * copy non-zero values to output buffer, which will be written
		 * in reverse order
		 */
		if (bufcopy(scratch,nx*ny,buffer,&b,bmax)) {
			/*
			 * quadtree is expanding data,
			 * change warning code and just fill buffer with bit-map
			 */
			write_bdirect(outfile,a,n,nqx,nqy,scratch,bit);
			goto bitplane_done;
		}
		/*
		 * do log2n reductions
		 */
		for (k = 1; k<log2n; k++) {
			qtree_reduce(scratch,ny,nx,ny,scratch);
			nx = (nx+1)>>1;
			ny = (ny+1)>>1;
			if (bufcopy(scratch,nx*ny,buffer,&b,bmax)) {
				write_bdirect(outfile,a,n,nqx,nqy,scratch,bit);
				goto bitplane_done;
			}
		}
		/*
		 * OK, we've got the code in buffer
		 * Write quadtree warning code, then
		 * write buffer in reverse order
		 */
		output_nybble(outfile,0xF);
		if (b==0) 
                {
		   if (bits_to_go>0) 
		   {
		      /* put out the last few bits */
		      output_nbits(outfile, bitbuffer & ((1<<bits_to_go)-1),
				   bits_to_go);
		   } 
		   else 
		   {
		      /* write a zero nybble if there are no 1's in array */
		      output_huffman(outfile,0);
		   }
		 } 
		 else
		 {
		    if (bits_to_go>0) 
		    {
		       /* put out the last few bits */
		       output_nbits(outfile, bitbuffer & ((1<<bits_to_go)-1),
				    bits_to_go);
		    }
		    for (i=b-1; i>=0; i--) 
		    {
		       output_nbits(outfile,buffer[i],8);
		    }
		 }
		 bitplane_done: ;
	}
	free(buffer);
	free(scratch);
}

/********************************************************************/

/* copy non-zero codes from array to buffer */
static int bufcopy(unsigned char a[], int n, unsigned char buffer[],
                   int *b, int bmax)
{
    int i;

    for (i = 0; i < n; i++)
    {
        if (a[i] != 0) 
        {
           /* add Huffman code for a[i] to buffer */
           bitbuffer |= code[a[i]] << bits_to_go;
           bits_to_go += ncode[a[i]];
           if (bits_to_go >= 8)
           {
               buffer[*b] = bitbuffer & 0xFF;
	       *b += 1;
	       /* return warning code if we fill buffer */
               if (*b >= bmax) return(1);
               bitbuffer >>= 8;
               bits_to_go -= 8;
	    }
	 }
     }
    return(0);
}

/********************************************************************/

/* Do first quadtree reduction step on bit BIT of array A.
 * Results put into B. 
 */
static void qtree_onebit(int a[], int n, int nx, int ny, 
                         unsigned char b[], int bit)
{
    int i, j, k;
    int b0, b1, b2, b3;
    int s10, s00;
 
    /* use selected bit to get amount to shift */
    b0 = 1<<bit;
    b1 = b0<<1;
    b2 = b0<<2;
    b3 = b0<<3;
    k = 0;  /* k is index of b[i/2,j/2]	*/
    for (i = 0; i<nx-1; i += 2) 
    {
        s00 = n*i;   /* s00 is index of a[i,j] */
        s10 = s00+n; /* s10 is index of a[i+1,j] */
        for (j = 0; j<ny-1; j += 2) 
        {
            b[k] = ( ( a[s10+1]     & b0)
		   | ((a[s10  ]<<1) & b1)
		   | ((a[s00+1]<<2) & b2)
		   | ((a[s00  ]<<3) & b3) ) >> bit;
            k += 1;
            s00 += 2;
            s10 += 2;
	}
        if (j < ny)
        {
	   /*
	    * row size is odd, do last element in row
	    * s00+1,s10+1 are off edge
	    */
            b[k] = ( ((a[s10  ]<<1) & b1)
		   | ((a[s00  ]<<3) & b3) ) >> bit;
            k += 1;
	}
    }
    if (i < nx) 
    {
	/*
	 * column size is odd, do last row
	 * s10,s10+1 are off edge
	 */
         s00 = n*i;
         for (j = 0; j<ny-1; j += 2)
         {
             b[k] = ( ((a[s00+1]<<2) & b2)
		    | ((a[s00  ]<<3) & b3) ) >> bit;
             k += 1;
             s00 += 2;
	 }
         if (j < ny) 
         {
	    /*
	     * both row and column size are odd, do corner element
	     * s00+1, s10, s10+1 are off edge
	     */
             b[k] = ( ((a[s00  ]<<3) & b3) ) >> bit;
             k += 1;
	 }
    }
}

/********************************************************************/

/*
 * do one quadtree reduction step on array a
 * results put into b (which may be the same as a)
 */
static void qtree_reduce(unsigned char a[], int n, int nx,
                         int ny, unsigned char b[])
{
    int i, j, k;
    int s10, s00;

    k = 0; /* k is index of b[i/2,j/2]	*/
    for (i = 0; i<nx-1; i += 2) 
    {
        s00 = n*i;      /* s00 is index of a[i,j] */
        s10 = s00+n;	/* s10 is index of a[i+1,j] */
        for (j = 0; j<ny-1; j += 2) 
        {
            b[k] = (a[s10+1] != 0)
		  | ( (a[s10  ] != 0) << 1)
		  | ( (a[s00+1] != 0) << 2)
		  | ( (a[s00  ] != 0) << 3);
            k += 1;
            s00 += 2;
            s10 += 2;
	}
        if (j < ny) 
        {
	   /*
	   * row size is odd, do last element in row
	   * s00+1,s10+1 are off edge
	   */
	    b[k] =  ( (a[s10  ] != 0) << 1)
		    | ( (a[s00  ] != 0) << 3);
	    k += 1;
	}
    }
    if (i < nx) 
    {
	   /*
	   * column size is odd, do last row
	   * s10,s10+1 are off edge
	   */
	    s00 = n*i;
	    for (j = 0; j<ny-1; j += 2) 
	    {
	        b[k] =  ( (a[s00+1] != 0) << 2) | ( (a[s00  ] != 0) << 3);
	        k += 1;
	        s00 += 2;
	    }
	    if (j < ny)
	    {
	        /*
		 * both row and column size are odd, do corner element
		 * s00+1, s10, s10+1 are off edge
		 */
	        b[k] = ( (a[s00  ] != 0) << 3);
	        k += 1;
	    }
    }
}

/********************************************************************/

static void write_bdirect(FILE *outfile, int a[], int n, int nqx,
                          int nqy, unsigned char scratch[], int bit)
{
int i;

	/* Write the direct bitmap warning code	 */
	output_nybble(outfile,0x0);

	/* Copy A to scratch array (again!), packing 4 bits/nybble */
	qtree_onebit(a,n,nqx,nqy,scratch,bit);

	/* write to outfile */
	for (i = 0; i < ((nqx+1)/2) * ((nqy+1)/2); i++) 
        {
	    output_nybble(outfile,scratch[i]);
	}
}

/********************************************************************/

void qtree_decode(FILE *infile, int a[], int n, int nqx, 
                  int nqy, int nbitplanes)
/* FILE *infile;
int a[];  a is 2-D array with dimensions (n,n) 
int n;	 length of full row in a 
int nqx;  partial length of row to decode 
int nqy;  partial length of column (<=n) 
int nbitplanes;	 number of bitplanes to decode */
{
    int log2n, k, bit, b, nqmax;
    int nx,ny,nfx,nfy,c;
    int nqx2, nqy2;
    unsigned char *scratch;

    /* log2n is log2 of max(nqx,nqy) rounded up to next power of 2 */
    nqmax = (nqx>nqy) ? nqx : nqy;
    log2n = (int) (log((float) nqmax)/log(2.0)+0.5);
    if (nqmax > (1<<log2n)) log2n += 1;

    /* allocate scratch array for working space */
    nqx2=(nqx+1)/2;
    nqy2=(nqy+1)/2;
    scratch = (unsigned char *) malloc(nqx2*nqy2);
    if (scratch == (unsigned char *) NULL) 
    {
        fprintf(stderr, "qtree_decode: insufficient memory\n");
        exit(-1);
    }
    /*
     * now decode each bit plane, starting at the top
     * A is assumed to be initialized to zero
     */
    for (bit = nbitplanes-1; bit >= 0; bit--) 
    {
    	/* Was bitplane was quadtree-coded or written directly? */
    	b = input_nybble(infile);
    	if (b == 0) 
        {
    	    /* bit map was written directly */
    	    read_bdirect(infile,a,n,nqx,nqy,scratch,bit);
    	}
        else if (b != 0xf) 
             {
    	        fprintf(stderr, "qtree_decode: bad format code %x\n",b);
    	        exit(-1);
    	     } 
             else 
             {
    	        /* bitmap was quadtree-coded, do log2n expansions 
                   read first code */
    		scratch[0] = input_huffman(infile);
    		/* do log2n expansions, reading codes from file as necessary */
    		nx = 1;
    		ny = 1;
    		nfx = nqx;
    		nfy = nqy;
    		c = 1<<log2n;
    		for (k = 1; k<log2n; k++) 
    		{
    	           /* this somewhat cryptic code generates the sequence
		    * n[k-1] = (n[k]+1)/2 where n[log2n]=nqx or nqy */
    		   c = c>>1;
    		   nx = nx<<1;
    		   ny = ny<<1;
    		   if (nfx <= c) { nx -= 1; } else { nfx -= c; }
    		   if (nfy <= c) { ny -= 1; } else { nfy -= c; }
    		   qtree_expand(infile,scratch,nx,ny,scratch);
		}
		/* copy last set of 4-bit codes to bitplane bit of array a */
                qtree_bitins(scratch,nqx,nqy,a,n,bit);
	     }
      }
      free(scratch);
}

/********************************************************************/

/*
 * do one quadtree expansion step on array a[(nqx+1)/2,(nqy+1)/2]
 * results put into b[nqx,nqy] (which may be the same as a)
 */
static void qtree_expand(FILE *infile, unsigned char a[],
                         int nx, int ny, unsigned char b[])
{
    int i;

    /* first copy a to b, expanding each 4-bit value */
    qtree_copy(a,nx,ny,b,ny);

    /* now read new 4-bit values into b for each non-zero element */
    for (i = nx*ny-1; i >= 0; i--) 
    {
        if (b[i] != 0) b[i] = input_huffman(infile);
    }
}

/********************************************************************/

/*
 * copy 4-bit values from a[(nx+1)/2,(ny+1)/2] to b[nx,ny], expanding
 * each value to 2x2 pixels
 * a,b may be same array
 */
static void qtree_copy(unsigned char a[], int nx, int ny, unsigned char b[],
                       int n)
{
    int i, j, k, nx2, ny2;
    int s00, s10;

    /* first copy 4-bit values to b start at end in case a,b are same array */
    nx2 = (nx+1)/2;
    ny2 = (ny+1)/2;
    k = ny2*(nx2-1)+ny2-1; /* k   is index of a[i,j] */
    for (i = nx2-1; i >= 0; i--) 
    {
        s00 = 2*(n*i+ny2-1); /* s00 is index of b[2*i,2*j] */
        for (j = ny2-1; j >= 0; j--) 
        {
            b[s00] = a[k];
            k -= 1;
            s00 -= 2;
	}
    }

    /* now expand each 2x2 block */
    for (i = 0; i<nx-1; i += 2) 
    {
        s00 = n*i; /* s00 is index of b[i,j] */
        s10 = s00+n; /* s10 is index of b[i+1,j]	*/
        for (j = 0; j<ny-1; j += 2) 
        {
            b[s10+1] =  b[s00]     & 1;
            b[s10  ] = (b[s00]>>1) & 1;
            b[s00+1] = (b[s00]>>2) & 1;
            b[s00  ] = (b[s00]>>3) & 1;
            s00 += 2;
            s10 += 2;
        }
        if (j < ny) 
        {
	  /* row size is odd, do last element in row s00+1, s10+1 are off edge*/
          b[s10  ] = (b[s00]>>1) & 1;
          b[s00  ] = (b[s00]>>3) & 1;
	}
    }

    if (i < nx) 
    {
	/* column size is odd, do last s10, s10+1 are off edge */
        s00 = n*i;
        for (j = 0; j<ny-1; j += 2) 
        {
            b[s00+1] = (b[s00]>>2) & 1;
            b[s00  ] = (b[s00]>>3) & 1;
            s00 += 2;
	}
        if (j < ny) 
        {
	   /* both row and column size are odd, do corner element
              s00+1, s10, s10+1 are off edge */
            b[s00  ] = (b[s00]>>3) & 1;
	}
    }
}

/********************************************************************/

/*
 * Copy 4-bit values from a[(nx+1)/2,(ny+1)/2] to b[nx,ny], expanding
 * each value to 2x2 pixels and inserting into bitplane BIT of B.
 * A,B may NOT be same array (it wouldn't make sense to be inserting
 * bits into the same array anyway.)
 */

static void qtree_bitins(unsigned char a[], int nx, int ny, 
                         int b[], int n, int bit)
{
    int i, j, k;
    int s00, s10;

    /* expand each 2x2 block */
    k = 0;	/* k   is index of a[i/2,j/2]	*/
    for (i = 0; i<nx-1; i += 2) 
    {
        s00 = n*i; /* s00 is index of b[i,j] */
        s10 = s00+n; /* s10 is index of b[i+1,j] */
        for (j = 0; j<ny-1; j += 2) 
        {
            b[s10+1] |= ( a[k]     & 1) << bit;
            b[s10  ] |= ((a[k]>>1) & 1) << bit;
            b[s00+1] |= ((a[k]>>2) & 1) << bit;
            b[s00  ] |= ((a[k]>>3) & 1) << bit;
            s00 += 2;
            s10 += 2;
            k += 1;
	}
        if (j < ny)
        {
            /* row size is odd, do last element in row
	     * s00+1, s10+1 are off edge */
            b[s10  ] |= ((a[k]>>1) & 1) << bit;
            b[s00  ] |= ((a[k]>>3) & 1) << bit;
            k += 1;
	}
    }

    if (i < nx) 
    {
	 /* column size is odd, do last row s10, s10+1 are off edge */
	 s00 = n*i;
	 for (j = 0; j<ny-1; j += 2) 
         {
	     b[s00+1] |= ((a[k]>>2) & 1) << bit;
	     b[s00  ] |= ((a[k]>>3) & 1) << bit;
	     s00 += 2;
	     k += 1;
         }
	 if (j < ny) 
         {
	    /* both row and column size are odd, do corner element
	     * s00+1, s10, s10+1 are off edge */
            b[s00  ] |= ((a[k]>>3) & 1) << bit;
            k += 1;
         }
    }
}

/********************************************************************/

static void read_bdirect(FILE *infile, int a[], int n, int nqx,
                         int nqy, unsigned char scratch[], int bit)
{
    int i;

    /* read bit image packed 4 pixels/nybble */
    for (i = 0; i < ((nqx+1)/2) * ((nqy+1)/2); i++) 
    {
        scratch[i] = input_nybble(infile);
    }

    /* insert in bitplane BIT of image A */
    qtree_bitins(scratch,nqx,nqy,a,n,bit);
}

/********************************************************************/

/*
 * Huffman decoding for fixed codes
 *
 * Coded values range from 0-15
 *
 * Huffman code values (hex):
 *
 *	3e, 00, 01, 08, 02, 09, 1a, 1b,
 *	03, 1c, 0a, 1d, 0b, 1e, 3f, 0c
 *
 * and number of bits in each code:
 *
 *	6,  3,  3,  4,  3,  4,  5,  5,
 *	3,  5,  4,  5,  4,  5,  6,  4
 */

static int input_huffman(FILE *infile)
{
    int c;

    /* get first 3 bits to start */
    c = input_nbits(infile,3);
    if (c < 4) 
    {
        /* this is all we need return 1,2,4,8 for c=0,1,2,3 */
        return(1<<c);
    }

    /* get the next bi */
    c = input_bit(infile) | (c<<1);
    if (c < 13) 
    {
	/* OK, 4 bits is enough */
        switch (c) 
        {
            case  8 : return(3);
            case  9 : return(5);
            case 10 : return(10);
            case 11 : return(12);
            case 12 : return(15);
	}
    }

    /* get yet another bit */
    c = input_bit(infile) | (c<<1);
    if (c < 31)
    {
	/* OK, 5 bits is enough */
        switch (c) 
        {
            case 26 : return(6);
            case 27 : return(7);
            case 28 : return(9);
            case 29 : return(11);
            case 30 : return(13);
	}
    }

    /* need the 6th bit */
    c = input_bit(infile) | (c<<1);
    if (c == 62)
    {
        return(0);
    } 
    else
    {
        return(14);
    }
}

/********************************************************************/

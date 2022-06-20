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
**    File:  IM_CompTool.cc
**
*******************************************************************************
**
**    DECRIPTION    ROUTINES for image compression
**    ---------- 
**
*******************************************************************************
**
** int readint(FILE *infile)
** 
** Read integer A one byte at a time from infile.
**
*******************************************************************************
**
** int myread(FILE *file, char buffer[], int n)
**
** read n bytes from file into buffer
** returns number of bytes read (=n) if successful, <=0 if not
**     
** this version is for VMS C: each read may return
** fewer than n bytes, so multiple reads may be needed
** to fill the buffer.
**     
**  I think this is unnecessary for Sun Unix, but it won't hurt
**  either, so I'll leave it in.   
**
*******************************************************************************
**
** void qread(FILE *infile, char *a, int n)
**
** read n bytes from file into buffer a
** 
*******************************************************************************
**
** void writeint(FILE *outfile, int a)
**
** write integer a to file
**
*******************************************************************************
**
** int mywrite(FILE *file, char buffer[], int n)
**
**  write n bytes from  buffer to file
**  returns number of bytes read (=n) if successful, <=0 if not
**
*******************************************************************************
**
** void qwrite(FILE *outfile, char *a, int n)
**
** write n bytes from  buffer to file
**
*****************************************************************************/

// static char sccsid[] = "@(#)IM_CompTool.cc 3.1 96/05/02 CEA 1995 @(#)";

#include"IM_Obj.h"
#include"IM_CompTool.h"
#include"IM_IOTools.h"

/* the following line comes from CFITSIO */
// static float testfloat = TESTFLOAT;  /* use to test floating pt format */

/* THE BIT BUFFER */
static int buffer_in; /* Bits waiting to be input	*/
static int bits_to_go_in;	/* Number of bits still in buffer_in */

int bitcount;
static int buffer_out; /* Bits buffer_outed for output	*/
static int bits_to_go_out;	/* Number of bits free in buffer_out */


/********************************************************************/

/* INITIALIZE BIT INPUT */

void start_inputing_bits()
{
    /* Buffer starts out with no bits in it */
    bits_to_go_in = 0;
}

/********************************************************************/

/* INITIALIZE FOR BIT OUTPUT */

void start_outputing_bits()
{
	buffer_out = 0;	/* Buffer is empty to start	*/
	bits_to_go_out = 8;	/* with	*/
	bitcount = 0;
}

/********************************************************************/

/* INPUT A BIT */

int input_bit(FILE *infile)
{
    if (bits_to_go_in == 0) 
    {	
        /* Read the next byte if no */
        buffer_in = getc(infile);	/* bits are left in buffer_in */
        if (buffer_in == EOF) 
        {
	    /* end of file is an error for this application */
            fprintf(stderr, "input_bit: unexpected end-of-file\n");
            exit(-1);
	}
        bits_to_go_in = 8;
    }
    /* Return the next bit */
    bits_to_go_in -= 1;
    return((buffer_in>>bits_to_go_in) & 1);
}

/********************************************************************/

/* OUTPUT A BIT */

void output_bit(FILE *outfile, int bit)
{
	buffer_out <<= 1;	/* Put bit at end of buffer_out */
	if (bit) buffer_out |= 1;
	bits_to_go_out -= 1;
	if (bits_to_go_out == 0) 
        {                               /* Output buffer_out if it is */
	 putc(buffer_out & 0xff,outfile);   /* now full */
	 bits_to_go_out = 8;
	 buffer_out = 0;
	}
	bitcount += 1;
}

/********************************************************************/
/* INPUT N BITS (N must be <= 8) */

int input_nbits(FILE *infile, int n)
{
    int c;

    if (bits_to_go_in < n)
    {
	/* need another byte's worth of bits */
        buffer_in <<= 8;
        c = getc(infile);
        if (c == EOF)
        {
	    /* end of file is an error for this application */
            fprintf(stderr, "input_nbits: unexpected end-of-file\n");
            exit(-1);
	}
        buffer_in |= c;
        bits_to_go_in += 8;
    }
    /*  now pick off the first n bits */
    bits_to_go_in -= n;
    return( (buffer_in>>bits_to_go_in) & ((1<<n)-1) );
}

/********************************************************************/

/* OUTPUT N BITS (N must be <= 8) */

void output_nbits(FILE *outfile, int bits, int n)
{
	/*
	 * insert bits at end of buffer_out
	 */
	buffer_out <<= n;
	buffer_out |= ( bits & ((1<<n)-1) );
	bits_to_go_out -= n;
	if (bits_to_go_out <= 0) {
		/*
		 * buffer_out full, put out top 8 bits
		 */
		putc((buffer_out>>(-bits_to_go_out)) & 0xff,outfile);
		bits_to_go_out += 8;
	}
	bitcount += n;
}

/********************************************************************/

/* FLUSH OUT THE LAST BITS */

void done_outputing_bits(FILE *outfile)
{
	if (bits_to_go_out < 8) {
		putc(buffer_out<<bits_to_go_out,outfile);
		/* count the garbage bits too */
		bitcount += bits_to_go_out;
	}
}

/********************************************************************/

int readint(FILE *infile)
{
    long int a;
    char *bufdata;

//	cout << "readint" << endl ;

    bufdata = new char [4]; 
    QFREAD(bufdata, 4, infile, "readint");

/* the 3 following lines come from CFITSIO */
#if BYTESWAPPED == 1
//	cout << "readint: BYTESWAPPED == True " << endl ;
    ffswap4((INT32BIT *) bufdata, 1); /* reverse order of bytes in each value */
#endif
    a = ((long int *)bufdata)[0];

    delete [] bufdata;
    return(a);
}

/********************************************************************/

float readfloat(FILE *infile)
{
    float a;
    char *bufdata;

    bufdata = new char [4]; 
    QFREAD(bufdata, 4, infile, "readfloat");

/* the n following lines come from CFITSIO */
#if BYTESWAPPED == 1
//	cout << "readfloat: BYTESWAPPED == True " << endl ;
    ffswap4((INT32BIT *) bufdata, 1); /* reverse order of bytes in each value */
#endif
    a = ((float *)bufdata)[0];

    delete [] bufdata;
    return(a);
}

/********************************************************************/

void writefloat(FILE *outfile, float a)
{
    char *bufdata;
    float *buff;

    bufdata = new char [4]; 
    buff = (float *) bufdata; 
    buff[0] = a;
#if BYTESWAPPED == 1
//	cout << "writefloat: BYTESWAPPED == True " << endl ;
     ffswap4((INT32BIT *) bufdata, 1); /* reverse order of bytes in each value */
#endif
   QFWRITE(bufdata,4, outfile, "writefloat");
   delete [] bufdata;
}

/********************************************************************/

void qread(FILE *infile, char *a, int n)
{
    if(myread(infile, a, n) != n) 
    {
    	fprintf(stderr, "Error when reading ...\n");
    	exit(-1);
    }
}

/********************************************************************/

int myread(FILE *file, char buffer[], int n)
{
    /*
     * read n bytes from file into buffer
     * returns number of bytes read (=n) if successful, <=0 if not
     *
     * this version is for VMS C: each read may return
     * fewer than n bytes, so multiple reads may be needed
     * to fill the buffer.
     *
     * I think this is unnecessary for Sun Unix, but it won't hurt
     * either, so I'll leave it in.
     */
    int nread, total;

    /* keep reading until we've read n bytes */
    total = 0;
    while ( (nread = fread(&buffer[total], 1, n-total, file)) > 0) 
    {
	total += nread;
	if (total==n) return(total);
    }

    /* end-of-file or error occurred if we got to here */
    return(nread);
}

/********************************************************************/

void writeint(FILE *outfile, int a)
{
    char *bufdata;
    long int *buff;

//	cout << "writeint" << endl ;

    bufdata = new char [4]; 
    buff = (long int *) bufdata; 
    buff[0] = a;
/* the 3 following lines come from CFITSIO */
#if BYTESWAPPED == 1
//	cout << "writeint: BYTESWAPPED == True " << endl ;
     ffswap4((INT32BIT *) bufdata, 1); /* reverse order of bytes in each value */
#endif
   QFWRITE(bufdata,4, outfile, "writeint");
   delete [] bufdata;
}

/********************************************************************/

void qwrite(FILE *outfile, char *a, int n)
{
	if(mywrite(outfile, a, n) != n) 
        {
		fprintf(stderr, "Error when writing ...");
		exit(-1);
	}
}

/********************************************************************/

int mywrite(FILE *file, char buffer[], int n)
{
    /*
     * read n bytes from file into buffer
     * returns number of bytes read (=n) if successful, <=0 if not
     */
    return ( fwrite(buffer, 1, n, file) );
}

/********************************************************************/


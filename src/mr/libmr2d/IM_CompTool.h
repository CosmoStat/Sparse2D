
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
**    File:  IM_CompTool.h
**
*******************************************************************************
**
**    DECRIPTION    Definition for image compression
**    ---------- 
**
*****************************************************************************/


void start_inputing_bits();
void start_outputing_bits();
int input_bit(FILE *infile);
void output_bit(FILE *outfile, int bit);
void output_nbits(FILE *outfile, int bits, int n);
void done_outputing_bits(FILE *outfile);

int  readint(FILE *infile);
void qread(FILE *infile, char *a, int n);
int  myread(FILE *file, char buffer[], int n);
float readfloat(FILE *infile);

void writeint(FILE *outfile, int a);
void qwrite(FILE *outfile, char *a, int n);
int  mywrite(FILE *file, char buffer[], int n);
void writefloat(FILE *outfile, float a);


int input_nbits(FILE *infile, int n);

void encode(FILE *outfile, int *a, int nx, int ny, float scale);
void dodecode(FILE *infile, int a[], int nx, int ny, unsigned char nbitplanes);
void decode(FILE *infile, int  **a, int  *nx, int  *ny, float  *scale);

void qtree_decode(FILE *infile, int a[], int n, int nqx, 
                   int nqy, int nbitplanes);
void qtree_encode(FILE *outfile, int a[], int n, int nqx, 
                  int nqy, int nbitplanes);


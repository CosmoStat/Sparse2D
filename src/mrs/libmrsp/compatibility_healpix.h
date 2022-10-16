#ifndef HEALPIX_COMPATIBILITY_H
#define HEALPIX_COMPATIBILITY_H

#include <string>
#include <vector>
#include "fitsio.h"
#include "arr.h"
#include "datatypes.h"

/*! \defgroup fitsgroup FITS-related functionality */
/*! \{ */

template<typename T> struct FITSUTIL {};

template<> struct FITSUTIL<signed char>
  { enum { DTYPE=TBYTE }; };
template<> struct FITSUTIL<short>
  { enum { DTYPE=TSHORT }; };
template<> struct FITSUTIL<int>
  { enum { DTYPE=TINT }; };
template<> struct FITSUTIL<long>
  { enum { DTYPE=TLONG }; };
template<> struct FITSUTIL<long long>
  { enum { DTYPE=TLONGLONG }; };
template<> struct FITSUTIL<float>
  { enum { DTYPE=TFLOAT }; };
template<> struct FITSUTIL<double>
  { enum { DTYPE=TDOUBLE }; };

/*! Converts a FITS type code (i.e. the data type of a FITS column) to the
    corresponding Planck type code. */
inline int ftc2type (int ftc)
  {
  switch (ftc)
    {
    case TLOGICAL : return PLANCK_BOOL;
    case TBYTE    : return PLANCK_INT8;
    case TSHORT   : return PLANCK_INT16;
    case TINT32BIT: return PLANCK_INT32;
    case TLONGLONG: return PLANCK_INT64;
    case TFLOAT   : return PLANCK_FLOAT32;
    case TDOUBLE  : return PLANCK_FLOAT64;
    case TSTRING  : return PLANCK_STRING;
        default: cerr << "ftc2type: unsupported component type" << endl;
            exit(-1);

    }
  }

/*! Converts a Planck type code to the corresponding FITS type code
    (i.e. the data type of a FITS column). */
inline int type2ftc (int type)
  {
  switch (type)
    {
    case PLANCK_BOOL   : return TLOGICAL;
    case PLANCK_INT8   : return TBYTE;
    case PLANCK_INT16  : return TSHORT;
    case PLANCK_INT32  : return TINT32BIT;
    case PLANCK_INT64  : return TLONGLONG;
    case PLANCK_FLOAT32: return TFLOAT;
    case PLANCK_FLOAT64: return TDOUBLE;
    case PLANCK_STRING : return TSTRING;
    default: cerr << "type2ftc: unsupported component type" << endl;
            exit(-1);
    }
  }

#endif

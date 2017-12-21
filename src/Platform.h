//============================================================================
//
//  This file is part of GFC, the GNSS FOUNDATION CLASS.
//
//  The GFC is free software; you can redistribute it and/or modify
//  it under the terms of the GNU Lesser General Public License as published
//  by the Free Software Foundation; either version 3.0 of the License, or
//  any later version.
//
//  The GFC is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with GFC; if not, write to the Free Software Foundation,
//  Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110, USA
//
//  Copyright 2015, lizhen
//
//============================================================================

/* This is where all platform specific includes, defines and crud should go.
   Unless, of course, it has to go somewhere else. :-)
*/

#ifndef GFC_PLATFORM_H
#define GFC_PLATFORM_H

typedef int               int_t;
typedef unsigned int        uint_t;
typedef long              long_t;
typedef float             float_t;
typedef double            double_t;
typedef long double       ldouble_t;


typedef signed char GINT8 ;


// used for microsoft visual studio

#ifdef _MSC_VER

#include <cstdlib>

#define HAVE_STRING_H 1
#define STDC_HEADERS  1

//// To get rid of 'stdint.h' for Microsoft visual studio
//#if (_MSC_VER < 1300 )
//    typedef signed char       int8_t;
//    typedef signed short      int16_t;
//    typedef signed int        int32_t;
//    typedef unsigned char     uint8_t;
//    typedef unsigned short    uint16_t;
//    typedef unsigned int      uint32_t;
//    typedef signed __int64    int64_t;
//    typedef unsigned __int64  uint64_t;
//#elif(_MSC_VER <= 1500)
//    typedef signed __int8     int8_t;
//    typedef signed __int16    int16_t;
//    typedef signed __int32    int32_t;
//    typedef unsigned __int8   uint8_t;
//    typedef unsigned __int16  uint16_t;
//    typedef unsigned __int32  uint32_t;
//    typedef signed __int64    int64_t;
//    typedef unsigned __int64  uint64_t;
//	
//	typedef int               int_t;
//	typedef unsigned int        uint_t;
//	typedef long              long_t;
//	typedef float             float_t;
//	typedef double            double_t;
//	typedef long double       ldouble_t;
//	
//	
//#else        
//// the other compiler
//#include <stdint.h>
//
//typedef int               int_t;
//typedef unsigned int        uint_t;
//typedef long              long_t;
//typedef float             float_t;
//typedef double            double_t;
//typedef long double       ldouble_t;
//
//
//#endif

//#include <sys/types.h>
#include <sys/timeb.h>

#elif defined __SUNPRO_CC

#include <sys/types.h>
#include <sys/timeb.h>

#else   

#include <stdint.h>

#endif  // _MSC_VER

#endif  // GFC_PLATFORM_H

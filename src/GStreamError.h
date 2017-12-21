#pragma ident "$Id: FFStreamError.hpp 3140 2012-06-18 15:03:02Z susancummins $"



/**
 * @file FFStreamError.hpp
 * Exceptions for FFStream
 */

#ifndef GFC_GSTREAMERROR_H
#define GFC_GSTREAMERROR_H

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







#include "GException.h"

namespace gfc
{
	  //用宏NEW_EXCEPTION_CLASS来定义一个新的类
      /// FFStreamError is an exception for when the file read doesn't
      /// match the specs for that file type.
      /// @ingroup exceptionclass
      /// @ingroup formattedfile
   NEW_GEXCEPTION_CLASS(GStreamError, gfc::GException);
}

#endif 

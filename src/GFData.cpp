
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



/**
 * @file GFData.cpp
 * Formatted File Data base class
 */

#include "GFData.h"
#include "GStream.h"

namespace gfc
{
   void GFData::putRecord(GStream& s) const
      throw(GStreamError, gfc::GStringException)
   { 
      s.tryGStreamPut(*this);
   }
   
   void GFData::getRecord(GStream& s)
      throw(GStreamError, gfc::GStringException)
   { 
      s.tryGStreamGet(*this);
   }
   
    // override and  friend function
   std::ostream& operator<<(std::ostream& o, const GFData& f)
         throw(GStreamError, gfc::GStringException)
   {
      GStream* ffs = dynamic_cast<GStream*>(&o);
      if (ffs)
      {
         f.putRecord(*ffs);
         return o;
      }
      else
      {
         GStreamError e("operator<< stream argument must be an GFStream");
         GFC_THROW(e);
      }
   }

   std::istream& operator>>(std::istream& i, GFData& f)
         throw(GStreamError, gfc::GStringException)
   {
      GStream* ffs = dynamic_cast<GStream*>(&i);
      if (ffs)
      {
         f.getRecord(*ffs);
         return i;
      }
      else
      {
         GStreamError e("operator<< stream argument must be an GFStream");
         GFC_THROW(e);
      }
   }
}

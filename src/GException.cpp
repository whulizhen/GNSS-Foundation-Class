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
 * @file Exception.cpp
 * Exceptions for all of GPSTK, including location information
 */
 
#include <sstream>
#include "GException.h"

using std::ostream;
using std::ostringstream;
using std::streambuf;
using std::string;
using std::endl;

namespace gfc
{
   	
   void GExceptionLocation::dump(ostream& s) const
      throw()
   { 
      s << getFileName() << ":" 
#ifdef __FUNCTION__
        << getFunctionName() << ":" 
#endif
        << getLineNumber(); 
   }

   GException::GException()
      throw()
   {
   }

   GException::GException(const string& errorText,
                        const unsigned long& errId,
                        const Severity& sever)
      throw()
   {
      text.push_back(errorText);
      errorId = errId;
      severity = sever;
   }

   GException::GException(const GException& e)
      throw()
         : errorId(e.errorId),
           locations(e.locations),
           severity(e.severity),
           text(e.text),
           streamBuffer(e.streamBuffer)
   {}

   GException& GException::operator=(const GException& e)
      throw()
   {
      errorId = e.errorId;
      locations = e.locations;
      severity = e.severity;
      text = e.text;
         // reuse existing stream objects, no matter.
         //streambuf(), ostream((streambuf*)this),
      streamBuffer = e.streamBuffer;

      return *this;
   }

   GException& GException::addLocation(
      const GExceptionLocation& location)
      throw()
   {
      locations.push_back(location);
      return *this;
   }

   const GExceptionLocation GException::getLocation(
      const size_t& index) const
      throw()
   {
      if ( (int_t)index < 0 || index>=getLocationCount())
      {
         return GExceptionLocation();
      }
      else
      {
         return locations[index];
      }
   }

   size_t GException::getLocationCount() const
      throw()
   {
      return locations.size();
   }

   GException& GException::addText(const string& errorText)
      throw()
   {
      text.push_back(errorText);
      return *this;
   }

   string GException::getText(const size_t& index) const
      throw()
   {
      if ( (int_t)index < 0 || index>=getTextCount())
      {
         string tmp;
         return tmp;
      }
      else
      {
         return text[index];
      }
   }

   size_t GException::getTextCount() const
      throw()
   {
      return text.size();
   }

   void GException::dump(ostream& s) const
      throw()
   {
      int i;
      for ( i=0; i<getTextCount(); i++ )
      {
         s << "text " << i << ":" << this->getText(i) << endl;
      }
      for ( i=0; i<getLocationCount(); i++)
      {
         s << "location " << i << ":" << getLocation(i).what() << endl;
      }
   }

   int GException::overflow(int c)
   {
      if (c == '\n' || !c)
      {
         if (streamBuffer.length() == 0)
         {
            return c;
         }
         addText(streamBuffer);
         streamBuffer = "";
         return c;
      }
      streamBuffer.append(1, (char)c);
      return c;
   }

   string GExceptionLocation::what() const
      throw()
   {
      ostringstream oss;
      this->dump(oss);
      return oss.str();
   }

   const char* GException::what() const
      throw()
   {
      ostringstream oss;
      this->dump(oss);
      return oss.str().c_str();
   }

    ostream& operator<<( ostream& s, 
                         const GException& e )
       throw()
    { 
       e.dump(s); 
       return s;
    }

    ostream& operator<<( ostream& s,
                         const GExceptionLocation& e )
       throw()
    {
       e.dump(s);
       return s;
    }

} // namespace gpstk


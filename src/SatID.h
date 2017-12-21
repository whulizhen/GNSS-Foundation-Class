//============================================================================
//
//  This file is part of GPSTk, the GPS Toolkit.
//
//  The GPSTk is free software; you can redistribute it and/or modify
//  it under the terms of the GNU Lesser General Public License as published
//  by the Free Software Foundation; either version 3.0 of the License, or
//  any later version.
//
//  The GPSTk is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//  GNU Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with GPSTk; if not, write to the Free Software Foundation,
//  Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110, USA
//  
//  Copyright 2004, The University of Texas at Austin
//
//============================================================================

//============================================================================
//
//This software developed by Applied Research Laboratories at the University of
//Texas at Austin, under contract to an agency or agencies within the U.S. 
//Department of Defense. The U.S. Government retains all rights to use,
//duplicate, distribute, disclose, or release this software. 
//
//Pursuant to DoD Directive 523024 
//
// DISTRIBUTION STATEMENT A: This software has been approved for public 
//                           release, distribution is unlimited.
//
//=============================================================================

#ifndef GFC_SATID_H
#define GFC_SATID_H

#include <iostream>
#include <iomanip>
#include <sstream>

//#include "gps_constants.hpp"

/**
 * @file SatID.hpp
 * gpstk::SatID - navigation system-independent representation of a satellite.
 */

namespace gfc
{
   // forward declarations
   class SatID;
	
   /// Satellite identifier consisting of a satellite number (PRN, etc.)
   /// and a satellite system
   class SatID
   {
   public:
      /// Supported satellite systems
      enum SatelliteSystem
      {
         sysGPS ,
         sysGalileo,
         sysGlonass,
         sysGeosync,
         sysLEO,
         sysTransit,
         sysBDS,
         sysQZSS,
         sysMixed,
		 sysINS,
		 sysGBAS,
		 sysSBAS,
         sysUserDefined,
         sysUnknown
      };

      /// empty constructor, creates an invalid object
      SatID() { id=-1; system=sysUnknown; }
	  	
      /// explicit constructor, no defaults
      /// @note if s is given a default value here,
      /// some compilers will silently cast int to SatID.
      SatID( int p, const SatelliteSystem& s) { id=p; system=s; }
	  	
      // operator=, copy constructor and destructor built by compiler
	  	
      /// Convenience method used by dump().
      static std::string convertSatelliteSystemToString(const SatelliteSystem& s)
      {
         switch(s)
         {
            case sysGPS:			return "GPS";           break;
            case sysGalileo:		return "Galileo";       break;
            case sysGlonass:		return "GLONASS";       break;
            case sysGeosync:		return "Geostationary"; break;
            case sysLEO:			return "LEO";           break;
            case sysTransit:		return "Transit";       break;
            case sysBDS:            return "BeiDou";        break;
            case sysQZSS:           return "QZSS";          break;
            case sysMixed:			return "Mixed";         break;
            case sysUserDefined:	return "UserDefined";   break;
            case sysUnknown:		return "Unknown";       break;
            default:                return "??";            break;
         };
      }

         /// Convenience output method.
      void dump( std::ostream& s ) const
      {
         s << convertSatelliteSystemToString(system) << " " << id;
      }

      /// operator == for SatID
      bool operator==(const SatID& right) const
      { return ((system == right.system) && (id == right.id)); }

      /// operator != for SatID
      bool operator!=(const SatID& right) const
      { return !(operator==(right)); }

      /// operator < for SatID : order by system, then number
      bool operator<(const SatID& right) const
      {
         if (system==right.system)
            return (id<right.id);
         return (system<right.system);
      }

      /// operator > for SatID
      bool operator>(const SatID& right) const
      {  return (!operator<(right) && !operator==(right)); }

      /// operator <= for SatID
      bool operator<=(const SatID& right) const
      { return (operator<(right) || operator==(right)); }

      /// operator >= for SatID
      bool operator>=(const SatID& right) const
      { return !(operator<(right)); }

      /// return true if this is a valid SatID
      /// @note assumes all id's are positive and less than 100;
      ///     plus GPS id's are less than or equal to MAX_PRN (32).
      /// @note this is not used internally in the gpstk library
      bool isValid() const
      {
         switch(system)
         {
            case sysGPS:      return ( id > 0 && id <= 32  );   //MAX_PRN_GPS
			case sysBDS:      return ( id > 0 && id <= 35  );   //MAX_PRN_BDS
            case sysGalileo:  return ( id > 0 && id <= 30  );   //MAX_PRN_GAL
            case sysGlonass:  return ( id > 0 && id <= 24  );   //MAX_PRN_GLS
            case sysGeosync:  return ( id > 0 && id <= 100 );   //MAX_PRN_GEO
            case sysLEO:      return ( id > 0 && id <= 100 );   //MAX_PRN_LEO
            case sysTransit:  return ( id > 0 && id <= 100 );   //MAX_PRN_TST
            default:          return ( id > 0 && id <= 100 );
         }
      }

      int id;                   ///< satellite identifier, e.g. PRN
      SatelliteSystem system;   ///< system for this satellite

   }; // class SatID

   namespace StringUtils
   {
      inline std::string asString(const SatID& p)
      {
         std::ostringstream oss;
         p.dump(oss);
         return oss.str();
      }
   }

      /// stream output for SatID
   inline std::ostream& operator<<(std::ostream& s, const SatID& p)
   {
      p.dump(s);
      return s;
   }

} // namespace gpstk

#endif

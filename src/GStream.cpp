#pragma ident "$Id: FFStream.cpp 3140 2012-06-18 15:03:02Z susancummins $"

/**
 * @file FFStream.cpp
 * Formatted File Stream base class
 */

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




#include "GStream.h"


namespace gfc
{
	   	
      /*
       * Overrides fstream:open so derived classes can make appropriate
       * internal changes (line count, header info, etc).
       */
   void GStream::open( const char* fn, std::ios::openmode mode )
   {
		
#ifdef _MSC_VER
      fstream::open(fn, mode);
#else
      std::fstream::open(fn, mode);
#endif
      filename = GString(fn);
      recordNumber = 0;
      
      std::ios::clear();
       
   }  // End of method 'GStream::open()'
    
    
    
      // A function to help debug GStreams
   void GStream::dumpState(std::ostream& s) const
   {

      s << "filename:" << filename
        << ", recordNumber:" << recordNumber;
      s << ", exceptions:";
		
       if (exceptions() & std::ios::badbit)  s << "bad ";
      if (exceptions() & std::ios::failbit) s << "fail ";
      if (exceptions() & std::ios::eofbit)  s << "eof ";
      if (exceptions() == 0) s << "none";
       
      s << ", rdstate:";
       
      if (rdstate() & std::ios::badbit)  s << "bad ";
      if (rdstate() & std::ios::failbit) s << "fail ";
      if (rdstate() & std::ios::eofbit)  s << "eof ";
      if (rdstate() == 0)  s << "none";
      s << std::endl;

   }  // End of method 'FFStream::dumpState()'



      // the crazy double try block is so that no gfc::GException throws
      // get masked, allowing all exception information (line numbers, text,
      // etc) to be retained.
   void GStream::tryGStreamGet(GFData& rec)
      throw(GStreamError, gfc::GStringException)
   {

         // Mark where we start in case there is an error.
      long initialPosition = tellg();
      unsigned long initialRecordNumber = recordNumber;
      clear();
       
      try
      {
         try
         {
            rec.reallyGetRecord(*this);
            recordNumber++;
         }
         catch (EndOfFile& e)
         {
            // EOF - do nothing - eof causes fail() to be set which
            // is handled by std::fstream
             e.addText("In record " + GString(recordNumber) );
             e.addText("In file " + filename);
             e.addLocation(FILE_LOCATION);
            mostRecentException = e;
         }
          catch (gfc::GException &e)
         {
             mostRecentException = GStreamError("gfc::GException thrown: " +
                                                GString(e.what()));
             mostRecentException.addText("In record " + GString(recordNumber));
            
            mostRecentException.addText("In file " + filename);
            mostRecentException.addLocation(FILE_LOCATION);
            clear();
            seekg(initialPosition);
            recordNumber = initialRecordNumber;
            setstate(std::ios::failbit);
            conditionalThrow();
         }
         catch (gfc::GStringException& e)
         {
            e.addText("In record " + GString(recordNumber));
            e.addText("In file " + filename);
            e.addLocation(FILE_LOCATION);
            mostRecentException = e;
            clear();
            seekg(initialPosition);
            recordNumber = initialRecordNumber;
            setstate(std::ios::failbit);
            conditionalThrow();
         }
            // catches some errors we can encounter
         catch (GStreamError& e)
         {
            e.addText("In record " + GString(recordNumber));
            e.addText("In file " + filename);
            e.addLocation(FILE_LOCATION);
            mostRecentException = e;
            clear();
            seekg(initialPosition);
            recordNumber = initialRecordNumber;
            setstate(std::ios::failbit);
            conditionalThrow();
         }
      }
         // this is if you throw an FFStream error in the above catch
         // block because the catch(...) below will mask it otherwise.
         // This also takes care of catching StringExceptions
      catch (gfc::GException &e)
      {
         GFC_RETHROW(e);
      }
      catch (std::ifstream::failure &e)
      {
            // setting failbit when catching FFStreamError can cause
            // this exception to be thrown. in this case, we don't want
            // to lose the exception info so only make a new exception
            // if this isn't a fail() case
         if (!fail())
         {
            mostRecentException = GStreamError("ifstream::failure thrown: " + GString(e.what()));
            mostRecentException.addText("In file " + filename);
            mostRecentException.addLocation(FILE_LOCATION);
         }
         conditionalThrow();
      }
    
       catch ( std::exception &e)
      {
         mostRecentException = GStreamError("gfc::GException thrown: " + GString(e.what()));
         mostRecentException.addText("In file " + filename);
         mostRecentException.addLocation(FILE_LOCATION);
         setstate(std::ios::failbit);
         conditionalThrow();
      }
      catch (...)
      {
         mostRecentException = GStreamError("Unknown Gexception thrown");
         mostRecentException.addText("In file " + filename);
         mostRecentException.addLocation(FILE_LOCATION);
         setstate(std::ios::failbit);
         conditionalThrow();
      }

   }  // End of method 'GStream::tryGStreamGet()'



      // the crazy double try block is so that no gpstk::Exception throws 
      // get masked, allowing all exception information (line numbers, text,
      // etc) to be retained.
   void GStream::tryGStreamPut(const GFData& rec)
      throw(GStreamError, gfc::GStringException)
   {
         // Mark where we start in case there is an error.
      long initialPosition = tellg();
      unsigned long initialRecordNumber = recordNumber;
      clear();
       
      try
      {
         try
         {
            rec.reallyPutRecord(*this);
            recordNumber++;
         }
         catch (gfc::GException &e)
         {
               // if this is a stream failure, don't mask it and let the
               // later catch block handle it
            if (dynamic_cast<std::ifstream::failure*>(&e))
               throw;
			
               // the catch(FFStreamError) below will add file information
               // to this exception
            mostRecentException = GStreamError("gfc::GException thrown: " + GString(e.what()));
            mostRecentException.addLocation(FILE_LOCATION);
            setstate(std::ios::failbit);
            conditionalThrow();
         }
         catch (gfc::GStringException& e)
         {
            e.addText("In record " + GString(recordNumber) );
            e.addText("In file " + filename);
            e.addLocation(FILE_LOCATION);
            mostRecentException = e;
            seekg(initialPosition);
            recordNumber = initialRecordNumber;
            setstate(std::ios::failbit);
            conditionalThrow();
         } 
            // catches some errors we can encounter
         catch (GStreamError& e)
         {
            e.addText("In record " + GString(recordNumber) );
            e.addText("In file " + filename);
            e.addLocation(FILE_LOCATION);
            mostRecentException = e;
            seekg(initialPosition);
            recordNumber = initialRecordNumber;
            setstate(std::ios::failbit);
            conditionalThrow();
         }
      }
         // this is if you throw an FFStream error in the above catch
         // block because the catch(...) below will mask it otherwise.
         // This also takes care of catching StringExceptions
      catch (gfc::GException &e)
      {
         GFC_RETHROW(e);
      }
      catch (std::ifstream::failure &e)
      {
            // setting failbit when catching FFStreamError can cause
            // this exception to be thrown. in this case, we don't want
            // to lose the exception info so only make a new exception
            // if this isn't a fail() case
         if (!fail())
         {
            mostRecentException = GStreamError("ifstream::failure thrown: " +
                                                GString(e.what()));
            mostRecentException.addText("In file " + filename);
            mostRecentException.addLocation(FILE_LOCATION);
         }
         conditionalThrow();
      }
    
       catch ( std::exception &e)
      {
         mostRecentException = GStreamError("std::exception thrown: " +
                                             GString(e.what()));
         mostRecentException.addText("In file " + filename);
         mostRecentException.addLocation(FILE_LOCATION);
         setstate(std::ios::failbit);
         conditionalThrow();
      }
      catch (...)
      {
         mostRecentException = GStreamError("Unknown exception thrown");
         mostRecentException.addText("In file " + filename);
         mostRecentException.addLocation(FILE_LOCATION);
         setstate(std::ios::failbit);
         conditionalThrow();
      }
       
   }  // End of method 'GStream::tryGStreamPut()'

}  // End of namespace gfc



/**
 * @file FFTextStream.hpp
 * An FFStream for text files
 */

#ifndef GFC_GTEXTSTREAM_H
#define GFC_GTEXTSTREAM_H


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
		
      /** @addtogroup formattedfile */
      //@{
		
      /**
       * An FFStream is meant for reading text.  This also includes an
       * internal line count and a read line method. When reading and
       * using the formattedGetLine() call, the lineNumber automatically
       * increments.  However, any other read and all write calls do not
       * update the line number - the derived class or programmer
       * needs to make sure that the reader or writer increments
       * lineNumber in these cases.
       */
   class GTextStream : public GStream
   {
   public:
         /// Destructor
      virtual ~GTextStream() {};
       
         /// Default constructor
      GTextStream()
            : lineNumber(0) {};
       
        /** Common constructor.
          *
          * @param fn file name.
          * @param mode file open mode (std::ios)
          */
      GTextStream( const char* fn,
                    std::ios::openmode mode=std::ios::in )
         : GStream(fn, mode), lineNumber(0)
      {};
       
       
         /** Common constructor.
          *
          * @param fn file name.
          * @param mode file open mode (std::ios)
          */
      GTextStream( const std::string& fn,
                    std::ios::openmode mode=std::ios::in )
         : GStream( fn.c_str(), mode ), lineNumber(0)
      {};
       
       
         /// Overrides open to reset the line number.
      virtual void open( const char* fn,
                         std::ios::openmode mode )
      { GStream::open(fn, mode); lineNumber = 0; };
       
       
         /// Overrides open to reset the line number.
      virtual void open( const std::string& fn,
                         std::ios::openmode mode )
      { open(fn.c_str(), mode); };
       
       
       
         /**
          * Like std::istream::getline but checks for EOF and removes '/r'.
          * Also increments lineNumber.  When \a expectEOF is true and EOF
          * is found, an gpstk::EndOfFile exception is thrown.  If
          * \a expectEOF is false and an EOF is encountered, an
          * gpstk::FFStreamError is thrown.
          * @param line is set to the value of the line read from the file.
          * @param expectEOF set true if finding EOF on this read is acceptable.
          * @throw EndOfFile if \a expectEOF is true and an EOF is encountered.
          * @throw FFStreamError if EOF is found and \a expectEOF is false
          * @throw gpstk::StringUtils::StringException when a string error occurs
          * or if any other error happens.
          * @warning There is a maximum line length of 1500 characters when
          * using this function.
          */
      inline void formattedGetLine( GString& line,
                                    const bool expectEOF = false )
         throw(EndOfFile, GStreamError, gfc::StringException);
			
			
   protected:


         /// calls FFStream::tryFFStreamGet and adds line number information
      virtual void tryGStreamGet(GFData& rec)
         throw(GStreamError, gfc::StringException)
      {

         unsigned int initialLineNumber = lineNumber;
		 	
         try
         {
            GStream::tryGStreamGet(rec);
         }
         catch(gfc::GException& e)
         {
             e.addText( GString("Near file line ") + GString(lineNumber) );
            lineNumber = initialLineNumber;
            mostRecentException = e;
            conditionalThrow();
         }

       };


         /// calls FFStream::tryFFStreamPut and adds line number information
      virtual void tryGStreamPut(const GFData& rec)
         throw(GStreamError, gfc::StringException)
      {
          
         unsigned int initialLineNumber = lineNumber;
          
         try
         {
            GStream::tryGStreamPut(rec);
         }
         catch(gfc::GException& e)
         {
             e.addText( std::string("Near file line ") + GString(lineNumber) );
            lineNumber = initialLineNumber;
            mostRecentException = e;
            conditionalThrow();
         }

      }
       
       
   private:
       /// The internal line count. When writing, make sure
       /// to increment this.
       unsigned int lineNumber;
       
       
   }; // End of class 'GTextStream'


    
      // inline function must get realized in the header file
      // the reason for checking ffs.eof() in the try AND catch block is
      // because if the user enabled exceptions on the stream with exceptions()
      // then eof could throw an exception, in which case we need to catch it
      // and rethrow an EOF or FFStream exception.  In any event, EndOfFile
      // gets thrown whenever there's an EOF and expectEOF is true
   void GTextStream::formattedGetLine( GString& line,
                                        const bool expectEOF )
         throw(EndOfFile, GStreamError, gfc::StringException)
   {
		
      try
      {
            // The following constant used to be 256, but with the change to
            // RINEX3 formats the possible length of a line increased
            // considerably. A RINEX3 observation file line for Galileo may
            // be 1277 characters long (taking into account all the possible
            // types of observations available, plus the end of line
            // characters), so this constant was conservatively set to
            // 1500 characters. Dagoberto Salazar.
         const int MAX_LINE_LENGTH = 1500;
		 char templine[MAX_LINE_LENGTH + 1]={0};
         getline(templine, MAX_LINE_LENGTH);
         lineNumber++;
            //check if line was longer than 256 characters, if so error
         if( fail() && !eof())
         {
            GStreamError err("Line too long");
            GFC_THROW(err);
         }
         line = templine;
          //gfc::StringUtils::stripTrailing(line, '\r');
          line.stripTrailing();
            // catch EOF when stream exceptions are disabled
         if ((gcount() == 0) && eof())
         {
            if (expectEOF)
            {
               EndOfFile err("EOF encountered");
               GFC_THROW(err);
            }
            else
            {
               GStreamError err("Unexpected EOF encountered");
               GFC_THROW(err);
            }
         }
      }
      catch(gfc::GException &e)
      {

            // catch EOF when exceptions are enabled
         if ( (gcount() == 0) && eof())
         {
            if (expectEOF)
            {
               EndOfFile err("EOF encountered");
               GFC_THROW(err);
            }
            else
            {
               GStreamError err("Unexpected EOF");
               GFC_THROW(err);
            }
         }
         else
         {
            GStreamError err("Critical file error: " +
                              GString(e.what()));
            GFC_THROW(err);
         }  // End of 'if ( (gcount() == 0) && eof())'

      }  // End of 'try-catch' block

   }  // End of method 'GTextStream::formattedGetLine()'

      //@}

}  // End of namespace gfc
#endif   // GFC_GTEXTSTREAM_H

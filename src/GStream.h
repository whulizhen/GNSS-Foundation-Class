#pragma ident "$Id: FFStream.h 3140 2012-06-18 15:03:02Z susancummins $"

/**
* @file FFStream.hpp
* Formatted File Stream, root class to provide formatted I/O
* operators ('<<' & '>>')
*/

#ifndef GFC_GSTREAM_H
#define GFC_GSTREAM_H

//============================================================================
//
//  This file is part of GFC, the GNSS FOUNDATION CLASS.
//
//  The GFC is free software; you can redistribute it and/or modify
//  it under the terms of the GNU Lesser General Public License as published
//  by the Free Software Foundation; either version 2.1 of the License, or
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



#include <iostream>
#include <fstream>
#include <string>
#include <typeinfo>

#include "GStreamError.h"
#include "GFData.h"

#include "GString.h"

#ifdef _MSC_VER
using namespace std;
#endif

namespace gfc
{
	/** @addtogroup formattedfile */
	//@{
    
	/// This gets thrown if a valid EOF occurs on formattedGetLine.
	/// @ingroup exceptionclass
	NEW_GEXCEPTION_CLASS(EndOfFile, gfc::GStreamError);
    
	/**
	* Formatted File Stream (GStream).
	* This is just a root class to provide the single point formatted i/o
	* operators (such as '<<' & '>>' ).
	*
	* As a special design consideration,
	* all exceptions thrown are based on gfc::GException - all
	* std::exception throws are rethrown as gfc::GException.
	* Furthermore, exceptions will not be thrown unless exceptions
	* are set to be thrown:
	* @code
	* fs.exceptions(std::fstream::failbit);
	* @endcode
	* where \c fs is the name of your file stream.
	* Then when an exception occurs, conditionalThrow() will throw the
	* last thrown exception.
	* Otherwise when an exception occurs, the stream sets
	* \c ios::fail and will not read any more.  Exceptions for this
	* class store the record number of the file for when the exception
	* occurred as well as the file name and any detailed information
	* about the error.  For gfc::GTextStream, the line number
	* of the file where the error was found is also recorded, allowing
	* for easy location of file problems.
	*
	* When operating on the file, recordNumber will automatically increment
	* with each read and write operation. When a file is opened with the
	* constructor or with open(), all internal GStream variables are
	* reset. Derived classes should make sure any of their internal
	* variables are reset when either of those function are called.
	*
	* Many file types have header data as part of the file format. When
	* reading the file, the reader is not required to explicitly read in
	* the header to access the data.  To facilitate this, each of these
	* stream classes has an internal header object that will store the
	* header. The stream keeps track of whether it read the
	* header or not, and reads the header if the internal state says
	* it hasn't been read.  When writing a file, the stream's
	* internal header is used for those formats which use header information
	* to determine what data is in the records.
	* See RinexObsHeader::reallyGetRecord() and
	* RinexObsData::reallyGetRecord()
	* for an example of this.
	*
	* \sa FFData for more information
	* \sa RinexObsData::reallyGetRecord() and
	*     RinexObsHeader::reallyGetRecord() for more information for files
	*     that read header data.
	*
	* @warning When using open(), the internal header data of the stream
	* is not guaranteed to be retained.
	*/
	class GStream : public std::fstream
	{
        
        /// FFData is a friend so it can access the try* functions.
        friend class GFData;
        
	public:
        
        /**
         * Default constructor
         */
        GStream()
        : recordNumber(0) {};
        
		/// Virtual destructor
		virtual ~GStream(void) {};
		
		/** Common constructor.
		*
		* @param fn file name.
		* @param mode file open mode (std::ios)
		*/
		GStream( const char* fn,
			std::ios::openmode mode = std::ios::in )
			:
#ifdef _MSC_VER
		fstream(fn, mode),
#else
		std::fstream(fn, mode),
#endif
			recordNumber(0), filename(fn)
		{
            // call the clear() function of father class fstream
            std::fstream::clear();
        }
        
        
		/** Common constructor.
		*
		* @param fn file name.
		* @param mode file open mode (std::ios)
		*/
		GStream( const GString& fn,
			std::ios::openmode mode=std::ios::in )
			:
#ifdef _MSC_VER
		fstream(fn.c_str(), mode),
#else
		std::fstream(fn.c_str(), mode),
#endif
			recordNumber(0), filename(fn)
		{
            std::fstream::clear();
        };
        
        
		/**
		* Overrides fstream:open so derived classes can make appropriate
		* internal changes (line count, header info, etc).
		*/
		virtual void open( const char* fn, std::ios::openmode mode );
        
        
		/**
		* Overrides fstream:open so derived classes can make appropriate
		* internal changes (line count, header info, etc).
		*/
		virtual void open( const GString& fn,std::ios::openmode mode )
		{ open( fn.c_str(), mode ); };
        
        
		/// A function to help debug GStreams
		void dumpState(std::ostream& s = std::cout) const;
        
        
		/**
		* Throws \a mostRecentException only if the stream is enabled
		* to throw exceptions when failbit is set.
		* You can set this behavior with the following line of code:
		* @code
		* ffstreamobject.exceptions(ifstream::failbit);
		* @endcode
		* where \a gstreamobject is the name of your stream object.
		*/
		inline void conditionalThrow(void) throw( GStreamError)
		{
            if ( std::ios::exceptions()  & std::fstream::failbit)
			{
				GFC_THROW(mostRecentException);
			}
		};
        
        
	protected:
		
		/// Encapsulates shared try/catch blocks for all file types
		/// to hide std::exception.
		virtual void tryGStreamGet( GFData& rec )
			throw(GStreamError, gfc::GStringException);
        
        
		/// Encapsulates shared try/catch blocks for all file types
		/// to hide std::exception.
		virtual void tryGStreamPut(const GFData& rec)
			throw(GStreamError, gfc::GStringException);
      
    private:
        
        ///@name Data members
        ///@{
        /// This stores the most recently thrown exception.
        GStreamError mostRecentException;
        
        /// keeps track of the number of records read
        unsigned int recordNumber;
        
        /// file name
        GString filename;
        
        //@}

	}; // End of class 'GStream'

	//@}

}  // End of namespace gfc

#endif   // GFC_GSTREAM_H


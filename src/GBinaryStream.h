
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
 * @file GBinaryStream.h
 * An GStream for binary file reading
 */

#ifndef GFC_GBINARYSTREAM_H
#define GFC_GBINARYSTREAM_H

#include "GStream.h"

namespace gfc
{
    /** @defgroup formattedfile Formatted file I/O */
    //@{
    
    /**
     * This is an FFStream that is required to be binary.  It also includes
     * functions for reading and writing binary file.  Otherwise, this
     * is the same as FFStream.
     */
    class GBinaryStream : public GStream
    {
    public:
        /// destructor
        virtual ~GBinaryStream() {};
        
        /// Default constructor
        GBinaryStream() {}
        
        /**
         * Constructor - opens the stream in binary mode if not set.
         * @param fn file name.
         * @param mode file open mode (std::ios)
         */
        GBinaryStream(const char* fn,
                       std::ios::openmode mode=std::ios::in|std::ios::binary)
        : GStream(fn, mode|std::ios::binary) {}
        
        /// Overrides open to ensure binary mode opens
        virtual void open(const char* fn, std::ios::openmode mode)
        { GStream::open(fn, mode|std::ios::binary); }
        
        /**
         * Reads a T-object directly from the stream
         * in binary form.
         * @throw FFStreamError when the size of the data read
         * from this stream doesn't match the size of a T-object.
         * @return a T-object
         */
        template <class T> T getData() throw(GStreamError, EndOfFile)
        {
            T data;
            getData((char*)&data, sizeof(T));
            return data;
        } // end of getData(GStream& strm)
        
        void getData(char* buff, size_t length) throw(GStreamError, EndOfFile)
        {
            try
            {
                read(buff, length);
            }
            catch(std::exception& exc)
            {
                if (gcount() != (std::streamsize)length && eof())
                {
                    EndOfFile err("EOF encountered");
                    GFC_THROW(err);
                }
                else
                {
                    GStreamError err(exc.what());
                    std::cout << err << std::endl;
                    GFC_THROW(err);
                }
            }
            catch(...)
                {
                    GStreamError err("Unknown exception");
                    GFC_THROW(err);
                }
        } // end of getData(char*, size_t))
                    
     /**
        * Writes a T-object directly from the stream
        * in binary form.
        * @param data the data to be written.
        * @throw FFStreamError when the size of the data written
        * to this stream doesn't match the size of a T-object.
        * @return a T-object
        */
        template <class T> void writeData(const T& data)
            throw(FFStreamError)
            {
                        //T temp = data;
                writeData((char*)&data, sizeof(T));
                return;
            } // end of writeData(FFStream& strm, const T& data)
                    
            void writeData(const char* buff, size_t length)
                throw(GStreamError)
                {
                    try
                    {
                        write(buff, length);
                    }
                    catch(std::exception& exc)
                    {
                        GStreamError err(exc.what());
                        GFC_THROW(err);
                    }
                    catch(...)
                    {
                        GStreamError err("Unknown exception");
                        GFC_THROW(err);
                    }
                        
                    if (fail() || bad())
                    {
                        GStreamError err("Error writing data");
                        GFC_THROW(err);
                    }
                   return;
                } // end of writeData(const char*, size_t)
                    
        }; // end of class GBinarayStream
        //@}
} // end of namespace gfc
#endif

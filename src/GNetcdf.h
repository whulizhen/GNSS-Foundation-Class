//
//  GNetcdf.hpp
//  GFC
//
//  Created by lizhen on 16/1/29.
//  Copyright © 2016年 lizhen. All rights reserved.
//

#ifndef GNetcdf_hpp
#define GNetcdf_hpp


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


#include <stdio.h>
#include "netcdf.h"
#include "GString.h"
#include "GTime.h"

namespace gfc
{
    
// for the process of netcdf data file, nc file
//
/*
 references:
   http://www.unidata.ucar.edu/software/netcdf/docs/group__datasets.html#ga555f117e74e9c4eed0d4423fedd27bfb
 
   http://www.unidata.ucar.edu/software/netcdf/docs/modules.html
 
*/
    
    
class GNetcdf
{
    
#define MAXCHAR_NAME 1000
#define MAXCHAR_VAL  1000
    
    struct ncAttribute
    {
        GString m_attName;
        int m_ID;
        int m_dataType; // the index in static variable dataType;
        GString m_value; // this value should be converted into different data types according to m_dataType
        ncAttribute()
        {
            m_ID = -1;
            m_dataType = -1;
        }
    };
    
    
    //  the definition of Dimension in nc data;
    struct ncDimension
    {
        GString m_dimName;
        size_t  m_dimVal;
        int m_ID;
        bool m_isUnlimited;
        double m_maxmum; // if the m_isUnlimited is true, then m_maxmum = -1;
        double m_minimum;
        double m_step;  // the step for this certain dimension
        std::vector<ncAttribute> m_attri; // usually, the dimension can appear in the Variables List, so we can get the attribute of the dimensions
        std::vector<double> m_data;  // store the data in this diminsion
        // when we read Variables information, we can know these attributes
        ncDimension()
        {
            m_dimVal = 0;
            m_ID = -1 ;
            m_isUnlimited = false;
            m_minimum = 0.0;
            m_maxmum = 0.0;
            m_step = 0.0;
        }
    };
    
    
    // this data structure just describle the Variable , it DOES NOT store the real data!!
    struct ncVariable
    {
        GString m_varName;
        std::vector<ncAttribute> m_attribute; // all the attributes belond to this variable
        int m_varID;  // the varID for searching
        int m_dataType; // most are in float, CANNOT process the datatype in string!!!!
        std::vector<int> m_dimen;  // the index of dimen in the dimension list, the last dimension varying fastest
                                // m_dimen[0] is the first dimension, m_dimen[1] is the second dimension ,......
    };
    
    
    
public:
    
    
    //NC_WRITE, NC_SHARE, or NC_WRITE|NC_SHARE
    enum NCMODE{ MODE_NOWRITE = NC_NOWRITE,MODE_WRITE = NC_WRITE, MODE_SHARE = NC_SHARE  };
    static GString DataType[13];
    //enum FileMode{};
    
    // construction
    GNetcdf()
    {
        m_fileID = 0;
        m_ndim = 0;
        m_nvar = 0;
        m_nattGLOBAL = 0;
    }
    
    GNetcdf(GString filename, int mode);
    
    //
    ~GNetcdf()
    {
        if( m_fileID != 0 )
        {
            nc_close(m_fileID);
            m_fileID = 0;
        }
    }
    
    void loadNCFile(GString filename, int mode);
    
    //get all the data of variable
    void getData( double* data, GString varName);
    
    //get one special data
    double getData( size_t* indexArray, GString varName);
    
    // data should be a one-dimension array, with the last dimension varying fastest
    void   getData( GString varName, double* data);
    
    void   getData( GString varName, size_t* start, size_t* count, double* data);
    
    
    
    std::vector<double> getDimensionData( gfc::GString dimName );
    
    //return value: datatype of the attribute
    // value: the actual data of the attribute
    GString getAttriValue( GString attriName,void* value );
    
    //double getValid_min(GString dimName);
    //double getValid_max(GString dimName);
    int varName2ID(GString varName, int* index);
    
    int getNdim() {return m_ndim;}
    int getNvar() {return m_nvar;}
    
private:
    
    void handle_error(int status);
    
    GString m_ncfileName;
    int m_fileID;   // netcdf file identifier for ecah file
    int m_ndim;     // the number of dimensions in this file
    int m_nvar;     // the number of variables in this file
    int m_nattGLOBAL;     // the number of global attitutes in this file(just for this file)
    GString  m_ncver;      // the version of the netcdf API
    std::vector<ncAttribute>   m_attriGLOBAL; // information on the global attributes
    std::vector<ncDimension> m_dimensions; // dimensions have attributes as well
    std::vector<ncVariable>    m_variables; // variables have their own attributes
    
};


}


#endif /* GNetcdf_hpp */

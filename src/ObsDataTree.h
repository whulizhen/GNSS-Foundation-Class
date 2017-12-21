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
//  Copyright 2015, lizhen, hpulizhen@163.com
//
//============================================================================

/**
 * @file DataTree.h
 * gfc::DataTree - Identifies types of values
 */

#ifndef GFC_OBSDATATREE_H
#define GFC_OBSDATATREE_H

#include "GSensor.h"
#include "GSourceID.h"
#include "GOBSData.h"
#include "GFData.h"

namespace gfc
{
    // class DataTree shoulb be a template class for the reuse of all the other data in different measuerement scenes, such as total stations , GNSS, INS etc.
    // here T1 means SourceID object and  T2 means SensorID object
    // the tree is just a abstract concept in logic. What we really need to pay attention to is the data storage structure!
    // actually, the tree just describle the relationship between the TREENODEs, the biggest problem is that every Node has different data structure
    //template< class CLASS_SOURCE, class CLASS_SENSOR, class CLASS_OBSTYPE >
    class ObsDataTree : public GFData
    {
        
    public:
        
        ObsDataTree( void );
        virtual ~ObsDataTree( void );
        
        //采用层次遍历构建多叉树，具体实现在子类中实现；这里作为类的接口
        virtual void buildObsTree();
        
        
    private:
        
        // the actually stored data
        //CLASS_SOURCE  m_source;  //as the root
        //std::list< CLASS_SENSOR > m_sensors; // as the second layer
        //std::list< CLASS_OBSTYPE >  m_data;     // as the third layer
        
    };
}

#endif
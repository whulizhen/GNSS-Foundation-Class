//
//  Source.h
//  GFC
//
//  Created by lizhen on 15/10/15.
//  Copyright © 2015年 lizhen. All rights reserved.
//

#ifndef GSource_h
#define GSource_h

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

#include "GString.h"

namespace gfc
{
    
    // class SourceID is the father class of all the classes related to stations
    // also, it is a kind of NODE in the observables TREE
    class GSourceID
    {
        
      public:
        
        GSourceID() {}
        
        GSourceID(GString sourcename)
        {
            
            m_sourceName = sourcename;
            
            m_sourceName.stripLeading_v();
            
            m_sourceName.stripTrailing_v();
            
        }
        
        void setSourceName(GString sourcename)
        {
            m_sourceName = sourcename;
            m_sourceName.stripLeading_v();
            
            m_sourceName.stripTrailing_v();
        }
        
        //GSourceID(char*  sourcename)  {m_sourceName = GString(sourcename);}
        
        virtual ~GSourceID() {}
        
        GString getSourceName() const
        { return m_sourceName; }
        
        
        GSourceID& operator= (const GSourceID& right)  //赋值重载
        {
            this->m_sourceName = right.m_sourceName;
            return *this;
        }
        
        GSourceID( const GSourceID& right )   //copy construction function
        {
            this->m_sourceName = right.m_sourceName;
        }
        
      private:
        
        GString  m_sourceName;
    };
    
    //these functions can not be member functions!!!
    bool   operator== (const GSourceID& left,const GSourceID& right) ;
    bool   operator!= (const GSourceID& left,const GSourceID& right) ;
    bool   operator>  (const GSourceID& left,const GSourceID& right) ;
    bool   operator<  (const GSourceID& left,const GSourceID& right) ;
    bool   operator>= (const GSourceID& left,const GSourceID& right) ;
    bool   operator<= (const GSourceID& left,const GSourceID& right) ;
    
    
}  // end of namespace gfc



#endif /* Source_h */

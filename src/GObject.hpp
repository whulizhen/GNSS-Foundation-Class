//
//  GObject.hpp
//  GFC
//  GObject is the father class of almost the GFC classes
//  Created by lizhen on 15/10/29.
//  Copyright © 2015年 lizhen. All rights reserved.
//

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

#ifndef GObject_hpp
#define GObject_hpp

#include <stdio.h>
#include "GString.h"
namespace gfc
{
    
    class GObject
    {
        
    public:
        
        GObject() { m_className = "GObject";}
        virtual ~GObject() {}
        virtual GString getClassName() const
        {
            return m_className;
        }
        
    private:
        GString m_className;
        
    };
    
    
    
}


#endif /* GObject_hpp */

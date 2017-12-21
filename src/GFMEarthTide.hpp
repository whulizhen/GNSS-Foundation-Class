//
//  GFMEarthTide.hpp
//  GFC
//
//  Created by lizhen on 20/07/2016.
//  Copyright © 2016 lizhen. All rights reserved.
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


#ifndef GFMEarthTide_hpp
#define GFMEarthTide_hpp

#include <stdio.h>

#include "GForceModel.hpp"

namespace gfc
{
    // solid earth tide correction
    class GFMEarthTide  : public GForceModel
    {
        
    public:
        
        GFMEarthTide()
        {
            setForceName("GFMEarthTide"); // should be N body
        }
        
        void doCompute(GVector& satpos_eci, GVector& sunpos_eci, GVector& moonpos_eci);
        
    private:
        
        
    };
    
    
} // end of namespace gfc



#endif /* GFMEarthTide_hpp */

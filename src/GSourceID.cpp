//
//  GSourceID.cpp
//  GFC
//
//  Created by lizhen on 24/08/2016.
//  Copyright © 2016 lizhen. All rights reserved.
//

#include <stdio.h>
#include "GSourceID.h"
namespace gfc
{
    bool  operator>=  (const GSourceID& left, const GSourceID& right)
    {
        return left.getSourceName() >= right.getSourceName();
    }
    
    bool  operator<=  (const GSourceID& left,const GSourceID& right)
    {
       return left.getSourceName() <= right.getSourceName();
    }
    
    bool  operator >  (const GSourceID& left,const GSourceID& right)
    {
        return left.getSourceName() > right.getSourceName();
    }
    
    bool  operator <  (const GSourceID& left,const GSourceID& right)
    {
        return left.getSourceName() < right.getSourceName();
    }
    
    bool   operator== (const GSourceID& left,const GSourceID& right)
    {
        return left.getSourceName() == right.getSourceName();
    }
    
    /*不等号重载*/
    bool operator!=(const GSourceID& left,const  gfc::GSourceID &right)
    {
        return !( left==right );
    }
    

    
    
} // end of namespace gfc
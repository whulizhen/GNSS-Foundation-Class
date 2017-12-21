//
//  GPOD.hpp
//  GFC
//
//  Created by lizhen on 01/04/2017.
//  Copyright © 2017 lizhen. All rights reserved.
//

#ifndef GPOD_hpp
#define GPOD_hpp

#include <stdio.h>

#include <vector>

#include "GNSSData.h"

#include "GOBSDataStorage.h"

#include "GRinex.hpp"

#include <fstream>

namespace gfc
{


class GPreciseOrbitDetermination
{
    
public:
    
    void loadObsdata(GString obsfilename);
    
    
private:
    
    // the raw obs data, multiple eppoches
    std::vector<GgnssStorage> rawdata;
    GRinex*  allobsfile;
    
    
};

}  // end of namespace gfc



#endif /* GPOD_hpp */

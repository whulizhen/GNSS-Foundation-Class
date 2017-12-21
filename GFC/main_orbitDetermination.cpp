//
//  main_orbitDetermination.cpp
//  GFC
//
//  Created by lizhen on 06/08/2017.
//  Copyright © 2017 lizhen. All rights reserved.
//

#include <stdio.h>
#include "GPreciseOrbitDetermination.hpp"
using namespace gfc;

int main()
{
    GPreciseOrbitDetermination pod;
    pod.loadObsdata("/Users/lizhen/myProject/GFC/GFC_Proj/Build/Products/Debug/obs/test.txt");
    
    return 0;
}

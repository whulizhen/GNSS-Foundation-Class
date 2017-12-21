//
//  GSLRValidation.hpp
//  GFC
//
//  Created by lizhen on 05/09/2016.
//  Copyright © 2016 lizhen. All rights reserved.
//

#ifndef GSLRValidation_hpp
#define GSLRValidation_hpp

#include <stdio.h>



#include "GSLRcrd.hpp"

#include "GIERS.hpp"

#include "GSpacecraft.hpp"


namespace gfc
{
    
    // this class is used to validate the orbits with SLR measurements
    class GSLRValidation
    {
        
        
        struct slrResInfo
        {
            double res;
            double obs;
            double geodis1;
            double geodis2;
            GString stationName;
            double ele;
            double beita; // degree
            double eps;   // the EPS angle
            double latitude;
            double temperature;
            double pressure;
            double earthTide;
            double staOffset;
            double satOffset;
            double trop;
            double rel_error;  //relativity
            
            double los[3];
            double dlra[3];
            double dxyz[3];
            double dtide[3];
        };
        
        
    public:
        
        void setOutputFile(GString filename)
        {
            m_outputFileName = filename;
        }
        
        void loadDataFile( GString filename );
        
        void setSpacecraft( GSpaceCraft* spaceVehicle );
        
        void computeResidual();
        void computeResidual2();
        
    private:
        
        GString  m_outputFileName;
        
        GslrStorage m_slrdata;
        
        GSLRcrd    m_datafile;
        
        GSpaceCraft* m_spacecraft; // storing the precise orbit of the spacecraft
        
        GslrStorage::iterator m_it;  // the iterator to the current spacecraft
        
    };

    
} // end of namespace gfc

#endif /* GSLRValidation_hpp */

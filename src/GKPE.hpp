//
//  GKPE.hpp
//  GFC
//
//  Created by lizhen on 13/07/2016.
//  Copyright © 2016 lizhen. All rights reserved.
//

#ifndef GKPE_hpp
#define GKPE_hpp

#include <stdio.h>

#include "GVector.hpp"
#include "GFCCONST.h"
#include <math.h>


namespace gfc
{
    /* the class for Keplerian Elements  */
    class GKeplerianElements
    {
        
    public:
        long double m_sma;  ///!< semi-major axis (km)
        long double m_ecc;  ///!< eccentricity of ellipse (ratio-->no units)
        long double m_inc;  ///!< inclination of orbital plane (rad)
        long double m_argp; ///!< argument of perigee (rad)
        long double m_raan; ///!< right ascension of ascending node (rad)
        long double m_tran; ///!< true anomaly (rad)
       // long double ecan; ///!< eccentric anomaly (rad)
       // long double mean; ///!< mean anomaly (rad)
        
        
        GKeplerianElements();
        GKeplerianElements(long double sma_in,long double ecc_in, long double inc_in, long double argp_in,
                                               long double raan_in,long double tran_in );
       
        
        /*convert from position and velocity to keplerian elements*/
        void PV2KP( GVector p, GVector v );
        void KP2PV(GVector& p, GVector& v);
        
        /*convert from keplerian elements to position and velocity */
        static void kpe2pv(GKeplerianElements& kpe,  GVector& p, GVector& v );
        static void propagate(GKeplerianElements& kpe,  double t, GVector& p, GVector& v);
        
    private:
        
       static long double get_ecan_from_tran(long double ecc, long double tran);
       static long double get_tran_from_ecan(long double ecan, long double ecc);
       static long double get_mean_from_ecan(long double ecan,long double ecc);
       static long double get_ecan_from_mean(long double mean, long double ecc);
        
    };

    
    
    
}  // end of namespace



#endif /* GKPE_hpp */

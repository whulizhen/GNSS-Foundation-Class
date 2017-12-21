//
//  GKPE.cpp
//  GFC
//
//  Created by lizhen on 13/07/2016.
//  Copyright © 2016 lizhen. All rights reserved.
//

#include "GKPE.hpp"

namespace gfc
{
    GKeplerianElements::GKeplerianElements()
    {
        m_sma = 0.0;  ///!< semi-major axis (km)
        m_ecc = 0.0;  ///!< eccentricity of ellipse (ratio-->no units)
        m_inc = 0.0;  ///!< inclination of orbital plane (rad)
        m_argp= 0.0; ///!< argument of perigee (rad)
        m_raan= 0.0; ///!< right ascension of ascending node (rad)
        m_tran= 0.0; ///!< true anomaly (rad)
    }
    
    GKeplerianElements::GKeplerianElements(long double sma_in,long double ecc_in, long double inc_in, long double argp_in,
                                           long double raan_in,long double tran_in )
    {
        m_sma = sma_in;  ///!< semi-major axis (km)
        m_ecc = ecc_in;  ///!< eccentricity of ellipse (ratio-->no units)
        m_inc = inc_in;  ///!< inclination of orbital plane (rad)
        m_argp= argp_in; ///!< argument of perigee (rad)
        m_raan= raan_in; ///!< right ascension of ascending node (rad)
        m_tran= tran_in; ///!< true anomaly (rad)
    }
    
    // predict the position of the satellite after time t
    void GKeplerianElements::propagate(GKeplerianElements& kpe,  double t, GVector& p, GVector& v)
    {
        GVector pos;
        long double GM = 398600.4415; //GCONST("GM");
        
        //the mean motion
        double n = sqrt(GM/kpe.m_sma)/ kpe.m_sma;
        
        //eccentric anomaly
        
        double ecan_t0 = get_ecan_from_tran(kpe.m_ecc, kpe.m_tran);
        
        // the mean anomaly at t0
        
        double M0 = get_mean_from_ecan(ecan_t0,kpe.m_ecc);
        
        double Mt = M0 + n*t;
        
        Mt = fmod(Mt, 2.0*M_PI);
        
        //because eccentricity is fixed
        long double ecan_t  = get_ecan_from_mean(Mt, kpe.m_ecc);
        
        double tran_t = get_tran_from_ecan(ecan_t, kpe.m_ecc);
        
        GKeplerianElements kpe_t = kpe;
        kpe_t.m_tran = tran_t;
        
        GKeplerianElements::kpe2pv(kpe_t, p, v);
        
    }
    
    
    void GKeplerianElements::kpe2pv(GKeplerianElements& kpe,  GVector& p, GVector& v )
    {
        long double GM = GCONST("GM");
        
        long double sqrt1me2 = sqrt(1.0 - kpe.m_ecc * kpe.m_ecc);
        
        long double ecan = get_ecan_from_tran(kpe.m_ecc, kpe.m_tran);
        
        long double cos_ecan = cos(ecan);
        long double sin_ecan = sin(ecan);
        
        // Compute the magnitude of the Gaussian vectors at the required point
        long double gaussX = kpe.m_sma * (cos_ecan - kpe.m_ecc);    // Magnitude of
        long double gaussY = kpe.m_sma * sqrt1me2 * sin_ecan; // Gaussian vectors
        
        long double XYdotcommon = sqrt(GM / kpe.m_sma) / (1.0 - kpe.m_ecc * cos_ecan);
        
        long double gaussXdot = -sin_ecan * XYdotcommon;           // Gaussian vel.
        long double gaussYdot = cos_ecan * sqrt1me2 * XYdotcommon; // components
        
        long double cos_inc = cos(kpe.m_inc);
        long double sin_inc = sin(kpe.m_inc);
        
        long double cos_argp = cos(kpe.m_argp);
        long double cos_raan = cos(kpe.m_raan);
        
        long double sin_argp = sin(kpe.m_argp);
        long double sin_raan = sin(kpe.m_raan);
        
        long double cc = cos_argp * cos_raan;
        long double cs = cos_argp * sin_raan;
        long double sc = sin_argp * cos_raan;
        long double ss = sin_argp * sin_raan;
        
        long double P[3]={0.0}, Q[3]={0.0}; // Components of the unit Gaussian vectors
        
        P[0] = cc - ss * cos_inc;
        P[1] = cs + sc * cos_inc;
        P[2] = sin_argp * sin_inc;
        
        Q[0] = -sc - cs * cos_inc;
        Q[1] = -ss + cc * cos_inc;
        Q[2] = cos_argp * sin_inc;
        
        p.x = gaussX * P[0] + gaussY * Q[0];
        p.y = gaussX * P[1] + gaussY * Q[1];
        p.z = gaussX * P[2] + gaussY * Q[2];
        
        v.x = gaussXdot * P[0] + gaussYdot * Q[0];
        v.y = gaussXdot * P[1] + gaussYdot * Q[1];
        v.z = gaussXdot * P[2] + gaussYdot * Q[2];
    }
    
    
    void GKeplerianElements::KP2PV(gfc::GVector &p, gfc::GVector &v)
    {
        GKeplerianElements::kpe2pv(*this, p, v);
    }
    
    /*convert from position and velocity to keplerian elements*/
    void GKeplerianElements::PV2KP( GVector p, GVector v )
    {
        double GM = GCONST("GM");
        double twopi = GCONST("2PI");
        
        long double  n[2];
        GVector e,h;
        // e = eccentricity or Lenz vector
        // h = specific angular momentum, independent of mass
        // n = node vector, 2D
        
        long double r, v2, rdotv;   // position mag, velocity mag squared, r dot v
        long double hx2hy2, e2, ne; // useful temporary variables
        long double b1, b2;         // used when calculating eccentricity vector
        long double h_mag, n_mag;   // magnitude of associated vectors
        
        r = p.norm(); // std::sqrt(x * x + y * y + z * z);
        v2 = v.norm2(); // u * u + v * v + w * w;
        rdotv = dotproduct(p, v);  //  dot_product();
        
        // compute the semi-major axis
        m_sma = r * GM / (GM + GM - v2 * r);
        
        h = crossproduct(p, v);
        hx2hy2 = h.x * h.x + h.y * h.y;
        h_mag = sqrt(hx2hy2 + h.z * h.z);

        /*
        // RAAN
        double myraan = atan2(h.x, -h.y);
        myraan = fmod(myraan, twopi);
        
        // argument of latitude = true anomaly + argument of perigee
        double myu = atan2 ( p.z*h_mag, -p.x*h.y+p.y*h.x );    // Arg. of latitude
        
        double ecosE = 1.0-r/m_sma;
        double esinE = dotproduct(p, v)/sqrt(GM*m_sma);
        e2 = ecosE*ecosE +esinE*esinE;
        double myecc = sqrt(e2);
        // eccentric anomaly
        double myE  = atan2(esinE,ecosE);                           // Eccentric anomaly
        double myM  = fmod(myE-esinE,twopi);                          // Mean anomaly
        double nu = atan2(sqrt(1.0-e2)*esinE, ecosE-e2);            // True anomaly(0--2pi)
        if(nu < 0.0)
        {
            nu += twopi;
        }
        double mu = fmod(myu-nu,twopi);                          // Arg. of perihelion
        */
        
        // compute the eccentricity or Lenz vector
        b1 = v2 / GM - 1.0 / r;
        b2 = rdotv / GM;
        
        //e = b1*p - b2*v;
        
        e.x = b1 * p.x - b2 * v.x;
        e.y = b1 * p.y - b2 * v.y;
        e.z = b1 * p.z - b2 * v.z;
        //e2 = e[0] * e[0] + e[1] * e[1] + e[2] * e[2];
        //ecc = sqrt(e2);
        e2= e.norm2();
        m_ecc = e.norm();
        
        // compute the node vector, [0, 0, 1] x h
        n[0] = -h.y;
        n[1] = h.x;
        n_mag = sqrt(hx2hy2);
        
        // compute the inclination
        m_inc = acos(h.z / h_mag);
        
        // compute the argument of perigee
        ne = sqrt(hx2hy2 * e2);
        
        if (ne > 0.0)
        {
            m_argp = acos((n[0] * e.x + n[1] * e.y) / ne);
            if (e.z < 0.0)
            {
                m_argp = twopi - m_argp;
            }
        }
        else
        {
            m_argp = 0; // Prevents NaN for exactly circular orbits
        }
        
        // compute the right ascension of the ascending node
        if ( n_mag > 0.0)
        {
            m_raan = acos(n[0] / n_mag);
            if( n[1] < 0.0 )
            {
                m_raan = twopi - m_raan;
            }
        }
        else
        {
            m_raan = 0; // Prevents NaN for exactly equatorial orbits
        }
        
        // compute the true anomaly
        m_tran = acos( dotproduct(e, p) / (m_ecc * r));
        if (rdotv < 0.0)
        {
            m_tran = twopi - m_tran;
        }
        
    }
    
    
    long double GKeplerianElements::get_tran_from_ecan(long double ecan, long double ecc)
    {
        
        long double sin_tran, cos_tran;
        
        // Not named correctly, these are both missing a division by (1 - ecc*cos(ecan))
        sin_tran = sqrt(1.0L - ecc * ecc) * sin(ecan);
        cos_tran = cos(ecan) - ecc;
        
        // But we don't do the division as it will cancel out in here anyway
        double tran = atan2(sin_tran, cos_tran);
        
        if (tran < 0.0L)
        {
            double twopi = GCONST("2PI");
            tran += twopi; // Bring E into range 0 to 2PI
        }
        
        return tran;
    }
    
    long double GKeplerianElements::get_ecan_from_tran(long double ecc, long double tran)
    {
        long double sin_ecan, cos_ecan;
        double twopi = GCONST("2PI");
        // Not named correctly, these are both missing a division by (1 + ecc*cos(tran))
        sin_ecan = sqrt(1.0 - ecc * ecc) * sin(tran);
        cos_ecan = cos(tran) + ecc;
        
        // But we don't do the division as it will cancel out in here anyway
        double ecan = atan2(sin_ecan, cos_ecan);
        
        if ( ecan < 0.0)
        {
            ecan += twopi; // Bring E into range 0 to 2PI
        }
        return ecan;
    }
    
    // get mean motion from eccentric anomaly
    long double GKeplerianElements::get_mean_from_ecan(long double ecan,long double ecc)
    {
        double mean = ecan - ecc * sin(ecan);
        return mean;
    }
    
    // get eccentric anomaly from mean motion
   long double GKeplerianElements::get_ecan_from_mean(long double mean, long double ecc)
    {
        long double E1=0.0L, E2=0.0L; // initial estimate, improved estimate
        
        //long double tol = 1.0E-15L; // maximum difference allowed
        long double tol = 0x1.0p-50L; // maximum difference allowed 2^(-50)
        
        int count = 0;
        E2 = mean;
        
        do
        {
            E1 = E2;
            E2 = E1 - (E1 - ecc*sin(E1) - mean)/(1.0L-ecc*cos(E1));
            count++;
        } while( (fabs(E2-E1)>tol) || (count > 5) );
        
        if(count > 5)
        {
            printf("warning: %lf\n", fabs(E1-E2));
        }
//        if ((fabs(mean) < 0.7) ||
//            (fabs(mean) > 5.6)) // if M < 40 degrees or M > 320 degrees
//        {
//            //apply Newton's method
//            do
//            {
//                E1 = E2;
//                E2 = E1 - (E1 - ecc * sin(E1) - mean) / (1.0 - ecc * cos(E1));
//            } while ( fabs(E1 - E2) > tol); // end of loop for Newton's method
//        }
//        else
//        {
//            //apply fixed point iteration
//            do
//            {
//                E1 = E2;
//                E2 = mean + ecc * sin(E1);
//            } while ( fabs(E1 - E2) > tol); // end of loop for fixed point iteration
//        }
        
        return E2;
    }
    
    
}

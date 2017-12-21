//
//  GFMEarthGravity.hpp
//  GFC
//
//  Created by lizhen on 09/04/2017.
//  Copyright © 2017 lizhen. All rights reserved.
//

#ifndef GFMEarthGravity_hpp
#define GFMEarthGravity_hpp

#include <stdio.h>

#include "GForceModel.hpp"
#include "GTidalCorrection.hpp"

namespace gfc
{
    /* the class for calculating the harmonic spherical gravity*/
    class GFMEarthGravity : public GForceModel
    {
        
    public:
        GFMEarthGravity() {setForceName("GFMGravity");}
        GFMEarthGravity(int n, int m) // n: degress m: order
        {
            setForceName("GFMGravity");
            
            m_desiredDegree = n;
            m_desiredOrder = m;
            
            //for the purpose of partial derivatives, it is always 2 degree higher
            m_nn = m_desiredDegree + 2;
            m_mm = m_nn;
            size = (m_nn+1)*(m_nn+2)/2;
            V.resize(size);
            W.resize(size);
            C.resize(size);
            S.resize(size);
            N.resize(size);
            
            solid_tide = false;
            time_variable = false;
            ocean_tide = false;
            
        }
        
        virtual ~ GFMEarthGravity() {}
        
        
        void getNormalFactor();
       
        void denormaliseCS(std::vector<double>& Cn, std::vector<double>& Sn);
        
        void loadEGM2008(gfc::GString filename);
        void loadJGM3(gfc::GString filename);
        
        void generateVW(GVector& pos_ecef);
        
        GVector calculateAcc(GVector& pos_ecef);
        
        void getPartialDerivatives(GVector& pos_ecef);
        
        void correctTideCS(GTime ctUTC,double ut1mutc);
        
        // must implement this function
        void doCompute(GTime ctUTC,double ut1mutc, GVector& pos_ecef  );
        
        bool time_variable;
        bool solid_tide;
        bool ocean_tide;
        bool polar_tide;
        
    private:
        
        
        //these variables are for solid earth tide corrections
        // ref IERS2010 techn Note36 Page80, equation 6.4
        //double refMJD = 46431.0; // reference MJD of the force model
        //double dotC20;
        //double dotC21;
        //double dotC30;
        //double dotC40;
        
        //double dotS21;
        
        double GM;
        double R;
        
        
        // the index for these coefficient is i*n+j, n 行m 列
        // index reference: Montenbruck, Page 67
        std::vector<double> V;
        std::vector<double> W;
        std::vector<double> C;
        std::vector<double> S;
        std::vector<double> N; // the denormalise coefficients
        
        // the copy of the denomalised coefficients
        std::vector<double> Cp;
        std::vector<double> Sp;
        
        int m_desiredDegree; // the current degree n
        int m_desiredOrder;  // the current order m
        
        int m_mm;  // for order
        int m_nn;  // for degree
        int size;  // the total size of V,W,C,S,and N vector
    };
    
    
    
} // end of namespace



#endif /* GFMEarthGravity_hpp */

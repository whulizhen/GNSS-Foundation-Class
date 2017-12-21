
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

//
//  GNSSDT.h
//  GFC
//
//  Created by lizhen on 15/10/18.
//  Copyright © 2015年 lizhen. All rights reserved.
//

#ifndef GNSSDATA_H
#define GNSSDATA_H

#include <map>
#include <fstream>
#include <algorithm>

#include "GOBSData.h"
#include "ObsDataTree.h"
#include "GSensor.h"
#include "GTime.h"
#include "GSourceID.h"

namespace gfc
{
    
    
    // GNSSSDT class must inherited from IDATATYPE
    class GgnssOBSType : public GOBSType
    {
        // how to implement the construction function ?
        
    public:
        
        GgnssOBSType()
        {
            m_trackingCode =' ';
            m_freq = NULL;
        }
        
//        // child should call the construction function of the Father
//        GgnssOBSType( GString myDataTypeName, GString trackingCode, GString freqName) : GOBSType(myDataTypeName )
//        {
//            m_trackingCode = trackingCode;
//            m_freqName = freqName;
//        }
        
        virtual ~GgnssOBSType() {}
        
        
        char m_trackingCode;  // tracking code
        
        const GString* m_freq;  // frequency number
        const GString* m_dataType; // the data type name
        
    private:
        
        static GString trackingCodeList; // all the possible tracking code for gnss observables，splitted by a space
        
    };
    
    
    // this class is the basic data structure of obs data tree
    // GNSSDATA is directly connected to Rinex related classes
    /*comments
     * sizeof(gnssOBSData) = 144 bytes, That is really large!!!
     * Assuming we can observe 20 satellites with 4 types observables(3 frequencies),
     * then the memory should be about 100Mb for 2880 epochs with intervals 30 seconds
     * If we store this data with std::vector, the memroy maybe doubled, That's 200Mb. The text file is about 30Mb
     * If we store this data with array, the memory should be about 6Mb.
     * Is there any method to reduce the memory ???
     * compare with gpstk, the sizeof class gnssRinex is 272 byes!!!
     * So it is not a good idea to store all the data in a whole day, It is better to process them in a data stream on behalf of the memory pressure !!!
     */
    class GgnssOBSData: public GOBSData
    {
        
    public:
        
        GgnssOBSData()
        {
            m_SSI = char(0);
            m_LLI = char(0);
        }
        
        //call the constructor of the father class
        GgnssOBSData(double dataValue, double dataStdDev,char ssi, char lli)
        {
            m_SSI = ssi;
            m_LLI = lli;
        }
        
        ~GgnssOBSData() { }
        
        void setSLI(char ssi, char lli) {m_SSI = ssi; m_LLI = lli;}
        void setSSI(char ssi)   {m_SSI = ssi;}
        void setLLI(char lli)   {m_LLI = lli;}
        void getSLI(char& ssi, char& lli) {ssi = m_SSI; lli = m_LLI;}
       
        //GgnssOBSType m_obstype;
        // tow more member variables except for data and std dev
        char m_SSI;  // Signal Strength Index
        char m_LLI;  // Lost of Lock Index
        
    };
    
    
    
    //template< class T, class OBSTYPE, class OBSDATA >
    class GgnssDataEpoch
    {
        
    public:
        
        GgnssDataEpoch()
        {
            int nsat = 100;
            int nobsdata = 50;
            m_indicatorlist.reserve(nsat);
            m_obstypelist.reserve(nobsdata*nsat);
            m_obsdatalist.reserve(nobsdata*nsat);
        };
        
        /// Copy constructor.
        //template<class OBSTYPE, class OBSDATA>
        GgnssDataEpoch( const GgnssDataEpoch & g )
        {
            m_indicatorlist = g.m_indicatorlist;
            m_obstypelist = g.m_obstypelist;
            m_obsdatalist = g.m_obsdatalist;
        }
        
        // data members
        GTime m_epoch; // epoch time, include the clock error
        
        std::vector< GSensorID > m_indicatorlist;  // satellite list
        
        // every satellite has the same obstype list
        
        // every obsdata has one obstype
        std::vector<  std::vector<GgnssOBSType> > m_obstypelist;
        
        std::vector<  std::vector<GgnssOBSData> > m_obsdatalist;
        
    };

    //template<class TA, class TB, class OBSTYPE, class OBSDATA >
    class GgnssStorage : public std::map<GSourceID, GgnssDataEpoch >
    {
        
    public:
        //remove TB
        void removeGSourceID(GSourceID& sourceID)
        {
            (*this).erase(sourceID);
        }
        
        
        void removeGSensorID( GSensorID&  sensorID )
        {
            std::map<GSourceID, GgnssDataEpoch >::iterator myit;
            
            for( myit = (*this).begin(); myit != (*this).end() ; myit++ )
            {
                std::vector<GSensorID>::iterator it =
                std::find(myit->second.m_indicatorlist.begin(), myit->second.m_indicatorlist.end(),sensorID);
                //find(myit->second.m_indicatorlist.begin(),myit->second.m_indicatorlist.end(),sensorID);  //  .find(sensorID);
                
                if( it == myit->second.m_indicatorlist.end() )
                {
                    printf("GOBSStorage: satellite does Not exist in removeGSensorID function\n");
                }
                else
                {
                    
                    int t =  it - myit->second.m_indicatorlist.begin();
                    
                    myit->second.m_indicatorlist.erase(it);
                    
                    
                    myit->second.m_obstypelist.erase(myit->second.m_obstypelist.begin()+t);
                    
                    myit->second.m_obsdatalist.erase(myit->second.m_obsdatalist.begin()+t);
                    
                }
                
            }
        }
        
        //remove the sensorID at sourceID
        void removeGSensor_GSource(GSensorID& sensorID, GSourceID& sourceID)
        {
            
        }
        
        
    };  //  end of class OBSMgr
    

    
    
} // end of namespace gfc


#endif /* GNSSDT_h */

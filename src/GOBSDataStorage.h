//
//  GOBSDataMgr.h
//  GFC
//
//  Created by lizhen on 23/08/2016.
//  Copyright © 2016 lizhen. All rights reserved.
//

#ifndef GOBSDataMgr_h
#define GOBSDataMgr_h

#include <algorithm>

#include "GSourceID.h"
#include "GTime.h"

namespace gfc
{
    
    
//    template< class T, class OBSTYPE, class OBSDATA >
//    class GDataEpoch
//    {
//        
//    public:
//        
//        GDataEpoch()
//        {
//            int nsat = 100;
//            int nobsdata = 50;
//            m_indicatorlist.reserve(nsat);
//            m_obstypelist.reserve(nobsdata*nsat);
//            m_obsdatalist.reserve(nobsdata*nsat);
//        };
//        
//        /// Copy constructor.
//        //template<class OBSTYPE, class OBSDATA>
//        GDataEpoch( const GDataEpoch< T, OBSTYPE,OBSDATA>& g )
//        {
//            m_indicatorlist = g.m_indicatorlist;
//            m_obstypelist = g.m_obstypelist;
//            m_obsdatalist = g.m_obsdatalist;
//        }
//        
//        // data members
//        GTime m_epoch; // epoch time, include the clock error
//        
//        std::vector< T > m_indicatorlist;  // satellite list
//        
//        // every satellite has the same obstype list
//        //std::vector< OBSTYPE > m_obstypelist;
//        
//        // every obsdata has one obstype
//        std::vector<  std::vector<OBSTYPE> > m_obstypelist;
//        
//        std::vector<  std::vector<OBSDATA> > m_obsdatalist;
//        
//    };
//    
//    
//    
//    template<class TA, class TB, class OBSTYPE, class OBSDATA >
//    class GOBSStorage : public std::map<TB, GDataEpoch< TA, OBSTYPE,OBSDATA> >
//    {
//        
//    public:
//        //remove TB
//        void removeGSourceID(TB& sourceID)
//        {
//            (*this).erase(sourceID);
//        }
//        
//        
//        void removeGSensorID( TA&  sensorID )
//        {
//            typename std::map<TB, GDataEpoch< TA, OBSTYPE,OBSDATA> >::iterator myit;
//            
//            for( myit = (*this).begin(); myit != (*this).end() ; myit++ )
//            {
//                typename std::vector<TA>::iterator it =
//                find(myit->second.m_indicatorlist.begin(),myit->second.m_indicatorlist.end(),sensorID);  //  .find(sensorID);
//                if( it == myit->second.m_indicatorlist.end() )
//                {
//                    printf("GOBSStorage: satellite does Not exist in removeGSensorID function\n");
//                }
//                else
//                {
//                    
//                    int t =  it - myit->second.m_indicatorlist.begin();
//                    
//                    myit->second.m_indicatorlist.erase(it);
//                   
//                    
//                    myit->second.m_obstypelist.erase(myit->second.m_obstypelist.begin()+t);
//                    
//                    myit->second.m_obsdatalist.erase(myit->second.m_obsdatalist.begin()+t);
//                
//                }
//                
//            }
//        }
//        
//        //remove the sensorID at sourceID
//        void removeGSensor_GSource(TA& sensorID, TB& sourceID)
//        {
//            
//        }
//        
//        
//    };  //  end of class OBSMgr
//    
    
    
    
}



#endif /* GOBSDataMgr_h */

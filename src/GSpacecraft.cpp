//
//  GSpacecraft.cpp
//  GFC
//
//  Created by lizhen on 16/4/14.
//  Copyright © 2016年 lizhen. All rights reserved.
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


#include "GSpacecraft.hpp"

namespace gfc
{
    
    //the declaration of the static variable
    std::map<GString,GSpaceCraft> GSpaceCraftMgr::gSpacecraft;
    int GSpaceCraftMgr::GPSMAXPRN = 32;
    int GSpaceCraftMgr::BDSMAXPRN = 35;
    int GSpaceCraftMgr::GLSMAXPRN = 30;
    int GSpaceCraftMgr::GALMAXPRN = 30;
    
    
    GSpaceCraft& GSpaceCraft::operator= (const GSpaceCraft& right)
    {
        this->m_scID = right.m_scID;
        this->m_state = right.m_state;
        this->m_pephStore = right.m_pephStore;
        this->m_bodyModel = right.m_bodyModel;
        return *this;
    }
    
    //copy constructor
    GSpaceCraft::GSpaceCraft(const GSpaceCraft& right)
    {
        this->m_scID = right.m_scID;
        this->m_state = right.m_state;
        m_pephStore = right.m_pephStore;
        m_bodyModel = right.m_bodyModel;
    }
    
    //put in the precise ephemeris
    void GSpaceCraft::pushPreciseEphemeris(GPreciseEphemeris pehp)
    {
        GTime t = pehp.getEpoch();
        m_pephStore[t] = pehp;
        //m_pephStore.insert(map<GTime, GPreciseEphemeris>::value_type(t,pehp) );
        
    }
    
    
    GPreciseEphemeris GSpaceCraft::getEphValue(int index)
    {
        std::map<GTime, GPreciseEphemeris>::iterator myit = m_pephStore.begin();
        for( int i = 0 ; i< index; i++)
        {
            myit++;
        }
        //myit = myit + index;
        return myit->second;
    }
    
    int GSpaceCraft::getNumofPreEph()
    {
        return m_pephStore.size();
    }
    
    
    
    int GSpaceCraft::getEphLocation(GTime& ttag, const int& nhalf,
                        std::map<GTime,GPreciseEphemeris>::const_iterator & it1,
                        std::map<GTime,GPreciseEphemeris>::const_iterator& it2,
                        bool& exactMatch,
                        bool exactReturn= true
                        )
    {
        
        int retval = 0;
        double myinterval = 900.0; // 9000 seconds
        bool checkDataGap = true;
        bool checkInterval = false;
        double gapInterval = 3600.0; // one hour is assumed a datagap
        double maxInterval = (2*nhalf-1)*myinterval;
        
        // cannot interpolate with one point
        if( m_pephStore.size() < 2 )
        {
            InvalidRequest e("Inadequate data (size < 2) for satellite "
                             + m_scID.getIDString() );
            //GFC_THROW(e);
            retval = 1;
            return retval;
        }
        
        // find the timetag in this table
        // NB. throw here if time systems do not match and are not "Any"
        // note it1 can always be the end of the map
        it1 = m_pephStore.find(ttag);
        
        // is it an exact match?
        exactMatch = (it1 != m_pephStore.end());
        
        // user must decide whether to return with exact value; e.g. without
        // velocity data, user needs the interval to compute v from x data
        if(exactMatch && exactReturn) return true;
        
        // lower_bound points to the first element with key >= ttag
        it1 = it2 = m_pephStore.lower_bound(ttag);
        
        // ttag is <= first time in table
        if( it1 == m_pephStore.begin() )
        {
            // at table begin but its an exact match && an interval of only 2
            if(exactMatch && nhalf==1)
            {
                ++(it2 = it1);
                return exactMatch;
            }
            
            InvalidRequest e("Inadequate data before(1) requested time for satellite "
                             + m_scID.getIDString() );
            retval = 2;
            return retval;
            GFC_THROW(e);
        }
        
        // move it1 down by one
        if( --it1 == m_pephStore.begin() )
        {
            // if an interval of only 2
            if(nhalf==1)
            {
                ++(it2 = it1);
                return exactMatch;
            }
            InvalidRequest e("Inadequate data before(2) requested time for satellite "
                             + m_scID.getIDString() );
            
            retval = 3;
            return retval;
            GFC_THROW(e);
        }
        
        if( checkDataGap && ( ( it2->first - it1->first ).toSeconds() ) > gapInterval)
        {
            InvalidRequest e("Gap at interpolation time for satellite "
                             + m_scID.getIDString() );
            printf("data gap found in precise ephemeris \n");
            
            retval = 4;
            return retval;
            GFC_THROW(e);
        }
        
        // now expand the interval to include 2*nhalf timesteps
        for(int k=0; k<nhalf-1; k++)
        {
            bool last(k==nhalf-2);        // true only on the last iteration
            // move left by one; if require full interval && out of room on left, fail
            if( --it1 == m_pephStore.begin() && !last)
            {
                InvalidRequest
                e("Inadequate data before(3) requested time for satellite "
                  + m_scID.getIDString());
                retval = 5;
                return retval;
                GFC_THROW(e);
            }
            
            if(++it2 == m_pephStore.end())
            {
                if(exactMatch && last && it1 != m_pephStore.begin())
                {
                    // exact match && at end of interval && with room to move down
                    it2--; it1--;  // move interval down by one
                }
                else
                {
                    InvalidRequest
                    e("Inadequate data after(2) requested time for satellite "
                      + m_scID.getIDString() );
                    
                    retval =6;
                    return retval;
                    
                    GFC_THROW(e);
                }
            }
            //LOG(INFO) << k << " expand right " << printTime(it2->first,"%F/%g");
        }

        // check that the interval is not too large
        if( checkInterval && ( (it2->first - it1->first ).toSeconds() ) > maxInterval)
        {
            InvalidRequest e("Interpolation interval too large for satellite "
                             + m_scID.getIDString() );
            
            retval = 7;
            return retval;
            GFC_THROW(e);
        }

        return retval;
        
    }
    
    /*
     
     
     */
    GPreciseEphemeris GSpaceCraft::getEphValue_test(GTime ct)
    {
        if(m_pephStore.size() == 0 )
        {
            printf("no precise ephemeris for satellite: %s\n", m_scID.getIDString().c_str());
            exit(1);
        }
        
        GPreciseEphemeris rec;
        
        std::map<GTime,GPreciseEphemeris>::const_iterator it1, it2,it3, kt;
        int NHALF = 5;
        int nhalf = NHALF;
        //check whether ct is between the time domain
        bool haveVelocity =false, exactMatch=false;
        GTime ct0; // the time for the first element
        it1 = m_pephStore.begin();
        it2 = m_pephStore.end(); it2--;  // it2 is the last element
        
        if( ct < it1->first || ct > it2->first )
        {
            return rec;
        }
        
        double velTest =  it1->second.m_u*it1->second.m_u + it1->second.m_v*it1->second.m_v + it1->second.m_w*it1->second.m_w;
        if( velTest < 1.0E-10 )
        {
            haveVelocity = false;
        }
        
        kt = m_pephStore.find(ct);
        
        // is it an exact match?
        exactMatch = (kt != m_pephStore.end());
        
        if(haveVelocity && exactMatch)
        {
            return kt->second;
        }
        
        if(exactMatch == true && kt== m_pephStore.begin())
        {
            return kt->second;
        }
        
        if( exactMatch == true && kt == it2 )
        {
            
            return kt->second;
        }
        //because kt++, here recovey to the original state
        
        // lower_bound points to the first element with key >= ct
        it1 = it3 = m_pephStore.lower_bound(ct);
        
        if(!exactMatch)
        {
            kt = it1;
        }
        
        
        //expand
        //int count_before = exactMatch==1?0:1;
        int count_before = 0;
        
        int count_after = exactMatch==1?0:1;
        
        while(1)
        {
            count_before++;
            it1--;
            if(it1==m_pephStore.begin())
            {
                break;
            }
            if(count_before>=nhalf)
            {
                break;
            }
        }
        
       
        while(1)
        {
            count_after++;
            it3++;
            if(it3==m_pephStore.end())
            {
                break;
            }
            if(count_after>=nhalf)
            {
                break;
            }
        }
        
        //get the minimum number in count_before and count_after
        nhalf = count_before<=count_after? count_before: count_after;
        
        if(nhalf < NHALF)
        {
            printf("%s: Lagrange Interpolation Accuracy Deficiency!\n", GTime::GTime2CivilTime(ct).TimeString().c_str());
            if(exactMatch)
            {
               return kt->second;
            }
            else
            {
                return rec;
            }
        }
        
        vector<double> times,P[3],V[3],A[3],sigP[3],sigV[3],sigA[3];
        //pull out the data
        int count =0;
        it1 = it3 = kt; // get the current position
        while(count++<nhalf)
        {
            it1--;
        }
        //the time at the begining of the series
        ct0 = it1->second.m_epoch;
        count = 0;
        int num = exactMatch?2*nhalf+1 : 2*nhalf;
        while ( count++ < num )
        {
            times.push_back((it1->second.m_epoch-ct0).toSeconds());
            it1->second.m_x;
            P[0].push_back(it1->second.m_x);
            P[1].push_back(it1->second.m_y);
            P[2].push_back(it1->second.m_z);
            V[0].push_back(it1->second.m_u);
            V[1].push_back(it1->second.m_v);
            V[2].push_back(it1->second.m_w);
            it1++;
        }
        
        double dt = (ct -ct0).toSeconds(), err;
        double pv_int[6]={0.0};
        //start interpolating
        if(haveVelocity) // get the velocity from interpolation
        {
            for(int i = 0 ; i< 3; i++)
            {
                pv_int[i] = GMath::LagrangeInterpolation(times,P[i],dt,err);
            }
            for(int i = 0 ; i< 3; i++)
            {
                pv_int[3+i] = GMath::LagrangeInterpolation(times,V[i],dt,err);
            }
        }
        else  // get the velocity from position
        {
            
            for(int i = 0 ; i< 3; i++ )
            {
                GMath::LagrangeInterpolation(times,P[i],dt,pv_int[i],pv_int[3+i] );
            }
            if(exactMatch) //keep position the same
            {
                pv_int[0] = kt->second.m_x;
                pv_int[1] = kt->second.m_y;
                pv_int[2] = kt->second.m_z;
            }
        }
            
        rec.setEpoch(ct);
        rec.setPos(pv_int);
        rec.setVel(pv_int+3);
        rec.setOK(true);
        
        return rec;
    }
    
    
    
    // the other way of getEphValue
    // ref: gpstk, PositionRecord::getValue
    GPreciseEphemeris GSpaceCraft::getEphValue(GTime ct)
    {
        
        if(m_pephStore.size() < 2 )
        {
            printf("precise ephemeris number not enough, number: %d\n", m_pephStore.size());
        }
        int Nhalf = 5;
        GPreciseEphemeris rec;
        bool isExact = false;
        int i = -1;
        bool haveVelocity=false, haveAcceleration = false;
        std::map<GTime,GPreciseEphemeris>::const_iterator it1, it2, kt;
        
        //check whether ct is between the time domain
        
        it1 = m_pephStore.begin();
        //it1++;
        it2 = m_pephStore.end();
        it2--;
        
        if( ct < it1->first || ct > it2->first )
        {
            return rec;
        }
        
        
        double velTest =  it1->second.m_u*it1->second.m_u + it1->second.m_v*it1->second.m_v + it1->second.m_w*it1->second.m_w;
        if( velTest < 1.0E-10 )
        {
            haveVelocity = false;
        }
        
        int retval = getEphLocation(ct, Nhalf, it1, it2, isExact);
        
        if(retval != 0 && isExact == false) // ephemeris is unavailable
        {
            return rec;
        }
        
        if( isExact && haveVelocity)
        {
            rec = it1->second;
            return rec;
        }
        
        
        // pull data out of the data table
        int n,Nlow(Nhalf-1),Nhi(Nhalf),Nmatch(Nhalf);
        GTime ct0(it1->first);
        vector<double> times,P[3],V[3],A[3],sigP[3],sigV[3],sigA[3];
        
        double pv[6]= {0.0};
        
        kt = it1; n=0;
        while(1)
        {
            //GTime tk(kt->first);
            
            // find index matching ttag
            if( isExact && fabs( (kt->first - ct).toSeconds() ) < 1.e-8 )
            {
                Nmatch = n;
            }
            
            times.push_back( ( kt->first - ct0).toSeconds() );          // sec
            
            pv[0] = kt->second.m_x; pv[1] = kt->second.m_y; pv[2] = kt->second.m_z;
            pv[3] = kt->second.m_u; pv[4] = kt->second.m_v; pv[5] = kt->second.m_w;
            
            //kt->second.getPV(pv);
            
            for(i=0; i<3; i++ )
            {
                P[i].push_back(pv[i] );
                V[i].push_back(pv[i+3]);
                //A[i].push_back(kt->second.Acc[i]);
                //sigP[i].push_back(kt->second.sigPos[i]);
                //sigV[i].push_back(kt->second.sigVel[i]);
                //sigA[i].push_back(kt->second.sigAcc[i]);
            }
            if(kt == it2) break;
            ++kt;
            ++n;
        };
        
        if(isExact && Nmatch == (int)(Nhalf-1)) { Nlow++; Nhi++; }
        
        // Lagrange interpolation
        //rec.sigAcc = rec.Acc = Triple(0,0,0);        // default
        double dt= (ct-ct0).toSeconds(), err;           // dt in seconds
        
        double pos_inter[3] = {0.0};
        double vel_inter[3] = {0.0};
        
        if(haveVelocity)
        {
            for(i=0; i<3; i++)
            {
                // interpolate the positions
                pos_inter[i] = GMath::LagrangeInterpolation(times,P[i],dt,err);
                if( haveAcceleration)
                {
                    // interpolate velocities and acclerations
                    vel_inter[i] = GMath::LagrangeInterpolation(times,V[i],dt,err);
                    //rec.Acc[i] = LagrangeInterpolation(times,A[i],dt,err);
                }
                else
                {
                    // interpolate velocities(dm/s) to get V and A
                    double inter_acc = 0.0;
                    GMath::LagrangeInterpolation(times,V[i],dt,vel_inter[i],inter_acc);
                    inter_acc *= 0.1;      // dm/s/s -> m/s/s
                }
                
                /*
                if(isExact)
                {
                    
                    rec.sigPos[i] = sigP[i][Nmatch];
                    rec.sigVel[i] = sigV[i][Nmatch];
                    if(haveAcceleration) rec.sigAcc[i] = sigA[i][Nmatch];
                }
                else
                {
                    // TD is this sigma related to 'err' in the Lagrange call?
                    rec.sigPos[i] = RSS(sigP[i][Nhi],sigP[i][Nlow]);
                    rec.sigVel[i] = RSS(sigV[i][Nhi],sigV[i][Nlow]);
                    if(haveAcceleration)
                        rec.sigAcc[i] = RSS(sigA[i][Nhi],sigA[i][Nlow]);
                }
                */
                
                // else Acc=sig_Acc=0   // TD can we do better?
            }
        }
        else
        {               // no V data - must interpolate position to get velocity
            for( i=0; i<3; i++ )
            {
                // interpolate positions(km) to get P and V
                GMath::LagrangeInterpolation(times,P[i],dt,pos_inter[i],vel_inter[i] );
                //rec.Vel[i] *= 10000.;         // km/sec -> dm/sec
                
                /*
                if(isExact)
                {
                    rec.sigPos[i] = sigP[i][Nmatch];
                }
                else {
                    rec.sigPos[i] = RSS(sigP[i][Nhi],sigP[i][Nlow]);
                }
                // TD
                rec.sigVel[i] = 0.0;
                */
            }
        }
        
        rec.setEpoch(ct);
        rec.setPos(pos_inter);
        rec.setVel(vel_inter);
        rec.setOK(true);
        
        return rec;
    }
    
    
    
    /*
     *
     * get the positon and velocity of time ct, GPST
     *
     */
    GPreciseEphemeris GSpaceCraft::getEphValue1( GTime ct )
    {
        
        GPreciseEphemeris rec;
        rec.setEpoch(ct);
        rec.setSCID(m_scID);
//        
//        const int Nhalf = 4; // 10 degree lagrange interplation
//        
//        GTime ep0 = m_pephStore[0].getEpoch();
//        GTime ep1 = m_pephStore[1].getEpoch();
//        GTime epEnd = m_pephStore[m_pephStore.size()-1].getEpoch();
//        
//        if(ct < ep0 || ct > epEnd)
//        {
//            return rec;
//        }
//        
//        double interval = (  ep1 - ep0 ).toSeconds();
//        
//        double dummyTime = (ct - ep0).toSeconds();
//        double dummy = dummyTime/interval;
//        int index0 = (int)(dummy);
//        
//        double pv[6]={0.0};
//        m_pephStore[index0].getPV(pv);
//        double test  = pv[3]*pv[3] + pv[4]*pv[4] + pv[5]*pv[5] ;
//        bool hasvelocity = false;
//        bool isExact = false;
//        if( fabs(test) > 1.0E-8 )
//        {
//            hasvelocity = true;
//        }
//        
//        if( fabs(dummy - index0 ) <1.0E-10 )  //find the exact epoch
//        {
//            if( ct == m_pephStore[index0].getEpoch() )
//            {
//                isExact = true;
//                if( hasvelocity )
//                {
//                    return m_pephStore[index0];
//                }
//            }
//            else
//            {
//                ObjectNotFound e("error in find the exact precise ephemeris");
//                //return rec;
//                GFC_THROW(e);
//            }
//        }
//        
//        
//        // cannot interpolate with one point
//        if( m_pephStore.size() < 2 )
//        {
//            InvalidRequest e("Inadequate precise ephemeris data (size < 2) for satellite "
//                             + m_scID.getIDString() );
//            return rec;
//            //GFC_THROW(e);
//        }
//
//        
//        int indexStart = index0; // the start index for interplation
//        int indexEnd   = index0; // the end   index for interplation
//        
//        if( isExact == true )
//        {
//            indexStart = index0 ;
//            indexEnd = index0+1;
//        }
//        else
//        {
//            indexStart = index0;
//            indexEnd = index0 + 1;
//        }
//        
//        GTime timeStart = m_pephStore[indexStart].getEpoch();
//        GTime timeEnd = m_pephStore[indexEnd].getEpoch();
//        
//        //from indexStart forward
//        for( int it = 0 ; it< Nhalf ; ++it )
//        {
//            //int i = indexStart - it;
//           
//            if( indexStart < 0 )
//            {
//                InvalidRequest e("Inadequate data before requested time for satellite "
//                                 + m_scID.getIDString() );
//                return rec;
//                //GFC_THROW(e);
//            }
//            GTime t = timeStart - interval*it;
//            double difTime =  ( m_pephStore[indexStart].getEpoch() - t ).toSeconds();
//            if( fabs(difTime) > 1.0E-8 )  // time is different, means data gap
//            {
//                InvalidRequest e(" data point is missing at interpolation time for satellite: "
//                                 + m_scID.getIDString() );
//                return rec;
//                //GFC_THROW(e);
//            }
//            
//             indexStart = indexStart - 1;
//        }
//        
//        //from the indexEnd backward
//        for( int it = 0 ; it< Nhalf ; it++ )
//        {
//            //int i = indexEnd + it;
//            if( indexEnd > m_pephStore.size() )
//            {
//                InvalidRequest e("Inadequate data before requested time for satellite "
//                                 + m_scID.getIDString() );
//                return rec;
//                //GFC_THROW(e);
//            }
//            GTime t = timeEnd + interval*it;
//            if( fabs( ( m_pephStore[indexEnd].getEpoch() - t ).toSeconds()) > 1.0E-8 )  // time is different, means data gap
//            {
//                InvalidRequest e(" data point is missing at interpolation time for satellite "
//                                 + m_scID.getIDString() );
//                return rec;
//                //GFC_THROW(e);
//            }
//            indexEnd = indexEnd + 1;
//        }
//        
//        
//        //if everything is good, then get all the data points
//        //std::vector<double> pX, pY, pZ,vX,vY,vZ;
//        std::vector<double> Times;
//        std::vector<double> P[3],V[3];
//        GTime time0 = m_pephStore[indexStart].getEpoch();
//        Times.resize(indexEnd-indexStart + 1);
//        
//        double interestedTime = (ct - time0).toSeconds(); // the current interesed time to interplate
//        
//        for( int i =0 ; i< 3; i++ )
//        {
//            P[i].resize(indexEnd - indexStart + 1);
//            V[i].resize(indexEnd - indexStart + 1);
//        }
//        
//        //pX.resize(indexEnd - indexStart + 1); pY.resize(indexEnd - indexStart + 1); pZ.resize(indexEnd - indexStart + 1);
//        //vX.resize(indexEnd - indexStart + 1); vY.resize(indexEnd - indexStart + 1); vZ.resize(indexEnd - indexStart + 1);
//        
//        for( int i = 0; i<indexEnd-indexStart+1 ; i++ )
//        {
//            double pv[6] = {0.0};
//            m_pephStore[i+indexStart].getPV(pv);
//            Times[i] = (m_pephStore[i+indexStart].getEpoch() - time0).toSeconds();
//            
//            if(hasvelocity == true)
//            {
//                for( int j = 0 ; j<  6 ; j++ )
//                {
//                    if( j < 3 )
//                    {
//                        P[j][i] = pv[j];
//                    }
//                    else
//                    {
//                        V[j-3][i] = pv[j];
//                    }
//                }
//            }
//            else if(hasvelocity == false)
//            {
//                for( int j = 0 ; j<  3 ; j++ )
//                {
//                    P[j][i] = pv[j];
//                }
//            }
//        }
//        
//        if( isExact == true && hasvelocity == false) // just need to interplate position to get velocity
//        {
//            double pos[3]={0.0}, vel[3]={0.0};
//            //LagrangeInterpolation(times,V[i],dt,rec.Vel[i],rec.Acc[i]);
//            for( int i = 0 ; i< 3; i++ )
//            {
//                //double p = 0.0, v = 0.0;
//                GMath::LagrangeInterpolation(Times, P[i], interestedTime,pos[i],vel[i]);
//            }
//            
//            rec.setPX(pv[0]);rec.setPY(pv[1]);rec.setPZ(pv[2]);
//            rec.setVel(vel);
//            rec.setOK(true);
//        }
//        else if(isExact == false)  // both need to interplate position and velocity
//        {
//            double pos[3]={0.0}, vel[3]={0.0};
//            //LagrangeInterpolation(times,V[i],dt,rec.Vel[i],rec.Acc[i]);
//            for( int i = 0 ; i< 3; i++ )
//            {
//                GMath::LagrangeInterpolation(Times, P[i], interestedTime,pos[i],vel[i]);
//            }
//            rec.setPos(pos);
//            rec.setVel(vel);
//            rec.setOK(true);
//        }
//        
       
        
        
        return rec;
    }
    
    
    void GSpaceCraft::bindSpacecraftModel(gfc::GSensorID scID)
    {
       
        GString svType = GSpacecraftModelMgr::sensorID2svType(scID);
        
        m_bodyModel = &GSpacecraftModelMgr::spacecraftModelStore[svType];
        
        int testc = 0;
    }
    
    //get the corrections for LRA or GNSS Antenna PCO
    // the final result is in ECEF because all the measurements are in ECEF
    void GSpaceCraft::getOffsetCorrection( int type, GVector& targetecef)
    {
        double distance_correction = 0.0;
        
        double offset_bfs[3] = {0.0}; // offset in BFS
        for( int i = 0 ; i< 3; i++ )
        {
            if( type == 0 )
            {
                offset_bfs[i] = m_bodyModel->m_offsetGNSS[i] - m_bodyModel->m_com[i];
            }
            else if(type == 1 )
            {
                offset_bfs[i] = m_bodyModel->m_offsetLRA[i] - m_bodyModel->m_com[i];
            }
        }
        
        //convet BFS offset to ECEF, because stationPos is usually in ECEF
        GVector target(offset_bfs[0],offset_bfs[1],offset_bfs[2]);
        
        targetecef =  m_state.attitude_ecef.convert2Target(target);
        
        //because target is much smaller comparing with the distance between satellite and station
        /*
           project the targetecef to line of sight vector
         
        | ->     -> |       | -> |     | -> |
        | a   -  b  |   =   | a  |  -  | b  |
        
         */
        
    }
    
    void GSpaceCraftMgr::Initializer()
    {
        
        GSensorID myid("ssTEST",1);
        GSpaceCraft mysat(myid);
        //set up the geometry information
        //mysat.bindSpacecraftModel(myid);
        gSpacecraft[myid.getIDString()]= mysat;
        
        myid.setData("ssLAGEOS", 52);  // LAGEOS-2 satellite, 51 for LAGEOS-1 satellite
        mysat.setSpaceCraftID(myid);
        gSpacecraft[myid.getIDString()] = mysat;
        
        myid.setData("ssGRACE", 1);  // GRACE-A satellite
        mysat.setSpaceCraftID(myid);
        gSpacecraft[myid.getIDString()] = mysat;
        
        myid.setData("ssGRACE", 2);  // GRACE-B satellite
        mysat.setSpaceCraftID(myid);
        gSpacecraft[myid.getIDString()] = mysat;
        
        
        for( int i = 1; i<= 32; i++ )  //GPS prn:32
        {
            GSensorID myid("ssGPS",i);
            
            GSpaceCraft mysat(myid);
            
            //set up the geometry information
            mysat.bindSpacecraftModel(myid);
            
            gSpacecraft[myid.getIDString()]= mysat;
            
        }
        for( int i = 1; i<= 32; i++ )  //GLS prn:32
        {
            GSensorID myid("ssGLO",i);
            GSpaceCraft mysat(myid);
            
            //set up the geometry information
            mysat.bindSpacecraftModel(myid);
            
            gSpacecraft[myid.getIDString()]= mysat;
        }
        for( int i = 1; i<= 35; i++ ) //BDS prn 35
        {
            GSensorID myid("ssBDS",i);
            GSpaceCraft mysat(myid);
            
            //set up the geometry information
            mysat.bindSpacecraftModel(myid);
            
            gSpacecraft[myid.getIDString()]= mysat;
        }
        
        for( int i = 1; i<= 30; i++ ) //GAL prn 35
        {
            GSensorID myid("ssGAL",i);
            GSpaceCraft mysat(myid);
            
            //set up the geometry information
            mysat.bindSpacecraftModel(myid);
            
            gSpacecraft[myid.getIDString()]= mysat;
        }
        
    }
    
    
    void GSpaceCraftMgr::loadGRACEEphemeris(GString ephemerisFile)
    {
        
        std::fstream infile(ephemerisFile.c_str());
        if ( infile.fail() )
        {
            std::cerr << "\nCould not open GRACE ephemeris file. Terminating...\n";
            exit(0);
        }
        const int MAX = 1024;
        char store[MAX]={0};
        long lineNumber = 0;
        int totalsat = 0;
        GPreciseEphemeris peph;
        while ( !infile.eof() )
        {
            
            infile.getline(store, MAX); //read past the file header
            GString reader(store);
            if(reader == "" )
            {
                break;
            }
            
            std::vector<GString> splitstr = reader.split();
            
            CivilTime ct;
            //default timesystem is GPS time
            ct.m_ts = GTimeSystem("tsGPS");
            ct.m_year = static_cast<int>(splitstr[0].asINT() ) ;
            ct.m_month = static_cast<int>( splitstr[1].asINT());
            ct.m_day = static_cast<int>(splitstr[2].asINT()) ;
            ct.m_hour = static_cast<int>(splitstr[3].asINT()) ;
            ct.m_minute = static_cast<int>(splitstr[4].asINT()) ;
            ct.m_second = splitstr[5].asDOUBLE();
            GTime epoch =  GTime::CivilTime2GTime(ct);
            peph.setEpoch(epoch);
            
            for(int i = 0 ; i< 2; i++)
            {
                GSensorID satid("ssGRACE",1);
                //the second line
                infile.getline(store, MAX); //read past the file header
                reader = store;
                // they are in ECEF frame
                splitstr = reader.split();
                if(splitstr[0] == "graceA")
                {
                    satid.setData("ssGRACE", 1);
                }
                else if(splitstr[0] == "graceB")
                {
                   satid.setData("ssGRACE", 2);
                }
                
                // the current research object is GPS prn 11, SVN 46
                
                peph.setSCID(satid);
                peph.setPX(splitstr[1].asDOUBLE()/1000.0);  //px  km
                peph.setPY(splitstr[2].asDOUBLE()/1000.0);  //px  km
                peph.setPZ(splitstr[3].asDOUBLE()/1000.0);  //px  km
                peph.setVX(splitstr[4].asDOUBLE()/1000.0);  //vx km
                peph.setVY(splitstr[5].asDOUBLE()/1000.0);  //vx km
                peph.setVZ(splitstr[6].asDOUBLE()/1000.0);  //vx km
                
                bool isexist = ( gfc::GSpaceCraftMgr::gSpacecraft.find(peph.getSensorID().getIDString()) != gfc::GSpaceCraftMgr::gSpacecraft.end() );
                if( isexist == true )
                {
                    gfc::GSpaceCraftMgr::gSpacecraft[peph.getSensorID().getIDString()].pushPreciseEphemeris(peph);
                }
                
            }
            
        }
        
    }
    
    void GSpaceCraftMgr::loadUCLEphemeris(GString ephemerisFile)
    {
        std::fstream infile(ephemerisFile.c_str());
        if ( infile.fail() )
        {
            std::cerr << "\nCould not open UCL ephemeris file. Terminating...\n";
            exit(0);
        }
        const int MAX = 1024;
        char store[MAX]={0};
        long lineNumber = 0;
        int totalsat = 0;
        GPreciseEphemeris peph;
        while ( !infile.eof() )
        {
            
            infile.getline(store, MAX); //read past the file header
            GString reader(store);
            if(reader == "" )
            {
                break;
            }
            std::vector<GString> splitstr = reader.split();
            
            // the current research object is GPS prn 11, SVN 46
            GSensorID satid("ssGPS",11);
            
            CivilTime ct;
            //default timesystem is UTC time
            ct.m_ts = GTimeSystem("tsUTC");
            ct.m_year = static_cast<int>(splitstr[1].asINT() ) ;
            ct.m_month = static_cast<int>( splitstr[2].asINT());
            ct.m_day = static_cast<int>(splitstr[3].asINT()) ;
            ct.m_hour = static_cast<int>(splitstr[4].asINT()) ;
            ct.m_minute = static_cast<int>(splitstr[5].asINT()) ;
            ct.m_second = splitstr[6].asDOUBLE();
            GTime epoch =  GTime::CivilTime2GTime(ct);
            
            // they are in ECI frame
            peph.setEpoch(epoch);
            peph.setSCID(satid);
            peph.setPX(splitstr[7].asDOUBLE());  //px  km
            peph.setPY(splitstr[8].asDOUBLE());  //px  km
            peph.setPZ(splitstr[9].asDOUBLE());  //px  km
            peph.setVX(splitstr[10].asDOUBLE());  //vx km
            peph.setVY(splitstr[11].asDOUBLE());  //vx km
            peph.setVZ(splitstr[12].asDOUBLE());  //vx km
            
            bool isexist = ( gfc::GSpaceCraftMgr::gSpacecraft.find(peph.getSensorID().getIDString()) != gfc::GSpaceCraftMgr::gSpacecraft.end() );
            if( isexist == true )
            {
                gfc::GSpaceCraftMgr::gSpacecraft[peph.getSensorID().getIDString()].pushPreciseEphemeris(peph);
            }
            
        }
        
    }
    
    void GSpaceCraftMgr::loadPreciseEphemeris( GString ephemerisFile)
    {
        //GTime t0;
        //GSensorID sc1;
        //m_peStorage[t0][sc1];
        std::fstream infile(ephemerisFile.c_str());
        if ( infile.fail() )
        {
            std::cerr << "\nCould not open sp3 precise ephemeris file:"<<ephemerisFile <<" Terminating...\n";
            exit(0);
        }
        const int MAX = 100;
        char store[MAX]={0};
        long lineNumber = 0;
        int totalsat = 0;
        bool firstCC = false;
        TimeSystem ts;
        
        GString reader;
        
        while ( !infile.eof() )
        {
            //infile.getline(store, MAX); //read past the file header
            //GString reader(store);
            
            std::getline(infile, reader);
            
            lineNumber++;
            
            if( lineNumber == 3 ) //lineNumber 3 has the total number of satellite
            {
                std::vector<GString> splitstr = reader.split();
                totalsat = static_cast<int>(splitstr[1].asINT());
            }
            
            if( reader[0] == '%' && reader[1] == 'c' && firstCC == false)
            {
                firstCC = true;
                std::vector<GString> tmpstr = reader.split();
                
                ts = GTimeSystem("tsGPS");
                
                if(tmpstr[3] == "GPS")
                {
                    ts = GTimeSystem("tsGPS");
                }
                else if(tmpstr[3] == "UTC")
                {
                    ts = GTimeSystem("tsUTC");
                }
                
            }
            
            
            if( reader[0] == '*' ) // start one epoch
            {
                std::vector<GString> splitstr = reader.split();
                CivilTime ct;
                //default timesystem is GPS time
                ct.m_ts = ts; //GTimeSystem("tsGPS");
                
                ct.m_year = static_cast<int>(splitstr[1].asINT() ) ;
                ct.m_month = static_cast<int>( splitstr[2].asINT());
                ct.m_day = static_cast<int>(splitstr[3].asINT()) ;
                ct.m_hour = static_cast<int>(splitstr[4].asINT()) ;
                ct.m_minute = static_cast<int>(splitstr[5].asINT()) ;
                ct.m_second = splitstr[6].asDOUBLE();
                GTime epoch =  GTime::CivilTime2GTime(ct);
                
                
                
                // this is a very important variable concerning the format of the sp3 file
                bool hasVelocity = false;
                
                GPreciseEphemeris peph;
                
                peph.setEpoch(epoch);
                
                for( int i = 0 ; i< totalsat; i++ )
                {
                    
                    std::getline(infile, reader);
                    peph.setData(reader);
                    
                    if( hasVelocity == true)
                    {
                        std::getline(infile, reader);
                        peph.setData(reader);
                    }
                    
                    if(peph.isOK() == false)
                    {
                        continue;
                    }
                    
                    //first check whether this satellite exist
                    bool isexist = ( gfc::GSpaceCraftMgr::gSpacecraft.find(peph.getSensorID().getIDString()) != gfc::GSpaceCraftMgr::gSpacecraft.end() );
                    
                    if( isexist == true )
                    {
                        gfc::GSpaceCraftMgr::gSpacecraft[peph.getSensorID().getIDString()].pushPreciseEphemeris(peph);
                    }
                    
                }
            }
        }
        
        //after loading, should check the repeated elements, keep the precise ephemeris unique
        // travel the whole map
 //       std::map< GString,GSpaceCraft>::iterator myiter = gfc::GSpaceCraftMgr::gSpacecraft.begin();
//        int size = gfc::GSpaceCraftMgr::gSpacecraft.size();
//        for(; myiter != gfc::GSpaceCraftMgr::gSpacecraft.end(); ++myiter)
//        {
//            myiter->second.checkUniqe();
//        }
        
        
    } // end of function: loadPreciseEphemeris
    
    
    
    // load the rinex3 broadcast ephemeris
    void GSpaceCraftMgr::loadBroadcastEphemeris( GString ephemerisFile)
    {
        
        
        
    }
        
    
    
    
} // end of namespace gfc

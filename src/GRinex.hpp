//
//  GRinex.hpp
//  GFC
//
//  Created by lizhen on 24/08/2016.
//  Copyright © 2016 lizhen. All rights reserved.
//

#ifndef GRinex_hpp
#define GRinex_hpp

#include <stdio.h>
#include <vector>
#include <fstream>
#include <utility>

#include "GString.h"
#include "GTime.h"

#include "GNSSData.h"

#include "GSensor.h"
#include "GSourceID.h"

#include "GOBSDataStorage.h"

namespace gfc
{
    
    // the class to proces the rinex file
    // read and write
    class GRinex
    {
        
        struct header
        {
            struct OBSTYPE
            {
                char sys;
                int  num;
                std::vector<GString> obscode;
                
                // the gnss obstype in the software
                std::vector<GgnssOBSType> gnssobsType;
                
            };
            // the definions for sat sys in Rinex 3
            enum SATSYS{G,R,E,J,C,I,S};
            
            double version;
            char filetype; // should be O
            char satsys;  // can be G, R, C, E, J,S,I,and M
            
            TimeSystem ts;
            
            GString PGM,RUNBY,DATE,markerName,markerNum;
            GString observer,agency,receiverNumber,receiverType,receiverVersion,antennaNum,antennaType;
            double approxPos[3];
            double antDel[3];
            
            std::vector<OBSTYPE> obstypes;
            int  sysIndex[7];  // G R E J C I S
            
            double interval;
            GTime timefirstObs, timelastObs;
            int leapsec;
            header()
            {
                version =0.0;
                filetype ='O';
                satsys = 'G';
                memset(sysIndex,-1,sizeof(int)*7);
                interval = 0.0;
                leapsec = 0.0;
            }
        };
        
    public:
        
        GRinex(){};
        
        void parseOfileHeader3();
        
        bool nextOBS( GgnssStorage& datastorage);
        
        GString satsys2SensorSys(char sys);
        int     satsys2Index(char sys);
        GString satsys2timeString(char sys);
        GString  timesysString(GString str);
        GgnssOBSType rinexObsCode2GobsType(char sys, GString obscode);
        
        /*
        GRinex& operator= (const GRinex& right)  //赋值重载
        {
            //this->m_timeSystemName = right.m_timeSystemName;
            this->h = right.h;
            this->m_station = right.m_station;
            this->rinexFile = right.rinexFile;
            this->stream = right.stream;
            return *this;
        }
        
        GRinex( const GRinex& right )   //拷贝构造函数
        {
            this->h = right.h;
            this->m_station = right.m_station;
            this->rinexFile = right.rinexFile;
            this->stream = right.stream;
        }
         */
        
        void setStation( GString fileName, GSourceID mystation )
        {
            m_station = mystation;
            
            rinexFile.open(fileName);
            
            if( rinexFile.is_open() == false )
            {
                printf("rinex file is unavailable\n");
                return;
            }
            
            parseOfileHeader3();
        }
        
    private:
        
        header  h; // the header struct
        
        GSourceID m_station; // the station indicator
        
        std::fstream  rinexFile;
        std::stringstream stream;
        
    };  // end of namespace gfc
    
    
    
} // end of namespace



#endif /* GRinex_hpp */

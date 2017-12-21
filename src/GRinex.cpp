//
//  GRinex.cpp
//  GFC
//
//  Created by lizhen on 24/08/2016.
//  Copyright © 2016 lizhen. All rights reserved.
//

#include "GRinex.hpp"

namespace gfc
{
    
    GString GRinex::satsys2timeString(char sys)
    {
        GString res;
        if(sys == 'G') {res = "tsGPS";}
        else if(sys == 'R') {res= "tsGLO";}
        else if(sys == 'E') {res= "tsGAL";}
        else if(sys == 'J') {res= "tsQZS";}
        else if(sys == 'C') {res= "tsBDS";}
        else if(sys == 'I') {res= "tsIRN";}
        else if(sys ==' ')  {res= "tsUNK";}
        else if(sys == 'M') {res= "tsUKN";}
        return res;
    } // end of function satsys2timeString
    
    // transfer from RINEX time system string to GFC time system string
    GString  GRinex::timesysString(GString str)
    {
        GString res;
        if(str == "GPS") {res = "tsGPS";}
        else if(str == "GLO") {res =  "tsGLO";}
        else if(str == "GAL") {res = "tsGAL";}
        else if(str == "QZS") {res = "tsQZS";}
        else if(str == "BDT") {res = "tsBDS";}
        else if(str == "IRN") {res = "tsIRN";}
        
        return res;
    }

    //convert from
    GgnssOBSType GRinex::rinexObsCode2GobsType(char sys, gfc::GString obscode)
    {
        GgnssOBSType gobstype;
        
        //tracking code
        gobstype.m_trackingCode = obscode[2];
        
        //determine the obs type
        switch (obscode[0])
        {
            // Code/Pseudorange
            case 'C':
                gobstype.m_dataType =  GOBSType::GetPointer("ot_RANGE");
                break;
            case 'L':
                gobstype.m_dataType =  GOBSType::GetPointer("ot_PHASE");
                break;
            case 'D':
                gobstype.m_dataType =  GOBSType::GetPointer("ot_DOPPLER");
                break;
                // carrier to noise ratio
            case 'S':
                gobstype.m_dataType =  GOBSType::GetPointer("ot_RSSI");
                break;
                
            default:
                break;
        }
        
        //determine the frequency according to the sat sys
        switch ( obscode[1] )
        {
            case '1':
                if( sys =='G' ) // gps
                {
                    gobstype.m_freq = GCarrierFreq::GetPointer("freq_L1");
                }
                else if(sys=='R') // GLONASS
                {
                    gobstype.m_freq = GCarrierFreq::GetPointer("freq_G1");
                }
                else if(sys=='E') //Galileo
                {
                    gobstype.m_freq = GCarrierFreq::GetPointer("freq_E1");
                }
                else if( sys =='C' ) // BDS
                {
                    gobstype.m_freq = GCarrierFreq::GetPointer("freq_B1");
                }

                break;
            case '2':
                if( sys =='G') // GPS
                {
                    gobstype.m_freq = GCarrierFreq::GetPointer("freq_L2");
                }
                else if(sys =='J') // QZSS
                {
                   // gobstype.m_freq = GCarrierFreq::GetPointer("freq_L2");
                }
                else if(sys =='R') // GLONASS
                {
                    gobstype.m_freq = GCarrierFreq::GetPointer("freq_G2");
                }
                break;
            case '3':
                
                break;
            case '4':
                
                break;
            case '5':
                if( sys =='G') // GPS
                {
                    gobstype.m_freq = GCarrierFreq::GetPointer("freq_L5");
                }
                else if(sys =='J') // QZSS
                {
                    // gobstype.m_freq = GCarrierFreq::GetPointer("freq_L2");
                }
                else if(sys=='S') // SBAS payload
                {
                    
                }
                else if(sys =='E') // GALILEO
                {
                    gobstype.m_freq = GCarrierFreq::GetPointer("freq_E5a");
                }
                else if( sys =='I' ) // IRNSS
                {
                    //gobstype.m_freq = GCarrierFreq::GetPointer("freq_B1");
                }
                break;
            case '6':
                if( sys =='E') // Galileo
                {
                    gobstype.m_freq = GCarrierFreq::GetPointer("freq_E6");
                }
                else if(sys =='J') // QZSS, LEX
                {
                    // gobstype.m_freq = GCarrierFreq::GetPointer("freq_L2");
                }
                else if(sys=='C') // BDS
                {
                    gobstype.m_freq = GCarrierFreq::GetPointer("freq_B3");
                }
                
                break;
            case '7':
                if( sys =='E') // Galileo
                {
                    gobstype.m_freq = GCarrierFreq::GetPointer("freq_E5b");
                }
                else if(sys=='C') // BDS
                {
                    gobstype.m_freq = GCarrierFreq::GetPointer("freq_B2");
                }
                
                break;
            case '8':
                // Galileo E5a+b frequency??
                if( sys == 'E' ) //galileo E5
                {
                    gobstype.m_freq = GCarrierFreq::GetPointer("freq_E5");
                }
                
                break;
            case '9':
                //IRNSS S frequency ?
                
                break;
                
            default:
                break;
        }
        
        return gobstype;
        
    } // end of function rinexObsCode2GobsType
    
    
    
    GString GRinex::satsys2SensorSys(char sys)
    {
        GString sensorsys;
        switch (sys)
        {
            case 'C':
                sensorsys = "ssBDS";
                break;
            case 'G':
                sensorsys = "ssGPS";
                break;
            case 'J':
                sensorsys = "ssQZS";
                break;
            case 'E':
                sensorsys = "ssGAL";
                break;
            case 'R':
                sensorsys = "ssGLO";
                break;
            case 'S':
                sensorsys = "ssSBS";
                break;
            case 'I':
                sensorsys = "ssIRN";
                break;
                
            default:
                break;
        }
        return sensorsys;
    }
    
    int GRinex::satsys2Index(char sys)
    {
        int index = -1;
        switch ( sys)
        {
            case 'G':
                index = header::G;
                break;
            case 'R':
                index = header::R;
                break;
            case 'E':
                index = header::E;
                break;
            case 'J':
                index = header::J;
                break;
            case 'C':
                index = header::C;
                break;
            case 'I':
                index = header::I;
                break;
            case 'S':
                index = header::S;
                break;
                
            default:
                break;
        }
        
        return index;
        
    } // end of function satsys2Index
    
   void GRinex::parseOfileHeader3()
    {
        while( ! rinexFile.eof() )
        {
            char pstr[81] = {0};
            rinexFile.getline(pstr, 81);
            
            if(strstr(pstr,"END OF HEADER") != NULL)
            {
                break;
            }
            
            if( strstr(pstr,"RINEX VERSION / TYPE") != NULL )
            {
                //format F9.2,11X,A1,19X,A1,19X
                sscanf(pstr, "%9lf%*11c%c%*19c%c%*19c%*s",&h.version,&h.filetype,&h.satsys );
            }
            else if(strstr(pstr,"PGM / RUN BY / DATE") != NULL)
            {
                //format A20,A20,A20
                char t1[21] = {0}, t2[21]={0},t3[21]={0};
                sscanf(pstr, "%20c%20c%20c%*s",t1,t2,t3);
                h.PGM = t1;h.RUNBY = t2;h.DATE = t3;
                // h.PGM.c_str(),h.RUNBY.c_str(),h.DATE.c_str() );
                
            }
            else if( strstr(pstr,"MARKER NAME") != NULL )
            {
                char t[61]={0};
                sscanf(pstr, "%60c%*s",t);
                h.markerName = t;
            }
            else if( strstr(pstr,"MARKER NUMBER") != NULL )
            {
                char t[21]={0};
                sscanf(pstr, "%20c%*s",t);
                h.markerNum = t;
            }
            else if( strstr(pstr,"OBSERVER / AGENCY") != NULL )
            {
                char t1[21]={0},t2[41]={0};
                sscanf(pstr, "%20c%40c%*s",t1,t2);
                h.observer = t1;
                h.agency   = t2;
            }
            else if( strstr(pstr,"REC # / TYPE / VERS") != NULL )
            {
                char t1[21]={0},t2[21]={0},t3[21]={0};
                sscanf(pstr, "%20c%20c%20c%*s",t1,t2,t3);
                h.receiverNumber = t1;
                h.receiverType   = t2;
                h.receiverVersion = t3;
            }
            else if( strstr(pstr,"ANT # / TYPE") != NULL )
            {
                char t1[21]={0},t2[21]={0};
                sscanf(pstr, "%20c%20c%*s",t1,t2);
                h.antennaNum = t1;
                h.antennaType   = t2;
                
            }
            else if( strstr(pstr,"APPROX POSITION XYZ") != NULL )
            {
                //char t1[20]={0},t2[20]={0},t3[20]={0};
                sscanf(pstr, "%14lf%14lf%14lf%*s",h.approxPos,h.approxPos+1,h.approxPos+2);
                int testc =0;
            }
            else if( strstr(pstr,"ANTENNA: DELTA H/E/N") != NULL )
            {
                //char t1[20]={0},t2[20]={0},t3[20]={0};
                sscanf(pstr, "%14lf%14lf%14lf%*s",h.antDel,h.antDel+1,h.antDel+2);
                int testc =0;
            }
            else if( strstr(pstr,"SYS / # / OBS TYPES") != NULL )
            {
                if( pstr[0] !=' ' )
                {
                    header::OBSTYPE ot;
                    ot.sys = pstr[0];
                    GString t(pstr);
                    ot.num  = t.substr(3,3).asINT();
                    int m1 = ot.num % 13;
                    int m2 = int(ot.num / 13);
                    if( m2 == 0 || (m2==1 && m1==0)) // less than 13
                    {
                        ot.obscode =  t.substr(6,54).split();
                    }
                    else // greater than 13
                    {
                        ot.obscode =  t.substr(6,54).split();
                        for( int i = 0 ; i< m2; i++ )
                        {
                            rinexFile.getline(pstr, 81);
                            GString t(pstr);
                            std::vector<GString> tt = t.substr(6,54).split();
                            for( int j = 0 ; j< tt.size(); j++)
                            {
                                ot.obscode.push_back(tt[j]);
                            }
                        }
                    }
                    
                    h.obstypes.push_back(ot);
                }
            }
            
            else if(strstr(pstr,"INTERVAL") != NULL)
            {
                sscanf(pstr, "%10lf%*s",&h.interval);
                int testc =0;
            }
            else if(strstr(pstr,"LEAP SECONDS") != NULL)
            {
                sscanf(pstr, "%6d%*6d%*6d%*6d%*3c%*s",&h.leapsec);
            }
            else if(strstr(pstr,"TIME OF FIRST OBS") != NULL)
            {
                CivilTime ct;
                char ts[4]={0};
                sscanf(pstr, "%6d%6d%6d%6d%6d%13lf%*5c%3c",&ct.m_year,&ct.m_month,&ct.m_day,
                                                           &ct.m_hour,&ct.m_minute,&ct.m_second,ts);
                ct.m_ts = timesysString(GString(ts));
                
                h.timefirstObs = GTime::CivilTime2GTime(ct);
                
            }
            else if(strstr(pstr,"TIME OF LAST OBS") != NULL)
            {
                CivilTime ct;
                char ts[4]={0};
                sscanf(pstr, "%6d%6d%6d%6d%6d%13lf%*5c%3c",&ct.m_year,&ct.m_month,&ct.m_day,
                       &ct.m_hour,&ct.m_minute,&ct.m_second,ts);
                ct.m_ts = timesysString(GString(ts));
                
                h.timelastObs = GTime::CivilTime2GTime(ct);
            }
        }
        
        //determine the time system in the record
        TimeSystem myts0("tsUKN");
        TimeSystem myts = satsys2timeString(h.satsys);
        if( myts == myts0 && h.timefirstObs.getTimeSystem() == myts0 &&h.timelastObs.getTimeSystem()==myts0)
        {
            printf("WARNING: no time system found! using the default GPS time!\n");
            h.ts = TimeSystem("tsGPS");
        }
        if( myts != myts0 )
        {
            h.ts = myts;
        }
        
        if(h.timefirstObs.getTimeSystem() != myts0)
        {
            h.ts = h.timefirstObs.getTimeSystem();
        }
        if(h.timelastObs.getTimeSystem() != myts0)
        {
            h.ts = h.timelastObs.getTimeSystem();
        }
        
        
        for(int i = 0 ; i< h.obstypes.size(); i++ )
        {
            
            h.sysIndex[satsys2Index(h.obstypes[i].sys)] = i;
            
            for( int j = 0 ; j< h.obstypes[i].obscode.size(); j++)
            {
                GgnssOBSType gtype = rinexObsCode2GobsType(h.obstypes[i].sys,
                                                           h.obstypes[i].obscode[j]);
                h.obstypes[i].gnssobsType.push_back(gtype);
                
            }
            
        }
        
        
        int testc =0;
    }
    
    
    bool GRinex::nextOBS(GgnssStorage& datastorage)
    {
        GString line;
        GSensorID sat;
        
        bool available = false;
        while( ! rinexFile.eof() )
        {
            if( available == true )
            {
                break;
            }
            
            getline( rinexFile, line ); // epoch time line
            
            if( line[0] == '>' && ( line[31] == '0' || line[31] == '1' ) ) // available data epoch
            {
                available = true;
                //one epoch data
                GgnssDataEpoch myepoch;
                
                stream.clear();
                stream.str(line);
                CivilTime ct;
                //stream.width(4);
                char c;
                int satnum = 0;
                stream >> c >> ct.m_year >> ct.m_month >> ct.m_day >> ct.m_hour
                >> ct.m_minute >> ct.m_second >>c >> satnum;
                ct.m_ts = h.ts;
                
                myepoch.m_epoch = GTime::CivilTime2GTime(ct);
                myepoch.m_obsdatalist.resize(satnum);
                myepoch.m_obstypelist.resize(satnum);
                for( int i = 0 ; i< satnum; i++ )
                {
                    getline( rinexFile, line );
                    
                    stream.clear();
                    
                    stream.str(line);
                    
                    int index = h.sysIndex[satsys2Index(line[0])];
                    char sys;
                    int prn = 0;
                    stream >> std::setw(1) >> sys >> std::setw(2) >> prn;
                    GString sensorsys = satsys2SensorSys(sys);
                    
                    myepoch.m_indicatorlist.push_back(GSensorID(sensorsys,prn) );
                    myepoch.m_obstypelist[i] = h.obstypes[index].gnssobsType;
                    myepoch.m_obsdatalist[i].resize(h.obstypes[index].num);
                    
                    int tag = 0;
                    for( int j = 0 ;j< h.obstypes[index].num; j++)
                    {
                        double value =0;
                        unsigned char lli =-1, ssi =-1;
                        
                        if(tag+3 >= line.length())
                        {
                            myepoch.m_obsdatalist[i][j].m_data_value = 0.0;
                            myepoch.m_obsdatalist[i][j].m_LLI = ' ';
                            myepoch.m_obsdatalist[i][j].m_SSI = ' ';
                            continue;
                        }
                        
                        myepoch.m_obsdatalist[i][j].m_data_value = line.substr(tag+3,14).asDOUBLE();
                        
                        if( tag + 17 >= line.length())
                        {
                            myepoch.m_obsdatalist[i][j].m_LLI = ' ';
                            myepoch.m_obsdatalist[i][j].m_SSI = ' ';
                        }
                        else
                        {
                            myepoch.m_obsdatalist[i][j].m_LLI = line[tag+17];
                            myepoch.m_obsdatalist[i][j].m_SSI = line[tag+18];
                        }
                        
                        tag += 16;
                    }
                    
                    // the satellite list
                    
                }
                
                // update return value
                
                datastorage[m_station] = myepoch;
                
                
            }
            
        }
        
        return available;
        
    } // end of function nextOBS
    
    
    
    
}  // end of namespace

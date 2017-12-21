//
//  GSLRcrd.hpp
//  GFC
//
//  Created by lizhen on 30/08/2016.
//  Copyright © 2016 lizhen. All rights reserved.
//

#ifndef GSLRcrd_hpp
#define GSLRcrd_hpp

#include <stdio.h>

#include <fstream>

#include <map>

#include "GString.h"

#include "GSourceID.h"

#include "GOBSData.h"

#include "GTime.h"

namespace gfc
{
    
    /*  Consolidated laser Ranging Data format (CRD)
     *       record and variable definitions for FORTRAN
     *       R. Ricklefs UT/CSR July 2007
     *  History:
     *  08/xx/07 - added H3 Target type (v0.27)
     *  11/26/07 - added H4 data quality alert
     *	       and #10 stop number
     *	       and #20 origin of values (v0.27) rlr.
     *  05/07/08 - Expand all configuration and data section character strings
     *             to allow up to 40 characters (plus NULL).
     *             - Added detector channel to normalpoint (11) and calibration (40)
     *             records.
     *             - Added field for 'crd' literal to H1.
     *             - Record '21' sky_clarity is not double rather than int.
     *               (v1.00 rlr)
     *
     * ----------------------------------------------------------------------
     */
    
    /* Common name abbreviations:
     *  app = applied
     *  cdp = Crustal Dynamics Project (old NASA name)
     *  CofM = Center of Mass
     *  corr = correction or corrected
     *  est = estimated
     *  ind = indicator
     *  num = number
     *  osc = oscillator
     *  off = offset
     *  PmM = peak minus mean
     *  sic = (Goddard) Satellite ID Code
     *  stn = station
     *  SC = spacecraft
     *  sys = system
     *  utc = Universal Time Coordinated
     *  xcv = receive
     */
    
    /* Ranging data header fields */
    /* H1 - format header */
    struct rh1 {
        char crd_literal[4];
        int format_version;
        int prod_year;
        int prod_mon;
        int prod_day;
        int prod_hour;
    };
    
    /* H2 - station header */
    struct rh2 {
        char stn_name[11];
        int cdp_pad_id;
        int cdp_sys_num;
        int cdp_occ_num;
        int stn_timescale;
    };
    
    /* H3 - spacecraft header */
    struct rh3 {
        char target_name[11];
        int ilrs_id;
        int sic;
        int norad;
        int SC_timescale;
        int target_type;
    };
    
    /* H4 - Session header */
    struct rh4 {
        int data_type;
        int start_year;
        int start_mon;
        int start_day;
        int start_hour;
        int start_min;
        int start_sec;
        int end_year;
        int end_mon;
        int end_day;
        int end_hour;
        int end_min;
        int end_sec;
        int data_release;
        int refraction_app_ind;
        int CofM_app_ind;
        int xcv_amp_app_ind;
        int stn_sysdelay_app_ind;
        int SC_sysdelay_app_ind;
        int range_type_ind;
        int data_qual_alert_ind;
    };
    
    /* Need indicators that these have been read? */
    /* H8 - End of Session footer */
    /* H9 - End of File footer */
    
    /* Ranging data configuration fields (1 of n) */
    /* C0 - System Configuration Record */
    struct rc0 {
        int detail_type;
        double xmit_wavelength;
        /**
         char sysconfig_id[5];
         char laserconfig_id[4];
         char detectorconfig_id[4];
         char timingconfig_id[4];
         char xponderconfig_id[4];
         **/
        char config_ids[10][41];
    };
    
    /* C1 - Laser Configuration Record */
    struct rc1 {
        int detail_type;
        char laser_config_id[41];
        char laser_type[41];
        double prim_wavelength;	/* Primary wavelength of laser */
        double nom_fire_rate;	/* Nominal fire rate of laser */
        double pulse_energy;
        double pulse_width;
        double beam_div;
        int pulses_in_semitrain;	/* for multi-pulse systems */
    };
    
    /* C2 - Detector Configuration Record */
    struct rc2 {
        int detail_type;
        char detector_config_id[41];
        char detector_type[41];
        double app_wavelength;
        double qe;			/* quantum efficiency (in %) */
        double voltage;
        double dark_count;
        char output_pulse_type[41];
        double output_pulse_width;
        double spectral_filter;
        double spectral_filter_xmission;	/* % transmission of filter */
        double spatial_filter;
        char signal_proc[41];	/* signal processing algorithm or pgm name */
    };
    
    /* C3 - Timing Configuration Record */
    struct rc3 {
        int detail_type;
        char timing_config_id[41];
        char time_source[41];
        char freq_source[41];
        char timer[41];
        char timer_serial_num[41];
        double epoch_delay_corr;
    };
    
    /* C4 - Transponder Configuration Record */
    struct rc4 {
        int detail_type;
        char xponder_config_id[41];
        long double est_stn_utc_offset;
        double est_stn_osc_drift;
        long double est_xponder_utc_offset;
        double est_xponder_osc_drift;
        long double xponder_clock_ref_time;
        int stn_off_drift_app_ind;
        int SC_off_drift_app_ind;
        int SC_time_simplified_ind;
    };
    
    /* Ranging data fields */
    /* Secofday: need int sec and int psec? */
    /* 10 - Range Record */
    struct rd10 {
        long double sec_of_day;
        long double time_of_flight;
        char sysconfig_id[41];
        int epoch_event;
        int filter_flag;
        int detector_channel;
        int stop_number;
        int xcv_amp;
    };
    
    /* 11 - Normal Point Record */
    struct rd11 {
        long double sec_of_day;
        long double time_of_flight;
        char sysconfig_id[41];
        int epoch_event;
        double np_window_length;
        int num_ranges;
        double bin_rms;
        double bin_skew;
        double bin_kurtosis;
        double bin_PmM;
        double return_rate;
        int detector_channel;
    };
    
    /* 12 - Range Supplement Record */
    struct rd12 {
        long double sec_of_day;
        char sysconfig_id[41];
        double refraction_corr;
        double target_CofM_corr;
        double nd_value;
        double time_bias;
    };
    
    /* 20 - Meteorological Record */
    struct rd20 {
        long double sec_of_day;
        double pressure;
        double temperature;
        double humidity;
        int value_origin;
    };
    
    /* 21 - Meteorological Supplement Record */
    struct rd21 {
        long double sec_of_day;
        double wind_speed;
        double wind_direction;
        char precip_type[41];
        int visibility;
        double sky_clarity;
        int atmospheric_seeing;
        int cloud_cover;
    };
    
    /* 30 - Pointing Angles Record */
    struct rd30 {
        long double sec_of_day;
        double azimuth;
        double elevation;
        int direction_ind;
        int angle_origin_ind;
        int refraction_corr_ind;
    };
    
    /* 40 - Calibration Record */
    struct rd40 {
        long double sec_of_day;
        int type_of_data;
        char sysconfig_id[41];
        int num_points_recorded;
        int num_points_used;
        double one_way_target_dist;
        double cal_sys_delay;
        double cal_delay_shift;
        double cal_rms;
        double cal_skew;
        double cal_kurtosis;
        double cal_PmM;
        int cal_type_ind;
        int cal_shift_type_ind;
        int detector_channel;
        
        rd40()
        {
            sec_of_day =0.0;
            type_of_data = -1;
            num_points_recorded =0;
            num_points_used =0;
            one_way_target_dist =0.0;
            cal_sys_delay =0.0;
            cal_delay_shift =0.0;
            cal_rms =0.0;
            cal_skew =0.0;
            cal_kurtosis =0.0;
            cal_PmM =0.0;
            cal_type_ind =0;
            cal_shift_type_ind =0;
            detector_channel =0;
        }
        
    };
    
    /* 50 - Session Statistics Record */
    struct rd50 {
        char sysconfig_id[41];
        double sess_rms;
        double sess_skew;
        double sess_kurtosis;
        double sess_PmM;
        int data_qual_ind;
    };
    
    /* 60 - Compatibility Record */
    struct rd60 {
        char sysconfig_id[41];
        int sys_change_ind;
        int sys_config_ind;
    };
    
    /* 9X - User Defined Record */
    struct rd9x {
        /**********
         Add userdefined record types and fields here
         **********/
    };
    
    /* 00 - Comment Record */
    struct rd00 {
        char comment[81];
    };
    
    
    struct GslrData : public GOBSData
    {
        GTime  obsepoch;
        int epoch_event;
        
        GslrData()
        {
            epoch_event = -1;
        }
        
        //call the constructor of the father class
        GslrData( double dataValue, double dataStdDev)
        {
            m_data_value = dataValue;
            m_std_dev = dataStdDev;
            
        }
    };
    
    
    class GSLRBlock : public GSourceID
    {
        
    public:
        
        
        //constructor
        GSLRBlock()
        {
            
        }
        
        /// Copy constructor.
        GSLRBlock( const GSLRBlock& g )
        {
            h1 = g.h1;h2 = g.h2;h3 = g.h3;h4 = g.h4;
            c0 = g.c0;c1 = g.c1;c2 = g.c2;c3 = g.c3; c4 = g.c4;
            
            //m_sourceName = g.m_sourceName;
            
            setSourceName(g.getSourceName());
            
            MeteorRecord = g.MeteorRecord;
            
            MeteorSupplement = g.MeteorSupplement;
            
            PointingAnglesRecord = g.PointingAnglesRecord;
            CalibrationRecord = g.CalibrationRecord;
            
            obstype = g.obstype;
            
            obsdata = g.obsdata;
        }
        
        //according to second of day , get the temperature and pressure
        void getTP( double secofday , double& temperature, double& pressure)
        {
            double t = 0.0, p =0.0;
            
            if( MeteorRecord.size() == 0 )
            {
               // call the Global Pressure and Temperature Model
                printf("**warning**: NO TEMPERATURE and PRESSURE data !\n");
            }
            else
            {
                for( int i = 0 ; i< MeteorRecord.size(); i++ )
                {
                    t += MeteorRecord[i].temperature;
                    p += MeteorRecord[i].pressure;
                }
                
                temperature = t / MeteorRecord.size();
                pressure = p / MeteorRecord.size();
            }
            
        }
        
        //data members
        
        rh1 h1;
        rh2 h2;
        rh3 h3;
        rh4 h4;
        rc0 c0;
        rc1 c1;
        rc2 c2;
        rc3 c3;
        rc4 c4;
        
        std::vector<rd20>  MeteorRecord;
        std::vector<rd21>  MeteorSupplement;
        
        std::vector<rd30>  PointingAnglesRecord ;
        std::vector<rd40>  CalibrationRecord ;
        
        //obstype must be otRange
        GOBSType obstype;
        
        std::vector<GslrData>  obsdata;
        
        
    };  // the class GSLRStation
    
    
    
    class GSLRStation
    {
        
    public:
        
        static void loadStationCoordinate( GString sinexFileName );
        static void loadStationEccentricity(GString sinexFileName);
        
        GTime  epoch;
        
        int code; // the slr code of this station
        // position of the station, unit: meter
        double staP[3];
        // velocity of the station, unit: meter/year
        double staV[3];
        // std of the position , unit: meter
        double stdP[3];
        // std of the velocity, unit: meter
        double stdV[3];
        
        // need to find this information for all stations
        double neu[3];  // the neu offset
        
        
        GSLRStation()
        {
            code = -1;
            memset(staP,0,sizeof(double)*3);
            memset(staV,0,sizeof(double)*3);
            memset(stdP,0,sizeof(double)*3);
            memset(stdV,0,sizeof(double)*3);
            memset(neu,0,sizeof(double)*3);
            
        }
        
        static std::map<int,GSLRStation>  stationInfo;
        
    private:
        
        
    };
    
    
    //class GgnssStorage : public std::map<GSourceID, GgnssDataEpoch >
    
   // typedef std::map<GSensorID,GSLRStation> GslrStorage ;
    
    class GslrStorage : public std::map<GSensorID,std::map<GSourceID,GSLRBlock> >
    {
        
    public:
        
        void push_data(GSLRBlock slrdata)
        {
            GSensorID sensorid = GSensorID::ilrsName2GSensorID(slrdata.h3.target_name);
            GSourceID sourceid(slrdata.h2.stn_name);
            
            slrdata.setSourceName(slrdata.h2.stn_name);
            
            std::map<GSensorID,std::map<GSourceID,GSLRBlock> >::iterator it;
            it = find(sensorid);
            if( it == end() )  // not exist
            {
                std::map<GSourceID,GSLRBlock> m;
                m[sourceid] = slrdata;
                (*this)[sensorid] = m;
            }
            else  // sensor exist
            {
                //check if sourceid exist
                std::map<GSourceID,GSLRBlock>::iterator myit =  it->second.find(sourceid);
                if( myit == it->second.end() )
                {
                    it->second[sourceid] = slrdata;
                }
                else  //exist
                {
                    //obsdata combination
                    //myit->second.obsdata = myit->second.obsdata +  slrdata.obsdata;
                    
                    myit->second.obsdata.insert(myit->second.obsdata.end(), slrdata.obsdata.begin(),slrdata.obsdata.end()) ;
                    
                }
            }
        }
        
        
    };
    
    
    // this class is to process the SLR crd file from ftp
    // ftp://edc.dgfi.tum.de/pub/slr/data/npt_crd/compassi3/2014
    //compassi3_20140330.npt
    class GSLRcrd
    {
        
        
        
        //reading function
        int read_h1 (char * str, struct rh1 *header);
        int read_h2 (char * str, struct rh2 *header);
        int read_h3 (char * str, struct rh3 *header);
        int read_h4 (char * str, struct rh4 *header);
        void read_h8 (char * str);
        void read_h9 (char * str);
        int  read_c0 (char * str, struct rc0 *config);
        int  read_c1 (char * str, struct rc1 *config);
        int  read_c2 (char * str, struct rc2 *config);
        int read_c3 (char * str, struct rc3 *config);
        int read_c4 (char * str, struct rc4 *config);
        int read_10 (char * str, struct rd10 *data_recd);
        int read_11 (char * str, struct rd11 *data_recd);
        int read_12 (char * str, struct rd12 *data_recd);
        int read_20 (char * str, struct rd20 *data_recd);
        int read_21 (char * str, struct rd21 *data_recd);
        int read_30 (char * str, struct rd30 *data_recd);
        int read_40 (char * str, struct rd40 *data_recd);
        int read_50 (char * str, struct rd50 *data_recd);
        int read_60 (char * str, struct rd60 *data_recd);
        int read_9x (char * str, struct rd9x *data_recd);
        int read_00 (char * str, struct rd00 *data_recd);
        
        //writing function
        void write_h1 (FILE * str_out, struct rh1 header);
        void write_h2 (FILE * str_out, struct rh2 header);
        void write_h3 (FILE * str_out, struct rh3 header);
        void write_h4 (FILE * str_out, struct rh4 header);
        void write_h8 (FILE * str_out);
        void write_h9 (FILE * str_out);
        void write_c0 (FILE * str_out, struct rc0 config);
        void write_c1 (FILE * str_out, struct rc1 config);
        void write_c2 (FILE * str_out, struct rc2 config);
        void write_c3 (FILE * str_out, struct rc3 config);
        void write_c4 (FILE * str_out, struct rc4 config);
        void write_10 (FILE * str_out, struct rd10 data_recd);
        void write_11 (FILE * str_out, struct rd11 data_recd);
        void write_12 (FILE * str_out, struct rd12 data_recd);
        void write_20 (FILE * str_out, struct rd20 data_recd);
        void write_21 (FILE * str_out, struct rd21 data_recd);
        void write_30 (FILE * str_out, struct rd30 data_recd);
        void write_40 (FILE * str_out, struct rd40 data_recd);
        void write_50 (FILE * str_out, struct rd50 data_recd);
        void write_60 (FILE * str_out, struct rd60 data_recd);
        void write_9x (FILE * str_out, struct rd9x data_recd);
        void write_00 (FILE * str_out, struct rd00 data_recd);
        
        
        
        char stro[256];
        
        std::fstream  slrcrdFile;
        
    public:
        
       void setdatabase(GString filename);
       void closedatbase();
        
       bool nextDataBlock(GslrStorage& slrstorage);
        
        
    };
    
    
    
} // end of namespace gfc


#endif /* GSLRcrd_hpp */

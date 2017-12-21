//
//  GSLRcrd.cpp
//  GFC
//
//  Created by lizhen on 30/08/2016.
//  Copyright © 2016 lizhen. All rights reserved.
//

#include "GSLRcrd.hpp"

namespace gfc
{
    
    std::map<int,GSLRStation>  GSLRStation::stationInfo;
    
    void GSLRStation::loadStationEccentricity(GString sinexFileName)
    {
        std::fstream  sinexFile;
        sinexFile.open(sinexFileName);
        
        if( sinexFile.is_open() == false )
        {
            printf("slr coordinate sinex file is unavailable\n");
            return;
        }
        
        bool start = false;
        GString buff;
        while( ! sinexFile.eof()  )
        {
            std::getline(sinexFile, buff);
            
            if( buff.find("-SITE/ECCENTRICITY") != -1 )
            {
                break;
            }
            
            //find the solution/estimate
            if(  buff.find("+SITE/ECCENTRICITY") != -1 )
            {
                std::getline(sinexFile, buff);
                std::getline(sinexFile, buff);
                start = true;
            }

            if( start == true)
            {
                std::vector<GString> stringVec = buff.split();
                int stationcode = stringVec[0].asINT();
                
                std::map<int, GSLRStation>::iterator myit = stationInfo.find(stationcode);
                
                if( myit != stationInfo.end())
                {
                    
                    // N E U
                    myit->second.neu[0] = stringVec[8].asDOUBLE() ;
                    myit->second.neu[1] = stringVec[9].asDOUBLE() ;
                    myit->second.neu[2] = stringVec[7].asDOUBLE() ;
                }
            }
            
        }
        
        
        sinexFile.close();
        
        
    }
    
    
    void GSLRStation::loadStationCoordinate( GString sinexFileName )
    {
        GSLRStation station;
        
        std::fstream  sinexFile;
        sinexFile.open(sinexFileName);
        
        if( sinexFile.is_open() == false )
        {
            printf("slr coordinate sinex file is unavailable\n");
            return;
        }
        
        GString buff;
        int index =-1, code = -1,soln =-1,s;
        char pt;
        GString ref_epoch,unit,type;
        bool start = false;
        int  test = 0;
        int num =0;
        while( ! sinexFile.eof()  )
        {
            std::getline(sinexFile, buff);
            
            if( buff.find("-SOLUTION/ESTIMATE") != -1 )
            {
                break;
            }

            //find the solution/estimate
            if(  buff.find("+SOLUTION/ESTIMATE") != -1 )
            {
                 std::getline(sinexFile, buff);
                 std::getline(sinexFile, buff);
                 start = true;
            }
            
            if(start == true)
            {
                
                std::vector<GString> stringVec = buff.split();
                
                station.code = stringVec[2].asINT();
                if(stringVec[1] == "STAX")
                {
                    station.staP[0] = stringVec[8].asDOUBLE();
                    station.stdP[0] = stringVec[9].asDOUBLE();
                    test++;
                }
                else if(stringVec[1] == "STAY")
                {
                    station.staP[1] = stringVec[8].asDOUBLE();
                    station.stdP[1] = stringVec[9].asDOUBLE();
                    test++;
                }
                else if(stringVec[1] == "STAZ")
                {
                    station.staP[2] = stringVec[8].asDOUBLE();
                    station.stdP[2] = stringVec[9].asDOUBLE();
                    test++;
                }
                else if(stringVec[1] == "VELX")
                {
                    station.staV[0] = stringVec[8].asDOUBLE();
                    station.stdV[0] = stringVec[9].asDOUBLE();
                    test++;
                }
                else if(stringVec[1] == "VELY")
                {
                    station.staV[1] = stringVec[8].asDOUBLE();
                    station.stdV[1] = stringVec[9].asDOUBLE();
                    test++;
                }
                else if(stringVec[1] == "VELZ")
                {
                    station.staV[2] = stringVec[8].asDOUBLE();
                    station.stdV[2] = stringVec[9].asDOUBLE();
                    test++;
                }
                
                if(test == 6)
                {
                   
                    
//                    if( station.code == 7090 )  //YARL station
//                    {
//                        station.neu[0] = -0.0096;  // North
//                        station.neu[1] = 0.0192;  // East
//                        station.neu[2] = 3.1813;  // Up
//                    }
                    
                    stationInfo[station.code] = station;
                    station.neu[0] = 0.0;
                    station.neu[1] = 0.0;
                    station.neu[2] = 0.0;
                    test =0;
                   
                    //printf("station Index: %d\n",num++);
                }
            }
            
            
        }
        
        sinexFile.close();
        
    }
    
    /* Ranging data header/footer records */
    /* H1 - format header */
    int GSLRcrd::read_h1 (char * str, struct rh1 *header)
    {
        char temp_crd_literal[4];
        int nstat;
        
        nstat= sscanf (str,
                       "%*s %s %d %d %d %d %d",
                       temp_crd_literal,
                       &header->format_version, &header->prod_year, &header->prod_mon,
                       &header->prod_day, &header->prod_hour);
        strncpy (header->crd_literal, temp_crd_literal, 3);	/* to preserve spaces */
        header->crd_literal[3]= '\0';
        return (nstat+1);
    }
    
    /* H2 - station header */
    int GSLRcrd::read_h2 (char * str, struct rh2 *header)
    {
        int nstat;
        
        nstat= sscanf (&str[14],
                       "%d %d %d %d",
                       &header->cdp_pad_id, &header->cdp_sys_num, &header->cdp_occ_num,
                       &header->stn_timescale);
        strncpy (header->stn_name, &str[3], 10);	/* to preserve spaces */
        header->stn_name[10]= '\0';
        /*  sscanf (str, "%s\n", stro);*/
        return (nstat+2);
    }
    
    /* H3 - spacecraft header */
    int GSLRcrd::read_h3 (char * str, struct rh3 *header)
    {
        int nstat;
        
        nstat= sscanf (&str[14],
                       "%d %d %d %d %d",
                       &header->ilrs_id, &header->sic,
                       &header->norad, &header->SC_timescale,&header->target_type);
        strncpy (header->target_name, &str[3], 10);	/* to preserve spaces */
        header->target_name[10]= '\0';
        return (nstat+2);
    }
    
    /* H4 - Session header */
    int GSLRcrd::read_h4 (char * str, struct rh4 *header)
    {
        int nstat;
        
        nstat= sscanf (str,
                       "%*s %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d",
                       &header->data_type,
                       &header->start_year, &header->start_mon, &header->start_day,
                       &header->start_hour, &header->start_min, &header->start_sec,
                       &header->end_year, &header->end_mon, &header->end_day,
                       &header->end_hour, &header->end_min, &header->end_sec,
                       &header->data_release, &header->refraction_app_ind,
                       &header->CofM_app_ind, &header->xcv_amp_app_ind,
                       &header->stn_sysdelay_app_ind,
                       &header->SC_sysdelay_app_ind,
                       &header->range_type_ind, &header->data_qual_alert_ind);
        return (nstat+1);
    }
    
    /* Need indicators that these have been read? */
    /* H8 - End of Session footer */
    void GSLRcrd::read_h8 (char * str)
    {
        sscanf (str, "H8");
    }
    
    /* H9 - End of File footer */
    void GSLRcrd::read_h9 (char * str)
    {
        sscanf (str, "H9");
    }
    
    /* Ranging data configuration records (1 of n) */
    /* C0 - System Configuration Record */
    int GSLRcrd::read_c0 (char * str, struct rc0 *config)
    {
        char temp_config_ids[10][256];
        int i;
        int nstat;
        
        for (i=0; i<10; i++)
        {
            temp_config_ids[i][0]= '\0';
        }
        
        nstat= sscanf (str,
                       "%*s %d %lf %s %s %s %s %s %s %s %s %s %s",
                       &config->detail_type, &config->xmit_wavelength,
                       &temp_config_ids[0], &temp_config_ids[1], &temp_config_ids[2],
                       &temp_config_ids[3], &temp_config_ids[4], &temp_config_ids[5],
                       &temp_config_ids[6], &temp_config_ids[7], &temp_config_ids[8],
                       &temp_config_ids[9]);
        for (i=0; i<10; i++)
        {
            strncpy (config->config_ids[i], temp_config_ids[i], 40);
            config->config_ids[i][40]= '\0';
        }
        
        return (nstat+1);
    }
    
    /* C1 - Laser Configuration Record */
    int GSLRcrd::read_c1 (char * str, struct rc1 *config)
    {
        char temp_laser_config_id[256];
        int nstat;
        
        nstat= sscanf (str,
                       "%*s %1d %s %s %lf %lf %lf %lf %lf %d",
                       &config->detail_type, &temp_laser_config_id, &config->laser_type,
                       &config->prim_wavelength, &config->nom_fire_rate, &config->pulse_energy,
                       &config->pulse_width, &config->beam_div, &config->pulses_in_semitrain);
        strncpy (config->laser_config_id, temp_laser_config_id, 40);
        config->laser_config_id[40]= '\0';
        return (nstat+1);
    }
    
    /* C2 - Detector Configuration Record */
    int GSLRcrd::read_c2 (char * str, struct rc2 *config)
    {
        char temp_detector_config_id[256];
        char temp_detector_type[256];
        char temp_output_pulse_type[256];
        char temp_signal_proc[256];
        int nstat;
        
        nstat= sscanf (str,
                       "%*s %d %s %s %lf %lf %lf %lf %s %lf %lf %lf %lf %s",
                       &config->detail_type, &temp_detector_config_id,
                       &temp_detector_type, &config->app_wavelength, &config->qe,
                       &config->voltage, &config->dark_count,
                       &temp_output_pulse_type, &config->output_pulse_width,
                       &config->spectral_filter, &config->spectral_filter_xmission,
                       &config->spatial_filter, &temp_signal_proc);
        strncpy (config->detector_config_id, temp_detector_config_id, 40);
        config->detector_config_id[40]= '\0';
        strncpy (config->detector_type, temp_detector_type, 40);
        config->detector_type[40]= '\0';
        strncpy (config->output_pulse_type, temp_output_pulse_type, 40);
        config->output_pulse_type[40]= '\0';
        strncpy (config->signal_proc, temp_signal_proc, 40);
        config->signal_proc[40]= '\0';
        return (nstat+1);
    }
    
    /* C3 - Timing Configuration Record */
    int GSLRcrd::read_c3 (char * str, struct rc3 *config)
    {
        char temp_str[5][256];
        char temp_timing_config_id[256];
        char temp_time_source[256];
        char temp_freq_source[256];
        char temp_timer[256];
        char temp_timer_serial_num[256];
        int nstat;
        
        nstat= sscanf (str,
                       "%*s %1d %s %s %s %s %s %lf",
                       &config->detail_type, &temp_timing_config_id,
                       &temp_time_source, &temp_freq_source, &temp_timer,
                       &temp_timer_serial_num, &config->epoch_delay_corr);
        strncpy (config->timing_config_id, temp_timing_config_id, 40);
        config->timing_config_id[40]= '\0';
        strncpy (config->time_source, temp_time_source, 40);
        config->time_source[40]= '\0';
        strncpy (config->freq_source, temp_freq_source, 40);
        config->freq_source[40]= '\0';
        strncpy (config->timer, temp_timer, 40);
        config->timer[40]= '\0';
        strncpy (config->timer_serial_num, temp_timer_serial_num, 40);
        config->timer_serial_num[40]= '\0';
        return (nstat+1);
    }
    
    /* C4 - Transponder Configuration Record */
    int GSLRcrd::read_c4 (char * str, struct rc4 *config)
    {
        char temp_xponder_config_id[256];
        int nstat;
        
        nstat= sscanf (str,
                       "%*s %d %s %Lf %lf %Lf %lf %Lf %d %d %d",
                       &config->detail_type, &config->xponder_config_id,
                       &config->est_stn_utc_offset, &config->est_stn_osc_drift,
                       &config->est_xponder_utc_offset, &config->est_xponder_osc_drift,
                       &config->xponder_clock_ref_time,
                       &config->stn_off_drift_app_ind,
                       &config->SC_off_drift_app_ind,
                       &config->SC_time_simplified_ind);
        strncpy (config->xponder_config_id, temp_xponder_config_id, 40);
        config->xponder_config_id[40]= '\0';
        return (nstat+1);
    }
    
    /* Ranging data records */
    /* Secofday: need int sec and int psec? */
    /* 10 - Range Record */
    int GSLRcrd::read_10 (char * str, struct rd10 *data_recd)
    {
        char temp_sysconfig_id[256];
        int nstat;
        
        nstat= sscanf (str,
                       "%*s %Lf %Lf %s %d %d %d %d %d",
                       &data_recd->sec_of_day, &data_recd->time_of_flight,
                       &temp_sysconfig_id, &data_recd->epoch_event,
                       &data_recd->filter_flag, &data_recd->detector_channel,
                       &data_recd->stop_number, &data_recd->xcv_amp);
        strncpy (data_recd->sysconfig_id, temp_sysconfig_id, 40);
        data_recd->sysconfig_id[40]= '\0';
        return (nstat+1);
    }
    
    /* 11 - Normal Point Record */
    int GSLRcrd::read_11 (char * str, struct rd11 *data_recd)
    {
        char temp_sysconfig_id[256];
        int nstat;
        
        nstat= sscanf (str,
                       "%*s %Lf %Lf %s %d %lf %d %lf %lf %lf %lf %lf %d",
                       &data_recd->sec_of_day, &data_recd->time_of_flight,
                       &temp_sysconfig_id, &data_recd->epoch_event,
                       &data_recd->np_window_length, &data_recd->num_ranges,
                       &data_recd->bin_rms, &data_recd->bin_skew, &data_recd->bin_kurtosis,
                       &data_recd->bin_PmM, &data_recd->return_rate,
                       &data_recd->detector_channel);
        strncpy (data_recd->sysconfig_id, temp_sysconfig_id, 40);
        data_recd->sysconfig_id[40]= '\0';
        return (nstat+1);
    }
    
    /* 12 - Range Supplement Record */
    int GSLRcrd::read_12 (char * str, struct rd12 *data_recd)
    {
        char temp_sysconfig_id[256];
        int nstat;
        
        nstat= sscanf (str,
                       "%*s %Lf %s %lf %lf %lf %lf",
                       &data_recd->sec_of_day, &temp_sysconfig_id,
                       &data_recd->refraction_corr, &data_recd->target_CofM_corr,
                       &data_recd->nd_value, &data_recd->time_bias);
        strncpy (data_recd->sysconfig_id, temp_sysconfig_id, 40);
        data_recd->sysconfig_id[40]= '\0';
        return (nstat+1);
    }
    
    /* 20 - Meteorological Record */
    int GSLRcrd::read_20 (char * str, struct rd20 *data_recd)
    {
        int nstat;
        
        nstat= sscanf (str,
                       "%*s %Lf %lf %lf %lf %d",
                       &data_recd->sec_of_day, &data_recd->pressure, &data_recd->temperature,
                       &data_recd->humidity, &data_recd->value_origin);
        return (nstat+1);
    }
    
    /* 21 - Meteorological Supplement Record */
    int GSLRcrd::read_21 (char * str, struct rd21 *data_recd)
    {
        char temp_precip_type[256];
        int nstat;
        
        nstat= sscanf (str,
                       "%*s %Lf %lf %lf %s %d %lf %d %d",
                       &data_recd->sec_of_day, &data_recd->wind_speed,
                       &data_recd->wind_direction, &temp_precip_type,
                       &data_recd->visibility, &data_recd->sky_clarity,
                       &data_recd->atmospheric_seeing, &data_recd->cloud_cover);
        strncpy (data_recd->precip_type, temp_precip_type, 40);
        data_recd->precip_type[40]= '\0';
        return (nstat+1);
    }
    
    /* 30 - Pointing Angles Record */
    int GSLRcrd::read_30 (char * str, struct rd30 *data_recd)
    {
        int nstat;
        
        nstat= sscanf (str,
                       "%*s %Lf %lf %lf %d %d %d",
                       &data_recd->sec_of_day, &data_recd->azimuth, &data_recd->elevation, 
                       &data_recd->direction_ind, &data_recd->angle_origin_ind, 
                       &data_recd->refraction_corr_ind);
        return (nstat+1);
    }
    
    /* 40 - Calibration Record */
    int GSLRcrd::read_40 (char * str, struct rd40 *data_recd)
    {
        char temp_sysconfig_id[256];
        int nstat;
        
        nstat= sscanf (str,
                       "%*s %Lf %d %s %d %d %lf %lf %lf %lf %lf %lf %lf %d %d %d",
                       &data_recd->sec_of_day, &data_recd->type_of_data, 
                       &temp_sysconfig_id, &data_recd->num_points_recorded, 
                       &data_recd->num_points_used, &data_recd->one_way_target_dist, 
                       &data_recd->cal_sys_delay, &data_recd->cal_delay_shift, 
                       &data_recd->cal_rms, &data_recd->cal_skew, &data_recd->cal_kurtosis,
                       &data_recd->cal_PmM, &data_recd->cal_type_ind, 
                       &data_recd->cal_shift_type_ind, &data_recd->detector_channel);
        strncpy (data_recd->sysconfig_id, temp_sysconfig_id, 40);
        data_recd->sysconfig_id[40]= '\0';
        return (nstat+1);
    }
    
    /* 50 - Session Statistics Record */
    int GSLRcrd::read_50 (char * str, struct rd50 *data_recd)
    {
        char temp_sysconfig_id[256];
        int nstat;
        
        nstat= sscanf (str,
                       "%*s %s %lf %lf %lf %lf %d",
                       &temp_sysconfig_id, &data_recd->sess_rms, &data_recd->sess_skew, 
                       &data_recd->sess_kurtosis, &data_recd->sess_PmM, 
                       &data_recd->data_qual_ind);
        strncpy (data_recd->sysconfig_id, temp_sysconfig_id, 40);
        data_recd->sysconfig_id[40]= '\0';
        return (nstat+1);
    }
    
    /* 60 - Compatibility Record */
    int GSLRcrd::read_60 (char * str, struct rd60 *data_recd)
    {
        char temp_sysconfig_id[256];
        int nstat;
        
        nstat= sscanf (str,
                       "%*s %s %d %d", 
                       &temp_sysconfig_id, &data_recd->sys_change_ind, 
                       &data_recd->sys_config_ind); 
        strncpy (data_recd->sysconfig_id, temp_sysconfig_id, 40);
        data_recd->sysconfig_id[40]= '\0';
        return (nstat+1);
    }
    
    /* 9X - User Defined Records 90-99 */
    int GSLRcrd::read_9x (char * str, struct rd9x *data_recd)
    {
        int nstat= 0;
        return (nstat);
    }
    
    /* 00 - Comment Record */
    int GSLRcrd::read_00 (char * str, struct rd00 *data_recd)
    {
        int nstat= 1;
        strncpy(str,data_recd->comment,80);
        data_recd->comment[80]= '\0';
        
        return (nstat);
    }

    
    /* Ranging data header/footer records */
    /* H1 - format header */
    void GSLRcrd::write_h1 (FILE * str_out, struct rh1 header)
    {
        sprintf (stro,
                 "h1     %2d %4d %2d %2d %2d\n",
                 header.format_version, header.prod_year, header.prod_mon,
                 header.prod_day, header.prod_hour);
        strncpy (&stro[3], header.crd_literal, 3);
        fputs (stro, str_out);
        /*  fprintf (str_out, "%s\n", stro);*/
    }
    
    /* H2 - station header */
    void GSLRcrd::write_h2 (FILE * str_out, struct rh2 header)
    {
        sprintf (stro,
                 "h2            %4d %2d %2d %2d",
                 header.cdp_pad_id, header.cdp_sys_num, header.cdp_occ_num,
                 header.stn_timescale);
        strncpy (&stro[3], header.stn_name, 10);
        fprintf (str_out, "%s\n", stro);
    }
    
    /* H3 - spacecraft header */
    void GSLRcrd::write_h3 (FILE * str_out, struct rh3 header)
    {
        sprintf (stro,
                 "h3             %07d %4d %8d %1d %1d", header.ilrs_id, header.sic,
                 header.norad, header.SC_timescale, header.target_type);
        strncpy (&stro[3], header.target_name, 10);
        fprintf (str_out, "%s\n", stro);
    }
    
    /* H4 - Session header */
    void GSLRcrd::write_h4 (FILE * str_out, struct rh4 header)
    {
        fprintf (str_out,
                 "h4 %2d %4d %2d %2d %2d %2d %2d %4d %2d %2d %2d %2d %2d %2d %1d %1d %1d %1d %1d %1d %1d\n",
                 header.data_type,
                 header.start_year, header.start_mon, header.start_day,
                 header.start_hour, header.start_min, header.start_sec,
                 header.end_year, header.end_mon, header.end_day,
                 header.end_hour, header.end_min, header.end_sec,
                 header.data_release, header.refraction_app_ind,
                 header.CofM_app_ind, header.xcv_amp_app_ind,
                 header.stn_sysdelay_app_ind,
                 header.SC_sysdelay_app_ind,
                 header.range_type_ind, header.data_qual_alert_ind);
    }
    
    /* Need indicators that these have been read? */
    /* H8 - End of Session footer */
    void GSLRcrd::write_h8 (FILE * str_out)
    {
        fprintf (str_out, "h8\n");
    }
    
    /* H9 - End of File footer */
    void GSLRcrd::write_h9 (FILE * str_out)
    {
        fprintf (str_out, "h9\n");
    }
    
    /* Ranging data configuration records (1 of n) */
    /* C0 - System Configuration Record */
    void GSLRcrd::write_c0 (FILE * str_out, struct rc0 config)
    {
        int i;
        
        fprintf (str_out,
                 "c0 %1d %.3f %-s",
                 config.detail_type, config.xmit_wavelength, config.config_ids[0]);
        for (i=1; i<10; i++)
        {
            if (config.config_ids[i][0] != '\0')
                fprintf (str_out, " %-s", config.config_ids[i]);
        }
        fprintf (str_out, "\n");
    }
    
    /* C1 - Laser Configuration Record */
    void GSLRcrd::write_c1 (FILE * str_out, struct rc1 config)
    {
        fprintf (str_out,
                 "c1 %d %-s %-s %.2f %.2f %.2f %.1f %.2f %d\n",
                 config.detail_type, config.laser_config_id, config.laser_type,
                 config.prim_wavelength, config.nom_fire_rate, config.pulse_energy,
                 config.pulse_width, config.beam_div, config.pulses_in_semitrain);
    }
    
    /* C2 - Detector Configuration Record */
    void GSLRcrd::write_c2 (FILE * str_out, struct rc2 config)
    {
        fprintf (str_out,
                 "c2 %d %-s %-s %.3f %.2f %.1f %.1f %-s %.1f %.2f %.1f %.1f %-s\n",
                 config.detail_type, config.detector_config_id, config.detector_type,
                 config.app_wavelength, config.qe, config.voltage, config.dark_count,
                 config.output_pulse_type, config.output_pulse_width,
                 config.spectral_filter, config.spectral_filter_xmission,
                 config.spatial_filter, config.signal_proc);
    }
    
    /* C3 - Timing Configuration Record */
    void GSLRcrd::write_c3 (FILE * str_out, struct rc3 config)
    {
        fprintf (str_out,
                 "c3 %d %-s %-s %-s %-s %-s %.1f\n",
                 config.detail_type, config.timing_config_id, config.time_source,
                 config.freq_source, config.timer, config.timer_serial_num,
                 config.epoch_delay_corr);
    }
    
    /* C4 - Transponder Configuration Record */
    void GSLRcrd::write_c4 (FILE * str_out, struct rc4 config)
    {
        fprintf (str_out,
                 "c4 %d %-s %.3Lf %.2f %.3Lf %.2f %.12Lf %d %d %d\n",
                 config.detail_type, config.xponder_config_id,
                 config.est_stn_utc_offset, config.est_stn_osc_drift,
                 config.est_xponder_utc_offset, config.est_xponder_osc_drift,
                 config.xponder_clock_ref_time,
                 config.stn_off_drift_app_ind,
                 config.SC_off_drift_app_ind,
                 config.SC_time_simplified_ind);
    }
    
    /* Ranging data records */
    /* 10 - Range Record */
    void GSLRcrd::write_10 (FILE * str_out, struct rd10 data_recd)
    {
        fprintf (str_out,
                 "10 %.12Lf %.12Lf %-s %d %d %d %d %d\n",
                 data_recd.sec_of_day, data_recd.time_of_flight,
                 data_recd.sysconfig_id, data_recd.epoch_event,
                 data_recd.filter_flag, data_recd.detector_channel,
                 data_recd.stop_number, data_recd.xcv_amp);
    }
    
    /* 11 - Normal Point Record */
    void GSLRcrd::write_11 (FILE * str_out, struct rd11 data_recd)
    {
        fprintf (str_out,
                 "11 %.12Lf %.12Lf %-s %d %.1f %d %.1f %.3f %.3f %.1f %.1f %d\n",
                 data_recd.sec_of_day, data_recd.time_of_flight,
                 data_recd.sysconfig_id, data_recd.epoch_event,
                 data_recd.np_window_length, data_recd.num_ranges,
                 data_recd.bin_rms, data_recd.bin_skew, data_recd.bin_kurtosis,
                 data_recd.bin_PmM, data_recd.return_rate,
                 data_recd.detector_channel);
    }
    
    /* 12 - Range Supplement Record */
    void GSLRcrd::write_12 (FILE * str_out, struct rd12 data_recd)
    {
        fprintf (str_out,
                 "12 %.7Lf %-s %.1f %.4f %.2f %.4f\n",
                 data_recd.sec_of_day, data_recd.sysconfig_id,
                 data_recd.refraction_corr, data_recd.target_CofM_corr,
                 data_recd.nd_value, data_recd.time_bias);
    }
    
    /* 20 - Meteorological Record */
    void GSLRcrd::write_20 (FILE * str_out, struct rd20 data_recd)
    {
        fprintf (str_out,
                 "20 %.3Lf %.2f %.2f %.0f %d\n",
                 data_recd.sec_of_day, data_recd.pressure, data_recd.temperature,
                 data_recd.humidity, data_recd.value_origin);
    }
    
    /* 21 - Meteorological Supplement Record */
    void GSLRcrd::write_21 (FILE * str_out, struct rd21 data_recd)
    {
        fprintf (str_out,
                 "21 %.3Lf %.1f %.1f %-s %d %.2f %d %d\n",
                 data_recd.sec_of_day, data_recd.wind_speed, data_recd.wind_direction,
                 data_recd.precip_type, data_recd.visibility, data_recd.sky_clarity,
                 data_recd.atmospheric_seeing, data_recd.cloud_cover);
    }
    
    /* 30 - Pointing Angles Record */
    void GSLRcrd::write_30 (FILE * str_out, struct rd30 data_recd)
    {
        fprintf (str_out,
                 "30 %.3Lf %.4f %.4f %d %d %d\n",
                 data_recd.sec_of_day, data_recd.azimuth, data_recd.elevation,
                 data_recd.direction_ind, data_recd.angle_origin_ind,
                 data_recd.refraction_corr_ind);
    }
    
    /* 40 - Calibration Record */
    void GSLRcrd::write_40 (FILE * str_out, struct rd40 data_recd)
    {
        /**
         fprintf (str_out,
         "40 %.7Lf %d %-s %d %d %.3f %.1f %.1f %.1f %.3f %.3f %.1f %d %d %d\n",
         data_recd.sec_of_day, data_recd.type_of_data, data_recd.sysconfig_id,
         data_recd.num_points_recorded, data_recd.num_points_used,
         data_recd.one_way_target_dist, data_recd.cal_sys_delay,
         data_recd.cal_delay_shift, data_recd.cal_rms, data_recd.cal_skew,
         data_recd.cal_kurtosis, data_recd.cal_PmM, data_recd.cal_type_ind,
         data_recd.cal_shift_type_ind, data_recd.detector_channel);
         **/
        /* Due to some problem w/ fedora 8 or crd_cal.... */
        fprintf (str_out,
                 "40 %.7Lf %1d %-s %d %d %.3f",
                 data_recd.sec_of_day, data_recd.type_of_data, data_recd.sysconfig_id,
                 data_recd.num_points_recorded, data_recd.num_points_used,
                 data_recd.one_way_target_dist);
        fprintf (str_out,
                 " %.1f %.1f %.1f %.3f %.3f %.1f %1d %1d %1d\n",
                 data_recd.cal_sys_delay,
                 data_recd.cal_delay_shift, data_recd.cal_rms, data_recd.cal_skew,
                 data_recd.cal_kurtosis, data_recd.cal_PmM, data_recd.cal_type_ind,
                 data_recd.cal_shift_type_ind, data_recd.detector_channel);
        
    }
    
    /* 50 - Session Statistics Record */
    void GSLRcrd::write_50 (FILE * str_out, struct rd50 data_recd)
    {
        fprintf (str_out,
                 "50 %-s %.1f %.3f %.3f %.1f %d\n",
                 data_recd.sysconfig_id, data_recd.sess_rms, data_recd.sess_skew, 
                 data_recd.sess_kurtosis, data_recd.sess_PmM, 
                 data_recd.data_qual_ind);
    }
    
    /* 60 - Compatibility Record */
    void GSLRcrd::write_60 (FILE * str_out, struct rd60 data_recd)
    {
        fprintf (str_out, 
                 "60 %-s %d %d\n", 
                 data_recd.sysconfig_id, data_recd.sys_change_ind, 
                 data_recd.sys_config_ind); 
    }
    
    /* 9X - User Defined Record */
    void GSLRcrd::write_9x (FILE * str_out, struct rd9x data_recd)
    {
        
    }
    
    /* 00 - Comment Record */
    void GSLRcrd::write_00 (FILE * str_out, struct rd00 data_recd)
    {
        fprintf (str_out,
                 "00 %-s\n", 
                 data_recd.comment);
    }

    
    void GSLRcrd::closedatbase()
    {
        slrcrdFile.close();
    }
    
    void GSLRcrd::setdatabase(gfc::GString filename)
    {
        slrcrdFile.open(filename);
        
        if( slrcrdFile.is_open() == false )
        {
            printf("slr crd file is unavailable\n");
            return;
        }
    }
    
    
    // the next obs
    bool GSLRcrd::nextDataBlock(GslrStorage& slrstorage)
    {
        char str[512] = {0};
        
        bool available = false;
        
        GSLRBlock slrstation;
        
        bool startRecord = false;
        GSensorID sensor;
        
        //station calibration data
        rd40 d40;
        
        while( ! slrcrdFile.eof() )
        {
            
            if(available == true)
            {
                break;
            }
            
            //getline( slrcrdFile, line ); // epoch time line
            slrcrdFile.getline(str, 512);
            
            if( strncmp( str, "H1", 2) == 0
               ||strncmp(str, "h1", 2) == 0
               )
            {
                read_h1(str, &slrstation.h1);
            }
            
            else if (strncmp (str, "H2", 2) == 0 ||
                     strncmp (str, "h2", 2) == 0)
            {
                read_h2 (str, &slrstation.h2);
            }
            else if (strncmp (str, "H3", 2) == 0 ||
                     strncmp (str, "h3", 2) == 0)
            {
                read_h3 (str, &slrstation.h3);
                sensor = GSensorID::ilrsName2GSensorID(slrstation.h3.target_name);
                if(sensor.Available() == true)
                {
                    startRecord = true;
                }
            }
            else if (strncmp (str, "H4", 2) == 0 ||
                     strncmp (str, "h4", 2) == 0)
            {
                read_h4 (str, &slrstation.h4);
            }
            
            // end of section
            else if (strncmp (str, "H8", 2) == 0 ||
                     strncmp (str, "h8", 2) == 0)
            {
                read_h8 (str);
                
              if( startRecord ==  true )
                {
                    //int sisze = slrstorage.size();
                    
                    // the end of data for this station, need to push_back in storage
                    //slrstorage[sensor] = slrstation;
                    slrstorage.push_data(slrstation);
                    
                    available = true;
                    
                    break;
                }
                else
                {
                    available = false;
                }
            }
            // end of file
            else if (strncmp (str, "H9", 2) == 0 ||
                     strncmp (str, "h9", 2) == 0)
            {
                read_h9 (str);
                break;
            }
            
            else if (strncmp (str, "C0", 2) == 0 ||
                     strncmp (str, "c0", 2) == 0)
            {
                read_c0 (str, &slrstation.c0);
            }
            else if (strncmp (str, "C1", 2) == 0 ||
                     strncmp (str, "c1", 2) == 0)
            {
                read_c1 (str, &slrstation.c1);
            }
            else if (strncmp (str, "C2", 2) == 0 ||
                     strncmp (str, "c2", 2) == 0)
            {
                read_c2 (str, &slrstation.c2);
            }
            else if (strncmp (str, "C3", 2) == 0 ||
                     strncmp (str, "c3", 2) == 0)
            {
                read_c3 (str, &slrstation.c3);
            }
            else if (strncmp (str, "C4", 2) == 0 ||
                     strncmp (str, "c4", 2) == 0)
            {
                read_c4 (str, &slrstation.c4);
            }
            else if (strncmp (str, "10", 2) == 0)
            {
                rd10 d10;
                read_10 (str, &d10);
                
                if( startRecord == true )
                {
                    GslrData slrdata;
                    
//                    double delt = d10.sec_of_day - d40.sec_of_day;
//                    //d10.sec_of_day += (d40.cal_sys_delay + delt*d40.cal_delay_shift)*1.0E-12;
//                    double bias = (d40.cal_sys_delay + delt*d40.cal_delay_shift)*1.0E-12;
//                    d10.sec_of_day += 2.0*bias;
                    
                    CivilTime ct;
                    ct.m_ts = TimeSystem("tsUTC");
                    ct.m_year = slrstation.h4.start_year;
                    ct.m_month = slrstation.h4.start_mon;
                    ct.m_day = slrstation.h4.start_day;
                    ct.m_hour = int(d10.sec_of_day/3600.0);
                    ct.m_minute = int((d10.sec_of_day - ct.m_hour * 3600.0)/60.0);
                    ct.m_second = d10.sec_of_day - ct.m_hour*3600.0 - ct.m_minute * 60.0;
                    
                    slrdata.obsepoch = GTime::CivilTime2GTime(ct);
                    
                    slrdata.m_data_value = d10.time_of_flight  ;
                    
                    slrdata.epoch_event = d10.epoch_event;
                    
                    slrstation.obsdata.push_back(slrdata);
                    
                }
                
            }
            else if (strncmp (str, "11", 2) == 0)
            {
                rd11 d11;
                read_11 (str, &d11);
                
                // data type convertion
                if(startRecord == true)
                {
                    GslrData slrdata;
                    
//                    double delt = d11.sec_of_day - d40.sec_of_day;
//                    //d11.sec_of_day += (d40.cal_sys_delay + delt*d40.cal_delay_shift)*1.0E-12;
//                    double bias = (d40.cal_sys_delay + delt*d40.cal_delay_shift)*1.0E-12;
//                    d11.sec_of_day += 2.0*bias;
//                    
                    
                    CivilTime ct;
                    ct.m_ts = TimeSystem("tsUTC");
                    ct.m_year = slrstation.h4.start_year;
                    ct.m_month = slrstation.h4.start_mon;
                    ct.m_day = slrstation.h4.start_day;
                    ct.m_hour = int(d11.sec_of_day/3600.0);
                    ct.m_minute = int((d11.sec_of_day - ct.m_hour * 3600.0)/60.0);
                    ct.m_second = d11.sec_of_day - ct.m_hour*3600.0 - ct.m_minute * 60.0;
                    
                    slrdata.obsepoch = GTime::CivilTime2GTime(ct);
                   
                    slrdata.m_data_value = d11.time_of_flight  ;
                    
                    slrdata.epoch_event = d11.epoch_event;
                    
                    slrstation.obsdata.push_back(slrdata);
                    
                }
                
                
            }
            else if (strncmp (str, "12", 2) == 0)
            {
                rd12  d12;
                read_12 (str, &d12);
                
            }
            else if (strncmp (str, "20", 2) == 0)
            {
                rd20 d20;
                read_20 (str, &d20);
                
                slrstation.MeteorRecord.push_back(d20);
                
            }
            else if (strncmp (str, "21", 2) == 0)
            {
                rd21 d21;
                read_21 (str, &d21);
                
            }
            else if (strncmp (str, "30", 2) == 0)
            {
                rd30 d30;
                
                read_30 (str, &d30);
                
            }
            else if (strncmp (str, "40", 2) == 0)
            {
                
                read_40 (str, &d40);
                
            }
            else if (strncmp (str, "50", 2) == 0)
            {
                rd50 d50;
                
                read_50 (str, &d50);
                
            }
            else if (strncmp (str, "60", 2) == 0)
            {
                rd60 d60;
                
                read_60 (str, &d60);
                
            }
            else if (strncmp (str, "00", 2) == 0)
            {
                rd00 d00;
                
                read_00 (str, &d00);
                
            }
            
            
        }
        
        
        return available;
    }
    
    
    
}  // end of namespace gfc

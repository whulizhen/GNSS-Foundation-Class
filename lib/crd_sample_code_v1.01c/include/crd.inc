C  Consolidated laser Ranging Data format (CRD)
C      record and variable definitions for FORTRAN
C      R. Ricklefs UT/CSR July 2007
C  History:
C  08/xx/07 - added H3 Target type (v0.27)
C  11/26/07 - added H4 data quality alert
C             and #10 stop number
C             and #20 origin of values (v0.27) rlr.
C  12/12/07 - commons were re-ordered double precision, integer, character
C            for efficiency and to keep the compiler happy. :-(
C  05/07/08 - Expand configuration and data record character fields to 
C             allow up to 40 characters.
C             - Added detector channel to normalpoint (11) and calibration (40)
C             records. 
C             - Added field for 'crd' literal to 'h1'.
C             - Record '21' sky_clarity is now double rather than int.
C               (v1.00 rlr)
C
C
C ----------------------------------------------------------------------
C
C  Header Records
C
C H1 - format header 
      COMMON /h1/ crd_literal, format_version, 
     &            prod_year, prod_mon, prod_day, prod_hour
      character*3 crd_literal
      integer format_version, 
     &            prod_year, prod_mon, prod_day, prod_hour

C H2 - station header
      COMMON /h2/ cdp_pad_id, cdp_sys_num, cdp_occ_num, stn_timescale,
     &            stn_name
      character*10 stn_name
      integer cdp_pad_id, cdp_sys_num, cdp_occ_num, stn_timescale

C H3 - spacecraft header
      COMMON /h3/ ilrs_id, sic, norad, SC_timescale, target_type,
     &            target_name
      character*10 target_name
      integer ilrs_id, sic, norad, SC_timescale, target_type

C H4 - Session header
      COMMON /h4/ data_type, start_year, start_mon, start_day,
     &            start_hour, start_min, start_sec,
     &            end_year, end_mon, end_day, 
     &            end_hour, end_min, end_sec,
     &            data_release, refraction_app_ind, CofM_app_ind,
     &            xcv_amp_app_ind, stn_sysdelay_app_ind,
     &            SC_sysdelay_app_ind, range_type_ind, 
     &            data_qual_alert_ind
      integer data_type, start_year, start_mon, start_day
      integer start_hour, start_min, start_sec
      integer end_year, end_mon, end_day, end_hour, end_min, end_sec
      integer data_release, refraction_app_ind, CofM_app_ind
      integer xcv_amp_app_ind, stn_sysdelay_app_ind
      integer SC_sysdelay_app_ind, range_type_ind, data_qual_alert_ind

C H8 - End of Session footer

C H9 - End of File footer

C
C  Configuration Records
C
C C0 - System Configuration Record
      COMMON /c0/ xmit_wavelength, c0_detail_type, config_ids
      integer c0_detail_type
      double precision xmit_wavelength
      character*40 config_ids(5)	! May be more later

C C1 - Laser Configuration Record
      COMMON /c1/ prim_wavelength, nom_fire_rate, pulse_energy, 
     &            pulse_width, beam_div, 
     &            c1_detail_type, pulses_in_semitrain,
     &            laser_type, laser_config_id
      integer c1_detail_type, pulses_in_semitrain
      character laser_config_id*40, laser_type*40
      double precision prim_wavelength, nom_fire_rate, pulse_energy
      double precision pulse_width, beam_div

C C2 - Detector Configuration Record
      COMMON /c2/ app_wavelength, qe, voltage, dark_count,
     &            output_pulse_width, spectral_filter,
     &            spectral_filter_xmission, spatial_filter, 
     &            c2_detail_type, signal_proc, 
     &            detector_config_id, detector_type, output_pulse_type
      integer c2_detail_type
      character detector_config_id*40, detector_type*40
      character output_pulse_type*40, signal_proc*40
      double precision app_wavelength, qe, voltage, dark_count
      double precision output_pulse_width, spectral_filter
      double precision spectral_filter_xmission, spatial_filter

C C3 - Timing Configuration Record
      COMMON /c3/ c3_detail_type, timing_config_id, time_source,
     &            freq_source, timer, timer_serial_num, epoch_delay_corr
      integer c3_detail_type
      character timing_config_id*40, time_source*40
      character freq_source*40, timer*40, timer_serial_num*40
      double precision epoch_delay_corr

C C4 - Transponder Configuration Record
      COMMON /c4/ 
     &            est_stn_utc_offset, est_stn_osc_drift,
     &            est_xponder_utc_offset, est_xponder_osc_drift,
     &            xponder_clock_ref_time, 
     &            c4_detail_type, stn_off_drift_app_ind,
     &            SC_off_drift_app_ind, SC_time_simplified_ind, 
     &            xponder_config_id
      integer c4_detail_type
      integer stn_off_drift_app_ind, SC_off_drift_app_ind
      integer SC_time_simplified_ind
      character xponder_config_id*40
      double precision est_stn_osc_drift, est_xponder_osc_drift
C The next 3 variables sould be "long double" or "real*16"
      double precision est_stn_utc_offset, est_xponder_utc_offset
      double precision xponder_clock_ref_time
C
C  Data Records
C
C 10 - Range Record
        COMMON /d10/ d10_sec_of_day, d10_time_of_flight,
     &               d10_epoch_event, filter_flag,
     &               d10_detector_channel, stop_number, xcv_amp, 
     &               d10_sysconfig_id
      character d10_sysconfig_id*40
      real*8 d10_sec_of_day, d10_time_of_flight
      integer d10_epoch_event, filter_flag, d10_detector_channel
      integer xcv_amp
      integer stop_number

C 11 - Normal Point Record
      COMMON /d11/ d11_sec_of_day, d11_time_of_flight,
     &             d11_epoch_event, np_window_length,
     &             bin_rms, bin_skew, bin_kurtosis, bin_PmM, 
     &             return_rate, d11_sysconfig_id, num_ranges, 
     &             d11_detector_channel 
      character*40 d11_sysconfig_id;
      double precision d11_sec_of_day, d11_time_of_flight
      double precision np_window_length
      double precision bin_rms, bin_skew, bin_kurtosis, bin_PmM
      double precision return_rate
      integer d11_epoch_event, num_ranges
      integer d11_detector_channel

C 12 - Range Supplement Record
      COMMON /d12/ d12_sec_of_day, refraction_corr, 
     &             target_CofM_corr, nd_value, time_bias, 
     &             d12_sysconfig_id
       character*40 d12_sysconfig_id
       double precision d12_sec_of_day, refraction_corr
       double precision target_CofM_corr, nd_value, time_bias

C 20 - Meteorological Record
      COMMON /d20/  d20_sec_of_day, pressure,
     &              temperature, humidity, value_origin
      real*8  d20_sec_of_day, pressure, temperature, humidity
      integer value_origin

C 21 - Meteorological Supplement Record
      COMMON /D21/ d21_sec_of_day, wind_speed, wind_direction,
     &             sky_clarity, visibility, atmospheric_seeing, 
     &             cloud_cover, precip_type
      character precip_type*40
      double precision d21_sec_of_day, wind_speed, wind_direction
      double precision sky_clarity
      integer visibility, atmospheric_seeing, cloud_cover

C 30 - Pointing Angles Record
      COMMON /D30/ d30_sec_of_day, azimuth, elevation, direction_ind,
     &             angle_origin_ind, refraction_corr_ind
      double precision d30_sec_of_day, azimuth, elevation
      integer direction_ind, angle_origin_ind, refraction_corr_ind

C 40 - Calibration Record
      COMMON /d40/ d40_sec_of_day, type_of_data, d40_sysconfig_id,
     &             num_points_recorded, num_point_used,
     &             one_way_target_dist, cal_sys_delay, cal_delay_shift,
     &             cal_rms, cal_skew, cal_kurtosis, cal_PmM, 
     &             cal_type_ind, cal_shift_type_ind, 
     &             d40_detector_channel
      character*40 d40_sysconfig_id
      real*8 d40_sec_of_day
      real*8 one_way_target_dist, cal_sys_delay, cal_delay_shift
      real*8 cal_rms, cal_skew, cal_kurtosis, cal_PmM
      integer type_of_data, num_points_recorded, num_point_used
      integer cal_type_ind, cal_shift_type_ind, d40_detector_channel

C 50 - Session Statistics Record
      COMMON /d50/ sess_rms, sess_skew, sess_kurtosis, sess_PmM, 
     &             data_qual_ind, d50_sysconfig_id
      character*40 d50_sysconfig_id
      double precision sess_rms, sess_skew, sess_kurtosis, sess_PmM
      integer data_qual_ind

C 60 - Compatibility Record
      COMMON /d60/ sys_change_ind,sys_config_ind, d60_sysconfig_id
      character*40 d60_sysconfig_id
      integer sys_change_ind,sys_config_ind

C 9X - User Defined Record

C 00 - Comment Record
      COMMON /d00/ comment
      character*80 comment

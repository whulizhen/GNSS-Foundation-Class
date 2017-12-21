C**-------------------------------------------------------------------------
C * Subroutines: write CRD data records to an output file
C *
C * Author: Randall Ricklefs / Univ of Texas / Center for Space Research
C *
C * History:
C *   July 06, 2007 - Initial version
C  05/07/08   - Expand configuration and data record character fields to
C               allow up to 40 characters.
C             - Added detector channel to normalpoint (11) and calibration (40)
C               records.
C             - Added field for 'crd' literal to 'h1'.
C             - Record '21' sky_clarity is not double rather than int.
C  06/24/08   - Record '11' np window length is now double rather than
C               int. (v1.00 rlr)
C  03/10/09   - Record H2 Epoch Timescale corrected from I1 to I2.
C               (v1.00a rlr).
C  03/10/09   - Record H3 changed to print leading zeros rather than
C               spaces for ilrs_id. (v1.00a rlr).
C *
C**-------------------------------------------------------------------------

C Ranging data header/footer records
C H1 - format header
      SUBROUTINE write_h1 (str)
      CHARACTER*512 str
      INCLUDE '../include/crd.inc'

      WRITE (str,1000,err=100) crd_literal, format_version,
     &           prod_year, prod_mon, prod_day, prod_hour
      RETURN

 100  WRITE(*,*) "Error writing CRD record type h1"
      STOP 1

 1000 FORMAT ("h1",1x,a3,1x,i2,1x,i4,1x,i2,1x,i2,1x,i2)
      END

C H2 - station header
      SUBROUTINE write_h2 (str)
      CHARACTER*512 str
      INCLUDE '../include/crd.inc'

      WRITE (str,1000,err=100) stn_name,
     &           cdp_pad_id, cdp_sys_num, cdp_occ_num, stn_timescale
      RETURN

 100  WRITE(*,*) "Error writing CRD record type h2"
      STOP 1

 1000 FORMAT ("h2",1x,a10,1x,i4,1x,i2,1x,i2,1x,i2)
      END

C H3 - spacecraft header
      SUBROUTINE write_h3 (str)
      CHARACTER*512 str
      INCLUDE '../include/crd.inc'

      WRITE (str,1000,err=100) target_name,
     &           ilrs_id, sic, norad, SC_timescale, target_type
      RETURN

 100  WRITE(*,*) "Error writing CRD record type h3"
      STOP 1

 1000 FORMAT ("h3",1x,a10,1x,i8.7,1x,i4,1x,i8,1x,i1,1x,i1)
      END

C H4 - Session header
      SUBROUTINE write_h4 (str)
      CHARACTER*512 str
      INCLUDE '../include/crd.inc'

      WRITE (str,1000,err=100)
     &           data_type, start_year, start_mon, start_day,
     &           start_hour, start_min, start_sec,
     &           end_year, end_mon, end_day, end_hour, end_min, end_sec,
     &           data_release, refraction_app_ind, CofM_app_ind,
     &           xcv_amp_app_ind, stn_sysdelay_app_ind,
     &           SC_sysdelay_app_ind, range_type_ind, 
     &           data_qual_alert_ind

      RETURN

 100  WRITE(*,*) "Error writing CRD record type h4"
      STOP 1

 1000 FORMAT ("h4",1x,i2,2(1x,i4,5(1x,i2)),1x,i2,7(1x,i1))
      END

C H8 - End of Session footer
      SUBROUTINE write_h8 (str)
      CHARACTER*512 str
      INCLUDE '../include/crd.inc'

      WRITE (str,1000,err=100)
      RETURN

 100  WRITE(*,*) "Error writing CRD record type h8"
      STOP 1

 1000 FORMAT ("h8")
      END

C H9 - End of File footer
      SUBROUTINE write_h9 (str)
      CHARACTER*512 str
      INCLUDE '../include/crd.inc'

      WRITE (str,1000,err=100)
      RETURN

 100  WRITE(*,*) "Error writing CRD record type h9"
      STOP 1

 1000 FORMAT ("h9")
      END

C Ranging data configuration records (1 of n)
C C0 - System Configuration Record
      SUBROUTINE write_c0 (str)
      CHARACTER*512 str
      INTEGER trimlen, ci_len(4), i        ! May be more than 4 later
      INCLUDE '../include/crd.inc'

      do i=1,4
        ci_len(i)= trimlen(config_ids(i))
      enddo
      
      WRITE (str,1000,err=100) c0_detail_type, xmit_wavelength, 
     &           (config_ids(i)(1:ci_len(i)),i=1,4)
      RETURN

 100  WRITE(*,*) "Error writing CRD record type c0"
      STOP 1

 1000 FORMAT ("c0 ",i1,1x,f8.3,4(1x,a))

      END

C C1 - Laser Configuration Record
      SUBROUTINE write_c1 (str)
      CHARACTER*512 str
      INTEGER trimlen
      INTEGER lci_len, lt_len
      INCLUDE '../include/crd.inc'

      lci_len= trimlen(laser_config_id)
      lt_len= trimlen(laser_type)
      WRITE (str,1000,err=100) 
     &           c1_detail_type, laser_config_id(1:lci_len), 
     &           laser_type(1:lt_len), prim_wavelength, nom_fire_rate, 
     &           pulse_energy, pulse_width, beam_div, 
     &           pulses_in_semitrain
      RETURN

 100  WRITE(*,*) "Error writing CRD record type c1"
      STOP 1

 1000 FORMAT ("c1 ",i1,1x,a,1x,a,1x,f10.2,1x,f10.2,1x,f10.2,1x,
     &          f6.1,1x,f5.2,1x,i4)
      END

C C2 - Detector Configdetector_type,uration Record
      SUBROUTINE write_c2 (str)
      CHARACTER*512 str
      INTEGER trimlen
      INTEGER dci_len, dt_len, opt_len, sp_len
      INCLUDE '../include/crd.inc'

      dci_len= trimlen(detector_config_id)
      dt_len= trimlen(detector_type)
      opt_len= trimlen(output_pulse_type)
      sp_len= trimlen(signal_proc)
      WRITE (str,1000,err=100) 
     &           c2_detail_type, detector_config_id(1:dci_len), 
     &           detector_type(1:dt_len), app_wavelength, qe, voltage, 
     &           dark_count, output_pulse_type(1:opt_len), 
     &           output_pulse_width, spectral_filter,
     &           spectral_filter_xmission, spatial_filter, 
     &           signal_proc(1:sp_len)
      RETURN

 100  WRITE(*,*) "Error writing CRD record type c2"
      STOP 1

 1000 FORMAT ("c2 ",i1,1x,a,1x,a,1x,f10.3,1x,f5.1,1x,f6.1,1x,f5.1,
     &          1x,a,1x,f5.1,1x,f5.2,1x,f5.1,1x,f5.2,1x,a)
      END

C C3 - Timing Configuratiming_config_id,tion Record
      SUBROUTINE write_c3 (str)
      CHARACTER*512 str
      INTEGER trimlen
      INTEGER tci_len, ts_len, fs_len, t_len, tsn_len
      INCLUDE '../include/crd.inc'

      tci_len= trimlen(timing_config_id)
      ts_len= trimlen(time_source)
      fs_len= trimlen(freq_source)
      t_len= trimlen(timer)
      tsn_len= trimlen(timer_serial_num)
      WRITE (str,1000,err=100) 
     &           c3_detail_type, timing_config_id(1:tci_len), 
     &           time_source(1:ts_len), freq_source(1:fs_len), 
     &           timer(1:t_len), timer_serial_num(1:tsn_len), 
     &           epoch_delay_corr
      RETURN

 100  WRITE(*,*) "Error writing CRD record type c3"
      STOP 1

 1000 FORMAT ("c3 ",i1,5(1x,a),1x,f6.1)
      END

C C4 - Transponder Configuration Record
      SUBROUTINE write_c4 (str)
      CHARACTER*512 str
      INTEGER trimlen, xci_len
      INCLUDE '../include/crd.inc'

      xci_len= trimlen(xponder_config_id)
      WRITE (str,1000,err=100)
     &           c4_detail_type, xponder_config_id(1:xci_len),
     &           est_stn_utc_offset, est_stn_osc_drift,
     &           est_xponder_utc_offset, est_xponder_osc_drift,
     &           xponder_clock_ref_time, stn_off_drift_app_ind,
     &           SC_off_drift_app_ind, SC_time_simplified_ind
      RETURN

 100  WRITE(*,*) "Error writing CRD record type c4"
      STOP 1

 1000 FORMAT ("c4 ",i1,1x,a,1x,f20.3,1x,f11.2,1x,f20.3,1x,f11.2,
     &          f20.12,3(1x,i1))
      END

C Ranging data records
C 10 - Range Record
      SUBROUTINE write_10 (str)
      CHARACTER*512 str
      INTEGER trimlen, sci_len
      INCLUDE '../include/crd.inc'

      sci_len= trimlen(d10_sysconfig_id)
      WRITE (str,1000,err=100) d10_sec_of_day,d10_time_of_flight,
     &           d10_sysconfig_id(1:sci_len), d10_epoch_event, 
     &           filter_flag, d10_detector_channel, stop_number, xcv_amp
      RETURN

 100  WRITE(*,*) "Error writing CRD record type 10"
      STOP 1

 1000 FORMAT ("10 ",f18.12,1x,f18.12,1x,a,1x,4(i1,1x),i5)
      END

C 11 - Normal Point Record
      SUBROUTINE write_11 (str)
      CHARACTER*512 str
      INTEGER trimlen, sci_len
      INCLUDE '../include/crd.inc'

      sci_len= trimlen(d11_sysconfig_id)
      WRITE (str,1000,err=100)  d11_sec_of_day, d11_time_of_flight,
     &           d11_sysconfig_id(1:sci_len), d11_epoch_event, 
     &           np_window_length, num_ranges, bin_rms, bin_skew, 
     &           bin_kurtosis, bin_PmM, return_rate, 
     &           d11_detector_channel

      RETURN

 100  WRITE(*,*) "Error writing CRD record type 11"
      STOP 1

 1000 FORMAT ("11 ",f18.12,1x,f18.12,1x,a,1x,i1,1x,f6.1,1x,i6,
     &            1x,f6.1,1x,f7.3,1x,f7.3,1x,f9.1,1x,f6.2,1x,i1)
      END

C 12 - Range Supplement Record
      SUBROUTINE write_12 (str)
      CHARACTER*512 str
      INTEGER trimlen, sci_len
      INCLUDE '../include/crd.inc'

      sci_len= trimlen(d12_sysconfig_id)
      WRITE (str,1000,err=100) d12_sec_of_day, 
     &           d12_sysconfig_id(1:sci_len), refraction_corr, 
     &           target_CofM_corr, nd_value, time_bias

      RETURN

 100  WRITE(*,*) "Error writing CRD record type 12"
      STOP 1

 1000 FORMAT ("12 ",f18.12,1x,a,1x,f6.1,1x,f6.4,1x,f5.2,1x,f8.4)
      END

C 20 - Meteorological Record
      SUBROUTINE write_20 (str)
      CHARACTER*512 str
      INCLUDE '../include/crd.inc'

      WRITE (str,1000,err=100) d20_sec_of_day, 
     &            pressure, temperature, humidity, value_origin

      RETURN

 100  WRITE(*,*) "Error writing CRD record type 20"
      STOP 1

 1000 FORMAT ("20 ",f9.3,1x,f7.2,1x,f6.2,1x,f4.0,1x,i1)
      END

C 21 - Meteorological Supplement Record
      SUBROUTINE write_21 (str)
      CHARACTER*512 str
      INTEGER trimlen, pt_len
      INCLUDE '../include/crd.inc'

      pt_len= trimlen(precip_type)
      WRITE (str,1000,err=100) d21_sec_of_day, wind_speed, 
     &           wind_direction, precip_type(1:pt_len), visibility, 
     &           sky_clarity, atmospheric_seeing, cloud_cover

      RETURN

 100  WRITE(*,*) "Error writing CRD record type 21"
      STOP 1

 1000 FORMAT ("21 ",f9.3,1x,f5.1,1x,f5.1,1x,a,1x,i3,1x,f4.2,1x,
     &        i2,1x,i2)
      END

C 30 - Pointing Angles Record
      SUBROUTINE write_30 (str)
      CHARACTER*512 str
      INCLUDE '../include/crd.inc'

      WRITE (str,1000,err=100) d30_sec_of_day, azimuth, elevation, 
     &           direction_ind, angle_origin_ind, refraction_corr_ind

      RETURN

 100  WRITE(*,*) "Error writing CRD record type 30"
      STOP 1

 1000 FORMAT ("30 ",f9.3,2(1x,f8.4),3(1x,i1))
      END

C 40 - Calibration Record
      SUBROUTINE write_40 (str)
      CHARACTER*512 str
      INTEGER trimlen, sci_len
      INCLUDE '../include/crd.inc'

      sci_len= trimlen(d40_sysconfig_id)
      WRITE (str,1000,err=100) d40_sec_of_day, 
     &           type_of_data, d40_sysconfig_id(1:sci_len),
     &           num_points_recorded, num_point_used,
     &           one_way_target_dist, cal_sys_delay, cal_delay_shift,
     &           cal_rms, cal_skew, cal_kurtosis, cal_PmM, cal_type_ind,
     &           cal_shift_type_ind, d40_detector_channel

      RETURN

 100  WRITE(*,*) "Error writing CRD record type 40"
      STOP 1

 1000 FORMAT ("40 ",f18.12,1x,i1,1x,a,1x,i8,1x,i8,1x,f7.3,1x,f10.1,
     &            1x,f8.1,1x,f6.1,1x,f7.3,1x,f7.3,1x,f6.1,3(1x,i1))
      END

C 50 - Session Statistics Record
      SUBROUTINE write_50 (str)
      CHARACTER*512 str
      INTEGER trimlen, sci_len
      INCLUDE '../include/crd.inc'

      sci_len= trimlen(d50_sysconfig_id)
      WRITE (str,1000,err=100)
     &           d50_sysconfig_id(1:sci_len), sess_rms, sess_skew,
     &           sess_kurtosis, sess_PmM, data_qual_ind

      RETURN

 100  WRITE(*,*) "Error writing CRD record type 50"
      STOP 1

 1000 FORMAT ("50 ",a,1x,f6.1,1x,f7.3,1x,f7.3,1x,f6.1,1x,i1)
      END

C 60 - Compatibility Record
      SUBROUTINE write_60 (str)
      CHARACTER*512 str
      INTEGER trimlen, sci_len
      INCLUDE '../include/crd.inc'

      sci_len= trimlen(d60_sysconfig_id)
      WRITE (str,1000,err=100) d60_sysconfig_id(1:sci_len),
     &           sys_change_ind, sys_config_ind

      RETURN

 100  WRITE(*,*) "Error writing CRD record type 60"
      STOP 1

 1000 FORMAT ("60 ",a4,1x,i1,1x,i1)
      END

C 9X - User Defined Records 90-99
CC      SUBROUTINE write_9x (str)
CC      CHARACTER*512 str
CC      INCLUDE '../include/crd.inc'
CC
CC      RETURN
CC
CC      END

C 00 - Comment Record
      SUBROUTINE write_00 (str)
      CHARACTER*512 str
      INTEGER trimlen, c_len
      INCLUDE '../include/crd.inc'

      c_len= trimlen(comment)
      WRITE (str,1000,err=100) comment(1:c_len)

      RETURN

 100  WRITE(*,*) "Error writing CRD record type 00"
      STOP 1

 1000 FORMAT ("00 ",a)
      END

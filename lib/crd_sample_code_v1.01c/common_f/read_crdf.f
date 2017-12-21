C**-------------------------------------------------------------------------
C * Subroutines: read CRD data record from an input string
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
C  09/29/08   - Initialize the variable in_arg in the C0 record. (v1.00a
C               rlr)
C  03/10/09   - Record H2 Epoch Timescale corrected from I1 to I2.
C               (v1.00a rlr).
C *
C**-------------------------------------------------------------------------

C Ranging data header/footer records
C H1 - format header
      SUBROUTINE read_h1 (str)
      CHARACTER*512 str
      INCLUDE '../include/crd.inc'

      READ (str(4:512),1000,err=100) crd_literal, format_version,
     &          prod_year, prod_mon, prod_day, prod_hour
      RETURN

 100  WRITE(*,*) "Error reading CRD record type h1"
      STOP 1

 1000 FORMAT (a3,1x,i2,1x,i4,1x,i2,1x,i2,1x,i2)
      END

C H2 - station header
      SUBROUTINE read_h2 (str)
      CHARACTER*512 str
      INCLUDE '../include/crd.inc'

      READ (str(4:512),1000,err=100) stn_name,
     &          cdp_pad_id, cdp_sys_num, cdp_occ_num, stn_timescale
      RETURN

 100  WRITE(*,*) "Error reading CRD record type h2"
      STOP 1

 1000 FORMAT (a10,1x,i4,1x,i2,1x,i2,1x,i2)
      END

C H3 - spacecraft header
      SUBROUTINE read_h3 (str)
      CHARACTER*512 str
      INCLUDE '../include/crd.inc'

      READ (str(4:512),1000,err=100) target_name,
     &          ilrs_id, sic, norad, SC_timescale, target_type
      RETURN

 100  WRITE(*,*) "Error reading CRD record type h3"
      STOP 1

 1000 FORMAT (a10,1x,i8,1x,i4,1x,i8,1x,i1,1x,i1)
      END

C H4 - Session header
      SUBROUTINE read_h4 (str)
      CHARACTER*512 str
      INCLUDE '../include/crd.inc'

      READ (str(4:512),1000,err=100) 
     &          data_type, start_year, start_mon, start_day,
     &          start_hour, start_min, start_sec,
     &          end_year, end_mon, end_day, end_hour, end_min, end_sec,
     &          data_release, refraction_app_ind, CofM_app_ind,
     &          xcv_amp_app_ind, stn_sysdelay_app_ind,
     &          SC_sysdelay_app_ind, range_type_ind, data_qual_alert_ind

      RETURN

 100  WRITE(*,*) "Error reading CRD record type h4"
      STOP 1

 1000 FORMAT (i2,2(1x,i4,5(1x,i2)),1x,i2,7(1x,i1))
      END

C H8 - End of Session footer
      SUBROUTINE read_h8 (str)
      CHARACTER*512 str
      END

C H9 - End of File footer
      SUBROUTINE read_h9 (str)
      CHARACTER*512 str
      END

C Ranging data configuration records (1 of n)
C C0 - System Configuration Record
      SUBROUTINE read_c0 (str)
      character*512 str
      character*256 temp_config_ids(5)
      INCLUDE '../include/crd.inc'

      integer n_arg
      logical in_arg

C  See how many parameters are on this line.
C  BE VERY CAREFUL that str is declared length 512 in calling program!
      n_arg= 0
      in_arg= .false.
      do ii= 1, len(str)
        if (str(ii:ii) .ne. " ") then
          if (.not.in_arg) then
             in_arg= .true.
           endif
         else
           if (in_arg .eqv. .true.) then
             in_arg= .false.
             n_arg= n_arg+ 1;
           endif
         endif
       enddo

C  Create default values
      do i=1,5
        temp_config_ids(i)= ""
      enddo

C  Read the correct number of parameters
      if (n_arg .eq. 4) then
        READ (str(3:120),*,err=100)
     &          c0_detail_type, xmit_wavelength, 
     &          temp_config_ids(1)
      elseif (n_arg .eq. 5) then
        READ (str(3:120),*,err=100)
     &          c0_detail_type, xmit_wavelength, 
     &          temp_config_ids(1),temp_config_ids(2)
      elseif (n_arg .eq. 6) then
        READ (str(3:120),*,err=100)
     &          c0_detail_type, xmit_wavelength, 
     &          temp_config_ids(1),temp_config_ids(2),
     &          temp_config_ids(3)
      elseif (n_arg .eq. 7) then
        READ (str(3:120),*,err=100)
     &          c0_detail_type, xmit_wavelength, 
     &          temp_config_ids(1),temp_config_ids(2),
     &          temp_config_ids(3),temp_config_ids(4)
      elseif (n_arg .eq. 8) then
        READ (str(3:120),*,err=100)
     &          c0_detail_type, xmit_wavelength, 
     &          temp_config_ids(1),temp_config_ids(2),
     &          temp_config_ids(3),temp_config_ids(4),
     &          temp_config_ids(5)
      else 
        goto 100
      endif

C In the future, there could be more than 5 config types.
      do i=1,5
        config_ids(i)= temp_config_ids(i)(1:40)
      enddo
      return

 100  write(*,*) "Error reading CRD record type c0"
      stop 1
      end

C C1 - Laser Configuration Record
      SUBROUTINE read_c1 (str)
      character*512 str
      character*256 temp_laser_config_id
      character*256 temp_laser_type
      INCLUDE '../include/crd.inc'

      READ (str(3:512),*,err=100) 
     &            c1_detail_type, temp_laser_config_id, temp_laser_type,
     &          prim_wavelength, nom_fire_rate, pulse_energy,
     &          pulse_width, beam_div, pulses_in_semitrain
      laser_config_id= temp_laser_config_id(1:40)
      laser_type= temp_laser_type(1:40)

      RETURN

 100      write(*,*) "Error reading CRD record type c1"
      stop 1
      end

C C2 - Detector Configuration Record
      SUBROUTINE read_c2 (str)
      character*512 str
      character*256 temp_detector_config_id, temp_detector_type
      character*256 temp_output_pulse_type, temp_signal_proc

      INCLUDE '../include/crd.inc'

      READ (str(3:512),*,err=100) 
     &          c2_detail_type, temp_detector_config_id, 
     &          temp_detector_type, app_wavelength, qe, voltage, 
     &          dark_count, temp_output_pulse_type, output_pulse_width,
     &          spectral_filter, spectral_filter_xmission, 
     &          spatial_filter, temp_signal_proc
      detector_config_id= temp_detector_config_id(1:40)
      detector_type= temp_detector_type(1:40)
      output_pulse_type= temp_output_pulse_type(1:40)
      signal_proc= temp_signal_proc(1:40)

      RETURN

 100      write(*,*) "Error reading CRD record type c2"
      stop 1
      end

C C3 - Timing Configuration Record
      SUBROUTINE read_c3 (str)
      character*512 str
      character*256 temp_timing_config_id, temp_time_source
      character*256 temp_freq_source, temp_timer, temp_timer_serial_num
      INCLUDE '../include/crd.inc'

      READ (str(3:512),*,err=100) 
     &          c3_detail_type, temp_timing_config_id, temp_time_source,
     &          temp_freq_source, temp_timer, temp_timer_serial_num, 
     &          epoch_delay_corr
      timing_config_id= temp_timing_config_id(1:40)
      time_source= temp_time_source(1:40)
      freq_source= temp_freq_source(1:40)
      timer= temp_timer(1:40)
      timer_serial_num= temp_timer_serial_num(1:40)

      RETURN

 100      write(*,*) "Error reading CRD record type c3"
      stop 1
      end

C C4 - Transponder Configuration Record
      SUBROUTINE read_c4 (str)
      character*512 str
      character*256 temp_xponder_config_id
      INCLUDE '../include/crd.inc'

      READ (str(3:512),*,err=100) 
     &          c4_detail_type, temp_xponder_config_id,
     &          est_stn_utc_offset, est_stn_osc_drift,
     &          est_xponder_utc_offset, est_xponder_osc_drift,
     &          xponder_clock_ref_time, stn_off_drift_app_ind,
     &          SC_off_drift_app_ind, SC_time_simplified_ind
      xponder_config_id= temp_xponder_config_id(1:40)

      RETURN

 100      write(*,*) "Error reading CRD record type c4"
      stop 1
      end

C Ranging data records
C 10 - Range Record
      SUBROUTINE read_10 (str)
      character*512 str
      character*256 temp_sysconfig_id
      INCLUDE '../include/crd.inc'

      READ (str(3:512),*,err=100) d10_sec_of_day,d10_time_of_flight,
     &          temp_sysconfig_id, d10_epoch_event, filter_flag, 
     &          d10_detector_channel, stop_number, xcv_amp
      d10_sysconfig_id= temp_sysconfig_id(1:40)

      return

 100      write(*,*) "Error reading CRD record type 10"
C      stop 1
      end

C 11 - Normal Point Record 
      SUBROUTINE read_11 (str)
      character*512 str
      character*256 temp_sysconfig_id
      INCLUDE '../include/crd.inc'

      READ (str(3:512),*,err=100) d11_sec_of_day, d11_time_of_flight,
     &          temp_sysconfig_id, d11_epoch_event, np_window_length,
     &          num_ranges, bin_rms, bin_skew, bin_kurtosis, bin_PmM,
     &          return_rate, d11_detector_channel
      d11_sysconfig_id= temp_sysconfig_id(1:40)

      return

 100      write(*,*) "Error reading CRD record type 11"
      stop 1
      end

C 12 - Range Supplement Record
      SUBROUTINE read_12 (str)
      character*512 str
      character*256 temp_sysconfig_id
      INCLUDE '../include/crd.inc'

      READ (str(3:512),*,err=100) d12_sec_of_day, temp_sysconfig_id,
     &          refraction_corr, target_CofM_corr, nd_value, 
     &          time_bias
      d12_sysconfig_id= temp_sysconfig_id(1:40)
      return

 100      write(*,*) "Error reading CRD record type 12"
      stop 1
      end

C 20 - Meteorological Record
      SUBROUTINE read_20 (str)
      character*512 str
      INCLUDE '../include/crd.inc'

      READ (str(3:512),*,err=100) d20_sec_of_day, pressure, 
     &          temperature, humidity, value_origin
      return

 100      write(*,*) "Error reading CRD record type 20"
      stop 1
      end

C 21 - Meteorological Supplement Record
      SUBROUTINE read_21 (str)
      character*512 str
      character*256 temp_precip_type
      INCLUDE '../include/crd.inc'

      READ (str(3:512),*,err=100) d21_sec_of_day, wind_speed, 
     &          wind_direction, temp_precip_type, visibility, 
     &          sky_clarity, atmospheric_seeing, cloud_cover
      precip_type= temp_precip_type(1:40);
      return

 100      write(*,*) "Error reading CRD record type 21"
      stop 1
      end

C 30 - Pointing Angles Record
      SUBROUTINE read_30 (str)
      character*512 str
      INCLUDE '../include/crd.inc'

      READ (str(3:512),*,err=100) d30_sec_of_day, azimuth, elevation, 
     &          direction_ind, angle_origin_ind, refraction_corr_ind
      return

 100      write(*,*) "Error reading CRD record type 30"
      stop 1
      end

C 40 - Calibration Record
      SUBROUTINE read_40 (str)
      character*512 str
      character*256 temp_sysconfig_id
      INCLUDE '../include/crd.inc'

      READ (str(3:512),*,err=100) d40_sec_of_day, 
     &          type_of_data, temp_sysconfig_id,
     &          num_points_recorded, num_point_used,
     &          one_way_target_dist, cal_sys_delay, cal_delay_shift,
     &          cal_rms, cal_skew, cal_kurtosis, cal_PmM, cal_type_ind,
     &          cal_shift_type_ind, d40_detector_channel
      d40_sysconfig_id= temp_sysconfig_id(1:40)

      return

 100      write(*,*) "Error reading CRD record type 40"
      stop 1
      end

C 50 - Session Statistics Record
      SUBROUTINE read_50 (str)
      character*512 str
      character*256 temp_sysconfig_id
      INCLUDE '../include/crd.inc'

      READ (str(3:512),*,err=100)
     &          temp_sysconfig_id, sess_rms, sess_skew,
     &          sess_kurtosis, sess_PmM, data_qual_ind
      d50_sysconfig_id= temp_sysconfig_id(1:40)

      return

 100      write(*,*) "Error reading CRD record type 50"
      stop 1
      end

C 60 - Compatibility Record
      SUBROUTINE read_60 (str)
      character*512 str
      character*256 temp_sysconfig_id
      INCLUDE '../include/crd.inc'

      READ (str(3:512),*,err=100) 
     &          temp_sysconfig_id, sys_change_ind, sys_config_ind
      d60_sysconfig_id= temp_sysconfig_id(1:40)

      return

 100      write(*,*) "Error reading CRD record type 60"
      stop 1
      end

C 9X - User Defined Records 90-99
CC      SUBROUTINE read_60 (str)
CC      character*512 str
CC      end

C 00 - Comment Record
      SUBROUTINE read_00 (str)
      character*512 str
      INCLUDE '../include/crd.inc'

      READ (str(3:83),*,err=100) comment

      return

 100      write(*,*) "Error reading CRD record type 00"
      stop 1
      end

#include <stdio.h>
#include <string.h>
#include "../include/crd.h"

/*-------------------------------------------------------------------------
 * Subroutines: write CRD data records to an output file.
 *
 * Author: Randall Ricklefs / Univ of Texas / Center for Space Research
 *
 * History:
 *   July 06, 2007 - Initial version
 *   05/07/08   - Expand all configuration and data section character strings
 *                to allow up to 40 characters (plus NULL). 
 *              - Added detector channel to normalpoint (11) and calibration
 *                (40) records. 
 *              - Added field for 'crd' literal to H1.
 *              - Record '21' sky_clarity is now double rather than int.
 *   06/24/08   - Record '11', np window length is now double rather than int.
 *                (v1.00 rlr)
 *   03/10/09   - Record H2 Epoch Timescale corrected from I1 to I2.
 *                (v1.00a rlr).
 *   03/11/09   - Record H3 changed to print leading zeros rather than spaces
 *                for ilrs_id. (v1.00a rlr)
 *
**-----------------------------------------------------------------------*/

char stro[256];

/* Ranging data header/footer records */
  /* H1 - format header */
write_h1 (FILE * str_out, struct rh1 header)
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
write_h2 (FILE * str_out, struct rh2 header)
{
  sprintf (stro,
           "h2            %4d %2d %2d %2d",
           header.cdp_pad_id, header.cdp_sys_num, header.cdp_occ_num, 
           header.stn_timescale);
  strncpy (&stro[3], header.stn_name, 10);
  fprintf (str_out, "%s\n", stro);
}

  /* H3 - spacecraft header */
write_h3 (FILE * str_out, struct rh3 header)
{
  sprintf (stro,
           "h3             %07d %4d %8d %1d %1d", header.ilrs_id, header.sic, 
	   header.norad, header.SC_timescale, header.target_type);
  strncpy (&stro[3], header.target_name, 10);
  fprintf (str_out, "%s\n", stro);
}

  /* H4 - Session header */
write_h4 (FILE * str_out, struct rh4 header)
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
write_h8 (FILE * str_out)
{
  fprintf (str_out, "h8\n");
}

  /* H9 - End of File footer */
write_h9 (FILE * str_out)
{
  fprintf (str_out, "h9\n");
}

/* Ranging data configuration records (1 of n) */
    /* C0 - System Configuration Record */
write_c0 (FILE * str_out, struct rc0 config)
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
write_c1 (FILE * str_out, struct rc1 config)
{
  fprintf (str_out, 
           "c1 %d %-s %-s %.2f %.2f %.2f %.1f %.2f %d\n", 
           config.detail_type, config.laser_config_id, config.laser_type,
	   config.prim_wavelength, config.nom_fire_rate, config.pulse_energy,
	   config.pulse_width, config.beam_div, config.pulses_in_semitrain);
}

    /* C2 - Detector Configuration Record */
write_c2 (FILE * str_out, struct rc2 config)
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
write_c3 (FILE * str_out, struct rc3 config)
{
  fprintf (str_out, 
           "c3 %d %-s %-s %-s %-s %-s %.1f\n",
           config.detail_type, config.timing_config_id, config.time_source,
	   config.freq_source, config.timer, config.timer_serial_num,
	   config.epoch_delay_corr);
}

    /* C4 - Transponder Configuration Record */
write_c4 (FILE * str_out, struct rc4 config)
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
write_10 (FILE * str_out, struct rd10 data_recd)
{
  fprintf (str_out,
           "10 %.12Lf %.12Lf %-s %d %d %d %d %d\n",
           data_recd.sec_of_day, data_recd.time_of_flight, 
	   data_recd.sysconfig_id, data_recd.epoch_event, 
	   data_recd.filter_flag, data_recd.detector_channel, 
	   data_recd.stop_number, data_recd.xcv_amp);
}

    /* 11 - Normal Point Record */
write_11 (FILE * str_out, struct rd11 data_recd)
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
write_12 (FILE * str_out, struct rd12 data_recd)
{
  fprintf (str_out,
           "12 %.7Lf %-s %.1f %.4f %.2f %.4f\n",
           data_recd.sec_of_day, data_recd.sysconfig_id,
	   data_recd.refraction_corr, data_recd.target_CofM_corr, 
	   data_recd.nd_value, data_recd.time_bias);
}

    /* 20 - Meteorological Record */
write_20 (FILE * str_out, struct rd20 data_recd)
{
  fprintf (str_out,
           "20 %.3Lf %.2f %.2f %.0f %d\n",
           data_recd.sec_of_day, data_recd.pressure, data_recd.temperature, 
	   data_recd.humidity, data_recd.value_origin);
}

    /* 21 - Meteorological Supplement Record */
write_21 (FILE * str_out, struct rd21 data_recd)
{
  fprintf (str_out,
           "21 %.3Lf %.1f %.1f %-s %d %.2f %d %d\n",
           data_recd.sec_of_day, data_recd.wind_speed, data_recd.wind_direction, 
	   data_recd.precip_type, data_recd.visibility, data_recd.sky_clarity,
	   data_recd.atmospheric_seeing, data_recd.cloud_cover);
}

    /* 30 - Pointing Angles Record */
write_30 (FILE * str_out, struct rd30 data_recd)
{
  fprintf (str_out,
           "30 %.3Lf %.4f %.4f %d %d %d\n",
           data_recd.sec_of_day, data_recd.azimuth, data_recd.elevation, 
	   data_recd.direction_ind, data_recd.angle_origin_ind, 
	   data_recd.refraction_corr_ind);
}

    /* 40 - Calibration Record */
write_40 (FILE * str_out, struct rd40 data_recd)
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
write_50 (FILE * str_out, struct rd50 data_recd)
{
  fprintf (str_out,
           "50 %-s %.1f %.3f %.3f %.1f %d\n",
           data_recd.sysconfig_id, data_recd.sess_rms, data_recd.sess_skew, 
	   data_recd.sess_kurtosis, data_recd.sess_PmM, 
	   data_recd.data_qual_ind);
}

    /* 60 - Compatibility Record */
write_60 (FILE * str_out, struct rd60 data_recd)
{
  fprintf (str_out, 
	 "60 %-s %d %d\n", 
	 data_recd.sysconfig_id, data_recd.sys_change_ind, 
         data_recd.sys_config_ind); 
}

    /* 9X - User Defined Record */
write_9x (FILE * str_out, struct rd9x data_recd)
{
}

    /* 00 - Comment Record */
write_00 (FILE * str_out, struct rd00 data_recd)
{
  fprintf (str_out,
	   "00 %-s\n", 
	   data_recd.comment);
}

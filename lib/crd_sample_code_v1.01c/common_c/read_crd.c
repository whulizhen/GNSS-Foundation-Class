#include <stdio.h>
#include <string.h>
#include "../include/crd.h"

/*-------------------------------------------------------------------------
 * Subroutines: read CRD data records from an input string.
 *
 * Author: Randall Ricklefs / Univ of Texas / Center for Space Research
 *
 * History:
 *   July 06, 2007 - Initial version
 *   June 24, 2008 - v1.00
 *
**-----------------------------------------------------------------------*/

char stro[256];

/* Ranging data header/footer records */
  /* H1 - format header */
read_h1 (char * str, struct rh1 *header)
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
read_h2 (char * str, struct rh2 *header)
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
read_h3 (char * str, struct rh3 *header)
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
read_h4 (char * str, struct rh4 *header)
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
read_h8 (char * str)
{
  sscanf (str, "H8");
}

  /* H9 - End of File footer */
read_h9 (char * str)
{
  sscanf (str, "H9");
}

/* Ranging data configuration records (1 of n) */
    /* C0 - System Configuration Record */
read_c0 (char * str, struct rc0 *config)
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
read_c1 (char * str, struct rc1 *config)
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
read_c2 (char * str, struct rc2 *config)
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
read_c3 (char * str, struct rc3 *config)
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
read_c4 (char * str, struct rc4 *config)
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
read_10 (char * str, struct rd10 *data_recd)
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
read_11 (char * str, struct rd11 *data_recd)
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
read_12 (char * str, struct rd12 *data_recd)
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
read_20 (char * str, struct rd20 *data_recd)
{
  int nstat;

  nstat= sscanf (str,
           "%*s %Lf %lf %lf %lf %d",
           &data_recd->sec_of_day, &data_recd->pressure, &data_recd->temperature, 
	   &data_recd->humidity, &data_recd->value_origin);
  return (nstat+1);
}

    /* 21 - Meteorological Supplement Record */
read_21 (char * str, struct rd21 *data_recd)
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
read_30 (char * str, struct rd30 *data_recd)
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
read_40 (char * str, struct rd40 *data_recd)
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
read_50 (char * str, struct rd50 *data_recd)
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
read_60 (char * str, struct rd60 *data_recd)
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
read_9x (char * str, struct rd9x *data_recd)
{
  int nstat= 0;
  return (nstat);
}

    /* 00 - Comment Record */
read_00 (char * str, struct rd00 *data_recd)
{
  int nstat= 1;
  strncpy(str,data_recd->comment,80);
  data_recd->comment[80]= '\0';
  
  return (nstat);
}

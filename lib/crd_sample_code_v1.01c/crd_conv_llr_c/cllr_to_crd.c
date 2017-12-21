#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#define FNLEN   256
#define FMT_VERSION 0

#include "crd.h"
#include "cospar_llr.h"

FILE *str_in, *str_out_ld;
fpos_t startpos;
fpos_t startpos2;
struct rllr0 llr0;
struct rllr1 llr1;
struct rllr2 llr2;
struct rllr3 llr3;
struct rllr4 llr4;
struct rllr5 llr5;
struct rllr6 llr6;
struct rllr7 llr7;
struct rh1 h1;
struct rh2 h2;
struct rh3 h3;
struct rh4 h4;
struct rc0 c0;
struct rc1 c1;
struct rc2 c2;
struct rc3 c3;
struct rc4 c4;
struct rd10 d10;
struct rd11 d11;
struct rd12 d12;
struct rd20 d20;
struct rd21 d21;
struct rd30 d30;
struct rd40 d40;
struct rd50 d50;
struct rd60 d60;
struct rd00 d00;

void get_sat_ids ();
void get_station_id ();
void get_sys_config ();
double mlrs_geo_calc();
void move_cllr_to_crd_hdr ();
void move_cllr5_to_crd_hdr ();
void move_cllr_to_crd_data ();
void setup_files ();
void sodtohms ();
void write_headers ();
void write_data ();
void write_end_of_data_block ();
void write_end_of_file ();

/*-------------------------------------------------------------------------
 * Program: cospar_llr_to_crd
 *
 * Purpose:
 * Converts old cospar LLR data to CRD format.
 *
 * Calling sequence:
 *   cospar_llr_to_crd -i cospar_filename -o crd_filename
 *
 * Note: The following cospar fields are not present in the CRD and could
 * be added in a new #61 lunar compatibility record. This will not be done at
 * at the current time, and may not be important enough to pursue. 
 *	Dark count,
 *	Lunar site count,
 *	Delay time bsae, and
 *	Shot-by-shot resolution.
 * The laser energy and laser pusle length fields are in the CRD but are not
 * yet written by this program. rlr.
 *
 * Author: Randall Ricklefs / Univ of Texas / Center for Space Research
 *
 * History:
 *   Nov 2, 2007 - Initial version
 *   Feb 2, 2010 - Because of round-off problems, change 1.e7 to
 *   		   (long double)1.e7 in second of days calulations and
 *   		   1.e13 to (long double)1.e13 in time of flight calcutions.
 *   		   1.01b. rlr.
 *
**-----------------------------------------------------------------------*/
int
main (argc, argv)
     int argc;
     char *argv[];
{
  int data_release = 0;
  int delta;
  int done = 0;
  int idummy;
  int first_pass = 1, p_first_pass = 1;
  int first_record = 0;
  int header_1_waiting= 0;
  int header_2_waiting= 0;
  int n_fr = 0;
  int new_pad= 0;
  int new_pass= 0;
  int new_target= 0;
  int next_new_pad = 1;
  int next_new_target = 1;
  int o_refl= -1, o_obsy= -1;
  int po_refl= -1, po_obsy= -1;
  int read_1= 0;
  int read_2= 0;
  int read_3 = 0;
  int read_4 = 0;
  int read_5 = 0;
  int reset_def = 1;
  int status;
  long filpos;
  long lclock;
  char istrdum[11];
  char str[256];
  double dttm, o_dttm= 0, po_dttm= 0;
  double ddelta;
  double jd, jdf;
  double esec_of_day;
  double ssec_of_day;
  struct tm gmt;

/*  Greetings */
  printf("Convert COSPAR Lunar Normalpoints to CRD Format: version 1.01b 02/02/2012\n");

/* get file names and open files  */
  setup_files (argc, argv);

/* get processing time */
  time (&lclock);
  gmt = *gmtime (&lclock);

/* Main data loop */
  while (done == 0)
    {
      fgetpos (str_in, &startpos);

    /* Phase 1 - read start and stop time of each data block;
     *	     ends when encountering new station, new reflector, or large
     *	     time jump. */
      while ((status = fgets (str, 256, str_in)) != NULL)
	{
	  /*printf("Prelim: [%s]\n",str);*/
	  if (strlen (str) <= 1)
	    continue;		/* Blank lines... */
       /* Run header */
	  else if (str[0] == '1')
	    {
	      if (read_1) break;
	      read_1= 1;
	    }
       /* Detail record */
	  else if (str[0] == '3')
	    {
	      read_llr3 (str, &llr3);
	      if (!read_3)
		{
		  read_3 = 1;
		  h4.start_year = llr3.year;
		  h4.start_mon = llr3.month;
		  h4.start_day = llr3.day;
		  h4.start_hour = llr3.hour;
		  h4.start_min = llr3.min;
		  h4.start_sec = llr3.sec/1.e7;

		  h4.end_year = llr3.year;
		  h4.end_mon = llr3.month;
		  h4.end_day = llr3.day;
		  h4.end_hour = llr3.hour;
		  h4.end_min = llr3.min;
		  h4.end_sec = llr3.sec/1.e7;
		}
	      else
		{
		  h4.end_year = llr3.year;
		  h4.end_mon = llr3.month;
		  h4.end_day = llr3.day;
		  h4.end_hour = llr3.hour;
		  h4.end_min = llr3.min;
		  h4.end_sec = llr3.sec/1.e7;
		}
	    }
       /* Normalpoint record */
	  else if (str[0] == '4' && !read_3)
	    {
	      read_llr4 (str, &llr4);
	      if (!read_4)
		{
		  read_4 = 1;
		  h4.start_year = llr4.year;
		  h4.start_mon = llr4.month;
		  h4.start_day = llr4.day;
		  h4.start_hour = llr4.hour;
		  h4.start_min = llr4.min;
		  h4.start_sec = llr4.sec/1.e7;

		  h4.end_year = llr4.year;
		  h4.end_mon = llr4.month;
		  h4.end_day = llr4.day;
		  h4.end_hour = llr4.hour;
		  h4.end_min = llr4.min;
		  h4.end_sec = llr4.sec/1.e7;
		}
	      else
		{
		  h4.end_year = llr4.year;
		  h4.end_mon = llr4.month;
		  h4.end_day = llr4.day;
		  h4.end_hour = llr4.hour;
		  h4.end_min = llr4.min;
		  h4.end_sec = llr4.sec/1.e7;
		}
	    }
       /* Mini-normalpoint record */
	  else if (str[0] == '5' && !read_3)
	    {
	      read_llr5 (str, &llr5);
	      if (!read_5)
		{
		  h4.start_year = llr5.year;
		  h4.start_mon = llr5.month;
		  h4.start_day = llr5.day;
		  h4.start_hour = llr5.hour;
		  h4.start_min = llr5.min;
		  h4.start_sec = llr5.sec/1.e7;

		  h4.end_year = llr5.year;
		  h4.end_mon = llr5.month;
		  h4.end_day = llr5.day;
		  h4.end_hour = llr5.hour;
		  h4.end_min = llr5.min;
		  h4.end_sec = llr5.sec/1.e7;
		}
	      else
		{
		  h4.end_year = llr5.year;
		  h4.end_mon = llr5.month;
		  h4.end_day = llr5.day;
		  h4.end_hour = llr5.hour;
		  h4.end_min = llr5.min;
		  h4.end_sec = llr5.sec/1.e7;
		}
	      if (po_refl != llr5.reflector) new_target= 1;
	      if (po_obsy != llr5.observatory_code) new_pad= 1;
	      grtojd (llr5.year, llr5.month, llr5.day, llr5.hour, llr5.min,
		0.e0, &jd, &jdf);
	      dttm= jd- 2440000.5+ jdf;
	      new_pass= fabs(dttm- o_dttm) > 1.0/24.0;
	      if (!p_first_pass &&
			(!read_5 || new_pass || new_pad || new_target))
		{
		  break;
                }
	      po_refl= llr5.reflector;
	      po_obsy= llr5.observatory_code;
	      po_dttm= dttm;
	      p_first_pass= 0;
	      read_5 = 1;
	    }
	}
      fsetpos (str_in, &startpos);
      read_1= 0;
      read_2= 0;
      read_3= 0;
      read_4= 0;
      read_5= 0;

/**printf("h4: %d %d %d %d %d %d   %d %d %d %d %d %d\n",h4.start_year,h4.start_mon,h4.start_day,h4.start_hour,h4.start_min,h4.start_sec,h4.end_year,h4.end_mon,h4.end_day,h4.end_hour,h4.end_min,h4.end_sec);**/

    /* Phase 2 - read and process data in block found in phase 1;
     *	     ends when encountering new station, new reflector, or large
     *	     time jump. */
      while ((status = fgets (str, 256, str_in)) != NULL)
	{
	  if (strlen (str) <= 1)
	    continue;		/* Blank lines... */
       /* System configuration record */
	  if (str[0] == '0')
	    {
	      read_llr0 (str, &llr0);
	    }
       /* Run header */
	  else if (str[0] == '1')
	    {
	      read_llr1 (str, &llr1);
	      if (o_refl != llr1.reflector) new_target= 1;
	      if (o_obsy != llr1.observatory_code) new_pad= 1;
	      if (read_1) 
                {
                  fsetpos (str_in, &startpos2);
                  break;
                }
	      read_1= 1;
	      header_1_waiting= 1;
	      o_refl= llr1.reflector;
	      o_obsy= llr1.observatory_code;
	    }
       /* Run sub-header */
	  else if (str[0] == '2')
	    {
	      read_2= 1;
	      header_2_waiting= 1;
	      read_llr2 (str, &llr2);
	    }
       /* Detail (range return) record */
	  else if (str[0] == '3')
	    {
	      read_llr3 (str, &llr3);
	      if (header_1_waiting)
		{
		  h4.data_type= 0;
	          move_cllr_to_crd_hdr (gmt);
		  write_headers (str_out_ld, first_pass, new_pad, new_target);
	          first_pass= 0;
	          new_target= 0;
	          new_pad= 0;
		  header_1_waiting= 0;
		}
	      move_cllr_to_crd_data (3);
	      write_data (str_out_ld, reset_def);
	      reset_def= 0;
	    }
       /* Normalpoint record */
	  else if (str[0] == '4')
	    {
	      read_llr4 (str, &llr4);
	      if (header_2_waiting)
		{
                  write_h8 (str_out_ld);
		  h4.data_type= 1;
	          move_cllr_to_crd_hdr (gmt);
		  write_headers (str_out_ld, first_pass, new_pad, new_target);
	          first_pass= 0;
	          new_target= 0;
	          new_pad= 0;
		  header_2_waiting= 0;
		  reset_def= 1;
		}
	      move_cllr_to_crd_data (4);
	      write_data (str_out_ld, reset_def);
	      reset_def= 0;
	    }
       /* Mini-normalpoint record */
	  else if (str[0] == '5')
	    {
	      read_llr5 (str, &llr5);
	      if (o_refl != llr5.reflector) new_target= 1;
	      if (o_obsy != llr5.observatory_code) new_pad= 1;
	      grtojd (llr5.year, llr5.month, llr5.day, llr5.hour, llr5.min,
		0.e0, &jd, &jdf);
	      dttm= jd- 2440000.5+ jdf;
	      new_pass= fabs(dttm- o_dttm) > 1.0/24.0;
	      if (!read_5 || new_pass || new_pad || new_target)
		{
		  h4.data_type= 1;
	          move_cllr5_to_crd_hdr (gmt);
		  write_headers (str_out_ld, first_pass, new_pad, new_target);
		  read_5= 1;
	          first_pass= 0;
	          new_target= 0;
	          new_pad= 0;
		  new_pass= 0;
		}
              d11.np_window_length = llr5.time_span;
	      move_cllr_to_crd_data (5);
	      write_data (str_out_ld, reset_def);
	      reset_def= 0;
	      o_refl= llr5.reflector;
	      o_obsy= llr5.observatory_code;
	      o_dttm= dttm;
	      if (read_5) break;
	    }
       /* Calibration record */
	  else if (str[0] == '6')
	    {
	      read_llr6 (str, &llr6);
	      /* ignore this line? */
	    }
       /* Comment record */
	  else if (str[0] == '7')
	    {
	      read_llr7 (str, &llr7);
	      sprintf(d00.comment, "%d %4d %2d %2d %2d %2d %2d ",
		llr7.laser_color, llr7.year, llr7.month, llr7.day,
		llr7.hour, llr7.min, llr7.source);
	      strncat(d00.comment, llr7.comment,56);
	      d00.comment[80]= '\0';
	      write_00 (str_out_ld, d00);
	    }
          fgetpos (str_in, &startpos2);
	}
      if (status == NULL) done= 1;
      read_1= 0;
      read_2= 0;
      read_5= 0;
      reset_def= 1;
      write_end_of_data_block (str_out_ld, h4.data_type);
    }

  /* clean up and go home... */
  write_end_of_file (str_out_ld);
  fclose (str_in);
  fclose (str_out_ld);
}

/*-------------------------------------------------------------------------
**
**      move_cllr_to_crd_hdr - Copy info from the input header records
**                        and elsewhere to the CRD header records.
**
**-----------------------------------------------------------------------*/
void
move_cllr_to_crd_hdr (struct tm gmt)
{
  char mon_name[4];
  int mon, day;
  int idummy, i;
  int cdp_pad_id= 0;

/* Format header */
  strcpy (h1.crd_literal,"CRD");
  h1.format_version = FMT_VERSION;
  h1.prod_year = gmt.tm_year + 1900;
  h1.prod_mon = gmt.tm_mon + 1;
  h1.prod_day = gmt.tm_mday;
  h1.prod_hour = gmt.tm_hour;

/* Station header */
  if (llr1.observatory_code == 71110)
    cdp_pad_id = 9999;
  else if (llr1.observatory_code == 71111)
    cdp_pad_id = 7086;
  else if (llr1.observatory_code == 71112)
    cdp_pad_id = 7080;
  get_station_id (cdp_pad_id, &h2.stn_name);
  h2.cdp_pad_id = cdp_pad_id;
  h2.cdp_sys_num = -1;
  h2.cdp_occ_num = -1;
  /* llr1.epoch_time_base; This is just UTC vs UT0, etc. Do this instead: */ 
  if (h4.start_year < 1985) /* What year should this be?? */
    h2.stn_timescale = 3;
  else 
    h2.stn_timescale = 4;

/* Target header */
  h3.ilrs_id = llr1.reflector + 100;
  get_sat_ids (0, &h3.ilrs_id, &h3.norad, &h3.sic, &h3.target_name, &idummy);
  h3.SC_timescale = 0;
  h3.target_type = 2;		/* lunar */
  if (llr1.data_file_version == ' ') h4.data_release= 0;
  /* from np detail recds; none for sed */
  else h4.data_release = tolower(llr1.data_file_version)- 'a';

  h4.refraction_app_ind = 0;
  h4.CofM_app_ind = 0;
  h4.xcv_amp_app_ind = 0;
  h4.stn_sysdelay_app_ind = 1;
  h4.SC_sysdelay_app_ind = 0;
  h4.range_type_ind = 2;

/* Configuration header */
  c0.detail_type = 0;
  c0.xmit_wavelength = llr1.laser_wavelength/10.;
  strcpy (c0.config_ids[0], "std");
  for (i = 1; i < 10; i++)
    {
      c0.config_ids[i][0] = '\0';
    }

/* Normalpoint record */
  if (h4.data_type == 1)
    d11.np_window_length = llr4.time_span;

/* Calibration record */
  d40.type_of_data = 0;
  strcpy (d40.sysconfig_id, c0.config_ids[0]);
  d40.num_points_recorded = -1;
  d40.num_points_used = -1;
  d40.one_way_target_dist = -1;
  d40.cal_delay_shift = 0;
  d40.cal_skew = -1;
  d40.cal_kurtosis = -1;
  d40.cal_PmM = -1;
  if (llr1.observatory_code == 71110)	/* Mcd 2.7m */
    {
      d40.cal_shift_type_ind = 3;
      d40.cal_type_ind = 0;
    }
  else if (llr1.observatory_code / 10 == 7111)	/* MLRS */
    {
      d40.cal_shift_type_ind = 3;
      d40.cal_type_ind = 3;
    }
  else
    {
      d40.cal_shift_type_ind = -1;
      d40.cal_type_ind = -1;
    }

/* Session (pass) statistical record */
  strcpy (d50.sysconfig_id, c0.config_ids[0]);
  d50.sess_skew = -1;
  d50.sess_kurtosis = -1;
  d50.sess_PmM = -1;
  if (llr1.data_quality_ind == ' ') d50.data_qual_ind= 0;
  else d50.data_qual_ind = tolower(llr1.data_quality_ind)- 'a'+ 1;

/* Compatibility record */
  strcpy (d60.sysconfig_id, c0.config_ids[0]);
  get_sys_config (h2.cdp_pad_id,
		  h4.start_year, h4.start_mon, h4.start_year, h4.start_hour,
		  &d60.sys_change_ind, &d60.sys_config_ind);
}

/*-------------------------------------------------------------------------
**
**      move_cllr5_to_crd_hdr - Copy info from the llr #5 (mini-normalpoint)
**                        records and elsewhere to the CRD header records.
**
**-----------------------------------------------------------------------*/
void
move_cllr5_to_crd_hdr (struct tm gmt)
{
  char mon_name[4];
  int mon, day;
  int idummy, i;
  int cdp_pad_id= 0;

/* Format header */
  h1.format_version = FMT_VERSION;
  h1.prod_year = gmt.tm_year + 1900;
  h1.prod_mon = gmt.tm_mon + 1;
  h1.prod_day = gmt.tm_mday;
  h1.prod_hour = gmt.tm_hour;

/* Station header */
  if (llr5.observatory_code == 71110)
    cdp_pad_id = 9999;
  else if (llr5.observatory_code == 71111)
    cdp_pad_id = 7086;
  else if (llr5.observatory_code == 71112)
    cdp_pad_id = 7080;
  get_station_id (cdp_pad_id, &h2.stn_name);
  h2.cdp_pad_id = cdp_pad_id;
  h2.cdp_sys_num = -1;
  h2.cdp_occ_num = -1;
  if (h4.start_year < 1985) /* What year should this be?? */
    h2.stn_timescale = 3;
  else 
    h2.stn_timescale = 4;

/* Target header */
  h3.ilrs_id = llr5.reflector + 100;
  get_sat_ids (0, &h3.ilrs_id, &h3.norad, &h3.sic, &h3.target_name, &idummy);
  h3.SC_timescale = 0;
  h3.target_type = 2;		/* lunar */
  if (llr5.data_file_version == ' ') h4.data_release= 0;
  /* from np detail recds; none for sed */
  else h4.data_release = tolower(llr5.data_file_version)- 'a';

  h4.refraction_app_ind = 0;
  h4.CofM_app_ind = 0;
  h4.xcv_amp_app_ind = 0;
  h4.stn_sysdelay_app_ind = 1;
  h4.SC_sysdelay_app_ind = 0;
  h4.range_type_ind = 2;

/* Configuration header */
  c0.detail_type = 0;
  c0.xmit_wavelength = llr5.laser_wavelength/10.;
  strcpy (c0.config_ids[0], "std");
  for (i = 1; i < 10; i++)
    {
      c0.config_ids[i][0] = '\0';
    }

/* Calibration record */
  d40.type_of_data = 0;
  strcpy (d40.sysconfig_id, c0.config_ids[0]);
  d40.num_points_recorded = -1;
  d40.num_points_used = -1;
  d40.one_way_target_dist = -1;
  d40.cal_sys_delay = 0;
  d40.cal_delay_shift = 0;
  d40.cal_rms = llr5.uncert_estimate/10.;
  d40.cal_skew = -1;
  d40.cal_kurtosis = -1;
  d40.cal_PmM = -1;
  if (llr5.observatory_code == 71110)	/* Mcd 2.7m */
    {
      d40.cal_shift_type_ind = -1;
      d40.cal_type_ind = -1;
    }
  else if (llr5.observatory_code / 10 == 7111)	/* MLRS */
    {
      d40.cal_shift_type_ind = 3;
      d40.cal_type_ind = 3;
    }
  else
    {
      d40.cal_shift_type_ind = -1;
      d40.cal_type_ind = -1;
    }

/* Session (pass) statistical record */
  strcpy (d50.sysconfig_id, c0.config_ids[0]);
  d50.sess_rms = llr5.uncert_estimate/10.;
  d50.sess_skew = -1;
  d50.sess_kurtosis = -1;
  d50.sess_PmM = -1;
  if (llr5.data_quality_ind == ' ') d50.data_qual_ind= 0;
  else d50.data_qual_ind = tolower(llr5.data_quality_ind)- 'a'+ 1;

/* Compatibility record */
  strcpy (d60.sysconfig_id, c0.config_ids[0]);
  get_sys_config (h2.cdp_pad_id,
		  h4.start_year, h4.start_mon, h4.start_year, h4.start_hour,
		  &d60.sys_change_ind, &d60.sys_config_ind);
}

/*-------------------------------------------------------------------------
**
**      move_cllr_to_crd_data - Copy info from the input llr records
**                       and elsewhere to the CRD data records.
**
**-----------------------------------------------------------------------*/
void
move_cllr_to_crd_data (int recd_type)
{
  long double sec_of_day;
  double geo_corr;
  char wind_dir[2];

/* Range Record */
  if (recd_type == 3)
    {
      sec_of_day= d10.sec_of_day = 
	llr3.hour * 3600 + llr3.min * 60 + (llr3.sec+ llr2.clock_offset)/ 
          (long double)1.e7;
      d10.time_of_flight = (llr3.time_of_flight- 
		llr3.electronic_delay- llr3.geometric_delay)/(long double)1.e13;
      strcpy (d10.sysconfig_id, c0.config_ids[0]);
      d10.epoch_event = 2;
      d10.filter_flag = llr3.filter_flag+1;
      d10.detector_channel = 0;
      d10.stop_number = 0;
      d10.xcv_amp = -1;
      d12.sec_of_day = sec_of_day;

  /* Calibration Record */
      d40.sec_of_day = sec_of_day;
      /* Geometric delay is wrapped into range and not included here, since
       * it varies with point angle. */
      d40.cal_sys_delay = llr3.electronic_delay/10.;
      d40.cal_rms = llr3.uncert_estimate/10.;
      d50.sess_rms = llr3.uncert_estimate/10.;

  /* Pointing Angle Record */
      d30.sec_of_day = sec_of_day;
      d30.direction_ind = 2;
      d30.angle_origin_ind = 2;
      d30.refraction_corr_ind = 1;	/* Assumption! */
      d30.azimuth = llr3.azimuth/1.e7;
      d30.elevation = llr3.elevation/1.e7;

  /* Meteorology Record */
      d20.sec_of_day = sec_of_day;
      d20.pressure = llr3.pressure/100.;
      d20.temperature = llr2.temperature/10.0+ 273.15;	/* deg C to deg K */
      d20.humidity = llr2.humidity;
      d21.sec_of_day = sec_of_day;
      d21.wind_speed= llr2.wind_speed*1000./3600.;	/* km/hr->m/sec */
      wind_dir[0]= tolower(llr2.wind_direction[0]);
      wind_dir[1]= tolower(llr2.wind_direction[1]);
      if (strncmp(wind_dir,"n ",2) == 0) d21.wind_direction= 0;
      else if (strncmp(wind_dir,"ne",2) == 0) d21.wind_direction= 45;
      else if (strncmp(wind_dir,"e ",2) == 0) d21.wind_direction= 90;
      else if (strncmp(wind_dir,"se",2) == 0) d21.wind_direction= 135;
      else if (strncmp(wind_dir,"s ",2) == 0) d21.wind_direction= 180;
      else if (strncmp(wind_dir,"sw",2) == 0) d21.wind_direction= 225;
      else if (strncmp(wind_dir,"w ",2) == 0) d21.wind_direction= 270;
      else if (strncmp(wind_dir,"nw",2) == 0) d21.wind_direction= 315;
      strcpy (d21.precip_type, "na");
      d21.visibility= -1;
      d21.sky_clarity= -1;
      d21.atmospheric_seeing= llr2.seeing/10.;
      d21.cloud_cover= -1;
    }
/* Normalpoint Record */
  else if (recd_type == 4)
    {
      sec_of_day = d11.sec_of_day =
	llr4.hour * 3600 + llr4.min * 60 + (llr4.sec+ llr2.clock_offset) / 
          (long double)1.e7;
      d11.time_of_flight = (llr4.time_of_flight- 
		llr4.electronic_delay- llr4.geometric_delay)/(long double)1.e13;
      strcpy (d11.sysconfig_id, c0.config_ids[0]);
      d11.epoch_event = 2;
      d11.num_ranges = llr4.number_of_returns;
      d11.bin_rms = llr4.uncert_estimate/10.;
      d40.cal_rms = llr4.uncert_estimate/10.;
      d11.bin_skew = -1;
      d11.bin_kurtosis = -1;
      d11.bin_PmM = -1;
      d11.return_rate = llr4.signal_to_noise/10.;
      d11.detector_channel= 0;
      d12.sec_of_day = sec_of_day;

  /* Calibration Record */
      d40.sec_of_day = sec_of_day;
      /* Geometric delay is wrapped into range and not included here, since
       * it varies with point angle. */
      d40.cal_sys_delay = llr4.electronic_delay/10.;
      d50.sess_rms = llr4.uncert_estimate/10.;

  /* Meteorology Record */
      d20.sec_of_day = sec_of_day;
      d20.pressure = llr4.pressure/100.;
      d20.temperature = llr2.temperature/10.0+ 273.15;	/* deg C to deg K */
      d20.humidity = llr2.humidity;
      d21.sec_of_day = sec_of_day;
    }
/* Mini-normalpoint Record */
  else if (recd_type == 5)
    {
      sec_of_day = d11.sec_of_day =
	llr5.hour * 3600 + llr5.min * 60 + llr5.sec / (long double)1.e7;
      d11.time_of_flight = llr5.time_of_flight/(long double)1.e13;
      strcpy (d11.sysconfig_id, c0.config_ids[0]);
      d11.epoch_event = 2;
      d11.num_ranges = llr5.number_of_returns;
      d11.bin_rms = llr5.uncert_estimate/10.;
      d11.bin_skew = -1;
      d11.bin_kurtosis = -1;
      d11.bin_PmM = -1;
      d11.return_rate = llr5.signal_to_noise/10.;
      d12.sec_of_day = sec_of_day;

  /* Calibration Record */
      d40.sec_of_day = sec_of_day;
      d40.cal_sys_delay = 0;
      d40.cal_rms = llr5.uncert_estimate/10.;
      d40.detector_channel= 0;
      d50.sess_rms = llr5.uncert_estimate/10.;

  /* Meteorology Record */
      d20.sec_of_day = sec_of_day;
      d20.pressure = llr5.pressure/100.;
      d20.temperature = llr5.temperature/10.0+ 273.15;	/* deg C to deg K */
      d20.humidity = llr5.humidity;
      d21.sec_of_day = sec_of_day;
    }

/* Range Supplement Record */
  strcpy (d12.sysconfig_id, c0.config_ids[0]);
  d12.refraction_corr = 0;
  d12.target_CofM_corr = 0;
  d12.nd_value = -1;
  d12.time_bias = -1;
}

/*-------------------------------------------------------------------------
**
**      sodtohms - Convert seconds of day to hour/minute/second
**
**-----------------------------------------------------------------------*/
void
sodtohms (long sod, int *h, int *m, int *s)
{
  *h = sod / 3600;
  *m = (sod - *h * 3600) / 60;
  *s = sod - *h * 3600 - *m * 60;
}

/*-------------------------------------------------------------------------
**
**      write_headers - Write CRD header records
**
**-----------------------------------------------------------------------*/
void
write_headers (FILE * str_out_crd, int first_pass, int new_pad,
	       int new_target)
{
  if (first_pass)
    write_h1 (str_out_crd, h1);
  if (first_pass || new_pad)
    write_h2 (str_out_crd, h2);
  if (first_pass || new_pad || new_target)
    write_h3 (str_out_crd, h3);
  write_h4 (str_out_crd, h4);
  write_c0 (str_out_crd, c0);
  write_60 (str_out_crd, d60);
}

/*-------------------------------------------------------------------------
**
**      write_data - Write CRD Data records
**
**-----------------------------------------------------------------------*/
void
write_data (FILE * str_out_crd, int reset_def)
{
  static double opressure = -999;
  static double otemperature = -999;
  static double ohumidity = -999;
  static double orefraction_corr = -999;
  static double otarget_CofM_corr = -999;
  static double ond_value = -999;
  static double otime_bias = -999;
  static double oazimuth = -999;
  static double oelevation = -999;
  static double ocal_delay_shift = -999;
  static int owind_speed= -1;
  static int owind_dir= -1;
  static int ovisibility= -1;
  static int osky_clarity= -1;
  static int oseeing= -1;
  static int ocloud_cover= -1;

  /* Starting new block, so reset all reference values */
  if (reset_def)
    {
      opressure = otemperature = ohumidity = -999;
      owind_speed= owind_dir= ovisibility= osky_clarity= oseeing= 
	ocloud_cover= -1;
      orefraction_corr = otarget_CofM_corr = ond_value = otime_bias = -999;
      oazimuth = oelevation = ocal_delay_shift = -999;
    }

  /* Normalpoint or range record? */
  if (h4.data_type != 1)
    write_10 (str_out_crd, d10);
  else
    write_11 (str_out_crd, d11);

  /* Met record */
  /* Write '20' first to keep shorter sod in time sequence */
  if (fabs (d20.pressure - opressure) > 0.009 || 
	fabs (d20.temperature - otemperature) > 0.09 || 
	fabs (d20.humidity - ohumidity) > 0.9)	
	/* For existing real data! Normally '5' */
    {
      write_20 (str_out_crd, d20);
      opressure = d20.pressure;
      otemperature = d20.temperature;
      ohumidity = d20.humidity;
    }

  /* Met supplement record */
  if (fabs (d21.wind_speed - owind_speed) > 1 || 
	fabs (d21.wind_direction - owind_dir) > 1 || 
	fabs (d21.visibility - ovisibility) > 1 || 
	fabs (d21.sky_clarity - osky_clarity) > 1 || 
	fabs (d21.atmospheric_seeing - oseeing) > 1 || 
	fabs (d21.cloud_cover - ocloud_cover) > 1)
    {
      write_21 (str_out_crd, d21);
      owind_speed= d21.wind_speed;
      owind_dir= d21.wind_direction;
      ovisibility= d21.visibility;
      osky_clarity= d21.sky_clarity;
      oseeing= d21.atmospheric_seeing;
      ocloud_cover= d21.cloud_cover;
    }

  /* Pointing angle record */
  if ((fabs (d30.azimuth - oazimuth) > 0.1 ||
       fabs (d30.elevation - oelevation) > 0.1) && h4.data_type != 1)
    {
      write_30 (str_out_crd, d30);
      oazimuth = d30.azimuth;
      oelevation = d30.elevation;
    }

  /* Range supplement record */
  if ((fabs (d12.refraction_corr - orefraction_corr) > 1 ||
       fabs (d12.target_CofM_corr - otarget_CofM_corr) > 1 ||
       fabs (d12.nd_value - ond_value) > 1 ||
       fabs (d12.time_bias - otime_bias) > 1) && h4.data_type != 1)
    {
      write_12 (str_out_crd, d12);
      orefraction_corr = d12.refraction_corr;
      otarget_CofM_corr = d12.target_CofM_corr;
      ond_value = d12.nd_value;
      otime_bias = d12.time_bias;
    }

  /* Calibration record */
  if (fabs (d40.cal_delay_shift - ocal_delay_shift) > 1)
    {
      write_40 (str_out_crd, d40);
      ocal_delay_shift = d40.cal_delay_shift;
    }
}

/*-------------------------------------------------------------------------
**
**      write_end_of_data_block - Write CRD End-of-data-block record
**
**-----------------------------------------------------------------------*/
void
write_end_of_data_block (FILE * str_out_crd, int data_type)
{
  if (data_type == 1)
    write_50 (str_out_crd, d50);
  write_h8 (str_out_crd);
}

/*-------------------------------------------------------------------------
**
**      write_end_of_file - Write CRD End-of-file record
**
**-----------------------------------------------------------------------*/
void
write_end_of_file (FILE * str_out_crd)
{
  write_h9 (str_out_crd);
}

/*-------------------------------------------------------------------------
**
**	setup_files		- Open input and output file names
**
**-----------------------------------------------------------------------*/
void
setup_files (int argc, char *argv[])
{
  char cospar_in[FNLEN] = { "" };
  char crd_out[FNLEN] = { "" };

  int found = 0, i;

/*  Get file names  */
  for (i = 1; i < argc; i += 2)
    {
      if (argv[i][0] == '-')
	{
	  if (argv[i][1] == 'i')
	    {
	      strcpy (cospar_in, argv[i + 1]);
	      printf ("cospar_in = [%s]\n", cospar_in);
	      found += 1;
	    }
	  else if (argv[i][1] == 'o')
	    {
	      if (found == 1)
		{
		  strcpy (crd_out, argv[i + 1]);
		  /*strcat (crd_out, ".frd"); */
		}
	      printf ("crd_out = [%s]\n", crd_out);
	      found += 2;
	    }
	}
    }
  if (found != 3)
    {
      printf ("Usage: cllr_to_crd -i cllrfile -o crdfile\n");
      exit (1);
    }

/*  open input COSPAR-format lunar data file  */
  if (strlen (cospar_in) > 0)
    {
      if ((str_in = fopen (cospar_in, "r")) == NULL)
	{
	  printf ("Could not open file %s\n", cospar_in);
	  exit (1);
	}
    }

/*  open output CRD file  */
  if ((str_out_ld = fopen (crd_out, "w")) == NULL)
    {
      printf ("Could not open file %s\n", crd_out);
      exit (1);
    }
}

/*-------------------------------------------------------------------------
**
**      get_sat_ids - Given a laser target's ilrs ID, norad ID, sic, or
**                    name, get the other 3 IDs.
**
**-----------------------------------------------------------------------*/
void
get_sat_ids (int mode, int *ilrs_id, int *norad_id, int *sic, char *target,
	     int *delta)
{
  FILE *sat_id_in;
  char *sat_id_file = "./targets.dat";
  char str[256], ttarget[11];
  int status, tilrs_id, tsic, tnorad_id = 0;
  int i, l;

  if (mode != 3)
    {
      for (i = 0; i < 10; i++)
	{
	  target[i] = ' ';
	}
      target[10] = '\0';
    }

  if ((sat_id_in = fopen (sat_id_file, "r")) == NULL)
    {
      printf ("Could not open file %s\n", sat_id_file);
      exit (1);
    }
  while ((status = fgets (str, 256, sat_id_in)) != NULL)
    {
      sscanf (str, "%s %*s %d %d %d %d", ttarget, &tsic, &tilrs_id,
	      &tnorad_id, delta);

      if (mode == 0 && tilrs_id == *ilrs_id)
	{
	  *norad_id = tnorad_id;
	  *sic = tsic;
	  strcpy (target, ttarget);
	  l = strlen (target);
	  if (l < 10)
	    target[l] = ' ';
	  target[10] = '\0';
	  break;
	}
      else if (mode == 1 && tsic == *sic)
	{
	  *ilrs_id = tilrs_id;
	  *norad_id = tnorad_id;
	  strcpy (target, ttarget);
	  l = strlen (target);
	  if (l < 10)
	    target[l] = ' ';
	  target[10] = '\0';
	  break;
	}
      else if (mode == 2 && tnorad_id == *norad_id)	/* useless now */
	{
	  *ilrs_id = tilrs_id;
	  *sic = tsic;
	  strcpy (target, ttarget);
	  l = strlen (target);
	  target[10] = '\0';
	  break;
	}
      else if (mode == 3 && strcmp (ttarget, target) == 0)
	{
	  *ilrs_id = tilrs_id;
	  *norad_id = tnorad_id;
	  *sic = tsic;
	  break;
	}
    }
  fclose (sat_id_in);
  /*printf ("sic %d delta %d\n", *sic, *delta); */
}

/*-------------------------------------------------------------------------
**
**       get_station_id - get the 10 character station id by using
**                        the marker id and sites.dat file.
**
**-----------------------------------------------------------------------*/
void
get_station_id (int marker, char *station)
{
  FILE *sttn_id_in;
  char *sttn_id_file = "./sites.dat";
  char str[256], tstation[11];
  int status, tmarker;
  int i;

  for (i = 0; i < 10; i++)
    {
      station[i] = ' ';
    }
  station[10] = '\0';
  if ((sttn_id_in = fopen (sttn_id_file, "r")) == NULL)
    {
      printf ("Could not open file %s\n", sttn_id_file);
      exit (1);
    }
  while ((status = fgets (str, 256, sttn_id_in)) != NULL)
    {
      sscanf (str, "%d %s", &tmarker, tstation);

      if (tmarker == marker)
	{
	  strncpy (station, tstation, 4);
	  break;
	}
    }
  fclose (sttn_id_in);

}

/*-------------------------------------------------------------------------
**
**       get_sys_config - Get the station's SCH and SCI based
**                        on date.
**
**-----------------------------------------------------------------------*/
void
get_sys_config (int pad, int year, int month, int day, int hour,
		int *sys_change_ind, int *sys_config_ind)
{
/* Currently a dummy routine... */
  *sys_change_ind = -1;
  *sys_config_ind = -1;
}

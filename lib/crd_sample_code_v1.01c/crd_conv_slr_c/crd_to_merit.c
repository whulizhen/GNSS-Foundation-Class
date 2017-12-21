#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#define FNLEN   256
#define FMT_VERSION 0

#include "crd.h"
#include "merit.h"

FILE *str_in, *str_out;
fpos_t startpos;
struct merit_fr fr;
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

void get_win_ind ();
void setup_files ();
double interp (double, double *, double *, int);

/*-------------------------------------------------------------------------
 * Program: crd_to_merit
 *
 * Purpose:
 * Converts ILRS CRD full rate data into full rate data in the old ILRS 
 * (Merit) format.
 *
 * Calling sequence:
 *   crd_to_merit -i crd_frd_filename -o cstg_frd_filename
 *
 * Author: Randall Ricklefs / Univ of Texas / Center for Space Research
 *
 * History:
 *   Oct 26, 2007 - Initial version
 *   jul 31, 2009 - Fixes for
 *                   1) incrementing day of year at midnight rollover
 *                   2) interpolating met data to time of '10' and '11'
 *                      records. Involved accumulating all met data in first
 *                      pass. rlr.
 *   Aug 18, 2009 - Fixes for
 *                   1) incrementing year and setting doy to 1 at midnight 
 *                      rollover on last day of the year, and
 *                   2) handling file without H9 to indicate EOF. rlr.
 *   May 14, 2015 - Handle case in which the start date in the header is
 *                  late in one day, but the data records start early in
 *                  the next day. rlr.
 *
**-----------------------------------------------------------------------*/

int
main (argc, argv)
     int argc;
     char *argv[];
{
  int data_type = 0;
  int data_release = 0;
  int done = 0;
  int dummy;
  int final_doy;
  int first_record = 0;
  int fr_ready = 0;
  int fnewday= 0;
  int n20 = 0;
  int n30 = 0;
  int n_fr = 0, n_sed = 0;
  int save_angle_origin_ind = -1;
  int status;
  int first_10= 1, first_11= 1;		// First #10 or #11 record of pass
  long lclock;
  char str[256];
  double sec_of_day;
  double ssec_of_day;
  double esec_of_day;
  double save_ang_sod[500];
  double save_azimuth[500];
  double save_elevation[500];
  double save_met_sod[500];
  double save_humidity[500];
  double save_pressure[500];
  double save_temperature[500];
  double header_sod= 0;
  long double d1xsec_of_day;
  long double old_d10sec_of_day=-1;
  long double old_d11sec_of_day=-1;
  struct tm gmt;

/*  Greetings */
  printf("Convert CRD Fullrate to Merit II Format: version 1.01a 05/14/2015\n");

/*  get file names and open files  */
  setup_files (argc, argv);

  fr.np_window_ind = 0;
  fr.num_ranges = 0;
  fr.sess_rms = 0.0;

  while (done == 0)
    {
/* Initialize variables for each pass */
      n20 = 0;
      n30 = 0;
      fnewday= 0;
      save_angle_origin_ind = -1;
      old_d10sec_of_day=-1;
      old_d11sec_of_day=-1;

/* read and process records critical to header... */
      fgetpos (str_in, &startpos);
      while ((status = (int)fgets (str, 256, str_in)) != (int)NULL)
	{
          if (status == (int)NULL) break;
	  if (isalpha (str[0]))
	    str[0] = tolower (str[0]);
	  if (strncmp (str, "h1", 2) == 0 || 
	      strncmp (str, "H1", 2) == 0)
	    {
	      read_h1 (str, &h1);
	      fr.format_version = 3;	/* latest format revision */
	    }
	  else if (strncmp (str, "h2", 2) == 0 ||
	           strncmp (str, "H2", 2) == 0)
	    {
	      read_h2 (str, &h2);
	      fr.cdp_pad_id = h2.cdp_pad_id;
	      fr.cdp_sys_num = h2.cdp_sys_num;
	      fr.cdp_occ_num = h2.cdp_occ_num;
	      fr.stn_timescale = h2.stn_timescale;
	    }
	  else if (strncmp (str, "h3", 2) == 0 ||
	           strncmp (str, "H3", 2) == 0)
	    {
	      read_h3 (str, &h3);
	      fr.ilrs_id = h3.ilrs_id;
	    }
	  else if (strncmp (str, "h4", 2) == 0 ||
	           strncmp (str, "H4", 2) == 0)
	    {
	      read_h4 (str, &h4);
	      fr.data_release = h4.data_release;
	      fr.year = h4.start_year % 100;
	      grtodoy (h4.start_year, h4.start_mon, h4.start_day, &fr.doy);
	      grtodoy (h4.start_year, 12, 31, &final_doy);
              header_sod= h4.start_hour*3600 + h4.start_min*60+ h4.start_sec;
              first_10= 1; 
              first_11= 1;		// First #10 or #11 record of pass

	      if (h4.refraction_app_ind == 0)
		fr.refraction_app_ind = 1;
	      else
		fr.refraction_app_ind = 0;

	      if (h4.CofM_app_ind == 0)
		fr.CofM_app_ind = 1;
	      else
		fr.CofM_app_ind = 0;

	      if (h4.xcv_amp_app_ind == 0)
		fr.xcv_amp_app_ind = 1;
	      else
		fr.xcv_amp_app_ind = 0;
	    }
	  else if (strncmp (str, "h8", 2) == 0 ||
	           strncmp (str, "H8", 2) == 0)
	    {
	      if (fr_ready)
		{
		  write_merit_fr (str_out, fr);
		  fr_ready = 0;
		}
	      read_h8 (str);
	      break;
	    }
	  else if (strncmp (str, "c0", 2) == 0 ||
	           strncmp (str, "C0", 2) == 0)
	    {
	      read_c0 (str, &c0);
	      fr.xmit_wavelength = c0.xmit_wavelength;
	      if (fr.xmit_wavelength < 1000)
		(fr.xmit_wavelength) *= 10;
	    }
	  else if (strncmp (str, "10", 2) == 0)
	    {
	      read_10 (str, &d10);
	      fr.epoch_event = d10.epoch_event;
	    }
	  else if (strncmp (str, "11", 2) == 0)
	    {
	      read_11 (str, &d11);
	      fr.epoch_event = d11.epoch_event;
	      /* For the moon, this should really take the max of all fr lengths */
	      get_win_ind ((int)d11.np_window_length, h3.ilrs_id,
			   &fr.np_window_ind, &dummy);
	    }
	  else if (strncmp (str, "12", 2) == 0)
	    {
	      read_12 (str, &d12);
	      fr.refraction_corr = d12.refraction_corr;
	      fr.target_CofM_corr = d12.target_CofM_corr;
	    }
	  else if (strncmp (str, "20", 2) == 0)
	    {

	      read_20 (str, &d20);
              /* New day ...  since last '20' or header; don't rely on 
               * stop time as it is '0' for some stations */
              save_met_sod[n20] = d20.sec_of_day;
              if ((n20 > 0 && save_met_sod[n20] < save_met_sod[0]) ||
                  (n20 == 0 && header_sod > save_met_sod[0]))
                save_met_sod[n20]+= 86400.e0;
	      save_humidity[n20] = d20.humidity;
	      save_pressure[n20] = d20.pressure * 10;
	      save_temperature[n20] = d20.temperature * 10;
	      if (n20 < 499) n20++;	/* Should scream, but >500 unlikely */
	    }
	  else if (strncmp (str, "30", 2) == 0)
	    {
	      read_30 (str, &d30);
              /* IF angle origin has changed within pass, use only first one */
              if (n30 == 0 || 
                 (n30 > 0 && save_angle_origin_ind == d30.angle_origin_ind))
                {
                  save_ang_sod[n30] = d30.sec_of_day;
                  if (n30 > 0 && save_ang_sod[n30] < save_ang_sod[0])
                    save_ang_sod[n30]+= 86400.e0;
	          save_azimuth[n30] = d30.azimuth * 1.e4;
	          save_elevation[n30] = d30.elevation * 1.e4;
	          fr.angle_origin_ind = d30.angle_origin_ind;
	          save_angle_origin_ind = d30.angle_origin_ind;
	          if (n30 < 499) n30++;	/* Should scream, but >500 unlikely */
                }
	    }
	  else if (strncmp (str, "40", 2) == 0)
	    {
	      read_40 (str, &d40);
	      fr.cal_sys_delay = d40.cal_sys_delay;
	      fr.cal_delay_shift = d40.cal_delay_shift;
	      fr.cal_rms = d40.cal_rms;
	      fr.cal_type_ind = d40.cal_type_ind - 2;
	      if (d40.cal_shift_type_ind == 3)
		fr.cal_type_ind += 5;
	    }
	  else if (strncmp (str, "50", 2) == 0)
	    {
	      read_50 (str, &d50);
	      fr.sess_rms = d50.sess_rms;
	    }
	  else if (strncmp (str, "60", 2) == 0)
	    {
	      read_60 (str, &d60);
	      fr.sys_change_ind = d60.sys_change_ind;
	      fr.sys_config_ind = d60.sys_config_ind;
	    }
	}

/* read and process the file... */
      fsetpos (str_in, &startpos);
      while ((status = (int)fgets (str, 256, str_in)) != (int)NULL)
	{
          if (status == (int)NULL) break;
	  if (isalpha (str[0]))
	    str[0] = tolower (str[0]);
	  if (strncmp (str, "h8", 2) == 0 ||
	      strncmp (str, "H8", 2) == 0)
	    {
	      read_h8 (str);
	      if (fr_ready)
		{
		  write_merit_fr (str_out, fr);
		  fr_ready = 0;
		}
	      break;
	    }
	  else if (strncmp (str, "h9", 2) == 0 ||
	           strncmp (str, "H9", 2) == 0)
	    {
	      read_h9 (str);
	      if (fr_ready)
		{
		  write_merit_fr (str_out, fr);
		  fr_ready = 0;
		}
	      fclose (str_in);
	      fclose (str_out);
	      done = 1;
	      break;
	    }
	  else if (strncmp (str, "c1", 2) == 0 ||
	           strncmp (str, "C1", 2) == 0)
	    {
	      read_c1 (str, &c1);
	    }
	  else if (strncmp (str, "c2", 2) == 0 ||
	           strncmp (str, "C2", 2) == 0)
	    {
	      read_c2 (str, &c2);
	    }
	  else if (strncmp (str, "c3", 2) == 0 ||
	           strncmp (str, "C3", 2) == 0)
	    {
	      read_c3 (str, &c3);
	    }
	  else if (strncmp (str, "c4", 2) == 0 ||
	           strncmp (str, "C4", 2) == 0)
	    {
	      read_c4 (str, &c4);
	    }
	  else if (strncmp (str, "00", 2) == 0)
	    {
	      read_00 (str, &d00);
	    }
	  else if (strncmp (str, "10", 2) == 0)
	    {
	      if (fr_ready)
		{
		  write_merit_fr (str_out, fr);
		  fr_ready = 0;
		}
	      read_10 (str, &d10);

	      /* Copy in data record info */
	      fr.sec_of_day = d10.sec_of_day * 1.e7;
	      d1xsec_of_day = d10.sec_of_day;
	      fr.time_of_flight = d10.time_of_flight;

              /* New day ...  since last '10' or header; don't rely on 
               * stop time as it is '0' for some stations */
	      if (d10.sec_of_day < old_d10sec_of_day ||
                  (first_10 && header_sod > d10.sec_of_day))
                {
		  fr.doy++;
                  fnewday= 1;
                  if (fr.doy > final_doy)
                    {
                      fr.doy= 1;
                      fr.year++;
                    }
                  first_10= 0;
                }
	      old_d10sec_of_day= d10.sec_of_day;

	      /* Assuming this isn't high-rep-rate! */
              /* Interpolate mets */
	      fr.time_of_flight *= 1.e12;
	      fr.humidity = interp((double)d10.sec_of_day+fnewday*86400.e0,
                 save_met_sod,save_humidity,n20) + 0.51;
	      fr.pressure = interp((double)d10.sec_of_day+fnewday*86400.e0,
                 save_met_sod,save_pressure,n20)+ 0.51;
	      fr.temperature = interp((double)d10.sec_of_day+fnewday*86400.e0,
                 save_met_sod,save_temperature,n20)+ 0.51;
              /*printf("pth: %d %lf %f %f %f\n",fnewday,(double)d10.sec_of_day+
                fnewday*86400.e0,fr.humidity,fr.pressure,fr.temperature);*/
	      fr.num_ranges = 0;

              /* Interpolate angles */
	      fr.azimuth = interp((double)d10.sec_of_day+fnewday*86400.e0,
                 save_ang_sod,save_azimuth,n30)+ 0.51;
	      fr.elevation = interp((double)d10.sec_of_day+fnewday*86400.e0,
                 save_ang_sod,save_elevation,n30)+ 0.51;

	      /* Assuming this isn't lunar! */
	      while (fr.num_ranges > 9999)
		{
		  fr.num_ranges /= 10;
		}
	      /* Write record */
	      fr_ready = 1;
	      /*write_merit_fr (str_out, fr); */
	    }
	  /* npts not likely, but just in case... */
	  else if (strncmp (str, "11", 2) == 0)
	    {
	      if (fr_ready)
		{
		  write_merit_fr (str_out, fr);
		  fr_ready = 0;
		}
	      read_11 (str, &d11);

	      /* Copy in data record info */
	      fr.sec_of_day = d11.sec_of_day * 1.e7;
	      d1xsec_of_day = d11.sec_of_day;

              /* New day ...  since last '10' or header; don't rely on 
               * stop time as it is '0' for some stations */
	      if (d10.sec_of_day < old_d10sec_of_day ||
                  (first_11 && header_sod > d10.sec_of_day))
	      //if (d11.sec_of_day < old_d11sec_of_day) 	/* New day */
                {
		  fr.doy++;
                  fnewday= 1;
                  if (fr.doy > final_doy)
                    {
                      fr.doy= 1;
                      fr.year++;
                    }
                  first_11= 0;
                }
	      old_d11sec_of_day= d11.sec_of_day;
	      fr.time_of_flight = d11.time_of_flight;
	      fr.humidity = interp((double)d11.sec_of_day+fnewday*86400.e0,
                 save_met_sod,save_humidity,n20) + 0.51;
	      fr.pressure = interp((double)d11.sec_of_day+fnewday*86400.e0,
                 save_met_sod,save_pressure,n20)+ 0.51;
	      fr.temperature = interp((double)d11.sec_of_day+fnewday*86400.e0,
                 save_met_sod,save_temperature,n20)+ 0.51;
              /*printf("pth: %d %lf %f %f %f\n",fnewday,(double)d11.sec_of_day+
                fnewday*86400.e0,fr.humidity,fr.pressure,fr.temperature);*/

	      /* Assuming this isn't high-rep-rate! */
	      fr.time_of_flight *= 1.e12;
	      fr.num_ranges = d11.num_ranges;

	      /* Assuming this isn't lunar! */
	      while (fr.num_ranges > 9999)
		{
		  fr.num_ranges /= 10;
		}
	      /* Write record */
	      fr_ready = 1;
	    }
	  else if (strncmp (str, "12", 2) == 0)
	    {
	      read_12 (str, &d12);
	      fr.refraction_corr= d12.refraction_corr;
	      fr.target_CofM_corr = d12.target_CofM_corr;
	    }
	  else if (strncmp (str, "21", 2) == 0)
	    {
	      read_21 (str, &d21);
	    }
	  else if (strncmp (str, "40", 2) == 0)
	    {
	      read_40 (str, &d40);
	    }
	  else if (strncmp (str, "50", 2) == 0)
	    {
	      read_50 (str, &d50);
	    }
	  else if (strncmp (str, "60", 2) == 0)
	    {
	      read_60 (str, &d60);
	    }
	  else if (strncmp (str, "9X", 1) == 0)
	    {
	      /* User-defined! */
	      /*read_00 (str, &d00); */
	    }
	}
      if (status == (int)NULL)
        {
	  if (fr_ready)
	    {
	      write_merit_fr (str_out, fr);
	      fr_ready = 0;
	    }
	  fclose (str_in);
	  fclose (str_out);
          done = 1;
        }
    }
}

/*-------------------------------------------------------------------------
**
**      get_win_ind             - Determine window indicator from the
**                                window length and ilrs id
**
**-----------------------------------------------------------------------*/
void
get_win_ind (int np_window_length, int ilrs_id, int *nwi, int *lnwi)
{

/* Lunar data: This test will not hold if other lunar reflectors are
   comissioned. */
  if (ilrs_id >= 100 && ilrs_id <= 104)
    {
      *nwi = 2;
      if (np_window_length <= 300)
	*lnwi = 1;
      if (np_window_length > 300 && np_window_length <= 600)
	*lnwi = 2;
      if (np_window_length > 600 && np_window_length <= 900)
	*lnwi = 3;
      if (np_window_length > 900 && np_window_length <= 1200)
	*lnwi = 4;
      if (np_window_length > 1200 && np_window_length <= 1500)
	*lnwi = 5;
      if (np_window_length > 1500 && np_window_length <= 1800)
	*lnwi = 6;
      if (np_window_length > 1800 && np_window_length <= 2100)
	*lnwi = 7;
      if (np_window_length > 2100 && np_window_length <= 2400)
	*lnwi = 8;
      if (np_window_length > 2400)
	*lnwi = 9;
      return;
    }

  if (np_window_length <= 5)
    *nwi = 1;
  if (np_window_length > 5 && np_window_length <= 15)
    *nwi = 3;
  if (np_window_length > 15 && np_window_length <= 20)
    *nwi = 4;
  if (np_window_length > 20 && np_window_length <= 30)
    *nwi = 5;
  if (np_window_length > 30 && np_window_length <= 60)
    *nwi = 6;
  if (np_window_length > 60 && np_window_length <= 120)
    *nwi = 7;
  if (np_window_length > 120 && np_window_length <= 180)
    *nwi = 8;
  if (np_window_length > 180 && np_window_length <= 300)
    *nwi = 9;
}

/*-------------------------------------------------------------------------
**
**	setup_files		- Open input and output file names
**
**-----------------------------------------------------------------------*/
void
setup_files (int argc, char *argv[])
{
  char crd_in[FNLEN] = { "" };
  char cstg_out[FNLEN] = { "" };

  int found = 0, i;

/*  Get file names  */
  for (i = 1; i < argc; i += 2)
    {
      if (argv[i][0] == '-')
	{
	  if (argv[i][1] == 'i')
	    {
	      strcpy (crd_in, argv[i + 1]);
	      printf ("crd_in = [%s]\n", crd_in);
	      found += 1;
	    }
	  else if (argv[i][1] == 'o')
	    {
	      if (found == 1)
		{
		  strcpy (cstg_out, argv[i + 1]);
		}
	      printf ("cstg_out = [%s]\n", cstg_out);
	      found += 2;
	    }
	}
    }
  if (found != 3)
    {
      printf ("Usage: crd_to_merit -i crd_frd_file -o merit_fr_file\n");
      exit (1);
    }

/*  open input crd file  */
  if (strlen (crd_in) > 0)
    {
      if ((str_in = fopen (crd_in, "r")) == NULL)
	{
	  printf ("Could not open file %s\n", crd_in);
	  exit (1);
	}
    }

/*  open output CSTG files  */
  if ((str_out = fopen (cstg_out, "w")) == NULL)
    {
      printf ("Could not open file %s\n", cstg_out);
      exit (1);
    }
}

/*-------------------------------------------------------------------------
**
**	interp		- Linear interpolation of y array
**			  07/31/09
**
**-----------------------------------------------------------------------*/
double
interp (double x, double *xarr, double *yarr, int nx)
{
  int i;

  if (nx == 1)
    return (yarr[0]);
  for (i=0;i<nx-1;i++)
    {
       if (x >= xarr[i] && x < xarr[i+1])
         return (yarr[i]+ (yarr[i+1]- yarr[i])*(x- xarr[i])/(xarr[i+1]-xarr[i]));
    }
  if (x < xarr[0])
    return(yarr[0]);
  if (x >= xarr[nx-1])
    return(yarr[nx-1]);
}

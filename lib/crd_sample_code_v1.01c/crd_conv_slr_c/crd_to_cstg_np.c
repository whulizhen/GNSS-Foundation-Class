#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#define FNLEN   256
#define FMT_VERSION 0

#include "crd.h"
#include "cstg.h"

FILE *str_in, *str_out;
fpos_t startpos;
struct cstg_hdr hdr;
struct cstg_sed sed;
struct cstg_np np;
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
 * Program: crd_to_cstg_np
 *
 * Purpose:
 * Converts ILRS CRD normalpoints into normalpoints in the old ilrs (cstg)
 * format
 *
 * Calling sequence:
 *   crd_to_cstg_np -i crd_np_filename -o cstg_np_filename
 *
 * Author: Randall Ricklefs / Univ of Texas / Center for Space Research
 *
 * History:
 *   July 2, 2007 - Initial version
 *   May 4, 2009  - Zero fill snr when return rate is set to -1 (n.a.). rlr.
 *   Aug 18, 2009 - Add intepolation for met values and 
 *                  properly handle file w/o ending H9. rlr.
 *
**-----------------------------------------------------------------------*/

int
main (argc, argv)
     int argc;
     char *argv[];
{
  int data_type = 0;
  int data_release = 0;
  int done= 0;
  int first_met= 1;
  int first_record = 0;
  int n_np = 0, n_sed = 0;
  int np_ready= 0;
  int fnewday= 0;
  int n20 = 0;
  int status;
  long lclock;
  char str[256];
  double sec_of_day;
  double ssec_of_day;
  double esec_of_day;
  double save_met_sod[500];
  double save_humidity[500];
  double save_pressure[500];
  double save_temperature[500];
  long double old_d11sec_of_day=-1;
  struct tm gmt;

/*  Greetings */
  printf("Convert CRD Normalpoints to CSTG Format: version 1.01 09/10/2009\n");

/*  get file names and open files  */
  setup_files (argc, argv);

  while (done == 0)
    {
/* Initialize variables for each pass */
      n20 = 0;
      fnewday= 0;
      old_d11sec_of_day=-1;

/* read and process records critical to header... */
      first_met= 1;
      fgetpos (str_in, &startpos);
      while ((status = (int)fgets (str, 256, str_in)) != (int)NULL)
	{
	  if (strncmp (str, "h1", 2) == 0 ||
	      strncmp (str, "H1", 2) == 0)
	    {
	      read_h1 (str, &h1);
	      hdr.format_version = 2;	/* latest format revision */
	    }
	  else if (strncmp (str, "h2", 2) == 0 ||
	           strncmp (str, "H2", 2) == 0)
	    {
	      read_h2 (str, &h2);
	      hdr.cdp_pad_id = h2.cdp_pad_id;
	      hdr.cdp_sys_num = h2.cdp_sys_num;
	      hdr.cdp_occ_num = h2.cdp_occ_num;
	      hdr.stn_timescale = h2.stn_timescale;
              /* No time scale > 9 in old format */
              if (hdr.stn_timescale > 9) hdr.stn_timescale= 0;
	    }
	  else if (strncmp (str, "h3", 2) == 0 ||
	           strncmp (str, "H3", 2) == 0)
	    {
	      read_h3 (str, &h3);
	      hdr.ilrs_id = h3.ilrs_id;
	    }
	  else if (strncmp (str, "h4", 2) == 0 ||
	           strncmp (str, "H4", 2) == 0)
	    {
	      read_h4 (str, &h4);
	      np.data_release = h4.data_release;
	      hdr.year = h4.start_year % 100;
	      grtodoy (h4.start_year, h4.start_mon, h4.start_day, &hdr.doy);
	    }
	  else if (strncmp (str, "h8", 2) == 0 ||
	           strncmp (str, "H8", 2) == 0)
	    {
	      if (np_ready)
                {
	          write_cstg_np (str_out, np);
		  np_ready= 0;
		}
	      read_h8 (str);

              /* Copy data to cstg header and write */
              hdr.checksum= 0;
              write_cstg_hdr (str_out, hdr, 1);
	      break;
	    }
	  else if (strncmp (str, "c0", 2) == 0 ||
	           strncmp (str, "C0", 2) == 0)
	    {
	      read_c0 (str, &c0);
	      hdr.xmit_wavelength = c0.xmit_wavelength;
	      if (hdr.xmit_wavelength < 1000)
		(hdr.xmit_wavelength) *= 10;
              hdr.xmit_wavelength= (int)(hdr.xmit_wavelength + 0.5);
	    }
	  else if (strncmp (str, "11", 2) == 0)
	    {
	      read_11 (str, &d11);
	      /* For the moon, this should really take the max of all np lengths */
	      get_win_ind ((int)d11.np_window_length, h3.ilrs_id,
			   &hdr.np_window_ind, &np.llr_np_window_ind);
	    }
          else if (strncmp (str, "20", 2) == 0)
            {
              read_20 (str, &d20);
              save_met_sod[n20] = d20.sec_of_day;
              if (n20 > 0 && save_met_sod[n20] < save_met_sod[0])
                save_met_sod[n20]+= 86400.e0;
              save_humidity[n20] = d20.humidity;
              save_pressure[n20] = d20.pressure * 10;
              save_temperature[n20] = d20.temperature * 10;
              if (n20 < 499) n20++;     /* Should scream, but >500 unlikely */
            }
	  else if (strncmp (str, "40", 2) == 0)
	    {
	      read_40 (str, &d40);
	      hdr.cal_sys_delay = (int)(d40.cal_sys_delay + 0.5);
	      hdr.cal_delay_shift = (int)(d40.cal_delay_shift + 0.5);
	      hdr.cal_rms = (int)(d40.cal_rms + 0.5);       /* Round off fix */
	      hdr.cal_type_ind = d40.cal_type_ind - 2;
	      if (d40.cal_shift_type_ind == 3)
		hdr.cal_type_ind += 5;
	    }
	  else if (strncmp (str, "50", 2) == 0)
	    {
	      read_50 (str, &d50);
	      hdr.sess_rms = (int)(d50.sess_rms + 0.5);     /* Round off fix */
	      hdr.data_qual_ind = d50.data_qual_ind;
	    }
	  else if (strncmp (str, "60", 2) == 0)
	    {
	      read_60 (str, &d60);
	      hdr.sys_change_ind = d60.sys_change_ind;
	      hdr.sys_config_ind = d60.sys_config_ind;
	    }
	}

/* read and process the file... */
      fsetpos (str_in, &startpos);
      while ((status = (int)fgets (str, 256, str_in)) != (int)NULL)
	{
	  if (strncmp (str, "h8", 2) == 0 ||
	      strncmp (str, "H8", 2) == 0)
	    {
	      read_h8 (str);
	      if (np_ready)
                {
	          write_cstg_np (str_out, np);
		  np_ready= 0;
		}
	      break;
	    }
	  else if (strncmp (str, "h9", 2) == 0 ||
	           strncmp (str, "H9", 2) == 0)
	    {
	      read_h9 (str);
	      if (np_ready)
                {
	          write_cstg_np (str_out, np);
		  np_ready= 0;
		}
	      fclose (str_in);
	      fclose (str_out);
	      done= 1;
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
	      read_10 (str, &d10);
	    }
	  else if (strncmp (str, "11", 2) == 0)
	    {
	      if (np_ready)
                {
	          write_cstg_np (str_out, np);
		  np_ready= 0;
		}
	      read_11 (str, &d11);

	      /* Copy in data record info */
	      np.sec_of_day = d11.sec_of_day * 1.e7;
	      np.time_of_flight = d11.time_of_flight;
              if (d11.sec_of_day < old_d11sec_of_day)   /* New day */
                {
                  fnewday= 1;
                }
              old_d11sec_of_day= d11.sec_of_day;


	      /* Assuming this isn't high-rep-rate! */
	      np.scale_or_tof_sec = 0;
	      if (np.time_of_flight >= 1.e0)
		{
		  np.scale_or_tof_sec = (int) np.time_of_flight;
		  np.time_of_flight -= np.scale_or_tof_sec;
		}
	      np.time_of_flight *= 1.e12;
	      np.num_ranges = d11.num_ranges;

              /* Interpolate mets */
              np.humidity = interp((double)d11.sec_of_day+fnewday*86400.e0,
                 save_met_sod,save_humidity,n20) + 0.51;
              np.pressure = interp((double)d11.sec_of_day+fnewday*86400.e0,
                 save_met_sod,save_pressure,n20)+ 0.51;
              np.temperature = interp((double)d11.sec_of_day+fnewday*86400.e0,
                 save_met_sod,save_temperature,n20)+ 0.51;

	      /* Assuming this isn't lunar! */
	      while (np.num_ranges > 9999)
		{
		  np.num_ranges /= 10;
		  np.scale_or_tof_sec++;
		}
	      np.bin_rms = (int)(d11.bin_rms + 0.5);       /* Round off fix */
	      np.snr = d11.return_rate;
              /* In case return_rate was set to -1 (na) or > 99. 05/22/09 */
              if (d11.return_rate < 0) np.snr= 0;
              if (d11.return_rate > 99) np.snr= 99;
	      np.checksum= 0;
	      /* Write record */
	      np_ready= 1;
	    }
	  else if (strncmp (str, "12", 2) == 0)
	    {
	      read_12 (str, &d12);
	    }
	  else if (strncmp (str, "21", 2) == 0)
	    {
	      read_21 (str, &d21);
	    }
	  else if (strncmp (str, "30", 2) == 0)
	    {
	      /* This is not expected for a np file! */
	      read_30 (str, &d30);
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
          if (np_ready)
            {
              write_cstg_np (str_out, np);
              np_ready = 0;
            }
          fclose (str_in);
          fclose (str_out);
          done = 1;
        }
    }
}

/*-------------------------------------------------------------------------
**
**	get_win_ind		- Determine window indicator from the
**				  window length and ilrs id
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
      printf ("Usage: crd_to_cstg_np -i crd_npt_file -o cstg_npt_file\n");
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
 * **
 * **      interp          - Linear interpolation of y array
 * **                        07/31/09
 * **
 * **-----------------------------------------------------------------------*/
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


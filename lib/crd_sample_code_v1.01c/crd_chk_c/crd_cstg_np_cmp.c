#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#define FNLEN   256
#define FMT_VERSION 0

#include "crd.h"
#include "cstg.h"

FILE *str_in1, *str_in2, *str_out1, *str_out2 = 0;
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

char pvern[] = "1.00";
int debug = 0;
int np_missing = 0, d11_missing = 0;

void check_hdr ();
int check_npt ();
void session_error_report ();
void get_win_ind ();
void setup_files ();
void end_it ();

/*-------------------------------------------------------------------------
 * Program: crd_cstg_np_cmp
 *
 * Purpose:
 * Compares CRD and old-style cstg normal point files
 *
 * Calling sequence:
 *   crd_cstg_np_cmp -1 crd_npt_file -2 cstg_npt_file -o report_file \
 *       [-s summary_file] [-v]
 *     -1 is the CRD-formatted normalpoint file
 *     -2 is the CSTG-formatted (old format) normalpoint file
 *     -o output report file
 *     -s one line summary is appended to this file
 *     -v prints program version number and stops.
 *
 * Author: Randall Ricklefs / Univ of Texas / Center for Space Research
 *
 * History:
 *   25 May 2008 - initial version
 *
**-----------------------------------------------------------------------*/

int
main (argc, argv)
     int argc;
     char *argv[];
{
  int check_status = 0;
  int crd_met_read = 0;		/* Has a meteorological record been read? */
  int cstg_npt_read = 0;	/* Has a cstg normalpoint record been read? */
  int read_sep = 0;		/* Has a cstg separator record been read? */
  int status;
  char str[256];
  char crd_in[FNLEN];

/*  get file names and open files  */
  setup_files (argc, argv, crd_in);

/* read and process all records... */
  while ((status = fgets (str, 256, str_in1)) != NULL)
    {
      if (debug)
	printf ("str_in1: [%s]\n", str);
      if (strncmp (str, "h1", 2) == 0 || strncmp (str, "H1", 2) == 0)
	{
          /* File header */
	  read_h1 (str, &h1);
	}
      else if (strncmp (str, "h2", 2) == 0 || strncmp (str, "H2", 2) == 0)
	{
          /* Station header */
	  read_h2 (str, &h2);
	}
      else if (strncmp (str, "h3", 2) == 0 || strncmp (str, "H3", 2) == 0)
	{
          /* Target header */
	  read_h3 (str, &h3);
	}
      else if (strncmp (str, "h4", 2) == 0 || strncmp (str, "H4", 2) == 0)
	{
          /* Session header */
	  read_h4 (str, &h4);
	}
      else if (strncmp (str, "h8", 2) == 0 || strncmp (str, "H8", 2) == 0)
	{
          /* End of CRD session data block. */
	  read_h8 (str);
	  if (debug)
	    printf ("read h8\n");

	  if (!read_sep)
	    {
	      /* read the rest of the cstg npt records */
	      do
		{
		  if ((status = fgets (str, 256, str_in2)) != NULL &&
		      strncmp (str, "99999", 5) != 0 &&
		      strncmp (str, "88888", 5) != 0)
		    d11_missing++;
		  if (debug)
		    printf ("np [%s] status %d d11_missing %d\n", str,
			    status, d11_missing);
		}
	      while (strncmp (str, "99999", 5) != 0 &&
		     strncmp (str, "88888", 5) != 0 && status != NULL);
	    }
	  if (debug)
	    printf ("np str at h8 = [%s]\n", str);
	  session_error_report (crd_in);
	  read_sep = 0;

	  /* Read rest of .np/.qld file and bypass sampled engineering 
	   * data. */
	  do
	    {
	      status = fgets (str, 256, str_in2);
	      if (debug)
		printf ("Bypassing %s\n", str);
	    }
	  while (strncmp (str, "99999", 5) != 0 &&
		 strncmp (str, "88888", 5) != 0 && status != NULL);
	}
      else if (strncmp (str, "h9", 2) == 0 || strncmp (str, "H9", 2) == 0)
	{
          /* End of CRD file. */
	  read_h9 (str);
	  if (debug)
	    printf ("read h9\n");
	  end_it (0);
	}
      else if (strncmp (str, "c0", 2) == 0 || strncmp (str, "C0", 2) == 0)
	{
          /* System configuration header */
	  read_c0 (str, &c0);

	  /* read 88888 or 99999 separator record */
	  do
	    {
	      if ((status = fgets (str, 256, str_in2)) == NULL)
		end_it (1);
	    }
	  while (strncmp (str, "99999", 5) != 0 &&
		 strncmp (str, "88888", 5) != 0);

	  /* read cstg np header */
	  if ((status = fgets (str, 256, str_in2)) == NULL)
	    end_it (1);
	  if (debug)
	    printf ("read header [%s]\n", str);
	  read_cstg_hdr (str, &hdr);
	  check_hdr ();
	}
      else if (strncmp (str, "11", 2) == 0)
	{
          /* Normal point record */
	  read_11 (str, &d11);

	  /* If the CRD normalpoint time is earlier that the cstg,
	   * skip the record */
	  do
	    {
	      if (check_status >= 0)	/* play catch-up */
		{
		  if ((status = fgets (str, 256, str_in2)) == NULL)
		    {
		      np_missing++;
		      if (debug)
			printf ("np missing (11)\n");
		      break;
		    }
		  if (strncmp (str, "99999", 5) == 0 ||
		      strncmp (str, "88888", 5) == 0)
		    {
		      read_sep = 1;
		      np_missing++;
		      break;
		    }
		  if (debug)
		    printf ("11 read np [%s]\n", str);
		  read_cstg_np (str, &np);
		  if (!cstg_npt_read && crd_met_read)
		    check_npt (20);
		  cstg_npt_read = 1;
		}
	    }
	  while ((check_status = check_npt (11)) > 0);
	}
      else if (strncmp (str, "20", 2) == 0)
	{
          /* Meteorological record */
	  read_20 (str, &d20);
	  if (cstg_npt_read)
	    check_npt (20);
	  crd_met_read = 1;
	}
      else if (strncmp (str, "40", 2) == 0)
	{
          /* Calibration record */
	  read_40 (str, &d40);
	  check_npt (40);
	}
      else if (strncmp (str, "50", 2) == 0)
	{
          /* Session summary record */
	  read_50 (str, &d50);
	  check_npt (50);
	}
      else if (strncmp (str, "60", 2) == 0)
	{
          /* Compatibility record */
	  read_60 (str, &d60);
	  check_npt (60);
	}
    }
}

double dpres[50];
double dtemp[50];
double dhum[50];
double dft[50];
double dr[50];
double dcsd[50];
double dcds[50];
double dbrms[50];
double dcrms[50];
double dsrms[50];
double toler = 0.15;
int ndr_err = 0;
int dnret[50];
int npth = 0, nfr = 0;
int ncal = 0;
int nsess = 0;
int ftbd11, ftbnp, dftb;
int nwin_ind_err = 0;
int crd_np_window_ind, crd_llr_np_window_ind;
int cti_err = 0;
int csti_err = 0;
int crd_cal_type_ind, crd_cal_shift_type_ind;
int sci_err = 0;
int scf_err = 0;
int dqi_err = 0;
int iid_err = 0;
int rectype_err = 0;
int year_err = 0;
int doy_err = 0;
int padid_err = 0;
int sysnum_err = 0;
int occnum_err = 0;
int wavelength_err = 0;
int timescale_err = 0;
int hdr_err = 0;
int fail = 0;
int marginal = 0;
int pass = 0;
int strongpass = 0;
int crd_doy;

/*-------------------------------------------------------------------------
**
**	check_hdr		- Check contents of cstg header against CRD
**
**-----------------------------------------------------------------------*/
void
check_hdr ()
{
  if (hdr.ilrs_id != h3.ilrs_id)
    iid_err = 1;
  if (hdr.year != h4.start_year)
    year_err = 1;
  grtodoy (h4.start_year % 100, h4.start_mon, h4.start_day, &crd_doy);
  if (hdr.doy != crd_doy)
    doy_err = 1;
  if (hdr.cdp_pad_id != h2.cdp_pad_id)
    padid_err = 1;
  if (hdr.cdp_sys_num != h2.cdp_sys_num)
    sysnum_err = 1;
  if (hdr.cdp_occ_num != h2.cdp_occ_num)
    occnum_err = 1;
  if (hdr.xmit_wavelength != c0.xmit_wavelength)
    wavelength_err = 1;
  if (hdr.stn_timescale != h2.stn_timescale)
    timescale_err = 1;
  hdr_err = iid_err + year_err + doy_err + padid_err + sysnum_err +
    occnum_err + wavelength_err + timescale_err;
}

/*-------------------------------------------------------------------------
**
**	check_npt		- Check contents of cstg normal point record
**                                against CRD records
**
**-----------------------------------------------------------------------*/
int
check_npt (int recd_type)
{
  double crd_cal_type_ind;
  double crd_cal_shift_type_ind;

  if (recd_type == 11)
    {

      /* Are both firing times in same normal point bin? */
      ftbd11 = d11.sec_of_day / d11.np_window_length;
      ftbnp = np.sec_of_day / d11.np_window_length;

/* kludge for MLRS lunar npt window length */
      if (h3.target_type == 2 &&
	  fabs (d11.sec_of_day - np.sec_of_day) < d11.np_window_length / 2.)
	ftbnp = ftbd11;

      if (debug)
	printf ("%Lf %Lf %f\n", d11.sec_of_day, np.sec_of_day,
		d11.np_window_length);
      if (debug)
	printf ("bins: %d %d\n", ftbd11, ftbnp);
      /* No */
      if (ftbd11 != ftbnp)
	{
	  dftb = ftbd11 - ftbnp;
	  if (debug)
	    printf ("dftb= %d\n", dftb);
	  if (dftb < 0)
	    {
	      np_missing++;
	      if (debug)
		printf ("np_missing (chknpt)\n");
	    }
	  else
	    {
	      d11_missing++;
	      if (debug)
		printf ("d11_missing (chknpt)\n");
	    }
	  return (dftb);	/* be sure to unread if necessary */
	}
      /* Yes */
      dft[nfr] = d11.sec_of_day - np.sec_of_day;
      if (debug)
	printf ("ft %.12Lf %.12Lf\n", d11.sec_of_day, np.sec_of_day);
      dr[nfr] =
	d11.time_of_flight - (np.scale_or_tof_sec + np.time_of_flight);
      if (debug)
	printf ("r %.12Lf %.12Lf\n", d11.time_of_flight, np.time_of_flight);
      get_win_ind ((int)d11.np_window_length, h3.ilrs_id,
		   &crd_np_window_ind, &crd_llr_np_window_ind);
      if (crd_np_window_ind != hdr.np_window_ind)
	nwin_ind_err++;
      dnret[nfr] = d11.num_ranges - np.num_ranges;
      dbrms[nfr] = d11.bin_rms - np.bin_rms;
      if (debug)
	printf ("n rms %f %f\n", d11.bin_rms, np.bin_rms);
      if (h4.data_release != np.data_release)
	ndr_err++;
      nfr++;

    }
  else if (recd_type == 20)
    {
      /* check dtime?? What if no np read? */
      if (debug)
	printf ("pth: %f %f    %f %f   %f %f\n", d20.pressure, np.pressure,
		d20.temperature, np.temperature, d20.humidity, np.humidity);
      dpres[npth] = d20.pressure - np.pressure;
      dtemp[npth] = d20.temperature - np.temperature;
      dhum[npth] = d20.humidity - np.humidity;
      npth++;
    }
  else if (recd_type == 40)
    {
      dcsd[ncal] = d40.cal_sys_delay - hdr.cal_sys_delay;
      dcds[ncal] = d40.cal_delay_shift - hdr.cal_delay_shift;
      dcrms[ncal] = d40.cal_rms - hdr.cal_rms;
      if (debug)
	printf ("cal: %f %f\n", d40.cal_rms, hdr.cal_rms);

      crd_cal_shift_type_ind = hdr.cal_type_ind > 4 ? 3 : 2;
      crd_cal_type_ind = hdr.cal_type_ind;
      if (crd_cal_type_ind > 4)
	crd_cal_type_ind -= 5;
      crd_cal_type_ind += 2;
      if (crd_cal_type_ind == 6)
	crd_cal_type_ind = 0;

      if (d40.cal_type_ind != crd_cal_type_ind)
	cti_err = 1;
      if (d40.cal_shift_type_ind != crd_cal_shift_type_ind)
	csti_err = 1;
      ncal++;
    }
  else if (recd_type == 50)
    {
      dsrms[nsess] = d50.sess_rms - hdr.sess_rms;
      if (debug)
	printf ("sess: %f %f\n", d50.sess_rms, hdr.sess_rms);
      if (hdr.data_qual_ind != d50.data_qual_ind)
	dqi_err = 1;
      nsess++;
    }
  else if (recd_type == 60)
    {
      if (d60.sys_change_ind != hdr.sys_change_ind)
	sci_err = 1;
      if (d60.sys_config_ind != hdr.sys_config_ind)
	scf_err = 1;
    }
  else
    {
      rectype_err++;
      fprintf (str_out1, "Unexpected record type: %d\n", recd_type);
    }
  return (0);
}

/*-------------------------------------------------------------------------
**
**	session_error_report		- Report on differences between
**					  input normal point files and
**                                        assign a "grade".
**
**-----------------------------------------------------------------------*/
void
session_error_report (char crd_in[FNLEN])
{
  int i;

  fprintf (str_out1, "CRD/CSTG Intercomparison Report for session\n");
  fprintf (str_out1, " Date: %4d/%02d/%02d (%3d) %02d:%02d:%02d UTC\n\
 Station: %10s %04d\n\
 Target: %10s %010d %04d %10d\n\n", h4.start_year, h4.start_mon, h4.start_day, hdr.doy, h4.start_hour, h4.start_min, h4.start_sec, h2.stn_name, h2.cdp_pad_id, h3.target_name, h3.ilrs_id, h3.sic, h3.norad);

  fprintf (str_out1, "The following disagreements were found between the\n\
  CRD and CSTG normalpoint files\n");

  if (year_err)
    fprintf (str_out1, " Year disagrees: %2d %2d\n", h4.start_year, hdr.year);
  if (doy_err)
    fprintf (str_out1, " Day of year disagrees: %3d %3d\n", crd_doy, hdr.doy);
  if (padid_err)
    fprintf (str_out1, " CDP pad ids disagree: %4d %4d\n", h2.cdp_pad_id,
	     hdr.cdp_pad_id);
  if (sysnum_err)
    fprintf (str_out1, " CDP system numbers disagree: %2d %2d\n",
	     h2.cdp_sys_num, hdr.cdp_sys_num);
  if (occnum_err)
    fprintf (str_out1, " CDP occupation numbers disagree: %2d %2d\n",
	     h2.cdp_occ_num, hdr.cdp_occ_num);
  if (ndr_err > 0)
    fprintf (str_out1,
	     " CRD Data release disagreed with CSTG normalpoints %d times\n",
	     ndr_err);
  if (wavelength_err)
    fprintf (str_out1, " Wavelengths disagree: %6.1f %6.1f\n",
	     c0.xmit_wavelength, hdr.xmit_wavelength);
  if (timescale_err)
    fprintf (str_out1, " Station time scales disagree: %2d %2d\n",
	     h2.stn_timescale, hdr.stn_timescale);
  if (nwin_ind_err)
    fprintf (str_out1, " Normalpoint window length disagree: %2d %2d\n",
	     hdr.np_window_ind, crd_np_window_ind);
  if (cti_err)
    fprintf (str_out1, " Calibration type indicators disagree: %2d %2d\n",
	     d40.cal_type_ind, crd_cal_type_ind);
  if (csti_err)
    fprintf (str_out1,
	     " Calibration shift type indicators disagree: %2d %2d\n",
	     d40.cal_shift_type_ind, crd_cal_shift_type_ind);
  if (sci_err)
    fprintf (str_out1, " Sytem change indicators disagree: %2d %2d\n",
	     d60.sys_change_ind, hdr.sys_change_ind);
  if (scf_err)
    fprintf (str_out1, " Sytem configuration indicators disagree: %2d %2d\n",
	     d60.sys_config_ind, hdr.sys_config_ind);
  if (dqi_err)
    fprintf (str_out1, " Data quality indicators disagree: %2d %2d\n",
	     d50.data_qual_ind, hdr.data_qual_ind);
  if (iid_err)
    fprintf (str_out1, " ILRS target IDs disagree: %10d %10d\n", h3.ilrs_id,
	     hdr.ilrs_id);
  if (rectype_err)
    fprintf (str_out1, " Unexpected record type %d times.\n", rectype_err);

  if (np_missing > 0)
    fprintf
      (str_out1,
       "There were %d CSTG normalpoint missing relative to the CRD file.\n",
       np_missing);
  if (np_missing > 3)
    fail++;
  else if (np_missing == 3)
    marginal++;
  else if (np_missing == 2)
    pass++;
  else
    strongpass++;
  if (d11_missing > 0)
    fprintf
      (str_out1,
       "There were %d CRD normalpoint missing relative to the CSTG file.\n",
       d11_missing);
  if (d11_missing > 3)
    fail++;
  else if (d11_missing == 3)
    marginal++;
  else if (d11_missing == 2)
    pass++;
  else
    strongpass++;

  int nok, nroundoff, nfurtheroff, nwayoff;
  nok = nroundoff = nfurtheroff = nwayoff = 0;
  for (i = 0; i < npth; i++)
    {
      if (fabs (dpres[i]) < 0.1)
	nok++;
      else if (fabs (dpres[i]) < 1.0)
	nroundoff++;
      else if (fabs (dpres[i]) < 10.0)
	nfurtheroff++;
      else
	nwayoff++;
    }
  if (nwayoff)
    fail++;
  else if (nfurtheroff)
    marginal++;
  else if (nroundoff)
    pass++;
  else
    strongpass++;
  fprintf (str_out1, " Of %d pressure measurements,\n\
   %d differed by < 0.1 mb;\n\
   %d differed by < 1.0 mb;\n\
   %d differed by < 10 mb; and\n\
   %d differed by more.\n", npth, nok, nroundoff, nfurtheroff, nwayoff);

  nok = nroundoff = nfurtheroff = nwayoff = 0;
  for (i = 0; i < npth; i++)
    {
      if (fabs (dtemp[i]) < 0.1)
	nok++;
      else if (fabs (dtemp[i]) < 1.0)
	nroundoff++;
      else if (fabs (dtemp[i]) < 10.0)
	nfurtheroff++;
      else
	nwayoff++;
    }
  if (nwayoff)
    fail++;
  else if (nfurtheroff)
    marginal++;
  else if (nroundoff)
    pass++;
  else
    strongpass++;
  fprintf (str_out1, " Of %d temperature measurements,\n\
   %d differed by < 0.1 K;\n\
   %d differed by < 1.0 K;\n\
   %d differed by < 10.0 K; and\n\
   %d differed by more.\n", npth, nok, nroundoff, nfurtheroff, nwayoff);

  nok = nroundoff = nfurtheroff = nwayoff = 0;
  for (i = 0; i < npth; i++)
    {
      if (fabs (dhum[i]) < 1)
	nok++;
      else if (fabs (dhum[i]) < 5)
	nroundoff++;
      else if (fabs (dhum[i]) < 10)
	nfurtheroff++;
      else
	nwayoff++;
    }
  if (nwayoff)
    fail++;
  else if (nfurtheroff)
    marginal++;
  else if (nroundoff)
    pass++;
  else
    strongpass++;
  fprintf (str_out1, " Of %d humidity measurements,\n\
   %d differed by < 1 %;\n\
   %d differed by < 5 %;\n\
   %d differed by < 10 %; and\n\
   %d differed by more.\n", npth, nok, nroundoff, nfurtheroff, nwayoff);

  nok = nroundoff = nfurtheroff = nwayoff = 0;
  for (i = 0; i < nfr; i++)
    {
      if (fabs (dft[i]) < 1.e-7)
	nok++;
      else if (fabs (dft[i]) < 5.e-7)
	nroundoff++;
      else if (fabs (dft[i]) < 1.e-6)
	nfurtheroff++;
      else
	nwayoff++;
    }
  /*if (nwayoff) */
  if (nwayoff > nfr * toler)
    fail++;
  else if (nfurtheroff)
    marginal++;
  else if (nroundoff)
    pass++;
  else
    strongpass++;
  fprintf (str_out1, " Of %d normal point seconds of day,\n\
   %d differed by < 0.1 psec ;\n\
   %d differed by < 500 nsec;\n\
   %d differed by < 1 microsec; and\n\
   %d differed by more.\n", nfr, nok, nroundoff, nfurtheroff, nwayoff);

  nok = nroundoff = nfurtheroff = nwayoff = 0;
  for (i = 0; i < nfr; i++)
    {
      if (fabs (dr[i]) < 1.e-12)
	nok++;
      else if (fabs (dr[i]) < 5.e-12)
	nroundoff++;
      else if (fabs (dr[i]) < 1.e-11)
	nfurtheroff++;
      else
	nwayoff++;
    }
  /*if (nwayoff) */
  if (nwayoff > nfr * toler)
    fail++;
  else if (nfurtheroff)
    marginal++;
  else if (nroundoff)
    pass++;
  else
    strongpass++;
  fprintf (str_out1, " Of %d normal point time of flight,\n\
   %d differed by < 1 psec ;\n\
   %d differed by < 5 psec;\n\
   %d differed by < 10 psec; and\n\
   %d differed by more.\n", nfr, nok, nroundoff, nfurtheroff, nwayoff);

  nok = nroundoff = nfurtheroff = nwayoff = 0;
  for (i = 0; i < nfr; i++)
    {
      if (fabs (dbrms[i]) < 1)
	nok++;
      else if (fabs (dbrms[i]) < 5.)
	nroundoff++;
      else if (fabs (dbrms[i]) < 10.)
	nfurtheroff++;
      else
	nwayoff++;
    }
  /*if (nwayoff) */
  if (nwayoff > nfr * toler)
    fail++;
  else if (nfurtheroff)
    marginal++;
  else if (nroundoff)
    pass++;
  else
    strongpass++;
  fprintf (str_out1, " Of %d normal point bin rms,\n\
   %d differed by < 1 psec ;\n\
   %d differed by < 5 psec;\n\
   %d differed by < 10 psec; and\n\
   %d differed by more.\n", nfr, nok, nroundoff, nfurtheroff, nwayoff);

  nok = nroundoff = nfurtheroff = nwayoff = 0;
  for (i = 0; i < nfr; i++)
    {
      if (abs (dnret[i]) < 1)
	nok++;
      else if (abs (dnret[i]) < 5)
	nroundoff++;
      else if (abs (dnret[i]) < 10)
	nfurtheroff++;
      else
	nwayoff++;
    }
  if (nwayoff)
    fail++;
  else if (nfurtheroff)
    marginal++;
  else if (nroundoff)
    pass++;
  else
    strongpass++;
  fprintf (str_out1, " Of %d normalpoints, the number of returns,\n\
   %d differed by < 1;\n\
   %d differed by < 5;\n\
   %d differed by < 10; and\n\
   %d differed by more.\n", nfr, nok, nroundoff, nfurtheroff, nwayoff);

  nok = nroundoff = nfurtheroff = nwayoff = 0;
  for (i = 0; i < ncal; i++)
    {
      if (fabs (dcsd[i]) < 1)
	nok++;
      else if (fabs (dcsd[i]) < 5.)
	nroundoff++;
      else if (fabs (dcsd[i]) < 10.)
	nfurtheroff++;
      else
	nwayoff++;
    }
  if (nwayoff)
    fail++;
  else if (nfurtheroff)
    marginal++;
  else if (nroundoff)
    pass++;
  else
    strongpass++;
  fprintf (str_out1, " Of %d calibration system delays,\n\
   %d differed by < 1 psec ;\n\
   %d differed by < 5 psec;\n\
   %d differed by < 10 psec; and\n\
   %d differed by more.\n", ncal, nok, nroundoff, nfurtheroff, nwayoff);

  nok = nroundoff = nfurtheroff = nwayoff = 0;
  for (i = 0; i < ncal; i++)
    {
      if (fabs (dcds[i]) < 1)
	nok++;
      else if (fabs (dcds[i]) < 5.)
	nroundoff++;
      else if (fabs (dcds[i]) < 10.)
	nfurtheroff++;
      else
	nwayoff++;
    }
  if (nwayoff)
    fail++;
  else if (nfurtheroff)
    marginal++;
  else if (nroundoff)
    pass++;
  else
    strongpass++;
  fprintf (str_out1, " Of %d calibration delay shifts,\n\
   %d differed by < 1 psec ;\n\
   %d differed by < 5 psec;\n\
   %d differed by < 10 psec; and\n\
   %d differed by more.\n", ncal, nok, nroundoff, nfurtheroff, nwayoff);

  nok = nroundoff = nfurtheroff = nwayoff = 0;
  for (i = 0; i < ncal; i++)
    {
      if (fabs (dcrms[i]) < 1)
	nok++;
      else if (fabs (dcrms[i]) < 5.)
	nroundoff++;
      else if (fabs (dcrms[i]) < 10.)
	nfurtheroff++;
      else
	nwayoff++;
    }
  if (nwayoff)
    fail++;
  else if (nfurtheroff)
    marginal++;
  else if (nroundoff)
    pass++;
  else
    strongpass++;
  fprintf (str_out1, " Of %d calibration rms,\n\
   %d differed by < 1 psec ;\n\
   %d differed by < 5 psec;\n\
   %d differed by < 10 psec; and\n\
   %d differed by more.\n", ncal, nok, nroundoff, nfurtheroff, nwayoff);

  nok = nroundoff = nfurtheroff = nwayoff = 0;
  for (i = 0; i < nsess; i++)
    {
      if (fabs (dsrms[i]) < 1.)
	nok++;
      else if (fabs (dsrms[i]) < 5.)
	nroundoff++;
      else if (fabs (dsrms[i]) < 10)
	nfurtheroff++;
      else
	nwayoff++;
    }
  if (nwayoff)
    fail++;
  else if (nfurtheroff)
    marginal++;
  else if (nroundoff)
    pass++;
  else
    strongpass++;
  fprintf (str_out1, " Of %d session rms,\n\
   %d differed by < 1 psec ;\n\
   %d differed by < 5 psec;\n\
   %d differed by < 10 psec; and\n\
   %d differed by more.\n", nsess, nok, nroundoff, nfurtheroff, nwayoff);

  if (fail || hdr_err)
    {
      fprintf (str_out1, "\n%s: failed\n", crd_in);
      if (str_out2)
	fprintf (str_out2, "%s: failed\n", crd_in);
    }
  else if (marginal)
    {
      fprintf (str_out1, "\n%s: marginal\n", crd_in);
      if (str_out2)
	fprintf (str_out2, "%s: marginal\n", crd_in);
    }
  else if (pass)
    {
      fprintf (str_out1, "\n%s: pass\n", crd_in);
      if (str_out2)
	fprintf (str_out2, "%s: pass\n", crd_in);
      return;
    }
  else if (strongpass)
    {
      fprintf (str_out1, "\n%s: strong pass\n", crd_in);
      if (str_out2)
	fprintf (str_out2, "%s: strong pass\n", crd_in);
      return;
    }
  else
    {
      fprintf (str_out1, "\n%s: internal program error\n", crd_in);
      if (str_out2)
	fprintf (str_out2, "%s: internal program error\n", crd_in);
    }
}

/*-------------------------------------------------------------------------
**
**	setup_files		- Open input and output file names
**
**-----------------------------------------------------------------------*/
void
get_win_ind (int np_window_length, int ilrs_id, int *nwi, int *lnwi)
{

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
setup_files (int argc, char *argv[], char crd_in[FNLEN])
{
  char cstg_in[FNLEN] = { "" };
  char report_out[FNLEN] = { "" };
  char summary_out[FNLEN] = { "" };

  int found = 0, i;

/*  Get file names  */
  for (i = 1; i < argc; i += 2)
    {
      if (argv[i][0] == '-')
	{
	  if (argv[i][1] == '1')
	    {
	      strcpy (crd_in, argv[i + 1]);
	      printf ("crd_in = [%s]\n", crd_in);
	      found += 1;
	    }
	  else if (argv[i][1] == '2')
	    {
	      strcpy (cstg_in, argv[i + 1]);
	      printf ("cstg_in = [%s]\n", cstg_in);
	      found += 2;
	    }
	  else if (argv[i][1] == 'o')
	    {
	      if (found >= 1)
		{
		  strcpy (report_out, argv[i + 1]);
		}
	      printf ("report_out = [%s]\n", report_out);
	      found += 3;
	    }
	  else if (argv[i][1] == 's')
	    {
	      if (found >= 1)
		{
		  strcpy (summary_out, argv[i + 1]);
		}
	      printf ("summary_out = [%s]\n", summary_out);
	      found += 4;
	    }
	  else if (argv[i][1] == 'v')
	    {
	      printf ("crd_cstg_np_cmp version: %s\n", pvern);
	      exit (1);
	    }
	}
    }
  if (found != 6 && found != 10)
    {
      printf
	("Usage: crd_cstg_np_cmp -1 crd_npt -2 cstg_npt -o cmp_report [-s cmp_summary] [-v]\n");
      exit (1);
    }

/*  open input crd file  */
  if ((str_in1 = fopen (crd_in, "r")) == NULL)
    {
      printf ("Could not open file %s\n", crd_in);
      exit (1);
    }

/*  open input cstg file  */
  if ((str_in2 = fopen (cstg_in, "r")) == NULL)
    {
      printf ("Could not open file %s\n", cstg_in);
      exit (1);
    }

/*  open output report file  */
  if ((str_out1 = fopen (report_out, "w")) == NULL)
    {
      printf ("Could not open file %s\n", report_out);
      exit (1);
    }

/*  open output summary file  */
  if (strlen (summary_out) > 0)
    {
      if ((str_out2 = fopen (summary_out, "a")) == NULL)
	{
	  printf ("Could not open file %s\n", summary_out);
	  exit (1);
	}
    }
}

/*-------------------------------------------------------------------------
**
**	end_it		- Program termination routine
**			  (Could just use exit(), but intend to add more
**                        information.)
**
**-----------------------------------------------------------------------*/
void
end_it (int n)
{
  exit (n);
}

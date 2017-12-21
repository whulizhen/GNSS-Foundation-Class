#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

#define FNLEN   256
#define MAXRS 10000
#define MAXPTH 10000
#define MAXFR 10000
#define VLIGHT 2.99792458e+08

#include "crd.h"
#include "merit.h"

FILE *str_in1, *str_in2, *str_out1, *str_out2;
struct merit_fr mrt;
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
int d10_missing = 0;
int mrt_missing = 0;

double check_frd ();
void session_error_report ();
void get_win_ind ();
void setup_files ();
void end_it ();

/*-------------------------------------------------------------------------
 * Program: crd_cstg_np_cmp
 *
 * Purpose:
 * Compares CRD and old-style MERIT-II full rate files
 *
 * Calling sequence:
 *   crd_merit_fr_cmp -1 crd_npt_file -2 merit_fr_file -o report_file \
 *       [-s summary_file] [-v]
 *     -1 is the CRD-formatted normalpoint file
 *     -2 is the merit-formatted (old format) full rate file
 *     -o output report file
 *     -s one line summary is appended to this file
 *     -v prints program version number and stops.
 *
 * Author: Randall Ricklefs / Univ of Texas / Center for Space Research
 *
 * History:
 *   28 May 2008 - initial version
 *
**-----------------------------------------------------------------------*/

int
main (argc, argv)
     int argc;
     char *argv[];
{
  char str[256];
  char crd_in[FNLEN];
  int crd_ang_read = 0;         /* Has a pointing angle record been read? */
  int crd_cal_read = 0;         /* Has a calibration record been read? */
  int crd_cfg_read = 0;         /* Has a sys configuration record been read? */
  int crd_met_read = 0;         /* Has a meteorological record been read? */
  int crd_rngsupp_read = 0;     /* Has a range supplement record been read? */
  int merit_fr_read = 0;        /* Has a merit full rate record been read? */
  int status;
  double check_status = 0;

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

	  /* read the rest of the merit full rate records */
	  do
	    {
	      if ((status = fgets (str, 256, str_in2)) != NULL)
		d10_missing++;
	      if (debug)
		printf ("mrt [%s] status %d d10_missing %d\n", str,
			status, d10_missing);
	    }
	  while (status != NULL);
	  if (debug)
	    printf ("read h8\n");
	  session_error_report (crd_in);
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
	}
      else if (strncmp (str, "10", 2) == 0)
	{
          /* Ranging record */
	  read_10 (str, &d10);

	  /* If the CRD full rate time is earlier than the merit,
	   * skip the record */
	  do
	    {
	      if (check_status >= 0)	/* play catch-up */
		{
		  if ((status = fgets (str, 256, str_in2)) == NULL)
		    {
		      mrt_missing++;
		      if (debug)
			printf ("mrt missing (10)\n");
		      break;
		    }
		  if (debug)
		    printf ("10 read mrt [%s]\n", str);
		  read_merit_fr (str, &mrt);
		  if (!merit_fr_read)
		    {
		      if (crd_rngsupp_read)
			check_frd (12);
		      if (crd_met_read)
			check_frd (20);
		      if (crd_ang_read)
			check_frd (30);
		      if (crd_cal_read)
			check_frd (40);
		      if (crd_cfg_read)
			check_frd (60);
		    }
		  merit_fr_read = 1;
		}
	    }
	  while ((check_status = check_frd (10)) > 0);
	}
      else if (strncmp (str, "12", 2) == 0)
	{
          /* Range supplement record */
	  read_12 (str, &d12);
	  if (merit_fr_read)
	    check_frd (12);
	  crd_rngsupp_read = 1;
	}
      else if (strncmp (str, "20", 2) == 0)
	{
          /* Meteorological record */
	  read_20 (str, &d20);
	  if (merit_fr_read)
	    check_frd (20);
	  crd_met_read = 1;
	}
      else if (strncmp (str, "30", 2) == 0)
	{
          /* Pointing angle record */
	  read_30 (str, &d30);
	  if (merit_fr_read)
	    check_frd (30);
	  crd_ang_read = 1;
	}
      else if (strncmp (str, "40", 2) == 0)
	{
          /* Calibration record */
	  read_40 (str, &d40);
	  if (merit_fr_read)
	    check_frd (40);
	  crd_cal_read = 1;
	}
      else if (strncmp (str, "50", 2) == 0)
	{
          /* Session summary record */
	  read_50 (str, &d50);
	  check_frd (50);
	}
      else if (strncmp (str, "60", 2) == 0)
	{
          /* Compatibility record */
	  read_60 (str, &d60);
	  if (merit_fr_read)
	    check_frd (60);
	  crd_cfg_read = 1;
	}
    }
}

double drefcorr[MAXRS];
double dCoMcorr[MAXRS];
double dpres[MAXPTH];
double dtemp[MAXPTH];
double dhum[MAXPTH];
double dft[MAXFR];
double dr[MAXFR];
double dxa[MAXFR];
double dcsd[50];
double dcds[50];
double dcrms[50];
double dsrms[50];
double toler = 0.15;
double dftb;
int ndr_err = 0;
int nrs = 0;
int npth = 0, nfr = 0;
int ncal = 0;
int nsess = 0;
int nwin_ind_err = 0;
int cti_err = 0;
int csti_err = 0;
int crd_cal_type_ind, crd_cal_shift_type_ind;
int sch_err = 0;
int sci_err = 0;
int iid_err = 0;
int rectype_err = 0;
int year_err = 0;
int doy_err = 0;
int padid_err = 0;
int sysnum_err = 0;
int occnum_err = 0;
int wavelength_err = 0;
int timescale_err = 0;
int aoi_err = 0;
int rai_err = 0;
int CoMai_err = 0;
int xaai_err = 0;
int ind_err = 0;
int fail = 0;
int marginal = 0;
int pass = 0;
int strongpass = 0;
int crd_doy;

/*-------------------------------------------------------------------------
**
**      check_frd               - Check contents of merit full rate record
**                                against CRD records
**
**-----------------------------------------------------------------------*/
double
check_frd (int recd_type)
{
  int crd_cal_type_ind;
  int crd_cal_shift_type_ind;
  int mrt_data_release;

  if (recd_type == 10)
    {

      if (mrt.ilrs_id != h3.ilrs_id)
	iid_err = 1;
      if (mrt.year != h4.start_year)
	year_err = 1;
      grtodoy (h4.start_year % 100, h4.start_mon, h4.start_day, &crd_doy);
      if (mrt.doy != crd_doy)
	doy_err = 1;
      if (mrt.cdp_pad_id != h2.cdp_pad_id)
	padid_err = 1;
      if (mrt.cdp_sys_num != h2.cdp_sys_num)
	sysnum_err = 1;
      if (mrt.cdp_occ_num != h2.cdp_occ_num)
	occnum_err = 1;
      if (mrt.xmit_wavelength != c0.xmit_wavelength)
	wavelength_err = 1;
      if (mrt.stn_timescale != h2.stn_timescale)
	timescale_err = 1;
      if (mrt.refraction_app_ind * h4.refraction_app_ind != 0 ||
	  mrt.refraction_app_ind + h4.refraction_app_ind != 1)
	rai_err = 1;
      if (mrt.CofM_app_ind * h4.CofM_app_ind != 0 ||
	  mrt.CofM_app_ind + h4.CofM_app_ind != 1)
	CoMai_err = 1;
      if (mrt.xcv_amp_app_ind * h4.xcv_amp_app_ind != 0 ||
	  mrt.xcv_amp_app_ind + h4.xcv_amp_app_ind != 1)
	xaai_err = 1;
      ind_err += rectype_err + iid_err + year_err + doy_err + padid_err +
	sysnum_err + occnum_err + wavelength_err + timescale_err +
	sch_err + sci_err + cti_err + csti_err + rai_err + CoMai_err +
	xaai_err;

      /* Are both firing times close to the same time? */
      if (debug)
	printf ("%Lf %Lf %d\n", d10.sec_of_day, mrt.sec_of_day);
      if (debug)
	printf ("bins: %f %f\n", d10.sec_of_day, mrt.sec_of_day);
      /* No */
      if (debug)
	printf ("sod %.13Lf %.13Lf dftb = %.13Lf\n",
		d10.sec_of_day, mrt.sec_of_day,
		d10.sec_of_day - mrt.sec_of_day);
      if (fabs (d10.sec_of_day - mrt.sec_of_day) > 1.e-5)
	{
	  dftb = d10.sec_of_day - mrt.sec_of_day;
	  if (debug)
	    printf ("dftb= %f\n", dftb);
	  if (dftb < 0)
	    {
	      mrt_missing++;
	      if (debug)
		printf ("mrt_missing (chkfrd)\n");
	    }
	  else
	    {
	      d10_missing++;
	      if (debug)
		printf ("d10_missing (chkfrd)\n");
	    }
	  return (dftb);
	}
      /* Yes */
      dft[nfr] = d10.sec_of_day - mrt.sec_of_day;
      if (debug)
	printf ("ft %.12Lf %.12Lf\n", d10.sec_of_day, mrt.sec_of_day);
      dr[nfr] = d10.time_of_flight - mrt.time_of_flight;
      if (debug)
	printf ("r %.12Lf %.12Lf\n", d10.time_of_flight, mrt.time_of_flight);
      dxa[nfr] = d10.xcv_amp - mrt.xcv_amp;
      if (h4.data_release != mrt.data_release - 1)
	ndr_err++;
      nfr++;

    }
  else if (recd_type == 12)
    {
      drefcorr[nrs] = d12.refraction_corr - mrt.refraction_corr;
      dCoMcorr[nrs] =
	d12.target_CofM_corr - mrt.target_CofM_corr * VLIGHT / 2.e9;
      nrs++;
    }
  else if (recd_type == 20)
    {
      /* check dtime?? What if no np read? */
      if (debug)
	printf ("pth: %f %f    %f %f   %f %f\n", d20.pressure, mrt.pressure,
		d20.temperature, mrt.temperature, d20.humidity, mrt.humidity);
      dpres[npth] = d20.pressure - mrt.pressure;
      dtemp[npth] = d20.temperature - mrt.temperature;
      dhum[npth] = d20.humidity - mrt.humidity;
      npth++;
    }
  else if (recd_type == 30)
    {
      if (mrt.angle_origin_ind != d30.angle_origin_ind)
	aoi_err = 1;
    }
  else if (recd_type == 40)
    {
      dcsd[ncal] = d40.cal_sys_delay - mrt.cal_sys_delay;
      dcds[ncal] = d40.cal_delay_shift - mrt.cal_delay_shift;
      dcrms[ncal] = d40.cal_rms - mrt.cal_rms;
      if (debug)
	printf ("cal: %f %f\n", d40.cal_rms, mrt.cal_rms);

      crd_cal_shift_type_ind = mrt.cal_type_ind > 4 ? 3 : 2;
      crd_cal_type_ind = mrt.cal_type_ind;
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
      dsrms[nsess] = d50.sess_rms - mrt.sess_rms;
      if (debug)
	printf ("sess: %f %f\n", d50.sess_rms, mrt.sess_rms);
      nsess++;
    }
  else if (recd_type == 60)
    {
      if (d60.sys_change_ind != mrt.sys_change_ind)
	sch_err = 1;
      if (d60.sys_config_ind != mrt.sys_config_ind)
	sci_err = 1;
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
**      session_error_report            - Report on differences between
**                                        input full rate files and
**                                        assign a "grade".
**
**-----------------------------------------------------------------------*/
void
session_error_report (char crd_in[FNLEN])
{
  int i;
  int nok, nroundoff, nfurtheroff, nwayoff;

  fprintf (str_out1,
	   "CRD/Merit-II Data Intercomparison Report for session\n");
  fprintf (str_out1, " Date: %4d/%02d/%02d (%3d) %02d:%02d:%02d UTC\n\
 Station: %10s %04d\n\
 Target: %10s %010d %04d %10d\n\n", h4.start_year, h4.start_mon, h4.start_day, mrt.doy, h4.start_hour, h4.start_min, h4.start_sec, h2.stn_name, h2.cdp_pad_id, h3.target_name, h3.ilrs_id, h3.sic, h3.norad);

  fprintf (str_out1, "The following disagreements were found between the\n\
  CRD and MERIT-II full rate files\n");

  if (year_err)
    fprintf (str_out1, " Year disagrees: %2d %2d\n", h4.start_year, mrt.year);
  if (doy_err)
    fprintf (str_out1, " Day of year disagrees: %3d %3d\n", crd_doy, mrt.doy);
  if (padid_err)
    fprintf (str_out1, " CDP pad ids disagree: %4d %4d\n", h2.cdp_pad_id,
	     mrt.cdp_pad_id);
  if (sysnum_err)
    fprintf (str_out1, " CDP system numbers disagree: %2d %2d\n",
	     h2.cdp_sys_num, mrt.cdp_sys_num);
  if (occnum_err)
    fprintf (str_out1, " CDP occupation numbers disagree: %2d %2d\n",
	     h2.cdp_occ_num, mrt.cdp_occ_num);
  if (ndr_err > 0)
    fprintf (str_out1,
	     " CRD Data release disagreed with MERIT full rate %d times\n",
	     ndr_err);
  if (wavelength_err)
    fprintf (str_out1, " Wavelengths disagree: %6.1f %6.1f\n",
	     c0.xmit_wavelength, mrt.xmit_wavelength);
  if (timescale_err)
    fprintf (str_out1, " Station time scales disagree: %2d %2d\n",
	     h2.stn_timescale, mrt.stn_timescale);
  if (cti_err)
    fprintf (str_out1, " Calibration type indicators disagree: %2d %2d\n",
	     d40.cal_type_ind, crd_cal_type_ind);
  if (csti_err)
    fprintf (str_out1,
	     " Calibration shift type indicators disagree: %2d %2d\n",
	     d40.cal_shift_type_ind, crd_cal_shift_type_ind);
  if (sch_err)
    fprintf (str_out1, " Sytem change indicators disagree: %2d %2d\n",
	     d60.sys_change_ind, mrt.sys_change_ind);
  if (sci_err)
    fprintf (str_out1, " Sytem configuration indicators disagree: %2d %2d\n",
	     d60.sys_config_ind, mrt.sys_config_ind);
  if (iid_err)
    fprintf (str_out1, " ILRS target IDs disagree: %10d %10d\n", h3.ilrs_id,
	     mrt.ilrs_id);
  if (rectype_err)
    fprintf (str_out1, " Unexpected record type %d times.\n", rectype_err);

  if (mrt_missing > 0)
    fprintf
      (str_out1,
       "There were %d MERIT full rate records missing relative to the CRD file.\n",
       mrt_missing);
  if (mrt_missing > 3)
    fail++;
  else if (mrt_missing == 3)
    marginal++;
  else if (mrt_missing == 2)
    pass++;
  else
    strongpass++;
  if (d10_missing > 0)
    fprintf
      (str_out1,
       "There were %d CRD full rate records missing relative to the MERIT file.\n",
       d10_missing);
  if (d10_missing > 3)
    fail++;
  else if (d10_missing == 3)
    marginal++;
  else if (d10_missing == 2)
    pass++;
  else
    strongpass++;

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
  fprintf (str_out1, " Of %d full rate seconds of day,\n\
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
  fprintf (str_out1, " Of %d full rate time of flight,\n\
   %d differed by < 1 psec ;\n\
   %d differed by < 5 psec;\n\
   %d differed by < 10 psec; and\n\
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

  if (h4.xcv_amp_app_ind)
    {
      nok = nroundoff = nfurtheroff = nwayoff = 0;
      for (i = 0; i < nfr; i++)
	{
	  if (fabs (dxa[i]) < 1.e-7)
	    nok++;
	  else if (fabs (dxa[i]) < 5.e-7)
	    nroundoff++;
	  else if (fabs (dxa[i]) < 1.e-6)
	    nfurtheroff++;
	  else
	    nwayoff++;
	}
      /*if (nwayoff)
         fail++;
         else if (nfurtheroff)
         marginal++;
         else if (nroundoff)
         pass++;
         else
         strongpass++; */
      fprintf (str_out1, " Of %d receiver amplitudes,\n\
   %d differed by < 1;\n\
   %d differed by < 5;\n\
   %d differed by < 10; and\n\
   %d differed by more.\n", nfr, nok, nroundoff, nfurtheroff, nwayoff);
    }

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

  /*if (refraction_corr_ind) */
  if (h4.refraction_app_ind)
    {
      nok = nroundoff = nfurtheroff = nwayoff = 0;
      for (i = 0; i < nrs; i++)
	{
	  if (fabs (drefcorr[i]) < 0.1)
	    nok++;
	  else if (fabs (drefcorr[i]) < 1.0)
	    nroundoff++;
	  else if (fabs (drefcorr[i]) < 10.0)
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
      fprintf (str_out1, " Of %d refraction corrections,\n\
   %d differed by < 0.1 ps;\n\
   %d differed by < 1.0 ps;\n\
   %d differed by < 10 ps; and\n\
   %d differed by more.\n", nrs, nok, nroundoff, nfurtheroff, nwayoff);
    }

  if (h4.CofM_app_ind)
    {
      nok = nroundoff = nfurtheroff = nwayoff = 0;
      for (i = 0; i < nrs; i++)
	{
	  if (fabs (dCoMcorr[i]) < 0.1)
	    nok++;
	  else if (fabs (dCoMcorr[i]) < 1.0)
	    nroundoff++;
	  else if (fabs (dCoMcorr[i]) < 10.0)
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
      fprintf (str_out1, " Of %d center of mass corrections,\n\
   %d differed by < 0.1 m;\n\
   %d differed by < 1.0 m;\n\
   %d differed by < 10 m; and\n\
   %d differed by more.\n", nrs, nok, nroundoff, nfurtheroff, nwayoff);
    }

  if (fail || ind_err)
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
setup_files (int argc, char *argv[], char crd_in[FNLEN])
{
  char merit_in[FNLEN] = { "" };
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
	      printf ("crd_in =   [%s]\n", crd_in);
	      found += 1;
	    }
	  else if (argv[i][1] == '2')
	    {
	      strcpy (merit_in, argv[i + 1]);
	      printf ("merit_in = [%s]\n", merit_in);
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
	("Usage: crd_merit_fr_cmp -1 crd_fr -2 merit_fr -o cmp_report [-s cmp_summary] [-v]\n");
      exit (1);
    }

/*  open input crd file  */
  if ((str_in1 = fopen (crd_in, "r")) == NULL)
    {
      printf ("Could not open file %s\n", crd_in);
      exit (1);
    }

/*  open input merit file  */
  if ((str_in2 = fopen (merit_in, "r")) == NULL)
    {
      printf ("Could not open file %s\n", merit_in);
      exit (1);
    }

/*  open output report file  */
  if ((str_out1 = fopen (report_out, "w")) == NULL)
    {
      printf ("Could not open file %s\n", report_out);
      exit (1);
    }

/*  open output summary file  */
  str_out2 = 0;
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
**      end_it          - Program termination routine
**                        (Could just use exit(), but intend to add more
**                        information.)
**
**-----------------------------------------------------------------------*/
void
end_it (int n)
{
  exit (n);
}

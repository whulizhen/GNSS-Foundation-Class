#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "crd.h"

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

char *strcasestr ();
int getifield ();
int read_c0 ();
int read_c1 ();
int read_c2 ();
int read_c3 ();
int read_c4 ();
int read_c0 ();
int read_10 ();
int read_11 ();
int read_12 ();
int read_20 ();
int read_21 ();
int read_30 ();
int read_40 ();
int read_50 ();
int read_60 ();

FILE *str_in;
int get_sat_ids ();

/*-------------------------------------------------------------------------
 * Program: crd_chk
 *
 * Purpose:
 * Reads and tests a CRD file for compliance with the format.
 *
 * Calling sequence:
 *   crd_chk [-b][-s] crd_filename
 *     -b gives a brief listing with errors only.
 *     -s write out-of-sequence records (NOT YET IMPLEMENTED. SET fout_of_seq)
 *     -x extended out-of-sequence check, below 1 msec 
 *		(NOT YET IMPLEMENTED. SET fex_seq)
 *
 * Author: Randall Ricklefs / Univ of Texas / Center for Space Research
 *
 * History:
 *   Nov 06, 2007 - Initial version
 *   Mar 10, 2009 - Correct column from which h4.data_qual_alert_ind is read. rlr.
 *   May 04, 2009 - Prevent "out of sequence" error when seconds of day wraps
 *                  around to zero. This effects all record types with second
 *                  of day fields. rlr.
 *   Nov 09, 2010 - Check sat label, norad, and ILRS IDs again data base. 
 *                  Tests are relative to ILRS ID, since it is in old 
 *                  normal point format. rlr.
 *   Jan 20, 2011 - Print sat ID error with the pass rather than at the end 
 *                  of processing. Check for out-of-sequence errors only within
 *                  the same record type. rlr.
 *   Jun 22, 2012 - Add "fail" code to provide non-zero return value when 
 *   		    program fails due to an error. Also, fail when there are
 *   		    fewer h8 records than h1 or h2 records. rlr.
 *
**-----------------------------------------------------------------------*/
main (argc, argv)
     int argc;
     char *argv[];
{
  int i;
  int fail = 0;
  int ftab = 0, nstat;
  int nh1 = 0, nh2 = 0, nh3 = 0, nh4 = 0, nh8 = 0, nh9 = 0, n00 = 0, n10 =
    0, n11 = 0, n12 = 0;
  int n20 = 0, n21 = 0, n30 = 0, n40 = 0, n50 = 0, n60 =
    0, n9x = 0;
  int nc0 = 0, nc1 = 0, nc2 = 0, nc3 = 0, nc4 = 0;
  int nh1_bl[10], nh1_blerr = 0, nh2_bl[25], nh2_blerr =
    0, nh3_bl[10], nh3_blerr = 0;
  int nh4_bl[5], nh4_blerr = 0, nh8_bl[5], nh8_blerr = 0;
  int nc0low = 0, nc1low = 0, nc2low = 0, nc3low = 0, nc4low = 0;
  int n10bad = 0, n11bad = 0, n12bad = 0, n20bad = 0, n21bad = 0, n22bad = 0;
  int n30bad = 0, n40bad = 0;
  int n50bad = 0, n60bad = 0, nc0bad = 0, nc2bad = 0, nc3bad = 0, nc4bad =
    0, nc1bad = 0;
  int n10low = 0, n11low = 0, n12low = 0, n20low = 0, n21low = 0, n30low = 0;
  int n40low = 0, n50low = 0, n60low = 0;
  int nseq = 0, tnseq = 0;
  int nc0dterr = 0, nc0xwerr = 0;
  int nc1dterr = 0, nc1pwerr = 0, nc1nfrerr = 0, nc1peerr = 0, nc1pwderr =
    0, nc1bderr = 0, nc1piserr = 0;
  int nc2dterr = 0, nc2awerr = 0, nc2qeerr = 0, nc2verr = 0, nc2dcerr =
    0, nc2opwerr = 0, nc2sferr = 0, nc2sfxerr = 0, nc2stferr = 0;
  int nc3dterr = 0, nc3edcerr = 0;
  int nc4dterr = 0, nc4sodaierr = 0, nc4scodaierr = 0, nc4sctsierr = 0;
  int nd10soderr = 0, nd10toferr = 0, nd10eeerr = 0, nd10fferr =
    0, nd10dcerr = 0, nd10snerr = 0, nd10xaerr = 0;
  int nd11soderr = 0, nd11toferr = 0, nd11eeerr = 0, nd11nwlerr =
    0, nd11nrerr = 0, nd11brerr = 0, nd11rrerr = 0;
  int nd12soderr = 0;
  int nd20soderr = 0, nd20preserr = 0, nd20temperr = 0, nd20humerr =
    0, nd20voerr = 0;
  int nd21soderr = 0, nd21wserr = 0, nd21wderr = 0, nd21viserr =
    0, nd21scerr = 0, nd21aserr = 0, nd21ccerr = 0;
  int nd30soderr = 0, nd30azerr = 0, nd30elerr = 0, nd30dierr =
    0, nd30aoierr = 0, nd30rcierr = 0;
  int nd40soderr = 0, nd40toderr = 0, nd40nprerr = 0, nd40npuerr =
    0, nd40owtderr = 0, nd40csderr = 0, nd40cdserr = 0, nd40crerr =
    0, nd40ctierr = 0, nd40cstierr = 0;
  int nd50srerr = 0, nd50dqierr = 0;
  int nd60scherr = 0, nd60scierr = 0;
  int fex_seq = 0;
  int fout_of_seq = 1;
  int red_flag_only = 0;
  int data_type, refraction_app_ind, CofM_app_ind, xcv_amp_app_ind;
  int stn_sysdelay_app_ind, SC_sysdelay_app_ind, range_type_ind, target_type;
  int data_qual_alert_ind;
  int tsic, tnorad, ttarget_type;
  int nh3sicerr = 0, nh3noraderr = 0, nh3ttypeerr = 0, nh3nameerr = 0;
  int nh3iderr = 0;
  long double last_10_sec_of_day = -1;
  long double last_11_sec_of_day = -1;
  long double last_12_sec_of_day = -1;
  long double last_20_sec_of_day = -1;
  long double last_21_sec_of_day = -1;
  long double last_30_sec_of_day = -1;
  long double last_40_sec_of_day = -1;
  long double test_sec_of_day = -1;
  char ltab[4];
  char crdstr[256];
  char str[256];
  char target_name[11];
  char ttarget_name[11];
  char stn_name[11];
  char infilename[256];
  char nvern[6] = "1.01b";
  char *yn_str[3] = { "no ", "yes", "err" };
  char *rti_str[6] = { "xmit ", "1-way", "2-way", " xcv ", "mixed", "error" };
  char *dqa_str[8] = { "good   ", "suspect", "poor   ", "error  " };
  char *dt_str[4] = { "frd", "npt", "qlk", "err" };
  char *tar_str[5] = { "err", "slr", "llr", "sxp", "axp" };

  if (argc == 3)
    {
      if (strcmp (argv[1], "-b") == 0)
	red_flag_only = 1;
      strcpy (infilename, argv[2]);
    }
  else if (argc == 2)
    {
      strcpy (infilename, argv[1]);
    }
  else
    {
      printf ("Usage: crd_chk [-b] filename\n");
    }

  if ((str_in = fopen (infilename, "r")) == NULL)
    {
      sprintf (str, "crd_chk: FATAL! could not open %s\n", infilename);
      perror (str);
      exit (-1);
    }

  printf ("CRD Check Version: %s\nCRD File name: %s\n\n", nvern, infilename);

  /* initialize all the arrays, in case there is nothing to read */

/** Read section **/
  while (fgets (crdstr, 256, str_in) != NULL)
    {
/*printf("crdstr=[%s]\n",crdstr);*/
      if (strncmp (crdstr, "h1", 2) == 0 || strncmp (crdstr, "H1", 2) == 0)
	{
	  if (strncmp (&crdstr[3], "CRD", 3) != 0
	      && strncmp (&crdstr[3], "crd", 3) != 0)
	    {
	      strncpy (ltab, &crdstr[3], 3);
	      ftab = 1;
	    }
	  getifield (crdstr, 7, 2, &h1.format_version);
	  getifield (crdstr, 10, 4, &h1.prod_year);
	  getifield (crdstr, 15, 2, &h1.prod_mon);
	  getifield (crdstr, 18, 2, &h1.prod_day);
	  getifield (crdstr, 21, 2, &h1.prod_hour);

	  i = 0;
	  if (crdstr[2] != ' ')
	    nh1_bl[i++] = 3;
	  if (crdstr[6] != ' ')
	    nh1_bl[i++] = 7;
	  if (crdstr[9] != ' ')
	    nh1_bl[i++] = 10;
	  if (crdstr[14] != ' ')
	    nh1_bl[i++] = 15;
	  if (crdstr[17] != ' ')
	    nh1_bl[i++] = 18;
	  if (crdstr[20] != ' ')
	    nh1_bl[i++] = 21;
	  nh1_blerr = i;
	  nh1++;
	}
      else
	if (strncmp (crdstr, "h2", 2) == 0 || strncmp (crdstr, "H2", 2) == 0)
	{
	  strncpy (h2.stn_name, &crdstr[3], 10);
	  /*for (i=0;i<10;i++)
	     if ( */
	  getifield (crdstr, 14, 4, &h2.cdp_pad_id);
	  getifield (crdstr, 19, 2, &h2.cdp_sys_num);
	  getifield (crdstr, 22, 2, &h2.cdp_occ_num);
	  getifield (crdstr, 25, 2, &h2.stn_timescale);

	  i = 0;
	  if (crdstr[2] != ' ')
	    nh2_bl[i++] = 3;
	  if (crdstr[13] != ' ')
	    nh2_bl[i++] = 13;
	  if (crdstr[18] != ' ')
	    nh2_bl[i++] = 18;
	  if (crdstr[21] != ' ')
	    nh2_bl[i++] = 22;
	  if (crdstr[24] != ' ')
	    nh2_bl[i++] = 25;
	  nh2_blerr = i;
	  nh2++;
	}
      else
	if (strncmp (crdstr, "h3", 2) == 0 || strncmp (crdstr, "H3", 2) == 0)
	{
	  strncpy (h3.target_name, &crdstr[3], 10);
	  h3.target_name[10] = '\0';
	  /*for (i=0;i<10;i++)
	     if ( */
	  getifield (crdstr, 14, 8, &h3.ilrs_id);
	  getifield (crdstr, 23, 4, &h3.sic);
	  getifield (crdstr, 28, 8, &h3.norad);
	  getifield (crdstr, 37, 1, &h3.SC_timescale);
	  getifield (crdstr, 39, 1, &h3.target_type);
	  if (h3.target_type < 1 || h3.target_type > 4)
	    target_type = 0;
	  else
	    target_type = h3.target_type;

	  /* Assume the ilrs id is correct and check everything else */
	    nh3iderr =
	      get_sat_ids (0, &h3.ilrs_id, &tnorad, &tsic, &ttarget_name,
			   &ttarget_type);
	  if (!nh3iderr)
	    {
	      if (strncasecmp (h3.target_name, ttarget_name, 10) != 0 &&
		  strcasestr (h3.target_name, "na") == 0)
		nh3nameerr = 1;
	      if (tsic != h3.sic)
		nh3sicerr = 1;
	      if (tnorad != h3.norad && h3.norad != -1)
		nh3noraderr = 1;
	      if (ttarget_type != h3.target_type)
		nh3ttypeerr = 1;
	    }

	  i = 0;
	  if (crdstr[2] != ' ')
	    nh3_bl[i++] = 3;
	  if (crdstr[13] != ' ')
	    nh3_bl[i++] = 14;
	  if (crdstr[22] != ' ')
	    nh3_bl[i++] = 23;
	  if (crdstr[27] != ' ')
	    nh3_bl[i++] = 28;
	  if (crdstr[36] != ' ')
	    nh3_bl[i++] = 37;
	  if (crdstr[38] != ' ')
	    nh3_bl[i++] = 39;
	  if (crdstr[38] != ' ')
	    nh3_bl[i++] = 39;
	  nh3_blerr = i;
	  nh3++;
	}
      else
	if (strncmp (crdstr, "h4", 2) == 0 || strncmp (crdstr, "H4", 2) == 0)
	{
	  getifield (crdstr, 3, 2, &h4.data_type);
	  getifield (crdstr, 6, 4, &h4.start_year);
	  getifield (crdstr, 11, 2, &h4.start_mon);
	  getifield (crdstr, 14, 2, &h4.start_day);
	  getifield (crdstr, 17, 2, &h4.start_hour);
	  getifield (crdstr, 20, 2, &h4.start_min);
	  getifield (crdstr, 23, 2, &h4.start_sec);
	  getifield (crdstr, 26, 4, &h4.end_year);
	  getifield (crdstr, 31, 2, &h4.end_mon);
	  getifield (crdstr, 34, 2, &h4.end_day);
	  getifield (crdstr, 37, 2, &h4.end_hour);
	  getifield (crdstr, 40, 2, &h4.end_min);
	  getifield (crdstr, 43, 2, &h4.end_sec);
	  getifield (crdstr, 46, 2, &h4.data_release);
	  getifield (crdstr, 49, 1, &h4.refraction_app_ind);
	  getifield (crdstr, 51, 1, &h4.CofM_app_ind);
	  getifield (crdstr, 53, 1, &h4.xcv_amp_app_ind);
	  getifield (crdstr, 55, 1, &h4.stn_sysdelay_app_ind);
	  getifield (crdstr, 57, 1, &h4.SC_sysdelay_app_ind);
	  getifield (crdstr, 59, 2, &h4.range_type_ind);
	  getifield (crdstr, 61, 1, &h4.data_qual_alert_ind);
	  if (h4.data_type < 0 || h4.data_type > 2)
	    data_type = 3;
	  else
	    data_type = h4.data_type;
	  if (h4.refraction_app_ind < 0 || h4.refraction_app_ind > 1)
	    refraction_app_ind = 2;
	  else
	    refraction_app_ind = h4.refraction_app_ind;
	  if (h4.CofM_app_ind < 0 || h4.CofM_app_ind > 1)
	    CofM_app_ind = 2;
	  else
	    CofM_app_ind = h4.CofM_app_ind;
	  if (h4.xcv_amp_app_ind < 0 || h4.xcv_amp_app_ind > 1)
	    xcv_amp_app_ind = 2;
	  else
	    xcv_amp_app_ind = h4.xcv_amp_app_ind;
	  if (h4.stn_sysdelay_app_ind < 0 || h4.stn_sysdelay_app_ind > 1)
	    stn_sysdelay_app_ind = 2;
	  else
	    stn_sysdelay_app_ind = h4.stn_sysdelay_app_ind;
	  if (h4.SC_sysdelay_app_ind < 0 || h4.SC_sysdelay_app_ind > 1)
	    SC_sysdelay_app_ind = 2;
	  else
	    SC_sysdelay_app_ind = h4.SC_sysdelay_app_ind;
	  if (h4.range_type_ind < 0 || h4.range_type_ind > 4)
	    range_type_ind = 5;
	  else
	    range_type_ind = h4.range_type_ind;
	  if (h4.data_qual_alert_ind < 0 || h4.data_qual_alert_ind > 2)
	    data_qual_alert_ind = 3;
	  else
	    data_qual_alert_ind = h4.data_qual_alert_ind;

	  i = 0;
	  if (crdstr[2] != ' ')
	    nh4_bl[i++] = 3;
	  if (crdstr[5] != ' ')
	    nh4_bl[i++] = 6;
	  if (crdstr[10] != ' ')
	    nh4_bl[i++] = 11;
	  if (crdstr[13] != ' ')
	    nh4_bl[i++] = 14;
	  if (crdstr[16] != ' ')
	    nh4_bl[i++] = 17;
	  if (crdstr[19] != ' ')
	    nh4_bl[i++] = 20;
	  if (crdstr[22] != ' ')
	    nh4_bl[i++] = 23;
	  if (crdstr[25] != ' ')
	    nh4_bl[i++] = 26;
	  if (crdstr[30] != ' ')
	    nh4_bl[i++] = 31;
	  if (crdstr[33] != ' ')
	    nh4_bl[i++] = 34;
	  if (crdstr[36] != ' ')
	    nh4_bl[i++] = 37;
	  if (crdstr[39] != ' ')
	    nh4_bl[i++] = 40;
	  if (crdstr[42] != ' ')
	    nh4_bl[i++] = 43;
	  if (crdstr[45] != ' ')
	    nh4_bl[i++] = 46;
	  if (crdstr[48] != ' ')
	    nh4_bl[i++] = 49;
	  if (crdstr[50] != ' ')
	    nh4_bl[i++] = 51;
	  if (crdstr[52] != ' ')
	    nh4_bl[i++] = 53;
	  if (crdstr[54] != ' ')
	    nh4_bl[i++] = 55;
	  if (crdstr[56] != ' ')
	    nh4_bl[i++] = 57;
	  if (crdstr[58] != ' ')
	    nh4_bl[i++] = 59;
	  if (crdstr[60] != ' ')
	    nh4_bl[i++] = 61;
	  nh4_blerr = i;
	  nh4++;

	  if (!red_flag_only)
	    {
	      strncpy (stn_name, h2.stn_name, 10);
	      strncpy (target_name, h3.target_name, 10);
	      target_name[10] = '\0';
	      stn_name[10] = '\0';
	      printf ("Station: Pad ID:%4d  Sys Num:%2d  Occ num: %1d  Name: [%s] time scale: %d\n\
Target: COSPAR:%8d    SIC:%4d    NORAD: %8d    Name: [%s]\n\
        time scale: %d    target type: %s\n\
Format version: %d    Data type: %s    Data release: %d\n\
Produced on %4d/%2d/%2d %2d hours\n\
Session Starting on %4d/%2d/%2d %2d:%2d:%2d \
Ending on %4d/%2d/%2d %2d:%2d:%2d UTC\n\
refraction applied: %s    CofM applied: %s    xcv amplitude applied: %s\n\
Station system delay applied: %s    Satellite system delay applied: %s\n\
Range type: %s    Date quality alert: %s\n\n",
		      h2.cdp_pad_id, h2.cdp_sys_num, h2.cdp_occ_num, h2.stn_name, h2.stn_timescale, h3.ilrs_id, h3.sic, h3.norad, h3.target_name, h3.SC_timescale, tar_str[target_type], h1.format_version, dt_str[h4.data_type], h4.data_release, h1.prod_year, h1.prod_mon, h1.prod_day, h1.prod_hour, h4.start_year, h4.start_mon, h4.start_day, h4.start_hour, h4.start_min, h4.start_sec, h4.end_year, h4.end_mon, h4.end_day, h4.end_hour, h4.end_min, h4.end_sec, yn_str[refraction_app_ind], yn_str[CofM_app_ind], yn_str[xcv_amp_app_ind], yn_str[stn_sysdelay_app_ind], yn_str[SC_sysdelay_app_ind], rti_str[range_type_ind],
		      dqa_str[data_qual_alert_ind]);
	    }
	  else
	    {
	      printf
		("Station: %4d ILRS ID: %8d Start Date: %4d/%2d/%2d %2d:%2d:%2d\n",
		 h2.cdp_pad_id, h3.ilrs_id, h4.start_year, h4.start_mon,
		 h4.start_day, h4.start_hour, h4.start_min, h4.start_sec);
	    }
	  for (i = 0; i < nh1_blerr; i++)
	    printf ("Error: H1 column %d is not blank but should be\n",
		    nh1_bl[i]);
	  for (i = 0; i < nh2_blerr; i++)
	    printf ("Error: H2 column %d is not blank but should be\n",
		    nh2_bl[i]);
	  for (i = 0; i < nh3_blerr; i++)
	    printf ("Error: H3 column %d is not blank but should be\n",
		    nh3_bl[i]);
	  for (i = 0; i < nh4_blerr; i++)
	    printf ("Error: H4 column %d is not blank but should be\n",
		    nh4_bl[i]);
	  for (i = 0; i < nh8_blerr; i++)
	    printf ("Error: H5 column %d is not blank but should be\n",
		    nh8_bl[i]);
	  nh1_blerr = nh2_blerr = nh3_blerr = nh4_blerr = nh8_blerr = 0;

	  if (ftab)
	    printf
	      ("Error: Header record 1 missing 'CRD' literal. Is [%s] instead.\n",
	       ltab);
	  ftab = 0;
	  if (nh3iderr == -1)
	    printf ("Error: Could not open targets.dat file\n");
	  if (nh3iderr == 1)
	    printf ("Error: ILRS ID %7d not found in targets.dat file\n",
		    h3.ilrs_id);
	  nh3iderr = 0;
	  if (nh3nameerr)
	    printf ("Error: Target name (%s) does not match official target name (%s)\n\
 based on ILRS ID (%07d)\n",
		    h3.target_name, ttarget_name, h3.ilrs_id);
	  nh3nameerr = 0;
	  if (nh3noraderr)
	    printf ("Error: Target NORAD ID (%05d) does not match official target NORAD ID (%05d)\n\
 based on ILRS ID (%07d)\n",
		    h3.norad, tnorad, h3.ilrs_id);
	  nh3noraderr = 0;
	  if (nh3sicerr)
	    printf ("Error: Target SIC (%04d) does not match official target SIC (%04d)\n\
 based on ILRS ID (%07d)\n",
		    h3.sic, tsic, h3.ilrs_id);
	  nh3sicerr = 0;
	  if (nh3ttypeerr)
	    printf ("Error: Target type (%d) does not match official target type (%d)\n\
 based on ILRS ID (%7d)\n",
		    h3.target_type, ttarget_type, h3.ilrs_id);
	  nh3ttypeerr = 0;
	  if (h3.target_type < 1 || h3.target_type > 4)
	    printf
	      ("Error: Target type (value = %d) beyond accepted values of 1 and 4\n",
	       h3.target_type);
/**
	  if ((h3.target_type == 3 || h3.target_type == 4) && nc4 == 0)
	    printf
	      ("Error: Tranponders require transponder configuration record C4\n");
	  nc4 = 0;
**/

	}
      else
	if (strncmp (crdstr, "h8", 2) == 0 || strncmp (crdstr, "H8", 2) == 0)
	{
	  nh8++;

	  if (nseq > 0)
	    {
	      printf ("Records not in sequence in %d cases\n", nseq);
	      tnseq++;
	    }
	  nseq = 0;
          /* 12/18/12 */
	  if ((h3.target_type == 3 || h3.target_type == 4) && nc4 == 0)
	    printf
	      ("Error: Transponders require transponder configuration record C4\n");
	  nc4 = 0;
	  last_10_sec_of_day = -1;
	  last_11_sec_of_day = -1;
	  last_12_sec_of_day = -1;
	  last_20_sec_of_day = -1;
	  last_21_sec_of_day = -1;
	  last_30_sec_of_day = -1;
	  last_40_sec_of_day = -1;
	  test_sec_of_day = -1;
	  h3.target_name[0] = '\0';
	  h3.ilrs_id = 0;
	  h3.norad = 0;
	  h3.sic = 0;
	  h3.target_type = 0;
	  h2.cdp_pad_id = 0;
	  h2.cdp_sys_num = 0;
	  h2.cdp_occ_num = 0;
	  h2.stn_name[0] = '\0';
	}
      else
	if (strncmp (crdstr, "h9", 2) == 0 || strncmp (crdstr, "H9", 2) == 0)
	{
	  nh9++;
	}
      else
	if (strncmp (crdstr, "c0", 2) == 0 || strncmp (crdstr, "C0", 2) == 0)
	{
	  nstat = read_c0 (crdstr, &c0);
	  nc0++;
	  if (nstat == 0)
	    nc0bad++;
	  /* May want to test for >= 4 and <= ... */
	  if (nstat > 0 && nstat < 4)
	    nc0low++;
	  if (c0.detail_type != 0)
	    nc0dterr++;
	  if ((fabs (c0.xmit_wavelength - 1064) > 1 &&
	       fabs (c0.xmit_wavelength - 532) > 1 &&
	       fabs (c0.xmit_wavelength - 266) > 1) || c0.xmit_wavelength < 0)
	    nc0xwerr++;
	}
      else
	if (strncmp (crdstr, "c1", 2) == 0 || strncmp (crdstr, "C1", 2) == 0)
	{
	  nstat = read_c1 (crdstr, &c1);
	  nc1++;
	  if (nstat == 0)
	    nc1bad++;
	  if (nstat > 0 && nstat < 10)
	    nc1low++;
	  if (c1.detail_type != 0)
	    nc1dterr++;
	  if (c1.prim_wavelength < -1)
	    nc1pwerr++;
	  if (c1.nom_fire_rate < -1)
	    nc1nfrerr++;
	  if (c1.pulse_energy < -1)
	    nc1peerr++;
	  if (c1.pulse_width < -1)
	    nc1pwderr++;
	  if (c1.beam_div < -1)
	    nc1bderr++;
	  if (c1.pulses_in_semitrain < -1)
	    nc1piserr++;
//printf("c1: %f %f %f %f %f\n", c1.prim_wavelength,c1.nom_fire_rate,c1.pulse_energy,c1.pulse_width,c1.beam_div);
	}
      else
	if (strncmp (crdstr, "c2", 2) == 0 || strncmp (crdstr, "C2", 2) == 0)
	{
	  nstat = read_c2 (crdstr, &c2);
	  nc2++;
	  if (nstat == 0)
	    nc2bad++;
	  if (nstat > 0 && nstat < 14)
	    nc2low++;
	  if (c2.detail_type != 0)
	    nc2dterr++;
	  if (c2.app_wavelength < -1)
	    nc2awerr++;
	  if (c2.qe < -1 || c2.qe > 100)
	    nc2qeerr++;
	  if (c2.voltage < -1.e6 || c2.voltage > 1e6)
	    nc2verr++;
	  if (c2.dark_count < -1 || c2.dark_count > 1e6)
	    nc2dcerr++;
	  if (c2.output_pulse_width < -1 || c2.output_pulse_width > 1e6)
	    nc2opwerr++;
	  if (c2.spectral_filter < -1 || c2.spectral_filter > 1e6)
	    nc2sferr++;
	  if (c2.spectral_filter_xmission < -1
	      || c2.spectral_filter_xmission > 100)
	    nc2sfxerr++;
	  if (c2.spatial_filter < -1 || c2.spatial_filter > 1e6)
	    nc2stferr++;
//printf("c2: %f %f %f %f %f %f %f %f\n",c2.app_wavelength,c2.qe,c2.voltage,c2.dark_count,c2.output_pulse_width,c2.spectral_filter,c2.spectral_filter_xmission,c2.spatial_filter);
	}
      else
	if (strncmp (crdstr, "c3", 2) == 0 || strncmp (crdstr, "C3", 2) == 0)
	{
	  nstat = read_c3 (crdstr, &c3);
	  nc3++;
	  if (nstat == 0)
	    nc3bad++;
	  if (nstat > 0 && nstat < 8)
	    nc3low++;
	  if (c3.detail_type != 0)
	    nc3dterr++;
	  //if (c3.epoch_delay_corr < -1)
	  if (c3.epoch_delay_corr < -1000000)
	    nc3edcerr++;
//printf("c3: %f\n", c3.epoch_delay_corr);
	}
      else
	if (strncmp (crdstr, "c4", 2) == 0 || strncmp (crdstr, "C4", 2) == 0)
	{
	  nstat = read_c4 (crdstr, &c4);
	  nc4++;
	  if (nstat == 0)
	    nc4bad++;
	  if (nstat > 0 && nstat < 11)
	    nc4low++;
	  if (c4.detail_type != 0)
	    nc4dterr++;
	  if (c4.stn_off_drift_app_ind < 0 || c4.stn_off_drift_app_ind > 3)
	    nc4sodaierr++;
	  if (c4.SC_off_drift_app_ind < 0 || c4.SC_off_drift_app_ind > 3)
	    nc4scodaierr++;
	  if (c4.SC_time_simplified_ind < 0 || c4.SC_time_simplified_ind > 1)
	    nc4sctsierr++;
	}
      else
	/* No choice but to use sscanf on free format records */
      if (strncmp (crdstr, "10", 2) == 0)
	{
	  nstat = read_10 (crdstr, &d10);
	  n10++;
	  if (nstat == 0)
	    n10bad++;
	  if (nstat > 0 && nstat < 9)
	    n10low++;
	  if (fex_seq)
	    test_sec_of_day = d10.sec_of_day;
	  else
	    test_sec_of_day =
	      (long) ((d10.sec_of_day + 0.0005) * 1000) / 1000.;
	  /* Want data in sequence, but wrap around from 86400 to 0 is OK */
	  if (test_sec_of_day < last_10_sec_of_day &&
	      (last_10_sec_of_day - test_sec_of_day) < 80000)
	    {
	      nseq++;
	      if (fout_of_seq)
		printf ("Error: out of seq: [%s]\n", crdstr);
	    }
	  last_10_sec_of_day = test_sec_of_day;
	  if (d10.sec_of_day < 0 || d10.sec_of_day > 86400)
	    nd10soderr++;
	  if (d10.time_of_flight < -1 || d10.time_of_flight > 10000)
	    nd10toferr++;
	  if (d10.epoch_event < 0 || d10.epoch_event > 6)
	    nd10eeerr++;
	  if (d10.filter_flag < 0 || d10.filter_flag > 2)
	    nd10fferr++;
	  if (d10.detector_channel < 0 /*|| d10.detector_channel > 100 */ )
	    nd10dcerr++;
	  if (d10.stop_number < 0 /*|| d10.stop_number > 100 */ )
	    nd10snerr++;
	  if (d10.xcv_amp < -1 /*|| d10.stop_number > 100 */ )
	    nd10xaerr++;
	}
      else if (strncmp (crdstr, "11", 2) == 0)
	{
	  nstat = read_11 (crdstr, &d11);
/**
          printf (
           "11 %.7Lf %.12Lf %-s %1d %4d %6d %.1f %.3f %.3f %.1f %.1f\n",
           d11.sec_of_day, d11.time_of_flight,
           d11.sysconfig_id, d11.epoch_event,
           d11.np_window_length, d11.num_ranges,
           d11.bin_rms, d11.bin_skew, d11.bin_kurtosis,
           d11.bin_PmM, d11.return_rate);
**/

	  n11++;
	  if (nstat == 0)
	    n11bad++;
	  if (nstat > 0 && nstat < 13)
	    n11low++;
	  if (fex_seq)
	    test_sec_of_day = d11.sec_of_day;
	  else
	    test_sec_of_day =
	      (long) ((d11.sec_of_day + 0.0005) * 1000) / 1000.;
	  if (test_sec_of_day < last_11_sec_of_day &&
	      (last_11_sec_of_day - test_sec_of_day) < 80000)
	    {
	      nseq++;
	      if (fout_of_seq)
		printf ("Error: out of seq: [%s]\n", crdstr);
	    }
	  last_11_sec_of_day = test_sec_of_day;
	  if (d11.sec_of_day < 0 || d11.sec_of_day > 86400)
	    nd11soderr++;
	  if (d11.epoch_event < 5 && (d11.time_of_flight < -1 ||
				      d11.time_of_flight > 10000))
	    nd11toferr++;
	  if (d11.epoch_event < 0 || d11.epoch_event > 6)
	    nd11eeerr++;
	  if (d11.np_window_length < 0 /*|| d11.np_window_length > 10000 */ )
	    nd11nwlerr++;
	  if (d11.num_ranges < 0 /*|| d11.num_ranges > 1.e6 */ )
	    nd11nrerr++;
	  if (d11.bin_rms < 0 /*|| d11.bin_rms > 1.e6 */ )
	    nd11brerr++;
	  if (d11.return_rate < -1 || d11.return_rate > 100)
	    nd11rrerr++;
	}
      else if (strncmp (crdstr, "12", 2) == 0)
	{
	  nstat = read_12 (crdstr, &d12);
	  n12++;
	  if (nstat == 0)
	    n12bad++;
	  if (nstat > 0 && nstat < 7)
	    n12low++;
	  if (fex_seq)
	    test_sec_of_day = d12.sec_of_day;
	  else
	    test_sec_of_day =
	      (long) ((d12.sec_of_day + 0.0005) * 1000) / 1000.;
	  if (test_sec_of_day < last_12_sec_of_day &&
	      (last_12_sec_of_day - test_sec_of_day) < 80000)
	    {
	      nseq++;
	      if (fout_of_seq)
		printf ("Error: out of seq: [%s]\n", crdstr);
	    }
	  last_12_sec_of_day = test_sec_of_day;
	  if (d12.sec_of_day < 0 || d12.sec_of_day > 86400)
	    nd12soderr++;
	}
      else if (strncmp (crdstr, "20", 2) == 0)
	{
	  nstat = read_20 (crdstr, &d20);
	  n20++;
	  if (nstat == 0)
	    n20bad++;
	  if (nstat > 0 && nstat < 6)
	    n20low++;
	  if (fex_seq)
	    test_sec_of_day = d20.sec_of_day;
	  else
	    test_sec_of_day =
	      (long) ((d20.sec_of_day + 0.0005) * 1000) / 1000.;
	  if (test_sec_of_day < last_20_sec_of_day &&
	      (last_20_sec_of_day - test_sec_of_day) < 80000)
	    {
	      nseq++;
	      if (fout_of_seq)
		printf ("Error: out of seq: [%s]\n", crdstr);
	    }
	  last_20_sec_of_day = test_sec_of_day;
	  if (d20.sec_of_day < 0 || d20.sec_of_day > 86400)
	    nd20soderr++;
	  if (d20.pressure < 500 || d20.pressure > 1200)
	    nd20preserr++;
	  if (d20.temperature < 200 || d20.temperature > 350)
	    nd20temperr++;
	  if (d20.humidity < 0 || d20.humidity > 100)
	    nd20humerr++;
	  if (d20.value_origin < -1 || d20.value_origin > 1)
	    nd20voerr++;
	}
      else if (strncmp (crdstr, "21", 2) == 0)
	{
	  nstat = read_21 (crdstr, &d21);
	  n21++;
	  if (nstat == 0)
	    n21bad++;
	  if (nstat > 0 && nstat < 9)
	    n21low++;
	  if (fex_seq)
	    test_sec_of_day = d21.sec_of_day;
	  else
	    test_sec_of_day =
	      (long) ((d21.sec_of_day + 0.0005) * 1000) / 1000.;
	  if (test_sec_of_day < last_21_sec_of_day &&
	      (last_21_sec_of_day - test_sec_of_day) < 80000)
	    {
	      nseq++;
	      if (fout_of_seq)
		printf ("Error: out of seq: [%s]\n", crdstr);
	    }
	  last_21_sec_of_day = test_sec_of_day;
	  if (d21.sec_of_day < 0 || d21.sec_of_day > 86400)
	    nd21soderr++;
	  if (d21.wind_speed < -1 || d21.wind_speed > 100)
	    nd21wserr++;
	  if (d21.wind_direction < -1 || d21.wind_direction > 360)
	    nd21wderr++;
	  if (d21.visibility < -1 || d21.visibility > 100)
	    nd21viserr++;
	  if (d21.sky_clarity < -1 || d21.sky_clarity > 100)
	    nd21scerr++;
	  if (d21.atmospheric_seeing < -1 || d21.atmospheric_seeing > 100)
	    nd21aserr++;
	  if (d21.cloud_cover < -1 || d21.cloud_cover > 100)
	    nd21ccerr++;
	}
      else if (strncmp (crdstr, "30", 2) == 0)
	{
	  nstat = read_30 (crdstr, &d30);
	  n30++;
	  if (nstat == 0)
	    n30bad++;
	  if (nstat > 0 && nstat < 7)
	    n30low++;
	  if (fex_seq)
	    test_sec_of_day = d30.sec_of_day;
	  else
	    test_sec_of_day =
	      (long) ((d30.sec_of_day + 0.0005) * 1000) / 1000.;
	  if (test_sec_of_day < last_30_sec_of_day &&
	      (last_30_sec_of_day - test_sec_of_day) < 80000)
	    {
	      nseq++;
	      if (fout_of_seq)
		printf ("Error: out of seq: [%s]\n", crdstr);
	    }
	  last_30_sec_of_day = test_sec_of_day;
	  if (d30.sec_of_day < 0 || d30.sec_of_day > 86400)
	    nd30soderr++;
	  if (d30.azimuth < -360 || d30.azimuth > 360)
	    nd30azerr++;	/* Allow bakwards wrap */
	  if (d30.elevation < -1 || d30.elevation > 180)
	    nd30elerr++;	/* allow dump mode? */
	  if (d30.direction_ind < -1 || d30.direction_ind > 2)
	    nd30dierr++;
	  if (d30.angle_origin_ind < 0 || d30.angle_origin_ind > 3)
	    nd30aoierr++;
	  if (d30.refraction_corr_ind < 0 || d30.refraction_corr_ind > 1)
	    nd30rcierr++;
	}
      else if (strncmp (crdstr, "40", 2) == 0)
	{
	  nstat = read_40 (crdstr, &d40);
	  n40++;
	  if (nstat == 0)
	    n40bad++;
	  if (nstat > 0 && nstat < 16)
	    n40low++;
	  if (fex_seq)
	    test_sec_of_day = d40.sec_of_day;
	  else
	    test_sec_of_day =
	      (long) ((d40.sec_of_day + 0.0005) * 1000) / 1000.;
	  if (test_sec_of_day < last_40_sec_of_day &&
	      (last_40_sec_of_day - test_sec_of_day) < 80000)
	    {
	      nseq++;
	      if (fout_of_seq)
		printf ("Error: out of seq: [%s]\n", crdstr);
	    }
	  last_40_sec_of_day = test_sec_of_day;
	  if (d40.sec_of_day < 0 || d40.sec_of_day > 86400)
	    nd40soderr++;
	  if (d40.type_of_data < 0 || d40.type_of_data > 5)
	    nd40toderr++;
	  if (d40.num_points_recorded < -1 || d40.num_points_recorded > 1.e8)
	    nd40nprerr++;
	  if (d40.num_points_used < -1 || d40.num_points_used > 1.e8)
	    nd40npuerr++;
	  if (d40.one_way_target_dist < -1.e4
	      || d40.one_way_target_dist > 1.e8)
	    nd40owtderr++;
	  if (d40.cal_sys_delay < -1.e4 || d40.cal_sys_delay > 1.e8)
	    nd40csderr++;
	  if (d40.cal_delay_shift < -1.e4 || d40.cal_delay_shift > 1.e8)
	    nd40cdserr++;
	  if (d40.cal_rms < -1 || d40.cal_rms > 1.e8)
	    nd40crerr++;
	  if (d40.cal_type_ind < 0 || d40.cal_type_ind > 5)
	    nd40ctierr++;
	  if (d40.cal_shift_type_ind < 0 || d40.cal_shift_type_ind > 4)
	    nd40cstierr++;
	}
      else if (strncmp (crdstr, "50", 2) == 0)
	{
	  nstat = read_50 (crdstr, &d50);
	  n50++;
	  if (nstat == 0)
	    n50bad++;
	  if (nstat > 0 && nstat < 7)
	    n50low++;
          // Changed to permit llr fr file to have 50 record. 
          //if (d50.sess_rms < 0 || d50.sess_rms > 1.e8)
	  if (((d50.sess_rms < 0 || d50.sess_rms > 1.e8) && h4.data_type != 0)||
	      (d50.sess_rms < -1 || d50.sess_rms > 1.e8) && h4.data_type == 0)
	    nd50srerr++;
	  if (d50.data_qual_ind < 0 || d50.data_qual_ind > 5)
	    nd50dqierr++;
	}
      else if (strncmp (crdstr, "60", 2) == 0)
	{
	  nstat = read_60 (crdstr, &d60);
	  n60++;
	  if (nstat == 0)
	    n60bad++;
	  if (nstat > 0 && nstat < 3)
	    n60low++;
	  if (d60.sys_change_ind < -1 || d60.sys_change_ind > 9)
	    nd60scherr++;
	  if (d60.sys_config_ind < -1 || d60.sys_config_ind > 9)
	    nd60scierr++;
	}
      else if (strncmp (crdstr, "9", 1) == 0)
	{
	  n9x++;
	}
      /* Comment */
      else if (strncmp (crdstr, "00", 2) == 0)
	{
	  n00++;
	}
    }

/* Recap */
  if (!red_flag_only)
    {
      printf ("\nNumber of records -\n");
      if (n10 > 0)
	printf (" 10 Range Record: %d\n", n10);
      if (n11 > 0)
	printf (" 11 Normal Point Record: %d\n", n11);
      if (n12 > 0)
	printf (" 12 Range Supplement Record: %d\n", n12);
      if (n20 > 0)
	printf (" 20 Meteorological Record: %d\n", n20);
      if (n21 > 0)
	printf (" 21 Met Supplement Record: %d\n", n21);
      if (n30 > 0)
	printf (" 30 Pointing Angles Record: %d\n", n30);
      if (n40 > 0)
	printf (" 40 Calibration Record: %d\n", n40);
      if (n50 > 0)
	printf (" 50 Session Statistics Record: %d\n", n50);
      if (n60 > 0)
	printf (" 60 Compatibility Record: %d\n", n60);
      if (n9x > 0)
	printf (" 9x User defined Records): %d\n", n9x);
      if (n00 > 0)
	printf (" 00 Comment Records: %d\n", n00);
      if (nc0 > 0)
	printf (" c0 System Configuration Records: %d\n", nc0);
      if (nc1 > 0)
	printf (" c1 Laser Configuration Records: %d\n", nc1);
      if (nc2 > 0)
	printf (" c2 Detector Configuration Records: %d\n", nc2);
      if (nc3 > 0)
	printf (" c3 Timing Configuration Records: %d\n", nc3);
      if (nc4 > 0)
	printf (" c4 Transponder Configuration Records: %d\n", nc4);

      printf ("\nErrors and Warnings -\n");
    }

/*  if (red_flag_only)
    printf ("CRD File name: %s \n\n", infilename);*/
  if (nh1 == 0)
    printf ("ERROR: No header H1\n");
  if (nh2 == 0)
    printf ("ERROR: No header H2\n");
  if (nh3 == 0)
    printf ("ERROR: No header H3\n");
  if (nh4 == 0)
    printf ("ERROR: No header H4\n");
  if (nc0 == 0)
    printf ("ERROR: No configuration record C0\n");
  if (n60 == 0 && nc1 == 0 && nc2 == 0 && nc3 == 0)
    printf ("ERROR: Must be configuration records C1-3 or 60 record\n");

  if (nh1 == 0 || nh2 == 0 || nh3 == 0 || nh4 == 0 || nc0 == 0 ||
     (n60 == 0 && nc1 == 0 && nc2 == 0 && nc3 == 0))
    fail = 1;

  if (!red_flag_only)
    {
      /*if (nc4 == 0) printf ("NOTE: No transponder configuration record C4\n"); */
    }
  if (nh8 == 0)
    printf ("ERROR: No end of session record H8\n");
  if (nh9 == 0)
    printf ("ERROR: No end of file record H9\n");
  if (nh1 > 1)
    printf ("Multi-pass file with: %d H1 Records\n", nh1);
  if (nh2 > 1)
    printf ("Multi-pass file with: %d H2 Records\n", nh2);
  if (nh3 > 1)
    printf ("Multi-pass file with: %d H3 Records\n", nh3);
  if (nh4 > 1)
    printf ("Multi-pass file with: %d H4 Records\n", nh4);
  if (nh8 > 1)
    printf ("Multi-pass file with: %d H8 Records\n", nh8);
  if (nh9 > 1)
    printf ("Multi-pass file with: %d H9 Records\n", nh9);
  if (nh1 > nh8 || nh2 > nh8)
    {
      printf ("Too many H1/H2 (%d/%d) or too few H8 Records (%d)\n", nh1, nh2, nh8);
      fail = 1;
    }
/**
  for (i = 0; i < nh1_blerr; i++)
    printf ("H1 column %d is not blank but should be\n", nh1_bl[i]);
  for (i = 0; i < nh2_blerr; i++)
    printf ("H2 column %d is not blank but should be\n", nh2_bl[i]);
  for (i = 0; i < nh3_blerr; i++)
    printf ("H3 column %d is not blank but should be\n", nh3_bl[i]);
  for (i = 0; i < nh4_blerr; i++)
    printf ("H4 column %d is not blank but should be\n", nh4_bl[i]);
  for (i = 0; i < nh8_blerr; i++)
    printf ("H5 column %d is not blank but should be\n", nh8_bl[i]);

  if (nh1_blerr > 0 || nh2_blerr > 0 || nh3_blerr > 0 || nh4_blerr > 0 ||
      nh8_blerr > 0)
    fail = 1;

  if (ftab)
    printf ("Header record 1 missing 'CRD' literal. Is [%s] instead.\n",
	    ltab);
  if (nh3iderr == -1)
    printf ("Could not open targets.dat file\n");
  if (nh3iderr == 1)
    printf ("ILRS ID %7d not found in targets.dat file\n",h3.ilrs_id);
  if (nh3nameerr)
    printf ("Target name (%s) does not match official target name (%s)\n\
 based on ILRS ID (%07d)\n",
            h3.target_name, ttarget_name, h3.ilrs_id);
  if (nh3noraderr)
    printf ("Target NORAD ID (%05d) does not match official target NORAD ID (%05d)\n\
 based on ILRS ID (%07d)\n",
            h3.norad, tnorad, h3.ilrs_id);
  if (nh3sicerr)
    printf ("Target SIC (%04d) does not match official target SIC (%04d)\n\
 based on ILRS ID (%07d)\n",
            h3.sic, tsic, h3.ilrs_id);
  if (nh3ttypeerr)
    printf ("Target type (%d) does not match official target type (%d)\n\
 based on ILRS ID (%7d)\n",
            h3.target_type, ttarget_type, h3.ilrs_id);
  if (h3.target_type < 1 || h3.target_type > 4)
    printf ("Target type (value = %d) beyond accepted values of 1 and 4\n",
	    h3.target_type);
  if ((h3.target_type == 3 || h3.target_type == 4) && nc4 == 0)
    printf ("Tranponders require transponder configuration record C4\n");
**/
  if (n10bad > 0)
    printf ("Bad input on %d records type 10\n", n10bad);
  if (n11bad > 0)
    printf ("Bad input on %d records type 11\n", n11bad);
  if (n12bad > 0)
    printf ("Bad input on %d records type 12\n", n12bad);
  if (n20bad > 0)
    printf ("Bad input on %d records type 20\n", n20bad);
  if (n21bad > 0)
    printf ("Bad input on %d records type 21\n", n21bad);
  if (n22bad > 0)
    printf ("Bad input on %d records type 22\n", n22bad);
  if (n30bad > 0)
    printf ("Bad input on %d records type 30\n", n30bad);
  if (n40bad > 0)
    printf ("Bad input on %d records type 40 ()\n", n40bad);
  if (n50bad > 0)
    printf ("Bad input on %d records type 50 ()\n", n50bad);
  if (n60bad > 0)
    printf ("Bad input on %d records type 60 ()\n", n60bad);
  if (nc0bad > 0)
    printf ("Bad input on %d records type c0 ()\n", nc0bad);
  if (nc1bad > 0)
    printf ("Bad input on %d records type c1 ()\n", nc1bad);
  if (nc2bad > 0)
    printf ("Bad input on %d records type c2 ()\n", nc2bad);
  if (nc3bad > 0)
    printf ("Bad input on %d records type c3 ()\n", nc3bad);
  if (nc4bad > 0)
    printf ("Bad input on %d records type c4 ()\n", nc4bad);
  if (n10low > 0)
    printf ("Too few fields on %d records type 10 ()\n", n10low);
  if (n11low > 0)
    printf ("Too few fields on %d records type 11 ()\n", n11low);
  if (n12low > 0)
    printf ("Too few fields on %d records type 12 ()\n", n12low);
  if (n20low > 0)
    printf ("Too few fields on %d records type 20 ()\n", n20low);
  if (n21low > 0)
    printf ("Too few fields on %d records type 21 ()\n", n21low);
  if (n30low > 0)
    printf ("Too few fields on %d records type 30 ()\n", n30low);
  if (n40low > 0)
    printf ("Too few fields on %d records type 40 ()\n", n40low);
  if (n50low > 0)
    printf ("Too few fields on %d records type 50 ()\n", n50low);
  if (n60low > 0)
    printf ("Too few fields on %d records type 60 ()\n", n60low);
  if (nc0low > 0)
    printf ("Too few fields on %d records type c0 ()\n", nc0low);
  if (nc1low > 0)
    printf ("Too few fields on %d records type c1 ()\n", nc1low);
  if (nc2low > 0)
    printf ("Too few fields on %d records type c2 ()\n", nc2low);
  if (nc3low > 0)
    printf ("Too few fields on %d records type c3 ()\n", nc3low);
  if (nc4low > 0)
    printf ("Too few fields on %d records type c4 ()\n", nc4low);

  //if (ftab || h3.target_type < 1 || h3.target_type > 4 ||
  //  ((h3.target_type == 3 || h3.target_type == 4) && nc4 == 0) ||
  if (ftab || target_type < 1 || target_type > 4 ||
    ((target_type == 3 || target_type == 4) && nc4 == 0) ||
    n10bad > 0 || n11bad > 0 || n12bad > 0 || n20bad > 0 || n21bad > 0 ||
    n22bad > 0 || n30bad > 0 || n40bad > 0 || n50bad > 0 || n60bad > 0 ||
    nc0bad > 0 || nc1bad > 0 || nc2bad > 0 || nc3bad > 0 || nc4bad > 0 ||
    n10low > 0 || n11low > 0 || n12low > 0 || n20low > 0 || n21low > 0 ||
    n30low > 0 || n50low > 0 || n60low > 0 || nc0low > 0 || nc1low > 0 ||
    nc2low > 0 || nc3low > 0 || nc4low > 0)
    fail = 1;
//printf("%d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d\n",
//ftab,target_type,nc4,n10bad,n11bad,n12bad,n20bad,n21bad,n22bad,n30bad,n40bad,n50bad,n60bad,nc0bad,nc1bad,nc2bad,nc3bad,nc4bad,n10low,n11low,n12low,n20low,n21low,n30low,n50low,n60low,nc0low,nc1low,nc2low,nc3low,nc4low);

  if (nc0dterr > 0)
    printf ("Record c0 detail type in error %d times\n", nc0dterr);
  if (nc0xwerr > 0)
    printf ("Record c0 transmit wavelength was non-standard %d times\n",
	    nc0xwerr);

  if (nc0dterr > 0 || nc0xwerr > 0) fail = 1;

  if (nc1dterr > 0)
    printf ("Record c1 detail type in error %d times\n", nc1dterr);
  if (nc1pwerr > 0)
    printf ("Record c1 primary wavelength in error %d times\n", nc1pwerr);
  if (nc1nfrerr > 0)
    printf ("Record c1 nominal fire rate in error %d times\n", nc1nfrerr);
  if (nc1peerr > 0)
    printf ("Record c1 pulse energy in error %d times\n", nc1peerr);
  if (nc1pwderr > 0)
    printf ("Record c1 pulse width in error %d times\n", nc1pwderr);
  if (nc1bderr > 0)
    printf ("Record c1 beam divergence in error %d times\n", nc1bderr);
  if (nc1piserr > 0)
    printf ("Record c1 pulses in semi-train in error %d times\n", nc1piserr);

  if (nc1dterr > 0 || nc1pwerr > 0 || nc1nfrerr > 0 || nc1peerr > 0 ||
    nc1pwderr > 0 || nc1bderr > 0 || nc1piserr > 0)
    fail = 1;

  if (nc2dterr > 0)
    printf ("Record c2 detail type in error %d times\n", nc2dterr);
  if (nc2awerr > 0)
    printf ("Record c2 applicable wavelength in error %d times\n", nc2awerr);
  if (nc2qeerr > 0)
    printf ("Record c2 quantum efficiency in error %d times\n", nc2qeerr);
  if (nc2verr > 0)
    printf ("Record c2 voltage in error %d times\n", nc2verr);
  if (nc2dcerr > 0)
    printf ("Record c2 dark count in error %d times\n", nc2dcerr);
  if (nc2opwerr > 0)
    printf ("Record c2 output pulse width in error %d times\n", nc2opwerr);
  if (nc2sferr > 0)
    printf ("Record c2 spectral filter in error %d times\n", nc2sferr);
  if (nc2sfxerr > 0)
    printf ("Record c2 spectral filter transmission in error %d times\n",
	    nc2sfxerr);
  if (nc2stferr > 0)
    printf ("Record c2 spatial filter in error %d times\n", nc2stferr);

  if (nc2dterr > 0 || nc2awerr > 0 || nc2qeerr > 0 || nc2verr > 0 ||
     nc2dcerr > 0 || nc2opwerr > 0 || nc2sferr > 0 || nc2sfxerr > 0 ||
     nc2stferr > 0)
     fail = 1;

  if (nc3dterr > 0)
    printf ("Record c3 detail type in error %d times\n", nc3dterr);
  if (nc3edcerr > 0)
    printf ("Record c3 epoch delay correction in error %d times\n",
	    nc3edcerr);

  if (nc3dterr > 0 || nc3edcerr > 0) fail = 1;

  if (nc4dterr > 0)
    printf ("Record c4 detail type in error %d times\n", nc4dterr);
  if (nc4sodaierr > 0)
    printf
      ("Record c4 station offset drift appied indicator in error %d times\n",
       nc4sodaierr);
  if (nc4scodaierr > 0)
    printf
      ("Record c4 spacecraft offset drift appied indicator in error %d times\n",
       nc4scodaierr);
  if (nc4sctsierr > 0)
    printf
      ("Record c4 spacecraft time simplified indicator in error %d times\n",
       nc4sctsierr);

  if (nc4dterr > 0 || nc4sodaierr > 0 || nc4scodaierr > 0 || nc4sctsierr > 0)
    fail = 1;

  if (nd10soderr > 0)
    printf ("Record 10 second of day in error %d times\n", nd10soderr);
  if (nd10toferr > 0)
    printf ("Record 10 time of flight in error %d times\n", nd10toferr);
  if (nd10eeerr > 0)
    printf ("Record 10 epoch event in error %d times\n", nd10eeerr);
  if (nd10fferr > 0)
    printf ("Record 10 filter flag in error %d times\n", nd10fferr);
  if (nd10dcerr > 0)
    printf ("Record 10 detector channel in error %d times\n", nd10dcerr);
  if (nd10snerr > 0)
    printf ("Record 10 stop number in error %d times\n", nd10snerr);
  if (nd10xaerr > 0)
    printf ("Record 10 xcv amplitude in error %d times\n", nd10xaerr);

  if (nd10soderr > 0 || nd10toferr > 0 || nd10eeerr > 0 || nd10fferr > 0 ||
    nd10dcerr > 0 || nd10snerr > 0 || nd10xaerr > 0)
    fail = 1;

  if (nd11soderr > 0)
    printf ("Record 11 second of day in error %d times\n", nd11soderr);
  if (nd11toferr > 0)
    printf ("Record 11 time of flight in error %d times\n", nd11toferr);
  if (nd11eeerr > 0)
    printf ("Record 11 epoch event in error %d times\n", nd11eeerr);
  if (nd11nwlerr > 0)
    printf ("Record 11 normalpoint window length in error %d times\n",
	    nd11nwlerr);
  if (nd11nrerr > 0)
    printf ("Record 11 number of ranges in error %d times\n", nd11nrerr);
  if (nd11brerr > 0)
    printf ("Record 11 bin rms in error %d times\n", nd11brerr);
  if (nd11rrerr > 0)
    printf ("Record 11 return rate in error %d times\n", nd11rrerr);

  if (nd11soderr > 0 || nd11toferr > 0 || nd11eeerr > 0 || nd11nwlerr > 0 ||
    nd11nrerr > 0 || nd11brerr > 0 || nd11rrerr > 0)
    fail = 1;

  if (nd12soderr > 0)
    printf ("Record 12 second of day in error %d times\n", nd12soderr);

  if (nd12soderr > 0) fail = 1;

  if (nd20soderr > 0)
    printf ("Record 20 second of day in error %d times\n", nd20soderr);
  if (nd20preserr > 0)
    printf ("Record 20 pressure in error %d times\n", nd20preserr);
  if (nd20temperr > 0)
    printf ("Record 20 temperature in error %d times\n", nd20temperr);
  if (nd20humerr > 0)
    printf ("Record 20 humidity in error %d times\n", nd20humerr);
  if (nd20voerr > 0)
    printf ("Record 20 origin of values in error %d times\n", nd20voerr);

  if (nd20soderr > 0 || nd20preserr > 0 || nd20temperr > 0 ||
    nd20humerr > 0 || nd20voerr > 0)
    fail = 1;

  if (nd21soderr > 0)
    printf ("Record 21 second of day in error %d times\n", nd21soderr);
  if (nd21wserr > 0)
    printf ("Record 21 wind speed in error %d times\n", nd21wserr);
  if (nd21wderr > 0)
    printf ("Record 21 wind direction in error %d times\n", nd21wderr);
  if (nd21viserr > 0)
    printf ("Record 21 visibility in error %d times\n", nd21viserr);
  if (nd21scerr > 0)
    printf ("Record 21 sky clarity in error %d times\n", nd21scerr);
  if (nd21aserr > 0)
    printf ("Record 21 atmospheric seeing in error %d times\n", nd21aserr);
  if (nd21ccerr > 0)
    printf ("Record 21 cloud coverage in error %d times\n", nd21ccerr);

  if (nd21soderr > 0 || nd21wserr > 0 || nd21wderr > 0 || nd21viserr > 0 ||
    nd21scerr > 0 || nd21aserr > 0 || nd21ccerr > 0)
    fail = 1;

  if (nd30soderr > 0)
    printf ("Record 30 second of day in error %d times\n", nd30soderr);
  if (nd30azerr > 0)
    printf ("Record 30 azimuth in error %d times\n", nd30azerr);
  if (nd30elerr > 0)
    printf ("Record 30 elevation in error %d times\n", nd30elerr);
  if (nd30dierr > 0)
    printf ("Record 30 direction indicator in error %d times\n", nd30dierr);
  if (nd30aoierr > 0)
    printf ("Record 30 angle origin indicator in error %d times\n",
	    nd30aoierr);
  if (nd30rcierr > 0)
    printf ("Record 30 refraction correction indicator in error %d times\n",
	    nd30rcierr);

  if (nd30soderr > 0 || nd30azerr > 0 || nd30elerr > 0 || nd30dierr > 0 ||
    nd30aoierr > 0 || nd30rcierr > 0)
    fail = 1;

  if (nd40soderr > 0)
    printf ("Record 40 second of day in error %d times\n", nd40soderr);
  if (nd40toderr > 0)
    printf ("Record 40 type of data in error %d times\n", nd40toderr);
  if (nd40nprerr > 0)
    printf ("Record 40 number of points recorded in error %d times\n",
	    nd40nprerr);
  if (nd40npuerr > 0)
    printf ("Record 40 number of points used in error %d times\n",
	    nd40npuerr);
  if (nd40owtderr > 0)
    printf ("Record 40 one-way target distance in error %d times\n",
	    nd40owtderr);
  if (nd40owtderr > 0)
    printf ("Record 40 one-way target distance in error %d times\n",
	    nd40owtderr);
  if (nd40csderr > 0)
    printf ("Record 40 calibration system delay in error %d times\n",
	    nd40csderr);
  if (nd40cdserr > 0)
    printf ("Record 40 calibration delay shift in error %d times\n",
	    nd40cdserr);
  if (nd40crerr > 0)
    printf ("Record 40 calibration RMS in error %d times\n", nd40crerr);
  if (nd40ctierr > 0)
    printf ("Record 40 calibration type indicator in error %d times\n",
	    nd40ctierr);
  if (nd40cstierr > 0)
    printf ("Record 40 calibration shift type indicator in error %d times\n",
	    nd40cstierr);

  if (nd40soderr > 0 || nd40toderr > 0 || nd40nprerr > 0 || nd40npuerr > 0 ||
    nd40owtderr > 0 || nd40owtderr > 0 || nd40csderr > 0 || nd40cdserr > 0 ||
    nd40crerr > 0 || nd40ctierr > 0 || nd40cstierr > 0)
    fail = 1;

  if (nd50srerr > 0)
    printf ("Record 50 session RMS in error %d times\n", nd50srerr);
  if (nd50dqierr > 0)
    printf ("Record 50 data quality indicator in error %d times\n",
	    nd50dqierr);

  if (nd50srerr > 0 || nd50dqierr > 0) fail = 1;

  if (nd60scherr > 0)
    printf ("Record 50 system change indicator in error %d times\n",
	    nd60scherr);
  if (nd60scierr > 0)
    printf ("Record 50 system configuration indicator in error %d times\n",
	    nd60scierr);

  if (nd60scherr > 0 || nd60scierr > 0) fail = 1;

  if (tnseq > 0)
    printf ("Records not in sequence in %d passes\n", tnseq);

  if (tnseq > 0 || nseq > 0) fail = 1;

  printf ("\n\n");

  if (fail == 0) exit (0);
  else exit (1);
}

int
get_sat_ids (int mode, int *cospar_id, int *norad_id, int *sic, char *target,
	     int *target_type)
{
  FILE *sat_id_in;
  //char *sat_id_file = "/data/lib/targets.dat";
  char *sat_id_file = "targets.dat";
  char str[256], ttarget[11];
  int status, tcospar_id, tsic, tnorad_id = 0, ttt, trti;
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
      return (-1);
      //printf ("Could not open file %s\n", sat_id_file);
      //exit (1);
    }
  while ((status = fgets (str, 256, sat_id_in)) != NULL)
    {
      sscanf (str, "%s %d %d %d %d %d", ttarget, &ttt, &tsic, &tcospar_id,
	      &tnorad_id, &trti);
      if (mode == 0 && tcospar_id == *cospar_id)
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
	  *cospar_id = tcospar_id;
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
	  *cospar_id = tcospar_id;
	  *sic = tsic;
	  strcpy (target, ttarget);
	  l = strlen (target);
	  target[10] = '\0';
	  break;
	}
      else if (mode == 3 && strcmp (ttarget, target) == 0)
	{
	  *cospar_id = tcospar_id;
	  *norad_id = tnorad_id;
	  *sic = tsic;
	  break;
	}
    }
  *target_type = ttt;
  /*printf ("%s %d %d %d %d %d\n", ttarget, ttt, tsic, tcospar_id,
     tnorad_id, trti); */
  fclose (sat_id_in);

  if ((mode == 0 && tcospar_id != *cospar_id) ||
      (mode == 1 && tsic != *sic) ||
      (mode == 2 && tnorad_id != *norad_id) ||
      (mode == 3 && strcmp (ttarget, target) != 0))
    {
      return (1);
    }
  else
    {
      return (0);
    }
}

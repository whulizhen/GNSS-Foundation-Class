#include <stdio.h>
#include "cstg.h"

static int first_sed= 1;
void check_sum ();

/*-------------------------------------------------------------------------
 * Subroutines:  write old ILRS (CSTG) normalpoint and sampled
 *               engineering data format records
 *
 * Author: Randall Ricklefs / Univ of Texas / Center for Space Research
 *
 * History:
 *   Nov 06, 2007 - Initial version
 *
**-----------------------------------------------------------------------*/

/*-------------------------------------------------------------------------
**
**      write_cstg_hdr - write header record for old ilrs (cstg)
**                       normalpoint and sampled engineering format
**
**-----------------------------------------------------------------------*/
write_cstg_hdr (FILE *str_out, struct cstg_hdr hdr, int file_type)
{
  char str[60];

  if (file_type == 0)
    fprintf(str_out, "88888\n");
  else if (file_type == 1)
    fprintf(str_out, "99999\n");
  sprintf(str,
    "%07d%02d%03d%04d%02d%02d%04d%08d%06d%04d%1d%1d%1d%1d%1d%04d%1d%02d%1d\n",
    hdr.ilrs_id, hdr.year, hdr.doy, hdr.cdp_pad_id, hdr.cdp_sys_num, 
    hdr.cdp_occ_num, (int)hdr.xmit_wavelength, (int)hdr.cal_sys_delay, 
    (int)hdr.cal_delay_shift, (int)hdr.cal_rms, hdr.np_window_ind, 
    hdr.stn_timescale, hdr.cal_type_ind, hdr.sys_change_ind, 
    hdr.sys_config_ind, (int)hdr.sess_rms, hdr.data_qual_ind, hdr.checksum, 
    hdr.format_version);
  check_sum(str,54);
  fputs(str, str_out);
}

/*-------------------------------------------------------------------------
**
**      write_cstg_np - write old ilrs (cstg) normalpoint record
**
**-----------------------------------------------------------------------*/
write_cstg_np (FILE *str_out, struct cstg_np np)
{
  char str[60];

  sprintf(str,
    "%012.0Lf%012.0Lf%07d%05d%04d%03d%04d%1d%1d%1d%02d%02d\n",
    np.sec_of_day, np.time_of_flight, (int)np.bin_rms, 
    (int)np.pressure, (int)np.temperature, (int)np.humidity,
    np.num_ranges, np.data_release, np.scale_or_tof_sec, np.llr_np_window_ind,
    (int)np.snr, np.checksum);
  check_sum(str,54);
  fputs(str, str_out);
}

/*-------------------------------------------------------------------------
**
**      write_cstg_sed - write old ilrs (cstg) sampled engineering record
**
**-----------------------------------------------------------------------*/
write_cstg_sed (FILE *str_out, struct cstg_sed sed)
{
  char str[70];

  sprintf(str,
    "%012.0Lf%012.0Lf%05d%04d%03d%08d%04d%01d%07d%06d%5s%02d\n",
    sed.sec_of_day, sed.time_of_flight,
    (int)sed.pressure, (int)sed.temperature, (int)sed.humidity,
    (int)sed.internal_burst_cal, sed.xcv_amp, 
    sed.angle_origin_ind, (int)sed.azimuth, (int)sed.elevation,
    "00000", sed.checksum);
  check_sum(str,69);
  fputs(str, str_out);
}

/*********************************************************************
*
*PURPOSE
*       DETERMINE THE CHECKSUM FOR ONE NORMAL POINT DATA RECORD
*       AND PLACE THAT SUM IN THE RECORD
*
*INPUT VARIABLE
*       LINE CHARACTER*(*)
*       LINE_LEN INTEGER*4
*OUTPUT VARIABLE
*       LINE CHARACTER*(*)
*
*       BY BRION CONKLIN BENDIX FIELD ENGINEERING DSG/OAS 4/90
*	Converted to 'c' by R. Ricklefs UT/CSR 3 July 2007
*
*********************************************************************/
void
check_sum(char *line,int line_len)
{
        char check[3];
        int i, sum;

        sum = 0;
        for (i=0;i<line_len-2;i++)
	  {
             if (line[i] == '1') sum+= 1;
             if (line[i] == '2') sum+= 2;
             if (line[i] == '3') sum+= 3;
             if (line[i] == '4') sum+= 4;
             if (line[i] == '5') sum+= 5;
             if (line[i] == '6') sum+= 6;
             if (line[i] == '7') sum+= 7;
             if (line[i] == '8') sum+= 8;
             if (line[i] == '9') sum+= 9;
          }
        sum%= 100;
        sprintf(check,"%02d", sum);
        line[line_len-2]= check[0];
        line[line_len-1]= check[1];
}

#include <stdio.h>
#include "cstg.h"

/*-------------------------------------------------------------------------
 * Subroutines: read old ilrs normalpoint and sampled enginnering data files
 *
 * Author: Randall Ricklefs / Univ of Texas / Center for Space Research
 *
 * History:
 *   Nov 06, 2007 - Initial version
 *   Feb 2, 2010 - Because of round-off problems, change 1.e7 to
 *                 (long double)1.e7 in second of days calulations.
 *                 1.01b. rlr.
 *
**-----------------------------------------------------------------------*/
/*-------------------------------------------------------------------------
**
**      read_cstg_hdr - read fields from old-format ilrs normalpoint and
**                      sampled engineering data header
**
**-----------------------------------------------------------------------*/
void
read_cstg_hdr (char *str, struct cstg_hdr *hdr)
{
  long syear;

  getifield (str, 0, 7, &hdr->ilrs_id);
  getifield (str, 7, 2, &syear);
  if (syear > 50)
    hdr->year = syear + 1900;
  else
    hdr->year = syear + 2000;
  getifield (str, 9, 3, &hdr->doy);
  getifield (str, 12, 4, &hdr->cdp_pad_id);
  getifield (str, 16, 2, &hdr->cdp_sys_num);
  getifield (str, 18, 2, &hdr->cdp_occ_num);
  getfield (str, 20, 4, &hdr->xmit_wavelength);
  if (hdr->xmit_wavelength < 1000 || hdr->xmit_wavelength > 2999)
    (hdr->xmit_wavelength) /= 10;
  getfield (str, 24, 8, &hdr->cal_sys_delay);
  getfield (str, 32, 6, &hdr->cal_delay_shift);
  getfield (str, 38, 4, &hdr->cal_rms);
  getifield (str, 42, 1, &hdr->np_window_ind);
  getifield (str, 43, 1, &hdr->stn_timescale);
  getifield (str, 44, 1, &hdr->cal_type_ind);
  getifield (str, 45, 1, &hdr->sys_change_ind);
  getifield (str, 46, 1, &hdr->sys_config_ind);
  getfield (str, 47, 4, &hdr->sess_rms);
  getifield (str, 51, 1, &hdr->data_qual_ind);
  getifield (str, 52, 2, &hdr->checksum);
  getifield (str, 54, 1, &hdr->format_version);
}

/*-------------------------------------------------------------------------
**
**      read_cstg_np - read fields from old-format ilrs normalpoint 
**                     data record
**
**-----------------------------------------------------------------------*/
void
read_cstg_np (char *str, struct cstg_np *np)
{
  double dtemp;

  getfield (str, 0, 12, &dtemp);
  np->sec_of_day= dtemp/(long double)1.e7;
  getfield (str, 12, 12, &dtemp);
  np->time_of_flight= dtemp/1.e12;
  getfield (str, 24, 7, &np->bin_rms);
  getfield (str, 31, 5, &np->pressure);
  np->pressure /= 10.;
  getfield (str, 36, 4, &np->temperature);
  (np->temperature) /= 10.;
  getfield (str, 40, 3, &np->humidity);
  getifield (str, 43, 4, &np->num_ranges);
  getifield (str, 47, 4, &np->data_release);
  getifield (str, 48, 1, &np->scale_or_tof_sec);
  getifield (str, 49, 1, &np->llr_np_window_ind);
  getifield (str, 50, 2, &np->snr);
  getifield (str, 52, 2, &np->checksum);
}

/*-------------------------------------------------------------------------
**
**      read_cstg_np - read fields from old-format ilrs sampled engineering 
**                     data record
**
**-----------------------------------------------------------------------*/
void
read_cstg_sampled_eng (char *str, struct cstg_sed *sed)
{
  double dtemp;

  getfield (str, 0, 12, &dtemp);
  sed->sec_of_day= dtemp/(long double)1.e7;
  getfield (str, 12, 12, &dtemp);
  sed->time_of_flight= dtemp/1.e12;
  getfield (str, 24, 5, &sed->pressure);
  sed->pressure /= 10.;
  getfield (str, 29, 4, &sed->temperature);
  sed->temperature /= 10.;
  getfield (str, 33, 3, &sed->humidity);
  getifield (str, 36, 8, &sed->internal_burst_cal);
  getifield (str, 44, 4, &sed->xcv_amp);
  getifield (str, 48, 1, &sed->angle_origin_ind);
  getfield (str, 49, 7, &sed->azimuth);
  sed->azimuth /= 1.e4;
  getfield (str, 56, 6, &sed->elevation);
  sed->elevation /= 1.e4;
  getifield (str, 67, 2, &sed->checksum);
}

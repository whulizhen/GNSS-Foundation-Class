#include <stdio.h>
#include "merit.h"

/*-------------------------------------------------------------------------
 * Subroutine: ILRS (MERIT II) fullrate reading routines.
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
**      read_merit_fr - Read fields from old ILRS fullrate data record.
**			Format was previsously known as Merit II.
**
**-----------------------------------------------------------------------*/
void
read_merit_fr (char *str, struct merit_fr *fr)
{
  long syear;
  double dtemp;

  getifield (str, 0, 7, &fr->ilrs_id);
  getifield (str, 7, 2, &syear);
  if (syear > 50)
    fr->year = syear+ 1900;
  else
    fr->year = syear+ 2000;
  getifield (str, 9, 3, &fr->doy);
  getfield (str, 12, 12, &dtemp);
  fr->sec_of_day = dtemp/(long double)1.e7;
  getifield (str, 24, 4, &fr->cdp_pad_id);
  getifield (str, 28, 2, &fr->cdp_sys_num);
  getifield (str, 30, 2, &fr->cdp_occ_num);
  getfield (str, 32, 7, &fr->azimuth);
  fr->azimuth /= 1.e4;
  getfield (str, 39, 6, &fr->elevation);
  fr->elevation /= 1.e4;
  getfield (str, 45, 12, &dtemp);
  fr->time_of_flight = dtemp/1.e12;
  getfield (str, 57, 7, &fr->sess_rms);
  getfield (str, 64, 4, &fr->xmit_wavelength);
  if (fr->xmit_wavelength < 1000 || fr->xmit_wavelength > 2999)
    fr->xmit_wavelength /= 10;
  getfield (str, 68, 5, &fr->pressure);
  fr->pressure /= 10.;
  getfield (str, 73, 4, &fr->temperature);
  fr->temperature /= 10.;
  getfield (str, 77, 3, &fr->humidity);
  getfield (str, 80, 5, &fr->refraction_corr);
  getfield (str, 85, 6, &fr->target_CofM_corr);
  getifield (str, 91, 5, &fr->xcv_amp);
  getfield (str, 96, 8, &fr->cal_sys_delay);
  getfield (str, 104, 6, &fr->cal_delay_shift);
  getfield (str, 110, 4, &fr->cal_rms);
  getifield (str, 114, 4, &fr->np_window_ind);
  getifield (str, 115, 4, &fr->num_ranges);
  getifield (str, 119, 1, &fr->epoch_event);
  getifield (str, 120, 1, &fr->stn_timescale);
  getifield (str, 121, 1, &fr->angle_origin_ind);
  getifield (str, 122, 1, &fr->refraction_app_ind);
  getifield (str, 123, 1, &fr->CofM_app_ind);
  getifield (str, 124, 1, &fr->xcv_amp_app_ind);
  getifield (str, 125, 1, &fr->cal_type_ind);
  getifield (str, 126, 1, &fr->sys_change_ind);
  getifield (str, 127, 1, &fr->sys_config_ind);
  getifield (str, 128, 1, &fr->format_version);
  getifield (str, 129, 1, &fr->data_release);
}

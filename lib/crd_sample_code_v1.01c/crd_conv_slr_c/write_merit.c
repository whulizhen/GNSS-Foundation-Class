#include <stdio.h>
#include "merit.h"

/*-------------------------------------------------------------------------
 * Subroutine:  write old ILRS (CSTG) fullrate data format records
 *
 * Author: Randall Ricklefs / Univ of Texas / Center for Space Research
 *
 * History:
 *   Nov 06, 2007 - Initial version
 *
**-----------------------------------------------------------------------*/
/*-------------------------------------------------------------------------
**
**      write_merit_fr - write old ILRS fullrate format records (previously
**			 known as the Merit II format)
**
**-----------------------------------------------------------------------*/
void
write_merit_fr (FILE *str_out, struct merit_fr fr)
{
  int syear;
  double dtemp;
  char str[256];

  sprintf (str,"%7d%02d%3d%12.0Lf%4d%2d%2d%7d%6d%12.0Lf%7d%4d%5d%4d%3d%5d%6d%5d%8d%6d%4d%1d%4d%1d%1d%1d%1d%1d%1d%1d%1d%1d%1d%1d\n",
    fr.ilrs_id,
    fr.year%100, fr.doy, fr.sec_of_day,
    fr.cdp_pad_id, fr.cdp_sys_num, fr.cdp_occ_num,
    (int)(fr.azimuth), (int)(fr.elevation),
    fr.time_of_flight,
    (int)(fr.sess_rms),
    (int)(fr.xmit_wavelength),
    (int)(fr.pressure), (int)(fr.temperature), (int)(fr.humidity),
    (int)fr.refraction_corr, (int)(fr.target_CofM_corr), fr.xcv_amp, 
    (int)fr.cal_sys_delay, (int)fr.cal_delay_shift, (int)fr.cal_rms, 
    fr.np_window_ind, fr.num_ranges,
    fr.epoch_event, fr.stn_timescale, fr.angle_origin_ind, 
    fr.refraction_app_ind, fr.CofM_app_ind, fr.xcv_amp_app_ind, 
    fr.cal_type_ind, fr.sys_change_ind, fr.sys_config_ind, 
    fr.format_version, fr.data_release);
  fputs(str,str_out);
}

#include <stdio.h>
#include <string.h>
#include "cospar_llr.h"

/*-------------------------------------------------------------------------
 * Subroutines: read old LLR COSPAR format - fullrate, normalpoint, and
 * 		mini-normalpoint data types.
 *
 * Author: Randall Ricklefs / Univ of Texas / Center for Space Research
 *
 * History:
 *   Nov 06, 2007 - Initial version
 *
**-----------------------------------------------------------------------*/

/* System Configuration Record */
void
read_llr0 (char *str, struct rllr0 *ld)
{

  getifield (str, 1, 1, &ld->laser_color);
  getifield (str, 2, 4, &ld->year);
  getifield (str, 6, 2, &ld->month);
  getifield (str, 8, 2, &ld->day);
  getifield (str, 10, 2, &ld->hour);
  getifield (str, 12, 2, &ld->min);
  getifield (str, 14, 2, &ld->subsys_changed);
  getifield (str, 16, 1, &ld->change_type);
  getifield (str, 17, 4, &ld->change_number);
  strncpy (ld->change_description,&str[21],68);
  ld->data_file_version= str[89];
}
 
/* Run Header */
void
read_llr1 (char *str, struct rllr1 *ld)
{

  getifield (str, 1, 1, &ld->laser_color);
  getifield (str, 2, 4, &ld->year);
  getifield (str, 6, 2, &ld->month);
  getifield (str, 8, 2, &ld->day);
  getifield (str, 10, 2, &ld->hour);
  getifield (str, 12, 2, &ld->min);
  getifield (str, 15, 3, &ld->doy);
  getifield (str, 19, 5, &ld->mjd);
  getifield (str, 24, 5, &ld->observatory_code);
  getifield (str, 29, 5, &ld->laser_wavelength);
  getifield (str, 34, 3, &ld->dark_count);
  getifield (str, 37, 3, &ld->lunar_site_count);
  getifield (str, 40, 3, &ld->star_count);
  strncpy (ld->star_name,&str[43],5);
  getifield (str, 48, 1, &ld->epoch_time_base);
  getifield (str, 49, 1, &ld->delay_time_base);
  getifield (str, 50, 2, &ld->reflector);
  ld->data_quality_ind= str[52];
  ld->detector= str[53];
  ld->data_file_version= str[89];
}
  
/* Run Sub-header */
void
read_llr2 (char *str, struct rllr2 *ld)
{

  getifield (str, 1, 1, &ld->laser_color);
  getifield (str, 2, 4, &ld->year);
  getifield (str, 6, 2, &ld->month);
  getifield (str, 8, 2, &ld->day);
  getifield (str, 10, 2, &ld->hour);
  getifield (str, 12, 2, &ld->min);
  getifield (str, 14, 2, &ld->laser_energy);
  getifield (str, 16, 6, &ld->laser_pulse_length);
  getifield (str, 22, 6, &ld->shot_by_shot_resolution);
  getifield (str, 28, 4, &ld->spectral_filter);
  getifield (str, 32, 4, &ld->spacial_filter);
  getifield (str, 36, 5, &ld->pmt_voltage);
  getifield (str, 41, 5, &ld->shots_out);
  getifield (str, 46, 3, &ld->seeing);
  getifield (str, 49, 4, &ld->temperature);
  getifield (str, 53, 2, &ld->humidity);
  getifield (str, 55, 2, &ld->wind_speed);
  strncpy (ld->wind_direction,&str[57],2);
  getifield (str, 59, 8, &ld->clock_offset);
  ld->data_file_version= str[89];
}
 
/* Detail Record */
void
read_llr3 (char *str, struct rllr3 *ld)
{

  long tof_u, tof_l;

  getifield (str, 1, 1, &ld->laser_color);
  getifield (str, 2, 4, &ld->year);
  getifield (str, 6, 2, &ld->month);
  getifield (str, 8, 2, &ld->day);
  getifield (str, 10, 2, &ld->hour);
  getifield (str, 12, 2, &ld->min);
  getifield (str, 14, 9, &ld->sec);
  getifield (str, 23, 9, &tof_u);
  getifield (str, 32, 5, &tof_l);
  ld->time_of_flight= tof_u*1.e5+ tof_l;
  getifield (str, 37, 1, &ld->vernier);
  getifield (str, 38, 9, &ld->electronic_delay);
  getifield (str, 47, 8, &ld->geometric_delay);
  getifield (str, 55, 6, &ld->uncert_estimate);
  getifield (str, 61, 7, &ld->clock_drift);
  getifield (str, 68, 6, &ld->pressure);
  getifield (str, 84, 1, &ld->filter_flag);
  ld->data_file_version= str[89];
  getifield (str, 90, 10, &ld->azimuth);
  getifield (str, 100, 10, &ld->elevation);
  getifield (str, 110, 9, &ld->range_residual);
}
 
/* Normal Point Record */
void
read_llr4 (char *str, struct rllr4 *ld)
{

  long tof_u, tof_l;

  getifield (str, 1, 1, &ld->laser_color);
  getifield (str, 2, 4, &ld->year);
  getifield (str, 6, 2, &ld->month);
  getifield (str, 8, 2, &ld->day);
  getifield (str, 10, 2, &ld->hour);
  getifield (str, 12, 2, &ld->min);
  getifield (str, 14, 9, &ld->sec);
  getifield (str, 23, 9, &tof_u);
  getifield (str, 32, 5, &tof_l);
  ld->time_of_flight= tof_u*1.e5+ tof_l;
  getifield (str, 37, 1, &ld->vernier);
  getifield (str, 38, 9, &ld->electronic_delay);
  getifield (str, 47, 8, &ld->geometric_delay);
  getifield (str, 55, 6, &ld->uncert_estimate);
  getifield (str, 61, 7, &ld->clock_drift);
  getifield (str, 68, 6, &ld->pressure);
  getifield (str, 74, 3, &ld->number_of_returns);
  getifield (str, 77, 3, &ld->signal_to_noise);
  getifield (str, 80, 4, &ld->time_span);
  getifield (str, 84, 1, &ld->filter_flag);
  ld->data_file_version= str[89];
  getifield (str, 90, 10, &ld->azimuth);
  getifield (str, 100, 10, &ld->elevation);
  getifield (str, 110, 9, &ld->range_residual);
}
 
/* Mini Normal Point Record */
void
read_llr5 (char *str, struct rllr5 *ld)
{

  long tof_u, tof_l;

  getifield (str, 1, 1, &ld->laser_color);
  getifield (str, 2, 4, &ld->year);
  getifield (str, 6, 2, &ld->month);
  getifield (str, 8, 2, &ld->day);
  getifield (str, 10, 2, &ld->hour);
  getifield (str, 12, 2, &ld->min);
  getifield (str, 14, 9, &ld->sec);
  getifield (str, 23, 9, &tof_u);
  getifield (str, 32, 5, &tof_l);
  ld->time_of_flight= tof_u*1.e5+ tof_l;
  getifield (str, 37, 1, &ld->reflector);
  getifield (str, 38, 5, &ld->observatory_code);
  getifield (str, 43, 3, &ld->number_of_returns);
  getifield (str, 46, 6, &ld->uncert_estimate);
  getifield (str, 52, 3, &ld->signal_to_noise);
  ld->data_quality_ind= str[55];
  getifield (str, 56, 6, &ld->pressure);
  getifield (str, 62, 4, &ld->temperature);
  getifield (str, 66, 2, &ld->humidity);
  getifield (str, 68, 5, &ld->laser_wavelength);
  ld->data_file_version= str[73];
  getifield (str, 74, 4, &ld->time_span);
  ld->detector= str[78];
}
 
/* Calibration Record */
void
read_llr6 (char *str, struct rllr6 *ld)
{
  int syear;

  getifield (str, 1, 1, &ld->laser_color);
  getifield (str, 2, 4, &ld->year);
  getifield (str, 6, 2, &ld->month);
  getifield (str, 8, 2, &ld->day);
  getifield (str, 10, 2, &ld->hour);
  getifield (str, 12, 2, &ld->min);
  getifield (str, 14, 9, &ld->sec);
  getifield (str, 23, 1, &ld->type_of_cal);
  getifield (str, 24, 6, &ld->calMcalaverage[0]);
  getifield (str, 30, 6, &ld->calMcalaverage[1]);
  getifield (str, 36, 6, &ld->calMcalaverage[2]);
  getifield (str, 42, 6, &ld->calMcalaverage[3]);
  getifield (str, 48, 6, &ld->calMcalaverage[4]);
  getifield (str, 54, 6, &ld->calMcalaverage[5]);
  getifield (str, 60, 6, &ld->calMcalaverage[6]);
  getifield (str, 66, 6, &ld->calMcalaverage[7]);
  getifield (str, 72, 6, &ld->calMcalaverage[8]);
  getifield (str, 78, 6, &ld->calMcalaverage[9]);
  ld->data_file_version= str[89];
}
 
/* Comment Record */
void
read_llr7 (char *str, struct rllr7 *ld)
{
  int syear;

  getifield (str, 1, 1, &ld->laser_color);
  getifield (str, 2, 4, &ld->year);
  getifield (str, 6, 2, &ld->month);
  getifield (str, 8, 2, &ld->day);
  getifield (str, 10, 2, &ld->hour);
  getifield (str, 12, 2, &ld->min);
  getifield (str, 14, 1, &ld->source);
  strncpy (ld->comment,&str[15],74);
  ld->comment[73]= '\0';
  ld->data_file_version= str[89];
printf("comment: [%s] %s]\n",str,ld->comment);
}
 

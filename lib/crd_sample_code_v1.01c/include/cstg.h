struct cstg_hdr
{
  int ilrs_id;
  int year;
  int doy;
  int cdp_pad_id;
  int cdp_sys_num;
  int cdp_occ_num;
  double xmit_wavelength;
  double cal_sys_delay;
  double cal_delay_shift;
  double cal_rms;
  int np_window_ind;
  int stn_timescale;
  int cal_type_ind;
  int sys_change_ind;
  int sys_config_ind;
  double sess_rms;
  int data_qual_ind;
  int checksum;
  int format_version;
};

struct cstg_np
{
  long double sec_of_day;
  long double time_of_flight;
  double bin_rms;
  double pressure;
  double temperature;
  double humidity;
  int num_ranges;
  int data_release;
  int scale_or_tof_sec;
  int llr_np_window_ind;
  double snr;
  int checksum;
};

struct cstg_sed
{
  long double sec_of_day;
  long double time_of_flight;
  double pressure;
  double temperature;
  double humidity;
  double internal_burst_cal;
  int xcv_amp;	/* really? */
  int angle_origin_ind;
  double azimuth;
  double elevation;
  char fill1[4];
  int checksum;
};

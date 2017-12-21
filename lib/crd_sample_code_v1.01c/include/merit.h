struct merit_fr
{
  int ilrs_id;
  int year;
  int doy;	/* Day of year */
  long double sec_of_day; 
  int cdp_pad_id;
  int cdp_sys_num;
  int cdp_occ_num;
  double azimuth;
  double elevation;
  long double time_of_flight;
  double sess_rms;
  double xmit_wavelength;
  double pressure;
  double temperature;
  double humidity;
  double refraction_corr;
  double target_CofM_corr;
  int xcv_amp;
  double cal_sys_delay;
  double cal_delay_shift;
  double cal_rms;
  int np_window_ind;
  int num_ranges;
  int epoch_event;
  int stn_timescale;
  int angle_origin_ind;
  int refraction_app_ind;
  int CofM_app_ind;
  int xcv_amp_app_ind;
  int cal_type_ind;
  int sys_change_ind;
  int sys_config_ind;
  int format_version;
  int data_release;
};

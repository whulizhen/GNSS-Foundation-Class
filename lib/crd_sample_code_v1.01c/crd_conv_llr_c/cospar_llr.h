/* Lunar Laser Ranging data header fields */
/* System Configuration Record */
  struct rllr0 {
    int laser_color;
    int year;
    int month;
    int day;
    int hour;
    int min;
    int subsys_changed;
    int change_type;
    int change_number;
    char change_description[69];
    char data_file_version;
  };

/* Run Header */
  struct rllr1 {
    int laser_color;
    int year;
    int month;
    int day;
    int hour;
    int min;
    int doy;
    int mjd;
    int observatory_code;
    int laser_wavelength;
    int dark_count;
    int lunar_site_count;
    int star_count;
    char star_name[6];
    int epoch_time_base;
    int delay_time_base;
    int reflector;
    int data_quality_ind;
    int detector;
    char data_file_version;
  };

/* Run Sub-header */
  struct rllr2 {
    int laser_color;
    int year;
    int month;
    int day;
    int hour;
    int min;
    int laser_energy;
    int laser_pulse_length;
    int shot_by_shot_resolution;
    int spectral_filter;
    int spacial_filter;
    int pmt_voltage;
    int shots_out;
    int seeing;
    int temperature;
    int humidity;
    int wind_speed;
    char wind_direction[3];
    int clock_offset;
    char data_file_version;
  };

/* Detail Record */
  struct rllr3 {
    int laser_color;
    int year;
    int month;
    int day;
    int hour;
    int min;
    long sec;
    double time_of_flight;
    int vernier;
    int electronic_delay;
    int geometric_delay;
    int uncert_estimate;
    int clock_drift;
    int pressure;
    int filter_flag;
    char data_file_version;
    long azimuth;
    long elevation;
    long range_residual;
  };

/* Normal Point Record */
  struct rllr4 {
    int laser_color;
    int year;
    int month;
    int day;
    int hour;
    int min;
    long sec;
    double time_of_flight;
    int vernier;
    int electronic_delay;
    int geometric_delay;
    int uncert_estimate;
    int clock_drift;
    int pressure;
    int number_of_returns;
    int signal_to_noise;
    int time_span;
    int filter_flag;
    char data_file_version;
    long azimuth;
    long elevation;
    long range_residual;
  };

/* Mini Normal Point Record */
  struct rllr5 {
    int laser_color;
    int year;
    int month;
    int day;
    int hour;
    int min;
    long sec;
    double time_of_flight;
    int reflector;
    int observatory_code;
    int number_of_returns;
    int uncert_estimate;
    int signal_to_noise;
    char data_quality_ind;
    int pressure;
    int temperature;
    int humidity;
    int laser_wavelength;
    char data_file_version;
    int time_span;
    int detector;
/**
    long azimuth;
    long elevation;
    long range_residual;
**/
  };

/* Calibration Record */
  struct rllr6 {
    int laser_color;
    int year;
    int month;
    int day;
    int hour;
    int min;
    int sec;
    int type_of_cal;
    int calMcalaverage[10];
    char data_file_version;
  };

/* Comment Record */
  struct rllr7 {
    int laser_color;
    int year;
    int month;
    int day;
    int hour;
    int min;
    int source;
    char comment[75];
    char data_file_version;
  };

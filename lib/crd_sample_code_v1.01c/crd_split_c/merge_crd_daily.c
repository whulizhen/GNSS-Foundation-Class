#include <stdio.h>
#ifdef BSD
#include <sys/dir.h>
#else
#include <dirent.h>
#endif
#include <sys/types.h>
#include <sys/stat.h>
#include <time.h>
#include <unistd.h>
#include "../include/crd.h"


/* The directory where the input CRD files reside: */
#define DATA_IN_DIR "merge_test_files/"

/* The directory where the output CRD files will be placed: */
#define DATA_DIR "merge_test_files/crd_merge/"

struct rh1 h1;
struct file_date
{
  int year;
  int mon;
  int day;
  int hour;
};

void get_pred_info ();
int compar_all (), compar_sic ();	/* programmer decides which to use */
double get_time ();
void get_sat_ids_short ();
void copy_file ();
void read_h1 ();
void write_h1 ();
void write_h9 ();

/*-------------------------------------------------------------------------
**
** Program: merge_crd_daily - 
**
** Purpose:
** Create daily CRD normal point (.npt), sampled engineering (quicklook, .qlk),
** and full rate (.frd) files from individual pass files. Input files names
** are expected to follow either the ILRS CRD station single pass prototcol
** or the NASA station/HTSI single pass naming convention. Output files are 
** written using the ILRS CRD data center naming convention. 
** See the ILRS CRD format document for the details.
**
** NOTE: a small change to the MJD test below could allow this program to 
** merge files by hour.
**
** Calling sequence:
**   merge_crd_daily
**
** All CRD .npt, and related files are taken from DATA_IN_DIR and written to
** DATA_DIR.
**
** Author: Randall Ricklefs / Univ of Texas / Center for Space Research
**
** History:
**   24 Jan 2008 - Initial version
**
**-----------------------------------------------------------------------*/
main ()
{

  DIR *dirp;			/*  directory pointer   */
#ifdef BSD
  struct direct *item;		/*  entry in directory  */
#else
  struct dirent *item;		/*  entry in directory  */
#endif
  struct stat buf;
  int i, ii;
  int first = 1, nf = 0, nitems = 0, nstr, mjd, last_mjd = -1;
  int doy, hour, minute, idum;
  long year, mon, day;
  int written_to_npt[200], written_to_qlk[200], written_to_frd[200];
  int std_naming, station_id, hr, rel;
  long ss, yy, ddd, hhmm, sic;
  long dt;
  char filename[256], filename_ext[256], files[256][200], satname[11];
  char passname[5000][70];	/* Selected passes' names */
  char satpassname[70], sys_str[256], str[120];
  char crd_lit[10], mon_name[4], testname[256], testname1[256];
  double dum, jdi, jdf, tdate, curr_date;
  FILE *str_in;
  FILE *str_out_npt[200], *str_out_qlk[200], *str_out_frd[200];
  struct file_date file_prod;


  /* Get the current system time for the 'H1' header. */
  get_time (&file_prod.year, &file_prod.mon, &file_prod.day, &file_prod.hour);

  /* Open the input directory */
  if ((dirp = opendir (DATA_IN_DIR)) == 0)
    {
      sprintf (sys_str, "Cannot open directory %s", DATA_IN_DIR);
      return;
    }

  /*  read one filename each time through the loop  */
  for (item = readdir (dirp); item != NULL; item = readdir (dirp))
    {
      if (strcmp (item->d_name, ".") == 0)
	continue;
      if (strcmp (item->d_name, "..") == 0)
	continue;

      /* Look at normalpoint file, only now. */
      if (strncmp (&item->d_name[strlen(item->d_name)-4], ".npt", 4) != 0 &&
	  strncmp (&item->d_name[strlen(item->d_name)-4], ".NPT") != 0)
	continue;

      /* If the normalpoint file is empty, we don't want to use the pass */
      strcpy (filename_ext, DATA_IN_DIR);
      strcat (filename_ext, item->d_name);
      stat (filename_ext, &buf);
      if (buf.st_size <= 0)
	continue;		/* silently ignore  empty npt files */

      /* What do we do when there are too many passes? */
      if (nitems >= 5000)
	continue;		/* silently ignore the rest */

      /* drop the extension but not the '.' */
      strncpy (testname, item->d_name, strlen (item->d_name) - 3);
      testname[strlen (item->d_name) - 3]= '\0';
      /* Standrd CRD station file naming convention, or NASA station? */
      std_naming = !(testname[0] == 's' && testname[3] == 'y');
      if (std_naming)
        {
	  for (ii=0; ii<strlen(testname); ii++)
	    {
	      if (testname[ii] == '_') testname1[ii]= ' ';
	      else if (testname[ii] == '.') testname1[ii]= ' ';
	      else testname1[ii]= testname[ii];
	    }
	  testname1[strlen(testname)]= '\0';
          sscanf(testname1,"%d %s %s %ld %d %d",
		&station_id,satname,crd_lit,&dt,&hr,&rel);

          sprintf (&passname[(nitems)++], "%8d %02d %s %04d %s", dt,
               hr, satname, station_id, testname);
	}
      else	/* NASA: sSSyYYdDDDtHHMM#SICX.ttt */
        {
	  getifield(testname, 1, 2, &ss);
	  getifield(testname, 4, 2, &yy);
	  getifield(testname, 7, 3, &ddd);
	  getifield(testname, 11, 4, &hhmm);
          if (strlen (testname) == 21)
	    getifield (testname, 16, 4, &sic);
          else if (strlen (testname) == 20)
	    getifield (testname, 16, 3, &sic);

          get_sat_ids_short (1, &idum, &idum, &sic, satname, &idum);
          doytogr (yy+100, ddd, &mon, &day, mon_name); 
          sprintf (&passname[(nitems)++], "%4d%02d%02d %02d %s %4d %s\0", 
	       yy+2000, mon,day, hhmm/100, satname, ss, testname);
	}
    }
  closedir (dirp);

  /* Sort qualifying passes using compar_all */
  qsort (passname, (size_t) nitems, (size_t) (70), compar_all);
   /** for (i = 0; i < nitems; i++)
    printf ("[%s]\n", passname[i]); **/ 

  /* For each pass, figure out which file is to be copied to and copy it */
  for (i = 0; i < nitems; i++)
    {
      getifield (passname[i], 0, 4, &year);
      getifield (passname[i], 4, 2, &mon);
      getifield (passname[i], 6, 2, &day);
      sscanf (passname[i], "%*ld %d %s %*d %s",&hour,satname,satpassname);

      /* Check time of pass for file name */
      grtojd (year - 1900, mon, day, hour, 0, 0.e0, &jdi, &jdf);
      mjd = jdi - 2400000.e0 + jdf - 0.5e0;
      /*printf ("mjd %d last_mjd %d\n", mjd, last_mjd);*/

      /* New day (or first pass). Close all existing files and start over. */
      if (mjd > last_mjd)
	{
	  /* first pass */
	  if (last_mjd > 0)
	    {
	      for (ii = 0; ii < nf; ii++)
		{
		  if (written_to_npt[ii])
		    write_h9 (str_out_npt[ii]);

		  if (written_to_qlk[ii])
		    write_h9 (str_out_qlk[ii]);

		  if (written_to_frd[ii])
		    write_h9 (str_out_frd[ii]);

		  fclose (str_out_npt[ii]);
		  fclose (str_out_qlk[ii]);
		  fclose (str_out_frd[ii]);
		}
	    }
	  last_mjd = mjd;
	  nf = 0;
	}

     /* Create the base filename and open all 3 related file names for this 
      * satellite and day */
      sprintf (filename, "%s_%4d%02d%02d", satname, year, mon, day);
      /**printf ("Data center file name: %d [%s]\n", nf, filename);**/

      /* Are the output files already open? */
      nstr = nf;
      for (ii = 0; ii < nf; ii++)
	{
	  if (strcmp (filename, files[ii]) == 0)
	    nstr = ii;
	}

      /* No match - open new ouput files */
      if (nstr == nf)
	{
          /* Normalpoint */
	  written_to_npt[nf] = 0;
	  strcpy (files[nf], filename);
	  strcpy (filename_ext, DATA_DIR);
	  strcat (filename_ext, filename);
	  strcat (filename_ext, ".npt");
	  if ((str_out_npt[nf] = fopen (filename_ext, "w")) == NULL)
	    {
	      printf ("Could not open file %s\n", filename_ext);
	      exit (1);
	    }

          /* Sampled engineering (Quicklook) */
	  written_to_qlk[nf] = 0;
	  strcpy (filename_ext, DATA_DIR);
	  strcat (filename_ext, filename);
	  strcat (filename_ext, ".qlk");
	  if ((str_out_qlk[nf] = fopen (filename_ext, "w")) == NULL)
	    {
	      printf ("Could not open file %s\n", filename_ext);
	      exit (1);
	    }

          /* Full rate */
	  written_to_frd[nf] = 0;
	  strcpy (filename_ext, DATA_DIR);
	  strcat (filename_ext, filename);
	  strcat (filename_ext, ".frd");
	  if ((str_out_frd[nf++] = fopen (filename_ext, "w")) == NULL)
	    {
	      printf ("Could not open file %s\n", filename_ext);
	      exit (1);
	    }
	  first = 1;
	}

      /* Open and copy each file type for this pass */
      /* Normalpoint */
      strcpy (filename_ext, DATA_IN_DIR);
      strcat (filename_ext, satpassname);
      strcat (filename_ext, "npt");
      if ((str_in = fopen (filename_ext, "r")) == NULL)
	{
	  printf ("Could not open file %s\n", filename_ext);
	}
      else
	{
	  copy_file (str_in, str_out_npt[nstr], first, file_prod, 
		&written_to_npt[nstr]);
	  fclose (str_in);
	}

      /* Sampled engineering (Quicklook) */
      strcpy (filename_ext, DATA_IN_DIR);
      strcat (filename_ext, satpassname);
      strcat (filename_ext, "qlk");
      if ((str_in = fopen (filename_ext, "r")) == NULL)
	{
	  printf ("Could not open file %s\n", filename_ext);
	}
      else
	{
	  copy_file (str_in, str_out_qlk[nstr], first, file_prod, 
		&written_to_qlk[nstr]);
	  fclose (str_in);
	}

      /* Full rate */
      strcpy (filename_ext, DATA_IN_DIR);
      strcat (filename_ext, satpassname);
      strcat (filename_ext, "frd");
      if ((str_in = fopen (filename_ext, "r")) == NULL)
	{
	  printf ("Could not open file %s\n", filename_ext);
	}
      else
	{
	  copy_file (str_in, str_out_frd[nstr], first, file_prod,
		&written_to_frd[nstr]);
	  fclose (str_in);
	}
      first = 0;
    }

  /* finish and close the last output files */
  for (ii = 0; ii < nf; ii++)
    {
      if (written_to_npt[ii])
	write_h9 (str_out_npt[ii]);

      if (written_to_qlk[ii])
	write_h9 (str_out_qlk[ii]);

      if (written_to_frd[ii])
	write_h9 (str_out_frd[ii]);

      fclose (str_out_npt[ii]);
      fclose (str_out_qlk[ii]);
      fclose (str_out_frd[ii]);
    }
  printf ("Merged %d passes!\n", nitems);
}

/*-------------------------------------------------------------------------
 * **
 * **      compar_all - Routine used by qsort to sort all files names in
 * **                    ascending date/time order.
 * **
 * **-----------------------------------------------------------------------*/
int
compar_all (in1, in2)
     char *in1;
     char *in2;
{
  return (strncmp (in1, in2, 50));
}

/*-------------------------------------------------------------------------
 * **
 * **      get_time - get the current data and time.
 * **
 * **-----------------------------------------------------------------------*/
double
get_time (int *year, int *mon, int *day, int *hour)
{
  struct tm gmt;
  long lclock;
  char mon_name[4];

  time (&lclock);
  gmt = *gmtime (&lclock);
  doytogr (gmt.tm_year, gmt.tm_yday + 1, mon, day, mon_name); 
  *year= gmt.tm_year+ 1900;
  *hour= gmt.tm_hour;
}

/*-------------------------------------------------------------------------
 * **
 * **      get_sat_ids_short - Given a laser target's ilrs ID, norad ID, sic,
 * **                    or name, get the other 3 IDs.
 * **
 * **      NOTE: this differs from get_sat_ids in that the target name
 * **		 does not have trailing spaces.
 * **
 * **-----------------------------------------------------------------------*/
void
get_sat_ids_short (int mode, int *cospar_id, int *norad_id, int *sic, 
	     char *target, int *delta)
{
  FILE *sat_id_in;
  char *sat_id_file = "/data/lib/targets.dat";
  char str[256], ttarget[11];
  int status, tcospar_id, tsic, tnorad_id = 0;
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
      printf ("Could not open file %s\n", sat_id_file);
      exit (1);
    }
  while ((status = fgets (str, 256, sat_id_in)) != NULL)
    {
      sscanf (str, "%s %*s %d %d %d %d", ttarget, &tsic, &tcospar_id,
	      &tnorad_id, delta);

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
	  /*l = strlen (target);
	     if (l < 10)
	     target[l] = ' ';
	     target[10] = '\0'; */
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
  fclose (sat_id_in);
}

/*-------------------------------------------------------------------------
 * **
 * **      copy_file - copies an input crd file to an output file,
 * **                    dropping or modifying headers as needed.
 * **
 * **-----------------------------------------------------------------------*/
void
copy_file (FILE * str_in, FILE * str_out, int first, 
	struct file_date file_prod, int *written) {
  char str[256];
  int status;

  while ((status = fgets (str, 256, str_in)) != NULL)
    {

      /* format/production info */
      if (strncmp (str, "h1", 2) == 0 || strncmp (str, "H1", 2) == 0)
	{
	  /* drop it after first pass */
	  if (first)
	    {
	      read_h1 (str, &h1);
	      h1.prod_year= file_prod.year;
	      h1.prod_mon= file_prod.mon;
	      h1.prod_day= file_prod.day;
	      h1.prod_hour= file_prod.hour;
	      write_h1 (str_out, h1);
	      *written = 1;
	    }
	}
      /* Station header */
      else if (strncmp (str, "h2", 2) == 0 || strncmp (str, "H2", 2) == 0)
	{
	  /* drop it after first pass */
	  if (first)
	    {
	      fputs (str, str_out);
	      *written = 1;
	    }
	}
      /* Satellite header */
      else if (strncmp (str, "h3", 2) == 0 || strncmp (str, "H3", 2) == 0)
	{
	  /* drop it after first pass */
	  if (first)
	    {
	      fputs (str, str_out);
	      *written = 1;
	    }
	}
      /* EOF header */
      else if (strncmp (str, "h9", 2) == 0 || strncmp (str, "H9", 2) == 0)
	{
	  /* drop it */
	}
      else
	{
	  fputs (str, str_out);
	  *written = 1;
	}
    }
}

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../include/crd.h"
struct rh1 h1;
struct rh2 h2;
struct rh3 h3;
struct rh4 h4;
struct rd00 d00;

/*-------------------------------------------------------------------------
 * Program: crd_split
 *
 * Purpose:
 * Split CRD data file into separate file for each pass and data type.
 * Output file names are generated automatically from the data file 
 * based on the protocols in section 5.1.1.1 (Station naming convention) in
 * the CRD format document.
 *
 * Calling sequence:
 *   crd_split crd_filename
 *
 * Author: Randall Ricklefs / Univ of Texas / Center for Space Research
 *
 * History:
 *   July 31, 2007 - Initial version
 *   Dec 5, 2008   - To support multiple segments files from the same pass, it is
 *                   necessary to include minutes on the file name. rlr. v1.00a
 *
**-----------------------------------------------------------------------*/
main (argc, argv)
     int argc;
     char *argv[];
{
  int h9_read= 1;	/* pretend we've already read h9 (eof) */
  int i;
  char comment [256];
  char out_name[256]= {};
  char str[256];
  char target_name[11];
  static char type_str[3][4]= {"frd\0","npt\0","qlk\0"};
  FILE *str_in, *str_out;

  if ((str_in = fopen (argv[1], "r")) == NULL)
    {
      printf ("Could not open file %s\n", argv[1]);
      printf("Usage: crd_split crd_filename\n");
      exit (1);
    }

  /* Copy and reformat data */
  while (fgets (str, 256, str_in) != NULL)
    {
      if (strncmp (str, "00", 2) == 0)
	{
          /* must read "00" so that we don't write to unopened file */
	  /* May need to expand this to multiple lines later. */
	  read_00 (str, &d00);
	}
      else if (strncmp (str, "h1", 2) == 0 ||
               strncmp (str, "H1", 2) == 0)
	{
	  read_h1 (str, &h1);
	}
      else if (strncmp (str, "h2", 2) == 0 ||
               strncmp (str, "H2", 2) == 0)
	{
	  read_h2 (str, &h2);
	}
      else if (strncmp (str, "h3", 2) == 0 ||
               strncmp (str, "H3", 2) == 0)
	{
	  read_h3 (str, &h3);
	}
      else if (strncmp (str, "h4", 2) == 0 ||
               strncmp (str, "H4", 2) == 0)
	{
	  read_h4 (str, &h4);

          /* Write eof to the open file and close it */
	  if (!h9_read)
            {
              write_h9 (str_out);
            }
          else
            {
	      h9_read= 0;
            }
	  if (strlen (out_name) > 0) fclose (str_out);

          /* open the file and write some headers */
          /* ssss_satellite_crd_yymmdd_hh_r.typ */
          strncpy(target_name,h3.target_name,10);
          for (i=9;i>=0;i--) 
            {
              if (isalnum(target_name[i]))
                {
                  target_name[i+1]= '\0';
                  break;
                }
            }
	  sprintf(out_name,"%04d_%-s_crd_%02d%02d%02d_%02d%02d_%1d.%-3s",
	    h2.cdp_pad_id, target_name, 
            h4.start_year%100, h4.start_mon, h4.start_day, h4.start_hour, 
            h4.start_min, h4.data_release, type_str[h4.data_type]);
	  printf("new file name:[%s]\n",out_name);
	  /* In case there are several pass_segment files within an hour,
	   * do an open for append rather than write */
          if ((str_out = fopen (out_name, "a")) == NULL)
            {
               printf ("Could not open file %s\n", out_name);
               exit (1);
            }
	  write_h1 (str_out, h1);
	  write_h2 (str_out, h2);
	  write_h3 (str_out, h3);
	  write_h4 (str_out, h4);
	}
      else if (strncmp (str, "h8", 2) == 0 ||
               strncmp (str, "H8", 2) == 0)
	{
	  read_h8 (str);
	  write_h8 (str_out);
	}
      else if (strncmp (str, "h9", 2) == 0 ||
               strncmp (str, "H9", 2) == 0)
	{
	  read_h9 (str);
	  write_h9 (str_out);
	  h9_read= 1;
	}
      else
	{
	  fputs (str, str_out);
	}

    }
}

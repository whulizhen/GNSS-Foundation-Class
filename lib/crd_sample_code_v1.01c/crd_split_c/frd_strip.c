#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../include/crd.h"
struct rd10 d10;
struct rd30 d30;

/*-------------------------------------------------------------------------
 * Program: frd_strip
 *
 * Purpose:
 * Remove station-specific records from full rate data, and compress some
 *      other records. This is the final step to create the fullrate file
 *	to send to the data center.
 *
 * Calling sequence:
 *   frd_strip input_crd_filename output_crd_filename
 *
 * Author: Randall Ricklefs / Univ of Texas / Center for Space Research
 *
 * History:
 *   July 31, 2007 - Initial version
 *
**-----------------------------------------------------------------------*/
main (argc, argv)
     int argc;
     char *argv[];
{
  char str[256];
  FILE *str_in, *str_out;

  if ((str_in = fopen (argv[1], "r")) == NULL)
    {
      printf ("Could not open file %s\n", argv[1]);
      exit (1);
    }
  if ((str_out = fopen (argv[2], "w")) == NULL)
    {
      printf ("Could not open file %s\n", argv[2]);
      exit (1);
    }

  /* Copy and reformat data */
  while (fgets (str, 256, str_in) != NULL)
    {
      if (strncmp (str, "9", 1) == 0)
	{
	  /* drop it */
	}
      else if (strncmp (str, "10", 2) == 0)
	{
	  /* squeeze out some spaces */
	  read_10 (str, &d10);
	  write_10 (str_out, d10);
	}
      else if (strncmp (str, "30", 2) == 0)
	{
	  /* drop the extra angle records */
	  read_30 (str, &d30);
	  if (d30.angle_origin_ind >= 0)
	    fputs (str, str_out);
	}
      else
	{
	  fputs (str, str_out);
	}

    }
}

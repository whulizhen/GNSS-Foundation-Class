#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "../include/crd.h"

struct rh1 h1;
struct rh2 h2;
struct rh3 h3;
struct rh4 h4;
struct rc0 c0;
struct rc1 c1;
struct rc2 c2;
struct rc3 c3;
struct rc4 c4;
struct rd10 d10;
struct rd11 d11;
struct rd12 d12;
struct rd20 d20;
struct rd21 d21;
struct rd30 d30;
struct rd40 d40;
struct rd50 d50;
struct rd60 d60;
struct rd00 d00;

/*  Test the consistency of the CRD c read/write routines. */
main (argc, argv)
     int argc;
     char *argv[];
{
  char str[512];
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
  printf("in: %s\nout: %s\n", argv[1], argv[2]);

  /* Copy and reformat data */
  while (fgets (str, 512, str_in) != NULL)
    {
      if (strncmp (str, "H1", 2) == 0 ||
          strncmp (str, "h1", 2) == 0)
        {
          read_h1 (str, &h1);
	  write_h1 (str_out, h1);
        }
      else if (strncmp (str, "H2", 2) == 0 ||
               strncmp (str, "h2", 2) == 0)
        {
          read_h2 (str, &h2);
	  write_h2 (str_out, h2);
        }
      else if (strncmp (str, "H3", 2) == 0 ||
               strncmp (str, "h3", 2) == 0)
        {
          read_h3 (str, &h3);
	  write_h3 (str_out, h3);
        }
      else if (strncmp (str, "H4", 2) == 0 ||
               strncmp (str, "h4", 2) == 0)
        {
          read_h4 (str, &h4);
	  write_h4 (str_out, h4);
        }
      else if (strncmp (str, "H8", 2) == 0 ||
               strncmp (str, "h8", 2) == 0)
        {
          read_h8 (str);
	  write_h8 (str_out);
        }
      else if (strncmp (str, "H9", 2) == 0 ||
               strncmp (str, "h9", 2) == 0)
        {
          read_h9 (str);
	  write_h9 (str_out);
        }
      else if (strncmp (str, "C0", 2) == 0 ||
               strncmp (str, "c0", 2) == 0)
        {
          read_c0 (str, &c0);
          write_c0 (str_out, c0);
        }
      else if (strncmp (str, "C1", 2) == 0 ||
               strncmp (str, "c1", 2) == 0)
        {
          read_c1 (str, &c1);
          write_c1 (str_out, c1);
        }
      else if (strncmp (str, "C2", 2) == 0 ||
               strncmp (str, "c2", 2) == 0)
        {
          read_c2 (str, &c2);
          write_c2 (str_out, c2);
        }
      else if (strncmp (str, "C3", 2) == 0 ||
               strncmp (str, "c3", 2) == 0)
        {
          read_c3 (str, &c3);
          write_c3 (str_out, c3);
        }
      else if (strncmp (str, "C4", 2) == 0 ||
               strncmp (str, "c4", 2) == 0)
        {
          read_c4 (str, &c4);
          write_c4 (str_out, c4);
        }
      else if (strncmp (str, "10", 2) == 0)
        {
          read_10 (str, &d10);
          write_10 (str_out, d10);
        }
      else if (strncmp (str, "11", 2) == 0)
        {
          read_11 (str, &d11);
          write_11 (str_out, d11);
        }
      else if (strncmp (str, "12", 2) == 0)
        {
          read_12 (str, &d12);
          write_12 (str_out, d12);
        }
      else if (strncmp (str, "20", 2) == 0)
        {
          read_20 (str, &d20);
          write_20 (str_out, d20);
        }
      else if (strncmp (str, "21", 2) == 0)
        {
          read_21 (str, &d21);
          write_21 (str_out, d21);
        }
      else if (strncmp (str, "30", 2) == 0)
        {
          read_30 (str, &d30);
          write_30 (str_out, d30);
        }
      else if (strncmp (str, "40", 2) == 0)
        {
          read_40 (str, &d40);
          write_40 (str_out, d40);
        }
      else if (strncmp (str, "50", 2) == 0)
        {
          read_50 (str, &d50);
          write_50 (str_out, d50);
        }
      else if (strncmp (str, "60", 2) == 0)
        {
          read_60 (str, &d60);
          write_60 (str_out, d60);
        }
      else if (strncmp (str, "00", 2) == 0)
        {
          read_00 (str, &d00);
          write_00 (str_out, d00);
        }
      /*fputs (str, str_out);*/

    }
}

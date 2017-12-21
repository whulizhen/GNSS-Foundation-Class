/*
C  CALCULATE THE GREGORIAN DATE FROM JULIAN DATE.
C
C  REVISIONS:
C	09/20/89 - Converted FORTRAN version to C. rlr.
C	06/09/93 - Eliminate round off that caused min=60 or sec=60. rlr.
C
C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
*/
jdtogr(jdint,jdf,year,month,day,hour,minute,second)

double	jdint,
	jdf,
	*second;
int	*year,
	*month,
	*day,
	*hour,
	*minute;
{
    double	jda1900,
		jdfc;
    int		tday;

    jdfc = jdf + (jdint- (long)(jdint)) + 0.5;
    jda1900 = (long)(jdint) - 2415079.0;
    while (jdfc >= 1.0) {
        jdfc  = jdfc - 1.0;
        jda1900= jda1900 + 1.0;
    }

    *hour   = jdfc*24.0+ 1.e-10;
    *minute = jdfc*1440.0 - *hour*60.0+ 1.e-8;
    *second = (double)((jdfc - *hour/24.0 - *minute/1440.0)*86400.0);
    *year   = ((4.0*jda1900)-1.0)/1461.0;
 
    tday   = ((4.0*jda1900)+3.0- (*year*1461.0))/4.0;
    *month  = ((5.0*tday)-3.0)/153.0;
    *day    = ((5.0*tday)+2.0 - (153.0*(*month)))/5.0;
 
    if (*month >= 10) {
        *month = *month - 9;
        *year  = *year + 1;
    } else {
	*month = *month + 3;
    }
 
}

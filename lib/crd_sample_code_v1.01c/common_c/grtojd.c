/*
C FOR A GIVEN GREGORIAN DATE & UNIVERSAL TIME COORDINATED (UTC) THIS
C ROUTINE CALCULATES THE EQUIVALENT JULIAN-DATE (RJD).
C
C INPUT:
C   YEAR   - (INTEGER) GREGORIAN YEAR LEAST SIG. 2 DIGITS
C   MONTH  - (INTEGER) GREGORIAN MONTH
C   DAY    - (INTEGER) GREGORIAN DAY
C   HOUR   - (INTEGER) UTC HOUR
C   MINUTE - (INTEGER) UTC MINUTE
C   SECOND - (D.P. REAL) UTC SECOND
C
C OUTPUT: 
C   JDINT  - DOUBLE PRECISION JULIAN DATE (WHOLE DAY)
C   JDF    - (FRACTIONAL DAY)
C   YEAR, MONTH, DAY, HOUR, MINUTE, SECOND - UNCHANGED
C
C  Revisions:
C	09/21/89 - Convert from FORTRAN to C. rlr.
C
*/
grtojd(year,month,day,hour,minute,second,jdint,jdf)

int	year,
	month,
	day,
	hour,
	minute;
double	*jdint,
	*jdf,
	second;

{ 
/* calculate # days since noon feburary 29, 1900 (julian date=2415078.0) */
 
    if (month <= 2) {
 
	*jdint = (long)(1461.0 * (year-1)/4.0)
		+ (long)((153.0 * (month+9) + 2.0)/5.0) + day;

    } else { 

	*jdint = (long)(1461.0 * year/4.0) 
		+ (long)((153.0 * (month-3) + 2.0)/5.0) + day;

    }
 
/* add fractional day and the jd for 2/29/1900 */
 
    *jdf = (hour + (minute + second/60.0)/60.0)/24.0+ 0.5;
    *jdint = *jdint + 2415078.0;
    while (*jdf >= 1.0) {
        *jdf= *jdf- 1.0;
        *jdint= *jdint+ 1.0;
    }
}

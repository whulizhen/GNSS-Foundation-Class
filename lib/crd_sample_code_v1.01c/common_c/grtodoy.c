grtodoy(year,month,day,doy)
int	year,
	month,
	day,
	*doy;
{
    long	jd_in,
		jd_jan0;

/* calculate # days since noon feburary 29, 1900 (julian date=2415078.0) */
 
    if (month <= 2) {
 
	jd_in = (long)(1461.0 * (year-1)/4.0)
		+ (long)((153.0 * (month+9) + 2.0)/5.0) + day;

    } else { 

	jd_in = (long)(1461.0 * year/4.0) 
		+ (long)((153.0 * (month-3) + 2.0)/5.0) + day;

    }

/* what is jd of jan 0? */
    jd_jan0 = (long)(1461.0 * (year-1)/4.0)+ 306.0;
    
/*  now doy is easy.  */
    *doy= jd_in- jd_jan0;
}
 

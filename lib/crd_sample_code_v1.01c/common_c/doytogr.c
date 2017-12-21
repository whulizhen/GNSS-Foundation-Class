doytogr(year,doy,month,day,mon_name)
    int	year,
	doy,
	*month,
	*day;
    char	mon_name[4];
{
    static char	*month_names[12] = {"Jan","Feb","Mar","Apr","May","Jun",
				    "Jul","Aug","Sep","Oct","Nov","Dec"};
    double	jda1900;
    int		i,
		tday,
		tyear;

    jda1900 = (long)(1461.0*(year-1)/4.0)+ 306.0+ doy;

    tyear   = ((4.0*jda1900)-1.0)/1461.0;
    tday    = ((4.0*jda1900)+3.0- (tyear*1461.0))/4.0;
    *month  = ((5.0*tday)-3.0)/153.0;
    *day    = ((5.0*tday)+2.0 - (153.0*(*month)))/5.0;
 
    if (*month >= 10) {
        *month = *month - 9;
    } else {
	*month = *month + 3;
    }

    for (i=0; (mon_name[i]= month_names[*month-1][i]) != 0; i++) {}
 
}

#include <ctype.h>
#include <math.h>
#include <stdio.h>
/*------------------------- test driver --------------------------------- 
main()
{
    double dbl;
    long lng;
    int i, j;
    char str[256];

    for (;;) {
	printf("Enter test value:");
	gets(str);
	i= getfield(str,0,20,&dbl);
	j= getifield(str,0,20,&lng);
	printf("Float conversion: %d, %g\n  Int conversion: %d, %d\n",
	    i,dbl,j,lng);
    }
}
------------------------------------------------------------------------*/
/*  getfield --  converts a string of characters to a double float point #.
**		 Routine checks for non-numeric characters and empty
**		 field, returning 0 for satifactory completion,
**		 -1 for syntax error and -2 for empty field.
*/
int
getfield(buffer, pos, fieldsize, value)
	char *buffer;		/*  string to parse		*/
	int pos;		/*  position in string to begin */
	int fieldsize;		/*  length of field to parse	*/
	double *value;		/*  value found in field	*/
{
	static char SccsId[] = "@(#)getfield.c	1.1\t12/06/90";
	int c, i, sawDecPoint, decDigits, sawMinus, numblanks;
	int sawExponent, sawExMinus, sawNumber;
	double intval, fracval, fracdiv, expval;
	
	fracdiv = 1;
	expval = 0.;
	intval = 0;
	fracval = 0;
	sawDecPoint = 0;
	sawExponent = 0;
	sawExMinus = 0;
	sawNumber= 0;
	sawMinus = 0;
	numblanks = 0;
	for (i=pos; i<fieldsize+pos; ++i) {
		c = buffer[i];
		if (c >= '0' && c <= '9') {
			if (sawExponent) {
				expval = expval * 10 + c - '0';
			} else if (sawDecPoint) {
				fracdiv *= 10;
				fracval = fracval * 10 + c - '0';
			} else {
				intval = intval * 10 + c - '0';
			}
			sawNumber= 1;
		} else if (c == '.') {
			if (sawDecPoint || sawExponent) return -1;
			sawDecPoint = 1;
		} else if (c == '-') {
			if (sawExponent) {
				if (sawExMinus) return -1;
				sawExMinus= 1;
			} else {
				if (sawMinus) return -1;
				sawMinus = 1;
			}
		} else if (c == 'e' || c == 'E' || c == 'd' || c == 'D') {
			if (sawExponent) return -1;
			if (!sawNumber) return -1;
			sawExponent = 1;
		} else if (isspace(c)) {
			numblanks ++ ;
		} else if (c == NULL) {
			numblanks += fieldsize- (i-pos);
			break;
		} else return -1;
	}
	if (sawDecPoint) {
		*value = (double)intval + (double)fracval / (double)fracdiv;
	} else {
		*value = (double)intval;
	}
	if (sawExMinus) expval= -expval;
	if (sawExponent) *value= *value * pow(10.,expval);
	if (sawMinus) *value = - *value;
	if (numblanks == fieldsize) return -2;	/* field all blanks */
	else return 0;
}

/*  getifield -- converts a string of characters to a long integer.
**		 Routine checks for non-numeric characters and empty
**		 field, returning 0 for satifactory completion,
**		 -1 for syntax error and -2 for empty field.
**
**  Note: routine DOES NOT CHECK FOR EVERFLOW (>2147483647= 2**31)
*/
int
getifield(buffer, pos, fieldsize, value)
	char *buffer;		/*  string to parse		*/
	int pos;		/*  position in string to begin */
	int fieldsize;		/*  length of field to parse	*/
	long *value;		/*  value found in field	*/
{
	int c, i, decDigits, sawMinus, numblanks;
	long intval;
	
	intval = 0;
	sawMinus = 0;
	numblanks = 0;
	for (i=pos; i<fieldsize+pos; ++i) {
		c = buffer[i];
		if (c >= '0' && c <= '9') {
			intval = intval * 10 + c - '0';
		} else if (c == '.') {
			return -1;
		} else if (c == '-') {
			if (sawMinus) return -1;
			sawMinus = 1;
		} else if (isspace(c)) {
			numblanks ++ ;
		} else if (c == NULL) {
			numblanks += fieldsize- (i-pos);
			break;
		} else return -1;
	}
	if (sawMinus) *value = - intval;
	else *value = intval;
	if (numblanks == fieldsize) return -2;	/* field all blanks */
	else return 0;
}

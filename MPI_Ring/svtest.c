/** @file 
 * \brief svtest.c - a very simple input validator  
 * \author Philip Papadopoulos
 * \email ppapadopoulos@ucsd.edu
 * \date 14 Jan 2018
*/
#include <stdio.h>
#include <stdlib.h>
#include "svalidate.h" 

/**
 * \fn int main(int argc, char *argv[])
 * \param argc - number of command-line arguments 
 * \param argv - string representations of arguments 
 */
int main (int argc, char * argv[])
{
	int i;
	for (i = 1; i < argc; i++)
	{
		char * str = argv[i];
		printf("=== input arg %d is '%s'\n",i,str);
		printf("    trimmed string is '%s'\n",str=trim(str));
		if ( isInteger(str) )
		{
			int ival;
			sscanf(str,"%d",&ival);
			printf("    int value is %d\n", ival);
		}
		if ( isFloat(str) )
		{
			float fval;
			sscanf(str,"%f",&fval);
			printf("    float value is %f\n", fval);
		}
	}
	exit(0);
}
// vim: ts=4

/** @file 
 * \brief header file for simple string validation functions 
 * \author Mingcheng Zhu
 * \email zhumc11@gmail.com
 * \date 19 Jan 2018
*/
#define FALSE 0
#define TRUE 1
/* Function prototypes  */
extern char * trim (char * input);
extern char * trimCopy( char *input);
extern int isInteger(char * input);
extern int isFloat(char * input);
extern int isDouble(char * input);

/*=====================================================================*
 *                                                                      *
 *                 Make index for multiple CMP Files.                   *
 *                                                                      *
 *======================================================================*
 *   Summary:                                                           *
 *      The program will scan CMP files and creates index files.        *
 *      Format of binary index files:                                   *
 *      typedef struct{                                                 *
 *         int   line;                                                  *
 *         int   cdp;                                                   *
 *         int   sx;                                                    *
 *         int   sy;                                                    *
 *         int   gx;                                                    *
 *         int   gy;                                                    *
 *         int   offset;                                                *
 *      }cmpFileIndex;                                                  *
 *======================================================================*
 *       Author: Feng Bo                                                *
 *       Date  : 2022/01/25.                                            *
 *======================================================================*/

//include
#include "fbopw3d.h"

int main( int argc , char *argv[] )
{
    int mystat  = 0;
    int verbose = 1;
    char    parFile[FILE_NAME_MAX_LENGTH]="";
    strcpy(parFile, argv[1]);

    int     nfile = -1;
	char	*Suffix_in="sgy";
	char	*Suffix_out="cmpIdx";
    char    **CMPFileNameList   = alloc2char(FILE_NAME_MAX_LENGTH, FILE_NUMBER_MAX);
    char    **CMPIndexNameList  = alloc2char(FILE_NAME_MAX_LENGTH, FILE_NUMBER_MAX);
    readpar_index(parFile, &nfile, CMPFileNameList, CMPIndexNameList, verbose);

mainProgStop:
    fprintf(stdout, "mystat = %d\n", mystat);

    free2char(CMPFileNameList);
    free2char(CMPIndexNameList);

	/*
	char	filename_in[FILE_NAME_MAX_LENGTH]="/home/data/sgyfile/LJ2017_cmp_Line2001.sgy";
	char	path_out[FILE_NAME_MAX_LENGTH]="/home/fb";
	char	filename_out[FILE_NAME_MAX_LENGTH]="";
	char	*Suffix_in="sgy";
	char	*Suffix_out="cmpIdx";
	replacePathSuffix(filename_in, Suffix_in, path_out, Suffix_out, filename_out, 1);
	printf("Index filename is %s\n", filename_out);
	 */

    return mystat;
}


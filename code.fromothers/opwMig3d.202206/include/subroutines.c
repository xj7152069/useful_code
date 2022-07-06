
#include "fbCommon.h"

void str_replace(char *str_src, int n, char * str_copy)
{
    int lenofstr;
    char *tmp;
    lenofstr = strlen(str_copy);
    //string to copy is shorter than string to find
    if(lenofstr < n)
    {
        tmp = str_src+n;
        while(*tmp)
        {
            *(tmp-(n-lenofstr)) = *tmp; //n-lenofstr, length of moving
            tmp++;
        }
        *(tmp-(n-lenofstr)) = *tmp; //move '\0'
    }
    else if(lenofstr > n) //string to copy longer than string to find
    {
        tmp = str_src;
        while(*tmp) tmp++;
        while( tmp>=(str_src+n) )
        {
            *(tmp+(lenofstr-n)) = *tmp;
            tmp--;
        }
    }
    strncpy(str_src, str_copy, lenofstr);
}

void    suffix_replace(char *str_source, char *str_find, char *str_to_copy)
{
    char *p;
    p = strstr(str_source, str_find);
    while(p){

        str_replace(p, strlen(str_find), str_to_copy);
        p = p + strlen(str_to_copy);
        p = strstr(p, str_find);
    }
}

int readpar(
        char    *parFile,
        int *fileNumber,
        char    **CMPFileNameList,
        char    **CMPIndexNameList,
        int verbose
        )
{
    char    tmp[FILE_NAME_MAX_LENGTH]="";
    char    Suffix_in[FILE_NAME_MAX_LENGTH]="";
    char    Suffix_out[FILE_NAME_MAX_LENGTH]="";

    FILE    *fp;
    if((fp=fopen(parFile,"r"))==NULL)
    {
        printf("Can not open the paramter file: %s to read!\n", parFile);
        return 1;
    }

    fscanf(fp,"%s%s", tmp, Suffix_in);  // line 01.
    fscanf(fp,"%s%s", tmp, Suffix_out); // line 02.
    fscanf(fp,"%s%d", tmp, fileNumber); // line 03.
    int nfile   = *fileNumber;
    for (int i = 0 ; i < nfile; ++i)
    {
        fscanf(fp,"%s", CMPFileNameList[i]);

        strcpy(CMPIndexNameList[i], CMPFileNameList[i]);

        suffix_replace(CMPIndexNameList[i], Suffix_in, Suffix_out);
    }

    if( verbose == 1 )
    {
        fprintf(stdout, "the input CMP-file suffix is: %s\n", Suffix_in);
        fprintf(stdout, "the output Index-file suffix is: %s\n", Suffix_out);
        fprintf(stdout, "number of CMP files is %d\n", nfile);

        fprintf(stdout, "Input Filename list:\n");
        for (int i = 0 ; i < nfile; ++i)
            fprintf(stdout, "the input CMP filename is: %s\n", CMPFileNameList[i]);
        fprintf(stdout, "\n");

        fprintf(stdout, "\nOutput Filename list:\n");
        for (int i = 0 ; i < nfile; ++i)
            fprintf(stdout, "the output Index filename is: %s\n", CMPIndexNameList[i]);
        fprintf(stdout, "\n");
    }

    fclose(fp);

    return 0;
}

int readPar_OPWD(
        int     myid,
        char    *parFile,
        int     *fileNumber,
        char    **CMPFileNameList,
        char    **CMPIndexNameList,
        char    *output_dir,
        long    *fileszieMaxByte,
        int     *cmpfileOffsetBeg,
        int     *conv,
        int     *endian,
        int     *ns,
        int     *dt_us,
        float   *Coor_S_R_Scale,
        int     *ntrMin_PWD,
        int     *ntau,
        int     *itaus,
        int     *itaue,
        int     *nphr,
        float   *phrmin,
        float   *phrmax,
        float   *dphr,
        int     *nline_loc,
        int     *firstCMPLineNo_PWD,
        int     *lastCMPLineNo_PWD,
        int     verbose
        )
{
    char    tmp[FILE_NAME_MAX_LENGTH]="";
    char    Suffix_in[FILE_NAME_MAX_LENGTH]="";
    char    Suffix_out[FILE_NAME_MAX_LENGTH]="";

    FILE    *fp;
    if((fp=fopen(parFile,"r"))==NULL)
    {
        printf("Can not open the paramter file: %s to read!\n", parFile);
        return 1;
    }

    fscanf(fp,"%s%s", tmp, Suffix_in);  // line 01.
    fscanf(fp,"%s%s", tmp, Suffix_out); // line 02.
    fscanf(fp,"%s%d", tmp, fileNumber); // line 03.
    int nfile   = *fileNumber;
    for (int i = 0 ; i < nfile; ++i)
    {
        fscanf(fp,"%s", CMPFileNameList[i]);

        strcpy(CMPIndexNameList[i], CMPFileNameList[i]);

        suffix_replace(CMPIndexNameList[i], Suffix_in, Suffix_out);
    }


    if( verbose == 1 && 0 == myid )
    {
        fprintf(stdout, "the input CMP-file suffix is: %s\n", Suffix_in);
        fprintf(stdout, "the output Index-file suffix is: %s\n", Suffix_out);
        fprintf(stdout, "number of CMP files is %d\n", nfile);

        fprintf(stdout, "Input Filename list:\n");
        for (int i = 0 ; i < nfile; ++i)
            fprintf(stdout, "the input CMP filename is: %s\n", CMPFileNameList[i]);
        fprintf(stdout, "\n");

        fprintf(stdout, "\nOutput Filename list:\n");
        for (int i = 0 ; i < nfile; ++i)
            fprintf(stdout, "the output Index filename is: %s\n", CMPIndexNameList[i]);
        fprintf(stdout, "\n");
    }

    char    par_dir[FILE_NAME_MAX_LENGTH]="";
    long    par_fileszieMaxByte;
    int     par_cmpfileOffsetBeg    = 3600; // =0, SU-format; =3600, SEGY-format.
    int     par_conv    = 1;                // = 0 ; assume data is in native format.
    int     par_endian  = 0;                // set =0 for little-endian machines(PC's,DEC,etc.).
    int     par_ns;
    int     par_dt_us;                      // unit of "dt_us" is microsecond (1E-6*second).
    float   par_Coor_S_R_Scale      = 1.0;  // set value of global variable.
    fscanf(fp,"%s%s", tmp, par_dir);
    fscanf(fp,"%s%ld",tmp, &par_fileszieMaxByte);
    fscanf(fp,"%s%d", tmp, &par_cmpfileOffsetBeg);
    fscanf(fp,"%s%d", tmp, &par_conv);
    fscanf(fp,"%s%d", tmp, &par_endian);
    fscanf(fp,"%s%d", tmp, &par_ns);
    fscanf(fp,"%s%d", tmp, &par_dt_us);
    fscanf(fp,"%s%f", tmp, &par_Coor_S_R_Scale);

    if( verbose == 1 && 0 == myid )
    {
        fprintf(stderr,"%s\n", par_dir);
        fprintf(stderr,"fileszieMaxByte=%ld\n",par_fileszieMaxByte);
        fprintf(stderr,"cmpfileOffsetBeg=%d\n", par_cmpfileOffsetBeg);
        fprintf(stderr,"conv=%d\n", par_conv);
        fprintf(stderr,"endian=%d\n", par_endian);
        fprintf(stderr,"ns=%d\n", par_ns);
        fprintf(stderr,"dt_us=%d\n", par_dt_us);
        fprintf(stderr,"Coor_S_R_Scale=%f\n", par_Coor_S_R_Scale);
    }

    strcpy(output_dir, par_dir);
    *fileszieMaxByte    = par_fileszieMaxByte;
    *cmpfileOffsetBeg   = par_cmpfileOffsetBeg;
    *conv   = par_conv;
    *endian = par_endian;
    *ns     = par_ns;
    *dt_us  = par_dt_us;
    *Coor_S_R_Scale = par_Coor_S_R_Scale;

    // parameters for tau-phr transform.
    int par_ntrMin_PWD, par_ntau, par_nphr;
    int par_itaus   = 0;
    int par_itaue   = par_ns;
    float   par_dphr, par_phrmin, par_phrmax;
    fscanf(fp,"%s%d", tmp, &par_ntrMin_PWD);
    fscanf(fp,"%s%d", tmp, &par_ntau);
    fscanf(fp,"%s%d", tmp, &par_nphr);
    fscanf(fp,"%s%f", tmp, &par_phrmax);
    par_phrmin  = 0.0 ;
    par_phrmax  = par_phrmax*1.0E-6 ;
    par_dphr    = (par_phrmax-par_phrmin)/(par_nphr-1);

    if( verbose == 1 && 0 == myid )
    {
        fprintf(stderr,"ntrMin_PWD=%d\n", par_ntrMin_PWD);
        fprintf(stderr,"ntau=%d\n", par_ntau);
        fprintf(stderr,"nphr=%d\n", par_nphr);
        fprintf(stderr,"phrmax_us_m=%f(us/m)\n", par_phrmax*1E6);
        fprintf(stderr,"phrmin_us_m=%f(us/m)\n", par_phrmin*1E6);
        fprintf(stderr,"par_dphr=%f(us/m)\n", par_dphr*1E6);
    }
    *ntrMin_PWD   = par_ntrMin_PWD;
    *ntau    = par_ntau;
    *itaus   = par_itaus;
    *itaue   = par_itaue;
    *nphr    = par_nphr;
    *phrmin  = par_phrmin;
    *phrmax  = par_phrmax;
    *dphr    = par_dphr;

    // parameters for parallel-calculation.
    int par_nline_loc;
    int par_firstCMPLineNo_PWD;
    int par_lastCMPLineNo_PWD;
    fscanf(fp,"%s%d", tmp, &par_nline_loc);
    fscanf(fp,"%s%d", tmp, &par_firstCMPLineNo_PWD);
    fscanf(fp,"%s%d", tmp, &par_lastCMPLineNo_PWD);

    fprintf(stderr,"nline_loc=%d\n", par_nline_loc);
    fprintf(stderr,"firstCMPLineNo_PWD=%d\n", par_firstCMPLineNo_PWD);
    fprintf(stderr,"lastCMPLineNo_PWD=%d\n", par_lastCMPLineNo_PWD);
    *nline_loc  = par_nline_loc;
    *firstCMPLineNo_PWD = par_firstCMPLineNo_PWD;
    *lastCMPLineNo_PWD  = par_lastCMPLineNo_PWD;

    fclose(fp);

    return 0;
}

void	replacePathSuffix(
	char	*filename_in,
	char	*Suffix_in,
	char	*path_out,
	char	*Suffix_out,
	char	*filename_out,
	int	verbose
)
{

	char	slash ='/';	// for Liunx system.
	char	*basename = strrchr(filename_in, slash);
	++ basename;
	char	filename_temp[FILE_NAME_MAX_LENGTH]="";

	if ( 1 == verbose)
		printf("Base name of %s is %s\n", filename_in, basename);
	sprintf(filename_temp, "%s/%s", path_out, basename);
	if ( 1 == verbose)
		printf("temp filename is %s\n", filename_temp);
	suffix_replace(filename_temp, Suffix_in, Suffix_out);
	if ( 1 == verbose)
		printf("Index filename is %s\n", filename_temp);
	strcpy(filename_out, filename_temp);
}

int readpar_index(
        char    *parFile,
        int *fileNumber,
        char    **CMPFileNameList,
        char    **CMPIndexNameList,
        int verbose
        )
{
    char    tmp[FILE_NAME_MAX_LENGTH]="";
    char    Suffix_in[FILE_NAME_MAX_LENGTH]="";
    char    Suffix_out[FILE_NAME_MAX_LENGTH]="";
    char    path_out[FILE_NAME_MAX_LENGTH]="";	// path for output index files.

    FILE    *fp;
    if((fp=fopen(parFile,"r"))==NULL)
    {
        printf("Can not open the paramter file: %s to read!\n", parFile);
        return 1;
    }

    fscanf(fp,"%s%s", tmp, Suffix_in);  // line 01.
    fscanf(fp,"%s%s", tmp, Suffix_out); // line 02.
    fscanf(fp,"%s%s", tmp, path_out);   // line 03.
    fscanf(fp,"%s%d", tmp, fileNumber); // line 04.
    int nfile   = *fileNumber;
    for (int i = 0 ; i < nfile; ++i)
    {
        fscanf(fp,"%s", CMPFileNameList[i]);

        strcpy(CMPIndexNameList[i], CMPFileNameList[i]);

        //suffix_replace(CMPIndexNameList[i], Suffix_in, Suffix_out);
        replacePathSuffix(CMPIndexNameList[i], Suffix_in,
		path_out, Suffix_out, CMPIndexNameList[i], 0);
    }

    if( verbose == 1 )
    {
        fprintf(stdout, "the input CMP-file suffix is: %s\n", Suffix_in);
        fprintf(stdout, "the output Index-file suffix is: %s\n", Suffix_out);
        fprintf(stdout, "number of CMP files is %d\n", nfile);

        fprintf(stdout, "Input Filename list:\n");
        for (int i = 0 ; i < nfile; ++i)
            fprintf(stdout, "the input CMP filename is: %s\n", CMPFileNameList[i]);
        fprintf(stdout, "\n");

        fprintf(stdout, "\nOutput Filename list:\n");
        for (int i = 0 ; i < nfile; ++i)
            fprintf(stdout, "the output Index filename is: %s\n", CMPIndexNameList[i]);
        fprintf(stdout, "\n");
    }

    fclose(fp);

    return 0;
}

int readPar_MakeIndex_OPWD(
        int     myid,
        char    *parFile,
        int     *fileNumber,
        char    **CMPFileNameList,
        char    **CMPIndexNameList,
        long    *fileszieMaxByte,
        int     *cmpfileOffsetBeg,
        int     *conv,
        int     *endian,
        int     *ns,
        int     *dt_us,
        int     OPWD_flag,
        char    *opwdata_dir,
        float   *Coor_S_R_Scale,
        int     *ntrMin_PWD,
        int     *ntau,
        int     *itaus,
        int     *itaue,
        int     *nphr,
        float   *phrmin,
        float   *phrmax,
        float   *dphr,
        int     *nline_loc,
        int     *firstCMPLineNo_PWD,
        int     *lastCMPLineNo_PWD,
        int     verbose
        )
{
    char    tmp[FILE_NAME_MAX_LENGTH]="";
    char    Suffix_in[FILE_NAME_MAX_LENGTH]="";
    char    Suffix_out[FILE_NAME_MAX_LENGTH]="";
    char    path_out[FILE_NAME_MAX_LENGTH]="";	// path for output index files.

    FILE    *fp;
    if((fp=fopen(parFile,"r"))==NULL)
    {
        printf("Can not open the paramter file: %s to read!\n", parFile);
        return 1;
    }

    fscanf(fp,"%s%s", tmp, Suffix_in);  // line 01.
    fscanf(fp,"%s%s", tmp, Suffix_out); // line 02.
    fscanf(fp,"%s%s", tmp, path_out);   // line 03.
    fscanf(fp,"%s%d", tmp, fileNumber); // line 04.
    int nfile   = *fileNumber;
    for (int i = 0 ; i < nfile; ++i)
    {
        fscanf(fp,"%s", CMPFileNameList[i]);

        strcpy(CMPIndexNameList[i], CMPFileNameList[i]);

        replacePathSuffix(CMPIndexNameList[i], Suffix_in,
		path_out, Suffix_out, CMPIndexNameList[i], 0);
    }

    if( verbose == 1 && 0 == myid )
    {
        fprintf(stdout, "the input CMP-file suffix is: %s\n", Suffix_in);
        fprintf(stdout, "the output Index-file suffix is: %s\n", Suffix_out);
        fprintf(stdout, "number of CMP files is %d\n", nfile);

        fprintf(stdout, "Input Filename list:\n");
        for (int i = 0 ; i < nfile; ++i)
            fprintf(stdout, "the input CMP filename is: %s\n", CMPFileNameList[i]);
        fprintf(stdout, "\n");

        fprintf(stdout, "\nOutput Filename list:\n");
        for (int i = 0 ; i < nfile; ++i)
            fprintf(stdout, "the output Index filename is: %s\n", CMPIndexNameList[i]);
        fprintf(stdout, "\n");
    }

    char    par_dir[FILE_NAME_MAX_LENGTH]="";
    long    par_fileszieMaxByte;
    int     par_cmpfileOffsetBeg    = 3600; // =0, SU-format; =3600, SEGY-format.
    int     par_conv    = 1;                // = 0 ; assume data is in native format.
    int     par_endian  = 0;                // set =0 for little-endian machines(PC's,DEC,etc.).
    int     par_ns;
    int     par_dt_us;                      // unit of "dt_us" is microsecond (1E-6*second).
    fscanf(fp,"%s%ld",tmp, &par_fileszieMaxByte);
    fscanf(fp,"%s%d", tmp, &par_cmpfileOffsetBeg);
    fscanf(fp,"%s%d", tmp, &par_conv);
    fscanf(fp,"%s%d", tmp, &par_endian);
    fscanf(fp,"%s%d", tmp, &par_ns);
    fscanf(fp,"%s%d", tmp, &par_dt_us);
    if( verbose == 1 && 0 == myid )
    {
        fprintf(stderr,"fileszieMaxByte=%ld\n",par_fileszieMaxByte);
        fprintf(stderr,"cmpfileOffsetBeg=%d\n", par_cmpfileOffsetBeg);
        fprintf(stderr,"conv=%d\n", par_conv);
        fprintf(stderr,"endian=%d\n", par_endian);
        fprintf(stderr,"ns=%d\n", par_ns);
        fprintf(stderr,"dt_us=%d\n", par_dt_us);
    }
    *fileszieMaxByte    = par_fileszieMaxByte;
    *cmpfileOffsetBeg   = par_cmpfileOffsetBeg;
    *conv   = par_conv;
    *endian = par_endian;
    *ns     = par_ns;
    *dt_us  = par_dt_us;
    if( 0 == OPWD_flag )	// the aboving parameters are used for making CMP index files.
        return  0;

    // the following parameters are used for OPWD.
    float   par_Coor_S_R_Scale      = 1.0;  // set value of global variable.
    fscanf(fp,"%s%s", tmp, par_dir);
    fscanf(fp,"%s%f", tmp, &par_Coor_S_R_Scale);
    if( verbose == 1 && 0 == myid )
    {
        fprintf(stderr,"OPW-data dir is %s\n", par_dir);
        fprintf(stderr,"Coor_S_R_Scale=%f\n", par_Coor_S_R_Scale);
    }
    strcpy(opwdata_dir, par_dir);
    *Coor_S_R_Scale = par_Coor_S_R_Scale;

    // parameters for tau-phr transform.
    int par_ntrMin_PWD, par_ntau, par_nphr;
    float   par_dphr, par_phrmin, par_phrmax;
    fscanf(fp,"%s%d", tmp, &par_ntrMin_PWD);
    fscanf(fp,"%s%d", tmp, &par_ntau);
    fscanf(fp,"%s%d", tmp, &par_nphr);
    fscanf(fp,"%s%f", tmp, &par_phrmax);
    par_phrmin  = 0.0 ;
    par_phrmax  = par_phrmax*1.0E-6 ;
    par_dphr    = (par_phrmax-par_phrmin)/(par_nphr-1);
    int par_itaus   = 0;
    int par_itaue   = par_ntau;

    if( verbose == 1 && 0 == myid )
    {
        fprintf(stderr,"ntrMin_PWD=%d\n", par_ntrMin_PWD);
        fprintf(stderr,"ntau=%d\n",  par_ntau);
        fprintf(stderr,"itaus=%d\n", par_itaus);
        fprintf(stderr,"itaue=%d\n", par_itaue);
        fprintf(stderr,"nphr=%d\n",  par_nphr);
        fprintf(stderr,"phrmax_us_m=%f(us/m)\n", par_phrmax*1E6);
        fprintf(stderr,"phrmin_us_m=%f(us/m)\n", par_phrmin*1E6);
        fprintf(stderr,"par_dphr=%f(us/m)\n",    par_dphr*1E6);
    }
    *ntrMin_PWD   = par_ntrMin_PWD;
    *ntau    = par_ntau;
    *itaus   = par_itaus;
    *itaue   = par_itaue;
    *nphr    = par_nphr;
    *phrmin  = par_phrmin;
    *phrmax  = par_phrmax;
    *dphr    = par_dphr;

    // parameters for parallel-calculation.
    int par_nline_loc;
    int par_firstCMPLineNo_PWD;
    int par_lastCMPLineNo_PWD;
    fscanf(fp,"%s%d", tmp, &par_nline_loc);
    fscanf(fp,"%s%d", tmp, &par_firstCMPLineNo_PWD);
    fscanf(fp,"%s%d", tmp, &par_lastCMPLineNo_PWD);

    fprintf(stderr,"nline_loc=%d\n", par_nline_loc);
    fprintf(stderr,"firstCMPLineNo_PWD=%d\n", par_firstCMPLineNo_PWD);
    fprintf(stderr,"lastCMPLineNo_PWD=%d\n", par_lastCMPLineNo_PWD);
    *nline_loc  = par_nline_loc;
    *firstCMPLineNo_PWD = par_firstCMPLineNo_PWD;
    *lastCMPLineNo_PWD  = par_lastCMPLineNo_PWD;

    fclose(fp);

    return 0;
}


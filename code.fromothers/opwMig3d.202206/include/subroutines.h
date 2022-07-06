// function prototpyes from multiple C-files.
// [BEGIN] in file: cmpFileIndex_nonstd.c //.
int open_cmp_gather(char *cmpGatherFilename);

int close_cmp_gather();

int init_cmp_gather_nonstd(
        // Input parameters.
        int fileOffsetBeg,      // =0, SU-format; =3600, SEGY-format.
        int conv,               // convert data to native format.
                                // = 0 ; assume data is in native format.
        int endian,             // set =0 for little-endian machines(PC's,DEC,etc.).
        long    ntr_cmpfile_max,
        // Output parameters.
        int *ncdp_first,
        int *ncdp_final,
        int *nline_first,
        int *nline_final,
        int *ntr_gather_max,
        long    *ntr_cmpfile_get,
        float   *Offset_X_Max,
        float   *Offset_Y_Max,
        float   *Offset_Max,
        float   *Offset_Min,
        cmpFileIndex *MycmpFileIndex,
        int verbose             // =1, print debug information.
        );

int write_cmpFileIndex(
        char    *CMPIndexName,
        long    ntr_cmpfile,
        cmpFileIndex    *MycmpFileIndex );

int read_cmpFileIndex(
        char    *CMPIndexName,
        long    *ntr_cmpfile,
        cmpFileIndex    *MycmpFileIndex );
// [END] in file: cmpFileIndex_nonstd.c //.

// [BEGIN] in file: subroutines.c //.
void    str_replace(char *str_src, int n, char * str_copy);

void    suffix_replace(char *str_source, char *str_find, char *str_to_copy);

int readpar(
        char    *parFile,
        int *fileNumber,
        char    **CMPFileNameList,
        char    **CMPIndexNameList,
        int verbose
        );
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
        );
void	replacePathSuffix(
	char	*filename_in,
	char	*Suffix_in,
	char	*path_out,
	char	*Suffix_out,
	char	*filename_out,
	int	verbose
	);
int readpar_index(
        char    *parFile,
        int *fileNumber,
        char    **CMPFileNameList,
        char    **CMPIndexNameList,
        int verbose
        );
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
        );
// [END] in file: subroutines.c //.

// [BEGIN] in file: sub_suvelan.c //.
int	weightedSemblance(
        // parameters for velocity scanning.
        float 	fv,		// first velocity.
        float	dv,		// velocity sampling interval.
        int	nv,		    // number of velocities.
        // parameters for the input data.
        float	**gather,	// the input CMP gather.
        int	ntr,		// trace number.
        float	*off,		// offset of each trace.
        int	ns,		// sampling number.
        float	dt,		// sampling interval (second).
        // parameters for the output data.
        float	**semb		// the output weighted-semblance spectrum.
        );
// [END] in file: sub_suvelan.c //.

// [BEGIN] in file: sinc_interpolation.c //.
float	sinc_interpolation ( float *trace , int lt , float dt , float time_s_r );
// [END] in file: sinc_interpolation.c //.

// [BEGIN] in file: tauph_transform.c //.
/*	Function Prototypes.		*/
float	Max_Float_1d(float *array , int n );
float	Min_Float_1d(float *array , int n );
int	Min_Int_1d(int *array , int n );
int	interpolation_taup_domain( float **in, float **out, int ns, int nph, int iphmin, float *vel, float dt, float dph);

int	tauph_transform_3d_to_2d( int ns, int ntau, int nphr, int ntr, int itaus, int itaue,
        float **tcmp, float **pdat, float *vrms, float *offset,
        float dt, float phrmin, float dphr);
int	tauph_transform_3d(int nphx, int nphy, int ns, int ntau, int ntr,
        int itaus, int itaue,
        float phxmin, float phymin, float dphx, float dphy, float dt,
        float **tcmp, float *vrms, float **header, float ***pdat);
int	tauph_transform_2d( int ns, int ntau, int nphx, int ntr, int itaus, int itaue,
        float **tcmp, float **pdat, float *vrms, float *offset,
        float dt, float phxmin, float dphx);
int	tauph_transform_3d_to_2d_mute( int ns, int ntau, int nphr, int ntr, int itaus, int itaue,
        float **tcmp, float **pdat, float *vrms, float *offset,
        float dt, float phrmin, float dphr);
int	tauph_transform_3d_to_2d_robust( int ns, int ntau, int nphr, int ntr, int itaus, int itaue,
        float **tcmp, float **pdat, float *vrms, float *offset,
        float dt, float phrmin, float dphr);
// [END] in file: tauph_transform.c //.

// [BEGIN] in file: local_tauph_transform.c //.
void postProcessingTauPspectrum_diag(
        int nphr,
        int ntau,
        float **taup2d_lslrt
        );

int postProcessingTauPspectrum(
        float **taup2d_frss,
        int nphr,
        int ntau,
        float **taup2d_lslrt
        );

int	slantStack_TD_2D_omp(
        float **trace,
        int ntrace,
        int ns,
        float dt,
        float *coor,
        float x0,
        float **tauppanel,
        int npsr,
        float psrmin,
        float dpsr
        );

int	slantStack_TD_2D(
        float **trace,
        int ntrace,
        int ns,
        float dt,
        float *coor,
        float x0,
        float **tauppanel,
        int npsr,
        float psrmin,
        float dpsr
        );

int	apply_LSLRT_to_CMPgather(
        int	ntrCMP,
        float   *offset,
        fbsegy_std  *suhdr,
        int	ns,
        float	dt_s,
        float	**gather,
        float	offsetWidth,
        int	nphr,
        float	phrmin,
        float	dphr,
        int	ntau,
        float	**taupSpectrum,
        int	reconstrutDataFlag,
        float	**gather_out,
        int	verbose
        );

int	apply_LocalLRT_to_CMPgather(
        int	ntrCMP,
        float   *offset,
        int	ns,
        float	dt_s,
        float	**gather,
        float	offsetWidth,
        int	nphr,
        float	phrmin,
        float	dphr,
        int	ntau,
        float	**taupSpectrum,
        int	verbose
        );
// [END] in file: local_tauph_transform.c //.

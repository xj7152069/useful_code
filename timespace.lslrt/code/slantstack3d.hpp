/************update in 2021.01.07**************/
/*
    
***********************************************/

#ifdef __cplusplus
extern "C" {
#endif
void Slantstack_CG_2D(float **recoverdata, float **tauppanel, float **trace,\
 float *coor, int ntrace, float x0, int ns, int ntau,float dt,\
 int npsr,float psrmin, float dpsr,\
 float factor_L2,float factor_L1, \
 int iterations_num, float residual_ratio, bool dorecover);

void apply_LocalLSLRT_to_CMPgather(int ntrCMP, float *offset_site, int ns, float dt_s,\
 float **gather, float offsetWidth, int nphr, float phrmin, float dphr, int ntau,\
 float **taup2d_lslrt, int mintrace, float factor_L2,float factor_L1, \
 int iterations_num, float residual_ratio);

void apply_GlobalLSLRT_to_CMPgather(int ntrCMP, float *offset, int ns, float dt_s,\
 float **gather, int nphr, float phrmin, float dphr, int ntau,float **taup2d_lslrt,\
 float factor_L2,float factor_L1, \
 int iterations_num, float residual_ratio);

void medianfilter(float **mat, int n1, int n2, int w1, int w2, float lmd);

#ifdef __cplusplus
}
#endif
// 2-D Slant-stack in the time - space domain.
/* discription: void Slantstack_CG_2D()
float **recoverdata;
the OUTPUT recover local seismic-gather, i.e. recoverdata[ntrace][ns].
This function will do the local LRT over local traces .
float **tauppanel; 
the OUTPUT tau-p spectrum, i.e. tauppanel[npsr][ns].
float **trace;
the INPUT local seismic-gather, i.e. trace[ntrace][ns].
float *coor;
coor[ntrace] stores the offset of each trace .
int ntrace;
the trace-number of the INPUT gather.
float x0;
Origin of spatial coordinates
int ns;
data time-sampling number.
int ntau;
tau-p spectrum time-sampling number.
float dt;
time-sampling interval, unit of dt should be second (s),
int npsr;
the ray-parameter number of tau-p spectrum.
float psrmin;
the minimum ray-parameter, Origin of ray-parameter(unit is s/m) .
float dpsr;
the sampling interval of ray-parameter， (unit is s/m) .
float factor_L2; 
L2, Tikhonov regularization parameter,default 1.0
float factor_L1; default 
L1, sparsity  regularization parameter, default 0.0
int iterations_num; default 45
the Maximum number of iterations of conjugate-gradient method
float residual_ratio; default 0.5
Conjugate gradient method iteration error, The smaller the more accurate
bool dorecover; default true
Whether to OUTPUT recover local seismic-gather
*/

/*discription: void medianfilter()
float **mat;
the INPUT seismic-gather, i.e. mat[ntrace][ns].
int n1;
Sampling of the first dimension of the INPUT seismic-gather.
int n2;
Sampling of the second dimension of the INPUT seismic-gather.
int w1;
The length of the filter window in the first dimension.
int w2;
The length of the filter window in the second dimension.
float lmd; default 3.0
Filter sensitivity, the smaller the value, the stronger the filtering effect;
*/



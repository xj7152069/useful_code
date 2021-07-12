#ifndef IO_FUNC_H
#define IO_FUNC_H
#include"ShengShen_head.h"


double profile_IO_func(struct PARAMETER *param, float **record);

double time_slice_IO_func(struct PARAMETER *param, float **time_slice_data, char name[]);

double PS_IO_func(struct PARAMETER *param, float **Vx,  float **Vz, int num_shot, int num_r);

#endif


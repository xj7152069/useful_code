#pragma once
#include"ShengShen_head.h"

//double rickerwavelet_func(int Lw, float freq, float dt, float *signal);
double rickerwavelet_func(struct PARAMETER *param, float *signal);

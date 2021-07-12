#ifndef MODEL_FUNC_H
#define MODEL_FUNC_H
#include "ShengShen_head.h"
#include "Stability_func.h"

double model_func(struct PARAMETER *param,float **V, float **V_p, float **V_s, float **DEN, float **lamda, float **miu);

#endif


#include "fbCommon.h"

float	sinc_interpolation ( float *trace , int lt , float dt , float time_s_r )
{

    float sinc,x,tdif,t_value,ddt;
    int ncurrent,it_start,it_final,lw,iw;

    lw=10;
    t_value=0.0;
    ddt = PI/dt;

    if(time_s_r>(lt*dt))
        return 0.0;
    else
    {
        ncurrent=(int)(time_s_r/dt+0.5);
        it_start=ncurrent-lw;
        if(it_start<1)
            it_start=1;
        it_final=ncurrent+lw;
        if(it_final>lt)
            it_final=lt;

        for(iw=it_start;iw<=it_final;iw++)
        {
            tdif=time_s_r-iw*dt;
            if(tdif==0.0)
            {
                //sinc=1.0;
                t_value+=*(trace+iw);  //*sinc;
            }
            else
            {
                //x=PI*tdif/dt;
                x = tdif * ddt;
                sinc=sin(x)/x;
                t_value+=*(trace+iw)*sinc;
            }
        }
        return t_value;
    }
}

#include "../lib/Core_func.h"

double core01_func(float **Vx_now, float **Vx_aft,
        float **Vxx_now, float **Vxx_aft,
        float **Vxz_now, float **Vxz_aft,
        float **Vz_now, float **Vz_aft,
        float **Vzx_now, float **Vzx_aft,
        float **Vzz_now, float **Vzz_aft,
        float **Txx_now, float **Txz_now,float **Tzz_now,
        float wavelet, float **DEN,float **absorbx,float **absorbz,
        struct PARAMETER *param)
{

    //center
    //Vx
    int threads_count=param->omp_ncores;
# pragma omp parallel for num_threads(threads_count)
    for(int i=param->PML;i<param->NX+param->PML; i++)
    {
        for(int j=param->PML; j<param->PML+param->NZ; j++)
        {
            Vxx_aft[i][j] = Vxx_now[i][j]
                +1.0/DEN[i][j]*(param->dt*1.0/param->dx)*(
                        param->C1*(Txx_now[i+0][j]-Txx_now[i-1][j])+
                        param->C2*(Txx_now[i+1][j]-Txx_now[i-2][j])+
                        param->C3*(Txx_now[i+2][j]-Txx_now[i-3][j])+
                        param->C4*(Txx_now[i+3][j]-Txx_now[i-4][j])+
                        param->C5*(Txx_now[i+4][j]-Txx_now[i-5][j])+
                        param->C6*(Txx_now[i+5][j]-Txx_now[i-6][j])+
                        param->C7*(Txx_now[i+6][j]-Txx_now[i-7][j])+
                        param->C8*(Txx_now[i+7][j]-Txx_now[i-8][j]));

            Vxz_aft[i][j] = Vxz_now[i][j]
                +1.0/DEN[i][j]*(param->dt*1.0/param->dz)*(
                        param->C1*(Txz_now[i][j+0]-Txz_now[i][j-1])+
                        param->C2*(Txz_now[i][j+1]-Txz_now[i][j-2])+
                        param->C3*(Txz_now[i][j+2]-Txz_now[i][j-3])+
                        param->C4*(Txz_now[i][j+3]-Txz_now[i][j-4])+
                        param->C5*(Txz_now[i][j+4]-Txz_now[i][j-5])+
                        param->C6*(Txz_now[i][j+5]-Txz_now[i][j-6])+
                        param->C7*(Txz_now[i][j+6]-Txz_now[i][j-7])+
                        param->C8*(Txz_now[i][j+7]-Txz_now[i][j-8]));
        }
    }

    //Vz
# pragma omp parallel for num_threads(threads_count)
    for(int i=param->PML;i<param->NX+param->PML; i++)
    {
        for(int j=param->PML; j<param->PML+param->NZ; j++)
        {
            Vzx_aft[i][j] = Vzx_now[i][j]
                +1.0/DEN[i][j]*(param->dt*1.0/param->dx)*(
                        param->C1*(Txz_now[i+1][j]-Txz_now[i-0][j])+
                        param->C2*(Txz_now[i+2][j]-Txz_now[i-1][j])+
                        param->C3*(Txz_now[i+3][j]-Txz_now[i-2][j])+
                        param->C4*(Txz_now[i+4][j]-Txz_now[i-3][j])+
                        param->C5*(Txz_now[i+5][j]-Txz_now[i-4][j])+
                        param->C6*(Txz_now[i+6][j]-Txz_now[i-5][j])+
                        param->C7*(Txz_now[i+7][j]-Txz_now[i-6][j])+
                        param->C8*(Txz_now[i+8][j]-Txz_now[i-7][j]));

            Vzz_aft[i][j] = Vzz_now[i][j]
                +1.0/DEN[i][j]*(param->dt*1.0/param->dz)*(
                        param->C1*(Tzz_now[i][j+1]-Tzz_now[i][j-0])+
                        param->C2*(Tzz_now[i][j+2]-Tzz_now[i][j-1])+
                        param->C3*(Tzz_now[i][j+3]-Tzz_now[i][j-2])+
                        param->C4*(Tzz_now[i][j+4]-Tzz_now[i][j-3])+
                        param->C5*(Tzz_now[i][j+5]-Tzz_now[i][j-4])+
                        param->C6*(Tzz_now[i][j+6]-Tzz_now[i][j-5])+
                        param->C7*(Tzz_now[i][j+7]-Tzz_now[i][j-6])+
                        param->C8*(Tzz_now[i][j+8]-Tzz_now[i][j-7]));
        }
    }


    //Left
    //Vx
# pragma omp parallel for num_threads(threads_count)
    for(int i=8; i<param->PML; i++)
    {
        for(int j=param->PML; j<param->NZ+ param->PML; j++)
        {
            Vxx_aft[i][j] = 1.0/(1.0+0.5*param->dt*absorbx[i][j])*(
                    (1.0-0.5*param->dt*absorbx[i][j])*Vxx_now[i][j]
                    +1.0/DEN[i][j]*(param->dt*1.0/param->dx)*(
                        param->C1*(Txx_now[i+0][j]-Txx_now[i-1][j])+
                        param->C2*(Txx_now[i+1][j]-Txx_now[i-2][j])+
                        param->C3*(Txx_now[i+2][j]-Txx_now[i-3][j])+
                        param->C4*(Txx_now[i+3][j]-Txx_now[i-4][j])+
                        param->C5*(Txx_now[i+4][j]-Txx_now[i-5][j])+
                        param->C6*(Txx_now[i+5][j]-Txx_now[i-6][j])+
                        param->C7*(Txx_now[i+6][j]-Txx_now[i-7][j])+
                        param->C8*(Txx_now[i+7][j]-Txx_now[i-8][j])));

            Vxz_aft[i][j] = Vxz_now[i][j]
                +1.0/DEN[i][j]*(param->dt*1.0/param->dz)*(
                        param->C1*(Txz_now[i][j+0]-Txz_now[i][j-1])+
                        param->C2*(Txz_now[i][j+1]-Txz_now[i][j-2])+
                        param->C3*(Txz_now[i][j+2]-Txz_now[i][j-3])+
                        param->C4*(Txz_now[i][j+3]-Txz_now[i][j-4])+
                        param->C5*(Txz_now[i][j+4]-Txz_now[i][j-5])+
                        param->C6*(Txz_now[i][j+5]-Txz_now[i][j-6])+
                        param->C7*(Txz_now[i][j+6]-Txz_now[i][j-7])+
                        param->C8*(Txz_now[i][j+7]-Txz_now[i][j-8]));
        }
    }
    //Vz
# pragma omp parallel for num_threads(threads_count)
    for(int i=8; i<param->PML; i++)
    {
        for(int j=param->PML; j<param->NZ+ param->PML; j++)
        {
            Vzx_aft[i][j] = 1.0/(1.0+0.5*param->dt*absorbx[i][j])*(
                    (1.0-0.5*param->dt*absorbx[i][j])*Vzx_now[i][j]
                    +1.0/DEN[i][j]*(param->dt*1.0/param->dx)*(
                        param->C1*(Txz_now[i+1][j]-Txz_now[i-0][j])+
                        param->C2*(Txz_now[i+2][j]-Txz_now[i-1][j])+
                        param->C3*(Txz_now[i+3][j]-Txz_now[i-2][j])+
                        param->C4*(Txz_now[i+4][j]-Txz_now[i-3][j])+
                        param->C5*(Txz_now[i+5][j]-Txz_now[i-4][j])+
                        param->C6*(Txz_now[i+6][j]-Txz_now[i-5][j])+
                        param->C7*(Txz_now[i+7][j]-Txz_now[i-6][j])+
                        param->C8*(Txz_now[i+8][j]-Txz_now[i-7][j])));

            Vzz_aft[i][j] = Vzz_now[i][j]
                +1.0/DEN[i][j]*(param->dt*1.0/param->dz)*(
                        param->C1*(Tzz_now[i][j+1]-Tzz_now[i][j-0])+
                        param->C2*(Tzz_now[i][j+2]-Tzz_now[i][j-1])+
                        param->C3*(Tzz_now[i][j+3]-Tzz_now[i][j-2])+
                        param->C4*(Tzz_now[i][j+4]-Tzz_now[i][j-3])+
                        param->C5*(Tzz_now[i][j+5]-Tzz_now[i][j-4])+
                        param->C6*(Tzz_now[i][j+6]-Tzz_now[i][j-5])+
                        param->C7*(Tzz_now[i][j+7]-Tzz_now[i][j-6])+
                        param->C8*(Tzz_now[i][j+8]-Tzz_now[i][j-7]));
        }
    }



    //right
    //Vx
# pragma omp parallel for num_threads(threads_count)
    for(int i=param->NX+param->PML;i<param->Nx-8; i++)
    {
        for(int j=param->PML; j<param->NZ+ param->PML; j++)
        {
            Vxx_aft[i][j] = 1.0/(1.0+0.5*param->dt*absorbx[i][j])*(
                    (1.0-0.5*param->dt*absorbx[i][j])*Vxx_now[i][j]
                    +1.0/DEN[i][j]*(param->dt*1.0/param->dx)*(
                        param->C1*(Txx_now[i+0][j]-Txx_now[i-1][j])+
                        param->C2*(Txx_now[i+1][j]-Txx_now[i-2][j])+
                        param->C3*(Txx_now[i+2][j]-Txx_now[i-3][j])+
                        param->C4*(Txx_now[i+3][j]-Txx_now[i-4][j])+
                        param->C5*(Txx_now[i+4][j]-Txx_now[i-5][j])+
                        param->C6*(Txx_now[i+5][j]-Txx_now[i-6][j])+
                        param->C7*(Txx_now[i+6][j]-Txx_now[i-7][j])+
                        param->C8*(Txx_now[i+7][j]-Txx_now[i-8][j])));

            Vxz_aft[i][j] = Vxz_now[i][j]
                +1.0/DEN[i][j]*(param->dt*1.0/param->dz)*(
                        param->C1*(Txz_now[i][j+0]-Txz_now[i][j-1])+
                        param->C2*(Txz_now[i][j+1]-Txz_now[i][j-2])+
                        param->C3*(Txz_now[i][j+2]-Txz_now[i][j-3])+
                        param->C4*(Txz_now[i][j+3]-Txz_now[i][j-4])+
                        param->C5*(Txz_now[i][j+4]-Txz_now[i][j-5])+
                        param->C6*(Txz_now[i][j+5]-Txz_now[i][j-6])+
                        param->C7*(Txz_now[i][j+6]-Txz_now[i][j-7])+
                        param->C8*(Txz_now[i][j+7]-Txz_now[i][j-8]));
        }
    }
    //Vz
# pragma omp parallel for num_threads(threads_count)
    for(int i=param->NX+param->PML;i<param->Nx-8; i++)
    {
        for(int j=param->PML; j<param->NZ+ param->PML; j++)
        {
            Vzx_aft[i][j] = 1.0/(1.0+0.5*param->dt*absorbx[i][j])*(
                    (1.0-0.5*param->dt*absorbx[i][j])*Vzx_now[i][j]
                    +1.0/DEN[i][j]*(param->dt*1.0/param->dx)*(
                        param->C1*(Txz_now[i+1][j]-Txz_now[i-0][j])+
                        param->C2*(Txz_now[i+2][j]-Txz_now[i-1][j])+
                        param->C3*(Txz_now[i+3][j]-Txz_now[i-2][j])+
                        param->C4*(Txz_now[i+4][j]-Txz_now[i-3][j])+
                        param->C5*(Txz_now[i+5][j]-Txz_now[i-4][j])+
                        param->C6*(Txz_now[i+6][j]-Txz_now[i-5][j])+
                        param->C7*(Txz_now[i+7][j]-Txz_now[i-6][j])+
                        param->C8*(Txz_now[i+8][j]-Txz_now[i-7][j])));

            Vzz_aft[i][j] = Vzz_now[i][j]
                +1.0/DEN[i][j]*(param->dt*1.0/param->dz)*(
                        param->C1*(Tzz_now[i][j+1]-Tzz_now[i][j-0])+
                        param->C2*(Tzz_now[i][j+2]-Tzz_now[i][j-1])+
                        param->C3*(Tzz_now[i][j+3]-Tzz_now[i][j-2])+
                        param->C4*(Tzz_now[i][j+4]-Tzz_now[i][j-3])+
                        param->C5*(Tzz_now[i][j+5]-Tzz_now[i][j-4])+
                        param->C6*(Tzz_now[i][j+6]-Tzz_now[i][j-5])+
                        param->C7*(Tzz_now[i][j+7]-Tzz_now[i][j-6])+
                        param->C8*(Tzz_now[i][j+8]-Tzz_now[i][j-7]));
        }
    }

    //Up
    //Vx
# pragma omp parallel for num_threads(threads_count)
    for(int i=param->PML;i<param->NX + param->PML; i++)
    {
        for(int j=8; j<param->PML; j++)
        {
            Vxx_aft[i][j] = Vxx_now[i][j]
                +1.0/DEN[i][j]*(param->dt*1.0/param->dx)*(
                        param->C1*(Txx_now[i+0][j]-Txx_now[i-1][j])+
                        param->C2*(Txx_now[i+1][j]-Txx_now[i-2][j])+
                        param->C3*(Txx_now[i+2][j]-Txx_now[i-3][j])+
                        param->C4*(Txx_now[i+3][j]-Txx_now[i-4][j])+
                        param->C5*(Txx_now[i+4][j]-Txx_now[i-5][j])+
                        param->C6*(Txx_now[i+5][j]-Txx_now[i-6][j])+
                        param->C7*(Txx_now[i+6][j]-Txx_now[i-7][j])+
                        param->C8*(Txx_now[i+7][j]-Txx_now[i-8][j]));

            Vxz_aft[i][j] = 1.0/(1.0+0.5*param->dt*absorbz[i][j])*(
                    (1.0-0.5*param->dt*absorbz[i][j])*Vxz_now[i][j]
                    +1.0/DEN[i][j]*(param->dt*1.0/param->dz)*(
                        param->C1*(Txz_now[i][j+0]-Txz_now[i][j-1])+
                        param->C2*(Txz_now[i][j+1]-Txz_now[i][j-2])+
                        param->C3*(Txz_now[i][j+2]-Txz_now[i][j-3])+
                        param->C4*(Txz_now[i][j+3]-Txz_now[i][j-4])+
                        param->C5*(Txz_now[i][j+4]-Txz_now[i][j-5])+
                        param->C6*(Txz_now[i][j+5]-Txz_now[i][j-6])+
                        param->C7*(Txz_now[i][j+6]-Txz_now[i][j-7])+
                        param->C8*(Txz_now[i][j+7]-Txz_now[i][j-8])));
        }
    }
    //Vz
# pragma omp parallel for num_threads(threads_count)
    for(int i=param->PML;i<param->NX + param->PML; i++)
    {
        for(int j=8; j<param->PML; j++)
        {
            Vzx_aft[i][j] = Vzx_now[i][j]
                +1.0/DEN[i][j]*(param->dt*1.0/param->dx)*(
                        param->C1*(Txz_now[i+1][j]-Txz_now[i-0][j])+
                        param->C2*(Txz_now[i+2][j]-Txz_now[i-1][j])+
                        param->C3*(Txz_now[i+3][j]-Txz_now[i-2][j])+
                        param->C4*(Txz_now[i+4][j]-Txz_now[i-3][j])+
                        param->C5*(Txz_now[i+5][j]-Txz_now[i-4][j])+
                        param->C6*(Txz_now[i+6][j]-Txz_now[i-5][j])+
                        param->C7*(Txz_now[i+7][j]-Txz_now[i-6][j])+
                        param->C8*(Txz_now[i+8][j]-Txz_now[i-7][j]));

            Vzz_aft[i][j] = 1.0/(1.0+0.5*param->dt*absorbz[i][j])*(
                    (1.0-0.5*param->dt*absorbz[i][j])*Vzz_now[i][j]
                    +1.0/DEN[i][j]*(param->dt*1.0/param->dz)*(
                        param->C1*(Tzz_now[i][j+1]-Tzz_now[i][j-0])+
                        param->C2*(Tzz_now[i][j+2]-Tzz_now[i][j-1])+
                        param->C3*(Tzz_now[i][j+3]-Tzz_now[i][j-2])+
                        param->C4*(Tzz_now[i][j+4]-Tzz_now[i][j-3])+
                        param->C5*(Tzz_now[i][j+5]-Tzz_now[i][j-4])+
                        param->C6*(Tzz_now[i][j+6]-Tzz_now[i][j-5])+
                        param->C7*(Tzz_now[i][j+7]-Tzz_now[i][j-6])+
                        param->C8*(Tzz_now[i][j+8]-Tzz_now[i][j-7])));
        }
    }

    //under
    //Vx
# pragma omp parallel for num_threads(threads_count)
    for(int i=8;i<param->Nx-8; i++)
    {
        for(int j=param->NZ+param->PML; j<param->Nz-8; j++)
        {
            Vxx_aft[i][j] = Vxx_now[i][j]
                +1.0/DEN[i][j]*(param->dt*1.0/param->dx)*(
                        param->C1*(Txx_now[i+0][j]-Txx_now[i-1][j])+
                        param->C2*(Txx_now[i+1][j]-Txx_now[i-2][j])+
                        param->C3*(Txx_now[i+2][j]-Txx_now[i-3][j])+
                        param->C4*(Txx_now[i+3][j]-Txx_now[i-4][j])+
                        param->C5*(Txx_now[i+4][j]-Txx_now[i-5][j])+
                        param->C6*(Txx_now[i+5][j]-Txx_now[i-6][j])+
                        param->C7*(Txx_now[i+6][j]-Txx_now[i-7][j])+
                        param->C8*(Txx_now[i+7][j]-Txx_now[i-8][j]));

            Vxz_aft[i][j] = 1.0/(1.0+0.5*param->dt*absorbz[i][j])*(
                    (1.0-0.5*param->dt*absorbz[i][j])*Vxz_now[i][j]
                    +1.0/DEN[i][j]*(param->dt*1.0/param->dz)*(
                        param->C1*(Txz_now[i][j+0]-Txz_now[i][j-1])+
                        param->C2*(Txz_now[i][j+1]-Txz_now[i][j-2])+
                        param->C3*(Txz_now[i][j+2]-Txz_now[i][j-3])+
                        param->C4*(Txz_now[i][j+3]-Txz_now[i][j-4])+
                        param->C5*(Txz_now[i][j+4]-Txz_now[i][j-5])+
                        param->C6*(Txz_now[i][j+5]-Txz_now[i][j-6])+
                        param->C7*(Txz_now[i][j+6]-Txz_now[i][j-7])+
                        param->C8*(Txz_now[i][j+7]-Txz_now[i][j-8])));
        }
    }
    //Vz
# pragma omp parallel for num_threads(threads_count)
    for(int i=8;i<param->Nx-8; i++)
    {
        for(int j=param->NZ+param->PML; j<param->Nz-8; j++)
        {
            Vzx_aft[i][j] = Vzx_now[i][j]
                +1.0/DEN[i][j]*(param->dt*1.0/param->dx)*(
                        param->C1*(Txz_now[i+1][j]-Txz_now[i-0][j])+
                        param->C2*(Txz_now[i+2][j]-Txz_now[i-1][j])+
                        param->C3*(Txz_now[i+3][j]-Txz_now[i-2][j])+
                        param->C4*(Txz_now[i+4][j]-Txz_now[i-3][j])+
                        param->C5*(Txz_now[i+5][j]-Txz_now[i-4][j])+
                        param->C6*(Txz_now[i+6][j]-Txz_now[i-5][j])+
                        param->C7*(Txz_now[i+7][j]-Txz_now[i-6][j])+
                        param->C8*(Txz_now[i+8][j]-Txz_now[i-7][j]));

            Vzz_aft[i][j] = 1.0/(1.0+0.5*param->dt*absorbz[i][j])*(
                    (1.0-0.5*param->dt*absorbz[i][j])*Vzz_now[i][j]
                    +1.0/DEN[i][j]*(param->dt*1.0/param->dz)*(
                        param->C1*(Tzz_now[i][j+1]-Tzz_now[i][j-0])+
                        param->C2*(Tzz_now[i][j+2]-Tzz_now[i][j-1])+
                        param->C3*(Tzz_now[i][j+3]-Tzz_now[i][j-2])+
                        param->C4*(Tzz_now[i][j+4]-Tzz_now[i][j-3])+
                        param->C5*(Tzz_now[i][j+5]-Tzz_now[i][j-4])+
                        param->C6*(Tzz_now[i][j+6]-Tzz_now[i][j-5])+
                        param->C7*(Tzz_now[i][j+7]-Tzz_now[i][j-6])+
                        param->C8*(Tzz_now[i][j+8]-Tzz_now[i][j-7])));
        }
    }


    //coner_left_Up
    //Vx
# pragma omp parallel for num_threads(threads_count)
    for(int i=8;i<param->PML; i++)
    {
        for(int j=8; j<param->PML; j++)
        {
            Vxx_aft[i][j] = 1.0/(1.0+0.5*param->dt*absorbx[i][j])*(
                    (1.0-0.5*param->dt*absorbx[i][j])*Vxx_now[i][j]
                    +1.0/DEN[i][j]*(param->dt*1.0/param->dx)*(
                        param->C1*(Txx_now[i+0][j]-Txx_now[i-1][j])+
                        param->C2*(Txx_now[i+1][j]-Txx_now[i-2][j])+
                        param->C3*(Txx_now[i+2][j]-Txx_now[i-3][j])+
                        param->C4*(Txx_now[i+3][j]-Txx_now[i-4][j])+
                        param->C5*(Txx_now[i+4][j]-Txx_now[i-5][j])+
                        param->C6*(Txx_now[i+5][j]-Txx_now[i-6][j])+
                        param->C7*(Txx_now[i+6][j]-Txx_now[i-7][j])+
                        param->C8*(Txx_now[i+7][j]-Txx_now[i-8][j])));

            Vxz_aft[i][j] = 1.0/(1.0+0.5*param->dt*absorbz[i][j])*(
                    (1.0-0.5*param->dt*absorbz[i][j])*Vxz_now[i][j]
                    +1.0/DEN[i][j]*(param->dt*1.0/param->dz)*(
                        param->C1*(Txz_now[i][j+0]-Txz_now[i][j-1])+
                        param->C2*(Txz_now[i][j+1]-Txz_now[i][j-2])+
                        param->C3*(Txz_now[i][j+2]-Txz_now[i][j-3])+
                        param->C4*(Txz_now[i][j+3]-Txz_now[i][j-4])+
                        param->C5*(Txz_now[i][j+4]-Txz_now[i][j-5])+
                        param->C6*(Txz_now[i][j+5]-Txz_now[i][j-6])+
                        param->C7*(Txz_now[i][j+6]-Txz_now[i][j-7])+
                        param->C8*(Txz_now[i][j+7]-Txz_now[i][j-8])));
        }
    }
    //Vz
# pragma omp parallel for num_threads(threads_count)
    for(int i=8;i<param->PML; i++)
    {
        for(int j=8; j<param->PML; j++)
        {
            Vzx_aft[i][j] = 1.0/(1.0+0.5*param->dt*absorbx[i][j])*(
                    (1.0-0.5*param->dt*absorbx[i][j])*Vzx_now[i][j]
                    +1.0/DEN[i][j]*(param->dt*1.0/param->dx)*(
                        param->C1*(Txz_now[i+1][j]-Txz_now[i-0][j])+
                        param->C2*(Txz_now[i+2][j]-Txz_now[i-1][j])+
                        param->C3*(Txz_now[i+3][j]-Txz_now[i-2][j])+
                        param->C4*(Txz_now[i+4][j]-Txz_now[i-3][j])+
                        param->C5*(Txz_now[i+5][j]-Txz_now[i-4][j])+
                        param->C6*(Txz_now[i+6][j]-Txz_now[i-5][j])+
                        param->C7*(Txz_now[i+7][j]-Txz_now[i-6][j])+
                        param->C8*(Txz_now[i+8][j]-Txz_now[i-7][j])));


            Vzz_aft[i][j] = 1.0/(1.0+0.5*param->dt*absorbz[i][j])*(
                    (1.0-0.5*param->dt*absorbz[i][j])*Vzz_now[i][j]
                    +1.0/DEN[i][j]*(param->dt*1.0/param->dz)*(
                        param->C1*(Tzz_now[i][j+1]-Tzz_now[i][j-0])+
                        param->C2*(Tzz_now[i][j+2]-Tzz_now[i][j-1])+
                        param->C3*(Tzz_now[i][j+3]-Tzz_now[i][j-2])+
                        param->C4*(Tzz_now[i][j+4]-Tzz_now[i][j-3])+
                        param->C5*(Tzz_now[i][j+5]-Tzz_now[i][j-4])+
                        param->C6*(Tzz_now[i][j+6]-Tzz_now[i][j-5])+
                        param->C7*(Tzz_now[i][j+7]-Tzz_now[i][j-6])+
                        param->C8*(Tzz_now[i][j+8]-Tzz_now[i][j-7])));
        }
    }

    //coner_right_Up
    //Vx
# pragma omp parallel for num_threads(threads_count)
    for(int i=param->NX+ param->PML;i<param->Nx-8; i++)
    {
        for(int j=8; j<param->PML; j++)
        {
            Vxx_aft[i][j] = 1.0/(1.0+0.5*param->dt*absorbx[i][j])*(
                    (1.0-0.5*param->dt*absorbx[i][j])*Vxx_now[i][j]
                    +1.0/DEN[i][j]*(param->dt*1.0/param->dx)*(
                        param->C1*(Txx_now[i+0][j]-Txx_now[i-1][j])+
                        param->C2*(Txx_now[i+1][j]-Txx_now[i-2][j])+
                        param->C3*(Txx_now[i+2][j]-Txx_now[i-3][j])+
                        param->C4*(Txx_now[i+3][j]-Txx_now[i-4][j])+
                        param->C5*(Txx_now[i+4][j]-Txx_now[i-5][j])+
                        param->C6*(Txx_now[i+5][j]-Txx_now[i-6][j])+
                        param->C7*(Txx_now[i+6][j]-Txx_now[i-7][j])+
                        param->C8*(Txx_now[i+7][j]-Txx_now[i-8][j])));

            Vxz_aft[i][j] = 1.0/(1.0+0.5*param->dt*absorbz[i][j])*(
                    (1.0-0.5*param->dt*absorbz[i][j])*Vxz_now[i][j]
                    +1.0/DEN[i][j]*(param->dt*1.0/param->dz)*(
                        param->C1*(Txz_now[i][j+0]-Txz_now[i][j-1])+
                        param->C2*(Txz_now[i][j+1]-Txz_now[i][j-2])+
                        param->C3*(Txz_now[i][j+2]-Txz_now[i][j-3])+
                        param->C4*(Txz_now[i][j+3]-Txz_now[i][j-4])+
                        param->C5*(Txz_now[i][j+4]-Txz_now[i][j-5])+
                        param->C6*(Txz_now[i][j+5]-Txz_now[i][j-6])+
                        param->C7*(Txz_now[i][j+6]-Txz_now[i][j-7])+
                        param->C8*(Txz_now[i][j+7]-Txz_now[i][j-8])));
        }
    }
    //Vz
# pragma omp parallel for num_threads(threads_count)
    for(int i=param->NX+ param->PML;i<param->Nx-8; i++)
    {
        for(int j=8; j<param->PML; j++)
        {
            Vzx_aft[i][j] = 1.0/(1.0+0.5*param->dt*absorbx[i][j])*(
                    (1.0-0.5*param->dt*absorbx[i][j])*Vzx_now[i][j]
                    +1.0/DEN[i][j]*(param->dt*1.0/param->dx)*(
                        param->C1*(Txz_now[i+1][j]-Txz_now[i-0][j])+
                        param->C2*(Txz_now[i+2][j]-Txz_now[i-1][j])+
                        param->C3*(Txz_now[i+3][j]-Txz_now[i-2][j])+
                        param->C4*(Txz_now[i+4][j]-Txz_now[i-3][j])+
                        param->C5*(Txz_now[i+5][j]-Txz_now[i-4][j])+
                        param->C6*(Txz_now[i+6][j]-Txz_now[i-5][j])+
                        param->C7*(Txz_now[i+7][j]-Txz_now[i-6][j])+
                        param->C8*(Txz_now[i+8][j]-Txz_now[i-7][j])));

            Vzz_aft[i][j] = 1.0/(1.0+0.5*param->dt*absorbz[i][j])*(
                    (1.0-0.5*param->dt*absorbz[i][j])*Vzz_now[i][j]
                    +1.0/DEN[i][j]*(param->dt*1.0/param->dz)*(
                        param->C1*(Tzz_now[i][j+1]-Tzz_now[i][j-0])+
                        param->C2*(Tzz_now[i][j+2]-Tzz_now[i][j-1])+
                        param->C3*(Tzz_now[i][j+3]-Tzz_now[i][j-2])+
                        param->C4*(Tzz_now[i][j+4]-Tzz_now[i][j-3])+
                        param->C5*(Tzz_now[i][j+5]-Tzz_now[i][j-4])+
                        param->C6*(Tzz_now[i][j+6]-Tzz_now[i][j-5])+
                        param->C7*(Tzz_now[i][j+7]-Tzz_now[i][j-6])+
                        param->C8*(Tzz_now[i][j+8]-Tzz_now[i][j-7])));
        }
    }

    //coner_right_Under
    //Vx
# pragma omp parallel for num_threads(threads_count)
    for(int i=param->NX+ param->PML;i<param->Nx-8; i++)
    {
        for(int j=param->NZ+ param->PML; j<param->Nz-8; j++)
        {
            Vxx_aft[i][j] = 1.0/(1.0+0.5*param->dt*absorbx[i][j])*(
                    (1.0-0.5*param->dt*absorbx[i][j])*Vxx_now[i][j]
                    +1.0/DEN[i][j]*(param->dt*1.0/param->dx)*(
                        param->C1*(Txx_now[i+0][j]-Txx_now[i-1][j])+
                        param->C2*(Txx_now[i+1][j]-Txx_now[i-2][j])+
                        param->C3*(Txx_now[i+2][j]-Txx_now[i-3][j])+
                        param->C4*(Txx_now[i+3][j]-Txx_now[i-4][j])+
                        param->C5*(Txx_now[i+4][j]-Txx_now[i-5][j])+
                        param->C6*(Txx_now[i+5][j]-Txx_now[i-6][j])+
                        param->C7*(Txx_now[i+6][j]-Txx_now[i-7][j])+
                        param->C8*(Txx_now[i+7][j]-Txx_now[i-8][j])));

            Vxz_aft[i][j] = 1.0/(1.0+0.5*param->dt*absorbz[i][j])*(
                    (1.0-0.5*param->dt*absorbz[i][j])*Vxz_now[i][j]
                    +1.0/DEN[i][j]*(param->dt*1.0/param->dz)*(
                        param->C1*(Txz_now[i][j+0]-Txz_now[i][j-1])+
                        param->C2*(Txz_now[i][j+1]-Txz_now[i][j-2])+
                        param->C3*(Txz_now[i][j+2]-Txz_now[i][j-3])+
                        param->C4*(Txz_now[i][j+3]-Txz_now[i][j-4])+
                        param->C5*(Txz_now[i][j+4]-Txz_now[i][j-5])+
                        param->C6*(Txz_now[i][j+5]-Txz_now[i][j-6])+
                        param->C7*(Txz_now[i][j+6]-Txz_now[i][j-7])+
                        param->C8*(Txz_now[i][j+7]-Txz_now[i][j-8])));
        }
    }
    //Vz
# pragma omp parallel for num_threads(threads_count)
    for(int i=param->NX+ param->PML;i<param->Nx-8; i++)
    {
        for(int j=param->NZ+ param->PML; j<param->Nz-8; j++)
        {
            Vzx_aft[i][j] = 1.0/(1.0+0.5*param->dt*absorbx[i][j])*(
                    (1.0-0.5*param->dt*absorbx[i][j])*Vzx_now[i][j]
                    +1.0/DEN[i][j]*(param->dt*1.0/param->dx)*(
                        param->C1*(Txz_now[i+1][j]-Txz_now[i-0][j])+
                        param->C2*(Txz_now[i+2][j]-Txz_now[i-1][j])+
                        param->C3*(Txz_now[i+3][j]-Txz_now[i-2][j])+
                        param->C4*(Txz_now[i+4][j]-Txz_now[i-3][j])+
                        param->C5*(Txz_now[i+5][j]-Txz_now[i-4][j])+
                        param->C6*(Txz_now[i+6][j]-Txz_now[i-5][j])+
                        param->C7*(Txz_now[i+7][j]-Txz_now[i-6][j])+
                        param->C8*(Txz_now[i+8][j]-Txz_now[i-7][j])));

            Vzz_aft[i][j] = 1.0/(1.0+0.5*param->dt*absorbz[i][j])*(
                    (1.0-0.5*param->dt*absorbz[i][j])*Vzz_now[i][j]
                    +1.0/DEN[i][j]*(param->dt*1.0/param->dz)*(
                        param->C1*(Tzz_now[i][j+1]-Tzz_now[i][j-0])+
                        param->C2*(Tzz_now[i][j+2]-Tzz_now[i][j-1])+
                        param->C3*(Tzz_now[i][j+3]-Tzz_now[i][j-2])+
                        param->C4*(Tzz_now[i][j+4]-Tzz_now[i][j-3])+
                        param->C5*(Tzz_now[i][j+5]-Tzz_now[i][j-4])+
                        param->C6*(Tzz_now[i][j+6]-Tzz_now[i][j-5])+
                        param->C7*(Tzz_now[i][j+7]-Tzz_now[i][j-6])+
                        param->C8*(Tzz_now[i][j+8]-Tzz_now[i][j-7])));
        }
    }

    //coner_left_Under
    //Vx
# pragma omp parallel for num_threads(threads_count)
    for(int i=8;i<param->PML; i++)
    {
        for(int j=param->NZ+ param->PML; j<param->Nz-8; j++)
        {
            Vxx_aft[i][j] = 1.0/(1.0+0.5*param->dt*absorbx[i][j])*(
                    (1.0-0.5*param->dt*absorbx[i][j])*Vxx_now[i][j]
                    +1.0/DEN[i][j]*(param->dt*1.0/param->dx)*(
                        param->C1*(Txx_now[i+0][j]-Txx_now[i-1][j])+
                        param->C2*(Txx_now[i+1][j]-Txx_now[i-2][j])+
                        param->C3*(Txx_now[i+2][j]-Txx_now[i-3][j])+
                        param->C4*(Txx_now[i+3][j]-Txx_now[i-4][j])+
                        param->C5*(Txx_now[i+4][j]-Txx_now[i-5][j])+
                        param->C6*(Txx_now[i+5][j]-Txx_now[i-6][j])+
                        param->C7*(Txx_now[i+6][j]-Txx_now[i-7][j])+
                        param->C8*(Txx_now[i+7][j]-Txx_now[i-8][j])));

            Vxz_aft[i][j] = 1.0/(1.0+0.5*param->dt*absorbz[i][j])*(
                    (1.0-0.5*param->dt*absorbz[i][j])*Vxz_now[i][j]
                    +1.0/DEN[i][j]*(param->dt*1.0/param->dz)*(
                        param->C1*(Txz_now[i][j+0]-Txz_now[i][j-1])+
                        param->C2*(Txz_now[i][j+1]-Txz_now[i][j-2])+
                        param->C3*(Txz_now[i][j+2]-Txz_now[i][j-3])+
                        param->C4*(Txz_now[i][j+3]-Txz_now[i][j-4])+
                        param->C5*(Txz_now[i][j+4]-Txz_now[i][j-5])+
                        param->C6*(Txz_now[i][j+5]-Txz_now[i][j-6])+
                        param->C7*(Txz_now[i][j+6]-Txz_now[i][j-7])+
                        param->C8*(Txz_now[i][j+7]-Txz_now[i][j-8])));
        }
    }
    //Vz
# pragma omp parallel for num_threads(threads_count)
    for(int i=8;i<param->PML; i++)
    {
        for(int j=param->NZ+ param->PML; j<param->Nz-8; j++)
        {
            Vzx_aft[i][j] = 1.0/(1.0+0.5*param->dt*absorbx[i][j])*(
                    (1.0-0.5*param->dt*absorbx[i][j])*Vzx_now[i][j]
                    +1.0/DEN[i][j]*(param->dt*1.0/param->dx)*(
                        param->C1*(Txz_now[i+1][j]-Txz_now[i-0][j])+
                        param->C2*(Txz_now[i+2][j]-Txz_now[i-1][j])+
                        param->C3*(Txz_now[i+3][j]-Txz_now[i-2][j])+
                        param->C4*(Txz_now[i+4][j]-Txz_now[i-3][j])+
                        param->C5*(Txz_now[i+5][j]-Txz_now[i-4][j])+
                        param->C6*(Txz_now[i+6][j]-Txz_now[i-5][j])+
                        param->C7*(Txz_now[i+7][j]-Txz_now[i-6][j])+
                        param->C8*(Txz_now[i+8][j]-Txz_now[i-7][j])));

            Vzz_aft[i][j] = 1.0/(1.0+0.5*param->dt*absorbz[i][j])*(
                    (1.0-0.5*param->dt*absorbz[i][j])*Vzz_now[i][j]
                    +1.0/DEN[i][j]*(param->dt*1.0/param->dz)*(
                        param->C1*(Tzz_now[i][j+1]-Tzz_now[i][j-0])+
                        param->C2*(Tzz_now[i][j+2]-Tzz_now[i][j-1])+
                        param->C3*(Tzz_now[i][j+3]-Tzz_now[i][j-2])+
                        param->C4*(Tzz_now[i][j+4]-Tzz_now[i][j-3])+
                        param->C5*(Tzz_now[i][j+5]-Tzz_now[i][j-4])+
                        param->C6*(Tzz_now[i][j+6]-Tzz_now[i][j-5])+
                        param->C7*(Tzz_now[i][j+7]-Tzz_now[i][j-6])+
                        param->C8*(Tzz_now[i][j+8]-Tzz_now[i][j-7])));
        }
    }


    //Vx\Vz
# pragma omp parallel for num_threads(threads_count)
    for(int i=0; i<param->Nx; i++)
    {
        for(int j=0; j<param->Nz; j++)
        {
            Vx_aft[i][j] = Vxx_aft[i][j] + Vxz_aft[i][j];
        }
    }
# pragma omp parallel for num_threads(threads_count)
    for(int i=0;i<param->Nx; i++)
    {
        for(int j=0; j<param->Nz; j++)
        {
            Vz_aft[i][j] = Vzx_aft[i][j] + Vzz_aft[i][j];
        }
    }

    //速度时间切片数值替换
# pragma omp parallel for num_threads(threads_count)
    for(int i=0; i<param->Nx; i++)
    {
        for(int j=0; j<param->Nz; j++)
        {
            Vx_now[i][j]  = Vx_aft[i][j];
            Vxx_now[i][j] = Vxx_aft[i][j];
            Vxz_now[i][j] = Vxz_aft[i][j];
            Vz_now[i][j]  = Vz_aft[i][j];
            Vzx_now[i][j] = Vzx_aft[i][j];
            Vzz_now[i][j] = Vzz_aft[i][j];
        }
    }

    return 0.0;
}


double core02_func(float **Txx_now, float **Txx_aft,
        float **Txxx_now, float **Txxx_aft,
        float **Txxz_now, float **Txxz_aft,
        float **Tzz_now, float **Tzz_aft,
        float **Tzzx_now, float **Tzzx_aft,
        float **Tzzz_now, float **Tzzz_aft,
        float **Txz_now, float **Txz_aft,
        float **Txzx_now, float **Txzx_aft,
        float **Txzz_now, float **Txzz_aft,
        float **Vx_now, float **Vz_now,
        float **lamda, float **miu, float wavelet,
        float **DEN,float **absorbx,float **absorbz,struct PARAMETER *param)
{
    int threads_count=param->omp_ncores;
    //Center
    //Txx
# pragma omp parallel for num_threads(threads_count)
    for(int i=param->PML;i<param->NX+param->PML; i++)
    {
        for(int j=param->PML; j<param->PML+param->NZ; j++)
        {
            Txxx_aft[i][j] = Txxx_now[i][j]
                +(lamda[i][j]+2.0*miu[i][j])*(param->dt*1.0/param->dx)*(
                        param->C1*(Vx_now[i+1][j]-Vx_now[i-0][j])+
                        param->C2*(Vx_now[i+2][j]-Vx_now[i-1][j])+
                        param->C3*(Vx_now[i+3][j]-Vx_now[i-2][j])+
                        param->C4*(Vx_now[i+4][j]-Vx_now[i-3][j])+
                        param->C5*(Vx_now[i+5][j]-Vx_now[i-4][j])+
                        param->C6*(Vx_now[i+6][j]-Vx_now[i-5][j])+
                        param->C7*(Vx_now[i+7][j]-Vx_now[i-6][j])+
                        param->C8*(Vx_now[i+8][j]-Vx_now[i-7][j]));

            Txxz_aft[i][j] = Txxz_now[i][j]
                +lamda[i][j]*(param->dt*1.0/param->dz)*(
                        param->C1*(Vz_now[i][j+0]-Vz_now[i][j-1])+
                        param->C2*(Vz_now[i][j+1]-Vz_now[i][j-2])+
                        param->C3*(Vz_now[i][j+2]-Vz_now[i][j-3])+
                        param->C4*(Vz_now[i][j+3]-Vz_now[i][j-4])+
                        param->C5*(Vz_now[i][j+4]-Vz_now[i][j-5])+
                        param->C6*(Vz_now[i][j+5]-Vz_now[i][j-6])+
                        param->C7*(Vz_now[i][j+6]-Vz_now[i][j-7])+
                        param->C8*(Vz_now[i][j+7]-Vz_now[i][j-8]));
        }
    }

    //Tzz
# pragma omp parallel for num_threads(threads_count)
    for(int i=param->PML;i<param->NX+param->PML; i++)
    {
        for(int j=param->PML; j<param->PML+param->NZ; j++)
        {
            Tzzx_aft[i][j] = Tzzx_now[i][j]
                +lamda[i][j]*(param->dt*1.0/param->dx)*(
                        param->C1*(Vx_now[i+1][j]-Vx_now[i-0][j])+
                        param->C2*(Vx_now[i+2][j]-Vx_now[i-1][j])+
                        param->C3*(Vx_now[i+3][j]-Vx_now[i-2][j])+
                        param->C4*(Vx_now[i+4][j]-Vx_now[i-3][j])+
                        param->C5*(Vx_now[i+5][j]-Vx_now[i-4][j])+
                        param->C6*(Vx_now[i+6][j]-Vx_now[i-5][j])+
                        param->C7*(Vx_now[i+7][j]-Vx_now[i-6][j])+
                        param->C8*(Vx_now[i+8][j]-Vx_now[i-7][j]));

            Tzzz_aft[i][j] = Tzzz_now[i][j]
                +(lamda[i][j]+2.0*miu[i][j])*(param->dt*1.0/param->dz)*(
                        param->C1*(Vz_now[i][j+0]-Vz_now[i][j-1])+
                        param->C2*(Vz_now[i][j+1]-Vz_now[i][j-2])+
                        param->C3*(Vz_now[i][j+2]-Vz_now[i][j-3])+
                        param->C4*(Vz_now[i][j+3]-Vz_now[i][j-4])+
                        param->C5*(Vz_now[i][j+4]-Vz_now[i][j-5])+
                        param->C6*(Vz_now[i][j+5]-Vz_now[i][j-6])+
                        param->C7*(Vz_now[i][j+6]-Vz_now[i][j-7])+
                        param->C8*(Vz_now[i][j+7]-Vz_now[i][j-8]));
        }
    }
    //Txz
# pragma omp parallel for num_threads(threads_count)
    for(int i=param->PML;i<param->NX+param->PML; i++)
    {
        for(int j=param->PML; j<param->PML+param->NZ; j++)
        {
            Txzx_aft[i][j] = Txzx_now[i][j]
                +miu[i][j]*(param->dt*1.0/param->dx)*(
                        param->C1*(Vz_now[i+0][j]-Vz_now[i-1][j])+
                        param->C2*(Vz_now[i+1][j]-Vz_now[i-2][j])+
                        param->C3*(Vz_now[i+2][j]-Vz_now[i-3][j])+
                        param->C4*(Vz_now[i+3][j]-Vz_now[i-4][j])+
                        param->C5*(Vz_now[i+4][j]-Vz_now[i-5][j])+
                        param->C6*(Vz_now[i+5][j]-Vz_now[i-6][j])+
                        param->C7*(Vz_now[i+6][j]-Vz_now[i-7][j])+
                        param->C8*(Vz_now[i+7][j]-Vz_now[i-8][j]));

            Txzz_aft[i][j] = Txzz_now[i][j]
                +miu[i][j]*(param->dt*1.0/param->dz)*(
                        param->C1*(Vx_now[i][j+1]-Vx_now[i][j-0])+
                        param->C2*(Vx_now[i][j+2]-Vx_now[i][j-1])+
                        param->C3*(Vx_now[i][j+3]-Vx_now[i][j-2])+
                        param->C4*(Vx_now[i][j+4]-Vx_now[i][j-3])+
                        param->C5*(Vx_now[i][j+5]-Vx_now[i][j-4])+
                        param->C6*(Vx_now[i][j+6]-Vx_now[i][j-5])+
                        param->C7*(Vx_now[i][j+7]-Vx_now[i][j-6])+
                        param->C8*(Vx_now[i][j+8]-Vx_now[i][j-7]));
        }
    }

    //left
    //Txx
# pragma omp parallel for num_threads(threads_count)
    for(int i=8; i<param->PML; i++)
    {
        for(int j=param->PML; j<param->NZ+ param->PML; j++)
        {
            Txxx_aft[i][j] = 1.0/(1.0+0.5*param->dt*absorbx[i][j])*(
                    (1.0-0.5*param->dt*absorbx[i][j])*Txxx_now[i][j]
                    +(lamda[i][j]+2.0*miu[i][j])*(param->dt*1.0/param->dx)*(
                        param->C1*(Vx_now[i+1][j]-Vx_now[i-0][j])+
                        param->C2*(Vx_now[i+2][j]-Vx_now[i-1][j])+
                        param->C3*(Vx_now[i+3][j]-Vx_now[i-2][j])+
                        param->C4*(Vx_now[i+4][j]-Vx_now[i-3][j])+
                        param->C5*(Vx_now[i+5][j]-Vx_now[i-4][j])+
                        param->C6*(Vx_now[i+6][j]-Vx_now[i-5][j])+
                        param->C7*(Vx_now[i+7][j]-Vx_now[i-6][j])+
                        param->C8*(Vx_now[i+8][j]-Vx_now[i-7][j])));

            Txxz_aft[i][j] = Txxz_now[i][j]
                +lamda[i][j]*(param->dt*1.0/param->dz)*(
                        param->C1*(Vz_now[i][j+0]-Vz_now[i][j-1])+
                        param->C2*(Vz_now[i][j+1]-Vz_now[i][j-2])+
                        param->C3*(Vz_now[i][j+2]-Vz_now[i][j-3])+
                        param->C4*(Vz_now[i][j+3]-Vz_now[i][j-4])+
                        param->C5*(Vz_now[i][j+4]-Vz_now[i][j-5])+
                        param->C6*(Vz_now[i][j+5]-Vz_now[i][j-6])+
                        param->C7*(Vz_now[i][j+6]-Vz_now[i][j-7])+
                        param->C8*(Vz_now[i][j+7]-Vz_now[i][j-8]));
        }
    }
    //Tzz
# pragma omp parallel for num_threads(threads_count)
    for(int i=8; i<param->PML; i++)
    {
        for(int j=param->PML; j<param->NZ+ param->PML; j++)
        {
            Tzzx_aft[i][j] = 1.0/(1.0+0.5*param->dt*absorbx[i][j])*(
                    (1.0-0.5*param->dt*absorbx[i][j])*Tzzx_now[i][j]
                    +lamda[i][j]*(param->dt*1.0/param->dx)*(
                        param->C1*(Vx_now[i+1][j]-Vx_now[i-0][j])+
                        param->C2*(Vx_now[i+2][j]-Vx_now[i-1][j])+
                        param->C3*(Vx_now[i+3][j]-Vx_now[i-2][j])+
                        param->C4*(Vx_now[i+4][j]-Vx_now[i-3][j])+
                        param->C5*(Vx_now[i+5][j]-Vx_now[i-4][j])+
                        param->C6*(Vx_now[i+6][j]-Vx_now[i-5][j])+
                        param->C7*(Vx_now[i+7][j]-Vx_now[i-6][j])+
                        param->C8*(Vx_now[i+8][j]-Vx_now[i-7][j])));

            Tzzz_aft[i][j] = Tzzz_now[i][j]
                +(lamda[i][j]+2.0*miu[i][j])*(param->dt*1.0/param->dz)*(
                        param->C1*(Vz_now[i][j+0]-Vz_now[i][j-1])+
                        param->C2*(Vz_now[i][j+1]-Vz_now[i][j-2])+
                        param->C3*(Vz_now[i][j+2]-Vz_now[i][j-3])+
                        param->C4*(Vz_now[i][j+3]-Vz_now[i][j-4])+
                        param->C5*(Vz_now[i][j+4]-Vz_now[i][j-5])+
                        param->C6*(Vz_now[i][j+5]-Vz_now[i][j-6])+
                        param->C7*(Vz_now[i][j+6]-Vz_now[i][j-7])+
                        param->C8*(Vz_now[i][j+7]-Vz_now[i][j-8]));
        }
    }
    //Txz
# pragma omp parallel for num_threads(threads_count)
    for(int i=8; i<param->PML; i++)
    {
        for(int j=param->PML; j<param->NZ+ param->PML; j++)
        {
            Txzx_aft[i][j] = 1.0/(1.0+0.5*param->dt*absorbx[i][j])*(
                    (1.0-0.5*param->dt*absorbx[i][j])*Txzx_now[i][j]
                    +miu[i][j]*(param->dt*1.0/param->dx)*(
                        param->C1*(Vz_now[i+0][j]-Vz_now[i-1][j])+
                        param->C2*(Vz_now[i+1][j]-Vz_now[i-2][j])+
                        param->C3*(Vz_now[i+2][j]-Vz_now[i-3][j])+
                        param->C4*(Vz_now[i+3][j]-Vz_now[i-4][j])+
                        param->C5*(Vz_now[i+4][j]-Vz_now[i-5][j])+
                        param->C6*(Vz_now[i+5][j]-Vz_now[i-6][j])+
                        param->C7*(Vz_now[i+6][j]-Vz_now[i-7][j])+
                        param->C8*(Vz_now[i+7][j]-Vz_now[i-8][j])));

            Txzz_aft[i][j] = Txzz_now[i][j]
                +miu[i][j]*(param->dt*1.0/param->dz)*(
                        param->C1*(Vx_now[i][j+1]-Vx_now[i][j-0])+
                        param->C2*(Vx_now[i][j+2]-Vx_now[i][j-1])+
                        param->C3*(Vx_now[i][j+3]-Vx_now[i][j-2])+
                        param->C4*(Vx_now[i][j+4]-Vx_now[i][j-3])+
                        param->C5*(Vx_now[i][j+5]-Vx_now[i][j-4])+
                        param->C6*(Vx_now[i][j+6]-Vx_now[i][j-5])+
                        param->C7*(Vx_now[i][j+7]-Vx_now[i][j-6])+
                        param->C8*(Vx_now[i][j+8]-Vx_now[i][j-7]));
        }
    }

    //right
    //Txx
# pragma omp parallel for num_threads(threads_count)
    for(int i=param->NX+param->PML;i<param->Nx-8; i++)
    {
        for(int j=param->PML; j<param->NZ+ param->PML; j++)
        {
            Txxx_aft[i][j] = 1.0/(1.0+0.5*param->dt*absorbx[i][j])*(
                    (1.0-0.5*param->dt*absorbx[i][j])*Txxx_now[i][j]
                    +(lamda[i][j]+2.0*miu[i][j])*(param->dt*1.0/param->dx)*(
                        param->C1*(Vx_now[i+1][j]-Vx_now[i-0][j])+
                        param->C2*(Vx_now[i+2][j]-Vx_now[i-1][j])+
                        param->C3*(Vx_now[i+3][j]-Vx_now[i-2][j])+
                        param->C4*(Vx_now[i+4][j]-Vx_now[i-3][j])+
                        param->C5*(Vx_now[i+5][j]-Vx_now[i-4][j])+
                        param->C6*(Vx_now[i+6][j]-Vx_now[i-5][j])+
                        param->C7*(Vx_now[i+7][j]-Vx_now[i-6][j])+
                        param->C8*(Vx_now[i+8][j]-Vx_now[i-7][j])));

            Txxz_aft[i][j] = Txxz_now[i][j]
                +lamda[i][j]*(param->dt*1.0/param->dz)*(
                        param->C1*(Vz_now[i][j+0]-Vz_now[i][j-1])+
                        param->C2*(Vz_now[i][j+1]-Vz_now[i][j-2])+
                        param->C3*(Vz_now[i][j+2]-Vz_now[i][j-3])+
                        param->C4*(Vz_now[i][j+3]-Vz_now[i][j-4])+
                        param->C5*(Vz_now[i][j+4]-Vz_now[i][j-5])+
                        param->C6*(Vz_now[i][j+5]-Vz_now[i][j-6])+
                        param->C7*(Vz_now[i][j+6]-Vz_now[i][j-7])+
                        param->C8*(Vz_now[i][j+7]-Vz_now[i][j-8]));
        }
    }
    //Tzz
# pragma omp parallel for num_threads(threads_count)
    for(int i=param->NX+param->PML;i<param->Nx-8; i++)
    {
        for(int j=param->PML; j<param->NZ+ param->PML; j++)
        {
            Tzzx_aft[i][j] = 1.0/(1.0+0.5*param->dt*absorbx[i][j])*(
                    (1.0-0.5*param->dt*absorbx[i][j])*Tzzx_now[i][j]
                    +lamda[i][j]*(param->dt*1.0/param->dx)*(
                        param->C1*(Vx_now[i+1][j]-Vx_now[i-0][j])+
                        param->C2*(Vx_now[i+2][j]-Vx_now[i-1][j])+
                        param->C3*(Vx_now[i+3][j]-Vx_now[i-2][j])+
                        param->C4*(Vx_now[i+4][j]-Vx_now[i-3][j])+
                        param->C5*(Vx_now[i+5][j]-Vx_now[i-4][j])+
                        param->C6*(Vx_now[i+6][j]-Vx_now[i-5][j])+
                        param->C7*(Vx_now[i+7][j]-Vx_now[i-6][j])+
                        param->C8*(Vx_now[i+8][j]-Vx_now[i-7][j])));

            Tzzz_aft[i][j] = Tzzz_now[i][j]
                +(lamda[i][j]+2.0*miu[i][j])*(param->dt*1.0/param->dz)*(
                        param->C1*(Vz_now[i][j+0]-Vz_now[i][j-1])+
                        param->C2*(Vz_now[i][j+1]-Vz_now[i][j-2])+
                        param->C3*(Vz_now[i][j+2]-Vz_now[i][j-3])+
                        param->C4*(Vz_now[i][j+3]-Vz_now[i][j-4])+
                        param->C5*(Vz_now[i][j+4]-Vz_now[i][j-5])+
                        param->C6*(Vz_now[i][j+5]-Vz_now[i][j-6])+
                        param->C7*(Vz_now[i][j+6]-Vz_now[i][j-7])+
                        param->C8*(Vz_now[i][j+7]-Vz_now[i][j-8]));
        }
    }
    //Txz
# pragma omp parallel for num_threads(threads_count)
    for(int i=param->NX+param->PML;i<param->Nx-8; i++)
    {
        for(int j=param->PML; j<param->NZ+ param->PML; j++)
        {
            Txzx_aft[i][j] = 1.0/(1.0+0.5*param->dt*absorbx[i][j])*(
                    (1.0-0.5*param->dt*absorbx[i][j])*Txzx_now[i][j]
                    +miu[i][j]*(param->dt*1.0/param->dx)*(
                        param->C1*(Vz_now[i+0][j]-Vz_now[i-1][j])+
                        param->C2*(Vz_now[i+1][j]-Vz_now[i-2][j])+
                        param->C3*(Vz_now[i+2][j]-Vz_now[i-3][j])+
                        param->C4*(Vz_now[i+3][j]-Vz_now[i-4][j])+
                        param->C5*(Vz_now[i+4][j]-Vz_now[i-5][j])+
                        param->C6*(Vz_now[i+5][j]-Vz_now[i-6][j])+
                        param->C7*(Vz_now[i+6][j]-Vz_now[i-7][j])+
                        param->C8*(Vz_now[i+7][j]-Vz_now[i-8][j])));

            Txzz_aft[i][j] = Txzz_now[i][j]
                +miu[i][j]*(param->dt*1.0/param->dz)*(
                        param->C1*(Vx_now[i][j+1]-Vx_now[i][j-0])+
                        param->C2*(Vx_now[i][j+2]-Vx_now[i][j-1])+
                        param->C3*(Vx_now[i][j+3]-Vx_now[i][j-2])+
                        param->C4*(Vx_now[i][j+4]-Vx_now[i][j-3])+
                        param->C5*(Vx_now[i][j+5]-Vx_now[i][j-4])+
                        param->C6*(Vx_now[i][j+6]-Vx_now[i][j-5])+
                        param->C7*(Vx_now[i][j+7]-Vx_now[i][j-6])+
                        param->C8*(Vx_now[i][j+8]-Vx_now[i][j-7]));
        }
    }
    //up
# pragma omp parallel for num_threads(threads_count)
    for(int i=param->PML;i<param->NX + param->PML; i++)
    {
        for(int j=8; j<param->PML; j++)
        {
            Txxx_aft[i][j] = Txxx_now[i][j]
                +(lamda[i][j]+2.0*miu[i][j])*(param->dt*1.0/param->dx)*(
                        param->C1*(Vx_now[i+1][j]-Vx_now[i-0][j])+
                        param->C2*(Vx_now[i+2][j]-Vx_now[i-1][j])+
                        param->C3*(Vx_now[i+3][j]-Vx_now[i-2][j])+
                        param->C4*(Vx_now[i+4][j]-Vx_now[i-3][j])+
                        param->C5*(Vx_now[i+5][j]-Vx_now[i-4][j])+
                        param->C6*(Vx_now[i+6][j]-Vx_now[i-5][j])+
                        param->C7*(Vx_now[i+7][j]-Vx_now[i-6][j])+
                        param->C8*(Vx_now[i+8][j]-Vx_now[i-7][j]));

            Txxz_aft[i][j] = 1.0/(1.0+0.5*param->dt*absorbz[i][j])*(
                    (1.0-0.5*param->dt*absorbz[i][j])*Txxz_now[i][j]
                    +lamda[i][j]*(param->dt*1.0/param->dz)*(
                        param->C1*(Vz_now[i][j+0]-Vz_now[i][j-1])+
                        param->C2*(Vz_now[i][j+1]-Vz_now[i][j-2])+
                        param->C3*(Vz_now[i][j+2]-Vz_now[i][j-3])+
                        param->C4*(Vz_now[i][j+3]-Vz_now[i][j-4])+
                        param->C5*(Vz_now[i][j+4]-Vz_now[i][j-5])+
                        param->C6*(Vz_now[i][j+5]-Vz_now[i][j-6])+
                        param->C7*(Vz_now[i][j+6]-Vz_now[i][j-7])+
                        param->C8*(Vz_now[i][j+7]-Vz_now[i][j-8])));
        }
    }
    //Tzz
# pragma omp parallel for num_threads(threads_count)
    for(int i=param->PML;i<param->NX + param->PML; i++)
    {
        for(int j=8; j<param->PML; j++)
        {
            Tzzx_aft[i][j] = Tzzx_now[i][j]
                +lamda[i][j]*(param->dt*1.0/param->dx)*(
                        param->C1*(Vx_now[i+1][j]-Vx_now[i-0][j])+
                        param->C2*(Vx_now[i+2][j]-Vx_now[i-1][j])+
                        param->C3*(Vx_now[i+3][j]-Vx_now[i-2][j])+
                        param->C4*(Vx_now[i+4][j]-Vx_now[i-3][j])+
                        param->C5*(Vx_now[i+5][j]-Vx_now[i-4][j])+
                        param->C6*(Vx_now[i+6][j]-Vx_now[i-5][j])+
                        param->C7*(Vx_now[i+7][j]-Vx_now[i-6][j])+
                        param->C8*(Vx_now[i+8][j]-Vx_now[i-7][j]));

            Tzzz_aft[i][j] = 1.0/(1.0+0.5*param->dt*absorbz[i][j])*(
                    (1.0-0.5*param->dt*absorbz[i][j])*Tzzz_now[i][j]
                    +(lamda[i][j]+2.0*miu[i][j])*(param->dt*1.0/param->dz)*(
                        param->C1*(Vz_now[i][j+0]-Vz_now[i][j-1])+
                        param->C2*(Vz_now[i][j+1]-Vz_now[i][j-2])+
                        param->C3*(Vz_now[i][j+2]-Vz_now[i][j-3])+
                        param->C4*(Vz_now[i][j+3]-Vz_now[i][j-4])+
                        param->C5*(Vz_now[i][j+4]-Vz_now[i][j-5])+
                        param->C6*(Vz_now[i][j+5]-Vz_now[i][j-6])+
                        param->C7*(Vz_now[i][j+6]-Vz_now[i][j-7])+
                        param->C8*(Vz_now[i][j+7]-Vz_now[i][j-8])));
        }
    }
    //Txz
# pragma omp parallel for num_threads(threads_count)
    for(int i=param->PML;i<param->NX + param->PML; i++)
    {
        for(int j=8; j<param->PML; j++)
        {
            Txzx_aft[i][j] = Txzx_now[i][j]
                +miu[i][j]*(param->dt*1.0/param->dx)*(
                        param->C1*(Vz_now[i+0][j]-Vz_now[i-1][j])+
                        param->C2*(Vz_now[i+1][j]-Vz_now[i-2][j])+
                        param->C3*(Vz_now[i+2][j]-Vz_now[i-3][j])+
                        param->C4*(Vz_now[i+3][j]-Vz_now[i-4][j])+
                        param->C5*(Vz_now[i+4][j]-Vz_now[i-5][j])+
                        param->C6*(Vz_now[i+5][j]-Vz_now[i-6][j])+
                        param->C7*(Vz_now[i+6][j]-Vz_now[i-7][j])+
                        param->C8*(Vz_now[i+7][j]-Vz_now[i-8][j]));

            Txzz_aft[i][j] = 1.0/(1.0+0.5*param->dt*absorbz[i][j])*(
                    (1.0-0.5*param->dt*absorbz[i][j])*Txzz_now[i][j]
                    +miu[i][j]*(param->dt*1.0/param->dz)*(
                        param->C1*(Vx_now[i][j+1]-Vx_now[i][j-0])+
                        param->C2*(Vx_now[i][j+2]-Vx_now[i][j-1])+
                        param->C3*(Vx_now[i][j+3]-Vx_now[i][j-2])+
                        param->C4*(Vx_now[i][j+4]-Vx_now[i][j-3])+
                        param->C5*(Vx_now[i][j+5]-Vx_now[i][j-4])+
                        param->C6*(Vx_now[i][j+6]-Vx_now[i][j-5])+
                        param->C7*(Vx_now[i][j+7]-Vx_now[i][j-6])+
                        param->C8*(Vx_now[i][j+8]-Vx_now[i][j-7])));
        }
    }

    //under
# pragma omp parallel for num_threads(threads_count)
    for(int i=param->PML;i<param->NX + param->PML; i++)
    {
        for(int j=param->NZ+param->PML;  j<param->Nz-8;  j++)
        {
            Txxx_aft[i][j] = Txxx_now[i][j]
                +(lamda[i][j]+2.0*miu[i][j])*(param->dt*1.0/param->dx)*(
                        param->C1*(Vx_now[i+1][j]-Vx_now[i-0][j])+
                        param->C2*(Vx_now[i+2][j]-Vx_now[i-1][j])+
                        param->C3*(Vx_now[i+3][j]-Vx_now[i-2][j])+
                        param->C4*(Vx_now[i+4][j]-Vx_now[i-3][j])+
                        param->C5*(Vx_now[i+5][j]-Vx_now[i-4][j])+
                        param->C6*(Vx_now[i+6][j]-Vx_now[i-5][j])+
                        param->C7*(Vx_now[i+7][j]-Vx_now[i-6][j])+
                        param->C8*(Vx_now[i+8][j]-Vx_now[i-7][j]));

            Txxz_aft[i][j] = 1.0/(1.0+0.5*param->dt*absorbz[i][j])*(
                    (1.0-0.5*param->dt*absorbz[i][j])*Txxz_now[i][j]
                    +lamda[i][j]*(param->dt*1.0/param->dz)*(
                        param->C1*(Vz_now[i][j+0]-Vz_now[i][j-1])+
                        param->C2*(Vz_now[i][j+1]-Vz_now[i][j-2])+
                        param->C3*(Vz_now[i][j+2]-Vz_now[i][j-3])+
                        param->C4*(Vz_now[i][j+3]-Vz_now[i][j-4])+
                        param->C5*(Vz_now[i][j+4]-Vz_now[i][j-5])+
                        param->C6*(Vz_now[i][j+5]-Vz_now[i][j-6])+
                        param->C7*(Vz_now[i][j+6]-Vz_now[i][j-7])+
                        param->C8*(Vz_now[i][j+7]-Vz_now[i][j-8])));
        }
    }
    //Tzz
# pragma omp parallel for num_threads(threads_count)
    for(int i=param->PML;i<param->NX + param->PML; i++)
    {
        for(int j=param->NZ+param->PML;  j<param->Nz-8;  j++)
        {
            Tzzx_aft[i][j] = Tzzx_now[i][j]
                +lamda[i][j]*(param->dt*1.0/param->dx)*(
                        param->C1*(Vx_now[i+1][j]-Vx_now[i-0][j])+
                        param->C2*(Vx_now[i+2][j]-Vx_now[i-1][j])+
                        param->C3*(Vx_now[i+3][j]-Vx_now[i-2][j])+
                        param->C4*(Vx_now[i+4][j]-Vx_now[i-3][j])+
                        param->C5*(Vx_now[i+5][j]-Vx_now[i-4][j])+
                        param->C6*(Vx_now[i+6][j]-Vx_now[i-5][j])+
                        param->C7*(Vx_now[i+7][j]-Vx_now[i-6][j])+
                        param->C8*(Vx_now[i+8][j]-Vx_now[i-7][j]));

            Tzzz_aft[i][j] = 1.0/(1.0+0.5*param->dt*absorbz[i][j])*(
                    (1.0-0.5*param->dt*absorbz[i][j])*Tzzz_now[i][j]
                    +(lamda[i][j]+2.0*miu[i][j])*(param->dt*1.0/param->dz)*(
                        param->C1*(Vz_now[i][j+0]-Vz_now[i][j-1])+
                        param->C2*(Vz_now[i][j+1]-Vz_now[i][j-2])+
                        param->C3*(Vz_now[i][j+2]-Vz_now[i][j-3])+
                        param->C4*(Vz_now[i][j+3]-Vz_now[i][j-4])+
                        param->C5*(Vz_now[i][j+4]-Vz_now[i][j-5])+
                        param->C6*(Vz_now[i][j+5]-Vz_now[i][j-6])+
                        param->C7*(Vz_now[i][j+6]-Vz_now[i][j-7])+
                        param->C8*(Vz_now[i][j+7]-Vz_now[i][j-8])));
        }
    }
    //Txz
# pragma omp parallel for num_threads(threads_count)
    for(int i=param->PML;i<param->NX + param->PML; i++)
    {
        for(int j=param->NZ+param->PML;  j<param->Nz-8;  j++)
        {
            Txzx_aft[i][j] = Txzx_now[i][j]
                +miu[i][j]*(param->dt*1.0/param->dx)*(
                        param->C1*(Vz_now[i+0][j]-Vz_now[i-1][j])+
                        param->C2*(Vz_now[i+1][j]-Vz_now[i-2][j])+
                        param->C3*(Vz_now[i+2][j]-Vz_now[i-3][j])+
                        param->C4*(Vz_now[i+3][j]-Vz_now[i-4][j])+
                        param->C5*(Vz_now[i+4][j]-Vz_now[i-5][j])+
                        param->C6*(Vz_now[i+5][j]-Vz_now[i-6][j])+
                        param->C7*(Vz_now[i+6][j]-Vz_now[i-7][j])+
                        param->C8*(Vz_now[i+7][j]-Vz_now[i-8][j]));

            Txzz_aft[i][j] = 1.0/(1.0+0.5*param->dt*absorbz[i][j])*(
                    (1.0-0.5*param->dt*absorbz[i][j])*Txzz_now[i][j]
                    +miu[i][j]*(param->dt*1.0/param->dz)*(
                        param->C1*(Vx_now[i][j+1]-Vx_now[i][j-0])+
                        param->C2*(Vx_now[i][j+2]-Vx_now[i][j-1])+
                        param->C3*(Vx_now[i][j+3]-Vx_now[i][j-2])+
                        param->C4*(Vx_now[i][j+4]-Vx_now[i][j-3])+
                        param->C5*(Vx_now[i][j+5]-Vx_now[i][j-4])+
                        param->C6*(Vx_now[i][j+6]-Vx_now[i][j-5])+
                        param->C7*(Vx_now[i][j+7]-Vx_now[i][j-6])+
                        param->C8*(Vx_now[i][j+8]-Vx_now[i][j-7])));
        }
    }

    //coner-left-up
    //Txx
# pragma omp parallel for num_threads(threads_count)
    for(int i=8; i<param->PML; i++)
    {
        for(int j=8; j<param->PML; j++)
        {
            Txxx_aft[i][j] = 1.0/(1.0+0.5*param->dt*absorbx[i][j])*(
                    (1.0-0.5*param->dt*absorbx[i][j])*Txxx_now[i][j]
                    +(lamda[i][j]+2.0*miu[i][j])*(param->dt*1.0/param->dx)*(
                        param->C1*(Vx_now[i+1][j]-Vx_now[i-0][j])+
                        param->C2*(Vx_now[i+2][j]-Vx_now[i-1][j])+
                        param->C3*(Vx_now[i+3][j]-Vx_now[i-2][j])+
                        param->C4*(Vx_now[i+4][j]-Vx_now[i-3][j])+
                        param->C5*(Vx_now[i+5][j]-Vx_now[i-4][j])+
                        param->C6*(Vx_now[i+6][j]-Vx_now[i-5][j])+
                        param->C7*(Vx_now[i+7][j]-Vx_now[i-6][j])+
                        param->C8*(Vx_now[i+8][j]-Vx_now[i-7][j])));

            Txxz_aft[i][j] = 1.0/(1.0+0.5*param->dt*absorbz[i][j])*(
                    (1.0-0.5*param->dt*absorbz[i][j])*Txxz_now[i][j]
                    +lamda[i][j]*(param->dt*1.0/param->dz)*(
                        param->C1*(Vz_now[i][j+0]-Vz_now[i][j-1])+
                        param->C2*(Vz_now[i][j+1]-Vz_now[i][j-2])+
                        param->C3*(Vz_now[i][j+2]-Vz_now[i][j-3])+
                        param->C4*(Vz_now[i][j+3]-Vz_now[i][j-4])+
                        param->C5*(Vz_now[i][j+4]-Vz_now[i][j-5])+
                        param->C6*(Vz_now[i][j+5]-Vz_now[i][j-6])+
                        param->C7*(Vz_now[i][j+6]-Vz_now[i][j-7])+
                        param->C8*(Vz_now[i][j+7]-Vz_now[i][j-8])));
        }
    }
    //Tzz
# pragma omp parallel for num_threads(threads_count)
    for(int i=8; i<param->PML; i++)
    {
        for(int j=8; j<param->PML; j++)
        {
            Tzzx_aft[i][j] = 1.0/(1.0+0.5*param->dt*absorbx[i][j])*(
                    (1.0-0.5*param->dt*absorbx[i][j])*Tzzx_now[i][j]
                    +lamda[i][j]*(param->dt*1.0/param->dx)*(
                        param->C1*(Vx_now[i+1][j]-Vx_now[i-0][j])+
                        param->C2*(Vx_now[i+2][j]-Vx_now[i-1][j])+
                        param->C3*(Vx_now[i+3][j]-Vx_now[i-2][j])+
                        param->C4*(Vx_now[i+4][j]-Vx_now[i-3][j])+
                        param->C5*(Vx_now[i+5][j]-Vx_now[i-4][j])+
                        param->C6*(Vx_now[i+6][j]-Vx_now[i-5][j])+
                        param->C7*(Vx_now[i+7][j]-Vx_now[i-6][j])+
                        param->C8*(Vx_now[i+8][j]-Vx_now[i-7][j])));

            Tzzz_aft[i][j] = 1.0/(1.0+0.5*param->dt*absorbz[i][j])*(
                    (1.0-0.5*param->dt*absorbz[i][j])*Tzzz_now[i][j]
                    +(lamda[i][j]+2.0*miu[i][j])*(param->dt*1.0/param->dz)*(
                        param->C1*(Vz_now[i][j+0]-Vz_now[i][j-1])+
                        param->C2*(Vz_now[i][j+1]-Vz_now[i][j-2])+
                        param->C3*(Vz_now[i][j+2]-Vz_now[i][j-3])+
                        param->C4*(Vz_now[i][j+3]-Vz_now[i][j-4])+
                        param->C5*(Vz_now[i][j+4]-Vz_now[i][j-5])+
                        param->C6*(Vz_now[i][j+5]-Vz_now[i][j-6])+
                        param->C7*(Vz_now[i][j+6]-Vz_now[i][j-7])+
                        param->C8*(Vz_now[i][j+7]-Vz_now[i][j-8])));
        }
    }
    //Txz
# pragma omp parallel for num_threads(threads_count)
    for(int i=8; i<param->PML; i++)
    {
        for(int j=8; j<param->PML; j++)
        {
            Txzx_aft[i][j] = 1.0/(1.0+0.5*param->dt*absorbx[i][j])*(
                    (1.0-0.5*param->dt*absorbx[i][j])*Txzx_now[i][j]
                    +miu[i][j]*(param->dt*1.0/param->dx)*(
                        param->C1*(Vz_now[i+0][j]-Vz_now[i-1][j])+
                        param->C2*(Vz_now[i+1][j]-Vz_now[i-2][j])+
                        param->C3*(Vz_now[i+2][j]-Vz_now[i-3][j])+
                        param->C4*(Vz_now[i+3][j]-Vz_now[i-4][j])+
                        param->C5*(Vz_now[i+4][j]-Vz_now[i-5][j])+
                        param->C6*(Vz_now[i+5][j]-Vz_now[i-6][j])+
                        param->C7*(Vz_now[i+6][j]-Vz_now[i-7][j])+
                        param->C8*(Vz_now[i+7][j]-Vz_now[i-8][j])));

            Txzz_aft[i][j] = 1.0/(1.0+0.5*param->dt*absorbz[i][j])*(
                    (1.0-0.5*param->dt*absorbz[i][j])*Txzz_now[i][j]
                    +miu[i][j]*(param->dt*1.0/param->dz)*(
                        param->C1*(Vx_now[i][j+1]-Vx_now[i][j-0])+
                        param->C2*(Vx_now[i][j+2]-Vx_now[i][j-1])+
                        param->C3*(Vx_now[i][j+3]-Vx_now[i][j-2])+
                        param->C4*(Vx_now[i][j+4]-Vx_now[i][j-3])+
                        param->C5*(Vx_now[i][j+5]-Vx_now[i][j-4])+
                        param->C6*(Vx_now[i][j+6]-Vx_now[i][j-5])+
                        param->C7*(Vx_now[i][j+7]-Vx_now[i][j-6])+
                        param->C8*(Vx_now[i][j+8]-Vx_now[i][j-7])));
        }
    }

    //coner-right-up
    //Txx
# pragma omp parallel for num_threads(threads_count)
    for(int i=param->PML+param->NX; i<param->Nx-8; i++)
    {
        for(int j=8; j<param->PML; j++)
        {
            Txxx_aft[i][j] = 1.0/(1.0+0.5*param->dt*absorbx[i][j])*(
                    (1.0-0.5*param->dt*absorbx[i][j])*Txxx_now[i][j]
                    +(lamda[i][j]+2.0*miu[i][j])*(param->dt*1.0/param->dx)*(
                        param->C1*(Vx_now[i+1][j]-Vx_now[i-0][j])+
                        param->C2*(Vx_now[i+2][j]-Vx_now[i-1][j])+
                        param->C3*(Vx_now[i+3][j]-Vx_now[i-2][j])+
                        param->C4*(Vx_now[i+4][j]-Vx_now[i-3][j])+
                        param->C5*(Vx_now[i+5][j]-Vx_now[i-4][j])+
                        param->C6*(Vx_now[i+6][j]-Vx_now[i-5][j])+
                        param->C7*(Vx_now[i+7][j]-Vx_now[i-6][j])+
                        param->C8*(Vx_now[i+8][j]-Vx_now[i-7][j])));

            Txxz_aft[i][j] = 1.0/(1.0+0.5*param->dt*absorbz[i][j])*(
                    (1.0-0.5*param->dt*absorbz[i][j])*Txxz_now[i][j]
                    +lamda[i][j]*(param->dt*1.0/param->dz)*(
                        param->C1*(Vz_now[i][j+0]-Vz_now[i][j-1])+
                        param->C2*(Vz_now[i][j+1]-Vz_now[i][j-2])+
                        param->C3*(Vz_now[i][j+2]-Vz_now[i][j-3])+
                        param->C4*(Vz_now[i][j+3]-Vz_now[i][j-4])+
                        param->C5*(Vz_now[i][j+4]-Vz_now[i][j-5])+
                        param->C6*(Vz_now[i][j+5]-Vz_now[i][j-6])+
                        param->C7*(Vz_now[i][j+6]-Vz_now[i][j-7])+
                        param->C8*(Vz_now[i][j+7]-Vz_now[i][j-8])));

        }
    }
    //Tzz
# pragma omp parallel for num_threads(threads_count)
    for(int i=param->PML+param->NX; i<param->Nx-8; i++)
    {
        for(int j=8; j<param->PML; j++)
        {
            Tzzx_aft[i][j] = 1.0/(1.0+0.5*param->dt*absorbx[i][j])*(
                    (1.0-0.5*param->dt*absorbx[i][j])*Tzzx_now[i][j]
                    +lamda[i][j]*(param->dt*1.0/param->dx)*(
                        param->C1*(Vx_now[i+1][j]-Vx_now[i-0][j])+
                        param->C2*(Vx_now[i+2][j]-Vx_now[i-1][j])+
                        param->C3*(Vx_now[i+3][j]-Vx_now[i-2][j])+
                        param->C4*(Vx_now[i+4][j]-Vx_now[i-3][j])+
                        param->C5*(Vx_now[i+5][j]-Vx_now[i-4][j])+
                        param->C6*(Vx_now[i+6][j]-Vx_now[i-5][j])+
                        param->C7*(Vx_now[i+7][j]-Vx_now[i-6][j])+
                        param->C8*(Vx_now[i+8][j]-Vx_now[i-7][j])));

            Tzzz_aft[i][j] = 1.0/(1.0+0.5*param->dt*absorbz[i][j])*(
                    (1.0-0.5*param->dt*absorbz[i][j])*Tzzz_now[i][j]
                    +(lamda[i][j]+2.0*miu[i][j])*(param->dt*1.0/param->dz)*(
                        param->C1*(Vz_now[i][j+0]-Vz_now[i][j-1])+
                        param->C2*(Vz_now[i][j+1]-Vz_now[i][j-2])+
                        param->C3*(Vz_now[i][j+2]-Vz_now[i][j-3])+
                        param->C4*(Vz_now[i][j+3]-Vz_now[i][j-4])+
                        param->C5*(Vz_now[i][j+4]-Vz_now[i][j-5])+
                        param->C6*(Vz_now[i][j+5]-Vz_now[i][j-6])+
                        param->C7*(Vz_now[i][j+6]-Vz_now[i][j-7])+
                        param->C8*(Vz_now[i][j+7]-Vz_now[i][j-8])));
        }
    }
    //Txz
# pragma omp parallel for num_threads(threads_count)
    for(int i=param->PML+param->NX; i<param->Nx-8; i++)
    {
        for(int j=8; j<param->PML; j++)
        {
            Txzx_aft[i][j] = 1.0/(1.0+0.5*param->dt*absorbx[i][j])*(
                    (1.0-0.5*param->dt*absorbx[i][j])*Txzx_now[i][j]
                    +miu[i][j]*(param->dt*1.0/param->dx)*(
                        param->C1*(Vz_now[i+0][j]-Vz_now[i-1][j])+
                        param->C2*(Vz_now[i+1][j]-Vz_now[i-2][j])+
                        param->C3*(Vz_now[i+2][j]-Vz_now[i-3][j])+
                        param->C4*(Vz_now[i+3][j]-Vz_now[i-4][j])+
                        param->C5*(Vz_now[i+4][j]-Vz_now[i-5][j])+
                        param->C6*(Vz_now[i+5][j]-Vz_now[i-6][j])+
                        param->C7*(Vz_now[i+6][j]-Vz_now[i-7][j])+
                        param->C8*(Vz_now[i+7][j]-Vz_now[i-8][j])));

            Txzz_aft[i][j] = 1.0/(1.0+0.5*param->dt*absorbz[i][j])*(
                    (1.0-0.5*param->dt*absorbz[i][j])*Txzz_now[i][j]
                    +miu[i][j]*(param->dt*1.0/param->dz)*(
                        param->C1*(Vx_now[i][j+1]-Vx_now[i][j-0])+
                        param->C2*(Vx_now[i][j+2]-Vx_now[i][j-1])+
                        param->C3*(Vx_now[i][j+3]-Vx_now[i][j-2])+
                        param->C4*(Vx_now[i][j+4]-Vx_now[i][j-3])+
                        param->C5*(Vx_now[i][j+5]-Vx_now[i][j-4])+
                        param->C6*(Vx_now[i][j+6]-Vx_now[i][j-5])+
                        param->C7*(Vx_now[i][j+7]-Vx_now[i][j-6])+
                        param->C8*(Vx_now[i][j+8]-Vx_now[i][j-7])));
        }
    }

    //coner-right-under
    //Txx
# pragma omp parallel for num_threads(threads_count)
    for(int i=param->PML+param->NX; i<param->Nx-8; i++)
    {
        for(int j=param->PML+param->NZ; j<param->Nz-8; j++)
        {
            Txxx_aft[i][j] = 1.0/(1.0+0.5*param->dt*absorbx[i][j])*(
                    (1.0-0.5*param->dt*absorbx[i][j])*Txxx_now[i][j]
                    +(lamda[i][j]+2.0*miu[i][j])*(param->dt*1.0/param->dx)*(
                        param->C1*(Vx_now[i+1][j]-Vx_now[i-0][j])+
                        param->C2*(Vx_now[i+2][j]-Vx_now[i-1][j])+
                        param->C3*(Vx_now[i+3][j]-Vx_now[i-2][j])+
                        param->C4*(Vx_now[i+4][j]-Vx_now[i-3][j])+
                        param->C5*(Vx_now[i+5][j]-Vx_now[i-4][j])+
                        param->C6*(Vx_now[i+6][j]-Vx_now[i-5][j])+
                        param->C7*(Vx_now[i+7][j]-Vx_now[i-6][j])+
                        param->C8*(Vx_now[i+8][j]-Vx_now[i-7][j])));

            Txxz_aft[i][j] = 1.0/(1.0+0.5*param->dt*absorbz[i][j])*(
                    (1.0-0.5*param->dt*absorbz[i][j])*Txxz_now[i][j]
                    +lamda[i][j]*(param->dt*1.0/param->dz)*(
                        param->C1*(Vz_now[i][j+0]-Vz_now[i][j-1])+
                        param->C2*(Vz_now[i][j+1]-Vz_now[i][j-2])+
                        param->C3*(Vz_now[i][j+2]-Vz_now[i][j-3])+
                        param->C4*(Vz_now[i][j+3]-Vz_now[i][j-4])+
                        param->C5*(Vz_now[i][j+4]-Vz_now[i][j-5])+
                        param->C6*(Vz_now[i][j+5]-Vz_now[i][j-6])+
                        param->C7*(Vz_now[i][j+6]-Vz_now[i][j-7])+
                        param->C8*(Vz_now[i][j+7]-Vz_now[i][j-8])));

        }
    }
    //Tzz
# pragma omp parallel for num_threads(threads_count)
    for(int i=param->PML+param->NX; i<param->Nx-8; i++)
    {
        for(int j=param->PML+param->NZ; j<param->Nz-8; j++)
        {
            Tzzx_aft[i][j] = 1.0/(1.0+0.5*param->dt*absorbx[i][j])*(
                    (1.0-0.5*param->dt*absorbx[i][j])*Tzzx_now[i][j]
                    +lamda[i][j]*(param->dt*1.0/param->dx)*(
                        param->C1*(Vx_now[i+1][j]-Vx_now[i-0][j])+
                        param->C2*(Vx_now[i+2][j]-Vx_now[i-1][j])+
                        param->C3*(Vx_now[i+3][j]-Vx_now[i-2][j])+
                        param->C4*(Vx_now[i+4][j]-Vx_now[i-3][j])+
                        param->C5*(Vx_now[i+5][j]-Vx_now[i-4][j])+
                        param->C6*(Vx_now[i+6][j]-Vx_now[i-5][j])+
                        param->C7*(Vx_now[i+7][j]-Vx_now[i-6][j])+
                        param->C8*(Vx_now[i+8][j]-Vx_now[i-7][j])));

            Tzzz_aft[i][j] = 1.0/(1.0+0.5*param->dt*absorbz[i][j])*(
                    (1.0-0.5*param->dt*absorbz[i][j])*Tzzz_now[i][j]
                    +(lamda[i][j]+2.0*miu[i][j])*(param->dt*1.0/param->dz)*(
                        param->C1*(Vz_now[i][j+0]-Vz_now[i][j-1])+
                        param->C2*(Vz_now[i][j+1]-Vz_now[i][j-2])+
                        param->C3*(Vz_now[i][j+2]-Vz_now[i][j-3])+
                        param->C4*(Vz_now[i][j+3]-Vz_now[i][j-4])+
                        param->C5*(Vz_now[i][j+4]-Vz_now[i][j-5])+
                        param->C6*(Vz_now[i][j+5]-Vz_now[i][j-6])+
                        param->C7*(Vz_now[i][j+6]-Vz_now[i][j-7])+
                        param->C8*(Vz_now[i][j+7]-Vz_now[i][j-8])));
        }
    }
    //Txz
# pragma omp parallel for num_threads(threads_count)
    for(int i=param->PML+param->NX; i<param->Nx-8; i++)
    {
        for(int j=param->PML+param->NZ; j<param->Nz-8; j++)
        {
            Txzx_aft[i][j] = 1.0/(1.0+0.5*param->dt*absorbx[i][j])*(
                    (1.0-0.5*param->dt*absorbx[i][j])*Txzx_now[i][j]
                    +miu[i][j]*(param->dt*1.0/param->dx)*(
                        param->C1*(Vz_now[i+0][j]-Vz_now[i-1][j])+
                        param->C2*(Vz_now[i+1][j]-Vz_now[i-2][j])+
                        param->C3*(Vz_now[i+2][j]-Vz_now[i-3][j])+
                        param->C4*(Vz_now[i+3][j]-Vz_now[i-4][j])+
                        param->C5*(Vz_now[i+4][j]-Vz_now[i-5][j])+
                        param->C6*(Vz_now[i+5][j]-Vz_now[i-6][j])+
                        param->C7*(Vz_now[i+6][j]-Vz_now[i-7][j])+
                        param->C8*(Vz_now[i+7][j]-Vz_now[i-8][j])));

            Txzz_aft[i][j] = 1.0/(1.0+0.5*param->dt*absorbz[i][j])*(
                    (1.0-0.5*param->dt*absorbz[i][j])*Txzz_now[i][j]
                    +miu[i][j]*(param->dt*1.0/param->dz)*(
                        param->C1*(Vx_now[i][j+1]-Vx_now[i][j-0])+
                        param->C2*(Vx_now[i][j+2]-Vx_now[i][j-1])+
                        param->C3*(Vx_now[i][j+3]-Vx_now[i][j-2])+
                        param->C4*(Vx_now[i][j+4]-Vx_now[i][j-3])+
                        param->C5*(Vx_now[i][j+5]-Vx_now[i][j-4])+
                        param->C6*(Vx_now[i][j+6]-Vx_now[i][j-5])+
                        param->C7*(Vx_now[i][j+7]-Vx_now[i][j-6])+
                        param->C8*(Vx_now[i][j+8]-Vx_now[i][j-7])));
        }
    }

    //coner-left-under
    //Txx
# pragma omp parallel for num_threads(threads_count)
    for(int i=8; i<param->PML; i++)
    {
        for(int j=param->PML+param->NZ; j<param->Nz-8; j++)
        {
            Txxx_aft[i][j] = 1.0/(1.0+0.5*param->dt*absorbx[i][j])*(
                    (1.0-0.5*param->dt*absorbx[i][j])*Txxx_now[i][j]
                    +(lamda[i][j]+2.0*miu[i][j])*(param->dt*1.0/param->dx)*(
                        param->C1*(Vx_now[i+1][j]-Vx_now[i-0][j])+
                        param->C2*(Vx_now[i+2][j]-Vx_now[i-1][j])+
                        param->C3*(Vx_now[i+3][j]-Vx_now[i-2][j])+
                        param->C4*(Vx_now[i+4][j]-Vx_now[i-3][j])+
                        param->C5*(Vx_now[i+5][j]-Vx_now[i-4][j])+
                        param->C6*(Vx_now[i+6][j]-Vx_now[i-5][j])+
                        param->C7*(Vx_now[i+7][j]-Vx_now[i-6][j])+
                        param->C8*(Vx_now[i+8][j]-Vx_now[i-7][j])));

            Txxz_aft[i][j] = 1.0/(1.0+0.5*param->dt*absorbz[i][j])*(
                    (1.0-0.5*param->dt*absorbz[i][j])*Txxz_now[i][j]
                    +lamda[i][j]*(param->dt*1.0/param->dz)*(
                        param->C1*(Vz_now[i][j+0]-Vz_now[i][j-1])+
                        param->C2*(Vz_now[i][j+1]-Vz_now[i][j-2])+
                        param->C3*(Vz_now[i][j+2]-Vz_now[i][j-3])+
                        param->C4*(Vz_now[i][j+3]-Vz_now[i][j-4])+
                        param->C5*(Vz_now[i][j+4]-Vz_now[i][j-5])+
                        param->C6*(Vz_now[i][j+5]-Vz_now[i][j-6])+
                        param->C7*(Vz_now[i][j+6]-Vz_now[i][j-7])+
                        param->C8*(Vz_now[i][j+7]-Vz_now[i][j-8])));
        }
    }
    //Tzz
# pragma omp parallel for num_threads(threads_count)
    for(int i=8; i<param->PML; i++)
    {
        for(int j=param->PML+param->NZ; j<param->Nz-8; j++)
        {
            Tzzx_aft[i][j] = 1.0/(1.0+0.5*param->dt*absorbx[i][j])*(
                    (1.0-0.5*param->dt*absorbx[i][j])*Tzzx_now[i][j]
                    +lamda[i][j]*(param->dt*1.0/param->dx)*(
                        param->C1*(Vx_now[i+1][j]-Vx_now[i-0][j])+
                        param->C2*(Vx_now[i+2][j]-Vx_now[i-1][j])+
                        param->C3*(Vx_now[i+3][j]-Vx_now[i-2][j])+
                        param->C4*(Vx_now[i+4][j]-Vx_now[i-3][j])+
                        param->C5*(Vx_now[i+5][j]-Vx_now[i-4][j])+
                        param->C6*(Vx_now[i+6][j]-Vx_now[i-5][j])+
                        param->C7*(Vx_now[i+7][j]-Vx_now[i-6][j])+
                        param->C8*(Vx_now[i+8][j]-Vx_now[i-7][j])));

            Tzzz_aft[i][j] = 1.0/(1.0+0.5*param->dt*absorbz[i][j])*(
                    (1.0-0.5*param->dt*absorbz[i][j])*Tzzz_now[i][j]
                    +(lamda[i][j]+2.0*miu[i][j])*(param->dt*1.0/param->dz)*(
                        param->C1*(Vz_now[i][j+0]-Vz_now[i][j-1])+
                        param->C2*(Vz_now[i][j+1]-Vz_now[i][j-2])+
                        param->C3*(Vz_now[i][j+2]-Vz_now[i][j-3])+
                        param->C4*(Vz_now[i][j+3]-Vz_now[i][j-4])+
                        param->C5*(Vz_now[i][j+4]-Vz_now[i][j-5])+
                        param->C6*(Vz_now[i][j+5]-Vz_now[i][j-6])+
                        param->C7*(Vz_now[i][j+6]-Vz_now[i][j-7])+
                        param->C8*(Vz_now[i][j+7]-Vz_now[i][j-8])));
        }
    }
    //Txz
# pragma omp parallel for num_threads(threads_count)
    for(int i=8; i<param->PML; i++)
    {
        for(int j=param->PML+param->NZ; j<param->Nz-8; j++)
        {
            Txzx_aft[i][j] = 1.0/(1.0+0.5*param->dt*absorbx[i][j])*(
                    (1.0-0.5*param->dt*absorbx[i][j])*Txzx_now[i][j]
                    +miu[i][j]*(param->dt*1.0/param->dx)*(
                        param->C1*(Vz_now[i+0][j]-Vz_now[i-1][j])+
                        param->C2*(Vz_now[i+1][j]-Vz_now[i-2][j])+
                        param->C3*(Vz_now[i+2][j]-Vz_now[i-3][j])+
                        param->C4*(Vz_now[i+3][j]-Vz_now[i-4][j])+
                        param->C5*(Vz_now[i+4][j]-Vz_now[i-5][j])+
                        param->C6*(Vz_now[i+5][j]-Vz_now[i-6][j])+
                        param->C7*(Vz_now[i+6][j]-Vz_now[i-7][j])+
                        param->C8*(Vz_now[i+7][j]-Vz_now[i-8][j])));

            Txzz_aft[i][j] = 1.0/(1.0+0.5*param->dt*absorbz[i][j])*(
                    (1.0-0.5*param->dt*absorbz[i][j])*Txzz_now[i][j]
                    +miu[i][j]*(param->dt*1.0/param->dz)*(
                        param->C1*(Vx_now[i][j+1]-Vx_now[i][j-0])+
                        param->C2*(Vx_now[i][j+2]-Vx_now[i][j-1])+
                        param->C3*(Vx_now[i][j+3]-Vx_now[i][j-2])+
                        param->C4*(Vx_now[i][j+4]-Vx_now[i][j-3])+
                        param->C5*(Vx_now[i][j+5]-Vx_now[i][j-4])+
                        param->C6*(Vx_now[i][j+6]-Vx_now[i][j-5])+
                        param->C7*(Vx_now[i][j+7]-Vx_now[i][j-6])+
                        param->C8*(Vx_now[i][j+8]-Vx_now[i][j-7])));
        }
    }

    //Txx、Tzz、Txz
# pragma omp parallel for num_threads(threads_count)
    for(int i=0;i<param->Nx; i++)
    {
        for(int j=0; j<param->Nz; j++)
        {
            Txx_aft[i][j] = Txxx_aft[i][j] + Txxz_aft[i][j];
        }
    }
# pragma omp parallel for num_threads(threads_count)
    for(int i=0;i<param->Nx; i++)
    {
        for(int j=0; j<param->Nz; j++)
        {
            Tzz_aft[i][j] = Tzzx_aft[i][j] + Tzzz_aft[i][j];
        }
    }
# pragma omp parallel for num_threads(threads_count)
    for(int i=0;i<param->Nx; i++)
    {
        for(int j=0; j<param->Nz; j++)
        {
            Txz_aft[i][j] = Txzx_aft[i][j] + Txzz_aft[i][j];
        }
    }

    //震源加载
    float gauss;
# pragma omp parallel for num_threads(threads_count)
    for(int i=0;i<1;i++)
    {
        for(int j=0;j<1;j++)
        {
            if(param->nx_location+i >= param->PML
            && param->nx_location+i <= param->PML+param->NX
            && param->nz_location+j >= param->PML
            && param->nz_location+j <= param->PML + param->NZ)
            {
                gauss = exp(-1.0*(i*i+j*j)*1.0/10.0);
                Txx_aft[param->nx_location+i][param->nz_location+j]
                    =gauss*wavelet+Txx_aft[param->nx_location+i][param->nz_location+j];
                Tzz_aft[param->nx_location+i][param->nz_location+j]
                    =gauss*wavelet+Tzz_aft[param->nx_location+i][param->nz_location+j];
                Txxx_aft[param->nx_location+i][param->nz_location+j]
                    = gauss*0.5*wavelet+Txxx_aft[param->nx_location+i][param->nz_location+j];
                Txxz_aft[param->nx_location+i][param->nz_location+j]
                    = gauss*0.5*wavelet+Txxz_aft[param->nx_location+i][param->nz_location+j];
                Tzzx_aft[param->nx_location+i][param->nz_location+j]
                    = gauss*0.5*wavelet+Tzzx_aft[param->nx_location+i][param->nz_location+j];
                Tzzz_aft[param->nx_location+i][param->nz_location+j]
                    = gauss*0.5*wavelet+Tzzz_aft[param->nx_location+i][param->nz_location+j];

                /*Txz_aft[param->nx_location+i][param->nz_location+j]*/
                    //= gauss*wavelet+Txz_aft[param->nx_location+i][param->nz_location+j];
                //Txzx_aft[param->nx_location+i][param->nz_location+j]
                    //= gauss*0.5*wavelet+Txzx_aft[param->nx_location+i][param->nz_location+j];
                //Txzz_aft[param->nx_location+i][param->nz_location+j]
                    //= gauss*0.5*wavelet+Txzz_aft[param->nx_location+i][param->nz_location+j];
            }
        }
    }


    //应力时间切片数值替换
# pragma omp parallel for num_threads(threads_count)
    for(int i=0; i<param->Nx; i++)
    {
        for(int j=0; j<param->Nz; j++)
        {
            Txx_now[i][j]  = Txx_aft[i][j];
            Txxx_now[i][j] = Txxx_aft[i][j];
            Txxz_now[i][j] = Txxz_aft[i][j];
            Txz_now[i][j]  = Txz_aft[i][j];
            Txzx_now[i][j] = Txzx_aft[i][j];
            Txzz_now[i][j] = Txzz_aft[i][j];
            Tzz_now[i][j]  = Tzz_aft[i][j];
            Tzzx_now[i][j] = Tzzx_aft[i][j];
            Tzzz_now[i][j] = Tzzz_aft[i][j];
        }
    }

    return 0.0;
}

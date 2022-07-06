#include "fbCommon.h"
#include "../main/cpp_code/slantstack3d.h"

// post-processing of the lslrt result.
void postProcessingTauPspectrum_diag(
        int nphr,
        int ntau,
        float **taup2d_lslrt
        )
{
    // set the cos taper along the ray-parameter axis.
    for ( int it = 0; it < ntau ; ++it )
    {
        int iphrend = (int)(nphr*(ntau-it-1)/ntau);
        for ( int iphr = iphrend ; iphr < nphr; ++iphr )
            taup2d_lslrt[iphr][it]	= 0.;

	    // apply taper function along phr axis.
        int taperLen    = 5;
        int iphrbeg   = iphrend - taperLen;
        if ( iphrbeg <= 0 ) iphrbeg   = 0;
        taperLen    = iphrend - iphrbeg + 1;
        if ( taperLen <= 2 )    continue;
        for ( int iphr = iphrbeg; iphr <= iphrend ; ++iphr )
        {
            float   x = 0.5*PI*(iphr-iphrbeg)/taperLen;
            taup2d_lslrt[iphr][it]  *= cos(x);
        }
    }
}

// post-processing of the lslrt result.
int postProcessingTauPspectrum(
        float **taup2d_frss,
        int nphr,
        int ntau,
        float **taup2d_lslrt
        )
{
    for ( int iphr = 0 ; iphr < nphr; ++iphr )
    {
        int itend   = ntau - 1;
        for ( int it = ntau-1 ; it >= 0 ; --it )
        {
            if ( 0 != taup2d_frss[iphr][it] )
            {
                itend   = it;
                goto    NextIphr;
            }
        }
NextIphr:
        ;
        if ( itend != (ntau-1) )
        {
            for ( int it = itend; it < ntau ; ++it )
                taup2d_lslrt[iphr][it]  = 0;

            // set the cos taper along the time axis.
            int taperLen    = 50;
            int itbeg   = itend - taperLen;
            if ( itbeg <= 0 )
                itbeg   = 0;
            taperLen    = itend - itbeg + 1;
            for ( int it = itbeg; it <= itend ; ++it )
            {
                float   x = 0.5*PI*(it-itbeg)/taperLen;
                taup2d_lslrt[iphr][it]  *= cos(x);
            }
        }
    }
    // set the cos taper along the ray-parameter axis.
    for ( int it = 0; it < ntau ; ++it )
    {
        int iphrend = nphr-1 ;
        for ( int iphr = nphr-1 ; iphr >= 0 ; --iphr )
        {
            if ( 0 != taup2d_lslrt[iphr][it] )
            {
                iphrend   = iphr;
                goto    NextIt;
            }
        }
NextIt:
        ;
        // further mute the LRT spectrum.
        //iphrend -= 20;
        iphrend -= 5;
        if ( iphrend <= 0 )
            iphrend   = 0;
        for ( int iphr = iphrend ; iphr < nphr; ++iphr )
            taup2d_lslrt[iphr][it]	= 0.;

	    // apply taper function along phr axis.
        int taperLen    = 10;
        int iphrbeg   = iphrend - taperLen;
        if ( iphrbeg <= 0 )
            iphrbeg   = 0;
        taperLen    = iphrend - iphrbeg + 1;
        for ( int iphr = iphrbeg; iphr <= iphrend ; ++iphr )
        {
            float   x = 0.5*PI*(iphr-iphrbeg)/taperLen;
            taup2d_lslrt[iphr][it]  *= cos(x);
        }
    }
}

// 2-D Linear Radon Transform in the time-space domain.
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
        )
{
    /*	discription.	*/
    /*	This function will do the local slant-stack over local traces.
     *
     *	float **trace;		the INPUT local seismic-gather, i.e. trace[ntrace][ns].
     *	int ntrace;		the trace-number of the INPUT gather.
     *	int ns;			the trace length.
     *	float dt; 		time-sampling interval, (unit of dt should be second),
     *
     *	float *coor;		coor[ntrace] stores the offset of each trace.
     *	float x0;		the beam-center coordinate.
     *
     *	float **tauppanel;	the OUTPUT tau-p spectrum, i.e. tauppanel[npsr][ns].
     *	int npsr;		the ray-parameter number of tau-p panel.
     *	float psrmin;		the minimum ray-parameter, (unit is s/m).
     *	float dpsr;		the interval of ray-parameter, (unit is s/m).
     *
     */
    /*
       for ( int itrace = 0 ; itrace < ntrace ; ++ itrace )
       {
       hilbert_transform_fftw( &trace[itrace][0], &trace[itrace][0], ns, 1 );
       }
       */

    memset( (void*)&tauppanel[0][0], 0, sizeof(float)*npsr*ns );

    // Local slant-stack of data from space-time window.
#pragma omp parallel for
    for ( int ipsr = 0 ; ipsr < npsr ; ++ ipsr )
    {
        float psr	= psrmin + dpsr * ipsr ;
        //printf("psr=%f\n", psr*1E6);

        for ( int it = 0 ; it < ns ; ++ it )
        {
            for ( int itrace = 0 ; itrace < ntrace ; ++ itrace )
            {
                float	xcurr	= coor[itrace];
                float	dis	=  xcurr - x0 ;

                float	tshift	= dis * psr ;
                int	itshift	= it + (int)(tshift/dt);
                if ( itshift >=0 && itshift < ns )
                    tauppanel[ipsr][it] += trace[itrace][itshift] ;
            }
        }
    }

    // Normalize.
    for ( int itrace = 0 ; itrace < npsr ; ++ itrace )
        for ( int it = 0 ; it < ns ; ++ it )
            tauppanel[itrace][it] /= ntrace;
    return 0;
}

// 2-D Linear Radon Transform in the time-space domain.
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
        )
{
    /*	discription.	*/
    /*	This function will do the local slant-stack over local traces.
     *
     *	float **trace;		the INPUT local seismic-gather, i.e. trace[ntrace][ns].
     *	int ntrace;		the trace-number of the INPUT gather.
     *	int ns;			the trace length.
     *	float dt; 		time-sampling interval, (unit of dt should be second),
     *
     *	float *coor;		coor[ntrace] stores the offset of each trace.
     *	float x0;		the beam-center coordinate.
     *
     *	float **tauppanel;	the OUTPUT tau-p spectrum, i.e. tauppanel[npsr][ns].
     *	int npsr;		the ray-parameter number of tau-p panel.
     *	float psrmin;		the minimum ray-parameter, (unit is s/m).
     *	float dpsr;		the interval of ray-parameter, (unit is s/m).
     *
     */
    /*
       for ( int itrace = 0 ; itrace < ntrace ; ++ itrace )
       {
       hilbert_transform_fftw( &trace[itrace][0], &trace[itrace][0], ns, 1 );
       }
       */

    memset( (void*)&tauppanel[0][0], 0, sizeof(float)*npsr*ns );

    // Local slant-stack of data from space-time window.
    for ( int ipsr = 0 ; ipsr < npsr ; ++ ipsr )
    {
        float psr	= psrmin + dpsr * ipsr ;
        //printf("psr=%f\n", psr*1E6);

        for ( int it = 0 ; it < ns ; ++ it )
        {
            for ( int itrace = 0 ; itrace < ntrace ; ++ itrace )
            {
                float	xcurr	= coor[itrace];
                float	dis	=  xcurr - x0 ;

                float	tshift	= dis * psr ;
                int	itshift	= it + (int)(tshift/dt);
                if ( itshift >=0 && itshift < ns )
                    tauppanel[ipsr][it] += trace[itrace][itshift] ;
            }
        }
    }

    // Normalize.
    for ( int itrace = 0 ; itrace < npsr ; ++ itrace )
        for ( int it = 0 ; it < ns ; ++ it )
            tauppanel[itrace][it] /= ntrace;
    return 0;
}


// apply the least-squares LRT to a single CMP gather.
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
        )
{
    int	 mystat	= 0;

    // set the sliding window over offset dim.
    float	offsetMin   = 0.;
    float	offsetMax   = 0.;
    for ( int itr = 0 ; itr < ntrCMP; ++itr )
    {
        if (offset[itr] > offsetMax)
            offsetMax	= offset[itr];
    }
    int Ngroup  = (int)( (offsetMax-offsetMin)/offsetWidth + 0.5);
    if ( 1 == verbose )
        fprintf(stdout, "ntrCMP=%d offsetMax=%f Nwindow=%d\n", ntrCMP, offsetMax, Ngroup);

    // apply the sliding spatial window along the offset direction.
    for ( int igroup = 0 ; igroup < Ngroup; ++igroup )
    {
        float   **gatherLoc    = alloc2float(ns, ntrCMP);
        float   *offsetLoc     = calloc(ntrCMP, FLOAT_SIZE_BYTES);
        fbsegy_std  *suhdrLoc  = calloc(ntrCMP, sizeof(fbsegy_std));
        memset( (void*)&gatherLoc[0][0], 0, FLOAT_SIZE_BYTES*1L*ns*ntrCMP);

        float	offset_beg = offsetMin + offsetWidth*igroup;
        float	offset_end = offset_beg + offsetWidth;
        if ( offset_end > offsetMax )
            offset_end = offsetMax;

        int jtr = 0;	// trace number within the spatial window: [offset_beg, offset_end]
        for ( int itr = 0 ; itr < ntrCMP; ++itr )
        {
            float   offsetNow   = offset[itr];
            if (offsetNow >= offset_beg && offsetNow <offset_end)
            {
                offsetLoc[jtr]  = offset[itr]*0.5;	// the half offset.
                suhdrLoc[jtr]   = suhdr[itr];
                suhdrLoc[jtr].offset    = (int)offsetNow;
                for ( int it = 0 ; it < ns; ++it )
                    gatherLoc[jtr][it]  = gather[itr][it];

                //fprintf(stdout, "jtr=%d offset=%f \n", jtr+1, offsetNow);
                ++jtr;
            }
        }
        int ntrLoc  = jtr;
        if ( 1 == verbose )
        fprintf(stdout, "igroup=%d ntrLoc=%d offset_beg=%f offset_end=%f\n",
                igroup+1, ntrLoc, offset_beg, offset_end);
        if ( 0 == ntrLoc )
            continue;
        //char    locaGatherFilename[FILE_NAME_MAX_LENGTH]="";
        //sprintf(locaGatherFilename,"./gatherLocOffset%f_%f.su", offset_beg, offset_end);
        //write_2d_float_wb_suhdr(gatherLoc, suhdrLoc, ntrLoc, ns, locaGatherFilename);

        float   x0  = 0;
        //fprintf(stdout, "ntr=%d x0=%f dt=%f(s)\n", ntrLoc, x0, dt_s);
        float   **taup2d_loc	= alloc2float(ns, nphr);
        /*
        // apply the classical LRT of the local traces.
        slantStack_TD_2D(gatherLoc, ntrLoc, ns, dt_s, offsetLoc, x0,
                taup2d_loc, nphr, phrmin, dphr);
        */
        // For XiangJian, add the time-domain Least-squares LRT code here.
        //void apply_GlobalLSLRT_to_CMPgather(int ntrCMP, float *offset, int ns, float dt_s,\
            float **gather, int nphr, float phrmin, float dphr, int ntau,float **taup2d_lslrt,\
            float factor_L2=1.0,float factor_L1=1.0, \
            int iterations_num=45, float residual_ratio=0.1)
        // Assume the "x0" has been removed from the coordinates "offsetLoc"!
        // A larger regularization "factor_L2" is helpful to the convergence in noisy data!
        float factor_L2=(1000.0), factor_L1=(0.0), iterations_num=(85), residual_ratio=(0.1);
        apply_GlobalLSLRT_to_CMPgather( ntrLoc, offsetLoc, ns, dt_s,\
            gatherLoc, nphr, phrmin, dphr, ntau,taup2d_loc,\
            factor_L2, factor_L1,iterations_num, residual_ratio);

        for ( int iphr = 0 ; iphr < nphr; ++iphr )
            for ( int it   = 0 ; it   < ntau; ++it   )
                taupSpectrum[iphr][it]	+= taup2d_loc[iphr][it];

        // free the allocated memory.
        free2float(taup2d_loc);
        free2float(gatherLoc);
        free(offsetLoc);
        free(suhdrLoc);
    }

    return mystat;
}

// apply the local-LRT to a single CMP gather.
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
        )
{
    int	 mystat	= 0;

    // set the sliding window over offset dim.
    float	offsetMin   = 0.;
    float	offsetMax   = 0.;
    for ( int itr = 0 ; itr < ntrCMP; ++itr )
    {
        if (offset[itr] > offsetMax)
            offsetMax	= offset[itr];
    }
    //2022.05.07, FengBo set the staring offset to avoid the residual surface waves.
    offsetMin   = 800.;
    int Ngroup  = (int)( (offsetMax-offsetMin)/offsetWidth + 0.5);
    if ( 1 == verbose )
        fprintf(stdout, "ntrCMP=%d offsetMax=%f Nwindow=%d\n", ntrCMP, offsetMax, Ngroup);

    // apply the sliding spatial window along the offset direction.
    for ( int igroup = 0 ; igroup < Ngroup; ++igroup )
    {
        float   **gatherLoc    = alloc2float(ns, ntrCMP);
        float   *offsetLoc     = calloc(ntrCMP, FLOAT_SIZE_BYTES);
        memset( (void*)&gatherLoc[0][0], 0, FLOAT_SIZE_BYTES*1L*ns*ntrCMP);

        float	offset_beg = offsetMin + offsetWidth*igroup;
        float	offset_end = offset_beg + offsetWidth;
        if ( offset_end > offsetMax )
            offset_end = offsetMax;

        int jtr = 0;	// trace number within the spatial window: [offset_beg, offset_end]
        for ( int itr = 0 ; itr < ntrCMP; ++itr )
        {
            float   offsetNow   = offset[itr];
            if (offsetNow >= offset_beg && offsetNow <offset_end)
            {
                offsetLoc[jtr]  = offset[itr]*0.5;	// the half offset.
                for ( int it = 0 ; it < ns; ++it )
                    gatherLoc[jtr][it]  = gather[itr][it];

                //fprintf(stdout, "jtr=%d offset=%f \n", jtr+1, offsetNow);
                ++jtr;
            }
        }
        int ntrLoc  = jtr;
        if ( 1 == verbose )
            fprintf(stdout, "igroup=%d ntrLoc=%d offset_beg=%f offset_end=%f\n",
                igroup+1, ntrLoc, offset_beg, offset_end);
        if ( 0 == ntrLoc )
            continue;

        char    locaGatherFilename[FILE_NAME_MAX_LENGTH]="";
        sprintf(locaGatherFilename,"./gatherLocOffset%f_%f.dat", offset_beg, offset_end);
        write_2d_float_wb(gatherLoc, ntrLoc, ns, locaGatherFilename);

        float   x0  = 0;
        //fprintf(stdout, "ntr=%d x0=%f dt=%f(s)\n", ntrLoc, x0, dt_s);
        float   **taup2d_loc	= alloc2float(ns, nphr);
        // apply the classical LRT of the local traces.
        slantStack_TD_2D(gatherLoc, ntrLoc, ns, dt_s, offsetLoc, x0,
                taup2d_loc, nphr, phrmin, dphr);

        for ( int iphr = 0 ; iphr < nphr; ++iphr )
            for ( int it   = 0 ; it   < ntau; ++it   )
                taupSpectrum[iphr][it]	+= taup2d_loc[iphr][it];

        // free the allocated memory.
        free2float(taup2d_loc);
        free2float(gatherLoc);
        free(offsetLoc);
    }

    return mystat;
}



#include "fbopw3d.h"

/*	tau-p transform functions.				*/
/*	School of Ocean and Earth Science, Tongji University.	*/
/*	Author: Feng Bo						*/
/*	Date:   2008.04						*/
/*	Update: 2008.08						*/
/*	Update: 2009.02						*/
/*	Update: 2009.06						*/
/*	Update: 2009.07						*/

/*	Function Prototypes.		*/
/*
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
*/

int	tauph_transform_3d_to_2d_robust( int ns, int ntau, int nphr, int ntr, int itaus, int itaue,
        float **tcmp, float **pdat, float *vrms, float *offset,
        float dt, float phrmin, float dphr)
    /*      Definition of some parameters in function interface.
            int	ns  ;		// Sampling number of the input CDP gather: tcmp[ntr][ns].
            int	ntau;		// Tau-domain sampling number of the output tau-ph data: pdat[nphr][ntau].
            int	nphr;		// The number of Ph-Rays.
            int	ntr ;		// Trace number within the current CDP gather.
            int	itaus ;		// The first Tau-domain sampling point number of the output tau-ph data: pdat[nphr][ntau], default is 0.
            int	itaue ;		// The last  Tau-domain sampling point number of the output tau-ph data: pdat[nphr][ntau], default is ntau.
            float	**tcmp;		// The current CDP gather is in buffer: tcmp[ntr][ns].
            float	**pdat;		// Tau-Ph transform of the current CDP gather is in buffer: pdat[nphr][ntau].
            float	*vrms ;		// The RMS velocity of the current CDP gather is in buffer: vrms[ns].
            float	*offset;	// The abs offset sqrt[(gx-sx)^2+(gy-sy)^2] of the current CDP gather is in buffer: offset[ntr].
            float	dt;		// Sampling rate of the input CDP gather. Warning: the unit of dt should be mili-second(ms).
            float	phrmin;		// The minimum value of Ph-Rays.
            float	dphr;		// The interval of Ph-Rays.
            */

{
    /*
       int	parameter_flag;
       parameter_flag = 1 ;	// Print the parameters on the screen.
       if ( parameter_flag == 1 )
       {
       printf(" ********  Parameters of 3d tau-ph data volume.  ******** \n");
       printf(" ns    = %d\n", ns   );
       printf(" ntau  = %d\n", ntau );
       printf(" nphr  = %d\n", nphr );
       printf(" ntr   = %d\n", ntr  );
       printf(" itaus = %d\n", itaus);
       printf(" itaue = %d\n", itaue);
       printf(" dt      = %f\n", dt     );
       printf(" phrmin  = %f\n", phrmin );
       printf(" dphr    = %f\n", dphr   );
       }
       */

    /*	Initializing the local variables.			*/
    dt 	*= 0.001 ;	//unit of dt is second(s).

    float	*verr=NULL ;
    verr    = (float *)calloc(ntau,FLOAT_SIZE_BYTES) ;

    float	velerrmin = 100.0;	/* The stack-velocity error range: minimum value. */
    float	velerrmax = 600.0 ;	/* The stack-velocity error range: maximum value. */

    /*	Initializing the output buffer: pdat[nphr][ntau].	*/
    zero2float(pdat,ntau,nphr);

    int	**vt0_taup=NULL;
    vt0_taup = alloc2int(ntau,nphr);
    zero2int(vt0_taup,ntau,nphr);

    /*	Search Mute Line.					*/
    int	*itmuteh=NULL;
    itmuteh = alloc1int(ntr);

    float	epsmin = 1.0E-30 ;
    int	it, itr, itmute ;
    for( itr = 0 ; itr < ntr ; ++ itr )
    {
        for( it = 0 ; it < ns ; ++ it )
        {
            if( fabs(tcmp[itr][it]) > epsmin )
            {
                itmuteh[itr] = it ;
                break ;
            }
        }
    }
    itmute  = Min_Int_1d(itmuteh,ntr) ;
    //write_1d_int_wb(itmuteh,ntr,"itmuteh.dat");

    float	offmax, offmin, velmax, velmin ;
    float	offmaxtmp, offmin2, offmax2, vvmax2, vvmin2 ;
    offmax  = Max_Float_1d(offset, ntr) ;
    offmin  = Min_Float_1d(offset, ntr) ;
    velmax  = Max_Float_1d(vrms,   ns) ;
    velmin  = Min_Float_1d(vrms,   ns) ;


    /*	Set the tau-axis range. tau=[t0min, t0max].		*/
    int	t0min, t0max ;
    t0min = MAX ( itmute, itaus ) ;
    t0max = MIN ( ns,     itaue ) ;

    int     it0, itau, iphr, iphrmin, iphrmax, itime ;
    float	phmax, phmin, phr, phr2, phrr ;
    float   t0, tau, time, vt0, vv2, hr ;
    float	rrmax, rrmin ;

    /*
       printf(" offmax = %f\n", offmax);
       printf(" offmin = %f\n", offmin);
       printf(" velmax = %f\n", velmax);
       printf(" velmin = %f\n", velmin);
       printf(" t0min  = %d\n", t0min);
       printf(" t0max  = %d\n", t0max);
       for( it0  = 0 ; it0 < ntau ; ++it0 )
       for( iphr = 0 ; iphr < nphr ; ++iphr )
       vt0_taup[iphr][it0] = 0 ;
       */

    /*	Scan the tau axis from t0min to t0max.			*/
    for( it0 = t0min ; it0 < t0max ; it0 += 1 )
    {
        t0	 = dt*it0 ;
        vt0	 = vrms[it0]*t0 ;
        verr[it0]= (velerrmax-velerrmin)*(it0-t0min)/(t0max-t0min) + velerrmin ;

        /*		Compute the range of Ph-value by Offset[min,max], t0, vrms(t0)&verror.	*/
        vv2	= vrms[it0]*vrms[it0];
        //vv2	= (vrms[it0]+verr[it0])*(vrms[it0]+verr[it0]) ;
        offmaxtmp = offmax;
        for( itr = 0 ; itr < ntr ; ++ itr )
        {
            time = sqrt(t0*t0 + (offset[itr]*offset[itr])/vv2);
            if( itmuteh[itr] > time/dt )
            {
                offmaxtmp = offset[itr] ;
                break ;
            }
        }

        /*		Ph_min=f(vt0max,hmin), Ph_max=f(vt0min,hmax).	*/
        offmax2 = offmaxtmp*offmaxtmp ;
        offmin2 = offmin*offmin ;
        velmax	= vrms[it0] + verr[it0] ;
        velmin	= vrms[it0] - verr[it0] ;
        vvmax2	= velmax*velmax ;
        vvmin2	= velmin*velmin ;

        phmax   = 2.0*offmaxtmp/(velmin*sqrt(vvmin2*t0*t0 + offmax2)) ;
        iphrmax = MIN( (int)(phmax/dphr), nphr ) ;
        phmin   = 2.0*offmin/(velmax*sqrt(vvmax2*t0*t0 + offmin2)) ;
        iphrmin = MIN( (int)(phmin/dphr), iphrmax ) ;

        //printf(" it=%4d iphrmax=%4d iphrmin=%4d\n", it0, iphrmax, iphrmin);
        /*		Scan the p axis from phmin to phmax.			*/
        vv2	= vrms[it0]*vrms[it0];
        for( iphr = iphrmin ; iphr < iphrmax ; ++ iphr )
        {
            phr  = phrmin + dphr*iphr ;
            phr2 = phr*phr ;

            /*			Compute the intercept time tau by t0, p, v(t0).		*/
            float s=1.0-phr2*vv2/4.0;
            if(s<0)	continue;
            //tau  = t0*sqrt(1.0-phr2*vv2/4.0) ;
            tau  = t0*sqrt(s);
            itau = (int)(tau/dt) + 1 ;
            vt0_taup[iphr][itau] = it0 ;

        }
    }
    //write_2d_int_wb(vt0_taup, nphr, ntau, "vt0.dat");

    /*	Finite-Range slant-stack begin.		*/
    for( itau = t0min ; itau < t0max ; ++ itau )
    {
        //printf(" itau      = %d\n", itau     );
        tau	= dt*itau ;

        for( iphr = 0 ; iphr < nphr ; ++ iphr )
        {
            if( 0 == vt0_taup[iphr][itau] )	continue ;

            phr	= phrmin + dphr*iphr ;
            phr2	= phr*phr ;

            it0	= vt0_taup[iphr][itau] ;
            velmax	= vrms[it0] + verr[it0] ;
            velmin	= vrms[it0] - verr[it0] ;

            vvmax2	= velmax*velmax ;
            vvmin2	= velmin*velmin ;
            rrmax	= phr*vvmax2*tau/(4.0-phr2*vvmax2);
            rrmin	= phr*vvmin2*tau/(4.0-phr2*vvmin2);

            for( itr = 0 ; itr < ntr ; ++ itr )
            {
                hr	= offset[itr]*0.5 ;
                if((hr>rrmin) && (hr<rrmax))
                {
                    time  = tau + phr*hr ;
                    itime = (int)(time/dt) ;
                    if( itime > 0 && itime < ns )
                        pdat[iphr][itau] += sinc_interpolation(tcmp[itr],ns,dt,time);
                }
            }
        }
    }
    //write_2d_int_wb(vt0_taup, nphr, ntau, "vt0_taup.dat") ;

    free(verr);
    free1int(itmuteh) ;
    free2int(vt0_taup) ;

    return 0 ;
}

int	tauph_transform_3d_to_2d_tmp( int ns, int ntau, int nphr, int ntr, int itaus, int itaue,
        float **tcmp, float **pdat, float *vrms, float *offset,
        float dt, float phrmin, float dphr)
    /*      Definition of some parameters in function interface.
            int	ns  ;		// Sampling number of the input CDP gather: tcmp[ntr][ns].
            int	ntau;		// Tau-domain sampling number of the output tau-ph data: pdat[nphr][ntau].
            int	nphr;		// The number of Ph-Rays.
            int	ntr ;		// Trace number within the current CDP gather.
            int	itaus ;		// The first Tau-domain sampling point number of the output tau-ph data: pdat[nphr][ntau], default is 0.
            int	itaue ;		// The last  Tau-domain sampling point number of the output tau-ph data: pdat[nphr][ntau], default is ntau.
            float	**tcmp;		// The current CDP gather is in buffer: tcmp[ntr][ns].
            float	**pdat;		// Tau-Ph transform of the current CDP gather is in buffer: pdat[nphr][ntau].
            float	*vrms ;		// The RMS velocity of the current CDP gather is in buffer: vrms[ns].
            float	*offset;	// The abs offset sqrt[(gx-sx)^2+(gy-sy)^2] of the current CDP gather is in buffer: offset[ntr].
            float	dt;		// Sampling rate of the input CDP gather. Warning: the unit of dt should be mili-second(ms).
            float	phrmin;		// The minimum value of Ph-Rays.
            float	dphr;		// The interval of Ph-Rays.
            */

{
    /*
       int	parameter_flag;
       parameter_flag = 1 ;	// Print the parameters on the screen.
       if ( parameter_flag == 1 )
       {
       printf(" ********  Parameters of 3d tau-ph data volume.  ******** \n");
       printf(" ns    = %d\n", ns   );
       printf(" ntau  = %d\n", ntau );
       printf(" nphr  = %d\n", nphr );
       printf(" ntr   = %d\n", ntr  );
       printf(" itaus = %d\n", itaus);
       printf(" itaue = %d\n", itaue);
       printf(" dt      = %f\n", dt     );
       printf(" phrmin  = %f\n", phrmin );
       printf(" dphr    = %f\n", dphr   );
       }
       */

    /*	Initializing the local variables.			*/
    float	epsmin = 1.0E-30 ;
    float	*velerr=NULL, velerrmin, velerrmax ;

    dt 	*= 0.001 ;	//unit of dt is second(s).

    velerrmin = 50 ;
    velerrmax = 300 ;
    velerr    = (float *)calloc(ntau,FLOAT_SIZE_BYTES) ;

    /*	Initializing the output buffer: pdat[nphr][ntau].	*/
    zero2float(pdat,ntau,nphr);
    int	*itmuteh=NULL;
    itmuteh = alloc1int(ntr);

    int	it, itr, itmute, it0, t0min, t0max ;
    float	t0, offmax, offmin, velmax, velmin, phmax, phmin, offmaxtmp ;

    int	**vt0_taup=NULL;
    vt0_taup = alloc2int(ntau,nphr);
    zero2int(vt0_taup,ntau,nphr);

    /*	Search Mute Line.					*/
    for( itr = 0 ; itr < ntr ; ++ itr )
    {
        for( it = 0 ; it < ns ; ++ it )
        {
            if( fabs(tcmp[itr][it]) > epsmin )
            {
                itmuteh[itr] = it ;
                break ;
            }
        }
    }
    itmute  = Min_Int_1d(itmuteh,ntr) ;

    offmax  = Max_Float_1d(offset, ntr) ;
    offmin  = Min_Float_1d(offset, ntr) ;
    velmax  = Max_Float_1d(vrms,   ns) ;
    velmin  = Min_Float_1d(vrms,   ns) ;

    /*	Set the tau-axis range. tau=[t0min, t0max].		*/
    t0min = MAX ( itmute, itaus ) ;
    t0max = MIN ( ns,     itaue ) ;

    int     itau, iphr, iphrmin, iphrmax, itime ;
    float   tau, time, vv2, hr, phr, phr2, phrr, vt0, vvtmp2 ;
    float   hhtmp, hhmax, hhmin, coef, offmin2, offmax2 ;
    float	rrmax, rrmin, vvmax2, vvmin2 ;
    float	velperc=0.10;

    printf(" offmax = %f\n", offmax);
    printf(" offmin = %f\n", offmin);
    printf(" velmax = %f\n", velmax);
    printf(" velmin = %f\n", velmin);
    printf(" t0min  = %d\n", t0min);
    printf(" t0max  = %d\n", t0max);
    /*
       for( it0  = 0 ; it0 < ntau ; ++it0 )
       for( iphr = 0 ; iphr < nphr ; ++iphr )
       vt0_taup[iphr][it0] = velmin ;
       */

    /*	Scan the tau axis from t0min to t0max.			*/
    for( it0 = t0min ; it0 < t0max ; it0 += 1 )
    {
        t0	= dt*it0 ;
        vt0	= vrms[it0]*t0 ;
        vv2	= vrms[it0]*vrms[it0];

        velerr[it0] = (velerrmax-velerrmin)*(it0-t0min)/(t0max-t0min) + velerrmin ;

        /*		Compute the range of Ph-value by Offset[min,max], t0, vrms(t0).	*/
        offmaxtmp = offmax;
        for( itr = 0 ; itr < ntr ; ++ itr )
        {
            time = sqrt(t0*t0 + (offset[itr]*offset[itr])/vv2);
            if( itmuteh[itr] > time )
            {
                offmaxtmp = offset[itr] ;
                break ;
            }
        }
        offmax2 = offmaxtmp*offmaxtmp ;
        offmin2 = offmin*offmin ;
        phmax   = 2.0*offmaxtmp/(vrms[it0]*sqrt(vt0*vt0 + offmax2)) ;
        iphrmax = MIN( (int)(phmax/dphr), nphr ) ;
        phmin   = 2.0*offmin/(vrms[it0]*sqrt(vt0*vt0 + offmin2)) ;
        iphrmin = MIN( (int)(phmin/dphr), iphrmax ) ;

        /*		Scan the p axis from phmin to phmax.			*/
        for( iphr = iphrmin ; iphr < iphrmax ; ++ iphr )
        {
            phr  = phrmin + dphr*iphr ;
            phr2 = phr*phr ;

            /*			Compute the intercept time tau by t0, p, v(t0).		*/
            tau  = t0*sqrt(1.0-phr2*vv2/4.0) ;
            itau = (int)(tau/dt) + 1 ;
            vt0_taup[iphr][itau] = it0 ;
            //vt0_taup[iphr][itau] = (int)vrms[it0] ;

        }
    }

    float	drr, drrs, drre, rrtmp ;
    drrs	= 50.0 ;
    drre	= 400.0 ;

    /*	Finite-Range slant-stack begin.		*/
    for( itau = t0min ; itau < t0max ; ++ itau )
    {
        //printf(" itau      = %d\n", itau     );
        tau	= dt*itau ;
        //velerr	= (velerrmax-velerrmin)*(it0-t0min)/(t0max-t0min) + velerrmin ;
        //drr	= (drre-drrs)*(itau-itaus)/(itaue-itaus) + drrs ;        //straight line relation

        for( iphr = 0 ; iphr < nphr ; ++ iphr )
        {
            if( 0 == vt0_taup[iphr][itau] )	continue ;

            phr	= phrmin + dphr*iphr ;
            phr2	= phr*phr ;

            it0	= vt0_taup[iphr][itau] ;
            velmax	= vrms[it0] + velerr[it0] ;
            velmin	= vrms[it0] - velerr[it0] ;
            //velmin	= vt0_taup[iphr][itau] - velerr ;
            //velmin	= vt0_taup[iphr][itau]*(1.0-velperc);

            vvmax2	= velmax*velmax ;
            vvmin2	= velmin*velmin ;
            rrmax	= phr*vvmax2*tau/(4.0-phr2*vvmax2);
            rrmin	= phr*vvmin2*tau/(4.0-phr2*vvmin2);

            /*
               vv2	= vt0_taup[iphr][itau]*vt0_taup[iphr][itau] ;
               rrtmp   = phr*vv2*tau/(4.0-phr2*vv2);
               rrtmp   = fabs(rrtmp);
               rrmin   = MAX( (rrtmp-drr), 0.0 );
               rrmax   = rrtmp + drr ;
               */

            for( itr = 0 ; itr < ntr ; ++ itr )
            {
                hr	= offset[itr]*0.5 ;
                if((hr>rrmin) && (hr<rrmax))
                {
                    time  = tau + phr*hr ;
                    itime = (int)(time/dt) ;
                    if( itime > 0 && itime < ns )
                        pdat[iphr][itau] += sinc_interpolation(tcmp[itr],ns,dt,time);
                    /*
                       if( itime > 0 && itime < ns )
                       pdat[iphr][itau] += tcmp[itr][itime] ;
                       */
                }
            }
        }
    }
    //write_2d_int_wb(vt0_taup, nphr, ntau, "vt0_taup.dat") ;

    free(velerr);
    free1int(itmuteh) ;
    free2int(vt0_taup) ;

    return 0 ;
}

int	tauph_transform_3d_to_2d_sa( int ns, int ntau, int nphr, int ntr, int itaus, int itaue,
        float **tcmp, float **pdat, float *vrms, float *offset,
        float dt, float phrmin, float dphr)
    /*      Definition of some parameters in function interface.
            int	ns  ;		// Sampling number of the input CDP gather: tcmp[ntr][ns].
            int	ntau;		// Tau-domain sampling number of the output tau-ph data: pdat[nphr][ntau].
            int	nphr;		// The number of Ph-Rays.
            int	ntr ;		// Trace number within the current CDP gather.
            int	itaus ;		// The first Tau-domain sampling point number of the output tau-ph data: pdat[nphr][ntau], default is 0.
            int	itaue ;		// The last  Tau-domain sampling point number of the output tau-ph data: pdat[nphr][ntau], default is ntau.
            float	**tcmp;		// The current CDP gather is in buffer: tcmp[ntr][ns].
            float	**pdat;		// Tau-Ph transform of the current CDP gather is in buffer: pdat[nphr][ntau].
            float	*vrms ;		// The RMS velocity of the current CDP gather is in buffer: vrms[ns].
            float	*offset;	// The abs offset sqrt[(gx-sx)^2+(gy-sy)^2] of the current CDP gather is in buffer: offset[ntr].
            float	dt;		// Sampling rate of the input CDP gather. Warning: the unit of dt should be mili-second(ms).
            float	phrmin;		// The minimum value of Ph-Rays.
            float	dphr;		// The interval of Ph-Rays.
            */

{
    /*
       int	parameter_flag;
       parameter_flag = 1 ;	// Print the parameters on the screen.
       if ( parameter_flag == 1 )
       {
       printf(" ********  Parameters of 3d tau-ph data volume.  ******** \n");
       printf(" ns    = %d\n", ns   );
       printf(" ntau  = %d\n", ntau );
       printf(" nphr  = %d\n", nphr );
       printf(" ntr   = %d\n", ntr  );
       printf(" itaus = %d\n", itaus);
       printf(" itaue = %d\n", itaue);
       printf(" dt      = %f\n", dt     );
       printf(" phrmin  = %f\n", phrmin );
       printf(" dphr    = %f\n", dphr   );
       }
       */

    /*	Initializing the local variables.			*/
    float	epsmin = 1.0E-30 ;
    float	velerr, velerrmin, velerrmax, dvelerr ;
    float	gvelmin, gvelmax ;	/*  Global velocity range in a CMP gather. */
    int	givel1, giveln ;

    dt 	*= 0.001  ;	//unit of dt is second(s).
    velerrmin = 100.0 ;
    velerrmax = 500.0 ;
    dvelerr   = 50.0  ;

    /*	Initializing the output buffer: pdat[nphr][ntau].	*/
    zero2float(pdat,ntau,nphr);
    int	*itmuteh=NULL;
    itmuteh = alloc1int(ntr);

    int	it, itr, itmute, it0, t0min, t0max ;
    float	t0, offmax, offmin, velmax, velmin, phmax, phmin, offmaxtmp ;

    /*	Search Mute Line.					*/
    for( itr = 0 ; itr < ntr ; ++ itr )
    {
        for( it = 0 ; it < ns ; ++ it )
        {
            if( fabs(tcmp[itr][it]) > epsmin )
            {
                itmuteh[itr] = it ;
                break ;
            }
        }
    }
    itmute  = Min_Int_1d(itmuteh,ntr) ;

    offmax  = Max_Float_1d(offset, ntr) ;
    offmin  = Min_Float_1d(offset, ntr) ;
    velmax  = Max_Float_1d(vrms,   ns) ;
    velmin  = Min_Float_1d(vrms,   ns) ;

    gvelmin = velmin - velerrmin ;
    gvelmax = velmax + velerrmax ;
    gvelmin = ((int)(gvelmin/dvelerr))*dvelerr ;
    gvelmax = ((int)(gvelmax/dvelerr))*dvelerr ;

    givel1    = (int)(gvelmin/dvelerr) ;
    giveln    = (int)(gvelmax/dvelerr) ;

    /*	Set the velocity scan number.	*/
    int	ivel, ivelid, nvel, ivel1, iveln ;
    float	vel_lmin, vel_lmax;	/* local velocity range at time t0.	*/

    nvel     = giveln - givel1 + 1 ;

    float	**spectrum=NULL ;
    spectrum = alloc2float(ns, nvel);
    zero2float(spectrum, ns, nvel);

    /*	Set the tau-axis range. tau=[t0min, t0max].		*/
    t0min = MAX ( itmute, itaus ) ;
    t0max = MIN ( ns,     itaue ) ;

    int     itau, iphr, iphrmin, iphrmax, itime ;
    float   tau, time, v0, vv2, hr, phr, phr2, phrr, vt0, vvtmp, vvtmp2 ;
    float   pdat_tmp, hhtmp, hhmax, hhmin, coef, offmin2, offmax2, phrmax ;

    phrmax = phrmin + dphr*(nphr-1) ;

    printf(" phrmax = %f\n", phrmax);
    printf(" phrmin = %f\n", phrmin);
    printf(" offmax = %f\n", offmax);
    printf(" offmin = %f\n", offmin);
    printf(" velmax = %f\n", velmax);
    printf(" velmin = %f\n", velmin);
    printf(" t0min  = %d\n", t0min);
    printf(" t0max  = %d\n", t0max);
    printf(" nvel   = %d\n", nvel );
    printf(" gvelmin=%f\n",  gvelmin);
    printf(" gvelmax=%f\n",  gvelmax);
    /*
    */

    float	taumin, taumax, tautmp ;
    float	sum1, sum2 ;
    int 	itaumin, itaumax ;

    /*	Scan the tau axis from t0min to t0max.			*/
    for( it0 = t0min ; it0 < t0max ; ++ it0 )
    {
        t0     = dt*it0 ;
        //vt0    = vrms[it0]*t0 ;
        //vv2    = vrms[it0]*vrms[it0];
        //printf(" it0 = %d \n", it0 );

        velerr   = (velerrmax-velerrmin)*(it0-t0min)/(t0max-t0min) + velerrmin ;
        vel_lmin = MAX( (vrms[it0]-velerr), gvelmin ) ;
        vel_lmax = MIN( (vrms[it0]+velerr), gvelmax ) ;
        ivel1    = MAX( (int)(vel_lmin/dvelerr), givel1 ) ;
        iveln    = MIN( (int)(vel_lmax/dvelerr), giveln ) ;

        //for( ivel = 0 ; ivel < nvel ; ++ ivel )
        for( ivel = ivel1 ; ivel <= iveln ; ++ ivel )
        {

            ivelid = ivel - givel1 ;
            v0     = dvelerr*ivel ;
            vt0    = v0*t0 ;
            vv2    = v0*v0 ;

            /*			Compute velocity spectrum of CMP gather.	*/
            /*
               v0     = 1500+dvelerr*ivel ;
               vv2    = v0*v0 ;
               spectrum[ivel][it0] = 0 ;
               vel_lmin = 0 ;
               vel_lmax = 0 ;
               for( itr = 0 ; itr < ntr ; ++ itr )
               {
               time  = sqrt(t0*t0+offset[itr]*offset[itr]/vv2) ;
               itime = (int)(time/dt) ;
               if( itime > 0 && itime < ns )
               {
            //spectrum[ivel][it0] += pow(tcmp[itr][itime],2) ;
            vel_lmin += tcmp[itr][itime] ;
            vel_lmax += pow(tcmp[itr][itime],2) ;
            }
            }
            spectrum[ivel][it0] = vel_lmin*vel_lmin/vel_lmax ;
            continue ;
            */

            /*			Compute the range of Ph-value by Offset[min,max], t0, vrms(t0).	*/
            offmaxtmp = offmax;
            for( itr = 0 ; itr < ntr ; ++ itr )
            {
                time = sqrt(t0*t0 + (offset[itr]*offset[itr])/vv2);
                if( itmuteh[itr] > time )
                {
                    offmaxtmp = offset[itr] ;
                    break ;
                }
            }
            offmax2 = offmaxtmp*offmaxtmp ;
            offmin2 = offmin*offmin ;

            //taumin  = t0*sqrt(1.0-offmax2/(t0*t0*vel_lmin*vel_lmin+offmax2)) ;
            //taumax  = t0*sqrt(1.0-offmin2/(t0*t0*vel_lmax*vel_lmax+offmin2)) ;
            taumin  = t0*sqrt(1.0-offmax2/(vt0*vt0+offmax2)) ;
            taumax  = t0*sqrt(1.0-offmin2/(vt0*vt0+offmin2)) ;
            //itaumin = MIN( (int)(taumin/dt), t0min );
            //itaumax = MIN( (int)(taumax/dt), (ns-1) );

            /*	Scan the tau axis from taumin to taumax.			*/
            spectrum[ivelid][it0] = 0 ;
            //for( itau = itaumin ; itau <= itaumax ; ++ itau )
            tau = taumax ;
            sum1 = 0 ;
            sum2 = 0 ;
            while( tau > taumin )
            {
                phr  = 2.0*sqrt(1.0-tau*tau/(t0*t0))/v0 ;
                phr2 = phr*phr ;

                if( phr > phrmax )	break ;
                /*	Compute the half-offset by tau, p, v(t0).		*/
                //vvtmp  = v0 - dvelerr;
                vvtmp  = v0 ;
                vvtmp2 = vvtmp*vvtmp ;
                hhmin  = phr*vvtmp2*tau/(4.0-phr2*vvtmp2);

                vvtmp  = v0 + dvelerr;
                vvtmp2 = vvtmp*vvtmp ;
                hhmax  = phr*vvtmp2*tau/(4.0-phr2*vvtmp2);

                pdat_tmp = 0 ;
                for( itr = 0 ; itr < ntr ; ++ itr )
                {
                    hr = offset[itr]*0.5 ;
                    if((hr>hhmin) && (hr<hhmax))
                    {
                        time  = tau + hr*phr ;
                        itime = (int)(time/dt) ;
                        if( itime > 0 && itime < ns )
                            pdat_tmp += sinc_interpolation(tcmp[itr],ns,dt,time);
                        //pdat_tmp += tcmp[itr][itime] ;
                    }
                }

                sum1 += pdat_tmp ;
                sum2 += pdat_tmp*pdat_tmp ;
                //spectrum[ivelid][it0] += (pdat_tmp*pdat_tmp) ;
                tau -= dt ;
            }
            if( sum2 < 1.0E-5 )
                spectrum[ivelid][it0] = 0 ;
            else
                spectrum[ivelid][it0] = sum1*sum1/sum2 ;
        }
    }

    write_2d_float_wb( spectrum, nvel, ns, "spectrum_taup.dat") ;

    free1int(itmuteh) ;
    free2float(spectrum);

    return 0 ;
}

int	tauph_transform_3d_to_2d_mute_modf( int ns, int ntau, int nphr, int ntr, int itaus, int itaue,
        float **tcmp, float **pdat, float *vrms, float *offset,
        float dt, float phrmin, float dphr)
    /*      Definition of some parameters in function interface.
            int	ns  ;		// Sampling number of the input CDP gather: tcmp[ntr][ns].
            int	ntau;		// Tau-domain sampling number of the output tau-ph data: pdat[nphr][ntau].
            int	nphr;		// The number of Ph-Rays.
            int	ntr ;		// Trace number within the current CDP gather.
            int	itaus ;		// The first Tau-domain sampling point number of the output tau-ph data: pdat[nphr][ntau], default is 0.
            int	itaue ;		// The last  Tau-domain sampling point number of the output tau-ph data: pdat[nphr][ntau], default is ntau.
            float	**tcmp;		// The current CDP gather is in buffer: tcmp[ntr][ns].
            float	**pdat;		// Tau-Ph transform of the current CDP gather is in buffer: pdat[nphr][ntau].
            float	*vrms ;		// The RMS velocity of the current CDP gather is in buffer: vrms[ns].
            float	*offset;	// The abs offset sqrt[(gx-sx)^2+(gy-sy)^2] of the current CDP gather is in buffer: offset[ntr].
            float	dt;		// Sampling rate of the input CDP gather. Warning: the unit of dt should be mili-second(ms).
            float	phrmin;		// The minimum value of Ph-Rays.
            float	dphr;		// The interval of Ph-Rays.
            */

{
    /*
       int	parameter_flag;
       parameter_flag = 1 ;	// Print the parameters on the screen.
       if ( parameter_flag == 1 )
       {
       printf(" ********  Parameters of 3d tau-ph data volume.  ******** \n");
       printf(" ns    = %d\n", ns   );
       printf(" ntau  = %d\n", ntau );
       printf(" nphr  = %d\n", nphr );
       printf(" ntr   = %d\n", ntr  );
       printf(" itaus = %d\n", itaus);
       printf(" itaue = %d\n", itaue);
       printf(" dt      = %f\n", dt     );
       printf(" phrmin  = %f\n", phrmin );
       printf(" dphr    = %f\n", dphr   );
       }
       */

    /*	Initializing the local variables.			*/
    float	epsmin = 1.0E-30 ;
    float	velerr, velerrmin, velerrmax ;
    dt 	*= 0.001 ;	//unit of dt is second(s).
    velerrmin = 50 ;
    velerrmax = 300 ;

    /*	Initializing the output buffer: pdat[nphr][ntau].	*/
    zero2float(pdat,ntau,nphr);
    int	*itmuteh=NULL;
    itmuteh = alloc1int(ntr);

    int	it, itr, itmute, it0, t0min, t0max ;
    float	t0, offmax, offmin, velmax, velmin, phmax, phmin, offmaxtmp ;

    /*	Search Mute Line.					*/
    for( itr = 0 ; itr < ntr ; ++ itr )
    {
        for( it = 0 ; it < ns ; ++ it )
        {
            if( fabs(tcmp[itr][it]) > epsmin )
            {
                itmuteh[itr] = it ;
                break ;
            }
        }
    }
    itmute  = Min_Int_1d(itmuteh,ntr) ;

    offmax  = Max_Float_1d(offset, ntr) ;
    offmin  = Min_Float_1d(offset, ntr) ;
    velmax  = Max_Float_1d(vrms,   ns) ;
    velmin  = Min_Float_1d(vrms,   ns) ;

    /*	Set the tau-axis range. tau=[t0min, t0max].		*/
    t0min = MAX ( itmute, itaus ) ;
    t0max = MIN ( ns,     itaue ) ;

    int     itau, iphr, iphrmin, iphrmax, itime ;
    float   tau, time, vv2, hr, phr, phr2, phrr, vt0, vvtmp2 ;
    float   hhtmp, hhmax, hhmin, coef, offmin2, offmax2 ;

    printf(" offmax = %f\n", offmax);
    printf(" offmin = %f\n", offmin);
    printf(" velmax = %f\n", velmax);
    printf(" velmin = %f\n", velmin);
    printf(" t0min  = %d\n", t0min);
    printf(" t0max  = %d\n", t0max);
    /*
    */

    /*	Scan the tau axis from t0min to t0max.			*/
    for( it0 = t0min ; it0 < t0max ; it0 += 1 )
    {
        t0     = dt*it0 ;
        vt0    = vrms[it0]*t0 ;
        vv2    = vrms[it0]*vrms[it0];
        velerr = (velerrmax-velerrmin)*(it0-t0min)/(t0max-t0min) + velerrmin ;
        //printf(" it0 = %d velerr = %f \n", it0, velerr);

        /*		Compute the range of Ph-value by Offset[min,max], t0, vrms(t0).	*/
        offmaxtmp = offmax;
        for( itr = 0 ; itr < ntr ; ++ itr )
        {
            time = sqrt(t0*t0 + (offset[itr]*offset[itr])/vv2);
            if( itmuteh[itr] > time )
            {
                offmaxtmp = offset[itr] ;
                break ;
            }
        }
        offmax2 = offmaxtmp*offmaxtmp ;
        offmin2 = offmin*offmin ;
        phmax   = 2.0*offmaxtmp/(vrms[it0]*sqrt(vt0*vt0 + offmax2)) ;
        iphrmax = MIN( (int)(phmax/dphr), nphr ) ;
        phmin   = 2.0*offmin/(vrms[it0]*sqrt(vt0*vt0 + offmin2)) ;
        iphrmin = MIN( (int)(phmin/dphr), iphrmax ) ;

        /*		Scan the p axis from phmin to phmax.			*/
        for( iphr = iphrmin ; iphr < iphrmax ; ++ iphr )
        {
            phr  = phrmin + dphr*iphr ;
            phr2 = phr*phr ;

            /*			Compute the intercept time tau by t0, p, v(t0).		*/
            //if( phr2 > 4.0/vv2 )	continue ;
            tau  = t0*sqrt(1.0-phr2*vv2/4.0) ;
            itau = (int)(tau/dt) + 1 ;
            tau  = dt*itau ;

            /*			Compute the half-offset by tau, p, v(t0).		*/
            vvtmp2  = (vrms[it0]-velerr)*(vrms[it0]-velerr) ;
            hhmin   = phr*vvtmp2*tau/(4.0-phr2*vvtmp2);

            vvtmp2  = (vrms[it0]+velerr)*(vrms[it0]+velerr) ;
            hhmax   = phr*vvtmp2*tau/(4.0-phr2*vvtmp2);

            //if( hhmin > hhmax ) printf("Error! hhmin > hhmax \n");
            for( itr = 0 ; itr < ntr ; ++ itr )
            {
                hr = offset[itr]*0.5 ;
                if((hr>hhmin) && (hr<hhmax))
                {
                    time  = tau + hr*phr ;
                    itime = (int)(time/dt) ;
                    if( itime > 0 && itime < ns )
                        pdat[iphr][itau] += sinc_interpolation(tcmp[itr],ns,dt,time);
                    /*
                       if( hr > hhtmp )
                       coef = cos(0.5*3.14*(hr-hhtmp)/(hhmax-hhtmp)) ;
                       if( hr < hhtmp )
                       coef = cos(0.5*3.14*(hr-hhtmp)/(hhmin-hhtmp)) ;
                       if( itime > 0 && itime < ns )
                       pdat[iphr][itau] += tcmp[itr][itime] ;
                       */
                }
            }
        }
    }

    /*	tau-p stack to P=0.	*/
    int	interp_flag = 1 ;
    int	itaup, iph, iph_tau ;
    float	taup, ph, tau0, tau01 ;
    float	ph_tau, intp_coe, iphmin ;

    iphmin = 0 ;
    if( interp_flag == 1 )
    {
        for( itau = t0min ; itau < t0max ; itau += 1 )
        {
            tau0    = dt*itau ;
            ph_tau  = 2.0/vrms[itau] ;
            iph_tau = ph_tau/dphr ;
            if( iph_tau > nphr )	iph_tau = nphr ;

            for( iph = 1 ; iph < iph_tau ; ++ iph )
            {
                ph   = dphr*iph ;
                taup = tau0*sqrt(1.0-ph*ph/(ph_tau*ph_tau));
                pdat[0][itau] += sinc_interpolation(pdat[iph], ns, dt, taup);
            }

            /*
               pdat[0][itau] /= (iph_tau-1.0);
               for( iph = 1 ; iph < 40 ; ++ iph )
               {
               ph    = dphr*iph ;
               taup  = tau0*sqrt(1.0-ph*ph/(ph_tau*ph_tau));
               itaup = (int)(taup/dt);
               taup  = dt*itaup;
               tau01 = taup/sqrt(1.0-ph*ph/(ph_tau*ph_tau));
               pdat[iph][itaup] = sinc_interpolation(pdat[0], ns, dt, tau01);
               }
               */
        }
    }

    free1int(itmuteh) ;

    return 0 ;
}

int	tauph_transform_3d_to_2d_mute( int ns, int ntau, int nphr, int ntr, int itaus, int itaue,
        float **tcmp, float **pdat, float *vrms, float *offset,
        float dt, float phrmin, float dphr)
    /*      Definition of some parameters in function interface.
            int	ns  ;		// Sampling number of the input CDP gather: tcmp[ntr][ns].
            int	ntau;		// Tau-domain sampling number of the output tau-ph data: pdat[nphr][ntau].
            int	nphr;		// The number of Ph-Rays.
            int	ntr ;		// Trace number within the current CDP gather.
            int	itaus ;		// The first Tau-domain sampling point number of the output tau-ph data: pdat[nphr][ntau], default is 0.
            int	itaue ;		// The last  Tau-domain sampling point number of the output tau-ph data: pdat[nphr][ntau], default is ntau.
            float	**tcmp;		// The current CDP gather is in buffer: tcmp[ntr][ns].
            float	**pdat;		// Tau-Ph transform of the current CDP gather is in buffer: pdat[nphr][ntau].
            float	*vrms ;		// The RMS velocity of the current CDP gather is in buffer: vrms[ns].
            float	*offset;	// The abs offset sqrt[(gx-sx)^2+(gy-sy)^2] of the current CDP gather is in buffer: offset[ntr].
            float	dt;		// Sampling rate of the input CDP gather. Warning: the unit of dt should be mili-second(ms).
            float	phrmin;		// The minimum value of Ph-Rays.
            float	dphr;		// The interval of Ph-Rays.
            */

{
    /*
       int	parameter_flag;
       parameter_flag = 1 ;	// Print the parameters on the screen.
       if ( parameter_flag == 1 )
       {
       printf(" ********  Parameters of 3d tau-ph data volume.  ******** \n");
       printf(" ns    = %d\n", ns   );
       printf(" ntau  = %d\n", ntau );
       printf(" nphr  = %d\n", nphr );
       printf(" ntr   = %d\n", ntr  );
       printf(" itaus = %d\n", itaus);
       printf(" itaue = %d\n", itaue);
       printf(" dt      = %f\n", dt     );
       printf(" phrmin  = %f\n", phrmin );
       printf(" dphr    = %f\n", dphr   );
       }
       */

    /*	Initializing the local variables.			*/
    float	epsmin = 1.0E-30 ;
    float	velerr, velerrmin, velerrmax ;
    dt 	*= 0.001 ;	//unit of dt is second(s).
    velerrmin = 50 ;
    velerrmax = 300 ;

    /*	Initializing the output buffer: pdat[nphr][ntau].	*/
    zero2float(pdat,ntau,nphr);
    int	*itmuteh=NULL;
    itmuteh = alloc1int(ntr);

    int	it, itr, itmute, it0, t0min, t0max ;
    float	t0, offmax, offmin, velmax, velmin, phmax, phmin, offmaxtmp ;

    /*	Search Mute Line.					*/
    for( itr = 0 ; itr < ntr ; ++ itr )
    {
        for( it = 0 ; it < ns ; ++ it )
        {
            if( fabs(tcmp[itr][it]) > epsmin )
            {
                itmuteh[itr] = it ;
                break ;
            }
        }
    }
    itmute  = Min_Int_1d(itmuteh,ntr) ;

    offmax  = Max_Float_1d(offset, ntr) ;
    offmin  = Min_Float_1d(offset, ntr) ;
    velmax  = Max_Float_1d(vrms,   ns) ;
    velmin  = Min_Float_1d(vrms,   ns) ;

    /*	Set the tau-axis range. tau=[t0min, t0max].		*/
    t0min = MAX ( itmute, itaus ) ;
    t0max = MIN ( ns,     itaue ) ;

    int     itau, iphr, iphrmin, iphrmax, itime ;
    float   tau, time, vv2, hr, phr, phr2, phrr, vt0, vvtmp2 ;
    float   hhtmp, hhmax, hhmin, coef, offmin2, offmax2 ;

    printf(" offmax = %f\n", offmax);
    printf(" offmin = %f\n", offmin);
    printf(" velmax = %f\n", velmax);
    printf(" velmin = %f\n", velmin);
    printf(" t0min  = %d\n", t0min);
    printf(" t0max  = %d\n", t0max);
    /*
    */

    /*	Scan the tau axis from t0min to t0max.			*/
    for( it0 = t0min ; it0 < t0max ; it0 += 1 )
    {
        t0     = dt*it0 ;
        vt0    = vrms[it0]*t0 ;
        vv2    = vrms[it0]*vrms[it0];
        velerr = (velerrmax-velerrmin)*(it0-t0min)/(t0max-t0min) + velerrmin ;
        //printf(" it0 = %d velerr = %f \n", it0, velerr);

        /*		Compute the range of Ph-value by Offset[min,max], t0, vrms(t0).	*/
        offmaxtmp = offmax;
        for( itr = 0 ; itr < ntr ; ++ itr )
        {
            time = sqrt(t0*t0 + (offset[itr]*offset[itr])/vv2);
            if( itmuteh[itr] > time )
            {
                offmaxtmp = offset[itr] ;
                break ;
            }
        }
        offmax2 = offmaxtmp*offmaxtmp ;
        offmin2 = offmin*offmin ;
        phmax   = 2.0*offmaxtmp/(vrms[it0]*sqrt(vt0*vt0 + offmax2)) ;
        iphrmax = MIN( (int)(phmax/dphr), nphr ) ;
        phmin   = 2.0*offmin/(vrms[it0]*sqrt(vt0*vt0 + offmin2)) ;
        iphrmin = MIN( (int)(phmin/dphr), iphrmax ) ;

        /*		Scan the p axis from phmin to phmax.			*/
        for( iphr = iphrmin ; iphr < iphrmax ; ++ iphr )
        {
            phr  = phrmin + dphr*iphr ;
            phr2 = phr*phr ;

            /*			Compute the intercept time tau by t0, p, v(t0).		*/
            //if( phr2 > 4.0/vv2 )	continue ;
            tau  = t0*sqrt(1.0-phr2*vv2/4.0) ;
            itau = (int)(tau/dt) + 1 ;
            tau  = dt*itau ;

            /*			Compute the half-offset by tau, p, v(t0).		*/
            vvtmp2  = (vrms[it0]-velerr)*(vrms[it0]-velerr) ;
            hhmin   = phr*vvtmp2*tau/(4.0-phr2*vvtmp2);

            vvtmp2  = (vrms[it0]+velerr)*(vrms[it0]+velerr) ;
            hhmax   = phr*vvtmp2*tau/(4.0-phr2*vvtmp2);

            //if( hhmin > hhmax ) printf("Error! hhmin > hhmax \n");
            for( itr = 0 ; itr < ntr ; ++ itr )
            {
                hr = offset[itr]*0.5 ;
                if((hr>hhmin) && (hr<hhmax))
                {
                    time  = tau + hr*phr ;
                    itime = (int)(time/dt) ;
                    if( itime > 0 && itime < ns )
                        pdat[iphr][itau] += sinc_interpolation(tcmp[itr],ns,dt,time);
                    /*
                       if( hr > hhtmp )
                       coef = cos(0.5*3.14*(hr-hhtmp)/(hhmax-hhtmp)) ;
                       if( hr < hhtmp )
                       coef = cos(0.5*3.14*(hr-hhtmp)/(hhmin-hhtmp)) ;
                       if( itime > 0 && itime < ns )
                       pdat[iphr][itau] += tcmp[itr][itime] ;
                       */
                }
            }
        }
    }

    /*	tau-p stack to P=0.	*/
    int	interp_flag = 1 ;
    int	itaup, iph, iph_tau ;
    float	taup, ph, tau0, tau01 ;
    float	ph_tau, intp_coe, iphmin ;

    iphmin = 0 ;
    if( interp_flag == 1 )
    {
        for( itau = t0min ; itau < t0max ; itau += 1 )
        {
            tau0    = dt*itau ;
            ph_tau  = 2.0/vrms[itau] ;
            iph_tau = ph_tau/dphr ;
            if( iph_tau > nphr )	iph_tau = nphr ;

            for( iph = 1 ; iph < iph_tau ; ++ iph )
            {
                ph   = dphr*iph ;
                taup = tau0*sqrt(1.0-ph*ph/(ph_tau*ph_tau));
                pdat[0][itau] += sinc_interpolation(pdat[iph], ns, dt, taup);
            }

            /*
               pdat[0][itau] /= (iph_tau-1.0);
               for( iph = 1 ; iph < 40 ; ++ iph )
               {
               ph    = dphr*iph ;
               taup  = tau0*sqrt(1.0-ph*ph/(ph_tau*ph_tau));
               itaup = (int)(taup/dt);
               taup  = dt*itaup;
               tau01 = taup/sqrt(1.0-ph*ph/(ph_tau*ph_tau));
               pdat[iph][itaup] = sinc_interpolation(pdat[0], ns, dt, tau01);
               }
               */
        }
    }

    free1int(itmuteh) ;

    return 0 ;
}

int	tauph_transform_3d_to_2d( int ns, int ntau, int nphr, int ntr, int itaus, int itaue,
        float **tcmp, float **pdat, float *vrms, float *offset,
        float dt, float phrmin, float dphr)
    /*      Definition of some parameters in function interface.
            int	ns  ;		// Sampling number of the input CDP gather: tcmp[ntr][ns].
            int	ntau;		// Tau-domain sampling number of the output tau-ph data: pdat[nphr][ntau].
            int	nphr;		// The number of Ph-Rays.
            int	ntr ;		// Trace number within the current CDP gather.
            int	itaus ;		// The first Tau-domain sampling point number of the output tau-ph data: pdat[nphr][ntau], default is 0.
            int	itaue ;		// The last  Tau-domain sampling point number of the output tau-ph data: pdat[nphr][ntau], default is ntau.
            float	**tcmp;		// The current CDP gather is in buffer: tcmp[ntr][ns].
            float	**pdat;		// Tau-Ph transform of the current CDP gather is in buffer: pdat[nphr][ntau].
            float	*vrms ;		// The RMS velocity of the current CDP gather is in buffer: vrms[ns].
            float	*offset;	// The abs offset sqrt[(gx-sx)^2+(gy-sy)^2] of the current CDP gather is in buffer: offset[ntr].
            float	dt;		// Sampling rate of the input CDP gather. Warning: the unit of dt should be mili-second(ms).
            float	phrmin;		// The minimum value of Ph-Rays.
            float	dphr;		// The interval of Ph-Rays.
            */

{
    /*
       int	parameter_flag;
       parameter_flag = 1 ;	// Print the parameters on the screen.
       if ( parameter_flag == 1 )
       {
       printf(" ********  Parameters of 2d tau-ph data volume.  ******** \n");
       printf(" ns    = %d\n", ns   );
       printf(" ntau  = %d\n", ntau );
       printf(" nphr  = %d\n", nphr );
       printf(" ntr   = %d\n", ntr  );
       printf(" itaus = %d\n", itaus);
       printf(" itaue = %d\n", itaue);
       printf(" dt      = %f\n", dt     );
       printf(" phrmin  = %f\n", phrmin );
       printf(" dphr    = %f\n", dphr   );
       }
       */

    /*	Initializing the output buffer: pdat[nphr][ntau].	*/
    zero2float(pdat,ntau,nphr);

    dt 	*= 0.001 ;	//unit of dt is second(s).

    int     itau, iphr, itr, itime ;
    float   tau, time, vv2, hr, phr, phr2, phrr, drr ;
    float   rrtmp, rrmax, rrmin ;

    float	drrs ;		// The minimum aperture(m) of Tau-domain transform.
    float	drre ;		// The maximum aperture(m) of Tau-domain transform.

    drr	= 250.0 ;
    drrs	= 50.0 ;
    drre	= 500.0 ;

    float	itaumute ;
    itaumute = 300.0 ;
    //float offmax=8000.0, pmaxtmp, tmptx;

    for( itau = itaus ; itau < itaue ; ++ itau )
    {
        //printf(" itau      = %d\n", itau     );
        tau = dt*itau ;
        vv2 = vrms[itau]*vrms[itau];
        drr = (drre-drrs)*(itau-itaus)/(itaue-itaus) + drrs ;        //straight line relation
        //if( itau < itaumute ) drr = 50.0 ;
        // drr     = 500.0 ;

        for( iphr = 0 ; iphr < nphr ; ++ iphr )
        {
            phr  = phrmin + dphr*iphr ;
            phr2 = phr*phr ;
            if( phr2 > 4.0/vv2 )    continue ;

            /* fengbo test. */
            /*
               tmptx = tau*0.5/offmax;
               pmaxtmp = 2.0*( sqrt(tmptx*tmptx + 1.0/vv2) - tmptx );
               if( fabs(phr) > 0.9*pmaxtmp )	continue ;
               */

            for( itr = 0 ; itr < ntr ; ++ itr )
            {
                hr      = offset[itr] ;
                phrr    = fabs(phr*hr*0.5) ;
                rrtmp   = 2.0*phr*vv2*tau/(4.0-phr2*vv2);
                rrtmp   = fabs(rrtmp);
                rrmin   = MAX( (rrtmp-drr), 0.0 );
                rrmax   = rrtmp + drr ;

                if((hr>rrmin) && (hr<rrmax))
                {
                    time  = tau + phrr ;
                    itime = (int)(time/dt) ;
                    if( itime > 0 && itime < ns )
                        pdat[iphr][itau] += sinc_interpolation(tcmp[itr],ns,dt,time);
                    /*
                       if( itime > 0 && itime < ns )
                       pdat[iphr][itau] += tcmp[itr][itime] ;
                       */
                }
            }
        }
    }

    /*    Mute Application.    */
    int	itlen ;		// The Tau-domain window length of attenuation area.
    itlen	= 50 ;
    for( iphr = 0 ; iphr < nphr ; ++ iphr )
    {
        for( itau = itaus ; itau < (itaus+itlen) ; ++ itau )
            pdat[iphr][itau] *= exp(5.0*(itau-itaus-itlen)/itlen);
    }

    return 0 ;
}



int	tauph_transform_3d(int nphx, int nphy, int ns, int ntau, int ntr,
        int itaus, int itaue,
        float phxmin, float phymin, float dphx, float dphy, float dt,
        float **tcmp, float *vrms, float **header, float ***pdat)
    /*      Definition of some parameters in function interface.
            int	nphx;		// The number of Ph-Rays in X component.
            int	nphy;		// The number of Ph-Rays in Y component.
            int	ns  ;		// Sampling number of the input CDP gather: tcmp[ntr][ns].
            int	ntau;		// Tau-domain sampling number of the output tau-ph data: pdat[nphr][ntau].
            int	ntr ;		// Trace number within the current CDP gather.
            int	itaus ;		// The first Tau-domain sampling point number of the output tau-ph data: pdat[nphr][ntau], default is 0.
            int	itaue ;		// The last  Tau-domain sampling point number of the output tau-ph data: pdat[nphr][ntau], default is ntau.
            float	**tcmp;		// The current CDP gather is in buffer: tcmp[ntr][ns].
            float	*vrms ;		// The RMS velocity of the current CDP gather is in buffer: vrms[ns].
            float	**header ;	// The Source-Receiver location is in buffer: header[ntr][4].
            float	***pdat;	// Tau-Ph transform of the current CDP gather is in buffer: pdat[nphy][nphx][ntau].
            float	dt;		// Sampling rate of the input CDP gather. Warning: the unit of dt should be mili-second(ms).
            float	phxmin;		// The minimum value of Ph-Rays in X component.
            float	phymin;		// The minimum value of Ph-Rays in Y component.
            float	dphx;		// The interval of Ph-Rays in X component.
            float	dphy;		// The interval of Ph-Rays in Y component.
            */
{
    int     parameter_flag;
    parameter_flag = 1 ;	// Print the parameters on the screen.
    if ( parameter_flag == 1 )
    {
        printf(" ********  Parameters of 3d tau-ph data volume.  ******** \n");
        printf(" nphx  = %d\n", nphx );
        printf(" nphy  = %d\n", nphy );
        printf(" ns    = %d\n", ns   );
        printf(" ntau  = %d\n", ntau );
        printf(" ntr   = %d\n", ntr  );
        printf(" itaus = %d\n", itaus);
        printf(" itaue = %d\n", itaue);
        printf(" dt      = %f\n", dt     );
        printf(" phxmin  = %f\n", phxmin );
        printf(" phymin  = %f\n", phymin );
        printf(" dphx    = %f\n", dphx   );
        printf(" dphy    = %f\n", dphy   );
    }

    /*      Initializing the output buffer: pdat[nphy][nphx][ntau].       */
    zero3float( pdat, ntau, nphx, nphy ) ;
    dt	*= 0.001 ;      //unit of dt is second(s).

    int	itau, iphx, iphy, itr, itime ;
    float	sx, gx, sy, gy, offsetx, offsety ;
    float	tau, time, dxx, dyy, vv2, temp ;
    float	phx, phx2, phy, phy2, phxx, phyy, phr, phr2 ;
    float	xxmin, xxmax, yymin, yymax, xxtmp, yytmp ;

    dxx	= 1000.0 ;
    dyy	= 600.0 ;

    for( itau = itaus ; itau < itaue ; ++ itau )
    {
        //printf(" itau      = %d\n", itau     );
        tau = dt*itau ;
        vv2 = vrms[itau]*vrms[itau] ;

        for( iphy = 0 ; iphy < nphy ; ++ iphy )
        {
            phy  = phymin + dphy*iphy ;
            phy2 = phy*phy ;

            for( iphx = 0 ; iphx < nphx ; ++ iphx )
            {
                phx  = phxmin + dphx*iphx ;
                phx2 = phx*phx ;

                phr2 = phy2 + phx2 ;
                if( phr2 > 4.0/vv2 )	continue ;

                for( itr = 0 ; itr < ntr ; ++ itr )
                {
                    sx = header[itr][0] ;
                    sy = header[itr][1] ;
                    gx = header[itr][2] ;
                    gy = header[itr][3] ;

                    offsetx = gx - sx ;
                    offsety = gy - sy ;
                    phxx = phx*offsetx*0.5 ;
                    phyy = phy*offsety*0.5 ;

                    if( phxx < 0.0 || phyy < 0.0 )	continue ;

                    if( offsetx >= 0.0 && offsety >= 0.0 )
                    {	/* phx>0 && phy>0 ; hx>0 && hy>0 ; */
                        yytmp = 2.0*phy*vv2*tau/(4.0-phy2*vv2);
                        yymin = MAX( (yytmp-dyy), 0.0 );
                        yymax = yytmp+dyy ;
                        xxtmp = 2.0*phx*vv2*tau/(4.0-phx2*vv2);
                        xxmin = MAX( (xxtmp-dxx), 0.0 );
                        xxmax = xxtmp+dxx ;
                    }
                    else if( offsetx >= 0.0 && offsety < 0.0 )
                    {	/* phx>0 && phy<0 ; hx>0 && hy<0 ; */
                        yytmp = 2.0*phy*vv2*tau/(4.0-phy2*vv2);
                        yymin = yytmp-dyy ;
                        yymax = MIN( (yytmp+dyy), 0.0 ) ;
                        xxtmp = 2.0*phx*vv2*tau/(4.0-phx2*vv2);
                        xxmin = MAX( (xxtmp-dxx), 0.0 );
                        xxmax = xxtmp+dxx ;
                    }
                    else if( offsetx < 0.0 && offsety >= 0.0 )
                    {	/* phx<0 && phy>0 ; hx<0 && hy>0 ; */
                        yytmp = 2.0*phy*vv2*tau/(4.0-phy2*vv2);
                        yymin = MAX( (yytmp-dyy), 0.0 );
                        yymax = yytmp+dyy ;
                        xxtmp = 2.0*phx*vv2*tau/(4.0-phx2*vv2);
                        xxmin = xxtmp-dxx ;
                        xxmax = MIN( (xxtmp+dxx), 0.0 ) ;
                    }
                    else if( offsetx < 0.0 && offsety < 0.0 )
                    {	/* phx<0 && phy<0 ; hx<0 && hy<0 ; */
                        yytmp = 2.0*phy*vv2*tau/(4.0-phy2*vv2);
                        yymin = yytmp-dyy ;
                        yymax = MIN( (yytmp+dyy), 0.0 ) ;
                        xxtmp = 2.0*phx*vv2*tau/(4.0-phx2*vv2);
                        xxmin = xxtmp-dxx ;
                        xxmax = MIN( (xxtmp+dxx), 0.0 ) ;
                    }

                    if( (offsetx>xxmin && offsetx<xxmax) &&( offsety>yymin && offsety<yymax) )
                    {
                        time  = tau + phxx + phyy ;
                        itime = (int)(time/dt) ;
                        if( itime > 0 && itime < ns )
                        {
                            temp  = sinc_interpolation(tcmp[itr],ns,dt,time);
                            pdat[iphy][iphx][itau] += temp ;
                        }
                        /*
                           itime = (int)(time/dt) ;
                           if( itime > 0 && itime < ns )
                           pdat[iphy][iphx][itau] += tcmp[itr][itime] ;
                           */
                    }
                }
            }
        }
    }
    return 0;
}

int	tauph_transform_2d( int ns, int ntau, int nphx, int ntr, int itaus, int itaue,
        float **tcmp, float **pdat, float *vrms, float *offset,
        float dt, float phxmin, float dphx)
    /*      Definition of some parameters in function interface.
            int	ns  ;		// Sampling number of the input CDP gather: tcmp[ntr][ns].
            int	ntau;		// Tau-domain sampling number of the output tau-ph data: pdat[nphr][ntau].
            int	nphx;		// The number of Ph-Rays in X component.
            int	ntr ;		// Trace number within the current CDP gather.
            int	itaus ;		// The first Tau-domain sampling point number of the output tau-ph data: pdat[nphr][ntau], default is 0.
            int	itaue ;		// The last  Tau-domain sampling point number of the output tau-ph data: pdat[nphr][ntau], default is ntau.
            float	**tcmp;		// The current CDP gather is in buffer: tcmp[ntr][ns].
            float	**pdat;		// Tau-Ph transform of the current CDP gather is in buffer: pdat[nphr][ntau].
            float	*vrms ;		// The RMS velocity of the current CDP gather is in buffer: vrms[ns].
            float	*offset;	// The abs offset sqrt[(gx-sx)^2+(gy-sy)^2] of the current CDP gather is in buffer: offset[ntr].
            float	dt;		// Sampling rate of the input CDP gather. Warning: the unit of dt should be mili-second(ms).
            float	phxmin;		// The minimum value of Ph-Rays in X component.
            float	dphx;		// The interval of Ph-Rays in X component.
            */
{
    int     parameter_flag;
    parameter_flag = 1 ;	// Print the parameters on the screen.
    if ( parameter_flag == 1 )
    {
        printf(" ********  Parameters of 3d tau-ph data volume.  ******** \n");
        printf(" nphx  = %d\n", nphx );
        printf(" ns    = %d\n", ns   );
        printf(" ntau  = %d\n", ntau );
        printf(" ntr   = %d\n", ntr  );
        printf(" itaus = %d\n", itaus);
        printf(" itaue = %d\n", itaue);
        printf(" dt      = %f\n", dt     );
        printf(" phxmin  = %f\n", phxmin );
        printf(" dphx    = %f\n", dphx   );
    }

    /*      Initializing the output buffer: pdat[nphr][ntau].       */
    zero2float( pdat, ntau, nphx );

    dt	*= 0.001 ;         //unit of dt is second(s).
    itaus	= 0   ;
    itaue	= ntau ;

    int	itau, iphx, itr, itime ;
    float	vrange, xxmin, xxmax, xxtmp ;
    float	tau, time, vv2, pp2, phx, phx2, phxx ;
    float	offx, dxx ;

    float	dxxs ;
    float	dxxe ;

    vrange	= 0.9 ;
    dxxs	= 50.0 ;
    dxxe	= 900.0 ;

    //float offmax=2880, pmaxtmp, tmptx;

    for( itau = itaus ; itau < itaue ; ++ itau )
    {
        //printf(" itau      = %d\n", itau     );
        tau = dt*itau ;
        vv2 = vrms[itau]*vrms[itau];
        dxx = (dxxe-dxxs)*(itau-itaus)/(itaue-itaus) + dxxs ;
        dxx = 300.0 ;

        for( iphx = 0 ; iphx < nphx ; ++ iphx )
        {
            phx  = phxmin + dphx*iphx ;
            phx2 = phx*phx ;
            if( phx2 > 4.0*vrange/vv2 )    continue ;

            /* fengbo test. */
            /*
               tmptx = tau*0.5/offmax;
               pmaxtmp = 2.0*( sqrt(tmptx*tmptx + 1.0/vv2) - tmptx );
               if( fabs(phx) > vrange*pmaxtmp )	continue ;
               */

            for( itr = 0 ; itr < ntr ; ++ itr )
            {
                offx = offset[itr] ;
                phxx = phx*offx*0.5 ;
                if( phxx < 0.0 )  continue ;

                if( offx >= 0.0 )
                {
                    xxtmp = 2.0*phx*vv2*tau/(4.0-phx2*vv2);
                    xxmin = MAX( (xxtmp-dxx), 0.0 );
                    xxmax = xxtmp+dxx ;
                }
                else if( offx < 0.0 )
                {
                    xxtmp = 2.0*phx*vv2*tau/(4.0-phx2*vv2);
                    xxmin = xxtmp-dxx ;
                    xxmax = MIN( (xxtmp+dxx), 0.0 ) ;
                }

                if ( ( offx > xxmin ) && ( offx < xxmax ) )
                {
                    time = tau + phxx ;
                    itime = (int)(time/dt) ;
                    if( ( itime > 0 ) && ( itime < ns ) )
                        pdat[iphx][itau] += sinc_interpolation(tcmp[itr],ns,dt,time);
                }
            }
        }

    }

    return 0 ;

}

int	interpolation_taup_domain( float **in, float **out, int ns, int nph, int iphmin, float *vel, float dt, float dph)
{
    float	taup, ph, tau0, tau01 ;
    int	itau, itaup, iph ;

    float	ph_tau, intp_coe ;
    int	iph_tau, nph_intp ;

    zero2float(out, ns, nph);

    intp_coe = 0.95 ;
    for( itau = 0 ; itau < ns ; ++ itau )
    {
        tau0    = dt*itau ;
        ph_tau  = 2.0/vel[itau] ;
        iph_tau = ph_tau/dph ;
        if( iph_tau > nph )	iph_tau = nph ;
        nph_intp = iphmin + 1.0*(nph-iphmin)*(ns-itau)/ns ;

        for( iph = 0 ; iph < iph_tau-1 ; ++ iph )
        {
            ph   = dph*iph ;
            taup = tau0*sqrt(1.0-ph*ph/(ph_tau*ph_tau));
            out[0][itau] += sinc_interpolation(in[iph], ns, dt, taup);
        }
        out[0][itau] /= (iph_tau-1);
        for( iph = 1 ; iph < MIN(iph_tau-1,nph_intp) ; ++ iph )
        {
            ph    = dph*iph ;
            taup  = tau0*sqrt(1.0-ph*ph/(ph_tau*ph_tau));
            itaup = (int)(taup/dt);
            taup  = dt*itaup;
            tau01 = taup/sqrt(1.0-ph*ph/(ph_tau*ph_tau));
            out[iph][itaup] = sinc_interpolation(out[0], ns, dt, tau01);
        }
    }

    return 0 ;
}

float Max_Float_1d(float *array , int n )
{
    int i ;
    float max;

    max = array[0] ;
    for ( i = 1 ; i < n ; ++i )
    {
        if ( array[i] > max )
        {
            max = array[i] ;
        }
    }

    return max ;
}

float Min_Float_1d(float *array , int n )
{
    int i ;
    float min;

    min = array[0] ;
    for ( i = 1 ; i < n ; ++i )
    {
        if ( array[i] < min )
        {
            min = array[i] ;
        }
    }

    return min ;
}

int Min_Int_1d(int *array , int n )
{
    int i ;
    int min;

    min = array[0] ;
    for ( i = 1 ; i < n ; ++i )
    {
        if ( array[i] < min )
        {
            min = array[i] ;
        }
    }

    return min ;
}

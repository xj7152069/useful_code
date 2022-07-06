#ifndef FBSEGY_SHENGLI_H
#define FBSEGY_SHENGLI_H

/* TYPEDEFS */
typedef struct {    /* segy - trace identification header */

        int tracl;      /* 1-4, trace sequence number within line */

        int tracr;      /* 5-8, trace sequence number within reel */

        int fldr;       /* 9-12, field record number */

    int tracf;  /* 13-16, trace number within field record */

    int ep;         /* 17-20, energy source point number */

    int cdp;    /* 21-24, CDP ensemble number */

    int cdpt;   /* 25-28, trace number within CDP ensemble */

    short trid; /* 29-30, trace identification code:
            1 = seismic data
            2 = dead
            3 = dummy
            4 = time break
            5 = uphole
            6 = sweep
            7 = timing
            8 = water break
            9---, N = optional use (N = 32,767)

            Following are CWP id flags:

             9 = autocorrelation

            10 = Fourier transformed - no packing
                 xr[0],xi[0], ..., xr[N-1],xi[N-1]

            11 = Fourier transformed - unpacked Nyquist
                 xr[0],xi[0],...,xr[N/2],xi[N/2]

            12 = Fourier transformed - packed Nyquist
                 even N:
                 xr[0],xr[N/2],xr[1],xi[1], ...,
                xr[N/2 -1],xi[N/2 -1]
                (note the exceptional second entry)
                 odd N:
                 xr[0],xr[(N-1)/2],xr[1],xi[1], ...,
                xr[(N-1)/2 -1],xi[(N-1)/2 -1],xi[(N-1)/2]
                (note the exceptional second & last entries)

            13 = Complex signal in the time domain
                 xr[0],xi[0], ..., xr[N-1],xi[N-1]

            14 = Fourier transformed - amplitude/phase
                 a[0],p[0], ..., a[N-1],p[N-1]

            15 = Complex time signal - amplitude/phase
                 a[0],p[0], ..., a[N-1],p[N-1]

            16 = Real part of complex trace from 0 to Nyquist

            17 = Imag part of complex trace from 0 to Nyquist

            18 = Amplitude of complex trace from 0 to Nyquist

            19 = Phase of complex trace from 0 to Nyquist

            21 = Wavenumber time domain (k-t)

            22 = Wavenumber frequency (k-omega)

            23 = Envelope of the complex time trace

            24 = Phase of the complex time trace

            25 = Frequency of the complex time trace

            30 = Depth-Range (z-x) traces

            43 = Seismic Data, Vertical Component

            44 = Seismic Data, Horizontal Component 1

            45 = Seismic Data, Horizontal Component 2

            46 = Seismic Data, Radial Component

            47 = Seismic Data, Transverse Component

            101 = Seismic data packed to bytes (by supack1)

            102 = Seismic data packed to 2 bytes (by supack2)
            */

    short nvs;  /* 31-32, number of vertically summed traces (see vscode
               in bhed structure) */

    short nhs;  /* 33-34, number of horizontally summed traces (see vscode
               in bhed structure) */

    short duse; /* 35-36, data use:
                1 = production
                2 = test */

    int offset; /* 37-40, distance from source point to receiver
               group (negative if opposite to direction
               in which the line was shot) */

    int gelev;  /* 41-44, receiver group elevation from sea level
               (above sea level is positive) */

    int selev;  /* 45-48, source elevation from sea level
               (above sea level is positive) */

    int sdepth; /* 49-52, source depth (positive) */

    int gdel;   /* 53-56, datum elevation at receiver group */

    int sdel;   /* 57-60, datum elevation at source */

    //int swdep;    /* 61-64, water depth at source */
    int line;   // line number, only for salt3d-segygather.//

    int gwdep;  /* 65-68, water depth at receiver group */

    short scalel;   /* 69-70, scale factor for previous 7 entries
               with value plus or minus 10 to the
               power 0, 1, 2, 3, or 4 (if positive,
               multiply, if negative divide) */

    short scalco;   /* 71-72, scale factor for next 4 entries
               with value plus or minus 10 to the
               power 0, 1, 2, 3, or 4 (if positive,
               multiply, if negative divide) */

    int  sx;    /* 73-76, X source coordinate */

    int  sy;    /* 77-80, Y source coordinate */

    int  gx;    /* 81-84, X group coordinate */

    int  gy;    /* 85-88, Y group coordinate */

    short counit;   /* 89-90, coordinate units code:
                for previous four entries
                1 = length (meters or feet)
                2 = seconds of arc (in this case, the
                X values are longitude and the Y values
                are latitude, a positive value designates
                the number of seconds east of Greenwich
                or north of the equator */

    short wevel;    /* 91-92, weathering velocity */

    short swevel;   /* 93-94, subweathering velocity */

    short sut;  /* 95-96, uphole time at source */

    short gut;  /* 97-98, uphole time at receiver group */

    short sstat;    /* 99-100, source static correction */

    short gstat;    /* 101-102, group static correction */

    short tstat;    /* 103-104, total static applied */

    short laga; /* 105-106, lag time A, time in ms between end of 240-
               byte trace identification header and time
               break, positive if time break occurs after
               end of header, time break is defined as
               the initiation pulse which maybe recorded
               on an auxiliary trace or as otherwise
               specified by the recording system */

    short lagb; /* 107-108, lag time B, time in ms between the time break
               and the initiation time of the energy source,
               may be positive or negative */

    short delrt;    /* 109-110, delay recording time, time in ms between
               initiation time of energy source and time
               when recording of data samples begins
               (for deep water work if recording does not
               start at zero time) */

    short muts; /* 111-112, mute time--start */

    short mute; /* 113-114, mute time--end */

    unsigned short ns;  /* 115-116, number of samples in this trace */

    unsigned short dt;  /* 117-118, sample interval; in micro-seconds */

    short gain; /* 119-120, gain type of field instruments code:
                1 = fixed
                2 = binary
                3 = floating point
                4 ---- N = optional use */

    short igc;  /* 121-122, instrument gain constant */

    short igi;  /* 123-124, instrument early or initial gain */

    short corr; /* 125-126, correlated:
                1 = no
                2 = yes */

    short sfs;  /* 127-128, sweep frequency at start */

    short sfe;  /* 129-130, sweep frequency at end */

    short slen; /* 131-132, sweep length in ms */

    short styp; /* 133-134, sweep type code:
                1 = linear
                2 = cos-squared
                3 = other */

    short stas; /* 135-136, sweep trace length at start in ms */

    short stae; /* 137-138, sweep trace length at end in ms */

    short tatyp;    /* 139-140, taper type: 1=linear, 2=cos^2, 3=other */

    short afilf;    /* 141-142, alias filter frequency if used */

    short afils;    /* 143-144, alias filter slope */

    short nofilf;   /* 145-146, notch filter frequency if used */

    short nofils;   /* 147-148, notch filter slope */

    short lcf;  /* 149-150, low cut frequency if used */

    short hcf;  /* 151-152, high cut frequncy if used */

    short lcs;  /* 153-154, low cut slope */

    short hcs;  /* 155-156, high cut slope */

    short year; /* 157-158, year data recorded */

    short day;  /* 159-160, day of year */

    short hour; /* 161-162, hour of day (24 hour clock) */

    short minute;   /* 163-164, minute of hour */

    short sec;  /* 165-166, second of minute */

    short timbas;   /* 167-168, time basis code:
                1 = local
                2 = GMT
                3 = other */

    short trwf; /* 169-170, trace weighting factor, defined as 1/2^N
               volts for the least sigificant bit */

    short grnors;   /* 171-172, geophone group number of roll switch
               position one */

    short grnofr;   /* 173-174, geophone group number of trace one within
               original field record */

    short grnlof;   /* 175-176, geophone group number of last trace within
               original field record */

    short gaps; /* 177-178, gap size (total number of groups dropped) */

    short otrav;    /* 179-180, overtravel taper code:
                1 = down (or behind)
                2 = up (or ahead) */

    /* local assignments */
    // Feng bo modified the following keywords for ShengLi data.
    int lineno; // 181-184, the LINE number. //

    int cdpno;  // 185-188, the CMP number. //

    //int   mx; // 189-192, x-coor. of mid-point. //
    int byte189;    // 189-192, x-coor. of mid-point. //

    //int   my; // 193-196, y-coor. of mid-point. //
    int byte193;    // 193-196, y-coor. of mid-point. //

    int byte197;// 197-200, not used. //

    //int   ntr;    /* 201-204, number of traces in a single CMP gather. */
    int byte201;    /* 201-204, number of traces in a single CMP gather. */

    short   unass[18];  /* unassigned--NOTE: last entry causes
               a break in the word alignment, if we REALLY
               want to maintain 240 bytes, the following
               entry should be an odd number of short/UINT2
               OR do the insertion above the "mark" keyword
               entry */
} fbsegy_nonstd;

#endif

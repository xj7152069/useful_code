#include <stdio.h>
#include <string.h>
#include <stdlib.h>

// SU道头 关键字 结构体
typedef struct {
    /* segy - trace identification header */
    int tracl;  /*	Trace sequence number within line numbers continue to increase if the
					same line continues across multiple SEG Y files.
					byte# 1-4	*/
    int tracr;  /*	Trace sequence number within SEG Y file each file starts with trace sequence one
					byte# 5-8	*/
    int fldr;   /*	Original field record number
					byte# 9-12	*/
    int tracf;  /*	Trace number within original field record
					byte# 13-16	*/
    int ep;     /*	energy source point number Used when more than one record occurs
					at the same effective surface location.
					byte# 17-20	*/
    int cdp;    /*	Ensemble number (i.e. CDP, CMP, CRP,...)
					byte# 21-24	*/
    int cdpt;   /*	trace number within the ensemble each ensemble starts with trace number one.
					byte# 25-28	*/
    short trid; /*	trace identification code:
					-1 = Other
					0 = Unknown
					1 = Seismic data
					2 = Dead
					3 = Dummy
					4 = Time break
					5 = Uphole
					6 = Sweep
					7 = Timing
					8 = Water break
					9 = Near-field gun signature
					10 = Far-field gun signature
					11 = Seismic pressure sensor
					12 = Multicomponent seismic sensor - Vertical component
					13 = Multicomponent seismic sensor - Cross-line component
					14 = Multicomponent seismic sensor - in-line component
					15 = Rotated multicomponent seismic sensor - Vertical component
					16 = Rotated multicomponent seismic sensor - Transverse component
					17 = Rotated multicomponent seismic sensor - Radial component
					18 = Vibrator reaction mass
					19 = Vibrator baseplate
					20 = Vibrator estimated ground force
					21 = Vibrator reference
					22 = Time-velocity pairs
					23 ... N = optional use (maximum N = 32,767)
					Following are CWP id flags:
						109 = autocorrelation
						110 = Fourier transformed - no packing xr[0],xi[0], ..., xr[N-1],xi[N-1]
						111 = Fourier transformed - unpacked Nyquist xr[0],xi[0],...,xr[N/2],xi[N/2]
						112 = Fourier transformed - packed Nyquist
					even N:
						xr[0],xr[N/2],xr[1],xi[1],...,xr[N/2 -1],xi[N/2 -1](note the exceptional second entry)
					odd N:
						xr[0],xr[(N-1)/2],xr[1],xi[1],...,xr[(N-1)/2 -1],xi[(N-1)/2 -1],xi[(N-1)/2]
						(note the exceptional second & last entries)
						113 = Complex signal in the time domain xr[0],xi[0], ..., xr[N-1],xi[N-1]
						114 = Fourier transformed - amplitude/phase a[0],p[0], ..., a[N-1],p[N-1]
						115 = Complex time signal - amplitude/phase a[0],p[0], ..., a[N-1],p[N-1]
						116 = Real part of complex trace from 0 to Nyquist
						117 = Imag part of complex trace from 0 to Nyquist
						118 = Amplitude of complex trace from 0 to Nyquist
						119 = Phase of complex trace from 0 to Nyquist
						121 = Wavenumber time domain (k-t)
						122 = Wavenumber frequency (k-omega)
						123 = Envelope of the complex time trace
						124 = Phase of the complex time trace
						125 = Frequency of the complex time trace
						130 = Depth-Range (z-x) traces
						201 = Seismic data packed to bytes (by supack1)
						202 = Seismic data packed to 2 bytes (by supack2)
						byte# 29-30	*/
    short nvs;	/*	Number of vertically summed traces yielding this trace. (1 is one trace,
					2 is two summed traces, etc.)
					byte# 31-32	*/
    short nhs;  /*	Number of horizontally summed traces yielding this trace. (1 is one trace
					2 is two summed traces, etc.)
					byte# 33-34	*/
    short duse; /*	Data use:
					1 = Production
					2 = Test
					byte# 35-36	*/
    int offset; /*	Distance from the center of the source point to the center of the receiver group
					(negative if opposite to direction in which the line was shot).
					byte# 37-40 */
    int gelev;  /*	Receiver group elevation from sea level (all elevations above the Vertical datum are
					positive and below are negative).
					byte# 41-44 */
    int selev;  /*	Surface elevation at source.
					byte# 45-48 */
    int sdepth; /*	Source depth below surface (a positive number).
					byte# 49-52 */
    int gdel;   /*	Datum elevation at receiver group.
					byte# 53-56	*/
    int sdel;   /*	Datum elevation at source.
					byte# 57-60	*/
    int swdep;  /*	Water depth at source.
					byte# 61-64 */
    int gwdep;  /*	Water depth at receiver group.
					byte# 65-68 */
    short scalel;   /*	Scalar to be applied to the previous 7 entries to give the real value.
						Scalar = 1, +10, +100, +1000, +10000.
						If positive, scalar is used as a multiplier, if negative, scalar is used as a divisor.
						byte# 69-70	*/
    short scalco;   /*	Scalar to be applied to the next 4 entries to give the real value.
						Scalar = 1, +10, +100, +1000, +10000.
						If positive, scalar is used as a multiplier, if negative, scalar is used as a divisor.
						byte# 71-72 */
    int  sx;    /*	Source coordinate - X
					byte# 73-76	*/
    int  sy;    /*	Source coordinate - Y
					byte# 77-80 */
    int  gx;    /*	Group coordinate - X
					byte# 81-84 */
    int  gy;    /*	Group coordinate - Y
					byte# 85-88 */
    short counit;   /*	Coordinate units: (for previous 4 entries and for the 7 entries before scalel)
						1 = Length (meters or feet)
						2 = Seconds of arc
						3 = Decimal degrees
						4 = Degrees, minutes, seconds (DMS)
						In case 2, the X values are longitude and the Y values are latitude, a positive value
						designates the number of seconds east of Greenwich or north of the equator
						In case 4, to encode +-DDDMMSS
						counit = +-DDD*10^4 + MM*10^2 + SS, with scalco = 1. To encode +-DDDMMSS.ss
						counit = +-DDD*10^6 + MM*10^4 + SS*10^2 with scalco = -100.
						byte# 89-90	*/
    short wevel;    /*	Weathering velocity.
						byte# 91-92 */
    short swevel;   /*	Subweathering velocity.
						byte# 93-94 */
    short sut;	  	/*	Uphole time at source in milliseconds.
						byte# 95-96	*/
    short gut;  	/*	Uphole time at receiver group in milliseconds.
						byte# 97-98	*/
    short sstat;    /*	Source static correction in milliseconds.
						byte# 99-100	*/
    short gstat;    /*	Group static correction in milliseconds.
						byte# 101-102	*/
    short tstat;    /*	Total static applied  in milliseconds. (Zero if no static has been applied.)
						byte# 103-104	*/
    short laga;		/*	Lag time A, time in ms between end of 240-byte trace identification header and time
						break, positive if time break occurs after end of header, time break is defined as
						the initiation pulse which maybe recorded on an auxiliary trace or as otherwise
						specified by the recording system
						byte# 105-106	*/
    short lagb; 	/*	lag time B, time in ms between the time break and the initiation time of the energy
						source, may be positive or negative
						byte# 107-108	*/
    short delrt;    /*	delay recording time, time in ms between initiation time of energy source and time
						when recording of data samples begins (for deep water work if recording does not
						start at zero time)
						byte# 109-110   */
    short muts;		/*	mute time--start
						byte# 111-112	*/
    short mute;		/*	mute time--end
						byte# 113-114	*/
    unsigned short ns;	/*	number of samples in this trace
							byte# 115-116	*/
    unsigned short dt;  /*	sample interval; in micro-seconds
							byte# 117-118	*/
    short gain;		/*	gain type of field instruments code:
						1 = fixed
						2 = binary
						3 = floating point
						4 ---- N = optional use
						byte# 119-120	*/
    short igc;		/*	instrument gain constant
						byte# 121-122	*/
    short igi;		/*	instrument early or initial gain
						byte# 123-124	*/
    short corr; 	/*	correlated:
						1 = no
						2 = yes
						byte# 125-126	*/
    short sfs;		/*	sweep frequency at start
						byte# 127-128	*/
    short sfe;		/*	sweep frequency at end
						byte# 129-130   */
    short slen;		/*	sweep length in ms
						byte# 131-132	*/
    short styp;		/*	sweep type code:
						1 = linear
						2 = cos-squared
						3 = other
						byte# 133-134	*/
    short stas;		/*	sweep trace length at start in ms
						byte# 135-136   */
    short stae;		/*	sweep trace length at end in ms
						byte# 137-138	*/
    short tatyp;    /*	taper type: 1=linear, 2=cos^2, 3=other
						byte# 139-140	*/
    short afilf;    /*	alias filter frequency if used
						byte# 141-142	*/
    short afils;    /*	alias filter slope
						byte# 143-144	*/
    short nofilf;   /*	notch filter frequency if used
						byte# 145-146	*/
    short nofils;   /*	notch filter slope
						byte# 147-148	*/
    short lcf;		/*	low cut frequency if used
						byte# 149-150	*/
    short hcf;		/*	high cut frequncy if used
						byte# 151-152	*/
    short lcs;		/*	low cut slope
						byte# 153-154	*/
    short hcs;		/*	high cut slope
						byte# 155-156	*/
    short year;		/*	year data recorded
						byte# 157-158	*/
	short day;		/*	day of year
						byte# 159-160	*/
    short hour; 	/*	hour of day (24 hour clock)
						byte# 161-162	*/
    short minute;   /*	minute of hour
						byte# 163-164	*/
    short sec;  	/*	second of minute
						byte# 165-166	*/
    short timbas;   /*	time basis code:
						1 = local
						2 = GMT
						3 = other
						byte# 167-168	*/
    short trwf;		/*	trace weighting factor, defined as 1/2^N volts for the least sigificant bit
						byte# 169-170	*/
    short grnors;   /*	geophone group number of roll switch position one
						byte# 171-172	*/
    short grnofr;   /*	geophone group number of trace one within original field record
						byte# 173-174	*/
    short grnlof;   /*	geophone group number of last trace within original field record
						byte# 175-176	*/
    short gaps; 	/*	gap size (total number of groups dropped)
						byte# 177-178	*/
    short otrav;    /*	overtravel taper code:
						1 = down (or behind)
						2 = up (or ahead)
						byte# 179-180	*/
    /* cwp local assignments */
    float d1;	/*	sample spacing for non-seismic data
					byte# 181-184	*/
    float f1;   /*	first sample location for non-seismic data
					byte# 185-188	*/
    int  iline;   /*	sample spacing between traces
					byte# 189-192	*/
    int xline;   /*	first trace location
					byte# 193-196	*/
    float ungpow;   /*	negative of power used for dynamic range compression
						byte# 197-200	*/
    float unscale;  /*	reciprocal of scaling factor to normalize range
						byte# 201-204	*/
    int ntr;	    /*	number of traces
						byte# 205-208	*/
    short mark;		/*	mark selected traces
						byte# 209-210	*/
    short shortpad; /*	alignment padding
						byte# 211-212	*/
    short unass[14];/*	unassigned--NOTE: last entry causes a break in the word alignment, if we REALLY
						want to maintain 240 bytes, the following entry should be an odd number of
						short/UINT2 OR do the insertion above the "mark" keyword entry
						byte# 213-240	*/
} SUHeaderDefine;

// Segy道头 关键字 结构体
typedef struct{
    // bhed - binary header(plus 3200 char Total 3600 byte)//
    char TextHead[3200];

    int jobid;  /* job identification number */

    int lino;   /* line number (only one line per reel) */

    int reno;   /* reel number */

    short ntrpr;    /* number of data traces per record */

    short nart; /* number of auxiliary traces per record */

    unsigned short hdt; /* sample interval in micro secs for this reel */

    unsigned short dto; /* same for original field recording */

    unsigned short hns; /* number of samples per trace for this reel */

    unsigned short nso; /* same for original field recording */

    short format;   /* data sample format code:
                       1 = floating point, 4 byte (32 bits)
                       2 = fixed point, 4 byte (32 bits)
                       3 = fixed point, 2 byte (16 bits)
                       4 = fixed point w/gain code, 4 byte (32 bits)
                       5 = IEEE floating point, 4 byte (32 bits)
                       8 = two's complement integer, 1 byte (8 bits)
                       */

    short fold; /* CDP fold expected per CDP ensemble */

    short tsort;    /* trace sorting code:
                       1 = as recorded (no sorting)
                       2 = CDP ensemble
                       3 = single fold continuous profile
                       4 = horizontally stacked */

    short vscode;   /* vertical sum code:
                       1 = no sum
                       2 = two sum ...
                       N = N sum (N = 32,767) */

    short hsfs; /* sweep frequency at start */

    short hsfe; /* sweep frequency at end */

    short hslen;    /* sweep length (ms) */

    short hstyp;    /* sweep type code:
                       1 = linear
                       2 = parabolic
                       3 = exponential
                       4 = other */

    short schn; /* trace number of sweep channel */

    short hstas;    /* sweep trace taper length at start if
                       tapered (the taper starts at zero time
                       and is effective for this length) */

    short hstae;    /* sweep trace taper length at end (the ending
                       taper starts at sweep length minus the taper
                       length at end) */

    short htatyp;   /* sweep trace taper type code:
                       1 = linear
                       2 = cos-squared
                       3 = other */

    short hcorr;    /* correlated data traces code:
                       1 = no
                       2 = yes */

    short bgrcv;    /* binary gain recovered code:
                       1 = yes
                       2 = no */

    short rcvm; /* amplitude recovery method code:
                   1 = none
                   2 = spherical divergence
                   3 = AGC
                   4 = other */

    short mfeet;    /* measurement system code:
                       1 = meters
                       2 = feet */

    short polyt;    /* impulse signal polarity code:
                       1 = increase in pressure or upward
                       geophone case movement gives
                       negative number on tape
                       2 = increase in pressure or upward
                       geophone case movement gives
                       positive number on tape */

    short vpol; /* vibratory polarity code:
                   code    seismic signal lags pilot by
                   1   337.5 to  22.5 degrees
                   2    22.5 to  67.5 degrees
                   3    67.5 to 112.5 degrees
                   4   112.5 to 157.5 degrees
                   5   157.5 to 202.5 degrees
                   6   202.5 to 247.5 degrees
                   7   247.5 to 292.5 degrees
                   8   293.5 to 337.5 degrees */

    short hunass[170];  /* unassigned */

}SegYHeaderDefine;

//---------------------------Sub Functions-----------------------------------//
void **alloc2 (size_t n1, size_t n2, size_t size);
float **alloc2float(size_t n1, size_t n2);
void free2float(float **p);
void free2 (void **p);
//---------------------------------------------------------------------------//
void swap_int_bytes (int *tni4);
void swap_float_bytes (float *tnf4);
void swap_short_bytes (short *tni2);
void swap_unsignedshort_bytes (unsigned short *tni2);
//---------------------------------------------------------------------------//
void seisDataBigEnd2LittleEnd(float **dat, int n1, int n2);
void su240ByteBigEnd2LittleEnd(SUHeaderDefine *parSU);
void segy400ByteBigEnd2LittleEnd(SegYHeaderDefine *parSegY);
void IBM2IEEE(float *fromf);
void seisDataIBM2IEEE(float **dat, int n1, int n2);
//---------------------------------------------------------------------------//

//--------------------------main function------------------------------------//
int main()
{
    //=======================================================================//
    // SegY or SU data :  Input & Output
    //=======================================================================//
    int flag_segy = 0;
    // flag_segy = 1 : *.segy
    //           = 0 : *.su

    //=====================//
    // fn
    //=====================//
    char fn_inpath[] = "../foward-8km-abs-50m-radon2.shot.su";
    char fn_outpath[] = "../foward-8km-abs-50m-radon2.shot.su.info";

    if(flag_segy == 1){printf("# Segy data\n");}
    else{printf("# SU data \n");}
    printf("# input  Path : %s\n",fn_inpath);
    printf("# output Path : %s\n",fn_outpath);

    //=====================//
    // open file
    //=====================//
    FILE *fp_in, *fp_out;
    if( (fp_in = fopen( fn_inpath, "rb" )) != NULL){}
    else{printf("ERROR : OPEN INPUT FILE\n");}
    if( (fp_out = fopen( fn_outpath, "w" )) != NULL){}
    else{printf("ERROR : OPEN OUTPUT FILE\n");}

    //=====================//
    // initialize
    //=====================//
    SegYHeaderDefine parSegY;
    SUHeaderDefine parSU;
    memset(&parSegY, 0 ,3600);
    memset(&parSU  , 0 ,240);

    //=====================//
    //read data segy header (3600 byte) + suheader (240 byte) + trace data
    //=====================//
    int flag_endian = 0;
    //           0 : little endian
    //           1 : big endian
    if(flag_segy == 1)
    {
        fread( &parSegY   , sizeof(parSegY),  1, fp_in);

        //for *.segy, judge its format
        //short format;   data sample format code:
        //  1 = floating point, 4 byte (32 bits)
        //  2 = fixed point, 4 byte (32 bits)
        //  3 = fixed point, 2 byte (16 bits)
        //  4 = fixed point w/gain code, 4 byte (32 bits)
        //  5 = IEEE floating point, 4 byte (32 bits)
        //  8 = two's complement integer, 1 byte (8 bits)
        if(parSegY.format == 1)
        {printf("ParSegy.format is 1 : floating point, 4 byte; Little End\n");}
        else if(parSegY.format == 2)
        {printf("ParSegy.format is 2 : fixed point, 4 byte; Little End\n");}
        else if(parSegY.format == 3)
        {printf("ParSegy.format is 3 : fixed point, 2 byte; Little End\n");}
        else if(parSegY.format == 4)
        {printf("ParSegy.format is 4 : fixed point w/gain code, 4 byte; Little End\n");}
        else if(parSegY.format == 5)
        {printf("ParSegy.format is 5 : IEEE floating point, 4 byte; Little End\n");}
        else if(parSegY.format == 8)
        {printf("ParSegy.format is 8 : two's complement integer, 1 byte; Little End\n");}
        else
        {
            // swap
            segy400ByteBigEnd2LittleEnd(&parSegY);

            if(parSegY.format == 1)
            {printf("ParSegy.format is 1 : floating point, 4 byte; BIG End\n"); flag_endian = 1;}
            else if(parSegY.format == 2)
            {printf("ParSegy.format is 2 : fixed point, 4 byte; BIG End\n");flag_endian = 1;}
            else if(parSegY.format == 3)
            {printf("ParSegy.format is 3 : fixed point, 2 byte; BIG End\n");flag_endian = 1;}
            else if(parSegY.format == 4)
            {printf("ParSegy.format is 4 : fixed point w/gain code, 4 byte; BIG End\n");flag_endian = 1;}
            else if(parSegY.format == 5)
            {printf("ParSegy.format is 5 : IEEE floating point, 4 byte; BIG End\n");flag_endian = 1;}
            else if(parSegY.format == 8)
            {printf("ParSegy.format is 8 : two's complement integer, 1 byte; BIG End\n");flag_endian = 1;}
            else{printf("ERROR : CAN NOT GET THE TYPE(FORMAT) OF *.SEGY\n");}
        }

    }
    // test
    if(0)
    {
        /*printf("parSegY.TextHead is  : %s \n",parSegY.TextHead);*/
        printf("    parSegY.format is: %d \n",parSegY.format);
        printf("    parSegY.ntrpr is : %d \n",parSegY.ntrpr);
        printf("    parSegY.hns is   : %d \n",parSegY.hns);
        /*printf("\n");*/
    }

    //alloc memory for dat (cmp or shot gather)
    int nx = 1; // 每次读入的道数

    int icounter = 0;
    while(!feof(fp_in))
    {
        if(icounter % 500 == 0)
        {
            int tracl = parSU.tracl;
            printf(" #Trace No.%d\n",icounter);
        }
        icounter ++;
        //for every trace
        //initialize
        memset(&parSU  , 0 ,240);

        //read
        fread( &parSU     , sizeof(parSU)  ,  1, fp_in);
        int nt = parSU.ns;
        /*printf("%d",nt);*/

        float **dat = alloc2float(nt,nx);

        fread( &dat[0][0], sizeof(float)  , nt, fp_in);

        //大小端转换
        if(flag_endian == 1)
        {
            su240ByteBigEnd2LittleEnd(&parSU);
            for(int i1=0; i1<nt; i1++){swap_float_bytes(&dat[i1][0]);}
        }

        //写出数据1 : 写出 *.su 数据
        if(0 && !feof(fp_in))
        {
            fwrite(&parSU,     sizeof(parSU), 1,  fp_out);
            fwrite(&dat[0][0], sizeof(float), nt, fp_out);
        }

        //写出数据2 ： 输出道头关键字信息
        if(1 && !feof(fp_in))
        {
            float sx     = (parSU.sx);
            float sy     = (parSU.sy);
            float selev  = (parSU.selev);
            float sdepth = (parSU.sdepth);
            float gx     = (parSU.gx);
            float gy     = (parSU.gy);
            float gelev  = (parSU.gelev);
            float offset = (parSU.offset);

            fwrite(&sx ,    sizeof(float) , 1, fp_out);
            fwrite(&sy ,    sizeof(float) , 1, fp_out);
            fwrite(&selev , sizeof(float) , 1, fp_out);
            fwrite(&sdepth ,sizeof(float) , 1, fp_out);
            fwrite(&gx ,    sizeof(float) , 1, fp_out);
            fwrite(&gy ,    sizeof(float) , 1, fp_out);
            fwrite(&gelev , sizeof(float) , 1, fp_out);
            fwrite(&offset ,sizeof(float) , 1, fp_out);
        }

        free2float(dat);

    }
    printf("# Total ntrace is : %d\n", icounter-1);

    //=====================//
    // test for su.txt
    //=====================//
    if(0)
    {
        /*printf("parSU.TextHead is  : %s \n",parSegY.TextHead);*/
        printf("    parSU.ns    is   : %d \n",parSU.ns);
        printf("    parSU.dt    is   : %d \n",parSU.dt);
        /*printf("\n");*/
    }

    //=====================//
    //  *.segy Output
    //=====================//
    fclose(fp_in);
    fclose(fp_out);


    return 0;
}


//---------------------------------------------------------------------------//
/* allocate a 2-d array */
void **alloc2 (size_t n1, size_t n2, size_t size)
{
    size_t i2;
    void **p;

    if ((p=(void**)malloc(n2*sizeof(void*)))==NULL)
        return NULL;

    if ((p[0]=(void*)malloc(n2*n1*size))==NULL) {
        free(p);
        return NULL;
    }
    for (i2=0; i2<n2; i2++)
        p[i2] = (char*)p[0]+size*n1*i2;
    return p;
}

/* allocate a 2-d array of floats */
/*  n1: fast dimension; n2: slow dimension */
float **alloc2float(size_t n1, size_t n2)
{
    return (float**)alloc2(n1,n2,sizeof(float));
}

/* free a 2-d array of floats */
void free2float(float **p)
{
    free2((void**)p);
}

/* free a 2-d array */
void free2 (void **p)
{
    free(p[0]);
    free(p);
}

//---------------------------------------------------------------------------//
void swap_int_bytes (int *tni4)
{
    *tni4 =
        (((*tni4 >> 24) & 0xff) | ((*tni4 & 0xff) << 24) |
         ((*tni4 >> 8) & 0xff00) | ((*tni4 & 0xff00) << 8));
}
//---------------------------------------------------------------------------//
void swap_float_bytes (float *tnf4)
{
    int *tni4=(int *)tnf4;
    *tni4=(((*tni4>>24)&0xff) | ((*tni4&0xff)<<24) |
            ((*tni4>>8)&0xff00) | ((*tni4&0xff00)<<8));
}
//---------------------------------------------------------------------------//
void swap_short_bytes (short *tni2)
{
    *tni2 = (((*tni2 >> 8) & 0xff) | ((*tni2 & 0xff) << 8));
}
//---------------------------------------------------------------------------//
void swap_unsignedshort_bytes (unsigned short *tni2)
{
    *tni2 = (((*tni2 >> 8) & 0xff) | ((*tni2 & 0xff) << 8));
}


//---------------------------------------------------------------------------//
void segy400ByteBigEnd2LittleEnd(SegYHeaderDefine *parSegY)
{
    // Big endina 2 Little endian for SegY 400 Byte
    swap_int_bytes(&(parSegY->jobid));
    swap_int_bytes(&(parSegY->lino));
    swap_int_bytes(&(parSegY->reno));
    swap_short_bytes(&(parSegY->ntrpr));
    swap_short_bytes(&(parSegY->nart));
    swap_unsignedshort_bytes(&(parSegY->hdt));
    swap_unsignedshort_bytes(&(parSegY->dto));
    swap_unsignedshort_bytes(&(parSegY->hns));
    swap_unsignedshort_bytes(&(parSegY->nso));
    swap_short_bytes(&(parSegY->format));
    swap_short_bytes(&(parSegY->fold));
    swap_short_bytes(&(parSegY->tsort));
    swap_short_bytes(&(parSegY->vscode));
    swap_short_bytes(&(parSegY->hsfs));
    swap_short_bytes(&(parSegY->hsfe));
    swap_short_bytes(&(parSegY->hslen));
    swap_short_bytes(&(parSegY->hstyp));
    swap_short_bytes(&(parSegY->schn));
    swap_short_bytes(&(parSegY->hstas));
    swap_short_bytes(&(parSegY->hstae));
    swap_short_bytes(&(parSegY->htatyp));
    swap_short_bytes(&(parSegY->hcorr));
    swap_short_bytes(&(parSegY->bgrcv));
    swap_short_bytes(&(parSegY->rcvm));
    swap_short_bytes(&(parSegY->mfeet));
    swap_short_bytes(&(parSegY->polyt));
    swap_short_bytes(&(parSegY->vpol));
    for(int i1=0; i1<170; i1++)
    {
        swap_short_bytes(&(parSegY->hunass[i1]));
    }
}

//---------------------------------------------------------------------------//
void su240ByteBigEnd2LittleEnd(SUHeaderDefine *parSU)
{
    // Big endina 2 Little endian for SU 240 Byte
    swap_int_bytes(&(parSU->tracl));
    swap_int_bytes(&(parSU->tracr));
    swap_int_bytes(&(parSU->fldr));
    swap_int_bytes(&(parSU->ep));
    swap_int_bytes(&(parSU->cdp));
    swap_int_bytes(&(parSU->cdpt));
    swap_short_bytes(&(parSU->trid));
    swap_short_bytes(&(parSU->nvs));
    swap_short_bytes(&(parSU->nhs));
    swap_short_bytes(&(parSU->duse));
    swap_int_bytes(&(parSU->offset));
    swap_int_bytes(&(parSU->gelev));
    swap_int_bytes(&(parSU->selev));
    swap_int_bytes(&(parSU->sdepth));
    swap_int_bytes(&(parSU->gdel));
    swap_int_bytes(&(parSU->sdel));
    swap_int_bytes(&(parSU->swdep));
    swap_int_bytes(&(parSU->gwdep));
    swap_short_bytes(&(parSU->scalel));
    swap_short_bytes(&(parSU->scalco));
    swap_int_bytes(&(parSU->sx));
    swap_int_bytes(&(parSU->sy));
    swap_int_bytes(&(parSU->gx));
    swap_int_bytes(&(parSU->gy));
    swap_short_bytes(&(parSU->counit));
    swap_short_bytes(&(parSU->wevel));
    swap_short_bytes(&(parSU->swevel));
    swap_short_bytes(&(parSU->sut));
    swap_short_bytes(&(parSU->gut));
    swap_short_bytes(&(parSU->sstat));
    swap_short_bytes(&(parSU->gstat));
    swap_short_bytes(&(parSU->tstat));
    swap_short_bytes(&(parSU->laga));
    swap_short_bytes(&(parSU->lagb));
    swap_short_bytes(&(parSU->delrt));
    swap_short_bytes(&(parSU->muts));
    swap_short_bytes(&(parSU->mute));
    swap_unsignedshort_bytes(&(parSU->ns));
    swap_unsignedshort_bytes(&(parSU->dt));
    swap_short_bytes(&(parSU->gain));
    swap_short_bytes(&(parSU->igc));
    swap_short_bytes(&(parSU->igi));
    swap_short_bytes(&(parSU->corr));
    swap_short_bytes(&(parSU->sfs));
    swap_short_bytes(&(parSU->sfe));
    swap_short_bytes(&(parSU->slen));
    swap_short_bytes(&(parSU->styp));
    swap_short_bytes(&(parSU->stas));
    swap_short_bytes(&(parSU->stae));
    swap_short_bytes(&(parSU->tatyp));
    swap_short_bytes(&(parSU->afilf));
    swap_short_bytes(&(parSU->afils));
    swap_short_bytes(&(parSU->nofilf));
    swap_short_bytes(&(parSU->nofils));
    swap_short_bytes(&(parSU->lcf));
    swap_short_bytes(&(parSU->hcf));
    swap_short_bytes(&(parSU->lcs));
    swap_short_bytes(&(parSU->hcs));
    swap_short_bytes(&(parSU->year));
    swap_short_bytes(&(parSU->day));
    swap_short_bytes(&(parSU->hour));
    swap_short_bytes(&(parSU->minute));
    swap_short_bytes(&(parSU->sec));
    swap_short_bytes(&(parSU->timbas));
    swap_short_bytes(&(parSU->trwf));
    swap_short_bytes(&(parSU->grnors));
    swap_short_bytes(&(parSU->grnofr));
    swap_short_bytes(&(parSU->grnlof));
    swap_short_bytes(&(parSU->gaps));
    swap_short_bytes(&(parSU->otrav));
    swap_float_bytes(&(parSU->d1));
    swap_float_bytes(&(parSU->f1));
    swap_int_bytes(&(parSU->iline));
    swap_int_bytes(&(parSU->xline));
    swap_float_bytes(&(parSU->ungpow));
    swap_float_bytes(&(parSU->unscale));
    swap_int_bytes(&(parSU->ntr));
    swap_short_bytes(&(parSU->shortpad));
    for(int i1=0; i1<14; i1++)
    {
        swap_short_bytes(&(parSU->unass[i1]));
    }
}

//---------------------------------------------------------------------------//
void seisDataBigEnd2LittleEnd(float **dat, int n1, int n2)
{
    // Big endina 2 Little endian for seismic data (single cmp or shot gather)
    for(int i2=0; i2<n2; i2++)
    {
        for(int i1=0; i1<n1; i1++)
        {
            swap_float_bytes(&(dat[i2][i1]));
        }
    }
}
//---------------------------------------------------------------------------//
void IBM2IEEE(float *fromf)
{
    int *from = (int *)fromf;
    // int fconv, fmant, i, t;
    int fconv, fmant, t;
    fconv = from[0];

    if (fconv)
    {
        fmant = 0x00ffffff & fconv;
        t = (int)((0x7f000000 & fconv) >> 22) - 130;
        while (!(fmant & 0x00800000))
        {
            --t;
            fmant <<= 1;
        }
        if (t > 254)
            fconv = (0x80000000 & fconv) | 0x7f7fffff;
        else if (t <= 0)
            fconv = 0;
        else
            fconv = (0x80000000 & fconv) | (t << 23) | (0x007fffff & fmant);
    }
    from[0] = fconv;
}

//---------------------------------------------------------------------------//
void seisDataIBM2IEEE(float **dat, int n1, int n2)
{
    for(int i2=0; i2<n2; i2++)
    {
        for(int i1=0; i1<n1; i1++)
        {
            IBM2IEEE(&(dat[i2][i1]));
        }
    }

}
//---------------------------------------------------------------------------//

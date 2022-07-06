#ifndef CMPFILEINDEX_H
#define CMPFILEINDEX_H

typedef struct  {
    int nline;          // the total line number in the project.
    int lineMin;
    int lineMax;
    int nfoldMax_proj;  // the maximum fold number in the project.
    int cdpMin_proj;
    int cdpMax_proj;
    int ntrace_proj;    // the total trace number in the project.
} ProjCMPInfo;  // statistics infomation of a project (may contain multiple CMP files).

typedef struct  {
    int line;           // the line number of each trace in a given CMP file.
    int cdp;            // the cdp  number of each trace in a given CMP file.
    int sx;
    int sy;
    int gx;
    int gy;
    int offset;         // offset=sqrt((gx-sx)*(gx-sx)+(gy-sy)*(gy-sy));
} cmpFileIndex; // statistics infomation of a CMP file (may contain multiple CMP lines).
// "cmpFileIndex" is created by scanning the header keywords of a given CMP file.

typedef struct  {
    int lineNo;         // the current line No.
    int ncdp_line;      // total cdp number with the current line.
    int ntrace_line;    // the total trace number of the current line.
    int cdpMin_line;    // the minimum cdp number of the current line.
    int cdpMax_line;    // the maximum cdp number of the current line.
    //int   nfoldMax_line;  // the maximum fold number of the current line.
    //int   *ntr_cdp;       // ntrace_cdp[ncdp], containing trace number in each CMP gather.
    // i.e. ntr_local = ntr_cdp[cdpNo-cdpMin];
} CMPLineInfo;  // statistics infomation of a CMP line.
// "CMPLineInfo" is created by scanning "cmpFileIndex" with a given line No. (lineNo)
//
typedef struct  {
    int ifile;          // the file index in multiple CMP files.
    int ntr_cdp;        // the trace number in the current CMP gather.
    long *itr_cdmpfile; // the trace index of the current CMP gather in the CMP file (ifile).
    // i.e. ntr_local = ntr_cdp[cdpNo-cdpMin];
} CMPInfo;  // statistics infomation of a CMP gather.

typedef struct  {
    int ntr_cdp;        // the trace number in the current CMP gather.
    int *sx;
    int *sy;
    int *gx;
    int *gy;
    int *offset;         // offset=sqrt((gx-sx)*(gx-sx)+(gy-sy)*(gy-sy));
    long *itr_cdmpfile;  // the trace index of the current CMP gather in the CMP file (ifile).
    float **gather; // the traces within a CMP gather: gather[ntr_cdp][ns].
    // i.e. ntr_local = ntr_cdp[cdpNo-cdpMin];
} CMPData;  // header infomation and traces in a CMP gather.

#endif  /* CMPFILEINDEX_H */

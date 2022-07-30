#include <string>
#include "segyhdr_getkeybyname2017.h"

void * GetSuKeyAdd(TrHdrStd *hdr,const std::string & key)
{
    return GetSuKeyAdd( hdr, key.c_str());
}

void * GetSuKeyAdd(TrHdrStd *hdr,const char *key)
{
    void *p = NULL;
    if      (strcmp(key,"tracl"   )==0)
    { p = &(hdr->tracl) ;}
    else if (strcmp(key,"tracr"   )==0) { p = &(hdr->tracr);}
    else if (strcmp(key,"fldr"    )==0) { p = &(hdr->fldr) ;}
    else if (strcmp(key,"tracf"   )==0) { p = &(hdr->tracf);}
    else if (strcmp(key,"ep"      )==0) { p = &(hdr->ep)  ;}
    else if (strcmp(key,"cdp"     )==0) { p = &(hdr->cdp) ;}
    else if (strcmp(key,"cdpt"    )==0) { p = &(hdr->cdpt);}
    else if (strcmp(key,"trid"    )==0) { p = &(hdr->trid);}
    else if (strcmp(key,"nvs"     )==0) { p = &(hdr->nvs) ;}
    else if (strcmp(key,"nhs"     )==0) { p = &hdr->nhs   ;}
    else if (strcmp(key,"duse"    )==0) { p = &hdr->duse  ;}
    else if (strcmp(key,"offset"  )==0) { p = &hdr->offset;}
    else if (strcmp(key,"gelev"   )==0) { p = &hdr->gelev ;}
    else if (strcmp(key,"selev"   )==0) { p = &hdr->selev ;}
    else if (strcmp(key,"sdepth"  )==0) { p = &hdr->sdepth;}
    else if (strcmp(key,"gdel"    )==0) { p = &hdr->gdel  ;}
    else if (strcmp(key,"sdel"    )==0) { p = &hdr->sdel  ;}
    else if (strcmp(key,"swdep"   )==0) { p = &hdr->swdep ;}
    else if (strcmp(key,"gwdep"   )==0) { p = &hdr->gwdep ;}
    else if (strcmp(key,"scalel"  )==0) { p = &hdr->scalel;}
    else if (strcmp(key,"scalco"  )==0) { p = &hdr->scalco;}
    else if (strcmp(key,"sx"      )==0) { p = &hdr->sx    ;}
    else if (strcmp(key,"sy"      )==0) { p = &hdr->sy    ;}
    else if (strcmp(key,"gx"      )==0) { p = &hdr->gx    ;}
    else if (strcmp(key,"gy"      )==0) { p = &hdr->gy    ;}
    else if (strcmp(key,"counit"  )==0) { p = &hdr->counit;}
    else if (strcmp(key,"wevel"   )==0) { p = &hdr->wevel ;}
    else if (strcmp(key,"swevel"  )==0) { p = &hdr->swevel;}
    else if (strcmp(key,"sut"     )==0) { p = &hdr->sut   ;}
    else if (strcmp(key,"gut"     )==0) { p = &hdr->gut   ;}
    else if (strcmp(key,"sstat"   )==0) { p = &hdr->sstat ;}
    else if (strcmp(key,"gstat"   )==0) { p = &hdr->gstat ;}
    else if (strcmp(key,"tstat"   )==0) { p = &hdr->tstat ;}
    else if (strcmp(key,"laga"    )==0) { p = &hdr->laga  ;}
    else if (strcmp(key,"lagb"    )==0) { p = &hdr->lagb  ;}
    else if (strcmp(key,"delrt"   )==0) { p = &hdr->delrt ;}
    else if (strcmp(key,"muts"    )==0) { p = &hdr->muts  ;}
    else if (strcmp(key,"mute"    )==0) { p = &hdr->mute  ;}
    else if (strcmp(key,"ns"      )==0) { p = &hdr->ns    ;}
    else if (strcmp(key,"dt"      )==0) { p = &hdr->dt    ;}
    else if (strcmp(key,"gain"    )==0) { p = &hdr->gain  ;}
    else if (strcmp(key,"igc"     )==0) { p = &hdr->igc   ;}
    else if (strcmp(key,"igi"     )==0) { p = &hdr->igi   ;}
    else if (strcmp(key,"corr"    )==0) { p = &hdr->corr  ;}
    else if (strcmp(key,"sfs"     )==0) { p = &hdr->sfs   ;}
    else if (strcmp(key,"sfe"     )==0) { p = &hdr->sfe   ;}
    else if (strcmp(key,"slen"    )==0) { p = &hdr->slen  ;}
    else if (strcmp(key,"styp"    )==0) { p = &hdr->styp  ;}
    else if (strcmp(key,"stas"    )==0) { p = &hdr->stas  ;}
    else if (strcmp(key,"stae"    )==0) { p = &hdr->stae  ;}
    else if (strcmp(key,"tatyp"   )==0) { p = &hdr->tatyp ;}
    else if (strcmp(key,"afilf"   )==0) { p = &hdr->afilf ;}
    else if (strcmp(key,"afils"   )==0) { p = &hdr->afils ;}
    else if (strcmp(key,"nofilf"  )==0) { p = &hdr->nofilf;}
    else if (strcmp(key,"nofils"  )==0) { p = &hdr->nofils;}
    else if (strcmp(key,"lcf"     )==0) { p = &hdr->lcf   ;}
    else if (strcmp(key,"hcf"     )==0) { p = &hdr->hcf   ;}
    else if (strcmp(key,"lcs"     )==0) { p = &hdr->lcs   ;}
    else if (strcmp(key,"hcs"     )==0) { p = &hdr->hcs   ;}
    else if (strcmp(key,"year"    )==0) { p = &hdr->year  ;}
    else if (strcmp(key,"day"     )==0) { p = &hdr->day   ;}
    else if (strcmp(key,"hour"    )==0) { p = &hdr->hour  ;}
    else if (strcmp(key,"minute"  )==0) { p = &hdr->minute;}
    else if (strcmp(key,"sec"     )==0) { p = &hdr->sec   ;}
    else if (strcmp(key,"timbas"  )==0) { p = &hdr->timbas;}
    else if (strcmp(key,"trwf"    )==0) { p = &hdr->trwf  ;}
    else if (strcmp(key,"grnors"  )==0) { p = &hdr->grnors;}
    else if (strcmp(key,"grnofr"  )==0) { p = &hdr->grnofr;}
    else if (strcmp(key,"grnlof"  )==0) { p = &hdr->grnlof;}
    else if (strcmp(key,"gaps"    )==0) { p = &hdr->gaps  ;}
    else if (strcmp(key,"otrav"   )==0) { p = &hdr->otrav ;}

    //else if (strcmp(key,"cmpline" )==0)
    //{ *p = &hdr->cmpline ;strcpy(Atype,"int");}

    //else if (strcmp(key,"indxshot" )==0)
    //{ *p = &hdr->indxshot;strcpy(Atype,"int");}

    //else if (strcmp(key,"recline" )==0)
    //{ *p = &hdr->recline ;strcpy(Atype,"int");}


    //////else if (strcmp(key,"d1"      )==0) { *p = &hdr->d1    ;strcpy(Atype,"float");}
    //////else if (strcmp(key,"f1"      )==0) { *p = &hdr->f1    ;strcpy(Atype,"float");}
    //////else if (strcmp(key,"d2"      )==0) { *p = &hdr->d2    ;strcpy(Atype,"float");}
    //////else if (strcmp(key,"f2"      )==0) { *p = &hdr->f2    ;strcpy(Atype,"float");}


    else
    {
        p=NULL;
        printf("%s ERROR! key=%s\n",__func__, key);
    }

    return p;
}

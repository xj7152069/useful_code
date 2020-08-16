/***************************************
*     Author:luofei
*     Date:2016-12-20 15:01
*     Filename:openCfile.c
*     Description:
*
*     Last modified:2016-12-20 15:01
****************************************/
#include"_lf.h"
#include"_lfFunC.h"

FILE *inC1,*inC2;

int cfile2_(float *in1,float *in2,int *nx, int *nz,int *itr,int *status)
{
	int nvx,nvz,itrace,len;
	char fn1[256],fn2[256],fn3[256],fn4[256];

	char par[256];
	FILE *fp;

	if(*status==1){
		if((fp=fopen("InvTomo.par","r"))==NULL)
		{
			printf("Open Parameter File Error!\n");
			return 0;
		}

		fgets(par,256,fp);
		fgets(par,256,fp);
		fgets(par,256,fp);
		fgets(par,256,fp);
		fgets(par,256,fp);
		fgets(par,256,fp);
		fgets(par,256,fp);
		fgets(fn1,256,fp);
		fclose(fp);
		if(fn1[strlen(fn1)-1]=='\n')
		{
			fn1[strlen(fn1)-1]='\0';
		}

		strcpy(fn2,fn1);

		strcat(fn1,"Ray.dat");
		strcat(fn2,"Fresnel.dat");
	}

//	strcpy(fn1,"/data3/luofei/FatRay/SR/Irregular/HaShan/Ray/Ray.dat");
//	strcpy(fn2,"/data3/luofei/FatRay/SR/Irregular/HaShan/Ray/Fresnel.dat");
//	printf("%s\n",fn1);fflush(stdout);	
//	printf("%d\n",*status);fflush(stdout);
	nvx    = *nx;
	nvz    = *nz;
	itrace = *itr;
	len    =(itrace-1)*nvx*nvz*4;
//	printf("itrace=%d len=%d\n",itrace,len);fflush(stdout);

	switch(*status){
		case 1:			
			inC1=fopen(fn1,"wb");
			inC2=fopen(fn2,"wb");
			break;
		case 2:
			fclose(inC1);
			fclose(inC2);
			break;
		case 3:
			fseek(inC1,len,SEEK_SET);
			fseek(inC2,len,SEEK_SET);
			fwrite(&in1[0],sizeof(float),nvx*nvz,inC1);
			fwrite(&in2[0],sizeof(float),nvx*nvz,inC2);
			break;
		default:
			break;
	}
	return 0;
}

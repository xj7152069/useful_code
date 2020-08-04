#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#include "/home/wcl/include/GeoSubroutineC/alloc_wcl.c"
#include "segy.h"

int main(int argv, char *argc)
{

	int i, iit;
	int nx, lt, nshot;
	int ix, it, ishot;

	FILE *fin, *fout;
	float **image;
	float ***image_tmp;
	float **image_stack;

	float *trace, *trace_hf;
	segy *su;
	segy head;
	
	nx=301;
	lt=2500;
	nshot=201;

/*
	image_tmp=alloc3float(nz, nx, nshot);
	zero3float(image_tmp, nz, nx, nshot);

	image=alloc2float(nz, nshot);
	zero2float(image, nz, nshot);

	image_stack=alloc2float(nz, nx);
	zero2float(image_stack, nz, nx);
	su=(segy *)malloc(sizeof(segy)* ntrace);
*/
	trace=alloc1float(lt);
	trace_hf=alloc1float(lt);

	srand((unsigned)time(NULL));

	for(i=0;i<10;i++)
		printf("%d, %d\n", i, rand()%100+1);



	fin=fopen("/data3/wcl/waxian_model/data/shot_waxian_mod_201s_301r_mute.su", "rb");
	fout=fopen("/data3/wcl/waxian_model/data/shot_waxian_mod_201s_301r_mute_wcl.su", "wb");

	for(ishot=0; ishot<nshot;ishot++)
		for(ix=0;ix<nx;ix++)
		{
			fread(&head, sizeof(segy), 1, fin);
			fread(&trace[0], sizeof(float), lt, fin);

			iit=rand()%25;
			printf("ishot=%d, ix=%d, iit=%d\n", ishot, ix, iit);
			for(it=0;it<lt;it++)
			{
				if(it+iit >=0 && it+iit<lt)
				{
					trace_hf[it]=trace[it+iit];
				}
			}

			fwrite(&head, sizeof(segy), 1, fout);
			fwrite(&trace_hf[0], sizeof(float), lt, fout);
		}
	fclose(fin);
	fclose(fout);
	
	free1float(trace);
	free1float(trace_hf);


/*
	fin=fopen("/data3/wcl/waxian_model/result/rtm/image_test111", "rb");
	for(ishot=0; ishot<nshot;ishot++)
		fread(&image_tmp[ishot][0][0], sizeof(float), nz*nx, fin);
	fclose(fin);

	ix=700;
	for(ishot=0; ishot<nshot;ishot++)
		for(iz=0;iz<nz;iz++)
			image[ishot][iz]=image_tmp[ishot][ix][iz];


	for(ix=0;ix<nx;ix++)
		for(iz=0;iz<nz;iz++)
			for(ishot=0;ishot<nshot;ishot++)
				image_stack[ix][iz]+=image_tmp[ishot][ix][iz];

	fout=fopen("/data3/wcl/waxian_model/result/rtm/image_test111_adjust", "wb");
	for(ix=0;ix<nx;ix++)
		for(ishot=0; ishot<nshot;ishot++)
			for(iz=0;iz<nz;iz++)
				fwrite(&image_tmp[ishot][ix][iz], sizeof(float), 1, fout);

	fclose(fout);

	fout=fopen("/data3/wcl/waxian_model/result/rtm/image_test111_700", "wb");
		for(ishot=0; ishot<nshot;ishot++)
			for(iz=0;iz<nz;iz++)
				fwrite(&image[ishot][iz], sizeof(float), 1, fout);
	fclose(fout);

	fout=fopen("/data3/wcl/waxian_model/result/rtm/image_test111_stack", "wb");
	for(ix=0;ix<nx;ix++)
		for(iz=0;iz<nz;iz++)
			fwrite(&image_stack[ix][iz], sizeof(float), 1, fout);
	fclose(fout);

	free3float(image_tmp);
	free2float(image);
	free2float(image_stack);
*/		

}





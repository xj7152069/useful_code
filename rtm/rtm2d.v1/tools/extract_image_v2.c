#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "/home/wcl/include/GeoSubroutineC/alloc_wcl.c"

int main(int argv, char *argc)
{


	int nx, nz, nshot;
	int ix, iz, ishot;

	FILE *fin, *fout;
	float **image;
	float ***image_tmp;
	float **image_stack;
	
//	nx=2579;
//	nz=795;
//	nshot=226;

	nx=2500;
	nz=301;
	nshot=201;

	image_tmp=alloc3float(nz, nx, nshot);
	zero3float(image_tmp, nz, nx, nshot);

	image=alloc2float(nz, nshot);
	zero2float(image, nz, nshot);


//	fin=fopen("/data3/wcl/waxian_model/result/rtm/image_test_wrong_data_sm10_adjust_ori_gather", "rb");
	fin=fopen("/data3/wcl/waxian_model/data/shot_waxian_mod_201s_301r_mute_wcl.dat", "rb");
	for(ishot=0; ishot<nshot;ishot++)
		for(ix=0;ix<nx;ix++)
			fread(&image_tmp[ishot][ix][0], sizeof(float), nz, fin);
	fclose(fin);

	ishot=0;


//	fout=fopen("/data3/wcl/waxian_model/result/rtm/image_test_wrong_data_sm10_adjust_ori_gather_900", "wb");
	fout=fopen("/data3/wcl/waxian_model/data/shot_waxian_mod_201s_301r_mute_wcl.dat_shot0", "wb");
		for(ix=0;ix<nx;ix++)
			for(iz=0;iz<nz;iz++)
				fwrite(&image_tmp[ishot][ix][iz], sizeof(float), 1, fout);
	fclose(fout);

/*	
	fout=fopen("/data3/wcl/waxian_model/result/rtm/image_test_wrong_datalap_stack", "wb");
	for(ix=0;ix<nx;ix++)
		for(iz=0;iz<nz;iz++)
			fwrite(&image_stack[ix][iz], sizeof(float), 1, fout);
	fclose(fout);
*/

	free3float(image_tmp);
	free2float(image);
		

}





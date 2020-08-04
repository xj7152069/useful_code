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

	nx=1801;
	nz=301;
	nshot=201;

	image_tmp=alloc3float(nz, nx, nshot);
	zero3float(image_tmp, nz, nx, nshot);

	image=alloc2float(nz, nshot);
	zero2float(image, nz, nshot);

	image_stack=alloc2float(nz, nx);
	zero2float(image_stack, nz, nx);



	fin=fopen("/data3/wcl/waxian_model/result/rtm/yzhou_image_test_wrong_data_sm10", "rb");
//	fin=fopen("/data3/wcl/nonlinear_beamforming/beijing/eage_image_process_", "rb");
	for(ishot=0; ishot<nshot;ishot++)
		fread(&image_tmp[ishot][0][0], sizeof(float), nz*nx, fin);
	fclose(fin);

	ix=1200;
	for(ishot=0; ishot<nshot;ishot++)
		for(iz=0;iz<nz;iz++)
			image[ishot][iz]=image_tmp[ishot][ix][iz];


	for(ix=0;ix<nx;ix++)
		for(iz=0;iz<nz;iz++)
			for(ishot=0;ishot<nshot;ishot++)
				image_stack[ix][iz]+=image_tmp[ishot][ix][iz];


	fout=fopen("/data3/wcl/waxian_model/result/rtm/yzhou_image_test_wrong_data_sm10_adjust", "wb");
	for(ix=0;ix<nx;ix++)
		for(ishot=0; ishot<nshot;ishot++)
			for(iz=0;iz<nz;iz++)
				fwrite(&image_tmp[ishot][ix][iz], sizeof(float), 1, fout);

	fclose(fout);


	fout=fopen("/data3/wcl/waxian_model/result/rtm/yzhou_image_test_wrong_data_sm10_1200", "wb");
		for(ishot=0; ishot<nshot;ishot++)
			for(iz=0;iz<nz;iz++)
				fwrite(&image[ishot][iz], sizeof(float), 1, fout);
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
	free2float(image_stack);
		

}





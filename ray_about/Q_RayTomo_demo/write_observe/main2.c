#include <stdio.h>
#include <stdlib.h>

#include "/home/wangpy/wpy/cwp_su/alloc.c"
#include "/home/wangpy/wpy/cwp_su/alloc.h"
#include "/home/wangpy/wpy/cwp_su/complex.c"
#include "/home/wangpy/wpy/cwp_su/complex.h"

int readpar( int *nx, int *ny, int *nz, int *nline, int *ns, int *nr,
		int *linex0, int *liney0, int *linez0, int *sx0, int *sy0, int *sz0, int *ofx0, int *ofy0, int *ofz0, 
		float *dx, float *dy, float *dz, float *dsx, float *dsy, float *dsz, float *dlinex, float *dliney, float *dlinez, 
		float *dofx, float *dofy, float *dofz, char *observefile, char *tobsfile, char *parafile);

int main(int argc, char *argv[])
{
	if(2 < argc)
	{
		printf("Parameter file ?\n");
		return 1;
	}

	int i, j, k, is, ir, iline, ilayer, zlayer;
	int nt, nx, ny, nz, nline, ns, nr, nof, nlayer;
	int flag_rp;
	int linex, liney, linez, sx, sy, sz, ofx, ofy, ofz, rx, ry, rz;
	int temp = 0;

	int linex0, liney0, linez0, sx0, sy0, sz0, ofx0, ofy0, ofz0;

	float dt, dx, dy, dz, dsx, dsy, dsz, dlinex, dliney, dlinez, dofx, dofy, dofz;
	float tobs;

	char observefile[256], tobsfile[256];
	
	FILE *fp_obf, *fp_tobs;

	flag_rp = readpar(&nx, &ny, &nz, &nline, &ns, &nr,
					&linex0, &liney0, &linez0, &sx0, &sy0, &sz0, &ofx0, &ofy0, &ofz0,
					&dx, &dy, &dz, &dsx, &dsy, &dsz, &dlinex, &dliney, &dlinez, 
					&dofx, &dofy, &dofz, observefile, tobsfile, argv[1]);
	
	fp_obf = fopen(observefile, "w");
	fp_tobs = fopen(tobsfile,"r");
	
	nlayer = 3;

	for(iline=0; iline<nline; iline++)
	{
		linex = linex0 + iline*dlinex;
		liney = liney0 + iline*dliney;
		linez = liney0 + iline*dlinez;

		fprintf(fp_obf, "LINE		%d		%d		%d\n", iline+1, linex, linez);
		
		for(ilayer=0; ilayer<nlayer; ilayer++)
		{
			if(ilayer==0)
				zlayer = 240*dz;
			if(ilayer==1)
				zlayer = 540*dz;
			if(ilayer==2)
				zlayer = 759*dz;

			for(is=0; is<ns; is++)
			{
				sx = sx0 + is*dsx;
				sy = sy0 + is*dsy;
				sz = sz0 + zlayer + is*dsz;

				fprintf(fp_obf, "SHOT		%d		%d		%d		%d		%d		%d		%d\n", ilayer*ns+is+1, sx, sy, sz, temp, temp, temp);

				for(ir=0; ir<nr; ir++)
				{
					ofx = ofx0 + ir*dofx;	
					ofy = ofy0 + ir*dofy;	
					ofz = ofz0 + (-zlayer) + ir*dofz;	
					
					rx = sx + ofx;
					ry = sy + ofy;
					rz = sz + ofz;

					fread(&tobs, sizeof(float), 1, fp_tobs);
//					tobs = 0;

					fprintf(fp_obf, "TRACE		%d		%d		%d		%d		%f\n", ilayer*ns*nr+is*nr+ir+1, rx, ry, rz, tobs*1000);
				}
			}
		}
	}

	fclose(fp_obf);
	fclose(fp_tobs);

	return 0;

}

int readpar( int *nx, int *ny, int *nz, int *nline, int *ns, int *nr,
		int *linex0, int *liney0, int *linez0, int *sx0, int *sy0, int *sz0, int *ofx0, int *ofy0, int *ofz0, 
		float *dx, float *dy, float *dz, float *dsx, float *dsy, float *dsz, float *dlinex, float *dliney, float *dlinez, 
		float *dofx, float *dofy, float *dofz, char *observefile, char *tobsfile, char *parafile)
{
	int i, j, k;
	char par_name[100];

	FILE *fp;
	fp = fopen(parafile, "r");
	
	fscanf(fp, "%s%s", par_name, observefile);
	fscanf(fp, "%s%s", par_name, tobsfile);

	fscanf(fp, "%s%d", par_name, nx);
	fscanf(fp, "%s%d", par_name, ny);
	fscanf(fp, "%s%d", par_name, nz);
	fscanf(fp, "%s%d", par_name, nline);
	fscanf(fp, "%s%d", par_name, linex0);
	fscanf(fp, "%s%d", par_name, liney0);
	fscanf(fp, "%s%d", par_name, linez0);
	fscanf(fp, "%s%d", par_name, ns);
	fscanf(fp, "%s%d", par_name, sx0);
	fscanf(fp, "%s%d", par_name, sy0);
	fscanf(fp, "%s%d", par_name, sz0);
	fscanf(fp, "%s%d", par_name, nr);
	fscanf(fp, "%s%d", par_name, ofx0);
	fscanf(fp, "%s%d", par_name, ofy0);
	fscanf(fp, "%s%d", par_name, ofz0);
	
	fscanf(fp, "%s%f", par_name, dx);
	fscanf(fp, "%s%f", par_name, dy);
	fscanf(fp, "%s%f", par_name, dz);
	fscanf(fp, "%s%f", par_name, dsx);
	fscanf(fp, "%s%f", par_name, dsy);
	fscanf(fp, "%s%f", par_name, dsz);
	fscanf(fp, "%s%f", par_name, dlinex);
	fscanf(fp, "%s%f", par_name, dliney);
	fscanf(fp, "%s%f", par_name, dlinez);
	fscanf(fp, "%s%f", par_name, dofx);
	fscanf(fp, "%s%f", par_name, dofy);
	fscanf(fp, "%s%f", par_name, dofz);
	
	fclose(fp);

	
	printf("Observefile = %s\n", observefile);
	printf("tobsfile = %s\n", tobsfile);
	printf("nx = %d\n", *nx);
	printf("ny = %d\n", *ny);
	printf("nz = %d\n", *nz);
	printf("nline = %d\n", *nline);
	printf("linex0 = %d\n", *linex0);
	printf("liney0 = %d\n", *liney0);
	printf("linez0 = %d\n", *linez0);
	printf("ns = %d\n", *ns);
	printf("sx0 = %d\n", *sx0);
	printf("sy0 = %d\n", *sy0);
	printf("sz0 = %d\n", *sz0);
	printf("nr = %d\n", *nr);
	printf("ofx0 = %d\n", *ofx0);
	printf("ofy0 = %d\n", *ofy0);
	printf("ofz0 = %d\n", *ofz0);
	printf("dx = %f\n", *dx);
	printf("dy = %f\n", *dy);
	printf("dz = %f\n", *dz);
	printf("dsx = %f\n", *dsx);
	printf("dsy = %f\n", *dsy);
	printf("dsz = %f\n", *dsz);
	printf("dlinex = %f\n", *dlinex);
	printf("dliney = %f\n", *dliney);
	printf("dlinez = %f\n", *dlinez);
	printf("dofx = %f\n", *dofx);
	printf("dofy = %f\n", *dofy);
	printf("dofz = %f\n", *dofz);
	
	return 0;
}





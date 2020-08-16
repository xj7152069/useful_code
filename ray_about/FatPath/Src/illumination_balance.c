#include "stdio.h"
#include "stdlib.h"
#include "string.h"
#include "math.h"

#include "include/alloc.h"

#include "./lsqr_packet_stanford/lsqr.c"
#include "./lsqr_packet_stanford/lsqr_mult_2.c"
#include "include/illumination_balance.c"

void illumination_(int *Nvx, int *Nvz, int *Nzz, int *N_rows, int *N_cols, 
					int *indx, int *jndx, double *path, double *rhs_vec)
{

	int nn, nvx, nvz, nzz, n_rows, n_cols, i_row, i_col, i;
	long tmpi, tmpj, t_val;
	float tmpp;

	nvx=*Nvx;
	nvz=*Nvz;
	nzz=*Nzz;
	n_rows=*N_rows;
	n_cols=*N_cols;

	/*-------------------------------------------------------------------------------------------*/
	//For illumination balance
	//ib ---- illumination balance
	int *i_val_ib;
	float *rhs_vec_ib;
	//For LSQR to solving illumination balance inversion problem
	lsqr_input *in_lsqr_ib;
	lsqr_output *out_lsqr_ib;
	lsqr_work *work_lsqr_ib;
	lsqr_func *func_lsqr_ib;
	matrix_a *prod_lsqr;
	matrix_a *prod_lsqr_ib;
	int n_rows_ib, n_cols_ib;
	float damp_val_ib = 0.0;
	float cond_lim_ib = 6000000.0;
	FILE *fp_lsqr_ib;
	char lsqr_print_FN_ib[256]="lsqr_ib.txt";

	/*-------------------------------------------------------------------------------------------*/

//	nn = 2 * (nvx + nvz);

	/*-------------------------------------------------------------------------------------------*/
	//For LSQR
	prod_lsqr = (matrix_a *) malloc (sizeof(matrix_a));
	prod_lsqr->n_rows = (long)n_rows;
	prod_lsqr->n_cols = (long)n_cols;

	prod_lsqr->m_row = alloc1matrix_row (prod_lsqr->n_rows);
	for (i_row=0; i_row<prod_lsqr->n_rows; i_row++)
	{
		prod_lsqr->m_row[i_row].n_val = 0;
//		prod_lsqr->m_row[i_row].m_col = alloc1matrix_col (nn);
	}

	//debug
	printf ("n_rows = %ld, n_cols = %ld\n", prod_lsqr->n_rows, prod_lsqr->n_cols);

	/*-------------------------------------------------------------------------------------------*/
	//Forward Convert Matrix
	for (i=0;i<nzz;i++)
	{
		tmpi=(long)indx[i]-1;
		prod_lsqr->m_row[tmpi].n_val++;
	}

	long *i_val;
	i_val = (long*)alloc1(n_rows,sizeof(long));

	for (i_row=0; i_row<prod_lsqr->n_rows; i_row++)
	{
		i_val[i_row]=0;
		nn = prod_lsqr->m_row[i_row].n_val;
		prod_lsqr->m_row[i_row].m_col = alloc1matrix_col(nn);
	}

	for (i=0;i<nzz;i++)
	{
		tmpi=(long)indx[i]-1;
		tmpj=(long)jndx[i]-1;
		tmpp=(float)path[i];

		t_val=i_val[tmpi];
		i_val[tmpi]++;

		prod_lsqr->m_row[tmpi].m_col[t_val].j_col = tmpj;
		prod_lsqr->m_row[tmpi].m_col[t_val].val = tmpp;

		indx[i]=0;
		jndx[i]=0;
		path[i]=0.0;

	}

	free(i_val);

	/*-------------------------------------------------------------------------------------------*/
	//For LSQR in illumination balance
	n_rows_ib = n_cols + 1;
	n_cols_ib = n_rows + 1;

	alloc_lsqr_mem (&in_lsqr_ib, &out_lsqr_ib, &work_lsqr_ib, &func_lsqr_ib, (long)n_rows_ib, (long)n_cols_ib);

	func_lsqr_ib->mat_vec_prod = mtx_mult;

	prod_lsqr_ib = (matrix_a *) malloc (sizeof(matrix_a));
	prod_lsqr_ib->n_rows = (long)n_rows_ib;
	prod_lsqr_ib->n_cols = (long)n_cols_ib;

	prod_lsqr_ib->m_row = alloc1matrix_row (prod_lsqr_ib->n_rows);
	for (i_row=0; i_row<prod_lsqr_ib->n_rows; i_row++)
	{
		prod_lsqr_ib->m_row[i_row].n_val = 0;
//		prod_lsqr_ib->m_row[i_row].m_col = alloc1matrix_col (nn);
	}

	//debug
	printf ("n_rows_ib = %ld, n_cols_ib = %ld\n", prod_lsqr_ib->n_rows, prod_lsqr_ib->n_cols);

	in_lsqr_ib->num_rows = (long)n_rows_ib;
	in_lsqr_ib->num_cols = (long)n_cols_ib;
	in_lsqr_ib->damp_val = (double)damp_val_ib;
	in_lsqr_ib->rel_mat_err = 5e-3;
	in_lsqr_ib->rel_rhs_err = 5e-3;
	in_lsqr_ib->cond_lim = cond_lim_ib;
	in_lsqr_ib->max_iter = (long)n_rows_ib + (long)n_cols_ib + 50;

	fp_lsqr_ib = NULL;
	fp_lsqr_ib = fopen (lsqr_print_FN_ib, "w");
	if (NULL == fp_lsqr_ib)
	{
		printf ("Open LSQR Print File Error.\n");
	}
	in_lsqr_ib->lsqr_fp_out = fp_lsqr_ib;

	/*-------------------------------------------------------------------------------------------*/

	i_val_ib = alloc1int (prod_lsqr->n_cols);
	rhs_vec_ib = alloc1float (prod_lsqr->n_rows);

	/*-------------------------------------------------------------------------------------------*/

	for (i_row=0; i_row<prod_lsqr->n_rows; i_row++)
		rhs_vec_ib[i_row] = (float)rhs_vec[i_row];

	printf ("\tSolving linear equations ...\n");
	illumination_balance (
		i_val_ib,
		rhs_vec_ib,
		prod_lsqr,			prod_lsqr_ib,
		in_lsqr_ib,
		out_lsqr_ib,
		work_lsqr_ib,
		func_lsqr_ib
	);
	printf ("\tSolving linear equations end...\n");

	for (i_row=0; i_row<prod_lsqr->n_rows; i_row++)
		rhs_vec[i_row] =(double)rhs_vec_ib[i_row];

	/*-------------------------------------------------------------------------------------------*/
	//Backward Convert Matrix
	i=0;
	for (i_row=0; i_row<prod_lsqr->n_rows; i_row++)
	{
		for (i_col=0;i_col<prod_lsqr->m_row[i_row].n_val;i_col++)
		{
			tmpi=i_row;
			tmpj=prod_lsqr->m_row[tmpi].m_col[i_col].j_col;
			tmpp=prod_lsqr->m_row[tmpi].m_col[i_col].val;

			indx[i]=(int)tmpi+1;
			jndx[i]=(int)tmpj+1;
			path[i]=(double)tmpp;
			i=i+1;
		}
		free1matrix_col (prod_lsqr->m_row[i_row].m_col);
	}
	
	free1int(i_val_ib);
	free1float(rhs_vec_ib);
	free_lsqr_mem (in_lsqr_ib, out_lsqr_ib, work_lsqr_ib, func_lsqr_ib);

	free1matrix_row(prod_lsqr->m_row);
	free1matrix_row(prod_lsqr_ib->m_row);

	free1matrix_a(prod_lsqr);
	free1matrix_a(prod_lsqr_ib);

	/*-------------------------------------------------------------------------------------------*/

}

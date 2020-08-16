
int illumination_balance (
	int *i_val_ib,
	float *rhs_vec_ib,
	matrix_a *kernel,			matrix_a *mtx_ib,
	lsqr_input *in_lsqr_ib,
	lsqr_output *out_lsqr_ib,
	lsqr_work *work_lsqr_ib,
	lsqr_func *func_lsqr_ib
)
{
	int i_row, i_col, i_val, itmp1;
	float r_ib, ftmp1;

	//----------------------------------------------------------------------------------------------
	//Count the n_val of each row in the new matrix
	//And allocate memory for each row
	for (i_row=0; i_row<kernel->n_rows; i_row++)
	{
		for (i_val=0; i_val<kernel->m_row[i_row].n_val; i_val++)
		{
			mtx_ib->m_row[kernel->m_row[i_row].m_col[i_val].j_col].n_val ++;
		}
	}

	for (i_row=0; i_row<mtx_ib->n_rows-1; i_row++)
	{
		mtx_ib->m_row[i_row].n_val ++;
		mtx_ib->m_row[i_row].m_col = alloc1matrix_col (mtx_ib->m_row[i_row].n_val);
	}
	mtx_ib->m_row[mtx_ib->n_rows-1].n_val = mtx_ib->n_cols - 1;
	mtx_ib->m_row[mtx_ib->n_rows-1].m_col = alloc1matrix_col (mtx_ib->m_row[mtx_ib->n_rows-1].n_val);
	//----------------------------------------------------------------------------------------------
	
	//----------------------------------------------------------------------------------------------
	//Assignment for the new matrix

	//Transpose of the kernel matrix
	memset ((void*)&i_val_ib[0], 0, sizeof(int)*kernel->n_cols);
	for (i_row=0; i_row<kernel->n_rows; i_row++)
	{
		for (i_val=0; i_val<kernel->m_row[i_row].n_val; i_val++)
		{
			itmp1 = kernel->m_row[i_row].m_col[i_val].j_col;
			mtx_ib->m_row[itmp1].m_col[i_val_ib[itmp1]].j_col = i_row;
			mtx_ib->m_row[itmp1].m_col[i_val_ib[itmp1]].val = kernel->m_row[i_row].m_col[i_val].val;
			i_val_ib[itmp1] ++;
		}
	}

	//The last collumn
	for (i_row=0; i_row<mtx_ib->n_rows-1; i_row++)
	{
		mtx_ib->m_row[i_row].m_col[mtx_ib->m_row[i_row].n_val-1].j_col = mtx_ib->n_cols - 1;
		mtx_ib->m_row[i_row].m_col[mtx_ib->m_row[i_row].n_val-1].val = -1.0;
	}

	//The last row
	r_ib = 0;
	for (i_row=0; i_row<kernel->n_rows; i_row++)
		for (i_val=0; i_val<kernel->m_row[i_row].n_val; i_val++)
			r_ib += kernel->m_row[i_row].m_col[i_val].val;
	ftmp1 = r_ib / (float)(kernel->n_rows * kernel->n_cols);
	mtx_ib->m_row[mtx_ib->n_rows-1].n_val = mtx_ib->n_cols - 1;
	for (i_val=0; i_val<mtx_ib->m_row[mtx_ib->n_rows-1].n_val; i_val++)
	{
		mtx_ib->m_row[mtx_ib->n_rows-1].m_col[i_val].j_col = i_val;
		mtx_ib->m_row[mtx_ib->n_rows-1].m_col[i_val].val = ftmp1;
	}
	//----------------------------------------------------------------------------------------------
	
	for (i_row=0; i_row<mtx_ib->n_rows-1; i_row++)
		in_lsqr_ib->rhs_vec->elements[i_row] = 0;

	in_lsqr_ib->rhs_vec->elements[mtx_ib->n_rows-1] = r_ib;

	for (i_col=0; i_col<mtx_ib->n_cols; i_col++)
		in_lsqr_ib->sol_vec->elements[i_col] = 0;

	lsqr (in_lsqr_ib, out_lsqr_ib, work_lsqr_ib, func_lsqr_ib, mtx_ib);

	for (i_row=0; i_row<kernel->n_rows; i_row++)
	{
		ftmp1 = in_lsqr_ib->sol_vec->elements[i_row];
		for (i_val=0; i_val<kernel->m_row[i_row].n_val; i_val++)
			kernel->m_row[i_row].m_col[i_val].val *= ftmp1;
		rhs_vec_ib[i_row] *= ftmp1;
	}
	
	for (i_row=0; i_row<mtx_ib->n_rows; i_row++)
		free1matrix_col (mtx_ib->m_row[i_row].m_col);

	return 0;
}

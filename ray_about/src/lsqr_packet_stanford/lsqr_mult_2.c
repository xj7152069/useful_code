typedef struct
{
	long j_col;
	float val;
}matrix_col;

typedef struct
{
	long n_val;
	matrix_col *m_col;
}matrix_row;

typedef struct
{
	long n_rows;
	long n_cols;
	matrix_row *m_row;
}matrix_a;

matrix_col *alloc1matrix_col (int n1)
{
	return (matrix_col *)alloc1 (n1, sizeof(matrix_col));
}

void free1matrix_col (matrix_col *p)
{
	free1 (p);
}

matrix_row *alloc1matrix_row (int n1)
{
	return (matrix_row *)alloc1 (n1, sizeof(matrix_row));
}

void free1matrix_row (matrix_row *p)
{
	free1 (p);
}

matrix_a *alloc1matrix_a (int n1)
{
	return (matrix_a *)alloc1 (n1, sizeof(matrix_a));
}

void free1matrix_a (matrix_a *p)
{
	free1 (p);
}

void mtx_mult (long mode, dvec *x, dvec *y, void *prod)
{
	int i, j;
	matrix_a *data;

	data = (matrix_a *) prod;

/*  Calculate Y = Y + A * X  */
	/*
	if (0 == mode)
	{
		for (i=0; i<data->n_rows; i++)
		{
			ftmp1 = 0;
			for (j=0; j<data->m_row[i].n_val; j++)
				ftmp1 += data->m_row[i].m_col[j].val * x->elements[data->m_row[i].m_col[j].j_col];
			y->elements[i] += ftmp1;
		}
	}
	*/
	if (0 == mode)
	{
		for (i=0; i<data->n_rows; i++)
		{
			for (j=0; j<data->m_row[i].n_val; j++)
				y->elements[i] += data->m_row[i].m_col[j].val * x->elements[data->m_row[i].m_col[j].j_col];
		}
	}
/*  Calculate X = X + A^T * Y  */
	/*
	else if (1 == mode)
	{
		for (i=0; i<data->n_cols; i++)
		{
			ftmp1 = 0;
			for (j=0; j<data->n_rows; j++)
			{
				for (k=0; k<data->m_row[j].n_val; k++)
				{
					if (i == data->m_row[j].m_col[k].j_col)
						break;
				}
				if (k == data->m_row[j].n_val)
					continue;
				ftmp1 += data->m_row[j].m_col[k].val * y->elements[j];
			}
			x->elements[i] += ftmp1;
		}
	}
	*/
	else if (1 == mode)
	{
		for (i=0; i<data->n_rows; i++)
		{
			for (j=0; j<data->m_row[i].n_val; j++)
			{
				x->elements[data->m_row[i].m_col[j].j_col] += data->m_row[i].m_col[j].val * y->elements[i];
			}
		}
	}
	else
	{
		printf ("Parameter Error in Matrix Multiplization Function.\n");
	}

	return;
}

int reading_kernel (matrix_a *kernel, char *kernel_file_name)
{
	int i_row, i_val;
	FILE *fp;

	fp = NULL;
	fp = fopen (kernel_file_name, "rb");
	if (NULL == fp)
	{
		printf ("Open file \"%s\" error.\n", kernel_file_name);
		return 1;
	}

	fread (&kernel->n_rows, sizeof(long), 1, fp);
	fread (&kernel->n_cols, sizeof(long), 1, fp);

	kernel->m_row = NULL;
	kernel->m_row = alloc1matrix_row (kernel->n_rows);
	if (NULL == kernel->m_row)
	{
		printf ("Memory allocation for \"m_row\" error.\n");
		return 1;
	}

	for (i_row=0; i_row<kernel->n_rows; i_row++)
	{
		fread (&kernel->m_row[i_row].n_val, sizeof(long), 1, fp);
		
		kernel->m_row[i_row].m_col = NULL;
		kernel->m_row[i_row].m_col = alloc1matrix_col (kernel->m_row[i_row].n_val);
		if (NULL == kernel->m_row[i_row].m_col)
		{
			printf ("Memory allocation for \"m_col\" error.\n");
			return 1;
		}

		for (i_val=0; i_val<kernel->m_row[i_row].n_val; i_val++)
		{
			fread (&kernel->m_row[i_row].m_col[i_val].j_col, sizeof(long), 1, fp);
			fread (&kernel->m_row[i_row].m_col[i_val].val, sizeof(float), 1, fp);
		}
	}

	fclose (fp);

	return 0;
}

int writing_kernel (matrix_a *kernel, char *kernel_file_name)
{
	int i_row, i_val;
	FILE *fp;

	fp = NULL;
	fp = fopen (kernel_file_name, "wb");
	if (NULL == fp)
	{
		printf ("Open file \"%s\" error.\n", kernel_file_name);
		return 1;
	}

	fwrite (&kernel->n_rows, sizeof(long), 1, fp);
	fwrite (&kernel->n_cols, sizeof(long), 1, fp);

	for (i_row=0; i_row<kernel->n_rows; i_row++)
	{
		fwrite (&kernel->m_row[i_row].n_val, sizeof(long), 1, fp);

		for (i_val=0; i_val<kernel->m_row[i_row].n_val; i_val++)
		{
			fwrite (&kernel->m_row[i_row].m_col[i_val].j_col, sizeof(long), 1, fp);
			fwrite (&kernel->m_row[i_row].m_col[i_val].val, sizeof(float), 1, fp);
		}
	}

	fclose (fp);

	return 0;
}

int free_kernel (matrix_a *kernel)
{
	int i_row;

	for (i_row=0; i_row<kernel->n_rows; i_row++)
	{
		free1matrix_col (kernel->m_row[i_row].m_col);
	}
	free1matrix_row (kernel->m_row);

	return 0;
}


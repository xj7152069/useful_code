typedef struct
{
	long n_rows;
	long n_cols;
	double **a;
}matrix_a;

void mtx_mult (long mode, dvec *x, dvec *y, void *prod)
{

	int i, j;
	double ftmp;
	matrix_a *data;

	data = (matrix_a *) prod;

/*  Calculate Y = Y + A * X  */
	if (0 == mode)
	{
		for (i=0; i<data->n_rows; i++)
		{
			ftmp = 0;
			for (j=0; j<data->n_cols; j++)
			{
				ftmp = ftmp + data->a[i][j] * x->elements[j];
			}
			y->elements[i] = y->elements[i] + ftmp;
		}
	}
/*  Calculate X = X + A^T * Y  */
	else if (1 == mode)
	{
		for (i=0; i<data->n_cols; i++)
		{
			ftmp = 0;
			for (j=0; j<data->n_rows; j++)
			{
				ftmp = ftmp + data->a[j][i] * y->elements[j];
			}
			x->elements[i] = x->elements[i] + ftmp;
		}
	}
	else
	{
		printf ("Parameter Error in Matrix Multiplization Function.\n");
	}

	return;
}

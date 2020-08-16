
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

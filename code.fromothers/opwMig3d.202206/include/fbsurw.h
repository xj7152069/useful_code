#ifndef FB_SURW_H
#define FB_SURW_H

void read_2d_float_rb_su ( float **trace , fbsegy_std *segy2d, int n1 , int n2 , char *fn );
void write_2d_float_wb_suhdr ( float **trace , fbsegy_std *segy2d, int n1 , int n2 , char *fn );
void write_3d_float_wb_suhdr ( float ***trace , fbsegy_std **segy2d, int n1 , int n2 , int n3, char *fn );

#endif

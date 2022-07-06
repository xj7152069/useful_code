  int spl_z(float *x, float *y, int n, float *dy, float t,
  	float *z, float *s)
  { int i,j;
    float h0,h1;
    //s=malloc(n*sizeof(float));

    for (j=0;j<=n-2;j++) s[j]=x[j+1]-x[j];
    if (t>=x[n-1]) i=n-2;
    else
      { i=0;
        while (t>x[i+1]) i=i+1;
      }
    h1=(x[i+1]-t)/s[i];
    h0=h1*h1;
    *z=(3.0*h0-2.0*h0*h1)*y[i];
    *z=*z+s[i]*(h0-h0*h1)*dy[i];
    h1=(t-x[i])/s[i];
    h0=h1*h1;
    *z=*z+(3.0*h0-2.0*h0*h1)*y[i+1];
    *z=*z-s[i]*(h0-h0*h1)*dy[i+1];

    //free(s);
    return 0;
  }

  int spl_dz(float *x, float *y, int n, float *dy, float t,
  	float *dz, float *s)
  { int i,j;
    float h0,h1;
    //s=malloc(n*sizeof(float));

    for (j=0;j<=n-2;j++) s[j]=x[j+1]-x[j];
    if (t>=x[n-1]) i=n-2;
    else
      { i=0;
        while (t>x[i+1]) i=i+1;
      }
    h1=(x[i+1]-t)/s[i];
    h0=h1*h1;
    *dz=6.0*(h0-h1)*y[i]/s[i];
    *dz=*dz+(3.0*h0-2.0*h1)*dy[i];
    h1=(t-x[i])/s[i];
    h0=h1*h1;
    *dz=*dz-6.0*(h0-h1)*y[i+1]/s[i];
    *dz=*dz+(3.0*h0-2.0*h1)*dy[i+1];

    //free(s);
    return 0;
  }

  int spl_ddz(float *x, float *y, int n, float *dy, float t,
  	float *ddz, float *s)
  { int i,j;
    float h0,h1;
    //s=malloc(n*sizeof(float));

    for (j=0;j<=n-2;j++) s[j]=x[j+1]-x[j];
    if (t>=x[n-1]) i=n-2;
    else
      { i=0;
        while (t>x[i+1]) i=i+1;
      }
    h1=(x[i+1]-t)/s[i];
    h0=h1*h1;
    *ddz=(6.0-12.0*h1)*y[i]/(s[i]*s[i]);
    *ddz=*ddz+(2.0-6.0*h1)*dy[i]/s[i];
    h1=(t-x[i])/s[i];
    h0=h1*h1;
    *ddz=*ddz+(6.0-12.0*h1)*y[i+1]/(s[i]*s[i]);
    *ddz=*ddz-(2.0-6.0*h1)*dy[i+1]/s[i];

    //free(s);
    return 0;
  }

  int spl_z_dz(float *x, float *y, int n, float *dy, float t,
  	float *z, float *dz, float *s)
  { int i,j;
    float h0,h1;
    //s=malloc(n*sizeof(float));

    for (j=0;j<=n-2;j++) s[j]=x[j+1]-x[j];
    if (t>=x[n-1]) i=n-2;
    else
      { i=0;
        while (t>x[i+1]) i=i+1;
      }
    h1=(x[i+1]-t)/s[i];
    h0=h1*h1;
    *z=(3.0*h0-2.0*h0*h1)*y[i];
    *z=*z+s[i]*(h0-h0*h1)*dy[i];
    *dz=6.0*(h0-h1)*y[i]/s[i];
    *dz=*dz+(3.0*h0-2.0*h1)*dy[i];
    h1=(t-x[i])/s[i];
    h0=h1*h1;
    *z=*z+(3.0*h0-2.0*h0*h1)*y[i+1];
    *z=*z-s[i]*(h0-h0*h1)*dy[i+1];
    *dz=*dz-6.0*(h0-h1)*y[i+1]/s[i];
    *dz=*dz+(3.0*h0-2.0*h1)*dy[i+1];

    //free(s);
    return 0;
  }

  int spl_z_dz_ddz(float *x, float *y, int n, float *dy, float t,
  	float *z, float *dz, float *ddz, float *s)
  { int i,j;
    float h0,h1;
    //s=malloc(n*sizeof(float));

    for (j=0;j<=n-2;j++) s[j]=x[j+1]-x[j];
    if (t>=x[n-1]) i=n-2;
    else
      { i=0;
        while (t>x[i+1]) i=i+1;
      }
    h1=(x[i+1]-t)/s[i];
    h0=h1*h1;
    *z=(3.0*h0-2.0*h0*h1)*y[i];
    *z=*z+s[i]*(h0-h0*h1)*dy[i];
    *dz=6.0*(h0-h1)*y[i]/s[i];
    *dz=*dz+(3.0*h0-2.0*h1)*dy[i];
	*ddz=(6.0-12.0*h1)*y[i]/(s[i]*s[i]);
	*ddz=ddz[j]+(2.0-6.0*h1)*dy[i]/s[i];
    h1=(t-x[i])/s[i];
    h0=h1*h1;
    *z=*z+(3.0*h0-2.0*h0*h1)*y[i+1];
    *z=*z-s[i]*(h0-h0*h1)*dy[i+1];
    *dz=*dz-6.0*(h0-h1)*y[i+1]/s[i];
    *dz=*dz+(3.0*h0-2.0*h1)*dy[i+1];
	*ddz=ddz[j]+(6.0-12.0*h1)*y[i+1]/(s[i]*s[i]);
	*ddz=ddz[j]-(2.0-6.0*h1)*dy[i+1]/s[i];

    //free(s);
    return 0;
  }

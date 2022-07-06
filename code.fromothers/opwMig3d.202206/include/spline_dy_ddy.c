  int spl_dy(float *x, float *y, int n, float *dy, float *s)
  { int j;
    float h0,h1,alpha,beta;
    s[0]=dy[0]; dy[0]=0.0;
    h0=x[1]-x[0];
    for (j=1;j<=n-2;j++)
      { h1=x[j+1]-x[j];
        alpha=h0/(h0+h1);
        beta=(1.0-alpha)*(y[j]-y[j-1])/h0;
        beta=3.0*(beta+alpha*(y[j+1]-y[j])/h1);
        dy[j]=-alpha/(2.0+(1.0-alpha)*dy[j-1]);
        s[j]=(beta-(1.0-alpha)*s[j-1]);
        s[j]=s[j]/(2.0+(1.0-alpha)*dy[j-1]);
        h0=h1;
      }
    for (j=n-2;j>=0;j--)
      dy[j]=dy[j]*dy[j+1]+s[j];
	/*
	for (j=0;j<=n-2;j++) s[j]=x[j+1]-x[j];
	for (j=0;j<=n-2;j++)
	{
		h1=s[j]*s[j];
	}
	h1=s[n-2]*s[n-2];
	*/

    return 0;
  }

  int spl_dy_ddy(float *x, float *y, int n, float *dy, float *ddy, float *s)
  { int j;
    float h0,h1,alpha,beta;
    //s=malloc(n*sizeof(float));
    s[0]=dy[0]; dy[0]=0.0;
    h0=x[1]-x[0];
    for (j=1;j<=n-2;j++)
      { h1=x[j+1]-x[j];
        alpha=h0/(h0+h1);
        beta=(1.0-alpha)*(y[j]-y[j-1])/h0;
        beta=3.0*(beta+alpha*(y[j+1]-y[j])/h1);
        dy[j]=-alpha/(2.0+(1.0-alpha)*dy[j-1]);
        s[j]=(beta-(1.0-alpha)*s[j-1]);
        s[j]=s[j]/(2.0+(1.0-alpha)*dy[j-1]);
        h0=h1;
      }
    for (j=n-2;j>=0;j--)
      dy[j]=dy[j]*dy[j+1]+s[j];
	for (j=0;j<=n-2;j++) s[j]=x[j+1]-x[j];
	for (j=0;j<=n-2;j++)
	{
		h1=s[j]*s[j];
		ddy[j]=6.0*(y[j+1]-y[j])/h1-2.0*(2.0*dy[j]+dy[j+1])/s[j];
	}
	h1=s[n-2]*s[n-2];
	ddy[n-1]=6.*(y[n-2]-y[n-1])/h1+2.*(2.*dy[n-1]+dy[n-2])/s[n-2];

    //free(s);
    return 0;
  }


/* This is code by Zhipeng Cao and Will Grissom who have given permission 
% for inclusion within this package. This code can be found also in their
% repo https://bitbucket.org/wgrissom/acptx/ Please cite appropriately.*/
        
#define PI 3.14159265

__global__
void kGradHess( double *vx,
		double *Sx,
		double *vy,
		double *Sy,
		double *vz,
		double *Sz,
		double *Ar,
		double *Ai,
		double *dr,
		double *di,
		double *x,
		double *y,
		double *z,
		int *m2v,
		const int ns,
		const int nd,
		const int nr)
{

  int j = blockIdx.x*blockDim.x + threadIdx.x;
  double ajdr, ajdi, sjl, tmp;
  double ajai_i, ajai_r, sji, tmp2;  
  int i,k;
    
  /* separable term */
  while(j < ns){
    for(i=0;i<nr;i++){
	
      ajdr = dr[j]*Ar[i*ns+j] + di[j]*Ai[i*ns+j];
      ajdi = dr[j]*Ai[i*ns+j] - di[j]*Ar[i*ns+j];
      sjl = atan2(ajdi,ajdr);
      
      switch(nd){
      case 1:
	tmp = 2*PI*ajdi*x[j];
	vx[i*ns+j] = -tmp;
	if(sjl != 0)
	  Sx[m2v[i*nr+i]*ns+j] = 2*PI*x[j]*tmp/sjl;
	break;
      case 2:
	tmp = 2*PI*ajdi*x[j];
	vx[i*ns+j] = -tmp;
	if(sjl != 0)
	  Sx[m2v[i*nr+i]*ns+j] = 2*PI*x[j]*tmp/sjl;
	
	tmp = 2*PI*ajdi*y[j];
	vy[i*ns+j] = -tmp;
	if(sjl != 0)
	  Sy[m2v[i*nr+i]*ns+j] = 2*PI*y[j]*tmp/sjl;
	break;
      case 3:
	tmp = 2*PI*ajdi*x[j];
	vx[i*ns+j] = -tmp;
	if(sjl != 0)
	  Sx[m2v[i*nr+i]*ns+j] = 2*PI*x[j]*tmp/sjl;
	
	tmp = 2*PI*ajdi*y[j];
	vy[i*ns+j] = -tmp;
	if(sjl != 0)
	  Sy[m2v[i*nr+i]*ns+j] = 2*PI*y[j]*tmp/sjl;
	
	tmp = 2*PI*ajdi*z[j];
	vz[i*ns+j] = -tmp;
	if(sjl != 0)
	  Sz[m2v[i*nr+i]*ns+j] = 2*PI*z[j]*tmp/sjl;
	break;
      }      
    }
    
    /* non-separable term */
    k = j;
    for(i=0;i<nr;i++){
      for(j=0;j<i;j++){  
	ajai_r = -Ar[i*ns+k]*Ar[j*ns+k] - Ai[i*ns+k]*Ai[j*ns+k];
	ajai_i = -Ar[i*ns+k]*Ai[j*ns+k] + Ai[i*ns+k]*Ar[j*ns+k];
	
	/* calculate phase angle */
	sji = atan2(ajai_i,ajai_r);
	
	switch(nd){
	case 1:
	  tmp = 2*PI*x[k]*ajai_i;
	  vx[i*ns+k] += tmp;
	  vx[j*ns+k] -= tmp;
	  
	  if(sji != 0)
	    {
	      tmp2 = 2*PI*x[k]*tmp/sji;
	      Sx[m2v[i*nr+i]*ns+k] += tmp2;
	      Sx[m2v[j*nr+i]*ns+k] = -tmp2;
	      Sx[m2v[j*nr+j]*ns+k] += tmp2;
	      /*Sx[m2v[i*nr+j]*ns+k] = -tmp2;*/
	    }
	  break;
	case 2:
	  tmp = 2*PI*x[k]*ajai_i;
	  vx[i*ns+k] += tmp;
	  vx[j*ns+k] -= tmp;
	  
	  if(sji != 0)
	    {
	      tmp2 = 2*PI*x[k]*tmp/sji;
	      Sx[m2v[i*nr+i]*ns+k] += tmp2;
	      Sx[m2v[j*nr+i]*ns+k] = -tmp2;
	      Sx[m2v[j*nr+j]*ns+k] += tmp2;
	      /*Sx[m2v[i*nr+j]*ns+k] = -tmp2;*/
	    }
	  
	  tmp = 2*PI*y[k]*ajai_i;
	  vy[i*ns+k] += tmp;
	  vy[j*ns+k] -= tmp;
	  
	  if(sji != 0)
	    {
	      tmp2 = 2*PI*y[k]*tmp/sji;
	      Sy[m2v[i*nr+i]*ns+k] += tmp2;
	      Sy[m2v[j*nr+i]*ns+k] = -tmp2;
	      Sy[m2v[j*nr+j]*ns+k] += tmp2;
	      /*Sy[m2v[i*nr+j]*ns+k] = -tmp2;*/
	    }
	  break;
	case 3:
	  tmp = 2*PI*x[k]*ajai_i;
	  vx[i*ns+k] += tmp;
	  vx[j*ns+k] -= tmp;
	  
	  if(sji != 0){
	    tmp2 = 2*PI*x[k]*tmp/sji;
	    Sx[m2v[i*nr+i]*ns+k] += tmp2;
	    Sx[m2v[j*nr+i]*ns+k] = -tmp2;
	    Sx[m2v[j*nr+j]*ns+k] += tmp2;
	    /*Sx[m2v[i*nr+j]*ns+k] -= tmp2;*/
	  }
	  
	  tmp = 2*PI*y[k]*ajai_i;
	  vy[i*ns+k] += tmp;
	  vy[j*ns+k] -= tmp;
	  
	  if(sji != 0){
	    tmp2 = 2*PI*y[k]*tmp/sji;
	    Sy[m2v[i*nr+i]*ns+k] += tmp2;
	    Sy[m2v[j*nr+i]*ns+k] = -tmp2;
	    Sy[m2v[j*nr+j]*ns+k] += tmp2;
	    /*Sy[m2v[i*nr+j]*ns+k] -= tmp2;*/
	  }
	  
	  tmp = 2*PI*z[k]*ajai_i;
	  vz[i*ns+k] += tmp;
	  vz[j*ns+k] -= tmp;
	  
	  if(sji != 0){
	    tmp2 = 2*PI*z[k]*tmp/sji;
	    Sz[m2v[i*nr+i]*ns+k] += tmp2;
	    Sz[m2v[j*nr+i]*ns+k] = -tmp2;
	    Sz[m2v[j*nr+j]*ns+k] += tmp2;
	    /*Sz[m2v[i*nr+j]*ns+k] -= tmp2;*/
	  }  
	  break;
	} /* nd switch */
	
      } /* nr inner loop */
    } /* nr outer loop */
    
    j = k;
    j += blockDim.x * gridDim.x;
    
  } /* while j < ns */

}



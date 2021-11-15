/* This is code by Zhipeng Cao and Will Grissom who have given permission 
% for inclusion within this package. This code can be found also in their
% repo https://bitbucket.org/wgrissom/acptx/ Please cite appropriately.*/

/***********************************************************************
* Program Name: calcvS.c
* Author:       Will Grissom
***********************************************************************/
/* #define NOTHREADS */

#ifndef NOTHREADS
#include <pthread.h>
#endif
#include <math.h>
#include "mex.h"

/* Input Arguments */

#define	A	prhs[0] /* excitation vectors */ 
#define D       prhs[1] /* desired patterns */ 
#define XX      prhs[2] /* encoding vectors */
#define NT      prhs[3] /* number of threads */

/* Output Arguments */

#define	VXOUT	plhs[0]
#define SXOUT   plhs[1]
#define	VYOUT	plhs[2]
#define SYOUT   plhs[3]
#define	VZOUT	plhs[4]
#define SZOUT   plhs[5]

#define PI 3.14159265

typedef struct {
  int m;
  int nrungs;
  double *Ar,*Ai;
  double *dr,*di;
  double *x,*y,*z;
  double *vx,*vy,*vz; /* each thread has its own copy of v, S; we sum these over space after all threads are done */
  double *Sx,*Sy,*Sz;
  int si,ei; /* start/end spatial indices */
} params;

void *docalcvS(void *arg)
{
    
  params *p = (params *)arg;

  /* separable term */
  int i = 0; int j = 0; int k = 0;
  int m = p->m;
  int nrungs = p->nrungs;
  double ajdr, ajdi, sjl, tmp;
  double ajai_i, ajai_r, sji, tmp2;
  
  for(i=0;i<nrungs;i++)
  {
      /* loop over space */
      for(j=p->si;j<p->ei;j++)
	{
	  ajdr = p->dr[j]*p->Ar[i*m+j] + p->di[j]*p->Ai[i*m+j];
	  ajdi = p->dr[j]*p->Ai[i*m+j] - p->di[j]*p->Ar[i*m+j];
	  sjl = atan2(ajdi,ajdr);
	  
	  tmp = 2*PI*ajdi*p->x[j];
	  p->vx[i] -= tmp;
	  if(sjl != 0)
	    p->Sx[i*nrungs+i] += 2*PI*p->x[j]*tmp/sjl;

	  tmp = 2*PI*ajdi*p->y[j];
	  p->vy[i] -= tmp;
	  if(sjl != 0)
	    p->Sy[i*nrungs+i] += 2*PI*p->y[j]*tmp/sjl;

	  tmp = 2*PI*ajdi*p->z[j];
	  p->vz[i] -= tmp;
	  if(sjl != 0)
	    p->Sz[i*nrungs+i] += 2*PI*p->z[j]*tmp/sjl;

	}
  }

  /* non-separable term */
  /* int k = 0; /* Zhipeng Commented out */
  /* double ajai_i, ajai_r, aiaj_i, aiaj_r, sji, sij, tmp2; /* Zhipeng Commented out */

  for(i=0;i<nrungs;i++)
    {
      for(j=0;j<i;j++)
	{  
	  /* loop over space */
	  for(k=p->si;k<p->ei;k++)
	    {
	      ajai_r = -p->Ar[i*m+k]*p->Ar[j*m+k] - p->Ai[i*m+k]*p->Ai[j*m+k];
	      ajai_i = -p->Ar[i*m+k]*p->Ai[j*m+k] + p->Ai[i*m+k]*p->Ar[j*m+k];
	      
	      /* calculate phase angle */
	      sji = atan2(ajai_i,ajai_r);
	      
	      tmp = 2*PI*p->x[k]*ajai_i;
	      p->vx[i] += tmp;
	      p->vx[j] -= tmp;

	      if(sji != 0)
		{
		  tmp2 = 2*PI*p->x[k]*tmp/sji;
		  p->Sx[i*nrungs+i] += tmp2;
		  p->Sx[j*nrungs+i] -= tmp2;
		  p->Sx[j*nrungs+j] += tmp2;
		  p->Sx[i*nrungs+j] -= tmp2;
		}

	      tmp = 2*PI*p->y[k]*ajai_i;
	      p->vy[i] += tmp;
	      p->vy[j] -= tmp;

	      if(sji != 0)
		{
		  tmp2 = 2*PI*p->y[k]*tmp/sji;
		  p->Sy[i*nrungs+i] += tmp2;
		  p->Sy[j*nrungs+i] -= tmp2;
		  p->Sy[j*nrungs+j] += tmp2;
		  p->Sy[i*nrungs+j] -= tmp2;
		}

	      tmp = 2*PI*p->z[k]*ajai_i;
	      p->vz[i] += tmp;
	      p->vz[j] -= tmp;

	      if(sji != 0)
		{
		  tmp2 = 2*PI*p->z[k]*tmp/sji;
		  p->Sz[i*nrungs+i] += tmp2;
		  p->Sz[j*nrungs+i] -= tmp2;
		  p->Sz[j*nrungs+j] += tmp2;
		  p->Sz[i*nrungs+j] -= tmp2;
		}


	    }
	}
    }

#ifndef NOTHREADS
  pthread_exit(NULL);
#endif
}

void mexFunction( int nlhs, mxArray *plhs[], 
		  int nrhs, const mxArray*prhs[] )
     
{ 
    double *Ar; 
    double *Ai;
    double *dr;
    double *di;
    double *vx,*vy,*vz;
    double *Sx,*Sy,*Sz;
    double *xx;
    int *nst; /* number of spatial locs in each thread */
    int m, nrungs;
    int nthreads;
    int i,j,k;

    int si;
    
#ifndef NOTHREADS
    pthread_t *threads;
#endif
    int rc;
    params *p;

    /* Check for proper number of arguments */
    
    if (nrhs != 4) { 
	mexErrMsgTxt("Four input arguments required."); 
    } else if (nlhs > 6) {
	mexErrMsgTxt("Too many output arguments."); 
    } 

    /* Assign pointers to the various parameters */ 
    Ar = mxGetPr(A);
    Ai = mxGetPi(A);
    xx = mxGetPr(XX);
    dr = mxGetPr(D);
    di = mxGetPi(D);

    m = mxGetM(A);
    nrungs = mxGetN(A);

    nthreads = mxGetScalar(NT);
    nst = mxCalloc(nthreads,sizeof(int));
    for(i = 0;i < nthreads;i++) nst[i] = (int)(m/nthreads);
    if(m - ((int)(m/nthreads))*nthreads > 0) /* need to add remainder to last thread */
      nst[nthreads-1] += m - ((int)(m/nthreads))*nthreads;    

    /* Create matrices for the return arguments */ 
    VXOUT = mxCreateDoubleMatrix(nrungs, 1, mxREAL); 
    SXOUT = mxCreateDoubleMatrix(nrungs, nrungs, mxREAL);
    vx = mxGetPr(VXOUT);
    Sx = mxGetPr(SXOUT);
    VYOUT = mxCreateDoubleMatrix(nrungs, 1, mxREAL); 
    SYOUT = mxCreateDoubleMatrix(nrungs, nrungs, mxREAL);
    vy = mxGetPr(VYOUT);
    Sy = mxGetPr(SYOUT);
    VZOUT = mxCreateDoubleMatrix(nrungs, 1, mxREAL); 
    SZOUT = mxCreateDoubleMatrix(nrungs, nrungs, mxREAL);
    vz = mxGetPr(VZOUT);
    Sz = mxGetPr(SZOUT);

    /* Build param structs for each thread */
    p = mxCalloc(nthreads,sizeof(params));
    si = 0;
    for(i = 0;i < nthreads;i++){
      
      /* scalars */
      p[i].m = m;
      p[i].nrungs = nrungs;
      
      /* spatial locs to work on */
      p[i].si = si;
      p[i].ei = si+nst[i];
      si += nst[i];

      /* pointers */
      p[i].Ar = Ar;
      p[i].Ai = Ai;
      p[i].dr = dr;
      p[i].di = di;
      p[i].x = xx;
      p[i].y = &xx[m];
      p[i].z = &xx[2*m];
      
      /* solutions */
      p[i].vx = mxCalloc(nrungs,sizeof(double));
      p[i].Sx = mxCalloc(nrungs*nrungs,sizeof(double));
      p[i].vy = mxCalloc(nrungs,sizeof(double));
      p[i].Sy = mxCalloc(nrungs*nrungs,sizeof(double));
      p[i].vz = mxCalloc(nrungs,sizeof(double));
      p[i].Sz = mxCalloc(nrungs*nrungs,sizeof(double));

    }
    
#ifndef NOTHREADS
    threads = mxCalloc(nthreads,sizeof(pthread_t));
    for(i = 0;i < nthreads;i++){
      rc = pthread_create(&threads[i], NULL, docalcvS, (void *)&p[i]);
      if (rc){
	mexErrMsgTxt("problem with return code from pthread_create()");
      }
    }

    /* wait for all threads to finish */
    for(i=0;i<nthreads;i++){
      pthread_join(threads[i],NULL);
    }
#else
    /* process them serially */
    for(i = 0;i < nthreads;i++) docalcvS((void *)&p[i]);
#endif

    /* sum results */
    for(i=0;i < nrungs;i++){
      for(j=0;j < nthreads;j++){
	vx[i] += p[j].vx[i];
	vy[i] += p[j].vy[i];
	vz[i] += p[j].vz[i];
      }
    }
    for(i=0;i < nrungs;i++){
      for(k=0;k < nrungs;k++){
	for(j=0;j < nthreads;j++){
	  Sx[i*nrungs+k] += p[j].Sx[i*nrungs+k];
	  Sy[i*nrungs+k] += p[j].Sy[i*nrungs+k];
	  Sz[i*nrungs+k] += p[j].Sz[i*nrungs+k];
	}
      }
    }    

    /* free memory */
    for(i=0;i < nthreads;i++){
      mxFree(p[i].vx);
      mxFree(p[i].Sx);
      mxFree(p[i].vy);
      mxFree(p[i].Sy);
      mxFree(p[i].vz);
      mxFree(p[i].Sz);
    }
    mxFree(p);
#ifndef NOTHREADS
    mxFree(threads);
#endif
    mxFree(nst);

    return;
    
}

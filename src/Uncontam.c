#include <R_ext/Applic.h>
#include <stdio.h>
#include <stdlib.h>
#include <R.h>
#include <Rmath.h>
#include <string.h>
#ifdef __APPLE__
#include <Accelerate/Accelerate.h>
#else
#include <R_ext/Lapack.h>
#include <R_ext/BLAS.h>
#endif
#include "functions.h"
#include "functionsC.h"
#define COMMENTS 0 

void densityU(int N, int p, int G, double *z, double **Sigma, 
               double **invSigma, double *mu,double *x,  
               double *logdet, double* PXgood){
  int g, n;
  double *delta = malloc(sizeof(double)*N*G);
  
  // compute delta
  for(g=0; g<G; g++ )
    mahalanobis(g, N, p, x, z, G, mu, invSigma[g], delta);
  
  // PXgood
  for(g=0; g<G; g++) {  
    for(n =0; n<N; n++){
      PXgood[n + N*g] = exp(- ((double)1.0/2)*delta[n + N*g]-((double)1.0/2.0)*logdet[g] - ((double)p/(double)2.0)*log(2*M_PI));
    }
  }
  free(delta);
}
double llikU(int N, int p, int G, double *z, double *prior, double **Sigma, 
                  double **invSigma, double *mu,double *x,  
                  double *logdet, int *lab, double* PXgood){
    int g, n; 
  // Compute llikelihood
  double rowSums=0.0;
  double llvalue=0.0;
  for(n =0; n<N; n++){
    for(g=0; g<G; g++) {
      if (lab[n] ==0){
        rowSums += prior[g] *PXgood[n + N*g];
      } else
        llvalue += 
          z[n+N*g] *(log(prior[g]) + log(PXgood[n + N*g]));
    }
    if (lab[n] ==0) llvalue += log(rowSums);
    rowSums = 0.0;
  }
  return(llvalue);
}
void estepU(int N, int p, int G, double *z, double *prior,double* PXgood, int *lab){
int g,n;
  double* znum  = malloc(sizeof(double)*G);
  double* zrowsum  = malloc(sizeof(double)*N);
  double  znumrowsum = 0;

  for(n =0; n<N; n++){ 
    znumrowsum = 0;
    zrowsum[n] = 0;
    for(g=0; g<G; g++) {
      znum[g] = prior[g]*(PXgood[n+g*N]);
      znumrowsum += znum[g];
    }
    for(g=0; g<G; g++){ 
      z[n+g*N] = znum[g]/znumrowsum;
      zrowsum[n] += z[n+g*N];
    }
  }  
  
  // z adjustment
  double colsum = 0;
  g=0;
  while (g<G) {
    colsum = 0;
    for(n =0; n<N; n++) colsum += z[n+g*N];
    g++;
    if (colsum <= p)
      for(n =0; n<N; n++) for(g=0; g<G; g++)
        z[n+g*N]=(z[n+g*N]+0.0000001)/(zrowsum[n]+G*0.0000001);
  } 
  
  for(n =0; n<N; n++) 
    if(lab[n] != 0){
      for(g=0; g<G; g++) z[n+g*N]=0;
      z[n+(lab[n]-1)*N]=1;
    }
free(znum);
free(zrowsum);
}

void loopU (int *NN, int *pp, int *GG, double *z, 
            double *sigmar, double *invsigmar, double *mu, double *mmtol, 
            int *mmmax, double *x, int *lab, char **covtype, 
            int *maxiter, double* threshold,double* prior,int* iteration,
            double *lllvalue,double *obslll, int* group){  
  int g, i, cc=0;
  int N = *NN;
  int p = *pp;
  int G = *GG;
  int exit = 0;
  double *D = (double*)malloc(sizeof(double)*p*p);
  double **Sigma      = malloc(sizeof(double*)*G);
  double **invSigma   = malloc(sizeof(double*)*G);
  double **sampcov    = malloc(sizeof(double*)*G);
  double *logdet      = malloc(sizeof(double)*G);
  double* PXgood = (double*)malloc(sizeof(double)*N*G);
  double loglik[] ={0,0,0};
  
  for(g=0; g < G; g++) {
    Sigma[g]        = malloc(sizeof(double)*p*p);
    invSigma[g]     = malloc(sizeof(double)*p*p);
    sampcov[g]      = malloc(sizeof(double)*p*p);
  }
  for(i=0; i<p*p; i++)
    D[i] = 0.0;
  for(i=0; i<p; i++)
    D[i*p + i] = (double)1.0;
  
  while(exit==0){
    cc++;
     
    get_pi( N, G, z, prior); // gets prior
   
    // gets Sigma, invSigma, logdet
    mstep(x, N, p, G, z, mu, sampcov, Sigma, invSigma, logdet,  *mmtol, *mmmax, D, covtype);
    
    // Returns loglik
    densityU( N, p, G, z, Sigma, invSigma, mu,x,logdet, PXgood);
    *lllvalue = llikU( N, p, G, z, prior, Sigma, invSigma, mu,x,logdet, lab, PXgood); 
    
    //Aitken's Stopping Criterion
    exit = stopcrit(G, *maxiter, cc, loglik, *lllvalue, *threshold);
    //E-Step 
    estepU(N, p, G, z, prior, PXgood, lab);
    Rprintf("*");
    //Rprintf("\n------------ iteration, %d \n",cc);
  }
  
  // Compute observed llikelihood
  double rowSums=0.0;
  int n;
  for(n =0; n<N; n++){
    for(g=0; g<G; g++) rowSums += prior[g] *PXgood[n + N*g];
    *obslll += log(rowSums);
    rowSums = 0.0;
  }
  
  get_group(G, N, z, group);
    // Prepare values to return
    for(g=0; g<G; g++){
    for(i=0; i<p*p; i++){
    sigmar[g*p*p +i] = Sigma[g][i];
    invsigmar[g*p*p +i] = invSigma[g][i];
    }
    }  
    *iteration = cc;
    
    
    for(g=0; g<G; g++) {
    free(sampcov[g]);
    free(Sigma[g]);
    free(invSigma[g]);
    }
    free(D);
    free(Sigma);
    free(invSigma); 
    free(sampcov);
    free(logdet);  
    free(PXgood);
}
void mstepU (int *NN, int *pp, int *GG, double *z, 
            double *sigmar, double *invsigmar, double *mu, double *mmtol, 
            int *mmmax, double *x, char **covtype, double *PXgood){
  int g, i; 
  int N = *NN;
  int p = *pp;
  int G = *GG;
  double *D = (double*)malloc(sizeof(double)*p*p);
  double **Sigma      = malloc(sizeof(double*)*G);
  double **invSigma   = malloc(sizeof(double*)*G);
  double **sampcov    = malloc(sizeof(double*)*G);
  double *logdet      = malloc(sizeof(double)*G);
  for(g=0; g < G; g++) {
    Sigma[g]        = malloc(sizeof(double)*p*p);
    invSigma[g]     = malloc(sizeof(double)*p*p);
    sampcov[g]      = malloc(sizeof(double)*p*p);
  }
  for(i=0; i<p*p; i++)
    D[i] = 0.0;
  for(i=0; i<p; i++)
    D[i*p + i] = (double)1.0;
    // gets Sigma, invSigma, logdet
    mstep(x, N, p, G, z, mu, sampcov, Sigma, invSigma, logdet,  *mmtol, *mmmax, D, covtype);
    densityU( N, p, G, z, Sigma, invSigma, mu,x,logdet, PXgood);
  // Prepare values to return
  for(g=0; g<G; g++){
    for(i=0; i<p*p; i++){
      sigmar[g*p*p +i] = Sigma[g][i];
      invsigmar[g*p*p +i] = invSigma[g][i];
    }
  }  
  for(g=0; g<G; g++) {
    free(sampcov[g]);
    free(Sigma[g]);
    free(invSigma[g]);
  }
  free(D);
  free(Sigma);
  free(invSigma); 
  free(sampcov);
  free(logdet);  
}

void dN(int *NN, int *pp, int *GG,double *x, double *mu, double *invSigmaR,
        double* pdf){
  int N = *NN;
  int p = *pp;
  int G = *GG;
  int g,i;
  double **invSigma   = malloc(sizeof(double*)*G);
  for(g=0; g<G; g++){
    invSigma[g]     = malloc(sizeof(double)*p*p);
    for(i=0; i<p*p; i++){
      invSigma[g][i] = invSigmaR[g*p*p +i];
    }
  } 
  get_PX(N, p, x, G, mu, invSigma, pdf);
  for(g=0; g<G; g++) {
    free(invSigma[g]);
  }
  free(invSigma);  
}


void RestepU(int* group, int *NN, int *pp, int *GG,double *x, double *mu, double *invSigmaR, 
             double *prior){
  int N = *NN;
  int p = *pp;
  int G = *GG;
  int i;
  double *pdf = (double*)malloc(sizeof(double)*N*G);
  double *z = (double*)malloc(sizeof(double)*N*G);
  int *lab = (int*)malloc(sizeof(int)*N);
  for(i=0; i<N; i++) lab[i]=0;
  dN(NN, pp, GG, x, mu, invSigmaR, pdf);
  estepU(N, p, G, z, prior, pdf, lab);
  get_group(G, N, z, group);
  
  free(pdf);
  free(z);
  free(lab);
}
void RllikelihoodU(double *lllvalue, int *NN, int *pp, int *GG,double *x, double *mu, 
                   double *invSigmaR, double *prior){
  int N = *NN;
  int p = *pp;
  int G = *GG;
  int i;
  int g, n;
  double *pdf = (double*)malloc(sizeof(double)*N*G);
  double **invSigma   = malloc(sizeof(double*)*G);
  for(g=0; g<G; g++){
    invSigma[g]     = malloc(sizeof(double)*p*p);
    for(i=0; i<p*p; i++){
      invSigma[g][i] = invSigmaR[g*p*p +i];
    }
  } 
  get_PX(N, p, x, G, mu, invSigma, pdf);

  // Compute llikelihood
  double rowSums=0.0;
  *lllvalue=0.0;
  for(n =0; n<N; n++){
    for(g=0; g<G; g++) rowSums += prior[g] * pdf[n + N*g];
    *lllvalue += log(rowSums);
    rowSums = 0.0;
  }
}   

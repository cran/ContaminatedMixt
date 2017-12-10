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

void get_group(int G, int N, double* z, int* group){
  int n,g,max_l=0;
  double max;
  for(n =0; n<N; n++){
  max=0;
  for(g=0; g<G; g++) 
    if(z[n + N*g] > max){
      max = z[n + N*g];
      max_l = g;
    } 
    group[n] = max_l+1;
  }
}


void mahalanobis2(int g, int N, int p, double *x, int G, double *mu, double *cov, double *delta)
{
  int i, j, n;
  double sum, insum;
  double *inv = (double *) malloc(sizeof(double)*p*p);
  for (n = 0; n < N; n++){ 
    sum = 0.0;
    for(j = 0; j < p; j++){
      insum = 0;
      for(i = 0; i < p; i++)
        insum += (x[n+ N*i] - mu[g+G*i])*cov[i+j*p];
      sum += insum*(x[n+ N*j] - mu [g+G*j]);
    }
    delta[n+N*g] = sum;
  }
  free(inv);
}


struct gmaxpar {
  int g;
  int N;
  double* z;
  double* v;
}; 
double gmax(double par, void *ex){
  double alpha = par;
  struct gmaxpar *ex1 = ex; 
  double* v = ex1 ->v;
  double* z = ex1 ->z;
  int g = ex1 ->g;
  int  N = ex1 -> N;
  int i;
  double sum =0;
  for(i=0; i < N; i++) 
    sum += z[i+g*N]*(v[i+g*N]*log(alpha)+(1-v[i+g*N])*log(1-alpha));
  return -sum;
}
void get_weights (int N, int p, int G,double* z, double* alphafix, double* alphamin,
                  double*v,double* eta, double* prior, double* alpha, double* fact){
  double* zv = malloc(sizeof(double)*N*G);
  double* colsums_z = malloc(sizeof(double)*G);
  double* colsums_zv = malloc(sizeof(double)*G);
  int g,i;
  
  //  prior    
  for(g = 0; g < G; g++){
    colsums_z[g] = 0;
    colsums_zv[g] = 0;
    for(i=0; i < N; i++) {
      colsums_z[g] += z[i+g*N];
      zv[i+g*N] = z[i+g*N]*v[i+g*N];
      colsums_zv[g] += zv[i+g*N];
    }
    prior[g] = colsums_z[g]/N;
  }
  
  // alpha
  if(alphafix[0] != -1)
    for(g = 0; g < G; g++)
      for(i=0; i < N; i++){ 
        alpha[g]= colsums_zv[g]/colsums_z[g];
        if (alpha[g] < alphamin[g]) alpha[g]= alphamin[g];
      }
  // fact    =correctionx   
  for(g = 0; g < G; g++)
    for(i=0; i < N; i++) 
      fact[i+g*N] = v[i+g*N] + ((double) 1.0 -v[i+g*N])/eta[g];
  
  free(zv);
  free(colsums_z);
  free(colsums_zv);
}

void CovarianceCN(int N, int p, int G, double *x, double *z, double *mu, int g, double *sampcov, double *fact, double *Wt, double *zfact) {
  int j,k,i;
  for(j=0; j < p; j++) {
    for(k=0; k < p; k++) {
      sampcov[j + k*p] = 0;
      for(i=0; i < N; i++) {
        sampcov[j + k*p] += Wt[i + N*g]* fact[i + N*g]*(x[i+N*j]-mu[g+G*j])*(x[i+N*k]-mu[g+G*k]);
      }
    }
  }
}

void get_zfact (int N, int p, int G, double *z, double *fact,double *Wt, double *zfact) {
  int g,i;
  double sum, sum_zfact;
  for(g = 0; g < G; g++){
    sum = 0;
    sum_zfact= 0;
    for(i=0; i < N; i++) {
      Wt[i + N*g] = z[i + N*g];
      sum  += Wt[i + N*g];
      
      zfact[i + N*g] = z[i + N*g] * fact[i + N*g];
      sum_zfact  += zfact[i + N*g];
    } 
    for(i=0; i < N; i++) {
      Wt[i + N*g] /= sum;
      zfact[i + N*g] /= sum_zfact;
    }
  }
}

void mstepC (double *x, int N, int p, int G, double *z, double *mu, double **sampcov,  double **Sigma, double **invSigma, double mtol, int mmax, double *D, char **covtype,double *fact, double *zfact)
{ int g;
  double *logdet      = malloc(sizeof(double)*G);
  double *pi = malloc(sizeof(double)*G);
  double *Wt = (double*)malloc(sizeof(double)*N*G);
  get_zfact (N, p, G, z, fact, Wt, zfact); 
  get_mu (p, G, N, x, zfact, mu);   
  for(g = 0; g < G; g++)
    CovarianceCN(N, p, G, x, z, mu, g, sampcov[g],fact, Wt, zfact);
  get_pi(N, G, z, pi);
  modeltype(p, pi, G, D, sampcov, Sigma, invSigma, logdet, p, mmax,covtype);
  free(Wt);
  free(pi);
  free(logdet);
}
void eta_max(int N, int p, int G, double* x, double *z, double* zfact, double *mu, 
             double** invSigma, double* v, double* eta){
  double *delta = malloc(sizeof(double)*N*G);
  double Aj = 0;
  double Bj = 0;
  int n,g = 0;
  
  for(g=0; g<G; g++) mahalanobis(g, N, p, x, zfact, G, mu, invSigma[g], delta);
  for(g=0; g<G; g++){
    for(n =0; n<N; n++){
      Aj += z[n+g*N]*(1.0-v[n+g*N]);
      Bj += z[n+g*N]*(1.0-v[n+g*N])*delta[n+N*g];
    }
    if(Bj/(p*Aj) > 1) eta[g] = Bj/(p*Aj);
    else eta[g] = 1;
  }
  free(delta);
}

void get_delta(int N, int p,double *x,int G,double *mu,
               double **invSigma, double*delta){
  int g;
  for(g=0; g<G; g++ )
  mahalanobis2(g, N, p, x, G, mu, invSigma[g], delta);
}
void get_PX(int N, int p,double *x,int G, double *mu, double **invSigma,double* PX){
  double *delta = malloc(sizeof(double)*N*G);
  double *logdet = malloc(sizeof(double)*G);
  int g,n;
  get_delta(N, p, x, G, mu, invSigma, delta);
  //Rprintf("delta %e \n",delta[0]);
  for(g=0; g<G; g++) {
    determinant(invSigma[g],p,p,&logdet[g]);
    logdet[g] = log(1/logdet[g]);
  }
  //Rprintf("logdet %e \n",logdet[0]);
  for(g=0; g<G; g++) {  
    for(n =0; n<N; n++){
      PX[n + N*g] = exp(- ((double)1.0/2)*delta[n + N*g]-((double)1.0/2.0)*logdet[g] - ((double)p/(double)2.0)*log(2*M_PI));
    }
  }
  free(logdet);
  free(delta);
}
void get_PXbad(int N, int p,double *x,int G, double *mu, double **invSigma,double *eta,double* PX){
  int i,g;
  double **invSigma_eta   = malloc(sizeof(double*)*G);
  for(g=0; g<G; g++){
    invSigma_eta[g]  = malloc(sizeof(double)*p*p);
    for(i=0; i<p*p; i++) invSigma_eta[g][i] = 1/eta[g]*invSigma[g][i];
  }
  // get PXbad  
  get_PX(N, p, x, G, mu, invSigma_eta, PX);
  
  for(g=0; g<G; g++) free(invSigma_eta[g]);
  free(invSigma_eta);
}
void density2( int N, int p, int G, double *z, double *prior, 
               double *eta, double **invSigma, double *mu,double *x,  
               double *fact, double *alpha, int *lab, 
               double* v, double* lllvalue, double*pdf){
  int g, n;
  double *PXgood = (double*)malloc(sizeof(double)*N*G);
  double *PXbad = (double*)malloc(sizeof(double)*N*G);
  double *Wt = (double*)malloc(sizeof(double)*N*G);
  double *zfact = (double*)malloc(sizeof(double)*N*G);
 
  get_zfact (N, p, G, z, fact, Wt, zfact);
  get_mu ( p,  G,  N, x, zfact, mu);
  get_PX(N, p, x, G, mu, invSigma, PXgood);
  get_PXbad(N, p, x, G, mu, invSigma,eta, PXbad);

// Compute pdf  
  for(g=0; g<G; g++) 
    for(n =0; n<N; n++)
      pdf[n + N*g] = alpha[g] * PXgood[n + N*g] + (1.0-alpha[g]) * PXbad[n + N*g];
  

// Compute llikelihood
  double rowSums=0.0;
  double llvalue=0.0;
  for(n =0; n<N; n++){
    for(g=0; g<G; g++) {
      if (lab[n] ==0){
        rowSums += prior[g] * pdf[n + N*g];
       } else
         llvalue += z[n+N*g] * log(prior[g] * pdf[n + N*g]);
    }
    if (lab[n] ==0) llvalue += log(rowSums);
    rowSums = 0.0;
  }
  *lllvalue = llvalue;

// update v
double a = 0;
double b = 0;
for(g=0; g<G; g++) {
  for(n =0; n<N; n++){
    a = PXgood[n + N*g]* (*alpha);
    b = PXbad[n + N*g]* (1-(*alpha));
    if( (a+b) != 0) v[n + N*g] = a/(a+b);
    else v[n + N*g] = 0;
  }
}

  free(zfact);
  free(Wt);
  free(PXgood);
  free(PXbad);
}


int stopcrit(int G, int maxiter, int cc, double loglik[3],double llvalue, 
             double threshold){    //Aitken's Stopping Criterion
  int exit=0;
  double a,b;
  
  if (cc==maxiter){
    exit = 1;
  } else{
    loglik[2] = loglik[1];loglik[1] = loglik[0];loglik[0]=llvalue;
    if (cc>2){
      if(fabs(loglik[1]-loglik[2])==0) exit = 1;
      else{
        a = (loglik[0]-loglik[1])/(loglik[1]-loglik[2]);
        b = loglik[1]+(1/(1-a)*(loglik[0]-loglik[1])); 
        if(fabs(b-loglik[0])< threshold) exit = 1;
      }
    }
  }
  return(exit);
}
void estepC(int N, int p, int G, double *z, double *prior,double* pdf, int *lab,
            double* alpha){
  double* znum  = malloc(sizeof(double)*G);
  double* zrowsum  = malloc(sizeof(double)*N);
  double  znumrowsum = 0;
  int n,g;
  for(n =0; n<N; n++){ 
    znumrowsum = 0;
    zrowsum[n] = 0;
    for(g=0; g<G; g++) {
      znum[g] = prior[g]*pdf[n+g*N];
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


void loopC (int *NN, int *pp, int *GG, double *z, 
            double *sigmar, double *invsigmar, double *mu, double *mmtol, 
            int *mmmax, double *x, int *lab, char **covtype, 
            int *maxiter, double* threshold,double* prior,int* iteration,
            double *lllvalue, double *obslll, int* group, double *v,
            double* eta, double* alpha, double* alphafix, double* alphamin){
  int g, i, cc=0;
  int N = *NN;
  int p = *pp;
  int G = *GG;
  int exit = 0;
  double *D = (double*)malloc(sizeof(double)*p*p);
  double *zfact = (double*)malloc(sizeof(double)*N*G);
  double *pdf = (double*)malloc(sizeof(double)*N*G);
  double *fact = (double*)malloc(sizeof(double)*N*G);
  double **Sigma      = malloc(sizeof(double*)*G);
  double **invSigma   = malloc(sizeof(double*)*G);
  double **sampcov    = malloc(sizeof(double*)*G);
  double loglik[] ={0,0,0};

  //FILE *f = fopen(*covtype, "a");
  //fprintf(f, "Some text: %s\n", *covtype);
  
  
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
    Rprintf("*");
    // get prior, alpha, fact 
    get_weights(N, p, G, z, alphafix, alphamin, v, eta, prior, alpha, fact);
    // get Sigma, invSigma, zfact
    mstepC(x, N, p, G, z, mu, sampcov, Sigma, invSigma, *mmtol, *mmmax, D, covtype,fact, zfact);  
    // Inflation parameters eta
    eta_max(N, p, G,x,z,zfact,mu,invSigma,v, eta);
    
    // updates v,lllvalue
    density2(N,p,G,z,prior,eta,invSigma,mu,x,fact,alpha,lab,v,lllvalue,pdf);
    //Aitken's Stopping Criterion
    exit = stopcrit(G, *maxiter, cc, loglik,*lllvalue, *threshold);
    
    //E-Step updates z
     estepC(N, p, G, z, prior, pdf, lab,alpha);
  }

  // Compute observed llikelihood
  double rowSums=0.0;
  *obslll = 0;
  int n;
  for(n =0; n<N; n++){
    for(g=0; g<G; g++) rowSums += prior[g] * pdf[n + N*g];
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
  free(zfact);
  free(pdf);
  free(fact);
  free(Sigma);
  free(invSigma); 
  free(sampcov);

  //fclose(f);
}
void RdCN(int *NN, int *pp, int *GG,double *x, double *mu, double *invSigmaR, double* eta,
            double* alpha, double* pdf){
  int N = *NN;
  int p = *pp;
  int G = *GG;
  int n,g,i;
  double **invSigma   = malloc(sizeof(double*)*G);
  double *PXgood = (double*)malloc(sizeof(double)*N*G);
  double *PXbad = (double*)malloc(sizeof(double)*N*G);

  for(g=0; g<G; g++){
    invSigma[g]     = malloc(sizeof(double)*p*p);
    for(i=0; i<p*p; i++){
      invSigma[g][i] = invSigmaR[g*p*p +i];
    }
  } 

  get_PX(N, p, x, G, mu, invSigma, PXgood);
  get_PXbad(N, p, x, G, mu, invSigma,eta, PXbad);
  for(g=0; g<G; g++) 
    for(n =0; n<N; n++)
      pdf[n + N*g] = alpha[g] * PXgood[n + N*g] + (1-alpha[g]) * PXbad[n + N*g];

  for(g=0; g<G; g++) {
    free(invSigma[g]);
  }
  free(invSigma);  
  free(PXgood);
  free(PXbad);
}

void RestepC(int* group, int *NN, int *pp, int *GG,double *x, double *mu, double *invSigmaR, double* eta,
             double* alpha,double *prior){
  int N = *NN;
  int p = *pp;
  int G = *GG;
  int i;
  double *pdf = (double*)malloc(sizeof(double)*N*G);
  double *z = (double*)malloc(sizeof(double)*N*G);
  int *lab = (int*)malloc(sizeof(int)*N);
  for(i=0; i<N; i++) lab[i]=0;
  RdCN(NN, pp, GG, x, mu, invSigmaR, eta, alpha, pdf);
  estepC(N, p, G, z, prior, pdf, lab,alpha);
  get_group(G, N, z, group);
  
  free(pdf);
  free(z);
  free(lab);
}

void RllikelihoodC(double *lllvalue, int *NN, int *pp, int *GG,double *x, double *mu, double *invSigmaR, double* eta,
             double* alpha,double *prior){
  int N = *NN;
  int p = *pp;
  int G = *GG;
  int i;
  int g, n;
  double *PXgood = (double*)malloc(sizeof(double)*N*G);
  double *PXbad = (double*)malloc(sizeof(double)*N*G);
  double *pdf = (double*)malloc(sizeof(double)*N*G);
  double **invSigma   = malloc(sizeof(double*)*G);
  for(g=0; g<G; g++){
    invSigma[g]     = malloc(sizeof(double)*p*p);
    for(i=0; i<p*p; i++){
      invSigma[g][i] = invSigmaR[g*p*p +i];
    }
  } 
  get_PX(N, p, x, G, mu, invSigma, PXgood);
  get_PXbad(N, p, x, G, mu, invSigma,eta, PXbad);
  
  // Compute pdf  
  for(g=0; g<G; g++) 
    for(n =0; n<N; n++)
      pdf[n + N*g] = alpha[g] * PXgood[n + N*g] + (1.0-alpha[g]) * PXbad[n + N*g];

  // Compute llikelihood
  double rowSums=0.0;
  *lllvalue=0.0;
  for(n =0; n<N; n++){
    for(g=0; g<G; g++) rowSums += prior[g] * pdf[n + N*g];
    *lllvalue += log(rowSums);
    rowSums = 0.0;
  }
}  

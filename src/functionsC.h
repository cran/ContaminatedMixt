int stopcrit(int G, int maxiter, int cc, double loglik[3],double llvalue, 
             double threshold);
double local_min2 ( double a, double b, double eps, double t, 
                    double f ( double x, void* ex ), double *x, void* ex );
double r8_epsilon ( void );
void get_weights (int N, int p, int G,double* z, double* alphafix,double* alphamin, double*v,double* eta, double* prior, double* alpha, double* fact);
void get_group(int G, int N, double* z, int* group);
void get_PX(int N, int p,double *x,int G, double *mu, double **invSigma,double* PX);
  

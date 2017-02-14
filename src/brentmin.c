# include <stdlib.h>
# include <stdio.h>
# include <math.h>
# include <time.h>
/*
This is a modified version of the LOCAL_MIN function published by John Burkardt
Original FORTRAN77 version by Richard Brent.
Our modifications add argument ex, that is used for passing additional parameters to f.
*/
/******************************************************************************/

double local_min2 ( double a, double b, double eps, double t, 
                   double f ( double x, void* ex ), double *x, void *ex )
  
/******************************************************************************/
/*
  Purpose:
    LOCAL_MIN seeks a local minimum of a function F(X) in an interval [A,B].
  Discussion:
    The method used is a combination of golden section search and
    successive parabolic interpolation.  Convergence is never much slower
    than that for a Fibonacci search.  If F has a continuous second
    derivative which is positive at the minimum (which is not at A or
    B), then convergence is superlinear, and usually of the order of
    about 1.324....
    The values EPS and T define a tolerance TOL = EPS * abs ( X ) + T.
    F is never evaluated at two points closer than TOL.  
    If F is a unimodal function and the computed values of F are always
    unimodal when separated by at least SQEPS * abs ( X ) + (T/3), then
    LOCAL_MIN approximates the abscissa of the global minimum of F on the 
    interval [A,B] with an error less than 3*SQEPS*abs(LOCAL_MIN)+T.  
    If F is not unimodal, then LOCAL_MIN may approximate a local, but 
    perhaps non-global, minimum to the same accuracy.
    Thanks to Jonathan Eggleston for pointing out a correction to the 
    golden section step, 01 July 2013.
  Licensing:
    This code is distributed under the GNU LGPL license. 
  Modified:
    01 July 2013
  Author:
    Original FORTRAN77 version by Richard Brent.
    C version by John Burkardt.
  Reference:
    Richard Brent,
    Algorithms for Minimization Without Derivatives,
    Dover, 2002,
    ISBN: 0-486-41998-3,
    LC: QA402.5.B74.
  Parameters:
    Input, double A, B, the endpoints of the interval.
    Input, double EPS, a positive relative error tolerance.
    EPS should be no smaller than twice the relative machine precision,
    and preferably not much less than the square root of the relative
    machine precision.
    Input, double T, a positive absolute error tolerance.
    Input, double F ( double x ), a user-supplied
    function whose local minimum is being sought.
    Output, double *X, the estimated value of an abscissa
    for which F attains a local minimum value in [A,B].
    Output, double LOCAL_MIN, the value F(X).
*/  
  
{
  double c;
  double d = 0;
  double e;
  double fu;
  double fv;
  double fw;
  double fx;
  double m;
  double p;
  double q;
  double r;
  double sa;
  double sb;
  double t2;
  double tol;
  double u;
  double v;
  double w;
  /*
  C is the square of the inverse of the golden ratio.
  */
  c = 0.5 * ( 3.0 - sqrt ( 5.0 ) );
  
  sa = a;
  sb = b;
  *x = sa + c * ( b - a );
  w = *x;
  v = w;
  e = 0.0;
  fx = f ( *x , ex);
  fw = fx;
  fv = fw;
  
  for ( ; ; )
  { 
    m = 0.5 * ( sa + sb ) ;
    tol = eps * fabs ( *x ) + t;
    t2 = 2.0 * tol;
    /*
    Check the stopping criterion.
    */
    if ( fabs ( *x - m ) <= t2 - 0.5 * ( sb - sa ) )
    {
      break;
    }
    /*
    Fit a parabola.
    */
    r = 0.0;
    q = r;
    p = q;
    
    if ( tol < fabs ( e ) )
    {
      r = ( *x - w ) * ( fx - fv );
      q = ( *x - v ) * ( fx - fw );
      p = ( *x - v ) * q - ( *x - w ) * r;
      q = 2.0 * ( q - r );
      if ( 0.0 < q )
      {
        p = - p;
      }
      q = fabs ( q );
      r = e;
      e = d;
    }
    
    if ( fabs ( p ) < fabs ( 0.5 * q * r ) && 
         q * ( sa - *x ) < p && 
         p < q * ( sb - *x ) )
    {
      /*
      Take the parabolic interpolation step.
      */
      d = p / q;
      u = *x + d;
      /*
      F must not be evaluated too close to A or B.
      */
      if ( ( u - sa ) < t2 || ( sb - u ) < t2 )
      {
        if ( *x < m )
        {
          d = tol;
        }
        else
        {
          d = - tol;
        }
      }
    }
    /*
    A golden-section step.
    */
    else
    {
      if ( *x < m )
      {
        e = sb - *x;
      }
      else
      {
        e = sa - *x;
      }
      d = c * e;
    }
    /*
    F must not be evaluated too close to X.
    */
    if ( tol <= fabs ( d ) )
    {
      u = *x + d;
    }
    else if ( 0.0 < d )
    {
      u = *x + tol;
    }
    else
    {
      u = *x - tol;
    }
    
    fu = f ( u, ex );
    /*
    Update A, B, V, W, and X.
    */
    if ( fu <= fx )
    {
      if ( u < *x )
      {
        sb = *x;
      }
      else
      {
        sa = *x;
      }
      v = w;
      fv = fw;
      w = *x;
      fw = fx;
      *x = u;
      fx = fu;
    }
    else
    {
      if ( u < *x )
      {
        sa = u;
      }
      else
      {
        sb = u;
      }
      
      if ( fu <= fw || w == *x )
      {
        v = w;
        fv = fw;
        w = u;
        fw = fu;
      }
      else if ( fu <= fv || v == *x || v == w )
      {
        v = u;
        fv = fu;
      }
    }
  }
  return fx;
}
/******************************************************************************/
/******************************************************************************/

double r8_epsilon ( )
  
/******************************************************************************/
  /*
  Purpose:
  
  R8_EPSILON returns the R8 round off unit.
  
  Discussion:
  
  R8_EPSILON is a number R which is a power of 2 with the property that,
  to the precision of the computer's arithmetic,
  1 < 1 + R
  but
  1 = ( 1 + R / 2 )
  
  Licensing:
  
  This code is distributed under the GNU LGPL license.
  
  Modified:
  
  01 September 2012
  
  Author:
  
  John Burkardt
  
  Parameters:
  
  Output, double R8_EPSILON, the R8 round-off unit.
  */
{
  const double value = 2.220446049250313E-016;
  
  return value;
}
/******************************************************************************/

/* Esempio
double eps;
double t;
eps = 10.0 *sqrt(r8_epsilon ( ));
t = eps;//
double X;
// other parameters to be passed to the objective function
struct gmaxpar ex;
ex.N = N;
ex.z = z;
ex.v = v;
for(g = 0; g < G; g++){
  ex.g = g;
  local_min2 (alphamin[g], 1,eps,t, gmax, &X, &ex );
  alpha[g] = X;   
}
*/

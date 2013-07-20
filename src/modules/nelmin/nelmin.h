#ifndef NELMIN_H
#define NELMIN_H


#include <stdlib.h>
//#include <vector.h>

//using namespace std;

// Nelder-Mead Minimization Algorithm ASA047
// from the Applied Statistics Algorithms available
// in STATLIB. Adapted from the C version by J. Burkhardt
// http://people.sc.fsu.edu/~jburkardt/c_src/asa047/asa047.html
//
/*
  Purpose:

    NELMIN minimizes a function using the Nelder-Mead algorithm.

  Discussion:

    This routine seeks the minimum value of a user-specified function.

    Simplex function minimisation procedure due to Nelder+Mead(1965),
    as implemented by O'Neill(1971, Appl.Statist. 20, 338-45), with
    subsequent comments by Chambers+Ertel(1974, 23, 250-1), Benyon(1976,
    25, 97) and Hill(1978, 27, 380-2)

    The function to be minimized must be defined by a function of
    the form

      function fn ( x, f )
      double fn
      double x(*)

    and the name of this subroutine must be declared EXTERNAL in the
    calling routine and passed as the argument FN.

    This routine does not include a termination test using the
    fitting of a quadratic surface.

  Licensing:

    This code is distributed under the GNU LGPL license. 

  Modified:

    28 October 2010

  Author:

    Original FORTRAN77 version by R ONeill.
    C version by John Burkardt.

  Reference:

    John Nelder, Roger Mead,
    A simplex method for function minimization,
    Computer Journal,
    Volume 7, 1965, pages 308-313.

    R ONeill,
    Algorithm AS 47:
    Function Minimization Using a Simplex Procedure,
    Applied Statistics,
    Volume 20, Number 3, 1971, pages 338-345.

  Parameters:

    Input, double FN ( double x[] ), the name of the routine which evaluates
    the function to be minimized.

    Input, int N, the number of variables.

    Input/output, double START[N].  On input, a starting point
    for the iteration.  On output, this data may have been overwritten.

    Output, double XMIN[N], the coordinates of the point which
    is estimated to minimize the function.

    Output, double YNEWLO, the minimum value of the function.

    Input, double REQMIN, the terminating limit for the variance
    of function values.

    Input, double STEP[N], determines the size and shape of the
    initial simplex.  The relative magnitudes of its elements should reflect
    the units of the variables.

    Input, int KONVGE, the convergence check is carried out 
    every KONVGE iterations.

    Input, int KCOUNT, the maximum number of function 
    evaluations.

    Output, int *ICOUNT, the number of function evaluations 
    used.

    Output, int *NUMRES, the number of restarts.

    Output, int *IFAULT, error indicator.
    0, no errors detected.
    1, REQMIN, N, or KONVGE has an illegal value.
    2, iteration terminated because KCOUNT was exceeded without convergence.
*/


//double objFunc(double *params, double *match_pairs, int N, double *poses_1, double *poses_2, int numPoses, double uHigh, double uLow, double u1) {
// int numPairs;
// int numPoses;

void nelmin ( int n, double start[], double xmin[], 
	      double *ynewlo, double reqmin, double step[], int konvge, int kcount, 
	      int *icount, int *numres, int *ifault,
	      int numPoses, int numPairs, double *d_matchPairs, double *d_offset, double *d_sum, double *poses_1, double *poses_2,
	      double uHigh, double uLow, double u1
	      )
{

  //double (*fn)(double*);

  const double ccoeff = 0.5;
  const double ecoeff = 2.0;
  const double eps = 0.001;
  const double rcoeff = 1.0;
  int ihi,ilo,jcount,l,nn;
  double del,dn,dnn;
  double rq,x,y2star,ylo,ystar,z;

  double *p;
  double *pstar;
  double *p2star;
  double *pbar;
  double *y;
  double offset[3];


  int i, j;

  //  Check the input parameters.
  if ( reqmin <= 0.0 ) { *ifault = 1; return; }
  if ( n < 1 ) { *ifault = 1; return; }
  if ( konvge < 1 ) { *ifault = 1; return; }

  //vector<double> p(n*(n+1));
  //vector<double> pstar(n);
  //vector<double> p2star(n);
  //vector<double> pbar(n);
  //vector<double> y(n+1);

  p = (double*) malloc(n*(n+1)*sizeof(double));
  pstar = (double*) malloc(n*sizeof(double));
  p2star = (double*) malloc(n*sizeof(double));
  pbar = (double*) malloc(n*sizeof(double));
  y = (double*) malloc((n+1)*sizeof(double));

  *icount = 0;
  *numres = 0;

  jcount = konvge; 
  dn = ( double ) ( n );
  nn = n + 1;
  dnn = ( double ) ( nn );
  del = 1.0;
  rq = reqmin * dn;

  //  Initial or restarted loop.
  for ( ; ; ) {

    for (i = 0; i < n; i++ ) {
      p[i+n*n] = start[i];
    }

    // CHECK1

    //y[n] = (*fn)( start );

    y[n] = objFunc(start, d_matchPairs, d_offset, d_sum, numPairs, poses_1, poses_2, numPoses, uHigh, uLow, u1, offset); 
	//printf("y[%d] = %lf\n", n, y[n]);

    *icount = *icount + 1;
    
    // CHECK2
   
    for (j = 0; j < n; j++ ) {
      x = start[j];
      start[j] = start[j] + step[j] * del;

      for (i = 0; i < n; i++ ) {
	p[i+j*n] = start[i];
      }

      //y[j] = (*fn)( start );
      y[j] = objFunc(start, d_matchPairs, d_offset, d_sum, numPairs, poses_1, poses_2, numPoses, uHigh, uLow, u1, offset);
	  //printf("y[%d] = %lf\n", j, y[j]);
      *icount = *icount + 1;
      start[j] = x;
    }

    // CHECK3

    //  The simplex construction is complete.
    //                    
    //  Find highest and lowest Y values.  YNEWLO = Y(IHI) indicates
    //  the vertex of the simplex to be replaced.
    ylo = y[0];
    ilo = 0;
    
    for (i = 1; i < nn; i++ ) {
      if ( y[i] < ylo ) {
	ylo = y[i]; ilo = i;
      }
    }

    //  Inner loop.
    for ( ; ; ) {
      if ( kcount <= *icount ) { break; }
      *ynewlo = y[0];

      ihi = 0;
      
      for (i = 1; i < nn; i++ ) {
        if ( *ynewlo < y[i] ) { *ynewlo = y[i]; ihi = i; }
      }
      //  Calculate PBAR, the centroid of the simplex vertices
      //  excepting the vertex with Y value YNEWLO.
      for (i = 0; i < n; i++ ) {
        z = 0.0;
        for (j = 0; j < nn; j++ ) { z = z + p[i+j*n]; }
        z = z - p[i+ihi*n];  
        pbar[i] = z / dn;
      }
      //  Reflection through the centroid.
      for (i = 0; i < n; i++ ) {
        pstar[i] = pbar[i] + rcoeff * ( pbar[i] - p[i+ihi*n] );
      }

      // CHECK4
      //ystar = (*fn)( &pstar[0] );
      ystar = objFunc(&pstar[0], d_matchPairs, d_offset, d_sum, numPairs, poses_1, poses_2, numPoses, uHigh, uLow, u1, offset); 
      *icount = *icount + 1;
      //  Successful reflection, so extension.
      if ( ystar < ylo ) {
        for (i = 0; i < n; i++ ) {
          p2star[i] = pbar[i] + ecoeff * ( pstar[i] - pbar[i] );
        }
        //y2star = (*fn)( &p2star[0] );
        y2star = objFunc(&p2star[0], d_matchPairs, d_offset, d_sum, numPairs, poses_1, poses_2, numPoses, uHigh, uLow, u1, offset); 
        *icount = *icount + 1;
	//  Check extension.
        if ( ystar < y2star ) {
          for (i = 0; i < n; i++ ) { p[i+ihi*n] = pstar[i]; }
          y[ihi] = ystar;
        } else { //  Retain extension or contraction.
          for (i = 0; i < n; i++ ) { p[i+ihi*n] = p2star[i]; }
          y[ihi] = y2star;
        }
      } else { //  No extension.
        l = 0;
        for (i = 0; i < nn; i++ ) {
	  if ( ystar < y[i] ) l += 1;
        }
	
        if ( 1 < l ) {
          for (i = 0; i < n; i++ ) { p[i+ihi*n] = pstar[i]; }
          y[ihi] = ystar;
        }
	//  Contraction on the Y(IHI) side of the centroid.
        else if ( l == 0 ) {
          for (i = 0; i < n; i++ ) {
            p2star[i] = pbar[i] + ccoeff * ( p[i+ihi*n] - pbar[i] );
          }
          //y2star = (*fn)( &p2star[0] );
          y2star = objFunc(&p2star[0], d_matchPairs, d_offset, d_sum, numPairs, poses_1, poses_2, numPoses, uHigh, uLow, u1, offset); 
          *icount = *icount + 1;
	  //  Contract the whole simplex.
          if ( y[ihi] < y2star ) {
            for (j = 0; j < nn; j++ ) {
              for (i = 0; i < n; i++ ) {
                p[i+j*n] = ( p[i+j*n] + p[i+ilo*n] ) * 0.5;
                xmin[i] = p[i+j*n];
              }
              //y[j] = (*fn)( xmin );
              y[j] = objFunc(xmin, d_matchPairs, d_offset, d_sum, numPairs, poses_1, poses_2, numPoses, uHigh, uLow, u1, offset); 
              *icount = *icount + 1;
            }
            ylo = y[0];
            ilo = 0;
	    
            for (i = 1; i < nn; i++ ) {
              if ( y[i] < ylo ) { ylo = y[i]; ilo = i; }
            }
            continue;
          }
	  //  Retain contraction.
          else {
            for (i = 0; i < n; i++ ) {
              p[i+ihi*n] = p2star[i];
            }
            y[ihi] = y2star;
          }
        }
	//  Contraction on the reflection side of the centroid.
        else if ( l == 1 ) {
          for (i = 0; i < n; i++ ) {
            p2star[i] = pbar[i] + ccoeff * ( pstar[i] - pbar[i] );
          }
          //y2star = (*fn)( &p2star[0] );
          y2star = objFunc(&p2star[0], d_matchPairs, d_offset, d_sum, numPairs, poses_1, poses_2, numPoses, uHigh, uLow, u1, offset); 
          *icount = *icount + 1;
	  //
	  //  Retain reflection?
	  //
          if ( y2star <= ystar ) {
            for (i = 0; i < n; i++ ) { p[i+ihi*n] = p2star[i]; }
            y[ihi] = y2star;
          }
          else {
            for (i = 0; i < n; i++ ) { p[i+ihi*n] = pstar[i]; }
            y[ihi] = ystar;
          }
        }
      }
      //  Check if YLO improved.
      if ( y[ihi] < ylo ) { ylo = y[ihi]; ilo = ihi; }
      jcount = jcount - 1;
      
      if ( 0 < jcount ) { continue; }
      //  Check to see if minimum reached.
      if ( *icount <= kcount ) {
        jcount = konvge;
	
        z = 0.0;
        for (i = 0; i < nn; i++ ) { z = z + y[i]; }
        x = z / dnn;
	
        z = 0.0;
        for (i = 0; i < nn; i++ ) {
          z = z + pow ( y[i] - x, 2 );
        }
	
        if ( z <= rq ) {break;}
      }
    }
    //  Factorial tests to check that YNEWLO is a local minimum.
    for (i = 0; i < n; i++ ) { xmin[i] = p[i+ilo*n]; }
    *ynewlo = y[ilo];
    
    if ( kcount < *icount ) { *ifault = 2; break; }

    *ifault = 0;

    for (i = 0; i < n; i++ ) {
      del = step[i] * eps;
      xmin[i] = xmin[i] + del;
      //z = (*fn)( xmin );
      z = objFunc(xmin, d_matchPairs, d_offset, d_sum, numPairs, poses_1, poses_2, numPoses, uHigh, uLow, u1, offset); 
      *icount = *icount + 1;
      if ( z < *ynewlo ) { *ifault = 2; break; }
      xmin[i] = xmin[i] - del - del;
      //z = (*fn)( xmin );
      z = objFunc(xmin, d_matchPairs, d_offset, d_sum, numPairs, poses_1, poses_2, numPoses, uHigh, uLow, u1, offset); 
      *icount = *icount + 1;
      if ( z < *ynewlo ) { *ifault = 2; break; }
      xmin[i] = xmin[i] + del;
    }
    
    if ( *ifault == 0 ) { break; }
    //  Restart the procedure.
    for (i = 0; i < n; i++ ) { start[i] = xmin[i]; }
    del = eps;
    *numres = *numres + 1;
  }
  return;
}




#endif

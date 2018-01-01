#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "NoFLA_UTV_WY_blk_var2.h"

#define max( a, b ) ( (a) > (b) ? (a) : (b) )
#define min( a, b ) ( (a) < (b) ? (a) : (b) )


#define PRINT_DATA


// ============================================================================
// Declaration of local prototypes.

static void matrix_generate( int m_A, int n_A, double * buff_A, int ldim_A );

static void print_double_matrix( char * name, int m_A, int n_A, 
                double * buff_A, int ldim_A );

static double check_utv( int m_A, int n_A, double * buff_A, int ldim_A,
                  double * buff_U, int ldim_U, double * buff_T, int ldim_T,
                  double * buff_V, int ldim_V );

static double check_ort( int m_Q, int n_Q, double * buff_Q, int ldim_Q );

static double check_dense_sv( int m_A, int n_A, double * buff_A, int ldim_A, 
                  double * buff_T, int ldim_T );

static void compute_svd( int m_A, int n_A, double * buff_A, int ldim_A,
                double * buff_sv );


// ============================================================================
// Declaration of external prototypes.

extern double dlange_();


// ============================================================================
int main( int argc, char *argv[] ) {
  int     nb_alg, pp, m_A, n_A, ldim_A, ldim_Acopy, ldim_U, ldim_V;
  double  * buff_A, * buff_Acopy, * buff_U, * buff_V;
  double  resid;

  // Some initializations.
  m_A    =  8;
  n_A    =  8;
  // We use a small block size to factorize small input matrices, but
  // larger blocksizes such as 64 should be used for larger matrices.
  nb_alg = 3;

  // Create matrices A, Acopy, U, and V.
  buff_A     = ( double * ) malloc( m_A * n_A * sizeof( double ) );
  ldim_A     = max( 1, m_A );

  buff_Acopy = ( double * ) malloc( m_A * n_A * sizeof( double ) );
  ldim_Acopy = max( 1, m_A );

  buff_U     = ( double * ) malloc( m_A * m_A * sizeof( double ) );
  ldim_U     = max( 1, m_A );

  buff_V     = ( double * ) malloc( n_A * n_A * sizeof( double ) );
  ldim_V     = max( 1, n_A );

  // Generate matrix A.
  matrix_generate( m_A, n_A, buff_A, ldim_A );
  // Get a copy of A into Acopy.
  dlacpy_( "All", & m_A, & n_A, buff_A, & ldim_A, buff_Acopy, & ldim_Acopy );

#ifdef PRINT_DATA
  print_double_matrix( "ai", m_A, n_A, buff_A, ldim_A );
#endif

  // Factorize matrix.
  printf( "%% Just before computing factorization.\n" );
  NoFLA_UTV_WY_blk_var2( m_A, n_A, buff_A, ldim_A, 
      1, m_A, m_A, buff_U, ldim_U, 
      1, n_A, n_A, buff_V, ldim_V, 
      nb_alg, 10, 2 );
  printf( "%% Just after computing factorization.\n" );

  // Print results.
#ifdef PRINT_DATA
  print_double_matrix( "af", m_A, n_A, buff_A, ldim_A );
  print_double_matrix( "uf", m_A, m_A, buff_U, ldim_U );
  print_double_matrix( "vf", n_A, n_A, buff_V, ldim_V );
#endif

  // Compute residuals.
  printf( "\n%% Computing residuals...\n" );

  resid = check_utv( m_A, n_A, buff_Acopy, ldim_Acopy, 
                   buff_U, ldim_U, buff_A, ldim_A, buff_V, ldim_V );
  printf( "Res. || A - U * T * V' ||_F / || A ||_F                = %15.5le\n",
          resid );

  resid = check_ort( m_A, m_A, buff_U, ldim_U );
  printf( "Res. || I - U' * U ||_F / || U ||_F                    = %15.5le\n",
          resid );

  resid = check_ort( n_A, n_A, buff_V, ldim_V );
  printf( "Res. || I - V' * V ||_F / || V ||_F                    = %15.5le\n",
          resid );

  resid = check_dense_sv( m_A, n_A, buff_Acopy, ldim_Acopy, buff_A, ldim_A );
  printf( "Res. || sin.val.A - sin.val.T ||_F / || sin.val.A ||_F = %15.5le\n",
          resid );

  // Free matrices.
  free( buff_A );
  free( buff_Acopy );
  free( buff_U );
  free( buff_V );

  printf( "%% End of Program\n\n" );

  return 0;
}

// ============================================================================
static void matrix_generate( int m_A, int n_A, double * buff_A, int ldim_A ) {
  int     i, j;

  srand( 10 );
  for ( j = 0; j < n_A; j++ ) {
    for ( i = 0; i < m_A; i++ ) {
      buff_A[ i + j * ldim_A ] = ( double ) rand() / ( double ) RAND_MAX;
    }
  }
}


// ============================================================================
static void matrix_generate2( int m_A, int n_A, double * buff_A, int ldim_A ) {
  int  i, j, num;

  //
  // Matrix with integer values.
  // ---------------------------
  //
  if( ( m_A > 0 )&&( n_A > 0 ) ) {
    num = 1;
    for ( j = 0; j < n_A; j++ ) {
      for ( i = ( j % m_A ); i < m_A; i++ ) {
        buff_A[ i + j * ldim_A ] = ( double ) num;
        num++;
      }
      for ( i = 0; i < ( j % m_A ); i++ ) {
        buff_A[ i + j * ldim_A ] = ( double ) num;
        num++;
      }
    }
    if( ( m_A > 0 )&&( n_A > 0 ) ) {
      buff_A[ 0 + 0 * ldim_A ] = 1.2;
    }
#if 0
    // Scale down matrix.
    if( num == 0.0 ) {
      rnum = 1.0;
    } else {
      rnum = 1.0 / num;
    }
    for ( j = 0; j < n_A; j++ ) {
      for ( i = 0; i < m_A; i++ ) {
        buff_A[ i + j * ldim_A ] *= rnum;
      }
    }
#endif
  }
}

// ============================================================================
static void print_double_matrix( char * name, int m_A, int n_A, 
                double * buff_A, int ldim_A ) {
  int  i, j;

  printf( "%s = [\n", name );
  for( i = 0; i < m_A; i++ ) {
    for( j = 0; j < n_A; j++ ) {
      printf( "%le ", buff_A[ i + j * ldim_A ] );
    }
    printf( "\n" );
  }
  printf( "];\n" );
}

// ============================================================================
static double check_utv( int m_A, int n_A, double * buff_A, int ldim_A,
                  double * buff_U, int ldim_U, double * buff_T, int ldim_T,
                  double * buff_V, int ldim_V ) {
// Compute: || U * T * V' - A ||, where T is m-by-m.
  double  * buff_UT, * buff_UTVt;
  int     ldim_UT, ldim_UTVt;
  double  dnorm_dif, dnorm_A, resid;
  double  d_one = 1.0, d_zero = 0.0, d_minus_one = -1.0;

  // Create auxiliary objects.
  //// FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, T, & UT );
  //// FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, T, & UTVt );
  //// FLA_Obj_create( FLA_Obj_datatype( A ), 1, 1, 0, 0, & nrm );
  buff_UT   = ( double * ) malloc( m_A * n_A * sizeof( double ) );
  buff_UTVt = ( double * ) malloc( m_A * n_A * sizeof( double ) );
  ldim_UT   = m_A;
  ldim_UTVt = m_A;

  // Compute: || U * T * V' - A ||.
  //// FLA_Gemm( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE, 
  ////           FLA_ONE, U, T, FLA_ZERO, UT );
  dgemm_( "No transpose", "No transpose", & m_A, & n_A, & m_A, 
          & d_one, buff_U, & ldim_U, buff_T, & ldim_T, 
          & d_zero, buff_UT, & ldim_UT );
  //// FLA_Copy( A, UTVt );
  dlacpy_( "All", & m_A, & n_A, buff_A, & ldim_A, buff_UTVt, & ldim_UTVt );
  //// FLA_Gemm( FLA_NO_TRANSPOSE, FLA_TRANSPOSE, 
  ////           FLA_ONE, UT, V, FLA_MINUS_ONE, UTVt );
  dgemm_( "No transpose", "Transpose", & m_A, & n_A, & n_A, 
          & d_one, buff_UT, & ldim_UT, buff_V, & ldim_V,
          & d_minus_one, buff_UTVt, & ldim_UTVt );
  //// FLA_Norm_frob( UTVt, nrm );
  //// dnorm_dif = *( ( double * ) FLA_Obj_buffer_at_view( nrm ) );
  dnorm_dif = dlange_( "Frobenius", & m_A, & n_A, buff_UTVt, & ldim_UTVt );
  //// printf( "  My norm of difference: %le\n", dnorm_dif );

  // Compute the norm of A.
  //// FLA_Norm_frob( A, nrm );
  //// dnorm_A = *( ( double * ) FLA_Obj_buffer_at_view( nrm ) );
  dnorm_A = dlange_( "Frobenius", & m_A, & n_A, buff_A, & ldim_A );
  //// printf( "  My norm of A: %le\n", dnorm_A );

  // Compute and return the quotient.
  if( dnorm_A == 0.0 ) {
    resid = dnorm_dif;
  } else {
    resid = dnorm_dif / dnorm_A;
  }
  //// printf( "resid: %le\n", resid );

  // Remove auxiliary objects.
  //// FLA_Obj_free( & UT );
  //// FLA_Obj_free( & UTVt );
  //// FLA_Obj_free( & nrm );
  free( buff_UT );
  free( buff_UTVt );

  return resid;
}

// ============================================================================
static double check_ort( int m_Q, int n_Q, double * buff_Q, int ldim_Q ) {
// Compute: || In - Q' * Q || / || Q ||, where Q is m-by-n.
  //// FLA_Obj Dif;
  double  * buff_Dif;
  int     ldim_Dif;
  double  dnorm_dif, dnorm_Q, resid;
  double  d_one = 1.0, d_zero = 0.0, d_minus_one = -1.0;
  
  // Create objects Dif and nrm.
  //// n_Q = FLA_Obj_width( Q );
  //// FLA_Obj_create( FLA_Obj_datatype( Q ), n_Q, n_Q, 0, 0, & Dif );
  //// FLA_Obj_create( FLA_Obj_datatype( Q ),   1,   1, 0, 0, & nrm );
  buff_Dif = ( double * ) malloc( n_Q * n_Q * sizeof( double ) );
  ldim_Dif = n_Q;

  // Compute the norm of the difference.
  //// FLA_Set_to_identity( Dif );
  dlaset_( "All", & n_Q, & n_Q, & d_zero, & d_one, buff_Dif, & ldim_Dif );
  //// FLA_Gemm( FLA_TRANSPOSE, FLA_NO_TRANSPOSE, 
  ////           FLA_ONE, Q, Q, FLA_MINUS_ONE, Dif );
  dgemm_( "Transpose", "No transpose", & n_Q, & n_Q, & m_Q, 
          & d_one, buff_Q, & ldim_Q, buff_Q, & ldim_Q,
          & d_minus_one, buff_Dif, & ldim_Dif );
  //// FLA_Norm_frob( Dif, nrm );
  //// dnorm_dif = *( ( double * ) FLA_Obj_buffer_at_view( nrm ) );
  dnorm_dif = dlange_( "Frobenius", & m_Q, & n_Q, buff_Dif, & ldim_Dif );
  //// printf( "  My norm of difference: %le\n", dnorm_dif );

  // Compute the norm of Q.
  //// FLA_Norm_frob( Q, nrm );
  //// dnorm_Q = *( ( double * ) FLA_Obj_buffer_at_view( nrm ) );
  dnorm_Q = dlange_( "Frobenius", & m_Q, & n_Q, buff_Q, & ldim_Q );
  //// printf( "  My norm of Q: %le\n", dnorm_Q );

  // Compute and return the quotient.
  if( dnorm_Q == 0.0 ) {
    resid = dnorm_dif;
  } else {
    resid = dnorm_dif / dnorm_Q;
  }

  // Remove objects.
  //// FLA_Obj_free( & Dif );
  //// FLA_Obj_free( & nrm );
  free( buff_Dif );

  return resid;
}                                                                               

// ============================================================================
static double check_dense_sv( int m_A, int n_A, double * buff_A, int ldim_A, 
                  double * buff_T, int ldim_T ) {
// Compute: || singular values of A - singular values of T ||, 
// where A is m-by-n.
  //// FLA_Obj Acopy, Tcopy, d1, d2;
  double  * buff_d1, * buff_d2;
  int     min_mn_A;
  double  dnorm_vs, dnorm_d, resid;
  int     i_one = 1;
  double  d_minus_one = -1.0;

  // Create objects d1 and d2.
  min_mn_A  = min( m_A, n_A );
  //// FLA_Obj_create( FLA_Obj_datatype( Acopy ), min_mn_A,   1, 0, 0, & d1 );
  //// FLA_Obj_create( FLA_Obj_datatype( Acopy ), min_mn_A,   1, 0, 0, & d2 );
  buff_d1 = ( double * ) malloc( min_mn_A * sizeof( double ) );
  buff_d2 = ( double * ) malloc( min_mn_A * sizeof( double ) );

  // Compute svd.
  compute_svd( m_A, n_A, buff_A, ldim_A, buff_d1 );
  compute_svd( m_A, n_A, buff_T, ldim_T, buff_d2 );
  //// FLA_Obj_show( " d1 = [ ", d1, "%le", " ];" );
  //// FLA_Obj_show( " d2 = [ ", d2, "%le", " ];" );

  // Compute residuals.
  //// FLA_Norm_frob( d1, res );
  //// norm_vs = *( ( double * ) FLA_Obj_buffer_at_view( res ) );
  dnorm_vs = dlange_( "Frobenius", & min_mn_A, & i_one, buff_d1, & min_mn_A );
  //// printf( "  My norm of d1: %le\n", dnorm_vs );

  //// FLA_Axpy( FLA_MINUS_ONE, d1, d2 );
  daxpy_( & min_mn_A, & d_minus_one, buff_d1, & i_one, buff_d2, & i_one );

  //// FLA_Norm_frob( d2, res );
  //// norm_d = *( ( double * ) FLA_Obj_buffer_at_view( res ) );
  dnorm_d = dlange_( "Frobenius", & min_mn_A, & i_one, buff_d2, & min_mn_A );
  //// printf( "  My norm of difference: %le\n", dnorm_d );

  if( dnorm_vs != 0.0 ) {
    resid = dnorm_d / dnorm_vs;
  } else {
    resid = dnorm_d;
  }

  // Remove object d1 and d2.
  //// FLA_Obj_free( & d1 );
  //// FLA_Obj_free( & d2 );
  free( buff_d1 );
  free( buff_d2 );

  return resid;
}

// ============================================================================
static void compute_svd( int m_A, int n_A, double * buff_A, int ldim_A,
                double * buff_sv ) {
// Compute singular values of A.
  //// FLA_Obj  Acopy, Work;
  double  * buff_Acopy, * buff_Work, * buff_U, * buff_VT, dwork;
  int     info, m_Acopy, n_Acopy, ldim_Acopy, ldim_U, ldim_VT, lwork;
  char    jobu, jobvt;

  //// printf( "compute_svd\n" );

  // Make a temporal copy of A into Acopy.
  //// FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, A, & Acopy );
  //// FLA_Copy( A, Acopy );                                  
  buff_Acopy = ( double * ) malloc( m_A * n_A * sizeof( double ) );
  m_Acopy    = m_A;
  n_Acopy    = n_A;
  ldim_Acopy = m_A;
  dlacpy_( "All", & m_A, & n_A, buff_A, & ldim_A, buff_Acopy, & ldim_Acopy );

  // Some initializations.
  jobu     = 'N';
  jobvt    = 'N';
  buff_U   = NULL;
  ldim_U   = m_Acopy;
  buff_VT  = NULL;
  ldim_VT  = n_Acopy;

  // Compute optimal length of workspace.
  lwork = -1;
  dgesvd_( & jobu, & jobvt, & m_Acopy, & n_Acopy, buff_Acopy, & ldim_Acopy,
           buff_sv, buff_U, & ldim_U, buff_VT, & ldim_VT,
           & dwork, & lwork, & info );
  if( info != 0 ) {
    fprintf( stderr, " *** Info after dgesvd_: %d \n", info );
  }
  // Add some additional space just in case.
  lwork = ( int ) dwork;
  //// printf( "  Optimal lwork: %d\n", lwork );

  // Create workspace.
  //// buff_Work = ( double * ) FLA_Obj_buffer_at_view( Work );
  buff_Work = ( double * ) malloc( lwork * sizeof( double ) );

  // Call to SUBROUTINE DGESVD( JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT,
  //                            WORK, LWORK, INFO )
  dgesvd_( & jobu, & jobvt, & m_Acopy, & n_Acopy, buff_Acopy, & ldim_Acopy,
           buff_sv, buff_U, & ldim_U, buff_VT, & ldim_VT,
           buff_Work, & lwork, & info );
  if( info != 0 ) {
    fprintf( stderr, " *** Info after dgesvd_: %d \n", info );
  }

  // Remove object Work.
  //// FLA_Obj_free( & Work );
  free( buff_Work );

  // Remove object Acopy.
  //// FLA_Obj_free( & Acopy );
  free( buff_Acopy );

  //// print_double_matrix( "sv", min( m_A, n_A ), 1, buff_sv, m_A );
}



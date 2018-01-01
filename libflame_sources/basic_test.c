#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "FLAME.h"
#include "FLA_UTV_UT_blk_var1.h"

#define PRINT_DATA


// ============================================================================
// Declaration of local prototypes.

static void matrix_generate( FLA_Obj A );

static int check_utv( FLA_Obj A, FLA_Obj U, FLA_Obj T, FLA_Obj V, 
               FLA_Obj resid );

static int check_ort( FLA_Obj Q, FLA_Obj resid );

static int check_dense_sv( FLA_Obj A, FLA_Obj T, FLA_Obj resid );

static int compute_svd( FLA_Obj A, FLA_Obj sv );


// ============================================================================
int main( int argc, char *argv[] ) {
  int      m_A, n_A, nb_alg;
  FLA_Obj  A, Acopy, U, V, resid;

  // Initialize FLAME.
  FLA_Init();

  // Some initializations.
  m_A     =  8;
  n_A     =  8;
  // We use a small block size to factorize small input matrices, but
  // larger blocksizes such as 64 should be used for larger matrices.
  nb_alg  = 3;

  // Create FLAME objects.
  FLA_Obj_create( FLA_DOUBLE, m_A, n_A, 0, 0, & A );
  FLA_Obj_create( FLA_DOUBLE, m_A, n_A, 0, 0, & Acopy );
  FLA_Obj_create( FLA_DOUBLE, m_A, m_A, 0, 0, & U );
  FLA_Obj_create( FLA_DOUBLE, n_A, n_A, 0, 0, & V );
  FLA_Obj_create( FLA_DOUBLE,   1,   1, 0, 0, & resid );

  // Generate matrix.
  //// FLA_Random_matrix( A );
  matrix_generate( A );
  FLA_Copy( A, Acopy );

  // Print initial data.
#ifdef PRINT_DATA
  FLA_Obj_show( " Ai = [ ", A, "%le", " ];" );
#endif

  // Factorize matrix.
  printf( "%% Just before computing factorization.\n" );
  FLA_UTV_UT_blk_var1( A, 1, U, 1, V, nb_alg, 10, 2 );
  printf( "%% Just after computing factorization.\n" );

  // Print results.
#ifdef PRINT_DATA
  FLA_Obj_show( " Af = [ ", A, "%le", " ];" );
  FLA_Obj_show( " Uf = [ ", U, "%le", " ];" );
  FLA_Obj_show( " Vf = [ ", V, "%le", " ];" );
#endif

  // Compute residuals.
  printf( "\n%% Computing residuals...\n" );

  check_utv( Acopy, U, A, V, resid );
  FLA_Obj_show( "Res. || A - U * T * V' ||_F / || A ||_F = [ ", 
      resid, "%15.5e", " ]; " );

  check_ort( U, resid );
  FLA_Obj_show( "Res. || I - U' * U ||_F / || U ||_F = [ ", 
      resid, "%15.5e", " ]; " );

  check_ort( V, resid );
  FLA_Obj_show( "Res. || I - V' * V ||_F / || V ||_F = [ ", 
      resid, "%15.5e", " ]; " );

  check_dense_sv( Acopy, A, resid );
  FLA_Obj_show( "Res. || sin.val.A - sin.val.T ||_F / || sin.val.A ||_F = [ ", 
      resid, "%15.5e", " ]; " );

  // Free objects.
  FLA_Obj_free( & A );
  FLA_Obj_free( & Acopy );
  FLA_Obj_free( & U );
  FLA_Obj_free( & V );
  FLA_Obj_free( & resid );

  // Finalize FLAME.
  printf( "%% End of Program\n\n" );
  FLA_Finalize();


  return 0;
}

// ============================================================================
static void matrix_generate( FLA_Obj A ) {
  double  * buff_A;
  int     m_A, n_A, ldim_A;
  int     i, j;

  buff_A = ( double * ) FLA_Obj_buffer_at_view( A );
  m_A    = FLA_Obj_length( A );
  n_A    = FLA_Obj_width ( A );
  ldim_A = FLA_Obj_col_stride( A );

  srand( 10 );
  for ( j = 0; j < n_A; j++ ) {
    for ( i = 0; i < m_A; i++ ) {
      buff_A[ i + j * ldim_A ] = ( double ) rand() / ( double ) RAND_MAX;
    }
  }
}

// ============================================================================
static void matrix_generate2( FLA_Obj A ) {
  double  * buff_A, rnum;
  int     m_A, n_A, ldim_A;
  int     i, j, num;

  buff_A = ( double * ) FLA_Obj_buffer_at_view( A );
  m_A    = FLA_Obj_length( A );
  n_A    = FLA_Obj_width ( A );
  ldim_A = FLA_Obj_col_stride( A );

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
    // Scale down matrix.
#if 0
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
int check_utv( FLA_Obj A, FLA_Obj U, FLA_Obj T, FLA_Obj V, FLA_Obj resid ) {
// Compute: || A - U * T * V' ||_F, where A is m-by-m.
  FLA_Obj  UT, UTVt, nrm;
  double   dnorm_dif, dnorm_A;

  // Create auxiliary objects.
  FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, T, & UT );
  FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, T, & UTVt );
  FLA_Obj_create( FLA_Obj_datatype( A ), 1, 1, 0, 0, & nrm );

  // Compute: || A - U * T * V' ||_F.
  FLA_Gemm( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE, 
            FLA_ONE, U, T, FLA_ZERO, UT );
  FLA_Copy( A, UTVt );
  FLA_Gemm( FLA_NO_TRANSPOSE, FLA_TRANSPOSE, 
            FLA_ONE, UT, V, FLA_MINUS_ONE, UTVt );
  FLA_Norm_frob( UTVt, nrm );
  FLA_Obj_extract_real_scalar( nrm, & dnorm_dif );

  // Compute: || A ||_F.
  FLA_Norm_frob( A, nrm );
  FLA_Obj_extract_real_scalar( nrm, & dnorm_A );

  // Compute and return the quotient.
  if( dnorm_A == 0.0 ) {
    *( ( double * ) FLA_Obj_buffer_at_view( resid ) ) = dnorm_dif;
  } else {
    *( ( double * ) FLA_Obj_buffer_at_view( resid ) ) = dnorm_dif / dnorm_A;
  }

  // Remove auxiliary objects.
  FLA_Obj_free( & UT );
  FLA_Obj_free( & UTVt );
  FLA_Obj_free( & nrm );

  return FLA_SUCCESS;
}

// ============================================================================
int check_ort( FLA_Obj Q, FLA_Obj resid ) {
// Compute: || In - Q' * Q || / || Q ||, where Q is m-by-n.
  FLA_Obj  Dif, nrm;
  double   dnorm_dif, dnorm_Q;
  int      n_Q; 
  
  // Create objects Dif and nrm.
  n_Q = FLA_Obj_width( Q );
  FLA_Obj_create( FLA_Obj_datatype( Q ), n_Q, n_Q, 0, 0, & Dif );
  FLA_Obj_create( FLA_Obj_datatype( Q ),   1,   1, 0, 0, & nrm );

  // Compute the norm of the difference.
  FLA_Set_to_identity( Dif );
  FLA_Gemm( FLA_TRANSPOSE, FLA_NO_TRANSPOSE, 
            FLA_ONE, Q, Q, FLA_MINUS_ONE, Dif );
  FLA_Norm_frob( Dif, nrm );
  FLA_Obj_extract_real_scalar( nrm, & dnorm_dif );
  //// printf( "   My norm of difference: %le\n", dnorm_dif );

  // Compute the norm of Q.
  FLA_Norm_frob( Q, nrm );
  FLA_Obj_extract_real_scalar( nrm, & dnorm_Q );
  //// printf( "   My norm of Q: %le\n", dnorm_Q );

  // Compute and return the quotient.
  if( dnorm_Q == 0.0 ) {
    *( ( double * ) FLA_Obj_buffer_at_view( resid ) ) = dnorm_dif;
  } else {
    *( ( double * ) FLA_Obj_buffer_at_view( resid ) ) = dnorm_dif / dnorm_Q;
  }

  // Remove objects.
  FLA_Obj_free( & Dif );
  FLA_Obj_free( & nrm );

  return FLA_SUCCESS;                                                           
}                                                                               

// ============================================================================
static int check_dense_sv( FLA_Obj A, FLA_Obj T, FLA_Obj resid ) {
// Compute: || singular values of A - singular values of T ||, 
// where A is m-by-n.
  FLA_Obj d1, d2, nrm;
  double  dnorm_vs, dnorm_dif;
  int     min_mn_A;

  // Create auxiliary objects.
  min_mn_A = FLA_Obj_min_dim( A );
  FLA_Obj_create( FLA_Obj_datatype( A ), min_mn_A, 1, 0, 0, & d1 );
  FLA_Obj_create( FLA_Obj_datatype( A ), min_mn_A, 1, 0, 0, & d2 );
  FLA_Obj_create( FLA_Obj_datatype( A ),        1, 1, 0, 0, & nrm );

  // Compute svd of A and T.
  compute_svd( A, d1 );
  compute_svd( T, d2 );
  //// FLA_Obj_show( " d1 = [ ", d1, "%le", " ];" );
  //// FLA_Obj_show( " d2 = [ ", d2, "%le", " ];" );

  // Compute residuals.
  FLA_Norm_frob( d1, nrm );
  FLA_Obj_extract_real_scalar( nrm, & dnorm_vs );

  FLA_Axpy( FLA_MINUS_ONE, d1, d2 );

  FLA_Norm_frob( d2, nrm );
  FLA_Obj_extract_real_scalar( nrm, & dnorm_dif );

  if( dnorm_vs == 0.0 ) {
    *( ( double * ) FLA_Obj_buffer_at_view( resid ) ) = dnorm_dif;
  } else {
    *( ( double * ) FLA_Obj_buffer_at_view( resid ) ) = dnorm_dif / dnorm_vs;
  }

  // Remove auxiliary objects.
  FLA_Obj_free( & d1 );
  FLA_Obj_free( & d2 );
  FLA_Obj_free( & nrm );

  return FLA_SUCCESS;
}

// ============================================================================
static int compute_svd( FLA_Obj A, FLA_Obj sv ) {
// Compute singular values of A.
  FLA_Obj  Acopy, Work;
  double   * buff_Acopy, * buff_Work, * buff_U, * buff_VT, * buff_sv, dwork;
  int      info, m_Acopy, n_Acopy, ldim_Acopy, ldim_U, ldim_VT, lwork;
  char     jobu, jobvt;

  //// printf( "compute_svd\n" );

  // Make a temporal copy of A into Acopy.
  FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, A, & Acopy );
  FLA_Copy( A, Acopy );                                  
  //// FLA_Obj_show( " A = [ ", A, "%le", " ] ; " );
  //// FLA_Obj_show( " Acopy = [ ", Acopy, "%le", " ] ; " ); 

  // Some initializations.
  buff_Acopy = ( double * ) FLA_Obj_buffer_at_view( Acopy );
  ldim_Acopy = FLA_Obj_col_stride( Acopy );
  m_Acopy    = FLA_Obj_length( Acopy );
  n_Acopy    = FLA_Obj_width( Acopy );
  buff_sv    = FLA_Obj_buffer_at_view( sv );

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
  lwork = ( int ) dwork;
  //// printf( "  Optimal lwork: %d\n", lwork );

  // Create workspace.
  FLA_Obj_create( FLA_Obj_datatype( A ), lwork, 1, 0, 0, & Work );
  buff_Work = ( double * ) FLA_Obj_buffer_at_view( Work );

  // Call to SUBROUTINE DGESVD( JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT,
  //                            WORK, LWORK, INFO )
  dgesvd_( & jobu, & jobvt, & m_Acopy, & n_Acopy, buff_Acopy, & ldim_Acopy,
           buff_sv, buff_U, & ldim_U, buff_VT, & ldim_VT,
           buff_Work, & lwork, & info );
  if( info != 0 ) {
    fprintf( stderr, " *** Info after dgesvd_: %d \n", info );
  }

  // Remove object Work.
  FLA_Obj_free( & Work );

  // Remove object Acopy.
  FLA_Obj_free( & Acopy );

  return FLA_SUCCESS;
}


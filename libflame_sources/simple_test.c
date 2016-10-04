#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "FLAME.h"

#define PRINT_DATA


// ============================================================================
// Declaration of local prototypes.

static void matrix_generate( FLA_Obj A );


// ============================================================================
int main( int argc, char *argv[] ) {
  int      m_A, n_A;
  FLA_Obj  A, U, V;

  // Initialize FLAME.
  FLA_Init();

  // Some initializations.
  m_A = 6;
  n_A = 6;

  // Create FLAME objects, and attach buffers.
  FLA_Obj_create( FLA_DOUBLE, m_A, n_A, 0, 0, & A );
  FLA_Obj_create( FLA_DOUBLE, m_A, m_A, 0, 0, & U );
  FLA_Obj_create( FLA_DOUBLE, n_A, n_A, 0, 0, & V );

  // Generate matrix.
  //// FLA_Random_matrix( A );
  matrix_generate( A );

  // Print initial data.
#ifdef PRINT_DATA
  FLA_Obj_show( " Ai = [ ", A, "%le", " ];" );
#endif

  // Factorize matrix.
  printf( "%% Just before computing factorization.\n" );
  // New factorization.
  // We use a small block size to factorize the small input matrix, but you
  // should use larger blocksizes such as 64 for larger matrices.
  //// FLA_UTV_UT_blk_var1( A, 1, U, 1, V, 64, 10, 2 );
  FLA_UTV_UT_blk_var1( A, 1, U, 1, V, 3, 10, 2 );
  printf( "%% Just after computing factorization.\n" );

  // Print results.
#ifdef PRINT_DATA
  FLA_Obj_show( " Af = [ ", A, "%le", " ];" );
  FLA_Obj_show( " Uf = [ ", U, "%le", " ];" );
  FLA_Obj_show( " Vf = [ ", V, "%le", " ];" );
#endif

  // Free objects.
  FLA_Obj_free( & A );
  FLA_Obj_free( & U );
  FLA_Obj_free( & V );

  // Finalize FLAME.
  printf( "%% End of Program\n" );
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
#if 0
#endif
  }
}


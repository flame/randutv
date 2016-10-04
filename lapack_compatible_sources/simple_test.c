#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#define max( a, b ) ( (a) > (b) ? (a) : (b) )
#define min( a, b ) ( (a) < (b) ? (a) : (b) )

#define PRINT_DATA


// ============================================================================
// Declaration of local prototypes.

static void matrix_generate( int m_A, int n_A, double * buff_A, int ldim_A );

static void print_double_matrix( char * name, int m_A, int n_A, 
                double * buff_A, int ldim_A );


// ============================================================================
int main( int argc, char *argv[] ) {
  int     nb_alg, pp, m_A, n_A, mn_A, ldim_A, ldim_U, ldim_V;
  double  * buff_A, * buff_U, * buff_V;

  // Create matrix A, matrix U, and matrix V.
  m_A      = 6;
  n_A      = 6;
  mn_A     = min( m_A, n_A );

  buff_A   = ( double * ) malloc( m_A * n_A * sizeof( double ) );
  ldim_A   = max( 1, m_A );

  buff_U   = ( double * ) malloc( m_A * m_A * sizeof( double ) );
  ldim_U   = max( 1, m_A );

  buff_V   = ( double * ) malloc( n_A * n_A * sizeof( double ) );
  ldim_V   = max( 1, n_A );

  // Generate matrix.
  matrix_generate( m_A, n_A, buff_A, ldim_A );

#ifdef PRINT_DATA
  print_double_matrix( "ai", m_A, n_A, buff_A, ldim_A );
#endif

  // Factorize matrix.
  printf( "%% Just before computing factorization.\n" );
  // New factorization.
  // We use a small block size to factorize the small input matrix, but you
  // should use larger blocksizes such as 64 for larger matrices.
  NoFLA_UTV_WY_blk_var2( m_A, n_A, buff_A, ldim_A, 
      1, m_A, m_A, buff_U, ldim_U, 
      1, n_A, n_A, buff_V, ldim_V, 
      3, 10, 2 );
      //// 64, 10, 2 );
  printf( "%% Just after computing factorization.\n" );

  // Print results.
#ifdef PRINT_DATA
  print_double_matrix( "af", m_A, n_A, buff_A, ldim_A );
  print_double_matrix( "uf", m_A, m_A, buff_U, ldim_U );
  print_double_matrix( "vf", n_A, n_A, buff_V, ldim_V );
#endif

  // Free matrices and vectors.
  free( buff_A );
  free( buff_U );
  free( buff_V );

  printf( "%% End of Program\n" );

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


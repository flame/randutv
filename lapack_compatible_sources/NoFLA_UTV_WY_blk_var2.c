/*
===============================================================================
Authors
===============================================================================

Per-Gunnar Martinsson
  Dept. of Applied Mathematics, 
  University of Colorado at Boulder, 
  526 UCB, Boulder, CO 80309-0526, USA

Gregorio Quintana-Orti
  Depto. de Ingenieria y Ciencia de Computadores, 
  Universitat Jaume I, 
  12.071 Castellon, Spain

Nathan Heavner
  Dept. of Applied Mathematics, 
  University of Colorado at Boulder, 
  526 UCB, Boulder, CO 80309-0526, USA

===============================================================================
Copyright
===============================================================================

Copyright (C) 2016, 
  Universitat Jaume I,
  University of Colorado at Boulder.

===============================================================================
Disclaimer
===============================================================================

This code is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY EXPRESSED OR IMPLIED.

*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "NoFLA_UTV_WY_blk_var2.h"


// ============================================================================
// Definition of macros.

#define max( a, b )  ( (a) > (b) ? (a) : (b) )
#define min( a, b )  ( (a) > (b) ? (b) : (a) )
#define dabs( a )    ( (a) >= 0.0 ? (a) : -(a) )


// ============================================================================
// Compilation declarations.

#undef PROFILE


// ============================================================================
// Declaration of local prototypes.

static int NoFLA_Set_to_identity( int m_A, int n_A, double * buff_A,
               int ldim_A );

static int NoFLA_Zero_strict_lower_triangular( int m_A, int n_A,
               double * buff_A, int ldim_A );

static int NoFLA_Zero( int m_A, int n_A, double * buff_A, int ldim_A );

static int NoFLA_Copy_vector_into_diagonal( double * v, int m_A, int n_A,
               double * buff_A, int ldim_A );

static int NoFLA_Multiply_BAB(
               char transa, char transb,
               int m_A, int n_A, double * buff_A, int ldim_A,
               int m_B, int n_B, double * buff_B, int ldim_B );

static int NoFLA_Multiply_BBA(
               char transa, char transb,
               int m_A, int n_A, double * buff_A, int ldim_A,
               int m_B, int n_B, double * buff_B, int ldim_B );

static int NoFLA_Compute_svd(
               int m_A, int n_A, double * buff_A, int ldim_A,
               int m_U, int n_U, double * buff_U, int ldim_U,
               int n_sv, double * buff_sv,
               int m_V, int n_V, double * buff_V, int ldim_V,
               int nb_alg );

static int NoFLA_Normal_random_matrix( int m_A, int n_A,
               double * buff_A, int ldim_A );

static double NoFLA_Normal_random_number( double mu, double sigma );

static int NoFLA_Apply_Q_WY_lhfc_blk_var2(
               int m_U, int n_U, double * buff_U, int ldim_U,
               int m_T, int n_T, double * buff_T, int ldim_T,
               int m_B, int n_B, double * buff_B, int ldim_B );

static int NoFLA_Apply_Q_WY_rnfc_blk_var2(
               int m_U, int n_U, double * buff_U, int ldim_U,
               int m_T, int n_T, double * buff_T, int ldim_T,
               int m_B, int n_B, double * buff_B, int ldim_B );

static int NoFLA_QRP_WY_unb_var2( int pivoting, int num_stages,
               int m_A, int n_A, double * buff_A, int ldim_A,
               double * buff_T, int ldim_T );

static int NoFLA_QRP_compute_norms(
               int m_A, int n_A, double * buff_A, int ldim_A,
               double * buff_d, double * buff_e );

static int NoFLA_QRP_downdate_partial_norms( int m_A, int n_A,
               double * buff_d,  int st_d,
               double * buff_e,  int st_e,
               double * buff_wt, int st_wt,
               double * buff_A,  int ldim_A );

static int NoFLA_QRP_pivot_G( int j_max_col,
               int m_G, double * buff_G, int ldim_G,
               double * buff_d, double * buff_e );


// ============================================================================
int NoFLA_UTV_WY_blk_var2(
        int m_A, int n_A, double * buff_A, int ldim_A,
        int build_u, int m_U, int n_U, double * buff_U, int ldim_U,
        int build_v, int m_V, int n_V, double * buff_V, int ldim_V,
        int nb_alg, int pp, int n_iter ) {
//
// randUTV: It computes the UTV factorization of matrix A.
//
// Main features:
//   * BLAS-3 based.
//   * Compact WY transformations are used instead of UT transformations.
//   * No use of libflame.
//
// Matrices A, U, and V must be stored in column-order.
//
// Arguments:
// ----------
// m_A:      Number of rows of matrix A.
// n_A:      Number of columns of matrix A.
// buff_A:   Address of data in matrix A. Matrix to be factorized.
// ldim_A:   Leading dimension of matrix A.
// build_u:  If build_u==1, matrix U is built.
// m_U:      Number of rows of matrix U.
// n_U:      Number of columns of matrix U.
// buff_U:   Address of data in matrix U.
// ldim_U:   Leading dimension of matrix U.
// build_v:  If build_v==1, matrix V is built.
// m_V:      Number of rows of matrix V.
// n_V:      Number of columns of matrix V.
// buff_V:   Address of data in matrix V.
// ldim_V:   Leading dimension of matrix V.
// nb_alg:   Block size. Usual values for nb_alg are 32, 64, etc.
// pp:       Oversampling size. Usual values for pp are 5, 10, etc.
// n_iter:   Number of "power" iterations. Usual values are 2.
//
// Final comments:
// ---------------
// This code has been created from a libflame code. Hence, you can find some
// commented calls to libflame routines. We have left them to make it easier
// to interpret the meaning of the C code.
//
  // Declaration of variables.
  double  d_one = 1.0, d_zero = 0.0;
  char    all = 'A', t = 'T', n = 'N';
  double  * buff_G, * buff_Y, * buff_S1, * buff_S2,
          * buff_SU, * buff_sv, * buff_SVT,
          * buff_SUtl, * buff_svl, * buff_SVTtl,
          * buff_A11, * buff_ABR, * buff_AB1, * buff_AB2,
          * buff_A01, * buff_A12,
          * buff_GBl, * buff_YBl, * buff_S1tl, * buff_S2tl,
          * buff_BR, * buff_C1, * buff_D1, * buff_CR, * buff_DR;
  int     i, j, bRow, mn_A;
  int     ldim_Y, ldim_G, ldim_S1, ldim_S2, ldim_SU, ldim_SVT, ldim_SA;
  int     m_YBl, n_YBl, m_GBl, n_GBl, m_CR, n_CR, m_AB1, n_AB1, m_AB2, n_AB2,
          m_WR, n_WR, m_XR, n_XR;
#ifdef PROFILE
  double  t1, t2, tt_by,
          tt_qr1_fact, tt_qr1_updt_a, tt_qr1_updt_v,
          tt_qr2_fact, tt_qr2_updt_a, tt_qr2_updt_u,
          tt_svd_fact, tt_svd_updt_a, tt_svd_updt_uv;
#endif

  // Executable Statements.
  //// printf( "%% NoFLA_UTV_WY_blk_var2.\n" );

  // Set seed for random generator.
  srand( 12 );

  // Check matrix dimensions.
  if( m_U != n_U ) {
    fprintf( stderr, "NoFLA_UTV_WY_blk_var2: Matrix U should be square.\n" ); 
    exit( -1 );
  }
  if( m_V != n_V ) {
    fprintf( stderr, "NoFLA_UTV_WY_blk_var2: Matrix V should be square.\n" ); 
    exit( -1 );
  }
  if( m_U != m_A ) {
    fprintf( stderr, "NoFLA_UTV_WY_blk_var2: Dims. of U and A do not match.\n");
    exit( -1 );
  }
  if( n_A != m_V ) {
    fprintf( stderr, "NoFLA_UTV_WY_blk_var2: Dims. of A and V do not match.\n");
    exit( -1 );
  }

#ifdef PROFILE
  tt_by          = 0.0;
  tt_qr1_fact    = 0.0;
  tt_qr1_updt_a  = 0.0;
  tt_qr1_updt_v  = 0.0;
  tt_qr2_fact    = 0.0;
  tt_qr2_updt_a  = 0.0;
  tt_qr2_updt_u  = 0.0;
  tt_svd_fact    = 0.0;
  tt_svd_updt_a  = 0.0;
  tt_svd_updt_uv = 0.0;
#endif

  // Create and initialize auxiliary objects.
  buff_G   = ( double * ) malloc( m_A * ( nb_alg + pp ) * sizeof( double ) );
  buff_Y   = ( double * ) malloc( n_A * ( nb_alg + pp ) * sizeof( double ) );
  buff_S1  = ( double * ) malloc( nb_alg * nb_alg * sizeof( double ) );
  buff_S2  = ( double * ) malloc( nb_alg * nb_alg * sizeof( double ) );
  buff_SU  = ( double * ) malloc( nb_alg * nb_alg * sizeof( double ) );
  buff_sv  = ( double * ) malloc( nb_alg * sizeof( double ) );
  buff_SVT = ( double * ) malloc( nb_alg * nb_alg * sizeof( double ) );

  // Some initializations.
  ldim_G     = m_A;
  ldim_Y     = n_A;
  ldim_S1    = nb_alg;
  ldim_S2    = nb_alg;
  ldim_SU    = nb_alg;
  ldim_SVT   = nb_alg;
  ldim_SA    = nb_alg;
  buff_SUtl  = & buff_SU [ 0 + 0 * ldim_SU ];
  buff_svl   = & buff_sv[ 0 ];
  buff_SVTtl = & buff_SVT[ 0 + 0 * ldim_SVT ];

  // %%% Initialize U and V and copy A onto T.
  // U = eye(m);
  // V = eye(n);
  // T = A;
  if( build_u == 1 ) {
    // MyFLA_Set_to_identity( U );
    NoFLA_Set_to_identity( m_U, n_U, buff_U, ldim_U );
  }
  if( build_v == 1 ) {
    // MyFLA_Set_to_identity( V );
    NoFLA_Set_to_identity( m_V, n_V, buff_V, ldim_V );
  }

  // Main Loop.
  mn_A = min( m_A, n_A );
  for( i = 0; i < mn_A; i += nb_alg ) {
    bRow = min( nb_alg, mn_A - i );

    // Some initializations for every iteration.
    m_YBl = n_A - i;
    n_YBl = bRow + pp;
    m_GBl = m_A - i;
    n_GBl = bRow + pp;
    m_CR  = m_A;
    n_CR  = n_A - i;
    m_AB1 = m_A - i;
    n_AB1 = bRow;
    m_AB2 = m_A - i;
    n_AB2 = n_A - i - bRow;
    m_WR  = m_U;
    n_WR  = n_U - i;
    m_XR  = m_V;
    n_XR  = n_V - i;

    buff_A11  = & buff_A[ i + i * ldim_A ];
    buff_ABR  = & buff_A[ i + i * ldim_A ];
    buff_AB1  = & buff_A[ i + i * ldim_A ];
    buff_AB2  = & buff_A[ i + min( i+bRow, n_A-1 ) * ldim_A ];
    buff_GBl  = & buff_G[ i + 0 * ldim_G ];
    buff_YBl  = & buff_Y[ i + 0 * ldim_Y ];
    buff_S1tl = & buff_S1[ 0 + 0 * ldim_S1 ];
    buff_S2tl = & buff_S2[ 0 + 0 * ldim_S2 ];
    buff_BR   = & buff_A[ 0 + i * ldim_A ];
    buff_A01  = & buff_A[ 0 + i * ldim_A ];
    buff_A12  = & buff_A[ i + min( i+bRow, n_A-1 ) * ldim_A ];
    buff_C1   = & buff_U[ 0 + i * ldim_U ];
    buff_D1   = & buff_V[ 0 + i * ldim_V ];
    buff_CR   = & buff_U[ 0 + i * ldim_U ];
    buff_DR   = & buff_V[ 0 + i * ldim_V ];

    // %%% Compute the "sampling" matrix Y.
    // Aloc = T([J2,I3],[J2,J3]);
    // Y    = Aloc'*randn(m-(j-1)*b,b+p);
#ifdef PROFILE
    t1 = FLA_Clock();
#endif
    // FLA_Part_1x2( GB,  & GBl, & None1,    bRow + pp, FLA_LEFT );
    // FLA_Part_1x2( YB,  & YBl, & None1,    bRow + pp, FLA_LEFT );
    // MyFLA_Normal_random_matrix( GBl );
    NoFLA_Normal_random_matrix( m_GBl, n_GBl, buff_GBl, ldim_G );
    // FLA_Gemm( FLA_TRANSPOSE, FLA_NO_TRANSPOSE,
    //           FLA_ONE, ABR, GBl, FLA_ZERO, YBl );
    dgemm_( & t, & n, & m_YBl, & n_YBl, & m_GBl,
            & d_one,  buff_ABR, & ldim_A,
                      buff_GBl, & ldim_G,
            & d_zero, buff_YBl, & ldim_Y );

    // %%% Perform "power iteration" if requested.
    // for i_iter = 1:n_iter
    //   Y = Aloc'*(Aloc*Y);
    // end
    for( j = 0; j < n_iter; j++ ) {
      // Reuse GBl.
      // FLA_Gemm( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE,
      //           FLA_ONE, ABR, YBl, FLA_ZERO, GBl );
      dgemm_( & n, & n, & m_GBl, & n_GBl, & m_YBl,
              & d_one,  buff_ABR, & ldim_A,
                        buff_YBl, & ldim_Y,
              & d_zero, buff_GBl, & ldim_G );
      // FLA_Gemm( FLA_TRANSPOSE, FLA_NO_TRANSPOSE,
      //           FLA_ONE, ABR, GBl, FLA_ZERO, YBl );
      dgemm_( & t, & n, & m_YBl, & n_YBl, & m_GBl,
              & d_one,  buff_ABR, & ldim_A,
                        buff_GBl, & ldim_G,
              & d_zero, buff_YBl, & ldim_Y );
    }
#ifdef PROFILE
    t2 = FLA_Clock();
    tt_by += ( t2 - t1 );
#endif

    // %%% Construct the local transform to be applied "from the left".
    // if (p > 0)
    //   [~,~,Jtmp] = qr(Y,0);
    //   [Vloc,~,~] = qr(Y(:,Jtmp(1:b)));
    // else
    //   [Vloc,~]   = LOCAL_nonpiv_QR(Y,b);
    // end
#ifdef PROFILE
    t1 = FLA_Clock();
#endif
    // FLA_Part_2x2( S1,  & S1tl,  & None1,
    //                    & None2, & None3,   bRow, bRow, FLA_TL );
    if( pp > 0 ) {
      // MyFLA_QRP_UT_unb( 1, bRow, YBl, S1tl );
      NoFLA_QRP_WY_unb_var2( 1, bRow, m_YBl, n_YBl, buff_YBl, ldim_Y,
          buff_S1tl, ldim_S1 );
    } else {
      // MyFLA_QRP_UT_unb( 0, bRow, YBl, S1tl );
      NoFLA_QRP_WY_unb_var2( 0, bRow, m_YBl, n_YBl, buff_YBl, ldim_Y,
          buff_S1tl, ldim_S1 );
    }
#ifdef PROFILE
    t2 = FLA_Clock();
    tt_qr1_fact += ( t2 - t1 );
#endif

    // %%% Apply the pivot matrix to rotate maximal mass into the "J2" column.
    // T(:,[J2,J3])  = T(:,[J2,J3])*Vloc;
#ifdef PROFILE
    t1 = FLA_Clock();
#endif
    // Update matrix A with transformations from the first QR.
    // FLA_Part_2x2( YBl,   & YBl1, & None1,
    //                      & YBl2, & None2,   bRow, bRow, FLA_TL );
    // MyFLA_Apply_Q_UT_rnfc_blk( YBl1, YBl2, S1tl, B1, B2 );
    NoFLA_Apply_Q_WY_rnfc_blk_var2(
        m_YBl, bRow, buff_YBl,  ldim_Y,
        bRow,  bRow, buff_S1tl, ldim_S1,
        m_CR,  n_CR, buff_BR,   ldim_A );

#ifdef PROFILE
    t2 = FLA_Clock();
    tt_qr1_updt_a += ( t2 - t1 );
#endif

    // Update matrix V with transformations from the first QR.
#ifdef PROFILE
    t1 = FLA_Clock();
#endif
    if( build_v == 1 ) {
      // MyFLA_Apply_Q_UT_rnfc_blk( YBl1, YBl2, S1tl, D1, D2 );
      NoFLA_Apply_Q_WY_rnfc_blk_var2(
          m_YBl, bRow, buff_YBl,  ldim_Y,
          bRow,  bRow, buff_S1tl, ldim_S1,
          m_XR,  n_XR, buff_DR,   ldim_V );
    }
#ifdef PROFILE
    t2 = FLA_Clock();
    tt_qr1_updt_v += ( t2 - t1 );
#endif

    // %%% Next determine the rotations to be applied "from the left".
    // [Uloc,Dloc]      = LOCAL_nonpiv_QR(T([J2,I3],J2));
    // Factorize [ A11; A21 ].
#ifdef PROFILE
    t1 = FLA_Clock();
#endif
    // FLA_Merge_2x1( A11,
    //                A21,   & AB1 );
    // FLA_Part_2x2( S2,  & S2tl,  & None1,
    //                    & None2, & None3,   bRow, bRow, FLA_TL );
    // MyFLA_QRP_UT_unb( 0, bRow, AB1, S2tl );
    NoFLA_QRP_WY_unb_var2( 0, bRow, m_AB1, n_AB1, buff_AB1, ldim_A,
        buff_S2tl, ldim_S2 );
#ifdef PROFILE
    t2 = FLA_Clock();
    tt_qr2_fact += ( t2 - t1 );
#endif

    // Update rest of matrix A with transformations from the second QR.
#ifdef PROFILE
    t1 = FLA_Clock();
#endif
    // MyFLA_Apply_Q_UT_lhfc_blk( A11, A21, S2tl, A12, A22 );
    NoFLA_Apply_Q_WY_lhfc_blk_var2(
        m_AB1, n_AB1, buff_AB1,  ldim_A,
        bRow,  bRow,  buff_S2tl, ldim_S2,
        m_AB2, n_AB2, buff_AB2,  ldim_A );
#ifdef PROFILE
    t2 = FLA_Clock();
    tt_qr2_updt_a += ( t2 - t1 );
#endif

    // Update matrix U with transformations from the second QR.
#ifdef PROFILE
    t1 = FLA_Clock();
#endif
    if( build_u == 1 ) {
      // MyFLA_Apply_Q_UT_rnfc_blk( A11, A21, S2tl, C1, C2 );
      NoFLA_Apply_Q_WY_rnfc_blk_var2(
          m_AB1, n_AB1, buff_AB1,  ldim_A,
          bRow,  bRow,  buff_S2tl, ldim_S2,
          m_WR,  n_WR,  buff_CR,   ldim_U );
    }
#ifdef PROFILE
    t2 = FLA_Clock();
    tt_qr2_updt_u += ( t2 - t1 );
#endif

    // Compute miniSVD.
    // [Utmp,Dtmp,Wloc] = svd(Dloc(1:b,:));
    // Dloc(1:b,:)      = Dtmp;
    // Uloc(:,1:b)      = Uloc(:,1:b)*Utmp;
    // Vloc(:,1:b)      = Vloc(:,1:b)*Wloc; % Update Vloc.
    // T([J2,I3],J2)    = Dloc;
    // T(J1,J2)         = T(J1,J2)*Wloc;
    // T([J2,I3],J3)    = Uloc'*T([J2,I3],J3);
    //
    // %%% Store away the ON matrices.
    // U(:,[J2,I3]) = U(:,[J2,I3])*Uloc;
    // V(:,[J2,J3]) = V(:,[J2,J3])*Vloc;
#ifdef PROFILE
    t1 = FLA_Clock();
#endif
    // FLA_Part_2x2( SU,  & SUtl,  & None1,
    //                    & None2, & None3,   bRow, bRow, FLA_TL );
    // FLA_Part_1x2( sv,  & svl,   & None1,   bRow, FLA_LEFT );
    // FLA_Part_2x2( SVT, & SVTtl, & None1,
    //                    & None2, & None3,   bRow, bRow, FLA_TL );
    // MyFLA_Zero_strict_lower_triangular( AB1 );
    NoFLA_Zero_strict_lower_triangular( m_AB1, n_AB1, buff_AB1, ldim_A );
    // MyFLA_Compute_svd( A11, SUtl, svl, SVTtl, bRow );
    NoFLA_Compute_svd(
        bRow, bRow, buff_A11, ldim_A,
        bRow, bRow, buff_SUtl, ldim_SU,
        bRow, buff_svl,
        bRow, bRow, buff_SVTtl, ldim_SVT,
        bRow );
    // MyFLA_Zero( A11 );
    NoFLA_Zero( bRow, bRow, buff_A11, ldim_A );
    // MyFLA_Copy_vector_into_diagonal( svl, A11 );
    NoFLA_Copy_vector_into_diagonal( buff_sv, bRow, bRow, buff_A11, ldim_A );
#ifdef PROFILE
    t2 = FLA_Clock();
    tt_svd_fact += ( t2 - t1 );
#endif

#ifdef PROFILE
    t1 = FLA_Clock();
#endif
    // Apply U of miniSVD to A.
    // FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, A12, & A12copy );
    // FLA_Copy( A12, A12copy );
    // FLA_Gemm( FLA_TRANSPOSE, FLA_NO_TRANSPOSE,
    //           FLA_ONE, SUtl, A12copy, FLA_ZERO, A12 );
    // FLA_Obj_free( & A12copy );
    NoFLA_Multiply_BAB( 't', 'n',
        bRow, bRow,           buff_SUtl, ldim_SU,
        bRow, n_A - i - bRow, buff_A12,  ldim_A );

    // Apply V of miniSVD to A.
    // FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, A01, & A01copy );
    // FLA_Copy( A01, A01copy );
    // FLA_Gemm( FLA_NO_TRANSPOSE, FLA_TRANSPOSE,
    //           FLA_ONE, A01copy, SVTtl, FLA_ZERO, A01 );
    // FLA_Obj_free( & A01copy );
    NoFLA_Multiply_BBA( 't', 'n',
        bRow, bRow, buff_SVTtl, ldim_SVT,
        i,    bRow, buff_A01,   ldim_A );
#ifdef PROFILE
    t2 = FLA_Clock();
    tt_svd_updt_a += ( t2 - t1 );
#endif

#ifdef PROFILE
    t1 = FLA_Clock();
#endif
    // Apply U of miniSVD to global U.
    if( build_u == 1 ) {
      // FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, C1, & C1copy );
      // FLA_Copy( C1, C1copy );
      // FLA_Gemm( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE,
      //           FLA_ONE, C1copy, SUtl, FLA_ZERO, C1 );
      // FLA_Obj_free( & C1copy );
      NoFLA_Multiply_BBA( 'n', 'n',
          bRow, bRow, buff_SUtl, ldim_SU,
          m_U,  bRow, buff_C1,   ldim_U );
    }

    // Apply V of miniSVD to global V.
    if( build_v == 1 ) {
      // FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, D1, & D1copy );
      // FLA_Copy( D1, D1copy );
      // FLA_Gemm( FLA_NO_TRANSPOSE, FLA_TRANSPOSE,
      //           FLA_ONE, D1copy, SVTtl, FLA_ZERO, D1 );
      // FLA_Obj_free( & D1copy );
      NoFLA_Multiply_BBA( 't', 'n',
          bRow, bRow, buff_SVTtl, ldim_SVT,
          m_V,  bRow, buff_D1,    ldim_V );
    }
#ifdef PROFILE
    t2 = FLA_Clock();
    tt_svd_updt_uv += ( t2 - t1 );
#endif

    // End of main loop.
  }

  // %%% Process the last block.
  // %%% At this point, no "block pivoting" is required.
  // J1               = 1:(n-b);
  // J2               = (b*(nstep-1)+1):n;
  // I2               = (b*(nstep-1)+1):m;
  // [Uloc,Dloc,Wloc] = svd(T(I2,J2));
  // T(I2,J2)         = Dloc;
  // T(J1,J2)         = T(J1,J2)*Wloc;
  // U( :,I2)         = U(:,I2)*Uloc;
  // V( :,J2)         = V(:,J2)*Wloc;
  //
  // return
  // The last block is processed inside the previous loop.

  // Remove auxiliary objects.
  free( buff_Y );
  free( buff_G );
  free( buff_S1 );
  free( buff_S2 );
  free( buff_SU );
  free( buff_sv );
  free( buff_SVT );

#ifdef PROFILE
  printf( "%% tt_build_y:     %le\n", tt_by );
  printf( "%% tt_qr1:         %le\n", tt_qr1_fact + tt_qr1_updt_a +
                                      tt_qr1_updt_v );
  printf( "%%     tt_qr1_fact:    %le\n", tt_qr1_fact );
  printf( "%%     tt_qr1_updt_a:  %le\n", tt_qr1_updt_a );
  printf( "%%     tt_qr1_updt_v:  %le\n", tt_qr1_updt_v );
  printf( "%% tt_qr2:         %le\n", tt_qr2_fact + tt_qr2_updt_a +
                                      tt_qr2_updt_u );
  printf( "%%     tt_qr2_fact:    %le\n", tt_qr2_fact );
  printf( "%%     tt_qr2_updt_a:  %le\n", tt_qr2_updt_a );
  printf( "%%     tt_qr2_updt_u:  %le\n", tt_qr2_updt_u );
  printf( "%% tt_svd:         %le\n", tt_svd_fact + tt_svd_updt_a +
                                      tt_svd_updt_uv);
  printf( "%%     tt_svd_fact:    %le\n", tt_svd_fact );
  printf( "%%     tt_svd_updt_a:  %le\n", tt_svd_updt_a );
  printf( "%%     tt_svd_updt_uv: %le\n", tt_svd_updt_uv );
  printf( "%% total_time:     %le\n",
          tt_by +
          tt_qr1_fact + tt_qr1_updt_a + tt_qr1_updt_v +
          tt_qr2_fact + tt_qr2_updt_a + tt_qr2_updt_u +
          tt_svd_fact + tt_svd_updt_a + tt_svd_updt_uv );
#endif

  return 0;
}

// ============================================================================
static int NoFLA_Set_to_identity( int m_A, int n_A, double * buff_A,
               int ldim_A ) {
// Set contents of object A to the identity matrix.
  int  i, j, mn_A;

  // Set the full matrix.
  for ( j = 0; j < n_A; j++ ) {
    for ( i = 0; i < m_A; i++ ) {
      buff_A[ i + j * ldim_A ] = 0.0;
    }
  }
  // Set the main diagonal.
  mn_A = min( m_A, n_A );
  for ( j = 0; j < mn_A; j++ ) {
    buff_A[ j + j * ldim_A ] = 1.0;
  }

  return 0;
}

// ============================================================================
static int NoFLA_Zero_strict_lower_triangular( int m_A, int n_A,
               double * buff_A, int ldim_A ) {
// Zero the strictly lower triangular part of matrix A.
  int  i, j;

  // Set the strictly lower triangular matrix.
  for ( j = 0; j < n_A; j++ ) {
    for ( i = j + 1; i < m_A; i++ ) {
      buff_A[ i + j * ldim_A ] = 0.0;
    }
  }

  return 0;
}

// ============================================================================
static int NoFLA_Zero( int m_A, int n_A, double * buff_A, int ldim_A ) {
// Set the contents of matrix A to zero.
  int  i, j;

  // Set the full matrix.
  for ( j = 0; j < n_A; j++ ) {
    for ( i = 0; i < m_A; i++ ) {
      buff_A[ i + j * ldim_A ] = 0.0;
    }
  }

  return 0;
}

// ============================================================================
static int NoFLA_Copy_vector_into_diagonal( double * v, int m_A, int n_A,
               double * buff_A, int ldim_A ) {
// Copy the contents of vector v into the diagonal of the matrix A.
  int  i, j, mn_A;

  // Copy vector into the diagonal.
  mn_A = min( m_A, n_A );
  for ( i = 0; i < mn_A; i++ ) {
    buff_A[ i + i * ldim_A ] = v[ i ];
  }

  return 0;
}

// ============================================================================
static int NoFLA_Multiply_BAB(
               char transa, char transb,
               int m_A, int n_A, double * buff_A, int ldim_A,
               int m_B, int n_B, double * buff_B, int ldim_B ) {
// Compute B := A * B.

  char    all = 'A';
  double  d_one = 1.0, d_zero = 0.0;
  double  * buff_Bcopy;
  int     ldim_Bcopy;

  if( ( m_B > 0 )&&( n_B > 0 ) ) {
    //// FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, B, & Bcopy );
    buff_Bcopy = ( double * ) malloc( m_B * n_B * sizeof( double ) );
    ldim_Bcopy = m_B;

    //// FLA_Copy( B, Bcopy );
    dlacpy_( & all, & m_B, & n_B, buff_B, & ldim_B,
                                  buff_Bcopy, & ldim_Bcopy );

    //// FLA_Gemm( transa, transb, FLA_ONE, A, Bcopy, FLA_ZERO, B );
    dgemm_( & transa, & transb, & m_B, & n_B, & m_B,
            & d_one,  buff_A, & ldim_A,
                      buff_Bcopy, & ldim_Bcopy,
            & d_zero, buff_B, & ldim_B );

    //// FLA_Obj_free( & Bcopy );
    free( buff_Bcopy );
  }

  return 0;
}

// ============================================================================
static int NoFLA_Multiply_BBA(
               char transa, char transb,
               int m_A, int n_A, double * buff_A, int ldim_A,
               int m_B, int n_B, double * buff_B, int ldim_B ) {
// Compute B := B * A.

  char    all = 'A';
  double  d_one = 1.0, d_zero = 0.0;
  double  * buff_Bcopy;
  int     ldim_Bcopy;

  if( ( m_B > 0 )&&( n_B > 0 ) ) {
    //// FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, B, & Bcopy );
    buff_Bcopy = ( double * ) malloc( m_B * n_B * sizeof( double ) );
    ldim_Bcopy = m_B;

    //// FLA_Copy( B, Bcopy );
    dlacpy_( & all, & m_B, & n_B, buff_B, & ldim_B,
                                  buff_Bcopy, & ldim_Bcopy );

    //// FLA_Gemm( transb, transa, FLA_ONE, Bcopy, A, FLA_ZERO, B );
    dgemm_( & transb, & transa, & m_B, & n_B, & n_B,
            & d_one,  buff_Bcopy, & ldim_Bcopy,
                      buff_A, & ldim_A,
            & d_zero, buff_B, & ldim_B );

    //// FLA_Obj_free( & Bcopy );
    free( buff_Bcopy );
  }

  return 0;
}

// ============================================================================
static int NoFLA_Compute_svd(
               int m_A, int n_A, double * buff_A, int ldim_A,
               int m_U, int n_U, double * buff_U, int ldim_U,
               int n_sv, double * buff_sv,
               int m_V, int n_V, double * buff_V, int ldim_V,
               int nb_alg ) {
// Compute:  U, and V of svd of A.
  char    all = 'A';
  double  * buff_Workspace;
  int     info, max_mn_A, min_mn_A, lwork;

  // Some initializations.
  max_mn_A = max( m_A, n_A );
  min_mn_A = min( m_A, n_A );

  // Create Workspace.
  // According to lapack's documentation,
  // workspace for dgebd2 should be: max( m, n ), and
  // workspace for dbdsqr should be: 2*n.
  // However, dgebd2 seems to need more.  So, workspace is increased.
  lwork  = max( 1,
                max( 3 * min_mn_A + max_mn_A, 5 * min_mn_A ) )
           + nb_alg * max_mn_A + 100000 + 10 * m_A + 10 * n_A;
  //// FLA_Obj_create( FLA_Obj_datatype( A ), lwork, 1, 0, 0, & Workspace );
  buff_Workspace = ( double * ) malloc( lwork * sizeof( double ) );
  //// printf( " lwork: %d\n ", lwork );

  // Call to SUBROUTINE DGESVD( JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT,
  //                            WORK, LWORK, INFO )
  dgesvd_( & all, & all, & m_A, & n_A,
           buff_A, & ldim_A, buff_sv,
           buff_U, & ldim_U, buff_V, & ldim_V,
           buff_Workspace, & lwork, & info );
  if( info != 0 ) {
    fprintf( stderr, " *** Info after dgesvd_f: %d \n", info );
  }

  // Remove object Work.
  //// FLA_Obj_free( & Workspace );
  free( buff_Workspace );

  return 0;
}

// ============================================================================
static int NoFLA_Normal_random_matrix( int m_A, int n_A,
               double * buff_A, int ldim_A ) {
//
// It generates a random matrix with normal distribution.
//
  int  i, j;

  // Main loop.
  for ( j = 0; j < n_A; j++ ) {
    for ( i = 0; i < m_A; i++ ) {
      buff_A[ i + j * ldim_A ] = NoFLA_Normal_random_number( 0.0, 1.0 );
    }
  }

  return 0;
}

// ============================================================================
static double NoFLA_Normal_random_number( double mu, double sigma ) {
  static int     alternate_calls = 0;
  static double  b1, b2;
  double         c1, c2, a, factor;

  // Quick return.
  if( alternate_calls == 1 ) {
    alternate_calls = ! alternate_calls;
    return( mu + sigma * b2 );
  }
  // Main loop.
  do {
    c1 = -1.0 + 2.0 * ( (double) rand() / RAND_MAX );
    c2 = -1.0 + 2.0 * ( (double) rand() / RAND_MAX );
    a = c1 * c1 + c2 * c2;
  } while ( ( a == 0 )||( a >= 1 ) );
  factor = sqrt( ( -2 * log( a ) ) / a );
  b1 = c1 * factor;
  b2 = c2 * factor;
  alternate_calls = ! alternate_calls;
  return( mu + sigma * b1 );
}

// ============================================================================
static int NoFLA_Apply_Q_WY_lhfc_blk_var2(
               int m_U, int n_U, double * buff_U, int ldim_U,
               int m_T, int n_T, double * buff_T, int ldim_T,
               int m_B, int n_B, double * buff_B, int ldim_B ) {
//
// It applies the transpose of a block transformation Q to a matrix B from
// the left:
//   B := Q' * B
// where:
//   Q = I - U * T' * U'.
//
  double  * buff_W;
  int     ldim_W;

  // Create auxiliary object.
  //// FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, B1, & W );
  buff_W = ( double * ) malloc( n_B * n_U * sizeof( double ) );
  ldim_W = max( 1, n_B );

  // Apply the block transformation.
  dlarfb_( "Left", "Transpose", "Forward", "Columnwise",
           & m_B, & n_B, & n_U, buff_U, & ldim_U, buff_T, & ldim_T,
           buff_B, & ldim_B, buff_W, & ldim_W );

  // Remove auxiliary object.
  //// FLA_Obj_free( & W );
  free( buff_W );

  return 0;
}

// ============================================================================
static int NoFLA_Apply_Q_WY_rnfc_blk_var2(
               int m_U, int n_U, double * buff_U, int ldim_U,
               int m_T, int n_T, double * buff_T, int ldim_T,
               int m_B, int n_B, double * buff_B, int ldim_B ) {
//
// It applies a block transformation Q to a matrix B from the right:
//   B = B * Q
// where:
//   Q = I - U * T' * U'.
//
  double  * buff_W;
  int     ldim_W;

  // Create auxiliary object.
  //// FLA_Obj_create_conf_to( FLA_TRANSPOSE, B1, & W );
  buff_W = ( double * ) malloc( m_B * n_U * sizeof( double ) );
  ldim_W = max( 1, m_B );

  // Apply the block transformation.
  dlarfb_( "Right", "No transpose", "Forward", "Columnwise",
           & m_B, & n_B, & n_U, buff_U, & ldim_U, buff_T, & ldim_T,
           buff_B, & ldim_B, buff_W, & ldim_W );

  // Remove auxiliary object.
  //// FLA_Obj_free( & W );
  free( buff_W );

  return 0;
}

// ============================================================================
static int NoFLA_QRP_WY_unb_var2( int pivoting, int num_stages,
               int m_A, int n_A, double * buff_A, int ldim_A,
               double * buff_T, int ldim_T ) {
//
// It computes an unblocked QR factorization of matrix A with or without
// pivoting. Matrices B and C are optionally pivoted, and matrix T is
// optionally built.
//
// Arguments:
// "pivoting": If pivoting==1, then QR factorization with pivoting is used.
// "numstages": It tells the number of columns that are factorized.
//   If "num_stages" is negative, the whole matrix A is factorized.
//   If "num_stages" is positive, only the first "num_stages" are factorized.
//
  int     j, mn_A, m_a21, m_A22, n_A22, n_dB, idx_max_col,
          i_one = 1, n_house_vector, m_rest;
  double  * buff_d, * buff_e, * buff_workspace, * buff_t, diag;
  int     idamax_();

  //// printf( "NoFLA_QRP_WY_unb_var2. pivoting: %d \n", pivoting );

  // Some initializations.
  mn_A    = min( m_A, n_A );

  // Set the number of stages, if needed.
  if( num_stages < 0 ) {
    num_stages = mn_A;
  }

  // Create auxiliary vectors.
  buff_d         = ( double * ) malloc( n_A * sizeof( double ) );
  buff_e         = ( double * ) malloc( n_A * sizeof( double ) );
  buff_workspace = ( double * ) malloc( n_A * sizeof( double ) );
  buff_t         = ( double * ) malloc( n_A * sizeof( double ) );

  if( pivoting == 1 ) {
    // Compute initial norms of A into d and e.
    NoFLA_QRP_compute_norms( m_A, n_A, buff_A, ldim_A, buff_d, buff_e );
  }

  // Main Loop.
  for( j = 0; j < num_stages; j++ ) {
    n_dB  = n_A - j;
    m_a21 = m_A - j - 1;
    m_A22 = m_A - j - 1;
    n_A22 = n_A - j - 1;

    if( pivoting == 1 ) {
      // Obtain the index of the column with largest 2-norm.
      idx_max_col = idamax_( & n_dB, & buff_d[ j ], & i_one ) - 1;

      // Swap columns of A, B, C, pivots, and norms vectors.
      NoFLA_QRP_pivot_G( idx_max_col,
          m_A, & buff_A[ 0 + j * ldim_A ], ldim_A,
          & buff_d[ j ],
          & buff_e[ j ] );
    }

    // Compute tau1 and u21 from alpha11 and a21 such that tau1 and u21
    // determine a Householder transform H such that applying H from the
    // left to the column vector consisting of alpha11 and a21 annihilates
    // the entries in a21 (and updates alpha11).
    n_house_vector = m_a21 + 1;
    dlarfg_( & n_house_vector,
             & buff_A[ j + j * ldim_A ],
             & buff_A[ min( m_A-1, j+1 ) + j * ldim_A ], & i_one,
             & buff_t[ j ] );

    // / a12t \ =  H / a12t \
    // \ A22  /      \ A22  /
    //
    // where H is formed from tau1 and u21.
    diag = buff_A[ j + j * ldim_A ];
    buff_A[ j + j * ldim_A ] = 1.0;
    m_rest = m_A22 + 1;
    dlarf_( "Left", & m_rest, & n_A22,
        & buff_A[ j + j * ldim_A ], & i_one,
        & buff_t[ j ],
        & buff_A[ j + ( j+1 ) * ldim_A ], & ldim_A,
        buff_workspace );
    buff_A[ j + j * ldim_A ] = diag;

    if( pivoting == 1 ) {
      // Update partial column norms.
      NoFLA_QRP_downdate_partial_norms( m_A22, n_A22,
          & buff_d[ j+1 ], 1,
          & buff_e[ j+1 ], 1,
          & buff_A[ j + ( j+1 ) * ldim_A ], ldim_A,
          & buff_A[ ( j+1 ) + min( n_A-1, ( j+1 ) ) * ldim_A ], ldim_A );
    }
  }

  // Build T.
  dlarft_( "Forward", "Columnwise", & m_A, & num_stages, buff_A, & ldim_A,
           buff_t, buff_T, & ldim_T );

  // Remove auxiliary vectors.
  free( buff_d );
  free( buff_e );
  free( buff_workspace );
  free( buff_t );

  return 0;
}

// ============================================================================
static int NoFLA_QRP_compute_norms(
               int m_A, int n_A, double * buff_A, int ldim_A,
               double * buff_d, double * buff_e ) {
//
// It computes the column norms of matrix A. The norms are stored into
// vectors d and e.
//
  int     j, i_one = 1;
  double  dnrm2_();

  // Main loop.
  for( j = 0; j < n_A; j++ ) {
    * buff_d = dnrm2_( & m_A, buff_A, & i_one );
    * buff_e = * buff_d;
    buff_A += ldim_A;
    buff_d++;
    buff_e++;
  }

  return 0;
}

// ============================================================================
static int NoFLA_QRP_downdate_partial_norms( int m_A, int n_A,
               double * buff_d,  int st_d,
               double * buff_e,  int st_e,
               double * buff_wt, int st_wt,
               double * buff_A,  int ldim_A ) {
//
// It updates (downdates) the column norms of matrix A. It uses Drmac's method.
//
  int     j, i_one = 1;
  double  * ptr_d, * ptr_e, * ptr_wt, * ptr_A;
  double  temp, temp2, temp5, tol3z;
  double  dnrm2_(), dlamch_();

  /*
*
*           Update partial column norms
*
          DO 30 J = I + 1, N
             IF( WORK( J ).NE.ZERO ) THEN
*
*                 NOTE: The following 4 lines follow from the analysis in
*                 Lapack Working Note 176.
*
                TEMP = ABS( A( I, J ) ) / WORK( J )
                TEMP = MAX( ZERO, ( ONE+TEMP )*( ONE-TEMP ) )
                TEMP2 = TEMP*( WORK( J ) / WORK( N+J ) )**2
                IF( TEMP2 .LE. TOL3Z ) THEN
                   IF( M-I.GT.0 ) THEN
                      WORK( J ) = DNRM2( M-I, A( I+1, J ), 1 )
                      WORK( N+J ) = WORK( J )
                   ELSE
                      WORK( J ) = ZERO
                      WORK( N+J ) = ZERO
                   END IF
                ELSE
                   WORK( J ) = WORK( J )*SQRT( TEMP )
                END IF
             END IF
 30       CONTINUE
  */

  // Some initializations.
  tol3z = sqrt( dlamch_( "Epsilon" ) );
  ptr_d  = buff_d;
  ptr_e  = buff_e;
  ptr_wt = buff_wt;
  ptr_A  = buff_A;

  // Main loop.
  for( j = 0; j < n_A; j++ ) {
    if( * ptr_d != 0.0 ) {
      temp = dabs( * ptr_wt ) / * ptr_d;
      temp = max( 0.0, ( 1.0 + temp ) * ( 1 - temp ) );
      temp5 = * ptr_d / * ptr_e;
      temp2 = temp * temp5 * temp5;
      if( temp2 <= tol3z ) {
        if( m_A > 0 ) {
          * ptr_d = dnrm2_( & m_A, ptr_A, & i_one );
          * ptr_e = *ptr_d;
        } else {
          * ptr_d = 0.0;
          * ptr_e = 0.0;
        }
      } else {
        * ptr_d = * ptr_d * sqrt( temp );
      }
    }
    ptr_A  += ldim_A;
    ptr_d  += st_d;
    ptr_e  += st_e;
    ptr_wt += st_wt;
  }

  return 0;
}


// ============================================================================
static int NoFLA_QRP_pivot_G( int j_max_col,
               int m_G, double * buff_G, int ldim_G,
               double * buff_d, double * buff_e ) {
//
// It pivots matrix G, pivot vector p, and norms vectors d and e.
// Matrices B and C are optionally pivoted.
//
  int     i_one = 1;
  double  * ptr_g1, * ptr_g2; //// , * ptr_b1, * ptr_b2, * ptr_c1, * ptr_c2;

  // Swap columns of G, pivots, and norms.
  if( j_max_col != 0 ) {

    // Swap full column 0 and column "j_max_col" of G.
    ptr_g1 = & buff_G[ 0 + 0         * ldim_G ];
    ptr_g2 = & buff_G[ 0 + j_max_col * ldim_G ];
    dswap_( & m_G, ptr_g1, & i_one, ptr_g2, & i_one );

    // Copy norms of column 0 to column "j_max_col".
    buff_d[ j_max_col ] = buff_d[ 0 ];
    buff_e[ j_max_col ] = buff_e[ 0 ];
  }

  return 0;
}



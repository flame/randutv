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
#include "FLAME.h"
#include "FLA_UTV_UT_blk_var1.h"


#undef PROFILE


// ============================================================================
// Declaration of local prototypes.

static FLA_Error MyFLA_Set_to_identity( FLA_Obj A );

static FLA_Error MyFLA_Zero_strict_lower_triangular( FLA_Obj A );

static FLA_Error MyFLA_Zero( FLA_Obj A );

static FLA_Error MyFLA_Copy_vector_into_diagonal( FLA_Obj v, FLA_Obj A );

static FLA_Error MyFLA_Normal_random_matrix( FLA_Obj A );

static double MyFLA_Normal_random_number( double mu, double sigma );

static FLA_Error MyFLA_Apply_Q_UT_rnfc_blk( FLA_Obj U11, FLA_Obj U21,
                     FLA_Obj T,
                     FLA_Obj B1, FLA_Obj B2 );

static FLA_Error MyFLA_Apply_Q_UT_lhfc_blk( FLA_Obj U11, FLA_Obj U21,
                     FLA_Obj T, 
                     FLA_Obj B1, FLA_Obj B2 );

static FLA_Error MyFLA_Compute_svd( FLA_Obj A, FLA_Obj U, FLA_Obj sv, 
                     FLA_Obj V, int nb_alg );

static FLA_Error MyFLA_QRP_UT_unb( int use_pivoting, int num_stages,
                     FLA_Obj A, FLA_Obj T );

static FLA_Error MyFLA_QRP_compute_norms( FLA_Obj A, FLA_Obj d, FLA_Obj e );

static FLA_Error MyFLA_Apply_H2_UT_l_opd( 
                     int m_u2_A2,
                     int n_a1t,
                     double * tau,
                     double * u2, int inc_u2,
                     double * a1t, int inc_a1t,
                     double * A2, int ldim_A2,
                     double * workspace );

static FLA_Error MyFLA_QRP_pivot_A( int j_max_col,
                     int m_A, double * buff_A, int ldim_A, 
                     double * buff_d, double * buff_e );

static FLA_Error MyFLA_QRP_update_partial_norms_opt( int m_A, int n_A,
                     double * buff_d,  int st_d,
                     double * buff_e,  int st_e,
                     double * buff_wt, int st_wt,
                     double * buff_A,  int ldim_A );


// ============================================================================
FLA_Error FLA_UTV_UT_blk_var1( FLA_Obj A, int build_u, FLA_Obj U, 
              int build_v, FLA_Obj V, int nb_alg, int pp, int n_iter ) {
//
// randUTV: It computes the UTV factorization of matrix A.
//
// Main features:
//   * BLAS-3 based.
//   * Use of Householder UT block transformations.
//   * Use of libflame.
//
// Arguments:
// ----------
// A:       (input)  Matrix to be factorized.
//          (output) Matrix factorized.
// build_u: (input)  If build_u==1, matrix U is built.
// U:       (output) Matrix U of the factorization.
// build_v: (input)  If build_v==1, matrix V is built.
// V:       (output) Matrix V of the factorization.
// nb_alg:  (input)  Block size. Usual values for nb_alg are 32, 64, etc.
// pp:      (input)  Oversampling size. Usual values for pp are 5, 10, etc.
// n_iter:  (input)  Number of "power" iterations. Usual values are 2.
//
  // Declaration of variables.
  FLA_Obj ATL, ATR,    A00, A01, A02,
          ABL, ABR,    A10, A11, A12,
                       A20, A21, A22;
  FLA_Obj UTL, UTR,    U00, U01, U02,
          UBL, UBR,    U10, U11, U12,
                       U20, U21, U22;
  FLA_Obj VTL, VTR,    V00, V01, V02,
          VBL, VBR,    V10, V11, V12,
                       V20, V21, V22;
  FLA_Obj YT,          Y0,
          YB,          Y1,
                       Y2;
  FLA_Obj GT,          G0,
          GB,          G1,
                       G2;
  FLA_Obj BL, BR,      B0,  B1,  B2;
  FLA_Obj CL, CR,      C0,  C1,  C2;
  FLA_Obj DL, DR,      D0,  D1,  D2;
  FLA_Obj AB1, Y, G, YBl, YBl1, YBl2, GBl, None1, None2, None3,
          S1, S1tl, S2, S2tl,
          SU, sv, SVT, SUtl, svl, SVTtl,
          C1copy, D1copy,
          A12copy, A01copy;
  int     bRow, m_A, n_A, dtype_A, j;
#ifdef PROFILE
  double  t1, t2, tt_by, 
          tt_qr1_fact, tt_qr1_updt_a, tt_qr1_updt_v, 
          tt_qr2_fact, tt_qr2_updt_a, tt_qr2_updt_u, 
          tt_svd_fact, tt_svd_updt_a, tt_svd_updt_uv;
#endif

  // Executable Statements.
  //// printf( "%% FLA_UTV_UT_blk_var1.\n" );

  // Set seed for random generator.
  srand( 12 );

  // Check matrix dimensions.
  if( FLA_Obj_length( U ) != FLA_Obj_width( U ) ) {
    FLA_Print_message( "FLA_UTV_UT_blk_var1: Matrix U should be square.", 
                       __FILE__, __LINE__ );
    FLA_Abort();
  }
  if( FLA_Obj_length( V ) != FLA_Obj_width( V ) ) {
    FLA_Print_message( "FLA_UTV_UT_blk_var1: Matrix V should be square.", 
                       __FILE__, __LINE__ );
    FLA_Abort();
  }
  if( FLA_Obj_length( U ) != FLA_Obj_length( A ) ) {
    FLA_Print_message( "FLA_UTV_UT_blk_var1: Dims. of U and A do not match.",
                       __FILE__, __LINE__ );
    FLA_Abort();
  }
  if( FLA_Obj_width( A ) != FLA_Obj_length( V ) ) {
    FLA_Print_message( "FLA_UTV_UT_blk_var1: Dims. of A and V do not match.",
                       __FILE__, __LINE__ );
    FLA_Abort();
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

  // Some initializations.
  m_A     = FLA_Obj_length( A );
  n_A     = FLA_Obj_width ( A );
  dtype_A = FLA_Obj_datatype( A );

  // %%% Initialize U and V and copy A onto T.
  // U = eye(m);
  // V = eye(n);
  // T = A;
  if( build_u == 1 ) {
    MyFLA_Set_to_identity( U );
  }
  if( build_v == 1 ) {
    MyFLA_Set_to_identity( V );
  }

  // Create and initialize auxiliary objects.
  FLA_Obj_create( dtype_A, n_A,    nb_alg + pp, 0, 0, & Y );
  FLA_Obj_create( dtype_A, m_A,    nb_alg + pp, 0, 0, & G );
  FLA_Obj_create( dtype_A, nb_alg, nb_alg,      0, 0, & S1 );
  FLA_Obj_create( dtype_A, nb_alg, nb_alg,      0, 0, & S2 );
  FLA_Obj_create( dtype_A, nb_alg, nb_alg,      0, 0, & SU );
  FLA_Obj_create( dtype_A, nb_alg, 1,           0, 0, & sv );
  FLA_Obj_create( dtype_A, nb_alg, nb_alg,      0, 0, & SVT );

  // Initial Partitioning.
  FLA_Part_2x2( A,    & ATL, & ATR,
                      & ABL, & ABR,    0, 0, FLA_TL );
  FLA_Part_2x1( Y,    & YT,
                      & YB,            0, FLA_TOP );
  FLA_Part_2x1( G,    & GT,
                      & GB,            0, FLA_TOP );
  FLA_Part_1x2( A,    & BL,  & BR,     0, FLA_LEFT );
  if( build_u == 1 ) {
    FLA_Part_2x2( U,    & UTL, & UTR,
                        & UBL, & UBR,    0, 0, FLA_TL );
    FLA_Part_1x2( U,    & CL,  & CR,     0, FLA_LEFT );
  }
  if( build_v == 1 ) {
    FLA_Part_2x2( V,    & VTL, & VTR,
                        & VBL, & VBR,    0, 0, FLA_TL );
    FLA_Part_1x2( V,    & DL,  & DR,     0, FLA_LEFT );
  }

  // Main Loop.
  while( ( FLA_Obj_length( ATL ) < FLA_Obj_length( A ) )&&
         ( FLA_Obj_width ( ATL ) < FLA_Obj_width ( A ) ) ) {
    bRow = min( FLA_Obj_min_dim( ABR ), nb_alg );

    // Iteration Initial Partitioning.
    FLA_Repart_2x2_to_3x3( ATL, /**/ ATR,       &A00, /**/ &A01, &A02,
                        /* ************* */   /* ******************** */
                                                &A10, /**/ &A11, &A12,
                           ABL, /**/ ABR,       &A20, /**/ &A21, &A22,
                           bRow, bRow, FLA_BR );
    FLA_Repart_2x1_to_3x1( YT,                  &Y0,
                        /* ** */              /* ** */
                                                &Y1,
                           YB,                  &Y2,        bRow, FLA_BOTTOM );
    FLA_Repart_2x1_to_3x1( GT,                  &G0,
                        /* ** */              /* ** */
                                                &G1,
                           GB,                  &G2,        bRow, FLA_BOTTOM );
    FLA_Repart_1x2_to_1x3( BL,  /**/ BR,        &B0, /**/ &B1, &B2,
                           bRow, FLA_RIGHT );
    if( build_u == 1 ) {
      FLA_Repart_2x2_to_3x3( UTL, /**/ UTR,       &U00, /**/ &U01, &U02,
                          /* ************* */   /* ******************** */
                                                  &U10, /**/ &U11, &U12,
                             UBL, /**/ UBR,       &U20, /**/ &U21, &U22,
                             bRow, bRow, FLA_BR );
      FLA_Repart_1x2_to_1x3( CL,  /**/ CR,        &C0, /**/ &C1, &C2,
                             bRow, FLA_RIGHT );
    }
    if( build_v == 1 ) {
      FLA_Repart_2x2_to_3x3( VTL, /**/ VTR,       &V00, /**/ &V01, &V02,
                          /* ************* */   /* ******************** */
                                                  &V10, /**/ &V11, &V12,
                             VBL, /**/ VBR,       &V20, /**/ &V21, &V22,
                             bRow, bRow, FLA_BR );
      FLA_Repart_1x2_to_1x3( DL,  /**/ DR,        &D0, /**/ &D1, &D2,
                             bRow, FLA_RIGHT );
    }

    // ------------------------------------------------------------------------

    // %%% Compute the "sampling" matrix Y.
    // Aloc = T([J2,I3],[J2,J3]);
    // Y    = Aloc'*randn(m-(j-1)*b,b+p);
#ifdef PROFILE
    t1 = FLA_Clock();
#endif
    FLA_Part_1x2( GB,  & GBl, & None1,    bRow + pp, FLA_LEFT );
    FLA_Part_1x2( YB,  & YBl, & None1,    bRow + pp, FLA_LEFT );
    MyFLA_Normal_random_matrix( GBl );
    FLA_Gemm( FLA_TRANSPOSE, FLA_NO_TRANSPOSE,
              FLA_ONE, ABR, GBl, FLA_ZERO, YBl );

    // %%% Perform "power iteration" if requested.
    // for i_iter = 1:n_iter
    //   Y = Aloc'*(Aloc*Y);
    // end
    for( j = 0; j < n_iter; j++ ) {
      // Reuse GBl.
      FLA_Gemm( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE,
                FLA_ONE, ABR, YBl, FLA_ZERO, GBl );
      FLA_Gemm( FLA_TRANSPOSE, FLA_NO_TRANSPOSE,
                FLA_ONE, ABR, GBl, FLA_ZERO, YBl );
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
    FLA_Part_2x2( S1,  & S1tl,  & None1,
                       & None2, & None3,   bRow, bRow, FLA_TL );
    if( pp > 0 ) {
      MyFLA_QRP_UT_unb( 1, bRow, YBl, S1tl );
    } else {
      MyFLA_QRP_UT_unb( 0, bRow, YBl, S1tl );
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
    FLA_Part_2x2( YBl,   & YBl1, & None1,
                         & YBl2, & None2,   bRow, bRow, FLA_TL );

    // Update matrix A with transformations from the first QR.
    MyFLA_Apply_Q_UT_rnfc_blk( YBl1, YBl2, S1tl, B1, B2 );
#ifdef PROFILE
    t2 = FLA_Clock();
    tt_qr1_updt_a += ( t2 - t1 );
#endif

    // Update matrix V with transformations from the first QR.
#ifdef PROFILE
    t1 = FLA_Clock();
#endif
    if( build_v == 1 ) {
      MyFLA_Apply_Q_UT_rnfc_blk( YBl1, YBl2, S1tl, D1, D2 );
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
    FLA_Merge_2x1( A11,
                   A21,   & AB1 );
    FLA_Part_2x2( S2,  & S2tl,  & None1,
                       & None2, & None3,   bRow, bRow, FLA_TL );
    MyFLA_QRP_UT_unb( 0, bRow, AB1, S2tl );
#ifdef PROFILE
    t2 = FLA_Clock();
    tt_qr2_fact += ( t2 - t1 );
#endif

    // Update rest of matrix A with transformations from the second QR.
#ifdef PROFILE
    t1 = FLA_Clock();
#endif
    MyFLA_Apply_Q_UT_lhfc_blk( A11, A21, S2tl, A12, A22 );
#ifdef PROFILE
    t2 = FLA_Clock();
    tt_qr2_updt_a += ( t2 - t1 );
#endif

    // Update matrix U with transformations from the second QR.
#ifdef PROFILE
    t1 = FLA_Clock();
#endif
    if( build_u == 1 ) {
      MyFLA_Apply_Q_UT_rnfc_blk( A11, A21, S2tl, C1, C2 );
    }
#ifdef PROFILE
    t2 = FLA_Clock();
    tt_qr2_updt_u += ( t2 - t1 );
#endif

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
    FLA_Part_2x2( SU,  & SUtl,  & None1,
                       & None2, & None3,   bRow, bRow, FLA_TL );
    FLA_Part_1x2( sv,  & svl,   & None1,   bRow, FLA_LEFT );
    FLA_Part_2x2( SVT, & SVTtl, & None1,
                       & None2, & None3,   bRow, bRow, FLA_TL );

    // Compute miniSVD.
    MyFLA_Zero_strict_lower_triangular( AB1 );
    MyFLA_Compute_svd( A11, SUtl, svl, SVTtl, bRow );
    MyFLA_Zero( A11 );
    MyFLA_Copy_vector_into_diagonal( svl, A11 );
#ifdef PROFILE
    t2 = FLA_Clock();
    tt_svd_fact += ( t2 - t1 );
#endif

#ifdef PROFILE
    t1 = FLA_Clock();
#endif
    // Apply U of miniSVD to A.
    FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, A12, & A12copy );
    FLA_Copy( A12, A12copy );
    FLA_Gemm( FLA_TRANSPOSE, FLA_NO_TRANSPOSE,
              FLA_ONE, SUtl, A12copy, FLA_ZERO, A12 );
    FLA_Obj_free( & A12copy );

    // Apply V of miniSVD to A.
    FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, A01, & A01copy );
    FLA_Copy( A01, A01copy );
    FLA_Gemm( FLA_NO_TRANSPOSE, FLA_TRANSPOSE,
              FLA_ONE, A01copy, SVTtl, FLA_ZERO, A01 );
    FLA_Obj_free( & A01copy );
#ifdef PROFILE
    t2 = FLA_Clock();
    tt_svd_updt_a += ( t2 - t1 );
#endif

#ifdef PROFILE
    t1 = FLA_Clock();
#endif
    // Apply U of miniSVD to global U.
    if( build_u == 1 ) {
      FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, C1, & C1copy );
      FLA_Copy( C1, C1copy );
      FLA_Gemm( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE,
                FLA_ONE, C1copy, SUtl, FLA_ZERO, C1 );
      FLA_Obj_free( & C1copy );
    }

    // Apply V of miniSVD to global V.
    if( build_v == 1 ) {
      FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, D1, & D1copy );
      FLA_Copy( D1, D1copy );
      FLA_Gemm( FLA_NO_TRANSPOSE, FLA_TRANSPOSE,
                FLA_ONE, D1copy, SVTtl, FLA_ZERO, D1 );
      FLA_Obj_free( & D1copy );
    }
#ifdef PROFILE
    t2 = FLA_Clock();
    tt_svd_updt_uv += ( t2 - t1 );
#endif

    // ------------------------------------------------------------------------
    // Iteration Final Repartitioning.

    FLA_Cont_with_3x3_to_2x2( &ATL, /**/ &ATR,       A00, A01, /**/ A02,
                                                     A10, A11, /**/ A12,
                            /* ************** */  /* ****************** */
                              &ABL, /**/ &ABR,       A20, A21, /**/ A22,
                              FLA_TL );
    FLA_Cont_with_3x1_to_2x1( &YT,                   Y0,
                                                     Y1,
                            /* ** */              /* ** */
                              &YB,                   Y2,     FLA_TOP );
    FLA_Cont_with_3x1_to_2x1( &GT,                   G0,
                                                     G1,
                            /* ** */              /* ** */
                              &GB,                   G2,     FLA_TOP );
    FLA_Cont_with_1x3_to_1x2( &BL,  /**/ &BR,        B0, B1, /**/ B2,
                              FLA_LEFT );
    if( build_u == 1 ) {
      FLA_Cont_with_3x3_to_2x2( &UTL, /**/ &UTR,       U00, U01, /**/ U02,
                                                       U10, U11, /**/ U12,
                              /* ************** */  /* ****************** */
                                &UBL, /**/ &UBR,       U20, U21, /**/ U22,
                                FLA_TL );
      FLA_Cont_with_1x3_to_1x2( &CL,  /**/ &CR,        C0, C1, /**/ C2,
                                FLA_LEFT );
    }
    if( build_v == 1 ) {
      FLA_Cont_with_3x3_to_2x2( &VTL, /**/ &VTR,       V00, V01, /**/ V02,
                                                       V10, V11, /**/ V12,
                              /* ************** */  /* ****************** */
                                &VBL, /**/ &VBR,       V20, V21, /**/ V22,
                                FLA_TL );
      FLA_Cont_with_1x3_to_1x2( &DL,  /**/ &DR,        D0, D1, /**/ D2,
                                FLA_LEFT );
    }
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

  // Remove auxiliary objects.
  FLA_Obj_free( & Y );
  FLA_Obj_free( & G );
  FLA_Obj_free( & S1 );
  FLA_Obj_free( & S2 );
  FLA_Obj_free( & SU );
  FLA_Obj_free( & sv );
  FLA_Obj_free( & SVT );

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

  return FLA_SUCCESS;
}

// ============================================================================
static FLA_Error MyFLA_Set_to_identity( FLA_Obj A ) {
// Set contents of object A to the identity matrix.

  switch( FLA_Obj_datatype( A ) ) {

  case FLA_DOUBLE:
    {
      int     m_A, n_A, rs_A, cs_A, i, j, mn_A;
      double  * buff_A;

      buff_A  = ( double * ) FLA_Obj_buffer_at_view( A );
      rs_A    = FLA_Obj_row_stride( A );
      cs_A    = FLA_Obj_col_stride( A );
      m_A     = FLA_Obj_length( A );
      n_A     = FLA_Obj_width( A );
      mn_A    = min( m_A, n_A );

      // Set the full matrix.
      for ( j = 0; j < n_A; j++ ) {
        for ( i = 0; i < m_A; i++ ) {
          buff_A[ i * rs_A + j * cs_A ] = 0.0;
        }
      }
      // Set the main diagonal.
      for ( j = 0; j < mn_A; j++ ) {
        buff_A[ j * rs_A + j * cs_A ] = 1.0;
      }
    }
    break;

  default:
    fprintf( stderr, "+++ ERROR in MyFLA_Set_to_identity:  " );
    fprintf( stderr, "Datatype not implemented: %d\n", FLA_Obj_datatype( A ) );
  }

  return FLA_SUCCESS;
}

// ============================================================================
static FLA_Error MyFLA_Zero_strict_lower_triangular( FLA_Obj A ) {
// Set the strictly lower triangular part of matrix A to zero.

  switch( FLA_Obj_datatype( A ) ) {

  case FLA_DOUBLE:
    {
      int     m_A, n_A, rs_A, cs_A, i, j;
      double  * buff_A;

      buff_A  = ( double * ) FLA_Obj_buffer_at_view( A );
      rs_A    = FLA_Obj_row_stride( A );
      cs_A    = FLA_Obj_col_stride( A );
      m_A     = FLA_Obj_length( A );
      n_A     = FLA_Obj_width( A );

      // Set the strictly lower triangular part matrix.
      for ( j = 0; j < n_A; j++ ) {
        for ( i = j + 1; i < m_A; i++ ) {
          buff_A[ i * rs_A + j * cs_A ] = 0.0;
        }
      }
    }
    break;

  default:
    fprintf( stderr, "+++ ERROR in MyFLA_Zero_strict_lower_triangular:  " );
    fprintf( stderr, "Datatype not implemented: %d\n", FLA_Obj_datatype( A ) );
  }

  return FLA_SUCCESS;
}

// ============================================================================
static FLA_Error MyFLA_Zero( FLA_Obj A ) {
// Set the contents of matrix A to zero.

  switch( FLA_Obj_datatype( A ) ) {

  case FLA_DOUBLE:
    {
      int     m_A, n_A, rs_A, cs_A, i, j;
      double  * buff_A;

      buff_A  = ( double * ) FLA_Obj_buffer_at_view( A );
      rs_A    = FLA_Obj_row_stride( A );
      cs_A    = FLA_Obj_col_stride( A );
      m_A     = FLA_Obj_length( A );
      n_A     = FLA_Obj_width( A );

      // Set the strictly lower triangular part matrix.
      for ( j = 0; j < n_A; j++ ) {
        for ( i = 0; i < m_A; i++ ) {
          buff_A[ i * rs_A + j * cs_A ] = 0.0;
        }
      }
    }
    break;

  default:
    fprintf( stderr, "+++ ERROR in MyFLA_Zero:  " );
    fprintf( stderr, "Datatype not implemented: %d\n", FLA_Obj_datatype( A ) );
  }

  return FLA_SUCCESS;
}

// ============================================================================
static FLA_Error MyFLA_Copy_vector_into_diagonal( FLA_Obj v, FLA_Obj A ) {
// Copy the contents of vector v into the diagonal of matrix A.

  switch( FLA_Obj_datatype( A ) ) {

  case FLA_DOUBLE:
    {
      int     mn_A, rs_A, cs_A, i, j;
      double  * buff_v, * buff_A;

      buff_v  = ( double * ) FLA_Obj_buffer_at_view( v );
      buff_A  = ( double * ) FLA_Obj_buffer_at_view( A );
      rs_A    = FLA_Obj_row_stride( A );
      cs_A    = FLA_Obj_col_stride( A );
      mn_A    = FLA_Obj_min_dim( A );

      // Copy the vector into the diagonal.
      for ( i = 0; i < mn_A; i++ ) {
        buff_A[ i * rs_A + i * cs_A ] = buff_v[ i ];
      }
    }
    break;

  default:
    fprintf( stderr, "+++ ERROR in MyFLA_Copy_vector_into_diagonal:  " );
    fprintf( stderr, "Datatype not implemented: %d\n", FLA_Obj_datatype( A ) );
  }

  return FLA_SUCCESS;
}

// ============================================================================
static FLA_Error MyFLA_Normal_random_matrix( FLA_Obj A ) {
// Set random matrix with normal distribution.

  switch( FLA_Obj_datatype( A ) ) {

  case FLA_DOUBLE:
    {
      double  * buff_A;
      int     m_A, n_A, ldim_A, i, j;

      // Some initializations.
      m_A     = FLA_Obj_length( A );
      n_A     = FLA_Obj_width ( A );
      buff_A  = ( double * ) FLA_Obj_buffer_at_view( A );
      ldim_A  = FLA_Obj_col_stride( A );

      // Main loop.
      for ( j = 0; j < n_A; j++ ) {
        for ( i = 0; i < m_A; i++ ) {
          buff_A[ i + j * ldim_A ] = MyFLA_Normal_random_number( 0.0, 1.0 );
        }
      }
    }
    break;

  default:
    fprintf( stderr, "+++ ERROR in MyFLA_Normal_random_matrix:  " );
    fprintf( stderr, "Datatype not implemented: %d\n", FLA_Obj_datatype( A ) );
  }

  return FLA_SUCCESS;
}

// ============================================================================
static double MyFLA_Normal_random_number( double mu, double sigma ) {
//
// It computes and returns a normal random number.
//
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
static FLA_Error MyFLA_Apply_Q_UT_rnfc_blk( FLA_Obj U11, FLA_Obj U21,
                     FLA_Obj T,
                     FLA_Obj B1, FLA_Obj B2 ) {
//
// Apply a unitary matrix Q to a matrix B from the right:
//   B = B * Q
// where:
//   B = [ B1; B2 ]
//   U = [ U11; U21 ]
//   Q = ( I - U * inv( T ) * U' ).
//
  FLA_Obj  W;

  //// printf( "MyFLA_Apply_Q_UT_rnfc_blk\n" );

  // Create auxiliary object.
  FLA_Obj_create_conf_to( FLA_TRANSPOSE, B1, & W );

  // W = B1^T;

  FLA_Copyt_external( FLA_TRANSPOSE, B1, W );

  // U11 = trilu( U11 );
  // U21 = U21;
  // Let W^T be conformal to B1.
  // W^T = ( B1 * U11 + B2 * U21 ) * inv( triu(T) );
  // W   = inv( triu(T)^T ) * ( U11^T * B1^T + U21^T * B2^T );

  FLA_Trmm_external( FLA_LEFT, FLA_LOWER_TRIANGULAR,
                     FLA_TRANSPOSE, FLA_UNIT_DIAG,
                     FLA_ONE, U11, W );

  FLA_Gemm_external( FLA_TRANSPOSE, FLA_TRANSPOSE,
                     FLA_ONE, U21, B2, FLA_ONE, W );

  FLA_Trsm_external( FLA_LEFT, FLA_UPPER_TRIANGULAR,
                     FLA_TRANSPOSE, FLA_NONUNIT_DIAG,
                     FLA_ONE, T, W );

  // B2 = B2 - W^T * U21';
  // B1 = B1 - W^T * U11';
  //    = B1 - ( conj(U11) * W )^T;

  FLA_Gemm_external( FLA_TRANSPOSE, FLA_CONJ_TRANSPOSE,
                     FLA_MINUS_ONE, W, U21, FLA_ONE, B2 );

  FLA_Trmm_external( FLA_LEFT, FLA_LOWER_TRIANGULAR,
                     FLA_CONJ_NO_TRANSPOSE, FLA_UNIT_DIAG,
                     FLA_MINUS_ONE, U11, W );

  FLA_Axpyt_external( FLA_TRANSPOSE, FLA_ONE, W, B1 );

  // Remove auxiliary object.
  FLA_Obj_free( & W );

  return FLA_SUCCESS;
}

// ============================================================================
static FLA_Error MyFLA_Apply_Q_UT_lhfc_blk( FLA_Obj U11, FLA_Obj U21,
                     FLA_Obj T, 
                     FLA_Obj B1, FLA_Obj B2 ) {
//
// Apply the conjugate-transpose of a unitary matrix Q to a matrix B from the
// left:
//   B := Q' B
// where:
//   B = [ B1; B2 ]
//   U = [ U11; U21 ]
//   Q = ( I - U inv(T) U' )'.
//
  FLA_Obj  W;

  //// printf( "MyFLA_Apply_Q_UT_lhfc_blk\n" );

  // Create auxiliary object.
  FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, B1, & W );

  // W = B1;

  FLA_Copyt_external( FLA_NO_TRANSPOSE, B1, W );

  // U11 = trilu( U11 );
  // U21 = U21;
  // W = inv( triu( T ) )' * ( U11' * B1 + U21' * B2 );

  FLA_Trmm_external( FLA_LEFT, FLA_LOWER_TRIANGULAR,
                     FLA_CONJ_TRANSPOSE, FLA_UNIT_DIAG,
                     FLA_ONE, U11, W );

  FLA_Gemm_external( FLA_CONJ_TRANSPOSE, FLA_NO_TRANSPOSE, 
                     FLA_ONE, U21, B2, FLA_ONE, W );

  FLA_Trsm_external( FLA_LEFT, FLA_UPPER_TRIANGULAR,
                     FLA_CONJ_TRANSPOSE, FLA_NONUNIT_DIAG,
                     FLA_ONE, T, W );

  // B2 = B2 - U21 * W;
  // B1 = B1 - U11 * W;

  FLA_Gemm_external( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE,
                     FLA_MINUS_ONE, U21, W, FLA_ONE, B2 );

  FLA_Trmm_external( FLA_LEFT, FLA_LOWER_TRIANGULAR,
                     FLA_NO_TRANSPOSE, FLA_UNIT_DIAG,
                     FLA_MINUS_ONE, U11, W );

  FLA_Axpyt_external( FLA_NO_TRANSPOSE, FLA_ONE, W, B1 );

  // Remove auxiliary object.
  FLA_Obj_free( & W );

  return FLA_SUCCESS;
}

// ============================================================================
static FLA_Error MyFLA_Compute_svd( FLA_Obj A, FLA_Obj U, FLA_Obj sv, 
                     FLA_Obj V, int nb_alg ) {
// Compute:  U, vs, and V of svd of A.
  FLA_Obj  Workspace;
  double   * buff_A, * buff_U, * buff_sv, * buff_V, * buff_Workspace;
  int      info, m_A, n_A, max_mn_A, min_mn_A, 
           ldim_A, ldim_U, ldim_V, lwork;
  char     all = 'A';

  // Some initializations.
  m_A      = FLA_Obj_length( A );
  n_A      = FLA_Obj_width( A );
  max_mn_A = max( m_A, n_A );
  min_mn_A = min( m_A, n_A );
  buff_A   = ( double * ) FLA_Obj_buffer_at_view( A );
  ldim_A   = FLA_Obj_col_stride( A );
  buff_U   = ( double * ) FLA_Obj_buffer_at_view( U );
  ldim_U   = FLA_Obj_col_stride( U );
  buff_sv  = ( double * ) FLA_Obj_buffer_at_view( sv );
  buff_V   = ( double * ) FLA_Obj_buffer_at_view( V );
  ldim_V   = FLA_Obj_col_stride( V );

  // Create object Workspace.
  // According to lapack's documentation,
  // workspace for dgebd2 should be: max( m, n ), and
  // workspace for dbdsqr should be: 2*n.
  // However, dgebd2 seems to need more.  So, workspace is increased.
  lwork  = max( 1, 
                max( 3 * min_mn_A + max_mn_A, 5 * min_mn_A ) ) 
           + nb_alg * max_mn_A + 100000 + 10 * m_A + 10 * n_A;
  FLA_Obj_create( FLA_Obj_datatype( A ), lwork, 1, 0, 0, & Workspace );
  buff_Workspace = ( double * ) FLA_Obj_buffer_at_view( Workspace );
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

  // Remove object Workspace.
  FLA_Obj_free( & Workspace );

  return FLA_SUCCESS;
}

// ============================================================================
static FLA_Error MyFLA_QRP_UT_unb( int use_pivoting, int num_stages,
                     FLA_Obj A, FLA_Obj T ) {
//
// Arguments:
// "use_pivoting" : 1=Matrix A is pivoted. 0=Matrix A is not pivoted.
//

  // Declaration of variables.
  FLA_Obj  t, d, e, workspace;
  int      j, m_A, n_A, mn_A, dtype_A, ldim_A,
           ldim_T, m_A20, n_A20, m_a21, m_A22, n_A22, n_dB, idx_max_col, 
           i_one = 1;
  double   * buff_A, * buff_t, * buff_T, 
           * buff_d, * buff_e, * buff_workspace, d_one = 1.0;

  int      idamax_();

  //// printf( "  %% MyFLA_QRP_UT_unb. \n" );

  // Some initializations.
  dtype_A = FLA_Obj_datatype( A );
  m_A     = FLA_Obj_length( A );
  n_A     = FLA_Obj_width ( A );
  mn_A    = min( m_A, n_A );
  ldim_A  = FLA_Obj_col_stride( A );
  buff_A  = ( double * ) FLA_Obj_buffer_at_view( A );

  buff_T  = ( double * ) FLA_Obj_buffer_at_view( T );
  ldim_T  = FLA_Obj_col_stride( T );

  FLA_Obj_create( dtype_A, n_A, 1, 0, 0, & t );
  FLA_Obj_create( dtype_A, n_A, 1, 0, 0, & d );
  FLA_Obj_create( dtype_A, n_A, 1, 0, 0, & e );
  FLA_Obj_create( dtype_A, n_A, 1, 0, 0, & workspace );

  buff_t         = ( double * ) FLA_Obj_buffer_at_view( t );
  buff_d         = ( double * ) FLA_Obj_buffer_at_view( d );
  buff_e         = ( double * ) FLA_Obj_buffer_at_view( e );
  buff_workspace = ( double * ) FLA_Obj_buffer_at_view( workspace );

  // Set the number of stages, if needed.
  if( num_stages < 0 ) {
    num_stages = mn_A;
  }

  // Compute initial norms of A into d and e.
  if( use_pivoting == 1 ) {
    MyFLA_QRP_compute_norms( A, d, e );
  }

  // Main Loop.
  for( j = 0; j < num_stages; j++ ) {
    n_dB  = n_A - j;
    m_a21 = m_A - j - 1;
    m_A22 = m_A - j - 1;
    n_A22 = n_A - j - 1;
    m_A20 = m_A - j - 1;
    n_A20 = j;

    // Obtain the index of the column with largest 2-norm.
    if( use_pivoting == 1 ) {
      //// FLA_Amax( dB, max_col );
      idx_max_col = idamax_( & n_dB, & buff_d[ j ], & i_one ) - 1;
    }

    // Swap columns of A, and norms.
    if( use_pivoting == 1 ) {
      MyFLA_QRP_pivot_A( idx_max_col,
          m_A, & buff_A[ 0 + j * ldim_A ], ldim_A,
          & buff_d[ j ],
          & buff_e[ j ] );
    }

    // Compute tau1 and u21 from alpha11 and a21 such that tau1 and u21
    // determine a Householder transform H such that applying H from the
    // left to the column vector consisting of alpha11 and a21 annihilates
    // the entries in a21 (and updates alpha11).
    //// FLA_Househ2( alpha11,
    ////              a21, tau1 );
    //// FLA_Househ2_UT( FLA_LEFT,
    ////                 alpha11,
    ////                 a21, tau1 );
    FLA_Househ2_UT_l_opd( m_a21,
                          & buff_A[ j + j * ldim_A ],
                          & buff_A[ ( j+1 ) + j * ldim_A ], i_one,
                          & buff_t[ j ] );

    // / a12t \ =  H / a12t \
    // \ A22  /      \ A22  /
    //
    // where H is formed from tau1 and u21.
    //// FLA_QR_Update_Rest_blk_b2( tau1, a21, a12t,
    ////                                       A22 );
    /// FLA_Apply_H2_UT( FLA_LEFT, tau1, a21, a12t,
    ///                                       A22 );
    MyFLA_Apply_H2_UT_l_opd(
        m_A22, //// m_u2_A2,
        n_A22, //// n_a1t,
        & buff_t[ j ], //// tau,
        & buff_A[ ( j+1 ) + j * ldim_A ], 1, //// u2,inc_u2,
        & buff_A[ j + ( j+1 ) * ldim_A ], ldim_A, //// a1t,inc_a1t,
        & buff_A[ ( j+1 ) + ( j+1 ) * ldim_A ], ldim_A,
        buff_workspace ); //// A2,rs_A2,cs_A2);

    // Build T.
    // rho11 = tau1;
    //// FLA_Copy( tau1, rho11 );
    buff_T[ j + j * ldim_T ] = buff_t[ j ];

    // t01 = a10t' + A20' * u21;
    //// FLA_Copyt_external( FLA_CONJ_TRANSPOSE, a10t, r01 );
    dcopy_( & j, & buff_A[ j + 0 * ldim_A ], & ldim_A,
                 & buff_T[ 0 + j * ldim_T ], & i_one );

    //// FLA_Gemv_external( FLA_CONJ_TRANSPOSE, FLA_ONE, A20, a21, 
    ////                    FLA_ONE, r01 );
    dgemv_( "Transpose", & m_A20, & n_A20,
        & d_one,
        & buff_A[ ( j+1 ) + 0 * ldim_A ], & ldim_A,
        & buff_A[ ( j+1 ) + j * ldim_A ], & i_one,
        & d_one,
        & buff_T[ 0 + j * ldim_T ], & i_one );

    // Update partial column norms.
    if( use_pivoting == 1 ) {
      //// MyFLA_QRP_update_partial_norms( dB, eB, a12t, A22 );
      MyFLA_QRP_update_partial_norms_opt( m_A22, n_A22, 
          & buff_d[ j+1 ], 1,
          & buff_e[ j+1 ], 1,
          & buff_A[ j + ( j+1 ) * ldim_A ], ldim_A,
          & buff_A[ ( j+1 ) + ( j+1 ) * ldim_A ], ldim_A );
    }
  }

  FLA_Obj_free( & t );
  FLA_Obj_free( & d );
  FLA_Obj_free( & e );
  FLA_Obj_free( & workspace );

  return FLA_SUCCESS;
}

// ============================================================================
static FLA_Error MyFLA_QRP_compute_norms( FLA_Obj A, FLA_Obj d, FLA_Obj e ) {

  switch( FLA_Obj_datatype( A ) ) {

  case FLA_DOUBLE: {
    int     m_A, n_A, ld_A, st_d, st_e, j, i_one = 1;
    double  * ptr_A, * ptr_d, * ptr_e, dnrm2_();

    m_A   = FLA_Obj_length( A );
    n_A   = FLA_Obj_width( A );
    ptr_A = ( double * ) FLA_Obj_buffer_at_view( A );
    ptr_d = ( double * ) FLA_Obj_buffer_at_view( d );
    ptr_e = ( double * ) FLA_Obj_buffer_at_view( e );
    ld_A  = FLA_Obj_col_stride( A );

    st_d = 0;
    if( FLA_Obj_length( d ) == 1 ) {
      st_d = FLA_Obj_col_stride( d );
    } else if( FLA_Obj_width( d ) == 1 ) {
      st_d = 1;
    } else {
      FLA_Print_message( "MyFLA_QRP_compute_norms: Object d is not a vector", 
                         __FILE__, __LINE__ );
      FLA_Abort();
    }
    st_e = 0;
    if( FLA_Obj_length( e ) == 1 ) {
      st_e = FLA_Obj_col_stride( e );
    } else if( FLA_Obj_width( e ) == 1 ) {
      st_e = 1;
    } else {
      FLA_Print_message( "MyFLA_QRP_compute_norms: Object e is not a vector", 
                         __FILE__, __LINE__ );
      FLA_Abort();
    }

    for( j = 0; j < n_A; j++ ) {
      *ptr_d = dnrm2_( & m_A, ptr_A, & i_one );
      *ptr_e = *ptr_d;
      ptr_A += ld_A;
      ptr_d += st_d;
      ptr_e += st_e;
    }

    break;
  }
  default:
    FLA_Print_message( "MyFLA_QRP_compute_norms: datatype not yet implemented", 
                       __FILE__, __LINE__ );
    FLA_Abort();
  }

  return FLA_SUCCESS;
}

// ============================================================================
static FLA_Error MyFLA_Apply_H2_UT_l_opd( 
                     int m_u2_A2,
                     int n_a1t,
                     double * tau,
                     double * u2, int inc_u2,
                     double * a1t, int inc_a1t,
                     double * A2, int ldim_A2, 
                     double * workspace ) {

  double  one_p       = 1.0;
  double  minus_one_p = -1.0;
  double  rtau;
  int     inc_w1t;

  // FLA_Obj w1t;
  double * w1t;

  // if ( FLA_Obj_has_zero_dim( a1t ) ) return FLA_SUCCESS;
  if ( n_a1t == 0 ) return FLA_SUCCESS;

  // FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, a1t, &w1t );
  //// w1t = ( double* ) FLA_malloc( n_a1t * sizeof( *a1t ) );
  w1t = workspace;
  inc_w1t = 1;

  // // w1t = a1t;
  // FLA_Copy_external( a1t, w1t );
  //// bl1_dcopyv( BLIS1_NO_CONJUGATE,
  ////             n_a1t,
  ////             a1t, inc_a1t, 
  ////             w1t, inc_w1t ); 
  dcopy_( & n_a1t,
          a1t, & inc_a1t, 
          w1t, & inc_w1t ); 

  // // w1t = w1t + u2' * A2;
  // // w1t = w1t + A2^T * conj(u2);
  // FLA_Gemvc_external( FLA_TRANSPOSE, FLA_CONJUGATE, FLA_ONE, A2, u2, 
  //     FLA_ONE, w1t );
  //// bl1_dgemv( BLIS1_TRANSPOSE,
  ////            BLIS1_CONJUGATE,
  ////            m_u2_A2,
  ////            n_a1t,
  ////            one_p,
  ////            A2, rs_A2, cs_A2,
  ////            u2, inc_u2,
  ////            one_p,
  ////            w1t, inc_w1t );
  dgemv_( "Transpose",
          & m_u2_A2,
          & n_a1t,
          & one_p,
          A2, & ldim_A2,
          u2, & inc_u2,
          & one_p,
          w1t, & inc_w1t );

  // // w1t = w1t / tau;
  // FLA_Inv_scalc_external( FLA_NO_CONJUGATE, tau, w1t );
  //// bl1_dinvscalv( BLIS1_NO_CONJUGATE,
  ////                n_a1t,
  ////                tau,
  ////                w1t, inc_w1t );
  if( * tau == 0.0 ) {
    fprintf( stderr, "ERROR in MyFLA_Apply_H2_UT_l_opd: Tau is zero.\n" );
  } else {
    rtau = 1.0 / ( * tau );
  }
  dscal_( & n_a1t,
          & rtau,
          w1t, & inc_w1t );

  // // a1t = - w1t + a1t;
  // FLA_Axpy_external( FLA_MINUS_ONE, w1t, a1t );
  //// bl1_daxpyv( BLIS1_NO_CONJUGATE,
  ////             n_a1t,
  ////             minus_one_p,
  ////             w1t, inc_w1t,
  ////             a1t, inc_a1t );
  daxpy_( & n_a1t,
          & minus_one_p,
          w1t, & inc_w1t,
          a1t, & inc_a1t );

  // // A2 = - u2 * w1t + A2;
  // FLA_Ger_external( FLA_MINUS_ONE, u2, w1t, A2 );
  //// bl1_dger( BLIS1_NO_CONJUGATE,
  ////           BLIS1_NO_CONJUGATE,
  ////           m_u2_A2,
  ////           n_a1t,
  ////           minus_one_p,
  ////           u2, inc_u2,
  ////           w1t, inc_w1t,
  ////           A2, rs_A2, cs_A2 );
  dger_( & m_u2_A2,
         & n_a1t,
         & minus_one_p,
         u2, & inc_u2,
         w1t, & inc_w1t,
         A2, & ldim_A2 );

  //// // FLA_Obj_free( &w1t );
  //// FLA_free( w1t );

  return FLA_SUCCESS;
}

// ============================================================================
static FLA_Error MyFLA_QRP_pivot_A( int j_max_col,
                     int m_A, double * buff_A, int ldim_A, 
                     double * buff_d, double * buff_e ) {

  // Declaration of variables.
  int     i_one = 1;
  double  * ptr_a1, * ptr_a2;

  //// j_max_col = *( ( int * ) FLA_Obj_buffer_at_view( max_col ) );
  //// printf( "j_max_col: %d\n", j_max_col );
  //// FLA_Obj_show( " max_col = [ ", max_col, "%d", " ];" );

  // Swap columns of G, pivots, and norms.
  if( j_max_col != 0 ) {

    // Swap full column 0 and column "j_max_col" of G.
    ptr_a1 = & buff_A[ 0 + 0         * ldim_A ];
    ptr_a2 = & buff_A[ 0 + j_max_col * ldim_A ];
    dswap_( & m_A, ptr_a1, & i_one, ptr_a2, & i_one );

    // Copy norms of column 0 to column "j_max_col".
    buff_d[ j_max_col ] = buff_d[ 0 ];
    buff_e[ j_max_col ] = buff_e[ 0 ];
  }

  return FLA_SUCCESS;
}

// ============================================================================
static FLA_Error MyFLA_QRP_update_partial_norms_opt( int m_A, int n_A,
                     double * buff_d,  int st_d,
                     double * buff_e,  int st_e,
                     double * buff_wt, int st_wt,
                     double * buff_A,  int ldim_A ) {
  //// printf( "  MyFLA_QRP_update_partial_norms_opt by Drmac. \n" );

  int     j, i_one = 1;
  double  * ptr_d, * ptr_e, * ptr_wt, * ptr_A;
  double  temp, temp2, temp5, tol3z;
  double  dnrm2_();

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
  tol3z = sqrt( FLA_Mach_params_opd( FLA_MACH_EPS ) );
  //// printf( "tol3z: %24.17le\n", tol3z );
  //// printf( "eps:   %24.17le\n", FLA_Mach_params_opd( FLA_MACH_EPS ) );
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
      //// printf( "temp2: %24.17lf   tol3z: %24.17lf\n", temp2, tol3z );
      if( temp2 <= tol3z ) {
        //// printf( "F. Cancel Drm %4d %14.6le %14.6le\n", 
        ////         j, * ptr_wt, * ptr_d );
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

  return FLA_SUCCESS;
}



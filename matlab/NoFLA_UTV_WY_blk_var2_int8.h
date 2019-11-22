// hacked version of randutv header
// ahb: hack to use pointer size (aka ptrdiff_t) MKL blas/lapack in MATLAB uses
#define int long int

int NoFLA_UTV_WY_blk_var2_int8(
        int m_A, int n_A, double * buff_A, int ldim_A,
        int build_u, int m_U, int n_U, double * buff_U, int ldim_U,
        int build_v, int m_V, int n_V, double * buff_V, int ldim_V,
        int nb_alg, int pp, int n_iter );

#undef int

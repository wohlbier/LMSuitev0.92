set (LIBSRC
     dcabs1.f
     izamax.f
     lsame.f
     xerbla.f
     zgemm.f
     zgeru.f
     zscal.f
     zswap.f
     ztrsm.f)

add_library (mini_blas_s ${LIBSRC})

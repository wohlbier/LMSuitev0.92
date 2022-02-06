set (LIBSRC
     zgesv.f
     zgetf2.f
     zgetrf.f
     zgetrs.f
     zlaswp.f)
add_library (mini_lapack_Z_s ${LIBSRC})

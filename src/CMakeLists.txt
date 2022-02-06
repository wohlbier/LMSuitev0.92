add_subdirectory (mini_blas)
add_subdirectory (mini_lapack)
add_subdirectory (mini_slatec)

if (${CMAKE_Fortran_COMPILER_ID} STREQUAL "GNU")
  set (CLEANOPS clean_ops_nix_gfortran.f90)
elseif (${CMAKE_Fortran_COMPILER_ID} STREQUAL "Intel")
  set (CLEANOPS clean_ops_nix_ifort.f90)
else ()
  message (FATAL_ERROR "Compiler not recognized.")
endif ()

set (SRC
     ${CLEANOPS}
     initialize.f90
     latte_dc_constant.f90
     main.f90
     make_frames.f90
     mini_blas
     mini_lapack
     mini_slatec
     movie_data_class.f90
     ode_functions.f90
     parameters.f90
     precision95.f90
     process_inputs.f90
     scan_data_class.f90
     solveTWT.f90
     TWTmodules.f90
     working_variables.f90)
add_executable (lmsuite ${SRC})
target_link_libraries (lmsuite
		       mini_slatec_s
		       mini_lapack_Z_s
		       mini_lapack_util_s
		       mini_blas_s)

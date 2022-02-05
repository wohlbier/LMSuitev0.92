!**************************************************************************!
! This file is part of LATTE/MUSE Numerical Suite (lmsuite)
!
! lmsuite is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 2 of the License, or
! (at your option) any later version.
!
! lmsuite is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with lmsuite; if not, write to the Free Software
! Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
!**************************************************************************!

subroutine initialize
  use working_variables
  use parameters, only : output_parameters

  print*, 'Initializing ... '
  print*, ''

  !! call clean outputs to clean up the outputs directory
  call clean_outputs

  !create frequency array
  call create_frequency_array(.true.)
  ! if only computing fundamental and harmonics for LOTS of harmonics, use
  ! routine below rather than one above. MUCH faster.
  !call create_harmonic_array(.true.)

  !for MUSE and S-MUSE compute the pair_matrix
  call make_pair_matrix
  !compute derived quantities
  call compute_derived_qtys
  !initialize dispersion arrays
  call initialize_dispersion_arrays

  !set logicals for plot types
  if (output_parameters % circuit_power_vs_z &
       .or. output_parameters % magnitude_vs_z &
       .or. output_parameters % phase_vs_z &
       .or. output_parameters % phase_differences_vs_z &
       .or. output_parameters % disk_orbits_vs_z &
       .or. output_parameters % conserved_quantity_vs_z &
       .or. output_parameters % beam_energy_vs_z &
       .or. output_parameters % dc_beam_vel_vs_z &
       .or. output_parameters % phase_space) then 
     vs_z_plot = .true.
  end if
  if (output_parameters % power_out_vs_freq &
       .or. output_parameters % gain_vs_freq &
       .or. output_parameters % phase_vs_freq) then
     vs_freq_plot = .true.
  end if
  if (output_parameters % power_out_vs_power_in &
       .or. output_parameters % gain_vs_power_in &
       .or. output_parameters % phase_out_vs_power_in) then
     vs_powin_plot = .true.
  end if

  !set logicals for runtime's
  if (interface_parameters % echo_runtime_one &
       .or. interface_parameters % echo_runtime_two &
       .or. interface_parameters % echo_runtime_three) then
     run_one = .true.
  end if
  if (interface_parameters % echo_runtime_two &
       .or. interface_parameters % echo_runtime_three) then
     run_two = .true.
  end if
  if (interface_parameters % echo_runtime_three) then
     run_three = .true.
  end if
end subroutine initialize

subroutine destroy
  use working_variables

  print*, 'De-allocating memory ... '
  print*, ''

  deallocate(fl)

  deallocate(vph_matrix)
  deallocate(K_matrix)
  deallocate(alpha_matrix)
  deallocate(scrf_matrix)
  deallocate(nscrf_matrix)

  deallocate(R_matrix)
  deallocate(L_matrix)
  deallocate(G_matrix)
  deallocate(C_matrix)
  deallocate(pC_matrix)

end subroutine destroy

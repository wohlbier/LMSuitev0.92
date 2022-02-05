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

module LATTE_module
  use precision95
  use parameters
  use working_variables, only : M
  use working_variables, only : N_disks

  implicit none
  integer :: n_eqns !(N) number of dependent fns whose sol'n desired
  complex(dp), dimension(:), allocatable :: y !vector of dep vars

contains
  subroutine init_latte
    !! ODE values that are model specific
    n_eqns = 3 * M + 2 * N_disks !system dimension
    allocate(y(n_eqns)) ! allocate y vector

    if (interface_parameters % echo_runtime_two) then
       print*, '  initializing latte'
       print "('   size of y vector is ', I7)", n_eqns
    end if
  end subroutine init_latte

  subroutine kill_latte
    if (interface_parameters % echo_runtime_two) then
       print*, '  de-allocating latte data structures'
    end if
    ! de-allocate y vector
    deallocate(y)
  end subroutine kill_latte
  !! commonTWT.i90 contains model independent subroutines
  include 'commonTWT.i90'
end module LATTE_module

module MUSE_module
  use precision95
  use parameters
  use working_variables, only : M

  implicit none
  integer :: n_eqns !(N) number of functions whose solutions are desired
  complex(dp), dimension(:), allocatable :: y

contains
  subroutine init_muse
    n_eqns =  5*(M+1)!system dimension
    allocate(y(n_eqns)) ! allocate y vector
  end subroutine init_muse
  subroutine kill_muse
    if (interface_parameters % echo_runtime_two) then
       print*, '  de-allocating muse data structures'
    end if
    ! de-allocate y vector
    deallocate(y)
  end subroutine kill_muse
  !! commonTWT.i90 contains model independent subroutines
  include 'commonTWT.i90'
end module MUSE_module

module SMUSE_module
  use precision95
  use parameters
  use working_variables, only : M

  implicit none
  integer :: n_eqns !(N) number of functions whose solutions are desired
  complex(dp), dimension(:), allocatable :: y

contains
  subroutine init_smuse
    n_eqns =  5*(M+1)!system dimension
    allocate(y(n_eqns)) ! allocate y vector
  end subroutine init_smuse
  subroutine kill_smuse
    if (interface_parameters % echo_runtime_two) then
       print*, '  de-allocating muse data structures'
    end if
    ! de-allocate y vector
    deallocate(y)
  end subroutine kill_smuse
  !! commonTWT.i90 contains model independent subroutines
  include 'commonTWT.i90'
end module SMUSE_module

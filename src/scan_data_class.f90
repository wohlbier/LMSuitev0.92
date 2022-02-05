!***************************************************************************!
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

!*** this is the scan data class module. the type definition defines the ***!
!*** class and the subroutine is a constructor. ****************************!
module scan_data_class
  use precision95
  !** no implicit variable naming **!
  implicit none

  !! define the scan data class
  type scan_data
     integer :: scanID
     logical :: two_parameter
     real(dp) :: min
     real(dp) :: max
     integer :: num_points
     logical :: distribute_log
     integer :: int_param_1
     integer :: int_param_2
     integer :: int_param_3
     real(dp) :: real_param_1
     real(dp) :: real_param_2
     logical :: logic_param_1
  end type scan_data

contains

  !! constructor for scan data class
  subroutine init_scan_data(this, scanID, two_parameter, min, max, &
       num_points, distribute_log, int_param_1, int_param_2, int_param_3, &
       real_param_1, real_param_2, logic_param_1)
    type (scan_data), intent (out) :: this
    integer, intent (in) :: scanID
    logical, intent (in) :: two_parameter
    real(dp), intent (in) :: min
    real(dp), intent (in) :: max
    integer, intent (in) :: num_points
    logical, intent (in) :: distribute_log
    integer, intent (in) :: int_param_1
    integer, intent (in) :: int_param_2
    integer, intent (in) :: int_param_3
    real(dp), intent (in) :: real_param_1
    real(dp), intent (in) :: real_param_2
    logical, intent (in) :: logic_param_1
    
    
    this % scanID = scanID
    this % two_parameter = two_parameter
    this % min = min
    this % max = max
    this % num_points = num_points
    this % distribute_log = distribute_log
    this % int_param_1 = int_param_1
    this % int_param_2 = int_param_2
    this % int_param_3 = int_param_3
    this % real_param_1 = real_param_1
    this % real_param_2 = real_param_2
    this % logic_param_1 = logic_param_1
  end subroutine init_scan_data

end module scan_data_class

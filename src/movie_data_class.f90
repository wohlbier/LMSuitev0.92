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

!*** this is the scan data class module. the type definition defines the ***!
!*** class and the subroutine is a constructor. ****************************!
module movie_data_class
  use precision95
  !** no implicit variable naming **!
  implicit none

  !! define the movie data class
  type movie_data
     integer :: movie_type
     real(dp) :: time
     integer :: num_frames
     real(dp) :: scale
  end type movie_data

contains

  !! constructor for movie data class
  subroutine init_movie_data(this, movie_type, time, num_frames, scale)
    type (movie_data), intent (out) :: this
    integer, intent (in) :: movie_type
    real(dp), intent(in) :: time
    integer, intent (in) :: num_frames
    real(dp), intent(in) :: scale

    this % movie_type = movie_type
    this % time = time
    this % num_frames = num_frames
    this % scale = scale

  end subroutine init_movie_data

end module movie_data_class

!**************************************************************************!
!******************* LATTE/MUSE Numerical Suite v0.92 *********************!
!******************** Copyright 2002 John G. Wohlbier *********************!
!**************************************************************************!

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

program main

  !** no implicit variable naming **!
  implicit none

  !** subroutines from external files **!
  external process_inputs
  external initialize
  external destroy
  external solveTWT

!*********************************************************!
!***************  load the input data  *******************!
!*********************************************************!
  
  call process_inputs

!*********************************************************!
!******  initialize plot and dispersion arrays  **********!
!*********************************************************!

  call initialize

!*********************************************************!
!********  run the code and generate the output  *********!
!*********************************************************!

  call solveTWT

!*********************************************************!
!*****************  free up the memory  ******************!
!*********************************************************!

  call destroy

!*********************************************************!
!*********************************************************!

end program main
